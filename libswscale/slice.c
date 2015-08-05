#include "swscale_internal.h"

static void free_lines(SwsSlice *s)
{
    int i;
    for (i = 0; i < 2; ++i) {
        int n = s->plane[i].available_lines;
        int j;
        for (j = 0; j < n; ++j) {
            av_freep(&s->plane[i].line[j]);
            if (s->is_ring)
               s->plane[i].line[j+n] = NULL;
        }
    }

    for (i = 0; i < 4; ++i)
        memset(s->plane[i].line, 0, sizeof(uint8_t*) * s->plane[i].available_lines * (s->is_ring ? 3 : 1));
    s->should_free_lines = 0;
}

/*
 slice lines contains extra bytes for vetorial code thus @size
 is the allocated memory size and @width is the number of pixels 
*/
static int alloc_lines(SwsSlice *s, int size, int width)
{
    int i;
    int idx[2] = {3, 2};

    s->should_free_lines = 1;
    s->width = width;

    for (i = 0; i < 2; ++i) {
        int n = s->plane[i].available_lines;
        int j;
        int ii = idx[i];

        av_assert0(n == s->plane[ii].available_lines);
        for (j = 0; j < n; ++j) {
            // chroma plane line U and V are expected to be contiguous in memory
            // by mmx vertical scaler code
            s->plane[i].line[j] = av_malloc(size * 2 + 32);
            if (!s->plane[i].line[j]) {
                free_lines(s);
                return AVERROR(ENOMEM);
            }
            s->plane[ii].line[j] = s->plane[i].line[j] + size + 16; 
            if (s->is_ring) {
               s->plane[i].line[j+n] = s->plane[i].line[j];
               s->plane[ii].line[j+n] = s->plane[ii].line[j];
            }
        }
    }

    return 0;
}

static int alloc_slice(SwsSlice *s, enum AVPixelFormat fmt, int lumLines, int chrLines, int h_sub_sample, int v_sub_sample, int ring)
{
    int i;
    int size[4] = { lumLines,
                    chrLines,
                    chrLines,
                    lumLines };

    s->h_chr_sub_sample = h_sub_sample;
    s->v_chr_sub_sample = v_sub_sample;
    s->fmt = fmt;
    s->is_ring = ring;
    s->should_free_lines = 0;

    for (i = 0; i < 4; ++i) {
        int n = size[i] * ( ring == 0 ? 1 : 3);
        s->plane[i].line = av_mallocz_array(sizeof(uint8_t*), n);
        if (!s->plane[i].line) 
            return AVERROR(ENOMEM);

        s->plane[i].tmp = ring ? s->plane[i].line + size[i] * 2 : NULL;
        s->plane[i].available_lines = size[i];
        s->plane[i].sliceY = 0;
        s->plane[i].sliceH = 0;
    }
    return 0;
}

static void free_slice(SwsSlice *s)
{
    int i;
    if (s) {
        if (s->should_free_lines)
            free_lines(s);
        for (i = 0; i < 4; ++i) {
            av_freep(&s->plane[i].line);
            s->plane[i].tmp = NULL;
        }
    }
}

int ff_rotate_slice(SwsSlice *s, int lum, int chr)
{
    int i;
    if (lum) {
        for (i = 0; i < 4; i+=3) {
            int n = s->plane[i].available_lines;
            int l = s->plane[i].sliceH;

            if (l+lum >= n * 2) {
                s->plane[i].sliceY += n;
                s->plane[i].sliceH -= n;
            }
        }
    }
    if (chr) {
        for (i = 1; i < 3; ++i) {
            int n = s->plane[i].available_lines;
            int l = s->plane[i].sliceH;

            if (l+chr >= n * 2) {
                s->plane[i].sliceY += n;
                s->plane[i].sliceH -= n;
            }
        }
    }
    return 0;
}

int ff_init_slice_from_src(SwsSlice * s, uint8_t *src[4], int stride[4], int srcW, int lumY, int lumH, int chrY, int chrH)
{
    int i = 0;

    const int start[4] = {lumY,
                    chrY,
                    chrY,
                    lumY};

    const int end[4] = {lumY +lumH,
                        chrY + chrH,
                        chrY + chrH,
                        lumY + lumH};

    s->width = srcW;

    for (i = 0; i < 4; ++i) {
        int j;
        int lines = end[i];
        lines = s->plane[i].available_lines < lines ? s->plane[i].available_lines : lines;

        if (end[i] > s->plane[i].sliceY+s->plane[i].sliceH) {
            if (start[i] <= s->plane[i].sliceY+1)
                s->plane[i].sliceY = FFMIN(start[i], s->plane[i].sliceY);
            else
                s->plane[i].sliceY = start[i];
            s->plane[i].sliceH = end[i] - s->plane[i].sliceY;
        } else {
            if (end[i] >= s->plane[i].sliceY)
                s->plane[i].sliceH = s->plane[i].sliceY + s->plane[i].sliceH - start[i];
            else
                s->plane[i].sliceH = end[i] - start[i];
            s->plane[i].sliceY = start[i];
        }

        for (j = start[i]; j < lines; j+= 1)
            s->plane[i].line[j] = src[i] + (start[i] + j) * stride[i];

    }

    return 0;
}

static int lum_h_scale(SwsContext *c, SwsFilterDescriptor *desc, int sliceY, int sliceH)
{
    FilterContext *instance = desc->instance;
    int srcW = desc->src->width;
    int dstW = desc->dst->width;
    int xInc = instance->xInc;

    int i;
    for (i = 0; i < sliceH; ++i) {
        uint8_t ** src = desc->src->plane[0].line;
        uint8_t ** dst = desc->dst->plane[0].line;
        int src_pos = sliceY+i - desc->src->plane[0].sliceY;
        int dst_pos = sliceY+i - desc->dst->plane[0].sliceY;

    
        if (c->hyscale_fast) {
            c->hyscale_fast(c, (int16_t*)dst[dst_pos], dstW, src[src_pos], srcW, xInc);
        } else {
            c->hyScale(c, (int16_t*)dst[dst_pos], dstW, (const uint8_t *)src[src_pos], instance->filter,
                       instance->filter_pos, instance->filter_size);
        }

        if (c->lumConvertRange)
            c->lumConvertRange((int16_t*)dst[dst_pos], dstW);

        desc->dst->plane[0].sliceH += 1;

        if (desc->alpha) {
            src = desc->src->plane[3].line;
            dst = desc->dst->plane[3].line;

            src_pos = sliceY+i - desc->src->plane[3].sliceY;
            dst_pos = sliceY+i - desc->dst->plane[3].sliceY;

            desc->dst->plane[3].sliceH += 1;

            if (c->hyscale_fast) {
                c->hyscale_fast(c, (int16_t*)dst[dst_pos], dstW, src[src_pos], srcW, xInc);
            } else {
                c->hyScale(c, (int16_t*)dst[dst_pos], dstW, (const uint8_t *)src[src_pos], instance->filter,
                            instance->filter_pos, instance->filter_size);
            }
        }
    }

    return sliceH;
}

static int lum_convert(SwsContext *c, SwsFilterDescriptor *desc, int sliceY, int sliceH)
{
    int srcW = desc->src->width;
    ColorContext * instance = desc->instance;
    uint32_t * pal = instance->pal;
    int i;

    desc->dst->plane[0].sliceY = sliceY;
    desc->dst->plane[0].sliceH = sliceH;
    desc->dst->plane[3].sliceY = sliceY;
    desc->dst->plane[3].sliceH = sliceH;

    for (i = 0; i < sliceH; ++i) {
        int sp0 = sliceY+i - desc->src->plane[0].sliceY;
        int sp1 = ((sliceY+i) >> desc->src->v_chr_sub_sample) - desc->src->plane[1].sliceY;
        const uint8_t * src[4] = { desc->src->plane[0].line[sp0],
                        desc->src->plane[1].line[sp1],
                        desc->src->plane[2].line[sp1],
                        desc->src->plane[3].line[sp0]};
        uint8_t * dst = desc->dst->plane[0].line[i];

        if (c->lumToYV12) {
            c->lumToYV12(dst, src[0], src[1], src[2], srcW, pal);
        } else if (c->readLumPlanar) {
            c->readLumPlanar(dst, src, srcW, c->input_rgb2yuv_table);
        } 
        
        
        if (desc->alpha) {
            dst = desc->dst->plane[3].line[i];
            if (c->alpToYV12) {
                c->alpToYV12(dst, src[3], src[1], src[2], srcW, pal);
            } else if (c->readAlpPlanar) {
                c->readAlpPlanar(dst, src, srcW, NULL);
            }
        }
    }

    return sliceH;
}

static int init_desc_fmt_convert(SwsFilterDescriptor *desc, SwsSlice * src, SwsSlice *dst, uint32_t *pal)
{
    ColorContext * li = av_malloc(sizeof(ColorContext));
    if (!li)
        return AVERROR(ENOMEM);
    li->pal = pal;
    desc->instance = li;

    desc->alpha = isALPHA(src->fmt) && isALPHA(dst->fmt);
    desc->src =src;
    desc->dst = dst;
    desc->process = &lum_convert;

    return 0;
}


static int init_desc_hscale(SwsFilterDescriptor *desc, SwsSlice *src, SwsSlice *dst, uint16_t *filter, int * filter_pos, int filter_size, int xInc)
{
    FilterContext *li = av_malloc(sizeof(FilterContext));
    if (!li)
        return AVERROR(ENOMEM);

    li->filter = filter;
    li->filter_pos = filter_pos;
    li->filter_size = filter_size;
    li->xInc = xInc;

    desc->instance = li;

    desc->alpha = isALPHA(src->fmt) && isALPHA(dst->fmt);
    desc->src = src;
    desc->dst = dst;

    desc->process = &lum_h_scale;

    return 0;
}

static int chr_h_scale(SwsContext *c, SwsFilterDescriptor *desc, int sliceY, int sliceH)
{
    FilterContext *instance = desc->instance;
    int srcW = FF_CEIL_RSHIFT(desc->src->width, desc->src->h_chr_sub_sample);
    int dstW = FF_CEIL_RSHIFT(desc->dst->width, desc->dst->h_chr_sub_sample);
    int xInc = instance->xInc;

    uint8_t ** src1 = desc->src->plane[1].line;
    uint8_t ** dst1 = desc->dst->plane[1].line;
    uint8_t ** src2 = desc->src->plane[2].line;
    uint8_t ** dst2 = desc->dst->plane[2].line;

    int src_pos1 = sliceY - desc->src->plane[1].sliceY;
    int dst_pos1 = sliceY - desc->dst->plane[1].sliceY;

    int src_pos2 = sliceY - desc->src->plane[2].sliceY;
    int dst_pos2 = sliceY - desc->dst->plane[2].sliceY;

    int i;
    for (i = 0; i < sliceH; ++i) {
        if (c->hcscale_fast) {
            c->hcscale_fast(c, (uint16_t*)dst1[dst_pos1+i], (uint16_t*)dst2[dst_pos2+i], dstW, src1[src_pos1+i], src2[src_pos2+i], srcW, xInc);
        } else {
            c->hcScale(c, (uint16_t*)dst1[dst_pos1+i], dstW, src1[src_pos1+i], instance->filter, instance->filter_pos, instance->filter_size);
            c->hcScale(c, (uint16_t*)dst2[dst_pos2+i], dstW, src2[src_pos2+i], instance->filter, instance->filter_pos, instance->filter_size);
        }

        if (c->chrConvertRange)
            c->chrConvertRange((uint16_t*)dst1[dst_pos1+i], (uint16_t*)dst2[dst_pos2+i], dstW);

        desc->dst->plane[1].sliceH += 1;
        desc->dst->plane[2].sliceH += 1;
    }
    return sliceH;
}

static int chr_convert(SwsContext *c, SwsFilterDescriptor *desc, int sliceY, int sliceH)
{
    int srcW = FF_CEIL_RSHIFT(desc->src->width, desc->src->h_chr_sub_sample);
    ColorContext * instance = desc->instance;
    uint32_t * pal = instance->pal;

    int sp0 = (sliceY - (desc->src->plane[0].sliceY >> desc->src->v_chr_sub_sample)) << desc->src->v_chr_sub_sample;
    int sp1 = sliceY - desc->src->plane[1].sliceY;

    int i;

    desc->dst->plane[1].sliceY = sliceY;
    desc->dst->plane[1].sliceH = sliceH;
    desc->dst->plane[2].sliceY = sliceY;
    desc->dst->plane[2].sliceH = sliceH;

    for (i = 0; i < sliceH; ++i) {
        const uint8_t * src[4] = { desc->src->plane[0].line[sp0+i],
                        desc->src->plane[1].line[sp1+i],
                        desc->src->plane[2].line[sp1+i],
                        desc->src->plane[3].line[sp0+i]};

        uint8_t * dst1 = desc->dst->plane[1].line[i];
        uint8_t * dst2 = desc->dst->plane[2].line[i];
        if (c->chrToYV12) {
            c->chrToYV12(dst1, dst2, src[0], src[1], src[2], srcW, pal);
        } else if (c->readChrPlanar) {
            c->readChrPlanar(dst1, dst2, src, srcW, c->input_rgb2yuv_table);
        }
    }
    return sliceH;
}

static int init_desc_cfmt_convert(SwsFilterDescriptor *desc, SwsSlice * src, SwsSlice *dst, uint32_t *pal)
{
    ColorContext * li = av_malloc(sizeof(ColorContext));
    if (!li)
        return AVERROR(ENOMEM);
    li->pal = pal;
    desc->instance = li;

    desc->src =src;
    desc->dst = dst;
    desc->process = &chr_convert;

    return 0;
}

static int init_desc_chscale(SwsFilterDescriptor *desc, SwsSlice *src, SwsSlice *dst, uint16_t *filter, int * filter_pos, int filter_size, int xInc)
{
    FilterContext *li = av_malloc(sizeof(FilterContext));
    if (!li)
        return AVERROR(ENOMEM);

    li->filter = filter;
    li->filter_pos = filter_pos;
    li->filter_size = filter_size;
    li->xInc = xInc;

    desc->instance = li;

    desc->alpha = isALPHA(src->fmt) && isALPHA(dst->fmt);
    desc->src = src;
    desc->dst = dst;

    desc->process = &chr_h_scale;

    return 0;
}

static void fill_ones(SwsSlice *s, int n, int is16bit)
{
    int i;
    for (i = 0; i < 4; ++i) {
        int j;
        int size = s->plane[i].available_lines;
        for (j = 0; j < size; ++j) {
            int k;
            int end = is16bit ? n>>1: n;
            // fill also one extra element
            end += 1;
            if (is16bit)
                for (k = 0; k < end; ++k)
                    ((int32_t*)(s->plane[i].line[j]))[k] = 1<<18;
            else
                for (k = 0; k < end; ++k)
                    ((int16_t*)(s->plane[i].line[j]))[k] = 1<<14;
        }   
    }
}

static int no_chr_scale(SwsContext *c, SwsFilterDescriptor *desc, int sliceY, int sliceH)
{
    desc->dst->plane[1].sliceY = sliceY + sliceH - desc->dst->plane[1].available_lines;
    desc->dst->plane[1].sliceH = desc->dst->plane[1].available_lines;
    desc->dst->plane[2].sliceY = sliceY + sliceH - desc->dst->plane[2].available_lines;
    desc->dst->plane[2].sliceH = desc->dst->plane[2].available_lines;
    return 0;
}

static int init_desc_no_chr(SwsFilterDescriptor *desc, SwsSlice * src, SwsSlice *dst)
{
    desc->src = src;
    desc->dst = dst;
    desc->alpha = 0;
    desc->instance = NULL;
    desc->process = &no_chr_scale;
    return 0;
}

#define FREE_FILTERS_ON_ERROR(err, ctx) if ((err) < 0) {            \
                                            ff_free_filters((ctx)); \
                                            return (err);           \
                                        }

int ff_init_filters(SwsContext * c)
{
    int i;
    int index;
    int num_ydesc;
    int num_cdesc;
    int need_lum_conv = c->lumToYV12 || c->readLumPlanar || c->alpToYV12 || c->readAlpPlanar;
    int need_chr_conv = c->chrToYV12 || c->readChrPlanar;
    int srcIdx, dstIdx;
    int dst_stride = FFALIGN(c->dstW * sizeof(int16_t) + 66, 16);

    uint32_t * pal = usePal(c->srcFormat) ? c->pal_yuv : (uint32_t*)c->input_rgb2yuv_table;
    int res = 0;

    if (c->dstBpc == 16)
        dst_stride <<= 1;

    num_ydesc = need_lum_conv ? 2 : 1;
    num_cdesc = need_chr_conv ? 2 : 1;

    c->numSlice = FFMAX(num_ydesc, num_cdesc) + 1;
    c->numDesc = num_ydesc + num_cdesc;
    c->descIndex[0] = num_ydesc;
    c->descIndex[1] = num_ydesc + num_cdesc;

    

    c->desc = av_mallocz_array(sizeof(SwsFilterDescriptor), c->numDesc);
    if (!c->desc)
        return AVERROR(ENOMEM);
    c->slice = av_mallocz_array(sizeof(SwsSlice), c->numSlice);


    res = alloc_slice(&c->slice[0], c->srcFormat, c->srcH, c->chrSrcH, c->chrSrcHSubSample, c->chrSrcVSubSample, 0);
    FREE_FILTERS_ON_ERROR(res, c);
    for (i = 1; i < c->numSlice-1; ++i) {
        res = alloc_slice(&c->slice[i], c->srcFormat, c->vLumFilterSize, c->vChrFilterSize, c->chrSrcHSubSample, c->chrSrcVSubSample, 0);
        FREE_FILTERS_ON_ERROR(res, c);
        res = alloc_lines(&c->slice[i], FFALIGN(c->srcW*2+78, 16), c->srcW);
        FREE_FILTERS_ON_ERROR(res, c);
    }
    res = alloc_slice(&c->slice[i], c->srcFormat, c->vLumFilterSize, c->vChrFilterSize, c->chrDstHSubSample, c->chrDstVSubSample, 1);
    FREE_FILTERS_ON_ERROR(res, c);
    res = alloc_lines(&c->slice[i], dst_stride, c->dstW);
    FREE_FILTERS_ON_ERROR(res, c);

    fill_ones(&c->slice[i], dst_stride>>1, c->dstBpc == 16);

    index = 0;
    srcIdx = 0;
    dstIdx = 1;

    if (need_lum_conv) {
        init_desc_fmt_convert(&c->desc[index], &c->slice[srcIdx], &c->slice[dstIdx], pal);
        c->desc[index].alpha = c->alpPixBuf != 0;
        ++index;
        srcIdx = dstIdx;
    }


    dstIdx = FFMAX(num_ydesc, num_cdesc);
    init_desc_hscale(&c->desc[index], &c->slice[index], &c->slice[dstIdx], c->hLumFilter, c->hLumFilterPos, c->hLumFilterSize, c->lumXInc);
    c->desc[index].alpha = c->alpPixBuf != 0;


    ++index;
    {
        srcIdx = 0;
        dstIdx = 1;
        if (need_chr_conv) {
            init_desc_cfmt_convert(&c->desc[index], &c->slice[srcIdx], &c->slice[dstIdx], pal);
            ++index;
            srcIdx = dstIdx;
        }

        dstIdx = FFMAX(num_ydesc, num_cdesc);
        if (c->needs_hcscale)
            init_desc_chscale(&c->desc[index], &c->slice[srcIdx], &c->slice[dstIdx], c->hChrFilter, c->hChrFilterPos, c->hChrFilterSize, c->chrXInc);
        else
            init_desc_no_chr(&c->desc[index], &c->slice[srcIdx], &c->slice[dstIdx]);
    }

    return 0;
}

int ff_free_filters(SwsContext *c)
{
    int i;
    if (c->desc) {
        for (i = 0; i < c->numDesc; ++i)
            av_freep(&c->desc[i].instance);
        av_freep(&c->desc);
    }

    if (c->slice) {
        for (i = 0; i < c->numSlice; ++i)
            free_slice(&c->slice[i]);
        av_freep(&c->slice);
    }
    return 0;
}






