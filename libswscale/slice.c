#include "swscale_internal.h"


static int alloc_slice(SwsSlice * s, enum AVPixelFormat fmt, int lumLines, int chrLines, int h_sub_sample, int v_sub_sample)
{
    int i;
    int err = 0;
 
    int size[4] = { lumLines,
                    chrLines,
                    chrLines,
                    lumLines };

    //s->width;
    s->h_chr_sub_sample = h_sub_sample;
    s->v_chr_sub_sample = v_sub_sample;
    s->fmt = fmt;

    for (i = 0; i < 4; ++i)
    {
        s->plane[i].line = av_malloc_array(sizeof(uint8_t*), size[i]);
        if (!s->plane[i].line) 
        {
            err = AVERROR(ENOMEM);
            break;
        }
        s->plane[i].available_lines = size[i];
        s->plane[i].sliceY = 0;
        s->plane[i].sliceH = 0;
    }

    if (err)
    {
        for (--i; i >= 0; --i)
            av_freep(&s->plane[i].line);
        return err;
    }
    return 1;
}

static void free_slice(SwsSlice *s)
{
    int i;
    for (i = 0; i < 4; ++i)
        av_freep(&s->plane[i].line);
}

int ff_init_slice_from_src(SwsSlice * s, uint8_t *src[4], int stride[4], int srcW, int lumY, int lumH, int chrY, int chrH, int skip)
{
    int i = 0;

    const int start[4] = {lumY,
                    chrY,
                    chrY,
                    lumY};
 
    const int stride1[4] = {stride[0],
                     stride[1] << skip,
                     stride[2] << skip,
                     stride[3]};
 
    const int height[4] = {lumH,
                    chrH,
                    chrH,
                    lumH};

    s->width = srcW;

    for (i = 0; i < 4; ++i)
    {
        int j;
        int lines = height[i];
        lines = s->plane[i].available_lines < lines ? s->plane[i].available_lines : lines;

        s->plane[i].sliceY = start[i];
        s->plane[i].sliceH = lines;

        for (j = 0; j < lines; j+= 1 << skip)
            s->plane[i].line[j] = src[i] + (start[i] + j) * stride1[i];

    }

    return 1;
}

int ff_init_slice_from_lp(SwsSlice *s, uint8_t ***linesPool, int dstW, int lumY, int lumH, int chrY, int chrH)
{
    int i;
    const int start[4] = {lumY,
                    chrY,
                    chrY,
                    lumY};

    const int height[4] = {lumH,
                    chrH,
                    chrH,
                    lumH};
    s->width = dstW;
    for (i = 0; i < 4; ++i)
    {
        int j;
        int lines = height[i];
        lines = s->plane[i].available_lines < lines ? s->plane[i].available_lines : lines;

        s->plane[i].sliceY = start[i];
        s->plane[i].sliceH = lines;

        for (j = 0; j < lines; ++j)
        {
            uint8_t * v = linesPool[i] ? linesPool[i][j] : NULL;
            s->plane[i].line[j] = v;
        }

    }
    return 1;
}

static int init_slice_1(SwsSlice *s, uint8_t *v, uint8_t *v2, int dstW, int sliceY, int sliceH)
{
    int i;
    uint8_t *ptr[4] = {v, v, v2, v2};
    s->width = dstW;
    for (i = 0; i < 4; ++i)
    {
        int j;
        int lines = s->plane[i].available_lines;

        s->plane[i].sliceY = sliceY;
        s->plane[i].sliceH = lines;

        for (j = 0; j < lines; ++j)
            s->plane[i].line[j] = ptr[i];

    }
    return 1;
}


static int lum_h_scale(SwsContext *c, SwsFilterDescriptor *desc, int sliceY, int sliceH)
{
    ScaleInstance *instance = desc->instance;
    int srcW = desc->src->width;
    int dstW = desc->dst->width;
    int xInc = instance->xInc;

    uint8_t ** src = desc->src->plane[0].line;
    uint8_t ** dst = desc->dst->plane[0].line;

    int src_pos = sliceY - desc->src->plane[0].sliceY;
    int dst_pos = sliceY - desc->dst->plane[0].sliceY;



    if (!c->hyscale_fast) {
        c->hyScale(c, (int16_t*)dst[dst_pos], dstW, (const uint8_t *)src[src_pos], instance->filter,
                   instance->filter_pos, instance->filter_size);
    } else { // fast bilinear upscale / crap downscale
        c->hyscale_fast(c, (int16_t*)dst[dst_pos], dstW, src[src_pos], srcW, xInc);
    }

    if (c->lumConvertRange)
        c->lumConvertRange((int16_t*)dst[dst_pos], dstW);


    if (desc->alpha)
    {
        src = desc->src->plane[3].line;
        dst = desc->dst->plane[3].line;

        src_pos = sliceY - desc->src->plane[3].sliceY;
        dst_pos = sliceY - desc->dst->plane[3].sliceY;



        if (!c->hyscale_fast) {
            c->hyScale(c, (int16_t*)dst[dst_pos], dstW, (const uint8_t *)src[src_pos], instance->filter,
                        instance->filter_pos, instance->filter_size);
        } else { // fast bilinear upscale / crap downscale
            c->hyscale_fast(c, (int16_t*)dst[dst_pos], dstW, src[src_pos], srcW, xInc);
        }
    }



    return 1;
}

static int lum_convert(SwsContext *c, SwsFilterDescriptor *desc, int sliceY, int sliceH)
{
    int srcW = desc->src->width;
    ConvertInstance * instance = desc->instance;
    uint32_t * pal = instance->pal;

    int sp = sliceY - desc->src->plane[0].sliceY;
    int dp = sliceY - desc->dst->plane[0].sliceY;

    const uint8_t * src[4] = { desc->src->plane[0].line[sp],
                        desc->src->plane[1].line[sp],
                        desc->src->plane[2].line[sp],
                        desc->src->plane[3].line[sp]};
    uint8_t * dst = desc->dst->plane[0].line[0/*dp*/];

    desc->dst->plane[0].sliceY = sliceY;
    desc->dst->plane[0].sliceH = sliceH;
    desc->dst->plane[3].sliceY = sliceY;
    desc->dst->plane[3].sliceH = sliceH;

    if (c->lumToYV12) {
        c->lumToYV12(dst, src[0], src[1], src[2], srcW, pal);
    } else if (c->readLumPlanar) {
        c->readLumPlanar(dst, src, srcW, c->input_rgb2yuv_table);
    } 
    
    
    if (desc->alpha)
    {
        dp = sliceY - desc->dst->plane[3].sliceY;
        dst = desc->dst->plane[3].line[dp];
        if (c->alpToYV12) {
            c->alpToYV12(dst, src[3], src[1], src[2], srcW, pal);
        } else if (c->readAlpPlanar) {
            c->readAlpPlanar(dst, src, srcW, NULL);
        }
    }

    return 1;
}

static int init_desc_fmt_convert(SwsFilterDescriptor *desc, SwsSlice * src, SwsSlice *dst, uint32_t *pal)
{
    ConvertInstance * li = av_malloc(sizeof(ConvertInstance));
    if (!li)
        return AVERROR(ENOMEM);
    li->pal = pal;
    desc->instance = li;

    desc->alpha = isALPHA(src->fmt) && isALPHA(dst->fmt);
    desc->src =src;
    desc->dst = dst;
    desc->process = &lum_convert;

    return 1;
}


static int init_desc_hscale(SwsFilterDescriptor *desc, SwsSlice *src, SwsSlice *dst, uint16_t *filter, int * filter_pos, int filter_size, int xInc)
{
    ScaleInstance *li = av_malloc(sizeof(ScaleInstance));
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

    return 1;
}

static int chr_h_scale(SwsContext *c, SwsFilterDescriptor *desc, int sliceY, int sliceH)
{
    ScaleInstance *instance = desc->instance;
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



    if (!c->hcscale_fast) {
        c->hcScale(c, (uint16_t*)dst1[dst_pos1], dstW, src1[src_pos1], instance->filter, instance->filter_pos, instance->filter_size);
        c->hcScale(c, (uint16_t*)dst2[dst_pos2], dstW, src2[src_pos2], instance->filter, instance->filter_pos, instance->filter_size);
    } else { // fast bilinear upscale / crap downscale
        c->hcscale_fast(c, (uint16_t*)dst1[dst_pos1], (uint16_t*)dst2[dst_pos2], dstW, src1[src_pos1], src2[src_pos2], srcW, xInc);
    }

    if (c->chrConvertRange)
        c->chrConvertRange((uint16_t*)dst1[dst_pos1], (uint16_t*)dst2[dst_pos2], dstW);

    return 1;
}

static int chr_convert(SwsContext *c, SwsFilterDescriptor *desc, int sliceY, int sliceH)
{
    int srcW = FF_CEIL_RSHIFT(desc->src->width, desc->src->h_chr_sub_sample);
    ConvertInstance * instance = desc->instance;
    uint32_t * pal = instance->pal;

    int sp = sliceY - desc->src->plane[1].sliceY;
    int dp = sliceY - desc->dst->plane[1].sliceY;

    const uint8_t * src[4] = { desc->src->plane[0].line[sp],
                        desc->src->plane[1].line[sp],
                        desc->src->plane[2].line[sp],
                        desc->src->plane[3].line[sp]};
    uint8_t * dst1 = desc->dst->plane[1].line[0/*dp*/];
    uint8_t * dst2 = desc->dst->plane[2].line[0/*dp*/];

    desc->dst->plane[1].sliceY = sliceY;
    desc->dst->plane[1].sliceH = sliceH;
    desc->dst->plane[2].sliceY = sliceY;
    desc->dst->plane[2].sliceH = sliceH;

    if (c->chrToYV12) {
        c->chrToYV12(dst1, dst2, src[0], src[1], src[2], srcW, pal);
    } else if (c->readChrPlanar) {
        c->readChrPlanar(dst1, dst2, src, srcW, c->input_rgb2yuv_table);
    }

    return 1;
}

static int init_desc_cfmt_convert(SwsFilterDescriptor *desc, SwsSlice * src, SwsSlice *dst, uint32_t *pal)
{
    ConvertInstance * li = av_malloc(sizeof(ConvertInstance));
    if (!li)
        return AVERROR(ENOMEM);
    li->pal = pal;
    desc->instance = li;

    desc->src =src;
    desc->dst = dst;
    desc->process = &chr_convert;

    return 1;
}

static int init_desc_chscale(SwsFilterDescriptor *desc, SwsSlice *src, SwsSlice *dst, uint16_t *filter, int * filter_pos, int filter_size, int xInc)
{
    ScaleInstance *li = av_malloc(sizeof(ScaleInstance));
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

    return 1;
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

    uint32_t * pal = usePal(c->srcFormat) ? c->pal_yuv : (uint32_t*)c->input_rgb2yuv_table;

    num_ydesc = need_lum_conv ? 2 : 1;
    num_cdesc = c->needs_hcscale ? (need_chr_conv ? 2 : 1) : 0;

    c->numSlice = FFMAX(num_ydesc, num_cdesc) + 1;
    c->numDesc = num_ydesc + num_cdesc;
    c->descIndex[0] = num_ydesc;
    c->descIndex[1] = num_ydesc + num_cdesc;

    

    c->desc = av_malloc_array(sizeof(SwsFilterDescriptor), c->numDesc);
    c->slice = av_malloc_array(sizeof(SwsSlice), c->numSlice);

    for (i = 0; i < c->numSlice-1; ++i)
        alloc_slice(&c->slice[i], c->srcFormat, c->vLumFilterSize, c->vChrFilterSize, c->chrSrcHSubSample, c->chrSrcVSubSample);
    alloc_slice(&c->slice[i], c->srcFormat, c->vLumFilterSize, c->vChrFilterSize, c->chrDstHSubSample, c->chrDstVSubSample);

    index = 0;
    srcIdx = 0;
    dstIdx = 1;

    // temp slice for color space conversion
    if (need_lum_conv || need_chr_conv)
        init_slice_1(&c->slice[dstIdx], c->formatConvBuffer, (c->formatConvBuffer + FFALIGN(c->srcW*2+78, 16)), c->srcW, 0, c->vLumFilterSize);

    if (need_lum_conv)
    {
        init_desc_fmt_convert(&c->desc[index], &c->slice[srcIdx], &c->slice[dstIdx], pal);
        c->desc[index].alpha = c->alpPixBuf != 0;
        ++index;
        srcIdx = dstIdx;
    }


    dstIdx = FFMAX(num_ydesc, num_cdesc);
    init_desc_hscale(&c->desc[index], &c->slice[index], &c->slice[dstIdx], c->hLumFilter, c->hLumFilterPos, c->hLumFilterSize, c->lumXInc);
    c->desc[index].alpha = c->alpPixBuf != 0;


    ++index;
    if (c->needs_hcscale)
    {
        srcIdx = 0;
        dstIdx = 1;
        if (need_chr_conv)
        {
            init_desc_cfmt_convert(&c->desc[index], &c->slice[srcIdx], &c->slice[dstIdx], pal);
            ++index;
            srcIdx = dstIdx;
        }

        dstIdx = FFMAX(num_ydesc, num_cdesc);
        init_desc_chscale(&c->desc[index], &c->slice[srcIdx], &c->slice[dstIdx], c->hChrFilter, c->hChrFilterPos, c->hChrFilterSize, c->chrXInc);
    }

    return 1;
}

int ff_free_filters(SwsContext *c)
{
    int i;
    for (i = 0; i < c->numDesc; ++i)
        av_freep(&c->desc->instance);

    av_freep(&c->desc);
    if (c->slice)
    {
        int i;
        for (i = 0; i < c->numSlice; ++i)
            free_slice(&c->slice[i]);
    }
    return 1;
}
