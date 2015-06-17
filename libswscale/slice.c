#include "swscale_internal.h"


static int alloc_slice(SwsSlice * s, enum AVPixelFormat fmt, int lines, int v_sub_sample, int h_sub_sample)
{
    int i;
    int err = 0;

    int size[4] = { lines,
                    FF_CEIL_RSHIFT(lines, v_sub_sample),
                    FF_CEIL_RSHIFT(lines, v_sub_sample),
                    lines };

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

int ff_init_slice_from_src(SwsSlice * s, uint8_t *src[4], int stride[4], int srcW, int sliceY, int sliceH, int skip)
{
    int i = 0;

    const int start[4] = {sliceY,
                    sliceY >> s->v_chr_sub_sample,
                    sliceY >> s->v_chr_sub_sample,
                    sliceY};

    const int stride1[4] = {stride[0],
                    stride[1] << skip,
                    stride[2] << skip,
                    stride[3]};

    const int height[4] = {sliceH,
                    FF_CEIL_RSHIFT(sliceH, s->v_chr_sub_sample),
                    FF_CEIL_RSHIFT(sliceH, s->v_chr_sub_sample),
                    sliceH};

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

int ff_init_slice_from_lp(SwsSlice *s, uint8_t ***linesPool, int dstW, int sliceY, int sliceH)
{
    int i;
    const int start[4] = {sliceY,
                    sliceY >> s->v_chr_sub_sample,
                    sliceY >> s->v_chr_sub_sample,
                    sliceY};
    const int height[4] = {sliceH,
                    FF_CEIL_RSHIFT(sliceH, s->v_chr_sub_sample),
                    FF_CEIL_RSHIFT(sliceH, s->v_chr_sub_sample),
                    sliceH};
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
    uint8_t *ptr[4] = {v, v, v, v2};
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
    LumScaleInstance *instance = desc->instance;
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
    LumConvertInstance * instance = desc->instance;
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
    LumConvertInstance * li = av_malloc(sizeof(LumConvertInstance));
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
    LumScaleInstance *li = av_malloc(sizeof(LumScaleInstance));
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

int ff_init_filters(SwsContext * c)
{
    int i;
    int need_convert = c->lumToYV12 || c->readLumPlanar || c->alpToYV12 || c->readAlpPlanar;

    c->numDesc = need_convert ? 2 : 1;
    c->desc = av_malloc_array(sizeof(SwsFilterDescriptor), c->numDesc);
    c->slice = av_malloc_array(sizeof(SwsSlice), c->numDesc + 1);

    for (i = 0; i < c->numDesc+1; ++i)
        alloc_slice(&c->slice[i], c->srcFormat, c->vLumFilterSize, 0, 0);

    i = 0;
    if (need_convert)
    {
        init_desc_fmt_convert(&c->desc[i], &c->slice[i], &c->slice[i+1], (uint32_t) usePal(c->srcFormat) ? c->pal_yuv : c->input_rgb2yuv_table);
        init_slice_1(&c->slice[i+1], c->formatConvBuffer, (c->formatConvBuffer + FFALIGN(c->srcW*2+78, 16)), c->srcW, 0, c->vLumFilterSize);
        c->desc[i].alpha = c->alpPixBuf != 0;
        ++i;
    }

    
    init_desc_hscale(&c->desc[i], &c->slice[i], &c->slice[i+1], c->hLumFilter, c->hLumFilterPos, c->hLumFilterSize, c->lumXInc);
    c->desc[i].alpha = c->alpPixBuf != 0;

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
        for (i = 0; i < c->numDesc+1; ++i)
            free_slice(&c->slice[i]);
    }
    return 1;
}
