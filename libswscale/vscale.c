/*
static int lum_p1_vscale(SwsContext *c, SwsFilterDescriptor *desc, int sliceY, int sliceH)
{
    int dstW = desc->dst->width;
    int sp = desc->src->plane[0].sliceH - 1;
    uint8_t *src = desc->src->plane[0].line[sp];
    
    int dp = sliceY - desc->dst->plane[0].sliceY;
    uint8_t *dst = desc->dst->plane[0].line[dp];
    
    c->yuv2plane1(src, dst, dstW, c->lumDither8, 0);
    if (desc->alpha) {
        int sp = desc->src->plane[3].sliceH - 1;
        uint8_t *src = desc->src->plane[3].line[sp];
        
        int dp = sliceY - desc->dst->plane[3].sliceY;
        uint8_t *dst = desc->dst->plane[3].line[dp];
        c->yuv2plane1(src, dst, dstW, c->lumDither8, 0);
    }

    return 1;
}

static int chr_p1_vscale(SwsContext *c, SwsFilterDescriptor *desc, int sliceY, int sliceH)
{
    const int chrSkipMask = (1 << desc->dst->v_chr_sub_sample) - 1;
    if (sliceY & chrSkipMask)
        return 0;
    else {
        int dstW = desc->dst->width;
        int sp1 = desc->src->plane[1].sliceH - 1;
        int sp2 = desc->src->plane[2].sliceH - 1;
        uint8_t *src1 = desc->src->plane[1].line[sp1];
        uint8_t *src2 = desc->src->plane[2].line[sp2];
        
        int dp1 = sliceY - desc->dst->plane[1].sliceY;
        int dp2 = sliceY - desc->dst->plane[2].sliceY;
        uint8_t *dst1 = desc->dst->plane[1].line[dp1];
        uint8_t *dst2 = desc->dst->plane[2].line[dp2];
        
        c->yuv2plane1(src1, dst1, dstW, c->chrDither8, 0);
        c->yuv2plane1(src2, dst2, dstW, c->chrDither8, 3);


        return 1;
    }
}
*/

static int lum_planar_vscale(SwsContext *c, SwsFilterDescriptor *desc, int sliceY, int sliceH)
{
    FilterContext *inst = desc->instance;
    int dstW = desc->dst->width;

    int sp = desc->src->plane[0].sliceH - inst->filter_size;
    int dp = sliceY - desc->dst->plane[0].sliceY;
    uint8_t **src = desc->src->plane[0].line + sp;
    uint8_t **dst = desc->dst->plane[0].line + dp;
    uint16_t *filter = inst->filter + sliceY * inst->filter_size;
    
    if (inst->filter_size == 1)
        c->yuv2plane1((const int16_t*)src[0], dst[0], dstW, c->lumDither8, 0);
    else
        c->yuv2planeX(filter, inst->filter_size, (const int16_t**)src, dst[0], dstW, c->lumDither8, 0);

    if (desc->alpha) {
        int sp = desc->src->plane[3].sliceH - inst->filter_size;
        int dp = sliceY - desc->dst->plane[3].sliceY;
        uint8_t **src = desc->src->plane[3].line + sp;
        uint8_t **dst = desc->dst->plane[3].line + dp;

        if (inst->filter_size == 1)
            c->yuv2plane1((const int16_t*)src[0], dst[0], dstW, c->lumDither8, 0);
        else
            c->yuv2planeX(filter, inst->filter_size, (const int16_t**)src, dst[0], dstW, c->lumDither8, 0);
    }

    return 1;
}

static int chr_planar_vscale(SwsContext *c, SwsFilterDescriptor *desc, int sliceY, int sliceH)
{
    const int chrSkipMask = (1 << desc->dst->v_chr_sub_sample) - 1;
    if (sliceY & chrSkipMask)
        return 0;
    else {
        FilterContext *inst = desc->instance;
        int dstW = FF_CEIL_RSHIFT(desc->dst->width, desc->dst->h_chr_sub_sample);
        int chrSliceY = sliceY >> desc->dst->v_chr_sub_sample;

        int sp1 = desc->src->plane[1].sliceH - inst->filter_size;
        int sp2 = desc->src->plane[2].sliceH - inst->filter_size;
        int dp1 = chrSliceY - desc->dst->plane[1].sliceY;
        int dp2 = chrSliceY - desc->dst->plane[2].sliceY;
        uint8_t **src1 = desc->src->plane[1].line + sp1;
        uint8_t **src2 = desc->src->plane[2].line + sp2;
        uint8_t **dst1 = desc->dst->plane[1].line + dp1;
        uint8_t **dst2 = desc->dst->plane[2].line + dp2;
        uint16_t *filter = inst->filter + sliceY * inst->filter_size;

        if (c->yuv2nv12cX) {
            c->yuv2nv12cX(c, filter, inst->filter_size, (const int16_t**)src1, (const int16_t**)src2, dst1[0], dstW);
        } else if (inst->filter_size == 1) {
            c->yuv2plane1((const int16_t*)src1[0], dst1[0], dstW, c->chrDither8, 0);
            c->yuv2plane1((const int16_t*)src2[0], dst2[0], dstW, c->chrDither8, 3);
        } else {
            c->yuv2planeX(filter, inst->filter_size, (const int16_t**)src1, dst1[0], dstW, c->chrDither8, 0);
            c->yuv2planeX(filter, inst->filter_size, (const int16_t**)src2, dst2[0], dstW, c->chrDither8, inst->isMMX ? (c->uv_offx2 >> 1) : 3);
        }
    }

    return 1;
}

static int packed_vscale(SwsContext *c, SwsFilterDescriptor *desc, int sliceY, int sliceH)
{
    FilterContext *inst = desc->instance;
    int dstW = desc->dst->width;
    int chrSliceY = sliceY >> desc->dst->v_chr_sub_sample;

    int lum_fsize = inst[0].filter_size;
    int chr_fsize = inst[1].filter_size;
    uint16_t *lum_filter = inst[0].filter;
    uint16_t *chr_filter = inst[1].filter;

    int sp0 = desc->src->plane[0].sliceH - lum_fsize;
    int sp1 = desc->src->plane[1].sliceH - chr_fsize;
    int sp2 = desc->src->plane[2].sliceH - chr_fsize;
    int sp3 = desc->src->plane[3].sliceH - lum_fsize;
    int dp = sliceY - desc->dst->plane[0].sliceY;
    uint8_t **src0 = desc->src->plane[0].line + sp0;
    uint8_t **src1 = desc->src->plane[1].line + sp1;
    uint8_t **src2 = desc->src->plane[2].line + sp2;
    uint8_t **src3 = desc->src->plane[3].line + sp3;
    uint8_t **dst = desc->dst->plane[0].line + dp;
    

    if (c->yuv2packed1 && lum_fsize == 1 && chr_fsize <= 2) { // unscaled RGB
        int chrAlpha = chr_fsize == 1 ? 0 : chr_filter[2 * sliceY + 1];
        c->yuv2packed1(c, (const int16_t*)*src0, (const int16_t**)src1, (const int16_t**)src2, (const int16_t*)(desc->alpha ? *src3 : NULL),  *dst, dstW, chrAlpha, sliceY);
    } else if (c->yuv2packed2 && lum_fsize == 2 && chr_fsize == 2) { // bilinear upscale RGB
        int lumAlpha = lum_filter[2 * sliceY + 1];
        int chrAlpha = chr_filter[2 * sliceY + 1];
        c->lumMmxFilter[2] =
        c->lumMmxFilter[3] = lum_filter[2 * sliceY]    * 0x10001;
        c->chrMmxFilter[2] =
        c->chrMmxFilter[3] = chr_filter[2 * chrSliceY] * 0x10001;
        c->yuv2packed2(c, (const int16_t**)src0, (const int16_t**)src1, (const int16_t**)src2, (const int16_t**)(desc->alpha ? src3 : NULL),
                    *dst, dstW, lumAlpha, chrAlpha, sliceY);
    } else { // general RGB
        c->yuv2packedX(c, lum_filter + sliceY * lum_fsize,
                    (const int16_t**)src0, lum_fsize, chr_filter + sliceY * chr_fsize,
                    (const int16_t**)src1, (const int16_t**)src2, chr_fsize, (const int16_t**)src3, *dst, dstW, sliceY);
    }
    return 1;
}

static int any_vscale(SwsContext *c, SwsFilterDescriptor *desc, int sliceY, int sliceH)
{
    FilterContext *inst = desc->instance;
    int dstW = desc->dst->width;
    int chrSliceY = sliceY >> desc->dst->v_chr_sub_sample;

    int lum_fsize = inst[0].filter_size;
    int chr_fsize = inst[1].filter_size;
    uint16_t *lum_filter = inst[0].filter;
    uint16_t *chr_filter = inst[1].filter;

    int sp0 = desc->src->plane[0].sliceH - lum_fsize;
    int sp1 = desc->src->plane[1].sliceH - chr_fsize;
    int sp2 = desc->src->plane[2].sliceH - chr_fsize;
    int sp3 = desc->src->plane[3].sliceH - lum_fsize;
    int dp0 = sliceY - desc->dst->plane[0].sliceY;
    int dp1 = chrSliceY - desc->dst->plane[1].sliceY;
    int dp2 = chrSliceY - desc->dst->plane[2].sliceY;
    int dp3 = sliceY - desc->dst->plane[3].sliceY;

    uint8_t **src0 = desc->src->plane[0].line + sp0;
    uint8_t **src1 = desc->src->plane[1].line + sp1;
    uint8_t **src2 = desc->src->plane[2].line + sp2;
    uint8_t **src3 = desc->src->plane[3].line + sp3;
    uint8_t *dst[4] = { desc->dst->plane[0].line[dp0],
                        desc->dst->plane[1].line[dp1],
                        desc->dst->plane[2].line[dp2],
                        desc->alpha ? desc->dst->plane[3].line[dp3] : NULL };


    av_assert1(!c->yuv2packed1 && !c->yuv2packed2);
    c->yuv2anyX(c, lum_filter + sliceY * lum_fsize,
             (const int16_t**)src0, lum_fsize, chr_filter + sliceY * chr_fsize,
             (const int16_t**)src1, (const int16_t**)src2, chr_fsize, (const int16_t**)src3, dst, dstW, sliceY);

    return 1;

}

void ff_set_desc_mmx(SwsContext *c, int enabled)
{
    if (isPlanarYUV(c->dstFormat) || (isGray(c->dstFormat) && !isALPHA(c->dstFormat))) {
        ((FilterContext *)c->desc[c->numDesc-1].instance)->isMMX = enabled;
        if (!isGray(c->dstFormat))
            ((FilterContext *)c->desc[c->numDesc-2].instance)->isMMX = enabled;
    } else {
        ((FilterContext *)c->desc[c->numDesc-1].instance)[0].isMMX = enabled;
        ((FilterContext *)c->desc[c->numDesc-1].instance)[1].isMMX = enabled;
    }
}

int ff_init_vscale(SwsContext *c, SwsFilterDescriptor *desc, SwsSlice *src, SwsSlice *dst);

int ff_init_vscale(SwsContext *c, SwsFilterDescriptor *desc, SwsSlice *src, SwsSlice *dst)
{
    FilterContext *lumCtx = NULL;
    FilterContext *chrCtx = NULL;

    if (isPlanarYUV(c->dstFormat) || (isGray(c->dstFormat) && !isALPHA(c->dstFormat))) {
        lumCtx = av_mallocz(sizeof(FilterContext));
        if (!lumCtx)
            return AVERROR(ENOMEM);
        

        desc[0].process = lum_planar_vscale;
        desc[0].instance = lumCtx;
        desc[0].src = src;
        desc[0].dst = dst;
        desc[0].alpha = c->alpPixBuf != 0;

        lumCtx->filter = c->vLumFilter;
        lumCtx->filter_size = c->vLumFilterSize;
        lumCtx->isMMX = c->use_mmx_vfilter;

        if (!isGray(c->dstFormat)) {
            chrCtx = av_mallocz(sizeof(FilterContext));
            if (!chrCtx)
                return AVERROR(ENOMEM);
            chrCtx->filter = c->vChrFilter;
            chrCtx->filter_size = c->vChrFilterSize;
            chrCtx->isMMX = c->use_mmx_vfilter;
            desc[1].process = chr_planar_vscale;
            desc[1].instance = chrCtx;
            desc[1].src = src;
            desc[1].dst = dst;
        }
    } else {
        lumCtx = av_mallocz_array(sizeof(FilterContext), 2);
        if (!lumCtx)
            return AVERROR(ENOMEM); 
        chrCtx = &lumCtx[1];

        desc[0].process = c->yuv2packedX ? packed_vscale : any_vscale;
        desc[0].instance = lumCtx;
        desc[0].src = src;
        desc[0].dst = dst;
        desc[0].alpha = c->alpPixBuf != 0;

        lumCtx->filter = c->vLumFilter;
        lumCtx->filter_size = c->vLumFilterSize;
        lumCtx->isMMX = c->use_mmx_vfilter;
        chrCtx->filter = c->vChrFilter;
        chrCtx->filter_size = c->vChrFilterSize;
        chrCtx->isMMX = c->use_mmx_vfilter;
    }
    return 0;
}


