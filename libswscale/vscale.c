static int lum_planar_vscale(SwsContext *c, SwsFilterDescriptor *desc, int sliceY, int sliceH)
{
    VScalerContext *inst = desc->instance;
    int dstW = desc->dst->width;

    int sp = desc->src->plane[0].sliceH - inst->filter_size;
    int dp = sliceY - desc->dst->plane[0].sliceY;
    uint8_t **src = desc->src->plane[0].line + sp;
    uint8_t **dst = desc->dst->plane[0].line + dp;
    uint16_t *filter = inst->filter[0] + (inst->isMMX ? 0 : sliceY * inst->filter_size);

    if (inst->filter_size == 1)
        ((yuv2planar1_fn)inst->pfn)((const int16_t*)src[0], dst[0], dstW, c->lumDither8, 0);
    else
        ((yuv2planarX_fn)inst->pfn)(filter, inst->filter_size, (const int16_t**)src, dst[0], dstW, c->lumDither8, 0);

    if (desc->alpha) {
        int sp = desc->src->plane[3].sliceH - inst->filter_size;
        int dp = sliceY - desc->dst->plane[3].sliceY;
        uint8_t **src = desc->src->plane[3].line + sp;
        uint8_t **dst = desc->dst->plane[3].line + dp;
        uint16_t *filter = inst->filter[1] + (inst->isMMX ? 0 : sliceY * inst->filter_size);

        if (inst->filter_size == 1)
            ((yuv2planar1_fn)inst->pfn)((const int16_t*)src[0], dst[0], dstW, c->lumDither8, 0);
        else
            ((yuv2planarX_fn)inst->pfn)(filter, inst->filter_size, (const int16_t**)src, dst[0], dstW, c->lumDither8, 0);
    }

    return 1;
}

static int chr_planar_vscale(SwsContext *c, SwsFilterDescriptor *desc, int sliceY, int sliceH)
{
    const int chrSkipMask = (1 << desc->dst->v_chr_sub_sample) - 1;
    if (sliceY & chrSkipMask)
        return 0;
    else {
        VScalerContext *inst = desc->instance;
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
        uint16_t *filter = inst->filter[0] + (inst->isMMX ? 0 : chrSliceY * inst->filter_size);

        if (c->yuv2nv12cX) {
            ((yuv2interleavedX_fn)inst->pfn)(c, filter, inst->filter_size, (const int16_t**)src1, (const int16_t**)src2, dst1[0], dstW);
        } else if (inst->filter_size == 1) {
            ((yuv2planar1_fn)inst->pfn)((const int16_t*)src1[0], dst1[0], dstW, c->chrDither8, 0);
            ((yuv2planar1_fn)inst->pfn)((const int16_t*)src2[0], dst2[0], dstW, c->chrDither8, 3);
        } else {
            ((yuv2planarX_fn)inst->pfn)(filter, inst->filter_size, (const int16_t**)src1, dst1[0], dstW, c->chrDither8, 0);
            ((yuv2planarX_fn)inst->pfn)(filter, inst->filter_size, (const int16_t**)src2, dst2[0], dstW, c->chrDither8, inst->isMMX ? (c->uv_offx2 >> 1) : 3);
        }
    }

    return 1;
}

static int packed_vscale(SwsContext *c, SwsFilterDescriptor *desc, int sliceY, int sliceH)
{
    VScalerContext *inst = desc->instance;
    int dstW = desc->dst->width;
    int chrSliceY = sliceY >> desc->dst->v_chr_sub_sample;

    int lum_fsize = inst[0].filter_size;
    int chr_fsize = inst[1].filter_size;
    uint16_t *lum_filter = inst[0].filter[0];
    uint16_t *chr_filter = inst[1].filter[0];

    int sp0 = desc->src->plane[0].sliceH - lum_fsize;
    int sp1 = desc->src->plane[1].sliceH - chr_fsize;
    int sp2 = desc->src->plane[2].sliceH - chr_fsize;
    int sp3 = desc->src->plane[3].sliceH - lum_fsize;
    int dp = sliceY - desc->dst->plane[0].sliceY;
    uint8_t **src0 = desc->src->plane[0].line + sp0;
    uint8_t **src1 = desc->src->plane[1].line + sp1;
    uint8_t **src2 = desc->src->plane[2].line + sp2;
    uint8_t **src3 = desc->alpha ? desc->src->plane[3].line + sp3 : NULL;
    uint8_t **dst = desc->dst->plane[0].line + dp;
    

    if (c->yuv2packed1 && lum_fsize == 1 && chr_fsize <= 2) { // unscaled RGB
        int chrAlpha = chr_fsize == 1 ? 0 : chr_filter[2 * sliceY + 1];
        ((yuv2packed1_fn)inst->pfn)(c, (const int16_t*)*src0, (const int16_t**)src1, (const int16_t**)src2, (const int16_t*)(desc->alpha ? *src3 : NULL),  *dst, dstW, chrAlpha, sliceY);
    } else if (c->yuv2packed2 && lum_fsize == 2 && chr_fsize == 2) { // bilinear upscale RGB
        int lumAlpha = lum_filter[2 * sliceY + 1];
        int chrAlpha = chr_filter[2 * sliceY + 1];
        c->lumMmxFilter[2] =
        c->lumMmxFilter[3] = lum_filter[2 * sliceY]    * 0x10001;
        c->chrMmxFilter[2] =
        c->chrMmxFilter[3] = chr_filter[2 * chrSliceY] * 0x10001;
        ((yuv2packed2_fn)inst->pfn)(c, (const int16_t**)src0, (const int16_t**)src1, (const int16_t**)src2, (const int16_t**)src3,
                    *dst, dstW, lumAlpha, chrAlpha, sliceY);
    } else { // general RGB
        ((yuv2packedX_fn)inst->pfn)(c, lum_filter + sliceY * lum_fsize,
                    (const int16_t**)src0, lum_fsize, chr_filter + sliceY * chr_fsize,
                    (const int16_t**)src1, (const int16_t**)src2, chr_fsize, (const int16_t**)src3, *dst, dstW, sliceY);
    }
    return 1;
}

static int any_vscale(SwsContext *c, SwsFilterDescriptor *desc, int sliceY, int sliceH)
{
    VScalerContext *inst = desc->instance;
    int dstW = desc->dst->width;
    int chrSliceY = sliceY >> desc->dst->v_chr_sub_sample;

    int lum_fsize = inst[0].filter_size;
    int chr_fsize = inst[1].filter_size;
    uint16_t *lum_filter = inst[0].filter[0];
    uint16_t *chr_filter = inst[1].filter[0];

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
    uint8_t **src3 = desc->alpha ? desc->src->plane[3].line + sp3 : NULL;
    uint8_t *dst[4] = { desc->dst->plane[0].line[dp0],
                        desc->dst->plane[1].line[dp1],
                        desc->dst->plane[2].line[dp2],
                        desc->alpha ? desc->dst->plane[3].line[dp3] : NULL };

    av_assert1(!c->yuv2packed1 && !c->yuv2packed2);
    ((yuv2anyX_fn)inst->pfn)(c, lum_filter + sliceY * lum_fsize,
             (const int16_t**)src0, lum_fsize, chr_filter + sliceY * chr_fsize,
             (const int16_t**)src1, (const int16_t**)src2, chr_fsize, (const int16_t**)src3, dst, dstW, sliceY);

    return 1;

}

void ff_set_desc_mmx(SwsContext *c, int enabled)
{
    int idx = c->numDesc-1;
    if (isPlanarYUV(c->dstFormat) || (isGray(c->dstFormat) && !isALPHA(c->dstFormat))) {
        if (!isGray(c->dstFormat)) {
            ((VScalerContext *)c->desc[idx].instance)->isMMX = enabled;
            ((VScalerContext *)c->desc[idx].instance)->filter[0] = enabled ? c->chrMmxFilter : c->vChrFilter;
            --idx;
        }
        ((VScalerContext *)c->desc[idx].instance)->isMMX = enabled;
        ((VScalerContext *)c->desc[idx].instance)->filter[0] = enabled ? c->lumMmxFilter : c->vLumFilter;
        ((VScalerContext *)c->desc[idx].instance)->filter[1] = enabled ? c->alpMmxFilter : c->vLumFilter;
    } else {
        ((VScalerContext *)c->desc[idx].instance)[0].isMMX = enabled;
        ((VScalerContext *)c->desc[idx].instance)[1].isMMX = enabled;
    }
}

int ff_init_vscale(SwsContext *c, SwsFilterDescriptor *desc, SwsSlice *src, SwsSlice *dst);

int ff_init_vscale(SwsContext *c, SwsFilterDescriptor *desc, SwsSlice *src, SwsSlice *dst)
{
    VScalerContext *lumCtx = NULL;
    VScalerContext *chrCtx = NULL;

    if (isPlanarYUV(c->dstFormat) || (isGray(c->dstFormat) && !isALPHA(c->dstFormat))) {
        lumCtx = av_mallocz(sizeof(FilterContext));
        if (!lumCtx)
            return AVERROR(ENOMEM);
        

        desc[0].process = lum_planar_vscale;
        desc[0].instance = lumCtx;
        desc[0].src = src;
        desc[0].dst = dst;
        desc[0].alpha = c->alpPixBuf != 0;

        lumCtx->filter[0] = c->use_mmx_vfilter ? c->lumMmxFilter : c->vLumFilter;
        lumCtx->filter[1] = c->use_mmx_vfilter ? c->alpMmxFilter : c->vLumFilter;
        lumCtx->filter_size = c->vLumFilterSize;
        lumCtx->isMMX = c->use_mmx_vfilter;

        if (!isGray(c->dstFormat)) {
            chrCtx = av_mallocz(sizeof(FilterContext));
            if (!chrCtx)
                return AVERROR(ENOMEM);
            chrCtx->filter[0] = c->use_mmx_vfilter ? c->chrMmxFilter : c->vChrFilter;
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

        lumCtx->filter[0] = c->vLumFilter;
        lumCtx->filter_size = c->vLumFilterSize;
        lumCtx->isMMX = c->use_mmx_vfilter;
        chrCtx->filter[0] = c->vChrFilter;
        chrCtx->filter_size = c->vChrFilterSize;
        chrCtx->isMMX = c->use_mmx_vfilter;
    }

    ff_init_pfn(c, c->yuv2plane1, c->yuv2planeX, c->yuv2nv12cX,
        c->yuv2packed1, c->yuv2packed2, c->yuv2packedX, c->yuv2anyX);
    return 0;
}

void ff_init_pfn(SwsContext *c,
    yuv2planar1_fn yuv2plane1,
    yuv2planarX_fn yuv2planeX,
    yuv2interleavedX_fn yuv2nv12cX,
    yuv2packed1_fn yuv2packed1,
    yuv2packed2_fn yuv2packed2,
    yuv2packedX_fn yuv2packedX,
    yuv2anyX_fn yuv2anyX)
{
    VScalerContext *lumCtx = NULL;
    VScalerContext *chrCtx = NULL;
    int idx = c->numDesc - 1;

    if (isPlanarYUV(c->dstFormat) || (isGray(c->dstFormat) && !isALPHA(c->dstFormat))) {
        if (!isGray(c->dstFormat)) {
            chrCtx = c->desc[idx].instance;
            --idx;
            if (yuv2nv12cX)               chrCtx->pfn = yuv2nv12cX;
            else if (c->vChrFilterSize == 1) chrCtx->pfn = yuv2plane1;
            else                             chrCtx->pfn = yuv2planeX;
        }

        lumCtx = c->desc[idx].instance;
        if (c->vLumFilterSize == 1) lumCtx->pfn = yuv2plane1;
        else                        lumCtx->pfn = yuv2planeX;

    } else {
        lumCtx = c->desc[idx].instance;

        if (yuv2packedX) {
            if (c->yuv2packed1 && c->vLumFilterSize == 1 && c->vChrFilterSize <= 2)
                lumCtx->pfn = yuv2packed1;
            else if (c->yuv2packed2 && c->vLumFilterSize == 2 && c->vChrFilterSize == 2)
                lumCtx->pfn = yuv2packed2;
            else
                lumCtx->pfn = yuv2packedX;
        } else
            lumCtx->pfn = yuv2anyX;
    }
}


