
/*!
 ************************************************************************
 * \file mv_search.h
 *
 * \brief
 *   array definition for motion search
 *
 * \author
 *    Inge Lille-Langoy               <inge.lille-langoy@telenor.com>   \n
 *    Alexis Michael Tourapis         <alexis.tourapis@dolby.com>       \n
 *
 ************************************************************************
 */

#ifndef _MV_SEARCH_H_
#define _MV_SEARCH_H_


extern void get_neighbors(Macroblock *currMB, PixelPos *block, int mb_x, int mb_y, int blockshape_x);
extern void set_access_method(int *access_method, MotionVector *blk, int min_x, int min_y, int max_x, int max_y);

extern void PrepareMEParams      (Slice *currSlice, MEBlock *mv_block, int ChromaMEEnable, int list, int ref);
extern void PrepareBiPredMEParams(Slice *currSlice, MEBlock *mv_block, int ChromaMEEnable, int list, int list_offset, int ref);

extern void Init_Motion_Search_Module  (ImageParameters *p_Img, InputParameters *p_Inp);
extern void Clear_Motion_Search_Module (ImageParameters *p_Img, InputParameters *p_Inp);

extern void  PartitionMotionSearch    (Macroblock *currMB, int, int, int*);
extern void  SubPartitionMotionSearch (Macroblock *currMB, int, int, int*);

extern void  Get_Direct_MV_Spatial_MBAFF  (Macroblock *currMB);
extern void  Get_Direct_MV_Spatial_Normal (Macroblock *currMB);
extern void  Get_Direct_MV_Temporal       (Macroblock *currMB);

extern void  FindSkipModeMotionVector     (Macroblock *currMB);

extern void init_ME_engine    (Macroblock *currMB);
extern int  BlockMotionSearch (Macroblock *currMB, MEBlock *mv_block, int,int, int*);
extern void init_mv_block     (Macroblock *currMB, MEBlock *mv_block, short blocktype, int list, char ref_idx, short mb_x, short mb_y);
extern void get_original_block(ImageParameters *p_Img, InputParameters *p_Inp, MEBlock *mv_block);
extern void free_mv_block     (InputParameters *p_Inp, MEBlock *mv_block);
extern void update_mv_block   (Macroblock *currMB, MEBlock *mv_block, int h, int v);
extern void get_search_range(MEBlock *mv_block, InputParameters *p_Inp, short ref, int blocktype);

static inline void add_mvs(MotionVector *mv0, const MotionVector *mv1)
{
  mv0->mv_x = (short) (mv0->mv_x + mv1->mv_x);
  mv0->mv_y = (short) (mv0->mv_y + mv1->mv_y);
}

static inline MotionVector add_MVs(MotionVector mv0, const MotionVector *mv1)
{
  mv0.mv_x = (short) (mv0.mv_x + mv1->mv_x);
  mv0.mv_y = (short) (mv0.mv_y + mv1->mv_y);
  
  return (mv0);
}

static inline MotionVector pad_MVs(MotionVector mv0, MEBlock *mv_block)
{
  mv0.mv_x = (short) (mv0.mv_x + mv_block->pos_x_padded);
  mv0.mv_y = (short) (mv0.mv_y + mv_block->pos_y_padded);
  
  return (mv0);
}

static inline int weight_cost(int lambda, int bits)
{  
#if (USE_RND_COST)
  return (rshift_rnd_sf((lambda) * (bits), LAMBDA_ACCURACY_BITS));
#else
  return (((lambda) * (bits)) >> LAMBDA_ACCURACY_BITS);
#endif
}

static inline int mv_cost(const ImageParameters *p_Img, int lambda, const MotionVector *mv, const MotionVector *pmv)
{
#if (USE_RND_COST)
  return (rshift_rnd_sf((lambda *(p_Img->mvbits[mv->mv_x - pmv->mv_x] + p_Img->mvbits[mv->mv_y - pmv->mv_y])), LAMBDA_ACCURACY_BITS));
#else
  return ((lambda *(p_Img->mvbits[mv->mv_x - pmv->mv_x] + p_Img->mvbits[mv->mv_y - pmv->mv_y]))>> LAMBDA_ACCURACY_BITS);
#endif
}

static inline int ref_cost(const ImageParameters *p_Img, int lambda, short ref, int list_offset)
{
  if (p_Img->listXsize[list_offset] <= 1)
    return 0;
  else
  {
#if (USE_RND_COST)    
    return (rshift_rnd_sf((lambda) * (p_Img->refbits[(ref)]), LAMBDA_ACCURACY_BITS));
#else
    return ((lambda *(p_Img->refbits[(ref)]))>> LAMBDA_ACCURACY_BITS);
#endif
  }
}

#endif

