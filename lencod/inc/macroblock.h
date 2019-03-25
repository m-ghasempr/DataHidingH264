
/*!
 ************************************************************************
 * \file
 *    macroblock.h
 *
 * \brief
 *    Arrays for macroblock processing
 *
 * \author
 *    Inge Lille-Langoy               <inge.lille-langoy@telenor.com>     \n
 *    Telenor Satellite Services                                          \n
 *    P.O.Box 6914 St.Olavs plass                                         \n
 *    N-0130 Oslo, Norway
 *
 ************************************************************************/

#ifndef _MACROBLOCK_H_
#define _MACROBLOCK_H_

#include "block.h"
#include "macroblock_p444.h"


static const short part_size[8][2] = 
{
  {4, 4},
  {4, 4},
  {4, 2},
  {2, 4},
  {2, 2},
  {2, 1},
  {1, 2},
  {1, 1}
};

static const short block_size[8][2] = 
{
  {16, 16},
  {16, 16},
  {16, 8},
  {8, 16},
  {8, 8},
  {8, 4},
  {4, 8},
  {4, 4}
};

extern void  next_macroblock  (Macroblock* currMB);
extern void  start_macroblock (Slice *currSlice, Macroblock** currMB, int mb_addr, Boolean mb_field);
extern void  reset_macroblock (Macroblock *currMB);
extern void  end_macroblock   (Macroblock* currMB, Boolean *end_of_slice, Boolean *recode_macroblock);
extern void  write_macroblock (Macroblock* currMB, int eos_bit);


extern int  luma_residual_coding_16x16 (Macroblock* currMB, int, int[2]);
extern int  luma_residual_coding_8x8   (Macroblock* currMB, int*, int64*, int, short, int[2], char *);
extern void luma_residual_coding       (Macroblock *currMB);
extern void chroma_residual_coding     (Macroblock *currMB);
extern void intra_chroma_RD_decision   (Macroblock *currMB, RD_PARAMS);
extern byte TransformDecision     (Macroblock *currMB, int, distblk*);
extern int  B8Mode2Value          (Slice *currSlice, short b8mode, short b8pdir);

extern int  writeMBLayerPSlice    (Macroblock *currMB, int rdopt, int *coeff_rate);
extern int  writeMBLayerBSlice    (Macroblock *currMB, int rdopt, int *coeff_rate);
extern int  writeMBLayerISlice    (Macroblock *currMB, int rdopt, int *coeff_rate);

extern void write_terminating_bit (Slice* currSlice, short bit);

extern int writeMotionInfo2NAL (Macroblock* currMB);
extern int  write_pslice_motion_info_to_NAL (Macroblock* currMB);
extern int  write_bslice_motion_info_to_NAL (Macroblock* currMB);
extern int  write_pslice_motion_info_to_NAL (Macroblock* currMB);
extern int  writeReferenceFrame   (Macroblock *currMB, int i, int j, int list_idx, int  ref);
extern int  writeMotionVector8x8  (Macroblock *currMB, int  i0, int  j0, int  i1, int  j1, int  refframe, int  list_idx, int  mv_mode, short bipred_me);

extern int  writeCoeff4x4_CABAC   (Macroblock *currMB, ColorPlane, int, int, int);
extern int  writeCoeff8x8_CABAC   (Macroblock* currMB, ColorPlane, int, int);
extern int  writeCoeff8x8         (Macroblock* currMB, ColorPlane, int, int, int);
extern int  writeCoeff16x16_CABAC (Macroblock* currMB, ColorPlane);
extern int  writeCoeff16x16_CAVLC (Macroblock* currMB, ColorPlane);

extern int  writeCoeff4x4_CAVLC_normal (Macroblock* currMB, int block_type, int b8, int b4, int param);
extern int  writeCoeff4x4_CAVLC_444    (Macroblock* currMB, int block_type, int b8, int b4, int param);

extern distblk distI16x16_sad  (Macroblock *currMB, imgpel **img_org, imgpel **pred_img, distblk min_cost);
extern distblk distI16x16_sse  (Macroblock *currMB, imgpel **img_org, imgpel **pred_img, distblk min_cost);
extern distblk distI16x16_satd (Macroblock *currMB, imgpel **img_org, imgpel **pred_img, distblk min_cost);
extern distblk find_sad_16x16_JM (Macroblock *currMB);

extern int   predict_nnz       (Macroblock *currMB, int block_type, int i,int j);
extern int   predict_nnz_chroma(Macroblock *currMB, int i,int j);

extern void SetLagrangianMultipliersOn (VideoParameters *p_Vid, InputParameters *p_Inp);
extern void SetLagrangianMultipliersOff(VideoParameters *p_Vid, InputParameters *p_Inp);

extern void set_modes_and_reframe (Macroblock* currMB, int b8, short* p_dir, int list_mode[2], char *list_ref_idx);
extern void set_modes_and_reframe_b_slice (Macroblock* currMB, int b8, short* p_dir, int list_mode[2], char *list_ref_idx);
extern void set_modes_and_reframe_p_slice (Macroblock* currMB, int b8, short* p_dir, int list_mode[2], char *list_ref_idx);


#endif

