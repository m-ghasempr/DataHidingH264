
/*!
 *************************************************************************************
 *
 * \file me_umhexsmp.h
 *
 * \brief
 *   Fast integer pixel and sub pixel motion estimation
 *   Improved and simplified from the original UMHexagonS algorithms
 *   See JVT-P021 for details
 *
 * \author
 *    Main contributors: (see contributors.h for copyright, address and affiliation details)
 *    - Zhibo Chen                      <chenzhibo@tsinghua.org.cn>
 *    - JianFeng Xu                     <fenax@video.mdc.tsinghua.edu.cn>
 *    - Wenfang Fu                      <fwf@video.mdc.tsinghua.edu.cn>
 *
 *    - Xiaoquan Yi                     <xyi@engr.scu.edu>
 *    - Jun Zhang                       <jzhang2@engr.scu.edu>
 *
 * \date
 *    6. Nov. 2006
 *************************************************************************************
 */

#ifndef _ME_UMHEXSMP_H_
#define _ME_UMHEXSMP_H_

#include "mbuffer.h"

struct umhex_smp_struct {
  uint16  SymmetricalCrossSearchThreshold1;
  uint16  SymmetricalCrossSearchThreshold2;
  uint16  ConvergeThreshold;
  uint16  SubPelThreshold1;
  uint16  SubPelThreshold3;

  byte  **SearchState;          //state for fractional pel search
  int  ***l0_cost;       //store SAD information needed for forward median and uplayer prediction
  int  ***l1_cost;       //store SAD information needed for backward median and uplayer prediction
  byte   *flag_intra;
  int     flag_intra_SAD;

  int     pred_SAD_uplayer;     // Up layer SAD prediction
  short   pred_MV_uplayer_X;    // Up layer MV predictor X-component
  short   pred_MV_uplayer_Y;    // Up layer MV predictor Y-component
};

typedef struct umhex_smp_struct UMHexSMPStruct;

extern void    smpUMHEX_init              (ImageParameters *p_Img);
extern int     smpUMHEX_get_mem           (ImageParameters *p_Img);
extern void    smpUMHEX_free_mem          (ImageParameters *p_Img);
extern void    smpUMHEX_decide_intrabk_SAD(Macroblock *currMB);
extern void    smpUMHEX_skip_intrabk_SAD  (Macroblock *currMB);
extern void    smpUMHEX_setup             (Macroblock *currMB, short, int, int, int, int, short ******);

extern int smpUMHEXBipredIntegerPelBlockMotionSearch (Macroblock *, int, MotionVector *, MotionVector *, MotionVector *, MotionVector *, MEBlock *, int, int, int);
extern int smpUMHEXIntegerPelBlockMotionSearch       (Macroblock *currMB, MotionVector *pred_mv, MEBlock *mv_block, int min_mcost, int lambda_factor);
extern int smpUMHEXSubPelBlockMotionSearch           (Macroblock *currMB, MotionVector *pred_mv, MEBlock *mv_block, int min_mcost, int lambda_factor);
extern int smpUMHEXFullSubPelBlockMotionSearch       (Macroblock *currMB, MotionVector *pred_mv, MEBlock *mv_block, int min_mcost, int lambda_factor);
extern int smpUMHEXSubPelBlockME                     (Macroblock *currMB, MotionVector *pred_mv, MEBlock *mv_block, int min_mcost, int *lambda);

#endif
