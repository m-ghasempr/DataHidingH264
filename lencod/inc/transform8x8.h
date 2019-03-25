
/*!
 ***************************************************************************
 *
 * \file transform8x8.h
 *
 * \brief
 *    prototypes of 8x8 transform functions
 *
 * \date
 *    9. October 2003
 *
 * \author
 *    Main contributors (see contributors.h for copyright, address and affiliation details)
 *    Yuri Vatis
 **************************************************************************/

#ifndef _TRANSFORM8X8_H_
#define _TRANSFORM8X8_H_
extern int    Mode_Decision_for_Intra8x8Macroblock       (Macroblock *currMB, int lambda, distblk *min_cost);
extern int    Mode_Decision_for_8x8IntraBlocks_JM_Low    (Macroblock *currMB, int b8, int lambda, distblk *min_cost);
extern int    Mode_Decision_for_8x8IntraBlocks_JM_High   (Macroblock *currMB, int b8, int lambda, distblk *min_cost);
extern int    Mode_Decision_for_8x8IntraBlocks_JM_Low444 (Macroblock *currMB, int b8, int lambda, distblk *min_cost);
extern int    Mode_Decision_for_8x8IntraBlocks_JM_High444(Macroblock *currMB, int b8, int lambda, distblk *min_cost);

extern distblk  rdcost_for_8x8_intra_blocks                (Macroblock *currMB, int *c_nz, int b8, int ipmode, int lambda, distblk min_rdcost, int mostProbableMode);
extern distblk  rdcost_for_8x8_intra_blocks_444            (Macroblock *currMB, int *c_nz, int b8, int ipmode, int lambda, distblk min_rdcost, int mostProbableMode);
extern void   compute_satd8x8_cost(VideoParameters *p_Vid, imgpel **cur_img, imgpel **mpr8x8, int pic_opix_x, distblk *cost, distblk min_cost);
extern void   compute_sse8x8_cost (VideoParameters *p_Vid, imgpel **cur_img, imgpel **mpr8x8, int pic_opix_x, distblk *cost, distblk min_cost);
extern void   compute_sad8x8_cost (VideoParameters *p_Vid, imgpel **cur_img, imgpel **mpr8x8, int pic_opix_x, distblk *cost, distblk min_cost);
extern void   compute_comp8x8_cost(VideoParameters *p_Vid, imgpel **cur_img, imgpel **mpr8x8, int pic_opix_x, distblk *cost, distblk min_cost);

extern int    dct_8x8       (Macroblock *currMB, ColorPlane pl, int b8, int *coeff_cost, int intra);
extern int    dct_8x8_cavlc (Macroblock *currMB, ColorPlane pl, int b8, int *coeff_cost, int intra);
extern int    dct_8x8_ls    (Macroblock *currMB, ColorPlane pl, int b8, int *coeff_cost, int intra);

#endif //_TRANSFORM8X8_H_
