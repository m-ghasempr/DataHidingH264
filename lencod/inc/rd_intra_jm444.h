/*!
 ***************************************************************************
 * \file
 *    rd_intra_jm444.h
 *
 * \author
 *    Alexis Michael Tourapis
 *
 * \date
 *    2 January 2008
 *
 * \brief
 *    Headerfile for JM rd based intra mode decision (High444)
 **************************************************************************
 */

#ifndef _RD_INTRA_JM444_H_
#define _RD_INTRA_JM444_H_

extern void Intra16x16_Mode_Decision444                (Macroblock *currMB);
extern int Mode_Decision_for_4x4IntraBlocks_JM_High444 (Macroblock *currMB, int  b8,  int  b4,  int  lambda,  distblk*  min_cost);
extern int Mode_Decision_for_4x4IntraBlocks_JM_Low444  (Macroblock *currMB, int  b8,  int  b4,  int  lambda,  distblk*  min_cost);
#endif

