/*!
 ***************************************************************************
 * \file
 *    rd_intra_jm.h
 *
 * \author
 *    Alexis Michael Tourapis
 *
 * \date
 *    2 January 2008
 *
 * \brief
 *    Headerfile for JM rd based intra mode decision
 **************************************************************************
 */

#ifndef _RD_INTRA_JM_H_
#define _RD_INTRA_JM_H_

extern void Intra16x16_Mode_Decision (Macroblock* currMB, int* i16mode, int is_cavlc);
extern int Mode_Decision_for_4x4IntraBlocks_JM_High (Slice *currSlice, Macroblock *currMB, int  b8,  int  b4,  double  lambda,  double*  min_cost, int cr_cbp[3], int is_cavlc);
extern int Mode_Decision_for_4x4IntraBlocks_JM_Low  (Slice *currSlice, Macroblock *currMB, int  b8,  int  b4,  double  lambda,  double*  min_cost, int cr_cbp[3], int is_cavlc);
#endif

