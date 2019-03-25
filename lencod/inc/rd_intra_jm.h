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

extern double min_rdcost_16x16                      (Macroblock *currMB, double lambda);
extern void Intra16x16_Mode_Decision_RDopt          (Macroblock *currMB, double lambda);
extern void Intra16x16_Mode_Decision_SAD            (Macroblock *currMB);
extern int Mode_Decision_for_4x4IntraBlocks_JM_High (Macroblock *currMB, int  b8,  int  b4,  double  lambda,  double*  min_cost);
extern int Mode_Decision_for_4x4IntraBlocks_JM_Low  (Macroblock *currMB, int  b8,  int  b4,  double  lambda,  double*  min_cost);

#endif

