
/*!
************************************************************************
* \file
*     me_epzs.h
*
* \author
*    Alexis Michael Tourapis        <alexis.tourapis@dolby.com>
*
* \date
*    11. August 2006
*
* \brief
*    Headerfile for EPZS motion estimation
**************************************************************************
*/


#ifndef _ME_EPZS_H_
#define _ME_EPZS_H_
#include "me_epzs_common.h"

// Functions
extern int EPZSPelBlockMotionSearch      (Macroblock *, MotionVector *, MEBlock *, int, int);
extern int EPZSPelBlockMotionSearchSubMB (Macroblock *, MotionVector *, MEBlock *, int, int);
extern int EPZSBiPredBlockMotionSearch   (Macroblock *, int, MotionVector *, MotionVector *, MotionVector *, MotionVector *, MEBlock *,  int, int, int);
extern int EPZSSubPelBlockMotionSearch   (Macroblock *, MotionVector *, MEBlock *mv_block, int, int*);
extern int EPZSSubPelBlockSearchBiPred   (Macroblock *,  MEBlock *mv_block, int list, 
                                         MotionVector *pred_mv1, MotionVector *pred_mv2, MotionVector *mv1, MotionVector *mv2, int min_mcost, int *lambda_factor);


#endif

