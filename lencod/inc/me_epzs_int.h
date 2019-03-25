
/*!
 ************************************************************************
 * \file
 *     me_epzs_int.h
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


#ifndef _ME_EPZS_INT_H_
#define _ME_EPZS_INT_H_

#include "me_epzs.h"

// Functions
extern distblk  EPZSIntPelBlockMotionSearch     (Macroblock *, MotionVector *, MEBlock *, distblk, int);
extern distblk  EPZSIntPelBlockMotionSearchSubMB(Macroblock *, MotionVector *, MEBlock *, distblk, int);
extern distblk  EPZSIntBiPredBlockMotionSearch  (Macroblock *, int, MotionVector *, MotionVector *, MotionVector *, MotionVector *, MEBlock *, int, distblk, int);
#endif

