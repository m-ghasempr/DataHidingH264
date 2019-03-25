
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
extern int  EPZSIntPelBlockMotionSearch     (Macroblock *, MotionVector *, MEBlock *, int, int);
extern int  EPZSIntPelBlockMotionSearchSubMB(Macroblock *, MotionVector *, MEBlock *, int, int);
extern int  EPZSIntBiPredBlockMotionSearch  (Macroblock *, int, MotionVector *, MotionVector *, MotionVector *, MotionVector *, MEBlock *, int, int, int);
#endif

