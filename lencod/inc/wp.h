/*!
 ***************************************************************************
 * \file
 *    wp.h
 *
 * \author
 *    Alexis Michael Tourapis
 *
 * \date
 *    22. February 2008
 *
 * \brief
 *    Headerfile for weighted prediction support
 **************************************************************************
 */

#ifndef _WP_H_
#define _WP_H_

#include "wp_lms.h"
#include "wp_mcprec.h"
#include "wp_mciter.h"

#define DEBUG_WP  0

void InitWP              (VideoParameters *p_Vid, InputParameters *p_Inp);

extern void EstimateWPBSliceAlg0(Slice *currSlice);
extern void EstimateWPPSliceAlg0(Slice *currSlice, int offset);
extern int  TestWPPSliceAlg0    (VideoParameters *p_Vid, int offset);
extern int  TestWPBSliceAlg0    (VideoParameters *p_Vid, int method);

extern double ComputeImgSum     (imgpel **CurrentImage, int height, int width);

#endif

