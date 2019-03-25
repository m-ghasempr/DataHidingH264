/*!
 ***************************************************************************
 * \file
 *    wp_lms.h
 *
 * \author
 *    Alexis Michael Tourapis
 *
 * \date
 *    22. February 2008
 *
 * \brief
 *    Headerfile for weighted prediction support using LMS
 **************************************************************************
 */

#ifndef _WP_LMS_H_
#define _WP_LMS_H_

extern void EstimateWPPSliceAlg1(Slice *currSlice, int offset);
extern void EstimateWPBSliceAlg1(Slice *currSlice);
extern int  TestWPPSliceAlg1    (Slice *currSlice, int offset);
extern int  TestWPBSliceAlg1    (Slice *currSlice, int method);

#endif

