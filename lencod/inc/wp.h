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

void InitWP(InputParameters *params);
void (*EstimateWPPSlice) (int offset);
void (*EstimateWPBSlice)(void);
int  (*TestWPPSlice)(int offset);
int  (*TestWPBSlice)(int method);

void estimate_weighting_factor_B_slice(void);
void estimate_weighting_factor_P_slice(int offset);
int  test_wp_P_slice(int offset);
int  test_wp_B_slice(int method);

double ComputeImgSum(imgpel **CurrentImage, int height, int width);
#endif

