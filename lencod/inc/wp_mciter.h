/*!
 ***************************************************************************
 * \file
 *    wp_mciter.h
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

#ifndef _WP_MCITERM_H_
#define _WP_MCITERM_H_


void EstimateWPBSliceAlg2(ImageParameters *img, InputParameters *params);
void EstimateWPPSliceAlg2(ImageParameters *img, InputParameters *params, int offset);
int  TestWPPSliceAlg2    (ImageParameters *img, InputParameters *params, int offset);
int  TestWPBSliceAlg2    (ImageParameters *img, InputParameters *params, int method);

void compute_offset();

#endif

