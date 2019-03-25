
/*!
 ************************************************************************
 * \file img_distortion.h
 *
 * \brief
 *    Distortion related definitions
 *
 * \author
 *    Main contributors (see contributors.h for copyright, address and affiliation details)
 *     - Alexis Michael Tourapis         <alexismt@ieee.org>
 *     - Woo-Shik Kim                    <wooshik.kim@usc.edu>
 *
 ************************************************************************
 */

#ifndef _IMG_DISTORTION_H_
#define _IMG_DISTORTION_H_

extern void accumulate_avslice(DistMetric metric[3], int slice_type, int frames);
extern void accumulate_average(DistMetric metric[3], int frames);
extern void find_distortion   (ImageParameters *p_Img, InputParameters *p_Inp, ImageData *imgData);
extern void select_img        (ImageParameters *p_Img, InputParameters *p_Inp, ImageStructure *imgSRC, ImageStructure *imgREF, ImageData *imgData);
extern void compute_distortion(ImageParameters *p_Img, InputParameters *p_Inp, ImageData *imgData);

#endif

