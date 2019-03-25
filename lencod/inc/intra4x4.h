
/*!
 ************************************************************************
 * \file intra4x4.h
 *
 * \brief
 *    intra 4x4 functions
 *
 * \author
 *    - Alexis Michael Tourapis         <alexismt@ieee.org>
 *
 ************************************************************************
 */

#ifndef _INTRA4x4_H_
#define _INTRA4x4_H_

extern void set_intrapred_4x4(Macroblock *currMB, ColorPlane pl, int img_x, int img_y, int *left_available, int *up_available, int *all_available);
extern void get_intrapred_4x4(Macroblock *currMB, ColorPlane pl, int i4x4_mode, int img_x, int img_y, int left_available, int up_available);

#endif

