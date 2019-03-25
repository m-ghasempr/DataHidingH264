
/*!
 ************************************************************************
 * \file image.h
 *
 * \brief
 *    prototypes for image.c
 *
 ************************************************************************
 */

#ifndef _IMAGE_H_
#define _IMAGE_H_

#include "mbuffer.h"

extern void calculate_frame_no(ImageParameters *p_Img, StorablePicture *p);
extern void find_snr          (ImageParameters *p_Img, StorablePicture *p, int *p_ref);
extern int  picture_order(ImageParameters *p_Img);

extern void decode_one_slice (Slice *currSlice);
extern int  read_new_slice(Slice *currSlice);
extern void exit_picture(ImageParameters *p_Img, StorablePicture **dec_picture);
extern int  decode_one_frame(ImageParameters *p_Img);

extern int  is_new_picture(StorablePicture *dec_picture, Slice *currSlice, OldSliceParams *p_old_slice);
extern void init_old_slice(OldSliceParams *p_old_slice);
// For 4:4:4 independent mode
extern void copy_dec_picture_JV( ImageParameters *p_Img, StorablePicture *dst, StorablePicture *src );

extern void frame_postprocessing(ImageParameters *p_Img);
extern void field_postprocessing(ImageParameters *p_Img);

#endif

