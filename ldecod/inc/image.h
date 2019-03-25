
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

extern StorablePicture *dec_picture;
extern StorablePicture *dec_picture_JV[MAX_PLANE];  //!< dec_picture to be used during 4:4:4 independent mode decoding

extern void calculate_frame_no(StorablePicture *p);
extern void find_snr(struct snr_par *snr, StorablePicture *p, int *p_ref);
extern int  picture_order(ImageParameters *img);

extern void decode_one_slice(ImageParameters *img, Slice *currSlice, struct inp_par *inp);
extern int  read_new_slice(Slice *currSlice);
extern void exit_picture(StorablePicture **dec_picture);
extern int  decode_one_frame(ImageParameters *img,struct inp_par *inp, struct snr_par *snr);

extern int  is_new_picture(StorablePicture *dec_picture, Slice *currSlice, OldSliceParams *p_old_slice);
extern void init_old_slice(OldSliceParams *p_old_slice);
// For 4:4:4 independent mode
extern void copy_dec_picture_JV( StorablePicture *dst, StorablePicture *src );

extern void frame_postprocessing(ImageParameters *img);
extern void field_postprocessing(ImageParameters *img);


#endif

