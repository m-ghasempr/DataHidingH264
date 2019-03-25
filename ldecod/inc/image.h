
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
extern StorablePicture *dec_picture_JV[MAX_PLANE];  //!< dec_picture to be used during 4:4:4 indpendent mode decoding

void calculate_frame_no(StorablePicture *p);
void find_snr(struct snr_par *snr, StorablePicture *p, int p_ref);
void get_block(int ref_frame, StorablePicture **list, int x_pos, int y_pos, struct img_par *img, int block[MB_BLOCK_SIZE][MB_BLOCK_SIZE]);
int  picture_order(struct img_par *img);

// For 4:4:4 independent mode
void copy_dec_picture_JV( StorablePicture *dst, StorablePicture *src );

#endif

