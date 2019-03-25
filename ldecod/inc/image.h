
/*!
 ************************************************************************
 * \file image.h
 *
 ************************************************************************
 */

#ifndef _IMAGE_H_
#define _IMAGE_H_

#include "mbuffer.h"

// this one is empty. keep it, maybe we will move some image.c function 
// declarations here

extern StorablePicture *dec_picture;

void find_snr(struct snr_par *snr, StorablePicture *p, FILE *p_ref);
void get_block(int ref_frame, StorablePicture **list, int x_pos, int y_pos, struct img_par *img, int block[BLOCK_SIZE][BLOCK_SIZE]);
int  picture_order(struct img_par *img);

#endif

