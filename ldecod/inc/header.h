/*!
 *************************************************************************************
 * \file header.h
 *
 * \brief
 *    Prototypes for header.c
 *************************************************************************************
 */

#ifndef _HEADER_H_
#define _HEADER_H_

extern int FirstPartOfSliceHeader(Slice *currSlice);
extern int RestOfSliceHeader     (Slice *currSlice);

extern void dec_ref_pic_marking(ImageParameters *p_Img, Bitstream *currStream);

extern void decode_poc(ImageParameters *p_Img);
extern int dumppoc(ImageParameters *p_Img);

#endif

