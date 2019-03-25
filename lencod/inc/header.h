
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

int SliceHeader(Slice* currSlice);
int Partition_BC_Header(int PartNo);
int writeERPS(SyntaxElement *sym, DataPartition *partition);
int dec_ref_pic_marking(Bitstream *bitstream, DecRefPicMarking_t *p_drpm, int idr_flag, int no_output_of_prior_pics_flag, int long_term_reference_flag );

#endif

