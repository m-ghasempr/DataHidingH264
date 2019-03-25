
/*!
 ************************************************************************
 * \file block.h
 *
 * \brief
 *    constant arrays for single block processing
 *
 * \author
 *    Inge Lille-Langoy               <inge.lille-langoy@telenor.com>    \n
 *    Telenor Satellite Services                                         \n
 *    P.O.Box 6914 St.Olavs plass                                        \n
 *    N-0130 Oslo, Norway
 *
 ************************************************************************
 */

#ifndef _BLOCK_H_
#define _BLOCK_H_

extern const byte SNGL_SCAN[16][2];
extern const byte FIELD_SCAN[16][2]; 
extern const byte FIELD_SCAN8x8[64][2];
extern const byte SNGL_SCAN8x8[64][2];
//! look up tables for FRExt_chroma support
extern const unsigned char subblk_offset_x[3][8][4];
extern const unsigned char subblk_offset_y[3][8][4];


#endif

