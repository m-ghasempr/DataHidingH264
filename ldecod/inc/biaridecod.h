/*
***********************************************************************
* COPYRIGHT AND WARRANTY INFORMATION
*
* Copyright 2001, International Telecommunications Union, Geneva
*
* DISCLAIMER OF WARRANTY
*
* These software programs are available to the user without any
* license fee or royalty on an "as is" basis. The ITU disclaims
* any and all warranties, whether express, implied, or
* statutory, including any implied warranties of merchantability
* or of fitness for a particular purpose.  In no event shall the
* contributor or the ITU be liable for any incidental, punitive, or
* consequential damages of any kind whatsoever arising from the
* use of these programs.
*
* This disclaimer of warranty extends to the user of these programs
* and user's customers, employees, agents, transferees, successors,
* and assigns.
*
* The ITU does not represent or warrant that the programs furnished
* hereunder are free of infringement of any third-party patents.
* Commercial implementations of ITU-T Recommendations, including
* shareware, may be subject to royalty fees to patent holders.
* Information regarding the ITU-T patent policy is available from
* the ITU Web site at http://www.itu.int.
*
* THIS IS NOT A GRANT OF PATENT RIGHTS - SEE THE ITU-T PATENT POLICY.
************************************************************************
*/

/*!
 ***************************************************************************
 * \file
 *    biaridecod.h
 *
 * \brief
 *    Headerfile for binary arithmetic decoder routines
 *
 * \author
 *    Detlev Marpe,
 *    Gabi Blättermann
 *    Copyright (C) 2000 HEINRICH HERTZ INSTITUTE All Rights Reserved.
 *
 * \date
 *    21. Oct 2000
 **************************************************************************
 */

#ifndef _BIARIDECOD_H_
#define _BIARIDECOD_H_


/************************************************************************
 * D e f i n i t i o n s
 ***********************************************************************
 */

void arideco_start_decoding(DecodingEnvironmentPtr eep, unsigned char *code_buffer, int firstbyte, int *code_len, int slice_type);
int  arideco_bits_read(DecodingEnvironmentPtr dep);
void arideco_done_decoding(DecodingEnvironmentPtr dep);
void biari_init_context (struct img_par *img, BiContextTypePtr ctx, const int* ini);
void rescale_cum_freq(BiContextTypePtr bi_ct);
unsigned int biari_decode_symbol(DecodingEnvironmentPtr dep, BiContextTypePtr bi_ct );
unsigned int biari_decode_symbol_eq_prob(DecodingEnvironmentPtr dep);
unsigned int biari_decode_final(DecodingEnvironmentPtr dep);

#endif  // BIARIDECOD_H_

