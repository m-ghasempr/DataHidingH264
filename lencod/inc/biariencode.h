
/*!
 ***************************************************************************
 * \file
 *    biariencode.h
 *
 * \brief
 *    Headerfile for binary arithmetic encoding routines
 *
 * \author
 *    - Detlev Marpe
 *    - Gabi Blaettermann
 *    - Gunnar Marten
 *
 *    Copyright (C) 2000 HEINRICH HERTZ INSTITUTE All Rights Reserved.
 *
 * \date
 *    21. Oct 2000
 **************************************************************************
 */


#ifndef _BIARIENCOD_H_
#define _BIARIENCOD_H_


/************************************************************************
 * D e f i n i t i o n s
 ***********************************************************************
 */

// some definitions to increase the readability of the source code

#define B_BITS         10    // Number of bits to represent the whole coding interval
#define BITS_TO_LOAD   16
#define MAX_BITS       26          //(B_BITS + BITS_TO_LOAD)
#define ONE            0x04000000  //(1 << MAX_BITS)
#define ONE_M1         0x03FFFFFF  //(ONE - 1)
#define HALF           0x0200      //(1 << (B_BITS-1))
#define QUARTER        0x0100      //(1 << (B_BITS-2))
#define MIN_BITS_TO_GO 0
#define B_LOAD_MASK    0xFFFF      // ((1<<BITS_TO_LOAD) - 1)

int get_pic_bin_count(void);
void reset_pic_bin_count(void);
void set_pic_bin_count(EncodingEnvironmentPtr eep);

void arienco_start_encoding(EncodingEnvironmentPtr eep, unsigned char *code_buffer, int *code_len);
void arienco_reset_EC(EncodingEnvironmentPtr eep);
void arienco_done_encoding(EncodingEnvironmentPtr eep);
void biari_init_context (BiContextTypePtr ctx, const char* ini);
void rescale_cum_freq(BiContextTypePtr bi_ct);
void biari_encode_symbol(EncodingEnvironmentPtr eep, signed short symbol, BiContextTypePtr bi_ct );
void biari_encode_symbol_eq_prob(EncodingEnvironmentPtr eep, signed short symbol);
void biari_encode_symbol_final(EncodingEnvironmentPtr eep, signed short symbol);

/*!
************************************************************************
* \brief
*    Returns the number of currently written bits
************************************************************************
*/
static inline int arienco_bits_written(EncodingEnvironmentPtr eep)
{
  return (((*eep->Ecodestrm_len) + eep->Epbuf + 1) << 3) + (eep->Echunks_outstanding * BITS_TO_LOAD) + BITS_TO_LOAD - eep->Ebits_to_go;
}

#endif  // BIARIENCOD_H

