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
 *************************************************************************************
 * \file biaridecod.c
 *
 * \brief
 *    binary arithmetic decoder routines
 * \date
 *    21. Oct 2000
 * \author
 *    Main contributors (see contributors.h for copyright, address and affiliation details)
 *    - Detlev Marpe                    <marpe@hhi.de>
 *    - Gabi Blaettermann               <blaetter@hhi.de>
 *************************************************************************************
 */

#include <stdlib.h>
#include <math.h>
#include "global.h"
#include "biaridecod.h"
#include "memalloc.h"

extern int symbolCount;

/************************************************************************
 * M a c r o s
 ************************************************************************
 */

#define get_byte(){                                                     \
            Dbuffer = Dcodestrm[(*Dcodestrm_len)++];        \
                        Dbits_to_go = 7;                                \
          }


/************************************************************************
 ************************************************************************
                      init / exit decoder
 ************************************************************************
 ************************************************************************/


/*!
 ************************************************************************
 * \brief
 *    Allocates memory for the DecodingEnvironment struct
 * \return DecodingContextPtr
 *    allocated memory
 ************************************************************************
 */
DecodingEnvironmentPtr arideco_create_decoding_environment()
{
  DecodingEnvironmentPtr dep;

  if ((dep = calloc(1,sizeof(DecodingEnvironment))) == NULL)
    no_mem_exit("arideco_create_decoding_environment: dep");
  return dep;
}


/*!
 ***********************************************************************
 * \brief
 *    Frees memory of the DecodingEnvironment struct
 ***********************************************************************
 */
void arideco_delete_decoding_environment(DecodingEnvironmentPtr dep)
{
  if (dep == NULL)
  {
    snprintf(errortext, ET_SIZE, "Error freeing dep (NULL pointer)");
    error (errortext, 200);
  }
  else
    free(dep);
}


/*!
 ************************************************************************
 * \brief
 *    Initializes the DecodingEnvironment for the arithmetic coder
 ************************************************************************
 */
void arideco_start_decoding(DecodingEnvironmentPtr dep, unsigned char *cpixcode,
                            int firstbyte, int *cpixcode_len, int slice_type )
{
  int i;
  int bits;
  bits = B_BITS;

  Dcodestrm = cpixcode;
  Dcodestrm_len = cpixcode_len;
  *Dcodestrm_len = firstbyte;

  Dbits_to_go = 0;
  dep->Dvalue = 0;

  for (i = 0; i < bits; i++)
  {
    if (--Dbits_to_go < 0)
      get_byte();
    dep->Dvalue <<= 1;
    if (Dbuffer & 0x80) 
      dep->Dvalue |= 0x01;
    Dbuffer <<= 1;
  }
  dep->Dlow = 0;

  /* Only necessary for new AC */
  dep->Drange = CACM99_HALF;

  if (slice_type == INTRA_IMG) //INTRA_SLICE
    dep->AC_next_state_MPS_64 = AC_next_state_MPS_64_INTRA;
  else
    dep->AC_next_state_MPS_64 = AC_next_state_MPS_64_INTER;
}


/*!
 ************************************************************************
 * \brief
 *    arideco_bits_read
 ************************************************************************
 */
int arideco_bits_read(DecodingEnvironmentPtr dep)
{
  return 8 * ((*Dcodestrm_len)-1) + (8 - Dbits_to_go) - CODE_VALUE_BITS;
}


/*!
 ************************************************************************
 * \brief
 *    arideco_done_decoding():
 ************************************************************************
 */
void arideco_done_decoding(DecodingEnvironmentPtr dep)
{
  (*Dcodestrm_len)++;
}



/*!
 ************************************************************************
 * \brief
 *    biari_decode_symbol():
 * \return
 *    the decoded symbol
 ************************************************************************
 */
unsigned int biari_decode_symbol(DecodingEnvironmentPtr dep, BiContextTypePtr bi_ct )
{
  int bit;
  register unsigned int range, value;
  
  range = dep->Drange;
  value = dep->Dvalue;

  /* scope  rLPS, rMPS  */
  {
    unsigned int rLPS, rMPS;
   
    rLPS = rLPS_table_64x4[bi_ct->state][(range-0x4001)>>12];
    rMPS = range - rLPS;

    if (value >= rMPS) 
    {
      value -= rMPS;
      range = rLPS;
      bit = !(bi_ct->MPS);
      if (!bi_ct->state)
        bi_ct->MPS = bi_ct->MPS ^ 1;               // switch LPS if necessary
      bi_ct->state = AC_next_state_LPS_64[bi_ct->state]; // next state
    } 
    else 
    {
      range = rMPS;
      bit = bi_ct->MPS;
      bi_ct->state = dep->AC_next_state_MPS_64[bi_ct->state]; // next state
    }
  }
  
  while (range <= CACM99_QUARTER)
  {
    /* Double range and value */
    range <<= 1;
    if (--Dbits_to_go < 0) 
    {
      get_byte();
    }
    /* Shift in next bit and add to value */
    value <<= 1;
    if (Dbuffer & 0x80)
      value |= 0x01;
    Dbuffer <<= 1;
  }
  
  dep->Drange = range;
  dep->Dvalue = value;

  return(bit);
}


/*!
 ************************************************************************
 * \brief
 *    biari_decode_symbol_eq_prob():
 * \return
 *    the decoded symbol
 ************************************************************************
 */
unsigned int biari_decode_symbol_eq_prob(DecodingEnvironmentPtr dep)
{
  int bit;
  register unsigned int value  = dep->Dvalue;
  register unsigned int range  = (dep->Drange>>1);

  
  if (value >= range) 
  {
    bit = 1;
    value -= range;
  } 
  else 
    bit = 0;

  if (--Dbits_to_go < 0) 
  {
    get_byte();
  }
  /* Shift in next bit and add to value */  
  value <<= 1;
  if (Dbuffer & 0x80)
    value |= 0x01;
  Dbuffer <<= 1;

  dep->Dvalue = value;
  dep->Drange = (range<<1);

  return(bit);
}



/*!
 ************************************************************************
 * \brief
 *    Initializes a given context with some pre-defined probabilities
 *    and a maximum symbol count for triggering the rescaling
 ************************************************************************
 */
void biari_init_context (struct img_par* img, BiContextTypePtr ctx, const int* ini)
{
  int pstate;

  pstate = ((ini[0]*(img->qp-SHIFT_QP))>>4) + ini[1]; 


  if (img->type==INTRA_IMG) pstate = min (max (27, pstate),  74);  // states from 0 to 23
  else                      pstate = min (max ( 0, pstate), 101);  // states from 0 to 50

  if (pstate>=51)
  {
    ctx->state  = pstate - 51;
    ctx->MPS    = 1;
  }
  else
  {
    ctx->state  = 50 - pstate;
    ctx->MPS    = 0;
  }
}

