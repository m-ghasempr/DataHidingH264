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
 * \file biariencode.c
 *
 * \brief
 *    Routines for binary arithmetic encoding
 * \author
 *    Main contributors (see contributors.h for copyright, address and affiliation details)
 *    - Detlev Marpe                    <marpe@hhi.de>
 *    - Gabi Blaettermann               <blaetter@hhi.de>
 *************************************************************************************
 */
#include <stdlib.h>
#include <math.h>
#include "global.h"
#include "biariencode.h"
#include <assert.h>


/*!
 ************************************************************************
 * Macro for writing bytes of code
 ***********************************************************************
 */
#ifdef NEW_CONSTRAINT_AC
#define put_byte() { \
  Ecodestrm[(*Ecodestrm_len)++] = Ebuffer; \
  Ebits_to_go = 8; \
  while (eep->C > 7) { \
  eep->C-=8; \
  eep->E++; \
  } \
                    }
#else
#define put_byte() { \
  Ecodestrm[(*Ecodestrm_len)++] = Ebuffer; \
  Ebits_to_go = 8; \
                    }
#endif

/*!
 ************************************************************************
 * \brief
 *    Allocates memory for the EncodingEnvironment struct
 ************************************************************************
 */
EncodingEnvironmentPtr arienco_create_encoding_environment()
{
  EncodingEnvironmentPtr eep;

  if ( (eep = (EncodingEnvironmentPtr) calloc(1,sizeof(EncodingEnvironment))) == NULL)
    no_mem_exit("arienco_create_encoding_environment: eep");

  return eep;
}



/*!
 ************************************************************************
 * \brief
 *    Frees memory of the EncodingEnvironment struct
 ************************************************************************
 */
void arienco_delete_encoding_environment(EncodingEnvironmentPtr eep)
{
  if (eep == NULL)
  {
    snprintf(errortext, ET_SIZE, "Error freeing eep (NULL pointer)");
    error (errortext, 200);
  }
  else
    free(eep);
}



/*!
 ************************************************************************
 * \brief
 *    Initializes the EncodingEnvironment for the arithmetic coder
 ************************************************************************
 */
void arienco_start_encoding(EncodingEnvironmentPtr eep,
                            unsigned char *code_buffer,
                            int *code_len, int *last_startcode, int slice_type )
{
  Elow = 0;
  Ebits_to_follow = 0;
  Ebuffer = 0;
  Ebits_to_go = 8;

  Ecodestrm = code_buffer;
  Ecodestrm_len = code_len;
  Ecodestrm_laststartcode = last_startcode;

  /* Only necessary for new AC */
  Erange = CACM99_HALF;

  if (slice_type == INTRA_IMG) //INTRA_SLICE
    eep->AC_next_state_MPS_64 = AC_next_state_MPS_64_INTRA;
  else
    eep->AC_next_state_MPS_64 = AC_next_state_MPS_64_INTER;

#ifdef NEW_CONSTRAINT_AC
  eep->C = 0;
  eep->B = *code_len;
  eep->E = 0;
#endif

}

/*!
 ************************************************************************
 * \brief
 *    Returns the number of currently written bits
 ************************************************************************
 */
int arienco_bits_written(EncodingEnvironmentPtr eep)
{
   return (8 * (*Ecodestrm_len-*Ecodestrm_laststartcode) + Ebits_to_follow + 8  - Ebits_to_go);
}

/*!
 ************************************************************************
 * \brief
 *    Determines the trailing bits needed for the termination of the 
 *    arihmetic coder
 ************************************************************************
 */
int get_trailing_bits(EncodingEnvironmentPtr eep)
{
  int nbits;
  unsigned int roundup, value, bits;
  
  for (nbits = 1; nbits <= B_BITS; nbits++)
  {
    roundup = (1 << (B_BITS - nbits)) - 1;
    bits = (Elow + roundup) >> (B_BITS - nbits);
    value = bits << (B_BITS - nbits);
    if (Elow <= value && value + roundup <= (Elow + (Erange - 1)) )
      break;
  }
  return (nbits);
}

/*!
 ************************************************************************
 * \brief
 *    Terminates the arithmetic coder and writes the trailing bits
 ************************************************************************
 */
void arienco_done_encoding(EncodingEnvironmentPtr eep)
{
  int nbits, i;
  unsigned int roundup, bits, value, bit;
  
  for (nbits = 1; nbits <= B_BITS; nbits++)
  {
    roundup = (1 << (B_BITS - nbits)) - 1;
    bits = (Elow + roundup) >> (B_BITS - nbits);
    value = bits << (B_BITS - nbits);
    if (Elow <= value && value + roundup <= (Elow + (Erange - 1)) )
      break;
  }
  /* output the nbits integer bits */
  for (i = 1; i <= nbits; i++)        
  {
    bit = (bits >> (nbits-i)) & 1;
    Ebuffer <<= 1;
    Ebuffer |= bit;
    if (--Ebits_to_go == 0) 
      put_byte();

    while (Ebits_to_follow > 0) 
    {
      Ebits_to_follow--;
      Ebuffer <<= 1;
      Ebuffer |= (!bit);
      if (--Ebits_to_go == 0) 
        put_byte();
    }
  }
#ifdef NEW_CONSTRAINT_AC

  while (eep->C > 7)
  { 
    eep->C-=8; 
    eep->E++; 
  }
  eep->B= 8*(*Ecodestrm_len - eep->B) + 8 - Ebits_to_go; // no of written bits
  eep->E= eep->E*8 + eep->C; // no of processed bins
  eep->E>>=2; // E/4

  if (eep->B >= 4*(img->current_mb_nr-img->currentSlice->start_mb_nr) ) // only if avg. # bits/MB >= 4
  { 
    if (eep->E > eep->B ) // check whether upper limit of R=4 is exceeded
    {       
      for (i = 0; i < eep->E - eep->B; i++)   // # stuffing bits: (E- R*B)/R     
      {
        Ebuffer <<= 1;  
        Ebuffer |= 0x01;
        if (--Ebits_to_go == 0) 
          put_byte();
  
      }
    } 
  } 

#endif

  if (Ebits_to_go != 8)
    // prevent SCE by filling lower bits with ones
    Ecodestrm[(*Ecodestrm_len)++] = (Ebuffer << Ebits_to_go) | ((1<<Ebits_to_go)-1);
 
}


/*!
 ************************************************************************
 * \brief
 *    Actually arithmetic encoding of one binary symbol by using
 *    the symbol counts of its associated context model
 ************************************************************************
 */
void biari_encode_symbol(EncodingEnvironmentPtr eep, signed short symbol, BiContextTypePtr bi_ct )
{
  unsigned long rLPS;
  unsigned long half, quarter;

  half = CACM99_HALF;
  quarter = CACM99_QUARTER;

  /* covers all cases where code does not bother to shift down symbol to be 
   * either 0 or 1, e.g. in some cases for cbp, mb_Type etc the code symply 
   * masks off the bit position and passes in the resulting value */

  if (symbol != 0) 
    symbol = 1;

  rLPS = rLPS_table_64x4[bi_ct->state][(Erange-0x4001)>>12];

  if (symbol != bi_ct->MPS) 
  {
    Elow += Erange - rLPS;
    Erange = rLPS;

    if (!bi_ct->state)
      bi_ct->MPS = bi_ct->MPS ^ 1;               // switch LPS if necessary
    bi_ct->state = AC_next_state_LPS_64[bi_ct->state]; // next state
  } 
  else 
  {
    Erange -= rLPS;
    bi_ct->state = eep->AC_next_state_MPS_64[bi_ct->state]; // next state
  }

  /* renormalise, as for arith_encode */
    
  while (Erange <= quarter)
  {
    if (Elow >= half)
    {
      /* BIT_PLUS_FOLLOW(1) */;
      Ebuffer <<= 1;
      Ebuffer |= 0x01;
      if (--Ebits_to_go == 0) 
        put_byte();
        
      while (Ebits_to_follow > 0) 
      {
        Ebits_to_follow--;
        Ebuffer <<= 1;
        if (--Ebits_to_go == 0) 
          put_byte();
      }
      Elow -= half;
    }
    else if (Elow+Erange <= half)
    {
      /* BIT_PLUS_FOLLOW(0); */
      Ebuffer <<= 1;
      if (--Ebits_to_go == 0) 
        put_byte();
        
      while (Ebits_to_follow > 0) 
      {
        Ebits_to_follow--;
        Ebuffer <<= 1;
        Ebuffer |= 0x01;
        if (--Ebits_to_go == 0) 
          put_byte();
      }
    }
    else
    {
      Ebits_to_follow++;
      Elow -= quarter;
    }
    Elow <<= 1;
    Erange <<= 1;
  }
#ifdef NEW_CONSTRAINT_AC
  eep->C++;
#endif

}




/*!
 ************************************************************************
 * \brief
 *    Arithmetic encoding of one binary symbol assuming 
 *    a fixed prob. distribution with p(symbol) = 0.5
 ************************************************************************
 */
void biari_encode_symbol_eq_prob(EncodingEnvironmentPtr eep, signed short symbol)
{
  unsigned int half = CACM99_HALF;
  unsigned int quarter = CACM99_QUARTER;

  
  Erange >>= 1; 
  if (symbol != 0)
  {
    Elow += Erange;
  }

  /* renormalise, as for biari_encode */ 
  if (Elow >= half)
  {
    /* BIT_PLUS_FOLLOW(1) */;
    Ebuffer <<= 1;
    Ebuffer |= 0x01;
    if (--Ebits_to_go == 0) 
      put_byte();
        
    while (Ebits_to_follow > 0) 
    {
      Ebits_to_follow--;
      Ebuffer <<= 1;
      if (--Ebits_to_go == 0) 
        put_byte();
    }
    Elow -= half;
  }
  else 
    if (Elow+Erange <= half)
    {
      /* BIT_PLUS_FOLLOW(0); */
      Ebuffer <<= 1;
      if (--Ebits_to_go == 0) 
          put_byte();
    
      while (Ebits_to_follow > 0) 
      {
        Ebits_to_follow--;
        Ebuffer <<= 1;
        Ebuffer |= 0x01;
        if (--Ebits_to_go == 0) 
          put_byte();     
      }
    }
    else
    {
      Ebits_to_follow++;
      Elow -= quarter;
    }

  Elow <<= 1;
  Erange <<=1;
#ifdef NEW_CONSTRAINT_AC
  eep->C++;
#endif
}


/*!
 ************************************************************************
 * \brief
 *    Initializes a given context with some pre-defined probabilities
 *    and a maximum symbol count for triggering the rescaling
 ************************************************************************
 */
void biari_init_context (BiContextTypePtr ctx, const int* ini)
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

