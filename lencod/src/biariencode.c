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

/*!
 ************************************************************************
 * Macro for writing bytes of code
 ***********************************************************************
 */
#define put_byte() { \
                     Ecodestrm[(*Ecodestrm_len)++] = Ebuffer; \
                     Ebits_to_go = 8; \
                    }

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
                            int *code_len )
{
  Elow = 0;
  Ebits_to_follow = 0;
  Ebuffer = 0;
  Ebits_to_go = 8;

  Ecodestrm = code_buffer;
  Ecodestrm_len = code_len;

  /* Only necessary for new AC */
  Erange = CACM99_HALF;			
}

/*!
 ************************************************************************
 * \brief
 *    Returns the number of currently written bits
 ************************************************************************
 */
int arienco_bits_written(EncodingEnvironmentPtr eep)
{
   return (8 * (*Ecodestrm_len) + Ebits_to_follow + 8 + 2 - Ebits_to_go);
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
    Ebuffer >>= 1;
    Ebuffer |= (bit<<7);
    if (--Ebits_to_go == 0) 
      put_byte();
				
    while (Ebits_to_follow > 0) 
    {
      Ebits_to_follow--;
      Ebuffer >>= 1;
      Ebuffer |= ((!bit)<<7);
      if (--Ebits_to_go == 0) 
        put_byte();
    }
  }
  if (Ebits_to_go != 8)
    Ecodestrm[(*Ecodestrm_len)++] = (Ebuffer >> Ebits_to_go);		
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
  int LPS;
  unsigned long cLPS, rLPS;
  unsigned long c0, c1;
  unsigned long half, quarter;
  unsigned long out_r;

  c0 = bi_ct->cum_freq[0]-bi_ct->cum_freq[1];
  c1 = bi_ct->cum_freq[1];
  half = CACM99_HALF;
  quarter = CACM99_QUARTER;

  /* covers all cases where code does not bother to shift down symbol to be 
   * either 0 or 1, e.g. in some cases for cbp, mb_Type etc the code symply 
	 * masks off the bit position and passes in the resulting value */
  if (symbol != 0) 
	  symbol = 1;
  
  /* From frequencies (c0 and c1) determine least probable symbol (LPS) and its
   * count (cLPS)				
   */
  if (c0 < c1) 		
  {				
    LPS = 0;		
    cLPS = c0;
  } 
  else 
  {
    LPS = 1;
    cLPS = c1;
  }

#if  AAC_FRAC_TABLE
  out_r = (Erange*(ARITH_CUM_FREQ_TABLE[bi_ct->cum_freq[0]]>>10))>>16;
#else
  out_r = Erange / (c0+c1);
#endif
  rLPS = out_r * cLPS;
  
  if (symbol == LPS) 
  {
    Elow += Erange - rLPS;
    Erange = rLPS;
  } 
  else 
  {
    Erange -= rLPS;
  }

  if (symbol != 0)
      bi_ct->cum_freq[1]++;

	if (++bi_ct->cum_freq[0] >= bi_ct->max_cum_freq) 
		rescale_cum_freq(bi_ct);
	

  /* renormalise, as for arith_encode */
  do 
  {
    while (Erange <= quarter)
    {
      if (Elow >= half)
      {
        /* BIT_PLUS_FOLLOW(1) */;
				Ebuffer >>= 1;
				Ebuffer |= 0x80;
				if (--Ebits_to_go == 0) 
					put_byte();
				
				while (Ebits_to_follow > 0) 
				{
					Ebits_to_follow--;
					Ebuffer >>= 1;
					if (--Ebits_to_go == 0) 
						put_byte();
				}
        Elow -= half;
      }
      else if (Elow+Erange <= half)
      {
        /* BIT_PLUS_FOLLOW(0); */
			  Ebuffer >>= 1;
			  if (--Ebits_to_go == 0) 
				  put_byte();
			  
			  while (Ebits_to_follow > 0) 
			  {
				  Ebits_to_follow--;
				  Ebuffer >>= 1;
				  Ebuffer |= 0x80;
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
  } while (0);
}


/*!
 ************************************************************************
 * \brief
 *    Initializes a given context with some pre-defined probabilities
 *    and a maximum symbol count for triggering the rescaling
 ************************************************************************
 */
void biari_init_context( BiContextTypePtr ctx, int ini_count_0, int ini_count_1, int max_cum_freq )
{
  ctx->in_use       = TRUE;
  ctx->max_cum_freq = max_cum_freq;


  ctx->cum_freq[1]  = ini_count_1;
  ctx->cum_freq[0]  = ini_count_0 + ini_count_1;
}

/*!
 ************************************************************************
 * \brief
 *    Copies the content (symbol counts) of a given context
 ************************************************************************
 */
void biari_copy_context( BiContextTypePtr ctx_orig, BiContextTypePtr ctx_dest )
{

  ctx_dest->in_use     =  ctx_orig->in_use;
  ctx_dest->max_cum_freq = ctx_orig->max_cum_freq;

  ctx_dest->cum_freq[1] = ctx_orig->cum_freq[1];
  ctx_dest->cum_freq[0] = ctx_orig->cum_freq[0];

}

/*!
 ************************************************************************
 * \brief
 *    Prints the content (symbol counts) of a given context model
 ************************************************************************
 */
void biari_print_context( BiContextTypePtr ctx )
{

  printf("0: %4d\t",ctx->cum_freq[0] - ctx->cum_freq[1]);
  printf("1: %4d",ctx->cum_freq[1]);

}


/*!
 ************************************************************************
 * \brief
 *    Rescales a given context model by halvening the symbol counts
 ************************************************************************
 */
void rescale_cum_freq( BiContextTypePtr   bi_ct)
{

  int old_cum_freq_of_one = bi_ct->cum_freq[1];

  bi_ct->cum_freq[1] = (bi_ct->cum_freq[1] + 1) >> 1;
  bi_ct->cum_freq[0] = bi_ct->cum_freq[1] +
                         ( ( bi_ct->cum_freq[0] - old_cum_freq_of_one + 1 ) >> 1);
}

