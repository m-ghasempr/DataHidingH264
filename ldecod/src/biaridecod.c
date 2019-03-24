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
                            int firstbyte, int *cpixcode_len )
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
    dep->Dvalue += dep->Dvalue  + (Dbuffer & 1);
    Dbuffer >>= 1;
  }
  dep->Dlow = 0;

  /* Only necessary for new AC */
  dep->Drange = CACM99_HALF;
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

  /* scope LPS, cLPS, rLPS, rMPS, and in_r */
  {
    int LPS;
    unsigned int cLPS, rLPS, rMPS;
    unsigned int in_r;

    /* scope c0 and c1 */
    {
	    register unsigned short c0, c1, one_over_c0;

      one_over_c0 = ARITH_CUM_FREQ_TABLE[bi_ct->cum_freq[0]]>>10;
	    c0 = bi_ct->cum_freq[0]-bi_ct->cum_freq[1];
	    c1 = bi_ct->cum_freq[1];

      /* Unit interval size */
#if  AAC_FRAC_TABLE
      in_r = ((range*one_over_c0)>>16);	  
#else
      in_r = range / (c0+c1);
#endif

      /* Determine LPS and cLPS. This should compile to conditional moves. */
	    LPS = (c0 <  c1) ? 0 : 1;
	    cLPS = (c0 < c1) ? c0 : c1;

    }

    /* Size of interval for LPS */
    rLPS = in_r * cLPS;
    rMPS = range - rLPS;
    /* Always set LPS interval at upper end of range. This always allocates the 
     * excess probabilty caused by the truncation during division to the MPS. 
     * Check if value lies in this upper end, if so, then bit is LPS otherwise 
     * bit is MPS */
    if (value >= rMPS) 
    {
      bit = LPS;
      value -= rMPS;
      range = rLPS;
    } 
    else 
    {
      bit = !(LPS);
      range -= rLPS;
    }
  }
  
  /* Increase model frequency counts appropriately */
  bi_ct->cum_freq[1] += bit;

  if (++bi_ct->cum_freq[0] >= bi_ct->max_cum_freq) 
		rescale_cum_freq(bi_ct);


  /* Renormalise when range <= 1/4 of max frequency count to maintain precise 
   * division of code space, minimizing loss of compression effectiveness. 
   * The idea is to prevent underflow, i.e. when the scaling operation 
   * maps different symbols of the model onto the same integer. The (low, high) 
   * pair can lie between:
   * - [0, 1/2) 
   * - [1/2, top) 
   * - [1/4, 3/4) 
   * The (low, high) pair can only become close together when they straddle 1/2 
   * (3rd case). We know unambiguously what to output if the (low,high) pair 
   * are, at the smallest:
   * - [1/4,1/2+epsilon)
   * - [1/2-epsilon,3/4)
   * This range is <= 1/4 of the max frequency count, hence we can (should) rescale 
   * at this time, i.e. when range <= 1/4
   */
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
    value += (Dbuffer & 1);
    Dbuffer >>= 1;   
  }
  
  dep->Drange = range;
  dep->Dvalue = value;

  return(bit);
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
 *    biari_copy_context():
 ************************************************************************
 */
void biari_copy_context( BiContextTypePtr ctx_orig, BiContextTypePtr ctx_dest )
{
  ctx_dest->in_use     =  ctx_orig->in_use;
  ctx_dest->max_cum_freq = ctx_orig->max_cum_freq;

  ctx_dest->cum_freq[1] = ctx_orig->cum_freq[1];
  ctx_dest->cum_freq[0] = ctx_orig->cum_freq[0];

  return;
}

/*!
 ***********************************************************************
 * \brief
 *    biari_print_context():
 ***********************************************************************
 */
void biari_print_context( BiContextTypePtr ctx )
{
  printf("0: %4d\t",ctx->cum_freq[0] - ctx->cum_freq[1]);
  printf("1: %4d",ctx->cum_freq[1]);

  return;
}


/*!
 ***********************************************************************
 * \brief
 *    Rescales a given context model by halvening the symbol counts
 *
 ***********************************************************************
 */
void rescale_cum_freq( BiContextTypePtr   bi_ct)
{
  int old_cum_freq_of_one = bi_ct->cum_freq[1];

  bi_ct->cum_freq[1] = (bi_ct->cum_freq[1] + 1) >> 1;
  bi_ct->cum_freq[0] = bi_ct->cum_freq[1] +
    ( ( bi_ct->cum_freq[0] - old_cum_freq_of_one + 1 ) >> 1);
}
