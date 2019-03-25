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

#define get_byte(){																					\
										Dbuffer = Dcodestrm[(*Dcodestrm_len)++];\
										Dbits_to_go = 7;												\
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
 *    allocates memory
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
 
	int value = 0;

	Dcodestrm = cpixcode;
  Dcodestrm_len = cpixcode_len;
  *Dcodestrm_len = firstbyte;
 
	{
		int i;
		Dbits_to_go = 0;
		for (i = 0; i < B_BITS -1 ; i++) // insertion of redundant bit
		{
			if (--Dbits_to_go < 0)
				get_byte();
			value = (value<<1)  | ((Dbuffer >> Dbits_to_go) & 0x01);
		}
	}
	dep->Drange = HALF-2;
  dep->Dvalue = value;

}


/*!
 ************************************************************************
 * \brief
 *    arideco_bits_read
 ************************************************************************
 */
int arideco_bits_read(DecodingEnvironmentPtr dep)
{
  return 8 * ((*Dcodestrm_len)-1) + (8 - Dbits_to_go) - 16;
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
  register unsigned int bit = bi_ct->MPS;
	register unsigned int value = dep->Dvalue;
  register unsigned int range = dep->Drange;
	register unsigned int rLPS = rLPS_table_64x4[bi_ct->state][(range>>6) & 0x03];

	range -= rLPS;

	if (value < range) /* MPS */ 
  	bi_ct->state = AC_next_state_MPS_64[bi_ct->state]; // next state
  else						  /* LPS */
  {
		value -= range;
		range = rLPS;
		bit = !bit;
    if (!bi_ct->state)			 // switch meaning of MPS if necessary	
			bi_ct->MPS ^= 0x01;              
		bi_ct->state = AC_next_state_LPS_64[bi_ct->state]; // next state 
	}
  
  while (range < QUARTER)
  {
    /* Double range */
    range <<= 1;
    if (--Dbits_to_go < 0) 
	    get_byte();   
    /* Shift in next bit and add to value */
		value = (value << 1) | ((Dbuffer >> Dbits_to_go) & 0x01);

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
  register unsigned int bit = 0;
	register unsigned int value  = (dep->Dvalue<<1);

	if (--Dbits_to_go < 0) 
		get_byte();	
	/* Shift in next bit and add to value */	
	value |= (Dbuffer >> Dbits_to_go) &  0x01;	
	if (value >= dep->Drange) 
	{
		bit = 1;
		value -= dep->Drange;
	}

	dep->Dvalue = value;

  return(bit);
}

/*!
 ************************************************************************
 * \brief
 *    biari_decode_symbol_final():
 * \return
 *    the decoded symbol
 ************************************************************************
 */
unsigned int biari_decode_final(DecodingEnvironmentPtr dep)
{
  register unsigned int value  = dep->Dvalue;
  register unsigned int range  = dep->Drange - 2;
		
	if (value >= range) 
	{
		return 1;
	}
	else
	{
		while (range < QUARTER)
		{
    /* Double range */
			range <<= 1;
			if (--Dbits_to_go < 0) 
				get_byte();   
			/* Shift in next bit and add to value */
			value = (value << 1) | ((Dbuffer >> Dbits_to_go) & 0x01);
		}	
		dep->Dvalue = value;
		dep->Drange = range;
		return 0;
	}
}



/*!
 ************************************************************************
 * \brief
 *    Initializes a given context with some pre-defined probability state
 ************************************************************************
 */
void biari_init_context (struct img_par* img, BiContextTypePtr ctx, const int* ini)
{
  int pstate;

  pstate = ((ini[0]*img->qp)>>4) + ini[1];
  pstate = min (max ( 1, pstate), 126);

  if ( pstate >= 64 )
  {
    ctx->state  = pstate - 64;
    ctx->MPS    = 1;
  }
  else
  {
    ctx->state  = 63 - pstate;
    ctx->MPS    = 0;
  }
}

