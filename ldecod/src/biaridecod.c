/*!
 *************************************************************************************
 * \file biaridecod.c
 *
 * \brief
 *   Binary arithmetic decoder routines.
 *
 *   This modified implementation of the M Coder is based on JVT-U084 
 *   with the choice of M_BITS = 16.
 *
 * \date
 *    21. Oct 2000
 * \author
 *    Main contributors (see contributors.h for copyright, address and affiliation details)
 *    - Detlev Marpe                    <marpe@hhi.de>
 *    - Gabi Blaettermann
 *    - Gunnar Marten
 *************************************************************************************
 */

#include "global.h"
#include "memalloc.h"
#include "biaridecod.h"


#define B_BITS    10      // Number of bits to represent the whole coding interval
#define HALF      0x0200  //(1 << (B_BITS-1))
#define QUARTER   0x0100  //(1 << (B_BITS-2))

/* Range table for  LPS */
static const byte rLPS_table_64x4[64][4]=
{
  { 128, 176, 208, 240},
  { 128, 167, 197, 227},
  { 128, 158, 187, 216},
  { 123, 150, 178, 205},
  { 116, 142, 169, 195},
  { 111, 135, 160, 185},
  { 105, 128, 152, 175},
  { 100, 122, 144, 166},
  {  95, 116, 137, 158},
  {  90, 110, 130, 150},
  {  85, 104, 123, 142},
  {  81,  99, 117, 135},
  {  77,  94, 111, 128},
  {  73,  89, 105, 122},
  {  69,  85, 100, 116},
  {  66,  80,  95, 110},
  {  62,  76,  90, 104},
  {  59,  72,  86,  99},
  {  56,  69,  81,  94},
  {  53,  65,  77,  89},
  {  51,  62,  73,  85},
  {  48,  59,  69,  80},
  {  46,  56,  66,  76},
  {  43,  53,  63,  72},
  {  41,  50,  59,  69},
  {  39,  48,  56,  65},
  {  37,  45,  54,  62},
  {  35,  43,  51,  59},
  {  33,  41,  48,  56},
  {  32,  39,  46,  53},
  {  30,  37,  43,  50},
  {  29,  35,  41,  48},
  {  27,  33,  39,  45},
  {  26,  31,  37,  43},
  {  24,  30,  35,  41},
  {  23,  28,  33,  39},
  {  22,  27,  32,  37},
  {  21,  26,  30,  35},
  {  20,  24,  29,  33},
  {  19,  23,  27,  31},
  {  18,  22,  26,  30},
  {  17,  21,  25,  28},
  {  16,  20,  23,  27},
  {  15,  19,  22,  25},
  {  14,  18,  21,  24},
  {  14,  17,  20,  23},
  {  13,  16,  19,  22},
  {  12,  15,  18,  21},
  {  12,  14,  17,  20},
  {  11,  14,  16,  19},
  {  11,  13,  15,  18},
  {  10,  12,  15,  17},
  {  10,  12,  14,  16},
  {   9,  11,  13,  15},
  {   9,  11,  12,  14},
  {   8,  10,  12,  14},
  {   8,   9,  11,  13},
  {   7,   9,  11,  12},
  {   7,   9,  10,  12},
  {   7,   8,  10,  11},
  {   6,   8,   9,  11},
  {   6,   7,   9,  10},
  {   6,   7,   8,   9},
  {   2,   2,   2,   2}
};


static const byte AC_next_state_MPS_64[64] =    
{
  1,2,3,4,5,6,7,8,9,10,
  11,12,13,14,15,16,17,18,19,20,
  21,22,23,24,25,26,27,28,29,30,
  31,32,33,34,35,36,37,38,39,40,
  41,42,43,44,45,46,47,48,49,50,
  51,52,53,54,55,56,57,58,59,60,
  61,62,62,63
};


static const byte AC_next_state_LPS_64[64] =    
{
  0, 0, 1, 2, 2, 4, 4, 5, 6, 7,
  8, 9, 9,11,11,12,13,13,15,15,
  16,16,18,18,19,19,21,21,22,22,
  23,24,24,25,26,26,27,27,28,29,
  29,30,30,30,31,32,32,33,33,33,
  34,34,35,35,35,36,36,36,37,37,
  37,38,38,63
};

static const byte renorm_table_32[32]={6,5,4,4,3,3,3,3,2,2,2,2,2,2,2,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};

#if TRACE
extern int symbolCount;
#endif

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
 *    finalize arithetic decoding():
 ************************************************************************
 */
void arideco_done_decoding(DecodingEnvironmentPtr dep)
{
  (*dep->Dcodestrm_len)++;
#if(TRACE==2)
  fprintf(p_trace, "done_decoding: %d\n", *dep->Dcodestrm_len);
#endif
}

/*!
 ************************************************************************
 * \brief
 *    read one byte from the bitstream
 ************************************************************************
 */
unsigned int getbyte(DecodingEnvironmentPtr dep)
{     
#if(TRACE==2)
  fprintf(p_trace, "get_byte: %d\n", (*dep->Dcodestrm_len));
#endif
  return dep->Dcodestrm[(*dep->Dcodestrm_len)++];
}

/*!
 ************************************************************************
 * \brief
 *    read two bytes from the bitstream
 ************************************************************************
 */
unsigned int getword(DecodingEnvironmentPtr dep)
{
  int d = *dep->Dcodestrm_len;
#if(TRACE==2)
  fprintf(p_trace, "get_byte: %d\n", d);
  fprintf(p_trace, "get_byte: %d\n", d + 1);
#endif
  *dep->Dcodestrm_len += 2;
  return ((dep->Dcodestrm[d]<<8) | dep->Dcodestrm[d+1]);
}
/*!
 ************************************************************************
 * \brief
 *    Initializes the DecodingEnvironment for the arithmetic coder
 ************************************************************************
 */
void arideco_start_decoding(DecodingEnvironmentPtr dep, unsigned char *code_buffer,
                            int firstbyte, int *code_len)
{

  dep->Dcodestrm      = code_buffer;
  dep->Dcodestrm_len  = code_len;
  *dep->Dcodestrm_len = firstbyte;

  dep->Dvalue = getbyte(dep);
  dep->Dvalue = (dep->Dvalue << 16) | getword(dep); // lookahead of 2 bytes: always make sure that bitstream buffer
                                        // contains 2 more bytes than actual bitstream
  dep->DbitsLeft = 15;
  dep->Drange = HALF-2;

#if (2==TRACE)
  fprintf(p_trace, "value: %d firstbyte: %d code_len: %d\n", dep->Dvalue >> dep->DbitsLeft, firstbyte, *code_len);
#endif
}


/*!
 ************************************************************************
 * \brief
 *    arideco_bits_read
 ************************************************************************
 */
int arideco_bits_read(DecodingEnvironmentPtr dep)
{ 
  int tmp = ((*dep->Dcodestrm_len) << 3) - dep->DbitsLeft;

#if (2==TRACE)
  fprintf(p_trace, "tmp: %d\n", tmp);
#endif
  return tmp;
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
  register unsigned int state = bi_ct->state;
  register unsigned int bit   = bi_ct->MPS;
  register unsigned int value = dep->Dvalue;
  register unsigned int range = dep->Drange;
  register int renorm = 1;
  register unsigned int rLPS  = rLPS_table_64x4[state][(range>>6) & 0x03];

  range -= rLPS;

  if(value < (range << dep->DbitsLeft))   //MPS
  {
    bi_ct->state = AC_next_state_MPS_64[state]; // next state 

    if( range >= QUARTER )
    {
      dep->Drange = range;
      return (bit);
    }
    else 
      range <<= 1;

  }
  else         // LPS 
  {
    value -= (range << dep->DbitsLeft);
    bit ^= 0x01;

    if (!state)          // switch meaning of MPS if necessary 
      bi_ct->MPS ^= 0x01; 

    bi_ct->state = AC_next_state_LPS_64[state]; // next state 

    renorm = renorm_table_32[(rLPS>>3) & 0x1F]; 
    range = (rLPS << renorm);
  }

  dep->Drange = range;
  dep->DbitsLeft -= renorm;
  if( dep->DbitsLeft > 0 )
  { 
    dep->Dvalue = value;
    return(bit);
  } 

  dep->Dvalue = (value << 16) | getword(dep);    // lookahead of 2 bytes: always make sure that bitstream buffer
                                                 // contains 2 more bytes than actual bitstream
  dep->DbitsLeft += 16;

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
  register int tmp_value;
  register int value = dep->Dvalue;

  if(--(dep->DbitsLeft) == 0)  
  {
    value = (value << 16) | getword( dep );  // lookahead of 2 bytes: always make sure that bitstream buffer
                                             // contains 2 more bytes than actual bitstream
    dep->DbitsLeft = 16;
  }
  tmp_value  = value - (dep->Drange << dep->DbitsLeft);

  if (tmp_value < 0)
  {
    dep->Dvalue = value;
    return 0;
  }
  else
  {
    dep->Dvalue = tmp_value;
    return 1;
  }
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
  register unsigned int range  = dep->Drange - 2;
  register int value  = dep->Dvalue;
  value -= (range << dep->DbitsLeft);

  if (value < 0) 
  {
    if( range >= QUARTER )
    {
      dep->Drange = range;
      return 0;
    }
    else 
    {   
      dep->Drange = (range << 1);
      if( --(dep->DbitsLeft) > 0 )
        return 0;
      else
      {
        dep->Dvalue = (dep->Dvalue << 16) | getword( dep ); // lookahead of 2 bytes: always make sure that bitstream buffer
                                                            // contains 2 more bytes than actual bitstream
        dep->DbitsLeft = 16;
        return 0;
      }
    }
  }
  else
  {
    return 1;
  }
}

/*!
 ************************************************************************
 * \brief
 *    Initializes a given context with some pre-defined probability state
 ************************************************************************
 */
void biari_init_context (int qp, BiContextTypePtr ctx, const char* ini)
{
  int pstate = ((ini[0]* qp )>>4) + ini[1];
  pstate = iClip3(1, 126, pstate);

  if ( pstate >= 64 )
  {
    ctx->state = (unsigned short) (pstate - 64);
    ctx->MPS   = 1;
  }
  else
  {
    ctx->state = (unsigned short) (63 - pstate);
    ctx->MPS   = 0;
  }
}

