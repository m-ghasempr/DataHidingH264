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
 * \file cabac.c
 *
 * \brief
 *    CABAC entropy coding routines
 *
 * \author
 *    Main contributors (see contributors.h for copyright, address and affiliation details)
 *    - Detlev Marpe                    <marpe@hhi.de>
 **************************************************************************************
 */

#include <stdlib.h>
#include <math.h>
#include <memory.h>
#include <assert.h> // for debugging
#include <string.h>
#include "cabac.h"
#include "memalloc.h"
#include "bitsbuf.h"
#include "header.h"
#include "elements.h"

extern const int BLOCK_STEP[8][2];
int symbolCount = 0;

/*!
 ************************************************************************
 * \brief
 *    Allocation of contexts models for the motion info
 *    used for arithmetic decoding
 *
 ************************************************************************
 */
MotionInfoContexts* create_contexts_MotionInfo(void)
{
  MotionInfoContexts *deco_ctx;

  deco_ctx = (MotionInfoContexts*) calloc(1, sizeof(MotionInfoContexts) );
  if( deco_ctx == NULL )
    no_mem_exit("create_contexts_MotionInfo: deco_ctx");

  return deco_ctx;
}


/*!
 ************************************************************************
 * \brief
 *    Allocates of contexts models for the texture info
 *    used for arithmetic decoding
 ************************************************************************
 */
TextureInfoContexts* create_contexts_TextureInfo(void)
{
  TextureInfoContexts *deco_ctx;

  deco_ctx = (TextureInfoContexts*) calloc(1, sizeof(TextureInfoContexts) );
  if( deco_ctx == NULL )
    no_mem_exit("create_contexts_TextureInfo: deco_ctx");

  return deco_ctx;
}


/*!
 ************************************************************************
 * \brief
 *    Initializes an array of contexts models with some pre-defined
 *    counts (ini_flag = 1) or with a flat histogram (ini_flag = 0)
 *
 ************************************************************************
 */
#define INIT_CTX(jj,ii,ctx,ini) {for(j=0;j<jj;j++)for(i=0;i<ii;i++)biari_init_context(img, (ctx)[j]+i,(ini)[j][i]);}


void init_contexts_MotionInfo (struct img_par *img, MotionInfoContexts *enco_ctx)
{
  int i,j;

  INIT_CTX (2, NUM_ABT_MODE_CTX, enco_ctx->ABT_mode_contexts, ABT_MODE_Ini);
  INIT_CTX (4, NUM_MB_TYPE_CTX,  enco_ctx->mb_type_contexts,  MB_TYPE_Ini );
  INIT_CTX (2, NUM_B8_TYPE_CTX,  enco_ctx->b8_type_contexts,  B8_TYPE_Ini );
  INIT_CTX (2, NUM_MV_RES_CTX,   enco_ctx->mv_res_contexts,   MV_RES_Ini  );
  INIT_CTX (2, NUM_REF_NO_CTX,   enco_ctx->ref_no_contexts,   REF_NO_Ini  );
  INIT_CTX (2, NUM_BWD_REF_NO_CTX, enco_ctx->bwd_ref_no_contexts, REF_NO_Ini);
  INIT_CTX (1, NUM_DELTA_QP_CTX, &enco_ctx->delta_qp_contexts, &DELTA_QP_Ini  );
  INIT_CTX (1, NUM_MB_AFF_CTX,   &enco_ctx->mb_aff_contexts,   &MB_AFF_Ini  );

  enco_ctx->slice_term_context.state = 63;
  enco_ctx->slice_term_context.MPS   =  0;
}

/*!
 ************************************************************************
 * \brief
 *    Initializes an array of contexts models with some pre-defined
 *    counts (ini_flag = 1) or with a flat histogram (ini_flag = 0)
 ************************************************************************
 */
void init_contexts_TextureInfo(struct img_par *img, TextureInfoContexts *enco_ctx)
{
  int i,j;
  int intra = (img->type==INTRA_IMG ? 1 : 0);

  INIT_CTX (3, NUM_CBP_CTX,  enco_ctx->cbp_contexts,    CBP_Ini[!intra]);
  INIT_CTX (9, NUM_IPR_CTX,  enco_ctx->ipr_contexts,    IPR_Ini   );
  INIT_CTX (1, NUM_CIPR_CTX, &enco_ctx->cipr_contexts,   &CIPR_Ini  ); //GB


  INIT_CTX (NUM_BLOCK_TYPES, NUM_BCBP_CTX,  enco_ctx->bcbp_contexts, BCBP_Ini[intra]);
  INIT_CTX (NUM_BLOCK_TYPES, NUM_MAP_CTX,   enco_ctx->map_contexts,  MAP_Ini [intra]);
  INIT_CTX (NUM_BLOCK_TYPES, NUM_LAST_CTX,  enco_ctx->last_contexts, LAST_Ini[intra]);
  INIT_CTX (NUM_BLOCK_TYPES, NUM_ONE_CTX,   enco_ctx->one_contexts,  ONE_Ini [intra]);
  INIT_CTX (NUM_BLOCK_TYPES, NUM_ABS_CTX,   enco_ctx->abs_contexts,  ABS_Ini [intra]);
}


/*!
 ************************************************************************
 * \brief
 *    Frees the memory of the contexts models
 *    used for arithmetic decoding of the motion info.
 ************************************************************************
 */
void delete_contexts_MotionInfo(MotionInfoContexts *deco_ctx)
{
  if( deco_ctx == NULL )
    return;

  free( deco_ctx );

  return;
}


/*!
 ************************************************************************
 * \brief
 *    Frees the memory of the contexts models
 *    used for arithmetic decoding of the texture info.
 ************************************************************************
 */
void delete_contexts_TextureInfo(TextureInfoContexts *deco_ctx)
{
  if( deco_ctx == NULL )
    return;

  free( deco_ctx );

  return;
}

void readFieldModeInfoFromBuffer_CABAC( SyntaxElement *se,
                                        struct inp_par *inp,
                                        struct img_par *img,
                                        DecodingEnvironmentPtr dep_dp)
{
  int a,b,act_ctx;
  MotionInfoContexts *ctx         = (img->currentSlice)->mot_ctx;
  Macroblock         *currMB      = &img->mb_data[img->map_mb_nr];//GB current_mb_nr];
  
  if (currMB->field_available[0] == NULL)
    b = 0;
  else
    b = currMB->field_available[0]->mb_field;
  if (currMB->field_available[1] == NULL)
    a = 0;
  else
    a = currMB->field_available[1]->mb_field;

  act_ctx = a + b;

  se->value1 = biari_decode_symbol (dep_dp, &ctx->mb_aff_contexts[act_ctx]);

#if TRACE
  fprintf(p_trace, "@%d %s\t\t%d\n",symbolCount++, se->tracestring, se->value1);
  fflush(p_trace);
#endif
}


int check_next_mb_and_get_field_mode_CABAC( SyntaxElement *se,
                                            struct img_par *img,
                                            struct inp_par *inp,
                                            DataPartition  *act_dp)
{
  BiContextTypePtr          mb_type_ctx_copy[4];
  BiContextTypePtr          mb_aff_ctx_copy;
  DecodingEnvironmentPtr    dep_dp_copy;

  int length;
  DecodingEnvironmentPtr    dep_dp = &(act_dp->de_cabac);

  int bframe = (img->type==B_IMG_1 || img->type==B_IMG_MULT);
  int skip   = 0;
  int field  = 0;
  int i;

  Macroblock *currMB;
  
  //get next MB
  img->current_mb_nr++;
  setRealMB_nr(img);
  currMB = &img->mb_data[img->map_mb_nr];
  currMB->slice_nr = img->current_slice_nr;
  CheckAvailabilityOfNeighborsForSkip(img);
  CheckAvailabilityOfNeighborsForAff(img);
  
  
  //create
  dep_dp_copy = (DecodingEnvironmentPtr) calloc(1, sizeof(DecodingEnvironment) );
  for (i=0;i<4;i++)
    mb_type_ctx_copy[i] = (BiContextTypePtr) calloc(NUM_MB_TYPE_CTX, sizeof(BiContextType) );
  mb_aff_ctx_copy = (BiContextTypePtr) calloc(NUM_MB_AFF_CTX, sizeof(BiContextType) );
  
  //copy
  memcpy(dep_dp_copy,dep_dp,sizeof(DecodingEnvironment));
  length = *(dep_dp_copy->Dcodestrm_len) = *(dep_dp->Dcodestrm_len);
  for (i=0;i<4;i++)
    memcpy(mb_type_ctx_copy[i], img->currentSlice->mot_ctx->mb_type_contexts[i],NUM_MB_TYPE_CTX*sizeof(BiContextType) );
  memcpy(mb_aff_ctx_copy, img->currentSlice->mot_ctx->mb_aff_contexts,NUM_MB_AFF_CTX*sizeof(BiContextType) );


  //check_next_mb
#if TRACE
  strncpy(se->tracestring, "Check MB skipflag", TRACESTRING_SIZE);
#endif
  readMB_skip_flagInfoFromBuffer_CABAC(se,inp,img,dep_dp);

  skip = (bframe)? (se->value1==0 && se->value2==0) : (se->value1==0);
  if (!skip)
  {
#if TRACE
    strncpy(se->tracestring, "Get Field mode", TRACESTRING_SIZE);
#endif
    readFieldModeInfoFromBuffer_CABAC( se,inp,img,dep_dp);
    field = se->value1;
  }

  //reset
  img->current_mb_nr--;
  setRealMB_nr(img);
  memcpy(dep_dp,dep_dp_copy,sizeof(DecodingEnvironment));
  *(dep_dp->Dcodestrm_len) = length;
  for (i=0;i<4;i++)
    memcpy(img->currentSlice->mot_ctx->mb_type_contexts[i],mb_type_ctx_copy[i], NUM_MB_TYPE_CTX*sizeof(BiContextType) );
  memcpy( img->currentSlice->mot_ctx->mb_aff_contexts,mb_aff_ctx_copy,NUM_MB_AFF_CTX*sizeof(BiContextType) );
  
  
  //delete
  free(dep_dp_copy);
  for (i=0;i<4;i++)
    free(mb_type_ctx_copy[i]);
  free(mb_aff_ctx_copy);
  
  return skip;
}




/*!
 ************************************************************************
 * \brief
 *    This function is used to arithmetically decode the motion
 *    vector data of a B-frame MB.
 ************************************************************************
 */
void readBiMVD2Buffer_CABAC( SyntaxElement *se,
                             struct inp_par *inp,
                             struct img_par *img,
                             DecodingEnvironmentPtr dep_dp)
{
  int i = img->subblock_x;
  int j = img->subblock_y;
  int a, b;
  int act_ctx;
  int act_sym;
  int mv_local_err;
  int mv_sign;
  int backward = se->value2 & 0x01;
  int k = (se->value2>>1); // MVD component

  MotionInfoContexts *ctx = img->currentSlice->mot_ctx;
  Macroblock *currMB = &img->mb_data[img->map_mb_nr];//GB current_mb_nr];


  if (j==0)
  {
    if (currMB->mb_available[0][1] == NULL)
      b = 0;
    else 
      b = absm((currMB->mb_available[0][1])->mvd[backward][BLOCK_SIZE-1][i][k]);
  }
  else
    b = absm(currMB->mvd[backward][j-1/*step_v*/][i][k]);
          
  if (i==0)
  {
    if (currMB->mb_available[1][0] == NULL)
      a = 0;
    else 
      a = absm((currMB->mb_available[1][0])->mvd[backward][j][BLOCK_SIZE-1][k]);
  }
  else
    a = absm(currMB->mvd[backward][j][i-1/*step_h*/][k]);

  if ((mv_local_err=a+b)<3)
    act_ctx = 5*k;
  else
  {
    if (mv_local_err>32)
      act_ctx=5*k+3;
    else
      act_ctx=5*k+2;
  }
  se->context = act_ctx;

  act_sym = biari_decode_symbol(dep_dp,&ctx->mv_res_contexts[0][act_ctx] );

  if (act_sym != 0)
  {
    act_ctx = 5*k+4;
    mv_sign = biari_decode_symbol_eq_prob(dep_dp);
    act_ctx=5*k;
    act_sym = unary_exp_golomb_mv_decode(dep_dp,ctx->mv_res_contexts[1]+act_ctx,3);
    act_sym++;

    if(mv_sign)
      act_sym = -act_sym;
  }
  se->value1 = act_sym;

#if TRACE
  fprintf(p_trace, "@%d      %s\t\t\t%d \n",symbolCount++, se->tracestring, se->value1);
  fflush(p_trace);
#endif
}





/*!
 ***************************************************************************
 * \brief
 *    This function is used to arithmetically encode the intra_block_modeABT
 ***************************************************************************
 */
void readABTIntraBlkModeInfo2Buffer_CABAC(SyntaxElement *se,
                                          struct inp_par *inp,
                                          struct img_par *img,
                                          DecodingEnvironmentPtr dep_dp)
{
  int                 act_sym = 0;
  int                 ftype   = (img->type==INTRA_IMG?0:1);
  MotionInfoContexts *ctx     = (img->currentSlice)->mot_ctx;

  if   (biari_decode_symbol (dep_dp, ctx->ABT_mode_contexts[ftype]+0))
  {
    if (biari_decode_symbol (dep_dp, ctx->ABT_mode_contexts[ftype]+2))   act_sym = 3;
    else                                                                 act_sym = 0;
  }
  else
  {
    if (biari_decode_symbol (dep_dp, ctx->ABT_mode_contexts[ftype]+1))   act_sym = 2;
    else                                                                 act_sym = 1;
  }

  se->value1 = act_sym;
}




/*!
 ************************************************************************
 * \brief
 *    This function is used to arithmetically decode the 8x8 block type.
 ************************************************************************
 */
void readB8_typeInfoFromBuffer_CABAC (SyntaxElement *se,
                                      struct inp_par *inp,
                                      struct img_par *img,
                                      DecodingEnvironmentPtr dep_dp)
{
  int act_sym = 0;
  int bframe  = (img->type==B_IMG_1 || img->type==B_IMG_MULT);

  MotionInfoContexts *ctx = (img->currentSlice)->mot_ctx;


  if (!bframe)
  {
    if (biari_decode_symbol (dep_dp, &ctx->b8_type_contexts[0][1]))
    {
      act_sym = 0;
    }
    else
    {
      if (biari_decode_symbol (dep_dp, &ctx->b8_type_contexts[0][2]))
      {
        act_sym = 4;
      }
      else
      {
        if (biari_decode_symbol (dep_dp, &ctx->b8_type_contexts[0][3]))
        {
          if (biari_decode_symbol (dep_dp, &ctx->b8_type_contexts[0][4])) act_sym = 2;
          else                                                            act_sym = 3;
        }
        else
        {
          act_sym = 1;
        }
      }
    }
  }
  else
  {
    if (biari_decode_symbol (dep_dp, &ctx->b8_type_contexts[1][0]))
    {
      if (biari_decode_symbol (dep_dp, &ctx->b8_type_contexts[1][1]))
      {
        if (biari_decode_symbol (dep_dp, &ctx->b8_type_contexts[1][2]))
        {
          act_sym=6;
          if (biari_decode_symbol (dep_dp, &ctx->b8_type_contexts[1][3])) act_sym+=4;
          if (biari_decode_symbol (dep_dp, &ctx->b8_type_contexts[1][3])) act_sym+=2;
          if (act_sym!=12)
          {
            if (biari_decode_symbol (dep_dp, &ctx->b8_type_contexts[1][3])) act_sym++;
          }
        }
        else
        {
          act_sym=2;
          if (biari_decode_symbol (dep_dp, &ctx->b8_type_contexts[1][3])) act_sym+=2;
          if (biari_decode_symbol (dep_dp, &ctx->b8_type_contexts[1][3])) act_sym+=1;
        }
      }
      else
      {
        if (biari_decode_symbol (dep_dp, &ctx->b8_type_contexts[1][3])) act_sym = 1;
        else                                                            act_sym = 0;
      }
      act_sym++;
    }
    else
    {
      act_sym= 0;
    }
  }
  se->value1 = act_sym;
//	if (act_sym == 13)				printf(" stop");
}

/*!
 ************************************************************************
 * \brief
 *    This function is used to arithmetically decode the macroblock
 *    type info of a given MB.
 ************************************************************************
 */
void readMB_skip_flagInfoFromBuffer_CABAC( SyntaxElement *se,
                                      struct inp_par *inp,
                                      struct img_par *img,
                                      DecodingEnvironmentPtr dep_dp)
{
  int a, b;
  int act_ctx;
  int bframe=(img->type==B_IMG_1 || img->type==B_IMG_MULT);
  MotionInfoContexts *ctx = (img->currentSlice)->mot_ctx;
  Macroblock *currMB = &img->mb_data[img->map_mb_nr];//GB current_mb_nr];


  if (bframe)
  {
    if (currMB->skip_mb_available[0][1] == NULL)
      b = 0;
    else
      b = (currMB->skip_mb_available[0][1]->mb_type==0 && currMB->skip_mb_available[0][1]->cbp==0 ? 0 : 1);
    if (currMB->skip_mb_available[1][0] == NULL)
      a = 0;
    else
      a = (currMB->skip_mb_available[1][0]->mb_type==0 && currMB->skip_mb_available[1][0]->cbp==0 ? 0 : 1);
    act_ctx = 7 + a + b;

    if (biari_decode_symbol (dep_dp, &ctx->mb_type_contexts[2][act_ctx]) == 0)
      se->value1 = se->value2 = 0;
    else
      se->value1 = se->value2 = 1;
  }
  else
  {
    if (currMB->skip_mb_available[0][1] == NULL)
      b = 0;
    else
      b = (( (currMB->skip_mb_available[0][1])->mb_type != 0) ? 1 : 0 );
    if (currMB->skip_mb_available[1][0] == NULL)
      a = 0;
    else
      a = (( (currMB->skip_mb_available[1][0])->mb_type != 0) ? 1 : 0 );
    act_ctx = a + b;

    if (biari_decode_symbol(dep_dp, &ctx->mb_type_contexts[1][act_ctx]) == 0)
      se->value1 = 0;
    else
      se->value1 = 1;
  }
#if TRACE
  fprintf(p_trace, "@%d %s\t\t%d\t%d %d\n",symbolCount++, se->tracestring, se->value1,a,b);
  fflush(p_trace);
#endif
  return;
}
/*!
 ************************************************************************
 * \brief
 *    This function is used to arithmetically decode the macroblock
 *    type info of a given MB.
 ************************************************************************
 */
void readMB_typeInfoFromBuffer_CABAC( SyntaxElement *se,
                                      struct inp_par *inp,
                                      struct img_par *img,
                                      DecodingEnvironmentPtr dep_dp)
{
  int a, b;
  int act_ctx;
  int act_sym;
  int bframe=(img->type==B_IMG_1 || img->type==B_IMG_MULT);
  int mode_sym;
  int ct = 0;
  int curr_mb_type;
  int useABT   = ((inp->abt==INTER_INTRA_ABT) || (inp->abt==INTER_ABT));


  MotionInfoContexts *ctx = (img->currentSlice)->mot_ctx;
  Macroblock *currMB = &img->mb_data[img->map_mb_nr];//GB current_mb_nr];

  if(img->type == INTRA_IMG)  // INTRA-frame
  {
    if (currMB->mb_available[0][1] == NULL)
      b = 0;
    else 
      b = (( (currMB->mb_available[0][1])->mb_type != I4MB) ? 1 : 0 );
    if (currMB->mb_available[1][0] == NULL)
      a = 0;
    else 
      a = (( (currMB->mb_available[1][0])->mb_type != I4MB) ? 1 : 0 );

    act_ctx = a + b;
    act_sym = biari_decode_symbol(dep_dp, ctx->mb_type_contexts[0] + act_ctx);
    se->context = act_ctx; // store context

    if (act_sym==0) // 4x4 Intra
    {
      curr_mb_type = act_sym;
    }
    else // 16x16 Intra
    {
      act_sym = 1;
      act_ctx = 4;
      mode_sym =  biari_decode_symbol(dep_dp, ctx->mb_type_contexts[0] + act_ctx ); // decoding of AC/no AC
      act_sym += mode_sym*12;
      act_ctx = 5;
      // decoding of cbp: 0,1,2
        mode_sym =  biari_decode_symbol(dep_dp, ctx->mb_type_contexts[0] + act_ctx );
      if (mode_sym!=0)
      {
        act_ctx=6;
        mode_sym = biari_decode_symbol(dep_dp, ctx->mb_type_contexts[0] + act_ctx );
        act_sym+=4;
        if (mode_sym!=0)
            act_sym+=4;
          }
        // decoding of I pred-mode: 0,1,2,3
        act_ctx = 7;
        mode_sym =  biari_decode_symbol(dep_dp, ctx->mb_type_contexts[0] + act_ctx );
        act_sym += mode_sym*2;
        act_ctx = 8;
        mode_sym =  biari_decode_symbol(dep_dp, ctx->mb_type_contexts[0] + act_ctx );
        act_sym += mode_sym;
        curr_mb_type = act_sym;
    }
  }
  else if(img->type == SI_IMG)  // SI-frame
  {
    // special ctx's for SI4MB
    if (currMB->mb_available[0][1] == NULL)
      b = 0;
    else 
      b = (( (currMB->mb_available[0][1])->mb_type != SI4MB) ? 1 : 0 );
    if (currMB->mb_available[1][0] == NULL)
      a = 0;
    else 
      a = (( (currMB->mb_available[1][0])->mb_type != SI4MB) ? 1 : 0 );

    act_ctx = a + b;
    act_sym = biari_decode_symbol(dep_dp, ctx->mb_type_contexts[4] + act_ctx);
    se->context = act_ctx; // store context

    if (act_sym==0) //  SI 4x4 Intra
    {
      curr_mb_type = 0;
    }
    else // analog INTRA_IMG
    {
      if (currMB->mb_available[0][1] == NULL)
        b = 0;
      else 
        b = (( (currMB->mb_available[0][1])->mb_type != I4MB) ? 1 : 0 );
      if (currMB->mb_available[1][0] == NULL)
        a = 0;
      else 
        a = (( (currMB->mb_available[1][0])->mb_type != I4MB) ? 1 : 0 );

      act_ctx = a + b;
      act_sym = biari_decode_symbol(dep_dp, ctx->mb_type_contexts[0] + act_ctx);
      se->context = act_ctx; // store context
      
      
      if (act_sym==0) // 4x4 Intra
      {
        curr_mb_type = 1;
      }
      else // 16x16 Intra
      {
        act_sym = 2;
        act_ctx = 4;
        mode_sym =  biari_decode_symbol(dep_dp, ctx->mb_type_contexts[0] + act_ctx ); // decoding of AC/no AC
        act_sym += mode_sym*12;
        act_ctx = 5;
        // decoding of cbp: 0,1,2
        mode_sym =  biari_decode_symbol(dep_dp, ctx->mb_type_contexts[0] + act_ctx );
        if (mode_sym!=0)
        {
          act_ctx=6;
          mode_sym = biari_decode_symbol(dep_dp, ctx->mb_type_contexts[0] + act_ctx );
          act_sym+=4;
          if (mode_sym!=0)
            act_sym+=4;
        }
        // decoding of I pred-mode: 0,1,2,3
        act_ctx = 7;
        mode_sym =  biari_decode_symbol(dep_dp, ctx->mb_type_contexts[0] + act_ctx );
        act_sym += mode_sym*2;
        act_ctx = 8;
        mode_sym =  biari_decode_symbol(dep_dp, ctx->mb_type_contexts[0] + act_ctx );
        act_sym += mode_sym;
        curr_mb_type = act_sym;
      }
    }
  }
  else
  {
    if (bframe)
    {
      ct = 1;
      if (currMB->mb_available[0][1] == NULL)
        b = 0;
			else
				b = (( (currMB->mb_available[0][1])->mb_type != 0) ? 1 : 0 );
			if (currMB->mb_available[1][0] == NULL)
				a = 0;
			else
				a = (( (currMB->mb_available[1][0])->mb_type != 0) ? 1 : 0 );

			act_ctx = a + b;

      if (biari_decode_symbol (dep_dp, &ctx->mb_type_contexts[2][act_ctx]))
      {
        if (biari_decode_symbol (dep_dp, &ctx->mb_type_contexts[2][4]))
        {
          if (biari_decode_symbol (dep_dp, &ctx->mb_type_contexts[2][5]))
          {
            act_sym=12;
            if (biari_decode_symbol (dep_dp, &ctx->mb_type_contexts[2][6])) act_sym+=8;
            if (biari_decode_symbol (dep_dp, &ctx->mb_type_contexts[2][6])) act_sym+=4;
            if (biari_decode_symbol (dep_dp, &ctx->mb_type_contexts[2][6])) act_sym+=2;

            if      (act_sym==24)  act_sym=11;
            else if (act_sym==26)  act_sym=22;
            else
            {
              if (act_sym==22)     act_sym=23;
              if ((!useABT) || (act_sym != 23))
                if (biari_decode_symbol (dep_dp, &ctx->mb_type_contexts[2][6])) act_sym+=1; 
            }
          }
          else
          {
            act_sym=3;
            if (biari_decode_symbol (dep_dp, &ctx->mb_type_contexts[2][6])) act_sym+=4;
            if (biari_decode_symbol (dep_dp, &ctx->mb_type_contexts[2][6])) act_sym+=2;
            if (biari_decode_symbol (dep_dp, &ctx->mb_type_contexts[2][6])) act_sym+=1;
          }
        }
        else
        {
          if (biari_decode_symbol (dep_dp, &ctx->mb_type_contexts[2][6])) act_sym=2;
          else                                                            act_sym=1;
        }
      }
      else
      {
        act_sym = 0;
      }
    }
    else // P-frame
    {
      {
        if (biari_decode_symbol(dep_dp, &ctx->mb_type_contexts[1][4] )) 
        {
          if (!useABT)
          {
            if (biari_decode_symbol(dep_dp, &ctx->mb_type_contexts[1][7] ))   act_sym = 7;
            else                                                              act_sym = 6;
          }
          else                                                                act_sym = 6;
        }
        else
        {
          if (biari_decode_symbol(dep_dp, &ctx->mb_type_contexts[1][5] ))
          {
            if (biari_decode_symbol(dep_dp, &ctx->mb_type_contexts[1][7] )) act_sym = 2;
            else                                                            act_sym = 3;
          }
          else
          {
            if (biari_decode_symbol(dep_dp, &ctx->mb_type_contexts[1][6] )) act_sym = 4;
            else                                                            act_sym = 1;
          }
        }
      }
    }

    if (act_sym<=6 || (((img->type == B_IMG_1 || img->type == B_IMG_MULT)?1:0) && act_sym<=23))
    {
      curr_mb_type = act_sym;
    }
    else  // additional info for 16x16 Intra-mode
    {
      act_ctx = 8;
      mode_sym =  biari_decode_symbol(dep_dp, ctx->mb_type_contexts[1] + act_ctx ); // decoding of AC/no AC
      act_sym += mode_sym*12;

      // decoding of cbp: 0,1,2
      act_ctx = 9;
      mode_sym = biari_decode_symbol(dep_dp, ctx->mb_type_contexts[1] + act_ctx );
      if (mode_sym != 0)
      {
        act_sym+=4;
        mode_sym = biari_decode_symbol(dep_dp, ctx->mb_type_contexts[1] + act_ctx );
        if (mode_sym != 0)
          act_sym+=4;
      }

      // decoding of I pred-mode: 0,1,2,3
      act_ctx = 10;
      mode_sym = biari_decode_symbol(dep_dp, ctx->mb_type_contexts[1] + act_ctx );
      act_sym += mode_sym*2;
      mode_sym = biari_decode_symbol(dep_dp, ctx->mb_type_contexts[1] + act_ctx );
      act_sym += mode_sym;
      curr_mb_type = act_sym;
    }
  }
  se->value1 = curr_mb_type;

//	if (curr_mb_type >= 23)				printf(" stopx");
#if TRACE
  fprintf(p_trace, "@%d%s\t\t\t%d\n",symbolCount++, se->tracestring, se->value1);
  fflush(p_trace);
#endif
}

/*!
 ************************************************************************
 * \brief
 *    This function is used to arithmetically decode a pair of
 *    intra prediction modes of a given MB.
 ************************************************************************
 */
void readIntraPredModeFromBuffer_CABAC( SyntaxElement *se,
                                        struct inp_par *inp,
                                        struct img_par *img,
                                        DecodingEnvironmentPtr dep_dp)
{
  TextureInfoContexts *ctx     = img->currentSlice->tex_ctx;

  se->value1 = (int) unary_bin_max_decode(dep_dp,ctx->ipr_contexts[0],1,8);
  se->value1--;
#if TRACE
  fprintf(p_trace, "@%d %s\t\t\t%d %d\n",symbolCount++, se->tracestring, se->value1,se->context);
  fflush(p_trace);
#endif
}

/*!
 ************************************************************************
 * \brief
 *    This function is used to arithmetically decode the reference
 *    parameter of a given MB.
 ************************************************************************
 */
void readRefFrameFromBuffer_CABAC(  SyntaxElement *se,
                                    struct inp_par *inp,
                                    struct img_par *img,
                                    DecodingEnvironmentPtr dep_dp)
{
  MotionInfoContexts *ctx = img->currentSlice->mot_ctx;
  Macroblock *currMB = &img->mb_data[img->map_mb_nr];//GB current_mb_nr];

  int   addctx = se->context;
  int   a, b;
  int   act_ctx;
  int   act_sym;
  int** refframe_array = ((img->type==B_IMG_1 || img->type==B_IMG_MULT) ? img->fw_refFrArr : refFrArr);
  int   block_y        = img->block_y;

  if( img->mb_frame_field_flag )
  {
    if( !img->mb_field )
    {
      refframe_array = ((img->type==B_IMG_1 || img->type==B_IMG_MULT) ? img->fw_refFrArr_frm : refFrArr_frm);
    }
    else if ( img->current_mb_nr % 2 )
    {
      refframe_array = ((img->type==B_IMG_1 || img->type==B_IMG_MULT) ? img->fw_refFrArr_bot : refFrArr_bot);
      block_y        = ( img->block_y - 4 ) / 2;
    }
    else
    {
      refframe_array = ((img->type==B_IMG_1 || img->type==B_IMG_MULT) ? img->fw_refFrArr_top : refFrArr_top);
      block_y        = img->block_y / 2;
    }
  }

  if (currMB->mb_available[0][1] == NULL)
    b = 0;
  else
    b = (refframe_array[block_y+img->subblock_y-1][img->block_x+img->subblock_x] > 0 ? 1 : 0);
  if (currMB->mb_available[1][0] == NULL)
    a = 0;
  else 
    a = (refframe_array[block_y+img->subblock_y][img->block_x+img->subblock_x-1] > 0 ? 1 : 0);

  act_ctx = a + 2*b;
  se->context = act_ctx; // store context

  act_sym = biari_decode_symbol(dep_dp,ctx->ref_no_contexts[addctx] + act_ctx );

  if (act_sym != 0)
  {
    act_ctx = 4;
    act_sym = unary_bin_decode(dep_dp,ctx->ref_no_contexts[addctx]+act_ctx,1);
    act_sym++;
  }
  se->value1 = act_sym;

#if TRACE
  fprintf(p_trace, "@%d %s\t\t\t%d \n",symbolCount++, se->tracestring, se->value1);
//  fprintf(p_trace," c: %d :%d \n",ctx->ref_no_contexts[addctx][act_ctx].cum_freq[0],ctx->ref_no_contexts[addctx][act_ctx].cum_freq[1]);
  fflush(p_trace);
  // commented out, does not compile. karll@real.com
#endif
}

/*!
 ************************************************************************
 * \brief
 *    This function is used to arithmetically decode the backward reference
 *    parameter of a given MB.
 ************************************************************************
 */
void readBwdRefFrameFromBuffer_CABAC(  SyntaxElement *se,
                                    struct inp_par *inp,
                                    struct img_par *img,
                                    DecodingEnvironmentPtr dep_dp)
{
  MotionInfoContexts *ctx = img->currentSlice->mot_ctx;
  Macroblock *currMB = &img->mb_data[img->map_mb_nr];//GB current_mb_nr];

  int   addctx = se->context;
  int   a, b;
  int   act_ctx;
  int   act_sym;
  int** refframe_array = img->bw_refFrArr;
  int   block_y        = img->block_y;
 
  if( img->mb_frame_field_flag )
  {
    if( !img->mb_field )
    {
      refframe_array = img->bw_refFrArr_frm;
    }
    else if ( img->current_mb_nr % 2 )
    {
      refframe_array = img->bw_refFrArr_bot;
      block_y        = ( img->block_y - 4 ) / 2;
    }
    else
    {
      refframe_array = img->bw_refFrArr_top;
      block_y        = img->block_y / 2;
    }
  }

#define REF_IDX(r) ((img->num_ref_pic_active_bwd>1 && (r)<2) ? (1-(r)):(r))

  if (currMB->mb_available[0][1] == NULL)
    b = 0;
  else
    b = (REF_IDX(refframe_array[block_y+img->subblock_y-1][img->block_x+img->subblock_x]) > 0 ? 1 : 0);
  if (currMB->mb_available[1][0] == NULL)
    a = 0;
  else 
    a = (REF_IDX(refframe_array[block_y+img->subblock_y][img->block_x+img->subblock_x-1]) > 0 ? 1 : 0);
#undef REF_IDX

  act_ctx = a + 2*b;
  se->context = act_ctx; // store context

  act_sym = biari_decode_symbol(dep_dp,ctx->bwd_ref_no_contexts[addctx] + act_ctx );

  if (act_sym != 0)
  {
    act_ctx = 4;
    act_sym = unary_bin_decode(dep_dp,ctx->bwd_ref_no_contexts[addctx]+act_ctx,1);
    act_sym++;
  }
  se->value1 = act_sym;

#if TRACE
//  fprintf(p_trace, "@%d%s\t\t\t%d",symbolCount++, se->tracestring, se->value1);
//  fprintf(p_trace," c: %d :%d \n",ctx->bwd_ref_no_contexts[addctx][act_ctx].cum_freq[0],ctx->bwd_ref_no_contexts[addctx][act_ctx].cum_freq[1]);
//  fflush(p_trace);
  // commented out, does not compile. karll@real.com
#endif
}


/*!
 ************************************************************************
 * \brief
 *    This function is used to arithmetically decode the motion
 *    vector data of a given MB.
 ************************************************************************
 */
void readMVDFromBuffer_CABAC(SyntaxElement *se,
                             struct inp_par *inp,
                             struct img_par *img,
                             DecodingEnvironmentPtr dep_dp)
{
  int i = img->subblock_x;
  int j = img->subblock_y;
  int a, b;
  int act_ctx;
  int act_sym;
  int mv_pred_res;
  int mv_local_err;
  int mv_sign;
  int k = se->value2; // MVD component

  MotionInfoContexts *ctx = img->currentSlice->mot_ctx;
  Macroblock *currMB = &img->mb_data[img->map_mb_nr];//GB current_mb_nr];

  if (j==0)
  {
    if (currMB->mb_available[0][1] == NULL)
      b = 0;
    else 
      b = absm((currMB->mb_available[0][1])->mvd[0][BLOCK_SIZE-1][i][k]);
  }
  else
    b = absm(currMB->mvd[0][j-1/*step_v*/][i][k]);
          
  if (i==0)
  {
    if (currMB->mb_available[1][0] == NULL)
      a = 0;
    else 
      a = absm((currMB->mb_available[1][0])->mvd[0][j][BLOCK_SIZE-1][k]);
  }
  else
    a = absm(currMB->mvd[0][j][i-1/*step_h*/][k]);

  if ((mv_local_err=a+b)<3)
    act_ctx = 5*k;
  else
  {
    if (mv_local_err>32)
      act_ctx=5*k+3;
    else
      act_ctx=5*k+2;
  }
  se->context = act_ctx;

  act_sym = biari_decode_symbol(dep_dp, &ctx->mv_res_contexts[0][act_ctx] );

  if (act_sym == 0)
  {
    mv_pred_res = 0;
  }
  else
  {
    act_ctx=5*k+4;
    mv_sign = biari_decode_symbol_eq_prob(dep_dp);
    act_ctx=5*k;
    act_sym = unary_exp_golomb_mv_decode(dep_dp,ctx->mv_res_contexts[1]+act_ctx,3);
    act_sym++;
    mv_pred_res = ((mv_sign != 0) ? (-act_sym) : act_sym);
  }
  se->value1 = mv_pred_res;

#if TRACE
  fprintf(p_trace, "@%d %s\t\t\t%d %d\n",symbolCount++, se->tracestring, se->value1,a+b);
  fflush(p_trace);
#endif
}

/*!
 ************************************************************************
 * \brief
 *    This function is used to arithmetically decode the delta qp
 *     of a given MB.
 ************************************************************************
 */
void readDquant_FromBuffer_CABAC(SyntaxElement *se,
                                struct inp_par *inp,
                                struct img_par *img,
                                DecodingEnvironmentPtr dep_dp)
{
  MotionInfoContexts *ctx = img->currentSlice->mot_ctx;
  Macroblock *currMB = &img->mb_data[img->map_mb_nr];//GB current_mb_nr];

  int act_ctx;
  int act_sym;
  int dquant;

  if (currMB->mb_available[1][0] == NULL)
    act_ctx = 0;
  else
    act_ctx = ( ((currMB->mb_available[1][0])->delta_quant != 0) ? 1 : 0);

  act_sym = biari_decode_symbol(dep_dp,ctx->delta_qp_contexts + act_ctx );
  if (act_sym != 0)
  {
    act_ctx = 2;
    act_sym = unary_bin_decode(dep_dp,ctx->delta_qp_contexts+act_ctx,1);
    act_sym++;
  }

  dquant = (act_sym+1)/2;
  if((act_sym & 0x01)==0)                           // lsb is signed bit
    dquant = -dquant;
  se->value1 = dquant;

#if TRACE
  fprintf(p_trace, "@%d%s\t\t\t%d\n",symbolCount++, se->tracestring, se->value1);
  fflush(p_trace);
#endif
}
/*!
 ************************************************************************
 * \brief
 *    This function is used to arithmetically decode the coded
 *    block pattern of a given MB.
 ************************************************************************
 */
void readCBPFromBuffer_CABAC(SyntaxElement *se,
                             struct inp_par *inp,
                             struct img_par *img,
                             DecodingEnvironmentPtr dep_dp)
{
  TextureInfoContexts *ctx = img->currentSlice->tex_ctx;
  Macroblock *currMB = &img->mb_data[img->map_mb_nr];//GB current_mb_nr];

  int mb_x, mb_y;
  int a, b;
  int curr_cbp_ctx, curr_cbp_idx;
  int cbp = 0;
  int cbp_bit;
  int mask;

  //  coding of luma part (bit by bit)
  for (mb_y=0; mb_y < 4; mb_y += 2)
  {
    for (mb_x=0; mb_x < 4; mb_x += 2)
    {
      if (currMB->b8mode[mb_y+(mb_x/2)]==IBLOCK)
        curr_cbp_idx = 0;
      else
        curr_cbp_idx = 1;

      if (mb_y == 0)
      {
        if (currMB->mb_available[0][1] == NULL)
          b = 0;
        else
          b = (( ((currMB->mb_available[0][1])->cbp & (1<<(2+mb_x/2))) == 0) ? 1 : 0);
      }
      else
        b = ( ((cbp & (1<<(mb_x/2))) == 0) ? 1: 0);

      if (mb_x == 0)
      {
        if (currMB->mb_available[1][0] == NULL)
          a = 0;
        else
          a = (( ((currMB->mb_available[1][0])->cbp & (1<<(mb_y+1))) == 0) ? 1 : 0);
      }
      else
        a = ( ((cbp & (1<<mb_y)) == 0) ? 1: 0);

      curr_cbp_ctx = a+2*b;
      mask = (1<<(mb_y+mb_x/2));
      cbp_bit = biari_decode_symbol(dep_dp, ctx->cbp_contexts[0] + curr_cbp_ctx );
      if (cbp_bit) cbp += mask;
    }
  }


  if ( se->type == SE_CBP_INTRA )
    curr_cbp_idx = 0;
  else
    curr_cbp_idx = 1;

  // coding of chroma part
  b = 0;
  if (currMB->mb_available[0][1] != NULL)
    b = ((currMB->mb_available[0][1])->cbp > 15) ? 1 : 0;

  a = 0;
  if (currMB->mb_available[1][0] != NULL)
    a = ((currMB->mb_available[1][0])->cbp > 15) ? 1 : 0;

  curr_cbp_ctx = a+2*b;
  cbp_bit = biari_decode_symbol(dep_dp, ctx->cbp_contexts[1] + curr_cbp_ctx );

  if (cbp_bit) // set the chroma bits
  {
    b = 0;
    if (currMB->mb_available[0][1] != NULL)
      if ((currMB->mb_available[0][1])->cbp > 15)
        b = (( ((currMB->mb_available[0][1])->cbp >> 4) == 2) ? 1 : 0);

    a = 0;
    if (currMB->mb_available[1][0] != NULL)
      if ((currMB->mb_available[1][0])->cbp > 15)
        a = (( ((currMB->mb_available[1][0])->cbp >> 4) == 2) ? 1 : 0);

    curr_cbp_ctx = a+2*b;
    cbp_bit = biari_decode_symbol(dep_dp, ctx->cbp_contexts[2] + curr_cbp_ctx );
    cbp += (cbp_bit == 1) ? 32 : 16;
  }

  se->value1 = cbp;

#if TRACE
  fprintf(p_trace, "@%d      %s\t\t\t%d\n",symbolCount++, se->tracestring, se->value1);
  fflush(p_trace);
#endif
}

/*!
 ************************************************************************
 * \brief
 *    This function is used to arithmetically decode the chroma
 *    intra prediction mode of a given MB.
 ************************************************************************
 */  //GB
void readCIPredMode_FromBuffer_CABAC(SyntaxElement *se,
                                     struct inp_par *inp,
                                     struct img_par *img,
                                     DecodingEnvironmentPtr dep_dp)
{

  TextureInfoContexts *ctx = img->currentSlice->tex_ctx;
  Macroblock          *currMB  = &img->mb_data[img->map_mb_nr];
  int                 act_ctx,a,b;
  int                 act_sym  = se->value1;

  if (currMB->mb_available[0][1] == NULL) b = 0;
  else  b = ( ((currMB->mb_available[0][1])->c_ipred_mode != 0) ? 1 : 0);

  if (currMB->mb_available[1][0] == NULL) a = 0;
  else  a = ( ((currMB->mb_available[1][0])->c_ipred_mode != 0) ? 1 : 0);

  act_ctx = a+b;

  act_sym = biari_decode_symbol(dep_dp, ctx->cipr_contexts + act_ctx );

  if (act_sym!=0) 
    act_sym = unary_bin_max_decode(dep_dp,ctx->cipr_contexts+3,0,2)+1;

//	act_sym = 0;
//	act_sym |= biari_decode_symbol_eq_prob(dep_dp);
//	act_sym |= (biari_decode_symbol_eq_prob(dep_dp)<<1);

  se->value1 = act_sym;

	//if (img->type==B_IMG_1 || img->type==B_IMG_MULT) 
	//	printf(" %d",act_sym);

#if TRACE
  fprintf(p_trace, "@%d%s\t\t\t%d\n",symbolCount++, se->tracestring, se->value1);
  fflush(p_trace);
#endif

}


static const int maxpos       [] = {16, 15, 64, 32, 32, 16,  4, 15};
static const int c1isdc       [] = { 1,  0,  1,  1,  1,  1,  1,  0};

static const int type2ctx_bcbp[] = { 0,  1,  2,  2,  3,  4,  5,  6}; // 7
static const int type2ctx_map [] = { 0,  1,  2,  3,  4,  5,  6,  7}; // 8
static const int type2ctx_last[] = { 0,  1,  2,  3,  4,  5,  6,  7}; // 8
static const int type2ctx_one [] = { 0,  1,  2,  3,  3,  4,  5,  6}; // 7
static const int type2ctx_abs [] = { 0,  1,  2,  3,  3,  4,  5,  6}; // 7




/*!
 ************************************************************************
 * \brief
 *    Read CBP4-BIT
 ************************************************************************
 */
int read_and_store_CBP_block_bit (Macroblock              *currMB,
                                  DecodingEnvironmentPtr  dep_dp,
                                  struct img_par          *img,
                                  int                     type)
{
#define BIT_SET(x,n)  ((int)(((x)&(1<<(n)))>>(n)))

  int y_ac        = (type==LUMA_16AC || type==LUMA_8x8 || type==LUMA_8x4 || type==LUMA_4x8 || type==LUMA_4x4);
  int y_dc        = (type==LUMA_16DC);
  int u_ac        = (type==CHROMA_AC && !img->is_v_block);
  int v_ac        = (type==CHROMA_AC &&  img->is_v_block);
  int u_dc        = (type==CHROMA_DC && !img->is_v_block);
  int v_dc        = (type==CHROMA_DC &&  img->is_v_block);
  int j           = (y_ac || u_ac || v_ac ? img->subblock_y : 0);
  int i           = (y_ac || u_ac || v_ac ? img->subblock_x : 0);
  int bit         = (y_dc ? 0 : y_ac ? 1+4*j+i : u_dc ? 17 : v_dc ? 18 : u_ac ? 19+2*j+i : 23+2*j+i);
  int ystep_back  = (y_ac ? 12 : u_ac || v_ac ? 2 : 0);
  int xstep_back  = (y_ac ?  3 : u_ac || v_ac ? 1 : 0);
  int ystep       = (y_ac ?  4 : u_ac || v_ac ? 2 : 0);
  int default_bit = (img->is_intra_block ? 1 : 0);
  int upper_bit   = default_bit;
  int left_bit    = default_bit;
  int cbp_bit     = 1;  // always one for 8x8 mode
  int ctx;


  if (type!=LUMA_8x8)
  {
    //--- get bits from neighbouring blocks ---
    if (j==0)
    {
      if (currMB->mb_available[0][1])
      {
        upper_bit = BIT_SET(currMB->mb_available[0][1]->cbp_bits,bit+ystep_back);
      }
    }
    else
    {
      upper_bit = BIT_SET(currMB->cbp_bits,bit-ystep);
    }
    if (i==0)
    {
      if (currMB->mb_available[1][0])
      {
        left_bit = BIT_SET(currMB->mb_available[1][0]->cbp_bits,bit+xstep_back);
      }
    }
    else
    {
      left_bit = BIT_SET(currMB->cbp_bits,bit-1);
    }
    ctx = 2*upper_bit+left_bit;


    //===== encode symbol =====
    cbp_bit = biari_decode_symbol (dep_dp, img->currentSlice->tex_ctx->bcbp_contexts[type2ctx_bcbp[type]] + ctx);
  }
  
  //--- set bits for current block ---
  if (cbp_bit)
  {
    if (type==LUMA_8x8)
    {
      currMB->cbp_bits   |= (1<< bit   );
      currMB->cbp_bits   |= (1<<(bit+1));
      currMB->cbp_bits   |= (1<<(bit+4));
      currMB->cbp_bits   |= (1<<(bit+5));
    }
    else if (type==LUMA_8x4)
    {
      currMB->cbp_bits   |= (1<< bit   );
      currMB->cbp_bits   |= (1<<(bit+1));
    }
    else if (type==LUMA_4x8)
    {
      currMB->cbp_bits   |= (1<< bit   );
      currMB->cbp_bits   |= (1<<(bit+4));
    }
    else
    {
      currMB->cbp_bits   |= (1<<bit);
    }
  }

  return cbp_bit;
}





//===== position -> ctx for MAP =====
//--- zig-zag scan ----
static const int  pos2ctx_map8x8 [] = { 0,  1,  2,  3,  4,  5,  5,  4,  4,  3,  3,  4,  4,  4,  5,  5,
                                        4,  4,  4,  4,  3,  3,  6,  7,  7,  7,  8,  9, 10,  9,  8,  7,
                                        7,  6, 11, 12, 13, 11,  6,  7,  8,  9, 14, 10,  9,  8,  6, 11,
                                       12, 13, 11,  6,  9, 14, 10,  9, 11, 12, 13, 11 ,14, 10, 12, 14}; // 15 CTX
static const int  pos2ctx_map8x4 [] = { 0,  1,  2,  3,  4,  5,  7,  8,  9, 10, 11,  9,  8,  6,  7,  8,
                                        9, 10, 11,  9,  8,  6, 12,  8,  9, 10, 11,  9, 13, 13, 14, 14}; // 15 CTX
static const int  pos2ctx_map4x4 [] = { 0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 14}; // 15 CTX
static const int* pos2ctx_map    [] = {pos2ctx_map4x4, pos2ctx_map4x4, pos2ctx_map8x8, pos2ctx_map8x4,
                                       pos2ctx_map8x4, pos2ctx_map4x4, pos2ctx_map4x4, pos2ctx_map4x4};
//--- interlace scan ----
static const int  pos2ctx_map8x8i[] = { 0,  1,  1,  2,  2,  3,  3,  4,  5,  6,  7,  7,  7,  8,  4,  5,
                                        6,  9, 10, 10,  8, 11, 12, 11,  9,  9, 10, 10,  8, 11, 12, 11,
                                        9,  9, 10, 10,  8, 11, 12, 11,  9,  9, 10, 10,  8, 13, 13,  9,
                                        9, 10, 10,  8, 13, 13,  9,  9, 10, 10, 14, 14, 14, 14, 14, 14}; // 15 CTX
static const int  pos2ctx_map8x4i[] = { 0,  1,  2,  3,  4,  5,  6,  3,  4,  5,  6,  3,  4,  7,  6,  8,
                                        9,  7,  6,  8,  9, 10, 11, 12, 12, 10, 11, 13, 13, 14, 14, 14}; // 15 CTX
static const int  pos2ctx_map4x8i[] = { 0,  1,  1,  1,  2,  3,  3,  4,  4,  4,  5,  6,  2,  7,  7,  8,
                                        8,  8,  5,  6,  9, 10, 10, 11, 11, 11, 12, 13, 13, 14, 14, 14}; // 15 CTX
static const int* pos2ctx_map_int[] = {pos2ctx_map4x4, pos2ctx_map4x4, pos2ctx_map8x8i,pos2ctx_map8x4i,
                                       pos2ctx_map4x8i,pos2ctx_map4x4, pos2ctx_map4x4, pos2ctx_map4x4};


//===== position -> ctx for LAST =====
static const int  pos2ctx_last8x8 [] = { 0,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,
                                         2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,
                                         3,  3,  3,  3,  3,  3,  3,  3,  4,  4,  4,  4,  4,  4,  4,  4,
                                         5,  5,  5,  5,  6,  6,  6,  6,  7,  7,  7,  7,  8,  8,  8,  8}; //  9 CTX
static const int  pos2ctx_last8x4 [] = { 0,  1,  1,  1,  1,  1,  1,  1,  2,  2,  2,  2,  2,  2,  2,  2,
                                         3,  3,  3,  3,  4,  4,  4,  4,  5,  5,  6,  6,  7,  7,  8,  8}; //  9 CTX

static const int  pos2ctx_last4x4 [] = { 0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15}; // 15 CTX
static const int* pos2ctx_last    [] = {pos2ctx_last4x4, pos2ctx_last4x4, pos2ctx_last8x8, pos2ctx_last8x4,
                                        pos2ctx_last8x4, pos2ctx_last4x4, pos2ctx_last4x4, pos2ctx_last4x4};





/*!
 ************************************************************************
 * \brief
 *    Read Significance MAP
 ************************************************************************
 */
int read_significance_map (Macroblock              *currMB,
                           DecodingEnvironmentPtr  dep_dp,
                           struct img_par          *img,
                           int                     type,
                           int                     coeff[])
{
  int   i, sig;
  int   coeff_ctr = 0;
  int   i0        = 0;
  int   i1        = maxpos[type]-1;

  if (!c1isdc[type])
  {
    i0++; i1++; coeff--;
  }

  for (i=i0; i<i1; i++) // if last coeff is reached, it has to be significant
  {
    //--- read significance symbol ---
    if (img->structure!=FRAME)
      sig = biari_decode_symbol   (dep_dp, img->currentSlice->tex_ctx->map_contexts [type2ctx_map [type]] + pos2ctx_map_int [type][i]);
    else
      sig = biari_decode_symbol   (dep_dp, img->currentSlice->tex_ctx->map_contexts [type2ctx_map [type]] + pos2ctx_map     [type][i]);
    if (sig)
    {
      coeff[i] = 1;
      coeff_ctr++;
      //--- read last coefficient symbol ---
      if (biari_decode_symbol (dep_dp, img->currentSlice->tex_ctx->last_contexts[type2ctx_last[type]] + pos2ctx_last[type][i]))
      {
        for (i++; i<i1+1; i++) coeff[i] = 0;
      }
    }
    else
    {
      coeff[i] = 0;
    }
  }
  //--- last coefficient must be significant if no last symbol was received ---
  if (i<i1+1)
  {
    coeff[i] = 1;
    coeff_ctr++;
  }

  return coeff_ctr;
}



/*!
 ************************************************************************
 * \brief
 *    Read Levels
 ************************************************************************
 */
void read_significant_coefficients (Macroblock              *currMB,
                                    DecodingEnvironmentPtr  dep_dp,
                                    struct img_par          *img,
                                    int                     type,
                                    int                     coeff[])
{
  int   i, ctx;
  int   c1 = 1;
  int   c2 = 0;

  for (i=maxpos[type]-1; i>=0; i--)
  {
    if (coeff[i]!=0)
    {
      ctx = min (c1,4);
      coeff[i] += biari_decode_symbol (dep_dp, img->currentSlice->tex_ctx->one_contexts[type2ctx_one[type]] + ctx);
      if (coeff[i]==2)
      {
        ctx = min (c2,4);
        coeff[i] += unary_exp_golomb_level_decode (dep_dp, img->currentSlice->tex_ctx->abs_contexts[type2ctx_abs[type]]+ctx);
        c1=0;
        c2++;
      }
      else if (c1)
      {
        c1++;
      }
      if (biari_decode_symbol_eq_prob(dep_dp))
      {
        coeff[i] *= -1;
      }
    }
  }
}


/*!
 ************************************************************************
 * \brief
 *    Read Block-Transform Coefficients
 ************************************************************************
 */
void readRunLevelFromBuffer_CABAC (SyntaxElement  *se,
                                   struct inp_par *inp,
                                   struct img_par *img,    
                                   DecodingEnvironmentPtr dep_dp)
{
  static int  coeff[64]; // one more for EOB
  static int  coeff_ctr = -1;
  static int  pos       =  0;

  Macroblock *currMB = &img->mb_data[img->map_mb_nr];//GB current_mb_nr];

  //--- read coefficients for whole block ---
  if (coeff_ctr < 0)
  {
    //===== decode CBP-BIT =====
    if ((coeff_ctr = read_and_store_CBP_block_bit (currMB, dep_dp, img, se->context)))
    {
      //===== decode significance map =====
      coeff_ctr = read_significance_map (currMB, dep_dp, img, se->context, coeff);

      //===== decode significant coefficients =====
      read_significant_coefficients     (currMB, dep_dp, img, se->context, coeff);
    }
  }

  //--- set run and level ---
  if (coeff_ctr)
  {
    //--- set run and level (coefficient) ---
    for (se->value2=0; coeff[pos]==0; pos++, se->value2++);
    se->value1=coeff[pos++];
  }
  else
  {
    //--- set run and level (EOB) ---
    se->value1 = se->value2 = 0;
  }
  //--- decrement coefficient counter and re-set position ---
  if (coeff_ctr-- == 0) pos=0;

#if TRACE
  fprintf(p_trace, "@%d%s\t\t\t%d\t%d\n",symbolCount++, se->tracestring, se->value1,se->value2);
  fflush(p_trace);
#endif
}



/*!
 ************************************************************************
 * \brief
 *    arithmetic decoding
 ************************************************************************
 */
int readSyntaxElement_CABAC(SyntaxElement *se, struct img_par *img, struct inp_par *inp, DataPartition *this_dataPart)
{
  int curr_len;
  DecodingEnvironmentPtr dep_dp = &(this_dataPart->de_cabac);

  curr_len = arideco_bits_read(dep_dp);

  // perform the actual decoding by calling the appropriate method
  se->reading(se, inp, img, dep_dp);

  return (se->len = (arideco_bits_read(dep_dp) - curr_len));
}

/*!
 ************************************************************************
 * \brief
 *    get slice and header
 ************************************************************************
 */
int readSliceCABAC(struct img_par *img, struct inp_par *inp)
{
  Slice *currSlice = img->currentSlice;
  Bitstream *currStream = currSlice->partArr[0].bitstream;
  unsigned char *code_buffer = currStream->streamBuffer;
  int *read_len = &(currStream->read_len);
  DecodingEnvironmentPtr dep = &((currSlice->partArr[0]).de_cabac);
  int current_header;
  int BitstreamLengthInBytes;
  int BitsUsedByHeader = 0, ByteStartPosition;
  int newframe = 0;   //WYK: Oct. 8, 2001, change the method to find a new frame
  int startcodeprefix_len; //Number of bytes taken by start code prefix

  currStream->frame_bitoffset =0;

  memset (code_buffer, 0xff, MAX_CODED_FRAME_SIZE);   // this prevents a buffer full with zeros
  BitstreamLengthInBytes = currStream->bitstream_length = GetOneSliceIntoSourceBitBuffer(img, inp, code_buffer, &startcodeprefix_len);

  // Here we are ready to interpret the picture and slice headers.  Since
  // SliceHeader() gets the data out of the UVLC's len/info
  // array, we need to convert the start of our slice to such a format.


  if (BitstreamLengthInBytes < startcodeprefix_len)
    return EOS;

  BitstreamLengthInBytes = currStream->bitstream_length = EBSPtoRBSP(code_buffer, currStream->bitstream_length, startcodeprefix_len);
  BitstreamLengthInBytes = currStream->bitstream_length = RBSPtoSODB(code_buffer, currStream->bitstream_length);

  // Now we have the bits between the current startcode (inclusive) and the
  // next start code in code_buffer.  Now decode the start codes and the headers
  currStream->frame_bitoffset += startcodeprefix_len * 8;
  BitsUsedByHeader+=SliceHeader(img, inp);

  //WYK: Oct. 8, 2001, change the method to find a new frame
  if(img->tr != img->tr_old)
    newframe = 1;
  else 
    newframe = 0;
// printf ("readSlice_CABAC: tr_old %d, tr %d, startcodeprefix_len %d, newframe %d\n", img->tr_old, img->tr, startcodeprefix_len, newframe);    
  img->tr_old = img->tr;
  // if the TR of current slice is not identical to the TR of previous received slice, we have a new frame
  if(newframe)
    current_header = SOP;
  else
    current_header = SOS;

  ByteStartPosition = currStream->frame_bitoffset/8;
  if ((currStream->frame_bitoffset)%8 != 0)
    ByteStartPosition++;
  arideco_start_decoding(dep, code_buffer, ByteStartPosition, read_len, img->type);

  currSlice->picture_id = img->tr;
  return current_header;

}

/*!
 ************************************************************************
 * \brief
 *    decoding of unary binarization using one or 2 distinct
 *    models for the first and all remaining bins; no terminating
 *    "0" for max_symbol
 ***********************************************************************
 */
unsigned int unary_bin_max_decode(DecodingEnvironmentPtr dep_dp,
                                  BiContextTypePtr ctx,
                                  int ctx_offset,
                                  unsigned int max_symbol)
{
  unsigned int l;
  unsigned int symbol;
  BiContextTypePtr ictx;

  symbol =  biari_decode_symbol(dep_dp, ctx );

  if (symbol==0)
    return 0;
  else
  {
    if (max_symbol == 1)
    return symbol;
    symbol=0;
    ictx=ctx+ctx_offset;
    do
    {
      l=biari_decode_symbol(dep_dp, ictx);
      symbol++;
    }
    while( (l!=0) && (symbol<max_symbol-1) );
    if ((l!=0) && (symbol==max_symbol-1))
      symbol++;
    return symbol;
  }
}


/*!
 ************************************************************************
 * \brief
 *    decoding of unary binarization using one or 2 distinct
 *    models for the first and all remaining bins
 ***********************************************************************
 */
unsigned int unary_bin_decode(DecodingEnvironmentPtr dep_dp,
                              BiContextTypePtr ctx,
                              int ctx_offset)
{
  unsigned int l;
  unsigned int symbol;
  BiContextTypePtr ictx;

  symbol = biari_decode_symbol(dep_dp, ctx );

  if (symbol==0)
    return 0;
  else
  {
    symbol=0;
    ictx=ctx+ctx_offset;
    do
    {
      l=biari_decode_symbol(dep_dp, ictx);
      symbol++;
    }
    while( l!=0 );
    return symbol;
  }
}


/*!
 ************************************************************************
 * \brief
 *    finding end of a slice in case this is not the end of a frame
 *
 * Unsure whether the "correction" below actually solves an off-by-one
 * problem or whether it introduces one in some cases :-(  Anyway,
 * with this change the bit stream format works with CABAC again.
 * StW, 8.7.02
 ************************************************************************
 */
int cabac_startcode_follows(struct img_par *img, struct inp_par *inp, int eos_bit)
{
  Slice         *currSlice  = img->currentSlice;
  int           *partMap    = assignSE2partition[currSlice->dp_mode];
  DataPartition *dP;
  unsigned int  bit;
  DecodingEnvironmentPtr dep_dp;
  
  if(img->type == B_IMG_1 || img->type == B_IMG_MULT) dP = &(currSlice->partArr[partMap[SE_BFRAME]]);
  else                                                dP = &(currSlice->partArr[partMap[SE_MBTYPE]]);
  dep_dp = &(dP->de_cabac);

  if( eos_bit )
  {
    bit = biari_decode_symbol (dep_dp, &(currSlice->mot_ctx->slice_term_context));
#if TRACE
//	strncpy(se->tracestring, "Decode Sliceterm", TRACESTRING_SIZE);
  fprintf(p_trace, "@%d %s\t\t%d\n",symbolCount++, "Decode Sliceterm", bit);
  fflush(p_trace);
#endif
  }
  else
  {
    bit = 0;
  }

  return (bit==1?1:0);
}




/*!
 ************************************************************************
 * \brief
 *    Exp Golomb binarization and decoding of a symbol without model
 *    k = exp golombparameter (k==0, uvlc)
 ************************************************************************
 */
unsigned int exp_golomb_decode( DecodingEnvironmentPtr dep_dp,
                                int k)
{
  unsigned int l;
  int symbol = 0;
  int binary_symbol = 0;

  do
  {
    l=biari_decode_symbol_eq_prob(dep_dp);
    if (l==1) 
    {
      symbol += (1<<k); 
      k++;
    }
  }
  while (l!=0);

  while (k--)                             //next binary part
    if (biari_decode_symbol_eq_prob(dep_dp)==1) 
      binary_symbol |= (1<<k);

  return (unsigned int) (symbol+binary_symbol);
}


/*!
 ************************************************************************
 * \brief
 *    Exp Golomb binarization and decoding of a symbol
 *    with prob. of 0.5
 ************************************************************************
 */
unsigned int exp_golomb_decode_eq_prob( DecodingEnvironmentPtr dep_dp,
                                        int k)
{
  unsigned int l;
  int symbol = 0;
  int binary_symbol = 0;

  do
  {
    l=biari_decode_symbol_eq_prob(dep_dp);
    if (l==1) 
    {
      symbol += (1<<k); 
      k++;
    }
  }
  while (l!=0);

  while (k--)                             //next binary part
    if (biari_decode_symbol_eq_prob(dep_dp)==1) 
      binary_symbol |= (1<<k);

  return (unsigned int) (symbol+binary_symbol);
}


/*!
 ************************************************************************
 * \brief
 *    Exp-Golomb decoding for LEVELS
 ***********************************************************************
 */
unsigned int unary_exp_golomb_level_decode( DecodingEnvironmentPtr dep_dp,
                                            BiContextTypePtr ctx)
{
  unsigned int l,k;
  unsigned int symbol;
  unsigned int exp_start = 13;

  symbol = biari_decode_symbol(dep_dp, ctx );

  if (symbol==0)
    return 0;
  else
  {
    symbol=0;
    k=1;
    do
    {
      l=biari_decode_symbol(dep_dp, ctx);
      symbol++;
      k++;
    }
    while((l!=0) && (k!=exp_start));
    if (l!=0)
      symbol += exp_golomb_decode_eq_prob(dep_dp,0)+1;
    return symbol;
  }
}




/*!
 ************************************************************************
 * \brief
 *    Exp-Golomb decoding for Motion Vectors
 ***********************************************************************
 */
unsigned int unary_exp_golomb_mv_decode(DecodingEnvironmentPtr dep_dp,
                                        BiContextTypePtr ctx,
                                        unsigned int max_bin)
{
  unsigned int l,k;
  unsigned int bin=1;
  unsigned int symbol;
  unsigned int exp_start = 8;

  BiContextTypePtr ictx=ctx;

  symbol = biari_decode_symbol(dep_dp, ictx );

  if (symbol==0)
    return 0;
  else
  {
    symbol=0;
    k=1;

    ictx++;
    do
    {
      l=biari_decode_symbol(dep_dp, ictx  );
      if ((++bin)==2) ictx++;
      if (bin==max_bin) ictx++;
      symbol++;
      k++;
    }
    while((l!=0) && (k!=exp_start));
    if (l!=0)
      symbol += exp_golomb_decode_eq_prob(dep_dp,3)+1;
    return symbol;
  }
}
