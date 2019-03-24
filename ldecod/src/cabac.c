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
#include <string.h>
#include "cabac.h"
#include "memalloc.h"
#include "elements.h"
#include "bitsbuf.h"
#include "header.h"

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
  int j;
  MotionInfoContexts *deco_ctx;

  deco_ctx = (MotionInfoContexts*) calloc(1, sizeof(MotionInfoContexts) );
  if( deco_ctx == NULL )
    no_mem_exit("create_contexts_MotionInfo: deco_ctx");

  for (j=0; j<2; j++)
  {
    deco_ctx->mb_type_contexts[j] = (BiContextTypePtr) malloc(NUM_MB_TYPE_CTX  * sizeof( BiContextType ) );

    if( deco_ctx->mb_type_contexts[j] == NULL )
      no_mem_exit("create_contexts_MotionInfo: deco_ctx->mb_type_contexts");

    deco_ctx->mv_res_contexts[j] = (BiContextTypePtr) malloc(NUM_MV_RES_CTX  * sizeof( BiContextType ) );

    if( deco_ctx->mv_res_contexts[j] == NULL )
      no_mem_exit("create_contexts_MotionInfo: deco_ctx->mv_res_contexts");
  }

  deco_ctx->ref_no_contexts = (BiContextTypePtr) malloc(NUM_REF_NO_CTX * sizeof( BiContextType ) );

  if( deco_ctx->ref_no_contexts == NULL )
    no_mem_exit("create_contexts_MotionInfo: deco_ctx->ref_no_contexts");

  deco_ctx->delta_qp_contexts = (BiContextTypePtr) malloc(NUM_DELTA_QP_CTX * sizeof( BiContextType ) );

  if( deco_ctx->delta_qp_contexts == NULL )
    no_mem_exit("create_contexts_MotionInfo: deco_ctx->delta_qp_contexts");
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

  int j,k;
  TextureInfoContexts *deco_ctx;

  deco_ctx = (TextureInfoContexts*) calloc(1, sizeof(TextureInfoContexts) );
  if( deco_ctx == NULL )
    no_mem_exit("create_contexts_TextureInfo: deco_ctx");

  for (j=0; j < 6; j++)
  {

    deco_ctx->ipr_contexts[j] = (BiContextTypePtr) malloc(NUM_IPR_CTX  * sizeof( BiContextType ) );

    if( deco_ctx->ipr_contexts[j] == NULL )
      no_mem_exit("create_contexts_TextureInfo: deco_ctx->ipr_contexts");

  }

  for (k=0; k<2; k++)
    for (j=0; j<3; j++)
    {

      deco_ctx->cbp_contexts[k][j] = (BiContextTypePtr) malloc(NUM_CBP_CTX  * sizeof( BiContextType ) );

      if( deco_ctx->cbp_contexts[k][j] == NULL )
        no_mem_exit("create_contexts_TextureInfo: deco_ctx->cbp_contexts");
    }


  for (j=0; j < NUM_TRANS_TYPE; j++)
  {

    deco_ctx->level_context[j] = (BiContextTypePtr) malloc(NUM_LEVEL_CTX  * sizeof( BiContextType ) );

    if( deco_ctx->level_context[j] == NULL )
      no_mem_exit("create_contexts_TextureInfo: deco_ctx->level_context");


    deco_ctx->run_context[j] = (BiContextTypePtr) malloc(NUM_RUN_CTX  * sizeof( BiContextType ) );

    if( deco_ctx->run_context[j] == NULL )
      no_mem_exit("create_contexts_TextureInfo: deco_ctx->run_context");

  }
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
void init_contexts_MotionInfo(struct img_par *img, MotionInfoContexts *deco_ctx, int ini_flag)
{

  int i,j;
  int scale_factor;
  int qp_factor;
  int ini[3];

  if(img->qp <= 10)
    qp_factor=0;
  else
    qp_factor=img->qp-10;

  if ( (img->width*img->height) <=  (IMG_WIDTH * IMG_HEIGHT) ) //  format <= QCIF
        scale_factor=1;
    else
        scale_factor=2;


  for (j=0; j<2; j++)
  {
    if (ini_flag)
    {
      for (i=0; i < NUM_MB_TYPE_CTX; i++)
      {
        ini[0] = MB_TYPE_Ini[j][i][0]+(MB_TYPE_Ini[j][i][3]*qp_factor)/10;
        ini[1] = MB_TYPE_Ini[j][i][1]+(MB_TYPE_Ini[j][i][4]*qp_factor)/10;
        ini[2] = MB_TYPE_Ini[j][i][2]*scale_factor;
        biari_init_context(deco_ctx->mb_type_contexts[j] + i,ini[0],ini[1],ini[2]);
      }
    }
    else
    {
      for (i=0; i < NUM_MB_TYPE_CTX; i++)
        biari_init_context(deco_ctx->mb_type_contexts[j] + i,1,1,100);
    }

    if (ini_flag)
    {
      for (i=0; i < NUM_MV_RES_CTX; i++)
        biari_init_context(deco_ctx->mv_res_contexts[j] + i,MV_RES_Ini[j][i][0]*scale_factor,MV_RES_Ini[j][i][1]*scale_factor,MV_RES_Ini[j][i][2]*scale_factor);
    }
    else
    {
      for (i=0; i < NUM_MV_RES_CTX; i++)
        biari_init_context(deco_ctx->mv_res_contexts[j] + i,1,1,1000);
    }
  }

  if (ini_flag)
  {
    for (i=0; i < NUM_REF_NO_CTX; i++)
      biari_init_context(deco_ctx->ref_no_contexts + i,REF_NO_Ini[i][0]*scale_factor,REF_NO_Ini[i][1]*scale_factor,REF_NO_Ini[i][2]*scale_factor);
  }
  else
  {
    for (i=0; i < NUM_REF_NO_CTX; i++)
      biari_init_context(deco_ctx->ref_no_contexts + i,1,1,1000);
  }

  if (ini_flag)
  {
    for (i=0; i < NUM_DELTA_QP_CTX; i++)
      biari_init_context(deco_ctx->delta_qp_contexts + i,DELTA_QP_Ini[i][0]*scale_factor,DELTA_QP_Ini[i][1]*scale_factor,DELTA_QP_Ini[i][2]*scale_factor);
  }
  else
  {
    for (i=0; i < NUM_DELTA_QP_CTX; i++)
      biari_init_context(deco_ctx->delta_qp_contexts + i,1,1,1000);
  }
}

/*!
 ************************************************************************
 * \brief
 *    Initializes an array of contexts models with some pre-defined
 *    counts (ini_flag = 1) or with a flat histogram (ini_flag = 0)
 *
 ************************************************************************
 */
void init_contexts_TextureInfo(struct img_par *img, TextureInfoContexts *deco_ctx, int ini_flag)
{

  int i,j,k;
  int scale_factor;
  int qp_factor;
  int ini[3];

  if(img->qp <= 10)
    qp_factor=0;
  else
    qp_factor=img->qp-10;

  if ( (img->width*img->height) <=  (IMG_WIDTH * IMG_HEIGHT) ) //  format <= QCIF
        scale_factor=1;
    else
        scale_factor=2;

  for (j=0; j < 6; j++)
  {
    if (ini_flag)
    {
      for (i=0; i < NUM_IPR_CTX; i++)
        biari_init_context(deco_ctx->ipr_contexts[j] + i,IPR_Ini[j][i][0]*scale_factor,IPR_Ini[j][i][1]*scale_factor,IPR_Ini[j][i][2]*scale_factor);
    }
    else
    {
      for (i=0; i < NUM_IPR_CTX; i++)
        biari_init_context(deco_ctx->ipr_contexts[j] + i,2,1,50);
    }
  }

  for (k=0; k<2; k++)
    for (j=0; j<3; j++)
    {
      if (ini_flag)
      {
        for (i=0; i < NUM_CBP_CTX; i++)
        {
          ini[0] = CBP_Ini[k][j][i][0]+(CBP_Ini[k][j][i][3]*qp_factor)/10;
          ini[1] = CBP_Ini[k][j][i][1]+(CBP_Ini[k][j][i][4]*qp_factor)/10;
          ini[2] = CBP_Ini[k][j][i][2]*scale_factor;
          biari_init_context(deco_ctx->cbp_contexts[k][j] + i,ini[0],ini[1],ini[2]);
        }
      }
      else
      {
        for (i=0; i < NUM_CBP_CTX; i++)
          biari_init_context(deco_ctx->cbp_contexts[k][j] + i,1,1,100);
      }
    }


  for (j=0; j < NUM_TRANS_TYPE; j++)
  {

    if (ini_flag)
    {
      for (i=0; i < NUM_LEVEL_CTX; i++)
      {
        ini[0] = (Level_Ini[j][i][0]+(Level_Ini[j][i][3]*qp_factor)/10)*scale_factor;
        ini[1] = (Level_Ini[j][i][1]+(Level_Ini[j][i][4]*qp_factor)/10)*scale_factor;
        ini[2] = Level_Ini[j][i][2]*scale_factor;
        biari_init_context(deco_ctx->level_context[j] + i,ini[0],ini[1],ini[2]);
      }
    }
    else
    {
      for (i=0; i < NUM_LEVEL_CTX; i++)
        biari_init_context(deco_ctx->level_context[j] + i,1,1,100);
    }

    if (ini_flag)
    {
      for (i=0; i < NUM_RUN_CTX; i++)
      {
        ini[0] = Run_Ini[j][i][0]*scale_factor;
        ini[1] = Run_Ini[j][i][1]*scale_factor;
        ini[2] = Run_Ini[j][i][2]*scale_factor;
        biari_init_context(deco_ctx->run_context[j] + i,ini[0],ini[1],ini[2]);
      }
    }
    else
    {
      for (i=0; i < NUM_RUN_CTX; i++)
        biari_init_context(deco_ctx->run_context[j] + i,1,1,100);
    }
  }
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

  int j;

  if( deco_ctx == NULL )
    return;

  for (j=0; j<2; j++)
  {
    if (deco_ctx->mb_type_contexts[j] != NULL)
      free(deco_ctx->mb_type_contexts[j] );

    if (deco_ctx->mv_res_contexts[j]  != NULL)
      free(deco_ctx->mv_res_contexts[j] );
  }

  if (deco_ctx->ref_no_contexts != NULL)
    free(deco_ctx->ref_no_contexts);

  if (deco_ctx->delta_qp_contexts != NULL)
    free(deco_ctx->delta_qp_contexts);

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

  int j,k;

  if( deco_ctx == NULL )
    return;

  for (j=0; j < 6; j++)
  {
    if (deco_ctx->ipr_contexts[j] != NULL)
      free(deco_ctx->ipr_contexts[j]);
  }

  for (k=0; k<2; k++)
    for (j=0; j<3; j++)
    {
      if (deco_ctx->cbp_contexts[k][j] != NULL)
        free(deco_ctx->cbp_contexts[k][j]);
    }

  for (j=0; j < NUM_TRANS_TYPE; j++)
  {
    if (deco_ctx->level_context[j] != NULL)
      free(deco_ctx->level_context[j]);

    if (deco_ctx->run_context[j] != NULL)
      free(deco_ctx->run_context[j]);
  }
  free( deco_ctx );

  return;

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
  int step_h, step_v;
  int i = img->subblock_x;
  int j = img->subblock_y;
  int a, b;
  int act_ctx, act_ctx1;
  int act_sym;
  int mv_local_err;
  int mv_sign;
  int backward = se->value2 & 0x01;
  int k = (se->value2>>1); // MVD component

  MotionInfoContexts *ctx = img->currentSlice->mot_ctx;
  Macroblock *currMB = &img->mb_data[img->current_mb_nr];

  if(backward == 0)
  {
    step_h=img->fw_blc_size_h/BLOCK_SIZE;  // horizontal stepsize
    step_v=img->fw_blc_size_v/BLOCK_SIZE;  // vertical stepsize
  }
  else
  {
    step_h=img->bw_blc_size_h/BLOCK_SIZE;  // horizontal stepsize
    step_v=img->bw_blc_size_v/BLOCK_SIZE;  // vertical stepsize
  }

  if (j==0)
  {
    if (currMB->mb_available[0][1] == NULL)
      b = 0;
    else
      b = absm((currMB->mb_available[0][1])->mvd[backward][BLOCK_SIZE-1][i][k]);
  }
  else
    b = absm(currMB->mvd[backward][j-step_v][i][k]);

  if (i==0)
  {
    if (currMB->mb_available[1][0] == NULL)
      a = 0;
    else
      a = absm((currMB->mb_available[1][0])->mvd[backward][j][BLOCK_SIZE-1][k]);
  }
  else
    a = absm(currMB->mvd[backward][j][i-step_h][k]);

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

  act_ctx1 = act_ctx;

  act_sym = biari_decode_symbol(dep_dp,&ctx->mv_res_contexts[0][act_ctx] );

  if (act_sym != 0)
  {
    act_ctx = 5*k+4;
    mv_sign = biari_decode_symbol(dep_dp, &ctx->mv_res_contexts[1][act_ctx] );
    act_ctx=5*k;
    act_sym = unary_mv_decode(dep_dp,ctx->mv_res_contexts[1]+act_ctx,3);
    act_sym++;

    if(mv_sign)
      act_sym = -act_sym;
  }
  se->value1 = act_sym;

#if TRACE
  fprintf(p_trace, "@%d      %s\t\t\t%3d \n",symbolCount++, se->tracestring, se->value1);
  fflush(p_trace);
#endif
}


/*!
 ************************************************************************
 * \brief
 *    This function is used to arithmetically decode the forward
 *    or backward bidirectional blocksize (for B frames only)
 ************************************************************************
 */
void readBiDirBlkSize2Buffer_CABAC( SyntaxElement *se,
                                    struct inp_par *inp,
                                    struct img_par *img,
                                    DecodingEnvironmentPtr dep_dp)
{
  int act_ctx;
  MotionInfoContexts *ctx = img->currentSlice->mot_ctx;

  act_ctx=4;

  // using the context models of mb_type
  se->value1 = unary_bin_max_decode(dep_dp,ctx->mb_type_contexts[1]+act_ctx,1,6);

#if TRACE
  fprintf(p_trace, "@%d%s\t\t\t%3d\n",symbolCount++, se->tracestring, se->value1);
  fflush(p_trace);
#endif
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
  int log_sym;
  int bin_sym;
  int mode_sym;
  int ct = 0;

  MotionInfoContexts *ctx = (img->currentSlice)->mot_ctx;
  Macroblock *currMB = &img->mb_data[img->current_mb_nr];

  int curr_mb_type;

  if(img->type == INTRA_IMG)  // INTRA-frame
  {
    if (currMB->mb_available[0][1] == NULL)
      b = 0;
    else
      b = (( (currMB->mb_available[0][1])->mb_type != 0) ? 1 : 0 );
    if (currMB->mb_available[1][0] == NULL)
      a = 0;
    else
      a = (( (currMB->mb_available[1][0])->mb_type != 0) ? 1 : 0 );

    act_ctx = a + b;
    act_sym = biari_decode_symbol(dep_dp, ctx->mb_type_contexts[0] + act_ctx);
    se->context = act_ctx; // store context

    if (act_sym==0) // 4x4 Intra
      curr_mb_type = act_sym;
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
  else
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

      // decode unary part
    log_sym = biari_decode_symbol(dep_dp, &ctx->mb_type_contexts[1][act_ctx] );

    if (log_sym != 0)
    {
      act_ctx=4;
      log_sym=0;
      do
      {
        act_sym = biari_decode_symbol(dep_dp, &ctx->mb_type_contexts[1][act_ctx] );
        if (log_sym==0) act_ctx=5;
        log_sym++;
      }
      while ( (act_sym!=0) && (log_sym<2+((img->type == B_IMG_1 || img->type == B_IMG_MULT)?1:0)));
      if( (act_sym!=0) && (log_sym==2+((img->type == B_IMG_1 || img->type == B_IMG_MULT)?1:0))  )
        log_sym++;
    }
    act_sym = (1<<log_sym);


    // decode binary part
    if (log_sym!=0)
    {
      act_ctx=6;
      if (log_sym==(3+((img->type == B_IMG_1 || img->type == B_IMG_MULT)?1:0)) ) log_sym=2; // only 2 LSBs are actually set for mode 7-9 (P-frame) or 15-17 (B-frame)
        bin_sym=0;
      do
      {
        log_sym--;
        bin_sym <<=1;
        bin_sym |=  (biari_decode_symbol(dep_dp,  &ctx->mb_type_contexts[1][act_ctx] ));
      }
      while (log_sym!=0);
      act_sym += bin_sym;
    }
    act_sym--;

    if (act_sym<=8 || (((img->type == B_IMG_1 || img->type == B_IMG_MULT)?1:0) && act_sym<=16))
      curr_mb_type = act_sym;
    else  // additional info for 16x16 Intra-mode
    {
      act_ctx = 7;
      mode_sym =  biari_decode_symbol(dep_dp, ctx->mb_type_contexts[1] + act_ctx ); // decoding of AC/no AC
      act_sym += mode_sym*12;

      // decoding of cbp: 0,1,2
      act_ctx = 8;
      mode_sym = biari_decode_symbol(dep_dp, ctx->mb_type_contexts[1] + act_ctx );
      if (mode_sym != 0)
      {
        act_sym+=4;
        mode_sym = biari_decode_symbol(dep_dp, ctx->mb_type_contexts[1] + act_ctx );
        if (mode_sym != 0)
          act_sym+=4;
      }

      // decoding of I pred-mode: 0,1,2,3
      act_ctx = 9;
      mode_sym = biari_decode_symbol(dep_dp, ctx->mb_type_contexts[1] + act_ctx );
      act_sym += mode_sym*2;
      mode_sym = biari_decode_symbol(dep_dp, ctx->mb_type_contexts[1] + act_ctx );
      act_sym += mode_sym;
      curr_mb_type = act_sym;
    }
  }
    se->value1 = curr_mb_type;

#if TRACE
    fprintf(p_trace, "@%d%s\t\t\t%3d\n",symbolCount++, se->tracestring, se->value1);
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
  static int prev_sym = 0;
  static int count = 0; // to detect a new row of intra prediction modes

  TextureInfoContexts *ctx = img->currentSlice->tex_ctx;

  if (count % 2 == 0)
    prev_sym = 0;

  se->value1 = unary_bin_max_decode(dep_dp,ctx->ipr_contexts[prev_sym],1,5);
  prev_sym = se->value1;
  se->value2 = unary_bin_max_decode(dep_dp,ctx->ipr_contexts[prev_sym],1,5);
  prev_sym = se->value2;


  if(++count == MB_BLOCK_SIZE/2) // all modes of one MB have been processed
    count=0;

#if TRACE
  fprintf(p_trace, "@%d%s\t\t\t%3d\n",symbolCount++, se->tracestring, se->value1);
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
  Macroblock *currMB = &img->mb_data[img->current_mb_nr];

  int a, b;
  int act_ctx;
  int act_sym;

  if (currMB->mb_available[0][1] == NULL)
    b = 0;
  else
      b = ( ((currMB->mb_available[0][1])->predframe_no != 0) ? 1 : 0);
  if (currMB->mb_available[1][0] == NULL)
    a = 0;
  else
      a = ( ((currMB->mb_available[1][0])->predframe_no != 0) ? 1 : 0);
  act_ctx = a + 2*b;
  se->context = act_ctx; // store context

  act_sym = biari_decode_symbol(dep_dp,ctx->ref_no_contexts + act_ctx );

  if (act_sym != 0)
  {
    act_ctx = 4;
    act_sym = unary_bin_decode(dep_dp,ctx->ref_no_contexts+act_ctx,1);
    act_sym++;
  }
  se->value1 = act_sym;

#if TRACE
  fprintf(p_trace, "@%d%s\t\t\t%3d",symbolCount++, se->tracestring, se->value1);
  fprintf(p_trace," c: %d :%d \n",ctx->ref_no_contexts[act_ctx].cum_freq[0],ctx->ref_no_contexts[act_ctx].cum_freq[1]);
  fflush(p_trace);
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
  int step_h, step_v;
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
  Macroblock *currMB = &img->mb_data[img->current_mb_nr];

  step_h=BLOCK_STEP[img->mb_mode][0];
  step_v=BLOCK_STEP[img->mb_mode][1];

  if (j==0)
  {
    if (currMB->mb_available[0][1] == NULL)
      b = 0;
    else
      b = absm((currMB->mb_available[0][1])->mvd[0][BLOCK_SIZE-1][i][k]);
  }
  else
    b = absm(currMB->mvd[0][j-step_v][i][k]);

  if (i==0)
  {
    if (currMB->mb_available[1][0] == NULL)
      a = 0;
    else
      a = absm((currMB->mb_available[1][0])->mvd[0][j][BLOCK_SIZE-1][k]);
  }
  else
    a = absm(currMB->mvd[0][j][i-step_h][k]);

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
        mv_pred_res = 0;
    else
    {
  
        act_ctx=5*k;
        act_sym = unary_mv_decode(dep_dp,ctx->mv_res_contexts[1]+act_ctx,3);
        act_sym++;
        act_ctx=5*k+4;
        mv_sign = biari_decode_symbol(dep_dp,&ctx->mv_res_contexts[1][act_ctx] );

        mv_pred_res = ((mv_sign != 0) ? (-act_sym) : act_sym);
    }
    se->value1 = mv_pred_res;

#if TRACE
    fprintf(p_trace, "@%d%s\t\t\t%3d\n",symbolCount++, se->tracestring, se->value1);
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
void readDquantFromBuffer_CABAC(SyntaxElement *se,
                                struct inp_par *inp,
                                struct img_par *img,
                                DecodingEnvironmentPtr dep_dp)
{
  MotionInfoContexts *ctx = img->currentSlice->mot_ctx;
  Macroblock *currMB = &img->mb_data[img->current_mb_nr];

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
  Macroblock *currMB = &img->mb_data[img->current_mb_nr];

  int mb_x, mb_y;
  int a, b;
  int curr_cbp_ctx, curr_cbp_idx;
  int cbp = 0;
  int cbp_bit;
    int mask;

  if ( se->type == SE_CBP_INTRA )
    curr_cbp_idx = 0;
  else
    curr_cbp_idx = 1;

  //  coding of luma part (bit by bit)
  for (mb_y=0; mb_y < 4; mb_y += 2)
  {
    for (mb_x=0; mb_x < 4; mb_x += 2)
    {

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
            cbp_bit = biari_decode_symbol(dep_dp, ctx->cbp_contexts[curr_cbp_idx][0] + curr_cbp_ctx );
            if (cbp_bit) cbp += mask;
    }
  }

  // coding of chroma part
  b = 0;
  if (currMB->mb_available[0][1] != NULL)
    b = ((currMB->mb_available[0][1])->cbp > 15) ? 1 : 0;

  a = 0;
  if (currMB->mb_available[1][0] != NULL)
    a = ((currMB->mb_available[1][0])->cbp > 15) ? 1 : 0;

  curr_cbp_ctx = a+2*b;
  cbp_bit = biari_decode_symbol(dep_dp, ctx->cbp_contexts[curr_cbp_idx][1] + curr_cbp_ctx );

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
    cbp_bit = biari_decode_symbol(dep_dp, ctx->cbp_contexts[curr_cbp_idx][2] + curr_cbp_ctx );
    cbp += (cbp_bit == 1) ? 32 : 16;
  }

  se->value1 = cbp;

#if TRACE
  fprintf(p_trace, "@%d      %s\t\t\t%3d\n",symbolCount++, se->tracestring, se->value1);
  fflush(p_trace);
#endif
}

/*!
 ************************************************************************
 * \brief
 *    This function is used to arithmetically decode level and
 *    run of a given MB.
 ************************************************************************
 */
void readRunLevelFromBuffer_CABAC(SyntaxElement *se,
                                  struct inp_par *inp,
                                  struct img_par *img,
                                  DecodingEnvironmentPtr dep_dp)
{
  int level;
  int run=0;
  const int curr_ctx_idx = se->context;
  int curr_level_ctx;
  int sign_of_level;
  int max_run;

  TextureInfoContexts *ctx = img->currentSlice->tex_ctx;
  // Macroblock *currMB = &img->mb_data[img->current_mb_nr];

  level = unary_level_decode(dep_dp,ctx->level_context[curr_ctx_idx]);

  if (level!=0)
  {
    curr_level_ctx = 3;
    sign_of_level = biari_decode_symbol(dep_dp, ctx->level_context[curr_ctx_idx] + curr_level_ctx );
    if (sign_of_level) level = (-1)*level;
      if (curr_ctx_idx != 0 && curr_ctx_idx != 6 && curr_ctx_idx != 5) // not double scan and not DC-chroma
        run = unary_bin_decode(dep_dp,ctx->run_context[curr_ctx_idx],1);
      else
      {
        max_run =  (curr_ctx_idx == 0) ? 7 : 3;  // if double scan max_run = 7; if DC-chroma max_run = 3;
        run = unary_bin_max_decode(dep_dp,ctx->run_context[curr_ctx_idx],1,max_run);
      }
  }
  se->value1 = level;
  se->value2 = run;


#if TRACE
  fprintf(p_trace, "@%d%s\t\t\t%3d \n",symbolCount++, se->tracestring, se->value1);
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
  int info;
  int BitsUsedByHeader = 0, ByteStartPosition;
  int newframe = 0;   //WYK: Oct. 8, 2001, change the method to find a new frame

  currStream->frame_bitoffset =0;

  memset (code_buffer, 0xff, MAX_CODED_FRAME_SIZE);   // this prevents a buffer full with zeros
  BitstreamLengthInBytes = currStream->bitstream_length = GetOneSliceIntoSourceBitBuffer(img, inp, code_buffer);

  // Here we are ready to interpret the picture and slice headers.  Since
  // SliceHeader() gets the data out of the UVLC's len/info
  // array, we need to convert the start of our slice to such a format.


  if (BitstreamLengthInBytes < 4)
    return EOS;

  // Now we have the bits between the current startcode (inclusive) and the
  // next start code in code_buffer.  Now decode the start codes and the headers
  if (31 != GetVLCSymbol (code_buffer, 0, &info, BitstreamLengthInBytes))
  {
    snprintf (errortext, ET_SIZE, "readSliceCABAC: Panic, expected start code symbol, found wrong len");
    error(errortext, 600);
  }
  currStream->frame_bitoffset +=31;
  BitsUsedByHeader+=SliceHeader(img, inp);

  //WYK: Oct. 8, 2001, change the method to find a new frame
  if(img->tr != img->tr_old)
    newframe = 1;
  else 
    newframe = 0;
  img->tr_old = img->tr;
    
  // if the TR of current slice is not identical to the TR of previous received slice, we have a new frame
  if(newframe)
    current_header = SOP;
  else
    current_header = SOS;

  ByteStartPosition = currStream->frame_bitoffset/8;
  if ((currStream->frame_bitoffset)%8 != 0)
    ByteStartPosition++;
  arideco_start_decoding(dep, code_buffer, ByteStartPosition, read_len);


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
 *    decoding of unary binarization of the absolute value
 *    of a level using 3 distinct models by separating the first,
 *    the second and all remaining bins
 ***********************************************************************
 */
unsigned int unary_level_decode(DecodingEnvironmentPtr dep_dp,
                                BiContextTypePtr ctx)
{
  unsigned int l;
  unsigned int symbol;
  int bin=1;
  BiContextTypePtr ictx=ctx;

  symbol = biari_decode_symbol(dep_dp, ictx );

  if (symbol==0)
    return 0;
  else
  {
    symbol=0;
    ictx++;
    do
    {
      l=biari_decode_symbol(dep_dp, ictx  );
      if ((++bin)==2) ictx++;
      symbol++;
    }
    while (l!=0);
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
 *    decoding of unary binarization of the absolute value of a
 *    mv component using 4 distinct models by separating the first,
 *    the second, intermediate and all remaining bins
 ***********************************************************************
 */
unsigned int unary_mv_decode(DecodingEnvironmentPtr dep_dp,
                             BiContextTypePtr ctx,
                             unsigned int max_bin)
{
  unsigned int l;
  unsigned int bin=1;
  unsigned int symbol;

  BiContextTypePtr ictx=ctx;

  symbol = biari_decode_symbol(dep_dp, ictx );

  if (symbol==0)
    return 0;
  else
  {
    symbol=0;
    ictx++;
    do
    {
      l=biari_decode_symbol(dep_dp, ictx  );
      if ((++bin)==2) ictx++;
      if (bin==max_bin) ictx++;
      symbol++;
    }
    while (l!=0);
    return symbol;
  }
}

/*!
 ************************************************************************
 * \brief
 *    finding end of a slice in case this is not the end of a frame
 ************************************************************************
 */
int cabac_startcode_follows(struct img_par *img, struct inp_par *inp)
{
  Slice *currSlice = img->currentSlice;
  if (img->current_mb_nr == currSlice->last_mb_nr)
    return TRUE;
  return FALSE;
}

