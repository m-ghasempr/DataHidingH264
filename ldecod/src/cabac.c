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

  for (j=0; j<3; j++)
  {
    deco_ctx->mb_type_contexts[j] = (BiContextTypePtr) malloc(NUM_MB_TYPE_CTX  * sizeof( BiContextType ) );
    if( deco_ctx->mb_type_contexts[j] == NULL )
      no_mem_exit("create_contexts_MotionInfo: deco_ctx->mb_type_contexts");
  }
  for (j=0; j<2; j++)
  {
    deco_ctx->b8_type_contexts[j] = (BiContextTypePtr) malloc(NUM_B8_TYPE_CTX * sizeof( BiContextType ) );
    if( deco_ctx->b8_type_contexts[j] == NULL ) 
      no_mem_exit("create_contexts_MotionInfo: deco_ctx->b8_type_contexts");

    deco_ctx->mv_res_contexts[j] = (BiContextTypePtr) malloc(NUM_MV_RES_CTX  * sizeof( BiContextType ) );
    if( deco_ctx->mv_res_contexts[j] == NULL )
      no_mem_exit("create_contexts_MotionInfo: deco_ctx->mv_res_contexts");

    deco_ctx->ref_no_contexts[j] = (BiContextTypePtr) malloc(NUM_REF_NO_CTX * sizeof( BiContextType ) );
    if( deco_ctx->ref_no_contexts[j] == NULL )
      no_mem_exit("create_contexts_MotionInfo: deco_ctx->ref_no_contexts");
  }

  deco_ctx->delta_qp_inter_contexts = (BiContextTypePtr) malloc(NUM_DELTA_QP_CTX * sizeof( BiContextType ) );
  if( deco_ctx->delta_qp_inter_contexts == NULL )
    no_mem_exit("create_contexts_MotionInfo: deco_ctx->delta_qp_inter_contexts");
  
  deco_ctx->delta_qp_intra_contexts = (BiContextTypePtr) malloc(NUM_DELTA_QP_CTX * sizeof( BiContextType ) );
  if( deco_ctx->delta_qp_intra_contexts == NULL )
    no_mem_exit("create_contexts_MotionInfo: deco_ctx->delta_qp_intra_contexts");
  
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
  const int max_ipr=9;

  deco_ctx = (TextureInfoContexts*) calloc(1, sizeof(TextureInfoContexts) );
  if( deco_ctx == NULL )
    no_mem_exit("create_contexts_TextureInfo: deco_ctx");

  for (j=0; j < max_ipr; j++)
  {
    deco_ctx->ipr_contexts[j] = (BiContextTypePtr) malloc(NUM_IPR_CTX  * sizeof( BiContextType ) );
    if( deco_ctx->ipr_contexts[j] == NULL )
      no_mem_exit("create_contexts_TextureInfo: deco_ctx->ipr_contexts");

  }

  for (k=0; k<2; k++)
  {
    for (j=0; j<3; j++)
    {
      deco_ctx->cbp_contexts[k][j] = (BiContextTypePtr) malloc(NUM_CBP_CTX  * sizeof( BiContextType ) );
      if( deco_ctx->cbp_contexts[k][j] == NULL )
        no_mem_exit("create_contexts_TextureInfo: deco_ctx->cbp_contexts");
    }
  }

  for (j=0; j < 4*NUM_TRANS_TYPE; j++)
  {
    deco_ctx->level_context[j] = (BiContextTypePtr) malloc(NUM_LEVEL_CTX  * sizeof( BiContextType ) );
    if( deco_ctx->level_context[j] == NULL )
      no_mem_exit("create_contexts_TextureInfo: deco_ctx->level_context");
  }

  for (j=0; j < 2*NUM_TRANS_TYPE; j++)
  {
    deco_ctx->run_context[j] = (BiContextTypePtr) malloc(NUM_RUN_CTX  * sizeof( BiContextType ) );
    if( deco_ctx->run_context[j] == NULL )
      no_mem_exit("create_contexts_TextureInfo: deco_ctx->run_context");
  }

  for (j=0; j < NUM_TRANS_TYPE; j++)
  {
    deco_ctx->coeff_count_context[j] = (BiContextTypePtr) malloc(NUM_COEFF_COUNT_CTX  * sizeof( BiContextType ) );
    if( deco_ctx->coeff_count_context[j] == NULL )
      no_mem_exit("create_contexts_TextureInfo: deco_ctx->coeff_count_context");
  }

  for (j=0; j < 2*NUM_TRANS_TYPE_ABT; j++)
  {
    deco_ctx->ABT_run_context[j] = (BiContextTypePtr) malloc(NUM_RUN_CTX_ABT  * sizeof( BiContextType ) );
    if( deco_ctx->ABT_run_context[j] == NULL )
      no_mem_exit("create_contexts_TextureInfo: deco_ctx->ABT_run_context");
  }

  for (j=0; j < NUM_TRANS_TYPE_ABT; j++)
  {
    deco_ctx->ABT_coeff_count_context[j] = (BiContextTypePtr) malloc(NUM_COEFF_COUNT_CTX_ABT  * sizeof( BiContextType ) );
    if( deco_ctx->ABT_coeff_count_context[j] == NULL )
      no_mem_exit("create_contexts_TextureInfo: deco_ctx->ABT_coeff_count_context");
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

  qp_factor=min(max(0,img->qp-10-SHIFT_QP),21);

  if ( (img->width*img->height) <=  (IMG_WIDTH * IMG_HEIGHT) ) //  format <= QCIF
        scale_factor=1;
    else
        scale_factor=2;


  for (j=0; j<3; j++)
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
  }

  for (j=0; j<2; j++)
  {
    if (ini_flag)
    {
      for (i=0; i < NUM_B8_TYPE_CTX; i++)
        biari_init_context(deco_ctx->b8_type_contexts[j] + i,B8_TYPE_Ini[j][i][0]*scale_factor,B8_TYPE_Ini[j][i][1]*scale_factor,B8_TYPE_Ini[j][i][2]*scale_factor);
    }
    else
    {
      for (i=0; i < NUM_B8_TYPE_CTX; i++)
        biari_init_context(deco_ctx->b8_type_contexts[j] + i,1,1,1000);
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

    if (ini_flag)
    {
      for (i=0; i < NUM_REF_NO_CTX; i++)
        biari_init_context(deco_ctx->ref_no_contexts[j] + i,REF_NO_Ini[i][0]*scale_factor,REF_NO_Ini[i][1]*scale_factor,REF_NO_Ini[i][2]*scale_factor);
    }
    else
    {
      for (i=0; i < NUM_REF_NO_CTX; i++)
        biari_init_context(deco_ctx->ref_no_contexts[j] + i,1,1,1000);
    }
  }

  if (ini_flag)
  {
    for (i=0; i < NUM_DELTA_QP_CTX; i++)
      biari_init_context(deco_ctx->delta_qp_inter_contexts + i,DELTA_QP_Ini[i][0]*scale_factor,DELTA_QP_Ini[i][1]*scale_factor,DELTA_QP_Ini[i][2]*scale_factor);
  }
  else
  {
    for (i=0; i < NUM_DELTA_QP_CTX; i++)
      biari_init_context(deco_ctx->delta_qp_inter_contexts + i,1,1,1000);
  }

  if (ini_flag)
  {
    for (i=0; i < NUM_DELTA_QP_CTX; i++)
      biari_init_context(deco_ctx->delta_qp_intra_contexts + i,DELTA_QP_Ini[i][0]*scale_factor,DELTA_QP_Ini[i][1]*scale_factor,DELTA_QP_Ini[i][2]*scale_factor);
  }
  else
  {
    for (i=0; i < NUM_DELTA_QP_CTX; i++)
      biari_init_context(deco_ctx->delta_qp_intra_contexts + i,1,1,1000);
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
  const int max_ipr=9;


  qp_factor=min(max(0,img->qp-10-SHIFT_QP),21);

  if ( (img->width*img->height) <=  (IMG_WIDTH * IMG_HEIGHT) ) //  format <= QCIF
        scale_factor=1;
    else
        scale_factor=2;

  for (j=0; j < max_ipr; j++)
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

  // other init mapping 
  qp_factor=min(max(0,28-(img->qp-SHIFT_QP)),24);

  for (j=0; j < 4*NUM_TRANS_TYPE; j++)
  {
    if (ini_flag)
    {
      for (i=0; i < NUM_LEVEL_CTX; i++)
      {
        ini[0] = (Level_Ini[j][i][0]+(Level_Ini[j][i][3]*qp_factor)/24)*scale_factor;
        ini[1] = (Level_Ini[j][i][1]+(Level_Ini[j][i][4]*qp_factor)/24)*scale_factor;
        ini[2] = Level_Ini[j][i][2]*scale_factor;
        biari_init_context(deco_ctx->level_context[j] + i,ini[0],ini[1],ini[2]);
      }
    }
    else
    {
      for (i=0; i < NUM_LEVEL_CTX; i++)
        biari_init_context(deco_ctx->level_context[j] + i,1,1,100);
    }
  }

  for (j=0; j < 2*NUM_TRANS_TYPE; j++)
  {
    if (ini_flag)
    {
      for (i=0; i < NUM_RUN_CTX; i++)
      {
        ini[0] = (Run_Ini[j][i][0]+(Run_Ini[j][i][3]*qp_factor)/24)*scale_factor;
        ini[1] = (Run_Ini[j][i][1]+(Run_Ini[j][i][4]*qp_factor)/24)*scale_factor;
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

  for (j=0; j < NUM_TRANS_TYPE; j++)
  {
    if (ini_flag)
    {
      for (i=0; i < NUM_COEFF_COUNT_CTX; i++)
      {
        ini[0] = (Coeff_Count_Ini[j][i][0]+(Coeff_Count_Ini[j][i][3]*qp_factor)/24)*scale_factor;
        ini[1] = (Coeff_Count_Ini[j][i][1]+(Coeff_Count_Ini[j][i][4]*qp_factor)/24)*scale_factor;
        ini[2] = Coeff_Count_Ini[j][i][2]*scale_factor;
        biari_init_context(deco_ctx->coeff_count_context[j] + i,ini[0],ini[1],ini[2]);
      }
    }
    else
    {
      for (i=0; i < NUM_COEFF_COUNT_CTX; i++)
        biari_init_context(deco_ctx->coeff_count_context[j] + i,1,1,100);
    }
  }

  for (j=0; j < 2*NUM_TRANS_TYPE_ABT; j++)
  {
    if (ini_flag)
    {
      for (i=0; i < NUM_RUN_CTX_ABT; i++)
      {
        ini[0] = (ABT_Run_Ini[j][i][0]+(ABT_Run_Ini[j][i][3]*qp_factor)/24)*scale_factor;
        ini[1] = (ABT_Run_Ini[j][i][1]+(ABT_Run_Ini[j][i][4]*qp_factor)/24)*scale_factor;
        ini[2] = ABT_Run_Ini[j][i][2]*scale_factor;
        biari_init_context(deco_ctx->ABT_run_context[j] + i,ini[0],ini[1],ini[2]);
      }
    }
    else
    {
      for (i=0; i < NUM_RUN_CTX_ABT; i++)
        biari_init_context(deco_ctx->ABT_run_context[j] + i,1,1,INICNT_ABT);
    }
  }

  for (j=0; j < NUM_TRANS_TYPE_ABT; j++)
  {
    if (ini_flag)
    {
      for (i=0; i < NUM_COEFF_COUNT_CTX_ABT; i++)
      {
        ini[0] = (ABT_Coeff_Count_Ini[j][i][0]+(ABT_Coeff_Count_Ini[j][i][3]*qp_factor)/24)*scale_factor;
        ini[1] = (ABT_Coeff_Count_Ini[j][i][1]+(ABT_Coeff_Count_Ini[j][i][4]*qp_factor)/24)*scale_factor;
        ini[2] = ABT_Coeff_Count_Ini[j][i][2]*scale_factor;
        biari_init_context(deco_ctx->ABT_coeff_count_context[j] + i,ini[0],ini[1],ini[2]);
      }
    }
    else
    {
      for (i=0; i < NUM_COEFF_COUNT_CTX_ABT; i++)
        biari_init_context(deco_ctx->ABT_coeff_count_context[j] + i,1,1,INICNT_ABT);
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

  for (j=0; j<3; j++)
  {
    if (deco_ctx->mb_type_contexts[j] != NULL)
      free(deco_ctx->mb_type_contexts[j]);
  }
  for (j=0; j<2; j++)
  {
    if (deco_ctx->b8_type_contexts[j] != NULL)
      free(deco_ctx->b8_type_contexts[j]);

    if (deco_ctx->mv_res_contexts[j]  != NULL)
      free(deco_ctx->mv_res_contexts [j]);

    if (deco_ctx->ref_no_contexts[j]  != NULL)
      free(deco_ctx->ref_no_contexts [j]);
  }

  if (deco_ctx->delta_qp_inter_contexts != NULL)
    free(deco_ctx->delta_qp_inter_contexts);

  if (deco_ctx->delta_qp_intra_contexts != NULL)
    free(deco_ctx->delta_qp_intra_contexts);

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

  static const int max_ipr=9;


  if( deco_ctx == NULL )
    return;

  for (j=0; j < max_ipr; j++)
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

  for (j=0; j < 4*NUM_TRANS_TYPE; j++)
  {
    if (deco_ctx->level_context[j] != NULL)
      free(deco_ctx->level_context[j]);
  }

  for (j=0; j < 2*NUM_TRANS_TYPE; j++)
  {
    if (deco_ctx->run_context[j] != NULL)
      free(deco_ctx->run_context[j]);
  }

  for (j=0; j < NUM_TRANS_TYPE; j++)
  {
    if (deco_ctx->coeff_count_context[j] != NULL)
      free(deco_ctx->coeff_count_context[j]);
  }

  for (j=0; j < 2*NUM_TRANS_TYPE_ABT; j++)
  {
    if (deco_ctx->ABT_run_context[j] != NULL)
      free(deco_ctx->ABT_run_context[j]);
  }

  for (j=0; j < NUM_TRANS_TYPE_ABT; j++)
  {
    if (deco_ctx->ABT_coeff_count_context[j] != NULL)
      free(deco_ctx->ABT_coeff_count_context[j]);
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
  Macroblock *currMB = &img->mb_data[img->current_mb_nr];


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
    mv_sign = biari_decode_symbol(dep_dp, &ctx->mv_res_contexts[1][act_ctx] );
    act_ctx=5*k;
    act_sym = unary_mv_decode(dep_dp,ctx->mv_res_contexts[1]+act_ctx,3);
    act_sym++;

    if(mv_sign)
      act_sym = -act_sym;
  }
  se->value1 = act_sym;

#if TRACE
  fprintf(p_trace, "@%d      %s\t\t\t%d ",symbolCount++, se->tracestring, se->value1);
  fflush(p_trace);
#endif
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

  MotionInfoContexts *ctx = (img->currentSlice)->mot_ctx;
  Macroblock *currMB = &img->mb_data[img->current_mb_nr];


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

    if (bframe)
    {
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
      if (biari_decode_symbol(dep_dp, &ctx->mb_type_contexts[1][act_ctx] ))
      {
        if (biari_decode_symbol(dep_dp, &ctx->mb_type_contexts[1][4] ))
        {
          if (biari_decode_symbol(dep_dp, &ctx->mb_type_contexts[1][7] ))   act_sym = 7;
          else                                                              act_sym = 6;
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
      else
      {
        act_sym = 0;
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

#if TRACE
  fprintf(p_trace, "@%d%s\t\t\t%d",symbolCount++, se->tracestring, se->value1);
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
  static const int    right[8] = {0, 0, 1, 1, 0, 0, 1, 1};
  Macroblock          *currMB  = &(img->mb_data[img->current_mb_nr]);
  TextureInfoContexts *ctx     = img->currentSlice->tex_ctx;
  int                 prev_sym;

  //--- first symbol ---
  if (right[se->context/2])             prev_sym = currMB->intra_pred_modes[se->context-3];
  else if (currMB->mb_available[1][0])  prev_sym = currMB->mb_available[1][0]->intra_pred_modes[se->context+5];
  else                                  prev_sym = 0;

  se->value1  = unary_bin_max_decode(dep_dp,ctx->ipr_contexts[prev_sym],1,8);

  //--- second symbol ---
  prev_sym = se->value1;
  se->value2  = unary_bin_max_decode(dep_dp,ctx->ipr_contexts[prev_sym],1,8);

#if TRACE
  fprintf(p_trace, "@%d%s\t\t\t%d\n",symbolCount++, se->tracestring, se->value1);
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

  int   addctx = se->context;
  int   a, b;
  int   act_ctx;
  int   act_sym;
  int** refframe_array = ((img->type==B_IMG_1 || img->type==B_IMG_MULT) ? img->fw_refFrArr : refFrArr);
  
  if (currMB->mb_available[0][1] == NULL)
    b = 0;
  else
    b = (refframe_array[img->block_y+img->subblock_y-1][img->block_x+img->subblock_x] > 0 ? 1 : 0);
  if (currMB->mb_available[1][0] == NULL)
    a = 0;
  else 
    a = (refframe_array[img->block_y+img->subblock_y][img->block_x+img->subblock_x-1] > 0 ? 1 : 0);

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
  fprintf(p_trace, "@%d%s\t\t\t%d",symbolCount++, se->tracestring, se->value1);
  fprintf(p_trace," c: %d :%d \n",ctx->ref_no_contexts[addctx][act_ctx].cum_freq[0],ctx->ref_no_contexts[addctx][act_ctx].cum_freq[1]);
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
    mv_sign = biari_decode_symbol(dep_dp,&ctx->mv_res_contexts[1][act_ctx] );

    act_ctx=5*k;
    act_sym = unary_mv_decode(dep_dp,ctx->mv_res_contexts[1]+act_ctx,3);
    act_sym++;
    mv_pred_res = ((mv_sign != 0) ? (-act_sym) : act_sym);
  }
  se->value1 = mv_pred_res;

#if TRACE
  fprintf(p_trace, "@%d%s\t\t\t%d\n",symbolCount++, se->tracestring, se->value1);
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
void readDquant_inter_FromBuffer_CABAC(SyntaxElement *se,
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

  act_sym = biari_decode_symbol(dep_dp,ctx->delta_qp_inter_contexts + act_ctx );
  if (act_sym != 0)
  {
    act_ctx = 2;
    act_sym = unary_bin_decode(dep_dp,ctx->delta_qp_inter_contexts+act_ctx,1);
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
 *    This function is used to arithmetically decode the delta qp
 *     of a given MB.
 ************************************************************************
 */
void readDquant_intra_FromBuffer_CABAC(SyntaxElement *se,
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

  act_sym = biari_decode_symbol(dep_dp,ctx->delta_qp_intra_contexts + act_ctx );
  if (act_sym != 0)
  {
    act_ctx = 2;
    act_sym = unary_bin_decode(dep_dp,ctx->delta_qp_intra_contexts+act_ctx,1);
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
      cbp_bit = biari_decode_symbol(dep_dp, ctx->cbp_contexts[curr_cbp_idx][0] + curr_cbp_ctx );
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
  fprintf(p_trace, "@%d      %s\t\t\t%d",symbolCount++, se->tracestring, se->value1);
  fflush(p_trace);
#endif
}

/*!
 ************************************************************************
 * \brief
 *    This function is used to arithmetically decode coeff_count,level and
 *    run of a given MB.
 ************************************************************************
 */
void readRunLevelFromBuffer_CABAC (SyntaxElement *se,
                                   struct inp_par *inp,
                                   struct img_par *img,    
                                   DecodingEnvironmentPtr dep_dp)
{
  int level = 0;
  int run   = 0;
  const int curr_ctx_idx = se->context;
  int curr_level_ctx, curr_run_ctx = 0;
  int sign_of_level;
  static int max_run;
  static int coeff_count=0;
  int act_ctx=0;
  int i,j, a,b;
  static int send_eob = 0;
  static int coeff_count_all = 0;
  int max_coeff[9] = {8,16,16,16,15,4,4,15,15};
  static int prev_sym[3] = {0,0,0};
  static int count[3] = {0,0,0};
  static int prevLevel = 0;
  int changed_ctx_idx;
  
  int b8, curr_abt_ctx_idx=-1; // =-1: avoid warnings
  static unsigned int max_ccnt_abt [6] = {63,32,16,63,32,16};
  static unsigned int max_coeff_abt[6] = {64,32,16,64,32,16};

  TextureInfoContexts *ctx = img->currentSlice->tex_ctx;
  Macroblock *currMB = &img->mb_data[img->current_mb_nr];

  if ((coeff_count == 0) && !send_eob)// get number of coeff.
  {
    switch(curr_ctx_idx)
    {
    case 0://max 8 coeffs, double_scan
    {
      //context determination
      act_ctx = (count[0] == 0) ? 0 : ((prev_sym[0] == 0) ? 1 : 2);
      coeff_count = biari_decode_symbol(dep_dp, ctx->coeff_count_context[curr_ctx_idx] + act_ctx );
      if (coeff_count != 0)
      {
        coeff_count = unary_bin_max_decode(dep_dp,ctx->coeff_count_context[curr_ctx_idx]+4,1,max_coeff[curr_ctx_idx]-1);
        coeff_count++;
      }
      prev_sym[0] = coeff_count;
      if(++count[0] == 2) count[0]=0;
      break;
    }
    case 1:
    case 2:
    {
      //context termination
      j = img->subblock_y;
      i = img->subblock_x;
      if (j==0)
        b = (currMB->mb_available[0][1] == NULL)?0:(((currMB->mb_available[0][1])->coeffs_count[BLOCK_SIZE-1][i]==0)?0:1);
      else
        b = (currMB->coeffs_count[j-1][i]==0)?0:1;
          
      if (i==0)
        a = (currMB->mb_available[1][0] == NULL)?0:(((currMB->mb_available[1][0])->coeffs_count[j][BLOCK_SIZE-1]==0)?0:1);
      else
        a = (currMB->coeffs_count[j][i-1]==0)?0:1;

      act_ctx = a+2*b;
      coeff_count = biari_decode_symbol(dep_dp, ctx->coeff_count_context[curr_ctx_idx] + act_ctx );
      if (coeff_count != 0)
      {
        coeff_count = unary_bin_max_decode(dep_dp,ctx->coeff_count_context[curr_ctx_idx]+4,1,max_coeff[curr_ctx_idx]-1);
        coeff_count++;
      }      
      currMB->coeffs_count[j][i] = coeff_count; 
      break;
    }
    case 4://max 15 coeffs // 16x16 luma ac
    {
      //context determination
      act_ctx = (count[2] == 0) ? 0 : ((prev_sym[2] == 0) ? 1 : 2);
      coeff_count = biari_decode_symbol(dep_dp, ctx->coeff_count_context[curr_ctx_idx] + act_ctx );
      if (coeff_count != 0)
      {
        coeff_count = unary_bin_max_decode(dep_dp,ctx->coeff_count_context[curr_ctx_idx]+4,1,max_coeff[curr_ctx_idx]-1);
        coeff_count++;
      }
      prev_sym[2] = coeff_count;
      if(++count[2] == 16) count[2]=0;
      break;
    }
    case 5:
    case 6://max 4 coeffs, chroma dc different ctx uv for first bin
    {
      act_ctx = (se->k == 0)?0:1;
      //act_ctx = (se->k < 4)?0:1;
      coeff_count = biari_decode_symbol(dep_dp, ctx->coeff_count_context[curr_ctx_idx] + act_ctx );
      if (coeff_count != 0)
      {
        coeff_count = unary_bin_max_decode(dep_dp,ctx->coeff_count_context[curr_ctx_idx]+4,1,max_coeff[curr_ctx_idx]-1);
        coeff_count++;
      }
      break;
    }
    case 7:
    case 8://max 15 coeffs, chroma ac different ctx uv for first bin
    {
      act_ctx = (se->k < 4)?0:1;
      coeff_count = biari_decode_symbol(dep_dp, ctx->coeff_count_context[curr_ctx_idx] + act_ctx );
      if (coeff_count != 0)
      {
        coeff_count = unary_bin_max_decode(dep_dp,ctx->coeff_count_context[curr_ctx_idx]+4,1,max_coeff[curr_ctx_idx]-1);
        coeff_count++;
      }
      break;
    }
    case 3:
    {
      act_ctx = 0;
      coeff_count = biari_decode_symbol(dep_dp, ctx->coeff_count_context[curr_ctx_idx] + act_ctx );
      if (coeff_count != 0)
      {
        coeff_count = unary_bin_max_decode(dep_dp,ctx->coeff_count_context[curr_ctx_idx]+4,1,max_coeff[curr_ctx_idx]-1);
        coeff_count++;
      }
      break;
    }
    case 9:  // 8x8 ABT luma inter
    case 10: // 4x8,8x4 ABT luma inter
    case 11: // 4x4 ABT luma inter
    case 12: // 8x8 ABT luma intra
    case 13: // 4x8,8x4 ABT luma intra
    case 14: // 4x4 ABT luma intra
    {
      //context termination
      j = img->subblock_y;
      i = img->subblock_x;
      b8 = ((j>1)?2:0) + ((i>1)?1:0);
      curr_abt_ctx_idx = curr_ctx_idx - NUM_TRANS_TYPE;

      if (j==0)
        b = (currMB->mb_available[0][1] == NULL)?0:(((currMB->mb_available[0][1])->coeffs_count[BLOCK_SIZE-1][i]==0)?0:1);
      else
        b = (currMB->coeffs_count[j-1][i]==0)?0:1;

      if (i==0)
        a = (currMB->mb_available[1][0] == NULL)?0:(((currMB->mb_available[1][0])->coeffs_count[j][BLOCK_SIZE-1]==0)?0:1);
      else
        a = (currMB->coeffs_count[j][i-1]==0)?0:1;

      act_ctx = a+2*b;

      coeff_count = unary_bin_max_decodeABT(dep_dp, ctx->ABT_coeff_count_context[curr_abt_ctx_idx]+act_ctx, 4-act_ctx,
                                            max_ccnt_abt[curr_abt_ctx_idx]);
      if ((curr_ctx_idx==9)||(curr_ctx_idx==12))
        coeff_count++;

      // copy coeff_count to all 4x4 positions belonging to the current subblock
      if ((curr_ctx_idx==9)||(curr_ctx_idx==12))
      {
        currMB->coeffs_count[j  ][i  ] =
        currMB->coeffs_count[j+1][i  ] =
        currMB->coeffs_count[j  ][i+1] =
        currMB->coeffs_count[j+1][i+1] = coeff_count;
      }
      else if ((curr_ctx_idx==10)||(curr_ctx_idx==13))
      {
        currMB->coeffs_count[j][i] = coeff_count;
        if (currMB->abt_mode[b8]==B8x4)
          currMB->coeffs_count[j  ][i+1] = coeff_count;
        else
          currMB->coeffs_count[j+1][i  ] = coeff_count;
      }
      else
      {
        currMB->coeffs_count[j][i] = coeff_count;
      }

      break;
    }
    default: printf("ERROR");
    }
    coeff_count_all = coeff_count;
    if (coeff_count == 0)
      send_eob = 1;
    else
      send_eob = 0;

    prevLevel = 0;
  }

  if (send_eob == 1)
  {
    se->value1 = 0;
    se->value2 =0;
    send_eob = 0;
#if TRACE
    fprintf(p_trace, "@%d%s\t\t\t%d\t%d\n",symbolCount++, se->tracestring, se->value1,se->value2);
    fflush(p_trace);
#endif
    return;
  }
  if (coeff_count != 0) //get the coeff. (run and level)
  {
    if (curr_ctx_idx < 9)
    {
      //determine run ctx
      switch (curr_ctx_idx)
      {
      case 1:
      case 2:
        curr_run_ctx = ((coeff_count_all) >= 4) ? 1 : 0;
        break;
      case 3:
        curr_run_ctx = ((coeff_count) >= 4) ? 1 : 0;
        break;
      case 4:
      case 7:
      case 8:
      case 0:
        curr_run_ctx = ((coeff_count) >= 3) ? 1 : 0;
        break;
      case 5:
      case 6:
        curr_run_ctx = ((coeff_count) >= 2) ? 1 : 0;
        break;
      }

      curr_run_ctx = 9*curr_run_ctx + curr_ctx_idx;
      // get run
      if (coeff_count == coeff_count_all)
        max_run = max_coeff[curr_ctx_idx]-coeff_count_all;
      if (max_run > 0)
        run = unary_bin_max_decode(dep_dp,ctx->run_context[curr_run_ctx],1,max_run);
      max_run -= run;
    }
    else
    {
      curr_abt_ctx_idx = curr_ctx_idx - NUM_TRANS_TYPE;
      curr_run_ctx = ((coeff_count_all) >= 4) ? 1 : 0;
      curr_run_ctx = NUM_TRANS_TYPE_ABT*curr_run_ctx + curr_abt_ctx_idx;        // 6->NUM_TRANS_TYPE_ABT to prevent selection of wrong context.
      assert( curr_run_ctx>=0);
      //send run
      if (coeff_count == coeff_count_all)
        max_run = max_coeff_abt[curr_abt_ctx_idx]-coeff_count_all;
      if ((max_run != 0) )
        run = unary_bin_max_decodeABT(dep_dp,ctx->ABT_run_context[curr_run_ctx],1,max_run);
      max_run -= run;
    }

    //get level
    if (curr_ctx_idx < 9)
      changed_ctx_idx = curr_ctx_idx;
    else
    {                        // use JM contexts for ABT level encoding:
      if (curr_ctx_idx < 12) // inter
        changed_ctx_idx = 1;
      else                   // intra
        changed_ctx_idx = 2;
    }
    changed_ctx_idx += 9*prevLevel;
    level = unary_level_decode(dep_dp,ctx->level_context[changed_ctx_idx]);
    level++;
    prevLevel     = level > 3 ? 3 : level;  
    //get sign
    curr_level_ctx = 3;
    sign_of_level = biari_decode_symbol(dep_dp, ctx->level_context[curr_ctx_idx] + curr_level_ctx);
    if (sign_of_level) level = (-1)*level;
    coeff_count--;
    if (coeff_count == 0)
    {
      send_eob = 1;
      prevLevel = 0;
    }
    else
      send_eob = 0;
  }

  se->value1 = level;
  se->value2 = run;

#if TRACE
  fprintf(p_trace, "@%d%s\t\t\t%d\t%d\n",symbolCount++, se->tracestring, se->value1, se->value2);
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

  if(inp->Encapsulated_NAL_Payload) 
  {
    BitstreamLengthInBytes = currStream->bitstream_length = EBSPtoRBSP(code_buffer, currStream->bitstream_length, startcodeprefix_len);
    BitstreamLengthInBytes = currStream->bitstream_length = RBSPtoSODB(code_buffer, currStream->bitstream_length);
  }

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
 *
 * Unsure whether the "correction" below actually solves an off-by-one
 * problem or whether it introduces one in some cases :-(  Anyway,
 * with this change the bit stream format works with CABAC again.
 * StW, 8.7.02
 ************************************************************************
 */
int cabac_startcode_follows(struct img_par *img, struct inp_par *inp)
{
  Slice *currSlice = img->currentSlice;
// printf ("cabac_startcode_follows returns %d\n", (img->current_mb_nr == currSlice->last_mb_nr));
//  if (img->current_mb_nr == currSlice->last_mb_nr)
  if (img->current_mb_nr > currSlice->last_mb_nr)
    return TRUE;
  return FALSE;
}

/*!
 ************************************************************************
 * \brief
 *    Unary binarization for values <16 with max 4 CTX; larger values use
 *    syncword and binary representation of value-16; no terminating "0"
 *    of sync for binary representation(finite symbol alphabet)
 *    Used for ABT.
 ************************************************************************
 */
unsigned int unary_bin_max_decodeABT(DecodingEnvironmentPtr dep_dp,
                                     BiContextTypePtr ctx,
                                     int ctx_offset,
                                     unsigned int max_symbol)
{
  unsigned int l, ibits=0, i=0;
  unsigned int symbol;
  int bin=1;
  BiContextTypePtr ictx=ctx;

  symbol =  biari_decode_symbol(dep_dp, ctx );

  if (symbol==0)
  {
    return 0;
  }
  else
  {
    symbol = 0;
    ictx = ctx + ctx_offset;
    do
    {
      if ((++bin)==17)
        break;
      l =  biari_decode_symbol(dep_dp, ictx );
      symbol++;
      if (bin==2) { ictx++; }
      if (bin==4) { ictx++; }
    }
    while (l>0);

    if (bin>16)
    {
      assert (max_symbol > 15);
      l=max_symbol-15;
      do { ibits++; }
      while ( (l>>=1) > 0);
      ictx = ctx + ctx_offset + 3;
      while (i<ibits)
      {
        l=biari_decode_symbol(dep_dp, ictx);
        symbol += l*(1<<i);
        i++;
      }
    }
    return symbol;
  }
}
