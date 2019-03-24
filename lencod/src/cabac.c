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
#include <assert.h>
#include "cabac.h"

/*!
 ************************************************************************
 * \brief
 *    Allocation of contexts models for the motion info
 *    used for arithmetic encoding
 ************************************************************************
 */
MotionInfoContexts* create_contexts_MotionInfo(void)
{
  int                 j;
  MotionInfoContexts* enco_ctx;

  enco_ctx = (MotionInfoContexts*) calloc(1, sizeof(MotionInfoContexts) );
  if( enco_ctx == NULL )
    no_mem_exit("create_contexts_MotionInfo: enco_ctx");

  for (j=0; j<3; j++)
  {
    enco_ctx->mb_type_contexts[j] = (BiContextTypePtr) malloc(NUM_MB_TYPE_CTX  * sizeof( BiContextType ) );
    if( enco_ctx->mb_type_contexts[j] == NULL )
      no_mem_exit("create_contexts_MotionInfo: enco_ctx->mb_type_contexts");
  }
  for (j=0; j<2; j++)
  {
    enco_ctx->b8_type_contexts[j] = (BiContextTypePtr) malloc(NUM_B8_TYPE_CTX * sizeof( BiContextType ) );      
    if( enco_ctx->b8_type_contexts[j] == NULL ) 
      no_mem_exit("create_contexts_MotionInfo: enco_ctx->b8_type_contexts");

    enco_ctx->mv_res_contexts[j] = (BiContextTypePtr) malloc(NUM_MV_RES_CTX * sizeof( BiContextType ) );
    if( enco_ctx->mv_res_contexts[j] == NULL )
      no_mem_exit("create_contexts_MotionInfo: enco_ctx->mv_res_contexts");

    enco_ctx->ref_no_contexts[j] = (BiContextTypePtr) malloc(NUM_REF_NO_CTX * sizeof( BiContextType ) );
    if( enco_ctx->ref_no_contexts[j] == NULL )
      no_mem_exit("create_contexts_MotionInfo: enco_ctx->ref_no_contexts");
  }

  enco_ctx->delta_qp_inter_contexts = (BiContextTypePtr) malloc(NUM_DELTA_QP_CTX * sizeof( BiContextType ) );
  if( enco_ctx->delta_qp_inter_contexts == NULL )
    no_mem_exit("create_contexts_MotionInfo: enco_ctx->delta_qp_inter_contexts");

  enco_ctx->delta_qp_intra_contexts = (BiContextTypePtr) malloc(NUM_DELTA_QP_CTX * sizeof( BiContextType ) );
  if( enco_ctx->delta_qp_intra_contexts == NULL )
    no_mem_exit("create_contexts_MotionInfo: enco_ctx->delta_qp_intra_contexts");

  return enco_ctx;
}


/*!
 ************************************************************************
 * \brief
 *    Allocates of contexts models for the texture info
 *    used for arithmetic encoding
 ************************************************************************
 */
TextureInfoContexts* create_contexts_TextureInfo(void)
{
  int                   j,k;
  TextureInfoContexts*  enco_ctx;
  static const int max_ipr=9;

 

  enco_ctx = (TextureInfoContexts*) calloc(1, sizeof(TextureInfoContexts) );
  if( enco_ctx == NULL )
    no_mem_exit("create_contexts_TextureInfo: enco_ctx");

  for (j=0; j < max_ipr; j++)
  {
    enco_ctx->ipr_contexts[j] = (BiContextTypePtr) malloc(NUM_IPR_CTX  * sizeof( BiContextType ) );
    if( enco_ctx->ipr_contexts[j] == NULL )
      no_mem_exit("create_contexts_TextureInfo: enco_ctx->ipr_contexts");
  }

  for (k=0; k<2; k++)
  {
    for (j=0; j<3; j++)
    {
      enco_ctx->cbp_contexts[k][j] = (BiContextTypePtr) malloc(NUM_CBP_CTX  * sizeof( BiContextType ) );
      if( enco_ctx->cbp_contexts[k][j] == NULL )
        no_mem_exit("create_contexts_TextureInfo: enco_ctx->cbp_contexts");
    }
  }

  for (j=0; j < 4*NUM_TRANS_TYPE; j++)
  {
    enco_ctx->level_context[j] = (BiContextTypePtr) malloc(NUM_LEVEL_CTX  * sizeof( BiContextType ) );
    if( enco_ctx->level_context[j] == NULL )
      no_mem_exit("create_contexts_TextureInfo: enco_ctx->level_context");
  }

  for (j=0; j < 2*NUM_TRANS_TYPE; j++)
  {
    enco_ctx->run_context[j] = (BiContextTypePtr) malloc(NUM_RUN_CTX  * sizeof( BiContextType ) );
    if( enco_ctx->run_context[j] == NULL )
      no_mem_exit("create_contexts_TextureInfo: enco_ctx->run_context");
  }
    
  for (j=0; j < NUM_TRANS_TYPE; j++)
  {
    enco_ctx->coeff_count_context[j] = (BiContextTypePtr) malloc(NUM_COEFF_COUNT_CTX  * sizeof( BiContextType ) );
    if( enco_ctx->coeff_count_context[j] == NULL )
      no_mem_exit("create_contexts_TextureInfo: enco_ctx->coeff_count_context");
  }

  for (j=0; j < 2*NUM_TRANS_TYPE_ABT; j++)
  {
    enco_ctx->ABT_run_context[j] = (BiContextTypePtr) malloc(NUM_RUN_CTX_ABT  * sizeof( BiContextType ) );
    if( enco_ctx->ABT_run_context[j] == NULL )
      no_mem_exit("create_contexts_TextureInfo: enco_ctx->ABT_run_context");
  }

  for (j=0; j < NUM_TRANS_TYPE_ABT; j++)
  {
    enco_ctx->ABT_coeff_count_context[j] = (BiContextTypePtr) malloc(NUM_COEFF_COUNT_CTX_ABT  * sizeof( BiContextType ) );
    if( enco_ctx->ABT_coeff_count_context[j] == NULL )
      no_mem_exit("create_contexts_TextureInfo: enco_ctx->ABT_coeff_count_context");
  }

  return enco_ctx;
}


/*!
 ************************************************************************
 * \brief
 *    Initializes an array of contexts models with some pre-defined
 *    counts (ini_flag = 1) or with a flat histogram (ini_flag = 0)
 ************************************************************************
 */
void init_contexts_MotionInfo(MotionInfoContexts *enco_ctx, int ini_flag)
{

  int i,j;
  int scale_factor;
  int qp_factor;
  int ini[3];

  if ( (img->width*img->height) <=  (IMG_WIDTH * IMG_HEIGHT) ) //  format <= QCIF
    scale_factor=1;
  else
    scale_factor=2;

  qp_factor=min(max(0,img->qp-10-SHIFT_QP),21);

  for (j=0; j<3; j++)
  {
    if (ini_flag)
    {
      for (i=0; i < NUM_MB_TYPE_CTX; i++)
      {
        ini[0] = MB_TYPE_Ini[j][i][0]+(MB_TYPE_Ini[j][i][3]*qp_factor)/10;
        ini[1] = MB_TYPE_Ini[j][i][1]+(MB_TYPE_Ini[j][i][4]*qp_factor)/10;
        ini[2] = MB_TYPE_Ini[j][i][2]*scale_factor;
        biari_init_context(enco_ctx->mb_type_contexts[j] + i,ini[0],ini[1],ini[2]);
      }
    }
    else
    {
      for (i=0; i < NUM_MB_TYPE_CTX; i++)
        biari_init_context(enco_ctx->mb_type_contexts[j] + i,1,1,100);
    }
  }
  for (j=0; j<2; j++)
  {
    if (ini_flag)
    {
      for (i=0; i < NUM_B8_TYPE_CTX; i++)
        biari_init_context(enco_ctx->b8_type_contexts[j] + i,B8_TYPE_Ini[j][i][0]*scale_factor,B8_TYPE_Ini[j][i][1]*scale_factor,B8_TYPE_Ini[j][i][2]*scale_factor);
    }
    else
    {
      for (i=0; i < NUM_B8_TYPE_CTX; i++)
        biari_init_context (enco_ctx->b8_type_contexts[j] + i, 1, 1, 1000);
    }

    if (ini_flag)
    {
      for (i=0; i < NUM_MV_RES_CTX; i++)
        biari_init_context(enco_ctx->mv_res_contexts[j] + i,MV_RES_Ini[j][i][0]*scale_factor,MV_RES_Ini[j][i][1]*scale_factor,MV_RES_Ini[j][i][2]*scale_factor);
    }
    else
    {
      for (i=0; i < NUM_MV_RES_CTX; i++)
        biari_init_context(enco_ctx->mv_res_contexts[j] + i,1,1,1000);
    }

    if (ini_flag)
    {
      for (i=0; i < NUM_REF_NO_CTX; i++)
        biari_init_context(enco_ctx->ref_no_contexts[j] + i,REF_NO_Ini[i][0]*scale_factor,REF_NO_Ini[i][1]*scale_factor,REF_NO_Ini[i][2]*scale_factor);
    }
    else
    {
      for (i=0; i < NUM_REF_NO_CTX; i++)
        biari_init_context(enco_ctx->ref_no_contexts[j] + i,1,1,1000);
    }
  }

  if (ini_flag)
  {
    for (i=0; i < NUM_DELTA_QP_CTX; i++)
      biari_init_context(enco_ctx->delta_qp_inter_contexts + i,DELTA_QP_Ini[i][0]*scale_factor,DELTA_QP_Ini[i][1]*scale_factor,DELTA_QP_Ini[i][2]*scale_factor);
  }
  else
  {
    for (i=0; i < NUM_DELTA_QP_CTX; i++)
      biari_init_context(enco_ctx->delta_qp_inter_contexts + i,1,1,1000);
  }

  if (ini_flag)
  {
    for (i=0; i < NUM_DELTA_QP_CTX; i++)
      biari_init_context(enco_ctx->delta_qp_intra_contexts + i,DELTA_QP_Ini[i][0]*scale_factor,DELTA_QP_Ini[i][1]*scale_factor,DELTA_QP_Ini[i][2]*scale_factor);
  }
  else
  {
    for (i=0; i < NUM_DELTA_QP_CTX; i++)
      biari_init_context(enco_ctx->delta_qp_intra_contexts + i,1,1,1000);
  }
}

/*!
 ************************************************************************
 * \brief
 *    Initializes an array of contexts models with some pre-defined
 *    counts (ini_flag = 1) or with a flat histogram (ini_flag = 0)
 ************************************************************************
 */
void init_contexts_TextureInfo(TextureInfoContexts *enco_ctx, int ini_flag)
{

  int i,j,k;
  int scale_factor;
  int qp_factor;
  int ini[3];
  const int max_ipr=9;

  if ( (img->width*img->height) <=  (IMG_WIDTH * IMG_HEIGHT) ) //  format <= QCIF
   scale_factor=1;
  else
   scale_factor=2;

  qp_factor=min(max(0,img->qp-10-SHIFT_QP),21);

  for (j=0; j < max_ipr; j++)
  {
    if (ini_flag)
    {
      for (i=0; i < NUM_IPR_CTX; i++)
        biari_init_context(enco_ctx->ipr_contexts[j] + i,IPR_Ini[j][i][0]*scale_factor,IPR_Ini[j][i][1]*scale_factor,IPR_Ini[j][i][2]*scale_factor);
    }
    else
    {
      for (i=0; i < NUM_IPR_CTX; i++)
        biari_init_context(enco_ctx->ipr_contexts[j] + i,2,1,50);
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
          biari_init_context(enco_ctx->cbp_contexts[k][j] + i,ini[0],ini[1],ini[2]);
        }
      }
      else
      {
        for (i=0; i < NUM_CBP_CTX; i++)
          biari_init_context(enco_ctx->cbp_contexts[k][j] + i,1,1,100);
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
        biari_init_context(enco_ctx->level_context[j] + i,ini[0],ini[1],ini[2]);
      }
    }
    else
    {
      for (i=0; i < NUM_LEVEL_CTX; i++)
        biari_init_context(enco_ctx->level_context[j] + i,1,1,100);
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
        biari_init_context(enco_ctx->run_context[j] + i,ini[0],ini[1],ini[2]);
      }
    }
    else
    {
      for (i=0; i < NUM_RUN_CTX; i++)
        biari_init_context(enco_ctx->run_context[j] + i,1,1,100);
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
        biari_init_context(enco_ctx->coeff_count_context[j] + i,ini[0],ini[1],ini[2]);
      }
    }
    else
    {
      for (i=0; i < NUM_COEFF_COUNT_CTX; i++)
        biari_init_context(enco_ctx->coeff_count_context[j] + i,1,1,100);
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
        biari_init_context(enco_ctx->ABT_run_context[j] + i,ini[0],ini[1],ini[2]);
      }
    }
    else
    {
      for (i=0; i < NUM_RUN_CTX_ABT; i++)
        biari_init_context(enco_ctx->ABT_run_context[j] + i,1,1,INICNT_ABT);
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
        biari_init_context(enco_ctx->ABT_coeff_count_context[j] + i,ini[0],ini[1],ini[2]);
      }
    }
    else
    {
      for (i=0; i < NUM_COEFF_COUNT_CTX_ABT; i++)
        biari_init_context(enco_ctx->ABT_coeff_count_context[j] + i,1,1,INICNT_ABT);
    }
  }
}


/*!
 ************************************************************************
 * \brief
 *    Frees the memory of the contexts models
 *    used for arithmetic encoding of the motion info.
 ************************************************************************
 */
void delete_contexts_MotionInfo(MotionInfoContexts *enco_ctx)
{
  int j;

  if( enco_ctx == NULL )
    return;

  for (j=0; j<3; j++)
  {
    if (enco_ctx->mb_type_contexts[j] != NULL)
      free(enco_ctx->mb_type_contexts[j] );
  }

  for (j=0; j<2; j++)
  {
    if (enco_ctx->b8_type_contexts[j] != NULL)
      free(enco_ctx->b8_type_contexts[j]);

    if (enco_ctx->mv_res_contexts[j]  != NULL)
      free(enco_ctx->mv_res_contexts [j]);

    if (enco_ctx->ref_no_contexts[j]  != NULL)
      free(enco_ctx->ref_no_contexts [j]);
  }

  if (enco_ctx->delta_qp_inter_contexts != NULL)
    free(enco_ctx->delta_qp_inter_contexts);

  if (enco_ctx->delta_qp_intra_contexts != NULL)
    free(enco_ctx->delta_qp_intra_contexts);

  free( enco_ctx );

  return;
}

/*!
 ************************************************************************
 * \brief
 *    Frees the memory of the contexts models
 *    used for arithmetic encoding of the texture info.
 ************************************************************************
 */
void delete_contexts_TextureInfo(TextureInfoContexts *enco_ctx)
{

  int j,k;

  static const int max_ipr=9;

  if( enco_ctx == NULL )
    return;

  for (j=0; j < max_ipr; j++)
  {
    if (enco_ctx->ipr_contexts[j] != NULL)
      free(enco_ctx->ipr_contexts[j]);
  }

  for (k=0; k<2; k++)
    for (j=0; j<3; j++)
    {
      if (enco_ctx->cbp_contexts[k][j] != NULL)
        free(enco_ctx->cbp_contexts[k][j]);
    }

  for (j=0; j < 4*NUM_TRANS_TYPE; j++)
  {
    if (enco_ctx->level_context[j] != NULL)
      free(enco_ctx->level_context[j]);
  }

  for (j=0; j < 2*NUM_TRANS_TYPE; j++)
  {
    if (enco_ctx->run_context[j] != NULL)
      free(enco_ctx->run_context[j]);
  }
      
  for (j=0; j < NUM_TRANS_TYPE; j++)
  {
    if (enco_ctx->coeff_count_context[j] != NULL)
      free(enco_ctx->coeff_count_context[j]);
  }

  for (j=0; j < 2*NUM_TRANS_TYPE_ABT; j++)
  {
    if (enco_ctx->ABT_run_context[j] != NULL)
      free(enco_ctx->ABT_run_context[j]);
  }

  for (j=0; j < NUM_TRANS_TYPE_ABT; j++)
  {
    if (enco_ctx->ABT_coeff_count_context[j] != NULL)
      free(enco_ctx->ABT_coeff_count_context[j]);
  }
  free( enco_ctx );

  return;

}


/*!
 **************************************************************************
 * \brief
 *    generates arithmetic code and passes the code to the buffer
 **************************************************************************
 */
int writeSyntaxElement_CABAC(SyntaxElement *se, DataPartition *this_dataPart)
{
  int curr_len;
  EncodingEnvironmentPtr eep_dp = &(this_dataPart->ee_cabac);

  curr_len = arienco_bits_written(eep_dp);

  // perform the actual coding by calling the appropriate method
  se->writing(se, eep_dp);

  if(se->type != SE_HEADER)
    this_dataPart->bitstream->write_flag = 1;

  return (se->len = (arienco_bits_written(eep_dp) - curr_len));
}


/*!
 ***************************************************************************
 * \brief
 *    This function is used to arithmetically encode the macroblock
 *    type info of a given MB.
 ***************************************************************************
 */
void writeMB_typeInfo2Buffer_CABAC(SyntaxElement *se, EncodingEnvironmentPtr eep_dp)
{
  int a, b;
  int act_ctx;
  int act_sym;
  int csym;
  int bframe   = (img->type==B_IMG);
  int mode_sym = 0;
  int mode16x16;

  MotionInfoContexts *ctx         = (img->currentSlice)->mot_ctx;
  Macroblock         *currMB      = &img->mb_data[img->current_mb_nr];
  int                curr_mb_type = se->value1;

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

    act_ctx     = a + b;
    act_sym     = curr_mb_type;
    se->context = act_ctx; // store context

    if (act_sym==0) // 4x4 Intra
    {
      biari_encode_symbol(eep_dp, 0, ctx->mb_type_contexts[0] + act_ctx );
    }
    else // 16x16 Intra
    {
      biari_encode_symbol(eep_dp, 1, ctx->mb_type_contexts[0] + act_ctx );
      mode_sym = act_sym-1; // Values in the range of 0...23
      act_ctx  = 4;
      act_sym  = mode_sym/12;
      biari_encode_symbol(eep_dp, (unsigned char) act_sym, ctx->mb_type_contexts[0] + act_ctx ); // coding of AC/no AC
      mode_sym = mode_sym % 12;
      act_sym  = mode_sym / 4; // coding of cbp: 0,1,2
      act_ctx  = 5;
      if (act_sym==0)
      {
        biari_encode_symbol(eep_dp, 0, ctx->mb_type_contexts[0] + act_ctx );
      }
      else
      {
        biari_encode_symbol(eep_dp, 1, ctx->mb_type_contexts[0] + act_ctx );
        act_ctx=6;
        if (act_sym==1)
        {
          biari_encode_symbol(eep_dp, 0, ctx->mb_type_contexts[0] + act_ctx );
        }
        else
        {
          biari_encode_symbol(eep_dp, 1, ctx->mb_type_contexts[0] + act_ctx );
        }
      }
      mode_sym = mode_sym % 4; // coding of I pred-mode: 0,1,2,3
      act_sym  = mode_sym/2;
      act_ctx  = 7;
      biari_encode_symbol(eep_dp, (unsigned char) act_sym, ctx->mb_type_contexts[0] + act_ctx );
      act_ctx  = 8;
      act_sym  = mode_sym%2;
      biari_encode_symbol(eep_dp, (unsigned char) act_sym, ctx->mb_type_contexts[0] + act_ctx );
    }
  }
  else // INTER
  {
    if (currMB->mb_available[0][1] == NULL)
      b = 0;
    else
      b = (( (currMB->mb_available[0][1])->mb_type != 0) ? 1 : 0 );
    if (currMB->mb_available[1][0] == NULL)
      a = 0;
    else
      a = (( (currMB->mb_available[1][0])->mb_type != 0) ? 1 : 0 );

    act_sym = curr_mb_type;

    if (act_sym>=(mode16x16=(bframe?24:7)))
    {
      mode_sym = act_sym-mode16x16;
      act_sym  = mode16x16; // 16x16 mode info
    }


    act_ctx = a + b;
    se->context = act_ctx; // store context

    if (!bframe)
    {
      switch (act_sym)
      {
      case 0:
        biari_encode_symbol (eep_dp, 0, &ctx->mb_type_contexts[1][act_ctx]);
        break;
      case 1:
        biari_encode_symbol (eep_dp, 1, &ctx->mb_type_contexts[1][act_ctx]);
        biari_encode_symbol (eep_dp, 0, &ctx->mb_type_contexts[1][4]);
        biari_encode_symbol (eep_dp, 0, &ctx->mb_type_contexts[1][5]);
        biari_encode_symbol (eep_dp, 0, &ctx->mb_type_contexts[1][6]);
        break;
      case 2:
        biari_encode_symbol (eep_dp, 1, &ctx->mb_type_contexts[1][act_ctx]);
        biari_encode_symbol (eep_dp, 0, &ctx->mb_type_contexts[1][4]);
        biari_encode_symbol (eep_dp, 1, &ctx->mb_type_contexts[1][5]);
        biari_encode_symbol (eep_dp, 1, &ctx->mb_type_contexts[1][7]);
        break;
      case 3:
        biari_encode_symbol (eep_dp, 1, &ctx->mb_type_contexts[1][act_ctx]);
        biari_encode_symbol (eep_dp, 0, &ctx->mb_type_contexts[1][4]);
        biari_encode_symbol (eep_dp, 1, &ctx->mb_type_contexts[1][5]);
        biari_encode_symbol (eep_dp, 0, &ctx->mb_type_contexts[1][7]);
        break;
      case 4:
      case 5:
        biari_encode_symbol (eep_dp, 1, &ctx->mb_type_contexts[1][act_ctx]);
        biari_encode_symbol (eep_dp, 0, &ctx->mb_type_contexts[1][4]);
        biari_encode_symbol (eep_dp, 0, &ctx->mb_type_contexts[1][5]);
        biari_encode_symbol (eep_dp, 1, &ctx->mb_type_contexts[1][6]);
        break;
      case 6:
        biari_encode_symbol (eep_dp, 1, &ctx->mb_type_contexts[1][act_ctx]);
        biari_encode_symbol (eep_dp, 1, &ctx->mb_type_contexts[1][4]);
        biari_encode_symbol (eep_dp, 0, &ctx->mb_type_contexts[1][7]);
        break;
      case 7:
        biari_encode_symbol (eep_dp, 1, &ctx->mb_type_contexts[1][act_ctx]);
        biari_encode_symbol (eep_dp, 1, &ctx->mb_type_contexts[1][4]);
        biari_encode_symbol (eep_dp, 1, &ctx->mb_type_contexts[1][7]);
        break;
      default:
        printf ("Unsupported MB-MODE in writeMB_typeInfo2Buffer_CABAC!\n");
        exit (1);
      }
    }
    else //===== B-FRAMES =====
    {
      if (act_sym==0)
      {
        biari_encode_symbol (eep_dp, 0, &ctx->mb_type_contexts[2][act_ctx]);
      }
      else if (act_sym<=2)
      {
        biari_encode_symbol (eep_dp, 1, &ctx->mb_type_contexts[2][act_ctx]);
        biari_encode_symbol (eep_dp, 0, &ctx->mb_type_contexts[2][4]);
        csym=act_sym-1;
        if (csym) biari_encode_symbol (eep_dp, 1, &ctx->mb_type_contexts[2][6]);
        else      biari_encode_symbol (eep_dp, 0, &ctx->mb_type_contexts[2][6]);
      }
      else if (act_sym<=10)
      {
        biari_encode_symbol (eep_dp, 1, &ctx->mb_type_contexts[2][act_ctx]);
        biari_encode_symbol (eep_dp, 1, &ctx->mb_type_contexts[2][4]);
        biari_encode_symbol (eep_dp, 0, &ctx->mb_type_contexts[2][5]);
        csym=(((act_sym-3)>>2)&0x01);
        if (csym) biari_encode_symbol (eep_dp, 1, &ctx->mb_type_contexts[2][6]);
        else      biari_encode_symbol (eep_dp, 0, &ctx->mb_type_contexts[2][6]);
        csym=(((act_sym-3)>>1)&0x01);
        if (csym) biari_encode_symbol (eep_dp, 1, &ctx->mb_type_contexts[2][6]);
        else      biari_encode_symbol (eep_dp, 0, &ctx->mb_type_contexts[2][6]);
        csym=((act_sym-3)&0x01);
        if (csym) biari_encode_symbol (eep_dp, 1, &ctx->mb_type_contexts[2][6]);
        else      biari_encode_symbol (eep_dp, 0, &ctx->mb_type_contexts[2][6]);
      }
      else if (act_sym==11 || act_sym==22)
      {
        biari_encode_symbol (eep_dp, 1, &ctx->mb_type_contexts[2][act_ctx]);
        biari_encode_symbol (eep_dp, 1, &ctx->mb_type_contexts[2][4]);
        biari_encode_symbol (eep_dp, 1, &ctx->mb_type_contexts[2][5]);
        biari_encode_symbol (eep_dp, 1, &ctx->mb_type_contexts[2][6]);
        biari_encode_symbol (eep_dp, 1, &ctx->mb_type_contexts[2][6]);
        if (act_sym==11)  biari_encode_symbol (eep_dp, 0, &ctx->mb_type_contexts[2][6]);
        else              biari_encode_symbol (eep_dp, 1, &ctx->mb_type_contexts[2][6]);
      }
      else
      {
        if (act_sym > 22) act_sym--;
        biari_encode_symbol (eep_dp, 1, &ctx->mb_type_contexts[2][act_ctx]);
        biari_encode_symbol (eep_dp, 1, &ctx->mb_type_contexts[2][4]);
        biari_encode_symbol (eep_dp, 1, &ctx->mb_type_contexts[2][5]);
        csym=(((act_sym-12)>>3)&0x01);
        if (csym) biari_encode_symbol (eep_dp, 1, &ctx->mb_type_contexts[2][6]);
        else      biari_encode_symbol (eep_dp, 0, &ctx->mb_type_contexts[2][6]);
        csym=(((act_sym-12)>>2)&0x01);
        if (csym) biari_encode_symbol (eep_dp, 1, &ctx->mb_type_contexts[2][6]);
        else      biari_encode_symbol (eep_dp, 0, &ctx->mb_type_contexts[2][6]);
        csym=(((act_sym-12)>>1)&0x01);
        if (csym) biari_encode_symbol (eep_dp, 1, &ctx->mb_type_contexts[2][6]);
        else      biari_encode_symbol (eep_dp, 0, &ctx->mb_type_contexts[2][6]);
        csym=((act_sym-12)&0x01);
        if (csym) biari_encode_symbol (eep_dp, 1, &ctx->mb_type_contexts[2][6]);
        else      biari_encode_symbol (eep_dp, 0, &ctx->mb_type_contexts[2][6]);
        if (act_sym >=22) act_sym++;
      }
    }

    if(act_sym==mode16x16) // additional info for 16x16 Intra-mode
    {
      act_ctx = 8;
      act_sym = mode_sym/12;
      biari_encode_symbol(eep_dp, (unsigned char) act_sym, ctx->mb_type_contexts[1] + act_ctx ); // coding of AC/no AC
      mode_sym = mode_sym % 12;

      act_sym = mode_sym / 4; // coding of cbp: 0,1,2
      act_ctx = 9;
      if (act_sym==0)
      {
        biari_encode_symbol(eep_dp, 0, ctx->mb_type_contexts[1] + act_ctx );
      }
      else
      {
        biari_encode_symbol(eep_dp, 1, ctx->mb_type_contexts[1] + act_ctx );
        if (act_sym==1)
        {
          biari_encode_symbol(eep_dp, 0, ctx->mb_type_contexts[1] + act_ctx );
        }
        else
        {
          biari_encode_symbol(eep_dp, 1, ctx->mb_type_contexts[1] + act_ctx );
        }
      }

      mode_sym = mode_sym % 4; // coding of I pred-mode: 0,1,2,3
      act_ctx  = 10;
      act_sym  = mode_sym/2;
      biari_encode_symbol(eep_dp, (unsigned char) act_sym, ctx->mb_type_contexts[1] + act_ctx );
      act_sym  = mode_sym%2;
      biari_encode_symbol(eep_dp, (unsigned char) act_sym, ctx->mb_type_contexts[1] + act_ctx );
    }
  }
}


/*!
 ***************************************************************************
 * \brief
 *    This function is used to arithmetically encode the 8x8 block
 *    type info
 ***************************************************************************
 */
void writeB8_typeInfo2Buffer_CABAC(SyntaxElement *se, EncodingEnvironmentPtr eep_dp)
{
  int act_ctx;
  int act_sym, csym;
  int bframe=(img->type==B_IMG);

  MotionInfoContexts *ctx    = (img->currentSlice)->mot_ctx;

  act_sym = se->value1;
  act_ctx = 0;

  if (!bframe)
  {
    switch (act_sym)
    {
    case 0:
      biari_encode_symbol (eep_dp, 1, &ctx->b8_type_contexts[0][1]);
      break;
    case 1:
      biari_encode_symbol (eep_dp, 0, &ctx->b8_type_contexts[0][1]);
      biari_encode_symbol (eep_dp, 0, &ctx->b8_type_contexts[0][2]);
      biari_encode_symbol (eep_dp, 0, &ctx->b8_type_contexts[0][3]);
      break;
    case 2:
      biari_encode_symbol (eep_dp, 0, &ctx->b8_type_contexts[0][1]);
      biari_encode_symbol (eep_dp, 0, &ctx->b8_type_contexts[0][2]);
      biari_encode_symbol (eep_dp, 1, &ctx->b8_type_contexts[0][3]);
      biari_encode_symbol (eep_dp, 1, &ctx->b8_type_contexts[0][4]);
      break;
    case 3:
      biari_encode_symbol (eep_dp, 0, &ctx->b8_type_contexts[0][1]);
      biari_encode_symbol (eep_dp, 0, &ctx->b8_type_contexts[0][2]);
      biari_encode_symbol (eep_dp, 1, &ctx->b8_type_contexts[0][3]);
      biari_encode_symbol (eep_dp, 0, &ctx->b8_type_contexts[0][4]);
      break;
    case 4:
      biari_encode_symbol (eep_dp, 0, &ctx->b8_type_contexts[0][1]);
      biari_encode_symbol (eep_dp, 1, &ctx->b8_type_contexts[0][2]);
      break;
    }
  }
  else //===== B-FRAME =====
  {
    if (act_sym==0)
    {
      biari_encode_symbol (eep_dp, 0, &ctx->b8_type_contexts[1][0]);
      return;
    }
    else
    {
      biari_encode_symbol (eep_dp, 1, &ctx->b8_type_contexts[1][0]);
      act_sym--;
    }
    if (act_sym<2)
    {
      biari_encode_symbol (eep_dp, 0, &ctx->b8_type_contexts[1][1]);
      if (act_sym==0)   biari_encode_symbol (eep_dp, 0, &ctx->b8_type_contexts[1][3]);
      else              biari_encode_symbol (eep_dp, 1, &ctx->b8_type_contexts[1][3]);
    }
    else if (act_sym<6)
    {
      biari_encode_symbol (eep_dp, 1, &ctx->b8_type_contexts[1][1]);
      biari_encode_symbol (eep_dp, 0, &ctx->b8_type_contexts[1][2]);
      csym=(((act_sym-2)>>1)&0x01);
      if (csym) biari_encode_symbol (eep_dp, 1, &ctx->b8_type_contexts[1][3]);
      else      biari_encode_symbol (eep_dp, 0, &ctx->b8_type_contexts[1][3]);
      csym=((act_sym-2)&0x01);
      if (csym) biari_encode_symbol (eep_dp, 1, &ctx->b8_type_contexts[1][3]);
      else      biari_encode_symbol (eep_dp, 0, &ctx->b8_type_contexts[1][3]);
    }
    else if (act_sym==12)
    {
      biari_encode_symbol (eep_dp, 1, &ctx->b8_type_contexts[1][1]);
      biari_encode_symbol (eep_dp, 1, &ctx->b8_type_contexts[1][2]);
      biari_encode_symbol (eep_dp, 1, &ctx->b8_type_contexts[1][3]);
      biari_encode_symbol (eep_dp, 1, &ctx->b8_type_contexts[1][3]);
    }
    else
    {
      biari_encode_symbol (eep_dp, 1, &ctx->b8_type_contexts[1][1]);
      biari_encode_symbol (eep_dp, 1, &ctx->b8_type_contexts[1][2]);
      csym=(((act_sym-6)>>2)&0x01);
      if (csym) biari_encode_symbol (eep_dp, 1, &ctx->b8_type_contexts[1][3]);
      else      biari_encode_symbol (eep_dp, 0, &ctx->b8_type_contexts[1][3]);
      csym=(((act_sym-6)>>1)&0x01);
      if (csym) biari_encode_symbol (eep_dp, 1, &ctx->b8_type_contexts[1][3]);
      else      biari_encode_symbol (eep_dp, 0, &ctx->b8_type_contexts[1][3]);
      csym=((act_sym-6)&0x01);
      if (csym) biari_encode_symbol (eep_dp, 1, &ctx->b8_type_contexts[1][3]);
      else      biari_encode_symbol (eep_dp, 0, &ctx->b8_type_contexts[1][3]);
    }
  }
}


/*!
 ****************************************************************************
 * \brief
 *    This function is used to arithmetically encode a pair of
 *    intra prediction modes of a given MB.
 ****************************************************************************
 */
void writeIntraPredMode2Buffer_CABAC(SyntaxElement *se, EncodingEnvironmentPtr eep_dp)
{
  static const int    right[8] = {0, 0, 1, 1, 0, 0, 1, 1};
  Macroblock          *currMB  = &(img->mb_data[img->current_mb_nr]);
  TextureInfoContexts *ctx     = img->currentSlice->tex_ctx;
  int                 prev_sym; 
 
  //--- first symbol ---
  if (right[se->context/2])             prev_sym = currMB->intra_pred_modes[se->context-3];
  else if (currMB->mb_available[1][0])  prev_sym = currMB->mb_available[1][0]->intra_pred_modes[se->context+5];
  else                                  prev_sym = 0;

  unary_bin_max_encode(eep_dp,(unsigned int) se->value1,ctx->ipr_contexts[prev_sym],1,8);

  //--- second symbol ---
  prev_sym = se->value1;
  unary_bin_max_encode(eep_dp,(unsigned int) se->value2,ctx->ipr_contexts[prev_sym],1,8);
}


/*!
 ****************************************************************************
 * \brief
 *    This function is used to arithmetically encode the reference
 *    parameter of a given MB.
 ****************************************************************************
 */
void writeRefFrame2Buffer_CABAC(SyntaxElement *se, EncodingEnvironmentPtr eep_dp)
{
  MotionInfoContexts  *ctx    = img->currentSlice->mot_ctx;
  Macroblock          *currMB = &img->mb_data[img->current_mb_nr];
  int                 addctx  = se->context;

  int   a, b;
  int   act_ctx;
  int   act_sym;
  int** refframe_array = (img->type==B_IMG ? fw_refFrArr : refFrArr);

  if (currMB->mb_available[0][1] == NULL)
    b = 0;
  else
    b = (refframe_array[img->block_y+img->subblock_y-1][img->block_x+img->subblock_x] > 0 ? 1 : 0);
  if (currMB->mb_available[1][0] == NULL)
    a = 0;
  else 
    a = (refframe_array[img->block_y+img->subblock_y][img->block_x+img->subblock_x-1] > 0 ? 1 : 0);

  act_ctx     = a + 2*b;
  se->context = act_ctx; // store context
  act_sym     = se->value1;

  if (act_sym==0)
  {
    biari_encode_symbol(eep_dp, 0, ctx->ref_no_contexts[addctx] + act_ctx );
  }
  else
  {
    biari_encode_symbol(eep_dp, 1, ctx->ref_no_contexts[addctx] + act_ctx);
    act_sym--;
    act_ctx=4;
    unary_bin_encode(eep_dp, act_sym,ctx->ref_no_contexts[addctx]+act_ctx,1);
  }
}


/*!
 ****************************************************************************
 * \brief
 *    This function is used to arithmetically encode the motion
 *    vector data of a given MB.
 ****************************************************************************
 */
void writeMVD2Buffer_CABAC(SyntaxElement *se, EncodingEnvironmentPtr eep_dp)
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

  MotionInfoContexts  *ctx         = img->currentSlice->mot_ctx;
  Macroblock          *currMB      = &img->mb_data[img->current_mb_nr];

  if (j==0)
  {
    if (currMB->mb_available[0][1] == NULL)
      b = 0;
    else 
      b = absm((currMB->mb_available[0][1])->mvd[0][BLOCK_SIZE-1][i][k]);
  }
  else
    b = absm(currMB->mvd[0][j-1][i][k]);
          
  if (i==0)
  {
    if (currMB->mb_available[1][0] == NULL)
      a = 0;
    else 
      a = absm((currMB->mb_available[1][0])->mvd[0][j][BLOCK_SIZE-1][k]);
  }
  else
    a = absm(currMB->mvd[0][j][i-1][k]);

  if ((mv_local_err=a+b)<3)
    act_ctx = 5*k;
  else
  {
    if (mv_local_err>32)
      act_ctx=5*k+3;
    else
      act_ctx=5*k+2;
  }
  mv_pred_res = se->value1;
  se->context = act_ctx;

  act_sym = absm(mv_pred_res);

  if (act_sym == 0)
    biari_encode_symbol(eep_dp, 0, &ctx->mv_res_contexts[0][act_ctx] );
  else
  {
    biari_encode_symbol(eep_dp, 1, &ctx->mv_res_contexts[0][act_ctx] );
    mv_sign = (mv_pred_res<0) ? 1: 0;
    act_ctx=5*k+4;
    biari_encode_symbol(eep_dp, (unsigned char) mv_sign, &ctx->mv_res_contexts[1][act_ctx] );
    act_sym--;
    act_ctx=5*k;
    unary_mv_encode(eep_dp,act_sym,ctx->mv_res_contexts[1]+act_ctx,3);
  }
}


/*!
 ****************************************************************************
 * \brief
 *    This function is used to arithmetically encode the coded
 *    block pattern of a given delta quant.
 ****************************************************************************
 */
void writeDquant_inter_CABAC(SyntaxElement *se, EncodingEnvironmentPtr eep_dp)
{
  MotionInfoContexts *ctx = img->currentSlice->mot_ctx;
  Macroblock *currMB = &img->mb_data[img->current_mb_nr];

  int act_ctx;
  int act_sym;
  int dquant = se->value1;
  int sign=0;

  if (dquant <= 0)
    sign = 1;
  act_sym = abs(dquant) << 1;

  act_sym += sign;
  act_sym --;

  if (currMB->mb_available[1][0] == NULL)
    act_ctx = 0;
  else
    act_ctx = ( ((currMB->mb_available[1][0])->delta_qp != 0) ? 1 : 0);

  if (act_sym==0)
  {
    biari_encode_symbol(eep_dp, 0, ctx->delta_qp_inter_contexts + act_ctx );
  }
  else
  {
    biari_encode_symbol(eep_dp, 1, ctx->delta_qp_inter_contexts + act_ctx);
    act_ctx=2;
    act_sym--;
    unary_bin_encode(eep_dp, act_sym,ctx->delta_qp_inter_contexts+act_ctx,1);
  }
}


/*!
 ****************************************************************************
 * \brief
 *    This function is used to arithmetically encode the coded
 *    block pattern of a given delta quant.
 ****************************************************************************
 */
void writeDquant_intra_CABAC(SyntaxElement *se, EncodingEnvironmentPtr eep_dp)
{
  MotionInfoContexts *ctx = img->currentSlice->mot_ctx;
  Macroblock *currMB = &img->mb_data[img->current_mb_nr];

  int act_ctx;
  int act_sym;
  int dquant = se->value1;
  int sign=0;

  if (dquant <= 0)
    sign = 1;
  act_sym = abs(dquant) << 1;

  act_sym += sign;
  act_sym --;

  if (currMB->mb_available[1][0] == NULL)
    act_ctx = 0;
  else
    act_ctx = ( ((currMB->mb_available[1][0])->delta_qp != 0) ? 1 : 0);

  if (act_sym==0)
  {
    biari_encode_symbol(eep_dp, 0, ctx->delta_qp_intra_contexts + act_ctx );
  }
  else
  {
    biari_encode_symbol(eep_dp, 1, ctx->delta_qp_intra_contexts + act_ctx);
    act_ctx=2;
    act_sym--;
    unary_bin_encode(eep_dp, act_sym,ctx->delta_qp_intra_contexts+act_ctx,1);
  }
}


/*!
 ****************************************************************************
 * \brief
 *    This function is used to arithmetically encode the motion
 *    vector data of a B-frame MB.
 ****************************************************************************
 */
void writeBiMVD2Buffer_CABAC(SyntaxElement *se, EncodingEnvironmentPtr eep_dp)
{
  int i = img->subblock_x;
  int j = img->subblock_y;
  int a, b;
  int act_ctx;
  int act_sym;
  int mv_pred_res;
  int mv_local_err;
  int mv_sign;
  int backward = se->value2 & 0x01;
  int k = (se->value2>>1); // MVD component

  MotionInfoContexts  *ctx    = img->currentSlice->mot_ctx;
  Macroblock          *currMB = &img->mb_data[img->current_mb_nr];

  if (j==0)
  {
    if (currMB->mb_available[0][1] == NULL)
      b = 0;
    else
      b = absm((currMB->mb_available[0][1])->mvd[backward][BLOCK_SIZE-1][i][k]);
  }
  else
    b = absm(currMB->mvd[backward][j-1][i][k]);

  if (i==0)
  {
    if (currMB->mb_available[1][0] == NULL)
      a = 0;
    else
      a = absm((currMB->mb_available[1][0])->mvd[backward][j][BLOCK_SIZE-1][k]);
  }
  else
    a = absm(currMB->mvd[backward][j][i-1][k]);

  if ((mv_local_err=a+b)<3)
    act_ctx = 5*k;
  else
  {
    if (mv_local_err>32)
      act_ctx=5*k+3;
    else
      act_ctx=5*k+2;
  }
  mv_pred_res = se->value1;
  se->context = act_ctx;

  act_sym = absm(mv_pred_res);

  if (act_sym == 0)
    biari_encode_symbol(eep_dp, 0, &ctx->mv_res_contexts[0][act_ctx] );
  else
  {
    biari_encode_symbol(eep_dp, 1, &ctx->mv_res_contexts[0][act_ctx] );
    mv_sign = (mv_pred_res<0) ? 1: 0;
    act_ctx=5*k+4;
    biari_encode_symbol(eep_dp, (unsigned char) mv_sign, &ctx->mv_res_contexts[1][act_ctx] );
    act_sym--;
    act_ctx=5*k;
    unary_mv_encode(eep_dp,act_sym,ctx->mv_res_contexts[1]+act_ctx,3);
  }
}


/*!
 ****************************************************************************
 * \brief
 *    This function is used to arithmetically encode the coded
 *    block pattern of an 8x8 block
 ****************************************************************************
 */
void writeCBP_BIT_CABAC (int b8, int bit, int cbp, Macroblock* currMB, int inter, EncodingEnvironmentPtr eep_dp)
{
  int a, b;

  //===== GET CONTEXT FOR CBP-BIT =====
  if (b8/2 == 0) // upper block is in upper macroblock
  {
    if (currMB->mb_available[0][1] == NULL)
      b = 0;
    else
      b = ((currMB->mb_available[0][1]->cbp & (1<<(b8+2))) == 0 ? 1 : 0);
  }
  else
    b   = ((cbp & (1<<(b8-2))) == 0 ? 1: 0);
  if (b8%2 == 0) // left block is in left macroblock
  {
    if (currMB->mb_available[1][0] == NULL)
      a = 0;
    else
      a = ((currMB->mb_available[1][0]->cbp & (1<<(b8+1))) == 0 ? 1 : 0);
  }
  else
    a   = ((cbp & (1<<(b8-1))) == 0 ? 1: 0);

  //===== WRITE BIT =====
  biari_encode_symbol (eep_dp, (unsigned char) bit,
                       img->currentSlice->tex_ctx->cbp_contexts[inter][0] + a+2*b);
}

/*!
 ****************************************************************************
 * \brief
 *    This function is used to arithmetically encode the coded
 *    block pattern of a macroblock
 ****************************************************************************
 */
void writeCBP2Buffer_CABAC(SyntaxElement *se, EncodingEnvironmentPtr eep_dp)
{
  TextureInfoContexts *ctx = img->currentSlice->tex_ctx;
  Macroblock *currMB = &img->mb_data[img->current_mb_nr];

  int a, b;
  int curr_cbp_ctx, curr_cbp_idx;
  int cbp = se->value1; // symbol to encode
  int cbp_bit;
  int b8;

  for (b8=0; b8<4; b8++)
  {
    curr_cbp_idx = (currMB->b8mode[b8] == IBLOCK ? 0 : 1);
    writeCBP_BIT_CABAC (b8, cbp&(1<<b8), cbp, currMB, curr_cbp_idx, eep_dp);
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
  cbp_bit = (cbp > 15 ) ? 1 : 0;
  biari_encode_symbol(eep_dp, (unsigned char) cbp_bit, ctx->cbp_contexts[curr_cbp_idx][1] + curr_cbp_ctx );

  if (cbp > 15)
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
    cbp_bit = ((cbp>>4) == 2) ? 1 : 0;
    biari_encode_symbol(eep_dp, (unsigned char) cbp_bit, ctx->cbp_contexts[curr_cbp_idx][2] + curr_cbp_ctx );
  }
}


/*!
 ****************************************************************************
 * \brief
 *    This function is used to arithmetically encode coeff_count, level and
 *    run of a given MB.
 ****************************************************************************
 */
void writeRunLevel2Buffer_CABAC(SyntaxElement *se, EncodingEnvironmentPtr eep_dp)
{
  const int level = se->value1;
  const int run = se->value2;
  int curr_ctx_idx = se->context;
  int sign_of_level;
  int max_run = 0;
  static int coeff_can[2][65]; // was initialized with zeros. init necessary? is only initialized at first block of whole seq.
  static int coeff_count=0;
  int i,j,a,b;
  int max_coeff[9] = {8,16,16,16,15,4,4,15,15};
  int curr_run_ctx = 0, curr_level_ctx;
  int act_ctx=0;
  static int count[3] = {0,0,0};
  static int prev_sym[3] = {0,0,0};
  int prevLevel = 0;
  int absLevel;
  int changed_ctx_idx=curr_ctx_idx;

  // ABT stuff
  int b8, curr_abt_ctx_idx=-1; // =-1: avoid warnings
  static unsigned int max_ccnt_abt  [6] = {63,32,16,63,32,16};
  static unsigned int max_coeff_abt [6] = {64,32,16,64,32,16};

  TextureInfoContexts *ctx = img->currentSlice->tex_ctx;
  Macroblock *currMB = &img->mb_data[img->current_mb_nr];

  if (level != 0)
  {
    coeff_can[0][coeff_count] = level;
    coeff_can[1][coeff_count] = run;
    coeff_count++;
  }
  else
  {
    //send coeff_count
    switch(curr_ctx_idx)
    {
    case 0://max 8 coeffs, double_scan
    {
      //context determination
      act_ctx = (count[0] == 0) ? 0 : ((prev_sym[0] == 0) ? 1 : 2);
      prev_sym[0] = coeff_count;
      if(++count[0] == 2) count[0]=0;
      
      if (coeff_count == 0)
        biari_encode_symbol(eep_dp, 0, ctx->coeff_count_context[curr_ctx_idx] + act_ctx);
      else
      {
        biari_encode_symbol(eep_dp, 1, ctx->coeff_count_context[curr_ctx_idx] + act_ctx);
        unary_bin_max_encode(eep_dp,(unsigned int) coeff_count-1,ctx->coeff_count_context[curr_ctx_idx]+4,1,max_coeff[curr_ctx_idx]-1);
      }
      break;
    }
    case 1: //max 16 coeffs //single scan
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
      if (coeff_count == 0)
        biari_encode_symbol(eep_dp, 0, ctx->coeff_count_context[curr_ctx_idx] + act_ctx);
      else
      {
        biari_encode_symbol(eep_dp, 1, ctx->coeff_count_context[curr_ctx_idx] + act_ctx);
        unary_bin_max_encode(eep_dp,(unsigned int) coeff_count-1,ctx->coeff_count_context[curr_ctx_idx]+4,1,max_coeff[curr_ctx_idx]-1);
      }
      currMB->coeffs_count[j][i] = coeff_count; 
      break;
    }
    case 4://max 15 coeffs // 16x16 luma ac
    {
      //context determination
      act_ctx = (count[2] == 0) ? 0 : ((prev_sym[2] == 0) ? 1 : 2);
      prev_sym[2] = coeff_count;
      if(++count[2] == 16) count[2]=0;

      if (coeff_count == 0)
        biari_encode_symbol(eep_dp, 0, ctx->coeff_count_context[curr_ctx_idx] + act_ctx);
      else
      {
        biari_encode_symbol(eep_dp, 1, ctx->coeff_count_context[curr_ctx_idx] + act_ctx);
        unary_bin_max_encode(eep_dp,(unsigned int) coeff_count-1,ctx->coeff_count_context[curr_ctx_idx]+4,1,max_coeff[curr_ctx_idx]-1);
      }
      break;
    }
    case 5:
    case 6://max 4 coeffs, chroma dc different ctx uv for first bin
    {
      act_ctx = (se->k == 0)?0:1;
      if (coeff_count == 0)
        biari_encode_symbol(eep_dp, 0, ctx->coeff_count_context[curr_ctx_idx] + act_ctx);
      else
      {
        biari_encode_symbol(eep_dp, 1, ctx->coeff_count_context[curr_ctx_idx] + act_ctx);
        unary_bin_max_encode(eep_dp,(unsigned int) coeff_count-1,ctx->coeff_count_context[curr_ctx_idx]+4,1,max_coeff[curr_ctx_idx]-1);
      }
      break;
    }
    case 7:
    case 8://max 15 coeffs, chroma ac different ctx uv for first bin
    {
      act_ctx = (se->k < 4)?0:1;
      if (coeff_count == 0)
        biari_encode_symbol(eep_dp, 0, ctx->coeff_count_context[curr_ctx_idx] + act_ctx);
      else
      {
        biari_encode_symbol(eep_dp, 1, ctx->coeff_count_context[curr_ctx_idx] + act_ctx);
        unary_bin_max_encode(eep_dp,(unsigned int) coeff_count-1,ctx->coeff_count_context[curr_ctx_idx]+4,1,max_coeff[curr_ctx_idx]-1);
      }
      break;
    }
    case 3:  //max 16 coeffs // 16x16 luma dc
    {
      act_ctx = 0;
      if (coeff_count == 0)
        biari_encode_symbol(eep_dp, 0, ctx->coeff_count_context[curr_ctx_idx] + act_ctx);
      else
      {
        biari_encode_symbol(eep_dp, 1, ctx->coeff_count_context[curr_ctx_idx] + act_ctx);
        unary_bin_max_encode(eep_dp,(unsigned int) coeff_count-1,ctx->coeff_count_context[curr_ctx_idx]+4,1,max_coeff[curr_ctx_idx]-1);
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
      if ((curr_ctx_idx==9)||(curr_ctx_idx==12))
        coeff_count--; // coeff_count can't be zero for 8x8 blocks - handled by CBP

      unary_bin_max_encodeABT(eep_dp,(unsigned int) coeff_count,
                              ctx->ABT_coeff_count_context[curr_abt_ctx_idx]+act_ctx, 4-act_ctx,
                              max_ccnt_abt[curr_abt_ctx_idx]);

      if ((curr_ctx_idx==9)||(curr_ctx_idx==12))
        coeff_count++; // save true coeff_count in currMB struct.

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
    prevLevel = 0;
    for(i=0; i<coeff_count;i++)
    {
      if (curr_ctx_idx < 9)
      {
        //determine run ctx
        switch (curr_ctx_idx)
        {
        case 1:
        case 2:
            curr_run_ctx = ((coeff_count) >= 4) ? 1 : 0;
            break;
        case 3:
            curr_run_ctx = ((coeff_count-i) >= 4) ? 1 : 0;
            break;
        case 4:
        case 7:
        case 8:
        case 0:
            curr_run_ctx = ((coeff_count-i) >= 3) ? 1 : 0;
            break;
        case 5:
        case 6:
            curr_run_ctx = ((coeff_count-i) >= 2) ? 1 : 0;
            break;
        }
        curr_run_ctx = 9*curr_run_ctx + curr_ctx_idx;
        //send run
        if (i==0)
          max_run = max_coeff[curr_ctx_idx]-coeff_count;
        if ((max_run != 0) )
          unary_bin_max_encode(eep_dp,(unsigned int) coeff_can[1][i],ctx->run_context[curr_run_ctx],1,max_run);
        max_run -= coeff_can[1][i];
      }
      else // ABT Run Coding
      {
        curr_abt_ctx_idx = curr_ctx_idx - NUM_TRANS_TYPE;
        curr_run_ctx = ((coeff_count) >= 4) ? 1 : 0;
        curr_run_ctx = NUM_TRANS_TYPE_ABT*curr_run_ctx + curr_abt_ctx_idx;        // 6->NUM_TRANS_TYPE_ABT to prevent selection of wrong context.
        //send run
        if (i==0)
          max_run = max_coeff_abt[curr_abt_ctx_idx]-coeff_count;
        assert(max_run>=0);

        if ((max_run != 0) )
          unary_bin_max_encodeABT(eep_dp,(unsigned int) coeff_can[1][i],ctx->ABT_run_context[curr_run_ctx],1,max_run);
        max_run -= coeff_can[1][i];
      }
      //determine level ctx
      if (curr_ctx_idx < 9)
        changed_ctx_idx = curr_ctx_idx;
      else
      {                        // use JM contexts for ABT level encoding:
        if (curr_ctx_idx < 12) // inter
          changed_ctx_idx = 1;
        else                   // intra
          changed_ctx_idx = 2;
      }
      absLevel    = absm(coeff_can[0][i]);
      changed_ctx_idx += 9*prevLevel;
      prevLevel     = absLevel > 3 ? 3 : absLevel;

      // send level
      unary_level_encode(eep_dp,(unsigned int) absm(coeff_can[0][i])-1,ctx->level_context[changed_ctx_idx]);
      sign_of_level = ((coeff_can[0][i] < 0) ? 1 : 0);
      curr_level_ctx = 3;
      biari_encode_symbol(eep_dp, (unsigned char) sign_of_level, ctx->level_context[curr_ctx_idx]+ curr_level_ctx );
    }
    coeff_count=0;
  }
}

/*!
 ************************************************************************
 * \brief
 *    Unary binarization and encoding of a symbol by using
 *    one or two distinct models for the first two and all
 *    remaining bins
*
************************************************************************/
void unary_bin_encode(EncodingEnvironmentPtr eep_dp,
                      unsigned int symbol,
                      BiContextTypePtr ctx,
                      int ctx_offset)
{
  unsigned int l;
  BiContextTypePtr ictx;

  if (symbol==0)
  {
    biari_encode_symbol(eep_dp, 0, ctx );
    return;
  }
  else
  {
    biari_encode_symbol(eep_dp, 1, ctx );
    l=symbol;
    ictx=ctx+ctx_offset;
    while ((--l)>0)
      biari_encode_symbol(eep_dp, 1, ictx);
    biari_encode_symbol(eep_dp, 0, ictx);
  }
  return;
}

/*!
 ************************************************************************
 * \brief
 *    Unary binarization and encoding of a symbol by using
 *    one or two distinct models for the first two and all
 *    remaining bins; no terminating "0" for max_symbol
 *    (finite symbol alphabet)
 ************************************************************************
 */
void unary_bin_max_encode(EncodingEnvironmentPtr eep_dp,
                          unsigned int symbol,
                          BiContextTypePtr ctx,
                          int ctx_offset,
                          unsigned int max_symbol)
{
  unsigned int l;
  BiContextTypePtr ictx;

  if (symbol==0)
  {
    biari_encode_symbol(eep_dp, 0, ctx );
    return;
  }
  else
  {
    biari_encode_symbol(eep_dp, 1, ctx );
    l=symbol;
    ictx=ctx+ctx_offset;
    while ((--l)>0)
      biari_encode_symbol(eep_dp, 1, ictx);
    if (symbol<max_symbol)
        biari_encode_symbol(eep_dp, 0, ictx);
  }
  return;
}

/*!
 ************************************************************************
 * \brief
 *    Unary binarization and encoding of a symbol by using
 *    three distinct models for the first, the second and all
 *    remaining bins
 ************************************************************************
 */
void unary_level_encode(EncodingEnvironmentPtr eep_dp,
                        unsigned int symbol,
                        BiContextTypePtr ctx)
{
  unsigned int l;
  int bin=1;
  BiContextTypePtr ictx=ctx;

  if (symbol==0)
  {
    biari_encode_symbol(eep_dp, 0, ictx );
    return;
  }
  else
  {
    biari_encode_symbol(eep_dp, 1, ictx );
    l=symbol;
    ictx++;
    while ((--l)>0)
    {
      biari_encode_symbol(eep_dp, 1, ictx  );
      if ((++bin)==2) ictx++;
    }
    biari_encode_symbol(eep_dp, 0, ictx  );
  }
  return;
}

/*!
 ************************************************************************
 * \brief
 *    Unary binarization and encoding of a symbol by using
 *    four distinct models for the first, the second, intermediate
 *    and all remaining bins
 ************************************************************************
 */
void unary_mv_encode(EncodingEnvironmentPtr eep_dp,
                        unsigned int symbol,
                        BiContextTypePtr ctx,
                        unsigned int max_bin)
{
  unsigned int l;
  unsigned int bin=1;
  BiContextTypePtr ictx=ctx;

  if (symbol==0)
  {
    biari_encode_symbol(eep_dp, 0, ictx );
    return;
  }
  else
  {
    biari_encode_symbol(eep_dp, 1, ictx );
    l=symbol;
    ictx++;
    while ((--l)>0)
    {
      biari_encode_symbol(eep_dp, 1, ictx  );
      if ((++bin)==2) ictx++;
      if (bin==max_bin) ictx++;
    }
    biari_encode_symbol(eep_dp, 0, ictx  );
  }
  return;
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
void unary_bin_max_encodeABT(EncodingEnvironmentPtr eep_dp,
                             unsigned int symbol,
                             BiContextTypePtr ctx,
                             int ctx_offset,
                             unsigned int max_symbol)
{
  unsigned int l, ibits=0;
  int bin=1;
  BiContextTypePtr ictx=ctx;

  assert (symbol<=max_symbol);

  if (symbol==0)
  {
    biari_encode_symbol(eep_dp, 0, ictx );
    return;
  }
  else
  {
    if (max_symbol > 15)
    {
      l=max_symbol-15;
      do { ibits++; }
      while ( (l>>=1) > 0);
    }

    biari_encode_symbol(eep_dp, 1, ictx );
    l=symbol;
    ictx = ctx + ctx_offset;
    while ((--l)>0)
    {
      if ((++bin)==17)
        break;
      biari_encode_symbol(eep_dp, 1, ictx  );
      if (bin==2) { ictx++; }
      if (bin==4) { ictx++; }
    }
    if (bin<16)
    {
      biari_encode_symbol(eep_dp, 0, ictx  );
    }
    if (symbol > 15)
    {
      l = symbol-15;
      ictx = ctx + ctx_offset + 3;
      do
      {
        if (l&1)
          biari_encode_symbol(eep_dp, 1, ictx  );
        else
          biari_encode_symbol(eep_dp, 0, ictx  );
        l>>=1;
      }
      while ((--ibits)>0);
    }
  }

  return;
}
