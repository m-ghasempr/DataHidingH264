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
 * \file context_ini.c
 *
 * \brief
 *    CABAC context initializations
 *
 * \author
 *    Main contributors (see contributors.h for copyright, address and affiliation details)
 *    - Detlev Marpe                    <marpe@hhi.de>
 *    - Heiko Schwarz                   <hschwarz@hhi.de>
 **************************************************************************************
 */

#define CONTEXT_INI_C

#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <string.h>

#include "defines.h"
#include "global.h"
#include "biaridecod.h"

#include "ctx_tables.h"




#define BIARI_CTX_INIT2(ii,jj,ctx,tab,num) \
{ \
  for (i=0; i<ii; i++) \
  for (j=0; j<jj; j++) \
  { \
    if      (img->type==I_SLICE)  biari_init_context (img, &(ctx[i][j]), &(tab ## _I[num][i][j][0])); \
    else                            biari_init_context (img, &(ctx[i][j]), &(tab ## _P[num][i][j][0])); \
  } \
}
#define BIARI_CTX_INIT1(jj,ctx,tab,num) \
{ \
  for (j=0; j<jj; j++) \
  { \
  if      (img->type==I_SLICE)  biari_init_context (img, &(ctx[j]), &(tab ## _I[num][0][j][0])); \
    else                            biari_init_context (img, &(ctx[j]), &(tab ## _P[num][0][j][0])); \
  } \
}


void
init_contexts (struct img_par* img)
{
  MotionInfoContexts*  mc = img->currentSlice->mot_ctx;
  TextureInfoContexts* tc = img->currentSlice->tex_ctx;
  int i, j;

  //printf("%d -", img->model_number);

  //--- motion coding contexts ---
  BIARI_CTX_INIT2 (3, NUM_MB_TYPE_CTX,   mc->mb_type_contexts,     INIT_MB_TYPE,    img->model_number);
  BIARI_CTX_INIT2 (2, NUM_B8_TYPE_CTX,   mc->b8_type_contexts,     INIT_B8_TYPE,    img->model_number);
  BIARI_CTX_INIT2 (2, NUM_MV_RES_CTX,    mc->mv_res_contexts,      INIT_MV_RES,     img->model_number);
  BIARI_CTX_INIT2 (2, NUM_REF_NO_CTX,    mc->ref_no_contexts,      INIT_REF_NO,     img->model_number);
  BIARI_CTX_INIT1 (   NUM_DELTA_QP_CTX,  mc->delta_qp_contexts,    INIT_DELTA_QP,   img->model_number);
  BIARI_CTX_INIT1 (   NUM_MB_AFF_CTX,    mc->mb_aff_contexts,      INIT_MB_AFF,     img->model_number);

  //--- texture coding contexts ---
  BIARI_CTX_INIT1 (                 NUM_IPR_CTX,  tc->ipr_contexts,     INIT_IPR,       img->model_number);
  BIARI_CTX_INIT1 (                 NUM_CIPR_CTX, tc->cipr_contexts,    INIT_CIPR,      img->model_number);
  BIARI_CTX_INIT2 (3,               NUM_CBP_CTX,  tc->cbp_contexts,     INIT_CBP,       img->model_number);
  BIARI_CTX_INIT2 (NUM_BLOCK_TYPES, NUM_BCBP_CTX, tc->bcbp_contexts,    INIT_BCBP,      img->model_number);
  BIARI_CTX_INIT2 (NUM_BLOCK_TYPES, NUM_MAP_CTX,  tc->map_contexts,     INIT_MAP,       img->model_number);
  BIARI_CTX_INIT2 (NUM_BLOCK_TYPES, NUM_LAST_CTX, tc->last_contexts,    INIT_LAST,      img->model_number);
  BIARI_CTX_INIT2 (NUM_BLOCK_TYPES, NUM_ONE_CTX,  tc->one_contexts,     INIT_ONE,       img->model_number);
  BIARI_CTX_INIT2 (NUM_BLOCK_TYPES, NUM_ABS_CTX,  tc->abs_contexts,     INIT_ABS,       img->model_number);
  BIARI_CTX_INIT2 (NUM_BLOCK_TYPES, NUM_MAP_CTX,  tc->fld_map_contexts, INIT_FLD_MAP,   img->model_number);
  BIARI_CTX_INIT2 (NUM_BLOCK_TYPES, NUM_LAST_CTX, tc->fld_last_contexts,INIT_FLD_LAST,  img->model_number);
}

