
/*!
 ***************************************************************************
 * \file ratectl.c
 *
 * \brief
 *    Rate Control algorithm
 *
 * \author
 *    Main contributors (see contributors.h for copyright, address and affiliation details)
 *     - Siwei Ma <swma@jdl.ac.cn>
 *     - Zhengguo LI<ezgli@lit.a-star.edu.sg>
 *
 * \date
 *   16 Jan. 2003
 **************************************************************************
 */

#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <memory.h>
#include <limits.h>

#include "global.h"
#include "ratectl.h"


/*!
 *************************************************************************************
 * \brief
 *    Update Rate Control Parameters
 *************************************************************************************
 */
void rc_store_mad(Macroblock *currMB)
{
  generic_RC->MADofMB[img->current_mb_nr] = ComputeMBMAD();

  if(input->basicunit < img->FrameSizeInMbs)
  {
    generic_RC->TotalMADBasicUnit += generic_RC->MADofMB[img->current_mb_nr];
  }  
}

/*!
 *************************************************************************************
 * \brief
 *    Update QP Parameters (in case of SKIP MBs or MBAFF)
 *************************************************************************************
 */

void update_qp_cbp(Macroblock *currMB, short best_mode)
{
  // delta_qp is present only for non-skipped macroblocks
  if ((currMB->cbp!=0 || best_mode == I16MB) && (best_mode != IPCM))
    currMB->prev_cbp = 1;
  else
  {
    currMB->delta_qp     = 0;
    currMB->qp           = currMB->prev_qp;
    
    set_chroma_qp(currMB);

    currMB->qp_scaled[0] = currMB->qp + img->bitdepth_luma_qp_scale - MIN_QP;
    currMB->qp_scaled[1] = currMB->qpc[0] + img->bitdepth_chroma_qp_scale;
    currMB->qp_scaled[2] = currMB->qpc[1] + img->bitdepth_chroma_qp_scale;

    select_dct(currMB);

    img->qp           = currMB->qp;
    currMB->prev_cbp  = 0;
  }

  if (input->MbInterlace)
  {
    // update rdopt buffered qps...
    rdopt->qp        = currMB->qp;
    rdopt->delta_qp  = currMB->delta_qp;
    rdopt->prev_cbp  = currMB->prev_cbp;
    memcpy(rdopt->qp_scaled, currMB->qp_scaled, 3 * sizeof(int));

    delta_qp_mbaff[currMB->mb_field][img->bot_MB] = currMB->delta_qp;
    qp_mbaff      [currMB->mb_field][img->bot_MB] = currMB->qp;
  }  
}

/*!
 *************************************************************************************
 * \brief
 *    map QP to Qstep
 *
 *************************************************************************************
*/
double QP2Qstep( int QP )
{
  int i;
  double Qstep;
  static const double QP2QSTEP[6] = { 0.625, 0.6875, 0.8125, 0.875, 1.0, 1.125 };

  Qstep = QP2QSTEP[QP % 6];
  for( i=0; i<(QP/6); i++)
    Qstep *= 2;

  return Qstep;
}


/*!
 *************************************************************************************
 * \brief
 *    map Qstep to QP
 *
 *************************************************************************************
*/
int Qstep2QP( double Qstep )
{
  int q_per = 0, q_rem = 0;

  //  assert( Qstep >= QP2Qstep(0) && Qstep <= QP2Qstep(51) );
  if( Qstep < QP2Qstep(0))
    return 0;
  else if (Qstep > QP2Qstep(51) )
    return 51;

  while( Qstep > QP2Qstep(5) )
  {
    Qstep /= 2.0;
    q_per += 1;
  }

  if (Qstep <= 0.65625)
  {
    Qstep = 0.625;
    q_rem = 0;
  }
  else if (Qstep <= 0.75)
  {
    Qstep = 0.6875;
    q_rem = 1;
  }
  else if (Qstep <= 0.84375)
  {
    Qstep = 0.8125;
    q_rem = 2;
  }
  else if (Qstep <= 0.9375)
  {
    Qstep = 0.875;
    q_rem = 3;
  }
  else if (Qstep <= 1.0625)
  {
    Qstep = 1.0;
    q_rem = 4;
  }
  else
  {
    Qstep = 1.125;
    q_rem = 5;
  }

  return (q_per * 6 + q_rem);
}

/*!
 ************************************************************************************
 * \brief
 *    calculate MAD for the current macroblock
 *
 * \return
 *    calculated MAD
 *
 *************************************************************************************
*/
int ComputeMBMAD()
{
  int k, l, sum = 0;

  for (k = 0; k < 16; k++)
    for (l = 0; l < 16; l++)
      sum += iabs(diffy[k][l]);

  return sum;
}

/*!
 *************************************************************************************
 * \brief
 *    Compute Frame MAD
 *
 *************************************************************************************
*/
double ComputeFrameMAD()
{
  int64 TotalMAD = 0;
  unsigned int i;
  for(i = 0; i < img->FrameSizeInMbs; i++)
    TotalMAD += generic_RC->MADofMB[i];
  return (double)TotalMAD / (256.0 * (double)img->FrameSizeInMbs);
}


/*!
 *************************************************************************************
 * \brief
 *    Copy JVT rate control objects
 *
 *************************************************************************************
*/
void rc_copy_generic( rc_generic *dst, rc_generic *src )
{
  /* buffer original addresses for which memory has been allocated */
  int *tmpMADofMB = dst->MADofMB;

  /* copy object */

  // This could be written as: *dst = *src;
  memcpy( (void *)dst, (void *)src, sizeof(rc_generic) );

  /* restore original addresses */
  dst->MADofMB = tmpMADofMB;

  /* copy MADs */
  memcpy( (void *)dst->MADofMB, (void *)src->MADofMB, img->FrameSizeInMbs * sizeof (int) );
}

/*!
 *************************************************************************************
 * \brief
 *    Dynamically allocate memory needed for generic rate control
 *
 *************************************************************************************
 */
void rc_alloc_generic( rc_generic **prc )
{
  *prc = (rc_generic *) malloc ( sizeof( rc_generic ) );
  if (NULL == *prc)
  {
    no_mem_exit("init_global_buffers: rc_alloc_generic");
  }
  (*prc)->MADofMB = (int *) calloc (img->FrameSizeInMbs, sizeof (int));
  if (NULL == (*prc)->MADofMB)
  {
    no_mem_exit("init_global_buffers: (*prc)->MADofMB");
  }
  (*prc)->FieldFrame = 1;
}


/*!
 *************************************************************************************
 * \brief
 *    Free memory needed for generic rate control
 *
 *************************************************************************************
 */
void rc_free_generic(rc_generic **prc)
{
  if (NULL!=(*prc)->MADofMB)
  {
    free ((*prc)->MADofMB);
    (*prc)->MADofMB = NULL;
  }
  if (NULL!=(*prc))
  {
    free ((*prc));
    (*prc) = NULL;
  }
}

/*!
 *************************************************************************************
 * \brief
 *    Initialize GOP Level Rate Control parameters
 *
 *************************************************************************************
 */
void rc_init_gop_params(void)
{
  int np, nb; 

  switch( input->RCUpdateMode )
  {
  case RC_MODE_1: case RC_MODE_3: 
    if (!(img->number) )
    {
      /* number of P frames */
      np = input->no_frames - 1;
      /* number of B frames */
      nb = np * input->successive_Bframe;

      rc_init_GOP(quadratic_RC, np, nb);
    }
    break;
  case RC_MODE_0: case RC_MODE_2:
    if (input->idr_period == 0)
    {
      /* number of P frames */
      np = input->no_frames - 1;
      /* number of B frames */
      nb = np * input->successive_Bframe;
    }
    else   
    {
      int M = input->successive_Bframe + 1;
      int N = M * input->idr_period;      
      int n = (img->number == 0) ? N - ( M - 1) : N;

      /* last GOP may contain less frames */
      if ((img->number / input->idr_period) >= (input->no_frames / input->idr_period))
      {
        if (img->number != 0)
          n = (input->no_frames - img->number) * (input->successive_Bframe + 1);
        else
          n = input->no_frames  + (input->no_frames - 1) * input->successive_Bframe;
      }

      /* number of P frames */
      np = (img->number == 0) ? 1 + ((n - 2) / M) : (n - 1) / M; 
      /* number of B frames */
      nb = n - np - 1;
    }
    rc_init_GOP(quadratic_RC, np, nb);
    break;
  default:
    break;
  }
}

/*!
 *************************************************************************************
 * \brief
 *    Initialize Frame Level Rate Control parameters
 *
 *************************************************************************************
 */

void rc_init_frame(int FrameNumberInFile)
{
  switch( input->RCUpdateMode )
  {
  case RC_MODE_0:  case RC_MODE_1:  case RC_MODE_2:  case RC_MODE_3:

  // update the number of MBs in the basic unit for MBAFF coding
  if( (input->MbInterlace) && (input->basicunit < img->FrameSizeInMbs) && (img->type == P_SLICE || (input->RCUpdateMode == RC_MODE_1 && img->number) ) )
    img->BasicUnit = input->basicunit << 1;
  else
    img->BasicUnit = input->basicunit;

    if ( input->RDPictureDecision )
    {    
      rc_copy_quadratic( quadratic_RC_init, quadratic_RC ); // store rate allocation quadratic...    
      rc_copy_generic( generic_RC_init, generic_RC ); // ...and generic model
    }
    rc_init_pict_ptr(quadratic_RC, 1,0,1, 1.0F);

    if( active_sps->frame_mbs_only_flag)
      generic_RC->TopFieldFlag=0;

    img->qp = updateQP(quadratic_RC, 0);
    break;
  default:
    break;
  }
}

/*!
 *************************************************************************************
 * \brief
 *    Initialize Sequence Level Rate Control parameters
 *
 *************************************************************************************
 */

void rc_init_sequence(void)
{
  switch( input->RCUpdateMode )
  {
  case RC_MODE_0:  case RC_MODE_1:  case RC_MODE_2:  case RC_MODE_3:
    rc_init_seq(quadratic_RC);
    break;
  default:
    break;
  }
}

void rc_store_slice_header_bits( int len )
{
  switch (input->RCUpdateMode)
  {
  case RC_MODE_0:  case RC_MODE_1:  case RC_MODE_2:  case RC_MODE_3:
    generic_RC->NumberofHeaderBits +=len;

    // basic unit layer rate control
    if(img->BasicUnit < img->FrameSizeInMbs)
      generic_RC->NumberofBasicUnitHeaderBits +=len;
    break;
  default:
    break;
  }
}

void update_qp_cbp_tmp(Macroblock *currMB, int cbp, int best_mode)
{
  if (((cbp!=0 || best_mode==I16MB) && (best_mode!=IPCM) ))
    currMB->prev_cbp = 1;
  else if ((cbp==0 && !input->RCEnable) || (best_mode==IPCM))
  {
    currMB->delta_qp  = 0;
    currMB->qp        = currMB->prev_qp;
    currMB->qp_scaled[0] = currMB->qp + img->bitdepth_luma_qp_scale - MIN_QP;
    currMB->qp_scaled[1] = currMB->qpc[0] + img->bitdepth_chroma_qp_scale;
    currMB->qp_scaled[2] = currMB->qpc[1] + img->bitdepth_chroma_qp_scale;

    select_dct(currMB);

    set_chroma_qp(currMB);
    img->qp           = currMB->qp;
    currMB->prev_cbp  = 0;
  }
}

/*!
*************************************************************************************
* \brief
*    Update Rate Control Difference
*************************************************************************************
*/
void rc_store_diff(int cpix_x, int cpix_y, imgpel prediction[16][16])
{
  int i, j;
  int *iDst;
  imgpel *Src1, *Src2;

  for(j = 0; j < MB_BLOCK_SIZE; j++)
  {
    iDst = diffy[j];
    Src1 = &pCurImg[cpix_y + j][cpix_x];
    Src2 = prediction[j];
    for (i = 0; i < MB_BLOCK_SIZE; i++)
    {
      iDst[i] = Src1[i] - Src2[i];
    }
  }
}

