
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
void update_rc(Macroblock *currMB, short best_mode)
{
  generic_RC->MADofMB[img->current_mb_nr] = calc_MAD();

  if(input->basicunit < img->FrameSizeInMbs)
  {
    generic_RC->TotalMADBasicUnit += generic_RC->MADofMB[img->current_mb_nr];
  }  
}

/*!
 *************************************************************************************
 * \brief
 *    Update QP Parameters (critical in case of SKIP MBs or MBAFF)
 *************************************************************************************
 */

void handle_qp(Macroblock *currMB, short best_mode)
{
  // delta_qp is present only for non-skipped macroblocks
  if ((currMB->cbp!=0 || best_mode==I16MB) && (best_mode!=IPCM))
    currMB->prev_cbp = 1;
  else
  {
    currMB->delta_qp = 0;
    currMB->qp = currMB->prev_qp;
    img->qp = currMB->qp;
    currMB->prev_cbp = 0;
  }

  if (input->MbInterlace)
  {
    // update rdopt buffered qps...
    rdopt->delta_qp = currMB->delta_qp;
    rdopt->qp = currMB->qp;
    rdopt->prev_cbp = currMB->prev_cbp;

    delta_qp_mbaff[currMB->mb_field][img->bot_MB] = currMB->delta_qp;
    qp_mbaff      [currMB->mb_field][img->bot_MB] = currMB->qp;
  }
  set_chroma_qp(currMB);
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
int calc_MAD()
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
void copy_rc_generic( rc_generic *dst, rc_generic *src )
{
  /* buffer original addresses for which memory has been allocated */
  int *tmpMADofMB = dst->MADofMB;

  /* copy object */
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
void generic_alloc( rc_generic **prc )
{
  *prc = (rc_generic *) malloc ( sizeof( rc_generic ) );
  if (NULL == *prc)
  {
    no_mem_exit("init_global_buffers: generic_alloc");
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
void generic_free(rc_generic **prc)
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
void init_GOP_rc()
{
  int M,N,n,np,nb;           //Rate control
  {
    if (img->type == I_SLICE && ((input->RCUpdateMode != RC_MODE_1 && input->RCUpdateMode != RC_MODE_3) || !(img->number) ) )
    {
      if (input->intra_period == 0)
      {
        n = input->no_frames + (input->no_frames - 1) * input->successive_Bframe;

        /* number of P frames */
        np = input->no_frames-1;

        /* number of B frames */
        nb = (input->no_frames - 1) * input->successive_Bframe;
      }
      else if ( input->RCUpdateMode != RC_MODE_1 && input->RCUpdateMode != RC_MODE_3 )
      {
        N = input->intra_period*(input->successive_Bframe+1);
        M = input->successive_Bframe+1;
        n = (img->number == 0) ? N - ( M - 1) : N;

        /* last GOP may contain less frames */
        if ((img->number / input->intra_period) >= (input->no_frames / input->intra_period))
        {
          if (img->number != 0)
            n = (input->no_frames - img->number) + (input->no_frames - img->number - 1) * input->successive_Bframe + input->successive_Bframe;
          else
            n = input->no_frames  + (input->no_frames - 1) * input->successive_Bframe;
        }

        /* number of P frames */
        if (img->number == 0)
          np = (n + 2 * (M - 1)) / M - 1; /* first GOP */
        else
          np = (n + (M - 1)) / M - 1;

        /* number of B frames */
        nb = n - np - 1;
      }
      else // applies RC to I and B slices
      {
        np = input->no_frames - 1; // includes I and P slices/frames except the very first IDR I_SLICE
        nb = (input->no_frames - 1) * input->successive_Bframe;
      }
      rc_init_GOP(quadratic_RC, np, nb);
    }
  }
}

/*!
 *************************************************************************************
 * \brief
 *    Initialize Frame Level Rate Control parameters
 *
 *************************************************************************************
 */
void init_frame_rc(int FrameNumberInFile)
{
  /*update the number of MBs in the basic unit for MB adaptive
  f/f coding*/
  if( (input->MbInterlace) && (input->basicunit < img->FrameSizeInMbs) && (img->type == P_SLICE || (input->RCUpdateMode == RC_MODE_1 && img->number) ) )
    img->BasicUnit = input->basicunit << 1;
  else
    img->BasicUnit = input->basicunit;

  if ( input->RDPictureDecision )
  {
    // store rate allocation quadratic...
    copy_rc_jvt( quadratic_RC_init, quadratic_RC );
    // ...and generic model
    copy_rc_generic( generic_RC_init, generic_RC );
  }

  rc_init_pict(quadratic_RC, 1,0,1, 1.0F);   

  img->qp  = updateQP(quadratic_RC, 0);

  //pic_type = img->type;
  //QP =0;

  if( active_sps->frame_mbs_only_flag)
    generic_RC->TopFieldFlag=0;
}
