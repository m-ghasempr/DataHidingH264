
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

#include <math.h>
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
  ImageParameters *p_Img = currMB->p_Img;
  InputParameters *p_Inp = currMB->p_Inp;
  RCGeneric *p_gen = p_Img->p_rc_gen;

  p_gen->MADofMB[currMB->mbAddrX] = ComputeMBMAD(currMB->p_slice->diffy);

  if(p_Inp->basicunit < p_Img->FrameSizeInMbs)
  {
    p_gen->TotalMADBasicUnit += p_gen->MADofMB[currMB->mbAddrX];
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
int Qstep2QP( double Qstep, int qp_offset )
{
  int q_per = 0, q_rem = 0;

  if( Qstep < QP2Qstep(MIN_QP))
    return MIN_QP;
  else if (Qstep > QP2Qstep(MAX_QP + qp_offset) )
    return (MAX_QP + qp_offset);

  while( Qstep > QP2Qstep(5) )
  {
    Qstep /= 2.0;
    q_per++;
  }

  if (Qstep <= 0.65625)
  {
    //Qstep = 0.625;
    q_rem = 0;
  }
  else if (Qstep <= 0.75)
  {
    //Qstep = 0.6875;
    q_rem = 1;
  }
  else if (Qstep <= 0.84375)
  {
    //Qstep = 0.8125;
    q_rem = 2;
  }
  else if (Qstep <= 0.9375)
  {
    //Qstep = 0.875;
    q_rem = 3;
  }
  else if (Qstep <= 1.0625)
  {
    //Qstep = 1.0;
    q_rem = 4;
  }
  else
  {
    //Qstep = 1.125;
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
int ComputeMBMAD(int diffy[16][16])
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
double ComputeFrameMAD(ImageParameters *p_Img)
{
  int64 TotalMAD = 0;
  unsigned int i;
  for(i = 0; i < p_Img->FrameSizeInMbs; i++)
    TotalMAD += p_Img->p_rc_gen->MADofMB[i];
  return (double)TotalMAD / (256.0 * (double)p_Img->FrameSizeInMbs);
}


/*!
 *************************************************************************************
 * \brief
 *    Copy JVT rate control objects
 *
 *************************************************************************************
*/
void rc_copy_generic( ImageParameters *p_Img, RCGeneric *dst, RCGeneric *src )
{
  /* buffer original addresses for which memory has been allocated */
  int *tmpMADofMB = dst->MADofMB;

  /* copy object */

  // This could be written as: *dst = *src;
  memcpy( (void *)dst, (void *)src, sizeof(RCGeneric) );

  /* restore original addresses */
  dst->MADofMB = tmpMADofMB;

  /* copy MADs */
  memcpy( (void *)dst->MADofMB, (void *)src->MADofMB, p_Img->FrameSizeInMbs * sizeof (int) );
}

/*!
 *************************************************************************************
 * \brief
 *    Dynamically allocate memory needed for generic rate control
 *
 *************************************************************************************
 */
void rc_alloc_generic( ImageParameters *p_Img, RCGeneric **p_quad )
{
  *p_quad = (RCGeneric *) malloc ( sizeof( RCGeneric ) );
  if (NULL == *p_quad)
  {
    no_mem_exit("rc_alloc_generic: rc_alloc_generic");
  }
  (*p_quad)->MADofMB = (int *) calloc (p_Img->FrameSizeInMbs, sizeof (int));
  if (NULL == (*p_quad)->MADofMB)
  {
    no_mem_exit("rc_alloc_generic: (*p_quad)->MADofMB");
  }
  (*p_quad)->FieldFrame = 1;
}


/*!
 *************************************************************************************
 * \brief
 *    Free memory needed for generic rate control
 *
 *************************************************************************************
 */
void rc_free_generic(RCGeneric **p_quad)
{
  if (NULL!=(*p_quad)->MADofMB)
  {
    free ((*p_quad)->MADofMB);
    (*p_quad)->MADofMB = NULL;
  }
  if (NULL!=(*p_quad))
  {
    free ((*p_quad));
    (*p_quad) = NULL;
  }
}

/*!
 *************************************************************************************
 * \brief
 *    Initialize GOP Level Rate Control parameters
 *
 *************************************************************************************
 */
void rc_init_gop_params(ImageParameters *p_Img, InputParameters *p_Inp)
{
  int np, nb;
  RCQuadratic *p_quad = p_Img->p_rc_quad;
  RCGeneric *p_gen = p_Img->p_rc_gen;

  switch( p_Inp->RCUpdateMode )
  {
  case RC_MODE_1: case RC_MODE_3: 
    if ( !(p_Img->number) )
    {
      /* number of P frames */
      np = p_Inp->no_frm_base - 1;
      /* number of B frames */
      nb = np * p_Inp->NumberBFrames;

      rc_init_GOP(p_Img, p_Inp, p_quad, p_gen, np, nb);
    }
    break;
  case RC_MODE_0: case RC_MODE_2:
    if (p_Inp->idr_period == 0)
    {
      if ( !(p_Img->number) )
      {
        /* number of P frames */
        np = p_Inp->no_frm_base - 1;
        /* number of B frames */
        nb = np * p_Inp->NumberBFrames;
        rc_init_GOP(p_Img, p_Inp, p_quad, p_gen, np, nb);
      }
    }
    else if ( (!p_Inp->adaptive_idr_period && ( p_Img->frm_number - p_Img->lastIDRnumber ) % p_Inp->idr_period == 0)
      || (p_Inp->adaptive_idr_period == 1 && ( p_Img->frm_number - imax(p_Img->lastIntraNumber, p_Img->lastIDRnumber) ) % p_Inp->idr_period == 0) )  
    {
      int M = p_Inp->NumberBFrames + 1;
      int N = M * p_Inp->idr_period;      
      int n = (p_Img->number == 0) ? N - ( M - 1) : N;

      /* last GOP may contain less frames */
      if ((p_Img->number / p_Inp->idr_period) >= (p_Inp->no_frm_base / p_Inp->idr_period))
      {
        if (p_Img->number != 0)
          n = (p_Inp->no_frm_base - p_Img->number) * (p_Inp->NumberBFrames + 1);
        else
          n = p_Inp->no_frm_base  + (p_Inp->no_frm_base - 1) * p_Inp->NumberBFrames;
      }

      /* number of P frames */
      np = (p_Img->number == 0) ? 1 + ((n - 2) / M) : (n - 1) / M; 
      /* number of B frames */
      nb = n - np - 1;
      rc_init_GOP(p_Img, p_Inp, p_quad, p_gen, np, nb);
    }
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

void rc_init_frame(ImageParameters *p_Img, InputParameters *p_Inp)
{
  switch( p_Inp->RCUpdateMode )
  {
  case RC_MODE_0:  case RC_MODE_1:  case RC_MODE_2:  case RC_MODE_3:

  // update the number of MBs in the basic unit for MBAFF coding
  if( (p_Inp->MbInterlace) && (p_Inp->basicunit < p_Img->FrameSizeInMbs) && (p_Img->type == P_SLICE || (p_Inp->RCUpdateMode == RC_MODE_1 && p_Img->number) ) )
    p_Img->BasicUnit = p_Inp->basicunit << 1;
  else
    p_Img->BasicUnit = p_Inp->basicunit;

    if ( p_Inp->RDPictureDecision )
    {    
      rc_copy_quadratic( p_Img, p_Inp, p_Img->p_rc_quad_init, p_Img->p_rc_quad ); // store rate allocation quadratic...    
      rc_copy_generic( p_Img, p_Img->p_rc_gen_init, p_Img->p_rc_gen ); // ...and generic model
    }
    p_Img->rc_init_pict_ptr(p_Img, p_Inp, p_Img->p_rc_quad, p_Img->p_rc_gen, 1,0,1, 1.0F);

    if( p_Img->active_sps->frame_mbs_only_flag)
      p_Img->p_rc_gen->TopFieldFlag=0;

    p_Img->qp = p_Img->updateQP(p_Img, p_Inp, p_Img->p_rc_quad, p_Img->p_rc_gen, 0) - p_Img->p_rc_quad->bitdepth_qp_scale;
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

void rc_init_sequence(ImageParameters *p_Img, InputParameters *p_Inp)
{
  switch( p_Inp->RCUpdateMode )
  {
  case RC_MODE_0:  case RC_MODE_1:  case RC_MODE_2:  case RC_MODE_3:
    rc_init_seq(p_Img, p_Inp, p_Img->p_rc_quad, p_Img->p_rc_gen);
    break;
  default:
    break;
  }
}

void rc_store_slice_header_bits( ImageParameters *p_Img, InputParameters *p_Inp, int len )
{
  switch (p_Inp->RCUpdateMode)
  {
  case RC_MODE_0:  case RC_MODE_1:  case RC_MODE_2:  case RC_MODE_3:
    p_Img->p_rc_gen->NumberofHeaderBits +=len;

    // basic unit layer rate control
    if(p_Img->BasicUnit < p_Img->FrameSizeInMbs)
      p_Img->p_rc_gen->NumberofBasicUnitHeaderBits +=len;
    break;
  default:
    break;
  }
}


/*!
*************************************************************************************
* \brief
*    Update Rate Control Difference
*************************************************************************************
*/
void rc_store_diff(int diffy[16][16], imgpel **p_curImg, int cpix_x,imgpel **prediction)
{
  int i, j;
  int *iDst;
  imgpel *Src1, *Src2;  

  for(j = 0; j < MB_BLOCK_SIZE; j++)
  {
    iDst = diffy[j];
    Src1 = &p_curImg[j][cpix_x];
    Src2 = prediction[j];
    for (i = 0; i < MB_BLOCK_SIZE; i++)
    {
      iDst[i] = Src1[i] - Src2[i];
    }
  }
}

