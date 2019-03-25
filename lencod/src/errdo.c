
/*!
 *************************************************************************************
 * \file errdo.c
 *
 * \brief
 *    Contains functions that estimate the distortion for the
 *    rate-distortion optimization with losses.
 * \date
 *    October 22nd, 2001
 *
 * \author
 *    Main contributors (see contributors.h for copyright, address and
 *    affiliation details)
 *    - Dimitrios Kontopodis                    <dkonto@eikon.tum.de>
 *    Code revamped July 2008 by:
 *    - Peshala Pahalawatta (ppaha@dolby.com)
 *    - Alexis Tourapis (atour@dolby.com)
 *	  Code modularized to support more distortion estimation algorithms in June 2009 by 
 *	  - Zhifeng Chen (zzchen@dolby.com)
 *************************************************************************************
 */

#include "global.h"
#include "memalloc.h"
#include "refbuf.h"
#include "image.h"
#include "md_common.h"
#include "errdo.h"
#include "errdo_dist_mhyp.h"

/*!
**************************************************************************************
* \brief 
*      Allocate memory for error resilient RDO.  
**************************************************************************************
*/
int allocate_errdo_mem(VideoParameters *p_Vid, InputParameters *p_Inp)
{
  int memory_size = 0;

  //allocate shared memory for all algorithms
  p_Vid->p_decs   = (Decoders *) malloc(sizeof(Decoders));
  memory_size += get_mem3Dint(&p_Vid->p_decs->res_img, MAX_PLANE, MB_BLOCK_SIZE, MB_BLOCK_SIZE);
  memory_size += get_mem3Dint(&p_Vid->p_decs->res_mb_best8x8, MAX_PLANE, MB_BLOCK_SIZE, MB_BLOCK_SIZE);

  p_Vid->p_decs->RCD_bestY_mb         = NULL;
  p_Vid->p_decs->RCD_bestY_b8x8         = NULL;
  p_Vid->p_decs->MVCD_bestY_mb         = NULL;
  p_Vid->p_decs->MVCD_bestY_b8x8         = NULL;
  p_Vid->p_decs->flag_bestY_mb         = NULL;
  p_Vid->p_decs->flag_bestY_b8x8         = NULL;
  p_Vid->p_decs->flag_wo_res         = NULL;
  p_Vid->p_decs->flag_wo_res_bestY_b8x8         = NULL;
  p_Vid->p_decs->trans_dist_bestY_mb         = NULL;
  p_Vid->p_decs->trans_dist_bestY_b8x8         = NULL;
  p_Vid->p_decs->trans_dist_wo_res         = NULL;   //it is used for P8x8, where residual may be set to 0
  p_Vid->p_decs->trans_dist_wo_res_bestY_b8x8   = NULL;   //it is used for P8x8, where residual may be set to 0
  p_Vid->p_decs->trans_err_bestY_mb         = NULL;
  p_Vid->p_decs->trans_err_bestY_b8x8         = NULL;
  p_Vid->p_decs->trans_err_wo_res         = NULL;   //it is used for P8x8, where residual may be set to 0
  p_Vid->p_decs->trans_err_wo_res_bestY_b8x8   = NULL;   //it is used for P8x8, where residual may be set to 0
  p_Vid->p_decs->dec_mb_pred         = NULL;
  p_Vid->p_decs->dec_mbY_best        = NULL;
  p_Vid->p_decs->dec_mb_pred_best8x8 = NULL;
  p_Vid->p_decs->dec_mbY_best8x8     = NULL;
  p_Vid->p_decs->first_moment_bestY_mb         = NULL;
  p_Vid->p_decs->first_moment_bestY_b8x8       = NULL;
  p_Vid->p_decs->first_moment_pred_bestY_b8x8       = NULL;
  p_Vid->p_decs->first_moment_pred       = NULL;
  p_Vid->p_decs->second_moment_bestY_mb        = NULL;
  p_Vid->p_decs->second_moment_bestY_b8x8      = NULL;
  p_Vid->p_decs->second_moment_pred_bestY_b8x8      = NULL;
  p_Vid->p_decs->second_moment_pred      = NULL;

  //Zhifeng 090630
  switch (p_Inp->de)
  {
  case RMPC:
    //allocate memory for rmpc algorithm
    memory_size += get_mem2Dint(&p_Vid->p_decs->trans_dist_bestY_mb, MB_BLOCK_SIZE, MB_BLOCK_SIZE);
    memory_size += get_mem3Dint(&p_Vid->p_decs->trans_dist_bestY_b8x8, 2, MB_BLOCK_SIZE, MB_BLOCK_SIZE);
    memory_size += get_mem2Dint(&p_Vid->p_decs->trans_dist_wo_res, MB_BLOCK_SIZE, MB_BLOCK_SIZE);
    memory_size += get_mem3Dint(&p_Vid->p_decs->trans_dist_wo_res_bestY_b8x8, 2, MB_BLOCK_SIZE, MB_BLOCK_SIZE);
#ifdef RMPC_PAPER
    memory_size += get_mem2Dint(&p_Vid->p_decs->RCD_bestY_mb, MB_BLOCK_SIZE, MB_BLOCK_SIZE);
    memory_size += get_mem3Dint(&p_Vid->p_decs->RCD_bestY_b8x8, 2, MB_BLOCK_SIZE, MB_BLOCK_SIZE);
    memory_size += get_mem2Dint(&p_Vid->p_decs->MVCD_bestY_mb, MB_BLOCK_SIZE, MB_BLOCK_SIZE);
    memory_size += get_mem3Dint(&p_Vid->p_decs->MVCD_bestY_b8x8, 2, MB_BLOCK_SIZE, MB_BLOCK_SIZE);
#endif  //#ifdef EMPC_PAPER
#ifdef RMPC_ERRDO
    memory_size += get_mem2D(&p_Vid->p_decs->flag_bestY_mb, MB_BLOCK_SIZE, MB_BLOCK_SIZE);
    memory_size += get_mem3D(&p_Vid->p_decs->flag_bestY_b8x8, 2, MB_BLOCK_SIZE, MB_BLOCK_SIZE);
    memory_size += get_mem2D(&p_Vid->p_decs->flag_wo_res, MB_BLOCK_SIZE, MB_BLOCK_SIZE);
    memory_size += get_mem3D(&p_Vid->p_decs->flag_wo_res_bestY_b8x8, 2, MB_BLOCK_SIZE, MB_BLOCK_SIZE);
#endif  //#ifdef RMPC_ERRDO
    break;
  case ERMPC_NN:
  case ERMPC_FAST:
  case EXTENDED_RMPC:
    //allocate memory for extended rmpc algorithm
    memory_size += get_mem2Dint(&p_Vid->p_decs->trans_dist_bestY_mb, MB_BLOCK_SIZE, MB_BLOCK_SIZE);
    memory_size += get_mem3Dint(&p_Vid->p_decs->trans_dist_bestY_b8x8, 2, MB_BLOCK_SIZE, MB_BLOCK_SIZE);
    memory_size += get_mem2Dint(&p_Vid->p_decs->trans_dist_wo_res, MB_BLOCK_SIZE, MB_BLOCK_SIZE);
    memory_size += get_mem3Dint(&p_Vid->p_decs->trans_dist_wo_res_bestY_b8x8, 2, MB_BLOCK_SIZE, MB_BLOCK_SIZE);
    memory_size += get_mem2Dint(&p_Vid->p_decs->trans_err_bestY_mb, MB_BLOCK_SIZE, MB_BLOCK_SIZE);
    memory_size += get_mem3Dint(&p_Vid->p_decs->trans_err_bestY_b8x8, 2, MB_BLOCK_SIZE, MB_BLOCK_SIZE);
    memory_size += get_mem2Dint(&p_Vid->p_decs->trans_err_wo_res, MB_BLOCK_SIZE, MB_BLOCK_SIZE);
    memory_size += get_mem3Dint(&p_Vid->p_decs->trans_err_wo_res_bestY_b8x8, 2, MB_BLOCK_SIZE, MB_BLOCK_SIZE);
    break;
  case LLN:
  case FAST_LLN:
    //allocate memory for lln algorithm
  memory_size += get_mem3Dpel(&p_Vid->p_decs->dec_mb_pred, p_Inp->NoOfDecoders, MB_BLOCK_SIZE, MB_BLOCK_SIZE);
  memory_size += get_mem3Dpel(&p_Vid->p_decs->dec_mbY_best, p_Inp->NoOfDecoders, MB_BLOCK_SIZE, MB_BLOCK_SIZE);
  memory_size += get_mem4Dpel(&p_Vid->p_decs->dec_mbY_best8x8, 2, p_Inp->NoOfDecoders, MB_BLOCK_SIZE, MB_BLOCK_SIZE);
  memory_size += get_mem4Dpel(&p_Vid->p_decs->dec_mb_pred_best8x8, 2, p_Inp->NoOfDecoders, MB_BLOCK_SIZE, MB_BLOCK_SIZE);
    break;
  case ROPE:
    //Zhifeng 090616
    //allocate memory for rope algorithm
    memory_size += get_mem2Dpel(&p_Vid->p_decs->first_moment_bestY_mb, MB_BLOCK_SIZE, MB_BLOCK_SIZE);
    memory_size += get_mem3Dpel(&p_Vid->p_decs->first_moment_bestY_b8x8, 2, MB_BLOCK_SIZE, MB_BLOCK_SIZE);
    memory_size += get_mem3Dpel(&p_Vid->p_decs->first_moment_pred_bestY_b8x8, 2, MB_BLOCK_SIZE, MB_BLOCK_SIZE);
    memory_size += get_mem2Dpel(&p_Vid->p_decs->first_moment_pred, MB_BLOCK_SIZE, MB_BLOCK_SIZE);
    memory_size += get_mem2Duint16(&p_Vid->p_decs->second_moment_bestY_mb, MB_BLOCK_SIZE, MB_BLOCK_SIZE);
    memory_size += get_mem3Duint16(&p_Vid->p_decs->second_moment_bestY_b8x8, 2, MB_BLOCK_SIZE, MB_BLOCK_SIZE);
    memory_size += get_mem3Duint16(&p_Vid->p_decs->second_moment_pred_bestY_b8x8, 2, MB_BLOCK_SIZE, MB_BLOCK_SIZE);
    memory_size += get_mem2Duint16(&p_Vid->p_decs->second_moment_pred, MB_BLOCK_SIZE, MB_BLOCK_SIZE);
    break;
  default:
    ;
  }
  return memory_size;
}

/*!
**************************************************************************************
* \brief 
*      free memory of error resilient RDO.  
**************************************************************************************
*/
void free_errdo_mem(VideoParameters *p_Vid)
{
  //for shared memory
  if (p_Vid->p_decs->res_img)
  {
    free_mem3Dint(p_Vid->p_decs->res_img);
    p_Vid->p_decs->res_img         = NULL;
  }
  if (p_Vid->p_decs->res_mb_best8x8)
  {
    free_mem3Dint(p_Vid->p_decs->res_mb_best8x8);
    p_Vid->p_decs->res_mb_best8x8         = NULL;
  }

  //for RMPC
  if (p_Vid->p_decs->RCD_bestY_mb)
  {
    free_mem2Dint(p_Vid->p_decs->RCD_bestY_mb);
    p_Vid->p_decs->RCD_bestY_mb         = NULL;
  }
  if (p_Vid->p_decs->RCD_bestY_b8x8)
  {
    free_mem3Dint(p_Vid->p_decs->RCD_bestY_b8x8);
    p_Vid->p_decs->RCD_bestY_b8x8         = NULL;
  }
  if (p_Vid->p_decs->MVCD_bestY_mb)
  {
    free_mem2Dint(p_Vid->p_decs->MVCD_bestY_mb);
    p_Vid->p_decs->MVCD_bestY_mb         = NULL;
  }
  if (p_Vid->p_decs->MVCD_bestY_b8x8)
  {
    free_mem3Dint(p_Vid->p_decs->MVCD_bestY_b8x8);
    p_Vid->p_decs->MVCD_bestY_b8x8         = NULL;
  }
  if (p_Vid->p_decs->flag_bestY_mb)
  {
    free_mem2D(p_Vid->p_decs->flag_bestY_mb);
    p_Vid->p_decs->flag_bestY_mb         = NULL;
  }
  if (p_Vid->p_decs->flag_bestY_b8x8)
  {
    free_mem3D(p_Vid->p_decs->flag_bestY_b8x8);
    p_Vid->p_decs->flag_bestY_b8x8         = NULL;
  }
  if (p_Vid->p_decs->flag_wo_res)
  {
    free_mem2D(p_Vid->p_decs->flag_wo_res);
    p_Vid->p_decs->flag_wo_res         = NULL;
  }
  if (p_Vid->p_decs->flag_wo_res_bestY_b8x8)
  {
    free_mem3D(p_Vid->p_decs->flag_wo_res_bestY_b8x8);
    p_Vid->p_decs->flag_wo_res_bestY_b8x8         = NULL;
  }
  if (p_Vid->p_decs->trans_dist_bestY_mb)
  {
    free_mem2Dint(p_Vid->p_decs->trans_dist_bestY_mb);
    p_Vid->p_decs->trans_dist_bestY_mb         = NULL;
  }
  if (p_Vid->p_decs->trans_dist_bestY_b8x8)
  {
    free_mem3Dint(p_Vid->p_decs->trans_dist_bestY_b8x8);
    p_Vid->p_decs->trans_dist_bestY_b8x8         = NULL;
  }
  if (p_Vid->p_decs->trans_dist_wo_res)
  {
    free_mem2Dint(p_Vid->p_decs->trans_dist_wo_res);
    p_Vid->p_decs->trans_dist_wo_res         = NULL;   //it is used for P8x8, where residual may be set to 0
  }
  if (p_Vid->p_decs->trans_dist_wo_res_bestY_b8x8)
  {
    free_mem3Dint(p_Vid->p_decs->trans_dist_wo_res_bestY_b8x8);
    p_Vid->p_decs->trans_dist_wo_res_bestY_b8x8   = NULL;   //it is used for P8x8, where residual may be set to 0
  }
  if (p_Vid->p_decs->trans_err_bestY_mb)
  {
    free_mem2Dint(p_Vid->p_decs->trans_err_bestY_mb);
    p_Vid->p_decs->trans_err_bestY_mb         = NULL;
  }
  if (p_Vid->p_decs->trans_err_bestY_b8x8)
  {
    free_mem3Dint(p_Vid->p_decs->trans_err_bestY_b8x8);
    p_Vid->p_decs->trans_err_bestY_b8x8         = NULL;
  }
  if (p_Vid->p_decs->trans_err_wo_res)
  {
    free_mem2Dint(p_Vid->p_decs->trans_err_wo_res);
    p_Vid->p_decs->trans_err_wo_res         = NULL;   //it is used for P8x8, where residual may be set to 0
  }
  if (p_Vid->p_decs->trans_err_wo_res_bestY_b8x8)
  {   
    free_mem3Dint(p_Vid->p_decs->trans_err_wo_res_bestY_b8x8);
    p_Vid->p_decs->trans_err_wo_res_bestY_b8x8   = NULL;   //it is used for P8x8, where residual may be set to 0
  }

  //for LLN
  if (p_Vid->p_decs->dec_mb_pred)
  {
    free_mem3Dpel(p_Vid->p_decs->dec_mb_pred);
    p_Vid->p_decs->dec_mb_pred         = NULL;
  }
  if (p_Vid->p_decs->dec_mbY_best)
  {
    free_mem3Dpel(p_Vid->p_decs->dec_mbY_best);
    p_Vid->p_decs->dec_mbY_best        = NULL;
  }
  if (p_Vid->p_decs->dec_mbY_best8x8)
  {   
    free_mem4Dpel(p_Vid->p_decs->dec_mbY_best8x8);
    p_Vid->p_decs->dec_mb_pred_best8x8 = NULL;
  }
  if (p_Vid->p_decs->dec_mb_pred_best8x8)
  {
    free_mem4Dpel(p_Vid->p_decs->dec_mb_pred_best8x8);
    p_Vid->p_decs->dec_mbY_best8x8     = NULL;
  }

  //for ROPE
  if (p_Vid->p_decs->first_moment_bestY_mb)
  {
    free_mem2Dpel(p_Vid->p_decs->first_moment_bestY_mb);
    p_Vid->p_decs->first_moment_bestY_mb         = NULL;
  }
  if (p_Vid->p_decs->first_moment_bestY_b8x8)
  {
    free_mem3Dpel(p_Vid->p_decs->first_moment_bestY_b8x8);
    p_Vid->p_decs->first_moment_bestY_b8x8       = NULL;
  }
  if (p_Vid->p_decs->first_moment_pred_bestY_b8x8)
  {
    free_mem3Dpel(p_Vid->p_decs->first_moment_pred_bestY_b8x8);
    p_Vid->p_decs->first_moment_pred_bestY_b8x8       = NULL;
  }
  if (p_Vid->p_decs->first_moment_pred)
  {
    free_mem2Dpel(p_Vid->p_decs->first_moment_pred);
    p_Vid->p_decs->first_moment_pred       = NULL;
  }
  if (p_Vid->p_decs->second_moment_bestY_mb)
  {
    free_mem2Duint16(p_Vid->p_decs->second_moment_bestY_mb);
    p_Vid->p_decs->second_moment_bestY_mb        = NULL;
  }
  if (p_Vid->p_decs->second_moment_bestY_b8x8)
  {
    free_mem3Duint16(p_Vid->p_decs->second_moment_bestY_b8x8);
    p_Vid->p_decs->second_moment_bestY_b8x8      = NULL;
  }
  if (p_Vid->p_decs->second_moment_pred_bestY_b8x8)
  {
    free_mem3Duint16(p_Vid->p_decs->second_moment_pred_bestY_b8x8);
    p_Vid->p_decs->second_moment_pred_bestY_b8x8      = NULL;
  }
  if (p_Vid->p_decs->second_moment_pred)
  {
    free_mem2Duint16(p_Vid->p_decs->second_moment_pred);
    p_Vid->p_decs->second_moment_pred      = NULL;
  }


  if ( p_Vid->p_decs != NULL )
  {
    free( p_Vid->p_decs );
    p_Vid->p_decs = NULL;
  }
}

/*!
 *************************************************************************************
* \brief 
 *    Initialize error concealment function
 *    (Currently only copy concealment is implemented. Can extend to other concealment
 *    types when available.)
 *
 *************************************************************************************
*/
void init_error_conceal(VideoParameters *p_Vid, int concealment_type)
{
  p_Vid->error_conceal_picture = copy_conceal_picture;
}


/******************************************************************************************
*
*  Finds reference picture with nearest POC to current picture to use for error concealment
*   
*******************************************************************************************
*/
StorablePicture* find_nearest_ref_picture(DecodedPictureBuffer *p_Dpb, int poc)
{
  unsigned int i;
  int min_poc_diff = 1000;
  int poc_diff;
  StorablePicture* refPic = NULL;

  for (i = 0; i < p_Dpb->ref_frames_in_buffer; i++)
  {
    if (p_Dpb->fs_ref[i]->is_used==3)
    {
      if ((p_Dpb->fs_ref[i]->frame->used_for_reference)&&(!p_Dpb->fs_ref[i]->frame->is_long_term))
      {
        poc_diff = iabs(p_Dpb->fs_ref[i]->frame->poc - poc);
        if (poc_diff < min_poc_diff)
        {
          refPic = p_Dpb->fs_ref[i]->frame;
          min_poc_diff = poc_diff;
        }
      }
    }        
  }
  return refPic;
}

/*!
 *************************************************************************************
* \brief 
 *    Gives the prediction residue for a block
 *************************************************************************************
*/
void errdo_compute_residue (Macroblock *currMB, imgpel **imgY, int **res_img, imgpel **mb_pred, int b8block, int block_size) 
{
  int i,j;
  int i0 = (b8block & 0x01)<<3,   i1 = i0+block_size;
  int j0 = (b8block >> 1)<<3,     j1 = j0+block_size;

  for (i = i0; i < i1; i++)
  {
    for (j = j0; j < j1; j++)
    {
      res_img[j][i] = (int)imgY[j][currMB->pix_x + i] - mb_pred[j][i];
    } 
  }
}


/*!
 *************************************************************************************
 * \brief
 *    Store the best 8x8 block of estimated distortion for errdo.
 *
 * \param 
 *
 * \note
 *	  For lln algorithm, we need to consider akip and direct?
 *
 *************************************************************************************
 */
void errdo_store_best_b8x8(Macroblock *currMB, int transform8x8, int block)
{
  VideoParameters *p_Vid = currMB->p_Vid;
  InputParameters *p_Inp = currMB->p_Inp;

  switch (p_Inp->de)
  {
  default:
    errdo_store_best_block_multihyp(p_Inp, p_Vid->p_decs->dec_mbY_best8x8[transform8x8], p_Vid->enc_picture->de_mem->p_dec_img[0], block, currMB->pix_x, currMB->pix_y, BLOCK_SIZE_8x8);
    errdo_store_best_block_multihyp(p_Inp, p_Vid->p_decs->dec_mb_pred_best8x8[transform8x8], p_Vid->p_decs->dec_mb_pred, block, 0, 0, BLOCK_SIZE_8x8); 
    break;
  }
}


/*!
 *************************************************************************************
 * \brief
 *    get the best 8x8 block of estimated distortion. Original code seems to have problem. Why not get p_Vid->p_decs->dec_mb_pred_best8x8?
 *	  But since there is errdo_get_best_P8x8() in SetCoeffAndReconstruction8x8(), errdo_get_best_b8x8 seems not necessary.
 *
 * \param 
 *
 * \note
 *	  For lln algorithm, we need to consider skip and direct?
 *    This function seems not necessary since we have errdo_get_best_P8x8
 *************************************************************************************
 */
void errdo_get_best_b8x8(Macroblock *currMB, int transform8x8, int block)
{
  VideoParameters *p_Vid = currMB->p_Vid;
  InputParameters *p_Inp = currMB->p_Inp;

  switch (p_Inp->de)
  {
  default:
    errdo_get_best_block_multihyp(currMB, p_Vid->enc_picture->de_mem->p_dec_img[0], p_Vid->p_decs->dec_mbY_best8x8[transform8x8], block, BLOCK_SIZE_8x8);
    break;
  }

}



/*!
 *************************************************************************************
* \brief 
 *    Store the best macroblock of estimated distortion for errdo.
*   
 * \param 
 *
 *************************************************************************************
*/
void errdo_store_best_MB(Macroblock *currMB)
{
  VideoParameters *p_Vid = currMB->p_Vid;
  InputParameters *p_Inp = currMB->p_Inp;

  switch (p_Inp->de)
  {
  default:
    errdo_store_best_block_multihyp(p_Inp, p_Vid->p_decs->dec_mbY_best, p_Vid->enc_picture->de_mem->p_dec_img[0], 0, currMB->pix_x, currMB->pix_y, MB_BLOCK_SIZE);
    break;
  }
}



/*!
 *************************************************************************************
 * \brief
 *    Get the best macroblock of estimated distortion for storable picture.
 *
 * \param 
 *
 *************************************************************************************
 */
void errdo_get_best_MB(Macroblock *currMB)
{   
  VideoParameters *p_Vid = currMB->p_Vid;
  InputParameters *p_Inp = currMB->p_Inp;

  switch (p_Inp->de)
  {
  default:
    errdo_get_best_block_multihyp(currMB, p_Vid->enc_picture->de_mem->p_dec_img[0], p_Vid->p_decs->dec_mbY_best, 0, MB_BLOCK_SIZE);
    break;
  }

}



/*!
 *************************************************************************************
 * \brief
 *    Store the best macroblock of estimated distortion for mode P8x8.
 *
 * \param 
 *
 *************************************************************************************
 */
void errdo_get_best_P8x8(Macroblock *currMB, int transform8x8)
{   
  VideoParameters *p_Vid = currMB->p_Vid;
  InputParameters *p_Inp = currMB->p_Inp;

  switch (p_Inp->de)
  {
  default:
    if (p_Vid->p_decs->rec_type == 0)
    {
      errdo_get_best_block_multihyp(currMB, p_Vid->enc_picture->de_mem->p_dec_img[0], p_Vid->p_decs->dec_mb_pred_best8x8[transform8x8], 0, MB_BLOCK_SIZE);
    }
    else
    {
      errdo_get_best_block_multihyp(currMB, p_Vid->enc_picture->de_mem->p_dec_img[0], p_Vid->p_decs->dec_mbY_best8x8[transform8x8], 0, MB_BLOCK_SIZE);
    }
  }
}


//Zhifeng 090611
/*!
 *************************************************************************************
 * \brief
 *    Initialize distortion estimation algorithm
 *    (Can extend to support more algorithms when available.)
 *  
 *************************************************************************************
 */
void init_distortion_estimation(VideoParameters *p_Vid, int de_algorithm)
{
  switch (de_algorithm)
  {
  default:
    p_Vid->estimate_distortion = errdo_distortion_estimation_multihyp;
    break;
  }
}


/*!
 *************************************************************************************
 * \brief
 *    allocate storable picture memory for errdo
*
 *************************************************************************************
*/
void errdo_alloc_storable_picture(StorablePicture *p, VideoParameters *p_Vid, InputParameters *p_Inp, int size_x, int size_y, int size_x_cr, int size_y_cr)
{
  Dist_Estm *s;
  int   dec, ndec, nplane;

  p->de_mem = (Dist_Estm *)malloc( sizeof(Dist_Estm) );
  s = p->de_mem;

  s->res_con_diff_Y   = NULL;
  s->res_con_diff_UV   = NULL;
  s->MV_con_diff_Y   = NULL;
  s->MV_con_diff_UV   = NULL;
  s->error_sign_flag_Y   = NULL;
  s->error_sign_flag_UV   = NULL;
  s->transmission_dist_Y   = NULL;
  s->transmission_dist_UV   = NULL;
  s->transmission_err_Y   = NULL;
  s->transmission_err_UV   = NULL;
  s->dec_imgY   = NULL;
  s->dec_imgUV  = NULL;
  s->mb_error_map = NULL;
  s->first_moment_Y   = NULL;
  s->first_moment_UV  = NULL;
  s->second_moment_Y   = NULL;
  s->second_moment_UV  = NULL;
  for (nplane = 0; nplane < 3; nplane++)
  {
    s->p_res_con_diff[nplane] = NULL;
    s->p_MV_con_diff[nplane] = NULL;
    s->p_error_sign_flag[nplane] = NULL;
    s->p_transmission_dist[nplane] = NULL;
    s->p_transmission_err[nplane] = NULL;
    s->p_dec_img[nplane] = NULL;
    s->p_first_moment[nplane] = NULL;
    s->p_second_moment[nplane] = NULL;
  }

  switch (p_Inp->de)
  {
  case RMPC:  //RMPC support slice data partitioning, so it keeps residual and MV in their respective memory
    get_mem2Dint(&(s->transmission_dist_Y), size_y, size_x);
    s->p_transmission_dist[0] = s->transmission_dist_Y;
    if (p_Vid->yuv_format != YUV400)
    {
      get_mem3Dint(&(s->transmission_dist_UV), 2, size_y_cr, size_x_cr);
      s->p_transmission_dist[1] = s->transmission_dist_UV[0];
      s->p_transmission_dist[2] = s->transmission_dist_UV[1];
    }
#ifdef RMPC_PAPER
    //allocate memory for the residual concealment difference, which is used in the reference paper to calculate propagation factor in the next frame
    get_mem2Dint(&(s->res_con_diff_Y), size_y, size_x);
    s->p_res_con_diff[0] = s->res_con_diff_Y;
    if (p_Vid->yuv_format != YUV400)
    {
      get_mem3Dint(&(s->res_con_diff_UV), 2, size_y_cr, size_x_cr);
      s->p_res_con_diff[1] = s->res_con_diff_UV[0];
      s->p_res_con_diff[2] = s->res_con_diff_UV[1];
    }

    //allocate memory for the MV concealment difference, which is used in the reference paper to calculate propagation factor in the next frame
    get_mem2Dint(&(s->MV_con_diff_Y), size_y, size_x);
    s->p_MV_con_diff[0] = s->MV_con_diff_Y;
    if (p_Vid->yuv_format != YUV400)
    {
      get_mem3Dint(&(s->MV_con_diff_UV), 2, size_y_cr, size_x_cr);
      s->p_MV_con_diff[1] = s->MV_con_diff_UV[0];
      s->p_MV_con_diff[2] = s->MV_con_diff_UV[1];
    }
#endif  //#ifdef RMPC_PAPER
#ifdef RMPC_ERRDO
    //allocate memory for the transmission error sign, which is used as an alternative to calculate propagation factor in the next frame
    get_mem2D(&(s->error_sign_flag_Y), size_y, size_x);
    s->p_error_sign_flag[0] = s->error_sign_flag_Y;
    if (p_Vid->yuv_format != YUV400)
    {
      get_mem3D(&(s->error_sign_flag_UV), 2, size_y_cr, size_x_cr);
      s->p_error_sign_flag[1] = s->error_sign_flag_UV[0];
      s->p_error_sign_flag[2] = s->error_sign_flag_UV[1];
    }
#endif  //#ifdef RMPC_ERRDO

    break;
  case ERMPC_NN:
  case ERMPC_FAST:
  case EXTENDED_RMPC: //Extended RMPC further considers interpolation filter, B-slice, and de-blocking filter
    get_mem2Dint(&(s->transmission_dist_Y), size_y, size_x);
    s->p_transmission_dist[0] = s->transmission_dist_Y;
    if (p_Vid->yuv_format != YUV400)
    {
      get_mem3Dint(&(s->transmission_dist_UV), 2, size_y_cr, size_x_cr);
      s->p_transmission_dist[1] = s->transmission_dist_UV[0];
      s->p_transmission_dist[2] = s->transmission_dist_UV[1];
    }
    get_mem2Dint(&(s->transmission_err_Y), size_y, size_x);
    s->p_transmission_err[0] = s->transmission_err_Y;
    if (p_Vid->yuv_format != YUV400)
    {
      get_mem3Dint(&(s->transmission_err_UV), 2, size_y_cr, size_x_cr);
      s->p_transmission_err[1] = s->transmission_err_UV[0];
      s->p_transmission_err[2] = s->transmission_err_UV[1];
    }
    break;

  case LLN:
  case FAST_LLN:
    ndec = p_Inp->NoOfDecoders;
    //check the consistent
    if (ndec == 0)
    {
      printf("ndec can not be zero for LLN and fast LLN algorithms, reset ndec to 30");
      ndec = 30;
    }
    get_mem3D(&(s->mb_error_map), ndec, size_y/MB_BLOCK_SIZE, size_x/MB_BLOCK_SIZE);
    get_mem3Dpel(&(s->dec_imgY), ndec, size_y, size_x);

    // This seems somewhat inefficient. Why not allocate array as [ndec][x] where x goes from 0 to 2?
    if ((s->p_dec_img[0] = (imgpel***)calloc(ndec,sizeof(imgpel**))) == NULL)
    {
      no_mem_exit("errdo.c: p_dec_img[0]");
    }

    if (p_Vid->yuv_format != YUV400)
    {
      get_mem4Dpel(&(s->dec_imgUV), ndec, 2, size_y_cr, size_x_cr);
      if ((s->p_dec_img[1] = (imgpel***)calloc(ndec,sizeof(imgpel**))) == NULL)
      {  
        no_mem_exit("errdo.c: p_dec_img[1]");
      }
      if ((s->p_dec_img[2] = (imgpel***)calloc(ndec,sizeof(imgpel**))) == NULL)
      {
        no_mem_exit("errdo.c: p_dec_img[2]");
      }
    }

    for (dec = 0; dec < ndec; dec++)
    {
      s->p_dec_img[0][dec] = s->dec_imgY[dec];
    }

    if (p_Vid->yuv_format != YUV400)
    {
      for (dec = 0; dec < ndec; dec++)
      {
        s->p_dec_img[1][dec] = s->dec_imgUV[dec][0];
        s->p_dec_img[2][dec] = s->dec_imgUV[dec][1];
      }
    }

    break;
  case ROPE:
    //allocate memory for the first moment
    get_mem2Dpel(&(s->first_moment_Y), size_y, size_x);
    s->p_first_moment[0] = s->first_moment_Y;

    if (p_Vid->yuv_format != YUV400)
    {
      get_mem3Dpel(&(s->first_moment_UV), 2, size_y_cr, size_x_cr);
      s->p_first_moment[1] = s->first_moment_UV[0];
      s->p_first_moment[2] = s->first_moment_UV[1];
    }

    //allocate memory for the second moment
    get_mem2Duint16(&(s->second_moment_Y), size_y, size_x);
    s->p_second_moment[0] = s->second_moment_Y;

    if (p_Vid->yuv_format != YUV400)
    {
      get_mem3Duint16(&(s->second_moment_UV), 2, size_y_cr, size_x_cr);
      s->p_second_moment[1] = s->second_moment_UV[0];
      s->p_second_moment[2] = s->second_moment_UV[1];
    }
    break;

  case LTI:
    break;

  default:
    ;
  }
}


/*!
 *************************************************************************************
 * \brief
 *    free storable picture memory for errdo
 *
 *************************************************************************************
 */
void errdo_free_storable_picture(StorablePicture* s)
{
  int nplane;
  Dist_Estm *p = s->de_mem;
  //free memory for RMPC and extended RMPC algorithms
  if (p->res_con_diff_Y)
  {
    free_mem2Dint(p->res_con_diff_Y);
    p->res_con_diff_Y   = NULL;
  }
  if (p->res_con_diff_UV)
  {
    free_mem3Dint(p->res_con_diff_UV);
    p->res_con_diff_UV   = NULL;
  }
  for (nplane = 0; nplane < 3; nplane++)
  {
    p->p_res_con_diff[nplane] = NULL;
  }

  if (p->MV_con_diff_Y)
  {
    free_mem2Dint(p->MV_con_diff_Y);
    p->MV_con_diff_Y   = NULL;
  }
  if (p->MV_con_diff_UV)
  {
    free_mem3Dint(p->MV_con_diff_UV);
    p->MV_con_diff_UV   = NULL;
  }
  for (nplane = 0; nplane < 3; nplane++)
  {
    p->p_MV_con_diff[nplane] = NULL;
  }

  if (p->error_sign_flag_Y)
  {
    free_mem2D(p->error_sign_flag_Y);
    p->error_sign_flag_Y   = NULL;
  } 
  if (p->error_sign_flag_UV)
  {
    free_mem3D(p->error_sign_flag_UV);
    p->error_sign_flag_UV   = NULL;
  }
  for (nplane = 0; nplane < 3; nplane++)
  {
    p->p_error_sign_flag[nplane] = NULL;
  }

  if (p->transmission_dist_Y)
  {
    free_mem2Dint(p->transmission_dist_Y);
    p->transmission_dist_Y   = NULL;
  }
  if (p->transmission_dist_UV)
  {
    free_mem3Dint(p->transmission_dist_UV);
    p->transmission_dist_UV   = NULL;
  }
  for (nplane = 0; nplane < 3; nplane++)
  {
    p->p_transmission_dist[nplane] = NULL;
  }

  if (p->transmission_err_Y)
  {
    free_mem2Dint(p->transmission_err_Y);
    p->transmission_err_Y   = NULL;
  }
  if (p->transmission_err_UV)
  {
    free_mem3Dint(p->transmission_err_UV);
    p->transmission_err_UV   = NULL;
  }
  for (nplane = 0; nplane < 3; nplane++)
  {
    p->p_transmission_err[nplane] = NULL;
  }

  //free memory for LLN and Faster LLN algorithms
  if (p->dec_imgY)
  {
    free_mem3Dpel(p->dec_imgY);
    p->dec_imgY   = NULL;
  }
  if (p->dec_imgUV)
  {
    free_mem4Dpel(p->dec_imgUV);
    p->dec_imgUV  = NULL;
  }
  for (nplane = 0; nplane < 3; nplane++)
  {
    if (p->p_dec_img[nplane])
    {  
      free(p->p_dec_img[nplane]);
      p->p_dec_img[nplane] = NULL;
    }
  }
  if (p->mb_error_map)
  {
    free_mem3D(p->mb_error_map);
    p->mb_error_map = NULL;
  }


  //free memory for ROPE and extended ROPE algorithms
  if (p->first_moment_Y)
  {
    free_mem2Dpel(p->first_moment_Y);
    p->first_moment_Y   = NULL;
  }
  if (p->first_moment_UV)
  {
    free_mem3Dpel(p->first_moment_UV);
    p->first_moment_UV  = NULL;
  }
  for (nplane = 0; nplane < 3; nplane++)
  {
    p->p_first_moment[nplane] = NULL;
  }
  if (p->second_moment_Y)
  {
    free_mem2Duint16(p->second_moment_Y);
    p->second_moment_Y   = NULL;
  }
  if (p->second_moment_UV)
  {
    free_mem3Duint16(p->second_moment_UV);
    p->second_moment_UV  = NULL;
  }
  for (nplane = 0; nplane < 3; nplane++)
  {
    p->p_second_moment[nplane] = NULL;
  }

  if (s->de_mem)
  {
    free(s->de_mem);
    s->de_mem   = NULL;
  }
}

