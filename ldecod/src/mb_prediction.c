/*!
 *************************************************************************************
 * \file mb_prediction.c
 *
 * \brief
 *    Macroblock prediction functions
 *
 * \author
 *    Main contributors (see contributors.h for copyright, address and affiliation details)
 *    - Alexis Michael Tourapis         <alexismt@ieee.org>
 *************************************************************************************
 */

#include "contributors.h"

#include "block.h"
#include "global.h"
#include "mbuffer.h"
#include "elements.h"
#include "errorconcealment.h"
#include "macroblock.h"
#include "fmo.h"
#include "cabac.h"
#include "vlc.h"
#include "image.h"
#include "mb_access.h"
#include "biaridecod.h"
#include "transform8x8.h"
#include "transform.h"
#include "mc_prediction.h"
#include "quant.h"
#include "intra4x4_pred.h"
#include "intra8x8_pred.h"
#include "intra16x16_pred.h"
#include "mv_prediction.h"
#include "mb_prediction.h"


int mb_pred_intra4x4(Macroblock *currMB, ColorPlane curr_plane, imgpel **currImg, StorablePicture *dec_picture)
{
  Slice *currSlice = currMB->p_Slice;
  int yuv = dec_picture->chroma_format_idc - 1;
  int i=0, j=0,k, j4=0,i4=0;  
  int j_pos, i_pos;
  int ioff,joff;
  int block8x8;   // needed for ABT
  currMB->itrans_4x4 = (currMB->is_lossless == FALSE) ? itrans4x4 : Inv_Residual_trans_4x4;    

  for (block8x8 = 0; block8x8 < 4; block8x8++)
  {
    for (k = block8x8 * 4; k < block8x8 * 4 + 4; k ++)
    {
      i =  (decode_block_scan[k] & 3);
      j = ((decode_block_scan[k] >> 2) & 3);

      ioff = (i << 2);
      joff = (j << 2);
      i4   = currMB->block_x + i;
      j4   = currMB->block_y + j;
      j_pos = j4 * BLOCK_SIZE;
      i_pos = i4 * BLOCK_SIZE;

      // PREDICTION
      //===== INTRA PREDICTION =====
      if (intrapred(currMB, curr_plane, ioff,joff,i4,j4) == SEARCH_SYNC)  /* make 4x4 prediction block mpr from given prediction p_Vid->mb_mode */
        return SEARCH_SYNC;                   /* bit error */
      // =============== 4x4 itrans ================
      // -------------------------------------------
      currMB->itrans_4x4  (currMB, curr_plane, ioff, joff);

      copy_image_data_4x4(&currImg[j_pos], &currSlice->mb_rec[curr_plane][joff], i_pos, ioff);
    }
  }

  // chroma decoding *******************************************************
  if ((dec_picture->chroma_format_idc != YUV400) && (dec_picture->chroma_format_idc != YUV444)) 
  {
    intra_cr_decoding(currMB, yuv);
  }

  return 1;
}


int mb_pred_intra16x16(Macroblock *currMB, ColorPlane curr_plane, imgpel **currImg, StorablePicture *dec_picture)
{
  int yuv = dec_picture->chroma_format_idc - 1;

  intrapred16x16(currMB, curr_plane, currMB->i16mode);
  currMB->ipmode_DPCM = (char) currMB->i16mode; //For residual DPCM
  // =============== 4x4 itrans ================
  // -------------------------------------------
  iMBtrans4x4(currMB, curr_plane, 0);

  // chroma decoding *******************************************************
  if ((dec_picture->chroma_format_idc != YUV400) && (dec_picture->chroma_format_idc != YUV444)) 
  {
    intra_cr_decoding(currMB, yuv);
  }
  return 1;
}

int mb_pred_intra8x8(Macroblock *currMB, ColorPlane curr_plane, imgpel **currImg, StorablePicture *dec_picture)
{
  Slice *currSlice = currMB->p_Slice;
  int yuv = dec_picture->chroma_format_idc - 1;

  int block8x8;   // needed for ABT
  currMB->itrans_8x8 = (currMB->is_lossless == FALSE) ? itrans8x8 : Inv_Residual_trans_8x8;

  for (block8x8 = 0; block8x8 < 4; block8x8++)
  {
    //=========== 8x8 BLOCK TYPE ============
    int ioff = (block8x8 & 0x01) << 3;
    int joff = (block8x8 >> 1  ) << 3;

    //PREDICTION
    intrapred8x8(currMB, curr_plane, ioff, joff);
    currMB->itrans_8x8  (currMB, curr_plane, ioff,joff);      // use DCT transform and make 8x8 block m7 from prediction block mpr

    copy_image_data_8x8(&currImg[currMB->pix_y + joff], &currSlice->mb_rec[curr_plane][joff], currMB->pix_x + ioff, ioff);
  }
  // chroma decoding *******************************************************
  if ((dec_picture->chroma_format_idc != YUV400) && (dec_picture->chroma_format_idc != YUV444)) 
  {
    intra_cr_decoding(currMB, yuv);
  }
  return 1;
}


static void set_chroma_vector(Macroblock *currMB, int *list_offset)
{
  Slice *currSlice = currMB->p_Slice;
  VideoParameters *p_Vid = currMB->p_Vid;

  if (!currSlice->mb_aff_frame_flag)
  {
    if(p_Vid->structure == TOP_FIELD)
    {
      int k,l;  
      for (l = LIST_0; l <= (LIST_1); l++)
      {
        for(k = 0; k < p_Vid->listXsize[l]; k++)
        {
          if(p_Vid->structure != p_Vid->listX[l][k]->structure)
            p_Vid->listX[l][k]->chroma_vector_adjustment = -2;
          else
            p_Vid->listX[l][k]->chroma_vector_adjustment= 0;
        }
      }
    }
    else if(p_Vid->structure == BOTTOM_FIELD)
    {
      int k,l;  
      for (l = LIST_0; l <= (LIST_1); l++)
      {
        for(k = 0; k < p_Vid->listXsize[l]; k++)
        {
          if (p_Vid->structure != p_Vid->listX[l][k]->structure)
            p_Vid->listX[l][k]->chroma_vector_adjustment = 2;
          else
            p_Vid->listX[l][k]->chroma_vector_adjustment= 0;
        }
      }
    }
    else
    {
      int k,l;  
      for (l = LIST_0; l <= (LIST_1); l++)
      {
        for(k = 0; k < p_Vid->listXsize[l]; k++)
        {
          p_Vid->listX[l][k]->chroma_vector_adjustment= 0;
        }
      }
    }
  }
  else
  {
    int mb_nr = (currMB->mbAddrX & 0x01);
    int k,l;  

    //////////////////////////
    // find out the correct list offsets
    if (currMB->mb_field)
    {
      *list_offset = mb_nr ? 4 : 2;

      for (l = LIST_0 + *list_offset; l <= (LIST_1 + *list_offset); l++)
      {
        for(k = 0; k < p_Vid->listXsize[l]; k++)
        {          
          if(mb_nr == 0 && p_Vid->listX[l][k]->structure == BOTTOM_FIELD)
            p_Vid->listX[l][k]->chroma_vector_adjustment = -2;
          else if(mb_nr == 1 && p_Vid->listX[l][k]->structure == TOP_FIELD)
            p_Vid->listX[l][k]->chroma_vector_adjustment = 2;
          else
            p_Vid->listX[l][k]->chroma_vector_adjustment= 0;
        }
      }
    }
    else
    {
      for (l = LIST_0; l <= (LIST_1); l++)
      {
        for(k = 0; k < p_Vid->listXsize[l]; k++)
        {
          p_Vid->listX[l][k]->chroma_vector_adjustment= 0;
        }
      }
    }
  }

  p_Vid->max_mb_vmv_r = (p_Vid->structure != FRAME || (currSlice->mb_aff_frame_flag && currMB->mb_field)) ? p_Vid->max_vmv_r >> 1 : p_Vid->max_vmv_r;
}

int mb_pred_skip(Macroblock *currMB, ColorPlane curr_plane, imgpel **currImg, StorablePicture *dec_picture)
{
  Slice *currSlice = currMB->p_Slice;
  VideoParameters *p_Vid = currMB->p_Vid;
  int curr_mb_field = ((currSlice->mb_aff_frame_flag)&&(currMB->mb_field));

  int list_offset = 0;

  set_chroma_vector(currMB, &list_offset);

  perform_mc(currMB, curr_plane, dec_picture, LIST_0, 0, 0, list_offset, MB_BLOCK_SIZE, MB_BLOCK_SIZE, curr_mb_field);

  copy_image_data_16x16(&currImg[currMB->pix_y], currSlice->mb_pred[curr_plane], currMB->pix_x, 0);

  if ((dec_picture->chroma_format_idc != YUV400) && (dec_picture->chroma_format_idc != YUV444)) 
  {

    copy_image_data(&dec_picture->imgUV[0][currMB->pix_c_y], currSlice->mb_pred[1], currMB->pix_c_x, 0, p_Vid->mb_size[1][0], p_Vid->mb_size[1][1]);
    copy_image_data(&dec_picture->imgUV[1][currMB->pix_c_y], currSlice->mb_pred[2], currMB->pix_c_x, 0, p_Vid->mb_size[1][0], p_Vid->mb_size[1][1]);
  }
  return 1;
}

int mb_pred_sp_skip(Macroblock *currMB, ColorPlane curr_plane, imgpel **currImg, StorablePicture *dec_picture)
{
  Slice *currSlice = currMB->p_Slice;
  int curr_mb_field = ((currSlice->mb_aff_frame_flag)&&(currMB->mb_field));
  int list_offset = 0;

  set_chroma_vector(currMB, &list_offset);

  perform_mc(currMB, curr_plane, dec_picture, LIST_0, 0, 0, list_offset, MB_BLOCK_SIZE, MB_BLOCK_SIZE, curr_mb_field);	
  iTransform(currMB, curr_plane, 1);
  return 1;
}

int mb_pred_p_inter8x8(Macroblock *currMB, ColorPlane curr_plane, imgpel **currImg, StorablePicture *dec_picture)
{
  int block8x8;   // needed for ABT
  int i=0, j=0,k;  

  Slice *currSlice = currMB->p_Slice;
  VideoParameters *p_Vid = currMB->p_Vid;
  int smb = p_Vid->type == SP_SLICE && IS_INTER(currMB);
  int curr_mb_field = ((currSlice->mb_aff_frame_flag)&&(currMB->mb_field));

  int list_offset = 0;

  set_chroma_vector(currMB, &list_offset);

  for (block8x8=0; block8x8<4; block8x8++)
  {
    int mv_mode  = currMB->b8mode[block8x8];
    int pred_dir = currMB->b8pdir[block8x8];

    int k_start = (block8x8 << 2);
    int k_inc = (mv_mode == SMB8x4) ? 2 : 1;
    int k_end = (mv_mode == SMB8x8) ? k_start + 1 : ((mv_mode == SMB4x4) ? k_start + 4 : k_start + k_inc + 1);

    int block_size_x = ( mv_mode == SMB8x4 || mv_mode == SMB8x8 ) ? SMB_BLOCK_SIZE : BLOCK_SIZE;
    int block_size_y = ( mv_mode == SMB4x8 || mv_mode == SMB8x8 ) ? SMB_BLOCK_SIZE : BLOCK_SIZE;

    for (k = k_start; k < k_end; k += k_inc)
    {
      i =  (decode_block_scan[k] & 3);
      j = ((decode_block_scan[k] >> 2) & 3);
      perform_mc(currMB, curr_plane, dec_picture, pred_dir, i, j, list_offset, block_size_x, block_size_y, curr_mb_field);
    }        
  }

  iTransform(currMB, curr_plane, smb); 
  return 1;
}

int mb_pred_p_inter16x16(Macroblock *currMB, ColorPlane curr_plane, imgpel **currImg, StorablePicture *dec_picture)
{
  int smb = (currMB->p_Vid->type == SP_SLICE);
  Slice *currSlice = currMB->p_Slice;
  int curr_mb_field = ((currSlice->mb_aff_frame_flag)&&(currMB->mb_field));
  int list_offset = 0;

  set_chroma_vector(currMB, &list_offset);

  perform_mc(currMB, curr_plane, dec_picture, currMB->b8pdir[0], 0, 0, list_offset, MB_BLOCK_SIZE, MB_BLOCK_SIZE, curr_mb_field);
  iTransform(currMB, curr_plane, smb);
  return 1;
}

int mb_pred_p_inter16x8(Macroblock *currMB, ColorPlane curr_plane, imgpel **currImg, StorablePicture *dec_picture)
{
  int smb = (currMB->p_Vid->type == SP_SLICE);
  Slice *currSlice = currMB->p_Slice;
  int curr_mb_field = ((currSlice->mb_aff_frame_flag)&&(currMB->mb_field));
  int list_offset = 0;

  set_chroma_vector(currMB, &list_offset);

  perform_mc(currMB, curr_plane, dec_picture, currMB->b8pdir[0], 0, 0, list_offset, MB_BLOCK_SIZE, BLOCK_SIZE_8x8, curr_mb_field);
  perform_mc(currMB, curr_plane, dec_picture, currMB->b8pdir[2], 0, 2, list_offset, MB_BLOCK_SIZE, BLOCK_SIZE_8x8, curr_mb_field);
  iTransform(currMB, curr_plane, smb); 
  return 1;
}

int mb_pred_p_inter8x16(Macroblock *currMB, ColorPlane curr_plane, imgpel **currImg, StorablePicture *dec_picture)
{
  int smb = (currMB->p_Vid->type == SP_SLICE);
  Slice *currSlice = currMB->p_Slice;
  int curr_mb_field = ((currSlice->mb_aff_frame_flag)&&(currMB->mb_field));
  int list_offset = 0;

  set_chroma_vector(currMB, &list_offset);

  perform_mc(currMB, curr_plane, dec_picture, currMB->b8pdir[0], 0, 0, list_offset, BLOCK_SIZE_8x8, MB_BLOCK_SIZE, curr_mb_field);
  perform_mc(currMB, curr_plane, dec_picture, currMB->b8pdir[1], 2, 0, list_offset, BLOCK_SIZE_8x8, MB_BLOCK_SIZE, curr_mb_field);
  iTransform(currMB, curr_plane, smb);
  return 1;
}

int mb_pred_b_dtemporal(Macroblock *currMB, ColorPlane curr_plane, imgpel **currImg, StorablePicture *dec_picture)
{
  short ref_idx;
  int refList;

  PicMotionParams *motion = &dec_picture->motion;
  int block_size_x, block_size_y;
  int k;
  int block8x8;   // needed for ABT
  Slice *currSlice = currMB->p_Slice;
  VideoParameters *p_Vid = currMB->p_Vid;
  int curr_mb_field = ((currSlice->mb_aff_frame_flag)&&(currMB->mb_field));

  MotionParams *colocated = &currSlice->p_colocated->frame;
  int list_offset = 0;

  set_chroma_vector(currMB, &list_offset);

  if (currMB->mb_field)
  {
    if(currMB->mbAddrX & 0x01)
    {
      colocated = &currSlice->p_colocated->bottom;
    }
    else
    {
      colocated = &currSlice->p_colocated->top;
    }
  }

  for (block8x8=0; block8x8<4; block8x8++)
  {      
    int pred_dir = currMB->b8pdir[block8x8];

    int k_start = (block8x8 << 2);
    int k_end = k_start;

    if (p_Vid->active_sps->direct_8x8_inference_flag)
    {
      block_size_x = SMB_BLOCK_SIZE;
      block_size_y = SMB_BLOCK_SIZE;
      k_end ++;
    }
    else
    {
      block_size_x = BLOCK_SIZE;
      block_size_y = BLOCK_SIZE;
      k_end += BLOCK_MULTIPLE;
    }

    for (k = k_start; k < k_start + BLOCK_MULTIPLE; k ++)
    {

      int i =  (decode_block_scan[k] & 3);
      int j = ((decode_block_scan[k] >> 2) & 3);
      int i4   = currMB->block_x + i;
      int j4   = currMB->block_y + j;
      int j6   = currMB->block_y_aff + j;
      assert (pred_dir<=2);

      refList = (colocated->ref_idx[LIST_0][j6][i4]== -1 ? LIST_1 : LIST_0);
      ref_idx =  colocated->ref_idx[refList][j6][i4];

      if(ref_idx==-1) // co-located is intra mode
      {
        memset( &motion->mv  [LIST_0][j4][i4][0], 0, 2* sizeof(short));
        memset( &motion->mv  [LIST_1][j4][i4][0], 0, 2* sizeof(short));

        motion->ref_idx[LIST_0][j4][i4] = 0;
        motion->ref_idx[LIST_1][j4][i4] = 0;
      }
      else // co-located skip or inter mode
      {
        int mapped_idx=0;
        int iref;

        for (iref=0;iref<imin(currSlice->num_ref_idx_l0_active,p_Vid->listXsize[LIST_0 + list_offset]);iref++)
        {
          if(p_Vid->structure==0 && curr_mb_field==0)
          {
            // If the current MB is a frame MB and the colocated is from a field picture,
            // then the colocated->ref_pic_id may have been generated from the wrong value of
            // frame_poc if it references it's complementary field, so test both POC values
            if(p_Vid->listX[0][iref]->top_poc*2 == colocated->ref_pic_id[refList][j6][i4] || p_Vid->listX[0][iref]->bottom_poc*2 == colocated->ref_pic_id[refList][j6][i4])
            {
              mapped_idx=iref;
              break;
            }
            else //! invalid index. Default to zero even though this case should not happen
              mapped_idx=INVALIDINDEX;
            continue;
          }

          if (dec_picture->ref_pic_num[p_Vid->current_slice_nr][LIST_0 + list_offset][iref]==colocated->ref_pic_id[refList][j6][i4])
          {
            mapped_idx=iref;
            break;
          }
          else //! invalid index. Default to zero even though this case should not happen
          {
            mapped_idx=INVALIDINDEX;
          }
        }
        if (INVALIDINDEX == mapped_idx)
        {
          error("temporal direct error: colocated block has ref that is unavailable",-1111);
        }
        else
        {
          int mv_scale = currSlice->mvscale[LIST_0 + list_offset][mapped_idx];

          //! In such case, an array is needed for each different reference.
          if (mv_scale == 9999 || p_Vid->listX[LIST_0+list_offset][mapped_idx]->is_long_term)
          {
            memcpy(&motion->mv  [LIST_0][j4][i4][0], &colocated->mv[refList][j6][i4][0], 2 * sizeof(short));
            memset(&motion->mv  [LIST_1][j4][i4][0], 0, 2 * sizeof(short));
          }
          else
          {
            motion->mv  [LIST_0][j4][i4][0]= (short) ((mv_scale * colocated->mv[refList][j6][i4][0] + 128 ) >> 8);
            motion->mv  [LIST_0][j4][i4][1]= (short) ((mv_scale * colocated->mv[refList][j6][i4][1] + 128 ) >> 8);

            motion->mv  [LIST_1][j4][i4][0]= (short) (motion->mv[LIST_0][j4][i4][0] - colocated->mv[refList][j6][i4][0]);
            motion->mv  [LIST_1][j4][i4][1]= (short) (motion->mv[LIST_0][j4][i4][1] - colocated->mv[refList][j6][i4][1]);
          }

          motion->ref_idx[LIST_0][j4][i4] = (char) mapped_idx; //p_Vid->listX[1][0]->ref_idx[refList][j4][i4];
          motion->ref_idx[LIST_1][j4][i4] = 0;
        }
      }
      // store reference picture ID determined by direct mode
      motion->ref_pic_id[LIST_0][j4][i4] = dec_picture->ref_pic_num[p_Vid->current_slice_nr][LIST_0 + list_offset][(short)motion->ref_idx[LIST_0][j4][i4]];
      motion->ref_pic_id[LIST_1][j4][i4] = dec_picture->ref_pic_num[p_Vid->current_slice_nr][LIST_1 + list_offset][(short)motion->ref_idx[LIST_1][j4][i4]];
    }
    for (k = k_start; k < k_end; k ++)
    {
      int i =  (decode_block_scan[k] & 3);
      int j = ((decode_block_scan[k] >> 2) & 3);
      perform_mc(currMB, curr_plane, dec_picture, pred_dir, i, j, list_offset, block_size_x, block_size_y, curr_mb_field);           
    }
  }

  if (currMB->cbp == 0)
  {
    copy_image_data_16x16(&currImg[currMB->pix_y], currSlice->mb_pred[curr_plane], currMB->pix_x, 0);

    if ((dec_picture->chroma_format_idc != YUV400) && (dec_picture->chroma_format_idc != YUV444)) 
    {
      copy_image_data(&dec_picture->imgUV[0][currMB->pix_c_y], currSlice->mb_pred[1], currMB->pix_c_x, 0, p_Vid->mb_size[1][0], p_Vid->mb_size[1][1]);
      copy_image_data(&dec_picture->imgUV[1][currMB->pix_c_y], currSlice->mb_pred[2], currMB->pix_c_x, 0, p_Vid->mb_size[1][0], p_Vid->mb_size[1][1]);
    }
  }
  else
    iTransform(currMB, curr_plane, 0); 
  return 1;
}

int mb_pred_b_dspatial(Macroblock *currMB, ColorPlane curr_plane, imgpel **currImg, StorablePicture *dec_picture)
{
  char l0_rFrame = -1, l1_rFrame = -1;
  PicMotionParams *motion = &dec_picture->motion;
  short pmvl0[2]={0,0}, pmvl1[2]={0,0};
  int k;
  int block8x8;
  Slice *currSlice = currMB->p_Slice;
  VideoParameters *p_Vid = currMB->p_Vid;
  int curr_mb_field = ((currSlice->mb_aff_frame_flag)&&(currMB->mb_field));

  MotionParams *colocated = &currSlice->p_colocated->frame;
  int list_offset = 0;
  int pred_dir = 0;
  
  Boolean has_zero_partitions = FALSE;
  int64 *ref_pic_num_l0, *ref_pic_num_l1;

  set_chroma_vector(currMB, &list_offset);

  if (currMB->mb_field)
  {
    if(currMB->mbAddrX & 0x01)
    {
      colocated = &currSlice->p_colocated->bottom;
    }
    else
    {
      colocated = &currSlice->p_colocated->top;
    }
  }

  prepare_direct_params(currMB, dec_picture, pmvl0, pmvl1, &l0_rFrame, &l1_rFrame);
  
  ref_pic_num_l0 = dec_picture->ref_pic_num[p_Vid->current_slice_nr][LIST_0 + list_offset];
  ref_pic_num_l1 = dec_picture->ref_pic_num[p_Vid->current_slice_nr][LIST_1 + list_offset];

  if (p_Vid->active_sps->direct_8x8_inference_flag)
  {
    int block_size_x = SMB_BLOCK_SIZE;
    int block_size_y = SMB_BLOCK_SIZE;
    for (block8x8 = 0; block8x8 < 4; block8x8++)
    {
      int k_start = (block8x8 << 2);

      for (k = k_start; k < k_start + BLOCK_MULTIPLE; k ++)
      {
        int i  =  (decode_block_scan[k] & 3);
        int j  = ((decode_block_scan[k] >> 2) & 3);
        int i4  = currMB->block_x + i;
        int j4  = currMB->block_y + j;
        int j6  = currMB->block_y_aff + j;

        //printf("%d %d\n", i, j);

        //===== DIRECT PREDICTION =====

        if (l0_rFrame >=0)
        {
          if (!l0_rFrame  && ((!colocated->moving_block[j6][i4]) && (!p_Vid->listX[LIST_1 + list_offset][0]->is_long_term)))
          {
            has_zero_partitions = TRUE;
            motion->mv     [LIST_0][j4][i4][0] = 0;
            motion->mv     [LIST_0][j4][i4][1] = 0;
            motion->ref_idx[LIST_0][j4][i4]    = 0;
          }
          else
          {
            has_zero_partitions = TRUE;
            motion->mv     [LIST_0][j4][i4][0] = pmvl0[0];
            motion->mv     [LIST_0][j4][i4][1] = pmvl0[1];
            motion->ref_idx[LIST_0][j4][i4]    = l0_rFrame;
          }
        }
        else
        {        
          motion->mv     [LIST_0][j4][i4][0] = 0;
          motion->mv     [LIST_0][j4][i4][1] = 0;
          motion->ref_idx[LIST_0][j4][i4]    = -1;
        }

        if (l1_rFrame >=0)
        {
          if  (l1_rFrame == 0 && ((!colocated->moving_block[j6][i4]) && (!p_Vid->listX[LIST_1 + list_offset][0]->is_long_term)))
          {
            has_zero_partitions = TRUE;
            motion->mv     [LIST_1][j4][i4][0] = 0;
            motion->mv     [LIST_1][j4][i4][1] = 0;
            motion->ref_idx[LIST_1][j4][i4]    = 0;
          }
          else
          {
            has_zero_partitions = TRUE;
            motion->mv     [LIST_1][j4][i4][0] = pmvl1[0];
            motion->mv     [LIST_1][j4][i4][1] = pmvl1[1];
            motion->ref_idx[LIST_1][j4][i4]    = l1_rFrame;
          }
        }
        else
        {
          motion->mv     [LIST_1][j4][i4][0] = 0;
          motion->mv     [LIST_1][j4][i4][1] = 0;
          motion->ref_idx[LIST_1][j4][i4]    = -1;
        }

        if (l1_rFrame == -1) 
            pred_dir = 0;
        else if (l0_rFrame == -1) 
          pred_dir = 1;
        else
          pred_dir = 2;

        if (l0_rFrame < 0 && l1_rFrame < 0)
        {
          motion->ref_idx[LIST_0][j4][i4] = 0;
          motion->ref_idx[LIST_1][j4][i4] = 0;
          pred_dir = 2;
        }

        motion->ref_pic_id[LIST_0][j4][i4] = ref_pic_num_l0[(short)motion->ref_idx[LIST_0][j4][i4]];
        motion->ref_pic_id[LIST_1][j4][i4] = ref_pic_num_l1[(short)motion->ref_idx[LIST_1][j4][i4]];
      }

      if (has_zero_partitions == TRUE)
      {
        int i =  (decode_block_scan[k_start] & 3);
        int j = ((decode_block_scan[k_start] >> 2) & 3);

        perform_mc(currMB, curr_plane, dec_picture, pred_dir, i, j, list_offset, block_size_x, block_size_y, curr_mb_field);           
      }
    }
  }
  else
  {
    int block_size_x = BLOCK_SIZE;
    int block_size_y = BLOCK_SIZE;
    for (block8x8 = 0; block8x8 < 4; block8x8++)
    {
      int k_start = (block8x8 << 2);
      int k_end = k_start + BLOCK_MULTIPLE;

      for (k = k_start; k < k_end; k ++)
      {
        int i  =  (decode_block_scan[k] & 3);
        int j  = ((decode_block_scan[k] >> 2) & 3);
        int i4  = currMB->block_x + i;
        int j4  = currMB->block_y + j;
        int j6  = currMB->block_y_aff + j;

        //===== DIRECT PREDICTION =====

        if (l0_rFrame >=0)
        {
          if (!l0_rFrame  && ((!colocated->moving_block[j6][i4]) && (!p_Vid->listX[LIST_1 + list_offset][0]->is_long_term)))
          {
            has_zero_partitions = TRUE;
            motion->mv     [LIST_0][j4][i4][0] = 0;
            motion->mv     [LIST_0][j4][i4][1] = 0;
            motion->ref_idx[LIST_0][j4][i4]    = 0;
          }
          else
          {
            has_zero_partitions = TRUE;
            motion->mv     [LIST_0][j4][i4][0] = pmvl0[0];
            motion->mv     [LIST_0][j4][i4][1] = pmvl0[1];
            motion->ref_idx[LIST_0][j4][i4]    = l0_rFrame;
          }
        }
        else
        {        
          motion->mv     [LIST_0][j4][i4][0] = 0;
          motion->mv     [LIST_0][j4][i4][1] = 0;
          motion->ref_idx[LIST_0][j4][i4]    = -1;
        }

        if (l1_rFrame >=0)
        {
          if  (l1_rFrame == 0 && ((!colocated->moving_block[j6][i4]) && (!p_Vid->listX[LIST_1 + list_offset][0]->is_long_term)))
          {
            has_zero_partitions = TRUE;
            motion->mv     [LIST_1][j4][i4][0] = 0;
            motion->mv     [LIST_1][j4][i4][1] = 0;
            motion->ref_idx[LIST_1][j4][i4]    = 0;
          }
          else
          {
            has_zero_partitions = TRUE;
            motion->mv     [LIST_1][j4][i4][0] = pmvl1[0];
            motion->mv     [LIST_1][j4][i4][1] = pmvl1[1];
            motion->ref_idx[LIST_1][j4][i4]    = l1_rFrame;
          }
        }
        else
        {
          motion->mv     [LIST_1][j4][i4][0] = 0;
          motion->mv     [LIST_1][j4][i4][1] = 0;
          motion->ref_idx[LIST_1][j4][i4]    = -1;
        }

        if (l0_rFrame < 0 && l1_rFrame < 0)
        {
          motion->ref_idx[LIST_0][j4][i4] = 0;
          motion->ref_idx[LIST_1][j4][i4] = 0;
        }

        if (l1_rFrame == -1) 
        {
          if (l0_rFrame == -1) 
            pred_dir = 2;
          else
            pred_dir = 0;
        }
        else if (l0_rFrame == -1) 
        {
          pred_dir = 1;
        }
        else                                               
          pred_dir = 2;

        motion->ref_pic_id[LIST_0][j4][i4] = ref_pic_num_l0[(short)motion->ref_idx[LIST_0][j4][i4]];
        motion->ref_pic_id[LIST_1][j4][i4] = ref_pic_num_l1[(short)motion->ref_idx[LIST_1][j4][i4]];
      }

      if (has_zero_partitions == TRUE)
      {
        for (k = k_start; k < k_end; k ++)
        {        
          int i =  (decode_block_scan[k] & 3);
          int j = ((decode_block_scan[k] >> 2) & 3);

          perform_mc(currMB, curr_plane, dec_picture, pred_dir, i, j, list_offset, block_size_x, block_size_y, curr_mb_field);           
        }
      }
    }
  }

  if (has_zero_partitions == FALSE)
  {
    perform_mc(currMB, curr_plane, dec_picture, pred_dir, 0, 0, list_offset, MB_BLOCK_SIZE, MB_BLOCK_SIZE, curr_mb_field);           
  }  

  if (currMB->cbp == 0)
  {
    copy_image_data_16x16(&currImg[currMB->pix_y], currSlice->mb_pred[curr_plane], currMB->pix_x, 0);

    if ((dec_picture->chroma_format_idc != YUV400) && (dec_picture->chroma_format_idc != YUV444)) 
    {
      copy_image_data(&dec_picture->imgUV[0][currMB->pix_c_y], currSlice->mb_pred[1], currMB->pix_c_x, 0, p_Vid->mb_size[1][0], p_Vid->mb_size[1][1]);
      copy_image_data(&dec_picture->imgUV[1][currMB->pix_c_y], currSlice->mb_pred[2], currMB->pix_c_x, 0, p_Vid->mb_size[1][0], p_Vid->mb_size[1][1]);
    }
  }
  else
    iTransform(currMB, curr_plane, 0); 
  return 1;
}

int mb_pred_b_inter8x8(Macroblock *currMB, ColorPlane curr_plane, imgpel **currImg, StorablePicture *dec_picture)
{
  short ref_idx;
  int refList;

  char l0_rFrame = -1, l1_rFrame = -1;
  PicMotionParams *motion = &dec_picture->motion;
  short pmvl0[2]={0,0}, pmvl1[2]={0,0};
  int block_size_x, block_size_y;
  int k;
  int block8x8;   // needed for ABT
  Slice *currSlice = currMB->p_Slice;
  VideoParameters *p_Vid = currMB->p_Vid;
  int curr_mb_field = ((currSlice->mb_aff_frame_flag)&&(currMB->mb_field));

  MotionParams *colocated = &currSlice->p_colocated->frame;
  int list_offset = 0;

  set_chroma_vector(currMB, &list_offset);

  if (currMB->mb_field)
  {
    if(currMB->mbAddrX & 0x01)
    {
      colocated = &currSlice->p_colocated->bottom;
    }
    else
    {
      colocated = &currSlice->p_colocated->top;
    }
  }

  // prepare direct modes
  if (currSlice->direct_spatial_mv_pred_flag && (!(currMB->b8mode[0] && currMB->b8mode[1] && currMB->b8mode[2] && currMB->b8mode[3])))
    prepare_direct_params(currMB, dec_picture, pmvl0, pmvl1, &l0_rFrame, &l1_rFrame);

  for (block8x8=0; block8x8<4; block8x8++)
  {
    int mv_mode  = currMB->b8mode[block8x8];
    int pred_dir = currMB->b8pdir[block8x8];
    
    if ( mv_mode != 0 )
    {
      int k_start = (block8x8 << 2);
      int k_inc = (mv_mode == SMB8x4) ? 2 : 1;
      int k_end = (mv_mode == SMB8x8) ? k_start + 1 : ((mv_mode == SMB4x4) ? k_start + 4 : k_start + k_inc + 1);

      block_size_x = ( mv_mode == SMB8x4 || mv_mode == SMB8x8 ) ? SMB_BLOCK_SIZE : BLOCK_SIZE;
      block_size_y = ( mv_mode == SMB4x8 || mv_mode == SMB8x8 ) ? SMB_BLOCK_SIZE : BLOCK_SIZE;

      for (k = k_start; k < k_end; k += k_inc)
      {
        int i =  (decode_block_scan[k] & 3);
        int j = ((decode_block_scan[k] >> 2) & 3);
        perform_mc(currMB, curr_plane, dec_picture, pred_dir, i, j, list_offset, block_size_x, block_size_y, curr_mb_field);
      }        
    }
    else
    {
      int k_start = (block8x8 << 2);
      int k_end = k_start;

      if (p_Vid->active_sps->direct_8x8_inference_flag)
      {
        block_size_x = SMB_BLOCK_SIZE;
        block_size_y = SMB_BLOCK_SIZE;
        k_end ++;
      }
      else
      {
        block_size_x = BLOCK_SIZE;
        block_size_y = BLOCK_SIZE;
        k_end += BLOCK_MULTIPLE;
      }

      // Prepare mvs (needed for deblocking and mv prediction
      if (currSlice->direct_spatial_mv_pred_flag)
      {
        int64 *ref_pic_num_l0 = dec_picture->ref_pic_num[p_Vid->current_slice_nr][LIST_0 + list_offset];
        int64 *ref_pic_num_l1 = dec_picture->ref_pic_num[p_Vid->current_slice_nr][LIST_1 + list_offset];

        for (k = k_start; k < k_start + BLOCK_MULTIPLE; k ++)
        {

          int i  =  (decode_block_scan[k] & 3);
          int j  = ((decode_block_scan[k] >> 2) & 3);
          int i4  = currMB->block_x + i;
          int j4  = currMB->block_y + j;
          int j6  = currMB->block_y_aff + j;

          assert (pred_dir<=2);


          //===== DIRECT PREDICTION =====

          if (l0_rFrame >=0)
          {
            if (!l0_rFrame  && ((!colocated->moving_block[j6][i4]) && (!p_Vid->listX[LIST_1 + list_offset][0]->is_long_term)))
            {
              motion->mv  [LIST_0][j4][i4][0] = 0;
              motion->mv  [LIST_0][j4][i4][1] = 0;
              motion->ref_idx[LIST_0][j4][i4] = 0;
            }
            else
            {
              motion->mv  [LIST_0][j4][i4][0] = pmvl0[0];
              motion->mv  [LIST_0][j4][i4][1] = pmvl0[1];
              motion->ref_idx[LIST_0][j4][i4] = l0_rFrame;
            }
          }
          else
          {
            motion->ref_idx[LIST_0][j4][i4] = -1;
            motion->mv  [LIST_0][j4][i4][0] = 0;
            motion->mv  [LIST_0][j4][i4][1] = 0;
          }

          if (l1_rFrame >=0)
          {
            if  (l1_rFrame==0 && ((!colocated->moving_block[j6][i4]) && (!p_Vid->listX[LIST_1 + list_offset][0]->is_long_term)))
            {

              motion->mv  [LIST_1][j4][i4][0] = 0;
              motion->mv  [LIST_1][j4][i4][1] = 0;
              motion->ref_idx[LIST_1][j4][i4] = l1_rFrame;
            }
            else
            {
              motion->mv  [LIST_1][j4][i4][0] = pmvl1[0];
              motion->mv  [LIST_1][j4][i4][1] = pmvl1[1];
              motion->ref_idx[LIST_1][j4][i4] = l1_rFrame;
            }
          }
          else
          {
            motion->mv  [LIST_1][j4][i4][0] = 0;
            motion->mv  [LIST_1][j4][i4][1] = 0;
            motion->ref_idx[LIST_1][j4][i4] = -1;
          }

          if (l0_rFrame < 0 && l1_rFrame < 0)
          {
            motion->ref_idx[LIST_0][j4][i4] = 0;
            motion->ref_idx[LIST_1][j4][i4] = 0;
          }

          if      (motion->ref_idx[LIST_1][j4][i4]==-1) 
          {
            pred_dir = 0;
            ref_idx  = (motion->ref_idx[LIST_0][j4][i4] != -1) ? motion->ref_idx[LIST_0][j4][i4] : 0;
          }
          else if (motion->ref_idx[LIST_0][j4][i4]==-1) 
          {
            pred_dir = 1;
            ref_idx  = (motion->ref_idx[LIST_1][j4][i4] != -1) ? motion->ref_idx[LIST_1][j4][i4] : 0;
          }
          else                                               
            pred_dir = 2;

          motion->ref_pic_id[LIST_0][j4][i4] = ref_pic_num_l0[(short)motion->ref_idx[LIST_0][j4][i4]];
          motion->ref_pic_id[LIST_1][j4][i4] = ref_pic_num_l1[(short)motion->ref_idx[LIST_1][j4][i4]];
        }
      }
      else
      {
        for (k = k_start; k < k_start + BLOCK_MULTIPLE; k ++)
        {

          int i =  (decode_block_scan[k] & 3);
          int j = ((decode_block_scan[k] >> 2) & 3);
          int i4   = currMB->block_x + i;
          int j4   = currMB->block_y + j;
          int j6   = currMB->block_y_aff + j;

          assert (pred_dir<=2);

          refList = (colocated->ref_idx[LIST_0][j6][i4]== -1 ? LIST_1 : LIST_0);
          ref_idx =  colocated->ref_idx[refList][j6][i4];

          if(ref_idx==-1) // co-located is intra mode
          {
            memset( &motion->mv  [LIST_0][j4][i4][0], 0, 2* sizeof(short));
            memset( &motion->mv  [LIST_1][j4][i4][0], 0, 2* sizeof(short));

            motion->ref_idx[LIST_0][j4][i4] = 0;
            motion->ref_idx[LIST_1][j4][i4] = 0;
          }
          else // co-located skip or inter mode
          {
            int mapped_idx=0;
            int iref;

            for (iref=0;iref<imin(currSlice->num_ref_idx_l0_active,p_Vid->listXsize[LIST_0 + list_offset]);iref++)
            {
              if(p_Vid->structure==0 && curr_mb_field==0)
              {
                // If the current MB is a frame MB and the colocated is from a field picture,
                // then the colocated->ref_pic_id may have been generated from the wrong value of
                // frame_poc if it references it's complementary field, so test both POC values
                if(p_Vid->listX[0][iref]->top_poc*2 == colocated->ref_pic_id[refList][j6][i4] || p_Vid->listX[0][iref]->bottom_poc*2 == colocated->ref_pic_id[refList][j6][i4])
                {
                  mapped_idx=iref;
                  break;
                }
                else //! invalid index. Default to zero even though this case should not happen
                  mapped_idx=INVALIDINDEX;
                continue;
              }

              if (dec_picture->ref_pic_num[p_Vid->current_slice_nr][LIST_0 + list_offset][iref]==colocated->ref_pic_id[refList][j6][i4])
              {
                mapped_idx=iref;
                break;
              }
              else //! invalid index. Default to zero even though this case should not happen
              {
                mapped_idx=INVALIDINDEX;
              }
            }
            if (INVALIDINDEX == mapped_idx)
            {
              error("temporal direct error: colocated block has ref that is unavailable",-1111);
            }
            else
            {
              int mv_scale = currSlice->mvscale[LIST_0 + list_offset][mapped_idx];

              //! In such case, an array is needed for each different reference.
              if (mv_scale == 9999 || p_Vid->listX[LIST_0+list_offset][mapped_idx]->is_long_term)
              {
                memcpy(&motion->mv  [LIST_0][j4][i4][0], &colocated->mv[refList][j6][i4][0], 2 * sizeof(short));
                memset(&motion->mv  [LIST_1][j4][i4][0], 0, 2 * sizeof(short));
              }
              else
              {
                motion->mv  [LIST_0][j4][i4][0]= (short) ((mv_scale * colocated->mv[refList][j6][i4][0] + 128 ) >> 8);
                motion->mv  [LIST_0][j4][i4][1]= (short) ((mv_scale * colocated->mv[refList][j6][i4][1] + 128 ) >> 8);

                motion->mv  [LIST_1][j4][i4][0]= (short) (motion->mv[LIST_0][j4][i4][0] - colocated->mv[refList][j6][i4][0]);
                motion->mv  [LIST_1][j4][i4][1]= (short) (motion->mv[LIST_0][j4][i4][1] - colocated->mv[refList][j6][i4][1]);
              }

              motion->ref_idx[LIST_0][j4][i4] = (char) mapped_idx; //p_Vid->listX[1][0]->ref_idx[refList][j4][i4];
              motion->ref_idx[LIST_1][j4][i4] = 0;
            }
          }
          // store reference picture ID determined by direct mode
          motion->ref_pic_id[LIST_0][j4][i4] = dec_picture->ref_pic_num[p_Vid->current_slice_nr][LIST_0 + list_offset][(short)motion->ref_idx[LIST_0][j4][i4]];
          motion->ref_pic_id[LIST_1][j4][i4] = dec_picture->ref_pic_num[p_Vid->current_slice_nr][LIST_1 + list_offset][(short)motion->ref_idx[LIST_1][j4][i4]];
        }
      }

      for (k = k_start; k < k_end; k ++)
      {
        int i =  (decode_block_scan[k] & 3);
        int j = ((decode_block_scan[k] >> 2) & 3);
        perform_mc(currMB, curr_plane, dec_picture, pred_dir, i, j, list_offset, block_size_x, block_size_y, curr_mb_field);           
      }
    }
  }

  iTransform(currMB, curr_plane, 0); 
  return 1;
}

/*!
 ************************************************************************
 * \brief
 *    Copy IPCM coefficients to decoded picture buffer and set parameters for this MB
 *    (for IPCM CABAC and IPCM CAVLC  28/11/2003)
 *
 * \author
 *    Dong Wang <Dong.Wang@bristol.ac.uk>
 ************************************************************************
 */

int mb_pred_ipcm(Macroblock *currMB)
{
  int i, j, k;
  Slice *currSlice = currMB->p_Slice;
  VideoParameters *p_Vid = currMB->p_Vid;
  StorablePicture *dec_picture = p_Vid->dec_picture;

  //Copy coefficients to decoded picture buffer
  //IPCM coefficients are stored in currSlice->cof which is set in function read_IPCM_coeffs_from_NAL()

  for(i = 0; i < MB_BLOCK_SIZE; ++i)
  {
    for(j = 0;j < MB_BLOCK_SIZE ; ++j)
    {
      dec_picture->imgY[currMB->pix_y + i][currMB->pix_x + j] = (imgpel) currSlice->cof[0][i][j];
    }
  }

  if ((dec_picture->chroma_format_idc != YUV400) && !IS_INDEPENDENT(p_Vid))
  {
    for (k = 0; k < 2; ++k)
    {
      for(i = 0; i < p_Vid->mb_cr_size_y; ++i)
      {
        for(j = 0;j < p_Vid->mb_cr_size_x; ++j)
        {
          dec_picture->imgUV[k][currMB->pix_c_y+i][currMB->pix_c_x + j] = (imgpel) currSlice->cof[k + 1][i][j];  
        }
      }
    }
  }

  // for deblocking filter
  update_qp(currMB, 0);

  // for CAVLC: Set the nz_coeff to 16.
  // These parameters are to be used in CAVLC decoding of neighbour blocks  
  memset(&p_Vid->nz_coeff[currMB->mbAddrX][0][0][0], 16, 3 * BLOCK_SIZE * BLOCK_SIZE * sizeof(byte));

  // for CABAC decoding of MB skip flag
  currMB->skip_flag = 0;

  //for deblocking filter CABAC
  currMB->cbp_blk[0] = 0xFFFF;

  //For CABAC decoding of Dquant
  currSlice->last_dquant = 0;
  return 1;
}

