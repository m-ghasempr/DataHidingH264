
/*!
 ***********************************************************************
 * \file image.c
 *
 * \brief
 *    Decode a Slice
 *
 * \author
 *    Main contributors (see contributors.h for copyright, address and affiliation details)
 *    - Inge Lille-Langoy               <inge.lille-langoy@telenor.com>
 *    - Rickard Sjoberg                 <rickard.sjoberg@era.ericsson.se>
 *    - Jani Lainema                    <jani.lainema@nokia.com>
 *    - Sebastian Purreiter             <sebastian.purreiter@mch.siemens.de>
 *    - Byeong-Moon Jeon                <jeonbm@lge.com>
 *    - Thomas Wedi                     <wedi@tnt.uni-hannover.de>
 *    - Gabi Blaettermann
 *    - Ye-Kui Wang                     <wyk@ieee.org>
 *    - Antti Hallapuro                 <antti.hallapuro@nokia.com>
 *    - Alexis Tourapis                 <alexismt@ieee.org>
 *    - Jill Boyce                      <jill.boyce@thomson.net>
 *    - Saurav K Bandyopadhyay          <saurav@ieee.org>
 *    - Zhenyu Wu                       <Zhenyu.Wu@thomson.net
 *    - Purvin Pandit                   <Purvin.Pandit@thomson.net>
 *
 ***********************************************************************
 */

#include "contributors.h"

#include <math.h>
#include <limits.h>

#include "global.h"
#include "image.h"
#include "fmo.h"
#include "annexb.h"
#include "nalu.h"
#include "parset.h"
#include "header.h"

#include "sei.h"
#include "output.h"
#include "mb_access.h"
#include "memalloc.h"
#include "macroblock.h"

#include "loopfilter.h"

#include "biaridecod.h"
#include "context_ini.h"
#include "cabac.h"
#include "vlc.h"
#include "quant.h"

#include "errorconcealment.h"
#include "erc_api.h"

extern int testEndian(void);

/*!
 ************************************************************************
 * \brief
 *    Initializes the parameters for a new picture
 ************************************************************************
 */
static void init_picture(ImageParameters *p_Img, Slice *currSlice, InputParameters *p_Inp)
{
  int i;
  int nplane;

  if (p_Img->dec_picture)
  {
    // this may only happen on slice loss
    exit_picture(p_Img, &p_Img->dec_picture);
  }
  if (p_Img->recovery_point)
    p_Img->recovery_frame_num = (p_Img->frame_num + p_Img->recovery_frame_cnt) % p_Img->MaxFrameNum;

  if (p_Img->idr_flag)
    p_Img->recovery_frame_num = p_Img->frame_num;

  if (p_Img->recovery_point == 0 &&
    p_Img->frame_num != p_Img->pre_frame_num &&
    p_Img->frame_num != (p_Img->pre_frame_num + 1) % p_Img->MaxFrameNum)
  {
    if (p_Img->active_sps->gaps_in_frame_num_value_allowed_flag == 0)
    {
      // picture error concealment
      if(p_Inp->conceal_mode !=0)
      {
        if((p_Img->frame_num) < ((p_Img->pre_frame_num + 1) % p_Img->MaxFrameNum))
        {
          /* Conceal lost IDR frames and any frames immediately
             following the IDR. Use frame copy for these since
             lists cannot be formed correctly for motion copy*/
          p_Img->conceal_mode = 1;
          p_Img->IDR_concealment_flag = 1;
          conceal_lost_frames(p_Img);
          //reset to original concealment mode for future drops
          p_Img->conceal_mode = p_Inp->conceal_mode;
        }
        else
        {
          //reset to original concealment mode for future drops
          p_Img->conceal_mode = p_Inp->conceal_mode;

          p_Img->IDR_concealment_flag = 0;
          conceal_lost_frames(p_Img);
        }
      }
      else
      {   /* Advanced Error Concealment would be called here to combat unintentional loss of pictures. */
        error("An unintentional loss of pictures occurs! Exit\n", 100);
      }
    }
    if(p_Img->conceal_mode == 0)
      fill_frame_num_gap(p_Img);
  }

  if(p_Img->nal_reference_idc)
  {
    p_Img->pre_frame_num = p_Img->frame_num;
  }

  p_Img->num_dec_mb = 0;

  //calculate POC
  decode_poc(p_Img);

  if (p_Img->recovery_frame_num == p_Img->frame_num &&
    p_Img->recovery_poc == 0x7fffffff)
    p_Img->recovery_poc = p_Img->framepoc;

  if(p_Img->nal_reference_idc)
    p_Img->last_ref_pic_poc = p_Img->framepoc;

  //  dumppoc (p_Img);

  if (p_Img->structure==FRAME ||p_Img->structure==TOP_FIELD)
  {
    gettime (&(p_Img->start_time));             // start time
  }

  p_Img->dec_picture = alloc_storable_picture (p_Img, (PictureStructure) p_Img->structure, p_Img->width, p_Img->height, p_Img->width_cr, p_Img->height_cr);
  p_Img->dec_picture->top_poc=p_Img->toppoc;
  p_Img->dec_picture->bottom_poc=p_Img->bottompoc;
  p_Img->dec_picture->frame_poc=p_Img->framepoc;
  p_Img->dec_picture->qp = p_Img->qp;
  p_Img->dec_picture->slice_qp_delta = currSlice->slice_qp_delta;
  p_Img->dec_picture->chroma_qp_offset[0] = p_Img->active_pps->chroma_qp_index_offset;
  p_Img->dec_picture->chroma_qp_offset[1] = p_Img->active_pps->second_chroma_qp_index_offset;

  // reset all variables of the error concealment instance before decoding of every frame.
  // here the third parameter should, if perfectly, be equal to the number of slices per frame.
  // using little value is ok, the code will allocate more memory if the slice number is larger
  ercReset(p_Img->erc_errorVar, p_Img->PicSizeInMbs, p_Img->PicSizeInMbs, p_Img->dec_picture->size_x);
  p_Img->erc_mvperMB = 0;

  switch (p_Img->structure )
  {
  case TOP_FIELD:
    {
      p_Img->dec_picture->poc=p_Img->toppoc;
      p_Img->number *= 2;
      break;
    }
  case BOTTOM_FIELD:
    {
      p_Img->dec_picture->poc=p_Img->bottompoc;
      p_Img->number = p_Img->number * 2 + 1;
      break;
    }
  case FRAME:
    {
      p_Img->dec_picture->poc=p_Img->framepoc;
      break;
    }
  default:
    error("p_Img->structure not initialized", 235);
  }

  p_Img->current_slice_nr=0;

  if (p_Img->type > SI_SLICE)
  {
    set_ec_flag(p_Img, SE_PTYPE);
    p_Img->type = P_SLICE;  // concealed element
  }

  // CAVLC init
  if (p_Img->active_pps->entropy_coding_mode_flag == CAVLC)
  {
    memset(&p_Img->nz_coeff[0][0][0][0], -1, p_Img->PicSizeInMbs * 48 *sizeof(byte)); // 3 * 4 * 4
  }

  if(p_Img->active_pps->constrained_intra_pred_flag)
  {
    for (i=0; i<(int)p_Img->PicSizeInMbs; ++i)
    {
      p_Img->intra_block[i] = 1;
    }
  }

  // Set the slice_nr member of each MB to -1, to ensure correct when packet loss occurs
  // TO set Macroblock Map (mark all MBs as 'have to be concealed')
  if( IS_INDEPENDENT(p_Img) )
  {
    for( nplane=0; nplane<MAX_PLANE; ++nplane )
    {      
      for(i=0; i<(int)p_Img->PicSizeInMbs; ++i)
      {
        p_Img->mb_data_JV[nplane][i].slice_nr = -1; 
        p_Img->mb_data_JV[nplane][i].ei_flag = 1;
        p_Img->mb_data_JV[nplane][i].dpl_flag = 0;
      }
    }
  }
  else
  {
    for(i=0; i<(int)p_Img->PicSizeInMbs; ++i)
    {
      p_Img->mb_data[i].slice_nr = -1; 
      p_Img->mb_data[i].ei_flag = 1;
      p_Img->mb_data[i].dpl_flag = 0;
    }
  }

  p_Img->dec_picture->slice_type = p_Img->type;
  p_Img->dec_picture->used_for_reference = (p_Img->nal_reference_idc != 0);
  p_Img->dec_picture->idr_flag = p_Img->idr_flag;
  p_Img->dec_picture->no_output_of_prior_pics_flag = p_Img->no_output_of_prior_pics_flag;
  p_Img->dec_picture->long_term_reference_flag = p_Img->long_term_reference_flag;
  p_Img->dec_picture->adaptive_ref_pic_buffering_flag = p_Img->adaptive_ref_pic_buffering_flag;

  p_Img->dec_picture->dec_ref_pic_marking_buffer = p_Img->dec_ref_pic_marking_buffer;
  p_Img->dec_ref_pic_marking_buffer = NULL;

  p_Img->dec_picture->MbaffFrameFlag = p_Img->MbaffFrameFlag;
  p_Img->dec_picture->PicWidthInMbs = p_Img->PicWidthInMbs;

  p_Img->get_mb_block_pos = p_Img->dec_picture->MbaffFrameFlag ? get_mb_block_pos_mbaff : get_mb_block_pos_normal;
  p_Img->getNeighbour = p_Img->dec_picture->MbaffFrameFlag ? getAffNeighbour : getNonAffNeighbour;

  p_Img->dec_picture->pic_num = p_Img->frame_num;
  p_Img->dec_picture->frame_num = p_Img->frame_num;

  p_Img->dec_picture->recovery_frame = (p_Img->frame_num == p_Img->recovery_frame_num);

  p_Img->dec_picture->coded_frame = (p_Img->structure==FRAME);

  p_Img->dec_picture->chroma_format_idc = p_Img->active_sps->chroma_format_idc;

  p_Img->dec_picture->frame_mbs_only_flag = p_Img->active_sps->frame_mbs_only_flag;
  p_Img->dec_picture->frame_cropping_flag = p_Img->active_sps->frame_cropping_flag;

  if (p_Img->dec_picture->frame_cropping_flag)
  {
    p_Img->dec_picture->frame_cropping_rect_left_offset   = p_Img->active_sps->frame_cropping_rect_left_offset;
    p_Img->dec_picture->frame_cropping_rect_right_offset  = p_Img->active_sps->frame_cropping_rect_right_offset;
    p_Img->dec_picture->frame_cropping_rect_top_offset    = p_Img->active_sps->frame_cropping_rect_top_offset;
    p_Img->dec_picture->frame_cropping_rect_bottom_offset = p_Img->active_sps->frame_cropping_rect_bottom_offset;
  }

#if (ENABLE_OUTPUT_TONEMAPPING)
  // store the necessary tone mapping sei into StorablePicture structure
  p_Img->dec_picture->seiHasTone_mapping = 0;

  if (p_Img->seiToneMapping->seiHasTone_mapping)
  {
    p_Img->dec_picture->seiHasTone_mapping    = 1;
    p_Img->dec_picture->tone_mapping_model_id = p_Img->seiToneMapping->model_id;
    p_Img->dec_picture->tonemapped_bit_depth  = p_Img->seiToneMapping->sei_bit_depth;
    p_Img->dec_picture->tone_mapping_lut      = malloc(sizeof(int)*(1<<p_Img->seiToneMapping->coded_data_bit_depth));
    if (NULL == p_Img->dec_picture->tone_mapping_lut)
    {
      no_mem_exit("init_picture: tone_mapping_lut");
    }
    memcpy(p_Img->dec_picture->tone_mapping_lut, p_Img->seiToneMapping->lut, sizeof(imgpel)*(1<<p_Img->seiToneMapping->coded_data_bit_depth));
    update_tone_mapping_sei(p_Img->seiToneMapping);
  }
#endif

  if( IS_INDEPENDENT(p_Img) )
  {
    p_Img->dec_picture_JV[0] = p_Img->dec_picture;
    p_Img->dec_picture_JV[1] = alloc_storable_picture (p_Img, (PictureStructure) p_Img->structure, p_Img->width, p_Img->height, p_Img->width_cr, p_Img->height_cr);
    copy_dec_picture_JV( p_Img, p_Img->dec_picture_JV[1], p_Img->dec_picture_JV[0] );
    p_Img->dec_picture_JV[2] = alloc_storable_picture (p_Img, (PictureStructure) p_Img->structure, p_Img->width, p_Img->height, p_Img->width_cr, p_Img->height_cr);
    copy_dec_picture_JV( p_Img, p_Img->dec_picture_JV[2], p_Img->dec_picture_JV[0] );
  }
}

void MbAffPostProc(ImageParameters *p_Img)
{

  imgpel temp[32][16];

  StorablePicture *dec_picture = p_Img->dec_picture;
  imgpel ** imgY  = dec_picture->imgY;
  imgpel ***imgUV = dec_picture->imgUV;

  short i, y, x0, y0, uv;
  for (i=0; i<(int)dec_picture->PicSizeInMbs; i+=2)
  {
    if (dec_picture->motion.mb_field[i])
    {
      get_mb_pos(p_Img, i, p_Img->mb_size[IS_LUMA], &x0, &y0);
      for (y=0; y<(2*MB_BLOCK_SIZE);++y)
        memcpy(temp[y], &imgY[y0+y][x0], MB_BLOCK_SIZE * sizeof(imgpel));

      for (y=0; y<MB_BLOCK_SIZE; ++y)
      {
         memcpy(&imgY[y0+(2*y  )][x0], temp[y              ], MB_BLOCK_SIZE * sizeof(imgpel));
         memcpy(&imgY[y0+(2*y+1)][x0], temp[y+MB_BLOCK_SIZE], MB_BLOCK_SIZE * sizeof(imgpel));
      }
       

      if (dec_picture->chroma_format_idc != YUV400)
      {
        x0 = x0 / (16/p_Img->mb_cr_size_x);
        y0 = y0 / (16/p_Img->mb_cr_size_y);

        for (uv=0; uv<2; ++uv)
        {
          for (y=0; y<(2*p_Img->mb_cr_size_y);++y)
            memcpy(temp[y], &imgUV[uv][y0+y][x0], p_Img->mb_cr_size_x * sizeof(imgpel));

          for (y=0; y<p_Img->mb_cr_size_y;++y)
          {
            memcpy(&imgUV[uv][y0+(2*y  )][x0], temp[y                  ], p_Img->mb_cr_size_x * sizeof(imgpel));
            memcpy(&imgUV[uv][y0+(2*y+1)][x0], temp[y+p_Img->mb_cr_size_y], p_Img->mb_cr_size_x * sizeof(imgpel));
          }
        }
      }
    }
  }
}

static void fill_wp_params(Slice *currSlice)
{
  ImageParameters *p_Img = currSlice->p_Img;
  
  int i, j, k;
  int comp;
  int log_weight_denom;
  int tb, td;
  int bframe = (p_Img->type==B_SLICE);
  int tx,DistScaleFactor;

  int max_l0_ref = currSlice->num_ref_idx_l0_active;
  int max_l1_ref = currSlice->num_ref_idx_l1_active;

  if (p_Img->active_pps->weighted_bipred_idc == 2 && bframe)
  {
    currSlice->luma_log2_weight_denom = 5;
    currSlice->chroma_log2_weight_denom = 5;
    currSlice->wp_round_luma   = 16;
    currSlice->wp_round_chroma = 16;

    for (i=0; i<MAX_REFERENCE_PICTURES; ++i)
    {
      for (comp=0; comp<3; ++comp)
      {
        log_weight_denom = (comp == 0) ? currSlice->luma_log2_weight_denom : currSlice->chroma_log2_weight_denom;
        currSlice->wp_weight[0][i][comp] = 1<<log_weight_denom;
        currSlice->wp_weight[1][i][comp] = 1<<log_weight_denom;
        currSlice->wp_offset[0][i][comp] = 0;
        currSlice->wp_offset[1][i][comp] = 0;
      }
    }
  }

  if (bframe)
  {
    for (i=0; i<max_l0_ref; ++i)
    {
      for (j=0; j<max_l1_ref; ++j)
      {
        for (comp = 0; comp<3; ++comp)
        {
          log_weight_denom = (comp == 0) ? currSlice->luma_log2_weight_denom : currSlice->chroma_log2_weight_denom;
          if (p_Img->active_pps->weighted_bipred_idc == 1)
          {
            currSlice->wbp_weight[0][i][j][comp] =  currSlice->wp_weight[0][i][comp];
            currSlice->wbp_weight[1][i][j][comp] =  currSlice->wp_weight[1][j][comp];
          }
          else if (p_Img->active_pps->weighted_bipred_idc == 2)
          {
            td = iClip3(-128,127,p_Img->listX[LIST_1][j]->poc - p_Img->listX[LIST_0][i]->poc);
            if (td == 0 || p_Img->listX[LIST_1][j]->is_long_term || p_Img->listX[LIST_0][i]->is_long_term)
            {
              currSlice->wbp_weight[0][i][j][comp] = 32;
              currSlice->wbp_weight[1][i][j][comp] = 32;
            }
            else
            {
              tb = iClip3(-128,127,p_Img->ThisPOC - p_Img->listX[LIST_0][i]->poc);

              tx = (16384 + iabs(td/2))/td;
              DistScaleFactor = iClip3(-1024, 1023, (tx*tb + 32 )>>6);

              currSlice->wbp_weight[1][i][j][comp] = DistScaleFactor >> 2;
              currSlice->wbp_weight[0][i][j][comp] = 64 - currSlice->wbp_weight[1][i][j][comp];
              if (currSlice->wbp_weight[1][i][j][comp] < -64 || currSlice->wbp_weight[1][i][j][comp] > 128)
              {
                currSlice->wbp_weight[0][i][j][comp] = 32;
                currSlice->wbp_weight[1][i][j][comp] = 32;
                currSlice->wp_offset[0][i][comp] = 0;
                currSlice->wp_offset[1][j][comp] = 0;
              }
            }
          }
        }
      }
   }
 }

  if (bframe && currSlice->MbaffFrameFlag)
  {
    for (i=0; i<2*max_l0_ref; ++i)
    {
      for (j=0; j<2*max_l1_ref; ++j)
      {
        for (comp = 0; comp<3; ++comp)
        {
          for (k=2; k<6; k+=2)
          {
            currSlice->wp_offset[k+0][i][comp] = currSlice->wp_offset[0][i>>1][comp];
            currSlice->wp_offset[k+1][j][comp] = currSlice->wp_offset[1][j>>1][comp];

            log_weight_denom = (comp == 0) ? currSlice->luma_log2_weight_denom : currSlice->chroma_log2_weight_denom;
            if (p_Img->active_pps->weighted_bipred_idc == 1)
            {
              currSlice->wbp_weight[k+0][i][j][comp] =  currSlice->wp_weight[0][i>>1][comp];
              currSlice->wbp_weight[k+1][i][j][comp] =  currSlice->wp_weight[1][j>>1][comp];
            }
            else if (p_Img->active_pps->weighted_bipred_idc == 2)
            {
              td = iClip3(-128,127,p_Img->listX[k+LIST_1][j]->poc - p_Img->listX[k+LIST_0][i]->poc);
              if (td == 0 || p_Img->listX[k+LIST_1][j]->is_long_term || p_Img->listX[k+LIST_0][i]->is_long_term)
              {
                currSlice->wbp_weight[k+0][i][j][comp] =   32;
                currSlice->wbp_weight[k+1][i][j][comp] =   32;
              }
              else
              {
                tb = iClip3(-128,127,((k==2)?p_Img->toppoc:p_Img->bottompoc) - p_Img->listX[k+LIST_0][i]->poc);

                tx = (16384 + iabs(td/2))/td;
                DistScaleFactor = iClip3(-1024, 1023, (tx*tb + 32 )>>6);

                currSlice->wbp_weight[k+1][i][j][comp] = DistScaleFactor >> 2;
                currSlice->wbp_weight[k+0][i][j][comp] = 64 - currSlice->wbp_weight[k+1][i][j][comp];
                if (currSlice->wbp_weight[k+1][i][j][comp] < -64 || currSlice->wbp_weight[k+1][i][j][comp] > 128)
                {
                  currSlice->wbp_weight[k+1][i][j][comp] = 32;
                  currSlice->wbp_weight[k+0][i][j][comp] = 32;
                  currSlice->wp_offset[k+0][i][comp] = 0;
                  currSlice->wp_offset[k+1][j][comp] = 0;
                }
              }
            }
          }
        }
      }
    }
  }
}

static void decode_slice(Slice *currSlice, int current_header)
{
  ImageParameters *p_Img = currSlice->p_Img;
  
  if (p_Img->active_pps->entropy_coding_mode_flag)
  {
    init_contexts  (currSlice);
    cabac_new_slice(currSlice);
  }

  if ( (p_Img->active_pps->weighted_bipred_idc > 0  && (currSlice->slice_type == B_SLICE)) || (p_Img->active_pps->weighted_pred_flag && currSlice->slice_type !=I_SLICE))
    fill_wp_params(currSlice);

  //printf("frame picture %d %d %d\n",p_Img->structure,p_Img->ThisPOC,currSlice->direct_spatial_mv_pred_flag);

  // decode main slice information
  if ((current_header == SOP || current_header == SOS) && currSlice->ei_flag == 0)
    decode_one_slice(currSlice);

  // setMB-Nr in case this slice was lost
  // if(currSlice->ei_flag)
  //   p_Img->current_mb_nr = currSlice->last_mb_nr + 1;
}


/*!
 ************************************************************************
 * \brief
 *    Error tracking: if current frame is lost or any reference frame of
 *                    current frame is lost, current frame is incorrect.
 ************************************************************************
 */
static void Error_tracking(ImageParameters *p_Img)
{
  Slice *currSlice = p_Img->currentSlice;
  int i;

  if(p_Img->redundant_pic_cnt == 0)
    {
      p_Img->Is_primary_correct = p_Img->Is_redundant_correct = 1;
    }

  if(p_Img->redundant_pic_cnt == 0 && p_Img->type != I_SLICE)
  {
    for(i=0;i<currSlice->num_ref_idx_l0_active;++i)
    {
      if(p_Img->ref_flag[i] == 0)  // any reference of primary slice is incorrect
      {
        p_Img->Is_primary_correct = 0; // primary slice is incorrect
      }
    }
  }
  else if(p_Img->redundant_pic_cnt != 0 && p_Img->type != I_SLICE)
  {
    if(p_Img->ref_flag[p_Img->redundant_slice_ref_idx] == 0)  // reference of redundant slice is incorrect
    {
      p_Img->Is_redundant_correct = 0;  // redundant slice is incorrect
    }
  }
}

/*!
 ***********************************************************************
 * \brief
 *    decodes one I- or P-frame
 *
 ***********************************************************************
 */

int decode_one_frame(ImageParameters *p_Img)
{

  InputParameters *p_Inp = p_Img->p_Inp;

  //SNRParameters   *snr   = p_Img->snr;
  int current_header;
  Slice *currSlice = p_Img->currentSlice;
  int i;

  currSlice->p_Img = p_Img;
  currSlice->p_Inp = p_Inp;
  p_Img->current_slice_nr = 0;
  p_Img->current_mb_nr = -4711;     // initialized to an impossible value for debugging -- correct value is taken from slice header
  currSlice->next_header = -8888; // initialized to an impossible value for debugging -- correct value is taken from slice header
  p_Img->num_dec_mb = 0;
  p_Img->newframe = 1;
  //currSlice->coeff[64]; // one more for EOB
  currSlice->coeff_ctr = -1;
  currSlice->pos       =  0;  


  while ((currSlice->next_header != EOS && currSlice->next_header != SOP))
  {
    current_header = read_new_slice(p_Img->currentSlice);

    // error tracking of primary and redundant slices.
    Error_tracking(p_Img);

    // If primary and redundant are received and primary is correct, discard the redundant
    // else, primary slice will be replaced with redundant slice.
    if(p_Img->frame_num == p_Img->previous_frame_num && p_Img->redundant_pic_cnt !=0
      && p_Img->Is_primary_correct !=0 && current_header != EOS)
    {
      continue;
    }

    // update reference flags and set current p_Img->ref_flag
    if(!(p_Img->redundant_pic_cnt != 0 && p_Img->previous_frame_num == p_Img->frame_num))
    {
      for(i=16;i>0;i--)
      {
        p_Img->ref_flag[i] = p_Img->ref_flag[i-1];
      }
    }
    p_Img->ref_flag[0] = p_Img->redundant_pic_cnt==0 ? p_Img->Is_primary_correct : p_Img->Is_redundant_correct;
    p_Img->previous_frame_num = p_Img->frame_num;

    if (current_header == EOS)
    {
      exit_picture(p_Img, &p_Img->dec_picture);
      return EOS;
    }

    if((p_Img->active_sps->chroma_format_idc==0)||(p_Img->active_sps->chroma_format_idc==3))
    {
      currSlice->linfo_cbp_intra = linfo_cbp_intra_other;
      currSlice->linfo_cbp_inter = linfo_cbp_inter_other;
    }
    else
    {
      currSlice->linfo_cbp_intra = linfo_cbp_intra_normal;
      currSlice->linfo_cbp_inter = linfo_cbp_inter_normal;
    }

    decode_slice(currSlice, current_header);

    p_Img->newframe = 0;
    ++(p_Img->current_slice_nr);
  }

  exit_picture(p_Img, &p_Img->dec_picture);

  return (SOP);
}


/*!
 ************************************************************************
 * \brief
 *    Convert file read buffer to source picture structure
 * \param imgX
 *    Pointer to image plane
 * \param buf
 *    Buffer for file output
 * \param size_x
 *    horizontal image size in pixel
 * \param size_y
 *    vertical image size in pixel
 * \param symbol_size_in_bytes
 *    number of bytes used per pel
 ************************************************************************
 */
void buffer2img (imgpel** imgX, unsigned char* buf, int size_x, int size_y, int symbol_size_in_bytes)
{
  int i,j;

  uint16 tmp16, ui16;
  unsigned long  tmp32, ui32;

  if (symbol_size_in_bytes> sizeof(imgpel))
  {
    error ("Source picture has higher bit depth than imgpel data type. \nPlease recompile with larger data type for imgpel.", 500);
  }

  if (( sizeof(char) == sizeof (imgpel)) && ( sizeof(char) == symbol_size_in_bytes))
  {
    // imgpel == pixel_in_file == 1 byte -> simple copy
    memcpy(&imgX[0][0], buf, size_x * size_y);
  }
  else
  {
    // sizeof (imgpel) > sizeof(char)
    if (testEndian())
    {
      // big endian
      switch (symbol_size_in_bytes)
      {
      case 1:
        {
          for(j = 0; j < size_y; ++j)
            for(i = 0; i < size_x; ++i)
            {
              imgX[j][i]= buf[i+j*size_x];
            }
          break;
        }
      case 2:
        {
          for(j=0;j<size_y;++j)
            for(i=0;i<size_x;++i)
            {
              memcpy(&tmp16, buf+((i+j*size_x)*2), 2);
              ui16  = (uint16) ((tmp16 >> 8) | ((tmp16&0xFF)<<8));
              imgX[j][i] = (imgpel) ui16;
            }
          break;
        }
      case 4:
        {
          for(j=0;j<size_y;++j)
            for(i=0;i<size_x;++i)
            {
              memcpy(&tmp32, buf+((i+j*size_x)*4), 4);
              ui32  = ((tmp32&0xFF00)<<8) | ((tmp32&0xFF)<<24) | ((tmp32&0xFF0000)>>8) | ((tmp32&0xFF000000)>>24);
              imgX[j][i] = (imgpel) ui32;
            }
        }
      default:
        {
           error ("reading only from formats of 8, 16 or 32 bit allowed on big endian architecture", 500);
           break;
        }
      }

    }
    else
    {
      // little endian
      if (symbol_size_in_bytes == 1)
      {
        for (j=0; j < size_y; ++j)
        {
          for (i=0; i < size_x; ++i)
          {
            imgX[j][i]=*(buf++);
          }
        }
      }
      else
      {
        for (j=0; j < size_y; ++j)
        {
          int jpos = j*size_x;
          for (i=0; i < size_x; ++i)
          {
            imgX[j][i]=0;
            memcpy(&(imgX[j][i]), buf +((i+jpos)*symbol_size_in_bytes), symbol_size_in_bytes);
          }
        }
      }

    }
  }
}


/*!
 ************************************************************************
 * \brief
 *    Calculate the value of frame_no
 ************************************************************************
*/
void calculate_frame_no(ImageParameters *p_Img, StorablePicture *p)
{
  InputParameters *p_Inp = p_Img->p_Inp;
  // calculate frame number
  int  psnrPOC = p_Img->active_sps->mb_adaptive_frame_field_flag ? p->poc /(p_Inp->poc_scale) : p->poc/(p_Inp->poc_scale);

  if (psnrPOC==0)// && p_Img->psnr_number)
    p_Img->idr_psnr_number = p_Img->number*p_Img->ref_poc_gap/(p_Inp->poc_scale);

  p_Img->psnr_number=imax(p_Img->psnr_number,p_Img->idr_psnr_number+psnrPOC);

  p_Img->frame_no = p_Img->idr_psnr_number + psnrPOC;
}


/*!
************************************************************************
* \brief
*    Find PSNR for all three components.Compare decoded frame with
*    the original sequence. Read p_Inp->jumpd frames to reflect frame skipping.
* \param p_Img
*      image encoding parameters for current picture
* \param p
*      picture to be compared
* \param p_ref
*      file pointer piont to reference YUV reference file
************************************************************************
*/
void find_snr(ImageParameters *p_Img, 
              StorablePicture *p,
              int *p_ref)
{
  InputParameters *p_Inp = p_Img->p_Inp;
  SNRParameters   *snr   = p_Img->snr;

  int i,j, k;
  int ret;
  int64 diff_comp[3] = {0};
  int64  status;
  int symbol_size_in_bytes = p_Img->pic_unit_bitsize_on_disk/8;
  int comp_size_x[3], comp_size_y[3];
  int64 framesize_in_bytes;

  unsigned int max_pix_value_sqd[3] = {iabs2(p_Img->max_imgpel_value_comp[0]),  iabs2(p_Img->max_imgpel_value_comp[1]), iabs2(p_Img->max_imgpel_value_comp[2])};

  Boolean rgb_output = (Boolean) (p_Img->active_sps->vui_seq_parameters.matrix_coefficients==0);
  unsigned char *buf;
  imgpel **cur_ref[3]  = {p_Img->imgY_ref, p->chroma_format_idc != YUV400 ? p_Img->imgUV_ref[0] : NULL, p->chroma_format_idc != YUV400 ? p_Img->imgUV_ref[1] : NULL};
  imgpel **cur_comp[3] = {p->imgY,  p->chroma_format_idc != YUV400 ? p->imgUV[0]  : NULL , p->chroma_format_idc!= YUV400 ? p->imgUV[1]  : NULL}; 
  // picture error concealment
  char yuv_types[4][6]= {"4:0:0","4:2:0","4:2:2","4:4:4"};
  
  comp_size_x[0] = p_Inp->source.width;
  comp_size_y[0] = p_Inp->source.height;
  comp_size_x[1] = comp_size_x[2] = p_Inp->source.width_cr;
  comp_size_y[1] = comp_size_y[2] = p_Inp->source.height_cr;

  framesize_in_bytes = (((int64) comp_size_x[0] * comp_size_y[0]) + ((int64) comp_size_x[1] * comp_size_y[1] ) * 2) * symbol_size_in_bytes;

  // KS: this buffer should actually be allocated only once, but this is still much faster than the previous version
  buf = malloc ( comp_size_x[0] * comp_size_y[0] * symbol_size_in_bytes );

  if (NULL == buf)
  {
    no_mem_exit("find_snr: buf");
  }

  status = lseek (*p_ref, framesize_in_bytes * p_Img->frame_no, SEEK_SET);
  if (status == -1)
  {
    fprintf(stderr, "Warning: Could not seek to frame number %d in reference file. Shown PSNR might be wrong.\n", p_Img->frame_no);
    free (buf);
    return;
  }

  if(rgb_output)
    lseek (*p_ref, framesize_in_bytes/3, SEEK_CUR);

  for (k = 0; k < ((p->chroma_format_idc != YUV400) ? 3 : 1); ++k)
  {

    if(rgb_output && k == 2)
      lseek (*p_ref, -framesize_in_bytes, SEEK_CUR);

    ret = read(*p_ref, buf, comp_size_x[k] * comp_size_y[k] * symbol_size_in_bytes);
    if (ret != comp_size_x[k] * comp_size_y[k] * symbol_size_in_bytes)
    {
      printf ("Warning: could not read from reconstructed file\n");
      memset (buf, 0, comp_size_x[k] * comp_size_y[k] * symbol_size_in_bytes);
      close(*p_ref);
      *p_ref = -1;
      break;
    }
    buffer2img(cur_ref[k], buf, comp_size_x[k], comp_size_y[k], symbol_size_in_bytes);

    for (j=0; j < comp_size_y[k]; ++j)
    {
      for (i=0; i < comp_size_x[k]; ++i)
      {
        diff_comp[k] += iabs2(cur_comp[k][j][i] - cur_ref[k][j][i]);
      }
    }

    // Collecting SNR statistics
    snr->snr[k] = psnr( max_pix_value_sqd[k], comp_size_x[k] * comp_size_y[k], (float) diff_comp[k]);   

    if (p_Img->number == 0) // first
    {
      snr->snra[k] = snr->snr[k];                                                        // keep luma snr for first frame
    }
    else
    {
      snr->snra[k] = (float)(snr->snra[k]*(snr->frame_ctr)+snr->snr[k])/(snr->frame_ctr + 1); // average snr chroma for all frames
    }
  }

   if(rgb_output)
     lseek (*p_ref, framesize_in_bytes*2/3, SEEK_CUR);

  free (buf);

  // picture error concealment
  if(p->concealed_pic)
  {
      fprintf(stdout,"%04d(P)  %8d %5d %5d %7.4f %7.4f %7.4f  %s %5d\n",
          p_Img->frame_no, p->frame_poc, p->pic_num, p->qp,
          snr->snr[0], snr->snr[1], snr->snr[2], yuv_types[p->chroma_format_idc], 0);
  }
}


void reorder_lists(Slice *currSlice)
{
  ImageParameters *p_Img = currSlice->p_Img;
  if ((currSlice->slice_type != I_SLICE)&&(currSlice->slice_type != SI_SLICE))
  {
    if (currSlice->ref_pic_list_reordering_flag_l0)
    {
      reorder_ref_pic_list(p_Img, p_Img->listX[0], &p_Img->listXsize[0],
                           currSlice->num_ref_idx_l0_active - 1,
                           currSlice->reordering_of_pic_nums_idc_l0,
                           currSlice->abs_diff_pic_num_minus1_l0,
                           currSlice->long_term_pic_idx_l0);
    }
    if (p_Img->no_reference_picture == p_Img->listX[0][currSlice->num_ref_idx_l0_active-1])
    {
      if (p_Img->non_conforming_stream)
        printf("RefPicList0[ num_ref_idx_l0_active_minus1 ] is equal to 'no reference picture'\n");
      else
        error("RefPicList0[ num_ref_idx_l0_active_minus1 ] is equal to 'no reference picture', invalid bitstream",500);
    }
    // that's a definition
    p_Img->listXsize[0] = currSlice->num_ref_idx_l0_active;
  }
  if (currSlice->slice_type == B_SLICE)
  {
    if (currSlice->ref_pic_list_reordering_flag_l1)
    {
      reorder_ref_pic_list(p_Img, p_Img->listX[1], &p_Img->listXsize[1],
                           currSlice->num_ref_idx_l1_active - 1,
                           currSlice->reordering_of_pic_nums_idc_l1,
                           currSlice->abs_diff_pic_num_minus1_l1,
                           currSlice->long_term_pic_idx_l1);
    }
    if (p_Img->no_reference_picture == p_Img->listX[1][currSlice->num_ref_idx_l1_active-1])
    {
      if (p_Img->non_conforming_stream)
        printf("RefPicList1[ num_ref_idx_l1_active_minus1 ] is equal to 'no reference picture'\n");
      else
        error("RefPicList1[ num_ref_idx_l1_active_minus1 ] is equal to 'no reference picture', invalid bitstream",500);
    }
    // that's a definition
    p_Img->listXsize[1] = currSlice->num_ref_idx_l1_active;
  }

  free_ref_pic_list_reordering_buffer(currSlice);
}


/*!
 ************************************************************************
 * \brief
 *    initialize ref_pic_num array
 ************************************************************************
 */
void set_ref_pic_num(ImageParameters *p_Img)
{
  int i,j;
  StorablePicture *dec_picture = p_Img->dec_picture;
  int slice_id=p_Img->current_slice_nr;

  for (i=0;i<p_Img->listXsize[LIST_0];++i)
  {
    dec_picture->ref_pic_num        [slice_id][LIST_0][i] = p_Img->listX[LIST_0][i]->poc * 2 + ((p_Img->listX[LIST_0][i]->structure==BOTTOM_FIELD)?1:0) ;
    dec_picture->frm_ref_pic_num    [slice_id][LIST_0][i] = p_Img->listX[LIST_0][i]->frame_poc * 2;
    dec_picture->top_ref_pic_num    [slice_id][LIST_0][i] = p_Img->listX[LIST_0][i]->top_poc * 2;
    dec_picture->bottom_ref_pic_num [slice_id][LIST_0][i] = p_Img->listX[LIST_0][i]->bottom_poc * 2 + 1;
    //printf("POCS %d %d %d %d ",p_Img->listX[LIST_0][i]->frame_poc,p_Img->listX[LIST_0][i]->bottom_poc,p_Img->listX[LIST_0][i]->top_poc,p_Img->listX[LIST_0][i]->poc);
    //printf("refid %d %d %d %d\n",(int) dec_picture->frm_ref_pic_num[LIST_0][i],(int) dec_picture->top_ref_pic_num[LIST_0][i],(int) dec_picture->bottom_ref_pic_num[LIST_0][i],(int) dec_picture->ref_pic_num[LIST_0][i]);
  }

  for (i=0;i<p_Img->listXsize[LIST_1];++i)
  {
    dec_picture->ref_pic_num        [slice_id][LIST_1][i] = p_Img->listX[LIST_1][i]->poc  *2 + ((p_Img->listX[LIST_1][i]->structure==BOTTOM_FIELD)?1:0);
    dec_picture->frm_ref_pic_num    [slice_id][LIST_1][i] = p_Img->listX[LIST_1][i]->frame_poc * 2;
    dec_picture->top_ref_pic_num    [slice_id][LIST_1][i] = p_Img->listX[LIST_1][i]->top_poc * 2;
    dec_picture->bottom_ref_pic_num [slice_id][LIST_1][i] = p_Img->listX[LIST_1][i]->bottom_poc * 2 + 1;
  }

  if (!p_Img->active_sps->frame_mbs_only_flag)
  {
    if (p_Img->structure==FRAME)
    {
      for (j=2;j<6;++j)
      {
        for (i=0;i<p_Img->listXsize[j];++i)
        {
          dec_picture->ref_pic_num        [slice_id][j][i] = p_Img->listX[j][i]->poc * 2 + ((p_Img->listX[j][i]->structure==BOTTOM_FIELD)?1:0);
          dec_picture->frm_ref_pic_num    [slice_id][j][i] = p_Img->listX[j][i]->frame_poc * 2 ;
          dec_picture->top_ref_pic_num    [slice_id][j][i] = p_Img->listX[j][i]->top_poc * 2 ;
          dec_picture->bottom_ref_pic_num [slice_id][j][i] = p_Img->listX[j][i]->bottom_poc * 2 + 1;
        }
      }
    }
  }
}


/*!
 ************************************************************************
 * \brief
 *    Reads new slice from bit_stream
 ************************************************************************
 */
int read_new_slice(Slice *currSlice)
{
  ImageParameters *p_Img = currSlice->p_Img;
  InputParameters *p_Inp = currSlice->p_Inp;
    
  NALU_t *nalu = AllocNALU(MAX_CODED_FRAME_SIZE);
  int current_header = 0;
  int BitsUsedByHeader;
  Bitstream *currStream;

  int slice_id_a, slice_id_b, slice_id_c;
  int redundant_pic_cnt_b, redundant_pic_cnt_c;  

  while (1)
  {
    if (0 == read_next_nalu(p_Img, nalu))
      return EOS;

process_nalu:

    switch (nalu->nal_unit_type)
    {
      case NALU_TYPE_SLICE:
      case NALU_TYPE_IDR:

        if (p_Img->recovery_point || nalu->nal_unit_type == NALU_TYPE_IDR)
        {
          if (p_Img->recovery_point_found == 0)
          {
            if (nalu->nal_unit_type != NALU_TYPE_IDR)
            {
              printf("Warning: Decoding does not start with an IDR picture.\n");
              p_Img->non_conforming_stream = 1;
            }
            else
              p_Img->non_conforming_stream = 0;
          }
          p_Img->recovery_point_found = 1;
        }

        if (p_Img->recovery_point_found == 0)
          break;

        p_Img->idr_flag = (nalu->nal_unit_type == NALU_TYPE_IDR);
        p_Img->nal_reference_idc = nalu->nal_reference_idc;
        currSlice->dp_mode = PAR_DP_1;
        currSlice->max_part_nr = 1;
        currSlice->ei_flag = 0;
        currStream = currSlice->partArr[0].bitstream;
        currStream->ei_flag = 0;
        currStream->frame_bitoffset = currStream->read_len = 0;
        memcpy (currStream->streamBuffer, &nalu->buf[1], nalu->len-1);
        currStream->code_len = currStream->bitstream_length = RBSPtoSODB(currStream->streamBuffer, nalu->len-1);

        // Some syntax of the Slice Header depends on the parameter set, which depends on
        // the parameter set ID of the SLice header.  Hence, read the pic_parameter_set_id
        // of the slice header first, then setup the active parameter sets, and then read
        // the rest of the slice header
        BitsUsedByHeader = FirstPartOfSliceHeader(currSlice);
        UseParameterSet (currSlice, currSlice->pic_parameter_set_id);
        BitsUsedByHeader+= RestOfSliceHeader (currSlice);

        fmo_init (p_Img);

        assign_quant_params (currSlice, p_Img->active_pps, p_Img->active_sps);

        // if primary slice is replaced with redundant slice, set the correct image type
        if(p_Img->redundant_pic_cnt && p_Img->Is_primary_correct==0 && p_Img->Is_redundant_correct)
        {
          p_Img->dec_picture->slice_type = p_Img->type;
        }

        if(is_new_picture(p_Img->dec_picture, currSlice, p_Img->old_slice))
        {
          init_picture(p_Img, currSlice, p_Inp);

          current_header = SOP;
          //check zero_byte if it is also the first NAL unit in the access unit
          CheckZeroByteVCL(p_Img, nalu);
        }
        else
          current_header = SOS;

        if (currSlice->slice_type == B_SLICE)
        {
          if( IS_INDEPENDENT(p_Img) )
          {
            int nplane = 0;
            for( nplane=0; nplane<MAX_PLANE; nplane++ )
            {
              if( NULL != currSlice->Co_located_JV[nplane] )
              {
                free_colocated(currSlice->Co_located_JV[nplane]);        
                currSlice->Co_located_JV[nplane] = NULL;
              }
              currSlice->Co_located_JV[nplane] = alloc_colocated (p_Img->width, p_Img->height, p_Img->active_sps->mb_adaptive_frame_field_flag);      
            }
          }
          else
          {
            if (NULL != currSlice->p_colocated)
            {
              free_colocated(currSlice->p_colocated);
              currSlice->p_colocated = NULL;
            }
            currSlice->p_colocated = alloc_colocated (p_Img->width, p_Img->height,p_Img->active_sps->mb_adaptive_frame_field_flag);
          }
        }

 
        init_lists(currSlice);
        reorder_lists (currSlice);

        if (p_Img->structure==FRAME)
        {
          init_mbaff_lists(p_Img);
        }

        // From here on, p_Img->active_sps, p_Img->active_pps and the slice header are valid
        if (currSlice->MbaffFrameFlag)
          p_Img->current_mb_nr = currSlice->start_mb_nr << 1;
        else
          p_Img->current_mb_nr = currSlice->start_mb_nr;

        if (p_Img->active_pps->entropy_coding_mode_flag)
        {
          int ByteStartPosition = currStream->frame_bitoffset/8;
          if (currStream->frame_bitoffset%8 != 0)
          {
            ++ByteStartPosition;
          }
          arideco_start_decoding (&currSlice->partArr[0].de_cabac, currStream->streamBuffer, ByteStartPosition, &currStream->read_len);
        }
        // printf ("read_new_slice: returning %s\n", current_header == SOP?"SOP":"SOS");
        FreeNALU(nalu);
        p_Img->recovery_point = 0;
        return current_header;
        break;
      case NALU_TYPE_DPA:
        // read DP_A
        currSlice->dpB_NotPresent =1; 
        currSlice->dpC_NotPresent =1; 

        p_Img->idr_flag          = FALSE;
        p_Img->nal_reference_idc = nalu->nal_reference_idc;
        currSlice->dp_mode     = PAR_DP_3;
        currSlice->max_part_nr = 3;
        currSlice->ei_flag     = 0;
        currStream             = currSlice->partArr[0].bitstream;
        currStream->ei_flag    = 0;
        currStream->frame_bitoffset = currStream->read_len = 0;
        memcpy (currStream->streamBuffer, &nalu->buf[1], nalu->len-1);
        currStream->code_len = currStream->bitstream_length = RBSPtoSODB(currStream->streamBuffer, nalu->len-1);

        BitsUsedByHeader     = FirstPartOfSliceHeader(currSlice);
        UseParameterSet (currSlice, currSlice->pic_parameter_set_id);
        BitsUsedByHeader    += RestOfSliceHeader (currSlice);

        fmo_init (p_Img);

        if(is_new_picture(p_Img->dec_picture, currSlice, p_Img->old_slice))
        {
          init_picture(p_Img, currSlice, p_Inp);
          current_header = SOP;
          CheckZeroByteVCL(p_Img, nalu);
        }
        else
          current_header = SOS;

        init_lists(currSlice);
        reorder_lists (currSlice);

        if (p_Img->structure==FRAME)
        {
          init_mbaff_lists(p_Img);
        }

        // From here on, p_Img->active_sps, p_Img->active_pps and the slice header are valid
        if (currSlice->MbaffFrameFlag)
          p_Img->current_mb_nr = currSlice->start_mb_nr << 1;
        else
          p_Img->current_mb_nr = currSlice->start_mb_nr;

        // Now I need to read the slice ID, which depends on the value of
        // redundant_pic_cnt_present_flag

        slice_id_a  = ue_v("NALU: DP_A slice_id", currStream);

        if (p_Img->active_pps->entropy_coding_mode_flag)
          error ("received data partition with CABAC, this is not allowed", 500);

        // continue with reading next DP
        if (0 == read_next_nalu(p_Img, nalu))
          return current_header;
        
        if ( NALU_TYPE_DPB == nalu->nal_unit_type)
        {
          // we got a DPB
          currStream             = currSlice->partArr[1].bitstream;
          currStream->ei_flag    = 0;
          currStream->frame_bitoffset = currStream->read_len = 0;

          memcpy (currStream->streamBuffer, &nalu->buf[1], nalu->len-1);
          currStream->code_len = currStream->bitstream_length = RBSPtoSODB(currStream->streamBuffer, nalu->len-1);

          slice_id_b  = ue_v("NALU: DP_B slice_id", currStream);

          currSlice->dpB_NotPresent = 0; 

          if ((slice_id_b != slice_id_a) || (nalu->lost_packets))
          {
            printf ("Waning: got a data partition B which does not match DP_A (DP loss!)\n");
            currSlice->dpB_NotPresent =1; 
            currSlice->dpC_NotPresent =1; 
          }
          else
          {
            if (p_Img->active_pps->redundant_pic_cnt_present_flag)
              redundant_pic_cnt_b = ue_v("NALU: DP_B redudant_pic_cnt", currStream);
            else
              redundant_pic_cnt_b = 0;

            // we're finished with DP_B, so let's continue with next DP
            if (0 == read_next_nalu(p_Img, nalu))
              return current_header;
          }
        }
        else
        {
          currSlice->dpB_NotPresent =1; 
        }

        // check if we got DP_C
        if ( NALU_TYPE_DPC == nalu->nal_unit_type)
        {
          currStream             = currSlice->partArr[2].bitstream;
          currStream->ei_flag    = 0;
          currStream->frame_bitoffset = currStream->read_len = 0;

          memcpy (currStream->streamBuffer, &nalu->buf[1], nalu->len-1);
          currStream->code_len = currStream->bitstream_length = RBSPtoSODB(currStream->streamBuffer, nalu->len-1);

          currSlice->dpC_NotPresent = 0;

          slice_id_c  = ue_v("NALU: DP_C slice_id", currStream);
          if ((slice_id_c != slice_id_a)|| (nalu->lost_packets))
          {
            printf ("Warning: got a data partition C which does not match DP_A(DP loss!)\n");
            //currSlice->dpB_NotPresent =1;
            currSlice->dpC_NotPresent =1;
          }

          if (p_Img->active_pps->redundant_pic_cnt_present_flag)
            redundant_pic_cnt_c = ue_v("NALU:SLICE_C redudand_pic_cnt", currStream);
          else
            redundant_pic_cnt_c = 0;
        }
        else
        {
          currSlice->dpC_NotPresent =1;
        }

        // check if we read anything else than the expected partitions
        if ((nalu->nal_unit_type != NALU_TYPE_DPB) && (nalu->nal_unit_type != NALU_TYPE_DPC))
        {
          // we have a NALI that we can't process here, so restart processing
          goto process_nalu;
          // yes, "goto" should not be used, but it's really the best way here before we restructure the decoding loop
          // (which should be taken care of anyway)
        }

        FreeNALU(nalu);
        return current_header;

        break;
      case NALU_TYPE_DPB:
        printf ("found data partition B without matching DP A, discarding\n");
        break;
      case NALU_TYPE_DPC:
        printf ("found data partition C without matching DP A, discarding\n");
        break;
      case NALU_TYPE_SEI:
        //printf ("read_new_slice: Found NALU_TYPE_SEI, len %d\n", nalu->len);
        InterpretSEIMessage(nalu->buf,nalu->len,p_Img);
        break;
      case NALU_TYPE_PPS:
        ProcessPPS(p_Img, nalu);
        break;
      case NALU_TYPE_SPS:
        ProcessSPS(p_Img, nalu);
        break;
      case NALU_TYPE_AUD:
//        printf ("read_new_slice: Found 'Access Unit Delimiter' NAL unit, len %d, ignored\n", nalu->len);
        break;
      case NALU_TYPE_EOSEQ:
//        printf ("read_new_slice: Found 'End of Sequence' NAL unit, len %d, ignored\n", nalu->len);
        break;
      case NALU_TYPE_EOSTREAM:
//        printf ("read_new_slice: Found 'End of Stream' NAL unit, len %d, ignored\n", nalu->len);
        break;
      case NALU_TYPE_FILL:
        printf ("read_new_slice: Found NALU_TYPE_FILL, len %d\n", (int) nalu->len);
        printf ("Skipping these filling bits, proceeding w/ next NALU\n");
        break;
      default:
        printf ("Found NALU type %d, len %d undefined, ignore NALU, moving on\n", (int) nalu->nal_unit_type, (int) nalu->len);
        break;
    }
  }
  FreeNALU(nalu);
        
  return  current_header;
}



/*!
 ************************************************************************
 * \brief
 *    finish decoding of a picture, conceal errors and store it
 *    into the DPB
 ************************************************************************
 */
void exit_picture(ImageParameters *p_Img, StorablePicture **dec_picture)
{
  InputParameters *p_Inp = p_Img->p_Inp;
  SNRParameters   *snr   = p_Img->snr;
  char yuv_types[4][6]= {"4:0:0","4:2:0","4:2:2","4:4:4"};
  int ercStartMB;
  int ercSegment;
  frame recfr;
  unsigned int i;
  int structure, frame_poc, slice_type, refpic, qp, pic_num, chroma_format_idc, is_idr;

  int64 tmp_time;                   // time used by decoding the last frame
  char   yuvFormat[10];
  int nplane;

  // return if the last picture has already been finished
  if (*dec_picture==NULL)
  {
    return;
  }

  recfr.p_Img = p_Img;
  recfr.yptr = &(*dec_picture)->imgY[0][0];
  if ((*dec_picture)->chroma_format_idc != YUV400)
  {
    recfr.uptr = &(*dec_picture)->imgUV[0][0][0];
    recfr.vptr = &(*dec_picture)->imgUV[1][0][0];
  }

  //! this is always true at the beginning of a picture
  ercStartMB = 0;
  ercSegment = 0;

  //! mark the start of the first segment
  if (!(*dec_picture)->MbaffFrameFlag)
  {
    ercStartSegment(0, ercSegment, 0 , p_Img->erc_errorVar);
    //! generate the segments according to the macroblock map
    for(i = 1; i<(*dec_picture)->PicSizeInMbs; ++i)
    {
      if(p_Img->mb_data[i].ei_flag != p_Img->mb_data[i-1].ei_flag)
      {
        ercStopSegment(i-1, ercSegment, 0, p_Img->erc_errorVar); //! stop current segment

        //! mark current segment as lost or OK
        if(p_Img->mb_data[i-1].ei_flag)
          ercMarkCurrSegmentLost((*dec_picture)->size_x, p_Img->erc_errorVar);
        else
          ercMarkCurrSegmentOK((*dec_picture)->size_x, p_Img->erc_errorVar);

        ++ercSegment;  //! next segment
        ercStartSegment(i, ercSegment, 0 , p_Img->erc_errorVar); //! start new segment
        ercStartMB = i;//! save start MB for this segment
      }
    }
    //! mark end of the last segment
    ercStopSegment((*dec_picture)->PicSizeInMbs-1, ercSegment, 0, p_Img->erc_errorVar);
    if(p_Img->mb_data[i-1].ei_flag)
      ercMarkCurrSegmentLost((*dec_picture)->size_x, p_Img->erc_errorVar);
    else
      ercMarkCurrSegmentOK((*dec_picture)->size_x, p_Img->erc_errorVar);

    //! call the right error concealment function depending on the frame type.
    p_Img->erc_mvperMB /= (*dec_picture)->PicSizeInMbs;

    p_Img->erc_img = p_Img;
    if((*dec_picture)->slice_type == I_SLICE || (*dec_picture)->slice_type == SI_SLICE) // I-frame
      ercConcealIntraFrame(p_Img, &recfr, (*dec_picture)->size_x, (*dec_picture)->size_y, p_Img->erc_errorVar);
    else
      ercConcealInterFrame(&recfr, p_Img->erc_object_list, (*dec_picture)->size_x, (*dec_picture)->size_y, p_Img->erc_errorVar, (*dec_picture)->chroma_format_idc);
  }

  //deblocking for frame or field
  if( IS_INDEPENDENT(p_Img) )
  {
    int colour_plane_id = p_Img->colour_plane_id;
    for( nplane=0; nplane<MAX_PLANE; ++nplane )
    {
      change_plane_JV( p_Img, nplane );
      DeblockPicture( p_Img, *dec_picture );
    }
    p_Img->colour_plane_id = colour_plane_id;
    make_frame_picture_JV(p_Img);
  }
  else
  {
    DeblockPicture( p_Img, *dec_picture );
  }

  if ((*dec_picture)->MbaffFrameFlag)
    MbAffPostProc(p_Img);

  if (p_Img->structure == FRAME)         // buffer mgt. for frame mode
    frame_postprocessing(p_Img);
  else
    field_postprocessing(p_Img);   // reset all interlaced variables

  structure  = (*dec_picture)->structure;
  slice_type = (*dec_picture)->slice_type;
  frame_poc  = (*dec_picture)->frame_poc;
  refpic     = (*dec_picture)->used_for_reference;
  qp         = (*dec_picture)->qp;
  pic_num    = (*dec_picture)->pic_num;
  is_idr     = (*dec_picture)->idr_flag;

  chroma_format_idc = (*dec_picture)->chroma_format_idc;

  store_picture_in_dpb(p_Img, *dec_picture);
  *dec_picture=NULL;

  if (p_Img->last_has_mmco_5)
  {
    p_Img->pre_frame_num = 0;
  }

  if (p_Inp->silent == FALSE)
  {
    if (structure==TOP_FIELD || structure==FRAME)
    {
      if(slice_type == I_SLICE && is_idr) // IDR picture
        strcpy(p_Img->cslice_type,"IDR");
      else if(slice_type == I_SLICE) // I picture
        strcpy(p_Img->cslice_type," I ");
      else if(slice_type == P_SLICE) // P pictures
        strcpy(p_Img->cslice_type," P ");
      else if(slice_type == SP_SLICE) // SP pictures
        strcpy(p_Img->cslice_type,"SP ");
      else if (slice_type == SI_SLICE)
        strcpy(p_Img->cslice_type,"SI ");
      else if(refpic) // stored B pictures
        strcpy(p_Img->cslice_type," B ");
      else // B pictures
        strcpy(p_Img->cslice_type," b ");
      if (structure==FRAME)
      {
        strncat(p_Img->cslice_type,")       ",8-strlen(p_Img->cslice_type));
      }
    }
    else if (structure==BOTTOM_FIELD)
    {
      if(slice_type == I_SLICE && is_idr) // IDR picture
        strncat(p_Img->cslice_type,"|IDR)",8-strlen(p_Img->cslice_type));
      else if(slice_type == I_SLICE) // I picture
        strncat(p_Img->cslice_type,"| I )",8-strlen(p_Img->cslice_type));
      else if(slice_type == P_SLICE) // P pictures
        strncat(p_Img->cslice_type,"| P )",8-strlen(p_Img->cslice_type));
      else if(slice_type == SP_SLICE) // SP pictures
        strncat(p_Img->cslice_type,"|SP )",8-strlen(p_Img->cslice_type));
      else if (slice_type == SI_SLICE)
        strncat(p_Img->cslice_type,"|SI )",8-strlen(p_Img->cslice_type));
      else if(refpic) // stored B pictures
        strncat(p_Img->cslice_type,"| B )",8-strlen(p_Img->cslice_type));
      else // B pictures
        strncat(p_Img->cslice_type,"| b )",8-strlen(p_Img->cslice_type));   
    }
  }

  if ((structure==FRAME)||structure==BOTTOM_FIELD)
  {
    gettime (&(p_Img->end_time));              // end time

    tmp_time  = timediff(&(p_Img->start_time), &(p_Img->end_time));
    p_Img->tot_time += tmp_time;
    tmp_time  = timenorm(tmp_time);

    sprintf(yuvFormat,"%s", yuv_types[chroma_format_idc]);
    
    if (p_Inp->silent == FALSE)
    {
      SNRParameters   *snr = p_Img->snr;
      if (p_Img->p_ref != -1)
        fprintf(stdout,"%05d(%s%5d %5d %5d %8.4f %8.4f %8.4f  %s %7d\n",
        p_Img->frame_no, p_Img->cslice_type, frame_poc, pic_num, qp, snr->snr[0], snr->snr[1], snr->snr[2], yuvFormat, (int) tmp_time);
    else
        fprintf(stdout,"%05d(%s%5d %5d %5d                             %s %7d\n",
        p_Img->frame_no, p_Img->cslice_type, frame_poc, pic_num, qp, yuvFormat, (int)tmp_time);
    }
    else
      fprintf(stdout,"Completed Decoding frame %05d.\r",snr->frame_ctr);

    fflush(stdout);

    if(slice_type == I_SLICE || slice_type == SI_SLICE || slice_type == P_SLICE || refpic)   // I or P pictures
      ++(p_Img->number);
    else
      ++(p_Img->Bframe_ctr);    // B pictures
    ++(snr->frame_ctr);

    ++(p_Img->g_nFrame);
  }

  p_Img->current_mb_nr = -4712;   // impossible value for debugging, StW
  p_Img->current_slice_nr = 0;
}

/*!
 ************************************************************************
 * \brief
 *    write the encoding mode and motion vectors of current
 *    MB to the buffer of the error concealment module.
 ************************************************************************
 */

void ercWriteMBMODEandMV(Macroblock *currMB)
{
  ImageParameters *p_Img = currMB->p_Img;
  int i, ii, jj, currMBNum = p_Img->current_mb_nr;
  StorablePicture *dec_picture = p_Img->dec_picture;
  int mbx = xPosMB(currMBNum, dec_picture->size_x), mby = yPosMB(currMBNum, dec_picture->size_x);
  objectBuffer_t *currRegion, *pRegion;
  short***  mv;

  currRegion = p_Img->erc_object_list + (currMBNum<<2);

  if(p_Img->type != B_SLICE) //non-B frame
  {
    for (i=0; i<4; ++i)
    {
      pRegion             = currRegion + i;
      pRegion->regionMode = (currMB->mb_type  ==I16MB  ? REGMODE_INTRA      :
                             currMB->b8mode[i]==IBLOCK ? REGMODE_INTRA_8x8  :
                             currMB->b8mode[i]==0      ? REGMODE_INTER_COPY :
                             currMB->b8mode[i]==1      ? REGMODE_INTER_PRED : REGMODE_INTER_PRED_8x8);
      if (currMB->b8mode[i]==0 || currMB->b8mode[i]==IBLOCK)  // INTRA OR COPY
      {
        pRegion->mv[0]    = 0;
        pRegion->mv[1]    = 0;
        pRegion->mv[2]    = 0;
      }
      else
      {
        ii              = 4*mbx + (i & 0x01)*2;// + BLOCK_SIZE;
        jj              = 4*mby + (i >> 1  )*2;
        if (currMB->b8mode[i]>=5 && currMB->b8mode[i]<=7)  // SMALL BLOCKS
        {
          pRegion->mv[0]  = (dec_picture->motion.mv[LIST_0][jj][ii][0] + dec_picture->motion.mv[LIST_0][jj][ii+1][0] + dec_picture->motion.mv[LIST_0][jj+1][ii][0] + dec_picture->motion.mv[LIST_0][jj+1][ii+1][0] + 2)/4;
          pRegion->mv[1]  = (dec_picture->motion.mv[LIST_0][jj][ii][1] + dec_picture->motion.mv[LIST_0][jj][ii+1][1] + dec_picture->motion.mv[LIST_0][jj+1][ii][1] + dec_picture->motion.mv[LIST_0][jj+1][ii+1][1] + 2)/4;
        }
        else // 16x16, 16x8, 8x16, 8x8
        {
          pRegion->mv[0]  = dec_picture->motion.mv[LIST_0][jj][ii][0];
          pRegion->mv[1]  = dec_picture->motion.mv[LIST_0][jj][ii][1];
//          pRegion->mv[0]  = dec_picture->motion.mv[LIST_0][4*mby+(i/2)*2][4*mbx+(i%2)*2+BLOCK_SIZE][0];
//          pRegion->mv[1]  = dec_picture->motion.mv[LIST_0][4*mby+(i/2)*2][4*mbx+(i%2)*2+BLOCK_SIZE][1];
        }
        p_Img->erc_mvperMB      += iabs(pRegion->mv[0]) + iabs(pRegion->mv[1]);
        pRegion->mv[2]    = dec_picture->motion.ref_idx[LIST_0][jj][ii];
      }
    }
  }
  else  //B-frame
  {
    for (i=0; i<4; ++i)
    {
      ii                  = 4*mbx + (i%2)*2;// + BLOCK_SIZE;
      jj                  = 4*mby + (i/2)*2;
      pRegion             = currRegion + i;
      pRegion->regionMode = (currMB->mb_type  ==I16MB  ? REGMODE_INTRA      :
                             currMB->b8mode[i]==IBLOCK ? REGMODE_INTRA_8x8  : REGMODE_INTER_PRED_8x8);
      if (currMB->mb_type==I16MB || currMB->b8mode[i]==IBLOCK)  // INTRA
      {
        pRegion->mv[0]    = 0;
        pRegion->mv[1]    = 0;
        pRegion->mv[2]    = 0;
      }
      else
      {
        int idx = (dec_picture->motion.ref_idx[0][jj][ii]<0)?1:0;
//        int idx = (currMB->b8mode[i]==0 && currMB->b8pdir[i]==2 ? LIST_0 : currMB->b8pdir[i]==1 ? LIST_1 : LIST_0);
//        int idx = currMB->b8pdir[i]==0 ? LIST_0 : LIST_1;
        mv                = dec_picture->motion.mv[idx];
        pRegion->mv[0]    = (mv[jj][ii][0] + mv[jj][ii+1][0] + mv[jj+1][ii][0] + mv[jj+1][ii+1][0] + 2)/4;
        pRegion->mv[1]    = (mv[jj][ii][1] + mv[jj][ii+1][1] + mv[jj+1][ii][1] + mv[jj+1][ii+1][1] + 2)/4;
        p_Img->erc_mvperMB      += iabs(pRegion->mv[0]) + iabs(pRegion->mv[1]);

        pRegion->mv[2]  = (dec_picture->motion.ref_idx[idx][jj][ii]);
/*
        if (currMB->b8pdir[i]==0 || (currMB->b8pdir[i]==2 && currMB->b8mode[i]!=0)) // forward or bidirect
        {
          pRegion->mv[2]  = (dec_picture->motion.ref_idx[LIST_0][jj][ii]);
          ///???? is it right, not only "p_Img->fw_refFrArr[jj][ii-4]"
        }
        else
        {
          pRegion->mv[2]  = (dec_picture->motion.ref_idx[LIST_1][jj][ii]);
//          pRegion->mv[2]  = 0;
        }
        */
      }
    }
  }
}

/*!
 ************************************************************************
 * \brief
 *    set defaults for old_slice
 *    NAL unit of a picture"
 ************************************************************************
 */
void init_old_slice(OldSliceParams *p_old_slice)
{
  p_old_slice->field_pic_flag = 0;

  p_old_slice->pps_id = INT_MAX;

  p_old_slice->frame_num = INT_MAX;

  p_old_slice->nal_ref_idc = INT_MAX;

  p_old_slice->idr_flag = FALSE;

  p_old_slice->pic_oder_cnt_lsb          = UINT_MAX;
  p_old_slice->delta_pic_oder_cnt_bottom = INT_MAX;

  p_old_slice->delta_pic_order_cnt[0] = INT_MAX;
  p_old_slice->delta_pic_order_cnt[1] = INT_MAX;
}

/*!
 ************************************************************************
 * \brief
 *    save slice parameters that are needed for checking of "first VCL
 *    NAL unit of a picture"
 ************************************************************************
 */
void exit_slice(Slice *currSlice, OldSliceParams *p_old_slice)
{
  ImageParameters *p_Img = currSlice->p_Img;

  if (currSlice->slice_type == B_SLICE)
  {
    if( IS_INDEPENDENT(p_Dec->p_Img) )
    {
      int nplane;
      for( nplane = 0; nplane < MAX_PLANE; ++nplane )
      {
        free_colocated(currSlice->Co_located_JV[nplane]);   
        currSlice->Co_located_JV[nplane] = NULL;
      }
    }
    else
    {
      free_colocated(currSlice->p_colocated);
      currSlice->p_colocated = NULL;
    }
  }

  p_old_slice->pps_id = currSlice->pic_parameter_set_id;

  p_old_slice->frame_num = p_Img->frame_num;

  p_old_slice->field_pic_flag = p_Img->field_pic_flag;

  if(p_Img->field_pic_flag)
  {
    p_old_slice->bottom_field_flag = p_Img->bottom_field_flag;
  }

  p_old_slice->nal_ref_idc   = p_Img->nal_reference_idc;

  p_old_slice->idr_flag = p_Img->idr_flag;
  if (p_Img->idr_flag)
  {
    p_old_slice->idr_pic_id = p_Img->idr_pic_id;
  }

  if (p_Img->active_sps->pic_order_cnt_type == 0)
  {
    p_old_slice->pic_oder_cnt_lsb          = p_Img->pic_order_cnt_lsb;
    p_old_slice->delta_pic_oder_cnt_bottom = p_Img->delta_pic_order_cnt_bottom;
  }

  if (p_Img->active_sps->pic_order_cnt_type == 1)
  {
    p_old_slice->delta_pic_order_cnt[0] = p_Img->delta_pic_order_cnt[0];
    p_old_slice->delta_pic_order_cnt[1] = p_Img->delta_pic_order_cnt[1];
  }
}

/*!
 ************************************************************************
 * \brief
 *    detect if current slice is "first VCL NAL unit of a picture"
 ************************************************************************
 */
int is_new_picture(StorablePicture *dec_picture, Slice *currSlice, OldSliceParams *p_old_slice)
{
  ImageParameters *p_Img = currSlice->p_Img;

  int result=0;

  result |= (NULL==dec_picture);

  result |= (p_old_slice->pps_id != currSlice->pic_parameter_set_id);

  result |= (p_old_slice->frame_num != p_Img->frame_num);

  result |= (p_old_slice->field_pic_flag != p_Img->field_pic_flag);

  if(p_Img->field_pic_flag && p_old_slice->field_pic_flag)
  {
    result |= (p_old_slice->bottom_field_flag != p_Img->bottom_field_flag);
  }

  result |= (p_old_slice->nal_ref_idc != p_Img->nal_reference_idc) && ((p_old_slice->nal_ref_idc == 0) || (p_Img->nal_reference_idc == 0));

  result |= ( p_old_slice->idr_flag != p_Img->idr_flag);

  if (p_Img->idr_flag && p_old_slice->idr_flag)
  {
    result |= (p_old_slice->idr_pic_id != p_Img->idr_pic_id);
  }

  if (p_Img->active_sps->pic_order_cnt_type == 0)
  {
    result |= (p_old_slice->pic_oder_cnt_lsb          != p_Img->pic_order_cnt_lsb);
    result |= (p_old_slice->delta_pic_oder_cnt_bottom != p_Img->delta_pic_order_cnt_bottom);
  }

  if (p_Img->active_sps->pic_order_cnt_type == 1)
  {
    result |= (p_old_slice->delta_pic_order_cnt[0] != p_Img->delta_pic_order_cnt[0]);
    result |= (p_old_slice->delta_pic_order_cnt[1] != p_Img->delta_pic_order_cnt[1]);
  }

  return result;
}



/*!
 ************************************************************************
 * \brief
 *    Prepare field and frame buffer after frame decoding
 ************************************************************************
 */
void frame_postprocessing(ImageParameters *p_Img)
{
}

/*!
 ************************************************************************
 * \brief
 *    Prepare field and frame buffer after field decoding
 ************************************************************************
 */
void field_postprocessing(ImageParameters *p_Img)
{
  p_Img->number /= 2;
}



/*!
 ************************************************************************
 * \brief
 *    copy StorablePicture *src -> StorablePicture *dst
 *    for 4:4:4 Independent mode
 ************************************************************************
 */
void copy_dec_picture_JV( ImageParameters *p_Img, StorablePicture *dst, StorablePicture *src )
{
  dst->top_poc              = src->top_poc;
  dst->bottom_poc           = src->bottom_poc;
  dst->frame_poc            = src->frame_poc;
  dst->qp                   = src->qp;
  dst->slice_qp_delta       = src->slice_qp_delta;
  dst->chroma_qp_offset[0]  = src->chroma_qp_offset[0];
  dst->chroma_qp_offset[1]  = src->chroma_qp_offset[1];

  dst->poc                  = src->poc;

  dst->slice_type           = src->slice_type;
  dst->used_for_reference   = src->used_for_reference;
  dst->idr_flag             = src->idr_flag;
  dst->no_output_of_prior_pics_flag = src->no_output_of_prior_pics_flag;
  dst->long_term_reference_flag = src->long_term_reference_flag;
  dst->adaptive_ref_pic_buffering_flag = src->adaptive_ref_pic_buffering_flag;

  dst->dec_ref_pic_marking_buffer = src->dec_ref_pic_marking_buffer;

  dst->MbaffFrameFlag       = src->MbaffFrameFlag;
  dst->PicWidthInMbs        = src->PicWidthInMbs;
  dst->pic_num              = src->pic_num;
  dst->frame_num            = src->frame_num;
  dst->recovery_frame       = src->recovery_frame;
  dst->coded_frame          = src->coded_frame;

  dst->chroma_format_idc    = src->chroma_format_idc;

  dst->frame_mbs_only_flag  = src->frame_mbs_only_flag;
  dst->frame_cropping_flag  = src->frame_cropping_flag;

  dst->frame_cropping_rect_left_offset   = src->frame_cropping_rect_left_offset;
  dst->frame_cropping_rect_right_offset  = src->frame_cropping_rect_right_offset;
  dst->frame_cropping_rect_top_offset    = src->frame_cropping_rect_top_offset;
  dst->frame_cropping_rect_bottom_offset = src->frame_cropping_rect_bottom_offset;

#if (ENABLE_OUTPUT_TONEMAPPING)
  // store the necessary tone mapping sei into StorablePicture structure
  dst->seiHasTone_mapping = src->seiHasTone_mapping;

  dst->seiHasTone_mapping    = src->seiHasTone_mapping;
  dst->tone_mapping_model_id = src->tone_mapping_model_id;
  dst->tonemapped_bit_depth  = src->tonemapped_bit_depth;
  if( src->tone_mapping_lut )
  {
    dst->tone_mapping_lut      = malloc(sizeof(int)*(1 << p_Img->seiToneMapping->coded_data_bit_depth));
    if (NULL == dst->tone_mapping_lut)
    {
      no_mem_exit("copy_dec_picture_JV: tone_mapping_lut");
    }
    memcpy(dst->tone_mapping_lut, src->tone_mapping_lut, sizeof(imgpel)*(1<<p_Img->seiToneMapping->coded_data_bit_depth));
  }
#endif
}


/*!
 ************************************************************************
 * \brief
 *    decodes one slice
 ************************************************************************
 */
void decode_one_slice(Slice *currSlice)
{
  ImageParameters *p_Img = currSlice->p_Img;
  Boolean end_of_slice = FALSE;
  Macroblock *currMB = NULL;
  p_Img->cod_counter=-1;

  setup_slice_methods(currSlice);

  if( IS_INDEPENDENT(p_Img) )
  {
    change_plane_JV( p_Img, p_Img->colour_plane_id );
  }

  set_ref_pic_num(p_Img);

  if (currSlice->slice_type == B_SLICE)
  {
    if( IS_INDEPENDENT(p_Img) )
    {
      compute_colocated_JV(currSlice, currSlice->p_colocated, p_Img->listX);
    }
    else
    {
      compute_colocated(currSlice, currSlice->p_colocated, p_Img->listX);
    }
  }

  //reset_ec_flags(p_Img);

  while (end_of_slice == FALSE) // loop over macroblocks
  {

#if TRACE
  fprintf(p_Dec->p_trace,"\n*********** POC: %i (I/P) MB: %i Slice: %i Type %d **********\n", p_Img->ThisPOC, p_Img->current_mb_nr, p_Img->current_slice_nr, currSlice->slice_type);
#endif

    // Initializes the current macroblock
    start_macroblock(currSlice, &currMB);
    // Get the syntax elements from the NAL
    currSlice->read_one_macroblock(currMB);
    decode_one_macroblock(currMB, p_Img->dec_picture);

    if(currSlice->MbaffFrameFlag && p_Img->dec_picture->motion.mb_field[p_Img->current_mb_nr])
    {
      currSlice->num_ref_idx_l0_active >>= 1;
      currSlice->num_ref_idx_l1_active >>= 1;
    }

    ercWriteMBMODEandMV(currMB);

    end_of_slice = exit_macroblock(currSlice, (!currSlice->MbaffFrameFlag||p_Img->current_mb_nr%2));
  }

  exit_slice(currSlice, p_Img->old_slice);
  //reset_ec_flags(p_Img);
}
