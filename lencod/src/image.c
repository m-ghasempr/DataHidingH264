
/*!
 *************************************************************************************
 * \file image.c
 *
 * \brief
 *    Code one image/slice
 *
 * \author
 *    Main contributors (see contributors.h for copyright, address and affiliation details)
 *     - Inge Lille-Langoy               <inge.lille-langoy@telenor.com>
 *     - Rickard Sjoberg                 <rickard.sjoberg@era.ericsson.se>
 *     - Jani Lainema                    <jani.lainema@nokia.com>
 *     - Sebastian Purreiter             <sebastian.purreiter@mch.siemens.de>
 *     - Byeong-Moon Jeon                <jeonbm@lge.com>
 *     - Yoon-Seong Soh                  <yunsung@lge.com>
 *     - Thomas Stockhammer              <stockhammer@ei.tum.de>
 *     - Detlev Marpe                    <marpe@hhi.de>
 *     - Guido Heising
 *     - Thomas Wedi                     <wedi@tnt.uni-hannover.de>
 *     - Ragip Kurceren                  <ragip.kurceren@nokia.com>
 *     - Antti Hallapuro                 <antti.hallapuro@nokia.com>
 *     - Alexis Michael Tourapis         <alexismt@ieee.org>
 *     - Athanasios Leontaris            <aleon@dolby.com>
 *************************************************************************************
 */
#include "contributors.h"

#include <math.h>
#include <time.h>

#include "global.h"

#include "filehandle.h"
#include "mbuffer.h"
#include "img_luma.h"
#include "img_chroma.h"
#include "img_distortion.h"
#include "intrarefresh.h"
#include "slice.h"
#include "fmo.h"
#include "sei.h"
#include "memalloc.h"
#include "nalu.h"
#include "ratectl.h"
#include "mb_access.h"
#include "context_ini.h"
#include "biariencode.h"
#include "enc_statistics.h"
#include "conformance.h"
#include "report.h"

#include "q_matrix.h"
#include "q_offsets.h"
#include "wp.h"
#include "input.h"
#include "image.h"
#include "errdo.h"
#include "img_process.h"
#include "rdopt.h"
#include "sei.h"

extern void DeblockFrame              (ImageParameters *p_Img, imgpel **, imgpel ***);

static void code_a_picture            (ImageParameters *p_Img, InputParameters *p_Inp, Picture *pic);
static void field_picture             (ImageParameters *p_Img, InputParameters *p_Inp, Picture *top, Picture *bottom);
static void prepare_enc_frame_picture (ImageParameters *p_Img, InputParameters *p_Inp, StorablePicture **stored_pic);
static void writeout_picture          (ImageParameters *p_Img, Picture *pic);
static byte picture_structure_decision(ImageParameters *p_Img, Picture *frame, Picture *top, Picture *bot);
static void distortion_fld            (ImageParameters *p_Img, InputParameters *p_Inp, Picture *field_pic, ImageData *imgData);


static void field_mode_buffer (ImageParameters *p_Img);
static void frame_mode_buffer (ImageParameters *p_Img, InputParameters *p_Inp);
static void init_frame        (ImageParameters *p_Img, InputParameters *p_Inp);
static void init_field        (ImageParameters *p_Img, InputParameters *p_Inp);
static void put_buffer_frame  (ImageParameters *p_Img);
static void put_buffer_top    (ImageParameters *p_Img);
static void put_buffer_bot    (ImageParameters *p_Img);

static void ReportFirstframe   (ImageParameters *p_Img, InputParameters *p_Inp, StatParameters *p_Stats, int64 tmp_time);
static void ReportI            (ImageParameters *p_Img, InputParameters *p_Inp, StatParameters *p_Stats, int64 tmp_time);
static void ReportP            (ImageParameters *p_Img, InputParameters *p_Inp, StatParameters *p_Stats, int64 tmp_time);
static void ReportB            (ImageParameters *p_Img, InputParameters *p_Inp, StatParameters *p_Stats, int64 tmp_time);
static void ReportNALNonVLCBits(ImageParameters *p_Img, InputParameters *p_Inp, StatParameters *p_Stats, int64 tmp_time);

extern void rd_picture_coding(ImageParameters *p_Img, InputParameters *p_Inp);

void MbAffPostProc(ImageParameters *p_Img)
{
  imgpel temp[32][16];

  StorablePicture *p_Enc_Pic = p_Img->enc_picture;
  imgpel ** imgY  = p_Enc_Pic->imgY;
  imgpel ***imgUV = p_Enc_Pic->imgUV;
  short i, y, x0, y0, uv;

  if (p_Img->yuv_format != YUV400)
  {
    for (i=0; i<(int)p_Img->PicSizeInMbs; i+=2)
    {
      if (p_Enc_Pic->motion.mb_field[i])
      {
        get_mb_pos(p_Img, i, p_Img->mb_size[IS_LUMA], &x0, &y0);
        for (y=0; y<(2*MB_BLOCK_SIZE);y++)
          memcpy(&temp[y],&imgY[y0+y][x0], MB_BLOCK_SIZE * sizeof(imgpel));

        for (y=0; y<MB_BLOCK_SIZE;y++)
        {
          memcpy(&imgY[y0+(2*y)][x0],temp[y], MB_BLOCK_SIZE * sizeof(imgpel));
          memcpy(&imgY[y0+(2*y + 1)][x0],temp[y+ MB_BLOCK_SIZE], MB_BLOCK_SIZE * sizeof(imgpel));
        }

        x0 = x0 / (16/p_Img->mb_cr_size_x);
        y0 = y0 / (16/p_Img->mb_cr_size_y);

        for (uv=0; uv<2; uv++)
        {
          for (y=0; y < (2 * p_Img->mb_cr_size_y); y++)
            memcpy(&temp[y],&imgUV[uv][y0+y][x0], p_Img->mb_cr_size_x * sizeof(imgpel));

          for (y=0; y<p_Img->mb_cr_size_y;y++)
          {
            memcpy(&imgUV[uv][y0+(2*y)][x0],temp[y], p_Img->mb_cr_size_x * sizeof(imgpel));
            memcpy(&imgUV[uv][y0+(2*y + 1)][x0],temp[y+ p_Img->mb_cr_size_y], p_Img->mb_cr_size_x * sizeof(imgpel));
          }
        }
      }
    }
  }
  else
  {
    for (i=0; i<(int)p_Img->PicSizeInMbs; i+=2)
    {
      if (p_Enc_Pic->motion.mb_field[i])
      {
        get_mb_pos(p_Img, i, p_Img->mb_size[IS_LUMA], &x0, &y0);
        for (y=0; y<(2*MB_BLOCK_SIZE);y++)
          memcpy(&temp[y],&imgY[y0+y][x0], MB_BLOCK_SIZE * sizeof(imgpel));

        for (y=0; y<MB_BLOCK_SIZE;y++)
        {
          memcpy(&imgY[y0+(2*y)][x0],temp[y], MB_BLOCK_SIZE * sizeof(imgpel));
          memcpy(&imgY[y0+(2*y + 1)][x0],temp[y+ MB_BLOCK_SIZE], MB_BLOCK_SIZE * sizeof(imgpel));
        }
      }
    }
  }
}

/*!
 ************************************************************************
 * \brief
 *    Sets slice type
 *
 ************************************************************************
 */
void set_slice_type(ImageParameters *p_Img, InputParameters *p_Inp, int slice_type)
{
  p_Img->type    = (short) slice_type;            // set slice type
  p_Img->RCMinQP = p_Inp->RCMinQP[p_Img->type];
  p_Img->RCMaxQP = p_Inp->RCMaxQP[p_Img->type];
}

static void code_a_plane(ImageParameters *p_Img, InputParameters *p_Inp)
{
  unsigned int NumberOfCodedMBs = 0;
  int SliceGroup = 0;
  // The slice_group_change_cycle can be changed here.
  // FmoInit() is called before coding each picture, frame or field
  p_Img->slice_group_change_cycle=1;
  FmoInit(p_Img, p_Img->active_pps, p_Img->active_sps);
  FmoStartPicture (p_Img);           //! picture level initialization of FMO

  CalculateQuant4x4Param (p_Img);
  CalculateOffset4x4Param(p_Img, p_Inp);

  if(p_Inp->Transform8x8Mode)
  {
    CalculateQuant8x8Param (p_Img);
    CalculateOffset8x8Param(p_Img, p_Inp);
  }

  reset_pic_bin_count(p_Img);
  p_Img->bytes_in_picture = 0;

  while (NumberOfCodedMBs < p_Img->PicSizeInMbs)       // loop over slices
  {
    // Encode one SLice Group
    while (!FmoSliceGroupCompletelyCoded (p_Img, SliceGroup))
    {
      // Encode the current slice
      if (!p_Img->MbaffFrameFlag)
        NumberOfCodedMBs += encode_one_slice (p_Img, p_Inp, SliceGroup, NumberOfCodedMBs);
      else
        NumberOfCodedMBs += encode_one_slice_MBAFF (p_Img, p_Inp, SliceGroup, NumberOfCodedMBs);

      FmoSetLastMacroblockInSlice (p_Img, p_Img->current_mb_nr);
      // Proceed to next slice
      p_Img->current_slice_nr++;
      p_Img->p_Stats->bit_slice = 0;
    }
    // Proceed to next SliceGroup
    SliceGroup++;
  }
  FmoEndPicture ();

  if ((p_Inp->SkipDeBlockNonRef == 0) || (p_Img->nal_reference_idc != 0))
    DeblockFrame (p_Img, p_Img->enc_picture->imgY, p_Img->enc_picture->imgUV); //comment out to disable deblocking filter 
}
/*!
 ************************************************************************
 * \brief
 *    Encodes a picture
 *
 *    This is the main picture coding loop.. It is called by all this
 *    frame and field coding stuff after the p_Img-> elements have been
 *    set up.  Not sure whether it is useful for MB-adaptive frame/field
 *    coding
 ************************************************************************
 */
static void code_a_picture(ImageParameters *p_Img, InputParameters *p_Inp, Picture *pic)
{
  int pl;

  p_Img->currentPicture = pic;
  p_Img->currentPicture->idr_flag = get_idr_flag(p_Img, p_Inp);

  if (p_Img->currentPicture->idr_flag && p_Inp->EnableIDRGOP && p_Img->frm_number)
  {
    p_Img->last_idr_number = p_Img->frm_no_in_file;
  }

  pic->no_slices = 0;

  RandomIntraNewPicture (p_Img);     //! Allocates forced INTRA MBs (even for fields!)

  if( IS_INDEPENDENT( p_Inp ) )
  {
    for( pl=0; pl<MAX_PLANE; pl++ )
    {
      p_Img->current_mb_nr = 0;
      p_Img->current_slice_nr = 0;
      p_Img->SumFrameQP = 0;
      p_Img->num_ref_idx_l0_active = 0;
      p_Img->num_ref_idx_l1_active = 0;

      p_Img->colour_plane_id = (char) pl;
      
      code_a_plane(p_Img, p_Inp);
    }
  }
  else
  {
    code_a_plane(p_Img, p_Inp);
  }

  if (p_Img->MbaffFrameFlag)
    MbAffPostProc(p_Img);
}

/*!
 ************************************************************************
 * \brief
 *    Determine whether picture is coded as IDR
 ************************************************************************
 */
byte get_idr_flag( ImageParameters *p_Img, InputParameters *p_Inp)
{
  int idr_flag;
  int idr_refresh;

  // currently this code only supports fixed enhancement layer distance
  if ( p_Inp->idr_period && !p_Inp->adaptive_idr_period )
  {
    idr_refresh = (( ( p_Img->frm_number - p_Img->lastIDRnumber ) % p_Inp->idr_period ) == 0);
  }
  else if ( p_Inp->idr_period && p_Inp->adaptive_idr_period == 1 )
  {
    idr_refresh = (( ( p_Img->frm_number - imax(p_Img->lastIntraNumber, p_Img->lastIDRnumber) ) % p_Inp->idr_period ) == 0);
  }
  else
  {
    idr_refresh = (p_Img->frm_number == 0);
  }

  idr_flag = ((!p_Img->gop_number) && (!(p_Img->structure==BOTTOM_FIELD)))
    || (idr_refresh && (p_Img->type == I_SLICE || p_Img->type==SI_SLICE)&& (!(p_Img->structure==BOTTOM_FIELD)));

  return (byte) idr_flag;
}

/*!
 ************************************************************************
 * \brief
 *    Update global stats
 ************************************************************************
 */
void update_global_stats(InputParameters *p_Inp, StatParameters *gl_stats, StatParameters *cur_stats)
{  
  int i, j, k;
  
  if (p_Inp->skip_gl_stats == 0)
  {
    for (i = 0; i < 4; i++)
    {
      gl_stats->intra_chroma_mode[i]    += cur_stats->intra_chroma_mode[i];
    }

    for (i = 0; i < 5; i++)
    {
      gl_stats->quant[i]                += cur_stats->quant[i];
      gl_stats->num_macroblocks[i]      += cur_stats->num_macroblocks[i];
      gl_stats->bit_use_mb_type [i]     += cur_stats->bit_use_mb_type[i];
      gl_stats->bit_use_header  [i]     += cur_stats->bit_use_header[i];
      gl_stats->tmp_bit_use_cbp [i]     += cur_stats->tmp_bit_use_cbp[i];
      gl_stats->bit_use_coeffC  [i]     += cur_stats->bit_use_coeffC[i];
      gl_stats->bit_use_coeff[0][i]     += cur_stats->bit_use_coeff[0][i];
      gl_stats->bit_use_coeff[1][i]     += cur_stats->bit_use_coeff[1][i]; 
      gl_stats->bit_use_coeff[2][i]     += cur_stats->bit_use_coeff[2][i]; 
      gl_stats->bit_use_delta_quant[i]  += cur_stats->bit_use_delta_quant[i];
      gl_stats->bit_use_stuffingBits[i] += cur_stats->bit_use_stuffingBits[i];

      for (k = 0; k < 2; k++)
        gl_stats->b8_mode_0_use[i][k] += cur_stats->b8_mode_0_use[i][k];

      for (j = 0; j < 15; j++)
      {
        gl_stats->mode_use[i][j]     += cur_stats->mode_use[i][j];
        gl_stats->bit_use_mode[i][j] += cur_stats->bit_use_mode[i][j];
        for (k = 0; k < 2; k++)
          gl_stats->mode_use_transform[i][j][k] += cur_stats->mode_use_transform[i][j][k];
      }
    }
  }
}

static void storeRedundantFrame(ImageParameters *p_Img)
{
  int j, k;
  if(p_Img->key_frame)
  {
    for(j=0; j<p_Img->height; j++)
    {
      memcpy(p_Img->imgY_tmp[j], p_Img->enc_frame_picture[0]->imgY[j], p_Img->width * sizeof(imgpel));
    }

    for (k = 0; k < 2; k++)
    {
      for(j=0; j<p_Img->height_cr; j++)
      {
        memcpy(p_Img->imgUV_tmp[k][j], p_Img->enc_frame_picture[0]->imgUV[k][j], p_Img->width_cr * sizeof(imgpel));
      }
    }
  }

  if(p_Img->redundant_coding)
  {
    for(j=0; j<p_Img->height; j++)
    {
      memcpy(p_Img->enc_frame_picture[0]->imgY[j], p_Img->imgY_tmp[j], p_Img->width * sizeof(imgpel));
    }
    for (k = 0; k < 2; k++)
    {
      for(j=0; j<p_Img->height_cr; j++)
      {
        memcpy(p_Img->enc_frame_picture[0]->imgUV[k][j], p_Img->imgUV_tmp[k][j], p_Img->width_cr * sizeof(imgpel));
      }
    }
  }
}
/*!
 ************************************************************************
 * \brief
 *    Free storable pictures
 ************************************************************************
 */
void free_pictures(ImageParameters *p_Img, InputParameters *p_Inp, int stored_pic)
{
  int i;
  for (i = 0; i < 6; i++)
  {
    if (i != stored_pic)
      free_storable_picture(p_Img, p_Inp, p_Img->enc_frame_picture[i]);
  }
}

/*!
 ************************************************************************
 * \brief
 *    Encodes one frame
 ************************************************************************
 */
int encode_one_frame (ImageParameters *p_Img, InputParameters *p_Inp)
{
  int prev_frame_no = 0; // POC200301
  int consecutive_non_reference_pictures = 0; // POC200301
  int        i;
  int   nplane;

  //Rate control
  int bits = 0;

  TIME_T start_time;
  TIME_T end_time;
  int64  tmp_time;

  p_Img->me_time = 0;
  p_Img->rd_pass = 0;

  if( IS_INDEPENDENT(p_Inp) )
  {
    for( nplane=0; nplane<MAX_PLANE; nplane++ ){
      p_Img->enc_frame_picture_JV[nplane] = NULL;
    }
  }

  for (i = 0; i < 6; i++)
    p_Img->enc_frame_picture[i]  = NULL;

  gettime(&start_time);          // start time in ms

  //Rate control
  p_Img->write_macroblock = FALSE;
  /*
  //Shankar Regunathan (Oct 2002)
  //Prepare Panscanrect SEI payload
  UpdatePanScanRectInfo (p_SEI);
  //Prepare Arbitrarydata SEI Payload
  UpdateUser_data_unregistered (p_SEI);
  //Prepare Registered data SEI Payload
  UpdateUser_data_registered_itu_t_t35 (p_SEI);
  //Prepare RandomAccess SEI Payload
  UpdateRandomAccess (p_Img);
  */


  put_buffer_frame (p_Img);      // sets the pointers to the frame structures
                               // (and not to one of the field structures)
  init_frame (p_Img, p_Inp);

  ReadOneFrame (p_Img, p_Inp, &p_Inp->input_file1, p_Img->frm_no_in_file, p_Inp->infile_header, &p_Inp->source, &p_Inp->output, p_Img->imgData0.frm_data); 
  PaddAutoCropBorders (p_Inp->output, p_Img->width, p_Img->height, p_Img->width_cr, p_Img->height_cr, p_Img->imgData0.frm_data);


  ProcessImage(p_Img, p_Inp);
  PaddAutoCropBorders (p_Inp->output, p_Img->width, p_Img->height, p_Img->width_cr, p_Img->height_cr, p_Img->imgData.frm_data);




  // Following code should consider optimal coding mode. Currently also does not support
  // multiple slices per frame.

  p_Img->p_Dist->frame_ctr++;

  if(p_Img->type == SP_SLICE)
  {
    if(p_Inp->sp2_frame_indicator)
    { // switching SP frame encoding
      p_Img->sp2_frame_indicator = TRUE;
      read_SP_coefficients(p_Img, p_Inp);
    }
  }
  else
  {
    p_Img->sp2_frame_indicator = FALSE;
  }
  if ( p_Inp->WPMCPrecision )
  {
    wpxInitWPXPasses(p_Img, p_Inp);
    p_Img->pWPX->curr_wp_rd_pass = p_Img->pWPX->wp_rd_passes;
    p_Img->pWPX->curr_wp_rd_pass->algorithm = WP_REGULAR;
  }

  if (p_Inp->PicInterlace == FIELD_CODING)
  {
    //Rate control
    if ( p_Inp->RCEnable && p_Inp->RCUpdateMode <= MAX_RC_MODE )
      p_Img->p_rc_gen->FieldControl = 1;

    p_Img->field_picture = 1;  // we encode fields
    field_picture (p_Img, p_Inp, p_Img->field_pic[0], p_Img->field_pic[1]);
    p_Img->fld_flag = TRUE;
  }
  else
  {
    int tmpFrameQP;
    int num_ref_idx_l0;
    int num_ref_idx_l1;

    //Rate control
    if ( p_Inp->RCEnable && p_Inp->RCUpdateMode <= MAX_RC_MODE )
      p_Img->p_rc_gen->FieldControl = 0;

    p_Img->field_picture = 0; // we encode a frame

    //Rate control
    if(p_Inp->RCEnable)
      rc_init_frame(p_Img, p_Inp);

    if (p_Inp->GenerateMultiplePPS)
      p_Img->active_pps = p_Img->PicParSet[0];

    frame_picture (p_Img, p_Inp, p_Img->frame_pic[0], &p_Img->imgData, 0);

    if(p_Inp->WPIterMC)
      p_Img->frameOffsetAvail = 1; 

    if ((p_Inp->RDPictureIntra || p_Img->type!=I_SLICE) && p_Inp->RDPictureDecision)
    {
      rd_picture_coding(p_Img, p_Inp);
    }

    tmpFrameQP = p_Img->SumFrameQP; // call it here since rd_picture_coding buffers it and may modify it
    num_ref_idx_l0 = p_Img->num_ref_idx_l0_active;
    num_ref_idx_l1 = p_Img->num_ref_idx_l1_active;



    if ((p_Img->type==SP_SLICE) && (p_Img->si_frame_indicator == FALSE) && (p_Inp->si_frame_indicator))
    {
      // once the picture has been encoded as a primary SP frame encode as an SI frame
      p_Img->si_frame_indicator = TRUE;
      frame_picture (p_Img, p_Inp, p_Img->frame_pic_si, &p_Img->imgData, 0);
    }

    if ((p_Img->type == SP_SLICE) && (p_Inp->sp_output_indicator))
    {
      // output the transformed and quantized coefficients (useful for switching SP frames)
      output_SP_coefficients(p_Img, p_Inp);
    }

    if (p_Inp->PicInterlace == ADAPTIVE_CODING)
    {
      //Rate control
      if ( p_Inp->RCEnable && p_Inp->RCUpdateMode <= MAX_RC_MODE )
        p_Img->p_rc_gen->FieldControl=1;
      p_Img->write_macroblock = FALSE;
      p_Img->bot_MB = FALSE;

      p_Img->field_picture = 1;  // we encode fields
      field_picture (p_Img, p_Inp, p_Img->field_pic[0], p_Img->field_pic[1]);

      if(p_Img->rd_pass == 0)
        p_Img->fld_flag = picture_structure_decision (p_Img, p_Img->frame_pic[0], p_Img->field_pic[0], p_Img->field_pic[1]);
      else if(p_Img->rd_pass == 1)
        p_Img->fld_flag = picture_structure_decision (p_Img, p_Img->frame_pic[1], p_Img->field_pic[0], p_Img->field_pic[1]);
      else
        p_Img->fld_flag = picture_structure_decision (p_Img, p_Img->frame_pic[2], p_Img->field_pic[0], p_Img->field_pic[1]);

      if ( p_Img->fld_flag )
      {
        tmpFrameQP = p_Img->SumFrameQP;
        num_ref_idx_l0 = p_Img->num_ref_idx_l0_active;
        num_ref_idx_l1 = p_Img->num_ref_idx_l1_active;
      }

      update_field_frame_contexts (p_Img, p_Img->fld_flag);

      //Rate control
      if ( p_Inp->RCEnable && p_Inp->RCUpdateMode <= MAX_RC_MODE )
        p_Img->p_rc_gen->FieldFrame = !(p_Img->fld_flag) ? 1 : 0;
    }
    else
      p_Img->fld_flag = FALSE;

    p_Img->SumFrameQP = tmpFrameQP;
    p_Img->num_ref_idx_l0_active = num_ref_idx_l0;
    p_Img->num_ref_idx_l1_active = num_ref_idx_l1;
  }

  p_Img->p_Stats->frame_counter++;
  p_Img->p_Stats->frame_ctr[p_Img->type]++;


  // Here, p_Img->structure may be either FRAME or BOTTOM FIELD depending on whether AFF coding is used
  // The picture structure decision changes really only the fld_flag

  if (p_Img->fld_flag)            // field mode (use field when fld_flag=1 only)
  {
    field_mode_buffer (p_Img);
    write_non_vcl_nalu(p_Img, p_Inp);
    writeout_picture  (p_Img, p_Img->field_pic[0]);
    writeout_picture  (p_Img, p_Img->field_pic[1]);
  }
  else                          //frame mode
  {
    frame_mode_buffer (p_Img, p_Inp);

    if (p_Img->type==SP_SLICE && p_Img->si_frame_indicator == 1)
    {
      write_non_vcl_nalu(p_Img, p_Inp);      
      writeout_picture  (p_Img, p_Img->frame_pic_si);
      p_Img->si_frame_indicator = FALSE;
    }
    else
    {
      write_non_vcl_nalu(p_Img, p_Inp);
      writeout_picture  (p_Img, p_Img->frame_pic[p_Img->rd_pass]);
    }
  }

  if (p_Img->frame_pic_si)
  {
    free_slice_list(p_Img->frame_pic_si);
  }

  for (i = 0; i < p_Img->frm_iter; i++)
  {
    if (p_Img->frame_pic[i])
    {
      free_slice_list(p_Img->frame_pic[i]);
    }
  }

  if (p_Img->field_pic)
  {
    for (i = 0; i < 2; i++)
    {
      if (p_Img->field_pic[i])
        free_slice_list(p_Img->field_pic[i]);
    }
  }

  /*
  // Tian Dong (Sept 2002)
  // in frame mode, the newly reconstructed frame has been inserted to the mem buffer
  // and it is time to prepare the spare picture SEI payload.
  if (p_Inp->InterlaceCodingOption == FRAME_CODING
  && p_Inp->SparePictureOption && p_Img->type != B_SLICE)
  CalculateSparePicture ();
  */

  //Rate control
  if ( p_Inp->RCEnable )
  {
    // we could add here a function pointer!
    bits = (int)( p_Img->p_Stats->bit_ctr - p_Img->p_Stats->bit_ctr_n )
      + (int)( p_Img->p_Stats->bit_ctr_filler_data - p_Img->p_Stats->bit_ctr_filler_data_n );

    if ( p_Inp->RCUpdateMode <= MAX_RC_MODE )
      p_Img->rc_update_pict_frame_ptr(p_Img, p_Inp, p_Img->p_rc_quad, p_Img->p_rc_gen, bits);
  }

  if (p_Inp->PicInterlace == FRAME_CODING)
  {
    if ((p_Inp->rdopt == 3) && (p_Img->nal_reference_idc != 0))
    {
      UpdateDecoders (p_Img, p_Inp, p_Img->enc_picture);      // simulate packet losses and move decoded image to reference buffers
    }

    if (p_Inp->RestrictRef)
      UpdatePixelMap (p_Img, p_Inp);
  }

  compute_distortion(p_Img, p_Inp, &p_Img->imgData);

  // redundant pictures: save reconstruction to calculate SNR and replace reference picture
  if(p_Inp->redundant_pic_flag)
  {
    storeRedundantFrame(p_Img);
  }

  if (p_Inp->PicInterlace == ADAPTIVE_CODING)
  {
    if (p_Img->fld_flag)
    {      
      update_global_stats(p_Inp, p_Img->p_Stats, &p_Img->enc_field_picture[0]->stats);
      update_global_stats(p_Inp, p_Img->p_Stats, &p_Img->enc_field_picture[1]->stats);
      // store bottom field
      store_picture_in_dpb (p_Img, p_Inp, p_Img->enc_field_picture[1], &p_Inp->output);
      free_storable_picture(p_Img, p_Inp, p_Img->enc_frame_picture[0]);
      free_storable_picture(p_Img, p_Inp, p_Img->enc_frame_picture[1]);
      free_storable_picture(p_Img, p_Inp, p_Img->enc_frame_picture[2]);
    }
    else
    {
      update_global_stats(p_Inp, p_Img->p_Stats, &p_Img->enc_frame_picture[p_Img->rd_pass]->stats);
      // replace top with frame
      if (p_Img->rd_pass==2)
      {
        replace_top_pic_with_frame(p_Img, p_Inp, p_Img->enc_frame_picture[2], &p_Inp->output);
        free_storable_picture     (p_Img, p_Inp, p_Img->enc_frame_picture[0]);
        free_storable_picture     (p_Img, p_Inp, p_Img->enc_frame_picture[1]);
      }
      else if (p_Img->rd_pass==1)
      {
        replace_top_pic_with_frame(p_Img, p_Inp, p_Img->enc_frame_picture[1], &p_Inp->output);
        free_storable_picture     (p_Img, p_Inp, p_Img->enc_frame_picture[0]);
        free_storable_picture     (p_Img, p_Inp, p_Img->enc_frame_picture[2]);
      }
      else
      {
        if(p_Inp->redundant_pic_flag==0 || (p_Img->key_frame==0))
        {
          replace_top_pic_with_frame(p_Img, p_Inp, p_Img->enc_frame_picture[0], &p_Inp->output);
          free_storable_picture     (p_Img, p_Inp, p_Img->enc_frame_picture[1]);
          free_storable_picture     (p_Img, p_Inp, p_Img->enc_frame_picture[2]);
        }
      }
      free_storable_picture(p_Img, p_Inp, p_Img->enc_field_picture[1]);      
    }
  }
  else
  {
    if (p_Img->fld_flag)
    {
      update_global_stats(p_Inp, p_Img->p_Stats, &p_Img->enc_field_picture[0]->stats);
      update_global_stats(p_Inp, p_Img->p_Stats, &p_Img->enc_field_picture[1]->stats);
      store_picture_in_dpb(p_Img, p_Inp, p_Img->enc_field_picture[1], &p_Inp->output);
    }
    else
    {      
      if ((p_Inp->redundant_pic_flag != 1) || (p_Img->key_frame == 0))
      {
        update_global_stats(p_Inp, p_Img->p_Stats, &p_Img->enc_frame_picture[p_Img->rd_pass]->stats);
        store_picture_in_dpb (p_Img, p_Inp, p_Img->enc_frame_picture[p_Img->rd_pass], &p_Inp->output);
        free_pictures(p_Img, p_Inp, p_Img->rd_pass);
      }
    }
  }

  p_Img->AverageFrameQP = isign(p_Img->SumFrameQP) * ((iabs(p_Img->SumFrameQP) + (int) (p_Img->FrameSizeInMbs >> 1))/ (int) p_Img->FrameSizeInMbs);  

  if ( p_Inp->RCEnable && p_Inp->RCUpdateMode <= MAX_RC_MODE && p_Img->type != B_SLICE && p_Inp->basicunit < p_Img->FrameSizeInMbs )
    p_Img->p_rc_quad->CurrLastQP = p_Img->AverageFrameQP + p_Img->p_rc_quad->bitdepth_qp_scale;

#ifdef _LEAKYBUCKET_
  // Store bits used for this frame and increment counter of no. of coded frames
  if (!p_Img->redundant_coding)
  {
    p_Img->Bit_Buffer[p_Img->total_frame_buffer++] = (long) (p_Img->p_Stats->bit_ctr - p_Img->p_Stats->bit_ctr_n)
      + (long)( p_Img->p_Stats->bit_ctr_filler_data - p_Img->p_Stats->bit_ctr_filler_data_n );
  }
#endif

  // POC200301: Verify that POC coding type 2 is not used if more than one consecutive
  // non-reference frame is requested or if decoding order is different from output order
  if (p_Img->pic_order_cnt_type == 2)
  {
    if (!p_Img->nal_reference_idc) consecutive_non_reference_pictures++;
    else consecutive_non_reference_pictures = 0;

    if (p_Img->frame_no < prev_frame_no || consecutive_non_reference_pictures>1)
      error("POC type 2 cannot be applied for the coding pattern where the encoding /decoding order of pictures are different from the output order.\n", -1);
    prev_frame_no = p_Img->frame_no;
  }

  gettime(&end_time);    // end time in ms
  tmp_time  = timediff(&start_time, &end_time);
  p_Img->tot_time += tmp_time;
  tmp_time  = timenorm(tmp_time);
  p_Img->me_time   = timenorm(p_Img->me_time);

  if (p_Img->p_Stats->bit_ctr_parametersets_n!=0 && p_Inp->Verbose != 3)
    ReportNALNonVLCBits(p_Img, p_Inp, p_Img->p_Stats, tmp_time);

  if (p_Img->frm_number == 0)
    ReportFirstframe(p_Img, p_Inp, p_Img->p_Stats, tmp_time);
  else
  {
    //Rate control
    if(p_Inp->RCEnable)
    {
      if ((!p_Inp->PicInterlace) && (!p_Inp->MbInterlace))
      {
        bits = (int) (p_Img->p_Stats->bit_ctr - p_Img->p_Stats->bit_ctr_n)
          + (int)( p_Img->p_Stats->bit_ctr_filler_data - p_Img->p_Stats->bit_ctr_filler_data_n );
      }
      else if ( p_Inp->RCUpdateMode <= MAX_RC_MODE )
      {
        bits = (int)(p_Img->p_Stats->bit_ctr - (p_Img->p_rc_quad->Pprev_bits))
          + (int)( p_Img->p_Stats->bit_ctr_filler_data - p_Img->p_Stats->bit_ctr_filler_data_n ); // used for rate control update
        p_Img->p_rc_quad->Pprev_bits = p_Img->p_Stats->bit_ctr + p_Img->p_Stats->bit_ctr_filler_data;
      }
    }

    p_Img->p_Stats->bit_counter[p_Img->type] += p_Img->p_Stats->bit_ctr - p_Img->p_Stats->bit_ctr_n;

    switch (p_Img->type)
    {
    case I_SLICE:
    case SI_SLICE:
      ReportI(p_Img, p_Inp, p_Img->p_Stats, tmp_time);
      break;
    case B_SLICE:
      ReportB(p_Img, p_Inp, p_Img->p_Stats, tmp_time);
      break;
    default:      // P
      ReportP(p_Img, p_Inp, p_Img->p_Stats, tmp_time);
    }
  }

  if (p_Inp->Verbose == 0)
  {
    //for (i = 0; i <= (p_Img->number & 0x0F); i++)
    //printf(".");
    //printf("                              \r");
    printf("Completed Encoding Frame %05d.\r", p_Img->frame_no);
  }
  // Flush output statistics
  fflush(stdout);

  //Rate control
  if(p_Inp->RCEnable)
    p_Img->rc_update_picture_ptr( p_Img, p_Inp, bits );

  // update bit counters
  p_Img->p_Stats->bit_ctr_n = p_Img->p_Stats->bit_ctr;
  p_Img->p_Stats->bit_ctr_parametersets_n = 0;
  p_Img->p_Stats->bit_ctr_filler_data_n = p_Img->p_Stats->bit_ctr_filler_data;

  if ( p_Img->type == I_SLICE && p_Img->nal_reference_idc)
  {
    //p_Img->lastINTRA = p_Img->frame_no;
    // Lets also handle the possibility of backward GOPs and hierarchical structures 
    if ( !(p_Img->b_frame_to_code) )
    {
      p_Img->lastINTRA       = imax(p_Img->lastINTRA, p_Img->frame_no);
      p_Img->lastIntraNumber = p_Img->frm_number;
    }
    if ( p_Img->currentPicture->idr_flag )
    {
      p_Img->lastIDRnumber = p_Img->frm_number;
    }
  }

  return ((p_Img->gop_number == 0)? 0 : 1);
}


/*!
 ************************************************************************
 * \brief
 *    This function write out a picture
 * \return
 *    0 if OK,                                                         \n
 *    1 in case of error
 *
 ************************************************************************
 */
static void writeout_picture(ImageParameters *p_Img, Picture *pic)
{
  int partition, slice;
  Slice  *currSlice;

  p_Img->currentPicture = pic;

  // loop over all slices of the picture
  for (slice=0; slice < pic->no_slices; slice++)
  {
    currSlice = pic->slices[slice];

    // loop over the partitions
    for (partition=0; partition < currSlice->max_part_nr; partition++)
    {
      // write only if the partition has content
      if (currSlice->partArr[partition].bitstream->write_flag )
      {
        p_Img->p_Stats->bit_ctr += p_Img->WriteNALU (p_Img, currSlice->partArr[partition].nal_unit);
      }
    }
  }
}


void copy_params(ImageParameters *p_Img, StorablePicture *enc_picture, seq_parameter_set_rbsp_t *active_sps)
{
  enc_picture->frame_mbs_only_flag = active_sps->frame_mbs_only_flag;
  enc_picture->frame_cropping_flag = active_sps->frame_cropping_flag;
  enc_picture->chroma_format_idc   = active_sps->chroma_format_idc;
  enc_picture->chroma_mask_mv_x    = p_Img->chroma_mask_mv_x;
  enc_picture->chroma_mask_mv_y    = p_Img->chroma_mask_mv_y;
  enc_picture->chroma_shift_y      = p_Img->chroma_shift_y;
  enc_picture->chroma_shift_x      = p_Img->chroma_shift_x;


  if (active_sps->frame_cropping_flag)
  {
    enc_picture->frame_cropping_rect_left_offset   = active_sps->frame_cropping_rect_left_offset;
    enc_picture->frame_cropping_rect_right_offset  = active_sps->frame_cropping_rect_right_offset;
    enc_picture->frame_cropping_rect_top_offset    = active_sps->frame_cropping_rect_top_offset;
    enc_picture->frame_cropping_rect_bottom_offset = active_sps->frame_cropping_rect_bottom_offset;
  }
  else
  {
    enc_picture->frame_cropping_rect_left_offset   = 0;
    enc_picture->frame_cropping_rect_right_offset  = 0;
    enc_picture->frame_cropping_rect_top_offset    = 0;
    enc_picture->frame_cropping_rect_bottom_offset = 0;
  }
}

/*!
 ************************************************************************
 * \brief
 *    Prepare and allocate an encoded frame picture structure
 ************************************************************************
 */
static void prepare_enc_frame_picture (ImageParameters *p_Img, InputParameters *p_Inp, StorablePicture **stored_pic)
{
  (*stored_pic)                 = alloc_storable_picture (p_Img, p_Inp, (PictureStructure) p_Img->structure, p_Img->width, p_Img->height, p_Img->width_cr, p_Img->height_cr);
  
  p_Img->ThisPOC                = p_Img->framepoc;
  (*stored_pic)->poc            = p_Img->framepoc;
  (*stored_pic)->top_poc        = p_Img->toppoc;
  (*stored_pic)->bottom_poc     = p_Img->bottompoc;
  (*stored_pic)->frame_poc      = p_Img->framepoc;
  (*stored_pic)->pic_num        = p_Img->frame_num;
  (*stored_pic)->frame_num      = p_Img->frame_num;
  (*stored_pic)->coded_frame    = 1;
  (*stored_pic)->MbaffFrameFlag = p_Img->MbaffFrameFlag = (Boolean) (p_Inp->MbInterlace != FRAME_CODING);
  
  p_Img->get_mb_block_pos           = p_Img->MbaffFrameFlag ? get_mb_block_pos_mbaff : get_mb_block_pos_normal;
  p_Img->getNeighbour        = p_Img->MbaffFrameFlag ? getAffNeighbour : getNonAffNeighbour;
  p_Img->enc_picture         = *stored_pic;

  copy_params(p_Img, p_Img->enc_picture, p_Img->active_sps);
}

static void calc_picture_bits(Picture *frame)
{
  int i, j;
  Slice *thisSlice = NULL;

  frame->bits_per_picture = 0;

  for ( i = 0; i < frame->no_slices; i++ )
  {
    thisSlice = frame->slices[i];

    for ( j = 0; j < thisSlice->max_part_nr; j++ )
      frame->bits_per_picture += 8 * ((thisSlice->partArr[j]).bitstream)->byte_pos;
  }
}
/*!
 ************************************************************************
 * \brief
 *    Encodes a frame picture
 ************************************************************************
 */
void frame_picture (ImageParameters *p_Img, InputParameters *p_Inp, Picture *frame, ImageData *imgData, int rd_pass)
{
  int nplane;
  p_Img->SumFrameQP = 0;
  p_Img->num_ref_idx_l0_active = 0;
  p_Img->num_ref_idx_l1_active = 0;
  p_Img->structure = FRAME;
  p_Img->PicSizeInMbs = p_Img->FrameSizeInMbs;
  //set mv limits to frame type
  update_mv_limits(p_Img, p_Inp, FALSE);


  if( IS_INDEPENDENT(p_Inp) )
  {
    for( nplane=0; nplane<MAX_PLANE; nplane++ )
    {
      prepare_enc_frame_picture( p_Img, p_Inp, &p_Img->enc_frame_picture_JV[nplane] );      
    }
  }
  else
  {
    prepare_enc_frame_picture( p_Img, p_Inp, &p_Img->enc_frame_picture[rd_pass] );
  }


  p_Img->fld_flag = FALSE;
  code_a_picture(p_Img, p_Inp, frame);

  if( IS_INDEPENDENT(p_Inp) )
  {
    make_frame_picture_JV(p_Img, p_Inp);
  }

  calc_picture_bits(frame);

  if (p_Img->structure==FRAME)
  {
    find_distortion (p_Img, p_Inp, imgData);
    frame->distortion = p_Img->p_Dist->metric[SSE];
  }
}


/*!
 ************************************************************************
 * \brief
 *    Encodes a field picture, consisting of top and bottom field
 ************************************************************************
 */
static void field_picture (ImageParameters *p_Img, InputParameters *p_Inp, Picture *top, Picture *bottom)
{
  //Rate control
  int old_pic_type;              // picture type of top field used for rate control
  int TopFieldBits;
  p_Img->SumFrameQP = 0;
  p_Img->num_ref_idx_l0_active = 0;
  p_Img->num_ref_idx_l1_active = 0;

  //set mv limits to field type
  update_mv_limits(p_Img, p_Inp, TRUE);

  //Rate control
  old_pic_type = p_Img->type;

  p_Img->number *= 2;
  p_Img->gop_number = (p_Img->number - p_Img->start_frame_no);
  p_Img->buf_cycle *= 2;
  p_Img->height    = (p_Inp->output.height + p_Img->auto_crop_bottom) >> 1;
  p_Img->height_cr = p_Img->height_cr_frame >> 1;
  p_Img->fld_flag  = TRUE;
  p_Img->PicSizeInMbs = p_Img->FrameSizeInMbs >> 1;
  // Top field

  p_Img->enc_field_picture[0]              = alloc_storable_picture (p_Img, p_Inp, (PictureStructure) p_Img->structure, p_Img->width, p_Img->height, p_Img->width_cr, p_Img->height_cr);
  p_Img->enc_field_picture[0]->poc         = p_Img->toppoc;
  p_Img->enc_field_picture[0]->frame_poc   = p_Img->toppoc;
  p_Img->enc_field_picture[0]->pic_num     = p_Img->frame_num;
  p_Img->enc_field_picture[0]->frame_num   = p_Img->frame_num;
  p_Img->enc_field_picture[0]->coded_frame = 0;
  p_Img->enc_field_picture[0]->MbaffFrameFlag = p_Img->MbaffFrameFlag = FALSE;
  p_Img->get_mb_block_pos = get_mb_block_pos_normal;
  p_Img->getNeighbour = getNonAffNeighbour;
  p_Img->ThisPOC = p_Img->toppoc;

  p_Img->structure = TOP_FIELD;
  p_Img->enc_picture = p_Img->enc_field_picture[0];
  copy_params(p_Img, p_Img->enc_picture, p_Img->active_sps);

  put_buffer_top (p_Img);
  init_field (p_Img, p_Inp);

  p_Img->fld_flag = TRUE;

  //Rate control
  if(p_Inp->RCEnable && p_Inp->RCUpdateMode <= MAX_RC_MODE)
    rc_init_top_field(p_Img, p_Inp);

  code_a_picture(p_Img, p_Inp, p_Img->field_pic[0]);
  p_Img->enc_picture->structure = TOP_FIELD;

  store_picture_in_dpb(p_Img, p_Inp, p_Img->enc_field_picture[0], &p_Inp->output);

  calc_picture_bits(top);

  //Rate control
  TopFieldBits=top->bits_per_picture;

  //  Bottom field
  p_Img->enc_field_picture[1]  = alloc_storable_picture (p_Img, p_Inp, (PictureStructure) p_Img->structure, p_Img->width, p_Img->height, p_Img->width_cr, p_Img->height_cr);
  p_Img->enc_field_picture[1]->poc=p_Img->bottompoc;
  p_Img->enc_field_picture[1]->frame_poc = p_Img->bottompoc;
  p_Img->enc_field_picture[1]->pic_num = p_Img->frame_num;
  p_Img->enc_field_picture[1]->frame_num = p_Img->frame_num;
  p_Img->enc_field_picture[1]->coded_frame = 0;
  p_Img->enc_field_picture[1]->MbaffFrameFlag = p_Img->MbaffFrameFlag = FALSE;
  p_Img->get_mb_block_pos = get_mb_block_pos_normal;
  p_Img->getNeighbour = getNonAffNeighbour;

  p_Img->ThisPOC = p_Img->bottompoc;
  p_Img->structure = BOTTOM_FIELD;
  p_Img->enc_picture = p_Img->enc_field_picture[1];
  copy_params(p_Img, p_Img->enc_picture, p_Img->active_sps);
  put_buffer_bot (p_Img);
  p_Img->number++;
  p_Img->gop_number = (p_Img->number - p_Img->start_frame_no);
  
  init_field (p_Img, p_Inp);

 if (p_Img->type == I_SLICE && p_Inp->IntraBottom!=1)
   set_slice_type(p_Img, p_Inp, (p_Inp->BRefPictures == 2) ? B_SLICE : P_SLICE);

  p_Img->fld_flag = TRUE;

  //Rate control
  if(p_Inp->RCEnable && p_Inp->RCUpdateMode <= MAX_RC_MODE)
    rc_init_bottom_field( p_Img, p_Inp, TopFieldBits );

  p_Img->enc_picture->structure = BOTTOM_FIELD;
  code_a_picture(p_Img, p_Inp, p_Img->field_pic[1]);

  calc_picture_bits(bottom);

  // the distortion for a field coded frame (consisting of top and bottom field)
  // lives in the top->distortion variables, the bottom-> are dummies
  distortion_fld (p_Img, p_Inp, top, &p_Img->imgData);
}

/*!
 ************************************************************************
 * \brief
 *    form frame picture from two field pictures
 ************************************************************************
 */
static void combine_field(ImageParameters *p_Img)
{
  int i, k;

  for (i = 0; i < (p_Img->height >> 1); i++)
  {
    memcpy(p_Img->imgY_com[i*2],     p_Img->enc_field_picture[0]->imgY[i], p_Img->width*sizeof(imgpel));     // top field
    memcpy(p_Img->imgY_com[i*2 + 1], p_Img->enc_field_picture[1]->imgY[i], p_Img->width*sizeof(imgpel)); // bottom field
  }

  if (p_Img->yuv_format != YUV400)
  {
    for (k = 0; k < 2; k++)
    {
      for (i = 0; i < (p_Img->height_cr >> 1); i++)
      {
        memcpy(p_Img->imgUV_com[k][i*2],     p_Img->enc_field_picture[0]->imgUV[k][i], p_Img->width_cr*sizeof(imgpel));
        memcpy(p_Img->imgUV_com[k][i*2 + 1], p_Img->enc_field_picture[1]->imgUV[k][i], p_Img->width_cr*sizeof(imgpel));
      }
    }
  }
}

/*!
 ************************************************************************
 * \brief
 *    Distortion Field
 ************************************************************************
 */
static void distortion_fld (ImageParameters *p_Img, InputParameters *p_Inp, Picture *field_pic, ImageData *imgData)
{

  p_Img->number /= 2;
  p_Img->gop_number = (p_Img->number - p_Img->start_frame_no);
  p_Img->buf_cycle /= 2;
  p_Img->height    = (p_Inp->output.height + p_Img->auto_crop_bottom);
  p_Img->height_cr = p_Img->height_cr_frame;

  combine_field (p_Img);

  p_Img->pCurImg   = imgData->frm_data[0];
  p_Img->pImgOrg[0] = imgData->frm_data[0];

  if (p_Inp->output.yuv_format != YUV400)
  {
    p_Img->pImgOrg[1] = imgData->frm_data[1];
    p_Img->pImgOrg[2] = imgData->frm_data[2];
  }

  find_distortion (p_Img, p_Inp, imgData);   // find snr from original frame picture
  field_pic->distortion = p_Img->p_Dist->metric[SSE];
}


/*!
 ************************************************************************
 * \brief
 *    RD decision of frame and field coding
 ************************************************************************
 */
static byte decide_fld_frame(float snr_frame_Y, float snr_field_Y, int bit_field, int bit_frame, double lambda_picture)
{
  double cost_frame, cost_field;

  cost_frame = bit_frame * lambda_picture + snr_frame_Y;
  cost_field = bit_field * lambda_picture + snr_field_Y;

  if (cost_field > cost_frame)
    return FALSE;
  else
    return TRUE;
}

/*!
 ************************************************************************
 * \brief
 *    Picture Structure Decision
 ************************************************************************
 */
static byte picture_structure_decision (ImageParameters *p_Img, Picture *frame, Picture *top, Picture *bot)
{
  double lambda_picture;
  float sse_frame, sse_field;
  int bit_frame, bit_field;

  lambda_picture = 0.68 * pow (2, p_Img->bitdepth_lambda_scale + ((p_Img->qp - SHIFT_QP) / 3.0)) * ((p_Img->type == B_SLICE) ? 1 : 1);

  sse_frame = frame->distortion.value[0] + frame->distortion.value[1] + frame->distortion.value[2];
  //! all distrortions of a field picture are accumulated in the top field
  sse_field = top->distortion.value[0] + top->distortion.value[1] + top->distortion.value[2];

  bit_field = top->bits_per_picture + bot->bits_per_picture;
  bit_frame = frame->bits_per_picture;
  return decide_fld_frame (sse_frame, sse_field, bit_field, bit_frame, lambda_picture);
}


/*!
 ************************************************************************
 * \brief
 *    Field Mode Buffer
 ************************************************************************
 */
static void field_mode_buffer (ImageParameters *p_Img)
{
  put_buffer_frame (p_Img);
}


/*!
 ************************************************************************
 * \brief
 *    Frame Mode Buffer
 ************************************************************************
 */
static void frame_mode_buffer (ImageParameters *p_Img, InputParameters *p_Inp)
{
  put_buffer_frame (p_Img);

  if ((p_Inp->PicInterlace != FRAME_CODING)||(p_Inp->MbInterlace != FRAME_CODING))
  {
    p_Img->height = p_Img->height / 2;
    p_Img->height_cr = p_Img->height_cr / 2;
    p_Img->number *= 2;
    p_Img->gop_number = (p_Img->number - p_Img->start_frame_no);

    put_buffer_top (p_Img);

    p_Img->number++;
    p_Img->gop_number = (p_Img->number - p_Img->start_frame_no);
    put_buffer_bot (p_Img);

    p_Img->number /= 2;         // reset the p_Img->number to field
    p_Img->gop_number = (p_Img->number - p_Img->start_frame_no);
    p_Img->height = (p_Inp->output.height + p_Img->auto_crop_bottom);
    p_Img->height_cr = p_Img->height_cr_frame;

    put_buffer_frame (p_Img);
  }
}


/*!
 ************************************************************************
 * \brief
 *    mmco initializations should go here
 ************************************************************************
 */
static void init_dec_ref_pic_marking_buffer(ImageParameters *p_Img)
{
  p_Img->dec_ref_pic_marking_buffer=NULL;
}


/*!
 ************************************************************************
 * \brief
 *    Initializes the parameters for a new frame
 ************************************************************************
 */
static void init_frame (ImageParameters *p_Img, InputParameters *p_Inp)
{
  int i, j;

  p_Img->current_mb_nr = 0;
  p_Img->current_slice_nr = 0;
  p_Img->p_Stats->bit_slice = 0;

  // The 'slice_nr' of each macroblock is set to -1 here, to guarantee the correct encoding
  // with FMO (if no FMO, encoding is correct without following assignment),
  // for which MBs may not be encoded with scan order
  if( IS_INDEPENDENT(p_Inp) )
  {
    for( j=0; j<MAX_PLANE; j++ ){
      for(i=0;i< ((int) (p_Img->FrameSizeInMbs));i++)
        p_Img->mb_data_JV[j][i].slice_nr=-1;
    }
  }
  else
  {
    for(i = 0; i < ((int) (p_Img->FrameSizeInMbs)); i++)
      p_Img->mb_data[i].slice_nr = -1;
  }

  if (p_Img->b_frame_to_code == 0)
  {
    //Rate control
    if(!p_Inp->RCEnable)                  // without using rate control
    {
      if (p_Img->type == I_SLICE)
      {
        //QP oscillation for secondary SP frames
        if ((p_Inp->qp2start > 0 && p_Img->frame_no >= p_Inp->qp2start && p_Inp->sp2_frame_indicator==0)||
          ( (p_Inp->qp2start > 0) && ( ( (p_Img->frame_no ) % (2*p_Inp->qp2start) ) >=p_Inp->qp2start ) && (p_Inp->sp2_frame_indicator==1)))
          p_Img->qp = p_Inp->qp[1][I_SLICE];
        else
          p_Img->qp = p_Inp->qp[0][I_SLICE];   // set quant. parameter for I-frame

        if (p_Img->redundant_coding)
        {
          //!KS: hard code qp increment
          p_Img->qp = imin(p_Img->qp + 5, 51);
        }
      }
      else
      {
        //QP oscillation for secondary SP frames
        if ((p_Inp->qp2start > 0 && p_Img->frame_no >= p_Inp->qp2start && p_Inp->sp2_frame_indicator==0)||
          ( (p_Inp->qp2start > 0) && ( ( (p_Img->frame_no ) % (2*p_Inp->qp2start) ) >=p_Inp->qp2start ) && (p_Inp->sp2_frame_indicator==1)))
          p_Img->qp = p_Inp->qp[1][P_SLICE] + (p_Img->nal_reference_idc ? 0 : p_Inp->DispPQPOffset);
        else
          p_Img->qp = p_Inp->qp[0][P_SLICE] + (p_Img->nal_reference_idc ? 0 : p_Inp->DispPQPOffset);

        if (p_Img->type == SP_SLICE)
        {
          if ( (p_Inp->qp2start > 0) && ( ( (p_Img->frame_no ) % (2*p_Inp->qp2start) ) >=p_Inp->qp2start ))
          {
            p_Img->qp   = p_Inp->qp[1][SP_SLICE];
            p_Img->qpsp = p_Inp->qpsp[1];
          }
          else
          {
            p_Img->qp   = p_Inp->qp[0][SP_SLICE];
            p_Img->qpsp = p_Inp->qpsp[0];
          }
        }
      }
    }

    p_Img->mb_y_intra = p_Img->mb_y_upd;  //  p_Img->mb_y_intra indicates which GOB to intra code for this frame

    if (p_Inp->intra_upd > 0) // if error robustness, find next GOB to update
    {
      p_Img->mb_y_upd = (p_Img->frm_number / p_Inp->intra_upd) % (p_Img->height / MB_BLOCK_SIZE);
    }
  }
  else
  {
    //Rate control
    if(!p_Inp->RCEnable && p_Inp->HierarchicalCoding == 0)                  // without using rate control
    {
      //QP oscillation for secondary SP frames
      if ((p_Inp->qp2start > 0 && p_Img->frame_no >= p_Inp->qp2start && p_Inp->sp2_frame_indicator==0)||
        ( (p_Inp->qp2start > 0) && ( ( (p_Img->frame_no ) % (2*p_Inp->qp2start) ) >=p_Inp->qp2start ) && (p_Inp->sp2_frame_indicator==1)))
      {
        p_Img->qp = p_Inp->qp[1][B_SLICE];
      }
      else
      {
        p_Img->qp = p_Inp->qp[0][B_SLICE];
      }

      if (p_Img->nal_reference_idc)
      {
        //QP oscillation for secondary SP frames
        if ((p_Inp->qp2start > 0 && p_Img->frame_no >= p_Inp->qp2start && p_Inp->sp2_frame_indicator==0)||
          ( (p_Inp->qp2start > 0) && ( ( (p_Img->frame_no ) % (2*p_Inp->qp2start) ) >=p_Inp->qp2start ) && (p_Inp->sp2_frame_indicator==1)))
        {
          p_Img->qp = iClip3(-p_Img->bitdepth_luma_qp_scale,51,p_Inp->qp[1][B_SLICE] + p_Inp->qpBRSOffset[1]);
        }
        else
        {
          p_Img->qp = iClip3(-p_Img->bitdepth_luma_qp_scale, 51, p_Inp->qp[0][B_SLICE] + p_Inp->qpBRSOffset[0]);
        }
      }
    }
    else if (!p_Inp->RCEnable && p_Inp->HierarchicalCoding !=0)
    {
      // Note that _CHANGE_QP_ does not work anymore for p_Img->gop_structure. Needs to be fixed
      p_Img->qp =  p_Img->gop_structure[p_Img->b_frame_to_code - 1].slice_qp;

    }
  }

  p_Img->no_output_of_prior_pics_flag = 0;
  p_Img->long_term_reference_flag = FALSE;

  init_dec_ref_pic_marking_buffer(p_Img);

  if(p_Inp->WPIterMC)
    p_Img->frameOffsetAvail = 0;

  // set parameters for direct mode and deblocking filter
  // currently selection is done at the frame level instead of slice level. This needs to be changed.  
  p_Img->direct_spatial_mv_pred_flag = (char) p_Inp->direct_spatial_mv_pred_flag;
  p_Img->DFDisableIdc                = (char) p_Inp->DFDisableIdc[p_Img->nal_reference_idc > 0][p_Img->type];
  p_Img->DFAlphaC0Offset             = (char) p_Inp->DFAlpha     [p_Img->nal_reference_idc > 0][p_Img->type];
  p_Img->DFBetaOffset                = (char) p_Inp->DFBeta      [p_Img->nal_reference_idc > 0][p_Img->type];

  p_Img->AdaptiveRounding            = p_Inp->AdaptiveRounding; 
}

/*!
 ************************************************************************
 * \brief
 *    Initializes the parameters for a new field
 ************************************************************************
 */
static void init_field (ImageParameters *p_Img, InputParameters *p_Inp)
{
  p_Img->current_mb_nr = 0;
  p_Img->current_slice_nr = 0;
  p_Img->p_Stats->bit_slice = 0;

  p_Inp->jumpd *= 2;
  p_Inp->NumberBFrames *= 2;
  p_Img->number /= 2;
  p_Img->gop_number = (p_Img->number - p_Img->start_frame_no);
  p_Img->buf_cycle /= 2;

  if (!p_Img->b_frame_to_code)
  {
      //Rate control
    if(!p_Inp->RCEnable)                  // without using rate control
    {
      if (p_Img->type == I_SLICE)
      {
        //QP oscillation for secondary SP frames
        if ((p_Inp->qp2start > 0 && p_Img->frame_no >= p_Inp->qp2start && p_Inp->sp2_frame_indicator==0)||
          ( (p_Inp->qp2start > 0) && ( ( (p_Img->frame_no ) % (2*p_Inp->qp2start) ) >=p_Inp->qp2start ) && (p_Inp->sp2_frame_indicator==1)))
          p_Img->qp = p_Inp->qp[1][I_SLICE];
        else
          p_Img->qp = p_Inp->qp[0][I_SLICE];   // set quant. parameter for I-frame
      }
      else
      {
        //QP oscillation for secondary SP frames
        if ((p_Inp->qp2start > 0 && p_Img->frame_no >= p_Inp->qp2start && p_Inp->sp2_frame_indicator==0)||
          ( (p_Inp->qp2start > 0) && ( ( (p_Img->frame_no ) % (2*p_Inp->qp2start) ) >=p_Inp->qp2start ) && (p_Inp->sp2_frame_indicator==1)))
          p_Img->qp = p_Inp->qp[1][P_SLICE] + (p_Img->nal_reference_idc ? 0 : p_Inp->DispPQPOffset);
        else
          p_Img->qp = p_Inp->qp[0][P_SLICE] + (p_Img->nal_reference_idc ? 0 : p_Inp->DispPQPOffset);

        if (p_Img->type == SP_SLICE)
        {
          if ( (p_Inp->qp2start > 0) && ( ( (p_Img->frame_no ) % (2*p_Inp->qp2start) ) >=p_Inp->qp2start ))
          {
            p_Img->qp = p_Inp->qp[1][SP_SLICE];
            p_Img->qpsp = p_Inp->qpsp[1];
          }
          else
          {
            p_Img->qp = p_Inp->qp[0][SP_SLICE];
            p_Img->qpsp = p_Inp->qpsp[0];
          }
        }
      }
    }
    p_Img->mb_y_intra = p_Img->mb_y_upd;  //  p_Img->mb_y_intra indicates which GOB to intra code for this frame

    if (p_Inp->intra_upd > 0) // if error robustness, find next GOB to update
    {
      p_Img->mb_y_upd =
        (p_Img->number / p_Inp->intra_upd) % (p_Img->width / MB_BLOCK_SIZE);
    }
  }
  else
  {
    //Rate control
    if(!p_Inp->RCEnable && p_Inp->HierarchicalCoding == 0)                  // without using rate control
    {
      //QP oscillation for secondary SP frames
      if ((p_Inp->qp2start > 0 && p_Img->frame_no >= p_Inp->qp2start && p_Inp->sp2_frame_indicator==0)||
        ( (p_Inp->qp2start > 0) && ( ( (p_Img->frame_no ) % (2*p_Inp->qp2start) ) >=p_Inp->qp2start ) && (p_Inp->sp2_frame_indicator==1)))
      {
        p_Img->qp = p_Inp->qp[1][B_SLICE];
      }
      else
        p_Img->qp = p_Inp->qp[0][B_SLICE];
      if (p_Img->nal_reference_idc)
      {
        //QP oscillation for secondary SP frames
        if ((p_Inp->qp2start > 0 && p_Img->frame_no >= p_Inp->qp2start && p_Inp->sp2_frame_indicator==0)||
          ( (p_Inp->qp2start > 0) && ( ( (p_Img->frame_no ) % (2*p_Inp->qp2start) ) >=p_Inp->qp2start ) && (p_Inp->sp2_frame_indicator==1)))
        {
          p_Img->qp = iClip3(-p_Img->bitdepth_luma_qp_scale,51,p_Inp->qp[1][B_SLICE] + p_Inp->qpBRSOffset[1]);
        }
        else
          p_Img->qp = iClip3(-p_Img->bitdepth_luma_qp_scale,51,p_Inp->qp[0][B_SLICE] + p_Inp->qpBRSOffset[0]);

      }
    }
    else if (!p_Inp->RCEnable && p_Inp->HierarchicalCoding != 0)
    {
      p_Img->qp =  p_Img->gop_structure[p_Img->b_frame_to_code - 1].slice_qp;
    }
  }
  p_Inp->jumpd /= 2;
  p_Inp->NumberBFrames /= 2;
  p_Img->buf_cycle *= 2;
  p_Img->number = 2 * p_Img->number + p_Img->fld_type;
  p_Img->gop_number = (p_Img->number - p_Img->start_frame_no);
}

/*!
 ************************************************************************
 * \brief
 *    Upsample 4 times, store them in out4x.  Color is simply copied
 *
 * \par Input:
 *    srcy, srcu, srcv, out4y, out4u, out4v
 *
 * \par Side Effects_
 *    Uses (writes) img4Y_tmp.  This should be moved to a static variable
 *    in this module
 ************************************************************************/
void UnifiedOneForthPix ( ImageParameters *p_Img, InputParameters *p_Inp, StorablePicture *s)
{
  int ypadded_size = s->size_y_padded;
  int xpadded_size = s->size_x_padded;

  // don't upsample twice
  if (s->imgY_sub)
    return;
  // Y component
  get_mem4Dpel (&(s->imgY_sub), 4, 4, ypadded_size, xpadded_size);
  if (NULL == s->imgY_sub)
    no_mem_exit("alloc_storable_picture: s->imgY_sub");
  s->p_img_sub[0] = s->imgY_sub;
  s->p_curr_img_sub = s->imgY_sub;
  s->p_curr_img = s->imgY;

  if ( p_Inp->ChromaMCBuffer || p_Img->P444_joined)
  {
    // UV components
    if ( p_Img->yuv_format != YUV400 )
    {
      if ( p_Img->yuv_format == YUV420 )
      {
        get_mem5Dpel (&(s->imgUV_sub), 2, 8, 8, ypadded_size>>1, xpadded_size>>1);
      }
      else if ( p_Img->yuv_format == YUV422 )
      {
        get_mem5Dpel (&(s->imgUV_sub), 2, 4, 8, ypadded_size, xpadded_size>>1);
      }
      else
      { // YUV444
        get_mem5Dpel (&(s->imgUV_sub), 2, 4, 4, ypadded_size, xpadded_size);
      }
      s->p_img_sub[1] = s->imgUV_sub[0];
      s->p_img_sub[2] = s->imgUV_sub[1];
    }
    else
    {
      s->p_img_sub[1] = NULL;
      s->p_img_sub[2] = NULL;
    }
  }
  else
  {
    s->p_img_sub[1] = NULL;
    s->p_img_sub[2] = NULL;
  }

  // derive the subpixel images for first component
  // No need to interpolate if intra only encoding
  //if (p_Inp->intra_period != 1)
  {
    getSubImagesLuma ( p_Img, p_Inp, s );

    // and the sub-images for U and V
    if ( (p_Img->yuv_format != YUV400) && (p_Inp->ChromaMCBuffer) )
    {
      if (p_Img->P444_joined)
      {
        //U
        select_plane(p_Img, PLANE_U);
        getSubImagesLuma (p_Img, p_Inp, s);
        //V
        select_plane(p_Img, PLANE_V);
        getSubImagesLuma (p_Img, p_Inp, s);
        //Y
        select_plane(p_Img, PLANE_Y);
      }
      else
        getSubImagesChroma( p_Img, p_Inp, s );
    }
  }
}

/*!
 ************************************************************************
 * \brief
 *    Upsample 4 times, store them in out4x.  Color is simply copied
 *    for 4:4:4 Independent mode
 *
 * \par Input:
 *    nplane
 *
 ************************************************************************/
void UnifiedOneForthPix_JV (ImageParameters *p_Img, InputParameters *p_Inp, int nplane, StorablePicture *s)
{
  int ypadded_size = s->size_y_padded;
  int xpadded_size = s->size_x_padded;

  if( nplane == 0 )
  {
    // don't upsample twice
    if (s->imgY_sub)
      return;
    // Y component
    get_mem4Dpel (&(s->imgY_sub), 4, 4, ypadded_size, xpadded_size);
    if (NULL == s->imgY_sub)
      no_mem_exit("UnifiedOneForthPix_JV: s->imgY_sub");

    // don't upsample twice
    if (s->imgUV_sub)
      return;
    // Y component
    get_mem5Dpel (&(s->imgUV_sub), 2, 4, 4, ypadded_size, xpadded_size);
    if (NULL == s->imgUV_sub)
      no_mem_exit("UnifiedOneForthPix_JV: s->imgUV_sub");

    s->p_img[0] = s->imgY;
    s->p_img[1] = s->imgUV[0];
    s->p_img[2] = s->imgUV[1];

    s->p_img_sub[0] = s->imgY_sub;
    s->p_img_sub[1] = s->imgUV_sub[0];
    s->p_img_sub[2] = s->imgUV_sub[1];
  }

  // derive the subpixel images for first component
  s->colour_plane_id = nplane;
  s->p_curr_img = s->p_img[nplane];
  s->p_curr_img_sub = s->p_img_sub[nplane];

  getSubImagesLuma ( p_Img, p_Inp, s );
}

  /*!
 ************************************************************************
 * \brief
 *    Just a placebo
 ************************************************************************
 */
Boolean dummy_slice_too_big (int bits_slice)
{
  return FALSE;
}


static void ReportSimple(ImageParameters *p_Img, char *pic_type, int cur_bits, DistMetric *metric, int tmp_time)
{
  printf ("%05d(%3s)%8d   %2d %7.3f %7.3f %7.3f %9d %7d    %3s    %d\n",
    p_Img->frame_no, pic_type, cur_bits, 
    p_Img->AverageFrameQP,
    metric->value[0], metric->value[1], metric->value[2], 
    tmp_time, (int) p_Img->me_time,
    p_Img->fld_flag ? "FLD" : "FRM", 
    p_Img->nal_reference_idc);
}

static void ReportVerbose(ImageParameters *p_Img, char *pic_type, int cur_bits, int wp_method, int lambda, DistMetric *mPSNR, int tmp_time, int direct_mode)
{
  printf ("%05d(%3s)%8d %1d %2d %2d %7.3f %7.3f %7.3f %9d %7d    %3s %5d %1d %2d %2d  %d   %d\n",
    p_Img->frame_no, pic_type, cur_bits, wp_method,
    p_Img->AverageFrameQP, lambda, 
    mPSNR->value[0], mPSNR->value[1], mPSNR->value[2],     
    tmp_time, (int) p_Img->me_time,
    p_Img->fld_flag ? "FLD" : "FRM", p_Img->intras, direct_mode,
    p_Img->num_ref_idx_l0_active, p_Img->num_ref_idx_l1_active, p_Img->rd_pass, p_Img->nal_reference_idc);
}

static void ReportVerboseNVB(ImageParameters *p_Img, char *pic_type, int cur_bits, int nvb_bits, int wp_method, int lambda, DistMetric *mPSNR, int tmp_time, int direct_mode)
{
  printf ("%05d(%3s)%8d %3d  %1d %2d %2d %7.3f %7.3f %7.3f %9d %7d    %3s %5d %1d %2d %2d  %d   %d\n",
    p_Img->frame_no, pic_type, cur_bits, nvb_bits, wp_method,
    p_Img->AverageFrameQP, lambda, 
    mPSNR->value[0], mPSNR->value[1], mPSNR->value[2],     
    tmp_time, (int) p_Img->me_time,
    p_Img->fld_flag ? "FLD" : "FRM", p_Img->intras, direct_mode,
    p_Img->num_ref_idx_l0_active, p_Img->num_ref_idx_l1_active, p_Img->rd_pass, p_Img->nal_reference_idc);
}

static void ReportVerboseFDN(ImageParameters *p_Img, char *pic_type, int cur_bits, int fdn_bits, int nvb_bits, int wp_method, int lambda, DistMetric *mPSNR, int tmp_time, int direct_mode)
{
  printf ("%05d(%3s)%8d %8d %3d  %1d %2d %2d %7.3f %7.3f %7.3f %9d %7d    %3s %5d %1d %2d %2d  %d   %d\n",
    p_Img->frame_no, pic_type, cur_bits, fdn_bits, nvb_bits, wp_method,
    p_Img->AverageFrameQP, lambda, 
    mPSNR->value[0], mPSNR->value[1], mPSNR->value[2],     
    tmp_time, (int) p_Img->me_time,
    p_Img->fld_flag ? "FLD" : "FRM", p_Img->intras, direct_mode,
    p_Img->num_ref_idx_l0_active, p_Img->num_ref_idx_l1_active, p_Img->rd_pass, p_Img->nal_reference_idc);
}

static void ReportVerboseSSIM(ImageParameters *p_Img, char *pic_type, int cur_bits, int wp_method, int lambda, DistMetric *mPSNR, DistMetric *mSSIM,int tmp_time, int direct_mode)
{
  printf ("%05d(%3s)%8d %1d %2d %2d %7.3f %7.3f %7.3f %7.4f %7.4f %7.4f %9d %7d    %3s %5d %1d %2d %2d  %d   %d\n",
    p_Img->frame_no, pic_type, cur_bits, wp_method,
    p_Img->AverageFrameQP, lambda, 
    mPSNR->value[0], mPSNR->value[1], mPSNR->value[2], 
    mSSIM->value[0], mSSIM->value[1], mSSIM->value[2], 
    tmp_time, (int) p_Img->me_time,
    p_Img->fld_flag ? "FLD" : "FRM", p_Img->intras, direct_mode,
    p_Img->num_ref_idx_l0_active, p_Img->num_ref_idx_l1_active,p_Img->rd_pass, p_Img->nal_reference_idc);
}

static void ReportVerboseNVBSSIM(ImageParameters *p_Img, char *pic_type, int cur_bits, int nvb_bits, int wp_method, int lambda, DistMetric *mPSNR, DistMetric *mSSIM,int tmp_time, int direct_mode)
{
  printf ("%05d(%3s)%8d %3d  %1d %2d %2d %7.3f %7.3f %7.3f %7.4f %7.4f %7.4f %9d %7d    %3s %5d %1d %2d %2d  %d   %d\n",
    p_Img->frame_no, pic_type, cur_bits, nvb_bits, wp_method,
    p_Img->AverageFrameQP, lambda, 
    mPSNR->value[0], mPSNR->value[1], mPSNR->value[2], 
    mSSIM->value[0], mSSIM->value[1], mSSIM->value[2], 
    tmp_time, (int) p_Img->me_time,
    p_Img->fld_flag ? "FLD" : "FRM", p_Img->intras, direct_mode,
    p_Img->num_ref_idx_l0_active, p_Img->num_ref_idx_l1_active, p_Img->rd_pass, p_Img->nal_reference_idc);
}

static void ReportVerboseFDNSSIM(ImageParameters *p_Img, char *pic_type, int cur_bits, int fdn_bits, int nvb_bits, int wp_method, int lambda, DistMetric *mPSNR, DistMetric *mSSIM,int tmp_time, int direct_mode)
{
  printf ("%05d(%3s)%8d %8d %3d  %1d %2d %2d %7.3f %7.3f %7.3f %7.4f %7.4f %7.4f %9d %7d    %3s %5d %1d %2d %2d  %d   %d\n",
    p_Img->frame_no, pic_type, cur_bits, fdn_bits, nvb_bits, wp_method,
    p_Img->AverageFrameQP, lambda, 
    mPSNR->value[0], mPSNR->value[1], mPSNR->value[2], 
    mSSIM->value[0], mSSIM->value[1], mSSIM->value[2], 
    tmp_time, (int) p_Img->me_time,
    p_Img->fld_flag ? "FLD" : "FRM", p_Img->intras, direct_mode,
    p_Img->num_ref_idx_l0_active, p_Img->num_ref_idx_l1_active, p_Img->rd_pass, p_Img->nal_reference_idc);
}


static void ReportNALNonVLCBits(ImageParameters *p_Img, InputParameters *p_Inp, StatParameters *p_Stats, int64 tmp_time)
{
  //! Need to add type (i.e. SPS, PPS, SEI etc).
  if (p_Inp->Verbose != 0)
    printf ("%05d(NVB)%8d \n", p_Img->frame_no, p_Stats->bit_ctr_parametersets_n);
}

static void ReportFirstframe(ImageParameters *p_Img, InputParameters *p_Inp, StatParameters *stats, int64 tmp_time)
{
  int cur_bits = (int)(stats->bit_ctr - stats->bit_ctr_n)
    + (int)(stats->bit_ctr_filler_data - stats->bit_ctr_filler_data_n);

  if (p_Inp->Verbose == 1)
  {
    ReportSimple(p_Img, "IDR", cur_bits, &p_Img->p_Dist->metric[PSNR], (int) tmp_time);
  }
  else if (p_Inp->Verbose == 2)
  {
    int lambda = (int) p_Img->lambda_me[I_SLICE][p_Img->AverageFrameQP][0];
    if (p_Inp->Distortion[SSIM] == 1)
      ReportVerboseSSIM(p_Img, "IDR", cur_bits, 0, lambda, &p_Img->p_Dist->metric[PSNR], &p_Img->p_Dist->metric[SSIM], (int) tmp_time, 0);
    else
      ReportVerbose(p_Img, "IDR", cur_bits, 0, lambda, &p_Img->p_Dist->metric[PSNR], (int) tmp_time, 0);
  }
  else if (p_Inp->Verbose == 3)
  {
    int lambda = (int) p_Img->lambda_me[I_SLICE][p_Img->AverageFrameQP][0];
    if (p_Inp->Distortion[SSIM] == 1)
      ReportVerboseNVBSSIM(p_Img, "IDR", cur_bits + stats->bit_ctr_parametersets_n, stats->bit_ctr_parametersets_n, 0, lambda, &p_Img->p_Dist->metric[PSNR], &p_Img->p_Dist->metric[SSIM], (int) tmp_time, 0);
    else
      ReportVerboseNVB(p_Img, "IDR", cur_bits + stats->bit_ctr_parametersets_n, stats->bit_ctr_parametersets_n, 0, lambda, &p_Img->p_Dist->metric[PSNR], (int) tmp_time, 0);
  }
  else if (p_Inp->Verbose == 4)
  {
    int lambda = (int) p_Img->lambda_me[I_SLICE][p_Img->AverageFrameQP][0];
    if (p_Inp->Distortion[SSIM] == 1)
      ReportVerboseFDNSSIM(p_Img, "IDR", cur_bits + stats->bit_ctr_parametersets_n, (int)(stats->bit_ctr_filler_data - stats->bit_ctr_filler_data_n), stats->bit_ctr_parametersets_n, 0, lambda, &p_Img->p_Dist->metric[PSNR], &p_Img->p_Dist->metric[SSIM], (int) tmp_time, 0);
    else
      ReportVerboseFDN(p_Img, "IDR", cur_bits + stats->bit_ctr_parametersets_n, (int)(stats->bit_ctr_filler_data - stats->bit_ctr_filler_data_n), stats->bit_ctr_parametersets_n, 0, lambda, &p_Img->p_Dist->metric[PSNR], (int) tmp_time, 0);
  }

  stats->bit_counter[I_SLICE] = stats->bit_ctr;
  stats->bit_ctr = 0;
}

static void ReportI(ImageParameters *p_Img, InputParameters *p_Inp, StatParameters *stats, int64 tmp_time)
{
  char pic_type[4];
  int  cur_bits = (int)(stats->bit_ctr - stats->bit_ctr_n)
    + (int)(stats->bit_ctr_filler_data - stats->bit_ctr_filler_data_n);

  if ((p_Inp->redundant_pic_flag == 0) || !p_Img->redundant_coding )
  {
    if (p_Img->currentPicture->idr_flag == TRUE)
      strcpy(pic_type,"IDR");
    else
      strcpy(pic_type," I ");
  }
  else
    strcpy(pic_type,"R");

  if (p_Inp->Verbose == 1)
  {
    ReportSimple(p_Img, pic_type, cur_bits, &p_Img->p_Dist->metric[PSNR], (int) tmp_time);
  }
  else if (p_Inp->Verbose == 2)
  {
    int lambda = (int) p_Img->lambda_me[I_SLICE][p_Img->AverageFrameQP][0];
    if (p_Inp->Distortion[SSIM] == 1)
    {
      ReportVerboseSSIM(p_Img, pic_type, cur_bits, 0, lambda, &p_Img->p_Dist->metric[PSNR], &p_Img->p_Dist->metric[SSIM], (int) tmp_time, 0);
    }
    else
    {
      ReportVerbose(p_Img, pic_type, cur_bits, 0, lambda, &p_Img->p_Dist->metric[PSNR], (int) tmp_time, 0);
    }
  }
  else if (p_Inp->Verbose == 3)
  {
    int lambda = (int) p_Img->lambda_me[I_SLICE][p_Img->AverageFrameQP][0];
    if (p_Inp->Distortion[SSIM] == 1)
    {
      ReportVerboseNVBSSIM(p_Img, pic_type, cur_bits + stats->bit_ctr_parametersets_n, stats->bit_ctr_parametersets_n, 0, lambda, &p_Img->p_Dist->metric[PSNR], &p_Img->p_Dist->metric[SSIM], (int) tmp_time, 0);
    }
    else
    {
      ReportVerboseNVB(p_Img, pic_type, cur_bits + stats->bit_ctr_parametersets_n, stats->bit_ctr_parametersets_n, 0, lambda, &p_Img->p_Dist->metric[PSNR], (int) tmp_time, 0);
    }
  }
  else if (p_Inp->Verbose == 4)
  {
    int lambda = (int) p_Img->lambda_me[I_SLICE][p_Img->AverageFrameQP][0];
    if (p_Inp->Distortion[SSIM] == 1)
    {
      ReportVerboseFDNSSIM(p_Img, pic_type, cur_bits + stats->bit_ctr_parametersets_n, (int)(stats->bit_ctr_filler_data - stats->bit_ctr_filler_data_n), stats->bit_ctr_parametersets_n, 0, lambda, &p_Img->p_Dist->metric[PSNR], &p_Img->p_Dist->metric[SSIM], (int) tmp_time, 0);
    }
    else
    {
      ReportVerboseFDN(p_Img, pic_type, cur_bits + stats->bit_ctr_parametersets_n, (int)(stats->bit_ctr_filler_data - stats->bit_ctr_filler_data_n), stats->bit_ctr_parametersets_n, 0, lambda, &p_Img->p_Dist->metric[PSNR], (int) tmp_time, 0);
    }
  }
}

static void ReportB(ImageParameters *p_Img, InputParameters *p_Inp, StatParameters *stats, int64 tmp_time)
{
  int cur_bits = (int)(stats->bit_ctr - stats->bit_ctr_n)
    + (int)(stats->bit_ctr_filler_data - stats->bit_ctr_filler_data_n);

  if (p_Inp->Verbose == 1)
  {
    ReportSimple(p_Img, " B ", cur_bits, &p_Img->p_Dist->metric[PSNR], (int) tmp_time);
  }
  else if (p_Inp->Verbose == 2)
  {
    int lambda = (int) p_Img->lambda_me[p_Img->nal_reference_idc ? 5 : B_SLICE][p_Img->AverageFrameQP][0];    
    if (p_Inp->Distortion[SSIM] == 1)
      ReportVerboseSSIM(p_Img, " B ", cur_bits, p_Img->active_pps->weighted_bipred_idc, lambda, &p_Img->p_Dist->metric[PSNR], &p_Img->p_Dist->metric[SSIM], (int) tmp_time, p_Img->direct_spatial_mv_pred_flag);
    else
      ReportVerbose(p_Img, " B ", cur_bits, p_Img->active_pps->weighted_bipred_idc, lambda, &p_Img->p_Dist->metric[PSNR], (int) tmp_time, p_Img->direct_spatial_mv_pred_flag);
  }
  else if (p_Inp->Verbose == 3)
  {
    int lambda = (int) p_Img->lambda_me[p_Img->nal_reference_idc ? 5 : B_SLICE][p_Img->AverageFrameQP][0];    
    if (p_Inp->Distortion[SSIM] == 1)
      ReportVerboseNVBSSIM(p_Img, " B ", cur_bits + stats->bit_ctr_parametersets_n, stats->bit_ctr_parametersets_n, p_Img->active_pps->weighted_bipred_idc, lambda, &p_Img->p_Dist->metric[PSNR], &p_Img->p_Dist->metric[SSIM], (int) tmp_time, p_Img->direct_spatial_mv_pred_flag);
    else
      ReportVerboseNVB(p_Img, " B ", cur_bits + stats->bit_ctr_parametersets_n, stats->bit_ctr_parametersets_n, p_Img->active_pps->weighted_bipred_idc, lambda, &p_Img->p_Dist->metric[PSNR], (int) tmp_time, p_Img->direct_spatial_mv_pred_flag);
  }
  else if (p_Inp->Verbose == 4)
  {
    int lambda = (int) p_Img->lambda_me[p_Img->nal_reference_idc ? 5 : B_SLICE][p_Img->AverageFrameQP][0];    
    if (p_Inp->Distortion[SSIM] == 1)
      ReportVerboseFDNSSIM(p_Img, " B ", cur_bits + stats->bit_ctr_parametersets_n, (int)(stats->bit_ctr_filler_data - stats->bit_ctr_filler_data_n), stats->bit_ctr_parametersets_n, p_Img->active_pps->weighted_bipred_idc, lambda, &p_Img->p_Dist->metric[PSNR], &p_Img->p_Dist->metric[SSIM], (int) tmp_time, p_Img->direct_spatial_mv_pred_flag);
    else
      ReportVerboseFDN(p_Img, " B ", cur_bits + stats->bit_ctr_parametersets_n, (int)(stats->bit_ctr_filler_data - stats->bit_ctr_filler_data_n), stats->bit_ctr_parametersets_n, p_Img->active_pps->weighted_bipred_idc, lambda, &p_Img->p_Dist->metric[PSNR], (int) tmp_time, p_Img->direct_spatial_mv_pred_flag);
  }
}

static void ReportP(ImageParameters *p_Img, InputParameters *p_Inp, StatParameters *stats, int64 tmp_time)
{
  char pic_type[4];
  int  cur_bits = (int)(stats->bit_ctr - stats->bit_ctr_n)
    + (int)(stats->bit_ctr_filler_data - stats->bit_ctr_filler_data_n);

  if (p_Img->type == SP_SLICE)
    strcpy(pic_type,"SP ");
  else if ((p_Inp->redundant_pic_flag == 0) || !p_Img->redundant_coding )
    strcpy(pic_type," P ");
  else
    strcpy(pic_type," R ");

  if (p_Inp->Verbose == 1)
  {
    ReportSimple(p_Img, pic_type, cur_bits, &p_Img->p_Dist->metric[PSNR], (int) tmp_time);
  }
  else if (p_Inp->Verbose == 2)
  {
    int lambda = (int) p_Img->lambda_me[P_SLICE][p_Img->AverageFrameQP][0];    
    if (p_Inp->Distortion[SSIM] == 1)
      ReportVerboseSSIM(p_Img, pic_type, cur_bits, p_Img->active_pps->weighted_pred_flag, lambda, &p_Img->p_Dist->metric[PSNR], &p_Img->p_Dist->metric[SSIM], (int) tmp_time, 0);
    else
      ReportVerbose(p_Img, pic_type, cur_bits, p_Img->active_pps->weighted_pred_flag, lambda, &p_Img->p_Dist->metric[PSNR], (int) tmp_time, 0);
  }
  else if (p_Inp->Verbose == 3)
  {
    int lambda = (int) p_Img->lambda_me[P_SLICE][p_Img->AverageFrameQP][0];    
    if (p_Inp->Distortion[SSIM] == 1)
      ReportVerboseNVBSSIM(p_Img, pic_type, cur_bits + stats->bit_ctr_parametersets_n, stats->bit_ctr_parametersets_n, p_Img->active_pps->weighted_pred_flag, lambda, &p_Img->p_Dist->metric[PSNR], &p_Img->p_Dist->metric[SSIM], (int) tmp_time, 0);
    else
      ReportVerboseNVB(p_Img, pic_type, cur_bits + stats->bit_ctr_parametersets_n, stats->bit_ctr_parametersets_n, p_Img->active_pps->weighted_pred_flag, lambda, &p_Img->p_Dist->metric[PSNR], (int) tmp_time, 0);
  }
  else if (p_Inp->Verbose == 4)
  {
    int lambda = (int) p_Img->lambda_me[P_SLICE][p_Img->AverageFrameQP][0];    
    if (p_Inp->Distortion[SSIM] == 1)
      ReportVerboseFDNSSIM(p_Img, pic_type, cur_bits + stats->bit_ctr_parametersets_n, (int)(stats->bit_ctr_filler_data - stats->bit_ctr_filler_data_n), stats->bit_ctr_parametersets_n, p_Img->active_pps->weighted_pred_flag, lambda, &p_Img->p_Dist->metric[PSNR], &p_Img->p_Dist->metric[SSIM], (int) tmp_time, 0);
    else
      ReportVerboseFDN(p_Img, pic_type, cur_bits + stats->bit_ctr_parametersets_n, (int)(stats->bit_ctr_filler_data - stats->bit_ctr_filler_data_n), stats->bit_ctr_parametersets_n, p_Img->active_pps->weighted_pred_flag, lambda, &p_Img->p_Dist->metric[PSNR], (int) tmp_time, 0);
  }
}

/*!
 ************************************************************************
 * \brief
 *    Calculates the absolute frame number in the source file out
 *    of various variables in p_Img-> and p_Inp->
 * \return
 *    frame number in the file to be read
 * \par side effects
 *    global variable frame_no updated -- dunno, for what this one is necessary
 ************************************************************************
 */
int CalculateFrameNumber(ImageParameters *p_Img, InputParameters *p_Inp)
{
  int frm_sign = (p_Img->gop_number && (p_Img->gop_number <= p_Inp->intra_delay)) ? -1 : 1;
  int delay    = (p_Img->gop_number <= p_Inp->intra_delay) ? p_Inp->intra_delay : 0;
  
  if (p_Img->b_frame_to_code)
  {
    if ((p_Img->gop_number && (p_Img->gop_number <= p_Inp->intra_delay)))
    {
      if (p_Inp->HierarchicalCoding)
        p_Img->frame_no = p_Img->start_tr_gop + (p_Inp->intra_delay - p_Img->gop_number) * p_Img->base_dist 
        +  (int) (p_Img->frame_interval * (double) (1 + p_Img->gop_structure[p_Img->b_frame_to_code - 1].display_no));
      else
        p_Img->frame_no = p_Img->start_tr_gop + (p_Inp->intra_delay - p_Img->gop_number) * p_Img->base_dist 
        + (int) (p_Img->frame_interval * (double) p_Img->b_frame_to_code);
    }
    else
    {
      if (p_Inp->HierarchicalCoding)
        p_Img->frame_no = p_Img->start_tr_gop + (p_Img->gop_number - 1) * p_Img->base_dist 
        + (int) (p_Img->frame_interval * (double) (1 + p_Img->gop_structure[p_Img->b_frame_to_code - 1].display_no));
      else
        p_Img->frame_no = p_Img->start_tr_gop + (p_Img->gop_number - 1) * p_Img->base_dist 
        + (int) (p_Img->frame_interval * (double) p_Img->b_frame_to_code);
    }
    //printf("frame_no %d %d %d\n",p_Img->frame_no,p_Img->gop_number - 1,(frm_sign * (p_Img->gop_number - 1) + p_Inp->intra_delay));
  }
  else
  {
    if (p_Inp->idr_period && p_Inp->EnableIDRGOP && p_Img->frm_number && 
      ((!p_Inp->adaptive_idr_period && (p_Img->frm_number - p_Img->lastIDRnumber ) % p_Inp->idr_period == 0)
      || (p_Inp->adaptive_idr_period == 1 && (p_Img->frm_number - imax(p_Img->lastIntraNumber, p_Img->lastIDRnumber) ) % p_Inp->idr_period == 0)) )
    {
      delay = p_Inp->intra_delay;
      p_Img->rewind_frame += p_Inp->NumberBFrames;
      p_Img->start_frame_no = p_Img->frm_number;
      p_Img->gop_number = (p_Img->number - p_Img->start_frame_no);
      p_Img->start_tr_gop = (p_Img->frm_number) * p_Img->base_dist - p_Img->rewind_frame;
    }

    p_Img->frame_no = p_Img->start_tr_gop + (frm_sign * p_Img->gop_number + delay) * p_Img->base_dist;

    if (p_Img->frame_no > p_Inp->last_frame)
      p_Img->frame_no = imin(p_Img->frame_no, (p_Inp->no_frames - 1) * (p_Inp->frame_skip + 1));
    //printf("%d %d %d %d %d\n", p_Img->frame_no, p_Img->start_tr_gop, p_Img->gop_number, delay, p_Inp->last_frame);
  }

  return p_Img->frame_no;
}

/*!
 ************************************************************************
 * \brief
 *    point to frame coding variables
 ************************************************************************
 */
static void put_buffer_frame(ImageParameters *p_Img)
{
  p_Img->pCurImg    = p_Img->imgData.frm_data[0];
  p_Img->pImgOrg[0] = p_Img->imgData.frm_data[0];    
  
  if (p_Img->yuv_format != YUV400)
  {
    p_Img->pImgOrg[1] = p_Img->imgData.frm_data[1];
    p_Img->pImgOrg[2] = p_Img->imgData.frm_data[2];
  }
}

/*!
 ************************************************************************
 * \brief
 *    point to top field coding variables
 ************************************************************************
 */
static void put_buffer_top(ImageParameters *p_Img)
{
  p_Img->fld_type = 0;

  p_Img->pCurImg    = p_Img->imgData.top_data[0];
  p_Img->pImgOrg[0] = p_Img->imgData.top_data[0];
  
  if (p_Img->yuv_format != YUV400)
  {
    p_Img->pImgOrg[1] = p_Img->imgData.top_data[1];
    p_Img->pImgOrg[2] = p_Img->imgData.top_data[2];
  }
}

/*!
 ************************************************************************
 * \brief
 *    point to bottom field coding variables
 ************************************************************************
 */
static void put_buffer_bot(ImageParameters *p_Img)
{
  p_Img->fld_type = 1;

  p_Img->pCurImg    = p_Img->imgData.bot_data[0];  
  p_Img->pImgOrg[0] = p_Img->imgData.bot_data[0];
  
  if (p_Img->yuv_format != YUV400)
  {
    p_Img->pImgOrg[1] = p_Img->imgData.bot_data[1];
    p_Img->pImgOrg[2] = p_Img->imgData.bot_data[2];
  }
}

/*!
*************************************************************************************
* Brief
*     Output SP frames coefficients
*************************************************************************************
*/
void output_SP_coefficients(ImageParameters *p_Img, InputParameters *p_Inp)
{
  int i,k;
  FILE *SP_coeff_file;
  int ret;
  if(p_Img->number_sp2_frames==0)
  {
    if ((SP_coeff_file = fopen(p_Inp->sp_output_filename,"wb")) == NULL)
    {
      printf ("Fatal: cannot open SP output file '%s', exit (-1)\n", p_Inp->sp_output_filename);
      exit (-1);
    }
    p_Img->number_sp2_frames++;
  }
  else
  {
    if ((SP_coeff_file = fopen(p_Inp->sp_output_filename,"ab")) == NULL)
    {
      printf ("Fatal: cannot open SP output file '%s', exit (-1)\n", p_Inp->sp_output_filename);
      exit (-1);
    }
  }

  for(i=0;i<p_Img->height;i++)
  {
    ret = fwrite(p_Img->lrec[i],sizeof(int),p_Img->width,SP_coeff_file);
    if (ret != p_Img->width)
    {
      error ("cannot write to SP output file", -1);
    }
  }

  for(k=0;k<2;k++)
  {
    for(i=0;i<p_Img->height_cr;i++)
    {
      ret = fwrite(p_Img->lrec_uv[k][i],sizeof(int),p_Img->width_cr,SP_coeff_file);
      if (ret != p_Img->width_cr)
      {
        error ("cannot write to SP output file", -1);
      }
    }
  }
  fclose(SP_coeff_file);
}

/*!
*************************************************************************************
* Brief
*     Read SP frames coefficients
*************************************************************************************
*/
void read_SP_coefficients(ImageParameters *p_Img, InputParameters *p_Inp)
{
  int i,k;
  FILE *SP_coeff_file;

  if ( (p_Inp->qp2start > 0) && ( ( (p_Img->frame_no ) % (2*p_Inp->qp2start) ) >=p_Inp->qp2start ))
  {
    if ((SP_coeff_file = fopen(p_Inp->sp2_input_filename1,"rb")) == NULL)
    {
      printf ("Fatal: cannot open SP input file '%s', exit (-1)\n", p_Inp->sp2_input_filename2);
      exit (-1);
    }
  }
  else
  {
    if ((SP_coeff_file = fopen(p_Inp->sp2_input_filename2,"rb")) == NULL)
    {
      printf ("Fatal: cannot open SP input file '%s', exit (-1)\n", p_Inp->sp2_input_filename1);
      exit (-1);
    }
  }

  if (0 != fseek (SP_coeff_file, p_Img->size * 3/2*p_Img->number_sp2_frames*sizeof(int), SEEK_SET))
  {
    printf ("Fatal: cannot seek in SP input file, exit (-1)\n");
    exit (-1);
  }
  p_Img->number_sp2_frames++;

  for(i=0;i<p_Img->height;i++)
  {
    if(p_Img->width!=(int)fread(p_Img->lrec[i],sizeof(int),p_Img->width,SP_coeff_file))
    {
      printf ("Fatal: cannot read in SP input file, exit (-1)\n");
      exit (-1);
    }
  }

  for(k=0;k<2;k++)
  {
    for(i=0;i<p_Img->height_cr;i++)
    {
      if(p_Img->width_cr!=(int)fread(p_Img->lrec_uv[k][i],sizeof(int),p_Img->width_cr,SP_coeff_file))
      {
        printf ("Fatal: cannot read in SP input file, exit (-1)\n");
        exit (-1);
      }
    }
  }
  fclose(SP_coeff_file);
}


/*!
*************************************************************************************
* Brief
*     Select appropriate image plane (for 444 coding)
*************************************************************************************
*/
void select_plane(ImageParameters *p_Img, ColorPlane color_plane)
{
  p_Img->pCurImg              = p_Img->pImgOrg[color_plane];
  p_Img->enc_picture->p_curr_img     = p_Img->enc_picture->p_img[color_plane];
  p_Img->enc_picture->p_curr_img_sub = p_Img->enc_picture->p_img_sub[color_plane];
  p_Img->max_imgpel_value     = (short) p_Img->max_pel_value_comp[color_plane];
  p_Img->dc_pred_value        = p_Img->dc_pred_value_comp[color_plane];
}

/*!
*************************************************************************************
* Brief
*     Is this picture the first access unit in this GOP?
*************************************************************************************
*/
static int is_gop_first_unit(ImageParameters *p_Img, InputParameters *p_Inp)
{
  return ( get_idr_flag(p_Img, p_Inp) || ( p_Img->type == I_SLICE && p_Inp->EnableOpenGOP ) );
}

/*!
*************************************************************************************
* Brief
*     AUD, SPS, PPS, and SEI messages
*************************************************************************************
*/
void write_non_vcl_nalu( ImageParameters *p_Img, InputParameters *p_Inp )
{
  // SPS + PPS
  if (p_Inp->ResendSPS == 3 && is_gop_first_unit(p_Img, p_Inp) && p_Img->number)
  {
    p_Img->p_Stats->bit_slice = rewrite_paramsets(p_Img, p_Inp);
  }

  if (p_Inp->ResendSPS == 2 && get_idr_flag(p_Img, p_Inp) && p_Img->number)
  {
    p_Img->p_Stats->bit_slice = rewrite_paramsets(p_Img, p_Inp);
  }
  if (p_Inp->ResendSPS == 1 && p_Img->type == I_SLICE && p_Img->frm_number != 0)
  {
    p_Img->p_Stats->bit_slice = rewrite_paramsets(p_Img, p_Inp);
  }
  // PPS
  if ( p_Inp->ResendPPS && p_Img->frm_number != 0
    && (p_Img->type != I_SLICE || p_Inp->ResendSPS != 1) 
    && (!(p_Inp->ResendSPS == 2 && get_idr_flag(p_Img, p_Inp)))
    && (!(p_Inp->ResendSPS == 3 && is_gop_first_unit(p_Img, p_Inp))) )
  {
    // Access Unit Delimiter NALU
    if ( p_Inp->SendAUD )
    {
      p_Img->p_Stats->bit_ctr_parametersets_n = Write_AUD_NALU(p_Img);
      p_Img->p_Stats->bit_ctr_parametersets_n += write_PPS(p_Img, p_Inp, 0, 0);
    }
    else
    {
      p_Img->p_Stats->bit_ctr_parametersets_n = write_PPS(p_Img, p_Inp, 0, 0);
    }
  }
  // Access Unit Delimiter NALU
  if ( p_Inp->SendAUD
    && (!(p_Inp->ResendPPS && p_Img->frm_number != 0))
    && (p_Img->type != I_SLICE || p_Inp->ResendSPS != 1) 
    && (!(p_Inp->ResendSPS == 2 && get_idr_flag(p_Img, p_Inp)))
    && (!(p_Inp->ResendSPS == 3 && is_gop_first_unit(p_Img, p_Inp))) )
  {
    p_Img->p_Stats->bit_ctr_parametersets_n += Write_AUD_NALU(p_Img);
  }

  UpdateSubseqInfo (p_Img, p_Inp, p_Img->layer);        // Tian Dong (Sept 2002)
  UpdateSceneInformation (p_Img->p_SEI, FALSE, 0, 0, -1); // JVT-D099, scene information SEI, nothing included by default

  //! Commented out by StW, needs fixing in SEI.h to keep the trace file clean
  //  PrepareAggregationSEIMessage (p_Img);

  // write tone mapping SEI message
  if (p_Inp->ToneMappingSEIPresentFlag)
  {
    UpdateToneMapping(p_Img->p_SEI);
  }

  PrepareAggregationSEIMessage(p_Img);
  
  p_Img->p_Stats->bit_ctr_parametersets_n += Write_SEI_NALU(p_Img, 0);
  // update seq NVB counter
  p_Img->p_Stats->bit_ctr_parametersets   += p_Img->p_Stats->bit_ctr_parametersets_n;
}


