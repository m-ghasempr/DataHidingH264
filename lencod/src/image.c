
/*!
 *************************************************************************************
 * \file image.c
 *
 * \brief
 *    Code one image/slice
 *
 * \author
 *    Main contributors (see contributors.h for copyright, address and affiliation details)
 *     - Inge Lille-Langøy               <inge.lille-langoy@telenor.com>
 *     - Rickard Sjoberg                 <rickard.sjoberg@era.ericsson.se>
 *     - Jani Lainema                    <jani.lainema@nokia.com>
 *     - Sebastian Purreiter             <sebastian.purreiter@mch.siemens.de>
 *     - Byeong-Moon Jeon                <jeonbm@lge.com>
 *     - Yoon-Seong Soh                  <yunsung@lge.com>
 *     - Thomas Stockhammer              <stockhammer@ei.tum.de>
 *     - Detlev Marpe                    <marpe@hhi.de>
 *     - Guido Heising                   <heising@hhi.de>
 *     - Thomas Wedi                     <wedi@tnt.uni-hannover.de>
 *     - Ragip Kurceren                  <ragip.kurceren@nokia.com>
 *     - Antti Hallapuro                 <antti.hallapuro@nokia.com>
 *************************************************************************************
 */
#include "contributors.h"

#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/timeb.h>
#include <string.h>
#include <memory.h>
#include <assert.h>

#include "global.h"
#include "image.h"
#include "refbuf.h"
#include "mbuffer.h"
#include "encodeiff.h"
#include "header.h"
#include "intrarefresh.h"
#include "fmo.h"
#include "sei.h"
#include "memalloc.h"
#include "nalu.h"

#define TO_SAVE 4711
#define FROM_SAVE 4712
static void copy_mv_to_or_from_save (int direction);

void code_a_picture(Picture *pic);
void frame_picture (Picture *frame);
void field_picture(Picture *top, Picture *bottom);

static int  writeout_picture(Picture *pic);

static int  picture_structure_decision(Picture *frame, Picture *top, Picture *bot);
static void distortion_fld (float *dis_fld_y, float *dis_fld_u, float *dis_fld_v);
static void find_snr();
static void find_distortion();

static void field_mode_buffer(int bit_field, float snr_field_y, float snr_field_u, float snr_field_v);
static void frame_mode_buffer (int bit_frame, float snr_frame_y, float snr_frame_u, float snr_frame_v);

static void init_frame();
static void init_field();

static void put_buffer_frame();
static void put_buffer_top();
static void put_buffer_bot();

static void interpolate_frame();
static void interpolate_frame_to_fb();

static void estimate_weighting_factor ();

static void store_field_MV(int frame_number);
static void store_direct_moving_flag (int frame_number);
static void copy_motion_vectors_MB(int bot_block); //!< For MB level field/frame coding tools

static void CopyFrameToOldImgOrgVariables (Sourceframe *sf);
static void CopyTopFieldToOldImgOrgVariables (Sourceframe *sf);
static void CopyBottomFieldToOldImgOrgVariables (Sourceframe *sf);
static Sourceframe *AllocSourceframe (int xs, int ys);
static void FreeSourceframe (Sourceframe *sf);
static void ReadOneFrame (int FrameNoInFile, int HeaderSize, int xs, int ys, Sourceframe *sf);
static void writeUnit(Bitstream* currStream);

#ifdef _ADAPT_LAST_GROUP_
int *last_P_no;
int *last_P_no_frm;
int *last_P_no_fld;
#endif

static void ReportFirstframe(int tmp_time);
static void ReportIntra(int tmp_time);
static void ReportSP(int tmp_time);
static void ReportBS(int tmp_time);
static void ReportP(int tmp_time);
static void ReportB(int tmp_time);

static int CalculateFrameNumber();  // Calculates the next frame number
static int FrameNumberInFile;       // The current frame number in the input file

static Sourceframe *srcframe;

StorablePicture *enc_picture;
StorablePicture *enc_frame_picture;
StorablePicture *enc_top_picture;
StorablePicture *enc_bottom_picture;


const int ONE_FOURTH_TAP[3][2] =
{
  {20,20},
  {-5,-4},
  { 1, 0},
};


/*!
 ************************************************************************
 * \brief
 *    Encodes a picture
 *
 *    This is the main picture coding loop.. It is called by all this
 *    frame and field coding stuff after the img-> elements have been
 *    set up.  Not sure whether it is useful for MB-adaptive frame/field
 *    coding
 ************************************************************************
 */
void code_a_picture(Picture *pic)
{
  int NumberOfCodedMBs = 0;
  int SliceGroup = 0;
  int i,j;

  img->currentPicture = pic;

  img->currentPicture->idr_flag = (!IMG_NUMBER) && (!(img->structure==BOTTOM_FIELD));

  pic->no_slices = 0;
  pic->distortion_u = pic->distortion_v = pic->distortion_y = 0.0;

  // restrict list 1 size
  if (img->structure==FRAME)
  {
    img->num_ref_idx_l0_active = img->buf_cycle;
    img->num_ref_idx_l1_active = (img->type==B_SLICE?1:0);
  }
  else
  {
    img->num_ref_idx_l0_active = 2*img->buf_cycle;
    img->num_ref_idx_l1_active = (img->type==B_SLICE?2:0);
  }

  // generate reference picture lists
  init_lists(img->type, img->structure);

	// update reference picture number index
	for (i=0;i<listXsize[LIST_0];i++)
  {
    enc_picture->ref_pic_num[LIST_0][i]=listX[LIST_0][i]->poc;
  }

  for (i=0;i<listXsize[LIST_1];i++)
  {
    enc_picture->ref_pic_num[LIST_1][i]=listX[LIST_1][i]->poc;
  }

  if (img->MbaffFrameFlag)
    for (j=2;j<6;j++)
      for (i=0;i<listXsize[j];i++)
      {
        enc_picture->ref_pic_num[j][i]=listX[j][i]->poc;        
      }

  // assign list 0 size from list size
  img->num_ref_idx_l0_active = listXsize[0];
  img->num_ref_idx_l1_active = listXsize[1];

  if (img->MbaffFrameFlag)
    init_mbaff_lists();

  RandomIntraNewPicture ();     //! Allocates forced INTRA MBs (even for fields!)
  FmoStartPicture ();           //! picture level initialization of FMO

  while (NumberOfCodedMBs < img->total_number_mb)       // loop over slices
  {
    // Encode one SLice Group
    while (!FmoSliceGroupCompletelyCoded (SliceGroup))
    {
      // Encode the current slice
      NumberOfCodedMBs += encode_one_slice (SliceGroup, pic);
      FmoSetLastMacroblockInSlice (img->current_mb_nr);
      // Proceed to next slice
      img->current_slice_nr++;
      stat->bit_slice = 0;
    }
    // Proceed to next SliceGroup
    SliceGroup++;
  }
  FmoEndPicture ();
  if (input->rdopt == 2 && (img->type != B_SLICE))
    for (j = 0; j < input->NoOfDecoders; j++)
      DeblockFrame (img, decs->decY_best[j], NULL);
  DeblockFrame (img, enc_picture->imgY, enc_picture->imgUV);
}



/*!
 ************************************************************************
 * \brief
 *    Encodes one frame
 ************************************************************************
 */
int encode_one_frame ()
{
  static int prev_frame_no = 0; // POC200301
  static int consecutive_non_reference_pictures = 0; // POC200301

#ifdef _LEAKYBUCKET_
  extern long Bit_Buffer[10000];
  extern unsigned long total_frame_buffer;
#endif

  time_t ltime1;
  time_t ltime2;

#ifdef WIN32
  struct _timeb tstruct1;
  struct _timeb tstruct2;
#else
  struct timeb tstruct1;
  struct timeb tstruct2;
#endif

  int tmp_time;
  int bits_frm = 0, bits_fld = 0;
  float dis_frm = 0, dis_frm_y = 0, dis_frm_u = 0, dis_frm_v = 0;
  float dis_fld = 0, dis_fld_y = 0, dis_fld_u = 0, dis_fld_v = 0;

#ifdef WIN32
  _ftime (&tstruct1);           // start time ms
#else
  ftime (&tstruct1);
#endif
  time (&ltime1);               // start time s

/*
  //Shankar Regunathan (Oct 2002)
  //Prepare Panscanrect SEI payload
  UpdatePanScanRectInfo ();
  //Prepare Arbitrarydata SEI Payload
  UpdateUser_data_unregistered ();
  //Prepare Registered data SEI Payload
  UpdateUser_data_registered_itu_t_t35 ();
  //Prepare RandomAccess SEI Payload
  UpdateRandomAccess ();
*/

  put_buffer_frame ();      // sets the pointers to the frame structures 
                            // (and not to one of the field structures)
  init_frame ();
  FrameNumberInFile = CalculateFrameNumber();

  srcframe = AllocSourceframe (img->width, img->height);
  ReadOneFrame (FrameNumberInFile, input->infile_header, img->width, img->height, srcframe);
  CopyFrameToOldImgOrgVariables (srcframe);

  if (img->type == B_SLICE)
    Bframe_ctr++;         // Bframe_ctr only used for statistics, should go to stat->

  if (input->InterlaceCodingOption == FIELD_CODING)
  {
    img->field_picture = 1;  // we encode fields
    field_picture (top_pic, bottom_pic);
    img->fld_flag = 1;
  }
  else
  {
    // For frame coding, turn MB level field/frame coding flag on
    if (input->InterlaceCodingOption >= MB_CODING)
      mb_adaptive = 1;
    img->field_picture = 0; // we encode a frame

    frame_picture (frame_pic);
    // For field coding, turn MB level field/frame coding flag off
    if (input->InterlaceCodingOption >= MB_CODING)
      mb_adaptive = 0;
    
    if (input->InterlaceCodingOption != FRAME_CODING)
    {
      img->field_picture = 1;  // we encode fields
      field_picture (top_pic, bottom_pic);
      
      //! Note: the distortion for a field coded picture is stored in the top field
      //! the distortion values in the bottom field are dummies
      dis_fld = top_pic->distortion_y + top_pic->distortion_u + top_pic->distortion_v;
      dis_frm = frame_pic->distortion_y + frame_pic->distortion_u + frame_pic->distortion_v;
      
      img->fld_flag = picture_structure_decision (frame_pic, top_pic, bottom_pic);
      update_field_frame_contexts (img->fld_flag);
    }
    else
      img->fld_flag = 0;
  }

  if (img->fld_flag)
    stat->bit_ctr_emulationprevention += stat->em_prev_bits_fld;
  else
    stat->bit_ctr_emulationprevention += stat->em_prev_bits_frm;

  if (img->type != B_SLICE)
  {
    img->pstruct_next_P = img->fld_flag;
  }

  // Here, img->structure may be either FRAME or BOTTOM FIELD depending on whether AFF coding is used
  // The picture structure decision changes really only the fld_flag

  if (img->fld_flag)            // field mode (use field when fld_flag=1 only)
  {
    field_mode_buffer (bits_fld, dis_fld_y, dis_fld_u, dis_fld_v);
    writeout_picture (top_pic);
    writeout_picture (bottom_pic);
  }
  else                          //frame mode
  {
    frame_mode_buffer (bits_frm, dis_frm_y, dis_frm_u, dis_frm_v);
    writeout_picture (frame_pic);
  }

  if (frame_pic)
    free_slice_list(frame_pic);
  if (top_pic)
    free_slice_list(top_pic);
  if (bottom_pic)
    free_slice_list(bottom_pic);

  /*
  // Tian Dong (Sept 2002)
  // in frame mode, the newly reconstructed frame has been inserted to the mem buffer
  // and it is time to prepare the spare picture SEI payload.
  if (input->InterlaceCodingOption == FRAME_CODING
      && input->SparePictureOption && img->type != B_SLICE)
    CalculateSparePicture ();
*/

  if (input->InterlaceCodingOption != FRAME_CODING)
    {
      store_field_MV (IMG_NUMBER);      // assume that img->number = frame_number
    }
  else
    store_direct_moving_flag (IMG_NUMBER);

  if (input->InterlaceCodingOption >= MB_CODING && img->number != 0 && input->successive_Bframe != 0 && img->type != B_SLICE)
  {
    copy_mv_to_or_from_save (FROM_SAVE);
  }
  if (input->InterlaceCodingOption == FRAME_CODING)
    {
      if (input->rdopt == 2 && img->type != B_SLICE)
        UpdateDecoders ();      // simulate packet losses and move decoded image to reference buffers

      if (input->RestrictRef)
        UpdatePixelMap ();
    }

  find_snr ();

  time (&ltime2);               // end time sec
#ifdef WIN32
  _ftime (&tstruct2);           // end time ms
#else
  ftime (&tstruct2);            // end time ms
#endif

  tmp_time = (ltime2 * 1000 + tstruct2.millitm) - (ltime1 * 1000 + tstruct1.millitm);
  tot_time = tot_time + tmp_time;

  // Write reconstructed images
  if (img->fld_flag)
  {
    store_picture_in_dpb(enc_top_picture);
    store_picture_in_dpb(enc_bottom_picture);
    enc_picture = enc_top_picture = enc_bottom_picture = NULL;
    if (enc_frame_picture)
    {
      free_storable_picture(enc_frame_picture);
      enc_frame_picture = NULL;
    }
  }
  else
  {
    store_picture_in_dpb(enc_frame_picture);
    enc_picture = enc_frame_picture = NULL;
    if (enc_top_picture)
    {
      free_storable_picture(enc_top_picture);
      enc_top_picture = NULL;
    }
    if (enc_bottom_picture)
    {
      free_storable_picture(enc_bottom_picture);
      enc_bottom_picture = NULL;
    }
  }


#ifdef _LEAKYBUCKET_
  // Store bits used for this frame and increment counter of no. of coded frames
  Bit_Buffer[total_frame_buffer] = stat->bit_ctr - stat->bit_ctr_n;
  total_frame_buffer++;
#endif

  // POC200301: Verify that POC coding type 2 is not used if more than one consecutive 
  // non-reference frame is requested or if decoding order is different from output order
  if (img->pic_order_cnt_type == 2)
  {
    if (!img->nal_reference_idc) consecutive_non_reference_pictures++;
    else consecutive_non_reference_pictures = 0;

    if (frame_no < prev_frame_no || consecutive_non_reference_pictures>1)
      error("POC type 2 cannot be applied for the coding pattern where the encoding /decoding order of pictures are different from the output order.\n", -1);
    prev_frame_no = frame_no;
  }

  if (IMG_NUMBER == 0)
    ReportFirstframe(tmp_time);
  else
  {
    switch (img->type)
    {
    case I_SLICE:
      stat->bit_ctr_P += stat->bit_ctr - stat->bit_ctr_n;
      ReportIntra(tmp_time);
      break;
    case SP_SLICE:
      stat->bit_ctr_P += stat->bit_ctr - stat->bit_ctr_n;
      ReportSP(tmp_time);
      break;
    case BS_IMG:
      stat->bit_ctr_P += stat->bit_ctr - stat->bit_ctr_n;
      ReportBS(tmp_time);
      break;
    case B_SLICE:
      stat->bit_ctr_B += stat->bit_ctr - stat->bit_ctr_n;
      ReportB(tmp_time);
      break;
    default:      // P, P_MULTPRED?
      stat->bit_ctr_P += stat->bit_ctr - stat->bit_ctr_n;
      ReportP(tmp_time);
    }
  }
  stat->bit_ctr_n = stat->bit_ctr;

  FreeSourceframe (srcframe);

  if (IMG_NUMBER == 0)
    return 0;
  else
    return 1;
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
static int writeout_picture(Picture *pic)
{
  Bitstream *currStream;
  int partition, slice;
  Slice *currSlice;

  img->currentPicture=pic;

  for (slice=0; slice<pic->no_slices; slice++)
  {
    currSlice = pic->slices[slice];
    for (partition=0; partition<currSlice->max_part_nr; partition++)
    {
      currStream = (currSlice->partArr[partition]).bitstream;
      assert (currStream->bits_to_go = 8);    //! should always be the case, the 
                                              //! byte alignment is done in terminate_slice
      writeUnit (currSlice->partArr[partition].bitstream);

    }           // partition loop
  }           // slice loop
  return 0;   
}


/*!
 ************************************************************************
 * \brief
 *    Encodes a frame picture
 ************************************************************************
 */
void frame_picture (Picture *frame)
{
  // if more than one B pictures, they will overwrite refFrArr_top, refFrArr_bot
  if (input->InterlaceCodingOption != FRAME_CODING && mb_adaptive && img->type == B_SLICE)
  {
    copy_mv_to_or_from_save (TO_SAVE);
  }

  img->structure = FRAME;
  enc_frame_picture  = alloc_storable_picture (img->structure, img->width, img->height, img->width_cr, img->height_cr);
  enc_frame_picture->poc=img->framepoc;
  enc_frame_picture->pic_num = img->frame_num;
  enc_frame_picture->coded_frame = 1;
  img->PicSizeInMbs = img->FrameSizeInMbs;

  enc_frame_picture->mb_adaptive_frame_field_flag = img->MbaffFrameFlag = (input->InterlaceCodingOption == MB_CODING);

  enc_picture=enc_frame_picture;


  stat->em_prev_bits_frm = 0;
  stat->em_prev_bits = &stat->em_prev_bits_frm;

  if (img->MbaffFrameFlag)
  {
    CopyTopFieldToOldImgOrgVariables (srcframe);
    CopyBottomFieldToOldImgOrgVariables (srcframe);
  }

  if (img->type != I_SLICE && (input->WeightedPrediction == 1 || (input->WeightedBiprediction > 0 && (img->type == B_SLICE))))
  {
    estimate_weighting_factor ();
  }

  img->fld_flag = 0;
  code_a_picture(frame);

  if (input->InterlaceCodingOption >= MB_CODING && mb_adaptive && img->number != 0 && input->successive_Bframe != 0 && img->type != B_SLICE)
  {
    copy_mv_to_or_from_save (TO_SAVE);
  }

  frame->bits_per_picture = 8 * ((((img->currentSlice)->partArr[0]).bitstream)->byte_pos);

  if (input->InterlaceCodingOption != FRAME_CODING)
    {
      find_distortion (snr, img);      
      frame->distortion_y = snr->snr_y;
      frame->distortion_u = snr->snr_u;
      frame->distortion_v = snr->snr_v;
    }
}


/*!
 ************************************************************************
 * \brief
 *    Encodes a field picture, consisting of top and bottom field
 ************************************************************************
 */
void field_picture (Picture *top, Picture *bottom)
{
  stat->em_prev_bits_fld = 0;
  stat->em_prev_bits = &stat->em_prev_bits_fld;
  img->number *= 2;
  img->buf_cycle *= 2;
  input->num_reference_frames = 2 * input->num_reference_frames + 1;
  img->height = input->img_height / 2;
  img->height_cr = input->img_height / 4;
  img->fld_flag = 1;
  img->PicSizeInMbs = img->FrameSizeInMbs/2;
  // Top field
  
//  img->bottom_field_flag = 0;
  put_buffer_top ();
  init_field ();
  if (img->type == B_SLICE)       //all I- and P-frames
    nextP_tr_fld--;

  CopyTopFieldToOldImgOrgVariables (srcframe);

  if (img->type != I_SLICE && (input->WeightedPrediction == 1 || (input->WeightedBiprediction > 0 && (img->type == B_SLICE || img->type == BS_IMG))))
  {
    estimate_weighting_factor ();
  }
  img->fld_flag = 1;
//  img->bottom_field_flag = 0;
 
  code_a_picture(top_pic);

  if (img->type != B_SLICE)       //all I- and P-frames
    interpolate_frame_to_fb ();

  top->bits_per_picture = 8 * ((((img->currentSlice)->partArr[0]).bitstream)->byte_pos);

  //  Bottom field
//  img->bottom_field_flag = 0;
  put_buffer_bot ();
  img->number++;

  init_field ();

  if (img->type == B_SLICE)       //all I- and P-frames
    nextP_tr_fld++;             //check once coding B field

  if (img->type == I_SLICE)
    img->type = P_SLICE;

  CopyBottomFieldToOldImgOrgVariables (srcframe);

  if (img->type != I_SLICE && (input->WeightedPrediction == 1 || (input->WeightedBiprediction > 0 && (img->type == B_SLICE || img->type == BS_IMG))))
  {
    estimate_weighting_factor ();
  }
  img->fld_flag = 1;
//  img->bottom_field_flag = 1;
  code_a_picture(bottom_pic);

  if (img->type != B_SLICE)       //all I- and P-frames
    interpolate_frame_to_fb (); 

  bottom->bits_per_picture = 8 * ((((img->currentSlice)->partArr[0]).bitstream)->byte_pos);

  // the distortion for a field coded frame (consisting of top and bottom field)
  // lives in the top->distortion varaibles, thye bottom-> are dummies
  distortion_fld (&top->distortion_y, &top->distortion_u, &top->distortion_v);
}


/*!
 ************************************************************************
 * \brief
 *    Distortion Field
 ************************************************************************
 */
static void distortion_fld (float *dis_fld_y, float *dis_fld_u, float *dis_fld_v)
{

  img->number /= 2;
  img->buf_cycle /= 2;
  img->height = input->img_height;
  img->height_cr = input->img_height / 2;
  img->total_number_mb =
    (img->width * img->height) / (MB_BLOCK_SIZE * MB_BLOCK_SIZE);
  input->num_reference_frames = (input->num_reference_frames - 1) / 2;
/*
  combine_field ();

  imgY = imgY_com;
  imgUV = imgUV_com;
  imgY_org = imgY_org_frm;
  imgUV_org = imgUV_org_frm;

  find_distortion (snr, img);   // find snr from original frame picture
*/
  *dis_fld_y = snr->snr_y;
  *dis_fld_u = snr->snr_u;
  *dis_fld_v = snr->snr_v;
}


/*!
 ************************************************************************
 * \brief
 *    Picture Structure Decision
 ************************************************************************
 */
static int picture_structure_decision (Picture *frame, Picture *top, Picture *bot)
{
  double lambda_picture;
  int spframe = (img->type == SP_SLICE);
  int bframe = (img->type == B_SLICE);
  float snr_frame, snr_field;
  int bit_frame, bit_field;

  lambda_picture = 0.85 * pow (2, (img->qp - SHIFT_QP) / 3.0) * (bframe
                                                                 || spframe ?
                                                                 4 : 1);

  snr_frame = frame->distortion_y + frame->distortion_u + frame->distortion_v;
  //! all distrortions of a field picture are accumulated in the top field
  snr_field = top->distortion_y + top->distortion_u + top->distortion_v;
  bit_field = top->bits_per_picture + bot->bits_per_picture;
  bit_frame = frame->bits_per_picture;

  return decide_fld_frame (snr_frame, snr_field, bit_field, bit_frame, lambda_picture);
}


/*!
 ************************************************************************
 * \brief
 *    Field Mode Buffer
 ************************************************************************
 */
static void field_mode_buffer (int bit_field, float snr_field_y, float snr_field_u, float snr_field_v)
{
  put_buffer_frame ();
  /*
  imgY = imgY_com;
  imgUV = imgUV_com;

  if (img->type != B_SLICE)       //all I- and P-frames 
    interpolate_frame_to_fb ();
*/
  input->no_fields += 1;

  snr->snr_y = snr_field_y;
  snr->snr_u = snr_field_u;
  snr->snr_v = snr_field_v;
}


/*!
 ************************************************************************
 * \brief
 *    Frame Mode Buffer
 ************************************************************************
 */
static void frame_mode_buffer (int bit_frame, float snr_frame_y, float snr_frame_u, float snr_frame_v)
{
  put_buffer_frame ();

  if (img->type != B_SLICE)       //all I- and P-frames
    interpolate_frame_to_fb ();

  if (input->InterlaceCodingOption != FRAME_CODING)
  {
    img->height = img->height / 2;
    img->height_cr = img->height_cr / 2;
    img->number *= 2;
    img->buf_cycle *= 2;
    
    put_buffer_top ();
    split_field_top ();
    
    if (img->type != B_SLICE)   //all I- and P-frames
    {
      interpolate_frame ();
    }
    img->number++;
    put_buffer_bot ();
    split_field_bot ();
    
    if (img->type != B_SLICE)   //all I- and P-frames
    {
      interpolate_frame ();
    }
    
    img->number /= 2;         // reset the img->number to field
    img->buf_cycle /= 2;
    img->height = input->img_height;
    img->height_cr = input->img_height / 2;
    img->total_number_mb =
      (img->width * img->height) / (MB_BLOCK_SIZE * MB_BLOCK_SIZE);
    
    snr->snr_y = snr_frame_y;
    snr->snr_u = snr_frame_u;
    snr->snr_v = snr_frame_v;
    put_buffer_frame ();
    
  }
}


/*!
 ************************************************************************
 * \brief
 *    mmco initializations should go here
 ************************************************************************
 */
static void init_dec_ref_pic_marking_buffer()
{
  img->dec_ref_pic_marking_buffer=NULL;
}


/*!
 ************************************************************************
 * \brief
 *    Initializes the parameters for a new frame
 ************************************************************************
 */
static void init_frame ()
{
  int i, j, k;
  int prevP_no, nextP_no;
  if (input->InterlaceCodingOption >= MB_CODING)
    img->structure = 3;
  else
    img->structure = 0;           //frame coding

  last_P_no = last_P_no_frm;

  img->current_mb_nr = 0;
  img->current_slice_nr = 0;
  stat->bit_slice = 0;

  img->mb_y = img->mb_x = 0;
  img->block_y = img->pix_y = img->pix_c_y = 0; 
  img->block_x = img->pix_x = img->block_c_x = img->pix_c_x = 0;

  if (img->type != B_SLICE)
  {
    img->tr = start_tr_in_this_IGOP + IMG_NUMBER * (input->jumpd + 1);
    
    img->imgtr_last_P_frm = img->imgtr_next_P_frm;
    img->imgtr_next_P_frm = img->tr;
    
#ifdef _ADAPT_LAST_GROUP_
    if (input->last_frame && img->number + 1 == input->no_frames)
      img->tr = input->last_frame;
#endif
    
    if (IMG_NUMBER != 0 && input->successive_Bframe != 0)     // B pictures to encode
      nextP_tr_frm = img->tr;
    
    if (img->type == I_SLICE)
      img->qp = input->qp0;   // set quant. parameter for I-frame
    else
    {
#ifdef _CHANGE_QP_
      if (input->qp2start > 0 && img->tr >= input->qp2start)
        img->qp = input->qpN2;
      else
#endif
        img->qp = input->qpN;
      
      if (img->type == SP_SLICE)
      {
        img->qp = input->qpsp;
        img->qpsp = input->qpsp_pred;
      }
      
    }

    img->mb_y_intra = img->mb_y_upd;  //  img->mb_y_intra indicates which GOB to intra code for this frame
    
    if (input->intra_upd > 0) // if error robustness, find next GOB to update
    {
      img->mb_y_upd = (IMG_NUMBER / input->intra_upd) % (img->height / MB_BLOCK_SIZE);
    }
  }
  else
  {
    img->p_interval = input->jumpd + 1;
    prevP_no = start_tr_in_this_IGOP + (IMG_NUMBER - 1) * img->p_interval;
    nextP_no = start_tr_in_this_IGOP + (IMG_NUMBER) * img->p_interval;
    
#ifdef _ADAPT_LAST_GROUP_
    last_P_no[0] = prevP_no;
    for (i = 1; i < img->buf_cycle; i++)
      last_P_no[i] = last_P_no[i - 1] - img->p_interval;
    
    if (input->last_frame && img->number + 1 == input->no_frames)
    {
      nextP_no = input->last_frame;
      img->p_interval = nextP_no - prevP_no;
    }
#endif
    
    img->b_interval =
      (int) ((float) (input->jumpd + 1) / (input->successive_Bframe + 1.0) +
      0.49999);
    
    img->tr = prevP_no + img->b_interval * img->b_frame_to_code;      // from prev_P

    if (img->tr >= nextP_no)
      img->tr = nextP_no - 1;
    
#ifdef _CHANGE_QP_
    if (input->qp2start > 0 && img->tr >= input->qp2start)
      img->qp = input->qpB2;
    else
#endif
      img->qp = input->qpB;

    // initialize arrays
    for (k = 0; k < 2; k++)
      for (i = 0; i < img->height / BLOCK_SIZE; i++)
        for (j = 0; j < img->width / BLOCK_SIZE + 4; j++)
        {
          tmp_fwMV[k][i][j] = 0;
          tmp_bwMV[k][i][j] = 0;
          dfMV[k][i][j] = 0;
          dbMV[k][i][j] = 0;
        }
    for (i = 0; i < img->height / BLOCK_SIZE; i++)
      for (j = 0; j < img->width / BLOCK_SIZE; j++)
      {
        fw_refFrArr[i][j] = bw_refFrArr[i][j] = -1;
      }

  }
  
  UpdateSubseqInfo (img->layer);        // Tian Dong (Sept 2002)
  UpdateSceneInformation (0, 0, 0, -1); // JVT-D099, scene information SEI, nothing included by default

  //! Commented out by StW, needs fixing in SEI.h to keep the trace file clean
  //  PrepareAggregationSEIMessage ();

  // JVT-D097
  if (img->type != B_SLICE && input->num_slice_groups_minus1 == 1 && input->FmoType > 3)
    {
      if (fmo_evlv_NewPeriod)
        FmoInitEvolvingMBAmap (input->FmoType, img->width / 16,
                               img->height / 16, MBAmap);

      FmoUpdateEvolvingMBAmap (input->FmoType, img->width / 16,
                               img->height / 16, MBAmap);
    }
  // End JVT-D097
  img->total_number_mb = (img->width * img->height) / (MB_BLOCK_SIZE * MB_BLOCK_SIZE);

  img->no_output_of_prior_pics_flag = 0;
  img->long_term_reference_flag = 0;

  init_dec_ref_pic_marking_buffer();
}

/*!
 ************************************************************************
 * \brief
 *    Initializes the parameters for a new field
 ************************************************************************
 */
static void init_field ()
{
  int i, j, k;
  int prevP_no, nextP_no;

  last_P_no = last_P_no_fld;

  //picture structure
  if (!img->fld_type)
    img->structure = 1;
  else
    img->structure = 2;

  img->current_mb_nr = 0;
  img->current_slice_nr = 0;
  stat->bit_slice = 0;

  input->jumpd *= 2;
  input->successive_Bframe *= 2;
  img->number /= 2;
  img->buf_cycle /= 2;

  img->mb_y = img->mb_x = 0;
  img->block_y = img->pix_y = img->pix_c_y = 0; // define vertical positions
  img->block_x = img->pix_x = img->block_c_x = img->pix_c_x = 0;        // define horizontal positions

  if (img->type != B_SLICE)
    {
      img->tr = img->number * (input->jumpd + 2) + img->fld_type;

      if (!img->fld_type)
        {
          img->imgtr_last_P_fld = img->imgtr_next_P_fld;
          img->imgtr_next_P_fld = img->tr;
        }

#ifdef _ADAPT_LAST_GROUP_
      if (input->last_frame && img->number + 1 == input->no_frames)
        img->tr = input->last_frame;
#endif
      if (img->number != 0 && input->successive_Bframe != 0)    // B pictures to encode
        nextP_tr_fld = img->tr;

      if (img->type == I_SLICE)
        img->qp = input->qp0;   // set quant. parameter for I-frame
      else
        {
#ifdef _CHANGE_QP_
          if (input->qp2start > 0 && img->tr >= input->qp2start)
            img->qp = input->qpN2;
          else
#endif
            img->qp = input->qpN;
          if (img->type == SP_SLICE)
            {
              img->qp = input->qpsp;
              img->qpsp = input->qpsp_pred;
            }

        }

      img->mb_y_intra = img->mb_y_upd;  //  img->mb_y_intra indicates which GOB to intra code for this frame

      if (input->intra_upd > 0) // if error robustness, find next GOB to update
        {
          img->mb_y_upd =
            (img->number / input->intra_upd) % (img->width / MB_BLOCK_SIZE);
        }
    }
  else
    {
      img->p_interval = input->jumpd + 2;
      prevP_no = (img->number - 1) * img->p_interval + img->fld_type;
      nextP_no = img->number * img->p_interval + img->fld_type;
#ifdef _ADAPT_LAST_GROUP_
      if (!img->fld_type)       // top field
        {
          last_P_no[0] = prevP_no + 1;
          last_P_no[1] = prevP_no;
          for (i = 1; i <= img->buf_cycle; i++)
            {
              last_P_no[2 * i] = last_P_no[2 * i - 2] - img->p_interval;
              last_P_no[2 * i + 1] = last_P_no[2 * i - 1] - img->p_interval;
            }
        }
      else                      // bottom field
        {
          last_P_no[0] = nextP_no - 1;
          last_P_no[1] = prevP_no;
          for (i = 1; i <= img->buf_cycle; i++)
            {
              last_P_no[2 * i] = last_P_no[2 * i - 2] - img->p_interval;
              last_P_no[2 * i + 1] = last_P_no[2 * i - 1] - img->p_interval;
            }
        }

      if (input->last_frame && img->number + 1 == input->no_frames)
        {
          nextP_no = input->last_frame;
          img->p_interval = nextP_no - prevP_no;
        }
#endif

      img->b_interval =
        (int) ((float) (input->jumpd + 1) / (input->successive_Bframe + 1.0) +
               0.49999);

      img->tr = prevP_no + (img->b_interval + 1) * img->b_frame_to_code;        // from prev_P
      if (img->tr >= nextP_no)
        img->tr = nextP_no - 1; // ?????

#ifdef _CHANGE_QP_
      if (input->qp2start > 0 && img->tr >= input->qp2start)
        img->qp = input->qpB2;
      else
#endif
        img->qp = input->qpB;

      // initialize arrays
      for (k = 0; k < 2; k++)
        for (i = 0; i < img->height / BLOCK_SIZE; i++)
          for (j = 0; j < img->width / BLOCK_SIZE + 4; j++)
            {
              tmp_fwMV[k][i][j] = 0;
              tmp_bwMV[k][i][j] = 0;
              dfMV[k][i][j] = 0;
              dbMV[k][i][j] = 0;
            }
      for (i = 0; i < img->height / BLOCK_SIZE; i++)
        for (j = 0; j < img->width / BLOCK_SIZE; j++)
          {
            fw_refFrArr[i][j] = bw_refFrArr[i][j] = -1;
          }

    }
  input->jumpd /= 2;
  input->successive_Bframe /= 2;
  img->buf_cycle *= 2;
  img->number = 2 * img->number + img->fld_type;
  img->total_number_mb = (img->width * img->height) / (MB_BLOCK_SIZE * MB_BLOCK_SIZE);
}


#define Clip(min,max,val) (((val)<(min))?(min):(((val)>(max))?(max):(val)))

/*!
 ************************************************************************
 * \brief
 *    Estimates reference picture weighting factors
 ************************************************************************
 */
static void estimate_weighting_factor ()
{
  int i, j, n;
  int x,z;
  int dc_org = 0;
  int index;
  int comp;
  int p0, pt;
  int fwd_ref[MAX_REFERENCE_PICTURES], bwd_ref[MAX_REFERENCE_PICTURES];

  int bframe = (img->type == B_SLICE) || (img->type == BS_IMG);
  int num_ref = min (img->number-((mref==mref_fld)&&img->fld_type&&bframe), img->buf_cycle);
  int dc_ref[MAX_REFERENCE_PICTURES];
  int log_weight_denom;
  int num_bwd_ref, num_fwd_ref;
  pel_t*  ref_pic;   
  pel_t*  ref_pic_w;   
  int default_weight;
  int default_weight_chroma;

  luma_log_weight_denom = 5;
  chroma_log_weight_denom = 5;
  wp_luma_round = 1 << (luma_log_weight_denom - 1);
  wp_chroma_round = 1 << (chroma_log_weight_denom - 1);
  default_weight = 1<<luma_log_weight_denom;
  default_weight_chroma = 1<<chroma_log_weight_denom;

  /* set all values to defaults */
  for (i = 0; i < 2; i++)
    for (j = 0; j < MAX_REFERENCE_PICTURES; j++)
      for (n = 0; n < 3; n++)
      {
        wp_weight[i][j][n] = default_weight;
        wp_offset[i][j][n] = 0;
      }

  for (i = 0; i < img->height; i++)
    for (j = 0; j < img->width; j++)
    {
      dc_org += imgY_org[i][j];
    }


 for (n = 0; n < num_ref; n++)
 {
   dc_ref[n] = 0;
   
   ref_pic       = img->type==B_SLICE? Refbuf11 [n] : Refbuf11[n];
   ref_pic_w       = img->type==B_SLICE? Refbuf11_w [n] : Refbuf11_w[n];
   
   // Y
   for (i = 0; i < img->height * img->width; i++)
   {
     dc_ref[n] += ref_pic[i];
   }
   
   if (dc_ref[n] != 0)
     weight[n][0] = (int) (default_weight * (double) dc_org / (double) dc_ref[n] + 0.5);
   else
     weight[n][0] = 2*default_weight;  // only used when reference picture is black
   
	  printf("dc_org = %d, dc_ref = %d, weight[%d] = %d\n",dc_org, dc_ref[n],n,weight[n][0]);
    
    /* for now always use default weight for chroma weight */
    weight[n][1] = default_weight_chroma;
    weight[n][2] = default_weight_chroma;



    /* store weighted reference pic for motion estimation */
    for (i = 0; i < img->height * img->width; i++)
    {
      ref_pic_w[i] = Clip (0, 255, ((int) ref_pic[i] * weight[n][0] + wp_luma_round) / default_weight);
    }
    for (i = 0; i < 4*(img->height + 2*IMG_PAD_SIZE) ; i++)
    {
      for (j = 0; j< 4*(img->width + 2*IMG_PAD_SIZE); j++)
      {
        mref_w[n][i][j] =   Clip (0, 255, ((int) mref[n][i][j] * weight[n][0] + wp_luma_round) / default_weight);
      }
    }
 }

  if ((img->type == P_SLICE)||(img->type == SP_SLICE))
  {
	  num_bwd_ref = 0;
	  num_fwd_ref = num_ref;
  }
  else
  {
    num_bwd_ref = (img->type == BS_IMG) ? num_ref : 1;
    num_fwd_ref = (img->type == BS_IMG) ? num_ref+1 : num_ref;
  }

//	printf("num_fwd_ref = %d num_bwd_ref = %d\n",num_fwd_ref,num_bwd_ref);

  {                             /* forward list */
    if ((img->type == P_SLICE || img->type == SP_SLICE) && input->WeightedPrediction)
    {
      for (index = 0; index < num_ref; index++)
      {
        wp_weight[0][index][0] = weight[index][0];
        wp_weight[0][index][1] = weight[index][1];
        wp_weight[0][index][2] = weight[index][2];
        // printf ("wp weight[%d] = %d  \n", index, wp_weight[0][index][0]);
      }
    }
    else if (img->type == BS_IMG && (input->WeightedBiprediction == 1))
    {
      for (index = 0; index < num_ref; index++)
      {
        wp_weight[0][index][0] = weight[index][0];
        wp_weight[0][index][1] = weight[index][1];
        wp_weight[0][index][2] = weight[index][2];
      }
      for (index = 0; index < num_ref; index++)
      {                     /* backward list */
        if (index == 0)
          n = 1;
        else if (index == 1)
          n = 0;
        else
          n = index;
      }
    }
    else if (img->type == B_SLICE && (input->WeightedBiprediction == 1))
    {
      for (index = 0; index < num_ref - 1; index++)
      {
        wp_weight[0][index][0] = weight[index + 1][0];
        wp_weight[0][index][1] = weight[index + 1][1];
        wp_weight[0][index][2] = weight[index + 1][2];
      }
      wp_weight[1][0][0] = weight[0][0];
      wp_weight[1][0][1] = weight[0][1];
      wp_weight[1][0][2] = weight[0][2];
    }
    else
    {
      for (index = 0; index < num_ref; index++)
      {
        wp_weight[0][index][0] = 1<<luma_log_weight_denom;
        wp_weight[0][index][1] = 1<<chroma_log_weight_denom;
        wp_weight[0][index][2] = 1<<chroma_log_weight_denom;
        wp_weight[1][index][0] = 1<<luma_log_weight_denom;
        wp_weight[1][index][1] = 1<<chroma_log_weight_denom;
        wp_weight[1][index][2] = 1<<chroma_log_weight_denom;
      }
    }

    if (input->WeightedBiprediction > 0 && (img->type == B_SLICE || img->type == BS_IMG))
    {
      if (img->type == BS_IMG )
      {
        for (index = 0; index < num_fwd_ref; index++)
        {
          fwd_ref[index] = index;
          if (index == 0)
            n = 1;
          else if (index == 1)
            n = 0;
          else
            n = index;
          bwd_ref[index] = n;
        }
      }
      else if (img->type == B_SLICE)
      {
        for (index = 0; index < num_fwd_ref; index++)
        {
          fwd_ref[index] = index+1;
        }
        bwd_ref[0] = 0; // only one possible backwards ref for traditional B picture in current software
      }
    }      


    if (img->type == B_SLICE || img->type == BS_IMG) // need to fill in wbp_weight values
    { 
      
      for (i = 0; i < num_fwd_ref; i++)
      {
        for (j = 0; j < num_bwd_ref; j++)
        {
          for (comp = 0; comp < 3; comp++)
          {
            log_weight_denom = (comp == 0) ? luma_log_weight_denom : chroma_log_weight_denom;
            if (input->WeightedBiprediction == 1)
            {
              wbp_weight[0][i][j][comp] = wp_weight[0][i][comp];
              wbp_weight[1][i][j][comp] = wp_weight[1][j][comp];
            }
            else if (input->WeightedBiprediction == 2)
            { // implicit mode
              pt = poc_distance (fwd_ref[i], bwd_ref[j]);
              p0 = poc_distance (fwd_ref[i], -1);
              if (pt == 0)
              {
                wbp_weight[1][i][j][comp] =  32 ;
                wbp_weight[0][i][j][comp] = 32;
              }	
              else
              {
                x = (16384 + (pt>>1))/pt;
                z = Clip(-1024, 1023, (x*p0 + 32 )>>6);
                wbp_weight[1][i][j][comp] = z>>2;
                if (wbp_weight[1][i][j][comp] < -64 || wbp_weight[1][i][j][comp] >128)
                  wbp_weight[1][i][j][comp] = 32;
                wbp_weight[0][i][j][comp] = 64 - wbp_weight[1][i][j][comp];
                
              }
              // if (comp == 0 )
              //   printf ("bpw weight[%d][%d] = %d  , %d \n", i, j, wbp_weight[0][i][j][0], wbp_weight[1][i][j][0]);
            }
          }
        }
      }
    }
  }
}

/*!
 ************************************************************************
 * \brief
 *    Choose interpolation method depending on MV-resolution
 ************************************************************************
 */
static void interpolate_frame_to_fb ()
{                             // write to mref[]
//  UnifiedOneForthPix (imgY, imgUV[0], imgUV[1], mref[0], mcef[0][0], mcef[0][1], Refbuf11[0]);
}

/*!
 ************************************************************************
 * \brief
 *    Choose interpolation method depending on MV-resolution
 ************************************************************************
 */
static void interpolate_frame ()
{                             // write to mref[]
//  UnifiedOneForthPix (imgY, imgUV[0], imgUV[1], mref[0], mcef[0][0], mcef[0][1], Refbuf11[0]);
}

static void GenerateFullPelRepresentation (pel_t ** Fourthpel,
                                           pel_t * Fullpel, int xsize,
                                           int ysize)
{
  int x, y;
  
  for (y = 0; y < ysize; y++)
    for (x = 0; x < xsize; x++)
      PutPel_11 (Fullpel, y, x, FastPelY_14 (Fourthpel, y * 4, x * 4));
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
void UnifiedOneForthPix (StorablePicture *s)
{
  int is;
  int i, j, j4;
  int ie2, je2, jj, maxy;
  
  byte **out4Y;
  byte  *ref11;
  byte  **imgY = s->imgY;
  
  // don't upsample twice
  if (s->imgY_ups || s->imgY_11)
    return;

  s->imgY_11 = malloc ((s->size_x * s->size_y) * sizeof (byte));
  if (NULL == s->imgY_11)
    no_mem_exit("alloc_storable_picture: s->imgY_11");
  
  get_mem2D (&(s->imgY_ups), (2*IMG_PAD_SIZE + s->size_y)*4, (2*IMG_PAD_SIZE + s->size_x)*4);

  out4Y = s->imgY_ups;
  ref11 = s->imgY_11;

  for (j = -IMG_PAD_SIZE; j < s->size_y + IMG_PAD_SIZE; j++)
  {
    for (i = -IMG_PAD_SIZE; i < s->size_x + IMG_PAD_SIZE; i++)
    {
      jj = max (0, min (s->size_y - 1, j));
      is =
              (ONE_FOURTH_TAP[0][0] *
               (imgY[jj][max (0, min (s->size_x - 1, i))] +
                imgY[jj][max (0, min (s->size_x - 1, i + 1))]) +
               ONE_FOURTH_TAP[1][0] *
               (imgY[jj][max (0, min (s->size_x - 1, i - 1))] +
                imgY[jj][max (0, min (s->size_x - 1, i + 2))]) +
               ONE_FOURTH_TAP[2][0] *
               (imgY[jj][max (0, min (s->size_x - 1, i - 2))] +
                imgY[jj][max (0, min (s->size_x - 1, i + 3))]));
            img4Y_tmp[j + IMG_PAD_SIZE][(i + IMG_PAD_SIZE) * 2] = imgY[jj][max (0, min (s->size_x - 1, i))] * 1024;    // 1/1 pix pos
            img4Y_tmp[j + IMG_PAD_SIZE][(i + IMG_PAD_SIZE) * 2 + 1] = is * 32;  // 1/2 pix pos
    }
  }
  
  for (i = 0; i < (s->size_x + 2 * IMG_PAD_SIZE) * 2; i++)
  {
    for (j = 0; j < s->size_y + 2 * IMG_PAD_SIZE; j++)
    {
      j4 = j * 4;
      maxy = s->size_y + 2 * IMG_PAD_SIZE - 1;
      // change for TML4, use 6 TAP vertical filter
      is =
              (ONE_FOURTH_TAP[0][0] *
               (img4Y_tmp[j][i] + img4Y_tmp[min (maxy, j + 1)][i]) +
               ONE_FOURTH_TAP[1][0] * (img4Y_tmp[max (0, j - 1)][i] +
                                       img4Y_tmp[min (maxy, j + 2)][i]) +
               ONE_FOURTH_TAP[2][0] * (img4Y_tmp[max (0, j - 2)][i] +
                                       img4Y_tmp[min (maxy, j + 3)][i])) / 32;
      
      PutPel_14 (out4Y, (j - IMG_PAD_SIZE) * 4, (i - IMG_PAD_SIZE * 2) * 2, (pel_t) max (0, min (255, (int) ((img4Y_tmp[j][i] + 512) / 1024))));  // 1/2 pix
      PutPel_14 (out4Y, (j - IMG_PAD_SIZE) * 4 + 2, (i - IMG_PAD_SIZE * 2) * 2, (pel_t) max (0, min (255, (int) ((is + 512) / 1024))));   // 1/2 pix
    }
  }
  
  /* 1/4 pix */
  /* luma */
  ie2 = (s->size_x + 2 * IMG_PAD_SIZE - 1) * 4;
  je2 = (s->size_y + 2 * IMG_PAD_SIZE - 1) * 4;
  
  for (j = 0; j < je2 + 4; j += 2)
    for (i = 0; i < ie2 + 3; i += 2)
    {
      /*  '-'  */
          PutPel_14 (out4Y, j - IMG_PAD_SIZE * 4, i - IMG_PAD_SIZE * 4 + 1,
                     (pel_t) (max
                              (0,
                               min (255,
                                    (int) (FastPelY_14
                                           (out4Y, j - IMG_PAD_SIZE * 4,
                                            i - IMG_PAD_SIZE * 4) +
                                           FastPelY_14 (out4Y,
                                                        j - IMG_PAD_SIZE * 4,
                                                        min (ie2 + 2,
                                                             i + 2) -
                                                        IMG_PAD_SIZE * 4)+1) /
                                    2))));
    }
    for (i = 0; i < ie2 + 4; i++)
    {
      for (j = 0; j < je2 + 3; j += 2)
      {
        if (i % 2 == 0)
        {
          /*  '|'  */
          PutPel_14 (out4Y, j - IMG_PAD_SIZE * 4 + 1,
                           i - IMG_PAD_SIZE * 4,
                           (pel_t) (max
                                    (0,
                                     min (255,
                                          (int) (FastPelY_14
                                                 (out4Y, j - IMG_PAD_SIZE * 4,
                                                  i - IMG_PAD_SIZE * 4) +
                                                 FastPelY_14 (out4Y,
                                                              min (je2 + 2,
                                                                   j + 2) -
                                                              IMG_PAD_SIZE *
                                                              4,
                                                              i -
                                                              IMG_PAD_SIZE *
                                                              4)+1) / 2))));
        }
        else if ((j % 4 == 0 && i % 4 == 1) || (j % 4 == 2 && i % 4 == 3))
        {
                /*  '/'  */
                PutPel_14 (out4Y, j - IMG_PAD_SIZE * 4 + 1,
                           i - IMG_PAD_SIZE * 4,
                           (pel_t) (max
                                    (0,
                                     min (255,
                                          (int) (FastPelY_14
                                                 (out4Y, j - IMG_PAD_SIZE * 4,
                                                  min (ie2 + 2,
                                                       i + 1) -
                                                  IMG_PAD_SIZE * 4) +
                                                 FastPelY_14 (out4Y,
                                                              min (je2 + 2,
                                                                   j + 2) -
                                                              IMG_PAD_SIZE *
                                                              4,
                                                              i -
                                                              IMG_PAD_SIZE *
                                                              4 - 1) + 1) / 2))));
        }
        else
        {
                /*  '\'  */
                PutPel_14 (out4Y, j - IMG_PAD_SIZE * 4 + 1,
                           i - IMG_PAD_SIZE * 4,
                           (pel_t) (max
                                    (0,
                                     min (255,
                                          (int) (FastPelY_14
                                                 (out4Y, j - IMG_PAD_SIZE * 4,
                                                  i - IMG_PAD_SIZE * 4 - 1) +
                                                 FastPelY_14 (out4Y,
                                                              min (je2 + 2,
                                                                   j + 2) -
                                                              IMG_PAD_SIZE *
                                                              4, min (ie2 + 2,
                                                                      i + 1) -
                                                              IMG_PAD_SIZE *
                                                              4) + 1) / 2))));
              }
          }
      }

    /*  Chroma: */
/*    for (j = 0; j < img->height_cr; j++)
      {
        memcpy (outU[j], imgU[j], img->width_cr);       // just copy 1/1 pix, interpolate "online" 
        memcpy (outV[j], imgV[j], img->width_cr);
      }
*/
    // Generate 1/1th pel representation (used for integer pel MV search)
    GenerateFullPelRepresentation (out4Y, ref11, s->size_x, s->size_y);

}


/*!
 ************************************************************************
 * \brief
 *    Find SNR for all three components
 ************************************************************************
 */
static void find_snr ()
{
  int i, j;
  int diff_y, diff_u, diff_v;
  int impix;
  
  //  Calculate  PSNR for Y, U and V.
  
  //     Luma.
  impix = img->height * img->width;
  
  diff_y = 0;
  for (i = 0; i < img->width; ++i)
  {
    for (j = 0; j < img->height; ++j)
    {
      diff_y += img->quad[imgY_org[j][i] - enc_picture->imgY[j][i]];
    }
  }
  
  //     Chroma.
  diff_u = 0;
  diff_v = 0;
  
  for (i = 0; i < img->width_cr; i++)
  {
    for (j = 0; j < img->height_cr; j++)
    {
      diff_u += img->quad[imgUV_org[0][j][i] - enc_picture->imgUV[0][j][i]];
      diff_v += img->quad[imgUV_org[1][j][i] - enc_picture->imgUV[1][j][i]];
    }
  }

  //  Collecting SNR statistics
  if (diff_y != 0)
  {
    snr->snr_y = (float) (10 * log10 (65025 * (float) impix / (float) diff_y));         // luma snr for current frame
    snr->snr_u = (float) (10 * log10 (65025 * (float) impix / (float) (4 * diff_u)));   // u croma snr for current frame, 1/4 of luma samples
    snr->snr_v = (float) (10 * log10 (65025 * (float) impix / (float) (4 * diff_v)));   // v croma snr for current frame, 1/4 of luma samples
  }
  
  if (img->number == 0)
  {
    snr->snr_y1 = (float) (10 * log10 (65025 * (float) impix / (float) diff_y));        // keep luma snr for first frame
    snr->snr_u1 = (float) (10 * log10 (65025 * (float) impix / (float) (4 * diff_u)));  // keep croma u snr for first frame
    snr->snr_v1 = (float) (10 * log10 (65025 * (float) impix / (float) (4 * diff_v)));  // keep croma v snr for first frame
    snr->snr_ya = snr->snr_y1;
    snr->snr_ua = snr->snr_u1;
    snr->snr_va = snr->snr_v1;
  }
  // B pictures
  else
  {
    snr->snr_ya = (float) (snr->snr_ya * (img->number + Bframe_ctr) + snr->snr_y) / (img->number + Bframe_ctr + 1); // average snr lume for all frames inc. first
    snr->snr_ua = (float) (snr->snr_ua * (img->number + Bframe_ctr) + snr->snr_u) / (img->number + Bframe_ctr + 1); // average snr u croma for all frames inc. first
    snr->snr_va = (float) (snr->snr_va * (img->number + Bframe_ctr) + snr->snr_v) / (img->number + Bframe_ctr + 1); // average snr v croma for all frames inc. first
  }
}

/*!
 ************************************************************************
 * \brief
 *    Find distortion for all three components
 ************************************************************************
 */
static void find_distortion ()
{
  int i, j;
  int diff_y, diff_u, diff_v;
  int impix;
  
  //  Calculate  PSNR for Y, U and V.
  
  //     Luma.
  impix = img->height * img->width;
  
  diff_y = 0;
  for (i = 0; i < img->width; ++i)
  {
    for (j = 0; j < img->height; ++j)
    {
      diff_y += img->quad[abs (imgY_org[j][i] - enc_picture->imgY[j][i])];
    }
  }
  
  //     Chroma.
  
  diff_u = 0;
  diff_v = 0;
  
  for (i = 0; i < img->width_cr; i++)
  {
    for (j = 0; j < img->height_cr; j++)
    {
      diff_u += img->quad[abs (imgUV_org[0][j][i] - enc_picture->imgUV[0][j][i])];
      diff_v += img->quad[abs (imgUV_org[1][j][i] - enc_picture->imgUV[1][j][i])];
    }
  }
  
  // Calculate real PSNR at find_snr_avg()
  snr->snr_y = (float) diff_y;
  snr->snr_u = (float) diff_u;
  snr->snr_v = (float) diff_v;
}


static void store_field_MV (int frame_number)
{
    int i, j;

    if (img->type != B_SLICE)     //all I- and P-frames
      {
        if (img->fld_flag)
          {
            for (i = 0; i < img->width / 4 + 4; i++)
              {
                for (j = 0; j < img->height / 8; j++)
                  {
                    tmp_mv_frm[0][2 * j][i] = tmp_mv_frm[0][2 * j + 1][i] =
                      tmp_mv_top[0][j][i];
                    tmp_mv_frm[0][2 * j][i] = tmp_mv_frm[0][2 * j + 1][i] = tmp_mv_top[0][j][i];        // ??
                    tmp_mv_frm[1][2 * j][i] = tmp_mv_frm[1][2 * j + 1][i] =
                      tmp_mv_top[1][j][i] * 2;
                    tmp_mv_frm[1][2 * j][i] = tmp_mv_frm[1][2 * j + 1][i] = tmp_mv_top[1][j][i] * 2;    // ??

                    if (input->direct_type
                        && (input->successive_Bframe != 0
                            || input->StoredBPictures > 0)
                        && (i < img->width / 4))
                      {

                        moving_block_frm[2 * j + 1][i] =
                          moving_block_frm[2 * j][i] =
                          ((refFrArr_top[j][i] != 0)
                           || (refFrArr_bot[j][i] != 0)
                           || (abs (tmp_mv_top[0][j][i + 4]) >> 1)
                           || (abs (tmp_mv_top[1][j][i + 4]) >> 1)
                           || (abs (tmp_mv_bot[0][j][i + 4]) >> 1)
                           || (abs (tmp_mv_bot[1][j][i + 4]) >> 1));

                        moving_block_top[j][i] = ((refFrArr_top[j][i] != 0)
                                                  ||
                                                  (abs
                                                   (tmp_mv_top[0][j][i + 4])
                                                   >> 1)
                                                  ||
                                                  (abs
                                                   (tmp_mv_top[1][j][i + 4])
                                                   >> 1));


                        moving_block_bot[j][i] = ((refFrArr_bot[j][i] != 0)
                                                  ||
                                                  (abs
                                                   (tmp_mv_bot[0][j][i + 4])
                                                   >> 1)
                                                  ||
                                                  (abs
                                                   (tmp_mv_bot[1][j][i + 4])
                                                   >> 1));
                      }

                    if ((i % 2 == 0) && (j % 2 == 0) && (i < img->width / 4))
                      {
                        if (refFrArr_top[j][i] == -1)
                          {
                            refFrArr_frm[2 * j][i] =
                              refFrArr_frm[2 * j + 1][i] = -1;
                            refFrArr_frm[2 * (j + 1)][i] =
                              refFrArr_frm[2 * (j + 1) + 1][i] = -1;
                            refFrArr_frm[2 * j][i + 1] =
                              refFrArr_frm[2 * j + 1][i + 1] = -1;
                            refFrArr_frm[2 * (j + 1)][i + 1] =
                              refFrArr_frm[2 * (j + 1) + 1][i + 1] = -1;
                          }
                        else
                          {
                            refFrArr_frm[2 * j][i] =
                              refFrArr_frm[2 * j + 1][i] =
                              (int) (refFrArr_top[j][i] / 2);
                            refFrArr_frm[2 * (j + 1)][i] =
                              refFrArr_frm[2 * (j + 1) + 1][i] =
                              (int) (refFrArr_top[j][i] / 2);
                            refFrArr_frm[2 * j][i + 1] =
                              refFrArr_frm[2 * j + 1][i + 1] =
                              (int) (refFrArr_top[j][i] / 2);
                            refFrArr_frm[2 * (j + 1)][i + 1] =
                              refFrArr_frm[2 * (j + 1) + 1][i + 1] =
                              (int) (refFrArr_top[j][i] / 2);
                          }
                      }
                  }
              }
          }
        else
          {
            for (i = 0; i < img->width / 4 + 4; i++)
              {
                for (j = 0; j < img->height / 8; j++)
                  {
                    tmp_mv_top[0][j][i] = tmp_mv_bot[0][j][i] =
                      (int) (tmp_mv_frm[0][2 * j][i]);
                    tmp_mv_top[1][j][i] = tmp_mv_bot[1][j][i] =
                      (int) ((tmp_mv_frm[1][2 * j][i]) / 2);
                    if (input->direct_type
                        && (input->successive_Bframe != 0
                            || input->StoredBPictures > 0)
                        && i < img->width / 4)
                      {
                        moving_block_top[j][i] = moving_block_bot[j][i] =
                          ((refFrArr_frm[2 * j][i] != 0)
                           || (refFrArr_frm[2 * j + 1][i] != 0)
                           || (abs (tmp_mv_frm[0][2 * j][i + 4]) >> 1)
                           || (abs (tmp_mv_frm[1][2 * j][i + 4]) >> 1)
                           || (abs (tmp_mv_frm[0][2 * j + 1][i + 4]) >> 1)
                           || (abs (tmp_mv_frm[1][2 * j + 1][i + 4]) >> 1));

                        moving_block_frm[2 * j][i] =
                          ((refFrArr_frm[2 * j][i] != 0)
                           || (abs (tmp_mv_frm[0][2 * j][i + 4]) >> 1)
                           || (abs (tmp_mv_frm[1][2 * j][i + 4]) >> 1));

                        moving_block_frm[2 * j + 1][i] =
                          ((refFrArr_frm[2 * j + 1][i] != 0)
                           || (abs (tmp_mv_frm[0][2 * j + 1][i + 4]) >> 1)
                           || (abs (tmp_mv_frm[1][2 * j + 1][i + 4]) >> 1));
                      }
                    if ((i % 2 == 0) && (j % 2 == 0) && (i < img->width / 4))
                      {
                        if (refFrArr_frm[2 * j][i] == -1)
                          {
                            refFrArr_top[j][i] = refFrArr_bot[j][i] = -1;
                            refFrArr_top[j + 1][i] = refFrArr_bot[j + 1][i] =
                              -1;
                            refFrArr_top[j][i + 1] = refFrArr_bot[j][i + 1] =
                              -1;
                            refFrArr_top[j + 1][i + 1] =
                              refFrArr_bot[j + 1][i + 1] = -1;
                          }
                        else
                          {
                            refFrArr_top[j][i] = refFrArr_bot[j][i] =
                              refFrArr_frm[2 * j][i] * 2;
                            refFrArr_top[j + 1][i] = refFrArr_bot[j + 1][i] =
                              refFrArr_frm[2 * j][i] * 2;
                            refFrArr_top[j][i + 1] = refFrArr_bot[j][i + 1] =
                              refFrArr_frm[2 * j][i] * 2;
                            refFrArr_top[j + 1][i + 1] =
                              refFrArr_bot[j + 1][i + 1] =
                              refFrArr_frm[2 * j][i] * 2;
                          }
                      }
                  }
              }
          }
      }
}

static void store_direct_moving_flag (int frame_number)
{
    int i, j;

    if (img->type != B_SLICE)     //all I- and P-frames
      for (i = 0; i < img->width / 4; i++)
        for (j = 0; j < img->height / 8; j++)
          if (input->direct_type
              && (input->successive_Bframe != 0
                  || input->StoredBPictures > 0))
            {
              moving_block_frm[2 * j][i] = ((refFrArr_frm[2 * j][i] != 0)
                                            ||
                                            (abs (tmp_mv_frm[0][2 * j][i + 4])
                                             >> 1)
                                            ||
                                            (abs (tmp_mv_frm[1][2 * j][i + 4])
                                             >> 1));
              moving_block_frm[2 * j + 1][i] =
                ((refFrArr_frm[2 * j + 1][i] != 0)
                 || (abs (tmp_mv_frm[0][2 * j + 1][i + 4]) >> 1)
                 || (abs (tmp_mv_frm[1][2 * j + 1][i + 4]) >> 1));
            }
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


/*! 
***************************************************************************
// For MB level field/frame coding
***************************************************************************
*/
void copy_rdopt_data (int bot_block)
{
  int mb_nr = img->current_mb_nr;
  Macroblock *currMB = &img->mb_data[mb_nr];
  int i, j, k, l;
  int **frefar =
    ((img->type == B_SLICE || img->type == BS_IMG) ? fw_refFrArr : refFrArr);
  int **brefar = bw_refFrArr;
  int **frefar_fld;
  int **brefar_fld = (bot_block ? bw_refFrArr_bot : bw_refFrArr_top);
  int bframe = (img->type == B_SLICE || img->type == BS_IMG);
  int mode;
  int offset_x, offset_y;
  int field_y;
  int block_y = (bot_block ? (img->block_y - 4) / 2 : (img->block_y / 2));

  frefar_fld = ((img->type == B_SLICE || img->type == BS_IMG) ? 
                (bot_block ? fw_refFrArr_bot : fw_refFrArr_top) : 
                (bot_block ? refFrArr_bot : refFrArr_top));

  mode = rdopt->mode;
  currMB->mb_type = rdopt->mb_type;   // copy mb_type 
  currMB->cbp = rdopt->cbp;   // copy cbp
  currMB->cbp_blk = rdopt->cbp_blk;   // copy cbp_blk
  img->i16offset = rdopt->i16offset;
  currMB->c_ipred_mode = rdopt->c_ipred_mode;

  if (img->type != B_SLICE && img->type != BS_IMG)
    field_mb[img->mb_y][img->mb_x] = MBPairIsField;
  
  
  if (img->field_mode)
  {
    if (img->top_field && (img->type != B_SLICE && img->type != BS_IMG)
      && mode == 0)
      TopFieldIsSkipped = 1;        //set top field MB skipped to 1 for skipped top field MBs
    else if (img->top_field && (img->type == B_SLICE || img->type == BS_IMG)
      && mode == 0)
      TopFieldIsSkipped = currMB->cbp ? 0 : 1;      // for direct mode check if cbp is zero (skipped) or not
    else if (img->top_field)        // don't change it for bottom field in case TopFieldIsSkipped=1
      TopFieldIsSkipped = 0;
/*    if(!bot_block && TopFieldIsSkipped)
      fprintf(stderr,"Skipping Top Field %d\n", img->current_mb_nr);
    else if(bot_block && mode == 0 && !TopFieldIsSkipped)
      fprintf(stderr,"Skipping Bottom Field %d\n", img->current_mb_nr);
    else if(bot_block && TopFieldIsSkipped && mode == 0)
      fprintf(stderr,"Skipping Field Pair %d\n", img->current_mb_nr);
*/
  }
  else
  {
    if (!bot_block && (img->type != B_SLICE && img->type != BS_IMG)
      && mode == 0)
      TopFrameIsSkipped = 1;        //set top frame MB skipped to 1 for skipped top field MBs
    else if (!bot_block && (img->type == B_SLICE || img->type == BS_IMG)
      && mode == 0)
      TopFrameIsSkipped = currMB->cbp ? 0 : 1;
    else if (!bot_block)
      TopFrameIsSkipped = 0;

/*    if(!bot_block && TopFrameIsSkipped)
      fprintf(stderr,"Skipping Top Frame %d\n", img->current_mb_nr);
    else if(bot_block && mode == 0 && !TopFrameIsSkipped)
      fprintf(stderr,"Skipping Bottom Frame %d\n", img->current_mb_nr);
    else if(bot_block && TopFrameIsSkipped && mode == 0)
      fprintf(stderr,"Skipping Frame Pair %d\n", img->current_mb_nr);
*/
  }
  

  for (i = 0; i < 6; i++)
    for (j = 0; j < 4; j++)
      for (k = 0; k < 2; k++)
        for (l = 0; l < 18; l++)
          img->cofAC[i][j][k][l] = rdopt->cofAC[i][j][k][l];

  for (i = 0; i < 3; i++)
    for (k = 0; k < 2; k++)
      for (l = 0; l < 18; l++)
        img->cofDC[i][k][l] = rdopt->cofDC[i][k][l];

  for (j = 0; j < 4; j++)
    for (i = 0; i < 4; i++)
    {
      frefar[img->block_y + j][img->block_x + i] = rdopt->frefar[j][i];     // copy refFrArr
      if (bframe)
        brefar[img->block_y + j][img->block_x + i] = rdopt->brefar[j][i];   // copy bw_refFrArr
    }

  if (!MBPairIsField)
  {
    //===== reconstruction values =====
    for (j = 0; j < 16; j++)
      for (i = 0; i < 16; i++)
      {
        enc_picture->imgY[img->pix_y + j][img->pix_x + i] = rdopt->rec_mbY[j][i];
      }


    for (j = 0; j < 8; j++)
      for (i = 0; i < 8; i++)
      {
        enc_picture->imgUV[0][img->pix_c_y + j][img->pix_c_x + i] = rdopt->rec_mbU[j][i];
        enc_picture->imgUV[1][img->pix_c_y + j][img->pix_c_x + i] = rdopt->rec_mbV[j][i];
      }
  }
  else
  {
    if (!bot_block)
      offset_y = img->pix_y;
    else
      offset_y = img->pix_y - MB_BLOCK_SIZE + 1;
    
  //===== reconstruction values =====
  for (j = 0; j < 16; j++)
    for (i = 0; i < 16; i++)
    {
      enc_picture->imgY[offset_y + (j << 1)][img->pix_x + i] = rdopt->rec_mbY[j][i];
    }

  if (!bot_block)
    offset_y = img->pix_c_y;
  else
    offset_y = img->pix_c_y + 1 - MB_BLOCK_SIZE / 2;
  
  for (j = 0; j < 8; j++)
    for (i = 0; i < 8; i++)
    {
      enc_picture->imgUV[0][offset_y + (j << 1)][img->pix_c_x + i] = rdopt->rec_mbU[j][i];
      enc_picture->imgUV[1][offset_y + (j << 1)][img->pix_c_x + i] = rdopt->rec_mbV[j][i];
    }
  }

  if (MBPairIsField)
  {
    if (!bot_block)
    {
      offset_x = img->pix_x;
      offset_y = img->pix_y >> 1;
      for (i = 0; i < 16; i++)
        for (j = 0; j < 16; j++)
          enc_top_picture->imgY[offset_y + i][offset_x + j] = rdopt->rec_mbY[i][j];
        
      offset_x = img->pix_c_x;
      offset_y = img->pix_c_y >> 1;

      for (i = 0; i < 8; i++)
        for (j = 0; j < 8; j++)
        {
          enc_top_picture->imgUV[0][offset_y + i][offset_x + j] = rdopt->rec_mbU[i][j];
          enc_top_picture->imgUV[1][offset_y + i][offset_x + j] = rdopt->rec_mbV[i][j];
        }
    }
    else
    {
      offset_x = img->pix_x;
      offset_y = (img->pix_y - MB_BLOCK_SIZE) >> 1;
      for (i = 0; i < 16; i++)
        for (j = 0; j < 16; j++)
          enc_bottom_picture->imgY[offset_y + i][offset_x + j] = rdopt->rec_mbY[i][j];
        
      offset_x = img->pix_c_x;
      offset_y = (img->pix_c_y - 8) >> 1;
      
      for (i = 0; i < 8; i++)
        for (j = 0; j < 8; j++)
        {
          enc_bottom_picture->imgUV[0][offset_y + i][offset_x + j] = rdopt->rec_mbU[i][j];
          enc_bottom_picture->imgUV[1][offset_y + i][offset_x + j] = rdopt->rec_mbV[i][j];
        }
    }
  }
  else
  {
    if (!bot_block)
    {
      offset_x = img->pix_x;
      offset_y = img->pix_y >> 1;
      for (i = 0; i < 8; i++)
        for (j = 0; j < 16; j++)
        {
          enc_top_picture->imgY[offset_y + i][offset_x + j] = rdopt->rec_mbY[2 * i][j];
          enc_bottom_picture->imgY[offset_y + i][offset_x + j] = rdopt->rec_mbY[2 * i + 1][j];
        }
        
      offset_x = img->pix_c_x;
      offset_y = img->pix_c_y >> 1;
        
      for (i = 0; i < 4; i++)
        for (j = 0; j < 8; j++)
        {
          enc_top_picture->imgUV[0][offset_y + i][offset_x + j] = rdopt->rec_mbU[2 * i][j];
          enc_top_picture->imgUV[1][offset_y + i][offset_x + j] = rdopt->rec_mbV[2 * i][j];
          enc_bottom_picture->imgUV[0][offset_y + i][offset_x + j] = rdopt->rec_mbU[2 * i + 1][j];
          enc_bottom_picture->imgUV[1][offset_y + i][offset_x + j] = rdopt->rec_mbV[2 * i + 1][j];
        }
    }
    else
    {
      offset_x = img->pix_x;
      offset_y = (img->pix_y - MB_BLOCK_SIZE) >> 1;
      offset_y += 8;
      
      for (i = 0; i < 8; i++)
        for (j = 0; j < 16; j++)
        {
          enc_top_picture->imgY[offset_y + i][offset_x + j] = rdopt->rec_mbY[2 * i][j];
          enc_bottom_picture->imgY[offset_y + i][offset_x + j] = rdopt->rec_mbY[2 * i + 1][j];
        }
        
      offset_x = img->pix_c_x;
      offset_y = (img->pix_c_y - 8) >> 1;
      offset_y += 4;
      
      for (i = 0; i < 4; i++)
        for (j = 0; j < 8; j++)
        {
          enc_top_picture->imgUV[0][offset_y + i][offset_x + j] = rdopt->rec_mbU[2 * i][j];
          enc_top_picture->imgUV[1][offset_y + i][offset_x + j] = rdopt->rec_mbV[2 * i][j];
          enc_bottom_picture->imgUV[0][offset_y + i][offset_x + j] = rdopt->rec_mbU[2 * i + 1][j];
          enc_bottom_picture->imgUV[1][offset_y + i][offset_x + j] = rdopt->rec_mbV[2 * i + 1][j];
        }
    }
  }

  for (i = 0; i < 4; i++)
  {
    currMB->b8mode[i] = rdopt->b8mode[i];
    currMB->b8pdir[i] = rdopt->b8pdir[i];
  }
  //==== reference frames =====
  for (j = 0; j < 4; j++)
    for (i = 0; i < 4; i++)
    {
      if (MBPairIsField)
      {
        frefar_fld[block_y + j][img->block_x + i] = rdopt->frefar[j][i];
        frefar[img->block_y + j][img->block_x + i] = rdopt->frefar[j][i] == -1 ? -1 : rdopt->frefar[j][i] / 2;    // Krit, add to match MBINTLC1 in decoder
      }
      else
      {
        frefar_fld[block_y + j][img->block_x + i] = rdopt->frefar[j][i] == -1 ? -1 : 2 * rdopt->frefar[j][i];
        frefar[img->block_y + j][img->block_x + i] = rdopt->frefar[j][i];
      }
    }
  if (bframe)
  {
    for (j = 0; j < 4; j++)
      for (i = 0; i < 4; i++)
      {
        if (MBPairIsField)
        {
          brefar_fld[block_y + j][img->block_x + i] = rdopt->brefar[j][i];
          brefar[img->block_y + j][img->block_x + i] = rdopt->brefar[j][i] == -1 ? -1 : rdopt->brefar[j][i] / 2;        // Krit, add to match MBINTLC1 in decoder
        }
        else
        {
          brefar_fld[block_y + j][img->block_x + i] = rdopt->brefar[j][i] == -1 ? -1 : 2 * rdopt->brefar[j][i];
          brefar[img->block_y + j][img->block_x + i] = rdopt->brefar[j][i];
        }
      }
  }
  if (img->type == BS_IMG)
  {
    for (j = 0; j < 4; j++)
      for (i = 0; i < 4; i++)
      {
        if (MBPairIsField)
        {
          if (bot_block)
          {
            refFrArr_bot[block_y + j][img->block_x + i] = rdopt->frefar[j][i];
          }
          else
          {
            refFrArr_top[block_y + j][img->block_x + i] = rdopt->frefar[j][i];
          }
          refFrArr[img->block_y + j][img->block_x + i] = rdopt->frefar[j][i] == -1 ? -1 : rdopt->frefar[j][i] / 2;
        }
        else
        {
          if (bot_block)
          {
            refFrArr_bot[block_y + j][img->block_x + i] = rdopt->frefar[j][i] == -1 ? -1 : 2 * rdopt->frefar[j][i];
          }
          else
          {
            refFrArr_top[block_y + j][img->block_x + i] = rdopt->frefar[j][i] == -1 ? -1 : 2 * rdopt->frefar[j][i];
          }
          refFrArr[img->block_y + j][img->block_x + i] = rdopt->frefar[j][i];
        }
      }
  }
  
  //==== intra prediction modes ====
  offset_y = MBPairIsField ? (bot_block ? (img->block_y-4)/2 : img->block_y/2) : img->block_y;
  field_y  = bot_block ? ((img->block_y - 4) / 2) : img->block_y / 2;

  if (mode == P8x8)
  {
    for (k = 0, j = 0; j < 4; j++)
      for (i = img->block_x; i < img->block_x + 4; i++, k++)
      {
//        img->ipredmode[i][img->block_y + j] = rdopt->ipredmode[i][offset_y + j];
        img->ipredmode[i][img->block_y + j] = rdopt->ipredmode[i][img->block_y + j];
        currMB->intra_pred_modes[k] = rdopt->intra_pred_modes[k];
      }
  }
  else if (mode != I4MB)
  {
    for (k = 0, j = 0; j < 4; j++)
      for (i = img->block_x; i < img->block_x + 4; i++, k++)
      {
        img->ipredmode[i][img->block_y + j] = DC_PRED;
        currMB->intra_pred_modes[k] = DC_PRED;
      }
  }
  else if (mode == I4MB)
  {
    for (k = 0, j = 0; j < 4; j++)
      for (i = img->block_x; i < img->block_x + 4; i++, k++)
      {
//        img->ipredmode[i][img->block_y + j] = rdopt->ipredmode[i][offset_y + j];
        img->ipredmode[i][img->block_y + j] = rdopt->ipredmode[i][img->block_y + j];
        currMB->intra_pred_modes[k] = rdopt->intra_pred_modes[k];
        
      }
      
  }

  // point the field motion vectors to either top or bottom field 
  // motion vectors
  copy_motion_vectors_MB (bot_block);
  
}                             // end of copy_rdopt_data
  
static void copy_motion_vectors_MB (int bot_block)
{
  int ***tmp_mv_fld, ***tmp_fwMV_fld, ***tmp_bwMV_fld;
  int *****mv_fld, *****all_mv_fld, *****all_bmv_fld;
  int *****p_fwMV_fld, *****p_bwMV_fld;
  int mode8, pdir8, l, by, bxr, bx;
  int **frefar, **refar;
  int by_f;
  int i, j, k, ref, ref_field, dref;
  int bframe = (img->type == B_SLICE || img->type == BS_IMG);
  Macroblock *currMB = &img->mb_data[img->current_mb_nr];
  int *****all_mv = img->all_mv;
  int *****all_bmv = img->all_bmv;
  int *****imgmv = img->mv;
  int *****p_fwMV = img->p_fwMV;
  int *****p_bwMV = img->p_bwMV;
  
  if (!bot_block)
  {
    tmp_mv_fld = tmp_mv_top;
    tmp_fwMV_fld = tmp_fwMV_top;
    tmp_bwMV_fld = tmp_bwMV_top;
    mv_fld = img->mv_top;
    all_mv_fld = img->all_mv_top;
    all_bmv_fld = img->all_bmv_top;
    p_fwMV_fld = img->p_fwMV_top;
    p_bwMV_fld = img->p_bwMV_top;
    frefar = fw_refFrArr_top;
    refar = refFrArr_top;
  }
  else
  {
    tmp_mv_fld = tmp_mv_bot;
    tmp_fwMV_fld = tmp_fwMV_bot;
    tmp_bwMV_fld = tmp_bwMV_bot;
    mv_fld = img->mv_bot;
    all_mv_fld = img->all_mv_bot;
    all_bmv_fld = img->all_bmv_bot;
    p_fwMV_fld = img->p_fwMV_bot;
    p_bwMV_fld = img->p_bwMV_bot;
    frefar = fw_refFrArr_bot;
    refar = refFrArr_bot;
  }

  if (MBPairIsField)
  {
    if (!bot_block)
    {
      all_mv = img->all_mv_top;
      all_bmv = img->all_bmv_top;
      imgmv = img->mv_top;
      p_fwMV = img->p_fwMV_top;
      p_bwMV = img->p_bwMV_top;
    }
    else
    {
      all_mv = img->all_mv_bot;
      all_bmv = img->all_bmv_bot;
      imgmv = img->mv_bot;
      p_fwMV = img->p_fwMV_bot;
      p_bwMV = img->p_bwMV_bot;
    }
    
  }

  if (!MBPairIsField)
  {
    for (j = 0; j < 4; j++)
      for (i = 0; i < 4; i++)
      {
        mode8 = currMB->b8mode[k = 2 * (j / 2) + (i / 2)];
        pdir8 = currMB->b8pdir[k];
        l = 2 * (j % 2) + (i % 2);
        by = img->block_y + j;
        bxr = img->block_x + i;
        bx = img->block_x + i + 4;
        
        if (!bot_block)
          by_f = img->block_y / 2 + j;
        else
          by_f = (img->block_y - 4) / 2 + j;
        
        ref = (bframe ? fw_refFrArr : refFrArr)[by][bxr];
        ref_field = (bframe ? frefar : refar)[by_f][bxr];
        
        if (!bframe)
        {
          if (mode8 != IBLOCK && ref != -1)
          {
            tmp_mv[0][by][bx] = rdopt->tmp_mv[0][j][i];
            tmp_mv[1][by][bx] = rdopt->tmp_mv[1][j][i];
          }
          else
          {
            tmp_mv[0][by][bx] = 0;
            tmp_mv[1][by][bx] = 0;
          }

          tmp_mv_fld[0][by_f][bx] = tmp_mv[0][by][bx];
          tmp_mv_fld[1][by_f][bx] = tmp_mv[1][by][bx] / 2;
        }
        else
        {
          if (pdir8 == -1)      // intra
          {
            tmp_fwMV[0][by][bx] = 0;
            tmp_fwMV[1][by][bx] = 0;
            tmp_bwMV[0][by][bx] = 0;
            tmp_bwMV[1][by][bx] = 0;
            dfMV[0][by][bx] = 0;
            dfMV[1][by][bx] = 0;
            dbMV[0][by][bx] = 0;
            dbMV[1][by][bx] = 0;
          }
          else if (pdir8 == 0)  // forward
          {
            tmp_fwMV[0][by][bx] = rdopt->tmp_fwMV[0][j][i];
            tmp_fwMV[1][by][bx] = rdopt->tmp_fwMV[1][j][i];
            tmp_bwMV[0][by][bx] = 0;
            tmp_bwMV[1][by][bx] = 0;
            dfMV[0][by][bx] = 0;
            dfMV[1][by][bx] = 0;
            dbMV[0][by][bx] = 0;
            dbMV[1][by][bx] = 0;
          }
          else if (pdir8 == 1)  // backward
          {
            tmp_fwMV[0][by][bx] = 0;
            tmp_fwMV[1][by][bx] = 0;
            tmp_bwMV[0][by][bx] = rdopt->tmp_bwMV[0][j][i];
            tmp_bwMV[1][by][bx] = rdopt->tmp_bwMV[1][j][i];
            dfMV[0][by][bx] = 0;
            dfMV[1][by][bx] = 0;
            dbMV[0][by][bx] = 0;
            dbMV[1][by][bx] = 0;
          }
          else if (mode8 != 0)  // bidirect
          {
            tmp_fwMV[0][by][bx] = rdopt->tmp_fwMV[0][j][i];
            tmp_fwMV[1][by][bx] = rdopt->tmp_fwMV[1][j][i];
            tmp_bwMV[0][by][bx] = rdopt->tmp_bwMV[0][j][i];
            tmp_bwMV[1][by][bx] = rdopt->tmp_bwMV[1][j][i];
            dfMV[0][by][bx] = 0;
            dfMV[1][by][bx] = 0;
            dbMV[0][by][bx] = 0;
            dbMV[1][by][bx] = 0;
          }
          else          // direct
          {
            dref = max (0, refFrArr[by][bxr]);
            tmp_fwMV[0][by][bx] = dfMV[0][by][bx] =
              rdopt->dfMV[0][j][i];
            tmp_fwMV[1][by][bx] = dfMV[1][by][bx] =
              rdopt->dfMV[1][j][i];
            tmp_bwMV[0][by][bx] = dbMV[0][by][bx] =
              rdopt->dbMV[0][j][i];
            tmp_bwMV[1][by][bx] = dbMV[1][by][bx] =
              rdopt->dbMV[1][j][i];
          }
          
          tmp_fwMV_fld[0][by_f][bx] = tmp_fwMV[0][by][bx];
          tmp_fwMV_fld[1][by_f][bx] = tmp_fwMV[1][by][bx] / 2;
          tmp_bwMV_fld[0][by_f][bx] = tmp_bwMV[0][by][bx];
          tmp_bwMV_fld[1][by_f][bx] = tmp_bwMV[1][by][bx] / 2;
          if (img->type == BS_IMG)
          {
            tmp_mv[0][by][bx] = tmp_fwMV[0][by][bx];
            tmp_mv[1][by][bx] = tmp_fwMV[1][by][bx];
            tmp_mv_fld[0][by_f][bx] = tmp_fwMV_fld[0][by_f][bx];
            tmp_mv_fld[1][by_f][bx] = tmp_fwMV_fld[1][by_f][bx];
          }
        }
      }
  }
  else
  {
    for (j = 0; j < 4; j++)
      for (i = 0; i < 4; i++)
      {
        mode8 = currMB->b8mode[k = 2 * (j / 2) + (i / 2)];
        pdir8 = currMB->b8pdir[k];
        l = 2 * (j % 2) + (i % 2);
        by = img->block_y + j;
        bxr = img->block_x + i;
        bx = img->block_x + i + 4;
        
        if (!bot_block)
          by_f = img->block_y / 2 + j;
        else
          by_f = (img->block_y - 4) / 2 + j;
        
        ref = (bframe ? fw_refFrArr : refFrArr)[by][bxr];
        ref_field = (bframe ? frefar : refar)[by_f][bxr];
        
        if (!bframe)
        {
          if (mode8 != IBLOCK && ref != -1)
          {
            tmp_mv_fld[0][by_f][bx] = rdopt->tmp_mv[0][j][i];
            tmp_mv_fld[1][by_f][bx] = rdopt->tmp_mv[1][j][i];
          }
          else
          {
            tmp_mv_fld[0][by_f][bx] = 0;
            tmp_mv_fld[1][by_f][bx] = 0;
          }
          
          tmp_mv[0][by][bx] = tmp_mv_fld[0][by_f][bx];
          tmp_mv[1][by][bx] = 2 * tmp_mv_fld[1][by_f][bx];
        }
        else
        {
          if (pdir8 == -1)      // intra
          {
            tmp_fwMV_fld[0][by_f][bx] = 0;
            tmp_fwMV_fld[1][by_f][bx] = 0;
            tmp_bwMV_fld[0][by_f][bx] = 0;
            tmp_bwMV_fld[1][by_f][bx] = 0;
            dfMV[0][by][bx] = 0;
            dfMV[1][by][bx] = 0;
            dbMV[0][by][bx] = 0;
            dbMV[1][by][bx] = 0;
          }
          else if (pdir8 == 0)  // forward
          {
            tmp_fwMV_fld[0][by_f][bx] = rdopt->tmp_fwMV[0][j][i];
            tmp_fwMV_fld[1][by_f][bx] = rdopt->tmp_fwMV[1][j][i];
            tmp_bwMV_fld[0][by_f][bx] = 0;
            tmp_bwMV_fld[1][by_f][bx] = 0;
            dfMV[0][by][bx] = 0;
            dfMV[1][by][bx] = 0;
            dbMV[0][by][bx] = 0;
            dbMV[1][by][bx] = 0;
          }
          else if (pdir8 == 1)  // backward
          {
            tmp_fwMV_fld[0][by_f][bx] = 0;
            tmp_fwMV_fld[1][by_f][bx] = 0;
            tmp_bwMV_fld[0][by_f][bx] = rdopt->tmp_bwMV[0][j][i];
            tmp_bwMV_fld[1][by_f][bx] = rdopt->tmp_bwMV[1][j][i];
            dfMV[0][by][bx] = 0;
            dfMV[1][by][bx] = 0;
            dbMV[0][by][bx] = 0;
            dbMV[1][by][bx] = 0;
          }
          else if (mode8 != 0)  // bidirect
          {
            tmp_fwMV_fld[0][by_f][bx] = rdopt->tmp_fwMV[0][j][i];
            tmp_fwMV_fld[1][by_f][bx] = rdopt->tmp_fwMV[1][j][i];
            tmp_bwMV_fld[0][by_f][bx] = rdopt->tmp_bwMV[0][j][i];
            tmp_bwMV_fld[1][by_f][bx] = rdopt->tmp_bwMV[1][j][i];
            dfMV[0][by][bx] = 0;
            dfMV[1][by][bx] = 0;
            dbMV[0][by][bx] = 0;
            dbMV[1][by][bx] = 0;
          }
          else          // direct
          {
            dref = max (0, refar[by_f][bxr]);
            tmp_fwMV_fld[0][by_f][bx] = dfMV[0][by][bx] =
              rdopt->dfMV[0][j][i];
            tmp_fwMV_fld[1][by_f][bx] = dfMV[1][by][bx] =
              rdopt->dfMV[1][j][i];
            tmp_bwMV_fld[0][by_f][bx] = dbMV[0][by][bx] =
              rdopt->dbMV[0][j][i];
            tmp_bwMV_fld[1][by_f][bx] = dbMV[1][by][bx] =
              rdopt->dbMV[1][j][i];
          }
          
          tmp_fwMV[0][by][bx] = tmp_fwMV_fld[0][by_f][bx];
          tmp_fwMV[1][by][bx] = tmp_fwMV_fld[1][by_f][bx] * 2;
          tmp_bwMV[0][by][bx] = tmp_bwMV_fld[0][by_f][bx];
          tmp_bwMV[1][by][bx] = tmp_bwMV_fld[1][by_f][bx] * 2;
          if (img->type == BS_IMG)
          {
            tmp_mv[0][by][bx] = tmp_fwMV[0][by][bx];
            tmp_mv[1][by][bx] = tmp_fwMV[1][by][bx];
            tmp_mv_fld[0][by_f][bx] = tmp_fwMV_fld[0][by_f][bx];
            tmp_mv_fld[1][by_f][bx] = tmp_fwMV_fld[1][by_f][bx];
          }
        }
      }
  }

  for (i = 0; i < 4; i++)
    for (j = 0; j < 4; j++)
      for (k = 0; k < img->buf_cycle; k++)
        for (l = 0; l < 9; l++)
        {
          all_mv[i][j][k][l][0] = rdopt->all_mv[i][j][k][l][0];
          all_bmv[i][j][k][l][0] = rdopt->all_bmv[i][j][k][l][0];
          p_fwMV[i][j][k][l][0] = rdopt->p_fwMV[i][j][k][l][0];
          p_bwMV[i][j][k][l][0] = rdopt->p_bwMV[i][j][k][l][0];
          imgmv[i][j][k][l][0] = rdopt->mv[i][j][k][l][0];
          
          all_mv[i][j][k][l][1] = rdopt->all_mv[i][j][k][l][1];
          all_bmv[i][j][k][l][1] = rdopt->all_bmv[i][j][k][l][1];
          p_fwMV[i][j][k][l][1] = rdopt->p_fwMV[i][j][k][l][1];
          p_bwMV[i][j][k][l][1] = rdopt->p_bwMV[i][j][k][l][1];
          imgmv[i][j][k][l][1] = rdopt->mv[i][j][k][l][1];
          
        }
  }
  

static void copy_mv_to_or_from_save (int direction)
{
  static int ***tmp_mv_top_save = NULL;
  static int ***tmp_mv_bot_save = NULL;
  static int **refFrArr_top_save = NULL;
  static int **refFrArr_bot_save = NULL;
  int i, j, dummy;
  const int height_field = img->height / 2;

  // Allocate memory if not yet done...
  if (tmp_mv_top_save == NULL)
    dummy = get_mem3Dint(&tmp_mv_top_save, 2, height_field/BLOCK_SIZE, img->width/BLOCK_SIZE+4);
  if (tmp_mv_bot_save == NULL)
    dummy = get_mem3Dint(&tmp_mv_bot_save, 2, height_field/BLOCK_SIZE, img->width/BLOCK_SIZE+4);
  if (refFrArr_top_save == NULL)
    dummy = get_mem2Dint(&refFrArr_top_save, height_field/BLOCK_SIZE, img->width/BLOCK_SIZE);
  if (refFrArr_bot_save == NULL)
    dummy = get_mem2Dint(&refFrArr_bot_save, height_field/BLOCK_SIZE, img->width/BLOCK_SIZE);

  // Copy
  if (direction == TO_SAVE)
  {
    for (i = 0; i < img->height / (2 * BLOCK_SIZE); i++)
      for (j = 0; j < img->width / BLOCK_SIZE; j++)
      {
        refFrArr_top_save[i][j] = refFrArr_top[i][j];
        refFrArr_bot_save[i][j] = refFrArr_bot[i][j];
      }
      if (input->successive_Bframe != 0 || input->StoredBPictures > 0)
        {
          for (i = 0; i < img->height / (2 * BLOCK_SIZE); i++)
            for (j = 0; j < img->width / BLOCK_SIZE + 4; j++)
            {
              tmp_mv_top_save[0][i][j] = tmp_mv_top[0][i][j];
              tmp_mv_bot_save[0][i][j] = tmp_mv_bot[0][i][j];
              tmp_mv_top_save[1][i][j] = tmp_mv_top[1][i][j];
              tmp_mv_bot_save[1][i][j] = tmp_mv_bot[1][i][j];
            }
        }

  }
  else if (direction == FROM_SAVE)
  {
    for (i = 0; i < img->height / (2 * BLOCK_SIZE); i++)
      for (j = 0; j < img->width / BLOCK_SIZE; j++)
      {
        refFrArr_top[i][j] = refFrArr_top_save[i][j];
        refFrArr_bot[i][j] = refFrArr_bot_save[i][j];
      }
    if (input->successive_Bframe != 0 || input->StoredBPictures > 0)
    {
      for (i = 0; i < img->height / (2 * BLOCK_SIZE); i++)
        for (j = 0; j < img->width / BLOCK_SIZE + 4; j++)
        {
          tmp_mv_top[0][i][j] = tmp_mv_top_save[0][i][j];
          tmp_mv_bot[0][i][j] = tmp_mv_bot_save[0][i][j];
          tmp_mv_top[1][i][j] = tmp_mv_top_save[1][i][j];
          tmp_mv_bot[1][i][j] = tmp_mv_bot_save[1][i][j];
        }
    }
  }
  else
  {
    printf ("copy_mv_tp_or_from_save: wrong driection %d\n", direction);
  }
}


static void ReportFirstframe(int tmp_time)
{
  printf ("%3d(I)  %8d %4d %7.4f %7.4f %7.4f  %5d       %3s \n",
          frame_no, stat->bit_ctr - stat->bit_ctr_n,
          img->qp, snr->snr_y, snr->snr_u, snr->snr_v, tmp_time,
          img->fld_flag ? "FLD" : "FRM");

  stat->bitr0 = stat->bitr;
  stat->bit_ctr_0 = stat->bit_ctr;
  stat->bit_ctr = 0;
}


static void ReportIntra(int tmp_time)
{

  printf ("%3d(I)  %8d %4d %7.4f %7.4f %7.4f  %5d       %3s \n",
          frame_no, stat->bit_ctr - stat->bit_ctr_n,
          img->qp, snr->snr_y, snr->snr_u, snr->snr_v, tmp_time,
          img->fld_flag ? "FLD" : "FRM");
}

static void ReportSP(int tmp_time)
{
  printf ("%3d(SP) %8d %4d %7.4f %7.4f %7.4f  %5d       %3s   %3d\n",
          frame_no, stat->bit_ctr - stat->bit_ctr_n, img->qp, snr->snr_y,
          snr->snr_u, snr->snr_v, tmp_time,
          img->fld_flag ? "FLD" : "FRM", intras);
}

static void ReportBS(int tmp_time)
{
  printf ("%3d(BS) %8d %4d %7.4f %7.4f %7.4f  %5d       %3s   %3d\n",
          frame_no, stat->bit_ctr - stat->bit_ctr_n, img->qp, snr->snr_y,
          snr->snr_u, snr->snr_v, tmp_time,
          img->fld_flag ? "FLD" : "FRM", intras);
}

static void ReportB(int tmp_time)
{
  printf ("%3d(B)  %8d %4d %7.4f %7.4f %7.4f  %5d       %3s \n",
          frame_no, stat->bit_ctr - stat->bit_ctr_n, img->qp,
          snr->snr_y, snr->snr_u, snr->snr_v, tmp_time,
          img->fld_flag ? "FLD" : "FRM");
}


static void ReportP(int tmp_time)
{            
  printf ("%3d(P)  %8d %4d %7.4f %7.4f %7.4f  %5d       %3s   %3d\n",
          frame_no, stat->bit_ctr - stat->bit_ctr_n, img->qp, snr->snr_y,
          snr->snr_u, snr->snr_v, tmp_time,
          img->fld_flag ? "FLD" : "FRM", intras);
}

/*!
 ************************************************************************
 * \brief
 *    Copies contents of a Sourceframe structure into the old-style
 *    variables imgY_org_frm and imgUV_org_frm.  No other side effects
 * \para
 *    sf: the source frame the frame is to be taken from
 ************************************************************************
 */

static void CopyFrameToOldImgOrgVariables (Sourceframe *sf)
{
  int x, y;

  for (y=0; y<sf->y_framesize; y++)
    for (x=0; x<sf->x_size; x++)
      imgY_org_frm [y][x] = sf->yf[y*sf->x_size+x];
  for (y=0; y<sf->y_framesize/2; y++)
    for (x=0; x<sf->x_size/2; x++)
    {
      imgUV_org_frm[0][y][x] = sf->uf[y*sf->x_size/2+x];
      imgUV_org_frm[1][y][x] = sf->vf[y*sf->x_size/2+x];
    }
}

/*!
 ************************************************************************
 * \brief
 *    Copies contents of a Sourceframe structure into the old-style
 *    variables imgY_org_top and imgUV_org_top.  No other side effects
 * \para
 *    sf: the source frame the field is to be taken from
 ************************************************************************
 */

static void CopyTopFieldToOldImgOrgVariables (Sourceframe *sf)
{
  int x, y;

  for (y=0; y<sf->y_fieldsize; y++)
    for (x=0; x<sf->x_size; x++)
      imgY_org_top [y][x] = sf->yt[y*sf->x_size+x];
  for (y=0; y<sf->y_fieldsize/2; y++)
    for (x=0;x<sf->x_size/2; x++)
    {
      imgUV_org_top[0][y][x] = sf->ut[y*sf->x_size/2+x];
      imgUV_org_top[1][y][x] = sf->vt[y*sf->x_size/2+x];
    }
}
/*!
 ************************************************************************
 * \brief
 *    Copies contents of a Sourceframe structure into the old-style
 *    variables imgY_org_bot and imgUV_org_bot.  No other side effects
 * \para
 *    sf: the source frame the field is to be taken from
 ************************************************************************
 */

static void CopyBottomFieldToOldImgOrgVariables (Sourceframe *sf)
{
  int x, y;

  for (y=0; y<sf->y_fieldsize; y++)
    for (x=0; x<sf->x_size; x++)
      imgY_org_bot [y][x] = sf->yb[y*sf->x_size+x];
  for (y=0; y<sf->y_fieldsize/2; y++)
    for (x=0;x<sf->x_size/2; x++)
    {
      imgUV_org_bot[0][y][x] = sf->ub[y*sf->x_size/2+x];
      imgUV_org_bot[1][y][x] = sf->vb[y*sf->x_size/2+x];
    }
}


/*!
 ************************************************************************
 * \brief
 *    Allocates Sourceframe structure
 * \para
 *    xs: horizontal size of frame in pixels
 *    ys: vertical size of frame in pixels, must be divisible by 2
 * \return
 *    pointer to initialized source frame structure
 ************************************************************************
 */

static Sourceframe *AllocSourceframe (int xs, int ys)
{
  Sourceframe *sf = NULL;
  const unsigned int bytes_y = xs*ys;
  const unsigned int bytes_uv = (xs*ys)/4;

  if ((sf = calloc (1, sizeof (Sourceframe))) == NULL) no_mem_exit ("ReadOneFrame: sf");
  if (sf->yf == NULL) if ((sf->yf = calloc (1, bytes_y)) == NULL) no_mem_exit ("ReadOneFrame: sf->yf");
  if (sf->yt == NULL) if ((sf->yt = calloc (1, bytes_y/2)) == NULL) no_mem_exit ("ReadOneFrame: sf->yt");
  if (sf->yb == NULL) if ((sf->yb = calloc (1, bytes_y/2)) == NULL) no_mem_exit ("ReadOneFrame: sf->yb");
  if (sf->uf == NULL) if ((sf->uf = calloc (1, bytes_uv)) == NULL) no_mem_exit ("ReadOneFrame: sf->uf");
  if (sf->ut == NULL) if ((sf->ut = calloc (1, bytes_uv/2)) == NULL) no_mem_exit ("ReadOneFrame: sf->ut");
  if (sf->ub == NULL) if ((sf->ub = calloc (1, bytes_uv/2)) == NULL) no_mem_exit ("ReadOneFrame: sf->ub");
  if (sf->vf == NULL) if ((sf->vf = calloc (1, bytes_uv)) == NULL) no_mem_exit ("ReadOneFrame: sf->vf");
  if (sf->vt == NULL) if ((sf->vt = calloc (1, bytes_uv/2)) == NULL) no_mem_exit ("ReadOneFrame: sf->vt");
  if (sf->vb == NULL) if ((sf->vb = calloc (1, bytes_uv/2)) == NULL) no_mem_exit ("ReadOneFrame: sf->vb");
  sf->x_size = xs;
  sf->y_framesize = ys;
  sf->y_fieldsize = ys/2;

  return sf;
}



/*!
 ************************************************************************
 * \brief
 *    Frees Sourceframe structure
 * \para
 *    pointer to Sourceframe previoously allocated with ALlocSourceframe()
 * \return
 *    none
 ************************************************************************
 */

static void FreeSourceframe (Sourceframe *sf)
{
  if (sf!=NULL) 
  {
    if (sf->yf != NULL) free (sf->yf);
    if (sf->yt != NULL) free (sf->yt);
    if (sf->yb != NULL) free (sf->yb);
    if (sf->uf != NULL) free (sf->uf);
    if (sf->ut != NULL) free (sf->ut);
    if (sf->ub != NULL) free (sf->ub);
    if (sf->vf != NULL) free (sf->vf);
    if (sf->vt != NULL) free (sf->vt);
    if (sf->vb != NULL) free (sf->vb);
    free (sf);
  }
}

/*!
 ************************************************************************
 * \brief
 *    Calculates the absolute frame number in the source file out
 *    of various variables in img-> and input->
 * \return
 *    frame number in the file to be read
 * \side effects
 *    global variable frame_no updated -- dunno, for what this one is necessary
 ************************************************************************
 */
static int CalculateFrameNumber()
{
  if (img->type == B_SLICE)
    frame_no = start_tr_in_this_IGOP + (IMG_NUMBER - 1) * (input->jumpd + 1) + img->b_interval * img->b_frame_to_code;
  else
    {
      frame_no = start_tr_in_this_IGOP + IMG_NUMBER * (input->jumpd + 1);
#ifdef _ADAPT_LAST_GROUP_
      if (input->last_frame && img->number + 1 == input->no_frames)
        frame_no = input->last_frame;
#endif
    }
  return frame_no;
}


/*!
 ************************************************************************
 * \brief
 *    Generate Field Component from Frame Components by copying
 * \para
 *    src: source frame component
 *    top: destination top field component
 *    bot: destination bottom field component
 *    xs: horizontal size of frame in pixels
 *    ys: vertical size of frame in pixels, must be divisible by 2
 ************************************************************************
 */
static void GenerateFieldComponent (char *src, char *top, char *bot, int xs, int ys)
{
  int fieldline;
  assert (ys % 2 == 0);

  for (fieldline = 0; fieldline < ys/2; fieldline++)
  {
    memcpy (&top[xs * fieldline], &src[xs * (fieldline * 2 + 0)], xs);
    memcpy (&bot[xs * fieldline], &src[xs * (fieldline * 2 + 1)], xs);
  }
}


/*!
 ************************************************************************
 * \brief
 *    Reads one new frame from file
 * \para
 *    FrameNoInFile: Frame number in the source file
 *    HeaderSize: Number of bytes in the source file to be skipped
 *    xs: horizontal size of frame in pixels, must be divisible by 16
 *    ys: vertical size of frame in pixels, must be divisible by 16 or
 *        32 in case of MB-adaptive frame/field coding
 *    sf: Sourceframe structure to which the frame is written
 ************************************************************************
 */
static void ReadOneFrame (int FrameNoInFile, int HeaderSize, int xs, int ys, Sourceframe *sf)
{
  int i;

  const unsigned int bytes_y = xs*ys;
  const unsigned int bytes_uv = (xs*ys)/4;
  const int framesize_in_bytes = bytes_y + 2*bytes_uv;

  assert (xs % MB_BLOCK_SIZE == 0);
  assert (ys % MB_BLOCK_SIZE == 0);
  assert (p_in != NULL);
  assert (sf != NULL);
  assert (sf->yf != NULL);

  assert (FrameNumberInFile == FrameNoInFile);
// printf ("ReadOneFrame: frame_no %d xs %d ys %d\n", FrameNoInFile, xs, ys);

  if (fseek (p_in, HeaderSize, SEEK_SET) != 0)
    error ("ReadOneFrame: cannot fseek to (Header size) in p_in", -1);

  // the reason for the following loop is to support source files bigger than
  // MAXINT.  In most operating systems, including Windows, it is possible to
  // fseek to file positions bigger than MAXINT by using this relative seeking
  // technique.  StW, 12/30/02
  for (i=0; i<FrameNoInFile; i++)
    if (fseek (p_in, framesize_in_bytes, SEEK_CUR) != 0) 
    {
      printf ("ReadOneFrame: cannot advance file pointer in p_in beyond frame %d, looping to picture zero\n", i);
      if (fseek (p_in, HeaderSize, SEEK_SET) != 0)
        error ("ReadOneFrame: cannot fseek to (Header size) in p_in", -1);
    }

  // Here we are at the correct position for the source frame in the file.  Now
  // read it.
  if (fread (sf->yf, 1, bytes_y, p_in) != bytes_y)
  {
    printf ("ReadOneFrame: cannot read %d bytes from input file, unexpected EOF?, exiting", bytes_y);
    exit (-1);
  }
  if (fread (sf->uf, 1, bytes_uv, p_in) != bytes_uv)
  {
    printf ("ReadOneFrame: cannot read %d bytes from input file, unexpected EOF?, exiting", bytes_uv);
    exit (-1);
  }
  if (fread (sf->vf, 1, bytes_uv, p_in) != bytes_uv)
  {
    printf ("ReadOneFrame: cannot read %d bytes from input file, unexpected EOF?, exiting", bytes_uv);
    exit (-1);
  }

  // Complete frame is read into sf->?f, now setup 
  // top and bottom field (sf->?t and sf->?b)

  GenerateFieldComponent (sf->yf, sf->yt, sf->yb, xs, ys);
  GenerateFieldComponent (sf->uf, sf->ut, sf->ub, xs/2, ys/2);
  GenerateFieldComponent (sf->vf, sf->vt, sf->vb, xs/2, ys/2);
}


/*!
 ************************************************************************
 * \brief
 *    point to frame coding variables 
 ************************************************************************
 */
static void put_buffer_frame()
{
  imgY_org = imgY_org_frm;
  imgUV_org = imgUV_org_frm;  
  tmp_mv = tmp_mv_frm;
  
	//mref = mref_frm;
  if (input->WeightedPrediction || input->WeightedBiprediction) 
    mref_w = mref_frm_w;
  //mcef = mcef_frm;  

  refFrArr = refFrArr_frm;
  fw_refFrArr = fw_refFrArr_frm;
  bw_refFrArr = bw_refFrArr_frm;

  //Refbuf11 = Refbuf11_frm;
  if (input->WeightedPrediction || input->WeightedBiprediction) 
        Refbuf11_w = Refbuf11_frm_w;
  if (input->direct_type && (input->successive_Bframe!=0 || input->StoredBPictures > 0))
    moving_block=moving_block_frm;
}

/*!
 ************************************************************************
 * \brief
 *    point to top field coding variables 
 ************************************************************************
 */
static void put_buffer_top()
{
  img->fld_type = 0;

  imgY_org = imgY_org_top;
  imgUV_org = imgUV_org_top;

  mref = mref_fld;
  if (input->WeightedPrediction || input->WeightedBiprediction) 
    mref_w = mref_fld_w;
  mcef = mcef_fld;

  Refbuf11 = Refbuf11_fld;  
  if (input->WeightedPrediction || input->WeightedBiprediction)
        Refbuf11_w = Refbuf11_fld_w;  
  tmp_mv = tmp_mv_top;
  refFrArr = refFrArr_top;
  fw_refFrArr = fw_refFrArr_top;
  bw_refFrArr = bw_refFrArr_top;

  if (input->direct_type && (input->successive_Bframe!=0 || input->StoredBPictures > 0))
    moving_block=moving_block_top;

}

/*!
 ************************************************************************
 * \brief
 *    point to bottom field coding variables 
 ************************************************************************
 */
static void put_buffer_bot()
{
  img->fld_type = 1;

  imgY_org = imgY_org_bot;
  imgUV_org = imgUV_org_bot;

  tmp_mv = tmp_mv_bot;
  refFrArr = refFrArr_bot;
  fw_refFrArr = fw_refFrArr_bot;
  bw_refFrArr = bw_refFrArr_bot;
  Refbuf11 = Refbuf11_fld;
  if (input->WeightedPrediction || input->WeightedBiprediction)
        Refbuf11_w = Refbuf11_fld_w;

  mref = mref_fld;
  if (input->WeightedPrediction || input->WeightedBiprediction) 
    mref_w = mref_fld_w;
  mcef = mcef_fld;
  if (input->direct_type && (input->successive_Bframe!=0 || input->StoredBPictures > 0))
    moving_block=moving_block_bot;
}

/*!
 ************************************************************************
 * \brief
 *    Writes a NAL unit of a partition or slice
 ************************************************************************
 */

static void writeUnit(Bitstream* currStream)
{
  NALU_t *nalu;
  assert (currStream->bits_to_go == 8);
  nalu = AllocNALU(img->width*img->height*4);
  nalu->startcodeprefix_len = 2+(img->current_mb_nr == 0?ZEROBYTES_SHORTSTARTCODE+1:ZEROBYTES_SHORTSTARTCODE);
//printf ("nalu->startcodeprefix_len %d\n", nalu->startcodeprefix_len);
  nalu->len = currStream->byte_pos +1;            // add one for the first byte of the NALU
//printf ("nalu->len %d\n", nalu->len);
  memcpy (&nalu->buf[1], currStream->streamBuffer, nalu->len-1);
  if (img->currentPicture->idr_flag)
  {
    nalu->nal_unit_type = NALU_TYPE_IDR;
    nalu->nal_reference_idc = NALU_PRIORITY_HIGHEST;
  }
  else if (img->type == B_SLICE)
  {
    nalu->nal_unit_type = NALU_TYPE_SLICE;
    nalu->nal_reference_idc = NALU_PRIORITY_DISPOSABLE;
  //  assert (img->nal_reference_idc != 0);
  }
  else   // non-disposable, non IDR slice
  {
    nalu->nal_unit_type = NALU_TYPE_SLICE;
    nalu->nal_reference_idc = NALU_PRIORITY_HIGH;
  }
  nalu->forbidden_bit = 0;
  stat->bit_ctr += WriteNALU (nalu);
  
  FreeNALU(nalu);
}


void get_block(int ref_frame, StorablePicture **list, int x_pos, int y_pos, int block[BLOCK_SIZE][BLOCK_SIZE])
{
  int dx, dy;
  int x, y;
  int i, j;
  int maxold_x,maxold_y;
  int result;
  int pres_x;
  int pres_y; 
  int tmp_res[4][9];
  static const int COEF[6] = {
			1, -5, 20, 20, -5, 1
		};

 
  dx = x_pos&3;
  dy = y_pos&3;
  x_pos = (x_pos-dx)/4;
  y_pos = (y_pos-dy)/4;

  maxold_x = img->width-1;
  maxold_y = img->height-1;

  if (enc_picture->mb_field[img->current_mb_nr])
    maxold_y = img->height/2 - 1;

if (dx == 0 && dy == 0) {  /* fullpel position */
    for (j = 0; j < BLOCK_SIZE; j++)
      for (i = 0; i < BLOCK_SIZE; i++)
        block[i][j] = list[ref_frame]->imgY[max(0,min(maxold_y,y_pos+j))][max(0,min(maxold_x,x_pos+i))];
  }
  else { /* other positions */

    if (dy == 0) { /* No vertical interpolation */

      for (j = 0; j < BLOCK_SIZE; j++) {
        for (i = 0; i < BLOCK_SIZE; i++) {
          for (result = 0, x = -2; x < 4; x++)
            result += list[ref_frame]->imgY[max(0,min(maxold_y,y_pos+j))][max(0,min(maxold_x,x_pos+i+x))]*COEF[x+2];
          block[i][j] = max(0, min(255, (result+16)/32));
        }
      }

      if ((dx&1) == 1) {
        for (j = 0; j < BLOCK_SIZE; j++)
          for (i = 0; i < BLOCK_SIZE; i++)
            block[i][j] = (block[i][j] + list[ref_frame]->imgY[max(0,min(maxold_y,y_pos+j))][max(0,min(maxold_x,x_pos+i+dx/2))] +1 )/2;
      }
    }
    else if (dx == 0) {  /* No horizontal interpolation */

      for (j = 0; j < BLOCK_SIZE; j++) {
        for (i = 0; i < BLOCK_SIZE; i++) {
          for (result = 0, y = -2; y < 4; y++)
            result += list[ref_frame]->imgY[max(0,min(maxold_y,y_pos+j+y))][max(0,min(maxold_x,x_pos+i))]*COEF[y+2];
          block[i][j] = max(0, min(255, (result+16)/32));
        }
      }

      if ((dy&1) == 1) {
        for (j = 0; j < BLOCK_SIZE; j++)
          for (i = 0; i < BLOCK_SIZE; i++)
           block[i][j] = (block[i][j] + list[ref_frame]->imgY[max(0,min(maxold_y,y_pos+j+dy/2))][max(0,min(maxold_x,x_pos+i))] +1 )/2;
      }
    }
    else if (dx == 2) {  /* Vertical & horizontal interpolation */

      for (j = -2; j < BLOCK_SIZE+3; j++) {
        for (i = 0; i < BLOCK_SIZE; i++)
          for (tmp_res[i][j+2] = 0, x = -2; x < 4; x++)
            tmp_res[i][j+2] += list[ref_frame]->imgY[max(0,min(maxold_y,y_pos+j))][max(0,min(maxold_x,x_pos+i+x))]*COEF[x+2];
      }

      for (j = 0; j < BLOCK_SIZE; j++) {
        for (i = 0; i < BLOCK_SIZE; i++) {
          for (result = 0, y = -2; y < 4; y++)
            result += tmp_res[i][j+y+2]*COEF[y+2];
          block[i][j] = max(0, min(255, (result+512)/1024));
        } 
      }

      if ((dy&1) == 1) {
        for (j = 0; j < BLOCK_SIZE; j++)
          for (i = 0; i < BLOCK_SIZE; i++)
            block[i][j] = (block[i][j] + max(0, min(255, (tmp_res[i][j+2+dy/2]+16)/32)) +1 )/2;
      }
    }
    else if (dy == 2) {  /* Horizontal & vertical interpolation */

      for (j = 0; j < BLOCK_SIZE; j++) {
        for (i = -2; i < BLOCK_SIZE+3; i++)
          for (tmp_res[j][i+2] = 0, y = -2; y < 4; y++)
            tmp_res[j][i+2] += list[ref_frame]->imgY[max(0,min(maxold_y,y_pos+j+y))][max(0,min(maxold_x,x_pos+i))]*COEF[y+2];
      }

      for (j = 0; j < BLOCK_SIZE; j++) {
        for (i = 0; i < BLOCK_SIZE; i++) {
          for (result = 0, x = -2; x < 4; x++)
            result += tmp_res[j][i+x+2]*COEF[x+2];
          block[i][j] = max(0, min(255, (result+512)/1024));
        }
      }

      if ((dx&1) == 1) {
        for (j = 0; j < BLOCK_SIZE; j++)
          for (i = 0; i < BLOCK_SIZE; i++)
            block[i][j] = (block[i][j] + max(0, min(255, (tmp_res[j][i+2+dx/2]+16)/32))+1)/2;
      }
    }
    else {  /* Diagonal interpolation */

      for (j = 0; j < BLOCK_SIZE; j++) {
        for (i = 0; i < BLOCK_SIZE; i++) {
          pres_y = dy == 1 ? y_pos+j : y_pos+j+1;
          pres_y = max(0,min(maxold_y,pres_y));
          for (result = 0, x = -2; x < 4; x++)
            result += list[ref_frame]->imgY[pres_y][max(0,min(maxold_x,x_pos+i+x))]*COEF[x+2];
          block[i][j] = max(0, min(255, (result+16)/32));
        }
      }

      for (j = 0; j < BLOCK_SIZE; j++) {
        for (i = 0; i < BLOCK_SIZE; i++) {
          pres_x = dx == 1 ? x_pos+i : x_pos+i+1;
          pres_x = max(0,min(maxold_x,pres_x));
          for (result = 0, y = -2; y < 4; y++)
            result += list[ref_frame]->imgY[max(0,min(maxold_y,y_pos+j+y))][pres_x]*COEF[y+2];
          block[i][j] = (block[i][j] + max(0, min(255, (result+16)/32)) +1 ) / 2;
        }
      }

    }
  }
}
