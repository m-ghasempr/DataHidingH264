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

#ifdef _ADAPT_LAST_GROUP_
int *last_P_no;
int *last_P_no_frm;
int *last_P_no_fld;
#endif

//! The followibng two variables are used for debug purposes.  They store
//! the status of the currStream data structure elements after the header
//! writing, and are used whether any MB bits were written during the
//! macroblock coding stage.  We need to worry about empty slices!
static int Byte_Pos_After_Header;
static int Bits_To_Go_After_Header;

static void UnifiedOneForthPix (pel_t **imgY, pel_t** imgU, pel_t **imgV,
        pel_t **out4Y, pel_t **outU, pel_t **outV, pel_t *ref11);

/*!
 ************************************************************************
 * \brief
 *    Encodes one frame
 ************************************************************************
 */
int encode_one_frame()
{
#ifdef _LEAKYBUCKET_
  extern long Bit_Buffer[10000];
  extern unsigned long total_frame_buffer;
#endif
  Bitstream *tmp_bitstream;
  int   i,j;

  time_t ltime1;   // for time measurement
  time_t ltime2;

#ifdef WIN32
  struct _timeb tstruct1;
  struct _timeb tstruct2;
#else
  struct timeb tstruct1;
  struct timeb tstruct2;
#endif

  int tmp_time;
  int  bits_frm=0, bits_fld=0;
  float dis_frm=0, dis_frm_y=0, dis_frm_u=0, dis_frm_v=0;
  float dis_fld=0, dis_fld_y=0, dis_fld_u=0, dis_fld_v=0;

#ifdef WIN32
  _ftime( &tstruct1 );    // start time ms
#else
  ftime(&tstruct1);
#endif
  time( &ltime1 );        // start time s


  if (input->InterlaceCodingOption==FIELD_CODING)
  {
    put_buffer_frame();  
    init_frame();
    read_one_new_frame();
    if (img->type==B_IMG)  Bframe_ctr++;

    field_picture(&bits_fld, &dis_fld_y, &dis_fld_u, &dis_fld_v);
    img->fld_flag = 1;
  }
  else 
  {
    
    // For frame coding, turn MB level field/frame coding flag on
    if(input->InterlaceCodingOption >= MB_CODING)
      mb_adaptive = 1;
    
    frame_picture(&bits_frm, &dis_frm_y, &dis_frm_u, &dis_frm_v);
    
    // For field coding, turn MB level field/frame coding flag off
    if(input->InterlaceCodingOption >= MB_CODING)
      mb_adaptive = 0;
    
    if(input->InterlaceCodingOption != FRAME_CODING)
    {
      field_picture(&bits_fld, &dis_fld_y, &dis_fld_u, &dis_fld_v);
      
      dis_fld = dis_fld_y + dis_fld_u + dis_fld_v;
      dis_frm = dis_frm_y + dis_frm_u + dis_frm_v;
      picture_structure_decision(bits_frm, bits_fld, dis_frm, dis_fld);
    }
    else 
      img->fld_flag = 0;
  }

  if (img->fld_flag)  stat->bit_ctr_emulationprevention += stat->em_prev_bits_fld;
  else                stat->bit_ctr_emulationprevention += stat->em_prev_bits_frm;



  if(img->type != B_IMG) {
    img->pstruct_next_P = img->fld_flag;
  }

  if (img->fld_flag)  // field mode (use field when fld_flag=1 only)
    field_mode_buffer(bits_fld, dis_fld_y, dis_fld_u, dis_fld_v);
  else  //frame mode
    frame_mode_buffer(bits_frm, dis_frm_y, dis_frm_u, dis_frm_v);
 
  if(input->of_mode==PAR_OF_26L)
    terminate_slice(1);

  if(input->InterlaceCodingOption != FRAME_CODING)
  {
    store_field_MV(IMG_NUMBER);  // assume that img->number = frame_number
    store_field_colB8mode();      // ABT: store B8-modes of collocated blocks for use in B-frames. mwi 020603
  }
  else
    store_direct_moving_flag(IMG_NUMBER);

  if(input->InterlaceCodingOption >= MB_CODING)
  {
    if (img->number!=0 && input->successive_Bframe!=0 && img->type != B_IMG)
    {
      for(i=0;i<img->height/(2*BLOCK_SIZE);i++)
        for(j=0;j<img->width/BLOCK_SIZE;j++)
        {
          refFrArr_top_save[i][j] = refFrArr_top[i][j];
          refFrArr_bot_save[i][j] = refFrArr_bot[i][j];
        }
      if (input->successive_Bframe!=0 || input->BipredictiveWeighting > 0)
      {
        for(i=0;i<img->height/(2*BLOCK_SIZE);i++)
          for(j=0;j<img->width/BLOCK_SIZE+4;j++)
          {
            tmp_mv_top_save[0][i][j] = tmp_mv_top[0][i][j];
            tmp_mv_bot_save[0][i][j] = tmp_mv_bot[0][i][j];
            tmp_mv_top_save[1][i][j] = tmp_mv_top[1][i][j];
            tmp_mv_bot_save[1][i][j] = tmp_mv_bot[1][i][j];
          }
      } 
      
    }
  }

  if(input->InterlaceCodingOption != FRAME_CODING)
  {
    if(!img->fld_flag) 
      tmp_bitstream = ((img->currentSlice)->partArr[0]).bitstream_fld;
    else
      tmp_bitstream = ((img->currentSlice)->partArr[0]).bitstream_frm;
        
    tmp_bitstream->stored_bits_to_go = 8;
    tmp_bitstream->stored_byte_pos = 0;
    tmp_bitstream->stored_byte_buf = 0;
  }

  if(input->InterlaceCodingOption == FRAME_CODING)
  {  
    if (input->rdopt==2 && img->type!=B_IMG)
      UpdateDecoders(); // simulate packet losses and move decoded image to reference buffers

    if (input->RestrictRef)
      UpdatePixelMap();
  }

  find_snr(snr,img);

  time(&ltime2);       // end time sec
#ifdef WIN32
  _ftime(&tstruct2);   // end time ms
#else
  ftime(&tstruct2);    // end time ms
#endif
  
  tmp_time=(ltime2*1000+tstruct2.millitm) - (ltime1*1000+tstruct1.millitm);
  tot_time=tot_time + tmp_time;


  // Write reconstructed images
  write_reconstructed_image();

#ifdef _LEAKYBUCKET_
  // Store bits used for this frame and increment counter of no. of coded frames
  Bit_Buffer[total_frame_buffer] = stat->bit_ctr - stat->bit_ctr_n;
  total_frame_buffer++;
#endif
  
  if(IMG_NUMBER == 0)
  {
    printf("%3d(I)  %8d %4d %7.4f %7.4f %7.4f  %5d       %3s \n",
        frame_no, stat->bit_ctr-stat->bit_ctr_n,
        img->qp, snr->snr_y, snr->snr_u, snr->snr_v, tmp_time, img->fld_flag?"FLD":"FRM");

    stat->bitr0=stat->bitr;
    stat->bit_ctr_0=stat->bit_ctr;
    stat->bit_ctr=0;
  }
  else
  {
    if (img->type == INTRA_IMG)
    {
      stat->bit_ctr_P += stat->bit_ctr-stat->bit_ctr_n;

      printf("%3d(I)  %8d %4d %7.4f %7.4f %7.4f  %5d       %3s \n",
        frame_no, stat->bit_ctr-stat->bit_ctr_n,
        img->qp, snr->snr_y, snr->snr_u, snr->snr_v, tmp_time, img->fld_flag?"FLD":"FRM");
    }
    else if (img->type != B_IMG)
    {
      stat->bit_ctr_P += stat->bit_ctr-stat->bit_ctr_n;
      if(img->types == SP_IMG)
        printf("%3d(SP) %8d %4d %7.4f %7.4f %7.4f  %5d       %3s   %3d\n",
          frame_no, stat->bit_ctr-stat->bit_ctr_n,
          img->qp, snr->snr_y, snr->snr_u, snr->snr_v, tmp_time, img->fld_flag?"FLD":"FRM", intras);
      else if(img->type == BS_IMG)
        printf("%3d(BS) %8d %4d %7.4f %7.4f %7.4f  %5d       %3s   %3d\n",
          frame_no, stat->bit_ctr-stat->bit_ctr_n,
          img->qp, snr->snr_y, snr->snr_u, snr->snr_v, tmp_time, img->fld_flag?"FLD":"FRM", intras);
      else
        printf("%3d(P)  %8d %4d %7.4f %7.4f %7.4f  %5d       %3s   %3d\n",
          frame_no, stat->bit_ctr-stat->bit_ctr_n,
          img->qp, snr->snr_y, snr->snr_u, snr->snr_v, tmp_time, img->fld_flag?"FLD":"FRM", intras);
    }
    else
    {
      stat->bit_ctr_B += stat->bit_ctr-stat->bit_ctr_n;


      printf("%3d(B)  %8d %4d %7.4f %7.4f %7.4f  %5d       %3s \n",
        frame_no, stat->bit_ctr-stat->bit_ctr_n, img->qp,
        snr->snr_y, snr->snr_u, snr->snr_v, tmp_time, img->fld_flag?"FLD":"FRM");
    }
  }
  stat->bit_ctr_n=stat->bit_ctr;

  // for Interim File Format output
  if ( input->of_mode == PAR_OF_IFF )
  {
    wrPictureInfo( box_pi.fpMeta );  // write the picture info to a temporary file
    freePictureInfo();
  }

  if(IMG_NUMBER == 0)
    return 0;
  else
    return 1;
}

void frame_picture(int *bits_frm, float *dis_frm_y, float *dis_frm_u, float *dis_frm_v)
{
  int i,j,k,l;

  int SliceGroup = 0;
  int NumberOfCodedMBs;

  if(input->InterlaceCodingOption != FRAME_CODING) 
  {
    if(mb_adaptive)
    {
     if (img->type == B_IMG)  // if more than one B pictures, they will overwrite refFrArr_top, refFrArr_bot
     {
       for(i=0;i<img->height/(2*BLOCK_SIZE);i++)
         for(j=0;j<img->width/BLOCK_SIZE;j++)
         {
           refFrArr_top[i][j] = refFrArr_top_save[i][j];  
           refFrArr_bot[i][j] = refFrArr_bot_save[i][j];
         }
        if (input->successive_Bframe!=0 || input->BipredictiveWeighting > 0)
        {
         for(i=0;i<img->height/(2*BLOCK_SIZE);i++)
           for(j=0;j<img->width/BLOCK_SIZE+4;j++)
           {
             tmp_mv_top[0][i][j] = tmp_mv_top_save[0][i][j];
             tmp_mv_bot[0][i][j] = tmp_mv_bot_save[0][i][j];
             tmp_mv_top[1][i][j] = tmp_mv_top_save[1][i][j];
             tmp_mv_bot[1][i][j] = tmp_mv_bot_save[1][i][j];

           }
        }

     }
    }
  }

  RandomIntraNewPicture();    //! Tells the RandomIntra that a new picture has started, no side effects
  put_buffer_frame();  
  stat->em_prev_bits_frm = 0;
  stat->em_prev_bits     = &stat->em_prev_bits_frm;

  // Initialize frame with all stat and img variables
  img->total_number_mb = (img->width * img->height)/(MB_BLOCK_SIZE*MB_BLOCK_SIZE);
  init_frame();

  init_mref(img);
  init_Refbuf(img);
  copy_mref(img);


  // Read one new frame
  read_one_new_frame();

  if (img->type == B_IMG)
    Bframe_ctr++;

  // for Interim File Format output
  if ( input->of_mode == PAR_OF_IFF )
  {
    box_pi.numPictures++;  // update header of AlternateTrackHeaderClump
    initPictureInfo();
  }

  for (i=0; i<img->currentSlice->max_part_nr; i++)
    ((img->currentSlice)->partArr[i]).bitstream = ((img->currentSlice)->partArr[i]).bitstream_frm; //temp

#if _DEBUG  //set temp_streamBuffer, limin 01/02/02
  tmp_streamBuffer = (img->currentSlice->partArr[0].bitstream_frm)->streamBuffer;
  pbits_to_go = &(img->currentSlice->partArr[0].bitstream_frm->bits_to_go);       
  pbyte_buf = &(img->currentSlice->partArr[0].bitstream_frm->byte_buf); 
  pbyte_pos = &(img->currentSlice->partArr[0].bitstream_frm->byte_pos); 
#endif

  // The difference betrween the following loop and the original code is
  // that the end of the picture is determined by summing up the number
  // of coded mBs, rather than rely on the property that the last coded
  // MB is the spatially last MB of the picture (which is not true in FMO)

  FmoStartPicture();   //! picture level initialization of FMO
  NumberOfCodedMBs = 0;
  SliceGroup = 0;

  while (NumberOfCodedMBs < img->total_number_mb) // loop over slices
  {
    // Encode one SLice Group
    while (!FmoSliceGroupCompletelyCoded (SliceGroup))
    {
      // Encode the current slice
      NumberOfCodedMBs += encode_one_slice(SliceGroup);
      FmoSetLastMacroblockInSlice (img->current_mb_nr);

    // Proceed to next slice
    img->current_slice_nr++;
    stat->bit_slice = 0;
    }
    // Proceed to next SliceGroup
    SliceGroup++;
  }
  FmoEndPicture();
  if (input->rdopt==2 && (img->type!=B_IMG) )
    for (j=0 ;j<input->NoOfDecoders; j++)  DeblockFrame(img, decs->decY_best[j], NULL ) ;

  // CAVLC init
  for (i=0;i < img->width/MB_BLOCK_SIZE; i++)
    for (j=0; j < img->height/MB_BLOCK_SIZE; j++)
      for (k=0;k<4;k++)
        for (l=0;l<6;l++)
          img->nz_coeff[i][j][k][l]=-1;


  if(input->InterlaceCodingOption >= MB_CODING)
  {
    if(mb_adaptive)
     if (img->number!=0 && input->successive_Bframe!=0 && img->type != B_IMG)
    {
       for(i=0;i<img->height/(2*BLOCK_SIZE);i++)
         for(j=0;j<img->width/BLOCK_SIZE;j++)
         {
           refFrArr_top_save[i][j] = refFrArr_top[i][j];
           refFrArr_bot_save[i][j] = refFrArr_bot[i][j];
         }
        if (input->successive_Bframe!=0 || input->BipredictiveWeighting > 0)
        {
           for(i=0;i<img->height/(2*BLOCK_SIZE);i++)
           for(j=0;j<img->width/BLOCK_SIZE+4;j++)
           {
             tmp_mv_top_save[0][i][j] = tmp_mv_top[0][i][j];
             tmp_mv_bot_save[0][i][j] = tmp_mv_bot[0][i][j];
             tmp_mv_top_save[1][i][j] = tmp_mv_top[1][i][j];
             tmp_mv_bot_save[1][i][j] = tmp_mv_bot[1][i][j];
           }
       }  

    }
  }

  DeblockFrame( img, imgY, imgUV ) ;

  *bits_frm = 8*((((img->currentSlice)->partArr[0]).bitstream)->byte_pos);

  if(input->InterlaceCodingOption != FRAME_CODING)
  {
    find_distortion(snr,img);       //distortion, not snr
    *dis_frm_y = snr->snr_y;
    *dis_frm_u = snr->snr_u;
    *dis_frm_v = snr->snr_v;
  }
}

void field_picture(int *bits_fld, float *dis_fld_y, float *dis_fld_u, float *dis_fld_v)
{
  stat->em_prev_bits_fld = 0;
  stat->em_prev_bits     = &stat->em_prev_bits_fld;

  top_field_picture(bits_fld);
  bottom_field_picture(bits_fld);
  distortion_fld(dis_fld_y, dis_fld_u, dis_fld_v);
}

void top_field_picture(int *bits_fld)
{
  int i,j,k,l;
  int NumberOfCodedMBs, SliceGroup;

  img->number *= 2;
  img->buf_cycle *= 2;
  input->no_multpred = 2*input->no_multpred+1;
  img->height = input->img_height/2;
  img->height_cr = input->img_height/4;
  img->total_number_mb = (img->width * img->height)/(MB_BLOCK_SIZE*MB_BLOCK_SIZE);

  put_buffer_top();

  // Initialize field with all stat and img variables
  for (i=0; i<img->currentSlice->max_part_nr; i++)
    ((img->currentSlice)->partArr[i]).bitstream = ((img->currentSlice)->partArr[i]).bitstream_fld; //temp

#if _DEBUG  //set temp_streamBuffer, limin 01/02/02
  tmp_streamBuffer = (img->currentSlice->partArr[0].bitstream_fld)->streamBuffer;
  pbits_to_go = &(img->currentSlice->partArr[0].bitstream_fld->bits_to_go);       
  pbyte_buf = &(img->currentSlice->partArr[0].bitstream_fld->byte_buf); 
  pbyte_pos = &(img->currentSlice->partArr[0].bitstream_fld->byte_pos); 
#endif

  init_field();

  if (img->type == B_IMG) //all I- and P-frames
    nextP_tr_fld--;

  read_one_new_field();  // top field

  // The difference betrween the following loop and the original code is
  // that the end of the picture is determined by summing up the number
  // of coded mBs, rather than rely on the property that the last coded
  // MB is the spatially last MB of the picture (which is not true in FMO)

  FmoStartPicture();   //! picture level initialization of FMO
  NumberOfCodedMBs = 0;
  SliceGroup = 0;

  while (NumberOfCodedMBs < img->total_number_mb) // loop over slices
  {
    // Encode one SLice Group
    while (!FmoSliceGroupCompletelyCoded (SliceGroup))
    {
      // Encode the current slice
      NumberOfCodedMBs += encode_one_slice(SliceGroup);
      FmoSetLastMacroblockInSlice (img->current_mb_nr);

    // Proceed to next slice
    img->current_slice_nr++;
    stat->bit_slice = 0;
    }
    // Proceed to next SliceGroup
    SliceGroup++;
  }
  FmoEndPicture();

  // CAVLC init
  for (i=0;i < img->width/MB_BLOCK_SIZE; i++)
    for (j=0; j < img->height/MB_BLOCK_SIZE; j++)
      for (k=0;k<4;k++)
        for (l=0;l<6;l++)
          img->nz_coeff[i][j][k][l]=-1;


  if (input->rdopt==2 && (img->type!=B_IMG) )
    for (j=0 ;j<input->NoOfDecoders; j++)  DeblockFrame(img, decs->decY_best[j], NULL ) ;


  DeblockFrame( img, imgY, imgUV ) ;
  if (img->type != B_IMG)               //all I- and P-frames
    interpolate_frame_to_fb();

  *bits_fld = 8*((((img->currentSlice)->partArr[0]).bitstream)->byte_pos);

}

void bottom_field_picture(int *bits_fld)
{
  int j; 

  Bitstream *tmp_bitstream;
  int bits_top, bit_bot;
  int NumberOfCodedMbs, SliceGroup;
  
  put_buffer_bot();

  img->number++;  
  bits_top = *bits_fld;

  init_field();

  if (img->type == B_IMG) //all I- and P-frames
    nextP_tr_fld++;       //check once coding B field

  if (img->type == INTRA_IMG)
    img->type = INTER_IMG;

  tmp_bitstream = ((img->currentSlice)->partArr[0]).bitstream;
  tmp_bitstream->stored_bits_to_go = tmp_bitstream->bits_to_go;
  tmp_bitstream->stored_byte_pos = tmp_bitstream->byte_pos;
  tmp_bitstream->stored_byte_buf = tmp_bitstream->byte_buf;

  read_one_new_field();  // bottom field, use frame number for reading from frame buffer

  // The difference betrween the following loop and the original code is
  // that the end of the picture is determined by summing up the number
  // of coded mBs, rather than rely on the property that the last coded
  // MB is the spatially last MB of the picture (which is not true in FMO)

  FmoStartPicture();   //! picture level initialization of FMO
  NumberOfCodedMbs = 0;
  SliceGroup = 0;

  while (NumberOfCodedMbs < img->total_number_mb) // loop over slices
  {
    // Encode one SLice Group
    while (!FmoSliceGroupCompletelyCoded (SliceGroup))
    {
      // Encode the current slice
      NumberOfCodedMbs += encode_one_slice(SliceGroup);
      FmoSetLastMacroblockInSlice (img->current_mb_nr);

    // Proceed to next slice
    img->current_slice_nr++;
    stat->bit_slice = 0;
    }
    // Proceed to next SliceGroup
    SliceGroup++;
  }
  FmoEndPicture();

  if (input->rdopt==2 && (img->type!=B_IMG) )
    for (j=0 ;j<input->NoOfDecoders; j++)  DeblockFrame(img, decs->decY_best[j], NULL ) ;


  DeblockFrame( img, imgY, imgUV ) ;
  if (img->type != B_IMG) //all I- and P-frames
    interpolate_frame_to_fb();

  *bits_fld = 8*((((img->currentSlice)->partArr[0]).bitstream)->byte_pos);
  bit_bot = *bits_fld - bits_top;
}

void distortion_fld(float *dis_fld_y, float *dis_fld_u, float *dis_fld_v)
{

  img->number /= 2;  
  img->buf_cycle /= 2;
  img->height = input->img_height;
  img->height_cr = input->img_height/2;
  img->total_number_mb = (img->width * img->height)/(MB_BLOCK_SIZE*MB_BLOCK_SIZE);
  input->no_multpred = (input->no_multpred-1)/2;

  combine_field();

  imgY = imgY_com; 
  imgUV = imgUV_com;
  imgY_org = imgY_org_frm;
  imgUV_org = imgUV_org_frm;
  
  find_distortion(snr,img);  // find snr from original frame picture
  
  *dis_fld_y = snr->snr_y;
  *dis_fld_u = snr->snr_u;
  *dis_fld_v = snr->snr_v;
}

void picture_structure_decision(int bit_frame, int bit_field, float snr_frame, float snr_field)
{
  double lambda_picture;
  int         spframe     = (img->type==INTER_IMG && img->types==SP_IMG);
  int         bframe      = (img->type==B_IMG);

  lambda_picture = 0.85 * pow (2, (img->qp-SHIFT_QP)/3.0) * (bframe||spframe?4:1); 

  img->fld_flag = decide_fld_frame(snr_frame, snr_field, bit_field, bit_frame, lambda_picture);    // modify fld_flag accordingly
}

void field_mode_buffer(int bit_field, float snr_field_y, float snr_field_u, float snr_field_v)
{

  int i;

  for (i=0; i<img->currentSlice->max_part_nr; i++)
    ((img->currentSlice)->partArr[i]).bitstream = ((img->currentSlice)->partArr[i]).bitstream_fld; //temp

  put_buffer_frame();
  imgY = imgY_com; 
  imgUV = imgUV_com;
  
  if (img->type != B_IMG) //all I- and P-frames   //    Krit: PLUS2, Change field buffer structure, 7/2
    interpolate_frame_to_fb();
    
  input->no_fields += 1;

  snr->snr_y = snr_field_y;
  snr->snr_u = snr_field_u;
  snr->snr_v = snr_field_v;
}

void frame_mode_buffer(int bit_frame, float snr_frame_y, float snr_frame_u, float snr_frame_v)
{
  int i;
  
  for (i=0; i<img->currentSlice->max_part_nr; i++)
    ((img->currentSlice)->partArr[i]).bitstream = ((img->currentSlice)->partArr[i]).bitstream_frm; //temp

  if(input->InterlaceCodingOption!=FRAME_CODING)
    if (input->symbol_mode == CABAC)
      img->current_mb_nr=img->current_mb_nr*2;

  put_buffer_frame();

  if (img->type != B_IMG) //all I- and P-frames
      interpolate_frame_to_fb();

  if(input->InterlaceCodingOption != FRAME_CODING)
  {
    img->height = img->height/2;
    img->height_cr = img->height_cr/2;
    img->number *= 2;
    img->buf_cycle *= 2;

    put_buffer_top();
    split_field_top();

    if (img->type != B_IMG) //all I- and P-frames
    {
      rotate_buffer();
      interpolate_frame();
    }
    img->number++;
    put_buffer_bot();
    split_field_bot();

    if (img->type != B_IMG) //all I- and P-frames
    {
      rotate_buffer();
      interpolate_frame();
    }
  
    img->number /= 2;  // reset the img->number to field
    img->buf_cycle /= 2;
    img->height = input->img_height;
    img->height_cr = input->img_height/2;
    img->total_number_mb = (img->width * img->height)/(MB_BLOCK_SIZE*MB_BLOCK_SIZE);

    snr->snr_y = snr_frame_y;
    snr->snr_u = snr_frame_u;
    snr->snr_v = snr_frame_v;
    put_buffer_frame();

  }
}


//! This is the old encode_one_slice.  It is left here to help those who
//! believe FMO has introduced a bug here to have a quick look at the old
//! code.  I (StW) do not believe that the new code is sriously buggy,
//! but one can never be sure
//!
//! PLEASE DO NOT simply comment in this function and comment out the
//! newer version below.  It will not only break FMO but also the rest
//! of the code.
//!
//! Karsten, please feel free to delete the old code and the comment
//! any time you deem appropriate (IMHO it is certainly irrelevant after 
//! Klagenfurt)
#if 0

/*!
 ************************************************************************
 * \brief
 *    Encodes one slice
 ************************************************************************
 */
void encode_one_slice(SyntaxElement *sym)
{
  Boolean end_of_slice = FALSE;
  Boolean recode_macroblock;
  int len;
  int short_used = 0, img_ref = 0;

  img->cod_counter=0;

  // Initializes the parameters of the current slice
  init_slice();
  if(input->Encapsulated_NAL_Payload) //stores the position after header. Needed for skipping byte stuffing in RTP Header
  {
    Bytes_After_Header = img->currentSlice->partArr[0].bitstream->byte_pos;
  }

  // Tian Dong: June 7, 2002 JVT-B042
  // When the pictures are put into different layers and subseq, not all the reference frames
  // in multi-frame buffer are valid for prediction. The acutual number of the valid reference
  // frames, fb->num_short_used, will be given by start_slice(sym).
  // Save the fb->short_used.
  if ( input->NumFramesInELSubSeq )
  {
    short_used = fb->short_used;
    img_ref = img->nb_references;
  }

  // Write slice or picture header
  len = start_slice(sym);

//  printf("short size, used, num-used: (%d,%d,%d)\n", fb->short_size, fb->short_used, fb->num_short_used);

  // Tian Dong: June 7, 2002 JVT-B042
  if ( input->NumFramesInELSubSeq )
  {
    fb->short_used = fb->num_short_used;
    img->nb_references = fb->short_used+fb->long_used;
  }

  Byte_Pos_After_Header = img->currentSlice->partArr[0].bitstream->byte_pos;
  Bits_To_Go_After_Header = img->currentSlice->partArr[0].bitstream->bits_to_go;

  if (input->of_mode==PAR_OF_RTP)
  {
    assert (Byte_Pos_After_Header > 0);     // there must be a header
    assert (Bits_To_Go_After_Header == 8);  // byte alignment must have been established
  }

  // Update statistics
  stat->bit_slice += len;
  stat->bit_use_header[img->type] += len;

  while (end_of_slice == FALSE) // loop over macroblocks
  {
    // recode_macroblock is used in slice mode two and three where
    // backing of one macroblock in the bitstream is possible
    recode_macroblock = FALSE;

    // Initializes the current macroblock
    start_macroblock();

    // Encode one macroblock
    encode_one_macroblock();

    // Pass the generated syntax elements to the NAL
    write_one_macroblock();

    // Terminate processing of the current macroblock
    terminate_macroblock(&end_of_slice, &recode_macroblock);

    if (recode_macroblock == FALSE)         // The final processing of the macroblock has been done
      proceed2nextMacroblock(); // Go to next macroblock

  }
  
  // Tian Dong: June 7, 2002 JVT-B042
  // Restore the short_used
  if ( input->NumFramesInELSubSeq )
  {
    fb->short_used = short_used;
    img->nb_references = img_ref;
  }

  terminate_slice(0);
}

#endif

/*!
 ************************************************************************
 * \brief
 *    Encodes one slice
 * \para
 *   returns the number of coded MBs in the SLice 
 ************************************************************************
 */

int encode_one_slice(int SliceGroupId)
{
  Boolean end_of_slice = FALSE;
  Boolean recode_macroblock;
  int len;
  int NumberOfCodedMBs = 0;
  int CurrentMbInScanOrder;
  int short_used = 0, img_ref = 0;
  int MBRowSize = img->width/MB_BLOCK_SIZE;
  double FrameRDCost, FieldRDCost;


  img->cod_counter=0;

  // Initializes the parameters of the current slice
  CurrentMbInScanOrder = FmoGetFirstMacroblockInSlice (SliceGroupId);
// printf ("\n\nEncode_one_slice: PictureID %d SliceGroupId %d  SliceID %d  FirstMB %d \n", img->tr, SliceGroupId, img->current_slice_nr, CurrentMbInScanOrder);

  set_MB_parameters (CurrentMbInScanOrder);
  init_slice(CurrentMbInScanOrder);
  if(input->Encapsulated_NAL_Payload) //stores the position after header. Needed for skipping byte stuffing in RTP Header
  {
    Bytes_After_Header = img->currentSlice->partArr[0].bitstream->byte_pos;
  }

  // Tian Dong: June 7, 2002 JVT-B042
  // When the pictures are put into different layers and subseq, not all the reference frames
  // in multi-frame buffer are valid for prediction. The acutual number of the valid reference
  // frames, fb->num_short_used, will be given by start_slice(sym).
  // Save the fb->short_used.
  if ( input->NumFramesInELSubSeq )
  {
    short_used = fb->short_used;
    img_ref = img->nb_references;
  }

  // Write slice or picture header
  len = start_slice();
//  printf("short size, used, num-used: (%d,%d,%d)\n", fb->short_size, fb->short_used, fb->num_short_used);

  // Tian Dong: June 7, 2002 JVT-B042
  if ( input->NumFramesInELSubSeq )
  {
    fb->short_used = fb->num_short_used;
    img->nb_references = fb->short_used+fb->long_used;
  }

  Byte_Pos_After_Header = img->currentSlice->partArr[0].bitstream->byte_pos;
  Bits_To_Go_After_Header = img->currentSlice->partArr[0].bitstream->bits_to_go;

  // Update statistics
  stat->bit_slice += len;
  stat->bit_use_header[img->type] += len;
// printf ("\n\n");

  while (end_of_slice == FALSE) // loop over macroblocks
  {
  if(input->InterlaceCodingOption < MB_CODING || mb_adaptive == 0)
  {
    // recode_macroblock is used in slice mode two and three where
    // backing of one macroblock in the bitstream is possible
    recode_macroblock = FALSE;
    set_MB_parameters (CurrentMbInScanOrder);


    // Initializes the current macroblock
    start_macroblock();

    // Encode one macroblock
    encode_one_macroblock();

    // Pass the generated syntax elements to the NAL
    write_one_macroblock();

    // Terminate processing of the current macroblock
    terminate_macroblock(&end_of_slice, &recode_macroblock);

    // printf ("encode_one_slice: mb %d,  slice %d,   bitbuf bytepos %d EOS %d\n", 
    //       img->current_mb_nr, img->current_slice_nr, 
    //       img->currentSlice->partArr[0].bitstream->byte_pos, end_of_slice);

    if (recode_macroblock == FALSE)         // The final processing of the macroblock has been done
    {
      CurrentMbInScanOrder = FmoGetNextMBNr (CurrentMbInScanOrder);
      if (CurrentMbInScanOrder == -1)   // end of slice
      {
    // printf ("FMO End of Slice Group detected, current MBs %d, force end of slice\n", NumberOfCodedMBs+1);
      end_of_slice = TRUE;
      }
      NumberOfCodedMBs++;       // only here we are sure that the coded MB is actually included in the slice
      proceed2nextMacroblock(CurrentMbInScanOrder); // Go to next macroblock
    } else
    {
      img->current_mb_nr--; /*KS*/
    }
  }
  else    // TBD -- Addition of FMO
  {
    // recode_macroblock is used in slice mode two and three where
    // backing of one macroblock in the bitstream is possible
    // code MB pair as frame MB 
    recode_macroblock = FALSE;
    img->field_mode   = 0;    // MB coded as frame
    img->top_field    = 0;    // Set top field to 0
    set_MB_parameters (CurrentMbInScanOrder);


    // Initializes the current macroblock
    start_macroblock();
    img->update_stats = 0;    // don't update any stats yet
    img->height     = input->img_height; // set image height as frame height
    rdopt       = &rddata_top_frame_mb; // store data in top frame MB 
    mcef        = mcef_frm; // set mcef to mcef_frm
    mref        =   mref_frm; // set mref to mref_frm
    Refbuf11      =   Refbuf11_frm; // set Refbuff to frm
    TopFrameIsSkipped =   0;
    WriteFrameFieldMBInHeader = 1;
    encode_one_macroblock();  // code the MB as frame
    field_mb[img->mb_y][img->mb_x] = 0;  // set the MB as field (for use in FindSkipMotionVector)
    FrameRDCost     =   rdopt->min_rdcost;
    //***   Top MB coded as frame MB ***//
    
    // go to the bottom MB in the MB pair
    CurrentMbInScanOrder =  img->current_mb_nr+MBRowSize;  
    set_MB_parameters (CurrentMbInScanOrder);
    start_macroblock();
    img->update_stats = 0;    // don't update any stats yet
    img->field_mode   = 0;    // MB coded as frame
    rdopt       = &rddata_bot_frame_mb; // store data in top frame MB
    WriteFrameFieldMBInHeader = TopFrameIsSkipped ? 1:0;
    field_mb[img->mb_y][img->mb_x] = 0;
    encode_one_macroblock();  // code the MB as frame
    field_mb[img->mb_y][img->mb_x] = 0;  // set the MB as field (for use in FindSkipMotionVector)
    FrameRDCost     +=   rdopt->min_rdcost;

    //***   Bottom MB coded as frame MB ***//


    // start coding the MB pair as a field MB pair
    CurrentMbInScanOrder -= MBRowSize;
    img->field_mode   =   1;    // MB coded as frame
    img->top_field    =   1;    // Set top field to 1
    set_MB_parameters (CurrentMbInScanOrder);
    img->buf_cycle <<=1;input->no_multpred <<= 1;
    img->num_ref_pic_active_fwd_minus1 <<= 1;img->num_ref_pic_active_fwd_minus1+=1;
    start_macroblock();



    img->height     =   input->img_height >> 1; // set image height as frame height
    rdopt       =   &rddata_top_field_mb; // store data in top frame MB 
    mcef        =   mcef_fld; // set mcef to mcef_frm
    mref        =   mref_mbfld; // set mref to mref_frm
    Refbuf11      =   Refbuf11_fld; // set Refbuff to frm
    TopFieldIsSkipped =   0;      // set the top field MB skipped flag to 0
    WriteFrameFieldMBInHeader = 1;
    encode_one_macroblock();  // code the MB as frame
    field_mb[img->mb_y][img->mb_x] = 1;  // set the MB as field (for use in FindSkipMotionVector)
    FieldRDCost     =   rdopt->min_rdcost;
    //***   Top MB coded as field MB ***//

    CurrentMbInScanOrder += MBRowSize;
    img->top_field    =   0;    // Set top field to 0
    set_MB_parameters (CurrentMbInScanOrder);
    start_macroblock();
    rdopt       =   &rddata_bot_field_mb; // store data in top frame MB 
    mcef        =   mcef_fld; // set mcef to mcef_frm
    mref        =   mref_mbfld; // set mref to mref_frm
    Refbuf11      =   Refbuf11_fld; // set Refbuff to frm
    WriteFrameFieldMBInHeader = TopFieldIsSkipped ? 1:0;
    encode_one_macroblock();  // code the MB as frame
    field_mb[img->mb_y][img->mb_x] = 1;  // set the MB as field (for use in FindSkipMotionVector)
    FieldRDCost     +=   rdopt->min_rdcost;
    //***   Bottom MB coded as field MB ***//

    // decide between frame/field MB pair
    if(FrameRDCost < FieldRDCost)
    {
      img->field_mode   =  0;
      img->buf_cycle    >>=1;
      input->no_multpred  >>= 1;
      MBPairIsField   = 0;
      mcef        = mcef_frm;
      mref        =   mref_frm;
      Refbuf11      =   Refbuf11_frm;
      img->height     = input->img_height; // set image height as frame height
      img->num_ref_pic_active_fwd_minus1-=1;
      img->num_ref_pic_active_fwd_minus1 >>= 1;


    }
    else
    {
      img->field_mode   =  1;
      MBPairIsField   = 1;
      mcef        = mcef_fld;
      mref        =   mref_mbfld;
      Refbuf11      =   Refbuf11_fld;
      img->height     = input->img_height/2; // set image height as frame height
    }

    if(MBPairIsField)
      img->top_field = 1;
    else
      img->top_field = 0;

    // go back to the Top MB in the MB pair
    CurrentMbInScanOrder -= MBRowSize;
    set_MB_parameters (CurrentMbInScanOrder);
    start_macroblock();   
    img->update_stats   = 1;    // Now update the stats
    rdopt         =   img->field_mode ? &rddata_top_field_mb:&rddata_top_frame_mb;
    copy_rdopt_data(0);       // copy the MB data for Top MB from the temp buffers
    WriteFrameFieldMBInHeader = 1;
    write_one_macroblock();     // write the Top MB data to the bitstream
    NumberOfCodedMBs++;       // only here we are sure that the coded MB is actually included in the slice
    terminate_macroblock(&end_of_slice, &recode_macroblock);  // done coding the Top MB 
      proceed2nextMacroblock(CurrentMbInScanOrder); // Go to next macroblock

    // go to the Bottom MB in the MB pair
    CurrentMbInScanOrder += MBRowSize;
    img->top_field      = 0;
    set_MB_parameters (CurrentMbInScanOrder);
    start_macroblock();   
    img->update_stats   = 1;    // Now update the stats
    rdopt         =   img->field_mode ? &rddata_bot_field_mb:&rddata_bot_frame_mb;
    copy_rdopt_data(1); // copy the MB data for Bottom MB from the temp buffers
    if(img->field_mode)
      WriteFrameFieldMBInHeader = TopFieldIsSkipped ? 1:0;
    else
      WriteFrameFieldMBInHeader = TopFrameIsSkipped ? 1:0;

    write_one_macroblock();     // write the Bottom MB data to the bitstream
    NumberOfCodedMBs++;       // only here we are sure that the coded MB is actually included in the slice
    terminate_macroblock(&end_of_slice, &recode_macroblock);  // done coding the Top MB 
      proceed2nextMacroblock(CurrentMbInScanOrder); // Go to next macroblock
    
    CurrentMbInScanOrder -= MBRowSize;

    if(MBPairIsField)       // if MB Pair was coded as field the buffer size variables back to frame mode
    {
      img->buf_cycle >>=1;
      input->no_multpred >>= 1;
      img->num_ref_pic_active_fwd_minus1-=1;
      img->num_ref_pic_active_fwd_minus1 >>= 1;
    }
    img->field_mode = img->top_field = 0; // reset to frame mode
    img->height     = input->img_height;  // reset the img->height  
    
    CurrentMbInScanOrder++;
    if(CurrentMbInScanOrder == img->total_number_mb-MBRowSize)
      end_of_slice = TRUE;    // just in case it does n't get set in terminate_macroblock

    if(CurrentMbInScanOrder%MBRowSize == 0)
      CurrentMbInScanOrder +=MBRowSize;
  }

  }

  // Tian Dong: June 7, 2002 JVT-B042
  // Restore the short_used
  if ( input->NumFramesInELSubSeq )
  {
    fb->short_used = short_used;
    img->nb_references = img_ref;
  }

  terminate_slice(0);
  return NumberOfCodedMBs;
}



/*!
 ************************************************************************
 * \brief
 *    Initializes the parameters for a new frame
 ************************************************************************
 */
void init_frame()
{
  int i,j,k;
  int prevP_no, nextP_no;
  if(input->InterlaceCodingOption >= MB_CODING)
      img->pstruct = 3;
  else
    img->pstruct = 0;    //frame coding

  last_P_no = last_P_no_frm;

  img->current_mb_nr=0;
  img->current_slice_nr=0;
  stat->bit_slice = 0;

  img->mb_y = img->mb_x = 0;
  img->block_y = img->pix_y = img->pix_c_y = 0;   // define vertical positions
  img->block_x = img->pix_x = img->block_c_x = img->pix_c_x = 0; // define horizontal positions

  // Tian Dong: June 7, 2002 JVT-B042
  // Initiate the actually number of ref frames that can be used.
  fb->num_short_used = fb->short_used;

  if(img->type != B_IMG)
  {
    img->refPicID_frm ++;
    img->refPicID = img->refPicID_frm;

//    img->tr=img->number*(input->jumpd+1);
    img->tr=start_tr_in_this_IGOP + IMG_NUMBER*(input->jumpd+1);

    img->imgtr_last_P_frm = img->imgtr_next_P_frm;
    img->imgtr_next_P_frm = img->tr;

    box_s.lastFrameNr = img->tr;  // serve as the duration of the segment.

#ifdef _ADAPT_LAST_GROUP_
    if (input->last_frame && img->number+1 == input->no_frames)
      img->tr=input->last_frame;
#endif

    if(IMG_NUMBER!=0 && input->successive_Bframe != 0)   // B pictures to encode
      nextP_tr_frm=img->tr;

    if (img->type == INTRA_IMG)
      img->qp = input->qp0;         // set quant. parameter for I-frame
    else
    {
#ifdef _CHANGE_QP_
      if (input->qp2start > 0 && img->tr >= input->qp2start)
        img->qp = input->qpN2;
      else
#endif
        img->qp = input->qpN;
      
      if (img->types==SP_IMG)
      {
        img->qp = input->qpsp;
        img->qpsp = input->qpsp_pred;
      }

    }

    img->mb_y_intra=img->mb_y_upd;   //  img->mb_y_intra indicates which GOB to intra code for this frame

    if (input->intra_upd > 0)          // if error robustness, find next GOB to update
    {
      img->mb_y_upd=(IMG_NUMBER/input->intra_upd) % (img->height/MB_BLOCK_SIZE);
    }
  }
  else
  {
    img->p_interval = input->jumpd+1;
//    prevP_no = (img->number-1)*img->p_interval;
//    nextP_no = img->number*img->p_interval;
    prevP_no = start_tr_in_this_IGOP + (IMG_NUMBER-1)*img->p_interval;
    nextP_no = start_tr_in_this_IGOP + (IMG_NUMBER)*img->p_interval;

#ifdef _ADAPT_LAST_GROUP_
    last_P_no[0] = prevP_no; for (i=1; i<img->buf_cycle; i++) last_P_no[i] = last_P_no[i-1]-img->p_interval;

    if (input->last_frame && img->number+1 == input->no_frames)
    {
      nextP_no        =input->last_frame;
      img->p_interval =nextP_no - prevP_no;
    }
#endif

    img->b_interval = (int)((float)(input->jumpd+1)/(input->successive_Bframe+1.0)+0.49999);

    img->tr= prevP_no+img->b_interval*img->b_frame_to_code; // from prev_P
    if(img->tr >= nextP_no)
      img->tr=nextP_no-1; 

#ifdef _CHANGE_QP_
    if (input->qp2start > 0 && img->tr >= input->qp2start)
      img->qp = input->qpB2;
    else
#endif
      img->qp = input->qpB;

    // initialize arrays
    for(k=0; k<2; k++)
      for(i=0; i<img->height/BLOCK_SIZE; i++)
        for(j=0; j<img->width/BLOCK_SIZE+4; j++)
        {
          tmp_fwMV[k][i][j]=0;
          tmp_bwMV[k][i][j]=0;
          dfMV[k][i][j]=0;
          dbMV[k][i][j]=0;
        }
    for(i=0; i<img->height/BLOCK_SIZE; i++)
      for(j=0; j<img->width/BLOCK_SIZE; j++)
      {
        fw_refFrArr[i][j]=bw_refFrArr[i][j]=-1;
      }

  }
}

/*!
 ************************************************************************
 * \brief
 *    Initializes the parameters for a new field
 ************************************************************************
 */
void init_field()
{
  int i,j,k;
  int prevP_no, nextP_no;

  last_P_no = last_P_no_fld;

  //picture structure
  if(!img->fld_type)
      img->pstruct = 1;
  else
    img->pstruct = 2;
  
  img->current_mb_nr=0;
  img->current_slice_nr=0;
  stat->bit_slice = 0;

  input->jumpd *= 2;
  input->successive_Bframe *= 2;
  img->number /= 2;
  img->buf_cycle /= 2;

#ifdef UMV
  img->mhor  = img->width *4-1;
  img->mvert = img->height*4-1;
#endif

  img->mb_y = img->mb_x = 0;
  img->block_y = img->pix_y = img->pix_c_y = 0;   // define vertical positions
  img->block_x = img->pix_x = img->block_c_x = img->pix_c_x = 0; // define horizontal positions

  if(img->type != B_IMG)
  {

    if(!img->fld_type)
      img->refPicID_fld++;   //increment by 1 for field 1 only.
    img->refPicID = img->refPicID_fld;

    img->tr=img->number*(input->jumpd+2)+img->fld_type;

    if(!img->fld_type) 
    {
      img->imgtr_last_P_fld = img->imgtr_next_P_fld;
      img->imgtr_next_P_fld = img->tr;
    }

#ifdef _ADAPT_LAST_GROUP_
    if (input->last_frame && img->number+1 == input->no_frames)
      img->tr=input->last_frame;
#endif
    if(img->number!=0 && input->successive_Bframe != 0)   // B pictures to encode
      nextP_tr_fld=img->tr;

    if (img->type == INTRA_IMG)
      img->qp = input->qp0;         // set quant. parameter for I-frame
    else
    {
#ifdef _CHANGE_QP_
      if (input->qp2start > 0 && img->tr >= input->qp2start)
        img->qp = input->qpN2;
      else
#endif
        img->qp = input->qpN;
      if (img->types==SP_IMG)
      {
        img->qp = input->qpsp;
        img->qpsp = input->qpsp_pred;
      }

    }

    img->mb_y_intra=img->mb_y_upd;   //  img->mb_y_intra indicates which GOB to intra code for this frame

    if (input->intra_upd > 0)          // if error robustness, find next GOB to update
    {
      img->mb_y_upd=(img->number/input->intra_upd) % (img->width/MB_BLOCK_SIZE);
    }
  }
  else
  {
    img->p_interval = input->jumpd+2;
    prevP_no = (img->number-1)*img->p_interval+img->fld_type;
    nextP_no = img->number*img->p_interval+img->fld_type;
#ifdef _ADAPT_LAST_GROUP_
    if (!img->fld_type)  // top field
    {
      last_P_no[0] = prevP_no + 1; 
      last_P_no[1] = prevP_no;
      for (i=1; i<=img->buf_cycle; i++) 
      {
        last_P_no[2*i] = last_P_no[2*i-2]-img->p_interval; 
        last_P_no[2*i+1] = last_P_no[2*i-1]-img->p_interval; 
      }
    }
    else    // bottom field
    {
      last_P_no[0] = nextP_no - 1; 
      last_P_no[1] = prevP_no;
      for (i=1; i<=img->buf_cycle; i++) 
      {
        last_P_no[2*i] = last_P_no[2*i-2]-img->p_interval; 
        last_P_no[2*i+1] = last_P_no[2*i-1]-img->p_interval; 
      }
    }

    if (input->last_frame && img->number+1 == input->no_frames)
      {
        nextP_no        =input->last_frame;
        img->p_interval =nextP_no - prevP_no;
      }
#endif

    img->b_interval = (int)((float)(input->jumpd+1)/(input->successive_Bframe+1.0)+0.49999);

    img->tr= prevP_no+(img->b_interval+1)*img->b_frame_to_code; // from prev_P
    if(img->tr >= nextP_no)
      img->tr=nextP_no-1; // ?????

#ifdef _CHANGE_QP_
    if (input->qp2start > 0 && img->tr >= input->qp2start)
      img->qp = input->qpB2;
    else
#endif
      img->qp = input->qpB;

    // initialize arrays
    for(k=0; k<2; k++)
      for(i=0; i<img->height/BLOCK_SIZE; i++)
        for(j=0; j<img->width/BLOCK_SIZE+4; j++)
        {
          tmp_fwMV[k][i][j]=0;
          tmp_bwMV[k][i][j]=0;
          dfMV[k][i][j]=0;
          dbMV[k][i][j]=0;
        }
    for(i=0; i<img->height/BLOCK_SIZE; i++)
      for(j=0; j<img->width/BLOCK_SIZE; j++)
      {
        fw_refFrArr[i][j]=bw_refFrArr[i][j]=-1;
      }

  }
  input->jumpd /= 2;
  input->successive_Bframe /= 2;
  img->buf_cycle *= 2;
  img->number = 2 * img->number + img->fld_type;
  if (img->type == BS_IMG)
  {
    img->num_ref_pic_active_fwd_minus1 = max(0,min(img->number-1,img->buf_cycle+img->fld_type-1));
    img->num_ref_pic_active_bwd_minus1 = max(0,min(img->number-1,3+img->fld_type));
  }
  else if (img->type == B_IMG)
  {
    img->num_ref_pic_active_fwd_minus1 = max(0,min(img->number-1,img->buf_cycle-1));
    img->num_ref_pic_active_bwd_minus1 = 0;
  }
  else
  {
    img->num_ref_pic_active_fwd_minus1 = max(0,min(img->number-1,img->buf_cycle+img->fld_type-1));
    img->num_ref_pic_active_bwd_minus1 = 0;
  }
}

/*!
 ************************************************************************
 * \brief
 *    Initializes the parameters for a new slice
 ************************************************************************
 */
void init_slice()
{

  int i;
  Slice *curr_slice = img->currentSlice;
  DataPartition *dataPart;
  Bitstream *currStream;

  curr_slice->picture_id = img->tr%256;
  curr_slice->qp = img->qp;
  curr_slice->start_mb_nr = img->current_mb_nr;
  curr_slice->slice_too_big = dummy_slice_too_big;

  for (i=0; i<curr_slice->max_part_nr; i++)
  {
    dataPart = &(curr_slice->partArr[i]);

    // in priciple it is possible to assign to each partition
    // a different entropy coding method
    if (input->symbol_mode == UVLC)
      dataPart->writeSyntaxElement = writeSyntaxElement_UVLC;
    else
      dataPart->writeSyntaxElement = writeSyntaxElement_CABAC;

    // A little hack until CABAC can handle non-byte aligned start positions   StW!
    // For UVLC, the stored_ positions in the bit buffer are necessary.  For CABAC,
    // the buffer is initialized to start at zero.

    if (input->symbol_mode == UVLC && (input->of_mode == PAR_OF_26L /*|| input->of_mode == PAR_OF_IFF */ ))    // Stw: added PAR_OF_26L check
    {
      currStream = dataPart->bitstream;
      currStream->bits_to_go  = currStream->stored_bits_to_go;
      currStream->byte_pos    = currStream->stored_byte_pos;
      currStream->byte_buf    = currStream->stored_byte_buf;
      currStream->tmp_byte_pos = currStream->stored_byte_pos; // temp store curr byte position
    } 
    else // anything but UVLC and PAR_OF_26L
    {   
      if (input->of_mode == PAR_OF_26L)
      {
        if((img->current_slice_nr==0)&&(img->pstruct!=2))
        {
          currStream = dataPart->bitstream;
          currStream->bits_to_go  = 8;
          currStream->byte_pos    = 0;
          currStream->byte_buf    = 0;
          currStream->tmp_byte_pos= 0;
        }
      }
      else
      {
        currStream = dataPart->bitstream;
        currStream->bits_to_go  = 8;
        currStream->byte_pos    = 0;
        currStream->byte_buf    = 0;
      }
    }
  }
}

/*!
 ************************************************************************
 * \brief
 *    Reads new frame from file and sets frame_no
 ************************************************************************
 */
void read_one_new_frame()
{
  int i, j, k, uv;
  int status; // frame_no;
  int frame_size = img->height*img->width*3/2;

  if(img->type == B_IMG)
//    frame_no = (img->number-1)*(input->jumpd+1)+img->b_interval*img->b_frame_to_code;
    frame_no = start_tr_in_this_IGOP + (IMG_NUMBER-1)*(input->jumpd+1)+img->b_interval*img->b_frame_to_code;
  else
  {
    frame_no = start_tr_in_this_IGOP + IMG_NUMBER*(input->jumpd+1);
    if ( input->of_mode == PAR_OF_IFF )
      box_ati.info[0].last_frame = frame_no;
#ifdef _ADAPT_LAST_GROUP_
      if (input->last_frame && img->number+1 == input->no_frames)
        frame_no=input->last_frame;
#endif
  }

  rewind (p_in);

  status  = fseek (p_in, frame_no * frame_size + input->infile_header, 0);

  if (status != 0)
  {
    snprintf(errortext, ET_SIZE, "Error in seeking frame no: %d\n", frame_no);
    error(errortext,1);
  }

  for (j=0; j < img->height; j++)
    for (i=0; i < img->width; i++)
      imgY_org[j][i]=fgetc(p_in);
  for (uv=0; uv < 2; uv++)
    for (j=0; j < img->height_cr ; j++)
      for (i=0; i < img->width_cr; i++)
        imgUV_org[uv][j][i]=fgetc(p_in);

  if(input->InterlaceCodingOption >= MB_CODING)
  {
    for (j=0; j < img->height/2; j++)
      for (k=2*j, i=0; i < img->width; i++)
        imgY_org_top[j][i]=imgY_org[k][i];
   for (uv=0; uv < 2; uv++)
     for (j=0; j < img->height_cr/2 ; j++)
       for (k=2*j, i=0; i < img->width_cr; i++)
         imgUV_org_top[uv][j][i]=imgUV_org[uv][k][i];

    for (j=0; j < img->height/2; j++)
      for (k=2*j+1, i=0; i < img->width; i++)
        imgY_org_bot[j][i]=imgY_org[k][i];
    for (uv=0; uv < 2; uv++)
      for (j=0; j < img->height_cr/2 ; j++)
        for (k=2*j+1, i=0; i < img->width_cr; i++)
          imgUV_org_bot[uv][j][i]=imgUV_org[uv][k][i];
  }
}

/*!
 ************************************************************************
 * \brief
 *    Reads new field from file and sets frame_no
 ************************************************************************
 */
void read_one_new_field()
{
  int i, j, uv;
  int status; // frame_no;
  int frame_size = img->height*img->width*3;  // (3/2) * 2
  int picture_no;

  picture_no = img->number / 2;
  if(img->type == B_IMG)
    frame_no = (picture_no-1)*(input->jumpd+1)+img->b_interval*img->b_frame_to_code;
  else
  {
    frame_no = picture_no*(input->jumpd+1);
#ifdef _ADAPT_LAST_GROUP_
      if (input->last_frame && picture_no+1 == input->no_frames)
        frame_no=input->last_frame;
#endif
  }

  rewind (p_in);
  status = fseek (p_in, frame_no * frame_size + input->infile_header, 0);

  if (status != 0)
  {
    snprintf(errortext, ET_SIZE, "Error in seeking frame no: %d\n", frame_no);
    error(errortext,1);
  }

  // img->height and img->height_cr must be adjusted for field before this function
  if (img->fld_type==0)  // top field
  {
    for (j=0; j < img->height; j++)
    {
      for (i=0; i < img->width; i++)
        imgY_org[j][i]=fgetc(p_in);
      for (i=0; i < img->width; i++)
        fgetc(p_in);
    }
    for (uv=0; uv < 2; uv++)
      for (j=0; j < img->height_cr; j++)
      {
        for (i=0; i < img->width_cr; i++)
          imgUV_org[uv][j][i]=fgetc(p_in);
        for (i=0; i < img->width_cr; i++)
          fgetc(p_in);
      }
  }
  else           // bottom field
  {
    for (j=0; j < img->height; j++)
    {
      for (i=0; i < img->width; i++)
        fgetc(p_in);
      for (i=0; i < img->width; i++)
        imgY_org[j][i]=fgetc(p_in);
    }
    for (uv=0; uv < 2; uv++)
      for (j=0; j < img->height_cr; j++)
      {
        for (i=0; i < img->width_cr; i++)
          fgetc(p_in);
        for (i=0; i < img->width_cr; i++)
          imgUV_org[uv][j][i]=fgetc(p_in);
      }
  }
}

/*!
 ************************************************************************
 * \brief
 *     Writes reconstructed image(s) to file
 *     This can be done more elegant!
 ************************************************************************
 */
void write_reconstructed_image()
{
  int i, j, k;
  int start=0, inc=1;

  if (p_dec != NULL)
  {
    if(img->type != B_IMG)
    {
      // write reconstructed image (IPPP)
      if(input->successive_Bframe==0)
      {
        for (i=start; i < img->height; i+=inc)
          for (j=0; j < img->width; j++)
            fputc(min(imgY[i][j],255),p_dec);

        for (k=0; k < 2; ++k)
          for (i=start; i < img->height/2; i+=inc)
            for (j=0; j < img->width/2; j++)
              fputc(min(imgUV[k][i][j],255),p_dec);
      }

      // write reconstructed image (IBPBP) : only intra written
      else if (IMG_NUMBER==0 && input->successive_Bframe!=0)
      {
        for (i=start; i < img->height; i+=inc)
          for (j=0; j < img->width; j++)
            fputc(min(imgY[i][j],255),p_dec);

        for (k=0; k < 2; ++k)
          for (i=start; i < img->height/2; i+=inc)
            for (j=0; j < img->width/2; j++)
              fputc(min(imgUV[k][i][j],255),p_dec);
      }

      // next P picture. This is saved with recon B picture after B picture coding
      if (IMG_NUMBER!=0 && input->successive_Bframe!=0)
      {
        for (i=start; i < img->height; i+=inc)
          for (j=0; j < img->width; j++)
            nextP_imgY[i][j]=imgY[i][j];
        for (k=0; k < 2; ++k)
          for (i=start; i < img->height/2; i+=inc)
            for (j=0; j < img->width/2; j++)
              nextP_imgUV[k][i][j]=imgUV[k][i][j];
      }
    }
    else
    {
      for (i=start; i < img->height; i+=inc)
        for (j=0; j < img->width; j++)
          fputc(min(imgY[i][j],255),p_dec);
      for (k=0; k < 2; ++k)
        for (i=start; i < img->height/2; i+=inc)
          for (j=0; j < img->width/2; j++)
            fputc(min(imgUV[k][i][j],255),p_dec);

      // If this is last B frame also store P frame
      if(img->b_frame_to_code == input->successive_Bframe)
      {
        // save P picture
        for (i=start; i < img->height; i+=inc)
          for (j=0; j < img->width; j++)
            fputc(min(nextP_imgY[i][j],255),p_dec);
        for (k=0; k < 2; ++k)
          for (i=start; i < img->height/2; i+=inc)
            for (j=0; j < img->width/2; j++)
              fputc(min(nextP_imgUV[k][i][j],255),p_dec);
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
void interpolate_frame_to_fb()
{                 // write to mref[]
  add_frame(img);
  init_mref(img);
  init_Refbuf(img);

  if(input->mv_res)
    oneeighthpix();
  else
    UnifiedOneForthPix(imgY, imgUV[0], imgUV[1],
               mref[0], mcef[0][0], mcef[0][1],
               Refbuf11[0]);
}

/*!
 ************************************************************************
 * \brief
 *    Choose interpolation method depending on MV-resolution
 ************************************************************************
 */
void interpolate_frame()
{                 // write to mref[]
  init_mref(img);
  init_Refbuf(img);

  if(input->mv_res)
    oneeighthpix();
  else
    UnifiedOneForthPix(imgY, imgUV[0], imgUV[1],
               mref[0], mcef[0][0], mcef[0][1],
               Refbuf11[0]);
}

static void GenerateFullPelRepresentation (pel_t **Fourthpel, pel_t *Fullpel, int xsize, int ysize)
{

  int x, y;

  for (y=0; y<ysize; y++)
    for (x=0; x<xsize; x++)
      PutPel_11 (Fullpel, y, x, FastPelY_14 (Fourthpel, y*4, x*4));
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
void UnifiedOneForthPix (pel_t **imgY, pel_t** imgU, pel_t **imgV,
                         pel_t **out4Y, pel_t **outU, pel_t **outV,
                         pel_t *ref11)
{
  int is;
  int i,j,j4;
  int ie2,je2,jj,maxy;


  for (j=-IMG_PAD_SIZE; j < img->height+IMG_PAD_SIZE; j++)
  {
    for (i=-IMG_PAD_SIZE; i < img->width+IMG_PAD_SIZE; i++)
    {
      jj = max(0,min(img->height-1,j));
      is=(ONE_FOURTH_TAP[0][0]*(imgY[jj][max(0,min(img->width-1,i  ))]+imgY[jj][max(0,min(img->width-1,i+1))])+
          ONE_FOURTH_TAP[1][0]*(imgY[jj][max(0,min(img->width-1,i-1))]+imgY[jj][max(0,min(img->width-1,i+2))])+
          ONE_FOURTH_TAP[2][0]*(imgY[jj][max(0,min(img->width-1,i-2))]+imgY[jj][max(0,min(img->width-1,i+3))]));
      img4Y_tmp[j+IMG_PAD_SIZE][(i+IMG_PAD_SIZE)*2  ]=imgY[jj][max(0,min(img->width-1,i))]*1024;    // 1/1 pix pos
      img4Y_tmp[j+IMG_PAD_SIZE][(i+IMG_PAD_SIZE)*2+1]=is*32;              // 1/2 pix pos
    }
  }
  
  for (i=0; i < (img->width+2*IMG_PAD_SIZE)*2; i++)
  {
    for (j=0; j < img->height+2*IMG_PAD_SIZE; j++)
    {
      j4=j*4;
      maxy = img->height+2*IMG_PAD_SIZE-1;
      // change for TML4, use 6 TAP vertical filter
      is=(ONE_FOURTH_TAP[0][0]*(img4Y_tmp[j         ][i]+img4Y_tmp[min(maxy,j+1)][i])+
          ONE_FOURTH_TAP[1][0]*(img4Y_tmp[max(0,j-1)][i]+img4Y_tmp[min(maxy,j+2)][i])+
          ONE_FOURTH_TAP[2][0]*(img4Y_tmp[max(0,j-2)][i]+img4Y_tmp[min(maxy,j+3)][i]))/32;

      PutPel_14 (out4Y, (j-IMG_PAD_SIZE)*4,   (i-IMG_PAD_SIZE*2)*2,(pel_t) max(0,min(255,(int)((img4Y_tmp[j][i]+512)/1024))));     // 1/2 pix
      PutPel_14 (out4Y, (j-IMG_PAD_SIZE)*4+2, (i-IMG_PAD_SIZE*2)*2,(pel_t) max(0,min(255,(int)((is+512)/1024)))); // 1/2 pix
    }
  }

  /* 1/4 pix */
  /* luma */
  ie2=(img->width+2*IMG_PAD_SIZE-1)*4;
  je2=(img->height+2*IMG_PAD_SIZE-1)*4;

  for (j=0;j<je2+4;j+=2)
    for (i=0;i<ie2+3;i+=2) {
      /*  '-'  */
      PutPel_14 (out4Y, j-IMG_PAD_SIZE*4, i-IMG_PAD_SIZE*4+1, (pel_t) (max(0,min(255,(int)(FastPelY_14(out4Y, j-IMG_PAD_SIZE*4, i-IMG_PAD_SIZE*4)+FastPelY_14(out4Y, j-IMG_PAD_SIZE*4, min(ie2+2,i+2)-IMG_PAD_SIZE*4))/2))));
    }
  for (i=0;i<ie2+4;i++)
  {
    for (j=0;j<je2+3;j+=2)
    {
      if( i%2 == 0 ) {
        /*  '|'  */
        PutPel_14 (out4Y, j-IMG_PAD_SIZE*4+1, i-IMG_PAD_SIZE*4, (pel_t)(max(0,min(255,(int)(FastPelY_14(out4Y, j-IMG_PAD_SIZE*4, i-IMG_PAD_SIZE*4)+FastPelY_14(out4Y, min(je2+2,j+2)-IMG_PAD_SIZE*4, i-IMG_PAD_SIZE*4))/2))));
      }
      else if( ((i&3) == 3)&&(((j+1)&3) == 3))
      {
        /* "funny posision" */
        PutPel_14 (out4Y, j-IMG_PAD_SIZE*4+1, i-IMG_PAD_SIZE*4, (pel_t) ((
          FastPelY_14 (out4Y, j-IMG_PAD_SIZE*4-2, i-IMG_PAD_SIZE*4-3) +
          FastPelY_14 (out4Y, min(je2,j+2)-IMG_PAD_SIZE*4, i-IMG_PAD_SIZE*4-3) +
          FastPelY_14 (out4Y, j-IMG_PAD_SIZE*4-2, min(ie2,i+1)-IMG_PAD_SIZE*4) +
          FastPelY_14 (out4Y, min(je2,j+2)-IMG_PAD_SIZE*4, min(ie2,i+1)-IMG_PAD_SIZE*4)
          + 2 )/4));
      }
      else if ((j%4 == 0 && i%4 == 1) || (j%4 == 2 && i%4 == 3)) {
        /*  '/'  */
        PutPel_14 (out4Y, j-IMG_PAD_SIZE*4+1, i-IMG_PAD_SIZE*4, (pel_t)(max(0,min(255,(int)(
                   FastPelY_14 (out4Y, j-IMG_PAD_SIZE*4, min(ie2+2,i+1)-IMG_PAD_SIZE*4) +
                   FastPelY_14(out4Y, min(je2+2,j+2)-IMG_PAD_SIZE*4, i-IMG_PAD_SIZE*4-1))/2))));
      }
      else {
        /*  '\'  */
        PutPel_14 (out4Y, j-IMG_PAD_SIZE*4+1, i-IMG_PAD_SIZE*4, (pel_t)(max(0,min(255,(int)(
          FastPelY_14 (out4Y, j-IMG_PAD_SIZE*4, i-IMG_PAD_SIZE*4-1) +
          FastPelY_14(out4Y, min(je2+2,j+2)-IMG_PAD_SIZE*4, min(ie2+2,i+1)-IMG_PAD_SIZE*4))/2))));
      }
    }
  }

  /*  Chroma: */
  for (j=0; j < img->height_cr; j++) {
    memcpy(outU[j],imgUV[0][j],img->width_cr); // just copy 1/1 pix, interpolate "online" 
    memcpy(outV[j],imgUV[1][j],img->width_cr);
  }

  // Generate 1/1th pel representation (used for integer pel MV search)
  GenerateFullPelRepresentation (out4Y, ref11, img->width, img->height);

}

/*!
 ************************************************************************
 * \brief
 *    Upsample 4 times for 1/8-pel estimation and store in buffer
 *    for multiple reference frames. 1/8-pel resolution is calculated
 *    during the motion estimation on the fly with bilinear interpolation.
 *
 ************************************************************************
 */
void oneeighthpix()
{
  static int h1[8] = {  -3, 12, -37, 229,  71, -21,  6, -1 };  
  static int h2[8] = {  -3, 12, -39, 158, 158, -39, 12, -3 };  
  static int h3[8] = {  -1,  6, -21,  71, 229, -37, 12, -3 };  

  int uv,x,y,y1,x4,y4,x4p;

  int nx_out, ny_out, nx_1, ny_1, maxy;
  int i0,i1,i2,i3;

  nx_out = 4*img->width;
  ny_out = 4*img->height;
  nx_1   = img->width-1;
  ny_1   = img->height-1;


  //horizontal filtering filtering
  for(y=-IMG_PAD_SIZE;y<img->height+IMG_PAD_SIZE;y++)
  {
    for(x=-IMG_PAD_SIZE;x<img->width+IMG_PAD_SIZE;x++)
    {
      y1 = max(0,min(ny_1,y));

      i0=(256*imgY[y1][max(0,min(nx_1,x))]);
      
      i1=(
        h1[0]* imgY[y1][max(0,min(nx_1,x-3))]  +
        h1[1]* imgY[y1][max(0,min(nx_1,x-2))]  +
        h1[2]* imgY[y1][max(0,min(nx_1,x-1))]  +
        h1[3]* imgY[y1][max(0,min(nx_1,x  ))]  +
        h1[4]* imgY[y1][max(0,min(nx_1,x+1))]  +
        h1[5]* imgY[y1][max(0,min(nx_1,x+2))]  +
        h1[6]* imgY[y1][max(0,min(nx_1,x+3))]  +                         
        h1[7]* imgY[y1][max(0,min(nx_1,x+4))] );
      
      
      i2=(
        h2[0]* imgY[y1][max(0,min(nx_1,x-3))]  +
        h2[1]* imgY[y1][max(0,min(nx_1,x-2))]  +
        h2[2]* imgY[y1][max(0,min(nx_1,x-1))]  +
        h2[3]* imgY[y1][max(0,min(nx_1,x  ))]  +
        h2[4]* imgY[y1][max(0,min(nx_1,x+1))]  +
        h2[5]* imgY[y1][max(0,min(nx_1,x+2))]  +
        h2[6]* imgY[y1][max(0,min(nx_1,x+3))]  +                         
        h2[7]* imgY[y1][max(0,min(nx_1,x+4))] );
      
      
      i3=(
        h3[0]* imgY[y1][max(0,min(nx_1,x-3))]  +
        h3[1]* imgY[y1][max(0,min(nx_1,x-2))]  +
        h3[2]* imgY[y1][max(0,min(nx_1,x-1))]  +
        h3[3]* imgY[y1][max(0,min(nx_1,x  ))]  +
        h3[4]* imgY[y1][max(0,min(nx_1,x+1))]  +
        h3[5]* imgY[y1][max(0,min(nx_1,x+2))]  +
        h3[6]* imgY[y1][max(0,min(nx_1,x+3))]  +                         
        h3[7]* imgY[y1][max(0,min(nx_1,x+4))] );
      
      x4=(x+IMG_PAD_SIZE)*4;

      img4Y_tmp[y+IMG_PAD_SIZE][x4  ] = i0;
      img4Y_tmp[y+IMG_PAD_SIZE][x4+1] = i1;
      img4Y_tmp[y+IMG_PAD_SIZE][x4+2] = i2;
      img4Y_tmp[y+IMG_PAD_SIZE][x4+3] = i3;
    }
  }

  maxy = img->height+2*IMG_PAD_SIZE-1;

  for(x4=0;x4<nx_out+2*IMG_PAD_SIZE*4;x4++)
  {
    for(y=0;y<=maxy;y++)
    {
      i0=(long int)(img4Y_tmp[y][x4]+256/2)/256;
      
      i1=(long int)( 
        h1[0]* img4Y_tmp[max(0   ,y-3)][x4]+
        h1[1]* img4Y_tmp[max(0   ,y-2)][x4]+
        h1[2]* img4Y_tmp[max(0   ,y-1)][x4]+
        h1[3]* img4Y_tmp[y][x4]            +
        h1[4]* img4Y_tmp[min(maxy,y+1)][x4]+
        h1[5]* img4Y_tmp[min(maxy,y+2)][x4]+
        h1[6]* img4Y_tmp[min(maxy,y+3)][x4]+ 
        h1[7]* img4Y_tmp[min(maxy,y+4)][x4]+ 256*256/2 ) / (256*256);
      
      i2=(long int)( 
        h2[0]* img4Y_tmp[max(0   ,y-3)][x4]+
        h2[1]* img4Y_tmp[max(0   ,y-2)][x4]+
        h2[2]* img4Y_tmp[max(0   ,y-1)][x4]+
        h2[3]* img4Y_tmp[y][x4]            +
        h2[4]* img4Y_tmp[min(maxy,y+1)][x4]+
        h2[5]* img4Y_tmp[min(maxy,y+2)][x4]+
        h2[6]* img4Y_tmp[min(maxy,y+3)][x4]+ 
        h2[7]* img4Y_tmp[min(maxy,y+4)][x4]+ 256*256/2 ) / (256*256);
      
      i3=(long int)( 
        h3[0]* img4Y_tmp[max(0   ,y-3)][x4]+
        h3[1]* img4Y_tmp[max(0   ,y-2)][x4]+
        h3[2]* img4Y_tmp[max(0   ,y-1)][x4]+
        h3[3]* img4Y_tmp[y][x4]            +
        h3[4]* img4Y_tmp[min(maxy,y+1)][x4]+
        h3[5]* img4Y_tmp[min(maxy,y+2)][x4]+
        h3[6]* img4Y_tmp[min(maxy,y+3)][x4]+ 
        h3[7]* img4Y_tmp[min(maxy,y+4)][x4]+ 256*256/2 ) / (256*256);
      
      y4  = (y-IMG_PAD_SIZE)*4;
      x4p = x4-IMG_PAD_SIZE*4;
  
      PutPel_14 (mref[0], y4,   x4p, (pel_t) max(0,min(255,i0)));   
      PutPel_14 (mref[0], y4+1, x4p, (pel_t) max(0,min(255,i1)));   
      PutPel_14 (mref[0], y4+2, x4p, (pel_t) max(0,min(255,i2)));
      PutPel_14 (mref[0], y4+3, x4p, (pel_t) max(0,min(255,i3)));   

    }
  }

  for(y=0;y<img->height;y++)
    for(x=0;x<img->width;x++)
      PutPel_11 (Refbuf11[0], y, x, FastPelY_14 (mref[0], y*4, x*4));

  for (uv=0; uv < 2; uv++)
    for (y=0; y < img->height_cr; y++)
      memcpy(mcef[0][uv][y],imgUV[uv][y],img->width_cr); // just copy 1/1 pix, interpolate "online"
  GenerateFullPelRepresentation (mref[0], Refbuf11[0], img->width, img->height);
  // Generate 1/1th pel representation (used for integer pel MV search)

}

/*!
 ************************************************************************
 * \brief
 *    Find SNR for all three components
 ************************************************************************
 */
void find_snr()
{
  int i,j;
  int diff_y,diff_u,diff_v;
  int impix;

  //  Calculate  PSNR for Y, U and V.

  //     Luma.
  impix = img->height*img->width;

  diff_y=0;
  for (i=0; i < img->width; ++i)
  {
    for (j=0; j < img->height; ++j)
    {
      diff_y += img->quad[imgY_org[j][i]-imgY[j][i]];
    }
  }

  //     Chroma.

  diff_u=0;
  diff_v=0;

  for (i=0; i < img->width_cr; i++)
  {
    for (j=0; j < img->height_cr; j++)
    {
      diff_u += img->quad[imgUV_org[0][j][i]-imgUV[0][j][i]];
      diff_v += img->quad[imgUV_org[1][j][i]-imgUV[1][j][i]];
    }
  }

  //  Collecting SNR statistics
  if (diff_y != 0)
  {
    snr->snr_y=(float)(10*log10(65025*(float)impix/(float)diff_y));        // luma snr for current frame
    snr->snr_u=(float)(10*log10(65025*(float)impix/(float)(4*diff_u)));    // u croma snr for current frame, 1/4 of luma samples
    snr->snr_v=(float)(10*log10(65025*(float)impix/(float)(4*diff_v)));    // v croma snr for current frame, 1/4 of luma samples
  }

  if (img->number == 0)
  {
    snr->snr_y1=(float)(10*log10(65025*(float)impix/(float)diff_y));       // keep luma snr for first frame
    snr->snr_u1=(float)(10*log10(65025*(float)impix/(float)(4*diff_u)));   // keep croma u snr for first frame
    snr->snr_v1=(float)(10*log10(65025*(float)impix/(float)(4*diff_v)));   // keep croma v snr for first frame
    snr->snr_ya=snr->snr_y1;
    snr->snr_ua=snr->snr_u1;
    snr->snr_va=snr->snr_v1;
  }
  // B pictures
  else
  {
    snr->snr_ya=(float)(snr->snr_ya*(img->number+Bframe_ctr)+snr->snr_y)/(img->number+Bframe_ctr+1);   // average snr lume for all frames inc. first
    snr->snr_ua=(float)(snr->snr_ua*(img->number+Bframe_ctr)+snr->snr_u)/(img->number+Bframe_ctr+1);   // average snr u croma for all frames inc. first
    snr->snr_va=(float)(snr->snr_va*(img->number+Bframe_ctr)+snr->snr_v)/(img->number+Bframe_ctr+1);   // average snr v croma for all frames inc. first
  }
}

/*!
 ************************************************************************
 * \brief
 *    Find distortion for all three components
 ************************************************************************
 */
void find_distortion()
{
  int i,j;
  int diff_y,diff_u,diff_v;
  int impix;

  //  Calculate  PSNR for Y, U and V.

  //     Luma.
  impix = img->height*img->width;

  diff_y=0;
  for (i=0; i < img->width; ++i)
  {
    for (j=0; j < img->height; ++j)
    {
      diff_y += img->quad[abs(imgY_org[j][i]-imgY[j][i])];
    }
  }

  //     Chroma.

  diff_u=0;
  diff_v=0;

  for (i=0; i < img->width_cr; i++)
  {
    for (j=0; j < img->height_cr; j++)
    {
      diff_u += img->quad[abs(imgUV_org[0][j][i]-imgUV[0][j][i])];
      diff_v += img->quad[abs(imgUV_org[1][j][i]-imgUV[1][j][i])];
    }
  }

  // Calculate real PSNR at find_snr_avg()
  snr->snr_y = (float)diff_y;
  snr->snr_u = (float)diff_u;
  snr->snr_v = (float)diff_v;
}


void rotate_buffer()
{
  Frame *f;

  f=fb->picbuf_short[1];
  fb->picbuf_short[1]=fb->picbuf_short[0];
  fb->picbuf_short[0]=f;

  mref[0]=fb->picbuf_short[0]->mref;
  mcef[0]=fb->picbuf_short[0]->mcef;
  mref[1]=fb->picbuf_short[1]->mref;
  mcef[1]=fb->picbuf_short[1]->mcef;
}

void store_field_MV(int frame_number)
{
  int i, j;  

  if (img->type != B_IMG) //all I- and P-frames
  {
    if (img->fld_flag)
    {
      for (i=0 ; i<img->width/4+4 ; i++)
      {
        for (j=0 ; j<img->height/8 ; j++)
        {
          tmp_mv_frm[0][2*j][i] = tmp_mv_frm[0][2*j+1][i] = tmp_mv_top[0][j][i];
          tmp_mv_frm[0][2*j][i] = tmp_mv_frm[0][2*j+1][i] = tmp_mv_top[0][j][i];    // ??
          tmp_mv_frm[1][2*j][i] = tmp_mv_frm[1][2*j+1][i] = tmp_mv_top[1][j][i]*2;
          tmp_mv_frm[1][2*j][i] = tmp_mv_frm[1][2*j+1][i] = tmp_mv_top[1][j][i]*2;  // ??
 
          if (input->direct_type && (input->successive_Bframe!=0 || input->BipredictiveWeighting > 0) && (i<img->width/4))
          {          

            moving_block_frm[2*j+1][i]=moving_block_frm[2*j][i]=
              ((refFrArr_top[j][i]!=0) || (refFrArr_bot[j][i]!=0) 
              || (abs(tmp_mv_top[0][j][i + 4])>>1) || (abs(tmp_mv_top[1][j][i+ 4])>>1) 
              || (abs(tmp_mv_bot[0][j][i + 4])>>1) || (abs(tmp_mv_bot[1][j][i+ 4])>>1));

            moving_block_top[j][i]=((refFrArr_top[j][i]!=0) 
              || (abs(tmp_mv_top[0][j][i+ 4])>>1) || (abs(tmp_mv_top[1][j][i+ 4])>>1));
            
            
            moving_block_bot[j][i]=((refFrArr_bot[j][i]!=0) 
              || (abs(tmp_mv_bot[0][j][i+ 4])>>1) || (abs(tmp_mv_bot[1][j][i+ 4])>>1));            
          }

          if ((i%2 == 0) && (j%2 == 0) && (i<img->width/4))
          {
            if (refFrArr_top[j][i] == -1)
            {
              refFrArr_frm[2*j][i] = refFrArr_frm[2*j+1][i] = -1;
              refFrArr_frm[2*(j+1)][i] = refFrArr_frm[2*(j+1)+1][i] = -1;
              refFrArr_frm[2*j][i+1] = refFrArr_frm[2*j+1][i+1] = -1;
              refFrArr_frm[2*(j+1)][i+1] = refFrArr_frm[2*(j+1)+1][i+1] = -1;
            }
            else
            {
              refFrArr_frm[2*j][i] = refFrArr_frm[2*j+1][i] = (int)(refFrArr_top[j][i]/2);
              refFrArr_frm[2*(j+1)][i] = refFrArr_frm[2*(j+1)+1][i] = (int)(refFrArr_top[j][i]/2);
              refFrArr_frm[2*j][i+1] = refFrArr_frm[2*j+1][i+1] = (int)(refFrArr_top[j][i]/2);
              refFrArr_frm[2*(j+1)][i+1] = refFrArr_frm[2*(j+1)+1][i+1] = (int)(refFrArr_top[j][i]/2);
            }
          }
        }
      }
    }
    else
    {
      for (i=0 ; i<img->width/4+4 ; i++)
      {
        for (j=0 ; j<img->height/8 ; j++)
        {
          tmp_mv_top[0][j][i] = tmp_mv_bot[0][j][i] = (int)(tmp_mv_frm[0][2*j][i]);
          tmp_mv_top[1][j][i] = tmp_mv_bot[1][j][i] = (int)((tmp_mv_frm[1][2*j][i])/2);
          if (input->direct_type && (input->successive_Bframe!=0 || input->BipredictiveWeighting > 0) && i<img->width/4)
          {
            moving_block_top[j][i]=moving_block_bot[j][i]=
              ((refFrArr_frm[2*j][i]!=0) || (refFrArr_frm[2*j + 1][i]!=0) 
              || (abs(tmp_mv_frm[0][2*j][i + 4])>>1) || (abs(tmp_mv_frm[1][2*j][i + 4])>>1) 
              || (abs(tmp_mv_frm[0][2*j+1][i + 4])>>1) || (abs(tmp_mv_frm[1][2*j+1][i + 4])>>1));
            
            moving_block_frm[2*j][i]=((refFrArr_frm[2*j][i]!=0) 
              || (abs(tmp_mv_frm[0][2*j][i+ 4])>>1) || (abs(tmp_mv_frm[1][2*j][i+ 4])>>1) );

            moving_block_frm[2*j+1][i]=((refFrArr_frm[2*j+1][i]!=0) 
              || (abs(tmp_mv_frm[0][2*j+1][i+ 4])>>1) || (abs(tmp_mv_frm[1][2*j+1][i+ 4])>>1));
          }
          if ((i%2 == 0) && (j%2 == 0) && (i<img->width/4))
          {
            if (refFrArr_frm[2*j][i] == -1)
            {
              refFrArr_top[j][i] = refFrArr_bot[j][i] = -1;
              refFrArr_top[j+1][i] = refFrArr_bot[j+1][i] = -1;
              refFrArr_top[j][i+1] = refFrArr_bot[j][i+1] = -1;
              refFrArr_top[j+1][i+1] = refFrArr_bot[j+1][i+1] = -1;
            }
            else
            {
              refFrArr_top[j][i] = refFrArr_bot[j][i] = refFrArr_frm[2*j][i]*2;
              refFrArr_top[j+1][i] = refFrArr_bot[j+1][i] = refFrArr_frm[2*j][i]*2;
              refFrArr_top[j][i+1] = refFrArr_bot[j][i+1] = refFrArr_frm[2*j][i]*2;
              refFrArr_top[j+1][i+1] = refFrArr_bot[j+1][i+1] = refFrArr_frm[2*j][i]*2;
            }
          }
        }
      }
    }
  }
}

void store_direct_moving_flag(int frame_number)
{
  int i, j;  

  if (img->type != B_IMG) //all I- and P-frames
    for (i=0 ; i<img->width/4 ; i++)
      for (j=0 ; j<img->height/8 ; j++)
        if (input->direct_type && (input->successive_Bframe!=0 || input->BipredictiveWeighting > 0))
        {            
          moving_block_frm[2*j][i]=((refFrArr_frm[2*j][i]!=0) 
            || (abs(tmp_mv_frm[0][2*j][i+ 4])>>1) || (abs(tmp_mv_frm[1][2*j][i+ 4])>>1) );          
          moving_block_frm[2*j+1][i]=((refFrArr_frm[2*j+1][i]!=0) 
            || (abs(tmp_mv_frm[0][2*j+1][i+ 4])>>1) || (abs(tmp_mv_frm[1][2*j+1][i+ 4])>>1));
        }      
}


// ABT
int  field2frame_mode(int fld_mode)
{
  int frm_mode;
  static const int field2frame_map[MAXMODE] = {-1, 1,1,3,3,4,6,6, -1, 7,1, -1, -1};

  frm_mode = field2frame_map[fld_mode];
  assert(frm_mode>0);

  return frm_mode;

}
int frame2field_mode(int frm_mode)
{
  int fld_mode;
  static const int frame2field_map[MAXMODE] = {-1, 2,5,4,5,5,7,7, -1, 7,2, -1, -1};

  fld_mode = frame2field_map[frm_mode];
  assert(fld_mode>0);

  return fld_mode;

}
void store_field_colB8mode()
{
  int i, j;

  if (input->abt)
  {
    if (img->type != B_IMG && img->type != BS_IMG) //all I- and P-frames
    {
      if (img->fld_flag)
      {
        for (i=0 ; i<img->width/B8_SIZE ; i++)
          for (j=0 ; j<img->height/MB_BLOCK_SIZE ; j++)
            colB8mode[FRAME][2*j][i] = colB8mode[FRAME][2*j+1][i] = field2frame_mode(colB8mode[TOP_FIELD][j][i]);
      }
      else
      {
        for (i=0 ; i<img->width/B8_SIZE ; i++)
          for (j=0 ; j<img->height/MB_BLOCK_SIZE ; j++)
            colB8mode[TOP_FIELD][j][i] = colB8mode[BOTTOM_FIELD][j][i] = frame2field_mode(colB8mode[FRAME][2*j][i]);
      }
    }
  }
}

/*!
 ************************************************************************
 * \brief
 *    Just a placebo
 ************************************************************************
 */
Boolean dummy_slice_too_big(int bits_slice)
{
  return FALSE;
}


/*!
 ************************************************************************
 * \brief
 *  
 ************************************************************************
 */

void assign_mem2mvs (RD_DATA *var)
{
  get_mem_mv(&(var->mv));
  get_mem_mv(&(var->all_mv));
  get_mem_mv(&(var->all_bmv));
  get_mem_mv(&(var->p_fwMV));
  get_mem_mv(&(var->p_bwMV));
}

/*! 
***************************************************************************
// For MB level field/frame coding
***************************************************************************
*/
void copy_rdopt_data(int bot_block)
{


  int mb_nr = img->current_mb_nr;
  Macroblock *currMB = &img->mb_data[mb_nr];
  int i,j,k,l;
  int     **frefar = ((img->type==B_IMG || img->type==BS_IMG) ? fw_refFrArr : refFrArr);
  int     **brefar = bw_refFrArr;
  int     **frefar_fld;
  int     **brefar_fld = (bot_block ? bw_refFrArr_bot:bw_refFrArr_top); 
  int     bframe   = (img->type==B_IMG || img->type==BS_IMG);
  int     mode;
  int     offset_x, offset_y;
  int     field_y;
  int     block_y = (bot_block ? (img->block_y-4)/2:(img->block_y/2));

  frefar_fld = ((img->type==B_IMG || img->type==BS_IMG) ? (bot_block ? fw_refFrArr_bot :fw_refFrArr_top) : 
                   (bot_block ? refFrArr_bot    :refFrArr_top));

  mode      = rdopt->mode;
  currMB->mb_type = rdopt->mb_type; // copy mb_type 
  currMB->cbp   = rdopt->cbp;   // copy cbp
  currMB->cbp_blk = rdopt->cbp_blk; // copy cbp_blk
  img->i16offset  = rdopt->i16offset;

  if(img->type != B_IMG && img->type != BS_IMG)
    field_mb[img->mb_y][img->mb_x] = MBPairIsField;

  
  if(img->field_mode)
  {
    if(img->top_field && (img->type != B_IMG && img->type != BS_IMG) && mode == 0)
      TopFieldIsSkipped=1;  //set top field MB skipped to 1 for skipped top field MBs
    else if(img->top_field && (img->type == B_IMG || img->type == BS_IMG) && mode == 0)
      TopFieldIsSkipped = currMB->cbp ? 0:1;  // for direct mode check if cbp is zero (skipped) or not
    else if(img->top_field)   // don't change it for bottom field in case TopFieldIsSkipped=1
      TopFieldIsSkipped=0;
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
    if(!bot_block && (img->type != B_IMG && img->type != BS_IMG) && mode == 0)
      TopFrameIsSkipped=1;  //set top frame MB skipped to 1 for skipped top field MBs
    else if(!bot_block && (img->type == B_IMG || img->type == BS_IMG) && mode == 0)
      TopFrameIsSkipped = currMB->cbp ? 0:1;
    else if(!bot_block) 
      TopFrameIsSkipped=0;

/*    if(!bot_block && TopFrameIsSkipped)
      fprintf(stderr,"Skipping Top Frame %d\n", img->current_mb_nr);
    else if(bot_block && mode == 0 && !TopFrameIsSkipped)
      fprintf(stderr,"Skipping Bottom Frame %d\n", img->current_mb_nr);
    else if(bot_block && TopFrameIsSkipped && mode == 0)
      fprintf(stderr,"Skipping Frame Pair %d\n", img->current_mb_nr);
*/
  }

  
  for(i=0;i<6;i++)
    for(j=0;j<4;j++)
      for(k=0;k<2;k++)
        for(l=0;l<18;l++)
          img->cofAC[i][j][k][l]=rdopt->cofAC[i][j][k][l];

  for(i=0;i<3;i++)
    for(k=0;k<2;k++)
      for(l=0;l<18;l++)
        img->cofDC[i][k][l] = rdopt->cofDC[i][k][l];
  
  for(j=0;j<4;j++)
    for(i=0;i<4;i++)
    {
      frefar[img->block_y + j][img->block_x + i] = rdopt->frefar[j][i]; // copy refFrArr
      if(bframe)
        brefar[img->block_y + j][img->block_x + i] = rdopt->brefar[j][i]; // copy bw_refFrArr
    }

  if(!MBPairIsField)
  {
    //===== reconstruction values =====
    for (j=0; j<16; j++)
      for (i=0; i<16; i++)
      {
        imgY[img->pix_y+j][img->pix_x+i] = rdopt->rec_mbY[j][i];    
      }


    for (j=0; j<8; j++)
      for (i=0; i<8; i++)
      {
        imgUV[0][img->pix_c_y+j][img->pix_c_x+i] = rdopt->rec_mbU[j][i];    
        imgUV[1][img->pix_c_y+j][img->pix_c_x+i] = rdopt->rec_mbV[j][i];    
      }
  }
  else 
  {
    if(!bot_block)
      offset_y = img->pix_y;
    else
      offset_y = img->pix_y-MB_BLOCK_SIZE+1;

    //===== reconstruction values =====
    for (j=0; j<16; j++)
      for (i=0; i<16; i++)
      {
        imgY[offset_y+(j<<1)][img->pix_x+i] = rdopt->rec_mbY[j][i];   
      }

    if(!bot_block)
      offset_y = img->pix_c_y;
    else
      offset_y = img->pix_c_y+1-MB_BLOCK_SIZE/2;

    for (j=0; j<8; j++)
      for (i=0; i<8; i++)
      {
        imgUV[0][offset_y+(j<<1)][img->pix_c_x+i] = rdopt->rec_mbU[j][i];   
        imgUV[1][offset_y+(j<<1)][img->pix_c_x+i] = rdopt->rec_mbV[j][i];   
      }
  }

  if(MBPairIsField)
  {
    if(!bot_block)
    {
       offset_x = img->pix_x;
       offset_y = img->pix_y >> 1;
       for(i=0;i<16;i++)
         for(j=0;j<16;j++)
           imgY_top[offset_y+i][offset_x+j] = rdopt->rec_mbY[i][j];

       offset_x = img->pix_c_x;
       offset_y = img->pix_c_y >> 1;

       for(i=0;i<8;i++)
         for(j=0;j<8;j++)
         {
           imgUV_top[0][offset_y+i][offset_x+j] = rdopt->rec_mbU[i][j];
           imgUV_top[1][offset_y+i][offset_x+j] = rdopt->rec_mbV[i][j];
         }
    }
    else
    {
       offset_x = img->pix_x;
       offset_y = (img->pix_y-MB_BLOCK_SIZE) >> 1;
       for(i=0;i<16;i++)
         for(j=0;j<16;j++)
           imgY_bot[offset_y+i][offset_x+j] = rdopt->rec_mbY[i][j];

       offset_x = img->pix_c_x;
       offset_y = (img->pix_c_y-8) >> 1;

       for(i=0;i<8;i++)
         for(j=0;j<8;j++)
         {
           imgUV_bot[0][offset_y+i][offset_x+j] = rdopt->rec_mbU[i][j];
           imgUV_bot[1][offset_y+i][offset_x+j] = rdopt->rec_mbV[i][j];
         }
    }
  }
  else
  {
    if(!bot_block)
    {
       offset_x = img->pix_x;
       offset_y = img->pix_y >> 1;
       for(i=0;i<8;i++)
         for(j=0;j<16;j++)
         {
           imgY_top[offset_y+i][offset_x+j] = rdopt->rec_mbY[2*i][j];
           imgY_bot[offset_y+i][offset_x+j] = rdopt->rec_mbY[2*i+1][j];
         }

       offset_x = img->pix_c_x;
       offset_y = img->pix_c_y >> 1;

       for(i=0;i<4;i++)
         for(j=0;j<8;j++)
         {
           imgUV_top[0][offset_y+i][offset_x+j] = rdopt->rec_mbU[2*i][j];
           imgUV_top[1][offset_y+i][offset_x+j] = rdopt->rec_mbV[2*i][j];
           imgUV_bot[0][offset_y+i][offset_x+j] = rdopt->rec_mbU[2*i+1][j];
           imgUV_bot[1][offset_y+i][offset_x+j] = rdopt->rec_mbV[2*i+1][j];
         }
    }
    else
    {
       offset_x = img->pix_x;
       offset_y = (img->pix_y-MB_BLOCK_SIZE) >> 1;
       offset_y += 8;

       for(i=0;i<8;i++)
         for(j=0;j<16;j++)
         {
           imgY_top[offset_y+i][offset_x+j] = rdopt->rec_mbY[2*i][j];
           imgY_bot[offset_y+i][offset_x+j] = rdopt->rec_mbY[2*i+1][j];
         }

       offset_x = img->pix_c_x;
       offset_y = (img->pix_c_y-8) >> 1;
       offset_y += 4;

       for(i=0;i<4;i++)
         for(j=0;j<8;j++)
         {
           imgUV_top[0][offset_y+i][offset_x+j] = rdopt->rec_mbU[2*i][j];
           imgUV_top[1][offset_y+i][offset_x+j] = rdopt->rec_mbV[2*i][j];
           imgUV_bot[0][offset_y+i][offset_x+j] = rdopt->rec_mbU[2*i+1][j];
           imgUV_bot[1][offset_y+i][offset_x+j] = rdopt->rec_mbV[2*i+1][j];
         }
    }
  }

  for (i=0; i<4; i++)
  {
    currMB->b8mode[i] = rdopt->b8mode[i];
    currMB->b8pdir[i] = rdopt->b8pdir[i];
  }
    //==== reference frames =====
  for (j=0; j<4; j++)
  for (i=0; i<4; i++)
  {
    if(MBPairIsField)
    {
      frefar_fld[block_y+j][img->block_x+i]   = rdopt->frefar[j][i];
      frefar[img->block_y+j][img->block_x+i]    = rdopt->frefar[j][i] == -1 ? -1: rdopt->frefar[j][i]/2;  // Krit, add to match MBINTLC1 in decoder
    }
    else
    {
      frefar_fld[block_y+j][img->block_x+i]   = rdopt->frefar[j][i] == -1 ? -1: 2*rdopt->frefar[j][i];
      frefar[img->block_y+j][img->block_x+i]    = rdopt->frefar[j][i];
    }
  }
  if (bframe)
  {
    for (j=0; j<4; j++)
    for (i=0; i<4; i++)
    {
      if(MBPairIsField)
      {
        brefar_fld[block_y+j][img->block_x+i]   = rdopt->brefar[j][i];
        brefar[img->block_y+j][img->block_x+i]    = rdopt->brefar[j][i] == -1 ? -1: rdopt->brefar[j][i]/2;  // Krit, add to match MBINTLC1 in decoder
      }
      else
      {
        brefar_fld[block_y+j][img->block_x+i]   = rdopt->brefar[j][i] == -1 ? -1: 2*rdopt->brefar[j][i];
        brefar[img->block_y+j][img->block_x+i]    = rdopt->brefar[j][i];
      }
    }
  }
  if (img->type==BS_IMG)
  {
    for (j=0; j<4; j++)
    for (i=0; i<4; i++)
    {
      if(MBPairIsField)
      {
        if (bot_block)
        {
          refFrArr_bot[block_y+j][img->block_x+i]   = rdopt->frefar[j][i];
        }
        else
        {
          refFrArr_top[block_y+j][img->block_x+i]   = rdopt->frefar[j][i];
        }
        refFrArr[img->block_y+j][img->block_x+i]    = rdopt->frefar[j][i] == -1 ? -1: rdopt->frefar[j][i]/2; 
      }
      else
      {
        if (bot_block)
        {
          refFrArr_bot[block_y+j][img->block_x+i]   = rdopt->frefar[j][i] == -1 ? -1: 2*rdopt->frefar[j][i];
        }
        else
        {
          refFrArr_top[block_y+j][img->block_x+i]   = rdopt->frefar[j][i] == -1 ? -1: 2*rdopt->frefar[j][i];
        }
        refFrArr[img->block_y+j][img->block_x+i]    = rdopt->frefar[j][i];
      }
    }
  }



    //==== intra prediction modes ====
  offset_y = MBPairIsField ? (bot_block ? (img->block_y-4)/2:img->block_y/2):img->block_y;
  field_y  = bot_block ? ((img->block_y-4)/2):img->block_y/2;

  if (mode==P8x8)
  {
    for (k=0, j=1; j<5; j++)
    for (     i=img->block_x+1; i<img->block_x+5; i++, k++)
    {
      img->ipredmode[i][img->block_y+j]=rdopt->ipredmode[i][offset_y+j];
      currMB->intra_pred_modes[k] = rdopt->intra_pred_modes[k];
      
      if(!bot_block)
        img->ipredmode_top[i][field_y+j]=rdopt->ipredmode[i][offset_y+j];
      else 
        img->ipredmode_bot[i][field_y+j]=rdopt->ipredmode[i][offset_y+j];
    }
  }
  else if (mode!=I4MB)
  {
    for (k=0, j=1; j<5; j++)
    for (     i=img->block_x+1; i<img->block_x+5; i++, k++)
    {
      img   ->ipredmode    [i][img->block_y+j] = 0;
      currMB->intra_pred_modes[k] = 0;

      if(!bot_block)
        img->ipredmode_top[i][field_y+j]=0;
      else 
        img->ipredmode_bot[i][field_y+j]=0;

    }
  }
  else if(mode == I4MB)
  {
    for (k=0, j=1; j<5; j++)
    for (i=img->block_x+1; i<img->block_x+5; i++, k++)
    {
       img->ipredmode[i][img->block_y+j] = rdopt->ipredmode[i][offset_y+j];
       currMB->intra_pred_modes[k]     = rdopt->intra_pred_modes[k];

      if(!bot_block)
        img->ipredmode_top[i][field_y+j]=rdopt->ipredmode[i][offset_y+j];
      else 
        img->ipredmode_bot[i][field_y+j]=rdopt->ipredmode[i][offset_y+j];
    }

  }

  if(!MBPairIsField && img->mb_y < 2) // to ensure that if the MB is coded as frame 
  {                 // invalid modes like TOP_PRED don't go into the 
    for (k=0, j=1; j<5; j++)    // zero row field MB's
    for (i=img->block_x+1; i<img->block_x+5; i++, k++)
    {
      if(!bot_block)
        img->ipredmode_top[i][field_y+j]=DC_PRED;
      else 
        img->ipredmode_bot[i][field_y+j]=DC_PRED;
    }

  }


   // point the field motion vectors to either top or bottom field 
   // motion vectors
   copy_motion_vectors_MB(bot_block);

} // end of copy_rdopt_data

 void copy_motion_vectors_MB(int bot_block)
 {
  int     ***tmp_mv_fld, ***tmp_fwMV_fld, ***tmp_bwMV_fld;
  int     *****mv_fld, *****all_mv_fld, *****all_bmv_fld;
  int     *****p_fwMV_fld, *****p_bwMV_fld;
  int     mode8,pdir8,l,by,bxr,bx;
  int     **frefar, **refar;
  int     by_f;
  int     i,j,k,ref,ref_field,dref;
  int     bframe = (img->type == B_IMG || img->type == BS_IMG);
  Macroblock* currMB = &img->mb_data[img->current_mb_nr];
  int *****all_mv = img->all_mv;
  int *****all_bmv = img->all_bmv;
  int *****imgmv   = img->mv;
  int *****p_fwMV  = img->p_fwMV;
  int *****p_bwMV  = img->p_bwMV;

  if(!bot_block)
  {
    tmp_mv_fld    = tmp_mv_top;
    tmp_fwMV_fld = tmp_fwMV_top;
    tmp_bwMV_fld = tmp_bwMV_top;
    mv_fld      = img->mv_top;
    all_mv_fld    = img->all_mv_top;
    all_bmv_fld   = img->all_bmv_top;
    p_fwMV_fld    = img->p_fwMV_top;
    p_bwMV_fld    = img->p_bwMV_top;
    frefar      = fw_refFrArr_top;
    refar     = refFrArr_top;
  }
  else 
  {
    tmp_mv_fld    = tmp_mv_bot;
    tmp_fwMV_fld = tmp_fwMV_bot;
    tmp_bwMV_fld = tmp_bwMV_bot;
    mv_fld      = img->mv_bot;
    all_mv_fld    = img->all_mv_bot;
    all_bmv_fld   = img->all_bmv_bot;
    p_fwMV_fld    = img->p_fwMV_bot;
    p_bwMV_fld    = img->p_bwMV_bot;
    frefar      = fw_refFrArr_bot;
    refar     = refFrArr_bot;
  }

  if(MBPairIsField)
  {
    if(!bot_block)
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

  if(!MBPairIsField)
  {
    for(j=0;j<4;j++)
      for(i=0;i<4;i++)
      {
        mode8 = currMB->b8mode[k=2*(j/2)+(i/2)];
        pdir8 = currMB->b8pdir[k];
        l     = 2*(j%2)+(i%2);
        by    = img->block_y+j;
        bxr   = img->block_x+i;
        bx    = img->block_x+i+4;

        if(!bot_block)
          by_f  = img->block_y/2 + j;
        else
          by_f  = (img->block_y - 4)/2 + j;

        ref   = (bframe?fw_refFrArr:refFrArr)[by][bxr];
        ref_field   = (bframe?frefar:refar)[by_f][bxr];

        if (!bframe)
        {
          if (mode8!=IBLOCK && ref != -1)
          {
            tmp_mv   [0][by][bx] = rdopt->tmp_mv[0][j][i];
            tmp_mv   [1][by][bx] = rdopt->tmp_mv[1][j][i];
          }
          else
          {
            tmp_mv   [0][by][bx] = 0;
            tmp_mv   [1][by][bx] = 0;
          }

          tmp_mv_fld [0][by_f][bx] = tmp_mv[0][by][bx];
          tmp_mv_fld [1][by_f][bx] = tmp_mv[1][by][bx]/2;
        }
        else
        {
          if (pdir8==-1) // intra
          {
            tmp_fwMV [0][by][bx] = 0;
            tmp_fwMV [1][by][bx] = 0;
            tmp_bwMV [0][by][bx] = 0;
            tmp_bwMV [1][by][bx] = 0;
            dfMV     [0][by][bx] = 0;
            dfMV     [1][by][bx] = 0;
            dbMV     [0][by][bx] = 0;
            dbMV     [1][by][bx] = 0;
          }
          else if (pdir8==0) // forward
          {
            tmp_fwMV [0][by][bx] = rdopt->tmp_fwMV[0][j][i];
            tmp_fwMV [1][by][bx] = rdopt->tmp_fwMV[1][j][i];
            tmp_bwMV [0][by][bx] = 0;
            tmp_bwMV [1][by][bx] = 0;
            dfMV     [0][by][bx] = 0;
            dfMV     [1][by][bx] = 0;
            dbMV     [0][by][bx] = 0;
            dbMV     [1][by][bx] = 0;
          }
          else if (pdir8==1) // backward
          {
            tmp_fwMV [0][by][bx] = 0;
            tmp_fwMV [1][by][bx] = 0;
            tmp_bwMV [0][by][bx] = rdopt->tmp_bwMV[0][j][i];
            tmp_bwMV [1][by][bx] = rdopt->tmp_bwMV[1][j][i];
            dfMV     [0][by][bx] = 0;
            dfMV     [1][by][bx] = 0;
            dbMV     [0][by][bx] = 0;
            dbMV     [1][by][bx] = 0;
          }
          else if (mode8!=0) // bidirect
          {
            tmp_fwMV [0][by][bx] = rdopt->tmp_fwMV[0][j][i];
            tmp_fwMV [1][by][bx] = rdopt->tmp_fwMV[1][j][i];
            tmp_bwMV [0][by][bx] = rdopt->tmp_bwMV[0][j][i];
            tmp_bwMV [1][by][bx] = rdopt->tmp_bwMV[1][j][i];
            dfMV     [0][by][bx] = 0;
            dfMV     [1][by][bx] = 0;
            dbMV     [0][by][bx] = 0;
            dbMV     [1][by][bx] = 0;
          }
          else // direct
          {
            dref = max(0,refFrArr[by][bxr]);
            tmp_fwMV [0][by][bx] = dfMV     [0][by][bx] = rdopt->dfMV[0][j][i];
            tmp_fwMV [1][by][bx] = dfMV     [1][by][bx] = rdopt->dfMV[1][j][i];
            tmp_bwMV [0][by][bx] = dbMV     [0][by][bx] = rdopt->dbMV[0][j][i];
            tmp_bwMV [1][by][bx] = dbMV     [1][by][bx] = rdopt->dbMV[1][j][i];
          }

          tmp_fwMV_fld [0][by_f][bx]   =  tmp_fwMV [0][by][bx];
          tmp_fwMV_fld [1][by_f][bx]   =  tmp_fwMV [1][by][bx]/2;
          tmp_bwMV_fld [0][by_f][bx]   =  tmp_bwMV [0][by][bx];
          tmp_bwMV_fld [1][by_f][bx]   =  tmp_bwMV [1][by][bx]/2;
          if (img->type==BS_IMG)
          {
            tmp_mv [0][by][bx]       = tmp_fwMV [0][by][bx];
            tmp_mv [1][by][bx]       = tmp_fwMV [1][by][bx];
            tmp_mv_fld [0][by_f][bx] = tmp_fwMV_fld [0][by_f][bx];
            tmp_mv_fld [1][by_f][bx] = tmp_fwMV_fld [1][by_f][bx];
          }
        }
      }
  }
  else
  {
    for(j=0;j<4;j++)
      for(i=0;i<4;i++)
      {
        mode8 = currMB->b8mode[k=2*(j/2)+(i/2)];
        pdir8 = currMB->b8pdir[k];
        l     = 2*(j%2)+(i%2);
        by    = img->block_y+j;
        bxr   = img->block_x+i;
        bx    = img->block_x+i+4;

        if(!bot_block)
          by_f  = img->block_y/2 + j;
        else
          by_f  = (img->block_y - 4)/2 + j;

        ref   = (bframe?fw_refFrArr:refFrArr)[by][bxr];
        ref_field   = (bframe?frefar:refar)[by_f][bxr];

        if (!bframe)
        {
          if (mode8!=IBLOCK && ref != -1)
          {
            tmp_mv_fld   [0][by_f][bx] = rdopt->tmp_mv[0][j][i];
            tmp_mv_fld   [1][by_f][bx] = rdopt->tmp_mv[1][j][i];
          }
          else
          {
            tmp_mv_fld   [0][by_f][bx] = 0;
            tmp_mv_fld   [1][by_f][bx] = 0;
          }

          tmp_mv [0][by][bx] = tmp_mv_fld[0][by_f][bx];
          tmp_mv [1][by][bx] = 2*tmp_mv_fld[1][by_f][bx];
        }
        else
        {
          if (pdir8==-1) // intra
          {
            tmp_fwMV_fld [0][by_f][bx] = 0;
            tmp_fwMV_fld [1][by_f][bx] = 0;
            tmp_bwMV_fld [0][by_f][bx] = 0;
            tmp_bwMV_fld [1][by_f][bx] = 0;
            dfMV     [0][by][bx] = 0;
            dfMV     [1][by][bx] = 0;
            dbMV     [0][by][bx] = 0;
            dbMV     [1][by][bx] = 0;
          }
          else if (pdir8==0) // forward
          {
            tmp_fwMV_fld [0][by_f][bx] = rdopt->tmp_fwMV[0][j][i];
            tmp_fwMV_fld [1][by_f][bx] = rdopt->tmp_fwMV[1][j][i];
            tmp_bwMV_fld [0][by_f][bx] = 0;
            tmp_bwMV_fld [1][by_f][bx] = 0;
            dfMV     [0][by][bx] = 0;
            dfMV     [1][by][bx] = 0;
            dbMV     [0][by][bx] = 0;
            dbMV     [1][by][bx] = 0;
          }
          else if (pdir8==1) // backward
          {
            tmp_fwMV_fld [0][by_f][bx] = 0;
            tmp_fwMV_fld [1][by_f][bx] = 0;
            tmp_bwMV_fld [0][by_f][bx] = rdopt->tmp_bwMV[0][j][i];
            tmp_bwMV_fld [1][by_f][bx] = rdopt->tmp_bwMV[1][j][i];
            dfMV     [0][by][bx] = 0;
            dfMV     [1][by][bx] = 0;
            dbMV     [0][by][bx] = 0;
            dbMV     [1][by][bx] = 0;
          }
          else if (mode8!=0) // bidirect
          {
            tmp_fwMV_fld [0][by_f][bx] = rdopt->tmp_fwMV[0][j][i];
            tmp_fwMV_fld [1][by_f][bx] = rdopt->tmp_fwMV[1][j][i];
            tmp_bwMV_fld [0][by_f][bx] = rdopt->tmp_bwMV[0][j][i];
            tmp_bwMV_fld [1][by_f][bx] = rdopt->tmp_bwMV[1][j][i];
            dfMV     [0][by][bx] = 0;
            dfMV     [1][by][bx] = 0;
            dbMV     [0][by][bx] = 0;
            dbMV     [1][by][bx] = 0;
          }
          else // direct
          {
            dref = max(0,refar[by_f][bxr]);
            tmp_fwMV_fld [0][by_f][bx] = dfMV[0][by][bx] = rdopt->dfMV[0][j][i];
            tmp_fwMV_fld [1][by_f][bx] = dfMV[1][by][bx] = rdopt->dfMV[1][j][i];
            tmp_bwMV_fld [0][by_f][bx] = dbMV[0][by][bx] = rdopt->dbMV[0][j][i];
            tmp_bwMV_fld [1][by_f][bx] = dbMV[1][by][bx] = rdopt->dbMV[1][j][i];
          }

          tmp_fwMV [0][by][bx] =  tmp_fwMV_fld [0][by_f][bx];
          tmp_fwMV [1][by][bx] =  tmp_fwMV_fld [1][by_f][bx]*2;
          tmp_bwMV [0][by][bx] =  tmp_bwMV_fld [0][by_f][bx];
          tmp_bwMV [1][by][bx] =  tmp_bwMV_fld [1][by_f][bx]*2;
          if (img->type==BS_IMG)
          {
            tmp_mv [0][by][bx]       = tmp_fwMV [0][by][bx];
            tmp_mv [1][by][bx]       = tmp_fwMV [1][by][bx];
            tmp_mv_fld [0][by_f][bx] = tmp_fwMV_fld [0][by_f][bx];
            tmp_mv_fld [1][by_f][bx] = tmp_fwMV_fld [1][by_f][bx];
          }
        }
      }
  }

  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      for(k=0;k<img->buf_cycle;k++)
        for(l=0;l<9;l++)
        {
          all_mv[i][j][k][l][0] = rdopt->all_mv[i][j][k][l][0];
          all_bmv[i][j][k][l][0] = rdopt->all_bmv[i][j][k][l][0];
          p_fwMV[i][j][k][l][0] = rdopt->p_fwMV[i][j][k][l][0];
          p_bwMV[i][j][k][l][0] = rdopt->p_bwMV[i][j][k][l][0];
          imgmv[i][j][k][l][0]  = rdopt->mv[i][j][k][l][0];
        
          all_mv[i][j][k][l][1] = rdopt->all_mv[i][j][k][l][1];
          all_bmv[i][j][k][l][1] =rdopt->all_bmv[i][j][k][l][1];
          p_fwMV[i][j][k][l][1] = rdopt->p_fwMV[i][j][k][l][1];
          p_bwMV[i][j][k][l][1] = rdopt->p_bwMV[i][j][k][l][1];
          imgmv[i][j][k][l][1]  = rdopt->mv[i][j][k][l][1];

        }
}


