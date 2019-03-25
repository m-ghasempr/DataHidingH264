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
 ***********************************************************************
 * \file image.c
 *
 * \brief
 *    Decode a Slice
 *
 * \author
 *    Main contributors (see contributors.h for copyright, address and affiliation details)
 *    - Inge Lille-Langøy               <inge.lille-langoy@telenor.com>
 *    - Rickard Sjoberg                 <rickard.sjoberg@era.ericsson.se>
 *    - Jani Lainema                    <jani.lainema@nokia.com>
 *    - Sebastian Purreiter             <sebastian.purreiter@mch.siemens.de>
 *    - Byeong-Moon Jeon                <jeonbm@lge.com>
 *    - Thomas Wedi                     <wedi@tnt.uni-hannover.de>
 *    - Gabi Blaettermann               <blaetter@hhi.de>
 *    - Ye-Kui Wang                     <wyk@ieee.org>
 *    - Antti Hallapuro                 <antti.hallapuro@nokia.com>
 ***********************************************************************
 */

#include "contributors.h"

#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <sys/timeb.h>
#include <string.h>
#include <assert.h>


#include "global.h"
#include "errorconcealment.h"
#include "image.h"
#include "mbuffer.h"
#include "fmo.h"
#include "nalu.h"
#include "parsetcommon.h"
#include "parset.h"
#include "header.h"
#include "rtp.h"
#include "sei.h"

#include "context_ini.h"


#include "erc_api.h"
extern objectBuffer_t *erc_object_list;
extern ercVariables_t *erc_errorVar;
extern frame erc_recfr;
extern int erc_mvperMB;
extern struct img_par *erc_img;

#ifdef _ADAPT_LAST_GROUP_
int *last_P_no;
int *last_P_no_frm;
int *last_P_no_fld;
#endif


/*!
 ***********************************************************************
 * \brief
 *    decodes one I- or P-frame
 *
 ***********************************************************************
 */

int decode_one_frame(struct img_par *img,struct inp_par *inp, struct snr_par *snr)
{
  int current_header;
  Slice *currSlice = img->currentSlice;
  int ercStartMB;
  int ercSegment;
  frame recfr;

  time_t ltime1;                  // for time measurement
  time_t ltime2;

#ifdef WIN32
  struct _timeb tstruct1;
  struct _timeb tstruct2;
#else
  struct timeb tstruct1;
  struct timeb tstruct2;
#endif

  int tmp_time;                   // time used by decoding the last frame


#ifdef WIN32
  _ftime (&tstruct1);             // start time ms
#else
  ftime (&tstruct1);              // start time ms
#endif
  time( &ltime1 );                // start time s

  FmoStartPicture();

  img->current_slice_nr = 0;
  img->current_mb_nr = -4711;     // initialized to an impossible value for debugging -- correct value is taken from slice header
  currSlice->next_header = -8888; // initialized to an impossible value for debugging -- correct value is taken from slice header
  currSlice->next_eiflag = 0;
  while ((currSlice->next_header != EOS && currSlice->next_header != SOP))
  {
// printf ("decode_one_frame: currSlice->next_header %d\n", currSlice->next_header);
    current_header = read_new_slice();

//    img->current_mb_nr = img->map_mb_nr = currSlice->start_mb_nr;//GB

// printf ("Now processing mb %d\n", img->current_mb_nr);
    if (current_header == EOS)
      return EOS;

    if (img->structure == FRAME)
      decode_frame_slice(img, inp, current_header);
    else 
      decode_field_slice(img, inp, current_header);

    img->current_slice_nr++;
  }
  FmoEndPicture();


  //deblocking for frame or top/bottom field
  DeblockFrame( img, imgY, imgUV ) ;
  if(img->structure != FRAME)       //if the previous pict is top or bottom field, 
  {
    FmoStartPicture();
    img->current_slice_nr = 0;
    currSlice->next_header = -8889;
    currSlice->next_eiflag = 0;
    while ((currSlice->next_header != EOS && currSlice->next_header != SOP))
    {   
      // read new slice
      current_header = read_new_slice();
//      img->current_mb_nr = currSlice->start_mb_nr;
      
      if (current_header == EOS)
        return EOS;
      
      decode_field_slice(img, inp, current_header);
      
      img->current_slice_nr++;  
    }
    FmoEndPicture();

    //deblocking bottom/top
    DeblockFrame( img, imgY, imgUV ) ;
  }

  recfr.yptr = &imgY[0][0];
  recfr.uptr = &imgUV[0][0][0];
  recfr.vptr = &imgUV[1][0][0];

  //! this is always true at the beginning of a frame
  ercStartMB = 0;
  ercSegment = 0;

/* !KS: This needs to be fixed for multiple slices
  
  //! mark the start of the first segment
  ercStartSegment(0, ercSegment, 0 , erc_errorVar);
  //! generate the segments according to the macroblock map
  for(i = 1; i<img->max_mb_nr; i++)
  {
    if(img->mb_data[i].ei_flag != img->mb_data[i-1].ei_flag)
    {
      ercStopSegment(i-1, ercSegment, 0, erc_errorVar); //! stop current segment
      
      //! mark current segment as lost or OK
      if(img->mb_data[i-1].ei_flag)
        ercMarkCurrSegmentLost(img->width, erc_errorVar);
      else
        ercMarkCurrSegmentOK(img->width, erc_errorVar);
      
      ercSegment++;  //! next segment
      ercStartSegment(i, ercSegment, 0 , erc_errorVar); //! start new segment
      ercStartMB = i;//! save start MB for this segment 
    }
  }
  //! mark end of the last segent
  ercStopSegment(img->max_mb_nr-1, ercSegment, 0, erc_errorVar);
  if(img->mb_data[i-1].ei_flag)
    ercMarkCurrSegmentLost(img->width, erc_errorVar);
  else
    ercMarkCurrSegmentOK(img->width, erc_errorVar);

  //! call the right error concealment function depending on the frame type.
  erc_mvperMB /= img->max_mb_nr;

  erc_img = img;
  if(img->type == INTRA_IMG || img->type == SI_IMG) // I-frame
    ercConcealIntraFrame(&recfr, img->width, img->height, erc_errorVar);
  else
    ercConcealInterFrame(&recfr, erc_object_list, img->width, img->height, erc_errorVar);
*/
  

  if (img->structure == FRAME)         // buffer mgt. for frame mode
    frame_postprocessing(img, inp);
  else
    field_postprocessing(img, inp);   // reset all interlaced variables

  post_poc( img );                    // POC200301

  if((img->type==B_IMG_1 || img->type==B_IMG_MULT) && !img->disposable_flag)
    copy_stored_B_motion_info(img);

  store_field_MV(img);

  if (p_ref)
    find_snr(snr,img,p_ref,inp->postfilter);      // if ref sequence exist

#ifdef WIN32
  _ftime (&tstruct2);   // end time ms
#else
  ftime (&tstruct2);    // end time ms
#endif
  time( &ltime2 );                                // end time sec
  tmp_time=(ltime2*1000+tstruct2.millitm) - (ltime1*1000+tstruct1.millitm);
  tot_time=tot_time + tmp_time;

  if(img->type == INTRA_IMG) // I picture
    fprintf(stdout,"%3d(I)  %3d %5d %7.4f %7.4f %7.4f %5d\n",
        frame_no, img->tr, img->qp,snr->snr_y,snr->snr_u,snr->snr_v,tmp_time);
  else if(img->type == INTER_IMG_1 || img->type == INTER_IMG_MULT) // P pictures
    fprintf(stdout,"%3d(P)  %3d %5d %7.4f %7.4f %7.4f %5d\n",
    frame_no, img->tr, img->qp,snr->snr_y,snr->snr_u,snr->snr_v,tmp_time);
  else if(img->type == SP_IMG_1 || img->type == SP_IMG_MULT) // SP pictures
    fprintf(stdout,"%3d(SP) %3d %5d %7.4f %7.4f %7.4f %5d\n",
    frame_no, img->tr, img->qp,snr->snr_y,snr->snr_u,snr->snr_v,tmp_time);
  else if (img->type == SI_IMG)
    fprintf(stdout,"%3d(SI) %3d %5d %7.4f %7.4f %7.4f %5d\n",
    frame_no, img->tr, img->qp,snr->snr_y,snr->snr_u,snr->snr_v,tmp_time);
  else if(!img->disposable_flag) // stored B pictures
    fprintf(stdout,"%3d(BS) %3d %5d %7.4f %7.4f %7.4f %5d\n",
        frame_no, img->tr, img->qp,snr->snr_y,snr->snr_u,snr->snr_v,tmp_time);
  else // B pictures
    fprintf(stdout,"%3d(B)  %3d %5d %7.4f %7.4f %7.4f %5d\n",
        frame_no, img->tr, img->qp,snr->snr_y,snr->snr_u,snr->snr_v,tmp_time);

  fflush(stdout);

  if(img->type == INTRA_IMG || img->type == INTER_IMG_1 || img->type == INTER_IMG_MULT) // I or P pictures
    copy_Pframe(img, inp->postfilter);  // imgY-->imgY_prev, imgUV-->imgUV_prev
  else if(img->type == SP_IMG_1 || img->type == SP_IMG_MULT || img->type == SI_IMG) // SP pictures
    copy_Pframe(img, inp->postfilter);  // imgY-->imgY_prev, imgUV-->imgUV_prev
  else if(!img->disposable_flag)  // stored B pictures
    copy_Pframe(img, inp->postfilter);  // imgY-->imgY_prev, imgUV-->imgUV_prev
  else // B pictures
    write_frame(img,inp->postfilter,p_out);         // write image to output YUV file
  
  //! TO 19.11.2001 Known Problem: for init_frame we have to know the picture type of the actual frame
  //! in case the first slice of the P-Frame following the I-Frame was lost we decode this P-Frame but 
  //! do not write it because it was assumed to be an I-Frame in init_frame. So we force the decoder to
  //! guess the right picture type. This is a hack a should be removed by the time there is a clean
  //! solution where we do not have to know the picture type for the function init_frame.
  if(img->type == INTRA_IMG)
    img->type = INTER_IMG_1;
  //! End TO 19.11.2001

  if(img->type <= INTRA_IMG || img->type >= SI_IMG || !img->disposable_flag)   // I or P pictures
    img->number++;
  else
    Bframe_ctr++;    // B pictures

  exit_frame(img, inp);

  if (img->structure != FRAME)
  {
    img->height /= 2;
    img->height_cr /= 2;
  }

  img->current_mb_nr = -4712;   // impossible value for debugging, StW
  img->current_slice_nr = 0;

  img->last_decoded_pic_id = img->tr; // JVT-D101

  return (SOP);
}


/*!
 ************************************************************************
 * \brief
 *    Find PSNR for all three components.Compare decoded frame with
 *    the original sequence. Read inp->jumpd frames to reflect frame skipping.
 ************************************************************************
 */
void find_snr(
  struct snr_par *snr,   //!< pointer to snr parameters
  struct img_par *img,   //!< pointer to image parameters
  FILE *p_ref,           //!< filestream to reference YUV file
  int postfilter)        //!< postfilterin on (=1) or off (=1)
{
  int i,j;
  int diff_y,diff_u,diff_v;
  int uv;
  int  status;
  Slice *currSlice = img->currentSlice;

#ifndef _ADAPT_LAST_GROUP_
  byte       diff;
#endif

#ifndef _ADAPT_LAST_GROUP_
  if(img->type<=INTRA_IMG || img->type >= SI_IMG || !img->disposable_flag) // I, P pictures
    frame_no=img->number*P_interval;
  else // B pictures
  {
    diff=nextP_tr-img->tr;
    frame_no=(img->number-1)*P_interval-diff;
  }
#else
  // TO 5.11.2001 We do have some problems finding the correct frame in the original sequence
  // if errors appear. In this case the method of using this p_frame_no, nextP_tr, prevP_tr
  // variables does not work. So I use the picture_id instead.

  // POC200301 The following modifications are done to make the decoder can get right frame_no
  // in case of more than one IDR pictures. I have not found any reasons to use something like
  // 256*modulo_ctr_xxx.

  //calculate frame number
    if (img->structure == FRAME)
      frame_no = currSlice->picture_id;
    else
      frame_no = currSlice->picture_id/2;
#endif

  rewind(p_ref);
  status = fseek (p_ref, frame_no*img->height*img->width*3/2, 0);
  if (status != 0)
  {
    snprintf(errortext, ET_SIZE, "Error in seeking img->tr: %d", img->tr);
    error(errortext, 500);
  }
  for (j=0; j < img->height; j++)
    for (i=0; i < img->width; i++)
      imgY_ref[j][i]=fgetc(p_ref);
  for (uv=0; uv < 2; uv++)
    for (j=0; j < img->height_cr ; j++)
      for (i=0; i < img->width_cr; i++)
        imgUV_ref[uv][j][i]=fgetc(p_ref);

  img->quad[0]=0;
  diff_y=0;
  for (j=0; j < img->height; ++j)
  {
    for (i=0; i < img->width; ++i)
    {
      diff_y += img->quad[abs(imgY[j][i]-imgY_ref[j][i])];
    }
  }

  // Chroma
  diff_u=0;
  diff_v=0;

  for (j=0; j < img->height_cr; ++j)
  {
    for (i=0; i < img->width_cr; ++i)
    {
      diff_u += img->quad[abs(imgUV_ref[0][j][i]-imgUV[0][j][i])];
      diff_v += img->quad[abs(imgUV_ref[1][j][i]-imgUV[1][j][i])];
    }
  }

  // Collecting SNR statistics
  if (diff_y != 0)
    snr->snr_y=(float)(10*log10(65025*(float)(img->width)*(img->height)/(float)diff_y));        // luma snr for current frame
  if (diff_u != 0)
    snr->snr_u=(float)(10*log10(65025*(float)(img->width)*(img->height)/(float)(4*diff_u)));    //  chroma snr for current frame
  if (diff_v != 0)
    snr->snr_v=(float)(10*log10(65025*(float)(img->width)*(img->height)/(float)(4*diff_v)));    //  chroma snr for current frame

  if (img->number == 0) // first
  {
    snr->snr_y1=(float)(10*log10(65025*(float)(img->width)*(img->height)/(float)diff_y));       // keep luma snr for first frame
    snr->snr_u1=(float)(10*log10(65025*(float)(img->width)*(img->height)/(float)(4*diff_u)));   // keep chroma snr for first frame
    snr->snr_v1=(float)(10*log10(65025*(float)(img->width)*(img->height)/(float)(4*diff_v)));   // keep chroma snr for first frame
    snr->snr_ya=snr->snr_y1;
    snr->snr_ua=snr->snr_u1;
    snr->snr_va=snr->snr_v1;

        if (diff_y == 0)   
      snr->snr_ya=50; // need to assign a reasonable large number so avg snr of entire sequece isn't infinite
        if (diff_u == 0)
      snr->snr_ua=50;
        if (diff_v == 0)
      snr->snr_va=50;

  }
  else
  {
    snr->snr_ya=(float)(snr->snr_ya*(img->number+Bframe_ctr)+snr->snr_y)/(img->number+Bframe_ctr+1); // average snr chroma for all frames
    snr->snr_ua=(float)(snr->snr_ua*(img->number+Bframe_ctr)+snr->snr_u)/(img->number+Bframe_ctr+1); // average snr luma for all frames
    snr->snr_va=(float)(snr->snr_va*(img->number+Bframe_ctr)+snr->snr_v)/(img->number+Bframe_ctr+1); // average snr luma for all frames
  }
}


/*!
 ************************************************************************
 * \brief
 *    Interpolation of 1/4 subpixel
 ************************************************************************
 */
void get_block(int ref_frame,int x_pos, int y_pos, struct img_par *img, int block[BLOCK_SIZE][BLOCK_SIZE])
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

  if (img->structure==FRAME && img->mb_field)
    maxold_y = img->height/2 - 1;

  if (dx == 0 && dy == 0) {  /* fullpel position */
    for (j = 0; j < BLOCK_SIZE; j++)
      for (i = 0; i < BLOCK_SIZE; i++)
        block[i][j] = mref[ref_frame][max(0,min(maxold_y,y_pos+j))][max(0,min(maxold_x,x_pos+i))];
  }
  else { /* other positions */

    if (dy == 0) { /* No vertical interpolation */

      for (j = 0; j < BLOCK_SIZE; j++) {
        for (i = 0; i < BLOCK_SIZE; i++) {
          for (result = 0, x = -2; x < 4; x++)
            result += mref[ref_frame][max(0,min(maxold_y,y_pos+j))][max(0,min(maxold_x,x_pos+i+x))]*COEF[x+2];
          block[i][j] = max(0, min(255, (result+16)/32));
        }
      }

      if ((dx&1) == 1) {
        for (j = 0; j < BLOCK_SIZE; j++)
          for (i = 0; i < BLOCK_SIZE; i++)
            block[i][j] = (block[i][j] + mref[ref_frame][max(0,min(maxold_y,y_pos+j))][max(0,min(maxold_x,x_pos+i+dx/2))])/2;
      }
    }
    else if (dx == 0) {  /* No horizontal interpolation */

      for (j = 0; j < BLOCK_SIZE; j++) {
        for (i = 0; i < BLOCK_SIZE; i++) {
          for (result = 0, y = -2; y < 4; y++)
            result += mref[ref_frame][max(0,min(maxold_y,y_pos+j+y))][max(0,min(maxold_x,x_pos+i))]*COEF[y+2];
          block[i][j] = max(0, min(255, (result+16)/32));
        }
      }

      if ((dy&1) == 1) {
        for (j = 0; j < BLOCK_SIZE; j++)
          for (i = 0; i < BLOCK_SIZE; i++)
           block[i][j] = (block[i][j] + mref[ref_frame][max(0,min(maxold_y,y_pos+j+dy/2))][max(0,min(maxold_x,x_pos+i))])/2;
      }
    }
    else if (dx == 2) {  /* Vertical & horizontal interpolation */

      for (j = -2; j < BLOCK_SIZE+3; j++) {
        for (i = 0; i < BLOCK_SIZE; i++)
          for (tmp_res[i][j+2] = 0, x = -2; x < 4; x++)
            tmp_res[i][j+2] += mref[ref_frame][max(0,min(maxold_y,y_pos+j))][max(0,min(maxold_x,x_pos+i+x))]*COEF[x+2];
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
            block[i][j] = (block[i][j] + max(0, min(255, (tmp_res[i][j+2+dy/2]+16)/32)))/2;
      }
    }
    else if (dy == 2) {  /* Horizontal & vertical interpolation */

      for (j = 0; j < BLOCK_SIZE; j++) {
        for (i = -2; i < BLOCK_SIZE+3; i++)
          for (tmp_res[j][i+2] = 0, y = -2; y < 4; y++)
            tmp_res[j][i+2] += mref[ref_frame][max(0,min(maxold_y,y_pos+j+y))][max(0,min(maxold_x,x_pos+i))]*COEF[y+2];
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
            block[i][j] = (block[i][j] + max(0, min(255, (tmp_res[j][i+2+dx/2]+16)/32)))/2;
      }
    }
    else {  /* Diagonal interpolation */

      for (j = 0; j < BLOCK_SIZE; j++) {
        for (i = 0; i < BLOCK_SIZE; i++) {
          pres_y = dy == 1 ? y_pos+j : y_pos+j+1;
          pres_y = max(0,min(maxold_y,pres_y));
          for (result = 0, x = -2; x < 4; x++)
            result += mref[ref_frame][pres_y][max(0,min(maxold_x,x_pos+i+x))]*COEF[x+2];
          block[i][j] = max(0, min(255, (result+16)/32));
        }
      }

      for (j = 0; j < BLOCK_SIZE; j++) {
        for (i = 0; i < BLOCK_SIZE; i++) {
          pres_x = dx == 1 ? x_pos+i : x_pos+i+1;
          pres_x = max(0,min(maxold_x,pres_x));
          for (result = 0, y = -2; y < 4; y++)
            result += mref[ref_frame][max(0,min(maxold_y,y_pos+j+y))][pres_x]*COEF[y+2];
          block[i][j] = (block[i][j] + max(0, min(255, (result+16)/32))) / 2;
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
int read_new_slice()
{
  NALU_t *nalu = AllocNALU(MAX_CODED_FRAME_SIZE);
  int current_header;
  int ret;
  int BitsUsedByHeader;
  Slice *currSlice = img->currentSlice;
  Bitstream *currStream;
  int newframe;
//  int i;

  while (1)
  {
    if (input->FileFormat == PAR_OF_ANNEXB)
      ret=GetAnnexbNALU (nalu);
    else
      ret=GetRTPNALU (nalu);

    NALUtoRBSP(nalu);
//    printf ("nalu->len %d\n", nalu->len);
    
    if (ret < 0)
      printf ("Error while getting the NALU in file format %s, exit\n", input->FileFormat==PAR_OF_ANNEXB?"Annex B":"RTP");
    if (ret == 0)
    {
//      printf ("read_new_slice: returning %s\n", "EOS");
      FreeNALU(nalu);
      return EOS;
    }

    // Got a NALU
    if (nalu->forbidden_bit)
    {
      printf ("Found NALU w/ forbidden_bit set, bit error?  Let's try...\n");
    }

    switch (nalu->nal_unit_type)
    {
      case NALU_TYPE_SLICE:
      case NALU_TYPE_IDR:
        img->idr_flag = (nalu->nal_unit_type == NALU_TYPE_IDR);
        img->disposable_flag = (nalu->nal_reference_idc == NALU_PRIORITY_DISPOSABLE);
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
        BitsUsedByHeader = FirstPartOfSliceHeader();
        UseParameterSet (currSlice->pic_parameter_set_id);
        BitsUsedByHeader+= RestOfSliceHeader ();
  
        if (img->currentSlice->structure!=0)
        {
          img->height /=2 ;
          img->height_cr /=2;
        }

        // FMO stuff
        //! note: this code does not yet support FMO.  Only the default MBAMap (all zeros,
        //! no scattering) is implemented.  Hence the assert earlier.
        if (img->currentSlice->structure!=0 && img->currentSlice->structure!=3)
        {
          FmoInit (img, input, active_sps->frame_width_in_mbs_minus1+1, (active_sps->frame_height_in_mbs_minus1+1)/2, NULL, 0);   // force a default MBAmap
        }
        else
        {
          FmoInit (img, input, active_sps->frame_width_in_mbs_minus1+1, active_sps->frame_height_in_mbs_minus1+1, NULL, 0);   // force a default MBAmap
        }


        // From here on, active_sps, active_pps and the slice header are valid
//        printf ("img->frame num %d, img->tr %d, img->disposable_flag %d\n", img->frame_num, img->tr, img->disposable_flag);
        img->current_mb_nr = img->map_mb_nr = currSlice->start_mb_nr;  //! potential interlace problem.  Is map_mb_nr correct?

        if (img->tr_old != img->tr)
        {
          newframe=1;
          img->tr_old = img->tr;
        }
        else
          newframe = 0;
        if (newframe)
          current_header = SOP;
        else
          current_header = SOS;

        if(img->structure != img->structure_old)        
          newframe |= 1;

        img->structure_old = img->structure; 
//! new stuff StW
        if(newframe)
          current_header = SOP;
        else
          current_header = SOS;

        if (active_pps->entropy_coding_mode == CABAC)
        {
          int ByteStartPosition = currStream->frame_bitoffset/8;
          if (currStream->frame_bitoffset%8 != 0) 
          {
//            printf ("Slice header ends NOT at a byte aligned position\n");
            ByteStartPosition++;
          }
//          else
//            printf ("SLice header ends at a byte aligned position\n");
// printf ("First CABAC Bytes to decode: %x %x %x %x %x\n", currStream->streamBuffer[ByteStartPosition],currStream->streamBuffer[ByteStartPosition+1],currStream->streamBuffer[ByteStartPosition+2],currStream->streamBuffer[ByteStartPosition+3],currStream->streamBuffer[ByteStartPosition+4],currStream->streamBuffer[ByteStartPosition+5]);

          arideco_start_decoding (&currSlice->partArr[0].de_cabac, currStream->streamBuffer, ByteStartPosition, &currStream->read_len, img->type);
        }
// printf ("read_new_slice: returning %s\n", current_header == SOP?"SOP":"SOS");
        FreeNALU(nalu);
        return current_header;
        break;
      case NALU_TYPE_DPA:
        //! The state machine here should follow the same ideas as the old readSliceRTP()
        //! basically:
        //! work on DPA (as above)
        //! read and process all following SEI/SPS/PPS/PD/Filler NALUs
        //! if next video NALU is dpB, 
        //!   then read and check whether it belongs to DPA, if yes, use it
        //! else
        //!   ;   // nothing
        //! read and process all following SEI/SPS/PPS/PD/Filler NALUs
        //! if next video NALU is dpC
        //!   then read and check whether it belongs to DPA (and DPB, if present), if yes, use it, done
        //! else
        //!   use the DPA (and the DPB if present)

        assert (0==1);
        break;
      case NALU_TYPE_DPB:
        printf ("read_new_slice: Found unexpected NALU_TYPE_DPB, len %d\n", nalu->len);
        printf ("ignoring... moving on with next NALU\n");
        //! Note: one could do something smarter here, e.g. checking the Slice ID
        //! in conjunction with redundant_pic_cnt to identify lost pictures
        assert (0==1);
        break;
      case NALU_TYPE_DPC:
        printf ("read_new_slice: Found NALU_TYPE_DPC, len %d\n", nalu->len);
        printf ("ignoring... moving on with next NALU\n");
        //! Note: one could do something smarter here, e.g. checking the Slice ID
        //! in conjunction with redundant_pic_cnt to identify lost pictures
        assert (0==1);
        break;
      case NALU_TYPE_SEI:
        printf ("read_new_slice: Found NALU_TYPE_SEI, len %d\n", nalu->len);
        InterpretSEIMessage(nalu->buf,nalu->len,img);
        break;
      case NALU_TYPE_PPS:
        ProcessPPS(nalu);
        break;

      case NALU_TYPE_SPS:
        ProcessSPS(nalu);
        break;
      case NALU_TYPE_PD:
        printf ("read_new_slice: Found NALU_TYPE_PD, len %d\n, ignored", nalu->len);
        break;
      case NALU_TYPE_FILL:
        printf ("read_new_slice: Found NALU_TYPE_FILL, len %d\n", nalu->len);
        printf ("Skipping these filling bits, proceeding w/ next NALU\n");
        break;
      default:
        printf ("Found NALU type %d, len %d undefined, ignore NALU, moving on\n", nalu->nal_unit_type, nalu->len);
    }
  }

  FreeNALU(nalu);

  return  current_header;
}


/*!
 ************************************************************************
 * \brief
 *    Initializes the parameters for a new frame
 ************************************************************************
 */
void init_frame(struct img_par *img, struct inp_par *inp)
{
  static int first_P = TRUE;
  int i,j,k,l;

// printf ("init_frame: img->tr %d, img->number %d, img->current_mb_nr %d\n", img->tr, img->number, img->current_mb_nr);
//  img->current_mb_nr=-4713;     // don't know why this should make sense.  
                            // First MB may be in a lost slcie, slices
                            // may be out-of-order... STW
  img->current_slice_nr=0;

//  img->mb_y = img->mb_x = 0;
//  img->block_y = img->pix_y = img->pix_c_y = 0; // define vertical positions
//  img->block_x = img->pix_x = img->pix_c_x = 0; // define horizontal positions

  last_P_no = last_P_no_frm;
  nextP_tr = nextP_tr_frm;

  //WYK: When entire non-B frames are lost, adjust the reference buffers
  //! TO 4.11.2001 Yes, but only for Bitstream mode! We do not loose anything in bitstream mode!
  //! Should remove this one time!

/* !KS removed refPicID from Header
#ifndef AFF //to be fixed later
  if(inp->FileFormat == PAR_OF_ANNEXB) //! TO 4.11.2001 just to make sure that this piece of code 
  {                              //! does not affect any other input mode where this refPicID is not supported
    j = img->refPicID-img->refPicID_old;
    if(j<0) j += 16;    // img->refPicID is 4 bit, wrapps at 15
    if(j > 1) //at least one non-B frame has been lost  
    {
      for(i=1; i<j; i++)  // j-1 non-B frames are lost
      {
        img->number++;
        copy2fb(img);
      }
    }
  }
#endif
*/  
  if (img->number == 0) // first picture
  {
    nextP_tr=prevP_tr=img->tr;
  }
  else if(img->type == INTRA_IMG || img->type == INTER_IMG_1 || img->type == INTER_IMG_MULT || img->type == SP_IMG_1 || img->type == SP_IMG_MULT || img->type == SI_IMG || !img->disposable_flag)
  {
#ifdef _ADAPT_LAST_GROUP_
    for (i = img->buf_cycle-1; i > 0; i--)
      last_P_no[i] = last_P_no[i-1];
    last_P_no[0] = nextP_tr;
#endif
    nextP_tr=img->tr;
    
    if(first_P) // first P picture
    {
      first_P = FALSE;
      P_interval=nextP_tr-prevP_tr; //! TO 4.11.2001 we get problems here in case the first P-Frame was lost
    }
    write_prev_Pframe(img, p_out);  // imgY_prev, imgUV_prev -> file
  }
  
  if (img->type > SI_IMG)
  {
    set_ec_flag(SE_PTYPE);
    img->type = INTER_IMG_1;  // concealed element
  }

  img->max_mb_nr = (img->width * img->height) / (MB_BLOCK_SIZE * MB_BLOCK_SIZE);

  // allocate memory for frame buffers
  if (img->number == 0) 
  {
    init_frame_buffers(inp, img); 
    init_global_buffers(inp, img); 
  }

  for(i=0;i<img->width/BLOCK_SIZE+1;i++)          // set edge to -1, indicate nothing to predict from
  {
    img->ipredmode[i+1][0]=-1;
    img->ipredmode[i+1][img->height/BLOCK_SIZE+1]=-1;
  }
  for(j=0;j<img->height/BLOCK_SIZE+1;j++)
  {
    img->ipredmode[0][j+1]=-1;
    img->ipredmode[img->width/BLOCK_SIZE+1][j+1]=-1;
  }

  for(i=0;i<img->width/BLOCK_SIZE+1;i++)          // set edge to -1, indicate nothing to predict from
    img->ipredmode_frm[i+1][0]=-1;
  for(j=0;j<img->height/BLOCK_SIZE+1;j++)
    img->ipredmode_frm[0][j+1]=-1;

  for(i=0;i<img->width/BLOCK_SIZE+1;i++)          // set edge to -1, indicate nothing to predict from
    img->ipredmode_top[i+1][0]=-1;
  for(j=0;j<(img->height /2)/BLOCK_SIZE+1;j++)
    img->ipredmode_top[0][j+1]=-1;

  for(i=0;i<img->width/BLOCK_SIZE+1;i++)          // set edge to -1, indicate nothing to predict from
    img->ipredmode_bot[i+1][0]=-1;
  for(j=0;j<(img->height /2)/BLOCK_SIZE+1;j++)
    img->ipredmode_bot[0][j+1]=-1;

  // CAVLC init
  for (i=0;i < img->width/MB_BLOCK_SIZE; i++)
    for (j=0; j < img->height/MB_BLOCK_SIZE; j++)
      for (k=0;k<4;k++)
        for (l=0;l<6;l++)
          img->nz_coeff[i][j][k][l]=-1;  // CAVLC

  if(img->constrained_intra_pred_flag)
  {
    for (i=0; i<img->width/MB_BLOCK_SIZE*img->height/MB_BLOCK_SIZE; i++)
    {
      img->intra_block[i][0] =img->intra_block[i][1] = img->intra_block[i][2] = img->intra_block[i][3] = 1;
    }
  }

  // WYK: Oct. 8, 2001. Set the slice_nr member of each MB to -1, to ensure correct when packet loss occurs
  // TO set Macroblock Map (mark all MBs as 'have to be concealed')
  for(i=0; i<img->max_mb_nr; i++)
  {
    img->mb_data[i].slice_nr = -1; 
    img->mb_data[i].ei_flag = 1;
  }

  fb = frm;

  imgY = imgY_frm;
  imgUV = imgUV_frm;

  mref = mref_frm;
  mcef = mcef_frm;

  img->mv = img->mv_frm;
  refFrArr = refFrArr_frm;
  moving_block = moving_block_frm;

  img->fw_refFrArr = img->fw_refFrArr_frm;
  img->bw_refFrArr = img->bw_refFrArr_frm;

//  printf("short size, used, (%d, %d)\n", frm->short_size, frm->short_used );

  // JVT-D097
  if (img->type!=B_IMG_1 && img->type!=B_IMG_MULT && 
      img->num_slice_groups_minus1 == 1 && img->mb_allocation_map_type > 3)
    FmoUpdateEvolvingMBAmap (img, inp, MBAmap);
  // End JVT-D097
}

/*!
 ************************************************************************
 * \brief
 *    exit a frame
 ************************************************************************
 */
void exit_frame(struct img_par *img, struct inp_par *inp)
{
  if(img->type==INTRA_IMG || img->type == INTER_IMG_1 || img->type == INTER_IMG_MULT || img->type == SP_IMG_1 || img->type == SP_IMG_MULT || img->type == SI_IMG || !img->disposable_flag)
    copy2fb(img);

  if (img->structure == FRAME)
  {
    fld->short_used = frm->short_used * 2;
    fld->short_size = frm->short_size * 2;
  }
}

/*!
 ************************************************************************
 * \brief
 *    write the encoding mode and motion vectors of current 
 *    MB to the buffer of the error concealment module.
 ************************************************************************
 */

void ercWriteMBMODEandMV(struct img_par *img,struct inp_par *inp)
{
  extern objectBuffer_t *erc_object_list;
  int i, ii, jj, currMBNum = img->current_mb_nr;
  int mbx = xPosMB(currMBNum,img->width), mby = yPosMB(currMBNum,img->width);
  objectBuffer_t *currRegion, *pRegion;
  Macroblock *currMB = &img->mb_data[currMBNum];
  int***  mv;

  currRegion = erc_object_list + (currMBNum<<2);

  if(img->type != B_IMG_1 && img->type != B_IMG_MULT) //non-B frame
  {
    for (i=0; i<4; i++)
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
        if (currMB->b8mode[i]>=5 && currMB->b8mode[i]<=7)  // SMALL BLOCKS
        {
          ii              = 4*mbx + (i%2)*2 + BLOCK_SIZE;
          jj              = 4*mby + (i/2)*2;
          pRegion->mv[0]  = (img->mv[ii][jj][0] + img->mv[ii+1][jj][0] + img->mv[ii][jj+1][0] + img->mv[ii+1][jj+1][0] + 2)/4;
          pRegion->mv[1]  = (img->mv[ii][jj][1] + img->mv[ii+1][jj][1] + img->mv[ii][jj+1][1] + img->mv[ii+1][jj+1][1] + 2)/4;
        }
        else // 16x16, 16x8, 8x16, 8x8
        {
          pRegion->mv[0]  = img->mv[4*mbx+(i%2)*2+BLOCK_SIZE][4*mby+(i/2)*2][0];
          pRegion->mv[1]  = img->mv[4*mbx+(i%2)*2+BLOCK_SIZE][4*mby+(i/2)*2][1];
        }
        erc_mvperMB      += mabs(pRegion->mv[0]) + mabs(pRegion->mv[1]);
        pRegion->mv[2]    = refFrArr[4*mby+(i/2)*2][4*mbx+(i%2)*2];
      }
    }
  }
  else  //B-frame
  {
    for (i=0; i<4; i++)
    {
      ii                  = 4*mbx + (i%2)*2 + BLOCK_SIZE;
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
        mv                = (currMB->b8mode[i]==0 && currMB->b8pdir[i]==2 ? img->dbMV : currMB->b8pdir[i]==1 ? img->bw_mv : img->fw_mv);
        pRegion->mv[0]    = (mv[ii][jj][0] + mv[ii+1][jj][0] + mv[ii][jj+1][0] + mv[ii+1][jj+1][0] + 2)/4;
        pRegion->mv[1]    = (mv[ii][jj][1] + mv[ii+1][jj][1] + mv[ii][jj+1][1] + mv[ii+1][jj+1][1] + 2)/4;
        erc_mvperMB      += mabs(pRegion->mv[0]) + mabs(pRegion->mv[1]);
        if (currMB->b8pdir[i]==0 || (currMB->b8pdir[i]==2 && currMB->b8mode[i]!=0)) // forward or bidirect
        {
          pRegion->mv[2]  = (img->fw_refFrArr[jj][ii-4]-1+img->buf_cycle) % img->buf_cycle;
          ///???? is it right, not only "img->fw_refFrArr[jj][ii-4]"
        }
        else
        {
          pRegion->mv[2]  = 0;
        }
      }
    }
  }
}


/*!
 ************************************************************************
 * \brief
 *    decodes one slice
 ************************************************************************
 */
void decode_one_slice(struct img_par *img,struct inp_par *inp)
{

  Boolean end_of_slice = FALSE;
  int read_flag;

  img->cod_counter=-1;

  reset_ec_flags();

  while (end_of_slice == FALSE) // loop over macroblocks
  {
    setMapMB_nr (img); //GB


#if TRACE
// Here was the slice nr from the img->mb_data used.  This slice number is only set after 
// the reconstruction of an MB and hence here not yet valid

//    fprintf(p_trace,"\n*********** Pic: %i (I/P) MB: %i Slice: %i Type %d **********\n", img->tr, img->map_mb_nr, img->mb_data[img->map_mb_nr].slice_nr, img->type);
  fprintf(p_trace,"\n*********** Pic: %i (I/P) MB: %i Slice: %i Type %d **********\n", img->tr, img->map_mb_nr, img->current_slice_nr, img->type);

#endif

    // Initializes the current macroblock
    start_macroblock(img,inp, img->current_mb_nr);
    // Get the syntax elements from the NAL
    read_flag = read_one_macroblock(img,inp);

    if (img->mb_frame_field_flag)
      init_super_macroblock(img,inp);

    // decode one macroblock
    if (img->mb_field)
      decode_super_macroblock(img,inp);
    else
      decode_one_macroblock(img,inp);


    if (img->mb_frame_field_flag)
      exit_super_macroblock(img,inp);

    if(img->mb_frame_field_flag && img->mb_field)
      img->num_ref_pic_active_fwd >>= 1;

    ercWriteMBMODEandMV(img,inp);

    end_of_slice=exit_macroblock(img,inp,(!img->mb_frame_field_flag||img->current_mb_nr%2));
  }
  reset_ec_flags();

  if(img->mb_frame_field_flag)
    img->buf_cycle = inp->buf_cycle+1; // reset the img->buf_cycle, otherwise free will cause problems

}


void decode_frame_slice(struct img_par *img,struct inp_par *inp, int current_header)
{
  Slice *currSlice = img->currentSlice;

  if (active_pps->entropy_coding_mode == CABAC)
  {
    init_contexts (img);
  }

  // init new frame
  if (current_header == SOP)
    init_frame(img, inp);

  // do reference frame buffer reordering
  reorder_mref(img);
  if ( (img->weighted_bipred_idc > 0  && (img->type == B_IMG_1 || img->type == B_IMG_MULT)) || (img->weighted_pred_flag && img->type !=INTRA_IMG))
    fill_wp_params(img);


  if (current_header == SOP)
  {
    if (img->number == 0) 
      ercInit(img->width, img->height, 1);
    // reset all variables of the error concealmnet instance before decoding of every frame.
    // here the third parameter should, if perfectly, be equal to the number of slices per frame.
    // using little value is ok, the code will alloc more memory if the slice number is larger
    ercReset(erc_errorVar, img->max_mb_nr, img->max_mb_nr, img->width);
    erc_mvperMB = 0;
  }
    

  // decode main slice information
  if ((current_header == SOP || current_header == SOS) && currSlice->ei_flag == 0)
    decode_one_slice(img,inp);
    
  // setMB-Nr in case this slice was lost
//  if(currSlice->ei_flag)  
//    img->current_mb_nr = currSlice->last_mb_nr + 1;

//! This code doesn't work with FMO or a slice-lossy environment!
//! StW NEEDS FIXING
  if(currSlice->next_eiflag && img->current_mb_nr != img->max_mb_nr)
    currSlice->next_header = SOS;
}



void decode_field_slice(struct img_par *img,struct inp_par *inp, int current_header)
{
  Slice *currSlice = img->currentSlice;

  if (active_pps->entropy_coding_mode == CABAC)
  {
    init_contexts (img);
  }
  
  // init new frame
  if (current_header == SOP)
  {
    if (img->structure == TOP_FIELD)
      init_top(img, inp); // set up field buffer in this function
    else
    {
      init_bottom(img, inp);
    }
  }
  
  // do reference frame buffer reordering
  reorder_mref(img);
  if ( (img->weighted_bipred_idc > 0  && (img->type == B_IMG_1 || img->type == B_IMG_MULT)) || (img->weighted_pred_flag && img->type !=INTRA_IMG))
    fill_wp_params(img);
  

  if (current_header == SOP)
  {
    if (img->number == 0) 
      ercInit(img->width, 2*img->height, 1);
    // reset all variables of the error concealmnet instance before decoding of every frame.
    // here the third parameter should, if perfectly, be equal to the number of slices per frame.
    // using little value is ok, the code will alloc more memory if the slice number is larger
    ercReset(erc_errorVar, img->max_mb_nr, img->max_mb_nr, img->width);
    erc_mvperMB = 0;
  }
  
  // decode main slice information
  if ((current_header == SOP || current_header == SOS) && currSlice->ei_flag == 0)
    decode_one_slice(img,inp);
  
  // setMB-Nr in case this slice was lost
//  if(currSlice->ei_flag)  
//    img->current_mb_nr = currSlice->last_mb_nr + 1;

//! This code doesn't work with FMO or a slice lossy environment or out-of-order slices
  if(currSlice->next_eiflag && img->current_mb_nr != img->max_mb_nr)
    currSlice->next_header = SOS;
}

/*!
 ************************************************************************
 * \brief
 *    Initializes the parameters for a new field
 ************************************************************************
 */
void init_top(struct img_par *img, struct inp_par *inp)
{
  static int first_P = TRUE;
  int i,j;

  img->buf_cycle *= 2;
  img->number *= 2;
//  img->current_mb_nr=-4714;   // impossible value, StW
  img->current_slice_nr=0;

  img->mb_y = img->mb_x = 0;
  img->block_y = img->pix_y = img->pix_c_y = 0; // define vertical positions
  img->block_x = img->pix_x = img->pix_c_x = 0; // define horizontal positions

  last_P_no = last_P_no_fld;
  nextP_tr = nextP_tr_fld;

  //WYK: When entire non-B frames are lost, adjust the reference buffers
  //! TO 4.11.2001 Yes, but only for Bitstream mode! We do not loose anything in bitstream mode!
  //! Should remove this one time!
/* !KS removed refPicID from Header
#ifndef AFF //to be fixed
  if(inp->FileFormat == PAR_OF_ANNEXB) //! TO 4.11.2001 just to make sure that this piece of code 
  {                              //! does not affect any other input mode where this refPicID is not supported
    j = img->refPicID-img->refPicID_old;
    if(j<0) j += 16;    // img->refPicID is 4 bit, wrapps at 15
    if(j > 1) //at least one non-B frame has been lost  
    {
      for(i=1; i<j; i++)  // j-1 non-B frames are lost
      {
        img->number++;
        copy2fb(img);
      }
    }
  }
#endif
*/
  if (img->number == 0) // first picture
  {
    nextP_tr=prevP_tr=img->tr;
  }
  else if(img->type == INTRA_IMG || img->type == INTER_IMG_1 || img->type == INTER_IMG_MULT || img->type == SP_IMG_1 || img->type == SP_IMG_MULT || img->type == SI_IMG || !img->disposable_flag)  // I or P pictures
  {
#ifdef _ADAPT_LAST_GROUP_

    for (i = img->buf_cycle; i > 0; i--)
      last_P_no[i] = last_P_no[i-1];
    last_P_no[0] = nextP_tr;
#endif
    nextP_tr=img->tr;

    if(first_P) // first P picture
    {
      first_P = FALSE;
      P_interval=nextP_tr-prevP_tr; //! TO 4.11.2001 we get problems here in case the first P-Frame was lost
    }
    write_prev_Pframe(img, p_out);  // imgY_prev, imgUV_prev -> file
  }
  
  if (img->type > SI_IMG)
  {
    set_ec_flag(SE_PTYPE);
    img->type = INTER_IMG_1;  // concealed element
  }

  img->max_mb_nr = (img->width * img->height) / (MB_BLOCK_SIZE * MB_BLOCK_SIZE);

  // allocate memory for frame buffers
  if (img->number == 0) 
  {
    init_frame_buffers(inp, img); 
    init_global_buffers(inp, img); 
    img->buf_cycle *= 2;
  }

  for(i=0;i<img->width/BLOCK_SIZE+1;i++)          // set edge to -1, indicate nothing to predict from
  {
    img->ipredmode[i+1][0]=-1;
    img->ipredmode[i+1][img->height/BLOCK_SIZE+1]=-1;
  }
  for(j=0;j<img->height/BLOCK_SIZE+1;j++)
  {
    img->ipredmode[0][j+1]=-1;
    img->ipredmode[img->width/BLOCK_SIZE+1][j+1]=-1;
  }

  if(img->constrained_intra_pred_flag)
  {
    for (i=0; i<img->width/MB_BLOCK_SIZE*img->height/MB_BLOCK_SIZE; i++)
    {
      img->intra_block[i][0] =img->intra_block[i][1] = img->intra_block[i][2] = img->intra_block[i][3] = 1;
    }
  }

  // WYK: Oct. 8, 2001. Set the slice_nr member of each MB to -1, to ensure correct when packet loss occurs
  for(i=0; i<img->max_mb_nr; i++)
    img->mb_data[i].slice_nr = -1;

  fb = fld;

  imgY = imgY_top;
  imgUV = imgUV_top;

  mref = mref_fld;
  mcef = mcef_fld;

  img->mv = img->mv_top;
  refFrArr = refFrArr_top;

  moving_block = moving_block_top;

  img->fw_refFrArr = img->fw_refFrArr_top;
  img->bw_refFrArr = img->bw_refFrArr_top;

  // JVT-D097
  if (img->type!=B_IMG_1 && img->type!=B_IMG_MULT && 
      img->num_slice_groups_minus1 == 1 && img->mb_allocation_map_type > 3)
    FmoUpdateEvolvingMBAmap (img, inp, MBAmap);
  // End JVT-D097
}

/*!
 ************************************************************************
 * \brief
 *    Initializes the parameters for a new field
 ************************************************************************
 */
void init_bottom(struct img_par *img, struct inp_par *inp)
{
  static int first_P = TRUE;
  int i,j;

  img->number++;
//  img->current_mb_nr=-4715;   // impossible value, StW
  img->current_slice_nr=0;
  img->buf_cycle *= 2;

  img->mb_y = img->mb_x = 0;
  img->block_y = img->pix_y = img->pix_c_y = 0; // define vertical positions
  img->block_x = img->pix_x = img->pix_c_x = 0; // define horizontal positions

  last_P_no = last_P_no_fld;

  //WYK: When entire non-B frames are lost, adjust the reference buffers
  //! TO 4.11.2001 Yes, but only for Bitstream mode! We do not loose anything in bitstream mode!
  //! Should remove this one time!
/* !KS removed refPicID from Header  
#ifndef AFF //to be fixed
  if(inp->FileFormat == PAR_OF_ANNEXB) //! TO 4.11.2001 just to make sure that this piece of code 
  {                              //! does not affect any other input mode where this refPicID is not supported
    j = img->refPicID-img->refPicID_old;
    if(j<0) j += 16;    // img->refPicID is 4 bit, wrapps at 15
    if(j > 1) //at least one non-B frame has been lost  
    {
      for(i=1; i<j; i++)  // j-1 non-B frames are lost
      {
        img->number++;
        copy2fb(img);
      }
    }
  }
#endif
*/
  if(img->type==INTRA_IMG || img->type == INTER_IMG_1 || img->type == INTER_IMG_MULT || img->type == SP_IMG_1 || img->type == SP_IMG_MULT || img->type == SI_IMG || !img->disposable_flag)
    copy2fb(img);       // trying to match exit_frame() in frame mode

  if (!img->mb_frame_field_flag)
  {
    if((img->type==B_IMG_1 || img->type==B_IMG_MULT) && !img->disposable_flag)
    {
      // copy motion information of stored B-picture for direct mode 
      for (i=0 ; i<img->width/4+4 ; i++)
      {
        for (j=0 ; j<img->height/4 ; j++)
        {
          img->mv_top[i][j][0] = img->fw_mv[i][j][0];
          img->mv_top[i][j][1] = img->fw_mv[i][j][1];        
          if (i<img->width/4)
          {
            refFrArr_top[j][i] = img->fw_refFrArr[j][i];
          }
        }
      }
    }
  }

  if (img->number == 0) // first picture
  {
    nextP_tr=prevP_tr=img->tr;
  }
  else if (img->type == INTRA_IMG || img->type == INTER_IMG_1 || img->type == INTER_IMG_MULT || img->type == SP_IMG_1 || img->type == SP_IMG_MULT || img->type == SI_IMG || !img->disposable_flag)  // I or P pictures
  {
#ifdef _ADAPT_LAST_GROUP_
    if (img->number==1)
    {
      for (i = img->buf_cycle; i > 0; i--)
        last_P_no[i] = last_P_no[i-1];
      last_P_no[0] = nextP_tr;
    }
#endif
    nextP_tr=img->tr;
    
    if(first_P) // first P picture
    {
      first_P = FALSE;
      P_interval=nextP_tr-prevP_tr; //! TO 4.11.2001 we get problems here in case the first P-Frame was lost
    }
  }

  if (img->type > SI_IMG)
  {
    set_ec_flag(SE_PTYPE);
    img->type = INTER_IMG_1;  // concealed element
  }

  img->max_mb_nr = (img->width * img->height) / (MB_BLOCK_SIZE * MB_BLOCK_SIZE);

  for(i=0;i<img->width/BLOCK_SIZE+1;i++)          // set edge to -1, indicate nothing to predict from
  {
    img->ipredmode[i+1][0]=-1;
    img->ipredmode[i+1][img->height/BLOCK_SIZE+1]=-1;
  }
  for(j=0;j<img->height/BLOCK_SIZE+1;j++)
  {
    img->ipredmode[0][j+1]=-1;
    img->ipredmode[img->width/BLOCK_SIZE+1][j+1]=-1;
  }

  if(img->constrained_intra_pred_flag)
  {
    for (i=0; i<img->width/MB_BLOCK_SIZE*img->height/MB_BLOCK_SIZE; i++)
    {
      img->intra_block[i][0] =img->intra_block[i][1] = img->intra_block[i][2] = img->intra_block[i][3] = 1;
    }
  }

  // WYK: Oct. 8, 2001. Set the slice_nr member of each MB to -1, to ensure correct when packet loss occurs
  for(i=0; i<img->max_mb_nr; i++)
    img->mb_data[i].slice_nr = -1;

  imgY = imgY_bot;
  imgUV = imgUV_bot;

  img->mv = img->mv_bot;
  moving_block = moving_block_bot;
  refFrArr = refFrArr_bot;
  mref = mref_fld;
  mcef = mcef_fld;

  img->fw_refFrArr = img->fw_refFrArr_bot;
  img->bw_refFrArr = img->bw_refFrArr_bot;

  // JVT-D097
  if (img->type!=B_IMG_1 && img->type!=B_IMG_MULT && 
      img->num_slice_groups_minus1 == 1 && img->mb_allocation_map_type > 3)
    FmoUpdateEvolvingMBAmap (img, inp, MBAmap);
  // End JVT-D097
}

/*!
 ************************************************************************
 * \brief
 *    Prepare field and frame buffer after frame decoding
 ************************************************************************
 */
void frame_postprocessing(struct img_par *img, struct inp_par *inp)
{
  int i;

  img->height = img->height/2;
  img->height_cr = img->height_cr/2;
  img->number *= 2;
  img->buf_cycle *= 2;

  fb = fld;
  mref = mref_fld;
  mcef = mcef_fld;
  imgY = imgY_top;
  imgUV = imgUV_top;
  if (img->type == INTRA_IMG || img->type == INTER_IMG_1 || img->type == INTER_IMG_MULT || img->type == SP_IMG_1 || img->type == SP_IMG_MULT || img->type == SI_IMG || !img->disposable_flag)  // I or P pictures
  {
    split_field_top(img);
    copy2fb(img);
  }

  img->number++;
  imgY = imgY_bot;
  imgUV = imgUV_bot;
  if (img->type == INTRA_IMG || img->type == INTER_IMG_1 || img->type == INTER_IMG_MULT || img->type == SP_IMG_1 || img->type == SP_IMG_MULT || img->type == SI_IMG || !img->disposable_flag)  // I or P pictures
  {
    split_field_bot(img);
    copy2fb(img);
  }

  fb = frm;
  mref = mref_frm;
  mcef = mcef_frm;
  imgY = imgY_frm;
  imgUV = imgUV_frm;
  
  img->height *= 2;
  img->height_cr *= 2;
  img->buf_cycle /= 2;
  img->number /= 2;

  if((img->number)&&(img->type==INTRA_IMG || img->type == INTER_IMG_1 || img->type == INTER_IMG_MULT || img->type == SP_IMG_1 || img->type == SP_IMG_MULT || img->type == SI_IMG || !img->disposable_flag))
  {
    for (i = img->buf_cycle; i > 2; i--)
    {
      last_P_no_fld[i] = last_P_no_fld[i-2];
      last_P_no_fld[i-1] = last_P_no_fld[i-3];
    }
    last_P_no_fld[0] = nextP_tr_fld+1;
    last_P_no_fld[1] = nextP_tr_fld;

    nextP_tr_fld = nextP_tr * 2;
    nextP_tr_frm = nextP_tr;
  }
}

/*!
 ************************************************************************
 * \brief
 *    Extract top field from a frame
 ************************************************************************
 */
void split_field_top(struct img_par *img)
{
  int i;

  for (i=0; i<img->height; i++)
  {
    memcpy(imgY[i], imgY_frm[i*2], img->width);
  }

  for (i=0; i<img->height_cr; i++)
  {
    memcpy(imgUV[0][i], imgUV_frm[0][i*2], img->width_cr);
    memcpy(imgUV[1][i], imgUV_frm[1][i*2], img->width_cr);
  }
}

/*!
 ************************************************************************
 * \brief
 *    Extract bottom field from a frame
 ************************************************************************
 */
void split_field_bot(struct img_par *img)
{
  int i;

  for (i=0; i<img->height; i++)
  {
    memcpy(imgY[i], imgY_frm[i*2 + 1], img->width);
  }

  for (i=0; i<img->height_cr; i++)
  {
    memcpy(imgUV[0][i], imgUV_frm[0][i*2 + 1], img->width_cr);
    memcpy(imgUV[1][i], imgUV_frm[1][i*2 + 1], img->width_cr);
  }
}

/*!
 ************************************************************************
 * \brief
 *    Prepare field and frame buffer after field decoding
 ************************************************************************
 */
void field_postprocessing(struct img_par *img, struct inp_par *inp)
{
  int i;

  combine_field(img);

  if(img->type==INTRA_IMG || img->type == INTER_IMG_1 || img->type == INTER_IMG_MULT || img->type == SP_IMG_1 || img->type == SP_IMG_MULT || img->type == SI_IMG || !img->disposable_flag)
  {
    img->number++;
    imgY = imgY_bot;
    imgUV = imgUV_bot;
    copy2fb(img);  // bottom field
    img->number--;
  }
 
  fb = frm;
  mref = mref_frm;
  imgY = imgY_frm;
  imgUV = imgUV_frm;

  if((img->number>1)&&(img->type==INTRA_IMG || img->type == INTER_IMG_1 || img->type == INTER_IMG_MULT || img->type == SP_IMG_1 || img->type == SP_IMG_MULT || img->type == SI_IMG || !img->disposable_flag))
  {
    for (i = img->buf_cycle-1; i > 0; i--)
      last_P_no_frm[i] = last_P_no_frm[i-1];
    last_P_no_frm[0] = nextP_tr_frm;

    nextP_tr_frm = nextP_tr / 2;
    nextP_tr_fld = nextP_tr;
  }

  img->height *= 2;
  img->height_cr *= 2;
  img->buf_cycle /= 2;
  img->number /= 2;
  img->max_mb_nr = (img->width * img->height) / (MB_BLOCK_SIZE * MB_BLOCK_SIZE);
}

/*!
 ************************************************************************
 * \brief
 *    Generate a frame from top and bottom fields
 ************************************************************************
 */
void combine_field(struct img_par *img)
{
  int i;

  for (i=0; i<img->height; i++)
  {
    memcpy(imgY_frm[i*2], imgY_top[i], img->width);     // top field
    memcpy(imgY_frm[i*2 + 1], imgY_bot[i], img->width); // bottom field
  }

  for (i=0; i<img->height_cr; i++)
  {
    memcpy(imgUV_frm[0][i*2], imgUV_top[0][i], img->width_cr);
    memcpy(imgUV_frm[0][i*2 + 1], imgUV_bot[0][i], img->width_cr);
    memcpy(imgUV_frm[1][i*2], imgUV_top[1][i], img->width_cr);
    memcpy(imgUV_frm[1][i*2 + 1], imgUV_bot[1][i], img->width_cr);
  }
}

/*!
 ************************************************************************
 * \brief
 *    Store information for use in B picture
 ************************************************************************
 */
void store_field_MV(struct img_par *img)
{
  int i, j;

  if(img->type==INTRA_IMG || img->type == INTER_IMG_1 || img->type == INTER_IMG_MULT || img->type == SP_IMG_1 || img->type == SP_IMG_MULT || img->type == SI_IMG || !img->disposable_flag)
  {
    if (img->structure != FRAME)
    {
      for (i=0 ; i<img->width/4+4 ; i++)
      {
        for (j=0 ; j<img->height/8 ; j++)
        {
          img->mv_frm[i][2*j][0] = img->mv_frm[i][2*j+1][0] = img->mv_top[i][j][0];
          img->mv_frm[i][2*j][0] = img->mv_frm[i][2*j+1][0] = img->mv_top[i][j][0];
          img->mv_frm[i][2*j][1] = img->mv_frm[i][2*j+1][1] = img->mv_top[i][j][1]*2;
          img->mv_frm[i][2*j][1] = img->mv_frm[i][2*j+1][1] = img->mv_top[i][j][1]*2;

          if (i<img->width/4)
          {          
            moving_block_frm[2*j+1][i]=moving_block_frm[2*j][i]=
              ((refFrArr_top[j][i]!=0) || (refFrArr_bot[j][i]!=0) 
              || (abs(img->mv_top[i+4][j][0])>>1) || (abs(img->mv_top[i+4][j][1])>>1) 
              || (abs(img->mv_bot[i+4][j][0])>>1) || (abs(img->mv_bot[i+4][j][1])>>1));

            moving_block_top[j][i]=((refFrArr_top[j][i]!=0) 
              || (abs(img->mv_top[i+4][j][0])>>1) || (abs(img->mv_top[i+4][j][1])>>1));

            moving_block_bot[j][i]=((refFrArr_bot[j][i]!=0) 
              || (abs(img->mv_bot[i+4][j][0])>>1) || (abs(img->mv_bot[i+4][j][1])>>1));
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
          img->mv_top[i][j][0] = img->mv_bot[i][j][0] = (int)(img->mv_frm[i][2*j][0]);
          img->mv_top[i][j][1] = img->mv_bot[i][j][1] = (int)((img->mv_frm[i][2*j][1])/2);

          if (i<img->width/4)
          {
            moving_block_top[j][i]=moving_block_bot[j][i]=
              ((refFrArr_frm[2*j][i]!=0) || (refFrArr_frm[2*j + 1][i]!=0) 
              || (abs(img->mv_frm[i+4][2*j][0])>>1) || (abs(img->mv_frm[i+4][2*j][1])>>1) 
              || (abs(img->mv_frm[i+4][2*j+1][0])>>1) || (abs(img->mv_frm[i+4][2*j+1][1])>>1));
            
            moving_block_frm[2*j][i]=((refFrArr_frm[2*j][i]!=0) 
              || (abs(img->mv_frm[i+4][2*j][0])>>1) || (abs(img->mv_frm[i+4][2*j][1])>>1) );

            moving_block_frm[2*j+1][i]=((refFrArr_frm[2*j+1][i]!=0) 
              || (abs(img->mv_frm[i+4][2*j+1][0])>>1) || (abs(img->mv_frm[i+4][2*j+1][1])>>1));
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
/*
void store_direct_moving_flag(struct img par *img)
{
  int i, j;  
    if(img->type==INTRA_IMG || img->type == INTER_IMG_1 || img->type == INTER_IMG_MULT || img->type == SP_IMG_1 || img->type == SP_IMG_MULT || img->type == SI_IMG)
    for (i=0 ; i<img->width/4 ; i++)
      for (j=0 ; j<img->height/8 ; j++)
        {                      
          moving_block_frm[2*j][i]=((refFrArr_frm[2*j][i]!=0) 
            || (abs(img->mv_frm[i+4][2*j][0])>>1) || (abs(img->mv_frm[i+4][2*j][1])>>1) );          
          moving_block_frm[2*j+1][i]=((refFrArr_frm[2*j+1][i]!=0) 
            || (abs(img->mv_frm[i+4][2*j+1][0])>>1) || (abs(img->mv_frm[i+4][2*j+1][1])>>1));
        }      
}
*/

void copy_stored_B_motion_info(struct img_par *img)
{
  int i,j;
  // copy motion information of stored B-picture for direct mode 
  if (img->structure != FRAME)
  {
    if (img->mb_frame_field_flag)
    {
      for (i=0 ; i<img->width/4+4 ; i++)
      {
        for (j=0 ; j<img->height/8 ; j++)
        {
          img->mv_top[i][j][0] = img->fw_mv_top[i][j][0];
          img->mv_top[i][j][1] = img->fw_mv_top[i][j][1];        
          img->mv_bot[i][j][0] = img->fw_mv_bot[i][j][0];
          img->mv_bot[i][j][1] = img->fw_mv_bot[i][j][1];        
          if (i<img->width/4)
          {
            refFrArr_top[j][i] = img->fw_refFrArr_top[j][i];
            refFrArr_bot[j][i] = img->fw_refFrArr_bot[j][i];
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
          img->mv_bot[i][j][0] = img->fw_mv[i][j][0];
          img->mv_bot[i][j][1] = img->fw_mv[i][j][1];        
          if (i<img->width/4)
          {
            refFrArr_bot[j][i] = img->fw_refFrArr[j][i];
          }
        }
      }
    }
  }
  else
  {
    for (i=0 ; i<img->width/4+4 ; i++)
    {
      for (j=0 ; j<img->height/4 ; j++)
      {
        img->mv_frm[i][j][0] = img->fw_mv[i][j][0];
        img->mv_frm[i][j][1] = img->fw_mv[i][j][1];        
        if (i<img->width/4)
        {
          refFrArr_frm[j][i] = img->fw_refFrArr[j][i];
        }
      }
    }
  }
}

void reset_wp_params(struct img_par *img)
{
  int i,comp;
  int log_weight_denom;

  for (i=0; i<MAX_REFERENCE_PICTURES; i++)
  {
    for (comp=0; comp<3; comp++)
    {
      log_weight_denom = (comp == 0) ? img->luma_log_weight_denom : img->chroma_log_weight_denom;
      img->wp_weight[0][i][comp] = 1<<log_weight_denom;
      img->wp_weight[1][i][comp] = 1<<log_weight_denom;
    }
  }
}


void fill_wp_params(struct img_par *img)
{
  int i, j, n;
  int comp;
  int log_weight_denom;
  int p0, pt;
//  int p1;
  int bframe = (img->type==B_IMG_1 || img->type==B_IMG_MULT);
  int fwd_ref[MAX_REFERENCE_PICTURES], bwd_ref[MAX_REFERENCE_PICTURES];
  int index;
  int max_bwd_ref, max_fwd_ref;
  int x,z;

  max_bwd_ref = img->num_ref_pic_active_bwd;
  max_fwd_ref = img->num_ref_pic_active_fwd;


if ((img->weighted_bipred_idc > 0) && (img->type == B_IMG_1 || img->type == B_IMG_MULT))
  {
    if (!img->disposable_flag )
    {
      for (index = 0; index < MAX_REFERENCE_PICTURES; index++)
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
    else 
    {
       for (index = 0; index < MAX_REFERENCE_PICTURES - 1; index++)
       {
         fwd_ref[index] = index+1;
       }
       bwd_ref[0] = 0; // only one possible backwards ref for traditional B picture in current software
    }
  }      

  if (img->weighted_bipred_idc == 2 && bframe)
  {
    img->luma_log_weight_denom = 5;
    img->chroma_log_weight_denom = 5;
    img->wp_round_luma = 16;
    img->wp_round_chroma = 16;

    for (i=0; i<MAX_REFERENCE_PICTURES; i++)
    {
      for (comp=0; comp<3; comp++)
      {
        log_weight_denom = (comp == 0) ? img->luma_log_weight_denom : img->chroma_log_weight_denom;
        img->wp_weight[0][i][comp] = 1<<log_weight_denom;
        img->wp_weight[1][i][comp] = 1<<log_weight_denom;
      }
        }
  }

  if (bframe)
  {

    for (i=0; i<max_fwd_ref; i++)
    {
      for (j=0; j<max_bwd_ref; j++)
      {
        for (comp = 0; comp<3; comp++)
        {
          log_weight_denom = (comp == 0) ? img->luma_log_weight_denom : img->chroma_log_weight_denom;
          if (img->weighted_bipred_idc == 1)
          {
            img->wbp_weight[0][i][j][comp] =  img->wp_weight[0][i][comp];
            img->wbp_weight[1][i][j][comp] =  img->wp_weight[1][j][comp];
          }
          else if (img->weighted_bipred_idc == 2)
          {
            pt = poc_distance (fwd_ref[i], bwd_ref[j]);
            if (pt == 0)
            {
              img->wbp_weight[0][i][j][comp] =   32;
              img->wbp_weight[1][i][j][comp] =   32;
            }
            else
            {
              p0 = poc_distance (fwd_ref[i], -1);
 //             p1 = poc_distance (-1, bwd_ref[j]);

			  x = (16384 + (pt>>1))/pt;
			  z = Clip3(-1024, 1023, (x*p0 + 32 )>>6);
              img->wbp_weight[1][i][j][comp] = z >> 2;
              img->wbp_weight[0][i][j][comp] = 64 - img->wbp_weight[1][i][j][comp];
			  if (img->wbp_weight[1][i][j][comp] < -64 || img->wbp_weight[1][i][j][comp] > 128)
			  {
				  img->wbp_weight[1][i][j][comp] = 32;
				  img->wbp_weight[0][i][j][comp] = 32;
			  }


               if (comp == 0)
                      printf ("bpw weight[%d][%d] = %d, %d\n", i,j,
                              img->wbp_weight[0][i][j][0], img->wbp_weight[1][i][j][0]);
            }
          }
        }
     }
   }
 }
}
