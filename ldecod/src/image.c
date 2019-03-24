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
 *    - Ye-Kui Wang                     <wangy@cs.tut.fi>
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
#include "decodeiff.h"
#include "fmo.h"


#if _ERROR_CONCEALMENT_
#include "erc_api.h"
extern objectBuffer_t *erc_object_list;
extern ercVariables_t *erc_errorVar;
extern frame erc_recfr;
extern int erc_mvperMB;
extern struct img_par *erc_img;
#endif

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
  int i;
  Slice *currSlice = img->currentSlice;
  extern FILE* bits;
  int Firstcall;

#if _ERROR_CONCEALMENT_
  int ercStartMB;
  int ercSegment;
  frame recfr;
#endif

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

  if ( inp->of_mode == PAR_OF_IFF ) // for Interim File Format
    rdPictureInfo( bits );

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
    // set the  corresponding read functions
    start_slice(img, inp);

    // read new slice
    Firstcall = (img->max_mb_nr == 0);    // A hack for the current byte string format.  In the current format,
                                          // the picture size is part of the slice header and unknown before
                                          // the first slice is read.  FMO needs to be initialized with the
                                          // picture format.  Hence, test the picture size for initialization
                                          // and initialize FMO in case of the first slice.
    current_header = read_new_slice(img, inp);
    if (Firstcall && (inp->of_mode == PAR_OF_26L || inp->of_mode == PAR_OF_IFF))
    {
      if (currSlice->structure)
      {
        FmoInit (img->width/16, img->height/8, NULL, 0);   // force a default MBAmap
      }
      else
      {
        FmoInit (img->width/16, img->height/16, NULL, 0);   // force a default MBAmap
      }
    }
    img->current_mb_nr = currSlice->start_mb_nr;
#ifdef _ABT_FLAG_IN_SLICE_HEADER_
    USEABT = currSlice->abt;
#endif

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
  
  if(img->structure != FRAME)		//if the previous pict is top or bottom field, 
  {
    FmoStartPicture();
    img->current_slice_nr = 0;
    currSlice->next_header = -8889;
    currSlice->next_eiflag = 0;
    while ((currSlice->next_header != EOS && currSlice->next_header != SOP))
    {	
      // set the  corresponding read functions
      start_slice(img, inp);
      
      // read new slice
      current_header = read_new_slice(img, inp);
      img->current_mb_nr = currSlice->start_mb_nr;
      
      if (current_header == EOS)
        return EOS;
      
      decode_field_slice(img, inp, current_header);
      
      img->current_slice_nr++;	
    }
    FmoEndPicture();
    //deblocking bottom/top
    DeblockFrame( img, imgY, imgUV ) ;
  }

#if _ERROR_CONCEALMENT_
  
  recfr.yptr = &imgY[0][0];
  recfr.uptr = &imgUV[0][0][0];
  recfr.vptr = &imgUV[1][0][0];

  //! this is always true at the beginning of a frame
  ercStartMB = 0;
  ercSegment = 0;
  
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
  if(img->type == INTRA_IMG) // I-frame
    ercConcealIntraFrame(&recfr, img->width, img->height, erc_errorVar);
  else
    ercConcealInterFrame(&recfr, erc_object_list, img->width, img->height, erc_errorVar);
#endif
  
  if (img->structure == FRAME)         // buffer mgt. for frame mode
    frame_postprocessing(img, inp);
  else
    field_postprocessing(img, inp);   // reset all interlaced variables
  store_field_MV(img);
  store_field_colB8mode(img);

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
  else // B pictures
    fprintf(stdout,"%3d(B)  %3d %5d %7.4f %7.4f %7.4f %5d\n",
        frame_no, img->tr, img->qp,snr->snr_y,snr->snr_u,snr->snr_v,tmp_time);

  fflush(stdout);

  if(img->type == INTRA_IMG || img->type == INTER_IMG_1 || img->type == INTER_IMG_MULT) // I or P pictures
    copy_Pframe(img, inp->postfilter);  // imgY-->imgY_prev, imgUV-->imgUV_prev
  else if(img->type == SP_IMG_1 || img->type == SP_IMG_MULT) // SP pictures
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

  if(img->type <= INTRA_IMG || img->type >= SP_IMG_1)   // I or P pictures
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
  static int modulo_ctr_frm=0,modulo_ctr_fld=0,pic_id_old_frm=0,pic_id_old_fld=0;
  static int modulo_ctr_frm_b=0,modulo_ctr_fld_b=0,pic_id_old_frm_b=0,pic_id_old_fld_b=0;
//  static int modulo_ctr = 0;
//  static int modulo_ctr_b = 0;
//  static int modulo_flag = 0;
 // static int modulo_flag_b = 0;
 // static int pic_id_old = 0, pic_id_old_b = 0;
  Slice *currSlice = img->currentSlice;

#ifndef _ADAPT_LAST_GROUP_
  byte       diff;
#endif

#ifndef _ADAPT_LAST_GROUP_
  if(img->type<=INTRA_IMG || img->type >= SP_IMG_1) // I, P pictures
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

  //calculate frame number
  if (img->type <= INTRA_IMG || img->type >= SP_IMG_1) 
  {
    if (img->structure == FRAME)
    { 
      if (currSlice->picture_id < pic_id_old_frm) 
        modulo_ctr_frm++;
      
      pic_id_old_frm = currSlice->picture_id;
      frame_no = currSlice->picture_id + (256*modulo_ctr_frm);
    }
    else
    {
      if (currSlice->picture_id < pic_id_old_fld) 
        modulo_ctr_fld++;
      
      pic_id_old_fld = currSlice->picture_id;
      frame_no = (currSlice->picture_id  + (256*modulo_ctr_fld))/2;
    }
  }
  else
  {
    if (img->structure == FRAME)
    { 
      if (currSlice->picture_id < pic_id_old_frm_b) 
        modulo_ctr_frm_b++;
      
      pic_id_old_frm_b = currSlice->picture_id;
      frame_no = currSlice->picture_id + (256*modulo_ctr_frm_b);
    }
    else
    {
      if (currSlice->picture_id < pic_id_old_fld_b) 
        modulo_ctr_fld_b++;
      
      pic_id_old_fld_b = currSlice->picture_id;
      frame_no = (currSlice->picture_id  + (256*modulo_ctr_fld_b))/2;
    }
  }
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
  {
    snr->snr_y=(float)(10*log10(65025*(float)(img->width)*(img->height)/(float)diff_y));        // luma snr for current frame
  }

  if (diff_u != 0)
  {
    snr->snr_u=(float)(10*log10(65025*(float)(img->width)*(img->height)/(float)(4*diff_u)));    //  chroma snr for current frame
    snr->snr_v=(float)(10*log10(65025*(float)(img->width)*(img->height)/(float)(4*diff_v)));    //  chroma snr for current frame
  }

  if (img->number == 0) // first
  {
    snr->snr_y1=(float)(10*log10(65025*(float)(img->width)*(img->height)/(float)diff_y));       // keep luma snr for first frame
    snr->snr_u1=(float)(10*log10(65025*(float)(img->width)*(img->height)/(float)(4*diff_u)));   // keep chroma snr for first frame
    snr->snr_v1=(float)(10*log10(65025*(float)(img->width)*(img->height)/(float)(4*diff_v)));   // keep chroma snr for first frame
    snr->snr_ya=snr->snr_y1;
    snr->snr_ua=snr->snr_u1;
    snr->snr_va=snr->snr_v1;
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
 *    Direct interpolation of a specific subpel position
 *
 * \author
 *    Thomas Wedi, 12.01.2001   <wedi@tnt.uni-hannover.de>
 *
 * \para Remarks:
 *    A further significant complexity reduction is possible
 *    if the direct interpolation is not performed on pixel basis
 *    but on block basis,
 *
 ************************************************************************
 */
void get_block(int ref_frame,int x_pos, int y_pos, struct img_par *img, int block[BLOCK_SIZE][BLOCK_SIZE])
{

  switch(img->mv_res)
  {
  case 0:

    get_quarterpel_block(ref_frame,x_pos,y_pos,img,block);
    break;

  case 1:

    get_eighthpel_block(ref_frame,x_pos,y_pos,img,block);
    break;

  default:

    snprintf(errortext, ET_SIZE, "wrong mv-resolution: %d",img->mv_res);
    error(errortext, 600);
    break;
  }
}


/*!
 ************************************************************************
 * \brief
 *    Interpolation of 1/4 subpixel
 ************************************************************************
 */
void get_quarterpel_block(int ref_frame,int x_pos, int y_pos, struct img_par *img, int block[BLOCK_SIZE][BLOCK_SIZE])
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

  if (dx == 0 && dy == 0) {  /* fullpel position */
    for (j = 0; j < BLOCK_SIZE; j++)
      for (i = 0; i < BLOCK_SIZE; i++)
        block[i][j] = mref[ref_frame][max(0,min(maxold_y,y_pos+j))][max(0,min(maxold_x,x_pos+i))];
  }
  else if (dx == 3 && dy == 3) {  /* funny position */
    for (j = 0; j < BLOCK_SIZE; j++)
      for (i = 0; i < BLOCK_SIZE; i++)
        block[i][j] = (
          mref[ref_frame][max(0,min(maxold_y,y_pos+j))  ][max(0,min(maxold_x,x_pos+i))  ]+
          mref[ref_frame][max(0,min(maxold_y,y_pos+j))  ][max(0,min(maxold_x,x_pos+i+1))]+
          mref[ref_frame][max(0,min(maxold_y,y_pos+j+1))][max(0,min(maxold_x,x_pos+i+1))]+
          mref[ref_frame][max(0,min(maxold_y,y_pos+j+1))][max(0,min(maxold_x,x_pos+i))  ]+2)/4;
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
 *    Interpolation of 1/8 subpixel
 ************************************************************************
 */

static void get_fullpel_block(int ref_frame,int x_pres, int y_pres, int max_x, int max_y, int block[BLOCK_SIZE][BLOCK_SIZE])
{
  int i, j;

  for (j = 0; j < BLOCK_SIZE; j++)
    for (i = 0; i < BLOCK_SIZE; i++)
      block[i][j] = mref[ref_frame][max(0,min(y_pres+j,max_y))][max(0,min(x_pres+i,max_x))];
}

static void interp_block_X(int ref_frame, int pres_x, int pres_y, const int *coefx, int max_x, int max_y, int block[BLOCK_SIZE][BLOCK_SIZE])
{
  int i, j, x;
  int result;

  for (j = 0; j < BLOCK_SIZE; j++) {
    for (i = 0; i < BLOCK_SIZE; i++) {
      for(result = 0, x = -3; x < 5; x++)
        result += mref[ref_frame][max(0,min(max_y,pres_y+j))][max(0,min(pres_x+i+x,max_x))]*coefx[x+3];
      block[i][j] = max(0, min(255, (result+128)/256));
    }
  }
}

static void interp_block_Y(int ref_frame, int pres_x, int pres_y, const int *coefy, int max_x, int max_y, int block[BLOCK_SIZE][BLOCK_SIZE])
{
  int i, j, y;
  int result;

  for (j = 0; j < BLOCK_SIZE; j++) {
    for (i = 0; i < BLOCK_SIZE; i++) {
      for(result = 0, y = -3; y < 5; y++)
        result += mref[ref_frame][max(0,min(pres_y+j+y,max_y))][max(0,min(max_x,pres_x+i))]*coefy[y+3];
      block[i][j] = max(0, min(255, (result+128)/256));
    }
  }
}

static void interp_block_3_1(int block[BLOCK_SIZE][BLOCK_SIZE],
                             int block1[BLOCK_SIZE][BLOCK_SIZE],
                             int block2[BLOCK_SIZE][BLOCK_SIZE])
{
  int i, j;

  for (j = 0; j < BLOCK_SIZE; j++)
    for (i = 0; i < BLOCK_SIZE; i++)
      block[i][j] = (3*block1[i][j] + block2[i][j] + 2) / 4;
}

static void average_block(int block[BLOCK_SIZE][BLOCK_SIZE], int block2[BLOCK_SIZE][BLOCK_SIZE])
{
  int i, j;

  for (j = 0; j < BLOCK_SIZE; j++)
    for (i = 0; i < BLOCK_SIZE; i++)
      block[i][j] = (block[i][j] + block2[i][j]) / 2;
}


void get_eighthpel_block(int ref_frame,int x_pos, int y_pos, struct img_par *img, int block[BLOCK_SIZE][BLOCK_SIZE])
{
  static const int COEF[3][8] = {
    {-3, 12, -37, 229,  71, -21,  6, -1},
    {-3, 12, -39, 158, 158, -39, 12, -3},
    {-1,  6, -21,  71, 229, -37, 12, -3}
  };

  int block2[BLOCK_SIZE][BLOCK_SIZE];
  int dx=0, x=0;
  int dy=0, y=0;
  int pres_x=0;
  int pres_y=0;
  int max_x=0,max_y=0;

  int tmp[4][4+8]; 
  int result;
  int i, j;


  dx = x_pos&7;
  dy = y_pos&7;
  pres_x = x_pos>>3;
  pres_y = y_pos>>3;
  max_x = img->width-1;
  max_y = img->height-1;


  /* choose filter depending on subpel position */
  if(dx==0 && dy==0) {                /* fullpel position */
    get_fullpel_block(ref_frame, pres_x, pres_y, max_x, max_y, block);
  }
  else if(dy == 0) {  /* Only horizontal interpolation */
    if (dx == 1)
      get_fullpel_block(ref_frame, pres_x, pres_y, max_x, max_y, block);
    else
      interp_block_X(ref_frame, pres_x, pres_y, COEF[dx/2-1], max_x, max_y, block);

    if ((dx&1) != 0) {
      if (dx == 7)
        get_fullpel_block(ref_frame, pres_x+1, pres_y, max_x, max_y, block2);
      else
        interp_block_X(ref_frame, pres_x, pres_y, COEF[dx/2], max_x, max_y, block2);

      average_block(block, block2);
    }
  }
  else if (dx == 0) {  /* Only vertical interpolation */
    if (dy == 1)
      get_fullpel_block(ref_frame, pres_x, pres_y, max_x, max_y, block);
    else
      interp_block_Y(ref_frame, pres_x, pres_y, COEF[dy/2-1], max_x, max_y, block);

    if ((dy&1) != 0) {
      if (dy == 7)
        get_fullpel_block(ref_frame, pres_x, pres_y+1, max_x, max_y, block2);
      else
        interp_block_Y(ref_frame, pres_x, pres_y, COEF[dy/2], max_x, max_y, block2);

      average_block(block, block2);
    }
  }
  else if ((dx&1) == 0) {  /* Horizontal 1/4-pel and vertical 1/8-pel interpolation */

    for (j = -3; j < BLOCK_SIZE+5; j++) {
      for (i = 0; i < BLOCK_SIZE; i++) {
        for (tmp[i][j+3] = 0, x = -3; x < 5; x++)
          tmp[i][j+3] += mref[ref_frame][max(0,min(pres_y+j,max_y))][max(0,min(pres_x+i+x,max_x))]*COEF[dx/2-1][x+3];
      }
    }

    if (dy == 1) {
      for (j = 0; j < BLOCK_SIZE; j++)
        for (i = 0; i < BLOCK_SIZE; i++)
          block[i][j] = max(0, min(255, (tmp[i][j+3]+128)/256));
    }
    else {
      for (j = 0; j < BLOCK_SIZE; j++) {
        for (i = 0; i < BLOCK_SIZE; i++) {
          for (result = 0, y = -3; y < 5; y++)
            result += tmp[i][j+y+3]*COEF[dy/2-1][y+3];
          block[i][j] = max(0, min(255, (result+32768)/65536));
        }
      }
    }

    if ((dy&1) != 0) {
      if (dy == 7) {
        for (j = 0; j < BLOCK_SIZE; j++)
          for (i = 0; i < BLOCK_SIZE; i++)
            block2[i][j] = max(0, min(255, (tmp[i][j+4]+128)/256));
      }
      else {
        for (j = 0; j < BLOCK_SIZE; j++) {
          for (i = 0; i < BLOCK_SIZE; i++) {
            for(result = 0, y = -3; y < 5; y++)
              result += tmp[i][j+y+3]*COEF[dy/2][y+3];
            block2[i][j] = max(0, min(255, (result+32768)/65536));
          }
        }
      }

      average_block(block, block2);
    }
  }
  else if ((dy&1) == 0) {  /* Vertical 1/4-pel and horizontal 1/8-pel interpolation */

    for (j = 0; j < BLOCK_SIZE; j++) {
      for (i = -3; i < BLOCK_SIZE+5; i++) {
        for (tmp[j][i+3] = 0, y = -3; y < 5; y++)
          tmp[j][i+3] += mref[ref_frame][max(0,min(pres_y+j+y,max_y))][max(0,min(pres_x+i,max_x))]*COEF[dy/2-1][y+3];
      }
    }

    if (dx == 1) {
      for (j = 0; j < BLOCK_SIZE; j++)
        for (i = 0; i < BLOCK_SIZE; i++)
          block[i][j] = max(0, min(255, (tmp[j][i+3]+128)/256));
    }
    else {
      for (j = 0; j < BLOCK_SIZE; j++) {
        for (i = 0; i < BLOCK_SIZE; i++) {
          for (result = 0, x = -3; x < 5; x++)
            result += tmp[j][i+x+3]*COEF[dx/2-1][x+3];
          block[i][j] = max(0, min(255, (result+32768)/65536));
        }
      }
    }

    if (dx == 7) {
      for (j = 0; j < BLOCK_SIZE; j++)
        for (i = 0; i < BLOCK_SIZE; i++)
          block2[i][j] = max(0, min(255, (tmp[j][i+4]+128)/256));
    }
    else {
      for (j = 0; j < BLOCK_SIZE; j++) {
        for (i = 0; i < BLOCK_SIZE; i++) {
          for (result = 0, x = -3; x < 5; x++)
            result += tmp[j][i+x+3]*COEF[dx/2][x+3];
          block2[i][j] = max(0, min(255, (result+32768)/65536));
        }
      }
    }

    average_block(block, block2);
  }
  else if ((dx == 1 || dx == 7) && (dy == 1 || dy == 7)) { /* Diagonal averaging */
    interp_block_X(ref_frame, pres_x, (y_pos+1)>>3, COEF[(dx+1)/4], max_x, max_y, block);
    interp_block_Y(ref_frame, (x_pos+1)>>3, pres_y, COEF[(dy+1)/4], max_x, max_y, block2);
    average_block(block, block2);
  }
  else if (dx == 1 || dx == 7 || dy == 1 || dy == 7) { /* Diagonal linear interpolation */

    interp_block_X(ref_frame, pres_x, (y_pos+3)>>3, COEF[1], max_x, max_y, block);
    interp_block_Y(ref_frame, (x_pos+3)>>3, pres_y, COEF[1], max_x, max_y, block2);

    if (dx == 1 || dx == 7)
      interp_block_3_1(block, block2, block);
    else
      interp_block_3_1(block, block, block2);
  }
  else { /* Diagonal interpolation using a full pixel and a center 1/2-pixel */

    for (j = -3; j < BLOCK_SIZE+5; j++) {
      for (i = 0; i < BLOCK_SIZE; i++) {
        for (tmp[i][j+3] = 0, x = -3; x < 5; x++)
          tmp[i][j+3] += mref[ref_frame][max(0,min(pres_y+j,max_y))][max(0,min(pres_x+i+x,max_x))]*COEF[1][x+3];
      }
    }

    for (j = 0; j < BLOCK_SIZE; j++) {
      for (i = 0; i < BLOCK_SIZE; i++) {
        for (result = 0, y = -3; y < 5; y++)
          result += tmp[i][j+y+3]*COEF[1][y+3];
        block[i][j] = max(0, min(255, (result+32768)/65536));
      }
    }

    get_fullpel_block(ref_frame, (x_pos+3)>>3, (y_pos+3)>>3, max_x, max_y, block2);
    interp_block_3_1(block, block, block2);
  }
}

/*!
 ************************************************************************
 * \brief
 *    Reads new slice (picture) from bit_stream
 ************************************************************************
 */
int read_new_slice(struct img_par *img, struct inp_par *inp)
{

    int current_header;
    Slice *currSlice = img->currentSlice;

    // read new slice
    current_header = currSlice->readSlice(img,inp);
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
  int i,j;

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
  
#ifndef AFF //to be fixed later
  if(inp->of_mode == PAR_OF_26L) //! TO 4.11.2001 just to make sure that this piece of code 
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
  
  if (img->number == 0) // first picture
  {
    nextP_tr=prevP_tr=img->tr;
  }
  else if(img->type == INTRA_IMG || img->type == INTER_IMG_1 || img->type == INTER_IMG_MULT || img->type == SP_IMG_1 || img->type == SP_IMG_MULT)  // I or P pictures
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
  
  if (img->type > SP_IMG_MULT)
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

  if(img->UseConstrainedIntraPred)
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

  img->fw_refFrArr = img->fw_refFrArr_frm;
  img->bw_refFrArr = img->bw_refFrArr_frm;

//  printf("short size, used, (%d, %d)\n", frm->short_size, frm->short_used );
}

/*!
 ************************************************************************
 * \brief
 *    exit a frame
 ************************************************************************
 */
void exit_frame(struct img_par *img, struct inp_par *inp)
{
  if(img->type==INTRA_IMG || img->type == INTER_IMG_1 || img->type == INTER_IMG_MULT || img->type == SP_IMG_1 || img->type == SP_IMG_MULT)
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
#if _ERROR_CONCEALMENT_

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
#endif

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

#if TRACE
    fprintf(p_trace,"\n*********** Pic: %i (I/P) MB: %i Slice: %i Type %d **********\n", img->tr, img->current_mb_nr, img->mb_data[img->current_mb_nr].slice_nr, img->type);
#endif

    // Initializes the current macroblock
    start_macroblock(img,inp, img->current_mb_nr);
// printf ("Decode_one_slice: after start_mb (%d)\n", img->current_mb_nr);
    // Get the syntax elements from the NAL
    read_flag = read_one_macroblock(img,inp);
// printf ("Decode_one_slice: after read_mb, currmb %d, currslice %d\n", img->current_mb_nr, img->current_slice_nr);

    // decode one macroblock
    switch(read_flag)
    {
    case DECODE_MB:
      decode_one_macroblock(img,inp);
      break;
    default:
        printf("need to trigger error concealment or something here\n ");
    }

#if _ERROR_CONCEALMENT_
    ercWriteMBMODEandMV(img,inp);
#endif

    end_of_slice=exit_macroblock(img,inp);
// printf (", exit-mb returns %d\n", end_of_slice);
  }
  reset_ec_flags();
}







void decode_frame_slice(struct img_par *img,struct inp_par *inp, int current_header)
{
  Slice *currSlice = img->currentSlice;

  if (inp->symbol_mode == CABAC)
  {
      init_contexts_MotionInfo(img, currSlice->mot_ctx, 1);
      init_contexts_TextureInfo(img, currSlice->tex_ctx, 1);
  }

  // init new frame
  if (current_header == SOP)
    init_frame(img, inp);

  // do reference frame buffer reordering
  reorder_mref(img);


#if _ERROR_CONCEALMENT_
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
  if(currSlice->ei_flag)  
    img->current_mb_nr = currSlice->last_mb_nr + 1;
#else
  // decode main slice information
  if (current_header == SOP || current_header == SOS)
    decode_one_slice(img,inp);
#endif

//! This code doesn't workj with FMO or a slice-lossy environment!

//  if(currSlice->next_eiflag && img->current_mb_nr != img->max_mb_nr)
//    currSlice->next_header = SOS;
}








void decode_field_slice(struct img_par *img,struct inp_par *inp, int current_header)
{
  Slice *currSlice = img->currentSlice;

  if (inp->symbol_mode == CABAC)
  {
    init_contexts_MotionInfo(img, currSlice->mot_ctx, 1);
    init_contexts_TextureInfo(img,currSlice->tex_ctx, 1);
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
  

#if _ERROR_CONCEALMENT_
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
#else

  // decode main slice information
  if (current_header == SOP || current_header == SOS)
    decode_one_slice(img,inp);

#endif

//! This code doesn't work with FMO or a slice lossy environment or out-of-order slices
//  if(currSlice->next_eiflag && img->current_mb_nr != img->max_mb_nr)
//    currSlice->next_header = SOS;
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

#ifndef AFF //to be fixed
  if(inp->of_mode == PAR_OF_26L) //! TO 4.11.2001 just to make sure that this piece of code 
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

  if (img->number == 0) // first picture
  {
    nextP_tr=prevP_tr=img->tr;
  }
  else if(img->type == INTRA_IMG || img->type == INTER_IMG_1 || img->type == INTER_IMG_MULT || img->type == SP_IMG_1 || img->type == SP_IMG_MULT)  // I or P pictures
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
  
  if (img->type > SP_IMG_MULT)
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

  if(img->UseConstrainedIntraPred)
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

  img->fw_refFrArr = img->fw_refFrArr_top;
  img->bw_refFrArr = img->bw_refFrArr_top;
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
  
#ifndef AFF //to be fixed
  if(inp->of_mode == PAR_OF_26L) //! TO 4.11.2001 just to make sure that this piece of code 
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

  if(img->type==INTRA_IMG || img->type == INTER_IMG_1 || img->type == INTER_IMG_MULT || img->type == SP_IMG_1 || img->type == SP_IMG_MULT)
    copy2fb(img);       // trying to match exit_frame() in frame mode

  if (img->number == 0) // first picture
  {
    nextP_tr=prevP_tr=img->tr;
  }
  else if (img->type == INTRA_IMG || img->type == INTER_IMG_1 || img->type == INTER_IMG_MULT || img->type == SP_IMG_1 || img->type == SP_IMG_MULT)  // I or P pictures
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

  if (img->type > SP_IMG_MULT)
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

  if(img->UseConstrainedIntraPred)
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
  refFrArr = refFrArr_bot;
  mref = mref_fld;
  mcef = mcef_fld;

  img->fw_refFrArr = img->fw_refFrArr_bot;
  img->bw_refFrArr = img->bw_refFrArr_bot;
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
  if (img->type == INTRA_IMG || img->type == INTER_IMG_1 || img->type == INTER_IMG_MULT || img->type == SP_IMG_1 || img->type == SP_IMG_MULT)  // I or P pictures
  {
    split_field_top(img);
    copy2fb(img);
  }

  img->number++;
  imgY = imgY_bot;
  imgUV = imgUV_bot;
  if (img->type == INTRA_IMG || img->type == INTER_IMG_1 || img->type == INTER_IMG_MULT || img->type == SP_IMG_1 || img->type == SP_IMG_MULT)  // I or P pictures
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

  if((img->number)&&(img->type==INTRA_IMG || img->type == INTER_IMG_1 || img->type == INTER_IMG_MULT || img->type == SP_IMG_1 || img->type == SP_IMG_MULT))
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

  if(img->type==INTRA_IMG || img->type == INTER_IMG_1 || img->type == INTER_IMG_MULT || img->type == SP_IMG_1 || img->type == SP_IMG_MULT)
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

  if((img->number>1)&&(img->type==INTRA_IMG || img->type == INTER_IMG_1 || img->type == INTER_IMG_MULT || img->type == SP_IMG_1 || img->type == SP_IMG_MULT))
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
  
  if(img->type==INTRA_IMG || img->type == INTER_IMG_1 || img->type == INTER_IMG_MULT || img->type == SP_IMG_1 || img->type == SP_IMG_MULT)
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

// ABT
int  field2frame_mode(int fld_mode)
{
  int frm_mode;
  static const int field2frame_map[MAXMODE] = {-1, 1,1,3,3,4,6,6, -1, 7,1, -1};

  frm_mode = field2frame_map[fld_mode];
  assert(frm_mode>0);

  return frm_mode;

}
int frame2field_mode(int frm_mode)
{
  int fld_mode;
  static const int frame2field_map[MAXMODE] = {-1, 2,5,4,5,5,7,7, -1, 7,2, -1};

  fld_mode = frame2field_map[frm_mode];
  assert(fld_mode>0);

  return fld_mode;

}
void store_field_colB8mode(struct img_par *img)
{
  int i, j;

  if (USEABT)
  {
    if(img->type==INTRA_IMG || img->type == INTER_IMG_1 || img->type == INTER_IMG_MULT || img->type == SP_IMG_1 || img->type == SP_IMG_MULT)
    {
      if (img->structure != FRAME)
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
