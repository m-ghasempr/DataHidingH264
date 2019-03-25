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
 *  \mainpage
 *     This is the H.26L encoder reference software. For detailed documentation
 *     see the comments in each file.
 *
 *  \author
 *     The main contributors are listed in contributors.h
 *
 *  \version
 *     JM 6.0a
 *
 *  \note
 *     tags are used for document system "doxygen"
 *     available at http://www.doxygen.org
 */
/*!
 *  \file
 *     lencod.c
 *  \brief
 *     TML encoder project main
 *  \author
 *   Main contributors (see contributors.h for copyright, address and affiliation details)
 *   - Inge Lille-Langøy               <inge.lille-langoy@telenor.com>
 *   - Rickard Sjoberg                 <rickard.sjoberg@era.ericsson.se>
 *   - Stephan Wenger                  <stewe@cs.tu-berlin.de>
 *   - Jani Lainema                    <jani.lainema@nokia.com>
 *   - Byeong-Moon Jeon                <jeonbm@lge.com>
 *   - Yoon-Seong Soh                  <yunsung@lge.com>
 *   - Thomas Stockhammer              <stockhammer@ei.tum.de>
 *   - Detlev Marpe                    <marpe@hhi.de>
 *   - Guido Heising                   <heising@hhi.de>
 *
 ***********************************************************************
 */

#include "contributors.h"

#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/timeb.h>
#include <stdlib.h>
#if defined WIN32
  #include <conio.h>
#endif
#include <assert.h>

#include "global.h"
#include "configfile.h"
#include "leaky_bucket.h"
#include "memalloc.h"
#include "mbuffer.h"
#include "encodeiff.h"
#include "intrarefresh.h"
#include "fmo.h"
#include "sei.h"

#define JM      "6"
#define VERSION "6.0a"

InputParameters inputs, *input = &inputs;
ImageParameters images, *img   = &images;
StatParameters  stats,  *stat  = &stats;
SNRParameters   snrs,   *snr   = &snrs;
Decoders decoders, *decs=&decoders;


#ifdef _ADAPT_LAST_GROUP_
int initial_Bframes = 0;
#endif

Boolean In2ndIGOP = FALSE;
int    start_frame_no_in_this_IGOP = 0;
int    start_tr_in_this_IGOP = 0;
int    FirstFrameIn2ndIGOP=0;


/*!
 ***********************************************************************
 * \brief
 *    Main function for encoder.
 * \param argc
 *    number of command line arguments
 * \param argv
 *    command line arguments
 * \return
 *    exit code
 ***********************************************************************
 */
void Init_Motion_Search_Module ();
void Clear_Motion_Search_Module ();

int main(int argc,char **argv)
{
  int image_type;
  int num_ref_pic_bak;
  
  p_dec = p_dec_u = p_dec_v = p_stat = p_log = p_datpart = p_trace = NULL;


//  isBigEndian = testEndian();

  Configure (argc, argv);
  init_img();
  AllocNalPayloadBuffer();

  frame_pic = malloc_picture();
  if (input->InterlaceCodingOption != FRAME_CODING)
  {
    top_pic = malloc_picture();
    bottom_pic = malloc_picture();
  }
  init_rdopt ();

  // allocate memory for frame buffers
  init_frame_buffers(input,img);
  init_global_buffers();

  Init_Motion_Search_Module ();

  information_init();

  // B pictures
  Bframe_ctr=0;
  img->refPicID=-1;  //WYK: Oct. 16, 2001
  img->refPicID_frm = -1;
  img->refPicID_fld = -1;
  tot_time=0;                 // time for total encoding session

#ifdef _ADAPT_LAST_GROUP_
  if (input->last_frame > 0)
    input->no_frames = 1 + (input->last_frame + input->jumpd) / (input->jumpd + 1);
  initial_Bframes = input->successive_Bframe;
#endif

  PatchInputNoFrames();
  img->pn = -1;    // Tian
  init_poc();

  // Write sequence header (with parameter sets)
  stat->bit_slice = start_sequence();
  start_frame_no_in_this_IGOP = 0;
  for (img->number=0; img->number < input->no_frames; img->number++)
  {
    img->disposable_flag = FALSE;
    img->pn = (img->pn+1) % (img->buf_cycle+1); 

    //much of this can go in init_frame() or init_field()?
    //poc for this frame or field
    img->toppoc = IMG_NUMBER*img->offset_for_ref_frame[0];
    if (input->InterlaceCodingOption==FRAME_CODING)
      img->bottompoc = img->toppoc;     //progressive
    else 
      img->bottompoc = img->toppoc+1;
    push_poc(img->toppoc, img->bottompoc, REFFRAME);               //push poc values into array
    
    //frame_num for this frame
    img->frame_num = IMG_NUMBER;
    img->idr_flag = !IMG_NUMBER;            //for 2nd field, idr is reset in SliceHeader() //5.0f
    
    //the following is sent in the slice header
    img->delta_pic_order_cnt[0]=0;
    
    SetImgType();

#ifdef _ADAPT_LAST_GROUP_
    if (input->successive_Bframe && input->last_frame && IMG_NUMBER+1 == input->no_frames)
    {                                           
      int bi = (int)((float)(input->jumpd+1)/(input->successive_Bframe+1.0)+0.499999);
      
      input->successive_Bframe = (input->last_frame-(img->number-1)*(input->jumpd+1))/bi-1;

      //about to code the last ref frame, adjust deltapoc         
      img->delta_pic_order_cnt[0]= -2*(initial_Bframes - input->successive_Bframe);
      img->toppoc += img->delta_pic_order_cnt[0];
      img->bottompoc += img->delta_pic_order_cnt[0];
      toprefpoc[0] = img->toppoc;                             //fix array values too
      bottomrefpoc[0] = img->bottompoc;
      
    }
#endif

    image_type = img->type;
    
    img->num_ref_pic_active_fwd_minus1 = max(0,img->nb_references - 1);
    num_ref_pic_bak          = img->num_ref_pic_active_fwd_minus1;
    if (input->StoredBPictures > 0)
    {
      if (input->InterlaceCodingOption>=MB_CODING)
      {
        img->num_ref_pic_active_bwd_minus1 = 1;
      }
      else
      {
        img->num_ref_pic_active_bwd_minus1 = 1;
      }
    }
    else
    {
      img->num_ref_pic_active_bwd_minus1 = 0;
    }
  
    // which layer the image belonged to?
    if ( IMG_NUMBER % (input->NumFramesInELSubSeq+1) == 0 )
      img->layer = 0;
    else
      img->layer = 1;



    encode_one_frame(); // encode one I- or P-frame


    img->nb_references += 1;
    img->nb_references = min(img->nb_references, img->buf_cycle +1); // Tian Dong. PLUS1, +1, June 7, 2002
    img->nb_references = min(img->nb_references, fb->short_used+fb->long_used);// Tian Dong. June 7, 2002

    if ((input->successive_Bframe != 0) && (IMG_NUMBER > 0)) // B-frame(s) to encode
    {
      img->type = B_IMG;            // set image type to B-frame
      img->num_ref_pic_active_bwd_minus1 = 0;

      img->types= INTER_IMG;

      if (input->NumFramesInELSubSeq == 0) 
        img->layer = 0;
      else 
        img->layer = 1;

      img->frame_num++;                 //increment frame_num once for B-frames
      img->idr_flag=0;                                        //B-frames cannot be idr  //5.0f
      img->disposable_flag = TRUE;      // These B frames are disposable (are the non-disposable B frames implemented???)

      for(img->b_frame_to_code=1; img->b_frame_to_code<=input->successive_Bframe; img->b_frame_to_code++)
      {
        img->num_ref_pic_active_fwd_minus1 = num_ref_pic_bak;       // coding field pictures would have changed this, so reset

        //! somewhere here the disposable flag was set -- B frames are always disposable in this encoder.
        //! This happens now in slice.c, terminate_slice, where the nal_reference_idc is set up
        //poc for this B frame
        img->toppoc = (IMG_NUMBER-1)*img->offset_for_ref_frame[0] + 2* img->b_frame_to_code;
        if (input->InterlaceCodingOption==FRAME_CODING)
          img->bottompoc = img->toppoc;     //progressive
        else 
          img->bottompoc = img->toppoc+1;
        
        push_poc(img->toppoc,img->bottompoc, NONREFFRAME);      //push pocs for B frame
        
        //the following is sent in the slice header
        img->delta_pic_order_cnt[0]= 2*(img->b_frame_to_code-1);
        
        encode_one_frame();  // encode one B-frame
      }
    }

    process_2nd_IGOP();
  }
  // terminate sequence
  terminate_sequence();

  fclose(p_in);
  if (p_dec)
    fclose(p_dec);
  if (p_trace)
    fclose(p_trace);

  Clear_Motion_Search_Module ();

  RandomIntraUninit();
  FmoUninit();
  
  // free structure for rd-opt. mode decision
  clear_rdopt ();

#ifdef _LEAKYBUCKET_
  calc_buffer();
#endif

  // report everything
  report();



  free_picture (frame_pic);
  if (top_pic)
    free_picture (top_pic);
  if (bottom_pic)
    free_picture (bottom_pic);

  // free allocated memory for frame buffers
  free_frame_buffers(input,img);
  free_global_buffers();

  // free image mem
  free_img ();
  FreeNalPayloadBuffer();
  return 0;
}


/*!
 ***********************************************************************
 * \brief
 *    Initializes the POC structure with appropriate parameters.
 * 
 ***********************************************************************
 */
void init_poc()
{
  int i;

  if(input->no_frames > (1<<15))error("too many frames",-998);
  if(input->no_multpred >= MAX_NO_POC_FRAMES)error("NumberReferenceFrames too large",-999);
  for(i=0; i<MAX_NO_POC_FRAMES; i++){toprefpoc[i] = bottomrefpoc[i] = 1<<29;}            //init with large 

  //the following should probably go in sequence parameters
  // frame poc's increase by 2, field poc's by 1

  // variable below gone for good, see defines.h LOG2_MAX_FRAME_NUM_MINUS4
  //  img->log2_max_frame_num_minus4 = 4;           // used to be 12, bug?;      
  img->pic_order_cnt_type=1;
  img->num_ref_frames_in_pic_order_cnt_cycle=1;
  img->delta_pic_order_always_zero_flag=0;
  img->offset_for_non_ref_pic =  -2*(input->successive_Bframe);

  if (input->InterlaceCodingOption==FRAME_CODING)
    img->offset_for_top_to_bottom_field=0;
  else    
    img->offset_for_top_to_bottom_field=1;
  img->offset_for_ref_frame[0] = 2*(input->successive_Bframe+1);

                                    //the following should probably go in picture parameters
  img->pic_order_present_flag=0;    //img->delta_pic_order_cnt[1] not sent
}

/*!
 ***********************************************************************
 * \brief
 *    Pushes values onto the POC ref arrays  toprefpoc[] & bottomrefpoc[] 
 * 
 ***********************************************************************
 */
void push_poc(unsigned int topvalue, unsigned int bottomvalue, unsigned int ref_frame_ind )
{
  int i;
  static int current_is_ref = 0;                //indicates if current top value is for a ref frame

         if(current_is_ref){
                 for(i=MAX_NO_POC_FRAMES-1; i>0; i--){          //move all the data down by one
                        toprefpoc[i] = toprefpoc[i-1] ; 
                        bottomrefpoc[i] = bottomrefpoc[i-1] ;
                 }
         }      
        
        toprefpoc[0] = topvalue;                //put new data
        bottomrefpoc[0] = bottomvalue;

        current_is_ref = ref_frame_ind;         //new data type

}


/*!
 ***********************************************************************
 * \brief
 *    Initializes the img->nz_coeff
 * \par Input:
 *    none
 * \par  Output:
 *    none
 * \ side effects
 *    sets omg->nz_coef[][][][] to -1
 ***********************************************************************
 */
void CAVLC_init()
{
  int i, j, k, l;

  for (i=0;i < img->width/MB_BLOCK_SIZE; i++)
    for (j=0; j < img->height/MB_BLOCK_SIZE; j++)
      for (k=0;k<4;k++)
        for (l=0;l<6;l++)
          img->nz_coeff[i][j][k][l]=-1;


}


/*!
 ***********************************************************************
 * \brief
 *    Initializes the Image structure with appropriate parameters.
 * \par Input:
 *    Input Parameters struct inp_par *inp
 * \par  Output:
 *    Image Parameters struct img_par *img
 ***********************************************************************
 */
void init_img()
{
  int i,j,size_x,size_y;

  img->no_multpred=input->no_multpred;
#ifdef _ADDITIONAL_REFERENCE_FRAME_
  img->buf_cycle = max (input->no_multpred, input->add_ref_frame+1);
#else
  img->buf_cycle = input->no_multpred;
#endif

  img->lindex=0;
  img->max_lindex=0;

  img->framerate=INIT_FRAME_RATE;   // The basic frame rate (of the original sequence)
  if(input->InterlaceCodingOption != FRAME_CODING) 
    img->buf_cycle *= 2;

  get_mem_mv (&(img->mv));
  get_mem_mv (&(img->p_fwMV));
  get_mem_mv (&(img->p_bwMV));
  get_mem_mv (&(img->all_mv));
  get_mem_mv (&(img->all_bmv));
  get_mem_mv (&(img->abp_all_dmv));

  get_mem_ACcoeff (&(img->cofAC));
  get_mem_DCcoeff (&(img->cofDC));

  if(input->InterlaceCodingOption >= MB_CODING) 
  {
    get_mem_mv (&(img->mv_top));
    get_mem_mv (&(img->p_fwMV_top));
    get_mem_mv (&(img->p_bwMV_top));
    get_mem_mv (&(img->all_mv_top));
    get_mem_mv (&(img->all_bmv_top));

    get_mem_mv (&(img->mv_bot));
    get_mem_mv (&(img->p_fwMV_bot));
    get_mem_mv (&(img->p_bwMV_bot));
    get_mem_mv (&(img->all_mv_bot));
    get_mem_mv (&(img->all_bmv_bot));

    get_mem_mv (&(rddata_top_frame_mb.mv));
    get_mem_mv (&(rddata_top_frame_mb.p_fwMV));
    get_mem_mv (&(rddata_top_frame_mb.p_bwMV));
    get_mem_mv (&(rddata_top_frame_mb.all_mv));
    get_mem_mv (&(rddata_top_frame_mb.all_bmv));

    get_mem_mv (&(rddata_bot_frame_mb.mv));
    get_mem_mv (&(rddata_bot_frame_mb.p_fwMV));
    get_mem_mv (&(rddata_bot_frame_mb.p_bwMV));
    get_mem_mv (&(rddata_bot_frame_mb.all_mv));
    get_mem_mv (&(rddata_bot_frame_mb.all_bmv));

    get_mem_mv (&(rddata_top_field_mb.mv));
    get_mem_mv (&(rddata_top_field_mb.p_fwMV));
    get_mem_mv (&(rddata_top_field_mb.p_bwMV));
    get_mem_mv (&(rddata_top_field_mb.all_mv));
    get_mem_mv (&(rddata_top_field_mb.all_bmv));

    get_mem_mv (&(rddata_bot_field_mb.mv));
    get_mem_mv (&(rddata_bot_field_mb.p_fwMV));
    get_mem_mv (&(rddata_bot_field_mb.p_bwMV));
    get_mem_mv (&(rddata_bot_field_mb.all_mv));
    get_mem_mv (&(rddata_bot_field_mb.all_bmv));
    get_mem_mv (&(img->abp_all_dmv_top));
    get_mem_mv (&(img->abp_all_dmv_bot));

    get_mem_ACcoeff (&(rddata_top_frame_mb.cofAC));
    get_mem_DCcoeff (&(rddata_top_frame_mb.cofDC));

    get_mem_ACcoeff (&(rddata_bot_frame_mb.cofAC));
    get_mem_DCcoeff (&(rddata_bot_frame_mb.cofDC));

    get_mem_ACcoeff (&(rddata_top_field_mb.cofAC));
    get_mem_DCcoeff (&(rddata_top_field_mb.cofDC));

    get_mem_ACcoeff (&(rddata_bot_field_mb.cofAC));
    get_mem_DCcoeff (&(rddata_bot_field_mb.cofDC));
  }

  if(input->InterlaceCodingOption != FRAME_CODING) 
    img->buf_cycle /= 2;

  if ((img->quad = (int*)calloc (511, sizeof(int))) == NULL)
    no_mem_exit ("init_img: img->quad");
  img->quad+=255;
  for (i=0; i < 256; ++i) // fix from TML1 / TML2 sw, truncation removed
  {
    img->quad[i]=img->quad[-i]=i*i;
  }

  img->width    = input->img_width;
  img->height   = input->img_height;
  img->width_cr = input->img_width/2;
  img->height_cr= input->img_height/2;

  if(((img->mb_data) = (Macroblock *) calloc((img->width/MB_BLOCK_SIZE) * (img->height/MB_BLOCK_SIZE),sizeof(Macroblock))) == NULL)
    no_mem_exit("init_img: img->mb_data");

  if(input->UseConstrainedIntraPred)
  {
    if(((img->intra_block) = (int**)calloc((j=(img->width/MB_BLOCK_SIZE) * (img->height/MB_BLOCK_SIZE)),sizeof(int))) == NULL)
      no_mem_exit("init_img: img->intra_block");
    for (i=0; i<j; i++)
    {
      if ((img->intra_block[i] = (int*)calloc(4, sizeof(int))) == NULL)
        no_mem_exit ("init_img: img->intra_block");
    }
  }

  // allocate memory for intra pred mode buffer for each block: img->ipredmode
  // int  img->ipredmode[90][74];
  size_x=img->width/BLOCK_SIZE+3;
  size_y=img->height/BLOCK_SIZE+3;
  get_mem2Dint(&(img->ipredmode), img->width/BLOCK_SIZE+3, img->height/BLOCK_SIZE+3);        //need two extra rows at right and bottom
  if(input->InterlaceCodingOption >= MB_CODING) 
  {
    get_mem2Dint(&(img->ipredmode_top), img->width/BLOCK_SIZE+3, img->height/BLOCK_SIZE+3);
    get_mem2Dint(&(img->ipredmode_bot), img->width/BLOCK_SIZE+3, img->height/BLOCK_SIZE+3);
    get_mem2Dint(&(rddata_top_frame_mb.ipredmode), img->width/BLOCK_SIZE+3, img->height/BLOCK_SIZE+3);
    get_mem2Dint(&(rddata_bot_frame_mb.ipredmode), img->width/BLOCK_SIZE+3, img->height/BLOCK_SIZE+3);
    get_mem2Dint(&(rddata_top_field_mb.ipredmode), img->width/BLOCK_SIZE+3, img->height/BLOCK_SIZE+3);
    get_mem2Dint(&(rddata_bot_field_mb.ipredmode), img->width/BLOCK_SIZE+3, img->height/BLOCK_SIZE+3);
  }
  // CAVLC mem
  if((img->nz_coeff = (int****)calloc(img->width/MB_BLOCK_SIZE,sizeof(int***))) == NULL)
    no_mem_exit("get_mem4global_buffers: nzcoeff");
  for(j=0;j<img->width/MB_BLOCK_SIZE;j++)
  {
    get_mem3Dint(&(img->nz_coeff[j]), img->height/MB_BLOCK_SIZE, 4, 6);
  }

  CAVLC_init();

  // Prediction mode is set to -1 outside the frame, indicating that no prediction can be made from this part
  for (i=0; i < img->width/BLOCK_SIZE+1; i++)
  {
    // img->ipredmode[i][0]=-1;
    // img->ipredmode[i][size_y-2]=-1;
    // img->ipredmode[i][size_y-1]=-1;
    img->ipredmode[i+1][0]=-1;
    img->ipredmode[i+1][img->height/BLOCK_SIZE+1]=-1;
  }
  for (j=0; j < img->height/BLOCK_SIZE+1; j++)
  {
    // img->ipredmode[0][j]=-1;
    // img->ipredmode[size_x-2][j]=-1;
    // img->ipredmode[size_x-1][j]=-1;
    img->ipredmode[0][j+1]=-1;
    img->ipredmode[img->width/BLOCK_SIZE+1][j+1]=-1;
  }

  if(input->InterlaceCodingOption >= MB_CODING)
  {
    for (i=0; i < img->width/BLOCK_SIZE+1; i++)
    {
    img->ipredmode_top[i+1][0]=-1;
    img->ipredmode_top[i+1][img->height/BLOCK_SIZE+1]=-1;
    img->ipredmode_bot[i+1][0]=-1;
    img->ipredmode_bot[i+1][img->height/BLOCK_SIZE+1]=-1;

    }
    for (j=0; j < img->height/BLOCK_SIZE+1; j++)
    {
    img->ipredmode_top[0][j+1]=-1;
    img->ipredmode_top[img->width/BLOCK_SIZE+1][j+1]=-1;
    img->ipredmode_bot[0][j+1]=-1;
    img->ipredmode_bot[img->width/BLOCK_SIZE+1][j+1]=-1;

    }
  }

  img->mb_y_upd=0;

  RandomIntraInit (img->width/16, img->height/16, input->RandomIntraMBRefresh);
  FmoInit (img->width/16, img->height/16, input->num_slice_groups_minus1, input->FmoType, NULL); 

  InitSEIMessages();  // Tian Dong (Sept 2002)

  // Initialize filtering parameters. If sending parameters, the offsets are 
  // multiplied by 2 since inputs are taken in "div 2" format.
  // If not sending paramters, all fields are cleared 
  if (input->LFSendParameters)
  {
    input->LFAlphaC0Offset <<= 1;
    input->LFBetaOffset <<= 1;
  }
  else
  {
    input->LFDisableIdc = 0;
    input->LFAlphaC0Offset = 0;
    input->LFBetaOffset = 0;
  }


}

/*!
 ***********************************************************************
 * \brief
 *    Free the Image structures
 * \par Input:
 *    Image Parameters struct img_par *img
 ***********************************************************************
 */
void free_img ()
{
  CloseSEIMessages(); // Tian Dong (Sept 2002)
  free_mem_mv (img->mv);
  free_mem_mv (img->p_fwMV);
  free_mem_mv (img->p_bwMV);
  free_mem_mv (img->all_mv);
  free_mem_mv (img->all_bmv);
  free_mem_mv (img->abp_all_dmv);
  if(input->InterlaceCodingOption >= MB_CODING) 
  {
    free_mem_mv (img->abp_all_dmv_top);
    free_mem_mv (img->abp_all_dmv_bot);
  }

  free_mem_ACcoeff (img->cofAC);
  free_mem_DCcoeff (img->cofDC);

  free (img->quad-255);
}



/*!
 ************************************************************************
 * \brief
 *    Allocates the picture structure along with its dependent
 *    data structures
 * \return
 *    Pointer to a Picture
 ************************************************************************
 */

Picture *malloc_picture()
{
  Picture *pic;
  if ((pic = calloc (1, sizeof (Picture))) == NULL) no_mem_exit ("malloc_picture: Picture structure");
  //! Note: slice structures are allocated as needed in code_a_picture
  return pic;
}

/*!
 ************************************************************************
 * \brief
 *    Frees a picture
 * \param
 *    pic: POinter to a Picture to be freed
 ************************************************************************
 */


void free_picture(Picture *pic)
{
  int i;

  if (pic != NULL)
  {
    for (i=0; i<pic->no_slices; i++)
      if (pic->slices[i] != NULL)
        free (pic->slices[i]);
    free (pic);
  }
}




/*!
 ************************************************************************
 * \brief
 *    Reports the gathered information to appropriate outputs
 * \par Input:
 *    struct inp_par *inp,                                            \n
 *    struct img_par *img,                                            \n
 *    struct stat_par *stat,                                          \n
 *    struct stat_par *stat                                           \n
 *
 * \par Output:
 *    None
 ************************************************************************
 */
void report()
{
  int bit_use[2][2] ;
  int i,j;
  char name[20];
  int bit_use_Bframe=0;
  int total_bits;
  float frame_rate;
  float mean_motion_info_bit_use[2];

#ifndef WIN32
  time_t now;
  struct tm *l_time;
  char string[1000];
#else
  char timebuf[128];
#endif
  bit_use[0][0]=1;
  bit_use[1][0]=max(1,input->no_frames-1);

  //  Accumulate bit usage for inter and intra frames
  bit_use[0][1]=bit_use[1][1]=0;


  for (i=0; i < 11; i++)
    bit_use[1][1] += stat->bit_use_mode_inter[0][i];


  for (j=0;j<2;j++)
  {
    bit_use[j][1]+=stat->bit_use_header[j];
    bit_use[j][1]+=stat->bit_use_mb_type[j];
    bit_use[j][1]+=stat->tmp_bit_use_cbp[j];
    bit_use[j][1]+=stat->bit_use_coeffY[j];
    bit_use[j][1]+=stat->bit_use_coeffC[j];
    bit_use[j][1]+=stat->bit_use_delta_quant[j];
    bit_use[j][1]+=stat->bit_use_stuffingBits[j];
  }

  // B pictures
  if(Bframe_ctr!=0)
  {
    bit_use_Bframe=0;
    for(i=0; i<11; i++)
      bit_use_Bframe += stat->bit_use_mode_inter[1][i]; 
    bit_use_Bframe += stat->bit_use_header[2];
    bit_use_Bframe += stat->bit_use_mb_type[2];
    bit_use_Bframe += stat->tmp_bit_use_cbp[2];
    bit_use_Bframe += stat->bit_use_coeffY[2];
    bit_use_Bframe += stat->bit_use_coeffC[2];
    bit_use_Bframe += stat->bit_use_delta_quant[2];
    bit_use_Bframe +=stat->bit_use_stuffingBits[2];

    stat->bitrate_P=(stat->bit_ctr_0+stat->bit_ctr_P)*(float)(img->framerate/(input->jumpd+1))/input->no_frames;
#ifdef _ADAPT_LAST_GROUP_
    stat->bitrate_B=(stat->bit_ctr_B)*(float)(img->framerate/(input->jumpd+1))*initial_Bframes/Bframe_ctr;
#else
    stat->bitrate_B=(stat->bit_ctr_B)*(float)(img->framerate/(input->jumpd+1))*input->successive_Bframe/Bframe_ctr;
#endif
  }
  else
  {
    if (input->no_frames > 1)
    {
      stat->bitrate=(bit_use[0][1]+bit_use[1][1])*(float)img->framerate/(input->no_frames*(input->jumpd+1));
    }
  }

  fprintf(stdout,"--------------------------------------------------------------------------\n");
  fprintf(stdout,   " Freq. for encoded bitstream       : %1.0f\n",(float)img->framerate/(float)(input->jumpd+1));
  if(input->hadamard)
    fprintf(stdout," Hadamard transform                : Used\n");
  else
    fprintf(stdout," Hadamard transform                : Not used\n");

  fprintf(stdout," Image format                      : %dx%d\n",input->img_width,input->img_height);

  if(input->intra_upd)
    fprintf(stdout," Error robustness                  : On\n");
  else
    fprintf(stdout," Error robustness                  : Off\n");
  fprintf(stdout,    " Search range                      : %d\n",input->search_range);


#ifdef _ADDITIONAL_REFERENCE_FRAME_
  if (input->add_ref_frame >= input->no_multpred)
  {
      fprintf(stdout,   " No of ref. frames used in P pred  : %d (+ no. %d)\n",input->no_multpred,input->add_ref_frame);
      if(input->successive_Bframe != 0)
        fprintf(stdout,   " No of ref. frames used in B pred  : %d (+ no. %d)\n",input->no_multpred,input->add_ref_frame);
  }
  else
#endif
  {
    fprintf(stdout,   " No of ref. frames used in P pred  : %d\n",input->no_multpred);
    if(input->successive_Bframe != 0)
      fprintf(stdout,   " No of ref. frames used in B pred  : %d\n",input->no_multpred);
  }
  fprintf(stdout,   " Total encoding time for the seq.  : %.3f sec \n",tot_time*0.001);

  // B pictures
  fprintf(stdout, " Sequence type                     :" );
  if(input->successive_Bframe==0 && input->StoredBPictures > 0)
    fprintf(stdout, " I-P-BS-BS (QP: I %d, P BS %d) \n",input->qp0, input->qpN);
  else if(input->successive_Bframe==1 && input->StoredBPictures > 0)
    fprintf(stdout, " I-B-P-BS-B-BS (QP: I %d, P BS %d, B %d) \n",input->qp0, input->qpN, input->qpB);
  else if(input->successive_Bframe==2 && input->StoredBPictures > 0)
    fprintf(stdout, " I-B-B-P-B-B-BS (QP: I %d, P BS %d, B %d) \n",input->qp0, input->qpN, input->qpB);
  else if(input->successive_Bframe==1)   fprintf(stdout, " IBPBP (QP: I %d, P %d, B %d) \n",
    input->qp0, input->qpN, input->qpB);
  else if(input->successive_Bframe==2) fprintf(stdout, " IBBPBBP (QP: I %d, P %d, B %d) \n",
    input->qp0, input->qpN, input->qpB);
  else if(input->successive_Bframe==0 && input->sp_periodicity==0) fprintf(stdout, " IPPP (QP: I %d, P %d) \n",   input->qp0, input->qpN);

  else fprintf(stdout, " I-P-P-SP-P (QP: I %d, P %d, SP (%d, %d)) \n",  input->qp0, input->qpN, input->qpsp, input->qpsp_pred);

  // report on entropy coding  method
  if (input->symbol_mode == UVLC)
    fprintf(stdout," Entropy coding method             : CAVLC\n");
  else
    fprintf(stdout," Entropy coding method             : CABAC\n");

#ifdef _FULL_SEARCH_RANGE_
  if (input->full_search == 2)
    fprintf(stdout," Search range restrictions         : none\n");
  else if (input->full_search == 1)
    fprintf(stdout," Search range restrictions         : older reference frames\n");
  else
    fprintf(stdout," Search range restrictions         : smaller blocks and older reference frames\n");
#endif

  if (input->rdopt)
    fprintf(stdout," RD-optimized mode decision        : used\n");
  else
    fprintf(stdout," RD-optimized mode decision        : not used\n");

  switch(input->partition_mode)
    {
    case PAR_DP_1:
      fprintf(stdout," Data Partitioning Mode            : 1 partition \n");
      break;
    case PAR_DP_3:
      fprintf(stdout," Data Partitioning Mode            : 3 partitions \n");
      break;
    default:
      fprintf(stdout," Data Partitioning Mode            : not supported\n");
      break;
    }


    switch(input->of_mode)
    {
    case PAR_OF_ANNEXB:
      fprintf(stdout," Output File Format                : H.26L Bit Stream File Format \n");
      break;
    case PAR_OF_RTP:
      fprintf(stdout," Output File Format                : RTP Packet File Format \n");
      break;
    default:
      fprintf(stdout," Output File Format                : not supported\n");
      break;
    }




  fprintf(stdout,"------------------ Average data all frames  ------------------------------\n");
  fprintf(stdout," SNR Y(dB)                         : %5.2f\n",snr->snr_ya);
  fprintf(stdout," SNR U(dB)                         : %5.2f\n",snr->snr_ua);
  fprintf(stdout," SNR V(dB)                         : %5.2f\n",snr->snr_va);

  if(Bframe_ctr!=0)
  {

    if(input->StoredBPictures > 0)
    {
      fprintf(stdout, " Total bits                        : %d (I %5d, P BS %5d, B %d) \n",
            total_bits=stat->bit_ctr_P + stat->bit_ctr_0 + stat->bit_ctr_B, stat->bit_ctr_0, stat->bit_ctr_P, stat->bit_ctr_B);
    }
    else
    {
      fprintf(stdout, " Total bits                        : %d (I %5d, P %5d, B %d) \n",
            total_bits=stat->bit_ctr_P + stat->bit_ctr_0 + stat->bit_ctr_B, stat->bit_ctr_0, stat->bit_ctr_P, stat->bit_ctr_B);
    }

    frame_rate = (float)(img->framerate *(input->successive_Bframe + 1)) / (float) (input->jumpd+1);
    stat->bitrate= ((float) total_bits * frame_rate)/((float) (input->no_frames + Bframe_ctr));

    fprintf(stdout, " Bit rate (kbit/s)  @ %2.2f Hz     : %5.2f\n", frame_rate, stat->bitrate/1000);

  }
  else if(input->StoredBPictures > 0)
  {
    fprintf(stdout, " Total bits                        : %d (I %5d, P BS %5d) \n",
    total_bits=stat->bit_ctr_P + stat->bit_ctr_0 , stat->bit_ctr_0, stat->bit_ctr_P);

    frame_rate = (float)img->framerate / ( (float) (input->jumpd + 1) );
    stat->bitrate= ((float) total_bits * frame_rate)/((float) input->no_frames );

    fprintf(stdout, " Bit rate (kbit/s)  @ %2.2f Hz     : %5.2f\n", frame_rate, stat->bitrate/1000);
  }
  else if (input->sp_periodicity==0)
  {
    fprintf(stdout, " Total bits                        : %d (I %5d, P %5d) \n",
    total_bits=stat->bit_ctr_P + stat->bit_ctr_0 , stat->bit_ctr_0, stat->bit_ctr_P);

    frame_rate = (float)img->framerate / ( (float) (input->jumpd + 1) );
    stat->bitrate= ((float) total_bits * frame_rate)/((float) input->no_frames );

    fprintf(stdout, " Bit rate (kbit/s)  @ %2.2f Hz     : %5.2f\n", frame_rate, stat->bitrate/1000);
  }else
  {
    fprintf(stdout, " Total bits                        : %d (I %5d, P %5d) \n",
    total_bits=stat->bit_ctr_P + stat->bit_ctr_0 , stat->bit_ctr_0, stat->bit_ctr_P);

    frame_rate = (float)img->framerate / ( (float) (input->jumpd + 1) );
    stat->bitrate= ((float) total_bits * frame_rate)/((float) input->no_frames );

    fprintf(stdout, " Bit rate (kbit/s)  @ %2.2f Hz     : %5.2f\n", frame_rate, stat->bitrate/1000);
  }

  fprintf(stdout, " Bits to avoid Startcode Emulation : %d \n", stat->bit_ctr_emulationprevention);
  fprintf(stdout, " Bits for parameter sets           : %d \n", stat->bit_ctr_parametersets);

  fprintf(stdout,"--------------------------------------------------------------------------\n");
  fprintf(stdout,"Exit JM %s encoder ver %s ", JM, VERSION);
  fprintf(stdout,"\n");

  // status file
  if ((p_stat=fopen("stat.dat","wt"))==0)
  {
    snprintf(errortext, ET_SIZE, "Error open file %s", "stat.dat");
    error(errortext, 500);
  }
  fprintf(p_stat," -------------------------------------------------------------- \n");
  fprintf(p_stat,"  This file contains statistics for the last encoded sequence   \n");
  fprintf(p_stat," -------------------------------------------------------------- \n");
  fprintf(p_stat,   " Sequence                     : %s\n",input->infile);
  fprintf(p_stat,   " No.of coded pictures         : %4d\n",input->no_frames+Bframe_ctr);
  fprintf(p_stat,   " Freq. for encoded bitstream  : %4.0f\n",frame_rate);

  // B pictures
  if(input->successive_Bframe != 0)
  {
    fprintf(p_stat,   " BaseLayer Bitrate(kb/s)      : %6.2f\n", stat->bitrate_P/1000);
    fprintf(p_stat,   " EnhancedLyaer Bitrate(kb/s)  : %6.2f\n", stat->bitrate_B/1000);
  }
  else
    fprintf(p_stat,   " Bitrate(kb/s)                : %6.2f\n", stat->bitrate/1000);

  if(input->hadamard)
    fprintf(p_stat," Hadamard transform           : Used\n");
  else
    fprintf(p_stat," Hadamard transform           : Not used\n");

  fprintf(p_stat,  " Image format                 : %dx%d\n",input->img_width,input->img_height);

  if(input->intra_upd)
    fprintf(p_stat," Error robustness             : On\n");
  else
    fprintf(p_stat," Error robustness             : Off\n");

  fprintf(p_stat,  " Search range                 : %d\n",input->search_range);

#ifdef _ADDITIONAL_REFERENCE_FRAME_
  if (input->add_ref_frame >= input->no_multpred)
    {
      fprintf(p_stat,   " No of frame used in P pred   : %d (+ no. %d)\n",input->no_multpred,input->add_ref_frame);
      if(input->successive_Bframe != 0)
        fprintf(p_stat, " No of frame used in B pred   : %d (+ no. %d)\n",input->no_multpred,input->add_ref_frame);
    }
  else
#endif
  {
    fprintf(p_stat,   " No of frame used in P pred   : %d\n",input->no_multpred);
    if(input->successive_Bframe != 0)
      fprintf(p_stat, " No of frame used in B pred   : %d\n",input->no_multpred);
  }
  if (input->symbol_mode == UVLC)
    fprintf(p_stat,   " Entropy coding method        : CAVLC\n");
  else
    fprintf(p_stat,   " Entropy coding method        : CABAC\n");

  if(input->InterlaceCodingOption >= MB_CODING)
    fprintf(p_stat, " MB Field Coding : On \n");

#ifdef _FULL_SEARCH_RANGE_
  if (input->full_search == 2)
    fprintf(p_stat," Search range restrictions    : none\n");
  else if (input->full_search == 1)
    fprintf(p_stat," Search range restrictions    : older reference frames\n");
  else
    fprintf(p_stat," Search range restrictions    : smaller blocks and older reference frames\n");
#endif
  if (input->rdopt)
    fprintf(p_stat," RD-optimized mode decision   : used\n");
  else
    fprintf(p_stat," RD-optimized mode decision   : not used\n");

  fprintf(p_stat," -------------------|---------------|---------------|\n");
  fprintf(p_stat,"     Item           |     Intra     |   All frames  |\n");
  fprintf(p_stat," -------------------|---------------|---------------|\n");
  fprintf(p_stat," SNR Y(dB)          |");
  fprintf(p_stat," %5.2f         |",snr->snr_y1);
  fprintf(p_stat," %5.2f         |\n",snr->snr_ya);
  fprintf(p_stat," SNR U/V (dB)       |");
  fprintf(p_stat," %5.2f/%5.2f   |",snr->snr_u1,snr->snr_v1);
  fprintf(p_stat," %5.2f/%5.2f   |\n",snr->snr_ua,snr->snr_va);

  // QUANT.
  fprintf(p_stat," Average quant      |");
  fprintf(p_stat," %5d         |",absm(input->qp0));
  fprintf(p_stat," %5.2f         |\n",(float)stat->quant1/max(1.0,(float)stat->quant0));

  // MODE
  fprintf(p_stat,"\n -------------------|---------------|\n");
  fprintf(p_stat,"   Intra            |   Mode used   |\n");
  fprintf(p_stat," -------------------|---------------|\n");

  fprintf(p_stat," Mode 0  intra old  | %5d         |\n",stat->mode_use_intra[I4MB]);
  fprintf(p_stat," Mode 1+ intra new  | %5d         |\n",stat->mode_use_intra[I16MB]);

  fprintf(p_stat,"\n -------------------|---------------|-----------------|\n");
  fprintf(p_stat,"   Inter            |   Mode used   | MotionInfo bits |\n");
  fprintf(p_stat," -------------------|---------------|-----------------|");
  fprintf(p_stat,"\n Mode  0  (copy)    | %5d         |    %8.2f     |",stat->mode_use_inter[0][0   ],(float)stat->bit_use_mode_inter[0][0   ]/(float)bit_use[1][0]);
  fprintf(p_stat,"\n Mode  1  (16x16)   | %5d         |    %8.2f     |",stat->mode_use_inter[0][1   ],(float)stat->bit_use_mode_inter[0][1   ]/(float)bit_use[1][0]);
  fprintf(p_stat,"\n Mode  2  (16x8)    | %5d         |    %8.2f     |",stat->mode_use_inter[0][2   ],(float)stat->bit_use_mode_inter[0][2   ]/(float)bit_use[1][0]);
  fprintf(p_stat,"\n Mode  3  (8x16)    | %5d         |    %8.2f     |",stat->mode_use_inter[0][3   ],(float)stat->bit_use_mode_inter[0][3   ]/(float)bit_use[1][0]);
  fprintf(p_stat,"\n Mode  4  (8x8)     | %5d         |    %8.2f     |",stat->mode_use_inter[0][P8x8],(float)stat->bit_use_mode_inter[0][P8x8]/(float)bit_use[1][0]);
  fprintf(p_stat,"\n Mode  5  intra old | %5d         |-----------------|",stat->mode_use_inter[0][I4MB]);
  fprintf(p_stat,"\n Mode  6+ intr.new  | %5d         |",stat->mode_use_inter[0][I16MB]);
  mean_motion_info_bit_use[0] = (float)(stat->bit_use_mode_inter[0][0] + stat->bit_use_mode_inter[0][1] + stat->bit_use_mode_inter[0][2] 
                                      + stat->bit_use_mode_inter[0][3] + stat->bit_use_mode_inter[0][P8x8])/(float) bit_use[1][0]; 

  // B pictures
  if(input->successive_Bframe!=0 && Bframe_ctr!=0)
  {
 
    fprintf(p_stat,"\n\n -------------------|---------------|-----------------|\n");
    fprintf(p_stat,"   B frame          |   Mode used   | MotionInfo bits |\n");
    fprintf(p_stat," -------------------|---------------|-----------------|");
    fprintf(p_stat,"\n Mode  0  (copy)    | %5d         |    %8.2f     |",stat->mode_use_inter[1][0   ],(float)stat->bit_use_mode_inter[1][0   ]/(float)Bframe_ctr);
    fprintf(p_stat,"\n Mode  1  (16x16)   | %5d         |    %8.2f     |",stat->mode_use_inter[1][1   ],(float)stat->bit_use_mode_inter[1][1   ]/(float)Bframe_ctr);
    fprintf(p_stat,"\n Mode  2  (16x8)    | %5d         |    %8.2f     |",stat->mode_use_inter[1][2   ],(float)stat->bit_use_mode_inter[1][2   ]/(float)Bframe_ctr);
    fprintf(p_stat,"\n Mode  3  (8x16)    | %5d         |    %8.2f     |",stat->mode_use_inter[1][3   ],(float)stat->bit_use_mode_inter[1][3   ]/(float)Bframe_ctr);
    fprintf(p_stat,"\n Mode  4  (8x8)     | %5d         |    %8.2f     |",stat->mode_use_inter[1][P8x8],(float)stat->bit_use_mode_inter[1][P8x8]/(float)Bframe_ctr);
    fprintf(p_stat,"\n Mode  5  intra old | %5d         |-----------------|",stat->mode_use_inter[1][I4MB]);
    fprintf(p_stat,"\n Mode  6+ intr.new  | %5d         |",stat->mode_use_inter[1][I16MB]);
    mean_motion_info_bit_use[1] = (float)(stat->bit_use_mode_inter[1][0] + stat->bit_use_mode_inter[1][1] + stat->bit_use_mode_inter[1][2] 
                                      + stat->bit_use_mode_inter[1][3] + stat->bit_use_mode_inter[1][P8x8])/(float) Bframe_ctr; 

  }

  fprintf(p_stat,"\n\n --------------------|----------------|----------------|----------------|\n");
  fprintf(p_stat,"  Bit usage:         |      Intra     |      Inter     |    B frame     |\n");
  fprintf(p_stat," --------------------|----------------|----------------|----------------|\n");

  fprintf(p_stat," Header              |");
  fprintf(p_stat," %10.2f     |",(float) stat->bit_use_header[0]/bit_use[0][0]);
  fprintf(p_stat," %10.2f     |",(float) stat->bit_use_header[1]/bit_use[1][0]);
  if(input->successive_Bframe!=0 && Bframe_ctr!=0)
    fprintf(p_stat," %10.2f     |",(float) stat->bit_use_header[2]/Bframe_ctr);
  else fprintf(p_stat," %10.2f     |", 0.);
  fprintf(p_stat,"\n");

  fprintf(p_stat," Mode                |");
  fprintf(p_stat," %10.2f     |",(float)stat->bit_use_mb_type[0]/bit_use[0][0]);
  fprintf(p_stat," %10.2f     |",(float)stat->bit_use_mb_type[1]/bit_use[1][0]);
  if(input->successive_Bframe!=0 && Bframe_ctr!=0)
    fprintf(p_stat," %10.2f     |",(float)stat->bit_use_mb_type[2]/Bframe_ctr);
  else fprintf(p_stat," %10.2f     |", 0.);
  fprintf(p_stat,"\n");

  fprintf(p_stat," Motion Info         |");
  fprintf(p_stat,"        ./.     |");
  fprintf(p_stat," %10.2f     |",mean_motion_info_bit_use[0]);
  if(input->successive_Bframe!=0 && Bframe_ctr!=0)
    fprintf(p_stat," %10.2f     |",mean_motion_info_bit_use[1]);
  else fprintf(p_stat," %10.2f     |", 0.);
  fprintf(p_stat,"\n");

  fprintf(p_stat," CBP Y/C             |");
  for (j=0; j < 2; j++)
  {
    fprintf(p_stat," %10.2f     |", (float)stat->tmp_bit_use_cbp[j]/bit_use[j][0]);
  }
  if(input->successive_Bframe!=0 && Bframe_ctr!=0)
    fprintf(p_stat," %10.2f     |", (float)stat->tmp_bit_use_cbp[2]/Bframe_ctr);
  else fprintf(p_stat," %10.2f     |", 0.);
  fprintf(p_stat,"\n");

  if(input->successive_Bframe!=0 && Bframe_ctr!=0)
    fprintf(p_stat," Coeffs. Y           | %10.2f     | %10.2f     | %10.2f     |\n",
      (float)stat->bit_use_coeffY[0]/bit_use[0][0], (float)stat->bit_use_coeffY[1]/bit_use[1][0], (float)stat->bit_use_coeffY[2]/Bframe_ctr);
  else
    fprintf(p_stat," Coeffs. Y           | %10.2f     | %10.2f     | %10.2f     |\n",
      (float)stat->bit_use_coeffY[0]/bit_use[0][0], (float)stat->bit_use_coeffY[1]/(float)bit_use[1][0], 0.);

  if(input->successive_Bframe!=0 && Bframe_ctr!=0)
    fprintf(p_stat," Coeffs. C           | %10.2f     | %10.2f     | %10.2f     |\n",
      (float)stat->bit_use_coeffC[0]/bit_use[0][0], (float)stat->bit_use_coeffC[1]/bit_use[1][0], (float)stat->bit_use_coeffC[2]/Bframe_ctr);
  else
    fprintf(p_stat," Coeffs. C           | %10.2f     | %10.2f     | %10.2f     |\n",
      (float)stat->bit_use_coeffC[0]/bit_use[0][0], (float)stat->bit_use_coeffC[1]/bit_use[1][0], 0.);

  if(input->successive_Bframe!=0 && Bframe_ctr!=0)
    fprintf(p_stat," Delta quant         | %10.2f     | %10.2f     | %10.2f     |\n",
      (float)stat->bit_use_delta_quant[0]/bit_use[0][0], (float)stat->bit_use_delta_quant[1]/bit_use[1][0], (float)stat->bit_use_delta_quant[2]/Bframe_ctr);
  else
    fprintf(p_stat," Delta quant         | %10.2f     | %10.2f     | %10.2f     |\n",
      (float)stat->bit_use_delta_quant[0]/bit_use[0][0], (float)stat->bit_use_delta_quant[1]/bit_use[1][0], 0.);

  if(input->successive_Bframe!=0 && Bframe_ctr!=0)
    fprintf(p_stat," Stuffing Bits       | %10.2f     | %10.2f     | %10.2f     |\n",
      (float)stat->bit_use_stuffingBits[0]/bit_use[0][0], (float)stat->bit_use_stuffingBits[1]/bit_use[1][0], (float)stat->bit_use_stuffingBits[2]/Bframe_ctr);
  else
    fprintf(p_stat," Stuffing Bits       | %10.2f     | %10.2f     | %10.2f     |\n",
      (float)stat->bit_use_stuffingBits[0]/bit_use[0][0], (float)stat->bit_use_stuffingBits[1]/bit_use[1][0], 0.);



  fprintf(p_stat," --------------------|----------------|----------------|----------------|\n");

  fprintf(p_stat," average bits/frame  |");
  for (i=0; i < 2; i++)
  {
    fprintf(p_stat," %10.2f     |", (float) bit_use[i][1]/(float) bit_use[i][0] );
  }
  if(input->successive_Bframe!=0 && Bframe_ctr!=0)
    fprintf(p_stat," %10.2f     |", (float) bit_use_Bframe/ (float) Bframe_ctr );
  else fprintf(p_stat," %10.2f     |", 0.);

  fprintf(p_stat,"\n");
  fprintf(p_stat," --------------------|----------------|----------------|----------------|\n");

  fclose(p_stat);

  // write to log file
  if ((p_log=fopen("log.dat","r"))==0)                      // check if file exist
  {
    if ((p_log=fopen("log.dat","a"))==NULL)            // append new statistic at the end
    {
      snprintf(errortext, ET_SIZE, "Error open file %s  \n","log.dat");
      error(errortext, 500);
    }
    else                                            // Create header for new log file
    {
      fprintf(p_log," ----------------------------------------------------------------------------------------------------------------------------------------------------------------- \n");
      fprintf(p_log,"|            Encoder statistics. This file is generated during first encoding session, new sessions will be appended                                              |\n");
      fprintf(p_log," ----------------------------------------------------------------------------------------------------------------------------------------------------------------- \n");
      fprintf(p_log,"| Date  | Time  |    Sequence        |#Img|Quant1|QuantN|Format|Hadamard|Search r|#Ref | Freq |Intra upd|SNRY 1|SNRU 1|SNRV 1|SNRY N|SNRU N|SNRV N|#Bitr P|#Bitr B|\n");
      fprintf(p_log," ----------------------------------------------------------------------------------------------------------------------------------------------------------------- \n");
    }
  }
  else
  {
    fclose (p_log);
    if ((p_log=fopen("log.dat","a"))==NULL)            // File exist,just open for appending
    {
      snprintf(errortext, ET_SIZE, "Error open file %s  \n","log.dat");
      error(errortext, 500);
    }
  }

#ifdef WIN32
  _strdate( timebuf );
  fprintf(p_log,"| %1.5s |",timebuf );

  _strtime( timebuf);
  fprintf(p_log," % 1.5s |",timebuf);
#else
  now = time ((time_t *) NULL); // Get the system time and put it into 'now' as 'calender time'
  time (&now);
  l_time = localtime (&now);
  strftime (string, sizeof string, "%d-%b-%Y", l_time);
  fprintf(p_log,"| %1.5s |",string );

  strftime (string, sizeof string, "%H:%M:%S", l_time);
  fprintf(p_log," %1.5s |",string );
#endif

  for (i=0;i<20;i++)
    name[i]=input->infile[i+max(0,strlen(input->infile)-20)]; // write last part of path, max 20 chars
  fprintf(p_log,"%20.20s|",name);

  fprintf(p_log,"%3d |",input->no_frames);
  fprintf(p_log,"  %2d  |",input->qp0);
  fprintf(p_log,"  %2d  |",input->qpN);


  fprintf(p_log,"%dx%d|",input->img_width,input->img_height);


  if (input->hadamard==1)
    fprintf(p_log,"   ON   |");
  else
    fprintf(p_log,"   OFF  |");

  fprintf(p_log,"   %2d   |",input->search_range );

  fprintf(p_log," %2d  |",input->no_multpred);

  fprintf(p_log," %2d  |",img->framerate/(input->jumpd+1));

  if (input->intra_upd==1)
    fprintf(p_log,"   ON    |");
  else
    fprintf(p_log,"   OFF   |");

  fprintf(p_log,"%5.3f|",snr->snr_y1);
  fprintf(p_log,"%5.3f|",snr->snr_u1);
  fprintf(p_log,"%5.3f|",snr->snr_v1);
  fprintf(p_log,"%5.3f|",snr->snr_ya);
  fprintf(p_log,"%5.3f|",snr->snr_ua);
  fprintf(p_log,"%5.3f|",snr->snr_va);
  if(input->successive_Bframe != 0)
  {
    fprintf(p_log,"%7.0f|",stat->bitrate_P);
    fprintf(p_log,"%7.0f|\n",stat->bitrate_B);
  }
  else
  {
    fprintf(p_log,"%7.0f|",stat->bitrate);
    fprintf(p_log,"%7.0f|\n",0.0);
  }

  fclose(p_log);

  p_log=fopen("data.txt","a");

  if(input->successive_Bframe != 0 && Bframe_ctr != 0) // B picture used
  {
    fprintf(p_log, "%3d %2d %2d %2.2f %2.2f %2.2f %5d "
          "%2.2f %2.2f %2.2f %5d "
        "%2.2f %2.2f %2.2f %5d %5d %.3f\n",
        input->no_frames, input->qp0, input->qpN,
        snr->snr_y1,
        snr->snr_u1,
        snr->snr_v1,
        stat->bit_ctr_0,
        0.0,
        0.0,
        0.0,
        0,
        snr->snr_ya,
        snr->snr_ua,
        snr->snr_va,
        (stat->bit_ctr_0+stat->bit_ctr)/(input->no_frames+Bframe_ctr),
        stat->bit_ctr_B/Bframe_ctr,
        (double)0.001*tot_time/(input->no_frames+Bframe_ctr));
  }
  else
  {
    fprintf(p_log, "%3d %2d %2d %2.2f %2.2f %2.2f %5d "
          "%2.2f %2.2f %2.2f %5d "
        "%2.2f %2.2f %2.2f %5d %5d %.3f\n",
        input->no_frames, input->qp0, input->qpN,
        snr->snr_y1,
        snr->snr_u1,
        snr->snr_v1,
        stat->bit_ctr_0,
        0.0,
        0.0,
        0.0,
        0,
        snr->snr_ya,
        snr->snr_ua,
        snr->snr_va,
        (stat->bit_ctr_0+stat->bit_ctr)/input->no_frames,
        0,
        (double)0.001*tot_time/input->no_frames);
  }

  fclose(p_log);

 
//  free(stat->mode_use_Bframe);
//  free(stat->bit_use_mode_Bframe);
}


/*!
 ************************************************************************
 * \brief
 *    Prints the header of the protocol.
 * \par Input:
 *    struct inp_par *inp
 * \par Output:
 *    none
 ************************************************************************
 */
void information_init()
{
  printf("--------------------------------------------------------------------------\n");
  printf(" Input YUV file                    : %s \n",input->infile);
  printf(" Output H.26L bitstream            : %s \n",input->outfile);
  if (p_dec != NULL)
    printf(" Output YUV file                   : %s \n",input->ReconFile);
  printf(" Output log file                   : log.dat \n");
  printf(" Output statistics file            : stat.dat \n");
  printf("--------------------------------------------------------------------------\n");
  printf(" Frame   Bit/pic   QP   SnrY    SnrU    SnrV    Time(ms) Frm/Fld IntraMBs\n");
}


/*!
 ************************************************************************
 * \brief
 *    Dynamic memory allocation of frame size related global buffers
 *    buffers are defined in global.h, allocated memory must be freed in
 *    void free_global_buffers()
 * \par Input:
 *    Input Parameters struct inp_par *inp,                            \n
 *    Image Parameters struct img_par *img
 * \return Number of allocated bytes
 ************************************************************************
 */
int init_global_buffers()
{
  int j,memory_size=0;
  int height_field = img->height/2;
#ifdef _ADAPT_LAST_GROUP_
  extern int *last_P_no_frm;
  extern int *last_P_no_fld;

  if ((last_P_no_frm = (int*)malloc(2*img->buf_cycle*sizeof(int))) == NULL)
    no_mem_exit("init_global_buffers: last_P_no");
  if(input->InterlaceCodingOption != FRAME_CODING)
    if ((last_P_no_fld = (int*)malloc(4*img->buf_cycle*sizeof(int))) == NULL)
      no_mem_exit("init_global_buffers: last_P_no");
#endif

  if(input->InterlaceCodingOption >= MB_CODING)
    memory_size += get_mem2Dint(&field_mb, img->height/MB_BLOCK_SIZE, img->width/MB_BLOCK_SIZE);


  // allocate memory for encoding frame buffers: imgY, imgUV
  // byte imgY[288][352];
  // byte imgUV[2][144][176];
  memory_size += get_mem2D(&imgY_frm, img->height, img->width);
  memory_size += get_mem3D(&imgUV_frm, 2, img->height_cr, img->width_cr);

  // allocate memory for reference frame buffers: imgY_org, imgUV_org
  // byte imgY_org[288][352];
  // byte imgUV_org[2][144][176];
  memory_size += get_mem2D(&imgY_org_frm, img->height, img->width);
  memory_size += get_mem3D(&imgUV_org_frm, 2, img->height_cr, img->width_cr);

  // allocate memory for temp P and B-frame motion vector buffer: tmp_mv, temp_mv_block
  // int tmp_mv[2][72][92];  ([2][72][88] should be enough)
  memory_size += get_mem3Dint(&tmp_mv_frm, 2, img->height/BLOCK_SIZE, img->width/BLOCK_SIZE+4);

  // allocate memory for reference frames of each block: refFrArr
  // int  refFrArr[72][88];
  memory_size += get_mem2Dint(&refFrArr_frm, img->height/BLOCK_SIZE, img->width/BLOCK_SIZE);

  memory_size += get_mem2Dint(&abp_type_FrArr, img->height/BLOCK_SIZE, img->width/BLOCK_SIZE);

  if(input->successive_Bframe!=0 || input->StoredBPictures > 0)
  {    
    // allocate memory for temp B-frame motion vector buffer: fw_refFrArr, bw_refFrArr
    // int ...refFrArr[72][88];
    memory_size += get_mem2Dint(&fw_refFrArr_frm, img->height/BLOCK_SIZE, img->width/BLOCK_SIZE);
    memory_size += get_mem2Dint(&bw_refFrArr_frm, img->height/BLOCK_SIZE, img->width/BLOCK_SIZE);

    // allocate memory for temp B-frame direct ref-frame buffer: fwdir_refFrArr, bwdir_refFrArr
    // int ...refFrArr[72][88];
    if(input->direct_type == DIR_SPATIAL)
    {
      memory_size += get_mem2Dint(&fwdir_refFrArr, img->height/BLOCK_SIZE, img->width/BLOCK_SIZE);
      memory_size += get_mem2Dint(&bwdir_refFrArr, img->height/BLOCK_SIZE, img->width/BLOCK_SIZE);
      // allocate memory for collocated motion stationarity - int could be replaced with boolean    
      memory_size += get_mem2Dint(&moving_block_frm, img->height/BLOCK_SIZE, img->width/BLOCK_SIZE);
    }

  if(input->InterlaceCodingOption >= MB_CODING)
  {
    memory_size += get_mem3Dint(&tmp_fwMV_top, 2, img->height/BLOCK_SIZE, img->width/BLOCK_SIZE+4);
    memory_size += get_mem3Dint(&tmp_fwMV_bot, 2, img->height/BLOCK_SIZE, img->width/BLOCK_SIZE+4);
    memory_size += get_mem3Dint(&tmp_bwMV_top, 2, img->height/BLOCK_SIZE, img->width/BLOCK_SIZE+4);
    memory_size += get_mem3Dint(&tmp_bwMV_bot, 2, img->height/BLOCK_SIZE, img->width/BLOCK_SIZE+4);
  }
  }


  // allocate memory for post filter frame buffers: imgY_pf, imgUV_pf
  // byte imgY_pf[288][352];
  // byte imgUV_pf[2][144][176];
  memory_size += get_mem2D(&imgY_pf, img->height, img->width);
  memory_size += get_mem3D(&imgUV_pf, 2, img->height_cr, img->width_cr);

  // allocate memory for B frame coding: nextP_imgY, nextP_imgUV
  // byte nextP_imgY[288][352];
  // byte nextP_imgUV[2][144][176];
  memory_size += get_mem2D(&nextP_imgY, img->height, img->width);
  memory_size += get_mem3D(&nextP_imgUV, 2, img->height_cr, img->width_cr);

  // allocate memory for multiple ref. frame buffers: mref, mcref
  //byte mref[MAX_MULT_PRED][1152][1408];  */   /* 1/4 pix luma
  //byte mcef[MAX_MULT_PRED][2][352][288]; */   /* pix chroma
  // rows and cols for croma component mcef[ref][croma][4x][4y] are switched
  // compared to luma mref[ref][4y][4x] for whatever reason
  // number of reference frames increased by one for next P-frame

  alloc_mref(img);

  alloc_Refbuf (img);

  if(input->successive_Bframe!=0 || input->StoredBPictures > 0)
  {
    // allocate memory for temp B-frame motion vector buffer: tmp_fwMV, tmp_bwMV, dfMV, dbMV
    // int ...MV[2][72][92];  ([2][72][88] should be enough)
    memory_size += get_mem3Dint(&tmp_fwMV, 2, img->height/BLOCK_SIZE, img->width/BLOCK_SIZE+4);
    memory_size += get_mem3Dint(&tmp_bwMV, 2, img->height/BLOCK_SIZE, img->width/BLOCK_SIZE+4);
    memory_size += get_mem3Dint(&dfMV,     2, img->height/BLOCK_SIZE, img->width/BLOCK_SIZE+4);
    memory_size += get_mem3Dint(&dbMV,     2, img->height/BLOCK_SIZE, img->width/BLOCK_SIZE+4);

  }

  // allocate memory for temp quarter pel luma frame buffer: img4Y_tmp
  // int img4Y_tmp[576][704];  (previously int imgY_tmp in global.h)
  memory_size += get_mem2Dint(&img4Y_tmp, img->height+2*IMG_PAD_SIZE, (img->width+2*IMG_PAD_SIZE)*4);

  if (input->rdopt==2)
  {
    memory_size += get_mem2Dint(&decs->resY, MB_BLOCK_SIZE, MB_BLOCK_SIZE);
    if ((decs->decref = (byte****) calloc(input->NoOfDecoders,sizeof(byte***))) == NULL) 
      no_mem_exit("init_global_buffers: decref");
    for (j=0 ; j<input->NoOfDecoders; j++)
    {
      memory_size += get_mem3D(&decs->decref[j], img->buf_cycle+1, img->height, img->width);
    }
    memory_size += get_mem2D(&decs->RefBlock, BLOCK_SIZE,BLOCK_SIZE);
    memory_size += get_mem3D(&decs->decY, input->NoOfDecoders, img->height, img->width);
    memory_size += get_mem3D(&decs->decY_best, input->NoOfDecoders, img->height, img->width);
    memory_size += get_mem2D(&decs->status_map, img->height/MB_BLOCK_SIZE,img->width/MB_BLOCK_SIZE);
    memory_size += get_mem2D(&decs->dec_mb_mode, img->width/MB_BLOCK_SIZE,img->height/MB_BLOCK_SIZE);
  }
  if (input->RestrictRef)
  {
    memory_size += get_mem2D(&pixel_map, img->height,img->width);
    memory_size += get_mem2D(&refresh_map, img->height/8,img->width/8);
  }

  if(input->InterlaceCodingOption != FRAME_CODING)
  {
    // allocate memory for encoding frame buffers: imgY, imgUV
    // byte imgY[288][352];
    // byte imgUV[2][144][176];
    memory_size += get_mem2D(&imgY_com, img->height, img->width);
    memory_size += get_mem3D(&imgUV_com, 2, img->height/2, img->width_cr);
    memory_size += get_mem2D(&imgY_top, height_field, img->width);
    memory_size += get_mem3D(&imgUV_top, 2, height_field/2, img->width_cr);
    memory_size += get_mem2D(&imgY_bot, height_field, img->width);
    memory_size += get_mem3D(&imgUV_bot, 2, height_field/2, img->width_cr);

    // allocate memory for reference frame buffers: imgY_org, imgUV_org
    // byte imgY_org[288][352];
    // byte imgUV_org[2][144][176];
    memory_size += get_mem2D(&imgY_org_top, height_field, img->width);
    memory_size += get_mem3D(&imgUV_org_top, 2, height_field/2, img->width_cr);
    memory_size += get_mem2D(&imgY_org_bot, height_field, img->width);
    memory_size += get_mem3D(&imgUV_org_bot, 2, height_field/2, img->width_cr);

    if(input->successive_Bframe!=0 || input->StoredBPictures > 0)
    {
      // allocate memory for temp B-frame motion vector buffer: fw_refFrArr, bw_refFrArr
      // int ...refFrArr[72][88];
      memory_size += get_mem2Dint(&fw_refFrArr_top, height_field/BLOCK_SIZE, img->width/BLOCK_SIZE);
      memory_size += get_mem2Dint(&bw_refFrArr_top, height_field/BLOCK_SIZE, img->width/BLOCK_SIZE);
      memory_size += get_mem2Dint(&fw_refFrArr_bot, height_field/BLOCK_SIZE, img->width/BLOCK_SIZE);
      memory_size += get_mem2Dint(&bw_refFrArr_bot, height_field/BLOCK_SIZE, img->width/BLOCK_SIZE);
      if(input->direct_type == DIR_SPATIAL)
      {
        // allocate memory for collocated motion stationarity - int could be replaced with boolean    
        memory_size += get_mem2Dint(&moving_block_top,  height_field/BLOCK_SIZE, img->width/BLOCK_SIZE);
        memory_size += get_mem2Dint(&moving_block_bot,  height_field/BLOCK_SIZE, img->width/BLOCK_SIZE);
      }
    }

    // allocate memory for temp P and B-frame motion vector buffer: tmp_mv, temp_mv_block
    // int tmp_mv[2][72][92];  ([2][72][88] should be enough)
    memory_size += get_mem3Dint(&tmp_mv_top, 2, height_field/BLOCK_SIZE, img->width/BLOCK_SIZE+4);
    memory_size += get_mem3Dint(&tmp_mv_bot, 2, height_field/BLOCK_SIZE, img->width/BLOCK_SIZE+4);

    // allocate memory for reference frames of each block: refFrArr
    // int  refFrArr[72][88];
    memory_size += get_mem2Dint(&refFrArr_top, height_field/BLOCK_SIZE, img->width/BLOCK_SIZE);
    memory_size += get_mem2Dint(&refFrArr_bot, height_field/BLOCK_SIZE, img->width/BLOCK_SIZE);
    memory_size += get_mem2Dint(&abp_type_FrArr_top, height_field/BLOCK_SIZE, img->width/BLOCK_SIZE);
    memory_size += get_mem2Dint(&abp_type_FrArr_bot, height_field/BLOCK_SIZE, img->width/BLOCK_SIZE);

  if(input->InterlaceCodingOption >= MB_CODING)
  {
//    memory_size += get_mem3Dint(&tmp_mv_top_save, 2, height_field/BLOCK_SIZE, img->width/BLOCK_SIZE+4);
//    memory_size += get_mem3Dint(&tmp_mv_bot_save, 2, height_field/BLOCK_SIZE, img->width/BLOCK_SIZE+4);
//    memory_size += get_mem2Dint(&refFrArr_top_save, height_field/BLOCK_SIZE, img->width/BLOCK_SIZE);
//    memory_size += get_mem2Dint(&refFrArr_bot_save, height_field/BLOCK_SIZE, img->width/BLOCK_SIZE);
  }

  }

  return (memory_size);
}

/*!
 ************************************************************************
 * \brief
 *    Free allocated memory of frame size related global buffers
 *    buffers are defined in global.h, allocated memory is allocated in
 *    int get_mem4global_buffers()
 * \par Input:
 *    Input Parameters struct inp_par *inp,                             \n
 *    Image Parameters struct img_par *img
 * \par Output:
 *    none
 ************************************************************************
 */
void free_global_buffers()
{
  int  i,j;

#ifdef _ADAPT_LAST_GROUP_
  extern int *last_P_no_frm;
  extern int *last_P_no_fld;
  free (last_P_no_frm);
  free (last_P_no_fld);
#endif

  free_mem2D(imgY_frm);
  free_mem3D(imgUV_frm,2);
  free_mem2D(imgY_org_frm);      // free ref frame buffers
  free_mem3D(imgUV_org_frm,2);
  free_mem3Dint(tmp_mv_frm,2);
  free_mem2Dint(refFrArr_frm);

  // free multiple ref frame buffers
  // number of reference frames increased by one for next P-frame
  free(mref_frm);
  if (input->WeightedPrediction || input->WeightedBiprediction)
    free(mref_frm_w);
  free(mcef_frm);

  free_mem2D(imgY_pf);       // free post filtering frame buffers
  free_mem3D(imgUV_pf,2);

  free_mem2D(nextP_imgY);    // free next frame buffers (for B frames)
  free_mem3D(nextP_imgUV,2);

  free (Refbuf11_frm);
  if (input->WeightedPrediction || input->WeightedBiprediction)
        free (Refbuf11_frm_w);

  // free multiple ref frame buffers
  // number of reference frames increased by one for next P-frame

  if(input->successive_Bframe!=0 || input->StoredBPictures > 0)
  {
    // free last P-frame buffers for B-frame coding
    free_mem3Dint(tmp_fwMV,2);
    free_mem3Dint(tmp_bwMV,2);
    free_mem3Dint(dfMV,2);
    free_mem3Dint(dbMV,2);
    free_mem2Dint(fw_refFrArr_frm);
    free_mem2Dint(bw_refFrArr_frm);    
    if (input->direct_type)
    {
      free_mem2Dint(moving_block_frm);
      free_mem2Dint(fwdir_refFrArr);
      free_mem2Dint(bwdir_refFrArr);
    }
  } // end if B frame
  free_mem2Dint(abp_type_FrArr);

  free_mem2Dint(img4Y_tmp);    // free temp quarter pel frame buffer

  // free mem, allocated in init_img()
  // free intra pred mode buffer for blocks
  free_mem2Dint(img->ipredmode);
  free(img->mb_data);

  if(input->UseConstrainedIntraPred)
  {
    j=(img->width/16)*(img->height/16);
    for (i=0; i<j; i++)
    {
      free (img->intra_block[i]);
    }
    free (img->intra_block);
  }

  if (input->rdopt==2)
  {
    free(decs->resY[0]);
    free(decs->resY);
    free(decs->RefBlock[0]);
    free(decs->RefBlock);
    for (j=0; j<input->NoOfDecoders; j++)
    {
      free(decs->decY[j][0]);
      free(decs->decY[j]);
      free(decs->decY_best[j][0]);
      free(decs->decY_best[j]);
      for (i=0; i<img->buf_cycle+1; i++)
      {
        free(decs->decref[j][i][0]);
        free(decs->decref[j][i]);
      }
      free(decs->decref[j]);
    }
    free(decs->decY);
    free(decs->decY_best);
    free(decs->decref);
    free(decs->status_map[0]);
    free(decs->status_map);
    free(decs->dec_mb_mode[0]);
    free(decs->dec_mb_mode);
  }
  if (input->RestrictRef)
  {
    free(pixel_map[0]);
    free(pixel_map);
    free(refresh_map[0]);
    free(refresh_map);
  }

  if(input->InterlaceCodingOption != FRAME_CODING)
  {
    free_mem2D(imgY_com);
    free_mem3D(imgUV_com,2);
    free_mem2D(imgY_top);
    free_mem3D(imgUV_top,2);
    free_mem2D(imgY_org_top);      // free ref frame buffers
    free_mem3D(imgUV_org_top,2);
    free_mem2D(imgY_bot);
    free_mem3D(imgUV_bot,2);
    free_mem2D(imgY_org_bot);      // free ref frame buffers
    free_mem3D(imgUV_org_bot,2);

    // free multiple ref frame buffers
    // number of reference frames increased by one for next P-frame
    free(mref_fld);
    if (input->WeightedPrediction || input->WeightedBiprediction)
      free(mref_fld_w);
    free(mcef_fld);

    if(input->successive_Bframe!=0 || input->StoredBPictures > 0)
    {
      // free last P-frame buffers for B-frame coding
      free_mem2Dint(fw_refFrArr_top);
      free_mem2Dint(bw_refFrArr_top);
      free_mem2Dint(fw_refFrArr_bot);
      free_mem2Dint(bw_refFrArr_bot);
      if (input->direct_type)
      {
        free_mem2Dint(moving_block_top);
        free_mem2Dint(moving_block_bot);
      }
    } // end if B frame

    free (Refbuf11_fld);
        if ( input->WeightedPrediction || input->WeightedBiprediction)
                free (Refbuf11_fld_w);

    free_mem3Dint(tmp_mv_top,2);
    free_mem3Dint(tmp_mv_bot,2);
    free_mem2Dint(refFrArr_top);
    free_mem2Dint(refFrArr_bot);
    free_mem2Dint(abp_type_FrArr_top);
    free_mem2Dint(abp_type_FrArr_bot);
  }

  // CAVLC mem free
  for(j=0;j<img->width/MB_BLOCK_SIZE;j++)
  for(i=0;i<img->height/MB_BLOCK_SIZE;i++)
  {
    if(img->nz_coeff[j][i][0] != NULL) free(img->nz_coeff[j][i][0]);
    if(img->nz_coeff[j][i]    != NULL) free(img->nz_coeff[j][i]);
  };
  if (img->nz_coeff !=NULL) free(img->nz_coeff );



}

/*!
 ************************************************************************
 * \brief
 *    Allocate memory for mv
 * \par Input:
 *    Image Parameters struct img_par *img                             \n
 *    int****** mv
 * \return memory size in bytes
 ************************************************************************
 */
int get_mem_mv (int****** mv)
{
  int i, j, k, l;

  if ((*mv = (int*****)calloc(4,sizeof(int****))) == NULL)
    no_mem_exit ("get_mem_mv: mv");
  for (i=0; i<4; i++)
  {
    if (((*mv)[i] = (int****)calloc(4,sizeof(int***))) == NULL)
      no_mem_exit ("get_mem_mv: mv");
    for (j=0; j<4; j++)
    {
      if (((*mv)[i][j] = (int***)calloc(img->buf_cycle,sizeof(int**))) == NULL)
        no_mem_exit ("get_mem_mv: mv");
      for (k=0; k<img->buf_cycle; k++)
      {
        if (((*mv)[i][j][k] = (int**)calloc(9,sizeof(int*))) == NULL)
          no_mem_exit ("get_mem_mv: mv");
        for (l=0; l<9; l++)
          if (((*mv)[i][j][k][l] = (int*)calloc(2,sizeof(int))) == NULL)
            no_mem_exit ("get_mem_mv: mv");
      }
    }
  }
  return 4*4*img->buf_cycle*9*2*sizeof(int);
}


/*!
 ************************************************************************
 * \brief
 *    Free memory from mv
 * \par Input:
 *    int****** mv
 ************************************************************************
 */
void free_mem_mv (int***** mv)
{
  int i, j, k, l;

  for (i=0; i<4; i++)
  {
    for (j=0; j<4; j++)
    {
      for (k=0; k<img->buf_cycle; k++)
      {
        for (l=0; l<9; l++)
          free (mv[i][j][k][l]);
        free (mv[i][j][k]);
      }
      free (mv[i][j]);
    }
    free (mv[i]);
  }
  free (mv);
}





/*!
 ************************************************************************
 * \brief
 *    Allocate memory for AC coefficients
 ************************************************************************
 */
int get_mem_ACcoeff (int***** cofAC)
{
  int i, j, k;

  if ((*cofAC = (int****)calloc (6, sizeof(int***))) == NULL)              no_mem_exit ("get_mem_ACcoeff: cofAC");
  for (k=0; k<6; k++)
  {
    if (((*cofAC)[k] = (int***)calloc (4, sizeof(int**))) == NULL)         no_mem_exit ("get_mem_ACcoeff: cofAC");
    for (j=0; j<4; j++)
    {
      if (((*cofAC)[k][j] = (int**)calloc (2, sizeof(int*))) == NULL)      no_mem_exit ("get_mem_ACcoeff: cofAC");
      for (i=0; i<2; i++)
      {
        if (((*cofAC)[k][j][i] = (int*)calloc (65, sizeof(int))) == NULL)  no_mem_exit ("get_mem_ACcoeff: cofAC"); // 18->65 for ABT
      }
    }
  }
  return 6*4*2*65*sizeof(int);// 18->65 for ABT
}

/*!
 ************************************************************************
 * \brief
 *    Allocate memory for DC coefficients
 ************************************************************************
 */
int get_mem_DCcoeff (int**** cofDC)
{
  int j, k;

  if ((*cofDC = (int***)calloc (3, sizeof(int**))) == NULL)           no_mem_exit ("get_mem_DCcoeff: cofDC");
  for (k=0; k<3; k++)
  {
    if (((*cofDC)[k] = (int**)calloc (2, sizeof(int*))) == NULL)      no_mem_exit ("get_mem_DCcoeff: cofDC");
    for (j=0; j<2; j++)
    {
      if (((*cofDC)[k][j] = (int*)calloc (65, sizeof(int))) == NULL)  no_mem_exit ("get_mem_DCcoeff: cofDC"); // 18->65 for ABT
    }
  }
  return 3*2*65*sizeof(int); // 18->65 for ABT
}


/*!
 ************************************************************************
 * \brief
 *    Free memory of AC coefficients
 ************************************************************************
 */
void free_mem_ACcoeff (int**** cofAC)
{
  int i, j, k;

  for (k=0; k<6; k++)
  {
    for (i=0; i<4; i++)
    {
      for (j=0; j<2; j++)
      {
        free (cofAC[k][i][j]);
      }
      free (cofAC[k][i]);
    }
    free (cofAC[k]);
  }
  free (cofAC);
}

/*!
 ************************************************************************
 * \brief
 *    Free memory of DC coefficients
 ************************************************************************
 */
void free_mem_DCcoeff (int*** cofDC)
{
  int i, j;

  for (j=0; j<3; j++)
  {
    for (i=0; i<2; i++)
    {
      free (cofDC[j][i]);
    }
    free (cofDC[j]);
  }
  free (cofDC);
}


/*!
 ************************************************************************
 * \brief
 *    form frame picture from two field pictures 
 ************************************************************************
 */
void combine_field()
{
  int i;

  for (i=0; i<img->height / 2; i++)
  {
    memcpy(imgY_com[i*2], imgY_top[i], img->width);     // top field
    memcpy(imgY_com[i*2 + 1], imgY_bot[i], img->width); // bottom field
  }

  for (i=0; i<img->height_cr / 2; i++)
  {
    memcpy(imgUV_com[0][i*2], imgUV_top[0][i], img->width_cr);
    memcpy(imgUV_com[0][i*2 + 1], imgUV_bot[0][i], img->width_cr);
    memcpy(imgUV_com[1][i*2], imgUV_top[1][i], img->width_cr);
    memcpy(imgUV_com[1][i*2 + 1], imgUV_bot[1][i], img->width_cr);
  }
}

/*!
 ************************************************************************
 * \brief
 *    extract top field from a frame 
 ************************************************************************
 */
void split_field_top()
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
 *    extract bottom field from a frame 
 ************************************************************************
 */
void split_field_bot()
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
 *    RD decision of frame and field coding 
 ************************************************************************
 */
int decide_fld_frame(float snr_frame_Y, float snr_field_Y, int bit_field, int bit_frame, double lambda_picture)
{
  double cost_frame, cost_field;

  cost_frame = bit_frame * lambda_picture + snr_frame_Y;
  cost_field = bit_field * lambda_picture + snr_field_Y;

  if (cost_field > cost_frame)
    return (0);
  else
    return (1);
}

/*!
 ************************************************************************
 * \brief
 *    Do some initializaiton work for encoding the 2nd IGOP
 ************************************************************************
 */
void process_2nd_IGOP()
{
  Boolean FirstIGOPFinished = FALSE;
  if ( img->number == input->no_frames-1 )
    FirstIGOPFinished = TRUE;
  if (input->NumFrameIn2ndIGOP==0) return;
  if (!FirstIGOPFinished || In2ndIGOP) return;
  In2ndIGOP = TRUE;

//  img->number = -1;
  img->pn = -1;
  start_frame_no_in_this_IGOP = input->no_frames;
  start_tr_in_this_IGOP = (input->no_frames-1)*(input->jumpd+1) +1;
  input->no_frames = input->no_frames + input->NumFrameIn2ndIGOP;

  reset_buffers();

  frm->picbuf_short[0]->used=0;
  frm->picbuf_short[0]->picID=-1;
  frm->picbuf_short[0]->lt_picID=-1;
  frm->short_used = 0;
  img->nb_references = 0;
}

/*!
 ************************************************************************
 * \brief
 *    SetImgType
 ************************************************************************
 */

void SetImgType()
{
  if (input->intra_period == 0)
  {
    if (IMG_NUMBER == 0)
    {
      img->type = INTRA_IMG;        // set image type for first image to I-frame
    }
    else
    {
      img->type = INTER_IMG;        // P-frame
      if(input->StoredBPictures > 0 && img->number > 1)
      {
        img->type = BS_IMG;
      }
      else if (input->sp_periodicity)
      {
        if ((IMG_NUMBER % input->sp_periodicity) ==0)
        {
          img->types=SP_IMG;
        }
        else img->types=INTER_IMG;
      }
    }
  }
  else
  {
    if ((IMG_NUMBER%input->intra_period) == 0)
    {
      img->type = INTRA_IMG;
    }
    else
    {
      img->type = INTER_IMG;        // P-frame
      if(input->StoredBPictures > 0 && img->number > 1)
      {
        img->type = BS_IMG;
      }
      else if (input->sp_periodicity)
      {
        if ((IMG_NUMBER % input->sp_periodicity) ==0)
            img->types=SP_IMG;
        else img->types=INTER_IMG;
      }
    }
  }

  if (frm->short_size==frm->short_used)
  {
    frm->picbuf_short[frm->short_size-1]->used = 0;
    frm->picbuf_short[frm->short_size-1]->picID = -1;
    frm->picbuf_short[frm->short_size-1]->lt_picID = -1;
    frm->short_used--;
    img->nb_references--;
  }
}


/*!
 ************************************************************************
 * \brief
 *    poc_distance()
 ************************************************************************
 */

int poc_distance( int refa, int refb)
{

  return toprefpoc[refb + 1] - toprefpoc[refa + 1];
}
 
