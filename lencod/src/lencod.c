/*!
 ***********************************************************************
 *  \mainpage
 *     This is the H.264/AVC encoder reference software. For detailed documentation
 *     see the comments in each file.
 *
 *     The JM software web site is located at:
 *     http://iphome.hhi.de/suehring/tml
 *
 *     For bug reporting and known issues see:
 *     https://ipbt.hhi.de
 *
 *  \author
 *     The main contributors are listed in contributors.h
 *
 *  \version
 *     JM 16.0 (FRExt)
 *
 *  \note
 *     tags are used for document system "doxygen"
 *     available at http://www.doxygen.org
 */
/*!
 *  \file
 *     lencod.c
 *  \brief
 *     H.264/AVC reference encoder project main()
 *  \author
 *   Main contributors (see contributors.h for copyright, address and affiliation details)
 *   - Inge Lille-Langoy               <inge.lille-langoy@telenor.com>
 *   - Rickard Sjoberg                 <rickard.sjoberg@era.ericsson.se>
 *   - Stephan Wenger                  <stewe@cs.tu-berlin.de>
 *   - Jani Lainema                    <jani.lainema@nokia.com>
 *   - Byeong-Moon Jeon                <jeonbm@lge.com>
 *   - Yoon-Seong Soh                  <yunsung@lge.com>
 *   - Thomas Stockhammer              <stockhammer@ei.tum.de>
 *   - Detlev Marpe                    <marpe@hhi.de>
 *   - Guido Heising
 *   - Valeri George                   <george@hhi.de>
 *   - Karsten Suehring                <suehring@hhi.de>
 *   - Alexis Michael Tourapis         <alexismt@ieee.org>
 ***********************************************************************
 */

#include "contributors.h"

#include <time.h>
#include <math.h>

#include "global.h"
#include "cconv_yuv2rgb.h"
#include "configfile.h"
#include "context_ini.h"
#include "explicit_gop.h"
#include "explicit_seq.h"
#include "filehandle.h"
#include "image.h"
#include "input.h"
#include "img_io.h"
#include "slice.h"
#include "intrarefresh.h"
#include "leaky_bucket.h"
#include "mc_prediction.h"
#include "memalloc.h"
#include "me_epzs_common.h"
#include "me_epzs_int.h"
#include "me_umhex.h"
#include "me_umhexsmp.h"
#include "output.h"
#include "parset.h"
#include "q_matrix.h"
#include "q_offsets.h"
#include "ratectl.h"
#include "report.h"
#include "rdoq.h"
#include "errdo.h"
#include "rdopt.h"
#include "wp_mcprec.h"
#include "mv_search.h"
#include "img_process.h"
#include "q_offsets.h"

static const int mb_width_cr[4] = {0,8, 8,16};
static const int mb_height_cr[4]= {0,8,16,16};

EncoderParams   *p_Enc = NULL;

static void SetLevelIndices(ImageParameters *p_Img);
static void chroma_mc_setup(ImageParameters *p_Img);

static int  init_global_buffers (ImageParameters *p_Img, InputParameters *p_Inp);
static void free_global_buffers (ImageParameters *p_Img, InputParameters *p_Inp);
static void free_img            (ImageParameters *p_Img, InputParameters *p_Inp);
static void free_params         (InputParameters *p_Inp);

static void init_img       (ImageParameters *p_Img, InputParameters *p_Inp);
static void init_poc       (ImageParameters *p_Img, InputParameters *p_Inp);
static void init_encoder   (ImageParameters *p_Img, InputParameters *p_Inp);
static void encode_sequence(ImageParameters *p_Img, InputParameters *p_Inp);

void init_stats (InputParameters *p_Inp, StatParameters *p_Stats)
{
  memset(p_Stats, 0, sizeof(StatParameters));
  p_Stats->NumberBFrames = p_Inp->NumberBFrames;
}

void init_dstats (DistortionParams *p_Dist)
{
  p_Dist->frame_ctr = 0;
  memset(p_Dist->metric, 0, TOTAL_DIST_TYPES * sizeof(DistMetric));
}

/*!
 ***********************************************************************
 * \brief
 *    Initialize encoding parameters.
 ***********************************************************************
 */
static void init_frame_params(ImageParameters *p_Img, InputParameters *p_Inp)
{
  int base_mul = 0;

  if (p_Inp->idr_period)
  {
    if (!p_Inp->adaptive_idr_period && ( p_Img->frm_number - p_Img->lastIDRnumber ) % p_Inp->idr_period == 0 )
      p_Img->nal_reference_idc = NALU_PRIORITY_HIGHEST;

    if (p_Inp->adaptive_idr_period == 1 && ( p_Img->frm_number - imax(p_Img->lastIntraNumber, p_Img->lastIDRnumber) ) % p_Inp->idr_period == 0 )
      p_Img->nal_reference_idc = NALU_PRIORITY_HIGHEST;
    else
      p_Img->nal_reference_idc = (p_Inp->DisposableP) ? (p_Img->frm_number + 1) & 0x01 : NALU_PRIORITY_LOW;

  }
  else
    p_Img->nal_reference_idc = (p_Img->frm_number && p_Inp->DisposableP) ? (p_Img->frm_number + 1) & 0x01 : NALU_PRIORITY_LOW;

  //much of this can go in init_frame() or init_field()?
  //poc for this frame or field
  if (p_Inp->idr_period)
  {
    if (!p_Inp->adaptive_idr_period)
      base_mul = ( p_Img->frm_number - p_Img->lastIDRnumber ) % p_Inp->idr_period;
    else if (p_Inp->adaptive_idr_period == 1)
      base_mul = (( p_Img->frm_number - imax(p_Img->lastIntraNumber, p_Img->lastIDRnumber) ) % p_Inp->idr_period == 0) ? 0 : ( p_Img->frm_number - p_Img->lastIDRnumber );
  }
  else 
    base_mul = ( p_Img->frm_number - p_Img->lastIDRnumber );

  if ((p_Img->frm_number - p_Img->lastIDRnumber) <= p_Inp->intra_delay)
  {    
    base_mul = -base_mul;
  }
  else
  {
    base_mul -= ( base_mul ? p_Inp->intra_delay :  0);    
  }

  p_Img->toppoc = base_mul * (2 * p_Img->base_dist);

  if ((p_Inp->PicInterlace==FRAME_CODING) && (p_Inp->MbInterlace==FRAME_CODING))
    p_Img->bottompoc = p_Img->toppoc;     //progressive
  else
    p_Img->bottompoc = p_Img->toppoc + 1;   //hard coded

  p_Img->framepoc = imin (p_Img->toppoc, p_Img->bottompoc);

  //the following is sent in the slice header
  p_Img->delta_pic_order_cnt[0] = 0;

  if ((p_Inp->BRefPictures == 1) && (p_Img->frm_number))
  {
    p_Img->delta_pic_order_cnt[0] = 2 * p_Inp->NumberBFrames;
  }  

  if (p_Inp->NumberBFrames && p_Inp->last_frame && ((p_Img->gop_number) + 1) == p_Inp->no_frm_base)
  {
    int bi = (int)((float)p_Img->base_dist / (p_Img->initial_Bframes + 1.0) + 0.499999);
    int new_bframes = ((p_Inp->last_frame - (p_Img->frm_number - 1) * p_Img->base_dist) / bi) - 1;

    //about to code the last ref frame, adjust delta poc
    p_Img->delta_pic_order_cnt[0]= -2*(p_Img->initial_Bframes - new_bframes);
    p_Img->toppoc    += p_Img->delta_pic_order_cnt[0];
    p_Img->bottompoc += p_Img->delta_pic_order_cnt[0];
    p_Img->framepoc   = imin (p_Img->toppoc, p_Img->bottompoc);
  }

  //frame_num for this frame
  if (p_Inp->idr_period && ((!p_Inp->adaptive_idr_period && ( p_Img->frm_number - p_Img->lastIDRnumber ) % p_Inp->idr_period == 0)
    || (p_Inp->adaptive_idr_period == 1 && ( p_Img->frm_number - imax(p_Img->lastIntraNumber, p_Img->lastIDRnumber) ) % p_Inp->idr_period == 0)) )
  {
    p_Img->frame_num = 0;
  }
}

/*!
 ***********************************************************************
 * \brief
 *    Allocate the Image structure
 * \par  Output:
 *    Image Parameters ImageParameters *p_Img
 ***********************************************************************
 */
static void alloc_img( ImageParameters **p_Img)
{
  if ((*p_Img = (ImageParameters *) calloc(1, sizeof(ImageParameters)))==NULL) 
    no_mem_exit("alloc_img: p_Img");
  if ((((*p_Img)->p_Dist)  = (DistortionParams *) calloc(1, sizeof(DistortionParams)))==NULL) 
    no_mem_exit("alloc_img: p_Dist");
  if ((((*p_Img)->p_Stats) = (StatParameters *) calloc(1, sizeof(StatParameters)))==NULL) 
    no_mem_exit("alloc_img: p_Stats");
  if (((*p_Img)->p_Dpb     = (DecodedPictureBuffer *) calloc(1, sizeof(DecodedPictureBuffer)))==NULL) 
    no_mem_exit("alloc_img: p_Dpb");
  if ((((*p_Img)->p_Quant)  = (QuantParameters *) calloc(1, sizeof(QuantParameters)))==NULL) 
    no_mem_exit("alloc_img: p_Quant");
  if ((((*p_Img)->p_QScale)  = (ScaleParameters *) calloc(1, sizeof(ScaleParameters)))==NULL) 
    no_mem_exit("alloc_img: p_QScale");
  if ((((*p_Img)->p_SEI)  = (SEIParameters *) calloc(1, sizeof(SEIParameters)))==NULL) 
    no_mem_exit("alloc_img: p_SEI");


  (*p_Img)->p_dec = -1;  
  (*p_Img)->p_log = NULL;
  (*p_Img)->f_annexb = NULL;
  // Init rtp related info
  (*p_Img)->f_rtp = NULL;
  (*p_Img)->CurrentRTPTimestamp = 0;         
  (*p_Img)->CurrentRTPSequenceNumber = 0;
}

/*!
 ***********************************************************************
 * \brief
 *    Allocate the Input structure
 * \par  Output:
 *    Input Parameters InputParameters *p_Img
 ***********************************************************************
 */
static void alloc_params( InputParameters **p_Inp )
{
  if ((*p_Inp = (InputParameters *) calloc(1, sizeof(InputParameters)))==NULL) 
    no_mem_exit("alloc_params: p_Inp");

  (*p_Inp)->top_left          = NULL;
  (*p_Inp)->bottom_right      = NULL;
  (*p_Inp)->slice_group_id    = NULL;
  (*p_Inp)->run_length_minus1 = NULL;
}


  /*!
 ***********************************************************************
 * \brief
 *    Allocate the Encoder Structure
 * \par  Output:
 *    Encoder Parameters
 ***********************************************************************
 */
static void alloc_encoder( EncoderParams **p_Enc)
{
  if ((*p_Enc = (EncoderParams *) calloc(1, sizeof(EncoderParams)))==NULL) 
    no_mem_exit("alloc_encoder: p_Enc");

  alloc_img(&((*p_Enc)->p_Img));
  alloc_params(&((*p_Enc)->p_Inp));
  (*p_Enc)->p_trace = NULL;
  (*p_Enc)->bufferSize = 0;
}

/*!
 ***********************************************************************
 * \brief
 *    Free the Encoder Structure
 ***********************************************************************
 */
static void free_encoder (EncoderParams *p_Enc)
{
  if ( p_Enc != NULL )
  {
    free( p_Enc );
  }
}

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
int main(int argc, char **argv)
{
  
  alloc_encoder(&p_Enc);

  Configure (p_Enc->p_Img, p_Enc->p_Inp, argc, argv);

  // init encoder
  init_encoder(p_Enc->p_Img, p_Enc->p_Inp);

  // encode sequence
  encode_sequence(p_Enc->p_Img, p_Enc->p_Inp);

  // terminate sequence
  free_encoder_memory(p_Enc->p_Img, p_Enc->p_Inp);

  free_params (p_Enc->p_Inp);  
  free_encoder(p_Enc);

  return 0;
}


/*!
 ***********************************************************************
 * \brief
 *    Initialize encoder
 ***********************************************************************
 */

static void init_encoder(ImageParameters *p_Img, InputParameters *p_Inp)
{
  p_Img->giRDOpt_B8OnlyFlag = FALSE;

  p_Img->p_log = NULL;

  p_Img->cabac_encoding = 0;

  p_Img->frame_statistic_start = 1;

  // Open Files
  OpenFiles(&p_Inp->input_file1);


  Init_QMatrix(p_Img, p_Inp);
  Init_QOffsetMatrix(p_Img, p_Inp);

  init_poc(p_Img, p_Inp);
  GenerateParameterSets(p_Img, p_Inp);
  SetLevelIndices(p_Img);

  init_img  (p_Img, p_Inp);

  if (p_Inp->rdopt == 3)
  {
    init_error_conceal(p_Img,p_Inp->ErrorConcealment); 
  }

#ifdef _LEAKYBUCKET_
  p_Img->initial_Bframes = 0;
  p_Img->Bit_Buffer = (long *)malloc((p_Inp->no_frames + 1) * sizeof(long));
  p_Img->total_frame_buffer = 0;
#endif

  // Prepare hierarchical coding structures. 
  // Code could be extended in the future to allow structure adaptation.
  if (p_Inp->HierarchicalCoding)
  {
    init_gop_structure(p_Img, p_Inp);
    if (p_Inp->NumberBFrames && p_Inp->HierarchicalCoding == 3)
      interpret_gop_structure(p_Img, p_Inp);
    else
      create_hierarchy(p_Img, p_Inp);
  }

  p_Img->p_Dpb->init_done = 0;

  init_dpb(p_Img, p_Inp, p_Img->p_Dpb);
  init_out_buffer(p_Img);
  init_stats (p_Inp, p_Img->p_Stats);
  init_dstats(p_Img->p_Dist);

  p_Img->enc_picture = NULL;

  init_global_buffers(p_Img, p_Inp);

  if ( p_Inp->WPMCPrecision )
  {
    wpxInitWPXPasses(p_Img, p_Inp);
  }

  Init_Motion_Search_Module (p_Img, p_Inp);
  information_init(p_Img, p_Inp, p_Img->p_Stats);

  if(p_Inp->DistortionYUVtoRGB)
    init_YUVtoRGB(p_Img, p_Inp);

  //Rate control
  if (p_Inp->RCEnable)
    rc_init_sequence(p_Img, p_Inp);

  p_Img->last_valid_reference = 0;
  p_Img->tot_time = 0;                 // time for total encoding session

  p_Img->initial_Bframes = p_Inp->NumberBFrames;

  PatchInputNoFrames(p_Inp);

  p_Img->type = I_SLICE;
  // Write sequence header (with parameter sets)
  p_Img->p_Stats->bit_ctr_filler_data = 0;
  p_Img->p_Stats->bit_ctr_filler_data_n = 0;
  p_Img->p_Stats->bit_ctr_parametersets = 0;
  p_Img->p_Stats->bit_slice = start_sequence(p_Img, p_Inp);

  if (p_Inp->UseRDOQuant)
    precalculate_unary_exp_golomb_level(p_Img);

   if (p_Inp->ExplicitSeqCoding)
     OpenExplicitSeqFile(p_Img, p_Inp);

  if ( p_Inp->ChromaMCBuffer )
    p_Img->OneComponentChromaPrediction4x4 = OneComponentChromaPrediction4x4_retrieve;
  else
    p_Img->OneComponentChromaPrediction4x4 = OneComponentChromaPrediction4x4_regenerate;

  p_Img->searchRange.min_x = -p_Inp->search_range << 2;
  p_Img->searchRange.max_x =  p_Inp->search_range << 2;
  p_Img->searchRange.min_y = -p_Inp->search_range << 2;
  p_Img->searchRange.max_y =  p_Inp->search_range << 2;

}

/*!
 ***********************************************************************
 * \brief
 *    Determine coding level a frame belongs to
 ***********************************************************************
 */
static int determine_coding_level(ImageParameters *p_Img, InputParameters *p_Inp, int curr_frame)
{
  int coding_level = 0;  

  if (curr_frame - p_Img->last_idr_number == 0)
    return coding_level;
  else
    coding_level  = (curr_frame - p_Img->last_idr_number - 1) % (1 + p_Inp->NumberBFrames);
  
  return coding_level;
}

/*!
 ************************************************************************
 * \brief
 *    Set the image type for I,P and SP pictures (not B!)
 ************************************************************************
 */
static void SetImgType(ImageParameters *p_Img, InputParameters *p_Inp, int gop_frame_num)
{
  if (gop_frame_num == 0)
  {
    int intra_refresh = (p_Inp->intra_period == 0) ? (p_Img->gop_number == 0) : (( ( p_Img->frm_number - p_Img->lastIntraNumber) % p_Inp->intra_period ) == 0);
    int idr_refresh;

    if ( p_Inp->idr_period && !p_Inp->adaptive_idr_period )
      idr_refresh = (( ( p_Img->frm_number - p_Img->lastIDRnumber  ) % p_Inp->idr_period   ) == 0);
    else if ( p_Inp->idr_period && p_Inp->adaptive_idr_period == 1 )
      idr_refresh = (( ( p_Img->frm_number - imax(p_Img->lastIntraNumber, p_Img->lastIDRnumber)  ) % p_Inp->idr_period   ) == 0);
    else
      idr_refresh = (p_Img->gop_number == 0);

    if (intra_refresh || idr_refresh)
    {
      set_slice_type( p_Img, p_Inp, I_SLICE );        // set image type for first image to I-frame
    }
    else
    {
      set_slice_type(p_Img, p_Inp, (p_Inp->sp_periodicity && ((p_Img->gop_number % p_Inp->sp_periodicity) == 0))
        ? SP_SLICE  : ((p_Inp->BRefPictures == 2) ? B_SLICE : P_SLICE) );
    }
  }
  else
  {
    if (p_Inp->HierarchicalCoding)
      set_slice_type( p_Img, p_Inp, p_Img->gop_structure[gop_frame_num - 1].slice_type);
    else
      set_slice_type( p_Img, p_Inp, ( p_Inp->PReplaceBSlice ) ? P_SLICE : B_SLICE);
  }
}

static void set_poc(ImageParameters *p_Img, InputParameters *p_Inp, double frame_to_code)
{
  p_Img->toppoc = (int) (p_Img->frame_interval * frame_to_code);
  if (p_Img->gop_number && (p_Img->gop_number <= p_Inp->intra_delay))
  {
    if(p_Inp->idr_period && ((!p_Inp->adaptive_idr_period && ( p_Img->frm_number - p_Img->lastIDRnumber ) % p_Inp->idr_period == 0)
      || (p_Inp->adaptive_idr_period == 1 && ( p_Img->frm_number - imax(p_Img->lastIntraNumber, p_Img->lastIDRnumber) ) % p_Inp->idr_period == 0)) )
      p_Img->toppoc = p_Img->toppoc;
    else
    {
      p_Img->toppoc += (-p_Img->frm_number + p_Img->lastIDRnumber) * p_Img->base_dist;
      p_Img->toppoc *= 2;
    }
  }
  else
  {
    if(p_Inp->idr_period && !p_Inp->adaptive_idr_period)
      p_Img->toppoc += ((( p_Img->frm_number - p_Img->lastIDRnumber - p_Inp->intra_delay) % p_Inp->idr_period ) - 1) * p_Img->base_dist;
    else
      p_Img->toppoc += (( p_Img->frm_number - p_Img->lastIDRnumber - p_Inp->intra_delay - 1 ) * p_Img->base_dist);
    p_Img->toppoc *= 2;
  }
}

/*!
************************************************************************
* \brief
*    Prepare first coding layer.
************************************************************************
*/
static void prepare_first_layer(ImageParameters *p_Img, InputParameters *p_Inp, int curr_frame_to_code)
{
  p_Img->number ++;
  p_Img->gop_number = (p_Img->number - p_Img->start_frame_no);
  p_Img->frm_number = p_Img->number;

  p_Img->frm_no_in_file = CalculateFrameNumber(p_Img, p_Inp);

  if (p_Inp->last_frame != 0 && (p_Inp->last_frame == p_Img->frm_no_in_file) && p_Inp->HierarchicalCoding)
  { 
    int numberBFrames = p_Inp->NumberBFrames;
    p_Inp->HierarchicalCoding = imin( 2, p_Inp->HierarchicalCoding);    
    p_Inp->NumberBFrames = p_Inp->no_frames - curr_frame_to_code - 1;    

    clear_gop_structure(p_Img);
    init_gop_structure(p_Img, p_Inp);
    create_hierarchy(p_Img, p_Inp);

    p_Inp->NumberBFrames = numberBFrames;
  }
  SetImgType(p_Img, p_Inp, 0);

  init_frame_params(p_Img, p_Inp);

  //Rate control
  if (p_Inp->RCEnable && p_Img->type == I_SLICE)
    rc_init_gop_params(p_Img, p_Inp);

  // which layer does the image belong to?
  p_Img->layer = (p_Img->gop_number % (p_Inp->NumFramesInELSubSeq + 1)) ? 0 : 1;
}

/*!
************************************************************************
* \brief
*    Prepare second coding layer.
************************************************************************
*/
static void prepare_second_layer(ImageParameters *p_Img, InputParameters *p_Inp, int enh_frame_to_code)
{
  p_Img->layer = (p_Inp->NumFramesInELSubSeq == 0) ? 0 : 1;      
  SetImgType(p_Img, p_Inp, enh_frame_to_code);

  if ((p_Img->gop_number > 0) && (p_Inp->EnableIDRGOP == 0 || p_Img->idr_gop_number)) // B-frame(s) to encode
  {
    if (p_Inp->HierarchicalCoding)
    {
      p_Img->nal_reference_idc = p_Img->gop_structure[enh_frame_to_code - 1].reference_idc;
      set_poc(p_Img, p_Inp, (double)(1 + p_Img->gop_structure[enh_frame_to_code - 1].display_no));

      if (p_Img->gop_number && (p_Img->gop_number <= p_Inp->intra_delay))
      {
        if (enh_frame_to_code == 1)
          p_Img->delta_pic_order_cnt[0] = p_Img->toppoc - 2*(p_Img->start_tr_gop  + (p_Inp->intra_delay - p_Img->gop_number)*(p_Img->base_dist));
        else
          p_Img->delta_pic_order_cnt[0] = p_Img->toppoc - 2*(p_Img->start_tr_gop  + (p_Inp->intra_delay - p_Img->gop_number)*(p_Img->base_dist) 
          + (int) (2.0 * p_Img->frame_interval * (double) (1 + p_Img->gop_structure[enh_frame_to_code - 2].display_no)));
      }
      else
      {
        if (enh_frame_to_code == 1)
          p_Img->delta_pic_order_cnt[0] = p_Img->toppoc - 2*(p_Img->start_tr_gop  + (p_Img->frm_number - p_Img->lastIDRnumber)*(p_Img->base_dist));
        else
          p_Img->delta_pic_order_cnt[0] = p_Img->toppoc - 2*(p_Img->start_tr_gop  + (p_Img->frm_number - p_Img->lastIDRnumber - 1)*(p_Img->base_dist) 
          + (int) (2.0 * p_Img->frame_interval * (double) (1+ p_Img->gop_structure[enh_frame_to_code - 2].display_no)));
      }
    }
    else
    {
      p_Img->nal_reference_idc = (p_Inp->BRefPictures == 1 ) ? NALU_PRIORITY_LOW : NALU_PRIORITY_DISPOSABLE;
      set_poc(p_Img, p_Inp, (double)enh_frame_to_code);

      //the following is sent in the slice header
      if (p_Inp->BRefPictures != 1)
      {
        p_Img->delta_pic_order_cnt[0]= 2 * (enh_frame_to_code - 1);
      }
      else
      {
        p_Img->delta_pic_order_cnt[0]= -2;
      }
    }

    p_Img->delta_pic_order_cnt[1]= 0;

    if ((p_Inp->PicInterlace==FRAME_CODING)&&(p_Inp->MbInterlace==FRAME_CODING))
      p_Img->bottompoc = p_Img->toppoc;
    else
      p_Img->bottompoc = p_Img->toppoc + 1;

    p_Img->framepoc = imin (p_Img->toppoc, p_Img->bottompoc);
    p_Img->frm_no_in_file = CalculateFrameNumber(p_Img, p_Inp);
  }
}

/*!
 ***********************************************************************
 * \brief
 *    Encode a sequence
 ***********************************************************************
 */
static void encode_sequence(ImageParameters *p_Img, InputParameters *p_Inp)
{
  int HierarchicalCoding = p_Inp->HierarchicalCoding;
  int NumberBFrames = p_Inp->NumberBFrames;
  int jumpd = p_Inp->jumpd;
  int curr_frame_to_code = 0;
  int enh_frame_to_code = 0;

  for (curr_frame_to_code = 0; curr_frame_to_code < p_Inp->no_frames; curr_frame_to_code++)
  {     
    // Update frame_num counter
    if (p_Img->last_ref_idc == 1)
    {
      p_Img->frame_num++;
      p_Img->frame_num %= p_Img->max_frame_num;
    }
    
    // Read explicit sequence coding information
    if (p_Inp->ExplicitSeqCoding)
    {
      ExpFrameInfo *info = &p_Img->expSeq->info[curr_frame_to_code % p_Img->expSeq->no_frames];
      ReadExplicitSeqFile(p_Img->expSeq, p_Img->expSFile, curr_frame_to_code);
      ExplicitUpdateImgParams (info, p_Img, p_Inp);
      p_Img->b_frame_to_code = 0;
    }
    else
    {
      enh_frame_to_code = determine_coding_level(p_Img, p_Inp, curr_frame_to_code);
      p_Img->b_frame_to_code = enh_frame_to_code;

      if (enh_frame_to_code == 0) 
        prepare_first_layer(p_Img, p_Inp, curr_frame_to_code);
      else 
      {
        prepare_second_layer(p_Img, p_Inp, enh_frame_to_code);
      }
    }

    // redundant frame initialization and allocation
    if (p_Inp->redundant_pic_flag)
    {
      Init_redundant_frame(p_Img, p_Inp);
      Set_redundant_frame(p_Img, p_Inp);
    }

    encode_one_frame(p_Img, p_Inp); // encode one frame;

    p_Img->last_ref_idc = p_Img->nal_reference_idc ? 1 : 0;

    // if key frame is encoded, encode one redundant frame
    if (p_Inp->redundant_pic_flag && p_Img->key_frame)
    {
      encode_one_redundant_frame(p_Img, p_Inp);
    }

    if (p_Img->type == I_SLICE && p_Inp->EnableOpenGOP)
      p_Img->last_valid_reference = p_Img->ThisPOC;

    if (p_Img->currentPicture->idr_flag)
    {
      p_Img->idr_gop_number = 0;
    }
    else
      p_Img->idr_gop_number ++;

    if (p_Inp->ReportFrameStats)
      report_frame_statistic(p_Img, p_Inp);

  }

  p_Inp->HierarchicalCoding = HierarchicalCoding;
  p_Inp->NumberBFrames      = NumberBFrames;
  p_Inp->jumpd = jumpd;
}


/*!
 ***********************************************************************
 * \brief
 *    Free memory allocated for the encoder
 ***********************************************************************
 */

void free_encoder_memory(ImageParameters *p_Img, InputParameters *p_Inp)
{
  terminate_sequence(p_Img, p_Inp);

  flush_dpb(p_Img, p_Inp, &p_Inp->output);

  CloseFiles(&p_Inp->input_file1);

  if (-1 != p_Img->p_dec)
    close(p_Img->p_dec);
  if (p_Enc->p_trace)
    fclose(p_Enc->p_trace);

  Clear_Motion_Search_Module (p_Img, p_Inp);

  RandomIntraUninit(p_Img);
  FmoUninit(p_Img);

  if (p_Inp->HierarchicalCoding)
    clear_gop_structure (p_Img);

#ifdef _LEAKYBUCKET_
  calc_buffer(p_Img, p_Inp);
#endif

  // report everything
  report(p_Img, p_Inp, p_Img->p_Stats);

#ifdef _LEAKYBUCKET_
  if (p_Img->Bit_Buffer != NULL)
  {
    free(p_Img->Bit_Buffer);
    p_Img->Bit_Buffer = NULL;
  }
#endif

  free_dpb(p_Img, p_Inp, p_Img->p_Dpb);

  uninit_out_buffer(p_Img, p_Inp);

  free_global_buffers(p_Img, p_Inp);

  FreeParameterSets(p_Img);

   if (p_Inp->ExplicitSeqCoding)
     CloseExplicitSeqFile(p_Img);

  // free image mem
  free_img (p_Img, p_Inp);
}

/*!
 ***********************************************************************
 * \brief
 *    Initializes the POC structure with appropriate parameters.
 *
 ***********************************************************************
 */
static void init_poc(ImageParameters *p_Img, InputParameters *p_Inp)
{
  //the following should probably go in sequence parameters
  // frame poc's increase by 2, field poc's by 1

  p_Img->pic_order_cnt_type=p_Inp->pic_order_cnt_type;

  p_Img->delta_pic_order_always_zero_flag = FALSE;
  p_Img->num_ref_frames_in_pic_order_cnt_cycle= 1;

  if (p_Inp->BRefPictures == 1)
  {
    p_Img->offset_for_non_ref_pic  =  0;
    p_Img->offset_for_ref_frame[0] =  2;
  }
  else
  {
    p_Img->offset_for_non_ref_pic  = -2*(p_Inp->NumberBFrames);
    p_Img->offset_for_ref_frame[0] =  2*(p_Inp->NumberBFrames + 1);
  }

  if ((p_Inp->PicInterlace==FRAME_CODING) && (p_Inp->MbInterlace==FRAME_CODING))
  {
    p_Img->offset_for_top_to_bottom_field = 0;
    p_Img->bottom_field_pic_order_in_frame_present_flag = FALSE;
    p_Img->delta_pic_order_cnt_bottom = 0;
  }
  else
  {
    p_Img->offset_for_top_to_bottom_field = 1;
    p_Img->bottom_field_pic_order_in_frame_present_flag = TRUE;
    p_Img->delta_pic_order_cnt_bottom = 1;
  }
}


/*!
 ***********************************************************************
 * \brief
 *    Initializes the Image structure with appropriate parameters.
 * \par Input:
 *    Input Parameters InputParameters *inp
 * \par  Output:
 *    Image Parameters ImageParameters *p_Img
 ***********************************************************************
 */
static void init_img( ImageParameters *p_Img, InputParameters *p_Inp)
{
  int i, j;
  int imgpel_abs_range;

  p_Img->number         = -1;
  p_Img->start_frame_no = 0;
  p_Img->gop_number     = (p_Img->number - p_Img->start_frame_no);
  p_Img->start_tr_gop   = 0;

  p_Img->last_idr_number = 0;
  // Color format
  p_Img->yuv_format  = p_Inp->output.yuv_format;
  p_Img->P444_joined = (p_Img->yuv_format == YUV444 && !IS_INDEPENDENT(p_Inp));  

  //pel bitdepth init
  p_Img->bitdepth_luma            = (short) p_Inp->output.bit_depth[0];
  p_Img->bitdepth_scale[0]        = 1 << (p_Img->bitdepth_luma - 8);
  p_Img->bitdepth_lambda_scale    = 2 * (p_Img->bitdepth_luma - 8);
  p_Img->bitdepth_luma_qp_scale   = 3 *  p_Img->bitdepth_lambda_scale;
  p_Img->dc_pred_value_comp[0]    =  (imgpel) (1<<(p_Img->bitdepth_luma - 1));
  p_Img->max_pel_value_comp[0] = (1<<p_Img->bitdepth_luma) - 1;
  p_Img->max_imgpel_value_comp_sq[0] = p_Img->max_pel_value_comp[0] * p_Img->max_pel_value_comp[0];

  p_Img->dc_pred_value            = p_Img->dc_pred_value_comp[0]; // set defaults
  p_Img->max_imgpel_value         = (short) p_Img->max_pel_value_comp[0];
  p_Img->mb_size[0][0]            = p_Img->mb_size[0][1] = MB_BLOCK_SIZE;

  // Initialization for RC QP parameters (could be placed in ratectl.c)
  p_Img->RCMinQP                = p_Inp->RCMinQP[P_SLICE];
  p_Img->RCMaxQP                = p_Inp->RCMaxQP[P_SLICE];

  p_Img->WalkAround = 0;
  p_Img->NumberOfMBs = 0;

  // Set current residue & prediction array pointers

  if (p_Img->active_sps->profile_idc == BASELINE || p_Img->active_sps->profile_idc == MAIN || p_Img->active_sps->profile_idc == EXTENDED)
    p_Img->min_IPCM_value = 1;  // See Annex A for restriction in pcm sample values for pre FRExt profiles
  else
    p_Img->min_IPCM_value = 0;

  if (p_Img->yuv_format != YUV400)
  {
    p_Img->bitdepth_chroma             = (short) p_Inp->output.bit_depth[1];
    p_Img->bitdepth_scale[1]           = 1 << (p_Img->bitdepth_chroma - 8);
    p_Img->dc_pred_value_comp[1]       = (imgpel) (1<<(p_Img->bitdepth_chroma - 1));
    p_Img->dc_pred_value_comp[2]       = p_Img->dc_pred_value_comp[1];
    p_Img->max_pel_value_comp[1]       = (1<<p_Img->bitdepth_chroma) - 1;
    p_Img->max_pel_value_comp[2]       = p_Img->max_pel_value_comp[1];
    p_Img->max_imgpel_value_comp_sq[1] = p_Img->max_pel_value_comp[1] * p_Img->max_pel_value_comp[1];
    p_Img->max_imgpel_value_comp_sq[2] = p_Img->max_pel_value_comp[2] * p_Img->max_pel_value_comp[2];
    p_Img->num_blk8x8_uv               = (1<<p_Img->yuv_format)&(~(0x1));
    p_Img->num_cdc_coeff               = p_Img->num_blk8x8_uv << 1;

    p_Img->mb_size[1][0] = p_Img->mb_size[2][0] = p_Img->mb_cr_size_x = (p_Img->yuv_format == YUV420 || p_Img->yuv_format == YUV422) ? 8 : 16;
    p_Img->mb_size[1][1] = p_Img->mb_size[2][1] = p_Img->mb_cr_size_y = (p_Img->yuv_format == YUV444 || p_Img->yuv_format == YUV422) ? 16 : 8;

    p_Img->bitdepth_chroma_qp_scale = 6*(p_Img->bitdepth_chroma - 8);

    p_Img->chroma_qp_offset[0] = p_Img->active_pps->cb_qp_index_offset;
    p_Img->chroma_qp_offset[1] = p_Img->active_pps->cr_qp_index_offset;
  }
  else
  {
    p_Img->bitdepth_chroma     = 0;
    p_Img->bitdepth_scale[1]   = 0;
    p_Img->max_pel_value_comp[1] = 0;
    p_Img->max_pel_value_comp[2] = p_Img->max_pel_value_comp[1];
    p_Img->max_imgpel_value_comp_sq[1] = p_Img->max_pel_value_comp[1] * p_Img->max_pel_value_comp[1];
    p_Img->max_imgpel_value_comp_sq[2] = p_Img->max_pel_value_comp[2] * p_Img->max_pel_value_comp[2];
    p_Img->num_blk8x8_uv       = 0;
    p_Img->num_cdc_coeff       = 0;
    p_Img->mb_size[1][0] = p_Img->mb_size[2][0] = p_Img->mb_cr_size_x = 0;
    p_Img->mb_size[1][1] = p_Img->mb_size[2][1] = p_Img->mb_cr_size_y = 0;

    p_Img->bitdepth_chroma_qp_scale = 0;
    p_Img->bitdepth_chroma_qp_scale = 0;

    p_Img->chroma_qp_offset[0] = 0;
    p_Img->chroma_qp_offset[1] = 0;
  }  

  p_Img->max_bitCount =  128 + 256 * p_Img->bitdepth_luma + 2 * p_Img->mb_cr_size_y * p_Img->mb_cr_size_x * p_Img->bitdepth_chroma;
  //p_Img->max_bitCount =  (128 + 256 * p_Img->bitdepth_luma + 2 *p_Img->mb_cr_size_y * p_Img->mb_cr_size_x * p_Img->bitdepth_chroma)*2;

  p_Img->max_qp_delta = (25 + (p_Img->bitdepth_luma_qp_scale>>1));
  p_Img->min_qp_delta = p_Img->max_qp_delta + 1;

  p_Img->num_ref_frames = p_Img->active_sps->num_ref_frames;
  p_Img->max_num_references   = p_Img->active_sps->frame_mbs_only_flag ? p_Img->active_sps->num_ref_frames : 2 * p_Img->active_sps->num_ref_frames;

  p_Img->buf_cycle = p_Inp->num_ref_frames;
  p_Img->base_dist = p_Inp->jumpd + 1;  

  // Intra/IDR related parameters
  p_Img->lastIntraNumber = 0;
  p_Img->lastINTRA       = 0;
  p_Img->lastIDRnumber   = 0;
  p_Img->last_ref_idc    = 0;
  p_Img->idr_refresh     = 0;
  p_Img->idr_gop_number  = 0;
  p_Img->rewind_frame    = 0;

  p_Img->DeblockCall     = 0;
  p_Img->framerate       = (float) p_Inp->output.frame_rate;   // The basic frame rate (of the original sequence)

  if (p_Inp->AdaptiveRounding)
  {
    if (p_Img->yuv_format != 0)
    {
      get_mem4Dint(&(p_Img->ARCofAdj4x4), 3, MAXMODE, MB_BLOCK_SIZE, MB_BLOCK_SIZE); //all modes
      get_mem4Dint(&(p_Img->ARCofAdj8x8), p_Img->P444_joined ? 3 : 1, MAXMODE, MB_BLOCK_SIZE, MB_BLOCK_SIZE); //modes 0, 1, 2, 3, P8x8
    }     
    else
    {
      get_mem4Dint(&(p_Img->ARCofAdj4x4), 1, MAXMODE, MB_BLOCK_SIZE, MB_BLOCK_SIZE); //all modes
      get_mem4Dint(&(p_Img->ARCofAdj8x8), 1, MAXMODE, MB_BLOCK_SIZE, MB_BLOCK_SIZE); //modes 0, 1, 2, 3, P8x8
    }
  }

  imgpel_abs_range = (imax(p_Img->max_pel_value_comp[0], p_Img->max_pel_value_comp[1]) + 1) * 2;

  p_Img->width         = (p_Inp->output.width  + p_Img->auto_crop_right);
  p_Img->height        = (p_Inp->output.height + p_Img->auto_crop_bottom);
  p_Img->width_blk     = p_Img->width  / BLOCK_SIZE;
  p_Img->height_blk    = p_Img->height / BLOCK_SIZE;
  p_Img->width_padded  = p_Img->width  + 2 * IMG_PAD_SIZE;
  p_Img->height_padded = p_Img->height + 2 * IMG_PAD_SIZE;

  if (p_Img->yuv_format != YUV400)
  {
    p_Img->width_cr = p_Img->width  * mb_width_cr [p_Img->yuv_format] / 16;
    p_Img->height_cr= p_Img->height * mb_height_cr[p_Img->yuv_format] / 16;
  }
  else
  {
    p_Img->width_cr = 0;
    p_Img->height_cr= 0;
  }

  p_Img->height_cr_frame = p_Img->height_cr;

  p_Img->size = p_Img->width * p_Img->height;
  p_Img->size_cr = p_Img->width_cr * p_Img->height_cr;

  p_Img->PicWidthInMbs    = p_Img->width  / MB_BLOCK_SIZE;
  p_Img->FrameHeightInMbs = p_Img->height / MB_BLOCK_SIZE;
  p_Img->FrameSizeInMbs   = p_Img->PicWidthInMbs * p_Img->FrameHeightInMbs;

  p_Img->PicHeightInMapUnits = ( p_Img->active_sps->frame_mbs_only_flag ? p_Img->FrameHeightInMbs : p_Img->FrameHeightInMbs >> 1 );

  if ((p_Img->b8x8info = (Block8x8Info *) calloc(1, sizeof(Block8x8Info))) == NULL)
     no_mem_exit("init_img: p_Img->block8x8info");

  if( IS_INDEPENDENT(p_Inp) )
  {
    for( i = 0; i < MAX_PLANE; i++ ){
      if ((p_Img->mb_data_JV[i] = (Macroblock *) calloc(p_Img->FrameSizeInMbs,sizeof(Macroblock))) == NULL)
        no_mem_exit("init_img: p_Img->mb_data_JV");
    }
    p_Img->mb_data = NULL;
  }
  else
  {
    if ((p_Img->mb_data = (Macroblock *) calloc(p_Img->FrameSizeInMbs, sizeof(Macroblock))) == NULL)
      no_mem_exit("init_img: p_Img->mb_data");
  }

  if (p_Inp->UseConstrainedIntraPred)
  {
    if ((p_Img->intra_block = (int*)calloc(p_Img->FrameSizeInMbs, sizeof(int))) == NULL)
      no_mem_exit("init_img: p_Img->intra_block");
  }

  if (p_Inp->CtxAdptLagrangeMult == 1)
  {
    if ((p_Img->mb16x16_cost_frame = (double*)calloc(p_Img->FrameSizeInMbs, sizeof(double))) == NULL)
    {
      no_mem_exit("init p_Img->mb16x16_cost_frame");
    }
  }
  get_mem2D((byte***)&(p_Img->ipredmode), p_Img->height_blk, p_Img->width_blk);        //need two extra rows at right and bottom
  get_mem2D((byte***)&(p_Img->ipredmode8x8), p_Img->height_blk, p_Img->width_blk);     // help storage for ipredmode 8x8, inserted by YV
  memset(&(p_Img->ipredmode[0][0])   , -1, p_Img->height_blk * p_Img->width_blk *sizeof(char));
  memset(&(p_Img->ipredmode8x8[0][0]), -1, p_Img->height_blk * p_Img->width_blk *sizeof(char));


  // CAVLC mem
  get_mem3Dint(&(p_Img->nz_coeff), p_Img->FrameSizeInMbs, 4, 4+p_Img->num_blk8x8_uv);
  
  get_mem2Dolm     (&(p_Img->lambda)   , 10, 52 + p_Img->bitdepth_luma_qp_scale, p_Img->bitdepth_luma_qp_scale);
  get_mem2Dodouble (&(p_Img->lambda_md), 10, 52 + p_Img->bitdepth_luma_qp_scale, p_Img->bitdepth_luma_qp_scale);
  get_mem3Dodouble (&(p_Img->lambda_me), 10, 52 + p_Img->bitdepth_luma_qp_scale, 3, p_Img->bitdepth_luma_qp_scale);
  get_mem3Doint    (&(p_Img->lambda_mf), 10, 52 + p_Img->bitdepth_luma_qp_scale, 3, p_Img->bitdepth_luma_qp_scale);

  if (p_Inp->CtxAdptLagrangeMult == 1)
  {
    get_mem2Dodouble(&(p_Img->lambda_mf_factor), 10, 52 + p_Img->bitdepth_luma_qp_scale, p_Img->bitdepth_luma_qp_scale);
  }

  p_Img->b_frame_to_code = 0;
  p_Img->GopLevels = (p_Inp->NumberBFrames) ? 1 : 0;
  p_Img->mb_y_upd  = 0;

  RandomIntraInit (p_Img, p_Img->PicWidthInMbs, p_Img->FrameHeightInMbs, p_Inp->RandomIntraMBRefresh);

  InitSEIMessages(p_Img, p_Inp); 

  initInput(p_Img, &p_Inp->source, &p_Inp->output);

  // Allocate I/O Frame memory
  AllocateFrameMemory(p_Img, p_Inp, &p_Inp->source);

  // Initialize filtering parameters. If sending parameters, the offsets are
  // multiplied by 2 since inputs are taken in "div 2" format.
  // If not sending parameters, all fields are cleared
  if (p_Inp->DFSendParameters)
  {
    for (j = 0; j < 2; j++)
    {
      for (i = 0; i < NUM_SLICE_TYPES; i++)
      {
        p_Inp->DFAlpha[j][i] <<= 1;
        p_Inp->DFBeta [j][i] <<= 1;
      }
    }
  }
  else
  {
    for (j = 0; j < 2; j++)
    {
      for (i = 0; i < NUM_SLICE_TYPES; i++)
      {
        p_Inp->DFDisableIdc[j][i] = 0;
        p_Inp->DFAlpha     [j][i] = 0;
        p_Inp->DFBeta      [j][i] = 0;
      }
    }
  }

  p_Img->ChromaArrayType = p_Inp->separate_colour_plane_flag ? 0 : p_Inp->output.yuv_format;
  p_Img->colour_plane_id = 0;

  if (p_Inp->RDPictureDecision)
    p_Img->frm_iter = 3;
  else
    p_Img->frm_iter = 1;

  p_Img->frame_interval = (double) (p_Inp->frame_skip + 1);

  p_Img->max_frame_num = 1 << (p_Img->log2_max_frame_num_minus4 + 4);
  p_Img->max_pic_order_cnt_lsb = 1 << (p_Img->log2_max_pic_order_cnt_lsb_minus4 + 4);

  create_context_memory (p_Img, p_Inp);
}


/*!
 ***********************************************************************
 * \brief
 *    Free the Image structures
 * \par Input:
 *    Image Parameters ImageParameters *p_Img
 ***********************************************************************
 */
static void free_img (ImageParameters *p_Img, InputParameters *p_Inp)
{
  // Delete Frame memory 
  DeleteFrameMemory(p_Img);

  CloseSEIMessages(p_Img, p_Inp); 

  free_context_memory (p_Img);

  if (p_Inp->AdaptiveRounding)
  {
    free_mem4Dint(p_Img->ARCofAdj4x4);
    free_mem4Dint(p_Img->ARCofAdj8x8);
  }


  free (p_Img->p_SEI);
  free (p_Img->p_QScale);
  free (p_Img->p_Quant);
  free (p_Img->p_Dpb);
  free (p_Img->p_Stats);
  free (p_Img->p_Dist);
  free (p_Img);
}


/*!
 ***********************************************************************
 * \brief
 *    Free the Input structures
 * \par Input:
 *    Input Parameters InputParameters *p_Inp
 ***********************************************************************
 */
static void free_params (InputParameters *p_Inp)
{
  if ( p_Inp != NULL )
  {
    if ( p_Inp->top_left != NULL )
      free( p_Inp->top_left );
    if ( p_Inp->bottom_right != NULL )
      free( p_Inp->bottom_right );
    if ( p_Inp->slice_group_id != NULL )
      free( p_Inp->slice_group_id );
    if ( p_Inp->run_length_minus1 != NULL )
      free( p_Inp->run_length_minus1 );
    free( p_Inp );
  }
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
  if (pic != NULL)
  {
    free_slice_list(pic);
    free (pic);
  }
}


/*!
 ************************************************************************
 * \brief
 *    memory allocation for original picture buffers
 ************************************************************************
 */
int init_orig_buffers(ImageParameters *p_Img, InputParameters *p_Inp, ImageData *imgData)
{
  int memory_size = 0;
  int nplane;

  // allocate memory for reference frame buffers: imgData->frm_data
  imgData->format           = p_Inp->output;
  imgData->format.width     = p_Img->width;    
  imgData->format.height    = p_Img->height;
  imgData->format.width_cr  = p_Img->width_cr;
  imgData->format.height_cr = p_Img->height_cr;
  imgData->format.yuv_format = p_Img->yuv_format;
  imgData->format.auto_crop_bottom = p_Img->auto_crop_bottom;
  imgData->format.auto_crop_right  = p_Img->auto_crop_right;
  imgData->format.auto_crop_bottom_cr = (p_Img->auto_crop_bottom * mb_height_cr [p_Img->yuv_format]) / MB_BLOCK_SIZE;
  imgData->format.auto_crop_right_cr  = (p_Img->auto_crop_right * mb_width_cr [p_Img->yuv_format]) / MB_BLOCK_SIZE;

  if( IS_INDEPENDENT(p_Inp) )
  {

    for( nplane=0; nplane<MAX_PLANE; nplane++ )
    {
      memory_size += get_mem2Dpel(&(imgData->frm_data[nplane]), p_Img->height, p_Img->width);
    }
  }
  else
  {
    //imgData->format = p_Inp->input_file1.format;    

    memory_size += get_mem2Dpel(&(imgData->frm_data[0]), p_Img->height, p_Img->width);

    if (p_Img->yuv_format != YUV400)
    {
      int i, j, k;
      memory_size += get_mem2Dpel(&(imgData->frm_data[1]), p_Img->height_cr, p_Img->width_cr);
      memory_size += get_mem2Dpel(&(imgData->frm_data[2]), p_Img->height_cr, p_Img->width_cr);

      if (sizeof(imgpel) == sizeof(unsigned char))
      {
        for (k = 1; k < 3; k++)
          memset(&(imgData->frm_data[k][0][0]), 128, p_Img->height_cr * p_Img->width_cr * sizeof(imgpel));
      }
      else
      {
        for (k = 1; k < 3; k++)
          for (j = 0; j < p_Img->height_cr; j++)
            for (i = 0; i < p_Img->width_cr; i++)
              imgData->frm_data[k][j][i] = 128;
      }
    }
  }

  if (!p_Img->active_sps->frame_mbs_only_flag)
  {
    // allocate memory for field reference frame buffers
    memory_size += init_top_bot_planes(imgData->frm_data[0], p_Img->height, &(imgData->top_data[0]), &(imgData->bot_data[0]));

    if (p_Img->yuv_format != YUV400)
    {

      memory_size += 4*(sizeof(imgpel**));

      memory_size += init_top_bot_planes(imgData->frm_data[1], p_Img->height_cr, &(imgData->top_data[1]), &(imgData->bot_data[1]));
      memory_size += init_top_bot_planes(imgData->frm_data[2], p_Img->height_cr, &(imgData->top_data[2]), &(imgData->bot_data[2]));
    }
  }
  return memory_size;
}

/*!
 ************************************************************************
 * \brief
 *    Dynamic memory allocation of frame size related global buffers
 *    buffers are defined in global.h, allocated memory must be freed in
 *    void free_global_buffers()
 * \par Input:
 *    Input Parameters InputParameters *inp,                            \n
 *    Image Parameters ImageParameters *p_Img
 * \return Number of allocated bytes
 ************************************************************************
 */
static int init_global_buffers(ImageParameters *p_Img, InputParameters *p_Inp)
{
  int j, memory_size=0;

  if ((p_Img->enc_frame_picture = (StorablePicture**)malloc(6 * sizeof(StorablePicture*))) == NULL)
    no_mem_exit("init_global_buffers: *p_Img->enc_frame_picture");

  for (j = 0; j < 6; j++)
    p_Img->enc_frame_picture[j] = NULL;

  if ((p_Img->enc_field_picture = (StorablePicture**)malloc(2 * sizeof(StorablePicture*))) == NULL)
    no_mem_exit("init_global_buffers: *p_Img->enc_field_picture");

  for (j = 0; j < 2; j++)
    p_Img->enc_field_picture[j] = NULL;

  if ((p_Img->frame_pic = (Picture**)malloc(p_Img->frm_iter * sizeof(Picture*))) == NULL)
    no_mem_exit("init_global_buffers: *p_Img->frame_pic");

  for (j = 0; j < p_Img->frm_iter; j++)
    p_Img->frame_pic[j] = malloc_picture();

  if (p_Inp->si_frame_indicator || p_Inp->sp_periodicity)
  {
    p_Img->si_frame_indicator = FALSE; //indicates whether the frame is SP or SI
    p_Img->number_sp2_frames=0;

    p_Img->frame_pic_si = malloc_picture();//picture buffer for the encoded SI picture
    //allocation of lrec and p_Img->lrec_uv for SI picture
    get_mem2Dint (&p_Img->lrec, p_Img->height, p_Img->width);
    get_mem3Dint (&p_Img->lrec_uv, 2, p_Img->height, p_Img->width);
  }

  // Allocate memory for field picture coding
  if (p_Inp->PicInterlace != FRAME_CODING)
  { 
    if ((p_Img->field_pic = (Picture**)malloc(2 * sizeof(Picture*))) == NULL)
      no_mem_exit("init_global_buffers: *p_Img->field_pic");

    for (j = 0; j < 2; j++)
      p_Img->field_pic[j] = malloc_picture();
  }

  // Init memory data for input & encoded images
  memory_size += init_orig_buffers(p_Img, p_Inp, &p_Img->imgData);
  memory_size += init_orig_buffers(p_Img, p_Inp, &p_Img->imgData0);
  
  memory_size += get_mem2Dshort(&PicPos, p_Img->FrameSizeInMbs + 1, 2);

  for (j = 0; j < (int) p_Img->FrameSizeInMbs + 1; j++)
  {
    PicPos[j][0] = (short) (j % p_Img->PicWidthInMbs);
    PicPos[j][1] = (short) (j / p_Img->PicWidthInMbs);
  }


  if (p_Inp->rdopt == 3)
  {
    memory_size += allocate_errdo_mem(p_Img, p_Inp);
  }

  if (p_Inp->RestrictRef)
  {
    memory_size += get_mem2D(&p_Img->pixel_map,   p_Img->height,   p_Img->width);
    memory_size += get_mem2D(&p_Img->refresh_map, p_Img->height >> 3, p_Img->width >> 3);
  }

  if (!p_Img->active_sps->frame_mbs_only_flag)
  {
    memory_size += get_mem2Dpel(&p_Img->imgY_com, p_Img->height, p_Img->width);

    if (p_Img->yuv_format != YUV400)
    {
      memory_size += get_mem3Dpel(&p_Img->imgUV_com, 2, p_Img->height_cr, p_Img->width_cr);
    }
  }

  // allocate and set memory relating to motion estimation
  if (!p_Inp->IntraProfile)
  {  
    if (p_Inp->SearchMode == UM_HEX)
    {
      if ((p_Img->p_UMHex = (UMHexStruct*)calloc(1, sizeof(UMHexStruct))) == NULL)
        no_mem_exit("init_mv_block: p_Img->p_UMHex");
      memory_size += UMHEX_get_mem(p_Img, p_Inp);
    }
    else if (p_Inp->SearchMode == UM_HEX_SIMPLE)
    {
      if ((p_Img->p_UMHexSMP = (UMHexSMPStruct*)calloc(1, sizeof(UMHexSMPStruct))) == NULL)
        no_mem_exit("init_mv_block: p_Img->p_UMHexSMP");

      smpUMHEX_init(p_Img);
      memory_size += smpUMHEX_get_mem(p_Img);
    }
    else if (p_Inp->SearchMode == EPZS)
    {
      memory_size += EPZSInit(p_Img);
    }

  }

  if (p_Inp->RCEnable)
    rc_allocate_memory(p_Img, p_Inp);

  if (p_Inp->redundant_pic_flag)
  {
    memory_size += get_mem2Dpel(&p_Img->imgY_tmp, p_Img->height, p_Img->width);
    memory_size += get_mem2Dpel(&p_Img->imgUV_tmp[0], p_Img->height_cr, p_Img->width_cr);
    memory_size += get_mem2Dpel(&p_Img->imgUV_tmp[1], p_Img->height_cr, p_Img->width_cr);
  }

  memory_size += get_mem2Dint (&p_Img->imgY_sub_tmp, p_Img->height_padded, p_Img->width_padded);

  if ( p_Inp->ChromaMCBuffer )
    chroma_mc_setup(p_Img);

  p_Img->padded_size_x       = (p_Img->width + 2 * IMG_PAD_SIZE);
  p_Img->padded_size_x_m8x8  = (p_Img->padded_size_x - BLOCK_SIZE_8x8);
  p_Img->padded_size_x_m4x4  = (p_Img->padded_size_x - BLOCK_SIZE);
  p_Img->cr_padded_size_x    = (p_Img->width_cr + 2 * p_Img->pad_size_uv_x);
  p_Img->cr_padded_size_x2   = (p_Img->cr_padded_size_x << 1);
  p_Img->cr_padded_size_x4   = (p_Img->cr_padded_size_x << 2);
  p_Img->cr_padded_size_x_m8 = (p_Img->cr_padded_size_x - 8);

  // RGB images for distortion calculation
  // Recommended to do this allocation (and de-allocation) in 
  // the appropriate file instead of here.
  if(p_Inp->DistortionYUVtoRGB)
  {
    memory_size += create_RGB_memory(p_Img);
  }

  p_Img->pWPX = NULL;
  if ( p_Inp->WPMCPrecision )
  {
    wpxInitWPXObject(p_Img);
  }

  memory_size += InitProcessImage( p_Img, p_Inp );

  return memory_size;
}


/*!
 ************************************************************************
 * \brief
 *    Free allocated memory of original picture buffers
 ************************************************************************
 */
void free_orig_planes(ImageParameters *p_Img, InputParameters *p_Inp, ImageData *imgData)
{
  if( IS_INDEPENDENT(p_Inp) )
  {
    int nplane;
    for( nplane=0; nplane<MAX_PLANE; nplane++ )
    {
      free_mem2Dpel(imgData->frm_data[nplane]);      // free ref frame buffers
    }
  }
  else
  {
    free_mem2Dpel(imgData->frm_data[0]);      // free ref frame buffers
    
    if (imgData->format.yuv_format != YUV400)
    {
      free_mem2Dpel(imgData->frm_data[1]);
      free_mem2Dpel(imgData->frm_data[2]);
    }
  }

  if (!p_Img->active_sps->frame_mbs_only_flag)
  {
    free_top_bot_planes(imgData->top_data[0], imgData->bot_data[0]);

    if (imgData->format.yuv_format != YUV400)
    {
      free_top_bot_planes(imgData->top_data[1], imgData->bot_data[1]);
      free_top_bot_planes(imgData->top_data[2], imgData->bot_data[2]);
    }
  }
}


/*!
 ************************************************************************
 * \brief
 *    Free allocated memory of frame size related global buffers
 *    buffers are defined in global.h, allocated memory is allocated in
 *    int get_mem4global_buffers()
 * \par Input:
 *    Input Parameters InputParameters *inp,                             \n
 *    Image Parameters ImageParameters *p_Img
 * \par Output:
 *    none
 ************************************************************************
 */
static void free_global_buffers(ImageParameters *p_Img, InputParameters *p_Inp)
{
  int  i,j;

  if (p_Img->enc_frame_picture)
    free (p_Img->enc_frame_picture);
  if (p_Img->frame_pic)
  {
    for (j = 0; j < p_Img->frm_iter; j++)
    {
      if (p_Img->frame_pic[j])
        free_picture (p_Img->frame_pic[j]);
    }
    free (p_Img->frame_pic);
  }

  if (p_Img->enc_field_picture)
    free (p_Img->enc_field_picture);
  if (p_Img->field_pic)
  {
    for (j = 0; j < 2; j++)
    {
      if (p_Img->field_pic[j])
        free_picture (p_Img->field_pic[j]);
    }
    free (p_Img->field_pic);
  }

  // Deallocation of SI picture related memory
  if (p_Inp->si_frame_indicator || p_Inp->sp_periodicity)
  {
    free_picture (p_Img->frame_pic_si);
    //deallocation of lrec and p_Img->lrec_uv for SI frames
    free_mem2Dint (p_Img->lrec);
    free_mem3Dint (p_Img->lrec_uv);
  }

  free_orig_planes(p_Img, p_Inp, &p_Img->imgData);
  free_orig_planes(p_Img, p_Inp, &p_Img->imgData0);

  // free lookup memory which helps avoid divides with PicWidthInMbs
  free_mem2Dshort(PicPos);
  // Free Qmatrices and offsets
  free_QMatrix(p_Img->p_Quant);
  free_QOffsets(p_Img->p_Quant, p_Inp);


  if ( p_Inp->WPMCPrecision )
  {
    wpxFreeWPXObject(p_Img);
  }

  if (p_Img->imgY_sub_tmp) // free temp quarter pel frame buffers
  {
    free_mem2Dint (p_Img->imgY_sub_tmp);
    p_Img->imgY_sub_tmp = NULL;
  }

  // free mem, allocated in init_img()
  // free intra pred mode buffer for blocks
  free_mem2D((byte**)p_Img->ipredmode);
  free_mem2D((byte**)p_Img->ipredmode8x8);
  free(p_Img->b8x8info);
  if( IS_INDEPENDENT(p_Inp) )
  {
    for( i=0; i<MAX_PLANE; i++ ){
      free(p_Img->mb_data_JV[i]);
    }
  }
  else
  {
    free(p_Img->mb_data);
  }

  if(p_Inp->UseConstrainedIntraPred)
  {
    free (p_Img->intra_block);
  }

  if (p_Inp->CtxAdptLagrangeMult == 1)
  {
    free(p_Img->mb16x16_cost_frame);
  }

  if (p_Inp->rdopt == 3)
  {
    free_errdo_mem(p_Img);
  }

  if (p_Inp->RestrictRef)
  {
    free(p_Img->pixel_map[0]);
    free(p_Img->pixel_map);
    free(p_Img->refresh_map[0]);
    free(p_Img->refresh_map);
  }

  if (!p_Img->active_sps->frame_mbs_only_flag)
  {
    free_mem2Dpel(p_Img->imgY_com);
    
    if (p_Img->yuv_format != YUV400)
    {
      free_mem3Dpel(p_Img->imgUV_com);
    }
  }

  free_mem3Dint(p_Img->nz_coeff);

  free_mem2Dolm     (p_Img->lambda, p_Img->bitdepth_luma_qp_scale);
  free_mem2Dodouble (p_Img->lambda_md, p_Img->bitdepth_luma_qp_scale);
  free_mem3Dodouble (p_Img->lambda_me, 10, 52 + p_Img->bitdepth_luma_qp_scale, p_Img->bitdepth_luma_qp_scale);
  free_mem3Doint    (p_Img->lambda_mf, 10, 52 + p_Img->bitdepth_luma_qp_scale, p_Img->bitdepth_luma_qp_scale);

  if (p_Inp->CtxAdptLagrangeMult == 1)
  {
    free_mem2Dodouble(p_Img->lambda_mf_factor, p_Img->bitdepth_luma_qp_scale);
  }

  if (!p_Inp->IntraProfile)
  {
    if (p_Inp->SearchMode == UM_HEX)
    {
      UMHEX_free_mem(p_Img, p_Inp);
    }
    else if (p_Inp->SearchMode == UM_HEX_SIMPLE)
    {
      smpUMHEX_free_mem(p_Img);
    }
    else if (p_Inp->SearchMode == EPZS)
    {
      EPZSDelete(p_Img);
    }
  }

  if (p_Inp->RCEnable)
    rc_free_memory(p_Img, p_Inp);

  if (p_Inp->redundant_pic_flag)
  {
    free_mem2Dpel(p_Img->imgY_tmp);
    free_mem2Dpel(p_Img->imgUV_tmp[0]);
    free_mem2Dpel(p_Img->imgUV_tmp[1]);
  }

  // Again process should be moved into cconv_yuv2rgb.c file for cleanliness
  // These should not be globals but instead only be visible through that code.
  if(p_Inp->DistortionYUVtoRGB)
  {
    delete_RGB_memory(p_Img);
  }

  ClearProcessImage( p_Img, p_Inp );
}


/*!
 ************************************************************************
 * \brief
 *    Allocate memory for AC coefficients
 ************************************************************************
 */
int get_mem_ACcoeff (ImageParameters *p_Img, int***** cofAC)
{
  int num_blk8x8 = BLOCK_SIZE + p_Img->num_blk8x8_uv;
  
  get_mem4Dint(cofAC, num_blk8x8, BLOCK_SIZE, 2, 65);

  return num_blk8x8 * BLOCK_SIZE * 2 * 65 * sizeof(int);// 18->65 for ABT
}

/*!
 ************************************************************************
 * \brief
 *    Allocate memory for AC coefficients
 ************************************************************************
 */
int get_mem_ACcoeff_new (int****** cofAC, int chroma)
{ 
  get_mem5Dint(cofAC, BLOCK_SIZE, chroma, BLOCK_SIZE, 2, 65);
  return chroma * BLOCK_SIZE * BLOCK_SIZE * 2 * 65 * sizeof(int);// 18->65 for ABT
}

/*!
 ************************************************************************
 * \brief
 *    Allocate memory for DC coefficients
 ************************************************************************
 */
int get_mem_DCcoeff (int**** cofDC)
{
  get_mem3Dint(cofDC, 3, 2, 18);
  return 3 * 2 * 18 * sizeof(int); 
}


/*!
 ************************************************************************
 * \brief
 *    Free memory of AC coefficients
 ************************************************************************
 */
void free_mem_ACcoeff (int**** cofAC)
{
  free_mem4Dint(cofAC);
}

/*!
 ************************************************************************
 * \brief
 *    Free memory of AC coefficients
 ************************************************************************
 */
void free_mem_ACcoeff_new (int***** cofAC)
{
  free_mem5Dint(cofAC);
}

/*!
 ************************************************************************
 * \brief
 *    Free memory of DC coefficients
 ************************************************************************
 */
void free_mem_DCcoeff (int*** cofDC)
{
  free_mem3Dint(cofDC);
}

/*!
 ************************************************************************
 * \brief
 *    Sets indices to appropriate level constraints, depending on 
 *    current level_idc
 ************************************************************************
 */
static void SetLevelIndices(ImageParameters *p_Img)
{
  switch(p_Img->active_sps->level_idc)
  {
  case 9:
    p_Img->LevelIndex=1;
    break;
  case 10:
    p_Img->LevelIndex=0;
    break;
  case 11:
    if (!IS_FREXT_PROFILE(p_Img->active_sps->profile_idc) && (p_Img->active_sps->constrained_set3_flag == 0))
      p_Img->LevelIndex=2;
    else
      p_Img->LevelIndex=1;
    break;
  case 12:
    p_Img->LevelIndex=3;
    break;
  case 13:
    p_Img->LevelIndex=4;
    break;
  case 20:
    p_Img->LevelIndex=5;
    break;
  case 21:
    p_Img->LevelIndex=6;
    break;
  case 22:
    p_Img->LevelIndex=7;
    break;
  case 30:
    p_Img->LevelIndex=8;
    break;
  case 31:
    p_Img->LevelIndex=9;
    break;
  case 32:
    p_Img->LevelIndex=10;
    break;
  case 40:
    p_Img->LevelIndex=11;
    break;
  case 41:
    p_Img->LevelIndex=12;
    break;
  case 42:
    if (!IS_FREXT_PROFILE(p_Img->active_sps->profile_idc))
      p_Img->LevelIndex=13;
    else
      p_Img->LevelIndex=14;
    break;
  case 50:
    p_Img->LevelIndex=15;
    break;
  case 51:
    p_Img->LevelIndex=16;
    break;
  default:
    fprintf ( stderr, "Warning: unknown LevelIDC, using maximum level 5.1 \n" );
    p_Img->LevelIndex=16;
    break;
  }
}

/*!
 ************************************************************************
 * \brief
 *    initialize key frames and corresponding redundant frames.
 ************************************************************************
 */
void Init_redundant_frame(ImageParameters *p_Img, InputParameters *p_Inp)
{
  if (p_Inp->redundant_pic_flag)
  {
    if (p_Inp->NumberBFrames)
    {
      error("B frame not supported when redundant picture used!", 100);
    }

    if (p_Inp->PicInterlace)
    {
      error("Interlace not supported when redundant picture used!", 100);
    }

    if (p_Inp->num_ref_frames < p_Inp->PrimaryGOPLength)
    {
      error("NumberReferenceFrames must be no less than PrimaryGOPLength", 100);
    }

    if ((1<<p_Inp->NumRedundantHierarchy) > p_Inp->PrimaryGOPLength)
    {
      error("PrimaryGOPLength must be greater than 2^NumRedundantHeirarchy", 100);
    }

    if (p_Inp->Verbose != 1)
    {
      error("Redundant slices not supported when Verbose != 1", 100);
    }
  }

  p_Img->key_frame = 0;
  p_Img->redundant_coding = 0;
  p_Img->redundant_pic_cnt = 0;
  p_Img->frameNuminGOP = p_Img->frm_number % p_Inp->PrimaryGOPLength;
  if (p_Img->frm_number == 0)
  {
    p_Img->frameNuminGOP = -1;
  }
}

/*!
 ************************************************************************
 * \brief
 *    allocate redundant frames in a primary GOP.
 ************************************************************************
 */
void Set_redundant_frame(ImageParameters *p_Img, InputParameters *p_Inp)
{
  int GOPlength = p_Inp->PrimaryGOPLength;

  //start frame of GOP
  if (p_Img->frameNuminGOP == 0)
  {
    p_Img->redundant_coding = 0;
    p_Img->key_frame = 1;
    p_Img->redundant_ref_idx = GOPlength;
  }

  //1/2 position
  if (p_Inp->NumRedundantHierarchy > 0)
  {
    if (p_Img->frameNuminGOP == GOPlength >> 1)
    {
      p_Img->redundant_coding = 0;
      p_Img->key_frame = 1;
      p_Img->redundant_ref_idx = GOPlength >> 1;
    }
  }

  //1/4, 3/4 position
  if (p_Inp->NumRedundantHierarchy > 1)
  {
    if (p_Img->frameNuminGOP == (GOPlength >> 2) || p_Img->frameNuminGOP == ((GOPlength*3) >> 2))
    {
      p_Img->redundant_coding = 0;
      p_Img->key_frame = 1;
      p_Img->redundant_ref_idx = GOPlength >> 2;
    }
  }

  //1/8, 3/8, 5/8, 7/8 position
  if (p_Inp->NumRedundantHierarchy > 2)
  {
    if (p_Img->frameNuminGOP == GOPlength >> 3 || p_Img->frameNuminGOP == ((GOPlength*3) >> 3)
      || p_Img->frameNuminGOP == ((GOPlength*5) >> 3) || p_Img->frameNuminGOP == ((GOPlength*7) & 0x03))
    {
      p_Img->redundant_coding = 0;
      p_Img->key_frame = 1;
      p_Img->redundant_ref_idx = GOPlength >> 3;
    }
  }

  //1/16, 3/16, 5/16, 7/16, 9/16, 11/16, 13/16 position
  if (p_Inp->NumRedundantHierarchy > 3)
  {
    if (p_Img->frameNuminGOP == (GOPlength >> 4) || p_Img->frameNuminGOP == ((GOPlength*3) >> 4)
      || p_Img->frameNuminGOP == ((GOPlength*5) >> 4) || p_Img->frameNuminGOP == ((GOPlength*7) >> 4)
      || p_Img->frameNuminGOP == ((GOPlength*9) >> 4) || p_Img->frameNuminGOP == ((GOPlength*11) >> 4)
      || p_Img->frameNuminGOP == ((GOPlength*13) >> 4))
    {
      p_Img->redundant_coding = 0;
      p_Img->key_frame = 1;
      p_Img->redundant_ref_idx = GOPlength >> 4;
    }
  }
}

/*!
 ************************************************************************
 * \brief
 *    encode one redundant frame.
 ************************************************************************
 */
void encode_one_redundant_frame(ImageParameters *p_Img, InputParameters *p_Inp)
{
  p_Img->key_frame = 0;
  p_Img->redundant_coding = 1;
  p_Img->redundant_pic_cnt = 1;

  if (!p_Img->currentPicture->idr_flag)
  {
    if (p_Img->type == I_SLICE)
    {
      set_slice_type( p_Img, p_Inp, P_SLICE );
    }
  }

  encode_one_frame(p_Img, p_Inp);
}

/*!
 ************************************************************************
 * \brief
 *    Setup Chroma MC Variables
 ************************************************************************
 */
static void chroma_mc_setup(ImageParameters *p_Img)
{
  // initialize global variables used for chroma interpolation and buffering
  if ( p_Img->yuv_format == YUV420 )
  {
    p_Img->pad_size_uv_x = IMG_PAD_SIZE >> 1;
    p_Img->pad_size_uv_y = IMG_PAD_SIZE >> 1;
    p_Img->chroma_mask_mv_y = 7;
    p_Img->chroma_mask_mv_x = 7;
    p_Img->chroma_shift_x = 3;
    p_Img->chroma_shift_y = 3;
  }
  else if ( p_Img->yuv_format == YUV422 )
  {
    p_Img->pad_size_uv_x = IMG_PAD_SIZE >> 1;
    p_Img->pad_size_uv_y = IMG_PAD_SIZE;
    p_Img->chroma_mask_mv_y = 3;
    p_Img->chroma_mask_mv_x = 7;
    p_Img->chroma_shift_y = 2;
    p_Img->chroma_shift_x = 3;
  }
  else
  { // YUV444
    p_Img->pad_size_uv_x = IMG_PAD_SIZE;
    p_Img->pad_size_uv_y = IMG_PAD_SIZE;
    p_Img->chroma_mask_mv_y = 3;
    p_Img->chroma_mask_mv_x = 3;
    p_Img->chroma_shift_y = 2;
    p_Img->chroma_shift_x = 2;
  }
  p_Img->shift_cr_y  = p_Img->chroma_shift_y - 2;
  p_Img->shift_cr_x  = p_Img->chroma_shift_x - 2;
}

