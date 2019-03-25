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
 *     JM 16.1 (FRExt)
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
 *   - Athanasios Leontaris            <aleon@dolby.com>
 ***********************************************************************
 */

#include "contributors.h"

#include <time.h>
#include <math.h>

#include "global.h"
#include "cconv_yuv2rgb.h"
#include "configfile.h"
#include "conformance.h"
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
#include "pred_struct.h"

static const int mb_width_cr[4] = {0, 8, 8,16};
static const int mb_height_cr[4]= {0, 8,16,16};

EncoderParams   *p_Enc = NULL;

static void SetLevelIndices     (VideoParameters *p_Vid);
static void chroma_mc_setup     (VideoParameters *p_Vid);
static void init_img            (VideoParameters *p_Vid);
static void init_encoder        (VideoParameters *p_Vid, InputParameters *p_Inp);
static int  init_global_buffers (VideoParameters *p_Vid, InputParameters *p_Inp);
static void free_global_buffers (VideoParameters *p_Vid, InputParameters *p_Inp);
static void free_img            (VideoParameters *p_Vid, InputParameters *p_Inp);
static void free_params         (InputParameters *p_Inp);

static void encode_sequence     (VideoParameters *p_Vid, InputParameters *p_Inp);

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
 *    Allocate the Image structure
 * \par  Output:
 *    Image Parameters VideoParameters *p_Vid
 ***********************************************************************
 */
static void alloc_img( VideoParameters **p_Vid)
{
  if ((*p_Vid = (VideoParameters *) calloc(1, sizeof(VideoParameters)))==NULL) 
    no_mem_exit("alloc_img: p_Vid");
  if ((((*p_Vid)->p_Dist)  = (DistortionParams *) calloc(1, sizeof(DistortionParams)))==NULL) 
    no_mem_exit("alloc_img: p_Dist");
  if ((((*p_Vid)->p_Stats) = (StatParameters *) calloc(1, sizeof(StatParameters)))==NULL) 
    no_mem_exit("alloc_img: p_Stats");
  if (((*p_Vid)->p_Dpb     = (DecodedPictureBuffer *) calloc(1, sizeof(DecodedPictureBuffer)))==NULL) 
    no_mem_exit("alloc_img: p_Dpb");
  if ((((*p_Vid)->p_Quant)  = (QuantParameters *) calloc(1, sizeof(QuantParameters)))==NULL) 
    no_mem_exit("alloc_img: p_Quant");
  if ((((*p_Vid)->p_QScale)  = (ScaleParameters *) calloc(1, sizeof(ScaleParameters)))==NULL) 
    no_mem_exit("alloc_img: p_QScale");
  if ((((*p_Vid)->p_SEI)  = (SEIParameters *) calloc(1, sizeof(SEIParameters)))==NULL) 
    no_mem_exit("alloc_img: p_SEI");


  (*p_Vid)->p_dec = -1;
  (*p_Vid)->p_log = NULL;
  (*p_Vid)->f_annexb = NULL;
  // Init rtp related info
  (*p_Vid)->f_rtp = NULL;
  (*p_Vid)->CurrentRTPTimestamp = 0;         
  (*p_Vid)->CurrentRTPSequenceNumber = 0;
}

/*!
 ***********************************************************************
 * \brief
 *    Allocate the Input structure
 * \par  Output:
 *    Input Parameters InputParameters *p_Vid
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

  alloc_img(&((*p_Enc)->p_Vid));
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

  Configure (p_Enc->p_Vid, p_Enc->p_Inp, argc, argv);

  // init encoder
  init_encoder(p_Enc->p_Vid, p_Enc->p_Inp);

  // encode sequence
  encode_sequence(p_Enc->p_Vid, p_Enc->p_Inp);

  // terminate sequence
  free_encoder_memory(p_Enc->p_Vid, p_Enc->p_Inp);

  free_params (p_Enc->p_Inp);  
  free_encoder(p_Enc);

  return 0;
}


/*!
 ************************************************************************
 * \brief
 *    calculate Ceil(Log2(uiVal))
 ************************************************************************
 */
unsigned CeilLog2( unsigned uiVal)
{
  unsigned uiTmp = uiVal-1;
  unsigned uiRet = 0;

  while( uiTmp != 0 )
  {
    uiTmp >>= 1;
    uiRet++;
  }
  return uiRet;
}

/*!
 ***********************************************************************
 * \brief
 *    Initialize encoder
 ***********************************************************************
 */

static void init_encoder(VideoParameters *p_Vid, InputParameters *p_Inp)
{
  p_Vid->p_Inp = p_Inp;
  p_Vid->giRDOpt_B8OnlyFlag = FALSE;
  p_Vid->p_log = NULL;

  p_Vid->cabac_encoding = 0;
  p_Vid->frame_statistic_start = 1;

  if (p_Inp->Log2MaxFNumMinus4 == -1)
  {    
    p_Vid->log2_max_frame_num_minus4 = iClip3(0,12, (int) (CeilLog2(p_Inp->no_frames) - 4)); // hack for now...
  }
  else  
    p_Vid->log2_max_frame_num_minus4 = p_Inp->Log2MaxFNumMinus4;

  if (p_Vid->log2_max_frame_num_minus4 == 0 && p_Inp->num_ref_frames == 16)
  {
    snprintf(errortext, ET_SIZE, " NumberReferenceFrames=%d and Log2MaxFNumMinus4=%d may lead to an invalid value of frame_num.", p_Inp->num_ref_frames, p_Inp-> Log2MaxFNumMinus4);
    error (errortext, 500);
  }

  // set proper p_Vid->log2_max_pic_order_cnt_lsb_minus4.
  if (p_Inp->Log2MaxPOCLsbMinus4 == - 1)
    p_Vid->log2_max_pic_order_cnt_lsb_minus4 = iClip3(0,12, (int) (CeilLog2( p_Inp->no_frames << 1 ) - 4)); // hack for now
  else
    p_Vid->log2_max_pic_order_cnt_lsb_minus4 = p_Inp->Log2MaxPOCLsbMinus4;

  if (((1<<(p_Vid->log2_max_pic_order_cnt_lsb_minus4 + 3)) < p_Inp->jumpd * 4) && p_Inp->Log2MaxPOCLsbMinus4 != -1)
    error("log2_max_pic_order_cnt_lsb_minus4 might not be sufficient for encoding. Increase value.",400);

  if (strlen (p_Inp->ReconFile) > 0 && (p_Vid->p_dec = open(p_Inp->ReconFile, OPENFLAGS_WRITE, OPEN_PERMISSIONS))==-1)
  {
    snprintf(errortext, ET_SIZE, "Error open file %s", p_Inp->ReconFile);
    error (errortext, 500);
  }

  if ((p_Inp->output.width & 0x0F) != 0)
  {
    p_Vid->auto_crop_right = 16 - (p_Inp->output.width & 0x0F);
  }
  else
  {
    p_Vid->auto_crop_right = 0;
  }

  if (p_Inp->PicInterlace || p_Inp->MbInterlace)
  {
    if ((p_Inp->output.height & 0x01) != 0)
    {
      error ("even number of lines required for interlaced coding", 500);
    }
    
    if ((p_Inp->output.height & 0x1F) != 0)
    {
      p_Vid->auto_crop_bottom = 32 - (p_Inp->output.height & 0x1F);
    }
    else
    {
      p_Vid->auto_crop_bottom=0;
    }
  }
  else
  {
    if ((p_Inp->output.height & 0x0F) != 0)
    {
      p_Vid->auto_crop_bottom = 16 - (p_Inp->output.height & 0x0F);
    }
    else
    {
      p_Vid->auto_crop_bottom = 0;
    }
  }
  if (p_Vid->auto_crop_bottom || p_Vid->auto_crop_right)
  {
    fprintf (stderr, "Warning: Automatic cropping activated: Coded frame Size: %dx%d\n", 
      p_Inp->output.width + p_Vid->auto_crop_right, p_Inp->output.height + p_Vid->auto_crop_bottom);
  }

  // read the slice group configuration file. Only for types 0, 2 or 6
  if ( 0 != p_Inp->num_slice_groups_minus1 )
  {
    read_slice_group_info(p_Vid, p_Inp);
  }
  //! the number of slice groups is forced to be 1 for slice group type 3-5
  if(p_Inp->num_slice_groups_minus1 > 0)
  {
    if( (p_Inp->slice_group_map_type >= 3) && (p_Inp->slice_group_map_type<=5) )
      p_Inp->num_slice_groups_minus1 = 1;
  }

  if(p_Inp->RCEnable)
  {
    if (p_Inp->basicunit == 0)
      p_Inp->basicunit = (p_Inp->output.height + p_Vid->auto_crop_bottom)*(p_Inp->output.width + p_Vid->auto_crop_right)/256;

    if ( ((p_Inp->output.height + p_Vid->auto_crop_bottom)*(p_Inp->output.width + p_Vid->auto_crop_right)/256) % p_Inp->basicunit != 0)
    {
      snprintf(errortext, ET_SIZE, "Frame size in macroblocks must be a multiple of BasicUnit.");
      error (errortext, 500);
    }
  }

  LevelCheck(p_Vid, p_Inp);

  // Open Files
  OpenFiles(&p_Inp->input_file1);


  Init_QMatrix(p_Vid, p_Inp);
  Init_QOffsetMatrix(p_Vid);

  init_poc(p_Vid);
  GenerateParameterSets(p_Vid);
  SetLevelIndices(p_Vid);

  init_img  (p_Vid);

  if (p_Inp->rdopt == 3)
  {
    init_error_conceal(p_Vid,p_Inp->ErrorConcealment); 
  }

#ifdef _LEAKYBUCKET_
  p_Vid->initial_Bframes = 0;
  p_Vid->Bit_Buffer = (long *)malloc((p_Inp->no_frames + 1) * sizeof(long));
  p_Vid->total_frame_buffer = 0;
#endif

  // Prepare hierarchical coding structures. 
  // Code could be extended in the future to allow structure adaptation.
  if (p_Inp->NumberBFrames && p_Inp->HierarchicalCoding == 3)
  {
    init_gop_structure(p_Vid, p_Inp);
    interpret_gop_structure(p_Vid, p_Inp);
  }

  p_Vid->p_Dpb->init_done = 0;

  init_dpb(p_Vid, p_Vid->p_Dpb);
  init_out_buffer(p_Vid);
  init_stats (p_Inp, p_Vid->p_Stats);
  init_dstats(p_Vid->p_Dist);

  p_Vid->enc_picture = NULL;

  init_global_buffers(p_Vid, p_Inp);

  if ( p_Inp->WPMCPrecision )
  {
    wpxInitWPXPasses(p_Vid, p_Inp);
  }

  Init_Motion_Search_Module (p_Vid, p_Inp);
  information_init(p_Vid, p_Inp, p_Vid->p_Stats);

  if(p_Inp->DistortionYUVtoRGB)
    init_YUVtoRGB(p_Vid, p_Inp);

  //Rate control
  if (p_Inp->RCEnable)
    rc_init_sequence(p_Vid, p_Inp);

  p_Vid->last_valid_reference = 0;
  p_Vid->tot_time = 0;                 // time for total encoding session
  p_Vid->last_bit_ctr_n = 0;

  p_Vid->initial_Bframes = p_Inp->NumberBFrames;  

  p_Vid->type = I_SLICE;
  // Write sequence header (with parameter sets)
  p_Vid->p_Stats->bit_ctr_filler_data = 0;
  p_Vid->p_Stats->bit_ctr_filler_data_n = 0;
  p_Vid->p_Stats->bit_ctr_parametersets = 0;
  p_Vid->p_Stats->bit_slice = start_sequence(p_Vid, p_Inp);

  if (p_Inp->UseRDOQuant)
    precalculate_unary_exp_golomb_level(p_Vid);

  if (p_Inp->ExplicitSeqCoding)
    OpenExplicitSeqFile(p_Vid, p_Inp);

  if ( p_Inp->ChromaMCBuffer )
    p_Vid->OneComponentChromaPrediction4x4 = OneComponentChromaPrediction4x4_retrieve;
  else
    p_Vid->OneComponentChromaPrediction4x4 = OneComponentChromaPrediction4x4_regenerate;

  p_Vid->searchRange.min_x = -p_Inp->search_range << 2;
  p_Vid->searchRange.max_x =  p_Inp->search_range << 2;
  p_Vid->searchRange.min_y = -p_Inp->search_range << 2;
  p_Vid->searchRange.max_y =  p_Inp->search_range << 2;
}

/*!
************************************************************************
* \brief
*    Prepare parameters for the current frame
************************************************************************
*/
static void prepare_frame_params(VideoParameters *p_Vid, InputParameters *p_Inp, int curr_frame_to_code)
{
  SeqStructure *p_seq_struct = p_Vid->p_pred;
  FrameUnitStruct *p_cur_frm = p_Vid->p_curr_frm_struct;

  if ( p_Inp->ExplicitSeqCoding )
  {
    ExpFrameInfo *info = &p_Vid->expSeq->info[curr_frame_to_code % p_Vid->expSeq->no_frames];
    ReadExplicitSeqFile(p_Vid->expSeq, p_Vid->expSFile, curr_frame_to_code);
    // override and overwrite whatever was in p_cur_frm
    populate_frame_explicit( info, p_Inp, p_cur_frm, p_seq_struct->max_num_slices );
  }

  // populate coding parameters for the current frame (frame_no, slice_type, nal_ref_idc, poc, reset frame_num)
  p_Vid->frame_no = p_cur_frm->frame_no;
  p_Vid->frm_no_in_file = (1 + p_Inp->frame_skip) * p_Vid->frame_no;
  if ( p_cur_frm->type == SI_SLICE )
  {
    set_slice_type( p_Vid, p_Inp, SP_SLICE );
  }
  else
  {
    set_slice_type( p_Vid, p_Inp, p_cur_frm->type ); // dominant slice type; perhaps it is better to use the type of the first slice?
  }
  p_Vid->nal_reference_idc = p_cur_frm->nal_ref_idc;
  switch( p_Vid->pic_order_cnt_type )
  {
  case 0:
    get_poc_type_zero( p_Vid, p_Inp, p_cur_frm );
    break;
  case 1:
    get_poc_type_one( p_Vid, p_Inp, p_cur_frm );
    break;
  default:
  case 2:
    get_poc_type_zero( p_Vid, p_Inp, p_cur_frm );
    break;
  }
  p_Vid->frame_num = ( p_cur_frm->idr_flag == 1 ) ? 0 : p_Vid->frame_num;
  
  //Rate control
  if (p_Inp->RCEnable && p_Vid->type == I_SLICE)
    rc_init_gop_params(p_Vid, p_Inp);

  // which layer does the image belong to?
  p_Vid->layer = ((p_Vid->curr_frm_idx - p_Vid->last_idr_code_order) % (p_Inp->NumFramesInELSubSeq + 1)) ? 0 : 1;  
}

/*!
 ***********************************************************************
 * \brief
 *    Encode a sequence
 ***********************************************************************
 */
static void encode_sequence(VideoParameters *p_Vid, InputParameters *p_Inp)
{
  int curr_frame_to_code;
  int frames_to_code = p_Inp->no_frames;
  int frame_num_bak = 0, frame_coded;
  SeqStructure *p_seq_struct = p_Vid->p_pred;
  FrameUnitStruct *p_frm = p_seq_struct->p_frm;

  for (curr_frame_to_code = 0; curr_frame_to_code < frames_to_code; curr_frame_to_code++)
  {     
    // determine whether to populate additional frames in the prediction structure
    if ( curr_frame_to_code >= p_Vid->p_pred->pop_start_frame )
    {
      populate_frm_struct( p_Vid, p_Inp, p_Vid->p_pred, p_Inp->FrmStructBufferLength, frames_to_code );
    }
    p_Vid->curr_frm_idx = curr_frame_to_code;
    p_Vid->p_curr_frm_struct = p_frm + ( p_Vid->curr_frm_idx % p_Vid->frm_struct_buffer ); // pointer to current frame structure
    p_Vid->number = curr_frame_to_code;

    if ( p_Vid->p_curr_frm_struct->frame_no >= p_Inp->no_frames )
    {
      continue;
    }

    // Update frame_num counter
    frame_num_bak = p_Vid->frame_num;
    if (p_Vid->last_ref_idc == 1)
    {      
      p_Vid->frame_num++;
      p_Vid->frame_num %= p_Vid->max_frame_num;
    }

    prepare_frame_params(p_Vid, p_Inp, curr_frame_to_code);

    // redundant frame initialization and allocation
    if (p_Inp->redundant_pic_flag)
    {
      Init_redundant_frame(p_Vid, p_Inp);
      Set_redundant_frame(p_Vid, p_Inp);
    }

    frame_coded = encode_one_frame(p_Vid, p_Inp); // encode one frame;
    if ( !frame_coded )
    {
      p_Vid->frame_num = frame_num_bak;
      continue;
    }

    p_Vid->last_ref_idc = p_Vid->nal_reference_idc ? 1 : 0;

    // if key frame is encoded, encode one redundant frame
    if (p_Inp->redundant_pic_flag && p_Vid->key_frame)
    {
      encode_one_redundant_frame(p_Vid, p_Inp);
    }

    if (p_Vid->type == I_SLICE && p_Inp->EnableOpenGOP)
      p_Vid->last_valid_reference = p_Vid->ThisPOC;

    if (p_Inp->ReportFrameStats)
      report_frame_statistic(p_Vid, p_Inp);

  }
}


/*!
 ***********************************************************************
 * \brief
 *    Free memory allocated for the encoder
 ***********************************************************************
 */

void free_encoder_memory(VideoParameters *p_Vid, InputParameters *p_Inp)
{
  terminate_sequence(p_Vid, p_Inp);

  flush_dpb(p_Vid, &p_Inp->output);

  CloseFiles(&p_Inp->input_file1);

  if (-1 != p_Vid->p_dec)
    close(p_Vid->p_dec);
  
  if (p_Enc->p_trace)
    fclose(p_Enc->p_trace);

  Clear_Motion_Search_Module (p_Vid, p_Inp);

  RandomIntraUninit(p_Vid);
  FmoUninit(p_Vid);

  if (p_Inp->NumberBFrames && p_Inp->HierarchicalCoding == 3)
  {
    clear_gop_structure (p_Vid);
  }

#ifdef _LEAKYBUCKET_
  calc_buffer(p_Vid, p_Inp);
#endif

  // report everything
  report(p_Vid, p_Inp, p_Vid->p_Stats);

#ifdef _LEAKYBUCKET_
  if (p_Vid->Bit_Buffer != NULL)
  {
    free(p_Vid->Bit_Buffer);
    p_Vid->Bit_Buffer = NULL;
  }
#endif

  free_dpb(p_Vid, p_Vid->p_Dpb);

  uninit_out_buffer(p_Vid);

  free_global_buffers(p_Vid, p_Inp);

  FreeParameterSets(p_Vid);

  if (p_Inp->ExplicitSeqCoding)
    CloseExplicitSeqFile(p_Vid);

  // free image mem
  free_img (p_Vid, p_Inp);
}

/*!
 ***********************************************************************
 * \brief
 *    Initializes the Image structure with appropriate parameters.
 * \par Input:
 *    Input Parameters InputParameters *inp
 * \par  Output:
 *    Image Parameters VideoParameters *p_Vid
 ***********************************************************************
 */
static void init_img( VideoParameters *p_Vid)
{
  InputParameters *p_Inp = p_Vid->p_Inp;
  int i, j;
  int imgpel_abs_range;

  p_Vid->number         = -1;

  p_Vid->last_idr_code_order = 0;
  p_Vid->last_idr_disp_order = 0;
  p_Vid->last_mmco_5_code_order = -1;
  p_Vid->last_mmco_5_disp_order = -1;
  // Color format
  p_Vid->yuv_format  = p_Inp->output.yuv_format;
  p_Vid->P444_joined = (p_Vid->yuv_format == YUV444 && !IS_INDEPENDENT(p_Inp));  

  //pel bitdepth init
  p_Vid->bitdepth_luma            = (short) p_Inp->output.bit_depth[0];
  p_Vid->bitdepth_scale[0]        = 1 << (p_Vid->bitdepth_luma - 8);
  p_Vid->bitdepth_lambda_scale    = 2 * (p_Vid->bitdepth_luma - 8);
  p_Vid->bitdepth_luma_qp_scale   = 3 *  p_Vid->bitdepth_lambda_scale;
  p_Vid->dc_pred_value_comp[0]    =  (imgpel) (1<<(p_Vid->bitdepth_luma - 1));
  p_Vid->max_pel_value_comp[0] = (1<<p_Vid->bitdepth_luma) - 1;
  p_Vid->max_imgpel_value_comp_sq[0] = p_Vid->max_pel_value_comp[0] * p_Vid->max_pel_value_comp[0];

  p_Vid->dc_pred_value            = p_Vid->dc_pred_value_comp[0]; // set defaults
  p_Vid->max_imgpel_value         = (short) p_Vid->max_pel_value_comp[0];
  p_Vid->mb_size[0][0]            = p_Vid->mb_size[0][1] = MB_BLOCK_SIZE;

  // Initialization for RC QP parameters (could be placed in ratectl.c)
  p_Vid->RCMinQP                = p_Inp->RCMinQP[P_SLICE];
  p_Vid->RCMaxQP                = p_Inp->RCMaxQP[P_SLICE];

  p_Vid->WalkAround = 0;
  p_Vid->NumberOfMBs = 0;

  // Set current residue & prediction array pointers

  if (p_Vid->active_sps->profile_idc == BASELINE || p_Vid->active_sps->profile_idc == MAIN || p_Vid->active_sps->profile_idc == EXTENDED)
    p_Vid->min_IPCM_value = 1;  // See Annex A for restriction in pcm sample values for pre FRExt profiles
  else
    p_Vid->min_IPCM_value = 0;

  if (p_Vid->yuv_format != YUV400)
  {
    p_Vid->bitdepth_chroma             = (short) p_Inp->output.bit_depth[1];
    p_Vid->bitdepth_scale[1]           = 1 << (p_Vid->bitdepth_chroma - 8);
    p_Vid->dc_pred_value_comp[1]       = (imgpel) (1<<(p_Vid->bitdepth_chroma - 1));
    p_Vid->dc_pred_value_comp[2]       = p_Vid->dc_pred_value_comp[1];
    p_Vid->max_pel_value_comp[1]       = (1<<p_Vid->bitdepth_chroma) - 1;
    p_Vid->max_pel_value_comp[2]       = p_Vid->max_pel_value_comp[1];
    p_Vid->max_imgpel_value_comp_sq[1] = p_Vid->max_pel_value_comp[1] * p_Vid->max_pel_value_comp[1];
    p_Vid->max_imgpel_value_comp_sq[2] = p_Vid->max_pel_value_comp[2] * p_Vid->max_pel_value_comp[2];
    p_Vid->num_blk8x8_uv               = (1<<p_Vid->yuv_format)&(~(0x1));
    p_Vid->num_cdc_coeff               = p_Vid->num_blk8x8_uv << 1;

    p_Vid->mb_size[1][0] = p_Vid->mb_size[2][0] = p_Vid->mb_cr_size_x = (p_Vid->yuv_format == YUV420 || p_Vid->yuv_format == YUV422) ? 8 : 16;
    p_Vid->mb_size[1][1] = p_Vid->mb_size[2][1] = p_Vid->mb_cr_size_y = (p_Vid->yuv_format == YUV444 || p_Vid->yuv_format == YUV422) ? 16 : 8;

    p_Vid->bitdepth_chroma_qp_scale = 6*(p_Vid->bitdepth_chroma - 8);

    p_Vid->chroma_qp_offset[0] = p_Vid->active_pps->cb_qp_index_offset;
    p_Vid->chroma_qp_offset[1] = p_Vid->active_pps->cr_qp_index_offset;
  }
  else
  {
    p_Vid->bitdepth_chroma     = 0;
    p_Vid->bitdepth_scale[1]   = 0;
    p_Vid->max_pel_value_comp[1] = 0;
    p_Vid->max_pel_value_comp[2] = p_Vid->max_pel_value_comp[1];
    p_Vid->max_imgpel_value_comp_sq[1] = p_Vid->max_pel_value_comp[1] * p_Vid->max_pel_value_comp[1];
    p_Vid->max_imgpel_value_comp_sq[2] = p_Vid->max_pel_value_comp[2] * p_Vid->max_pel_value_comp[2];
    p_Vid->num_blk8x8_uv       = 0;
    p_Vid->num_cdc_coeff       = 0;
    p_Vid->mb_size[1][0] = p_Vid->mb_size[2][0] = p_Vid->mb_cr_size_x = 0;
    p_Vid->mb_size[1][1] = p_Vid->mb_size[2][1] = p_Vid->mb_cr_size_y = 0;

    p_Vid->bitdepth_chroma_qp_scale = 0;
    p_Vid->bitdepth_chroma_qp_scale = 0;

    p_Vid->chroma_qp_offset[0] = 0;
    p_Vid->chroma_qp_offset[1] = 0;
  }  

  p_Vid->max_bitCount =  128 + 256 * p_Vid->bitdepth_luma + 2 * p_Vid->mb_cr_size_y * p_Vid->mb_cr_size_x * p_Vid->bitdepth_chroma;
  //p_Vid->max_bitCount =  (128 + 256 * p_Vid->bitdepth_luma + 2 *p_Vid->mb_cr_size_y * p_Vid->mb_cr_size_x * p_Vid->bitdepth_chroma)*2;

  p_Vid->max_qp_delta = (25 + (p_Vid->bitdepth_luma_qp_scale>>1));
  p_Vid->min_qp_delta = p_Vid->max_qp_delta + 1;

  p_Vid->num_ref_frames = p_Vid->active_sps->num_ref_frames;
  p_Vid->max_num_references = p_Vid->active_sps->frame_mbs_only_flag ? p_Vid->active_sps->num_ref_frames : 2 * p_Vid->active_sps->num_ref_frames;

  p_Vid->base_dist = p_Inp->jumpd + 1;  

  // Intra/IDR related parameters
  p_Vid->lastIntraNumber = 0;
  p_Vid->lastINTRA       = 0;
  p_Vid->last_ref_idc    = 0;
  p_Vid->idr_refresh     = 0;

  p_Vid->DeblockCall     = 0;
  p_Vid->framerate       = (float) p_Inp->output.frame_rate;   // The basic frame rate (of the original sequence)

  if (p_Inp->AdaptiveRounding)
  {
    if (p_Vid->yuv_format != 0)
    {
      get_mem4Dint(&(p_Vid->ARCofAdj4x4), 3, MAXMODE, MB_BLOCK_SIZE, MB_BLOCK_SIZE); //all modes
      get_mem4Dint(&(p_Vid->ARCofAdj8x8), p_Vid->P444_joined ? 3 : 1, MAXMODE, MB_BLOCK_SIZE, MB_BLOCK_SIZE); //modes 0, 1, 2, 3, P8x8
    }     
    else
    {
      get_mem4Dint(&(p_Vid->ARCofAdj4x4), 1, MAXMODE, MB_BLOCK_SIZE, MB_BLOCK_SIZE); //all modes
      get_mem4Dint(&(p_Vid->ARCofAdj8x8), 1, MAXMODE, MB_BLOCK_SIZE, MB_BLOCK_SIZE); //modes 0, 1, 2, 3, P8x8
    }
  }

  imgpel_abs_range = (imax(p_Vid->max_pel_value_comp[0], p_Vid->max_pel_value_comp[1]) + 1) * 2;

  p_Vid->width         = (p_Inp->output.width  + p_Vid->auto_crop_right);
  p_Vid->height        = (p_Inp->output.height + p_Vid->auto_crop_bottom);
  p_Vid->width_blk     = p_Vid->width  / BLOCK_SIZE;
  p_Vid->height_blk    = p_Vid->height / BLOCK_SIZE;
  p_Vid->width_padded  = p_Vid->width  + 2 * IMG_PAD_SIZE;
  p_Vid->height_padded = p_Vid->height + 2 * IMG_PAD_SIZE;

  if (p_Vid->yuv_format != YUV400)
  {
    p_Vid->width_cr = p_Vid->width  * mb_width_cr [p_Vid->yuv_format] / 16;
    p_Vid->height_cr= p_Vid->height * mb_height_cr[p_Vid->yuv_format] / 16;
  }
  else
  {
    p_Vid->width_cr = 0;
    p_Vid->height_cr= 0;
  }

  p_Vid->height_cr_frame = p_Vid->height_cr;

  p_Vid->size = p_Vid->width * p_Vid->height;
  p_Vid->size_cr = p_Vid->width_cr * p_Vid->height_cr;

  p_Vid->PicWidthInMbs    = p_Vid->width  / MB_BLOCK_SIZE;
  p_Vid->FrameHeightInMbs = p_Vid->height / MB_BLOCK_SIZE;
  p_Vid->FrameSizeInMbs   = p_Vid->PicWidthInMbs * p_Vid->FrameHeightInMbs;

  p_Vid->PicHeightInMapUnits = ( p_Vid->active_sps->frame_mbs_only_flag ? p_Vid->FrameHeightInMbs : p_Vid->FrameHeightInMbs >> 1 );

  if ((p_Vid->b8x8info = (Block8x8Info *) calloc(1, sizeof(Block8x8Info))) == NULL)
     no_mem_exit("init_img: p_Vid->block8x8info");

  if( IS_INDEPENDENT(p_Inp) )
  {
    for( i = 0; i < MAX_PLANE; i++ ){
      if ((p_Vid->mb_data_JV[i] = (Macroblock *) calloc(p_Vid->FrameSizeInMbs,sizeof(Macroblock))) == NULL)
        no_mem_exit("init_img: p_Vid->mb_data_JV");
    }
    p_Vid->mb_data = NULL;
  }
  else
  {
    if ((p_Vid->mb_data = (Macroblock *) calloc(p_Vid->FrameSizeInMbs, sizeof(Macroblock))) == NULL)
      no_mem_exit("init_img: p_Vid->mb_data");
  }

  if (p_Inp->UseConstrainedIntraPred)
  {
    if ((p_Vid->intra_block = (int*)calloc(p_Vid->FrameSizeInMbs, sizeof(int))) == NULL)
      no_mem_exit("init_img: p_Vid->intra_block");
  }

  if (p_Inp->CtxAdptLagrangeMult == 1)
  {
    if ((p_Vid->mb16x16_cost_frame = (double*)calloc(p_Vid->FrameSizeInMbs, sizeof(double))) == NULL)
    {
      no_mem_exit("init p_Vid->mb16x16_cost_frame");
    }
  }
  get_mem2D((byte***)&(p_Vid->ipredmode), p_Vid->height_blk, p_Vid->width_blk);        //need two extra rows at right and bottom
  get_mem2D((byte***)&(p_Vid->ipredmode8x8), p_Vid->height_blk, p_Vid->width_blk);     // help storage for ipredmode 8x8, inserted by YV
  memset(&(p_Vid->ipredmode[0][0])   , -1, p_Vid->height_blk * p_Vid->width_blk *sizeof(char));
  memset(&(p_Vid->ipredmode8x8[0][0]), -1, p_Vid->height_blk * p_Vid->width_blk *sizeof(char));


  // CAVLC mem
  get_mem3Dint(&(p_Vid->nz_coeff), p_Vid->FrameSizeInMbs, 4, 4+p_Vid->num_blk8x8_uv);
  
  get_mem2Dolm     (&(p_Vid->lambda)   , 10, 52 + p_Vid->bitdepth_luma_qp_scale, p_Vid->bitdepth_luma_qp_scale);
  get_mem2Dodouble (&(p_Vid->lambda_md), 10, 52 + p_Vid->bitdepth_luma_qp_scale, p_Vid->bitdepth_luma_qp_scale);
  get_mem3Dodouble (&(p_Vid->lambda_me), 10, 52 + p_Vid->bitdepth_luma_qp_scale, 3, p_Vid->bitdepth_luma_qp_scale);
  get_mem3Doint    (&(p_Vid->lambda_mf), 10, 52 + p_Vid->bitdepth_luma_qp_scale, 3, p_Vid->bitdepth_luma_qp_scale);

  if (p_Inp->CtxAdptLagrangeMult == 1)
  {
    get_mem2Dodouble(&(p_Vid->lambda_mf_factor), 10, 52 + p_Vid->bitdepth_luma_qp_scale, p_Vid->bitdepth_luma_qp_scale);
  }

  p_Vid->mb_y_upd  = 0;

  RandomIntraInit (p_Vid, p_Vid->PicWidthInMbs, p_Vid->FrameHeightInMbs, p_Inp->RandomIntraMBRefresh);

  InitSEIMessages(p_Vid, p_Inp); 

  initInput(p_Vid, &p_Inp->source, &p_Inp->output);

  // Allocate I/O Frame memory
  AllocateFrameMemory(p_Vid, p_Inp, &p_Inp->source);

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

  p_Vid->ChromaArrayType = p_Inp->separate_colour_plane_flag ? 0 : p_Inp->output.yuv_format;
  p_Vid->colour_plane_id = 0;

  if (p_Inp->RDPictureDecision)
    p_Vid->frm_iter = 3;
  else
    p_Vid->frm_iter = 1;

  p_Vid->max_frame_num = 1 << (p_Vid->log2_max_frame_num_minus4 + 4);
  p_Vid->max_pic_order_cnt_lsb = 1 << (p_Vid->log2_max_pic_order_cnt_lsb_minus4 + 4);

  p_Vid->prev_frame_no = 0; // POC200301
  p_Vid->consecutive_non_reference_pictures = 0; // POC200301

  p_Vid->p_Inp = p_Inp;

  create_context_memory (p_Vid, p_Inp);
}


/*!
 ***********************************************************************
 * \brief
 *    Free the Image structures
 * \par Input:
 *    Image Parameters VideoParameters *p_Vid
 ***********************************************************************
 */
static void free_img (VideoParameters *p_Vid, InputParameters *p_Inp)
{
  // Delete Frame memory 
  DeleteFrameMemory(p_Vid);

  CloseSEIMessages(p_Vid, p_Inp); 

  free_context_memory (p_Vid);

  if (p_Inp->AdaptiveRounding)
  {
    free_mem4Dint(p_Vid->ARCofAdj4x4);
    free_mem4Dint(p_Vid->ARCofAdj8x8);
  }


  free (p_Vid->p_SEI);
  free (p_Vid->p_QScale);
  free (p_Vid->p_Quant);
  free (p_Vid->p_Dpb);
  free (p_Vid->p_Stats);
  free (p_Vid->p_Dist);
  free (p_Vid);
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
int init_orig_buffers(VideoParameters *p_Vid, ImageData *imgData)
{
  InputParameters *p_Inp = p_Vid->p_Inp;
  int memory_size = 0;
  int nplane;

  // allocate memory for reference frame buffers: imgData->frm_data
  imgData->format           = p_Inp->output;
  imgData->format.width     = p_Vid->width;    
  imgData->format.height    = p_Vid->height;
  imgData->format.width_cr  = p_Vid->width_cr;
  imgData->format.height_cr = p_Vid->height_cr;
  imgData->format.yuv_format = p_Vid->yuv_format;
  imgData->format.auto_crop_bottom = p_Vid->auto_crop_bottom;
  imgData->format.auto_crop_right  = p_Vid->auto_crop_right;
  imgData->format.auto_crop_bottom_cr = (p_Vid->auto_crop_bottom * mb_height_cr [p_Vid->yuv_format]) / MB_BLOCK_SIZE;
  imgData->format.auto_crop_right_cr  = (p_Vid->auto_crop_right * mb_width_cr [p_Vid->yuv_format]) / MB_BLOCK_SIZE;

  if( IS_INDEPENDENT(p_Inp) )
  {

    for( nplane=0; nplane<MAX_PLANE; nplane++ )
    {
      memory_size += get_mem2Dpel(&(imgData->frm_data[nplane]), p_Vid->height, p_Vid->width);
    }
  }
  else
  {
    //imgData->format = p_Inp->input_file1.format;    

    memory_size += get_mem2Dpel(&(imgData->frm_data[0]), p_Vid->height, p_Vid->width);

    if (p_Vid->yuv_format != YUV400)
    {
      int i, j, k;
      memory_size += get_mem2Dpel(&(imgData->frm_data[1]), p_Vid->height_cr, p_Vid->width_cr);
      memory_size += get_mem2Dpel(&(imgData->frm_data[2]), p_Vid->height_cr, p_Vid->width_cr);

      if (sizeof(imgpel) == sizeof(unsigned char))
      {
        for (k = 1; k < 3; k++)
          memset(&(imgData->frm_data[k][0][0]), 128, p_Vid->height_cr * p_Vid->width_cr * sizeof(imgpel));
      }
      else
      {
        for (k = 1; k < 3; k++)
          for (j = 0; j < p_Vid->height_cr; j++)
            for (i = 0; i < p_Vid->width_cr; i++)
              imgData->frm_data[k][j][i] = 128;
      }
    }
  }

  if (!p_Vid->active_sps->frame_mbs_only_flag)
  {
    // allocate memory for field reference frame buffers
    memory_size += init_top_bot_planes(imgData->frm_data[0], p_Vid->height, &(imgData->top_data[0]), &(imgData->bot_data[0]));

    if (p_Vid->yuv_format != YUV400)
    {

      memory_size += 4*(sizeof(imgpel**));

      memory_size += init_top_bot_planes(imgData->frm_data[1], p_Vid->height_cr, &(imgData->top_data[1]), &(imgData->bot_data[1]));
      memory_size += init_top_bot_planes(imgData->frm_data[2], p_Vid->height_cr, &(imgData->top_data[2]), &(imgData->bot_data[2]));
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
 *    Image Parameters VideoParameters *p_Vid
 * \return Number of allocated bytes
 ************************************************************************
 */
static int init_global_buffers(VideoParameters *p_Vid, InputParameters *p_Inp)
{
  int j, memory_size=0;

  if ((p_Vid->enc_frame_picture = (StorablePicture**)malloc(6 * sizeof(StorablePicture*))) == NULL)
    no_mem_exit("init_global_buffers: *p_Vid->enc_frame_picture");

  for (j = 0; j < 6; j++)
    p_Vid->enc_frame_picture[j] = NULL;

  if ((p_Vid->enc_field_picture = (StorablePicture**)malloc(2 * sizeof(StorablePicture*))) == NULL)
    no_mem_exit("init_global_buffers: *p_Vid->enc_field_picture");

  for (j = 0; j < 2; j++)
    p_Vid->enc_field_picture[j] = NULL;

  if ((p_Vid->frame_pic = (Picture**)malloc(p_Vid->frm_iter * sizeof(Picture*))) == NULL)
    no_mem_exit("init_global_buffers: *p_Vid->frame_pic");

  for (j = 0; j < p_Vid->frm_iter; j++)
    p_Vid->frame_pic[j] = malloc_picture();

  if (p_Inp->si_frame_indicator || p_Inp->sp_periodicity)
  {
    p_Vid->number_sp2_frames=0;
    p_Vid->frame_pic_si = malloc_picture();//picture buffer for the encoded SI picture
    //allocation of lrec and p_Vid->lrec_uv for SI picture
    get_mem2Dint (&p_Vid->lrec, p_Vid->height, p_Vid->width);
    get_mem3Dint (&p_Vid->lrec_uv, 2, p_Vid->height, p_Vid->width);
  }

  // Allocate memory for field picture coding
  if (p_Inp->PicInterlace != FRAME_CODING)
  { 
    if ((p_Vid->field_pic = (Picture**)malloc(2 * sizeof(Picture*))) == NULL)
      no_mem_exit("init_global_buffers: *p_Vid->field_pic");

    for (j = 0; j < 2; j++)
      p_Vid->field_pic[j] = malloc_picture();
  }

  // Init memory data for input & encoded images
  memory_size += init_orig_buffers(p_Vid, &p_Vid->imgData);
  memory_size += init_orig_buffers(p_Vid, &p_Vid->imgData0);
  
  memory_size += get_mem2Dshort(&PicPos, p_Vid->FrameSizeInMbs + 1, 2);

  for (j = 0; j < (int) p_Vid->FrameSizeInMbs + 1; j++)
  {
    PicPos[j][0] = (short) (j % p_Vid->PicWidthInMbs);
    PicPos[j][1] = (short) (j / p_Vid->PicWidthInMbs);
  }


  if (p_Inp->rdopt == 3)
  {
    memory_size += allocate_errdo_mem(p_Vid, p_Inp);
  }

  if (p_Inp->RestrictRef)
  {
    memory_size += get_mem2D(&p_Vid->pixel_map,   p_Vid->height,   p_Vid->width);
    memory_size += get_mem2D(&p_Vid->refresh_map, p_Vid->height >> 3, p_Vid->width >> 3);
  }

  if (!p_Vid->active_sps->frame_mbs_only_flag)
  {
    memory_size += get_mem2Dpel(&p_Vid->imgY_com, p_Vid->height, p_Vid->width);

    if (p_Vid->yuv_format != YUV400)
    {
      memory_size += get_mem3Dpel(&p_Vid->imgUV_com, 2, p_Vid->height_cr, p_Vid->width_cr);
    }
  }

  // allocate and set memory relating to motion estimation
  if (!p_Inp->IntraProfile)
  {  
    if (p_Inp->SearchMode == UM_HEX)
    {
      if ((p_Vid->p_UMHex = (UMHexStruct*)calloc(1, sizeof(UMHexStruct))) == NULL)
        no_mem_exit("init_mv_block: p_Vid->p_UMHex");
      memory_size += UMHEX_get_mem(p_Vid, p_Inp);
    }
    else if (p_Inp->SearchMode == UM_HEX_SIMPLE)
    {
      if ((p_Vid->p_UMHexSMP = (UMHexSMPStruct*)calloc(1, sizeof(UMHexSMPStruct))) == NULL)
        no_mem_exit("init_mv_block: p_Vid->p_UMHexSMP");

      smpUMHEX_init(p_Vid);
      memory_size += smpUMHEX_get_mem(p_Vid);
    }
    else if (p_Inp->SearchMode == EPZS)
    {
      memory_size += EPZSInit(p_Vid);
    }

  }

  if (p_Inp->RCEnable)
    rc_allocate_memory(p_Vid, p_Inp);

  if (p_Inp->redundant_pic_flag)
  {
    memory_size += get_mem2Dpel(&p_Vid->imgY_tmp, p_Vid->height, p_Vid->width);
    memory_size += get_mem2Dpel(&p_Vid->imgUV_tmp[0], p_Vid->height_cr, p_Vid->width_cr);
    memory_size += get_mem2Dpel(&p_Vid->imgUV_tmp[1], p_Vid->height_cr, p_Vid->width_cr);
  }

  memory_size += get_mem2Dint (&p_Vid->imgY_sub_tmp, p_Vid->height_padded, p_Vid->width_padded);

  if ( p_Inp->ChromaMCBuffer )
    chroma_mc_setup(p_Vid);

  p_Vid->padded_size_x       = (p_Vid->width + 2 * IMG_PAD_SIZE);
  p_Vid->padded_size_x_m8x8  = (p_Vid->padded_size_x - BLOCK_SIZE_8x8);
  p_Vid->padded_size_x_m4x4  = (p_Vid->padded_size_x - BLOCK_SIZE);
  p_Vid->cr_padded_size_x    = (p_Vid->width_cr + 2 * p_Vid->pad_size_uv_x);
  p_Vid->cr_padded_size_x2   = (p_Vid->cr_padded_size_x << 1);
  p_Vid->cr_padded_size_x4   = (p_Vid->cr_padded_size_x << 2);
  p_Vid->cr_padded_size_x_m8 = (p_Vid->cr_padded_size_x - 8);

  // RGB images for distortion calculation
  // Recommended to do this allocation (and de-allocation) in 
  // the appropriate file instead of here.
  if(p_Inp->DistortionYUVtoRGB)
  {
    memory_size += create_RGB_memory(p_Vid);
  }

  p_Vid->pWPX = NULL;
  if ( p_Inp->WPMCPrecision )
  {
    wpxInitWPXObject(p_Vid);
  }

  memory_size += InitProcessImage( p_Vid, p_Inp );

  p_Vid->p_pred = init_seq_structure( p_Vid, p_Inp, &memory_size );

  return memory_size;
}


/*!
 ************************************************************************
 * \brief
 *    Free allocated memory of original picture buffers
 ************************************************************************
 */
void free_orig_planes(VideoParameters *p_Vid, ImageData *imgData)
{
  if( IS_INDEPENDENT(p_Vid->p_Inp) )
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

  if (!p_Vid->active_sps->frame_mbs_only_flag)
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
 *    Image Parameters VideoParameters *p_Vid
 * \par Output:
 *    none
 ************************************************************************
 */
static void free_global_buffers(VideoParameters *p_Vid, InputParameters *p_Inp)
{
  int  i,j;

  if (p_Vid->enc_frame_picture)
    free (p_Vid->enc_frame_picture);
  if (p_Vid->frame_pic)
  {
    for (j = 0; j < p_Vid->frm_iter; j++)
    {
      if (p_Vid->frame_pic[j])
        free_picture (p_Vid->frame_pic[j]);
    }
    free (p_Vid->frame_pic);
  }

  if (p_Vid->enc_field_picture)
    free (p_Vid->enc_field_picture);
  if (p_Vid->field_pic)
  {
    for (j = 0; j < 2; j++)
    {
      if (p_Vid->field_pic[j])
        free_picture (p_Vid->field_pic[j]);
    }
    free (p_Vid->field_pic);
  }

  // Deallocation of SI picture related memory
  if (p_Inp->si_frame_indicator || p_Inp->sp_periodicity)
  {
    free_picture (p_Vid->frame_pic_si);
    //deallocation of lrec and p_Vid->lrec_uv for SI frames
    free_mem2Dint (p_Vid->lrec);
    free_mem3Dint (p_Vid->lrec_uv);
  }

  free_orig_planes(p_Vid, &p_Vid->imgData);
  free_orig_planes(p_Vid, &p_Vid->imgData0);

  // free lookup memory which helps avoid divides with PicWidthInMbs
  free_mem2Dshort(PicPos);
  // Free Qmatrices and offsets
  free_QMatrix(p_Vid->p_Quant);
  free_QOffsets(p_Vid->p_Quant, p_Inp);


  if ( p_Inp->WPMCPrecision )
  {
    wpxFreeWPXObject(p_Vid);
  }

  if (p_Vid->imgY_sub_tmp) // free temp quarter pel frame buffers
  {
    free_mem2Dint (p_Vid->imgY_sub_tmp);
    p_Vid->imgY_sub_tmp = NULL;
  }

  // free mem, allocated in init_img()
  // free intra pred mode buffer for blocks
  free_mem2D((byte**)p_Vid->ipredmode);
  free_mem2D((byte**)p_Vid->ipredmode8x8);

  if(p_Vid->ipredmode4x4_line)
  {
    free_mem2D((byte**)p_Vid->ipredmode4x4_line);
    p_Vid->ipredmode4x4_line=NULL;
  }
  if(p_Vid->ipredmode8x8_line)
  {
    free_mem2D((byte **)p_Vid->ipredmode8x8_line);
    p_Vid->ipredmode8x8_line = NULL;
  }

  if( IS_INDEPENDENT(p_Inp) )
  {
    for( i=0; i<MAX_PLANE; i++ ){
      free(p_Vid->mb_data_JV[i]);
    }
  }
  else
  {
    free(p_Vid->mb_data);
  }

  if(p_Inp->UseConstrainedIntraPred)
  {
    free (p_Vid->intra_block);
  }

  if (p_Inp->CtxAdptLagrangeMult == 1)
  {
    free(p_Vid->mb16x16_cost_frame);
  }

  if (p_Inp->rdopt == 3)
  {
    free_errdo_mem(p_Vid);
  }

  if (p_Inp->RestrictRef)
  {
    free(p_Vid->pixel_map[0]);
    free(p_Vid->pixel_map);
    free(p_Vid->refresh_map[0]);
    free(p_Vid->refresh_map);
  }

  if (!p_Vid->active_sps->frame_mbs_only_flag)
  {
    free_mem2Dpel(p_Vid->imgY_com);
    
    if (p_Vid->yuv_format != YUV400)
    {
      free_mem3Dpel(p_Vid->imgUV_com);
    }
  }

  free_mem3Dint(p_Vid->nz_coeff);

  free_mem2Dolm     (p_Vid->lambda, p_Vid->bitdepth_luma_qp_scale);
  free_mem2Dodouble (p_Vid->lambda_md, p_Vid->bitdepth_luma_qp_scale);
  free_mem3Dodouble (p_Vid->lambda_me, 10, 52 + p_Vid->bitdepth_luma_qp_scale, p_Vid->bitdepth_luma_qp_scale);
  free_mem3Doint    (p_Vid->lambda_mf, 10, 52 + p_Vid->bitdepth_luma_qp_scale, p_Vid->bitdepth_luma_qp_scale);

  if (p_Inp->CtxAdptLagrangeMult == 1)
  {
    free_mem2Dodouble(p_Vid->lambda_mf_factor, p_Vid->bitdepth_luma_qp_scale);
  }

  if (!p_Inp->IntraProfile)
  {
    if (p_Inp->SearchMode == UM_HEX)
    {
      UMHEX_free_mem(p_Vid, p_Inp);
    }
    else if (p_Inp->SearchMode == UM_HEX_SIMPLE)
    {
      smpUMHEX_free_mem(p_Vid);
    }
    else if (p_Inp->SearchMode == EPZS)
    {
      EPZSDelete(p_Vid);
    }
  }

  if (p_Inp->RCEnable)
    rc_free_memory(p_Vid, p_Inp);

  if (p_Inp->redundant_pic_flag)
  {
    free_mem2Dpel(p_Vid->imgY_tmp);
    free_mem2Dpel(p_Vid->imgUV_tmp[0]);
    free_mem2Dpel(p_Vid->imgUV_tmp[1]);
  }

  // Again process should be moved into cconv_yuv2rgb.c file for cleanliness
  // These should not be globals but instead only be visible through that code.
  if(p_Inp->DistortionYUVtoRGB)
  {
    delete_RGB_memory(p_Vid);
  }

  ClearProcessImage( p_Vid, p_Inp );

  free_seq_structure( p_Vid->p_pred );
}


/*!
 ************************************************************************
 * \brief
 *    Allocate memory for AC coefficients
 ************************************************************************
 */
int get_mem_ACcoeff (VideoParameters *p_Vid, int***** cofAC)
{
  int num_blk8x8 = BLOCK_SIZE + p_Vid->num_blk8x8_uv;
  
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
static void SetLevelIndices(VideoParameters *p_Vid)
{
  switch(p_Vid->active_sps->level_idc)
  {
  case 9:
    p_Vid->LevelIndex=1;
    break;
  case 10:
    p_Vid->LevelIndex=0;
    break;
  case 11:
    if (!IS_FREXT_PROFILE(p_Vid->active_sps->profile_idc) && (p_Vid->active_sps->constrained_set3_flag == 0))
      p_Vid->LevelIndex=2;
    else
      p_Vid->LevelIndex=1;
    break;
  case 12:
    p_Vid->LevelIndex=3;
    break;
  case 13:
    p_Vid->LevelIndex=4;
    break;
  case 20:
    p_Vid->LevelIndex=5;
    break;
  case 21:
    p_Vid->LevelIndex=6;
    break;
  case 22:
    p_Vid->LevelIndex=7;
    break;
  case 30:
    p_Vid->LevelIndex=8;
    break;
  case 31:
    p_Vid->LevelIndex=9;
    break;
  case 32:
    p_Vid->LevelIndex=10;
    break;
  case 40:
    p_Vid->LevelIndex=11;
    break;
  case 41:
    p_Vid->LevelIndex=12;
    break;
  case 42:
    if (!IS_FREXT_PROFILE(p_Vid->active_sps->profile_idc))
      p_Vid->LevelIndex=13;
    else
      p_Vid->LevelIndex=14;
    break;
  case 50:
    p_Vid->LevelIndex=15;
    break;
  case 51:
    p_Vid->LevelIndex=16;
    break;
  default:
    fprintf ( stderr, "Warning: unknown LevelIDC, using maximum level 5.1 \n" );
    p_Vid->LevelIndex=16;
    break;
  }
}

/*!
 ************************************************************************
 * \brief
 *    initialize key frames and corresponding redundant frames.
 ************************************************************************
 */
void Init_redundant_frame(VideoParameters *p_Vid, InputParameters *p_Inp)
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

  p_Vid->key_frame = 0;
  p_Vid->redundant_coding = 0;
  p_Vid->redundant_pic_cnt = 0;
  p_Vid->frameNuminGOP = p_Vid->curr_frm_idx % p_Inp->PrimaryGOPLength;
  if (p_Vid->curr_frm_idx == 0)
  {
    p_Vid->frameNuminGOP = -1;
  }
}

/*!
 ************************************************************************
 * \brief
 *    allocate redundant frames in a primary GOP.
 ************************************************************************
 */
void Set_redundant_frame(VideoParameters *p_Vid, InputParameters *p_Inp)
{
  int GOPlength = p_Inp->PrimaryGOPLength;

  //start frame of GOP
  if (p_Vid->frameNuminGOP == 0)
  {
    p_Vid->redundant_coding = 0;
    p_Vid->key_frame = 1;
    p_Vid->redundant_ref_idx = GOPlength;
  }

  //1/2 position
  if (p_Inp->NumRedundantHierarchy > 0)
  {
    if (p_Vid->frameNuminGOP == GOPlength >> 1)
    {
      p_Vid->redundant_coding = 0;
      p_Vid->key_frame = 1;
      p_Vid->redundant_ref_idx = GOPlength >> 1;
    }
  }

  //1/4, 3/4 position
  if (p_Inp->NumRedundantHierarchy > 1)
  {
    if (p_Vid->frameNuminGOP == (GOPlength >> 2) || p_Vid->frameNuminGOP == ((GOPlength*3) >> 2))
    {
      p_Vid->redundant_coding = 0;
      p_Vid->key_frame = 1;
      p_Vid->redundant_ref_idx = GOPlength >> 2;
    }
  }

  //1/8, 3/8, 5/8, 7/8 position
  if (p_Inp->NumRedundantHierarchy > 2)
  {
    if (p_Vid->frameNuminGOP == GOPlength >> 3 || p_Vid->frameNuminGOP == ((GOPlength*3) >> 3)
      || p_Vid->frameNuminGOP == ((GOPlength*5) >> 3) || p_Vid->frameNuminGOP == ((GOPlength*7) & 0x03))
    {
      p_Vid->redundant_coding = 0;
      p_Vid->key_frame = 1;
      p_Vid->redundant_ref_idx = GOPlength >> 3;
    }
  }

  //1/16, 3/16, 5/16, 7/16, 9/16, 11/16, 13/16 position
  if (p_Inp->NumRedundantHierarchy > 3)
  {
    if (p_Vid->frameNuminGOP == (GOPlength >> 4) || p_Vid->frameNuminGOP == ((GOPlength*3) >> 4)
      || p_Vid->frameNuminGOP == ((GOPlength*5) >> 4) || p_Vid->frameNuminGOP == ((GOPlength*7) >> 4)
      || p_Vid->frameNuminGOP == ((GOPlength*9) >> 4) || p_Vid->frameNuminGOP == ((GOPlength*11) >> 4)
      || p_Vid->frameNuminGOP == ((GOPlength*13) >> 4))
    {
      p_Vid->redundant_coding = 0;
      p_Vid->key_frame = 1;
      p_Vid->redundant_ref_idx = GOPlength >> 4;
    }
  }
}

/*!
 ************************************************************************
 * \brief
 *    encode one redundant frame.
 ************************************************************************
 */
void encode_one_redundant_frame(VideoParameters *p_Vid, InputParameters *p_Inp)
{
  p_Vid->key_frame = 0;
  p_Vid->redundant_coding = 1;
  p_Vid->redundant_pic_cnt = 1;

  if (!p_Vid->currentPicture->idr_flag)
  {
    if (p_Vid->type == I_SLICE)
    {
      set_slice_type( p_Vid, p_Inp, P_SLICE );
    }
  }

  encode_one_frame(p_Vid, p_Inp);
}

/*!
 ************************************************************************
 * \brief
 *    Setup Chroma MC Variables
 ************************************************************************
 */
static void chroma_mc_setup(VideoParameters *p_Vid)
{
  // initialize global variables used for chroma interpolation and buffering
  if ( p_Vid->yuv_format == YUV420 )
  {
    p_Vid->pad_size_uv_x = IMG_PAD_SIZE >> 1;
    p_Vid->pad_size_uv_y = IMG_PAD_SIZE >> 1;
    p_Vid->chroma_mask_mv_y = 7;
    p_Vid->chroma_mask_mv_x = 7;
    p_Vid->chroma_shift_x = 3;
    p_Vid->chroma_shift_y = 3;
  }
  else if ( p_Vid->yuv_format == YUV422 )
  {
    p_Vid->pad_size_uv_x = IMG_PAD_SIZE >> 1;
    p_Vid->pad_size_uv_y = IMG_PAD_SIZE;
    p_Vid->chroma_mask_mv_y = 3;
    p_Vid->chroma_mask_mv_x = 7;
    p_Vid->chroma_shift_y = 2;
    p_Vid->chroma_shift_x = 3;
  }
  else
  { // YUV444
    p_Vid->pad_size_uv_x = IMG_PAD_SIZE;
    p_Vid->pad_size_uv_y = IMG_PAD_SIZE;
    p_Vid->chroma_mask_mv_y = 3;
    p_Vid->chroma_mask_mv_x = 3;
    p_Vid->chroma_shift_y = 2;
    p_Vid->chroma_shift_x = 2;
  }
  p_Vid->shift_cr_y  = p_Vid->chroma_shift_y - 2;
  p_Vid->shift_cr_x  = p_Vid->chroma_shift_x - 2;
}

