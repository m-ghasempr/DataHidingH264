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
 *     JM 15.1 (FRExt)
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
#include <sys/timeb.h>
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
#include "memalloc.h"
#include "me_epzs.h"
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
#include "mv-search.h"

InputParameters  *params = NULL;
ImageParameters  *img = NULL;


int initial_Bframes;

int    cabac_encoding;
int    frame_statistic_start;
extern ColocatedParams *Co_located;
extern ColocatedParams *Co_located_JV[MAX_PLANE];  //!< Co_located to be used during 4:4:4 independent mode encoding
extern double *mb16x16_cost_frame;

extern int  allocate_block_mem(void);
extern void free_block_mem(void);
extern int allocate_mb_mem(void);
extern void free_mb_mem(void);

static void SetLevelIndices(void);
static void chroma_mc_setup(void);
static int  get_mem_mv  (short*******);
static int  get_mem_bipred_mv (short********);
static void free_mem_bipred_mv (short*******);
static void free_mem_mv (short******);
static void init_img( ImageParameters *img, InputParameters *params);
static void init_poc(void);

static void init_encoder(ImageParameters *img);
static void encode_sequence(ImageParameters *img);


void init_stats (StatParameters *stats)
{
  memset(stats, 0, sizeof(StatParameters));
  stats->NumberBFrames = params->NumberBFrames;
}

void init_dstats (DistortionParams *dist)
{
  dist->frame_ctr = 0;
  memset(dist->metric, 0, TOTAL_DIST_TYPES * sizeof(DistMetric));
}

/*!
 ***********************************************************************
 * \brief
 *    Initialize encoding parameters.
 ***********************************************************************
 */
static void init_frame_params(ImageParameters *p_img)
{
  int base_mul = 0;

  if (params->idr_period)
  {
    if (!params->adaptive_idr_period && ( p_img->frm_number - p_img->lastIDRnumber ) % params->idr_period == 0 )
      p_img->nal_reference_idc = NALU_PRIORITY_HIGHEST;

    if (params->adaptive_idr_period == 1 && ( p_img->frm_number - imax(p_img->lastIntraNumber, p_img->lastIDRnumber) ) % params->idr_period == 0 )
      p_img->nal_reference_idc = NALU_PRIORITY_HIGHEST;
    else
      p_img->nal_reference_idc = (params->DisposableP) ? (p_img->frm_number + 1)% 2 : NALU_PRIORITY_LOW;
  }
  else
    p_img->nal_reference_idc = (p_img->frm_number && params->DisposableP) ? (p_img->frm_number + 1)% 2 : NALU_PRIORITY_LOW;

  //much of this can go in init_frame() or init_field()?
  //poc for this frame or field
  if (params->idr_period)
  {
    if (!params->adaptive_idr_period)
      base_mul = ( p_img->frm_number - p_img->lastIDRnumber ) % params->idr_period;
    else if (params->adaptive_idr_period == 1)
      base_mul = (( p_img->frm_number - imax(p_img->lastIntraNumber, p_img->lastIDRnumber) ) % params->idr_period == 0) ? 0 : ( p_img->frm_number - p_img->lastIDRnumber );
  }
  else 
    base_mul = ( p_img->frm_number - p_img->lastIDRnumber );

  if ((p_img->frm_number - p_img->lastIDRnumber) <= params->intra_delay)
  {    
    base_mul = -base_mul;
  }
  else
  {
    base_mul -= ( base_mul ? params->intra_delay :  0);    
  }

  p_img->toppoc = base_mul * (2 * p_img->base_dist);

  if ((params->PicInterlace==FRAME_CODING) && (params->MbInterlace==FRAME_CODING))
    p_img->bottompoc = p_img->toppoc;     //progressive
  else
    p_img->bottompoc = p_img->toppoc + 1;   //hard coded

  p_img->framepoc = imin (p_img->toppoc, p_img->bottompoc);

  //the following is sent in the slice header
  p_img->delta_pic_order_cnt[0] = 0;

  if ((params->BRefPictures == 1) && (p_img->frm_number))
  {
    p_img->delta_pic_order_cnt[0] = 2 * params->NumberBFrames;
  }  

  if (params->NumberBFrames && params->last_frame && ((p_img->gop_number) + 1) == params->no_frm_base)
  {
    int bi = (int)((float)p_img->base_dist / (initial_Bframes + 1.0) + 0.499999);
    int new_bframes = ((params->last_frame - (p_img->frm_number - 1) * p_img->base_dist) / bi) - 1;

    //about to code the last ref frame, adjust delta poc
    p_img->delta_pic_order_cnt[0]= -2*(initial_Bframes - new_bframes);
    p_img->toppoc    += p_img->delta_pic_order_cnt[0];
    p_img->bottompoc += p_img->delta_pic_order_cnt[0];
    p_img->framepoc   = imin (p_img->toppoc, p_img->bottompoc);
  }

  //frame_num for this frame
  if (params->idr_period && ((!params->adaptive_idr_period && ( p_img->frm_number - p_img->lastIDRnumber ) % params->idr_period == 0)
    || (params->adaptive_idr_period == 1 && ( p_img->frm_number - imax(p_img->lastIntraNumber, p_img->lastIDRnumber) ) % params->idr_period == 0)) )
  {
    p_img->frame_num = 0;
  }
}

/*!
 ***********************************************************************
 * \brief
 *    Allocate the Image structure
 * \par  Output:
 *    Image Parameters ImageParameters *img
 ***********************************************************************
 */
static void alloc_img( ImageParameters **p_img)
{
  if ((*p_img = (ImageParameters *) calloc(1, sizeof(ImageParameters)))==NULL) 
    no_mem_exit("alloc_img: p_img");

  get_mem3Dint(&((*p_img)->mb_rres),   MAX_PLANE, MB_BLOCK_SIZE, MB_BLOCK_SIZE);
  get_mem3Dint(&((*p_img)->mb_ores),   MAX_PLANE, MB_BLOCK_SIZE, MB_BLOCK_SIZE);
  get_mem4Dpel(&((*p_img)->mpr_4x4),   MAX_PLANE, 9, MB_BLOCK_SIZE, MB_BLOCK_SIZE);
  get_mem4Dpel(&((*p_img)->mpr_8x8),   MAX_PLANE, 9, MB_BLOCK_SIZE, MB_BLOCK_SIZE);
  get_mem4Dpel(&((*p_img)->mpr_16x16), MAX_PLANE, 5, MB_BLOCK_SIZE, MB_BLOCK_SIZE);
  get_mem3Dpel(&((*p_img)->mb_pred),   MAX_PLANE, MB_BLOCK_SIZE, MB_BLOCK_SIZE);
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
  p_dec = -1;

  alloc_img( &img);
  allocate_block_mem();
  allocate_mb_mem();  

  params = (InputParameters *) malloc(sizeof(InputParameters));

  Configure (argc, argv);

  // init encoder
  init_encoder(img);

  // encode sequence
  encode_sequence(img);

  // terminate sequence
  free_mb_mem();
  free_block_mem();
  free_encoder_memory(img);

  if ( params != NULL )
    free( params );

  return 0;
}


/*!
 ***********************************************************************
 * \brief
 *    Initialize encoder
 ***********************************************************************
 */

static void init_encoder(ImageParameters *img)
{
  giRDOpt_B8OnlyFlag = 0;

  p_log = NULL;

  cabac_encoding = 0;

  frame_statistic_start = 1;

  // Open Files
  OpenFiles(&params->input_file1);

  Init_QMatrix();
  Init_QOffsetMatrix();

  init_poc();
  GenerateParameterSets();
  SetLevelIndices();

  init_img(img, params);
  init_rdopt(params);


  if (params->rdopt == 3)
  {
    init_error_conceal(params->ErrorConcealment); 
  }

#ifdef _LEAKYBUCKET_
  initial_Bframes = 0;
  Bit_Buffer = (long *)malloc((params->no_frames + 1) * sizeof(long));
  total_frame_buffer = 0;
#endif

  // Prepare hierarchical coding structures. 
  // Code could be extended in the future to allow structure adaptation.
  if (params->HierarchicalCoding)
  {
    init_gop_structure(params);
    if (params->NumberBFrames && params->HierarchicalCoding == 3)
      interpret_gop_structure();
    else
      create_hierarchy(params);
  }

  dpb.init_done = 0;
  init_dpb();
  init_out_buffer();
  init_stats (stats);
  init_dstats(dist);

  enc_picture = NULL;

  init_global_buffers();

  if ( params->WPMCPrecision )
  {
    wpxInitWPXPasses(params);
  }

  create_context_memory ();
  Init_Motion_Search_Module ();
  information_init(img, params, stats);

  if(params->DistortionYUVtoRGB)
    init_YUVtoRGB();

  //Rate control
  if (params->RCEnable)
    rc_init_sequence();

  img->last_valid_reference = 0;
  tot_time = 0;                 // time for total encoding session

  initial_Bframes = params->NumberBFrames;

  PatchInputNoFrames();

  img->type = I_SLICE;
  // Write sequence header (with parameter sets)
  stats->bit_ctr_parametersets = 0;
  stats->bit_slice = start_sequence();

  if (params->UseRDOQuant)
    precalculate_unary_exp_golomb_level();

   if (params->ExplicitSeqCoding)
     OpenExplicitSeqFile(params);
}

/*!
 ***********************************************************************
 * \brief
 *    Determine coding level a frame belongs to
 ***********************************************************************
 */
static int determine_coding_level(ImageParameters *p_img, InputParameters *pparams, int curr_frame)
{
  int coding_level = 0;  

  if (curr_frame - p_img->last_idr_number == 0)
    return coding_level;
  else
    coding_level  = (curr_frame - p_img->last_idr_number - 1) % (1 + pparams->NumberBFrames);
  
  return coding_level;
}

/*!
 ************************************************************************
 * \brief
 *    Set the image type for I,P and SP pictures (not B!)
 ************************************************************************
 */
static void SetImgType(ImageParameters *img, int gop_frame_num)
{
  if (gop_frame_num == 0)
  {
    int intra_refresh = (params->intra_period == 0) ? (img->gop_number == 0) : (( ( img->frm_number - img->lastIntraNumber) % params->intra_period ) == 0);
    int idr_refresh;

    if ( params->idr_period && !params->adaptive_idr_period )
      idr_refresh = (( ( img->frm_number - img->lastIDRnumber  ) % params->idr_period   ) == 0);
    else if ( params->idr_period && params->adaptive_idr_period == 1 )
      idr_refresh = (( ( img->frm_number - imax(img->lastIntraNumber, img->lastIDRnumber)  ) % params->idr_period   ) == 0);
    else
      idr_refresh = (img->gop_number == 0);

    if (intra_refresh || idr_refresh)
    {
      set_slice_type( img, I_SLICE );        // set image type for first image to I-frame
    }
    else
    {
      set_slice_type(img, (params->sp_periodicity && ((img->gop_number % params->sp_periodicity) == 0))
        ? SP_SLICE  : ((params->BRefPictures == 2) ? B_SLICE : P_SLICE) );
    }
  }
  else
  {
    if (params->HierarchicalCoding)
      set_slice_type( img, gop_structure[gop_frame_num - 1].slice_type);
    else
      set_slice_type( img, ( params->PReplaceBSlice ) ? P_SLICE : B_SLICE);
  }
}

static void set_poc(ImageParameters *img, double frame_to_code)
{
  img->toppoc = (int) (img->frame_interval * frame_to_code);
  if (img->gop_number && (img->gop_number <= params->intra_delay))
  {
    if(params->idr_period && ((!params->adaptive_idr_period && ( img->frm_number - img->lastIDRnumber ) % params->idr_period == 0)
      || (params->adaptive_idr_period == 1 && ( img->frm_number - imax(img->lastIntraNumber, img->lastIDRnumber) ) % params->idr_period == 0)) )
      img->toppoc = img->toppoc;
    else
    {
      img->toppoc += (-img->frm_number + img->lastIDRnumber) * img->base_dist;
      img->toppoc *= 2;
    }
  }
  else
  {
    if(params->idr_period && !params->adaptive_idr_period)
      img->toppoc += ((( img->frm_number - img->lastIDRnumber - params->intra_delay) % params->idr_period ) - 1) * img->base_dist;
    else
      img->toppoc += (( img->frm_number - img->lastIDRnumber - params->intra_delay - 1 ) * img->base_dist);
    img->toppoc *= 2;
  }
}

/*!
************************************************************************
* \brief
*    Prepare first coding layer.
************************************************************************
*/
static void prepare_first_layer(ImageParameters *p_img, int curr_frame_to_code)
{
  p_img->number ++;
  p_img->gop_number = (p_img->number - p_img->start_frame_no);
  p_img->frm_number = p_img->number;

  p_img->frm_no_in_file = CalculateFrameNumber(p_img);

  if (params->last_frame != 0 && (params->last_frame == p_img->frm_no_in_file) && params->HierarchicalCoding)
  { 
    int numberBFrames = params->NumberBFrames;
    params->HierarchicalCoding = imin( 2, params->HierarchicalCoding);    
    params->NumberBFrames = params->no_frames - curr_frame_to_code - 1;    

    clear_gop_structure();
    init_gop_structure(params);
    create_hierarchy(params);

    params->NumberBFrames = numberBFrames;
  }
  SetImgType(p_img, 0);

  init_frame_params(p_img);

  //Rate control
  if (params->RCEnable && p_img->type == I_SLICE)
    rc_init_gop_params();

  // which layer does the image belong to?
  p_img->layer = (p_img->gop_number % (params->NumFramesInELSubSeq + 1)) ? 0 : 1;
}

/*!
************************************************************************
* \brief
*    Prepare second coding layer.
************************************************************************
*/
static void prepare_second_layer(ImageParameters *p_img, int enh_frame_to_code)
{
  p_img->layer = (params->NumFramesInELSubSeq == 0) ? 0 : 1;      
  SetImgType(p_img, enh_frame_to_code);

  if ((p_img->gop_number > 0) && (params->EnableIDRGOP == 0 || p_img->idr_gop_number)) // B-frame(s) to encode
  {
    if (params->HierarchicalCoding)
    {
      p_img->nal_reference_idc = gop_structure[enh_frame_to_code - 1].reference_idc;
      set_poc(p_img, (double)(1 + gop_structure[enh_frame_to_code - 1].display_no));

      if (p_img->gop_number && (p_img->gop_number <= params->intra_delay))
      {
        if (enh_frame_to_code == 1)
          p_img->delta_pic_order_cnt[0] = p_img->toppoc - 2*(p_img->start_tr_gop  + (params->intra_delay - p_img->gop_number)*(p_img->base_dist));
        else
          p_img->delta_pic_order_cnt[0] = p_img->toppoc - 2*(p_img->start_tr_gop  + (params->intra_delay - p_img->gop_number)*(p_img->base_dist) 
          + (int) (2.0 * p_img->frame_interval * (double) (1 + gop_structure[enh_frame_to_code - 2].display_no)));
      }
      else
      {
        if (enh_frame_to_code == 1)
          p_img->delta_pic_order_cnt[0] = p_img->toppoc - 2*(p_img->start_tr_gop  + (p_img->frm_number - p_img->lastIDRnumber)*(p_img->base_dist));
        else
          p_img->delta_pic_order_cnt[0] = p_img->toppoc - 2*(p_img->start_tr_gop  + (p_img->frm_number - p_img->lastIDRnumber - 1)*(p_img->base_dist) 
          + (int) (2.0 * p_img->frame_interval * (double) (1+ gop_structure[enh_frame_to_code - 2].display_no)));
        }
    }
    else
    {
      p_img->nal_reference_idc = (params->BRefPictures == 1 ) ? NALU_PRIORITY_LOW : NALU_PRIORITY_DISPOSABLE;
      set_poc(p_img, (double)enh_frame_to_code);

      //the following is sent in the slice header
      if (params->BRefPictures != 1)
      {
        p_img->delta_pic_order_cnt[0]= 2 * (enh_frame_to_code - 1);
      }
      else
      {
        p_img->delta_pic_order_cnt[0]= -2;
      }
    }

    p_img->delta_pic_order_cnt[1]= 0;

    if ((params->PicInterlace==FRAME_CODING)&&(params->MbInterlace==FRAME_CODING))
      p_img->bottompoc = p_img->toppoc;
    else
      p_img->bottompoc = p_img->toppoc + 1;

    p_img->framepoc = imin (p_img->toppoc, p_img->bottompoc);
    p_img->frm_no_in_file = CalculateFrameNumber(p_img);
  }
}

/*!
 ***********************************************************************
 * \brief
 *    Encode a sequence
 ***********************************************************************
 */
static void encode_sequence(ImageParameters *img)
{
  int HierarchicalCoding = params->HierarchicalCoding;
  int NumberBFrames = params->NumberBFrames;
  int jumpd = params->jumpd;
  int curr_frame_to_code = 0;
  int enh_frame_to_code = 0;

  for (curr_frame_to_code = 0; curr_frame_to_code < params->no_frames; curr_frame_to_code++)
  {     
    // Update frame_num counter
    if (img->last_ref_idc == 1)
    {
      img->frame_num++;
      img->frame_num %= max_frame_num;
    }
   
    // Read explicit sequence coding information
    if (params->ExplicitSeqCoding)
    {
      ExpFrameInfo *info = &expSeq->info[curr_frame_to_code % expSeq->no_frames];
      ReadExplicitSeqFile(expSeq, curr_frame_to_code);
      ExplicitUpdateImgParams(info, img);
      img->b_frame_to_code = 0;
    }
    else
    {
      enh_frame_to_code = determine_coding_level(img, params, curr_frame_to_code);
      img->b_frame_to_code = enh_frame_to_code;

      if (enh_frame_to_code == 0) 
        prepare_first_layer(img, curr_frame_to_code);
      else 
      {
        prepare_second_layer(img, enh_frame_to_code);
      }
    }

    // redundant frame initialization and allocation
    if (params->redundant_pic_flag)
    {
      Init_redundant_frame();
      Set_redundant_frame();
    }

    encode_one_frame(img); // encode one frame;

    img->last_ref_idc = img->nal_reference_idc ? 1 : 0;

    // if key frame is encoded, encode one redundant frame
    if (params->redundant_pic_flag && key_frame)
    {
      encode_one_redundant_frame();
    }

    if (img->type == I_SLICE && params->EnableOpenGOP)
      img->last_valid_reference = img->ThisPOC;

    if (img->currentPicture->idr_flag)
    {
      img->idr_gop_number = 0;
    }
    else
      img->idr_gop_number ++;

    if (params->ReportFrameStats)
      report_frame_statistic();
  }

  params->HierarchicalCoding = HierarchicalCoding;
  params->NumberBFrames      = NumberBFrames;
  params->jumpd = jumpd;
}


/*!
 ***********************************************************************
 * \brief
 *    Free memory allocated for the encoder
 ***********************************************************************
 */

void free_encoder_memory(ImageParameters *img)
{
  int nplane;
  terminate_sequence();

  flush_dpb();

  CloseFiles(&params->input_file1);

  if (-1 != p_dec)
    close(p_dec);
  if (p_trace)
    fclose(p_trace);

  Clear_Motion_Search_Module ();

  RandomIntraUninit();
  FmoUninit();

  if (params->HierarchicalCoding)
    clear_gop_structure ();

  // free structure for rd-opt. mode decision
  clear_rdopt (params);

#ifdef _LEAKYBUCKET_
  calc_buffer();
#endif

  // report everything
  report(img, params, stats);

#ifdef _LEAKYBUCKET_
  if (Bit_Buffer != NULL)
  {
    free(Bit_Buffer);
    Bit_Buffer = NULL;
  }
#endif

  free_dpb();

  if( IS_INDEPENDENT(params) )
  {
    for( nplane=0; nplane<MAX_PLANE; nplane++ )
    {
      free_colocated(Co_located_JV[nplane]);
    }
  }
  else
  {
    free_colocated(Co_located);
  }

  uninit_out_buffer();

  free_global_buffers();

  // free image mem
  free_img ();
  free_context_memory ();
  FreeParameterSets();
   if (params->ExplicitSeqCoding)
     CloseExplicitSeqFile();
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
  //the following should probably go in sequence parameters
  // frame poc's increase by 2, field poc's by 1

  img->pic_order_cnt_type=params->pic_order_cnt_type;

  img->delta_pic_order_always_zero_flag = FALSE;
  img->num_ref_frames_in_pic_order_cnt_cycle= 1;

  if (params->BRefPictures == 1)
  {
    img->offset_for_non_ref_pic  =  0;
    img->offset_for_ref_frame[0] =  2;
  }
  else
  {
    img->offset_for_non_ref_pic  = -2*(params->NumberBFrames);
    img->offset_for_ref_frame[0] =  2*(params->NumberBFrames + 1);
  }

  if ((params->PicInterlace==FRAME_CODING) && (params->MbInterlace==FRAME_CODING))
  {
    img->offset_for_top_to_bottom_field = 0;
    img->pic_order_present_flag = FALSE;
    img->delta_pic_order_cnt_bottom = 0;
  }
  else
  {
    img->offset_for_top_to_bottom_field = 1;
    img->pic_order_present_flag = TRUE;
    img->delta_pic_order_cnt_bottom = 1;
  }
}


static inline int alloc_rddata(ImageParameters *img, RD_DATA *rd_data)
{
  int alloc_size = 0;

  alloc_size += get_mem2Dpel(&(rd_data->rec_mbY), MB_BLOCK_SIZE, MB_BLOCK_SIZE);
  alloc_size += get_mem3Dpel(&(rd_data->rec_mb_cr), 2, MB_BLOCK_SIZE, MB_BLOCK_SIZE);

  alloc_size += get_mem_ACcoeff (&(rd_data->cofAC));
  alloc_size += get_mem_DCcoeff (&(rd_data->cofDC));  

  if (!params->IntraProfile)
  {          
    alloc_size += get_mem_mv (&(rd_data->all_mv));
  }
  
  alloc_size += get_mem2D((byte***)&(rd_data->ipredmode), img->height_blk, img->width_blk);
  alloc_size += get_mem3D((byte****)&(rd_data->refar), 2, 4, 4);

  return alloc_size;
}

static inline void free_rddata(RD_DATA *rd_data)
{
  free_mem3D((byte***) rd_data->refar);
  free_mem2D((byte**) rd_data->ipredmode);

  if (!params->IntraProfile)
  {          
    free_mem_mv (rd_data->all_mv);
  }
  free_mem_DCcoeff (rd_data->cofDC);
  free_mem_ACcoeff (rd_data->cofAC);  

  free_mem3Dpel(rd_data->rec_mb_cr);
  free_mem2Dpel(rd_data->rec_mbY);
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
void CAVLC_init(void)
{
  memset(&img->nz_coeff[0][0][0], 0, img->PicSizeInMbs * 4 * (4 + img->num_blk8x8_uv)* sizeof(int));
}


/*!
 ***********************************************************************
 * \brief
 *    Initializes the Image structure with appropriate parameters.
 * \par Input:
 *    Input Parameters struct inp_par *inp
 * \par  Output:
 *    Image Parameters ImageParameters *img
 ***********************************************************************
 */
static void init_img( ImageParameters *img, InputParameters *params)
{
  int i, j;
  int byte_abs_range;

  static int mb_width_cr[4] = {0,8, 8,16};
  static int mb_height_cr[4]= {0,8,16,16};

  img->number         = -1;
  img->start_frame_no = 0;
  img->gop_number     = (img->number - img->start_frame_no);
  img->start_tr_gop   = 0;

  img->last_idr_number = 0;
  // Color format
  img->yuv_format  = params->output.yuv_format;
  img->P444_joined = (img->yuv_format == YUV444 && !IS_INDEPENDENT(params));  

  //pel bitdepth init
  img->bitdepth_luma            = params->output.bit_depth[0];
  img->bitdepth_scale[0]        = 1 << (img->bitdepth_luma - 8);
  img->bitdepth_lambda_scale    = 2 * (img->bitdepth_luma - 8);
  img->bitdepth_luma_qp_scale   = 3 *  img->bitdepth_lambda_scale;
  img->dc_pred_value_comp[0]    =  1<<(img->bitdepth_luma - 1);
  img->max_imgpel_value_comp[0] = (1<<img->bitdepth_luma) - 1;
  img->max_imgpel_value_comp_sq[0] = img->max_imgpel_value_comp[0] * img->max_imgpel_value_comp[0];

  img->dc_pred_value            = img->dc_pred_value_comp[0]; // set defaults
  img->max_imgpel_value         = img->max_imgpel_value_comp[0];
  img->mb_size[0][0]            = img->mb_size[0][1] = MB_BLOCK_SIZE;

  // Initialization for RC QP parameters (could be placed in ratectl.c)
  img->RCMinQP                = params->RCMinQP[P_SLICE];
  img->RCMaxQP                = params->RCMaxQP[P_SLICE];

  // Set current residue & prediction array pointers
  img->curr_res = img->mb_rres[0];
  img->curr_prd = img->mb_pred[0];
  if (active_sps->profile_idc == 66 || active_sps->profile_idc == 77 || active_sps->profile_idc == 88)
    img->min_IPCM_value = 1;  // See Annex A for restriction in pcm sample values for pre FRExt profiles
  else
    img->min_IPCM_value = 0;

  if (img->yuv_format != YUV400)
  {
    img->bitdepth_chroma          = params->output.bit_depth[1];
    img->bitdepth_scale[1]        = 1 << (img->bitdepth_chroma - 8);
    img->dc_pred_value_comp[1]    = 1<<(img->bitdepth_chroma - 1);
    img->dc_pred_value_comp[2]    = img->dc_pred_value_comp[1];
    img->max_imgpel_value_comp[1] = (1<<img->bitdepth_chroma) - 1;
    img->max_imgpel_value_comp[2] = img->max_imgpel_value_comp[1];
    img->max_imgpel_value_comp_sq[1] = img->max_imgpel_value_comp[1] * img->max_imgpel_value_comp[1];
    img->max_imgpel_value_comp_sq[2] = img->max_imgpel_value_comp[2] * img->max_imgpel_value_comp[2];
    img->num_blk8x8_uv            = (1<<img->yuv_format)&(~(0x1));
    img->num_cdc_coeff            = img->num_blk8x8_uv << 1;

    img->mb_size[1][0] = img->mb_size[2][0] = img->mb_cr_size_x = (img->yuv_format == YUV420 || img->yuv_format == YUV422) ? 8 : 16;
    img->mb_size[1][1] = img->mb_size[2][1] = img->mb_cr_size_y = (img->yuv_format == YUV444 || img->yuv_format == YUV422) ? 16 : 8;

    img->bitdepth_chroma_qp_scale = 6*(img->bitdepth_chroma - 8);

    img->chroma_qp_offset[0] = active_pps->cb_qp_index_offset;
    img->chroma_qp_offset[1] = active_pps->cr_qp_index_offset;
  }
  else
  {
    img->bitdepth_chroma     = 0;
    img->bitdepth_scale[1]   = 0;
    img->max_imgpel_value_comp[1] = 0;
    img->max_imgpel_value_comp[2] = img->max_imgpel_value_comp[1];
    img->max_imgpel_value_comp_sq[1] = img->max_imgpel_value_comp[1] * img->max_imgpel_value_comp[1];
    img->max_imgpel_value_comp_sq[2] = img->max_imgpel_value_comp[2] * img->max_imgpel_value_comp[2];
    img->num_blk8x8_uv       = 0;
    img->num_cdc_coeff       = 0;
    img->mb_size[1][0] = img->mb_size[2][0] = img->mb_cr_size_x = 0;
    img->mb_size[1][1] = img->mb_size[2][1] = img->mb_cr_size_y = 0;

    img->bitdepth_chroma_qp_scale = 0;
    img->bitdepth_chroma_qp_scale = 0;

    img->chroma_qp_offset[0] = 0;
    img->chroma_qp_offset[1] = 0;
  }  

  //img->pic_unit_size_on_disk = (imax(img->bitdepth_luma , img->bitdepth_chroma) > 8) ? 16 : 8;
  img->pic_unit_size_on_disk = (imax(params->source.bit_depth[0], params->source.bit_depth[1]) > 8) ? 16 : 8;
  img->out_unit_size_on_disk = (imax(params->output.bit_depth[0], params->output.bit_depth[1]) > 8) ? 16 : 8;

  img->max_bitCount =  128 + 256 * img->bitdepth_luma + 2 * img->mb_cr_size_y * img->mb_cr_size_x * img->bitdepth_chroma;
  //img->max_bitCount =  (128 + 256 * img->bitdepth_luma + 2 *img->mb_cr_size_y * img->mb_cr_size_x * img->bitdepth_chroma)*2;

  img->max_qp_delta = (25 + (img->bitdepth_luma_qp_scale>>1));
  img->min_qp_delta = img->max_qp_delta + 1;

  img->num_ref_frames = active_sps->num_ref_frames;
  img->max_num_references   = active_sps->frame_mbs_only_flag ? active_sps->num_ref_frames : 2 * active_sps->num_ref_frames;

  img->buf_cycle = params->num_ref_frames;
  img->base_dist = params->jumpd + 1;  

  // Intra/IDR related parameters
  img->lastIntraNumber = 0;
  img->lastINTRA       = 0;
  img->lastIDRnumber   = 0;
  img->last_ref_idc    = 0;
  img->idr_refresh     = 0;
  img->idr_gop_number  = 0;
  img->rewind_frame    = 0;

  img->DeblockCall     = 0;
  img->framerate       = (float) params->output.frame_rate;   // The basic frame rate (of the original sequence)

  // Allocate proper memory space for different parameters (i.e. MVs, coefficients, etc)
  if (!params->IntraProfile)
  {
    get_mem_mv (&(img->all_mv));

    if (params->BiPredMotionEstimation)
    {
      get_mem_bipred_mv(&(img->bipred_mv)); 
    }
  }

  get_mem_ACcoeff (&(img->cofAC));
  get_mem_DCcoeff (&(img->cofDC));

  if (params->AdaptiveRounding)
  {
    if (img->yuv_format != 0)
    {
      get_mem4Dint(&(img->ARCofAdj4x4), 3, MAXMODE, MB_BLOCK_SIZE, MB_BLOCK_SIZE); //all modes
      get_mem4Dint(&(img->ARCofAdj8x8), img->P444_joined ? 3 : 1, MAXMODE, MB_BLOCK_SIZE, MB_BLOCK_SIZE); //modes 0, 1, 2, 3, P8x8
    }     
    else
    {
      get_mem4Dint(&(img->ARCofAdj4x4), 1, MAXMODE, MB_BLOCK_SIZE, MB_BLOCK_SIZE); //all modes
      get_mem4Dint(&(img->ARCofAdj8x8), 1, MAXMODE, MB_BLOCK_SIZE, MB_BLOCK_SIZE); //modes 0, 1, 2, 3, P8x8
    }
  }

  byte_abs_range = (imax(img->max_imgpel_value_comp[0], img->max_imgpel_value_comp[1]) + 1) * 2;

  if ((img->quad = (int*)calloc (byte_abs_range, sizeof(int))) == NULL)
    no_mem_exit ("init_img: img->quad");

  img->quad += byte_abs_range/2;
  for (i=0; i < byte_abs_range/2; ++i)
  {
    img->quad[i] = img->quad[-i] = i * i;
  }

  img->width         = (params->output.width  + img->auto_crop_right);
  img->height        = (params->output.height + img->auto_crop_bottom);
  img->width_blk     = img->width  / BLOCK_SIZE;
  img->height_blk    = img->height / BLOCK_SIZE;
  img->width_padded  = img->width  + 2 * IMG_PAD_SIZE;
  img->height_padded = img->height + 2 * IMG_PAD_SIZE;

  if (img->yuv_format != YUV400)
  {
    img->width_cr = img->width  * mb_width_cr [img->yuv_format] / 16;
    img->height_cr= img->height * mb_height_cr[img->yuv_format] / 16;
  }
  else
  {
    img->width_cr = 0;
    img->height_cr= 0;
  }

  img->height_cr_frame = img->height_cr;

  img->size = img->width * img->height;
  img->size_cr = img->width_cr * img->height_cr;

  img->PicWidthInMbs    = img->width  / MB_BLOCK_SIZE;
  img->FrameHeightInMbs = img->height / MB_BLOCK_SIZE;
  img->FrameSizeInMbs   = img->PicWidthInMbs * img->FrameHeightInMbs;

  img->PicHeightInMapUnits = ( active_sps->frame_mbs_only_flag ? img->FrameHeightInMbs : img->FrameHeightInMbs/2 );

  if(params->MbInterlace)
  {
    alloc_rddata(img, &rddata_top_frame_mb);
    alloc_rddata(img, &rddata_bot_frame_mb);
    if ( params->MbInterlace != FRAME_MB_PAIR_CODING )
    {
      alloc_rddata(img, &rddata_top_field_mb);
      alloc_rddata(img, &rddata_bot_field_mb);
    }
  }

  if (params->UseRDOQuant && params->RDOQ_QP_Num > 1)
  {
    alloc_rddata(img, &rddata_trellis_curr);
    alloc_rddata(img, &rddata_trellis_best);

    if (!params->IntraProfile)
    {          
      if (params->Transform8x8Mode && params->RDOQ_CP_MV)
      {
        get_mem5Dshort(&tmp_mv8, 2, img->max_num_references, 4, 4, 2);
        get_mem3Dint  (&motion_cost8, 2, img->max_num_references, 4);
      }
    }
  }

  if ((img->b8x8info = (Block8x8Info *) calloc(1, sizeof(Block8x8Info))) == NULL)
     no_mem_exit("init_img: img->block8x8info");

  if( IS_INDEPENDENT(params) )
  {
    for( i = 0; i < MAX_PLANE; i++ ){
      if ((img->mb_data_JV[i] = (Macroblock *) calloc(img->FrameSizeInMbs,sizeof(Macroblock))) == NULL)
        no_mem_exit("init_img: img->mb_data_JV");
    }
    img->mb_data = NULL;
  }
  else
  {
    if ((img->mb_data = (Macroblock *) calloc(img->FrameSizeInMbs, sizeof(Macroblock))) == NULL)
      no_mem_exit("init_img: img->mb_data");
  }

  if (params->UseConstrainedIntraPred)
  {
    if ((img->intra_block = (int*)calloc(img->FrameSizeInMbs, sizeof(int))) == NULL)
      no_mem_exit("init_img: img->intra_block");
  }

  if (params->CtxAdptLagrangeMult == 1)
  {
    if ((mb16x16_cost_frame = (double*)calloc(img->FrameSizeInMbs, sizeof(double))) == NULL)
    {
      no_mem_exit("init mb16x16_cost_frame");
    }
  }
  get_mem2D((byte***)&(img->ipredmode), img->height_blk, img->width_blk);        //need two extra rows at right and bottom
  get_mem2D((byte***)&(img->ipredmode8x8), img->height_blk, img->width_blk);     // help storage for ipredmode 8x8, inserted by YV
  memset(&(img->ipredmode[0][0])   , -1, img->height_blk * img->width_blk *sizeof(char));
  memset(&(img->ipredmode8x8[0][0]), -1, img->height_blk * img->width_blk *sizeof(char));


  if (params->MbInterlace)
  {
    get_mem2D((byte***)&(rddata_bot_frame_mb.ipredmode), img->height_blk, img->width_blk);
    get_mem2D((byte***)&(rddata_top_field_mb.ipredmode), img->height_blk, img->width_blk);
    get_mem2D((byte***)&(rddata_bot_field_mb.ipredmode), img->height_blk, img->width_blk);
  }

  // CAVLC mem
  get_mem3Dint(&(img->nz_coeff), img->FrameSizeInMbs, 4, 4+img->num_blk8x8_uv);
  
  get_mem2Dolm     (&(img->lambda), 10, 52 + img->bitdepth_luma_qp_scale, img->bitdepth_luma_qp_scale);
  get_mem2Dodouble (&(img->lambda_md), 10, 52 + img->bitdepth_luma_qp_scale, img->bitdepth_luma_qp_scale);
  get_mem3Dodouble (&(img->lambda_me), 10, 52 + img->bitdepth_luma_qp_scale, 3, img->bitdepth_luma_qp_scale);
  get_mem3Doint    (&(img->lambda_mf), 10, 52 + img->bitdepth_luma_qp_scale, 3, img->bitdepth_luma_qp_scale);

  if (params->CtxAdptLagrangeMult == 1)
  {
    get_mem2Dodouble(&(img->lambda_mf_factor), 10, 52 + img->bitdepth_luma_qp_scale, img->bitdepth_luma_qp_scale);
  }

  img->b_frame_to_code = 0;
  img->GopLevels = (params->NumberBFrames) ? 1 : 0;
  img->mb_y_upd=0;

  RandomIntraInit (img->PicWidthInMbs, img->FrameHeightInMbs, params->RandomIntraMBRefresh);

  InitSEIMessages(); 

  initInput(&params->source, &params->output);

  // Allocate I/O Frame memory
  AllocateFrameMemory(img, params, params->source.size);

  // Initialize filtering parameters. If sending parameters, the offsets are
  // multiplied by 2 since inputs are taken in "div 2" format.
  // If not sending parameters, all fields are cleared
  if (params->DFSendParameters)
  {
    for (j = 0; j < 2; j++)
    {
      for (i = 0; i < NUM_SLICE_TYPES; i++)
      {
        params->DFAlpha[j][i] <<= 1;
        params->DFBeta [j][i] <<= 1;
      }
    }
  }
  else
  {
    for (j = 0; j < 2; j++)
    {
      for (i = 0; i < NUM_SLICE_TYPES; i++)
      {
        params->DFDisableIdc[j][i] = 0;
        params->DFAlpha     [j][i] = 0;
        params->DFBeta      [j][i] = 0;
      }
    }
  }

  img->ChromaArrayType = params->separate_colour_plane_flag ? 0 : params->output.yuv_format;
  img->colour_plane_id = 0;

  if (params->RDPictureDecision)
    img->frm_iter = 3;
  else
    img->frm_iter = 1;

  img->frame_interval = (double) (params->frame_skip + 1);
}


/*!
 ***********************************************************************
 * \brief
 *    Free the Image structures
 * \par Input:
 *    Image Parameters ImageParameters *img
 ***********************************************************************
 */
void free_img ()
{
  // Delete Frame memory 
  DeleteFrameMemory();

  CloseSEIMessages(); 
  if (!params->IntraProfile)
  {
    free_mem_mv (img->all_mv);

    if (params->BiPredMotionEstimation)
    {
      free_mem_bipred_mv (img->bipred_mv);
    }
  }

  free_mem_ACcoeff (img->cofAC);
  free_mem_DCcoeff (img->cofDC);

  if (params->AdaptiveRounding)
  {
    free_mem4Dint(img->ARCofAdj4x4);
    free_mem4Dint(img->ARCofAdj8x8);
  }


  if (params->MbInterlace)
  {
    free_rddata(&rddata_top_frame_mb);
    free_rddata(&rddata_bot_frame_mb);

    if ( params->MbInterlace != FRAME_MB_PAIR_CODING )
    {
      free_rddata(&rddata_top_field_mb);
      free_rddata(&rddata_bot_field_mb);
    }
  }

  if (params->UseRDOQuant && params->RDOQ_QP_Num > 1)
  {
    free_rddata(&rddata_trellis_curr);
    free_rddata(&rddata_trellis_best);

    if (!params->IntraProfile)
    {    
      if (params->Transform8x8Mode && params->RDOQ_CP_MV)
      {
        free_mem5Dshort(tmp_mv8);
        free_mem3Dint(motion_cost8);
      }
    }
  }

  free (img->quad - (imax(img->max_imgpel_value_comp[0],img->max_imgpel_value_comp[1]) + 1));
  
  free_mem3Dpel(img->mb_pred  );
  free_mem4Dpel(img->mpr_16x16);
  free_mem4Dpel(img->mpr_8x8  );
  free_mem4Dpel(img->mpr_4x4  );
  free_mem3Dint(img->mb_ores  );
  free_mem3Dint(img->mb_rres  );
  free (img);
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
int init_orig_buffers(ImageData *imgData)
{
  int memory_size = 0;
  int nplane;

  // allocate memory for reference frame buffers: imgData->frm_data
  if( IS_INDEPENDENT(params) )
  {
    for( nplane=0; nplane<MAX_PLANE; nplane++ )
    {
      imgData->format = params->output;
      imgData->format.width = img->width;
      imgData->format.height = img->height;
      imgData->format.width_cr = img->width_cr;
      imgData->format.height_cr = img->height_cr;
      memory_size += get_mem2Dpel(&(imgData->frm_data[nplane]), img->height, img->width);
    }
  }
  else
  {
    //imgData->format = params->input_file1.format;
    imgData->format           = params->output;
    imgData->format.width     = img->width;    
    imgData->format.height    = img->height;
    imgData->format.width_cr  = img->width_cr;
    imgData->format.height_cr = img->height_cr;
    memory_size += get_mem2Dpel(&(imgData->frm_data[0]), img->height, img->width);
    
    if (img->yuv_format != YUV400)
    {
      int i, j, k;
      memory_size += get_mem2Dpel(&(imgData->frm_data[1]), img->height_cr, img->width_cr);
      memory_size += get_mem2Dpel(&(imgData->frm_data[2]), img->height_cr, img->width_cr);

      if (sizeof(imgpel) == sizeof(unsigned char))
      {
        for (k = 1; k < 3; k++)
          memset(&(imgData->frm_data[k][0][0]), 128, img->height_cr * img->width_cr * sizeof(imgpel));
      }
      else
      {
        for (k = 1; k < 3; k++)
          for (j = 0; j < img->height_cr; j++)
            for (i = 0; i < img->width_cr; i++)
              imgData->frm_data[k][j][i] = 128;
      }
    }
  }

  if (!active_sps->frame_mbs_only_flag)
  {
    // allocate memory for field reference frame buffers
    memory_size += init_top_bot_planes(imgData->frm_data[0], img->height, img->width, &(imgData->top_data[0]), &(imgData->bot_data[0]));

    if (img->yuv_format != YUV400)
    {

      memory_size += 4*(sizeof(imgpel**));

      memory_size += init_top_bot_planes(imgData->frm_data[1], img->height_cr, img->width_cr, &(imgData->top_data[1]), &(imgData->bot_data[1]));
      memory_size += init_top_bot_planes(imgData->frm_data[2], img->height_cr, img->width_cr, &(imgData->top_data[2]), &(imgData->bot_data[2]));
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
 *    Input Parameters struct inp_par *inp,                            \n
 *    Image Parameters ImageParameters *img
 * \return Number of allocated bytes
 ************************************************************************
 */
int init_global_buffers(void)
{
  int j, memory_size=0;

  if ((enc_frame_picture = (StorablePicture**)malloc(6 * sizeof(StorablePicture*))) == NULL)
    no_mem_exit("init_global_buffers: *enc_frame_picture");

  for (j = 0; j < 6; j++)
    enc_frame_picture[j] = NULL;

  if ((enc_field_picture = (StorablePicture**)malloc(2 * sizeof(StorablePicture*))) == NULL)
    no_mem_exit("init_global_buffers: *enc_field_picture");

  for (j = 0; j < 2; j++)
    enc_field_picture[j] = NULL;

  if ((frame_pic = (Picture**)malloc(img->frm_iter * sizeof(Picture*))) == NULL)
    no_mem_exit("init_global_buffers: *frame_pic");

  for (j = 0; j < img->frm_iter; j++)
    frame_pic[j] = malloc_picture();

  if (params->si_frame_indicator || params->sp_periodicity)
  {
    si_frame_indicator=0; //indicates whether the frame is SP or SI
    number_sp2_frames=0;

    frame_pic_si = malloc_picture();//picture buffer for the encoded SI picture
    //allocation of lrec and lrec_uv for SI picture
    get_mem2Dint (&lrec, img->height, img->width);
    get_mem3Dint (&lrec_uv, 2, img->height, img->width);
  }

  // Allocate memory for field picture coding
  if (params->PicInterlace != FRAME_CODING)
  { 
    if ((field_pic = (Picture**)malloc(2 * sizeof(Picture*))) == NULL)
      no_mem_exit("init_global_buffers: *field_pic");

    for (j = 0; j < 2; j++)
      field_pic[j] = malloc_picture();
  }

  // Init memory data for input & encoded images
  memory_size += init_orig_buffers(&imgData);
  memory_size += init_orig_buffers(&imgData1);
  
  memory_size += get_mem2Dint(&PicPos, img->FrameSizeInMbs + 1, 2);

  for (j = 0; j < (int) img->FrameSizeInMbs + 1; j++)
  {
    PicPos[j][0] = (j % img->PicWidthInMbs);
    PicPos[j][1] = (j / img->PicWidthInMbs);
  }

  if (params->WeightedPrediction || params->WeightedBiprediction || params->GenerateMultiplePPS)
  {
    // Currently only use up to 32 references. Need to use different indicator such as maximum num of references in list
    memory_size += get_mem3Dint(&wp_weight, 6, MAX_REFERENCE_PICTURES, 3);
    memory_size += get_mem3Dint(&wp_offset, 6, MAX_REFERENCE_PICTURES, 3);

    memory_size += get_mem4Dint(&wbp_weight, 6, MAX_REFERENCE_PICTURES, MAX_REFERENCE_PICTURES, 3);
  }

  // allocate memory for reference frames of each block: refFrArr

  if ((params->NumberBFrames != 0) || (params->BRefPictures > 0) || (params->ProfileIDC != 66))
  {
    memory_size += get_mem3D((byte ****)(void*)&direct_ref_idx, 2, img->height_blk, img->width_blk);
    memory_size += get_mem2D((byte ***) (void*)&direct_pdir      , img->height_blk, img->width_blk);
  }

  if (params->rdopt == 3)
  {
    memory_size += allocate_errdo_mem(params);
  }

  if (params->RestrictRef)
  {
    memory_size += get_mem2D(&pixel_map,   img->height,   img->width);
    memory_size += get_mem2D(&refresh_map, img->height/8, img->width/8);
  }

  if (!active_sps->frame_mbs_only_flag)
  {
    memory_size += get_mem2Dpel(&imgY_com, img->height, img->width);

    if (img->yuv_format != YUV400)
    {
      memory_size += get_mem3Dpel(&imgUV_com, 2, img->height_cr, img->width_cr);
    }
  }

  // allocate and set memory relating to motion estimation
  if (!params->IntraProfile)
  {  
    if (params->SearchMode == UM_HEX)
    {
      memory_size += UMHEX_get_mem();
    }
    else if (params->SearchMode == UM_HEX_SIMPLE)
    {
      smpUMHEX_init();
      memory_size += smpUMHEX_get_mem();
    }
    else if (params->SearchMode == EPZS)
      memory_size += EPZSInit(params, img);
  }

  if (params->RCEnable)
    rc_allocate_memory();

  if (params->redundant_pic_flag)
  {
    memory_size += get_mem2Dpel(&imgY_tmp, img->height, img->width);
    memory_size += get_mem2Dpel(&imgUV_tmp[0], img->height_cr, img->width_cr);
    memory_size += get_mem2Dpel(&imgUV_tmp[1], img->height_cr, img->width_cr);
  }

  memory_size += get_mem2Dint (&imgY_sub_tmp, img->height_padded, img->width_padded);

  if ( params->ChromaMCBuffer )
    chroma_mc_setup();

  img_padded_size_x       = (img->width + 2 * IMG_PAD_SIZE);
  img_padded_size_x2      = (img_padded_size_x << 1);
  img_padded_size_x4      = (img_padded_size_x << 2);
  img_padded_size_x_m8    = (img_padded_size_x - 8);
  img_padded_size_x_m8x8  = (img_padded_size_x - BLOCK_SIZE_8x8);
  img_padded_size_x_m4x4  = (img_padded_size_x - BLOCK_SIZE);
  img_cr_padded_size_x    = (img->width_cr + 2 * img_pad_size_uv_x);
  img_cr_padded_size_x2   = (img_cr_padded_size_x << 1);
  img_cr_padded_size_x4   = (img_cr_padded_size_x << 2);
  img_cr_padded_size_x_m8 = (img_cr_padded_size_x - 8);

  // RGB images for distortion calculation
  // Recommended to do this allocation (and de-allocation) in 
  // the appropriate file instead of here.
  if(params->DistortionYUVtoRGB)
  {
    memory_size += create_RGB_memory(img);
  }

  pWPX = NULL;
  if ( params->WPMCPrecision )
  {
    wpxInitWPXObject();
  }

  return memory_size;
}


/*!
 ************************************************************************
 * \brief
 *    Free allocated memory of original picture buffers
 ************************************************************************
 */
void free_orig_planes(ImageData *imgData)
{
  if( IS_INDEPENDENT(params) )
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
    
    if (img->yuv_format != YUV400)
    {
      free_mem2Dpel(imgData->frm_data[1]);
      free_mem2Dpel(imgData->frm_data[2]);
    }
  }

  if (!active_sps->frame_mbs_only_flag)
  {
    free_top_bot_planes(imgData->top_data[0], imgData->bot_data[0]);

    if (img->yuv_format != YUV400)
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
 *    Input Parameters struct inp_par *inp,                             \n
 *    Image Parameters ImageParameters *img
 * \par Output:
 *    none
 ************************************************************************
 */
void free_global_buffers(void)
{
  int  i,j;

  if (enc_frame_picture)
    free (enc_frame_picture);
  if (frame_pic)
  {
    for (j = 0; j < img->frm_iter; j++)
    {
      if (frame_pic[j])
        free_picture (frame_pic[j]);
    }
    free (frame_pic);
  }

  if (enc_field_picture)
    free (enc_field_picture);
  if (field_pic)
  {
    for (j = 0; j < 2; j++)
    {
      if (field_pic[j])
        free_picture (field_pic[j]);
    }
    free (field_pic);
  }

  // Deallocation of SI picture related memory
  if (params->si_frame_indicator || params->sp_periodicity)
  {
    free_picture (frame_pic_si);
    //deallocation of lrec and lrec_uv for SI frames
    free_mem2Dint (lrec);
    free_mem3Dint (lrec_uv);
  }

  free_orig_planes(&imgData);
  free_orig_planes(&imgData1);

  // free lookup memory which helps avoid divides with PicWidthInMbs
  free_mem2Dint(PicPos);
  // Free Qmatrices and offsets
  free_QMatrix();
  free_QOffsets();

  if ( params->WPMCPrecision )
  {
    wpxFreeWPXObject();
  }

  if (params->WeightedPrediction || params->WeightedBiprediction || params->GenerateMultiplePPS)
  {
    free_mem3Dint(wp_weight );
    free_mem3Dint(wp_offset );
    free_mem4Dint(wbp_weight);
  }

  if ((stats->NumberBFrames != 0) || (params->BRefPictures > 0)||(params->ProfileIDC != 66))
  {
    free_mem3D((byte ***)direct_ref_idx);
    free_mem2D((byte **) direct_pdir);
  } // end if B frame

  if (imgY_sub_tmp) // free temp quarter pel frame buffers
  {
    free_mem2Dint (imgY_sub_tmp);
    imgY_sub_tmp=NULL;
  }

  // free mem, allocated in init_img()
  // free intra pred mode buffer for blocks
  free_mem2D((byte**)img->ipredmode);
  free_mem2D((byte**)img->ipredmode8x8);
  free(img->b8x8info);

  if( IS_INDEPENDENT(params) )
  {
    for( i=0; i<MAX_PLANE; i++ ){
      free(img->mb_data_JV[i]);
    }
  }
  else
  {
    free(img->mb_data);
  }

  if(params->UseConstrainedIntraPred)
  {
    free (img->intra_block);
  }

  if (params->CtxAdptLagrangeMult == 1)
  {
    free(mb16x16_cost_frame);
  }

  if (params->rdopt == 3)
  {
    free_errdo_mem();
  }

  if (params->RestrictRef)
  {
    free(pixel_map[0]);
    free(pixel_map);
    free(refresh_map[0]);
    free(refresh_map);
  }

  if (!active_sps->frame_mbs_only_flag)
  {
    free_mem2Dpel(imgY_com);
    
    if (img->yuv_format != YUV400)
    {
      free_mem3Dpel(imgUV_com);
    }
  }

  free_mem3Dint(img->nz_coeff);

  free_mem2Dolm     (img->lambda, img->bitdepth_luma_qp_scale);
  free_mem2Dodouble (img->lambda_md, img->bitdepth_luma_qp_scale);
  free_mem3Dodouble (img->lambda_me, 10, 52 + img->bitdepth_luma_qp_scale, img->bitdepth_luma_qp_scale);
  free_mem3Doint    (img->lambda_mf, 10, 52 + img->bitdepth_luma_qp_scale, img->bitdepth_luma_qp_scale);

  if (params->CtxAdptLagrangeMult == 1)
  {
    free_mem2Dodouble(img->lambda_mf_factor, img->bitdepth_luma_qp_scale);
  }

  if (!params->IntraProfile)
  {
    if (params->SearchMode == UM_HEX)
    {
      UMHEX_free_mem();
    }
    else if (params->SearchMode == UM_HEX_SIMPLE)
    {
      smpUMHEX_free_mem();
    }
    else if (params->SearchMode == EPZS)
    {
      EPZSDelete(params);
    }
  }


  if (params->RCEnable)
    rc_free_memory();

  if (params->redundant_pic_flag)
  {
    free_mem2Dpel(imgY_tmp);
    free_mem2Dpel(imgUV_tmp[0]);
    free_mem2Dpel(imgUV_tmp[1]);
  }

  // Again process should be moved into cconv_yuv2rgb.c file for cleanliness
  // These should not be globals but instead only be visible through that code.
  if(params->DistortionYUVtoRGB)
  {
    delete_RGB_memory();
  }
}

/*!
 ************************************************************************
 * \brief
 *    Allocate memory for mv
 * \par Input:
 *    Image Parameters ImageParameters *img                             \n
 *    int****** mv
 * \return memory size in bytes
 ************************************************************************
 */
static int get_mem_mv (short ******* mv)
{
  // LIST, reference, block_type, block_y, block_x, component
  get_mem6Dshort(mv, 2, img->max_num_references, 9, 4, 4, 2);

  return 2 * img->max_num_references * 9 * 4 * 4 * 2 * sizeof(short);
}

/*!
 ************************************************************************
 * \brief
 *    Allocate memory for bipredictive mv 
 * \par Input:
 *    Image Parameters ImageParameters *img                             \n
 *    int****** mv
 * \return memory size in bytes
 ************************************************************************
 */
static int get_mem_bipred_mv (short******** bipred_mv) 
{
  get_mem7Dshort(bipred_mv, 2, 2, img->max_num_references, 9, 4, 4, 2);
  
  return 2 * 2 * img->max_num_references * 9 * 4 * 4 * 2 * sizeof(short);
}

/*!
 ************************************************************************
 * \brief
 *    Free memory from mv
 * \par Input:
 *    int****** mv
 ************************************************************************
 */
static void free_mem_mv (short****** mv)
{
  free_mem6Dshort(mv);
}


/*!
 ************************************************************************
 * \brief
 *    Free memory from mv
 * \par Input:
 *    int****** mv
 ************************************************************************
 */
static void free_mem_bipred_mv (short******* bipred_mv) 
{
  free_mem7Dshort(bipred_mv);
}

/*!
 ************************************************************************
 * \brief
 *    Allocate memory for AC coefficients
 ************************************************************************
 */
int get_mem_ACcoeff (int***** cofAC)
{
  int num_blk8x8 = BLOCK_SIZE + img->num_blk8x8_uv;
  
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
static void SetLevelIndices(void)
{
  switch(active_sps->level_idc)
  {
  case 9:
    img->LevelIndex=1;
    break;
  case 10:
    img->LevelIndex=0;
    break;
  case 11:
    if (!IS_FREXT_PROFILE(active_sps->profile_idc) && (active_sps->constrained_set3_flag == 0))
      img->LevelIndex=2;
    else
      img->LevelIndex=1;
    break;
  case 12:
    img->LevelIndex=3;
    break;
  case 13:
    img->LevelIndex=4;
    break;
  case 20:
    img->LevelIndex=5;
    break;
  case 21:
    img->LevelIndex=6;
    break;
  case 22:
    img->LevelIndex=7;
    break;
  case 30:
    img->LevelIndex=8;
    break;
  case 31:
    img->LevelIndex=9;
    break;
  case 32:
    img->LevelIndex=10;
    break;
  case 40:
    img->LevelIndex=11;
    break;
  case 41:
    img->LevelIndex=12;
    break;
  case 42:
    if (!IS_FREXT_PROFILE(active_sps->profile_idc))
      img->LevelIndex=13;
    else
      img->LevelIndex=14;
    break;
  case 50:
    img->LevelIndex=15;
    break;
  case 51:
    img->LevelIndex=16;
    break;
  default:
    fprintf ( stderr, "Warning: unknown LevelIDC, using maximum level 5.1 \n" );
    img->LevelIndex=16;
    break;
  }
}

/*!
 ************************************************************************
 * \brief
 *    initialize key frames and corresponding redundant frames.
 ************************************************************************
 */
void Init_redundant_frame()
{
  if (params->redundant_pic_flag)
  {
    if (params->NumberBFrames)
    {
      error("B frame not supported when redundant picture used!", 100);
    }

    if (params->PicInterlace)
    {
      error("Interlace not supported when redundant picture used!", 100);
    }

    if (params->num_ref_frames < params->PrimaryGOPLength)
    {
      error("NumberReferenceFrames must be no less than PrimaryGOPLength", 100);
    }

    if ((1<<params->NumRedundantHierarchy) > params->PrimaryGOPLength)
    {
      error("PrimaryGOPLength must be greater than 2^NumRedundantHeirarchy", 100);
    }

    if (params->Verbose != 1)
    {
      error("Redundant slices not supported when Verbose != 1", 100);
    }
  }

  key_frame = 0;
  redundant_coding = 0;
  img->redundant_pic_cnt = 0;
  frameNuminGOP = img->frm_number % params->PrimaryGOPLength;
  if (img->frm_number == 0)
  {
    frameNuminGOP = -1;
  }
}

/*!
 ************************************************************************
 * \brief
 *    allocate redundant frames in a primary GOP.
 ************************************************************************
 */
void Set_redundant_frame()
{
  int GOPlength = params->PrimaryGOPLength;

  //start frame of GOP
  if (frameNuminGOP == 0)
  {
    redundant_coding = 0;
    key_frame = 1;
    redundant_ref_idx = GOPlength;
  }

  //1/2 position
  if (params->NumRedundantHierarchy > 0)
  {
    if (frameNuminGOP == GOPlength/2)
    {
      redundant_coding = 0;
      key_frame = 1;
      redundant_ref_idx = GOPlength/2;
    }
  }

  //1/4, 3/4 position
  if (params->NumRedundantHierarchy > 1)
  {
    if (frameNuminGOP == GOPlength/4 || frameNuminGOP == GOPlength*3/4)
    {
      redundant_coding = 0;
      key_frame = 1;
      redundant_ref_idx = GOPlength/4;
    }
  }

  //1/8, 3/8, 5/8, 7/8 position
  if (params->NumRedundantHierarchy > 2)
  {
    if (frameNuminGOP == GOPlength/8 || frameNuminGOP == GOPlength*3/8
      || frameNuminGOP == GOPlength*5/8 || frameNuminGOP == GOPlength*7/8)
    {
      redundant_coding = 0;
      key_frame = 1;
      redundant_ref_idx = GOPlength/8;
    }
  }

  //1/16, 3/16, 5/16, 7/16, 9/16, 11/16, 13/16 position
  if (params->NumRedundantHierarchy > 3)
  {
    if (frameNuminGOP == GOPlength/16 || frameNuminGOP == GOPlength*3/16
      || frameNuminGOP == GOPlength*5/16 || frameNuminGOP == GOPlength*7/16
      || frameNuminGOP == GOPlength*9/16 || frameNuminGOP == GOPlength*11/16
      || frameNuminGOP == GOPlength*13/16)
    {
      redundant_coding = 0;
      key_frame = 1;
      redundant_ref_idx = GOPlength/16;
    }
  }
}

/*!
 ************************************************************************
 * \brief
 *    encode on redundant frame.
 ************************************************************************
 */
void encode_one_redundant_frame()
{
  key_frame = 0;
  redundant_coding = 1;
  img->redundant_pic_cnt = 1;

  if (!img->currentPicture->idr_flag)
  {
    if (img->type == I_SLICE)
    {
      set_slice_type( img, P_SLICE );
    }
  }

  encode_one_frame(img);
}

/*!
 ************************************************************************
 * \brief
 *    Setup Chroma MC Variables
 ************************************************************************
 */
static void chroma_mc_setup(void)
{
  // initialize global variables used for chroma interpolation and buffering
  if ( img->yuv_format == YUV420 )
  {
    img_pad_size_uv_x = IMG_PAD_SIZE >> 1;
    img_pad_size_uv_y = IMG_PAD_SIZE >> 1;
    chroma_mask_mv_y = 7;
    chroma_mask_mv_x = 7;
    chroma_shift_x = 3;
    chroma_shift_y = 3;
  }
  else if ( img->yuv_format == YUV422 )
  {
    img_pad_size_uv_x = IMG_PAD_SIZE >> 1;
    img_pad_size_uv_y = IMG_PAD_SIZE;
    chroma_mask_mv_y = 3;
    chroma_mask_mv_x = 7;
    chroma_shift_y = 2;
    chroma_shift_x = 3;
  }
  else
  { // YUV444
    img_pad_size_uv_x = IMG_PAD_SIZE;
    img_pad_size_uv_y = IMG_PAD_SIZE;
    chroma_mask_mv_y = 3;
    chroma_mask_mv_x = 3;
    chroma_shift_y = 2;
    chroma_shift_x = 2;
  }
  shift_cr_y  = chroma_shift_y - 2;
  shift_cr_x  = chroma_shift_x - 2;
}

