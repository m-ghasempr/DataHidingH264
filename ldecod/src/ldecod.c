
/*!
 ***********************************************************************
 *  \mainpage
 *     This is the H.264/AVC decoder reference software. For detailed documentation
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
 *     ldecod.c
 *  \brief
 *     H.264/AVC reference decoder project main()
 *  \author
 *     Main contributors (see contributors.h for copyright, address and affiliation details)
 *     - Inge Lille-Langøy       <inge.lille-langoy@telenor.com>
 *     - Rickard Sjoberg         <rickard.sjoberg@era.ericsson.se>
 *     - Stephan Wenger          <stewe@cs.tu-berlin.de>
 *     - Jani Lainema            <jani.lainema@nokia.com>
 *     - Sebastian Purreiter     <sebastian.purreiter@mch.siemens.de>
 *     - Byeong-Moon Jeon        <jeonbm@lge.com>
 *     - Gabi Blaettermann
 *     - Ye-Kui Wang             <wyk@ieee.org>
 *     - Valeri George           <george@hhi.de>
 *     - Karsten Suehring        <suehring@hhi.de>
 *
 ***********************************************************************
 */

#include "contributors.h"

#include <sys/stat.h>

#include "global.h"
#include "annexb.h"
#include "image.h"
#include "memalloc.h"
#include "mc_prediction.h"
#include "mbuffer.h"
#include "leaky_bucket.h"
#include "fmo.h"
#include "output.h"
#include "cabac.h"
#include "parset.h"
#include "sei.h"
#include "erc_api.h"
#include "quant.h"
#include "block.h"
#include "nalu.h"
#include "img_io.h"

#define LOGFILE     "log.dec"
#define DATADECFILE "dataDec.txt"
#define TRACEFILE   "trace_dec.txt"

// Decoder definition. This should be the only global variable in the entire
// software. Global variables should be avoided.
DecoderParams  *p_Dec;

// Prototypes of static functions
static void init_conf   (ImageParameters *p_Img, InputParameters *p_Inp, char *config_filename);
static void Report      (ImageParameters *p_Img);
static void init        (ImageParameters *p_Img);
static void malloc_slice(InputParameters *p_Inp, ImageParameters *p_Img);
static void free_slice  (Slice *currSlice);

/*!
 ************************************************************************
 * \brief
 *    Error handling procedure. Print error message to stderr and exit
 *    with supplied code.
 * \param text
 *    Error message
 * \param code
 *    Exit code
 ************************************************************************
 */
void error(char *text, int code)
{
  fprintf(stderr, "%s\n", text);
  flush_dpb(p_Dec->p_Img);
  exit(code);
}


/*!
 ***********************************************************************
 * \brief
 *   print help message and exit
 ***********************************************************************
 */
void JMDecHelpExit (void)
{
  fprintf( stderr, "\n   ldecod [-h] {[defdec.cfg] | {[-p pocScale][-i bitstream.264]...[-o output.yuv] [-r reference.yuv] [-uv]}}\n\n"
    "## Parameters\n\n"

    "## Options\n"
    "   -h  :  prints function usage\n"
    "       :  parse <defdec.cfg> for decoder operation.\n"
    "   -i  :  Input file name. \n"
    "   -o  :  Output file name. If not specified default output is set as test_dec.yuv\n\n"
    "   -r  :  Reference file name. If not specified default output is set as test_rec.yuv\n\n"
    "   -p  :  Poc Scale. \n"
    "   -uv :  write chroma components for monochrome streams(4:2:0)\n"
    "   -lp :  By default the deblocking filter for High Intra-Only profile is off \n\t  regardless of the flags in the bitstream. In the presence of\n\t  this option, the loop filter usage will then be determined \n\t  by the flags and parameters in the bitstream.\n\n" 

    "## Supported video file formats\n"
    "   Input : .264 -> H.264 bitstream files. \n"
    "   Output: .yuv -> RAW file. Format depends on bitstream information. \n\n"

    "## Examples of usage:\n"
    "   ldecod\n"
    "   ldecod  -h\n"
    "   ldecod  default.cfg\n"
    "   ldecod  -i bitstream.264 -o output.yuv -r reference.yuv\n");

  exit(-1);
}


void Configure(ImageParameters *p_Img, InputParameters *p_Inp, int ac, char *av[])
{
  int CLcount = 1;
  char *config_filename=NULL;
  
  p_Img->p_Inp = p_Inp;

  strcpy(p_Inp->infile,"test.264");      //! set default bitstream name
  strcpy(p_Inp->outfile,"test_dec.yuv"); //! set default output file name
  strcpy(p_Inp->reffile,"test_rec.yuv"); //! set default reference file name

  p_Inp->FileFormat = PAR_OF_ANNEXB;
  p_Inp->ref_offset=0;
  p_Inp->poc_scale=2;
  p_Inp->silent = FALSE;
  p_Inp->intra_profile_deblocking = 0;

#ifdef _LEAKYBUCKET_
  p_Inp->R_decoder=500000;          //! Decoder rate
  p_Inp->B_decoder=104000;          //! Decoder buffer size
  p_Inp->F_decoder=73000;           //! Decoder initial delay
  strcpy(p_Inp->LeakyBucketParamFile,"leakybucketparam.cfg");    // file where Leaky Bucket parameters (computed by encoder) are stored
#endif

  if (ac==2)
  {
    if (0 == strncmp (av[1], "-h", 2))
    {
      JMDecHelpExit();
    }
    else if (0 == strncmp (av[1], "-s", 2))
    {
      p_Inp->silent = TRUE;
    }
    else
    {
      config_filename=av[1];
      init_conf(p_Img, p_Inp, av[1]);
    }
    CLcount=2;
  }

  if (ac>=3)
  {
    if (0 == strncmp (av[1], "-i", 2))
    {
      strcpy(p_Inp->infile,av[2]);
      CLcount = 3;
    }
    if (0 == strncmp (av[1], "-h", 2))
    {
      JMDecHelpExit();
    }
    if (0 == strncmp (av[1], "-s", 2))
    {
      p_Inp->silent = TRUE;
    }
  }

  // Parse the command line

  while (CLcount < ac)
  {
    if (0 == strncmp (av[CLcount], "-h", 2))
    {
      JMDecHelpExit();
    }
    else if (0 == strncmp (av[CLcount], "-s", 2))
    {
      p_Inp->silent = TRUE;
      ++CLcount;
    }
    else if (0 == strncmp (av[CLcount], "-i", 2))  //! Input file
    {
      strcpy(p_Inp->infile,av[CLcount+1]);
      CLcount += 2;
    }
    else if (0 == strncmp (av[CLcount], "-o", 2))  //! Output File
    {
      strcpy(p_Inp->outfile,av[CLcount+1]);
      CLcount += 2;
    }
    else if (0 == strncmp (av[CLcount], "-r", 2))  //! Reference File
    {
      strcpy(p_Inp->reffile,av[CLcount+1]);
      CLcount += 2;
    }
    else if (0 == strncmp (av[CLcount], "-p", 2))  //! Poc Scale
    {
      sscanf (av[CLcount+1], "%d", &p_Inp->poc_scale);
      CLcount += 2;
    }
    else if (0 == strncmp (av[CLcount], "-uv", 3))  //! indicate UV writing for 4:0:0
    {
      p_Inp->write_uv = 1;
      ++CLcount;
    }
    else if (0 == strncmp (av[CLcount], "-lp", 3))  
    {
      p_Inp->intra_profile_deblocking = 1;
      ++CLcount;
    }
    else
    {
      snprintf(errortext, ET_SIZE, "Invalid syntax. Use ldecod -h for proper usage");
      error(errortext, 300);
    }
  }

#if TRACE
  if ((p_Dec->p_trace = fopen(TRACEFILE,"w"))==0)             // append new statistic at the end
  {
    snprintf(errortext, ET_SIZE, "Error open file %s!",TRACEFILE);
    error(errortext,500);
  }
#endif

  if ((p_Img->p_out = open(p_Inp->outfile, OPENFLAGS_WRITE, OPEN_PERMISSIONS))==-1)
  {
    snprintf(errortext, ET_SIZE, "Error open file %s ",p_Inp->outfile);
    error(errortext,500);
  }

  fprintf(stdout,"----------------------------- JM %s %s -----------------------------\n", VERSION, EXT_VERSION);
  fprintf(stdout," Decoder config file                    : %s \n",config_filename);
  fprintf(stdout,"--------------------------------------------------------------------------\n");
  fprintf(stdout," Input H.264 bitstream                  : %s \n",p_Inp->infile);
  fprintf(stdout," Output decoded YUV                     : %s \n",p_Inp->outfile);
  fprintf(stdout," Output status file                     : %s \n",LOGFILE);


  if ((p_Img->p_ref = open(p_Inp->reffile,OPENFLAGS_READ))==-1)
  {
    fprintf(stdout," Input reference file                   : %s does not exist \n",p_Inp->reffile);
    fprintf(stdout,"                                          SNR values are not available\n");
  }
  else
    fprintf(stdout," Input reference file                   : %s \n",p_Inp->reffile);

  fprintf(stdout,"--------------------------------------------------------------------------\n");
#ifdef _LEAKYBUCKET_
  fprintf(stdout," Rate_decoder        : %8ld \n",p_Inp->R_decoder);
  fprintf(stdout," B_decoder           : %8ld \n",p_Inp->B_decoder);
  fprintf(stdout," F_decoder           : %8ld \n",p_Inp->F_decoder);
  fprintf(stdout," LeakyBucketParamFile: %s \n",p_Inp->LeakyBucketParamFile); // Leaky Bucket Param file
  calc_buffer(p_Inp);
  fprintf(stdout,"--------------------------------------------------------------------------\n");
#endif
  if (!p_Inp->silent)
  {
    fprintf(stdout,"POC must = frame# or field# for SNRs to be correct\n");
    fprintf(stdout,"--------------------------------------------------------------------------\n");
    fprintf(stdout,"  Frame          POC  Pic#   QP    SnrY     SnrU     SnrV   Y:U:V Time(ms)\n");
    fprintf(stdout,"--------------------------------------------------------------------------\n");
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
  if ((*p_Img   =  (ImageParameters *) calloc(1, sizeof(ImageParameters)))==NULL) 
    no_mem_exit("alloc_img: p_Img");

  if (((*p_Img)->old_slice = (OldSliceParams *) calloc(1, sizeof(OldSliceParams)))==NULL) 
    no_mem_exit("alloc_img: p_Img->old_slice");

  if (((*p_Img)->snr =  (SNRParameters *)calloc(1, sizeof(SNRParameters)))==NULL) 
    no_mem_exit("alloc_img: p_Img->snr");  

  if (((*p_Img)->p_Dpb =  (DecodedPictureBuffer*)calloc(1, sizeof(DecodedPictureBuffer)))==NULL) 
    no_mem_exit("alloc_img: p_Img->p_Dpb");  

  (*p_Img)->p_Dpb->init_done = 0;
  
  (*p_Img)->global_init_done = 0;

#if (ENABLE_OUTPUT_TONEMAPPING)  
  if (((*p_Img)->seiToneMapping =  (ToneMappingSEI*)calloc(1, sizeof(ToneMappingSEI)))==NULL) 
    no_mem_exit("alloc_img: (*p_Img)->seiToneMapping");  
#endif

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
}

  /*!
 ***********************************************************************
 * \brief
 *    Allocate the Decoder Structure
 * \par  Output:
 *    Decoder Parameters
 ***********************************************************************
 */
static void alloc_decoder( DecoderParams **p_Dec)
{
  if ((*p_Dec = (DecoderParams *) calloc(1, sizeof(DecoderParams)))==NULL) 
    no_mem_exit("alloc_decoder: p_Dec");

  alloc_img(&((*p_Dec)->p_Img));
  alloc_params(&((*p_Dec)->p_Inp));
  (*p_Dec)->p_trace = NULL;
  (*p_Dec)->bufferSize = 0;
  (*p_Dec)->bitcounter = 0;
}

/*!
 ***********************************************************************
 * \brief
 *    Free the Image structure
 * \par  Input:
 *    Image Parameters ImageParameters *p_Img
 ***********************************************************************
 */
static void free_img( ImageParameters *p_Img)
{
  //free_mem3Dint(p_Img->fcf    ); 
  if (p_Img != NULL)
  {
    free_annex_b (p_Img);
#if (ENABLE_OUTPUT_TONEMAPPING)  
    if (p_Img->seiToneMapping != NULL)
    {
      free (p_Img->seiToneMapping);
      p_Img->seiToneMapping = NULL;
    }
#endif

    if (p_Img->bitsfile != NULL)
    {
      free (p_Img->bitsfile);
      p_Img->bitsfile = NULL;
    }

    if (p_Img->p_Dpb != NULL)
    {
      free (p_Img->p_Dpb);
      p_Img->p_Dpb = NULL;
    }
    if (p_Img->snr != NULL)
    {
      free (p_Img->snr);
      p_Img->snr = NULL;
    }
    if (p_Img->old_slice != NULL)
    {
      free (p_Img->old_slice);
      p_Img->old_slice = NULL;
    }

    free (p_Img);
    p_Img = NULL;
  }
}
/*!
 ***********************************************************************
 * \brief
 *    main function for TML decoder
 ***********************************************************************
 */
int main(int argc, char **argv)
{  
  alloc_decoder(&p_Dec);

  Configure (p_Dec->p_Img, p_Dec->p_Inp, argc, argv);

  initBitsFile(p_Dec->p_Img, p_Dec->p_Inp->FileFormat);

  p_Dec->p_Img->bitsfile->OpenBitsFile(p_Dec->p_Img, p_Dec->p_Inp->infile);
  
  // Allocate Slice data struct
  malloc_slice(p_Dec->p_Inp, p_Dec->p_Img);
  init_old_slice(p_Dec->p_Img->old_slice);

  init(p_Dec->p_Img);
 
  init_out_buffer(p_Dec->p_Img);  

  while (decode_one_frame(p_Dec->p_Img) != EOS)
    ;

  Report(p_Dec->p_Img);
  free_slice(p_Dec->p_Img->currentSlice);
  FmoFinit(p_Dec->p_Img);

  free_global_buffers(p_Dec->p_Img);
  flush_dpb(p_Dec->p_Img);

#if (PAIR_FIELDS_IN_OUTPUT)
  flush_pending_output(p_Dec->p_Img, p_Dec->p_Img->p_out);
#endif

  p_Dec->p_Img->bitsfile->CloseBitsFile(p_Dec->p_Img);

  close(p_Dec->p_Img->p_out);

  if (p_Dec->p_Img->p_ref != -1)
    close(p_Dec->p_Img->p_ref);

#if TRACE
  fclose(p_Dec->p_trace);
#endif

  ercClose(p_Dec->p_Img, p_Dec->p_Img->erc_errorVar);

  CleanUpPPS(p_Dec->p_Img);
  free_dpb(p_Dec->p_Img);
  uninit_out_buffer(p_Dec->p_Img);

  free (p_Dec->p_Inp);
  free_img (p_Dec->p_Img);
  free(p_Dec);

  return 0;
}


/*!
 ***********************************************************************
 * \brief
 *    Initilize some arrays
 ***********************************************************************
 */
static void init(ImageParameters *p_Img)  //!< image parameters
{
  int i;
  InputParameters *p_Inp = p_Img->p_Inp;
  p_Img->oldFrameSizeInMbs = -1;

  p_Img->imgY_ref  = NULL;
  p_Img->imgUV_ref = NULL;

  p_Img->recovery_point = 0;
  p_Img->recovery_point_found = 0;
  p_Img->recovery_poc = 0x7fffffff; /* set to a max value */

  p_Img->idr_psnr_number = p_Inp->ref_offset;
  p_Img->psnr_number=0;

  p_Img->number = 0;
  p_Img->type = I_SLICE;

  p_Img->dec_ref_pic_marking_buffer = NULL;

  p_Img->g_nFrame = 0;
  // B pictures
  p_Img->Bframe_ctr = p_Img->snr->frame_ctr = 0;

  // time for total decoding session
  p_Img->tot_time = 0;

  p_Img->dec_picture = NULL;
  // reference flag initialization
  for(i=0;i<17;++i)
  {
    p_Img->ref_flag[i] = 1;
  }

  p_Img->MbToSliceGroupMap = NULL;
  p_Img->MapUnitToSliceGroupMap = NULL;

  p_Img->LastAccessUnitExists  = 0;
  p_Img->NALUCount = 0;


  p_Img->out_buffer = NULL;
  p_Img->pending_output = NULL;
  p_Img->pending_output_state = FRAME;
  p_Img->recovery_flag = 0;


#if (ENABLE_OUTPUT_TONEMAPPING)
  init_tone_mapping_sei(p_Img->seiToneMapping);
#endif

}

/*!
 ***********************************************************************
 * \brief
 *    Initialize FREXT variables
 ***********************************************************************
 */
void init_frext(ImageParameters *p_Img)  //!< image parameters
{
  //pel bitdepth init
  p_Img->bitdepth_luma_qp_scale   = 6 * (p_Img->bitdepth_luma - 8);

  if(p_Img->bitdepth_luma > p_Img->bitdepth_chroma || p_Img->active_sps->chroma_format_idc == YUV400)
    p_Img->pic_unit_bitsize_on_disk = (p_Img->bitdepth_luma > 8)? 16:8;
  else
    p_Img->pic_unit_bitsize_on_disk = (p_Img->bitdepth_chroma > 8)? 16:8;
  p_Img->dc_pred_value_comp[0]    = 1<<(p_Img->bitdepth_luma - 1);
  p_Img->max_imgpel_value_comp[0] = (1<<p_Img->bitdepth_luma) - 1;
  p_Img->mb_size[0][0] = p_Img->mb_size[0][1] = MB_BLOCK_SIZE;

  if (p_Img->active_sps->chroma_format_idc != YUV400)
  {
    //for chrominance part
    p_Img->bitdepth_chroma_qp_scale = 6 * (p_Img->bitdepth_chroma - 8);
    p_Img->dc_pred_value_comp[1]    = (1 << (p_Img->bitdepth_chroma - 1));
    p_Img->dc_pred_value_comp[2]    = p_Img->dc_pred_value_comp[1];
    p_Img->max_imgpel_value_comp[1] = (1 << p_Img->bitdepth_chroma) - 1;
    p_Img->max_imgpel_value_comp[2] = (1 << p_Img->bitdepth_chroma) - 1;
    p_Img->num_blk8x8_uv = (1 << p_Img->active_sps->chroma_format_idc) & (~(0x1));
    p_Img->num_uv_blocks = (p_Img->num_blk8x8_uv >> 1);
    p_Img->num_cdc_coeff = (p_Img->num_blk8x8_uv << 1);
    p_Img->mb_size[1][0] = p_Img->mb_size[2][0] = p_Img->mb_cr_size_x  = (p_Img->active_sps->chroma_format_idc==YUV420 || p_Img->active_sps->chroma_format_idc==YUV422)?  8 : 16;
    p_Img->mb_size[1][1] = p_Img->mb_size[2][1] = p_Img->mb_cr_size_y  = (p_Img->active_sps->chroma_format_idc==YUV444 || p_Img->active_sps->chroma_format_idc==YUV422)? 16 :  8;
  }
  else
  {
    p_Img->bitdepth_chroma_qp_scale = 0;
    p_Img->max_imgpel_value_comp[1] = 0;
    p_Img->max_imgpel_value_comp[2] = 0;
    p_Img->num_blk8x8_uv = 0;
    p_Img->num_uv_blocks = 0;
    p_Img->num_cdc_coeff = 0;
    p_Img->mb_size[1][0] = p_Img->mb_size[2][0] = p_Img->mb_cr_size_x  = 0;
    p_Img->mb_size[1][1] = p_Img->mb_size[2][1] = p_Img->mb_cr_size_y  = 0;
  }
  p_Img->mb_size_blk[0][0] = p_Img->mb_size_blk[0][1] = p_Img->mb_size[0][0] >> 2;
  p_Img->mb_size_blk[1][0] = p_Img->mb_size_blk[2][0] = p_Img->mb_size[1][0] >> 2;
  p_Img->mb_size_blk[1][1] = p_Img->mb_size_blk[2][1] = p_Img->mb_size[1][1] >> 2;

  p_Img->mb_size_shift[0][0] = p_Img->mb_size_shift[0][1] = CeilLog2_sf (p_Img->mb_size[0][0]);
  p_Img->mb_size_shift[1][0] = p_Img->mb_size_shift[2][0] = CeilLog2_sf (p_Img->mb_size[1][0]);
  p_Img->mb_size_shift[1][1] = p_Img->mb_size_shift[2][1] = CeilLog2_sf (p_Img->mb_size[1][1]);
}

/*!
************************************************************************
* \brief
*    exit with error message if reading from config file failed
************************************************************************
*/
static inline void conf_read_check (int val, int expected)
{
  if (val != expected)
  {
    error ("init_conf: error reading from config file", 500);
  }
}

/*!
 ************************************************************************
 * \brief
 *    Read parameters from configuration file
 *
 * \par Input:
 *    Name of configuration filename
 *
 * \par Output
 *    none
 ************************************************************************
 */
static void init_conf(ImageParameters *p_Img, InputParameters *p_Inp, char *config_filename)
{
  FILE *fd;
  int NAL_mode;

  // picture error concealment
  long int temp;
  char tempval[100];

  // read the decoder configuration file
  if((fd=fopen(config_filename,"r")) == NULL)
  {
    snprintf(errortext, ET_SIZE, "Error: Control file %s not found\n",config_filename);
    error(errortext, 300);
  }

  conf_read_check (fscanf(fd,"%s",p_Inp->infile), 1);                // H.264 compressed input bitstream
  conf_read_check (fscanf(fd,"%*[^\n]"), 0);

  conf_read_check (fscanf(fd,"%s",p_Inp->outfile), 1);               // RAW (YUV/RGB) output file
  conf_read_check (fscanf(fd,"%*[^\n]"), 0);

  conf_read_check (fscanf(fd,"%s",p_Inp->reffile), 1);               // reference file
  conf_read_check (fscanf(fd,"%*[^\n]"), 0);

  conf_read_check (fscanf(fd,"%d",&(p_Inp->write_uv)), 1);           // write UV in YUV 4:0:0 mode
  conf_read_check (fscanf(fd,"%*[^\n]"), 0);

  conf_read_check (fscanf(fd,"%d",&(NAL_mode)), 1);                // NAL mode
  conf_read_check (fscanf(fd,"%*[^\n]"), 0);

  switch(NAL_mode)
  {
  case 0:
    p_Inp->FileFormat = PAR_OF_ANNEXB;
    break;
  case 1:
    p_Inp->FileFormat = PAR_OF_RTP;
    break;
  default:
    snprintf(errortext, ET_SIZE, "NAL mode %i is not supported", NAL_mode);
    error(errortext,400);
  }

  conf_read_check (fscanf(fd,"%d,",&p_Inp->ref_offset), 1);   // offset used for SNR computation
  conf_read_check (fscanf(fd,"%*[^\n]"), 0);

  conf_read_check (fscanf(fd,"%d,",&p_Inp->poc_scale), 1);   // offset used for SNR computation
  conf_read_check (fscanf(fd,"%*[^\n]"), 0);


  if (p_Inp->poc_scale < 1 || p_Inp->poc_scale > 10)
  {
    snprintf(errortext, ET_SIZE, "Poc Scale is %d. It has to be within range 1 to 10",p_Inp->poc_scale);
    error(errortext,1);
  }

  p_Inp->write_uv=1;

  // picture error concealment
  p_Img->conceal_mode = p_Inp->conceal_mode = 0;
  p_Img->ref_poc_gap = p_Inp->ref_poc_gap = 2;
  p_Img->poc_gap = p_Inp->poc_gap = 2;

#ifdef _LEAKYBUCKET_
  conf_read_check (fscanf(fd,"%ld,",&p_Inp->R_decoder), 1);             // Decoder rate
  conf_read_check (fscanf(fd, "%*[^\n]"), 0);
  conf_read_check (fscanf(fd,"%ld,",&p_Inp->B_decoder), 1);             // Decoder buffer size
  conf_read_check (fscanf(fd, "%*[^\n]"), 0);
  conf_read_check (fscanf(fd,"%ld,",&p_Inp->F_decoder), 1);             // Decoder initial delay
  conf_read_check (fscanf(fd, "%*[^\n]"), 0);
  conf_read_check (fscanf(fd,"%s",p_Inp->LeakyBucketParamFile), 1);    // file where Leaky Bucket params (computed by encoder) are stored
  conf_read_check (fscanf(fd,"%*[^\n]"), 0);
#endif

  /* since error concealment parameters are added at the end of
  decoder conf file we need to read the leakybucket params to get to
  those parameters */
#ifndef _LEAKYBUCKET_
  conf_read_check (fscanf(fd,"%ld,",&temp), 1);
  conf_read_check (fscanf(fd, "%*[^\n]"), 0);
  conf_read_check (fscanf(fd,"%ld,",&temp), 1);
  conf_read_check (fscanf(fd, "%*[^\n]"), 0);
  conf_read_check (fscanf(fd,"%ld,",&temp), 1);
  conf_read_check (fscanf(fd, "%*[^\n]"), 0);
  conf_read_check (fscanf(fd,"%s",tempval), 1);
  conf_read_check (fscanf(fd,"%*[^\n]"), 0);
#endif

  conf_read_check (fscanf(fd,"%d",&p_Inp->conceal_mode), 1);   // Mode of Error Concealment
  conf_read_check (fscanf(fd,"%*[^\n]"), 0);
  p_Img->conceal_mode = p_Inp->conceal_mode;
  conf_read_check (fscanf(fd,"%d",&p_Inp->ref_poc_gap), 1);   // POC gap depending on pattern
  conf_read_check (fscanf(fd,"%*[^\n]"), 0);
  p_Img->ref_poc_gap = p_Inp->ref_poc_gap;
  conf_read_check (fscanf(fd,"%d",&p_Inp->poc_gap), 1);   // POC gap between consecutive frames in display order
  conf_read_check (fscanf(fd,"%*[^\n]"), 0);
  p_Img->poc_gap = p_Inp->poc_gap;
  conf_read_check (fscanf(fd,"%d,",&p_Inp->silent), 1);     // use silent decode mode
  conf_read_check (fscanf(fd,"%*[^\n]"), 0);
  conf_read_check (fscanf(fd,"%d,",&p_Inp->intra_profile_deblocking), 1);     // use deblocking filter in intra only profile
  conf_read_check (fscanf(fd,"%*[^\n]"), 0);

  fclose (fd);
}

/*!
 ************************************************************************
 * \brief
 *    Reports the gathered information to appropriate outputs
 *
 * \par Input:
 *    InputParameters *p_Inp,
 *    ImageParameters *p_Img,
 *    struct snr_par *stat
 *
 * \par Output:
 *    None
 ************************************************************************
 */
static void Report(ImageParameters *p_Img)
{
  pic_parameter_set_rbsp_t *active_pps = p_Img->active_pps;
  InputParameters *p_Inp = p_Img->p_Inp;
  SNRParameters   *snr   = p_Img->snr;
#define OUTSTRING_SIZE 255
  char string[OUTSTRING_SIZE];
  FILE *p_log;
  static const char yuv_formats[4][4]= { {"400"}, {"420"}, {"422"}, {"444"} };

#ifndef WIN32
  time_t  now;
  struct tm *l_time;
#else
  char timebuf[128];
#endif

  // normalize time
  p_Img->tot_time  = timenorm(p_Img->tot_time);

  if (p_Inp->silent == FALSE)
  {
    fprintf(stdout,"-------------------- Average SNR all frames ------------------------------\n");
    fprintf(stdout," SNR Y(dB)           : %5.2f\n",snr->snra[0]);
    fprintf(stdout," SNR U(dB)           : %5.2f\n",snr->snra[1]);
    fprintf(stdout," SNR V(dB)           : %5.2f\n",snr->snra[2]);
    fprintf(stdout," Total decoding time : %.3f sec (%.3f fps)\n",p_Img->tot_time*0.001,(snr->frame_ctr ) * 1000.0 / p_Img->tot_time);
    fprintf(stdout,"--------------------------------------------------------------------------\n");
    fprintf(stdout," Exit JM %s decoder, ver %s ",JM, VERSION);
    fprintf(stdout,"\n");
  }
  else
  {
    fprintf(stdout,"\n----------------------- Decoding Completed -------------------------------\n");
    fprintf(stdout," Total decoding time : %.3f sec (%.3f fps)\n",p_Img->tot_time*0.001, (snr->frame_ctr) * 1000.0 / p_Img->tot_time);
    fprintf(stdout,"--------------------------------------------------------------------------\n");
    fprintf(stdout," Exit JM %s decoder, ver %s ",JM, VERSION);
    fprintf(stdout,"\n");
  }

  // write to log file

  snprintf(string, OUTSTRING_SIZE, "%s", LOGFILE);

  if ((p_log=fopen(string,"r"))==0)                    // check if file exist
  {
    if ((p_log=fopen(string,"a"))==0)
    {
      snprintf(errortext, ET_SIZE, "Error open file %s for appending",string);
      error(errortext, 500);
    }
    else                                              // Create header to new file
    {
      fprintf(p_log," -------------------------------------------------------------------------------------------------------------------\n");
      fprintf(p_log,"|  Decoder statistics. This file is made first time, later runs are appended               |\n");
      fprintf(p_log," ------------------------------------------------------------------------------------------------------------------- \n");
      fprintf(p_log,"|   ver  | Date  | Time  |    Sequence        |#Img| Format  | YUV |Coding|SNRY 1|SNRU 1|SNRV 1|SNRY N|SNRU N|SNRV N|\n");
      fprintf(p_log," -------------------------------------------------------------------------------------------------------------------\n");
    }
  }
  else
  {
    fclose(p_log);
    p_log=fopen(string,"a");                    // File exist,just open for appending
  }

  fprintf(p_log,"|%s/%-4s", VERSION, EXT_VERSION);

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
  fprintf(p_log,"| %1.5s |",string );
#endif

  fprintf(p_log,"%20.20s|",p_Inp->infile);

  fprintf(p_log,"%3d |",p_Img->number);
  fprintf(p_log,"%4dx%-4d|", p_Img->width, p_Img->height);
  fprintf(p_log," %s |", &(yuv_formats[p_Img->yuv_format][0]));

  if (active_pps)
  {
    if (active_pps->entropy_coding_mode_flag == CAVLC)
      fprintf(p_log," CAVLC|");
    else
      fprintf(p_log," CABAC|");
  }

  fprintf(p_log,"%6.3f|",snr->snr1[0]);
  fprintf(p_log,"%6.3f|",snr->snr1[1]);
  fprintf(p_log,"%6.3f|",snr->snr1[2]);
  fprintf(p_log,"%6.3f|",snr->snra[0]);
  fprintf(p_log,"%6.3f|",snr->snra[1]);
  fprintf(p_log,"%6.3f|",snr->snra[2]);
  fprintf(p_log,"\n");
  fclose(p_log);

  snprintf(string, OUTSTRING_SIZE,"%s", DATADECFILE);
  p_log=fopen(string,"a");

  if(p_Img->Bframe_ctr != 0) // B picture used
  {
    fprintf(p_log, "%3d %2d %2d %2.2f %2.2f %2.2f %5d "
      "%2.2f %2.2f %2.2f %5d "
      "%2.2f %2.2f %2.2f %5d %.3f\n",
      p_Img->number, 0, p_Img->qp,
      snr->snr1[0],
      snr->snr1[1],
      snr->snr1[2],
      0,
      0.0,
      0.0,
      0.0,
      0,
      snr->snra[0],
      snr->snra[1],
      snr->snra[2],
      0,
      (double)0.001*p_Img->tot_time/(p_Img->number + p_Img->Bframe_ctr - 1));
  }
  else
  {
    fprintf(p_log, "%3d %2d %2d %2.2f %2.2f %2.2f %5d "
      "%2.2f %2.2f %2.2f %5d "
      "%2.2f %2.2f %2.2f %5d %.3f\n",
      p_Img->number, 0, p_Img->qp,
      snr->snr1[0],
      snr->snr1[1],
      snr->snr1[2],
      0,
      0.0,
      0.0,
      0.0,
      0,
      snr->snra[0],
      snr->snra[1],
      snr->snra[2],
      0,
      (double)0.001*p_Img->tot_time/p_Img->number);
  }
  fclose(p_log);
}

/*!
 ************************************************************************
 * \brief
 *    Allocates a stand-alone partition structure.  Structure should
 *    be freed by FreePartition();
 *    data structures
 *
 * \par Input:
 *    n: number of partitions in the array
 * \par return
 *    pointer to DataPartition Structure, zero-initialized
 ************************************************************************
 */

DataPartition *AllocPartition(int n)
{
  DataPartition *partArr, *dataPart;
  int i;

  partArr = (DataPartition *) calloc(n, sizeof(DataPartition));
  if (partArr == NULL)
  {
    snprintf(errortext, ET_SIZE, "AllocPartition: Memory allocation for Data Partition failed");
    error(errortext, 100);
  }

  for (i=0; i<n; ++i) // loop over all data partitions
  {
    dataPart = &(partArr[i]);
    dataPart->bitstream = (Bitstream *) calloc(1, sizeof(Bitstream));
    if (dataPart->bitstream == NULL)
    {
      snprintf(errortext, ET_SIZE, "AllocPartition: Memory allocation for Bitstream failed");
      error(errortext, 100);
    }
    dataPart->bitstream->streamBuffer = (byte *) calloc(MAX_CODED_FRAME_SIZE, sizeof(byte));
    if (dataPart->bitstream->streamBuffer == NULL)
    {
      snprintf(errortext, ET_SIZE, "AllocPartition: Memory allocation for streamBuffer failed");
      error(errortext, 100);
    }
  }
  return partArr;
}




/*!
 ************************************************************************
 * \brief
 *    Frees a partition structure (array).
 *
 * \par Input:
 *    Partition to be freed, size of partition Array (Number of Partitions)
 *
 * \par return
 *    None
 *
 * \note
 *    n must be the same as for the corresponding call of AllocPartition
 ************************************************************************
 */


void FreePartition (DataPartition *dp, int n)
{
  int i;

  assert (dp != NULL);
  assert (dp->bitstream != NULL);
  assert (dp->bitstream->streamBuffer != NULL);
  for (i=0; i<n; ++i)
  {
    free (dp[i].bitstream->streamBuffer);
    free (dp[i].bitstream);
  }
  free (dp);
}


/*!
 ************************************************************************
 * \brief
 *    Allocates the slice structure along with its dependent
 *    data structures
 *
 * \par Input:
 *    Input Parameters InputParameters *p_Inp,  ImageParameters *p_Img
 ************************************************************************
 */
static void malloc_slice(InputParameters *p_Inp, ImageParameters *p_Img)
{
  int memory_size = 0;
  Slice *currSlice;

  p_Img->currentSlice = (Slice *) calloc(1, sizeof(Slice));
  if ( (currSlice = p_Img->currentSlice) == NULL)
  {
    snprintf(errortext, ET_SIZE, "Memory allocation for Slice datastruct in NAL-mode %d failed", p_Inp->FileFormat);
    error(errortext,100);
  }
  //  p_Img->currentSlice->rmpni_buffer=NULL;
  //! you don't know whether we do CABAC hre, hence initialize CABAC anyway
  // if (p_Inp->symbol_mode == CABAC)

  // create all context models
  currSlice->mot_ctx = create_contexts_MotionInfo();
  currSlice->tex_ctx = create_contexts_TextureInfo();


  currSlice->max_part_nr = 3;  //! assume data partitioning (worst case) for the following mallocs()
  currSlice->partArr = AllocPartition(currSlice->max_part_nr);
  currSlice->p_colocated = NULL;

  memory_size += get_mem3Dint(&(currSlice->wp_weight), 2, MAX_REFERENCE_PICTURES, 3);
  memory_size += get_mem3Dint(&(currSlice->wp_offset), 6, MAX_REFERENCE_PICTURES, 3);
  memory_size += get_mem4Dint(&(currSlice->wbp_weight), 6, MAX_REFERENCE_PICTURES, MAX_REFERENCE_PICTURES, 3);

  memory_size += get_mem3Dpel(&(currSlice->mb_pred), MAX_PLANE, MB_BLOCK_SIZE, MB_BLOCK_SIZE);
  memory_size += get_mem3Dpel(&(currSlice->mb_rec ), MAX_PLANE, MB_BLOCK_SIZE, MB_BLOCK_SIZE);
  memory_size += get_mem3Dint(&(currSlice->mb_rres), MAX_PLANE, MB_BLOCK_SIZE, MB_BLOCK_SIZE);
  memory_size += get_mem3Dint(&(currSlice->cof    ), MAX_PLANE, MB_BLOCK_SIZE, MB_BLOCK_SIZE);
  //  memory_size += get_mem3Dint(&(currSlice->fcf    ), MAX_PLANE, MB_BLOCK_SIZE, MB_BLOCK_SIZE);

  allocate_pred_mem(currSlice);
}


/*!
 ************************************************************************
 * \brief
 *    Memory frees of the Slice structure and of its dependent
 *    data structures
 *
 * \par Input:
 *    Input Parameters InputParameters *p_Inp,  ImageParameters *p_Img
 ************************************************************************
 */
static void free_slice(Slice *currSlice)
{
  free_pred_mem(currSlice);

  free_mem3Dint(currSlice->cof    );
  free_mem3Dint(currSlice->mb_rres);
  free_mem3Dpel(currSlice->mb_rec );
  free_mem3Dpel(currSlice->mb_pred);


  free_mem3Dint(currSlice->wp_weight );
  free_mem3Dint(currSlice->wp_offset );
  free_mem4Dint(currSlice->wbp_weight);

  FreePartition (currSlice->partArr, 3);
  // if (p_Inp->symbol_mode == CABAC)
  if (1)
  {
    // delete all context models
    delete_contexts_MotionInfo(currSlice->mot_ctx);
    delete_contexts_TextureInfo(currSlice->tex_ctx);
  }
  free(currSlice);

  currSlice = NULL;
}

/*!
 ************************************************************************
 * \brief
 *    Dynamic memory allocation of frame size related global buffers
 *    buffers are defined in global.h, allocated memory must be freed in
 *    void free_global_buffers()
 *
 *  \par Input:
 *    Input Parameters InputParameters *p_Inp, Image Parameters ImageParameters *p_Img
 *
 *  \par Output:
 *     Number of allocated bytes
 ***********************************************************************
 */
int init_global_buffers(ImageParameters *p_Img)
{
  int memory_size=0;
  int i;

  if (p_Img->global_init_done)
  {
    free_global_buffers(p_Img);
  }

  // allocate memory for reference frame in find_snr
  memory_size += get_mem2Dpel(&p_Img->imgY_ref, p_Img->height, p_Img->width);

  if (p_Img->active_sps->chroma_format_idc != YUV400)
    memory_size += get_mem3Dpel(&p_Img->imgUV_ref, 2, p_Img->height_cr, p_Img->width_cr);
  else
    p_Img->imgUV_ref=NULL;

  // allocate memory in structure p_Img
  if( IS_INDEPENDENT(p_Img) )
  {
    for( i=0; i<MAX_PLANE; ++i )
    {
      if(((p_Img->mb_data_JV[i]) = (Macroblock *) calloc(p_Img->FrameSizeInMbs, sizeof(Macroblock))) == NULL)
        no_mem_exit("init_global_buffers: p_Img->mb_data");
    }
    p_Img->mb_data = NULL;
  }
  else
  {
    if(((p_Img->mb_data) = (Macroblock *) calloc(p_Img->FrameSizeInMbs, sizeof(Macroblock))) == NULL)
      no_mem_exit("init_global_buffers: p_Img->mb_data");
  }

  if(((p_Img->intra_block) = (int*)calloc(p_Img->FrameSizeInMbs, sizeof(int))) == NULL)
    no_mem_exit("init_global_buffers: p_Img->intra_block");

  memory_size += get_mem2Dint(&PicPos,p_Img->FrameSizeInMbs + 1,2);  //! Helper array to access macroblock positions. We add 1 to also consider last MB.

  for (i = 0; i < (int) p_Img->FrameSizeInMbs + 1;++i)
  {
    PicPos[i][0] = (i % p_Img->PicWidthInMbs);
    PicPos[i][1] = (i / p_Img->PicWidthInMbs);
  }

  memory_size += get_mem2D(&(p_Img->ipredmode), 4*p_Img->FrameHeightInMbs, 4*p_Img->PicWidthInMbs);

  // CAVLC mem
  memory_size += get_mem4D(&(p_Img->nz_coeff), p_Img->FrameSizeInMbs, 3, BLOCK_SIZE, BLOCK_SIZE);

  memory_size += get_mem2Dint(&(p_Img->siblock), p_Img->FrameHeightInMbs, p_Img->PicWidthInMbs);

  init_qp_process(p_Img);

  p_Img->global_init_done = 1;

  p_Img->oldFrameSizeInMbs = p_Img->FrameSizeInMbs;

  return (memory_size);
}

/*!
 ************************************************************************
 * \brief
 *    Free allocated memory of frame size related global buffers
 *    buffers are defined in global.h, allocated memory is allocated in
 *    int init_global_buffers()
 *
 * \par Input:
 *    Input Parameters InputParameters *p_Inp, Image Parameters ImageParameters *p_Img
 *
 * \par Output:
 *    none
 *
 ************************************************************************
 */
void free_global_buffers(ImageParameters *p_Img)
{  
  free_mem2Dpel (p_Img->imgY_ref);

  if (p_Img->imgUV_ref)
    free_mem3Dpel (p_Img->imgUV_ref);

  // CAVLC free mem
  free_mem4D(p_Img->nz_coeff);

  free_mem2Dint(p_Img->siblock);

  // free mem, allocated for structure p_Img
  if (p_Img->mb_data != NULL)
    free(p_Img->mb_data);

  free_mem2Dint(PicPos);

  free (p_Img->intra_block);
  free_mem2D(p_Img->ipredmode);

  free_qp_matrices(p_Img);

  p_Img->global_init_done = 0;

}

void report_stats_on_error(void)
{
  //free_encoder_memory(p_Img);
  exit (-1);
}
