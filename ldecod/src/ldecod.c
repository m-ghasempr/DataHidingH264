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
 *     This is the H.26L decoder reference software. For detailed documentation
 *     see the comments in each file.
 *
 *  \author
 *     The main contributors are listed in contributors.h
 *
 *  \version
 *     JM 5.0c
 *
 *  \note
 *     tags are used for document system "doxygen"
 *     available at http://www.doxygen.org
 *
 *  \par Limitations:
 *     Using different NAL's the assignment of partition-id to containing
 *     syntax elements may got lost, if this information is not transmitted.
 *     The same has to be stated for the partitionlength if partitions are
 *     merged by the NAL.
 *  \par
 *     The presented solution in Q15-K-16 solves both of this problems as the
 *     departitioner parses the bitstream before decoding. Due to syntax element
 *     dependencies both, partition bounds and partitionlength information can
 *     be parsed by the departitioner.
 *
 *  \par Handling partition information in external file:
 *     As the TML is still a work in progress, it makes sense to handle this
 *     information for simplification in an external file, here called partition
 *     information file, which can be found by the extension .dp extending the
 *     original encoded H.26L bitstream. In this file partition-ids followed by its
 *     partitionlength is written. Instead of parsing the bitstream we get the
 *     partition information now out of this file.
 *     This information is assumed to be never sent over transmission channels
 *     (simulation scenarios) as it's information we allways get using a
 *     "real" departitioner before decoding
 *
 *  \par Extension of Interim File Format:
 *     Therefore a convention has to be made within the interim file format.
 *     The underlying NAL has to take care of fulfilling these conventions.
 *     All partitions have to be bytealigned to be readable by the decoder,
 *     So if the NAL-encoder merges partitions, >>this is only possible to use the
 *     VLC structure of the H.26L bitstream<<, this bitaligned structure has to be
 *     broken up by the NAL-decoder. In this case the NAL-decoder is responsable to
 *     read the partitionlength information from the partition information file.
 *     Partitionlosses are signaled with a partition of zero length containing no
 *     syntax elements.
 *
 */
/*!
 *  \file
 *     ldecod.c
 *  \brief
 *     TML decoder project main()
 *  \author
 *     Main contributors (see contributors.h for copyright, address and affiliation details)
 *     - Inge Lille-Langøy       <inge.lille-langoy@telenor.com>
 *     - Rickard Sjoberg         <rickard.sjoberg@era.ericsson.se>
 *     - Stephan Wenger          <stewe@cs.tu-berlin.de>
 *     - Jani Lainema            <jani.lainema@nokia.com>
 *     - Sebastian Purreiter     <sebastian.purreiter@mch.siemens.de>
 *     - Byeong-Moon Jeon        <jeonbm@lge.com>
 *     - Gabi Blaettermann       <blaetter@hhi.de>
 *     - Ye-Kui Wang             <wyk@ieee.org>
 *
 ***********************************************************************
 */

#include "contributors.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <sys/timeb.h>

#if defined WIN32
  #include <conio.h>
#endif

#include "global.h"
#include "bitsbuf.h"
#include "rtp.h"
#include "memalloc.h"
#include "mbuffer.h"
#include "leaky_bucket.h"
#include "decodeiff.h"
#include "fmo.h"

#if _ERROR_CONCEALMENT_
#include "erc_api.h"
#endif

#define JM          "5"
#define VERSION     "5.0c"

#define LOGFILE     "log.dec"
#define DATADECFILE "dataDec.txt"
#define TRACEFILE   "trace_dec.txt"

#if _ERROR_CONCEALMENT_
extern objectBuffer_t *erc_object_list;
extern ercVariables_t *erc_errorVar;
#endif


/*!
 ***********************************************************************
 * \brief
 *    main function for TML decoder
 ***********************************************************************
 */
int main(int argc, char **argv)
{
  extern FILE* bits;
  struct inp_par    *inp;         // input parameters from input configuration file
  struct snr_par    *snr;         // statistics
  struct img_par    *img;         // image parameters

    // allocate memory for the structures
  if ((inp =  (struct inp_par *)calloc(1, sizeof(struct inp_par)))==NULL) no_mem_exit("main: inp");
  if ((snr =  (struct snr_par *)calloc(1, sizeof(struct snr_par)))==NULL) no_mem_exit("main: snr");
  if ((img =  (struct img_par *)calloc(1, sizeof(struct img_par)))==NULL) no_mem_exit("main: img");

  // Read Configuration File
  if (argc != 2)
  {
    snprintf(errortext, ET_SIZE, "Usage: %s <config.dat> \n\t<config.dat> defines decoder parameters",argv[0]);
    error(errortext, 300);
  }

  isBigEndian = testEndian();

  // Initializes Configuration Parameters with configuration file
  init_conf(inp, argv[1]);

  img->UseConstrainedIntraPred = inp->UseConstrainedIntraPred;

  if (inp->of_mode == PAR_OF_RTP)
  {
    extern FILE *bits;
    // Read the first RTP packet conmtaining a header packet, and set the initial parameter set
    RTPSequenceHeader (img, inp, bits);
  }
  if (inp->of_mode == PAR_OF_IFF)
  {
    // Read the first boxes, and set the initial parameter set
    if ( -1 == IFFSequenceHeader( img, inp, bits ) )
    {
      snprintf(errortext, ET_SIZE, "Error: The input file is not in the Interim File Format\n");
      error(errortext, 300);
    }
  }
// printf ("In main: some pictrue information: %d x %d, with %d reference frames %d\n", img->height, img->width, img->buf_cycle, inp->buf_cycle);
#ifndef _ABT_FLAG_IN_SLICE_HEADER_
  USEABT = inp->abt; // set global ABT flag. (0=OFF, 1=Inter, 2=Inter&Intra)
#endif

  // Allocate Slice data struct
  malloc_slice(inp,img);

  init(img);
  img->number=0;
  img->type = INTRA_IMG;
  img->tr_old = -1; // WYK: Oct. 8, 2001, for detection of a new frame
  img->refPicID = -1; // WYK: for detection of a new non-B frame
  img->imgtr_last_P = 0;
  img->imgtr_next_P = 0;
  img->mmco_buffer=NULL;
  img->last_decoded_pic_id = -1; // JVT-D101

  // B pictures
  Bframe_ctr=0;

  // time for total decoding session
  tot_time = 0;
  if ( inp->of_mode == PAR_OF_IFF )
    while ( parse_one_box(img, inp, snr, bits) != -1 );
  else
    while (decode_one_frame(img, inp, snr) != EOS);

  // B PICTURE : save the last P picture
  write_prev_Pframe(img, p_out);

  report(inp, img, snr);

  free_slice(inp,img);

  FmoFinit();

  free_frame_buffers(inp, img);
  free_global_buffers(inp, img);

  terminateInterimFile();
  CloseBitstreamFile();

  fclose(p_out);
  if (p_ref)
    fclose(p_ref);
#if TRACE
  fclose(p_trace);
#endif

#if _ERROR_CONCEALMENT_
  ercClose(erc_errorVar);
#endif

  free (inp);
  free (snr);
  free (img);
  
  //while( !kbhit() ); 
  return 0;
}


/*!
 ***********************************************************************
 * \brief
 *    Initilize some arrays
 ***********************************************************************
 */
void init(struct img_par *img)  //!< image parameters
{
  int i;

  // initilize quad matrix used in snr routine
  for (i=0; i <  256; i++)
  {
    img->quad[i]=i*i; // fix from TML 1, truncation removed
  }
}

/*!
 ************************************************************************
 * \brief
 *    Read input from configuration file
 *
 * \par Input:
 *    Name of configuration filename
 *
 * \par Output
 *    none
 ************************************************************************
 */
void init_conf(struct inp_par *inp,
               char *config_filename)
{
  FILE *fd;
  int NAL_mode;

  // read the decoder configuration file
  if((fd=fopen(config_filename,"r")) == NULL)
  {
    snprintf(errortext, ET_SIZE, "Error: Control file %s not found\n",config_filename);
    error(errortext, 300);
  }

  fscanf(fd,"%s",inp->infile);                // H.26L compressed input bitsream
  fscanf(fd,"%*[^\n]");

  fscanf(fd,"%s",inp->outfile);               // YUV 4:2:2 input format
  fscanf(fd,"%*[^\n]");

  fscanf(fd,"%s",inp->reffile);               // reference file
  fscanf(fd,"%*[^\n]");

    // Symbol mode
  fscanf(fd,"%d,",&inp->symbol_mode);        // 0: UVLC 1: CABAC, may be overwritten ni case of RTP NAL
  fscanf(fd,"%*[^\n]");
  if (inp->symbol_mode != UVLC && inp->symbol_mode != CABAC)
  {
    snprintf(errortext, ET_SIZE, "Unsupported symbol mode=%d, use UVLC=0 or CABAC=1",inp->symbol_mode);
    error(errortext,1);
  }

  // UseConstrainedIntraPred
  fscanf(fd,"%d,",&inp->UseConstrainedIntraPred);        // 0: UsePred   1: ConstrainPred, may be overwritten in case of RTP NAL
  fscanf(fd,"%*[^\n]");
  if(inp->UseConstrainedIntraPred != 0 && inp->UseConstrainedIntraPred != 1)
  {
    snprintf(errortext, ET_SIZE, "Unsupported value=%d on constrained intra pred",inp->UseConstrainedIntraPred);
    error(errortext,1);
  }

  // Frame buffer size
  fscanf(fd,"%d,",&inp->buf_cycle);   // may be overwritten in case of RTP NAL
  fscanf(fd,"%*[^\n]");
  if (inp->buf_cycle < 1)
  {
    snprintf(errortext, ET_SIZE, "Frame Buffer Size is %d. It has to be at least 1",inp->buf_cycle);
    error(errortext,1);
  }

  fscanf(fd,"%d",&(NAL_mode));                // NAL mode
    fscanf(fd,"%*[^\n]");

  switch(NAL_mode)
  {
  case 0:
    inp->of_mode = PAR_OF_26L;
    // Note: Data Partitioning in 26L File Format not yet supported
    inp->partition_mode = PAR_DP_1;
    break;
  case 1:
    inp->of_mode = PAR_OF_RTP;
    inp->partition_mode = PAR_DP_3;         // DP_3 forces malloc_slcie to reserve memory
                                            // for three partitions.  In the RTP NAL, it can
                                            // be chanegd on a slice basis whether to use
                                            // one or three partitions
    break;
  case 2:
    inp->of_mode = PAR_OF_IFF;
    inp->partition_mode = PAR_DP_1;
    break;
  default:
    snprintf(errortext, ET_SIZE, "NAL mode %i is not supported", NAL_mode);
    error(errortext,400);
  }

#ifdef _LEAKYBUCKET_
  fscanf(fd,"%ld,",&inp->R_decoder);             // Decoder rate
  fscanf(fd, "%*[^\n]");
  fscanf(fd,"%ld,",&inp->B_decoder);             // Decoder buffer size
  fscanf(fd, "%*[^\n]");
  fscanf(fd,"%ld,",&inp->F_decoder);             // Decoder initial delay
  fscanf(fd, "%*[^\n]"); 
  fscanf(fd,"%s",inp->LeakyBucketParamFile);    // file where Leaky Bucket params (computed by encoder) are stored
  fscanf(fd,"%*[^\n]");
#endif

#ifndef _ABT_FLAG_IN_SLICE_HEADER_
  fscanf(fd,"%d,",&inp->abt); // Adaptive Block Transforms ABT
  fscanf(fd,"%*[^\n]");
  if ( (inp->abt < 0) || (inp->abt > 2) )
    {
      snprintf(errortext, ET_SIZE, "ABT Mode  %d not defined.",inp->abt);
      error(errortext,1);
    }
#endif

  // Loop Filter parameters flag
  fscanf(fd,"%d,",&inp->LFParametersFlag);   // 0: No Params  1: Read Filter Params, may be overwritten in case of RTP NAL
  fscanf(fd,"%*[^\n]");
  if(inp->LFParametersFlag != 0 && inp->LFParametersFlag != 1)
  {
    snprintf(errortext, ET_SIZE, "Unsupported value=%d on loop filter parameters flag",inp->LFParametersFlag);
    error(errortext,1);
  }


  fclose (fd);


#if TRACE
  if ((p_trace=fopen(TRACEFILE,"w"))==0)             // append new statistic at the end
  {
    snprintf(errortext, ET_SIZE, "Error open file %s!",TRACEFILE);
    error(errortext,500);
  }
#endif


  if (OpenBitstreamFile (inp->infile) < 0)
  {
    snprintf (errortext, ET_SIZE, "Cannot open bitstream file '%s'", inp->infile);
    error(errortext,500);
  }
  if ((p_out=fopen(inp->outfile,"wb"))==0)
  {
    snprintf(errortext, ET_SIZE, "Error open file %s ",inp->outfile);
    error(errortext,500);
  }

  fprintf(stdout,"--------------------------------------------------------------------------\n");
  fprintf(stdout," Decoder config file                    : %s \n",config_filename);
  fprintf(stdout,"--------------------------------------------------------------------------\n");
  fprintf(stdout," Input H.26L bitstream                  : %s \n",inp->infile);
  fprintf(stdout," Output decoded YUV 4:2:0               : %s \n",inp->outfile);
  fprintf(stdout," Output status file                     : %s \n",LOGFILE);
  if ((p_ref=fopen(inp->reffile,"rb"))==0)
  {
    fprintf(stdout," Input reference file                   : %s does not exist \n",inp->reffile);
    fprintf(stdout,"                                          SNR values are not available\n");
  }
  else
    fprintf(stdout," Input reference file                   : %s \n",inp->reffile);
#ifndef _ABT_FLAG_IN_SLICE_HEADER_
  if ( inp->abt )
  {
    fprintf(stdout," Adaptive Block Transforms              : Used ");
    if (inp->abt==INTER_ABT)
      fprintf(stdout,"(Inter only)\n");
    else
      fprintf(stdout,"(Inter and Intra)\n");
  }
  else
    fprintf(stdout," Adaptive Block Transforms              : Not used \n");
#endif
  fprintf(stdout,"--------------------------------------------------------------------------\n");
#ifdef _LEAKYBUCKET_
  fprintf(stdout," Rate_decoder        : %8ld \n",inp->R_decoder);
  fprintf(stdout," B_decoder           : %8ld \n",inp->B_decoder);
  fprintf(stdout," F_decoder           : %8ld \n",inp->F_decoder);
  fprintf(stdout," LeakyBucketParamFile: %s \n",inp->LeakyBucketParamFile); // Leaky Bucket Param file
  calc_buffer(inp);
  fprintf(stdout,"--------------------------------------------------------------------------\n");
#endif
  fprintf(stdout,"Frame    TR    QP  SnrY    SnrU    SnrV   Time(ms)\n");
}

/*!
 ************************************************************************
 * \brief
 *    Reports the gathered information to appropriate outputs
 *
 * \par Input:
 *    struct inp_par *inp,
 *    struct img_par *img,
 *    struct snr_par *stat
 *
 * \par Output:
 *    None
 ************************************************************************
 */
void report(struct inp_par *inp, struct img_par *img, struct snr_par *snr)
{
  #define OUTSTRING_SIZE 255
  char string[OUTSTRING_SIZE];
  FILE *p_log;

#ifndef WIN32
  time_t  now;
  struct tm *l_time;
#else
  char timebuf[128];
#endif

  fprintf(stdout,"-------------------- Average SNR all frames ------------------------------\n");
  fprintf(stdout," SNR Y(dB)           : %5.2f\n",snr->snr_ya);
  fprintf(stdout," SNR U(dB)           : %5.2f\n",snr->snr_ua);
  fprintf(stdout," SNR V(dB)           : %5.2f\n",snr->snr_va);
  fprintf(stdout," Total decoding time : %.3f sec \n",tot_time*0.001);
  fprintf(stdout,"--------------------------------------------------------------------------\n");
  fprintf(stdout," Exit JM %s decoder, ver %s ",JM,VERSION);
#if ( INI_CTX == 0 )
  fprintf(stdout,"No CABAC Initialization. ");
  fprintf(stdout," ABT_max_count %d ",INICNT_ABT);
#endif
  fprintf(stdout,"\n");
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
      fprintf(p_log," ------------------------------------------------------------------------------------------\n");
      fprintf(p_log,"|  Decoder statistics. This file is made first time, later runs are appended               |\n");
      fprintf(p_log," ------------------------------------------------------------------------------------------ \n");
      fprintf(p_log,"| Date  | Time  |    Sequence        |#Img|Format|SNRY 1|SNRU 1|SNRV 1|SNRY N|SNRU N|SNRV N|\n");
      fprintf(p_log," ------------------------------------------------------------------------------------------\n");
    }
  }
  else
  { 
    fclose(p_log);
    p_log=fopen(string,"a");                    // File exist,just open for appending
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
  fprintf(p_log,"| %1.5s |",string );
#endif

  fprintf(p_log,"%20.20s|",inp->infile);

  fprintf(p_log,"%3d |",img->number);

  fprintf(p_log,"%6.3f|",snr->snr_y1);
  fprintf(p_log,"%6.3f|",snr->snr_u1);
  fprintf(p_log,"%6.3f|",snr->snr_v1);
  fprintf(p_log,"%6.3f|",snr->snr_ya);
  fprintf(p_log,"%6.3f|",snr->snr_ua);
  fprintf(p_log,"%6.3f|\n",snr->snr_va);

  fclose(p_log);

  snprintf(string, OUTSTRING_SIZE,"%s", DATADECFILE);
  p_log=fopen(string,"a");

  if(Bframe_ctr != 0) // B picture used
  {
    fprintf(p_log, "%3d %2d %2d %2.2f %2.2f %2.2f %5d "
      "%2.2f %2.2f %2.2f %5d "
      "%2.2f %2.2f %2.2f %5d %.3f\n",
      img->number, 0, img->qp,
      snr->snr_y1,
      snr->snr_u1,
      snr->snr_v1,
      0,
      0.0,
      0.0,
      0.0,
      0,
      snr->snr_ya,
      snr->snr_ua,
      snr->snr_va,
      0,
      (double)0.001*tot_time/(img->number+Bframe_ctr-1));
  }
  else
  {
    fprintf(p_log, "%3d %2d %2d %2.2f %2.2f %2.2f %5d "
      "%2.2f %2.2f %2.2f %5d "
      "%2.2f %2.2f %2.2f %5d %.3f\n",
      img->number, 0, img->qp,
      snr->snr_y1,
      snr->snr_u1,
      snr->snr_v1,
      0,
      0.0,
      0.0,
      0.0,
      0,
      snr->snr_ya,
      snr->snr_ua,
      snr->snr_va,
      0,
      (double)0.001*tot_time/img->number);
  }
  fclose(p_log);
}

/*!
 ************************************************************************
 * \brief
 *    Allocates the slice structure along with its dependent
 *    data structures
 *
 * \par Input:
 *    Input Parameters struct inp_par *inp,  struct img_par *img
 ************************************************************************
 */
void malloc_slice(struct inp_par *inp, struct img_par *img)
{
  int i;
  DataPartition *dataPart;
  Slice *currSlice;
  const int buffer_size = MAX_CODED_FRAME_SIZE; // picture size unknown at this time, this value is to check

  switch(inp->of_mode) // init depending on NAL mode
  {
    case PAR_OF_IFF:
      // Current File Format
      img->currentSlice = (Slice *) calloc(1, sizeof(Slice));
      if ( (currSlice = img->currentSlice) == NULL)
      {
        snprintf(errortext, ET_SIZE, "Memory allocation for Slice datastruct in NAL-mode %d failed", inp->of_mode);
        error(errortext,100);
      }
      if (inp->symbol_mode == CABAC)
      {
        // create all context models
        currSlice->mot_ctx = create_contexts_MotionInfo();
        currSlice->tex_ctx = create_contexts_TextureInfo();
      }

      switch(inp->partition_mode)
      {
      case PAR_DP_1:
        currSlice->max_part_nr = 1;
        break;
      case PAR_DP_3:
        error("Data Partitioning Mode 3 in 26L-Format not supported",1);
        break;
      default:
        error("Data Partitioning Mode not supported!",1);
        break;
      }

      currSlice->partArr = (DataPartition *) calloc(1, sizeof(DataPartition));
      if (currSlice->partArr == NULL)
      {
        snprintf(errortext, ET_SIZE, "Memory allocation for Data Partition datastruct in NAL-mode %d failed", inp->of_mode);
        error(errortext, 100);
      }
      dataPart = currSlice->partArr;
      dataPart->bitstream = (Bitstream *) calloc(1, sizeof(Bitstream));
      if (dataPart->bitstream == NULL)
      {
        snprintf(errortext, ET_SIZE, "Memory allocation for Bitstream datastruct in NAL-mode %d failed", inp->of_mode);
        error(errortext, 100);
      }
      dataPart->bitstream->streamBuffer = (byte *) calloc(buffer_size, sizeof(byte));
      if (dataPart->bitstream->streamBuffer == NULL)
      {
        snprintf(errortext, ET_SIZE, "Memory allocation for bitstream buffer in NAL-mode %d failed", inp->of_mode);
        error(errortext, 100);
      }
      return;
    
    case PAR_OF_26L:
      // Current File Format
      img->currentSlice = (Slice *) calloc(1, sizeof(Slice));
      if ( (currSlice = img->currentSlice) == NULL)
      {
        snprintf(errortext, ET_SIZE, "Memory allocation for Slice datastruct in NAL-mode %d failed", inp->of_mode);
        error(errortext,100);
      }

      img->currentSlice->rmpni_buffer=NULL;

      if (inp->symbol_mode == CABAC)
      {
        // create all context models
        currSlice->mot_ctx = create_contexts_MotionInfo();
        currSlice->tex_ctx = create_contexts_TextureInfo();
      }

      switch(inp->partition_mode)
      {
      case PAR_DP_1:
        currSlice->max_part_nr = 1;
        break;
      case PAR_DP_3:
        error("Data Partitioning Mode 3 in 26L-Format not supported",1);
        break;
      default:
        error("Data Partitioning Mode not supported!",1);
        break;
      }

      currSlice->partArr = (DataPartition *) calloc(1, sizeof(DataPartition));
      if (currSlice->partArr == NULL)
      {
        snprintf(errortext, ET_SIZE, "Memory allocation for Data Partition datastruct in NAL-mode %d failed", inp->of_mode);
        error(errortext, 100);
      }
      dataPart = currSlice->partArr;
      dataPart->bitstream = (Bitstream *) calloc(1, sizeof(Bitstream));
      if (dataPart->bitstream == NULL)
      {
        snprintf(errortext, ET_SIZE, "Memory allocation for Bitstream datastruct in NAL-mode %d failed", inp->of_mode);
        error(errortext, 100);
      }
      dataPart->bitstream->streamBuffer = (byte *) calloc(buffer_size, sizeof(byte));
      if (dataPart->bitstream->streamBuffer == NULL)
      {
        snprintf(errortext, ET_SIZE, "Memory allocation for bitstream buffer in NAL-mode %d failed", inp->of_mode);
        error(errortext, 100);
      }
      return;
    
    case PAR_OF_RTP:
      img->currentSlice = (Slice *) calloc(1, sizeof(Slice));
      if ( (currSlice = img->currentSlice) == NULL)
      {
        snprintf(errortext, ET_SIZE, "Memory allocation for Slice datastruct in NAL-mode %d failed", inp->of_mode);
        error(errortext, 100);
      }

      if (inp->symbol_mode == CABAC)
      {
        // create all context models
        currSlice->mot_ctx = create_contexts_MotionInfo();
        currSlice->tex_ctx = create_contexts_TextureInfo();
      }

      switch(inp->partition_mode)
      {
      case PAR_DP_1:
        currSlice->max_part_nr = 1;
        break;
      case PAR_DP_3:
        currSlice->max_part_nr = 3;
        break;
      default:
        error("Data Partitioning Mode not supported!",1);
        break;
      }

      currSlice->partArr = (DataPartition *) calloc(3, sizeof(DataPartition));
      if (currSlice->partArr == NULL)
      {
        snprintf(errortext, ET_SIZE, "Memory allocation for Data Partition datastruct in NAL-mode %d failed", inp->of_mode);
        error(errortext, 100);
      }

      for (i=0; i<3; i++) // loop over all data partitions
      {
        dataPart = &(currSlice->partArr[i]);
        dataPart->bitstream = (Bitstream *) calloc(1, sizeof(Bitstream));
        if (dataPart->bitstream == NULL)
        {
          snprintf(errortext, ET_SIZE, "Memory allocation for Bitstream datastruct in NAL-mode %d failed", inp->of_mode);
          error(errortext, 100);
        }
        dataPart->bitstream->streamBuffer = (byte *) calloc(buffer_size, sizeof(byte));
        if (dataPart->bitstream->streamBuffer == NULL)
        {
          snprintf(errortext, ET_SIZE, "Memory allocation for bitstream buffer in NAL-mode %d failed", inp->of_mode);
          error(errortext, 100);
        }
      }
      return;

    default:
      snprintf(errortext, ET_SIZE, "Output File Mode %d not supported", inp->of_mode);
      error(errortext, 600);
  }

}

/*!
 ************************************************************************
 * \brief
 *    Memory frees of the Slice structure and of its dependent
 *    data structures
 *
 * \par Input:
 *    Input Parameters struct inp_par *inp,  struct img_par *img
 ************************************************************************
 */
void free_slice(struct inp_par *inp, struct img_par *img)
{
  int i;
  DataPartition *dataPart;
  Slice *currSlice = img->currentSlice;

  switch(inp->of_mode) // init depending on NAL mode
  {
    case PAR_OF_IFF:
      dataPart = currSlice->partArr;  // only one active data partition
      if (dataPart->bitstream->streamBuffer != NULL)
        free(dataPart->bitstream->streamBuffer);
      if (dataPart->bitstream != NULL)
        free(dataPart->bitstream);
      if (currSlice->partArr != NULL)
        free(currSlice->partArr);
      if (inp->symbol_mode == CABAC)
      {
        // delete all context models
        delete_contexts_MotionInfo(currSlice->mot_ctx);
        delete_contexts_TextureInfo(currSlice->tex_ctx);
      }
      if (currSlice != NULL)
        free(img->currentSlice);
      break;
    case PAR_OF_26L:
      // Current File Format
      dataPart = currSlice->partArr;  // only one active data partition
      if (dataPart->bitstream->streamBuffer != NULL)
        free(dataPart->bitstream->streamBuffer);
      if (dataPart->bitstream != NULL)
        free(dataPart->bitstream);
      if (currSlice->partArr != NULL)
        free(currSlice->partArr);
      if (inp->symbol_mode == CABAC)
      {
        // delete all context models
        delete_contexts_MotionInfo(currSlice->mot_ctx);
        delete_contexts_TextureInfo(currSlice->tex_ctx);
      }
      if (currSlice != NULL)
        free(img->currentSlice);
      break;
    case PAR_OF_RTP:
      // RTP File Format.
      // Here, mallocSLice is always called with 3 partitions, although sometimes only one is used
      for (i=0; i<3; i++) // loop over all data partitions
      {
        dataPart = &(currSlice->partArr[i]);
        if (dataPart->bitstream->streamBuffer != NULL)
          free(dataPart->bitstream->streamBuffer);
        if (dataPart->bitstream != NULL)
          free(dataPart->bitstream);
      }
      if (currSlice->partArr != NULL)
        free(currSlice->partArr);
      if (inp->symbol_mode == CABAC)
      {
        // delete all context models
        delete_contexts_MotionInfo(currSlice->mot_ctx);
        delete_contexts_TextureInfo(currSlice->tex_ctx);
      }
      if (currSlice != NULL)
        free(img->currentSlice);
      break;
    default:
      snprintf(errortext, ET_SIZE,  "Output File Mode %d not supported", inp->of_mode);
      error(errortext, 400);
  }

}

/*!
 ************************************************************************
 * \brief
 *    Dynamic memory allocation of frame size related global buffers
 *    buffers are defined in global.h, allocated memory must be freed in
 *    void free_global_buffers()
 *
 *  \par Input:
 *    Input Parameters struct inp_par *inp, Image Parameters struct img_par *img
 *
 *  \par Output:
 *     Number of allocated bytes
 ***********************************************************************
 */
int init_global_buffers(struct inp_par *inp, struct img_par *img)
{
  int i,j;

  int memory_size=0;
#ifdef _ADAPT_LAST_GROUP_
  extern int *last_P_no_frm;
  extern int *last_P_no_fld;
#endif

  img->buf_cycle = inp->buf_cycle+1;

  img->buf_cycle *= 2;

  if (img->structure != FRAME)
  {
    img->height *= 2;         // set height to frame (twice of field) for normal variables
    img->height_cr *= 2;      // set height to frame (twice of field) for normal variables
  }

#ifdef _ADAPT_LAST_GROUP_
  if ((last_P_no_frm = (int*)malloc(2*img->buf_cycle*sizeof(int))) == NULL)
    no_mem_exit("get_mem4global_buffers: last_P_no_frm");
  if ((last_P_no_fld = (int*)malloc(2*img->buf_cycle*sizeof(int))) == NULL)
    no_mem_exit("get_mem4global_buffers: last_P_no_fld");
#endif

  // allocate memory for encoding frame buffers: imgY, imgUV
  // byte imgY[288][352];
  // byte imgUV[2][144][176];
  memory_size += get_mem2D(&imgY_frm, img->height, img->width);    // processing memory in frame mode
  memory_size += get_mem3D(&imgUV_frm, 2, img->height_cr, img->width_cr); // processing memory in frame mode

  memory_size += get_mem2D(&imgY_top, img->height/2, img->width);    // processing memory in field mode
  memory_size += get_mem3D(&imgUV_top, 2, img->height_cr/2, img->width_cr);  // processing memory in field mode

  memory_size += get_mem2D(&imgY_bot, img->height/2, img->width);    // processing memory in field mode
  memory_size += get_mem3D(&imgUV_bot, 2, img->height_cr/2, img->width_cr);  // processing memory in field mode

  // allocate memory for multiple ref. frame buffers: mref, mcref
  // rows and cols for croma component mcef[ref][croma][4x][4y] are switched
  // compared to luma mref[ref][4y][4x] for whatever reason
  // number of reference frames increased by one for next P-frame
  alloc_mref(img);
  
  // allocate memory for imgY_prev
  memory_size += get_mem2D(&imgY_prev, img->height, img->width);
  memory_size += get_mem3D(&imgUV_prev, 2, img->height_cr, img->width_cr);

  // allocate memory for reference frames of each block: refFrArr
  // int  refFrArr[72][88];
  memory_size += get_mem2Dint(&refFrArr_frm, img->height/BLOCK_SIZE, img->width/BLOCK_SIZE);
  memory_size += get_mem2Dint(&refFrArr_top, img->height/BLOCK_SIZE, img->width/BLOCK_SIZE);
  memory_size += get_mem2Dint(&refFrArr_bot, img->height/BLOCK_SIZE, img->width/BLOCK_SIZE);
  
  // allocate memory for collocated motion stationarity - int could be replaced with boolean    
  memory_size += get_mem2Dint(&moving_block_frm, img->height/BLOCK_SIZE, img->width/BLOCK_SIZE);
  memory_size += get_mem2Dint(&moving_block_top,  img->height/BLOCK_SIZE, img->width/BLOCK_SIZE);
  memory_size += get_mem2Dint(&moving_block_bot,  img->height/BLOCK_SIZE, img->width/BLOCK_SIZE);

  // allocate memory for reference frame in find_snr
  // byte imgY_ref[288][352];
  // byte imgUV_ref[2][144][176];
  memory_size += get_mem2D(&imgY_ref, img->height, img->width);
  memory_size += get_mem3D(&imgUV_ref, 2, img->height_cr, img->width_cr);

  // allocate memory in structure img
  if(((img->mb_data) = (Macroblock *) calloc((img->width/MB_BLOCK_SIZE) * (img->height/MB_BLOCK_SIZE),sizeof(Macroblock))) == NULL)
    no_mem_exit("init_global_buffers: img->mb_data");
  if(img->UseConstrainedIntraPred)
  {
    if(((img->intra_block) = (int**)calloc((j=(img->width/MB_BLOCK_SIZE) * (img->height/MB_BLOCK_SIZE)),sizeof(int))) == NULL)
      no_mem_exit("init_global_buffers: img->intra_block");
    for (i=0; i<j; i++)
    {
      if ((img->intra_block[i] = (int*)calloc(4, sizeof(int))) == NULL)
        no_mem_exit ("init_global_buffers: img->intra_block");
    }
  }
  // img => int mv[92][72][3]
  memory_size += get_mem3Dint(&(img->mv_frm),img->width/BLOCK_SIZE +4, img->height/BLOCK_SIZE,3);
  memory_size += get_mem3Dint(&(img->mv_top),img->width/BLOCK_SIZE +4, img->height/BLOCK_SIZE,3);
  memory_size += get_mem3Dint(&(img->mv_bot),img->width/BLOCK_SIZE +4, img->height/BLOCK_SIZE,3);
  // img => int ipredmode[90][74]
  memory_size += get_mem2Dint(&(img->ipredmode),img->width/BLOCK_SIZE +2 , img->height/BLOCK_SIZE +2);
  // int dfMV[92][72][3];
  memory_size += get_mem3Dint(&(img->dfMV),img->width/BLOCK_SIZE +4, img->height/BLOCK_SIZE,3);
  // int dbMV[92][72][3];
  memory_size += get_mem3Dint(&(img->dbMV),img->width/BLOCK_SIZE +4, img->height/BLOCK_SIZE,3);
  // int fw_refFrArr[72][88];
  memory_size += get_mem2Dint(&(img->fw_refFrArr_frm),img->height/BLOCK_SIZE,img->width/BLOCK_SIZE);
  // int bw_refFrArr[72][88];
  memory_size += get_mem2Dint(&(img->bw_refFrArr_frm),img->height/BLOCK_SIZE,img->width/BLOCK_SIZE);
  // int fw_mv[92][72][3];
  memory_size += get_mem3Dint(&(img->fw_mv),img->width/BLOCK_SIZE +4, img->height/BLOCK_SIZE,3);
  // int bw_mv[92][72][3];
  memory_size += get_mem3Dint(&(img->bw_mv),img->width/BLOCK_SIZE +4, img->height/BLOCK_SIZE,3);

  memory_size += get_mem3Dint(&colB8mode,3,img->height/B8_SIZE, img->width/B8_SIZE); // collocated ABT block mode
  
  // int fw_refFrArr[72][88];
  memory_size += get_mem2Dint(&(img->fw_refFrArr_top),img->height/BLOCK_SIZE,img->width/BLOCK_SIZE);
  // int bw_refFrArr[72][88];
  memory_size += get_mem2Dint(&(img->bw_refFrArr_top),img->height/BLOCK_SIZE,img->width/BLOCK_SIZE);
  // int fw_refFrArr[72][88];
  memory_size += get_mem2Dint(&(img->fw_refFrArr_bot),img->height/BLOCK_SIZE,img->width/BLOCK_SIZE);
  // int bw_refFrArr[72][88];
  memory_size += get_mem2Dint(&(img->bw_refFrArr_bot),img->height/BLOCK_SIZE,img->width/BLOCK_SIZE);

  memory_size += get_mem2Dint(&(img->ipredmode_frm),img->width/BLOCK_SIZE +2 , img->height/BLOCK_SIZE +2);
  memory_size += get_mem2Dint(&(img->ipredmode_top),img->width/BLOCK_SIZE +2 , (img->height /2)/BLOCK_SIZE +2);
  memory_size += get_mem2Dint(&(img->ipredmode_bot),img->width/BLOCK_SIZE +2 , (img->height /2)/BLOCK_SIZE +2);

  memory_size += get_mem3Dint(&(img->fw_mv_frm),img->width/BLOCK_SIZE +4, img->height/BLOCK_SIZE,3);
  memory_size += get_mem3Dint(&(img->fw_mv_top),img->width/BLOCK_SIZE +4, (img->height/2)/BLOCK_SIZE,3);
  memory_size += get_mem3Dint(&(img->fw_mv_bot),img->width/BLOCK_SIZE +4, (img->height/2)/BLOCK_SIZE,3);

  memory_size += get_mem3Dint(&(img->bw_mv_frm),img->width/BLOCK_SIZE +4, img->height/BLOCK_SIZE,3);
  memory_size += get_mem3Dint(&(img->bw_mv_top),img->width/BLOCK_SIZE +4, (img->height/2)/BLOCK_SIZE,3);
  memory_size += get_mem3Dint(&(img->bw_mv_bot),img->width/BLOCK_SIZE +4, (img->height/2)/BLOCK_SIZE,3);

  memory_size += get_mem3Dint(&(img->dfMV_top),img->width/BLOCK_SIZE +4, img->height/BLOCK_SIZE,3);
  memory_size += get_mem3Dint(&(img->dbMV_top),img->width/BLOCK_SIZE +4, img->height/BLOCK_SIZE,3);

  memory_size += get_mem3Dint(&(img->dfMV_bot),img->width/BLOCK_SIZE +4, img->height/BLOCK_SIZE,3);
  memory_size += get_mem3Dint(&(img->dbMV_bot),img->width/BLOCK_SIZE +4, img->height/BLOCK_SIZE,3);

  memory_size += get_mem2Dint(&(img->field_anchor),img->height/BLOCK_SIZE,img->width/BLOCK_SIZE);
  memory_size += get_mem2Dint(&(field_mb), img->height/MB_BLOCK_SIZE, img->width/MB_BLOCK_SIZE);


  // CAVLC mem
  if((img->nz_coeff = (int****)calloc(img->width/MB_BLOCK_SIZE,sizeof(int***))) == NULL)
    no_mem_exit("get_mem4global_buffers: nzcoeff");
  for(j=0;j<img->width/MB_BLOCK_SIZE;j++)
  {
    memory_size += get_mem3Dint(&(img->nz_coeff[j]), img->height/MB_BLOCK_SIZE, 4, 6);
  }

  memory_size += get_mem2Dint(&(img->siblock),img->width/MB_BLOCK_SIZE  , img->height/MB_BLOCK_SIZE);

  if (img->structure != FRAME)
  {
    img->height /= 2;      // reset height for normal variables
    img->height_cr /= 2;   // reset height for normal variables
  }
  
  img->buf_cycle = inp->buf_cycle+1;

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
 *    Input Parameters struct inp_par *inp, Image Parameters struct img_par *img
 *
 * \par Output:
 *    none
 *
 ************************************************************************
 */
void free_global_buffers(struct inp_par *inp, struct img_par *img)
{
  int  i,j;
#ifdef _ADAPT_LAST_GROUP_
  extern int *last_P_no_frm;
  extern int *last_P_no_fld;
  free (last_P_no_frm);
  free (last_P_no_fld);
#endif

  free_mem2D(imgY_frm);
  free_mem2D(imgY_top);
  free_mem2D(imgY_bot);
  free_mem3D(imgUV_frm,2);
  free_mem3D(imgUV_top,2);
  free_mem3D(imgUV_bot,2);
  free_mem2D(imgY_prev);
  free_mem3D(imgUV_prev,2);

  // free multiple ref frame buffers
  free (mref_frm);
  free (mcef_frm);

  free (mref_fld);
  free (mcef_fld);

  free_mem2Dint(refFrArr_frm);
  free_mem2Dint(refFrArr_top);
  free_mem2Dint(refFrArr_bot);

  free_mem2Dint(moving_block_frm);
  free_mem2Dint(moving_block_top);
  free_mem2Dint(moving_block_bot);

  free_mem2D (imgY_ref);
  free_mem3D (imgUV_ref,2);
//  free_mem2D (imgY_tmp);
//  free_mem3D (imgUV_tmp,2);

  // CAVLC free mem
  for(j=0;j<img->width/MB_BLOCK_SIZE;j++)
  for(i=0;i<img->height/MB_BLOCK_SIZE;i++)
  {
    if(img->nz_coeff[j][i][0] != NULL) free(img->nz_coeff[j][i][0]);
    if(img->nz_coeff[j][i]    != NULL) free(img->nz_coeff[j][i]);
  };
  if (img->nz_coeff !=NULL) free(img->nz_coeff );
  free_mem2Dint(img->siblock);

  // free mem, allocated for structure img
  if (img->mb_data       != NULL) free(img->mb_data);

  if(img->UseConstrainedIntraPred)
  {
    j = (img->width/16)*(img->height/16);
    for (i=0; i<j; i++)
    {
      free (img->intra_block[i]);
    }
    free (img->intra_block);
  }
  free_mem3Dint(img->mv_frm,img->width/BLOCK_SIZE + 4);
  free_mem3Dint(img->mv_top,img->width/BLOCK_SIZE + 4);
  free_mem3Dint(img->mv_bot,img->width/BLOCK_SIZE + 4);

  free_mem2Dint (img->ipredmode);

  free_mem3Dint(img->dfMV,img->width/BLOCK_SIZE + 4);
  free_mem3Dint(img->dbMV,img->width/BLOCK_SIZE + 4);

  free_mem2Dint(img->fw_refFrArr_frm);
  free_mem2Dint(img->bw_refFrArr_frm);
  free_mem2Dint(img->fw_refFrArr_top);
  free_mem2Dint(img->bw_refFrArr_top);
  free_mem2Dint(img->fw_refFrArr_bot);
  free_mem2Dint(img->bw_refFrArr_bot);

  free_mem3Dint(img->fw_mv,img->width/BLOCK_SIZE + 4);
  free_mem3Dint(img->bw_mv,img->width/BLOCK_SIZE + 4);

  free_mem3Dint(colB8mode,3); // collocated ABT block mode

  free_mem2Dint (img->ipredmode_frm);
  free_mem2Dint (img->ipredmode_top);
  free_mem2Dint (img->ipredmode_bot);

  free_mem3Dint(img->fw_mv_frm,img->width/BLOCK_SIZE + 4);
  free_mem3Dint(img->fw_mv_top,img->width/BLOCK_SIZE + 4);
  free_mem3Dint(img->fw_mv_bot,img->width/BLOCK_SIZE + 4);

  free_mem3Dint(img->bw_mv_frm,img->width/BLOCK_SIZE + 4);
  free_mem3Dint(img->bw_mv_top,img->width/BLOCK_SIZE + 4);
  free_mem3Dint(img->bw_mv_bot,img->width/BLOCK_SIZE + 4);

  free_mem3Dint(img->dfMV_top,img->width/BLOCK_SIZE + 4);
  free_mem3Dint(img->dbMV_top,img->width/BLOCK_SIZE + 4);

  free_mem3Dint(img->dfMV_bot,img->width/BLOCK_SIZE + 4);
  free_mem3Dint(img->dbMV_bot,img->width/BLOCK_SIZE + 4);

  free_mem2Dint(img->field_anchor);
  free_mem2Dint(field_mb);

}

