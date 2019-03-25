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
 **************************************************************************************
 * \file
 *    slice.c
 * \brief
 *    generate the slice header, setup the bit buffer for slices,
 *    and generates the slice NALU(s)

 * \author
 *    Main contributors (see contributors.h for copyright, address and affiliation details)
 *      - Thomas Stockhammer            <stockhammer@ei.tum.de>
 *      - Detlev Marpe                  <marpe@hhi.de>
 *      - Stephan Wenger                <stewe@cs.tu-berlin.de>
 ***************************************************************************************
 */

#include "contributors.h"

#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/timeb.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include "global.h"
#include "header.h"
#include "rtp.h"
#include "nalu.h"
#include "annexb.h"
#include "parset.h"
#include "fmo.h"
#include "vlc.h"
#include "image.h"

//! Local declarations

static Slice *malloc_slice();
static void  free_slice(Slice *slice);
static void  init_slice();


/*!
 ************************************************************************
 * \brief
 *    init_ref_pic_list_reordering initializations should go here
 ************************************************************************
 */
void init_ref_pic_list_reordering()
{
  Slice* currSlice = img->currentSlice;

  currSlice->ref_pic_list_reordering_flag_l0 = 0;
  currSlice->ref_pic_list_reordering_flag_l1 = 0;
}


/*!
 ************************************************************************
 *  \brief
 *     This function generates the slice (and partition) header(s) 
 *
 *  \return number of bits used for the slice (and partition) header(s)
 *
 *  \par Side effects:
 *      Adds slice/partition header symbols to the symbol buffer
 *      increments Picture->no_slices, allocates memory for the
 *      slice, sets img->currSlice
 ************************************************************************
*/
int start_slice()
{
  EncodingEnvironmentPtr eep;
  Slice *currSlice = img->currentSlice;
  Bitstream *currStream;
  int header_len = 0;
  int i;
  int NumberOfPartitions = (input->partition_mode == PAR_DP_1?1:3);

  init_ref_pic_list_reordering();

  RTPUpdateTimestamp (img->tr);   // this has no side effects, just leave it for all NALs

  for (i=0; i<NumberOfPartitions; i++)
  {
    currStream = (currSlice->partArr[i]).bitstream;

    currStream->write_flag = 0;
    if (i==0)     // First partition
      header_len += SliceHeader (0);
    else          // Second/Third partition
      header_len += Partition_BC_Header(i);
     
    //! Initialize CABAC
    if (input->symbol_mode == CABAC)
    {
      eep = &((currSlice->partArr[i]).ee_cabac);
      if (currStream->bits_to_go != 8)
        header_len+=currStream->bits_to_go;
      writeVlcByteAlign(currStream);
      arienco_start_encoding(eep, currStream->streamBuffer, &(currStream->byte_pos)/*, &(currStream->last_startcode)*/,img->type);
    } else 
    {
      // Initialize CA-VLC
      CAVLC_init();
    }
  }
  if(input->symbol_mode == CABAC)
  {
    init_contexts();
  }
  return header_len;
}



/*!
 ************************************************************************
 * \brief
 *    This function terminates a slice (but doesn't write it out), 
 *    the old terminate_slice (0)
 * \return
 *    0 if OK,                                                         \n
 *    1 in case of error
 *
 ************************************************************************
 */
int terminate_slice()
{
  int bytes_written;
  Bitstream *currStream;
  Slice *currSlice = img->currentSlice;
  EncodingEnvironmentPtr eep;
  int i;
  int byte_pos_before_startcode_emu_prevention;

  if (input->symbol_mode == CABAC)
    write_terminating_bit (1);      // only once, not for all partitions
  
  for (i=0; i<currSlice->max_part_nr; i++)
  {
    currStream = (currSlice->partArr[i]).bitstream;
    if (input->symbol_mode == UVLC)
    {
      SODBtoRBSP(currStream);
      byte_pos_before_startcode_emu_prevention = currStream->byte_pos;
      currStream->byte_pos = RBSPtoEBSP(currStream->streamBuffer, 0 , currStream->byte_pos, 0);
      *(stat->em_prev_bits) += (currStream->byte_pos - byte_pos_before_startcode_emu_prevention) * 8;
    }
    else     // CABAC
    {
      eep = &((currSlice->partArr[i]).ee_cabac);
      // terminate the arithmetic code
      arienco_done_encoding(eep);
      currStream->bits_to_go = eep->Ebits_to_go;
      currStream->byte_buf = 0;
      bytes_written = currStream->byte_pos;
      byte_pos_before_startcode_emu_prevention= currStream->byte_pos;
      currStream->byte_pos = RBSPtoEBSP(currStream->streamBuffer, 0, currStream->byte_pos, eep->E);
      *(stat->em_prev_bits) += (currStream->byte_pos - byte_pos_before_startcode_emu_prevention) * 8;
    }           // CABAC
  }           // partition loop
  if( input->symbol_mode == CABAC )
  {
    store_contexts();
  }

  return 0;   
}


void SetStateVariablesForFrameMode()
{
  mcef = mcef_frm;      // set mcef to mcef_frm
  mref = mref_frm;      // set mref to mref_frm
  mref_w = mref_frm_w;  // set mref to mref_frm
  Refbuf11 = Refbuf11_frm;      // set Refbuff to frm
  Refbuf11_w = Refbuf11_frm_w;  // set Refbuff to frm
  img->height = input->img_height;      // set image height as frame height
}


void SetStateVariablesForFieldMode()
{
  mcef = mcef_fld;      // set mcef to mcef_frm
  mref = mref_mbfld;    // set mref to mref_frm
  mref_w = mref_mbfld_w;        // set mref to mref_frm
  Refbuf11 = Refbuf11_fld;      // set Refbuff to frm
  Refbuf11_w = Refbuf11_fld_w;  // set Refbuff to frm
}


/*!
 ************************************************************************
 * \brief
 *    Encodes one slice
 * \para
 *   returns the number of coded MBs in the SLice 
 ************************************************************************
 */
int encode_one_slice (int SliceGroupId, Picture *pic)
{
  Boolean end_of_slice = FALSE;
  Boolean recode_macroblock;
  int len, i, j;
  int NumberOfCodedMBs = 0;
  int CurrentMbInScanOrder;
  int MBRowSize = img->width / MB_BLOCK_SIZE;
  double FrameRDCost, FieldRDCost;

  img->cod_counter = 0;

  CurrentMbInScanOrder = FmoGetFirstMacroblockInSlice (SliceGroupId);
// printf ("\n\nEncode_one_slice: PictureID %d SliceGroupId %d  SliceID %d  FirstMB %d \n", img->tr, SliceGroupId, img->current_slice_nr, CurrentMbInScanOrder);

  set_MB_parameters (CurrentMbInScanOrder);
  init_slice (CurrentMbInScanOrder);
  Bytes_After_Header = img->currentSlice->partArr[0].bitstream->byte_pos;

  if (input->symbol_mode==CABAC)
  {
    SetCtxModelNumber ();
  }

/*
  // Tian Dong: June 7, 2002 JVT-B042
  // When the pictures are put into different layers and subseq, not all the reference frames
  // in multi-frame buffer are valid for prediction. The acutual number of the valid reference
  // frames, fb->num_short_used, will be given by start_slice(sym).
  // Save the fb->short_used.
  if (input->NumFramesInELSubSeq)
    {
      short_used = fb->short_used;
      img_ref = img->nb_references;
    }
*/

  len = start_slice ();

//  printf("short size, used, num-used: (%d,%d,%d)\n", fb->short_size, fb->short_used, fb->num_short_used);

/*
  // Tian Dong: June 7, 2002 JVT-B042
  if (input->NumFramesInELSubSeq)
    {
      fb->short_used = fb->num_short_used;
      img->nb_references = fb->short_used + fb->long_used;
    }
*/
  // Update statistics
  stat->bit_slice += len;
  stat->bit_use_header[img->type] += len;
// printf ("\n\n");

  while (end_of_slice == FALSE) // loop over macroblocks
  {
    if (!img->MbaffFrameFlag)
    {
      recode_macroblock = FALSE;
      set_MB_parameters (CurrentMbInScanOrder);
      start_macroblock ();
      encode_one_macroblock ();
      write_one_macroblock (1);
      terminate_macroblock (&end_of_slice, &recode_macroblock);

// printf ("encode_one_slice: mb %d,  slice %d,   bitbuf bytepos %d EOS %d\n", 
//       img->current_mb_nr, img->current_slice_nr, 
//       img->currentSlice->partArr[0].bitstream->byte_pos, end_of_slice);

      if (recode_macroblock == FALSE)       // The final processing of the macroblock has been done
      {
        CurrentMbInScanOrder = FmoGetNextMBNr (CurrentMbInScanOrder);
        if (CurrentMbInScanOrder == -1)   // end of slice
        {
// printf ("FMO End of Slice Group detected, current MBs %d, force end of slice\n", NumberOfCodedMBs+1);
          end_of_slice = TRUE;
        }
        NumberOfCodedMBs++;       // only here we are sure that the coded MB is actually included in the slice
        proceed2nextMacroblock (CurrentMbInScanOrder);
      }
      else
      {
        //! The statement below breaks obviously FMO.  I believe the correct statement would be
        //! img->current_mb_nr = CurrentMbInScanOrder;  
        //! It's now tested with an assert() (take it out if I'm wrong) and should be changed
        //! as soon as someone works on FMO
        
        img->current_mb_nr--;/*KS*/  
        assert (img->current_mb_nr == CurrentMbInScanOrder);
      }
    }
    else                      // TBD -- Addition of FMO
    {

//! This following ugly code breaks slices, at least for a slice mode that accumulates a certain
//! number of bits into one slice.  
//! The suggested algorithm is as follows:
//!
//! SaveState (Bitstream, stats,  etc. etc.);
//! BitsForThisMBPairInFrameMode = CodeMB (Upper, FRAME_MODE) + CodeMB (Lower, FRAME_MODE);
//! DistortionForThisMBPairInFrameMode = CalculateDistortion(Upper) + CalculateDistortion (Lower);
//! RestoreState();
//! BitsForThisMBPairInFieldMode = CodeMB (Upper, FIELD_MODE) + CodeMB (Lower, FIELD_MODE);
//! DistortionForThisMBPairInFrameMode = CalculateDistortion(Upper) + CalculateDistortion (Lower);
//! FrameFieldMode = Decision (...)
//! RestoreState()
//! if (FrameFieldMode == FRAME) {
//!   CodeMB (Upper, FRAME); CodeMB (Lower, FRAME);
//! } else {
//!   CodeMB (Upper FIELD); CodeMB (Lower, FIELD);
//! }
//!
//! Open questions/issues:
//!   1. CABAC/CA-VLC state:  It seems that the CABAC/CA_VLC states are changed during the
//!      dummy encoding processes (for the R-D based selection), but that they are never
//!      reset, once the selection is made.  I believe that this breaks the MB-adaptive
//!      frame/field coding.  The necessary code for the state saves is readily available
//!      in macroblock.c, start_macroblock() and terminate_macroblock() (this code needs
//!      to be double checked that it works with CA-VLC as well
//!   2. would it be an option to allocate Bitstreams with zero data in them (or copy the
//!      already generated bitstream) for the "test coding"?  


          // code MB pair as frame MB 
      recode_macroblock = FALSE;
      img->field_mode = 0;  // MB coded as frame
      img->top_field = 0;   // Set top field to 0

      set_MB_parameters (CurrentMbInScanOrder);
      start_macroblock ();
      
      SetStateVariablesForFrameMode();
      rdopt = &rddata_top_frame_mb; // store data in top frame MB 
      TopFrameIsSkipped = 0;
      WriteFrameFieldMBInHeader = 1;
      encode_one_macroblock ();     // code the MB as frame
      field_mb[img->mb_y][img->mb_x] = 0;   // set the MB as field (for use in FindSkipMotionVector)
      FrameRDCost = rdopt->min_rdcost;
      //***   Top MB coded as frame MB ***//
      
      // go to the bottom MB in the MB pair
      CurrentMbInScanOrder = img->current_mb_nr + MBRowSize;
      img->field_mode = 0;  // MB coded as frame  //GB
      set_MB_parameters (CurrentMbInScanOrder);
      start_macroblock ();
      rdopt = &rddata_bot_frame_mb; // store data in top frame MB
      WriteFrameFieldMBInHeader = TopFrameIsSkipped ? 1 : 0;
      field_mb[img->mb_y][img->mb_x] = 0;
      encode_one_macroblock ();     // code the MB as frame
      field_mb[img->mb_y][img->mb_x] = 0;   // set the MB as field (for use in FindSkipMotionVector)
      FrameRDCost += rdopt->min_rdcost;

      //***   Bottom MB coded as frame MB ***//
      
      
      // start coding the MB pair as a field MB pair
      CurrentMbInScanOrder -= MBRowSize;                //! FMO problem and generally dirty, just to go back like this
      img->field_mode = 1;  // MB coded as frame
      img->top_field = 1;   // Set top field to 1
      set_MB_parameters (CurrentMbInScanOrder);
      img->buf_cycle <<= 1;
      input->num_reference_frames <<= 1;
      img->num_ref_idx_l0_active <<= 1;
      img->num_ref_idx_l0_active += 1;
      start_macroblock ();
      
      
      
      img->height = input->img_height >> 1; 
      rdopt = &rddata_top_field_mb; // store data in top frame MB 
      SetStateVariablesForFieldMode();
      TopFieldIsSkipped = 0;        // set the top field MB skipped flag to 0
      WriteFrameFieldMBInHeader = 1;
      encode_one_macroblock ();     // code the MB as frame
      field_mb[img->mb_y][img->mb_x] = 1;   // set the MB as field (for use in FindSkipMotionVector)
      FieldRDCost = rdopt->min_rdcost;
      //***   Top MB coded as field MB ***//
      
      CurrentMbInScanOrder += MBRowSize;
      img->top_field = 0;   // Set top field to 0
      set_MB_parameters (CurrentMbInScanOrder);
      start_macroblock ();
      rdopt = &rddata_bot_field_mb; // store data in top frame MB 
      SetStateVariablesForFieldMode();
      WriteFrameFieldMBInHeader = TopFieldIsSkipped ? 1 : 0;
      encode_one_macroblock ();     // code the MB as frame
      field_mb[img->mb_y][img->mb_x] = 1;   // set the MB as field (for use in FindSkipMotionVector)
      FieldRDCost += rdopt->min_rdcost;
      //***   Bottom MB coded as field MB ***//
      
      // decide between frame/field MB pair
      if (FrameRDCost < FieldRDCost)
      {
        img->field_mode = 0;
        img->buf_cycle >>= 1;
        input->num_reference_frames >>= 1;
        MBPairIsField = 0;
        SetStateVariablesForFrameMode();
        img->num_ref_idx_l0_active -= 1;
        img->num_ref_idx_l0_active >>= 1;
        
      }
      else
      {
        img->field_mode = 1;
        MBPairIsField = 1;
        SetStateVariablesForFieldMode();
        img->height = input->img_height / 2;      // set image height as frame height
      }
      
      
      if (MBPairIsField)
        img->top_field = 1;
      else
        img->top_field = 0;
      
      // go back to the Top MB in the MB pair
      CurrentMbInScanOrder -= MBRowSize;
      set_MB_parameters (CurrentMbInScanOrder);
      start_macroblock ();
      for (i=0;i<4;i++)
        for (j=0;j<4;j++)
        {
          img->field_anchor[img->block_y+j][img->block_x+i] = img->field_mode;
        }
      rdopt =  img->field_mode ? &rddata_top_field_mb : &rddata_top_frame_mb;
      copy_rdopt_data (0);  // copy the MB data for Top MB from the temp buffers
      WriteFrameFieldMBInHeader = 1;
      write_one_macroblock (1);     // write the Top MB data to the bitstream
      NumberOfCodedMBs++;   // only here we are sure that the coded MB is actually included in the slice
      terminate_macroblock (&end_of_slice, &recode_macroblock);     // done coding the Top MB 
      proceed2nextMacroblock (CurrentMbInScanOrder);        // Go to next macroblock
      
      // go to the Bottom MB in the MB pair
      CurrentMbInScanOrder += MBRowSize;
      img->top_field = 0;
      set_MB_parameters (CurrentMbInScanOrder);
      start_macroblock ();

      for (i=0;i<4;i++)
        for (j=0;j<4;j++)
        {
          img->field_anchor[img->block_y+j][img->block_x+i] = img->field_mode;
        }
        rdopt = img->field_mode ? &rddata_bot_field_mb : &rddata_bot_frame_mb;
        copy_rdopt_data (1);  // copy the MB data for Bottom MB from the temp buffers
        if (img->field_mode)
          WriteFrameFieldMBInHeader = TopFieldIsSkipped ? 1 : 0;
        else
            WriteFrameFieldMBInHeader = TopFrameIsSkipped ? 1 : 0;

        write_one_macroblock (0);     // write the Bottom MB data to the bitstream
        NumberOfCodedMBs++;   // only here we are sure that the coded MB is actually included in the slice
        terminate_macroblock (&end_of_slice, &recode_macroblock);     // done coding the Top MB 
        proceed2nextMacroblock (CurrentMbInScanOrder);        // Go to next macroblock
        
        CurrentMbInScanOrder -= MBRowSize;
        
        if (MBPairIsField)    // if MB Pair was coded as field the buffer size variables back to frame mode
        {
          img->buf_cycle >>= 1;
          input->num_reference_frames >>= 1;
          img->num_ref_idx_l0_active -= 1;
          img->num_ref_idx_l0_active >>= 1;
        }
        img->field_mode = img->top_field = 0; // reset to frame mode
        img->height = input->img_height;      // reset the img->height  
        
        CurrentMbInScanOrder++;     //! Breaks FMO
        if (CurrentMbInScanOrder == img->total_number_mb - MBRowSize)
          end_of_slice = TRUE;        // just in case it does n't get set in terminate_macroblock
        
        if (CurrentMbInScanOrder % MBRowSize == 0)    //! Breaks FMO
          CurrentMbInScanOrder += MBRowSize;
        }
        
    }
/*
  // Tian Dong: June 7, 2002 JVT-B042
  // Restore the short_used
  if (input->NumFramesInELSubSeq)
    {
      fb->short_used = short_used;
      img->nb_references = img_ref;
    }
*/
  terminate_slice ();
  return NumberOfCodedMBs;
}



/*!
 ************************************************************************
 * \brief
 *    Initializes the parameters for a new slice and
 *     allocates the memory for the coded slice in the Picture structure
 *  \par Side effects:
 *      Adds slice/partition header symbols to the symbol buffer
 *      increments Picture->no_slices, allocates memory for the
 *      slice, sets img->currSlice
 ************************************************************************
 */
static void init_slice ()
{
  int i;
  Picture *currPic = img->currentPicture;
  Slice *curr_slice;
  DataPartition *dataPart;
  Bitstream *currStream;

  // Allocate new Slice in the current Picture, and set img->currentSlice
  assert (currPic != NULL);
  currPic->no_slices++;
  if (currPic->no_slices >= MAXSLICEPERPICTURE)
    error ("Too many slices per picture, incease MAXLSICESPERPICTURE in global.h, exitus", -1);
  currPic->slices[currPic->no_slices-1] = malloc_slice();
  curr_slice = currPic->slices[currPic->no_slices-1];
  img->currentSlice = curr_slice;

  curr_slice->picture_id = img->tr % 256;
  curr_slice->qp = img->qp;
  curr_slice->start_mb_nr = img->current_mb_nr;
  curr_slice->slice_too_big = dummy_slice_too_big;

  for (i = 0; i < curr_slice->max_part_nr; i++)
  {
    dataPart = &(curr_slice->partArr[i]);
    if (input->symbol_mode == UVLC)
      dataPart->writeSyntaxElement = writeSyntaxElement_UVLC;
    else
      dataPart->writeSyntaxElement = writeSyntaxElement_CABAC;
    
    currStream = dataPart->bitstream;
    currStream->bits_to_go = 8;
    currStream->byte_pos = 0;
    currStream->byte_buf = 0;
  }
}


/*!
 ************************************************************************
 * \brief
 *    Allocates a slice structure along with its dependentdata structures
 * \return
 *    Pointer to a Slice
 ************************************************************************
 */
static Slice *malloc_slice()
{
  int i;
  DataPartition *dataPart;
  Slice *slice;
  const int buffer_size = (img->width * img->height * 4); // AH 190202: There can be data expansion with 
                                                          // low QP values. So, we make sure that buffer 
                                                          // does not everflow. 4 is probably safe multiplier.

  if ((slice = (Slice *) calloc(1, sizeof(Slice))) == NULL) no_mem_exit ("malloc_slie: slice structure");

  if (input->symbol_mode == CABAC)
    {
      // create all context models
      slice->mot_ctx = create_contexts_MotionInfo();
      slice->tex_ctx = create_contexts_TextureInfo();
    }

  slice->max_part_nr = input->partition_mode==0?1:3;

  slice->num_mb = 0;          // no coded MBs so far

  if ((slice->partArr = (DataPartition *) calloc(slice->max_part_nr, sizeof(DataPartition))) == NULL) no_mem_exit ("malloc_slice: partArr");
  for (i=0; i<slice->max_part_nr; i++) // loop over all data partitions
  {
    dataPart = &(slice->partArr[i]);
    if ((dataPart->bitstream = (Bitstream *) calloc(1, sizeof(Bitstream))) == NULL) no_mem_exit ("malloc_slice: Bitstream");
    if ((dataPart->bitstream->streamBuffer = (byte *) calloc(buffer_size, sizeof(byte))) == NULL) no_mem_exit ("malloc_slice: StreamBuffer");
    // Initialize storage of bitstream parameters
  }
  return slice;
}


/*!
 ************************************************************************
 * \brief
 *    Memory frees of all Slice structures and of its dependent
 *    data structures
 * \par Input:
 *    Image Parameters struct struct img_par *img
 ************************************************************************
 */
void free_slice_list(Picture *currPic)
{
  int i;

  for (i=0; i<currPic->no_slices; i++)
  {
    free_slice (currPic->slices[i]);
    currPic->slices[i]=NULL;
  }
}


/*!
 ************************************************************************
 * \brief
 *    Memory frees of the Slice structure and of its dependent
 *    data structures
 * \para slice:
 *    Slice to be freed
 ************************************************************************
 */
static void free_slice(Slice *slice)
{
  int i;
  DataPartition *dataPart;

  if (slice != NULL)
  {
    for (i=0; i<slice->max_part_nr; i++) // loop over all data partitions
    {
      dataPart = &(slice->partArr[i]);
      if (dataPart != NULL)
      {
        if (dataPart->bitstream->streamBuffer != NULL)
          free(dataPart->bitstream->streamBuffer);
        if (dataPart->bitstream != NULL)
          free(dataPart->bitstream);
      }
    }
    if (slice->partArr != NULL)
      free(slice->partArr);
    if (input->symbol_mode == CABAC)
    {
      delete_contexts_MotionInfo(slice->mot_ctx);
      delete_contexts_TextureInfo(slice->tex_ctx);
    }
    //free(img->currentSlice);
    free(slice);
  }
}


/*!
 ************************************************************************
 * \brief
 *    This function set the value of a bit in a bitstream to 1
 ************************************************************************
 */
void modify_redundant_pic_cnt(unsigned char *buffer)
{
  unsigned char tmp = 1 << (rpc_bits_to_go-1);
  buffer[rpc_bytes_to_go] |= tmp;
}

