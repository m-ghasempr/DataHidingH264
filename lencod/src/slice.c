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
#include "encodeiff.h"
#include "sei.h"
#include "nalu.h"
#include "annexb.h"
#include "parset.h"
#include "fmo.h"
#include "vlc.h"

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
      init_contexts_MotionInfo (currSlice->mot_ctx);
      init_contexts_TextureInfo(currSlice->tex_ctx);
    } else 
    {
      // Initialize CA-VLC
      CAVLC_init();
    }
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
  //! is this correct considering that below in the partirtion loop the 
  //! terminating bit is written again???  StW 020103
  
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
//      write_terminating_bit (1);
      eep = &((currSlice->partArr[i]).ee_cabac);
      // terminate the arithmetic code
      arienco_done_encoding(eep);
      currStream->bits_to_go = eep->Ebits_to_go;
      currStream->byte_buf = 0;
      bytes_written = currStream->byte_pos /*-currStream->tmp_byte_pos*/;
      byte_pos_before_startcode_emu_prevention= currStream->byte_pos;
      currStream->byte_pos = RBSPtoEBSP(currStream->streamBuffer, 0, currStream->byte_pos, eep->E);
      *(stat->em_prev_bits) += (currStream->byte_pos - byte_pos_before_startcode_emu_prevention) * 8;
//        currStream->tmp_byte_pos = currStream->byte_pos;
    }           // CABAC
  }           // partition loop
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
  int len;
  int NumberOfCodedMBs = 0;
  int CurrentMbInScanOrder;
  int short_used = 0, img_ref = 0;
  int MBRowSize = img->width / MB_BLOCK_SIZE;
  double FrameRDCost, FieldRDCost;
  int LastCodedMB = -1;

  img->cod_counter = 0;
  CurrentMbInScanOrder = FmoGetFirstMacroblockInSlice (SliceGroupId);
// printf ("\n\nEncode_one_slice: PictureID %d SliceGroupId %d  SliceID %d  FirstMB %d \n", img->tr, SliceGroupId, img->current_slice_nr, CurrentMbInScanOrder);

  set_MB_parameters (CurrentMbInScanOrder);
  init_slice (CurrentMbInScanOrder);
  Bytes_After_Header = img->currentSlice->partArr[0].bitstream->byte_pos;

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
      if (input->InterlaceCodingOption < MB_CODING || mb_adaptive == 0)
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
//          img->update_stats = 0;        // don't update any stats yet                 //! This variable seems never be used
          
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
//          img->update_stats = 0;        // don't update any stats yet                 //! This variable seems never be used
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
          input->no_multpred <<= 1;
          img->num_ref_pic_active_fwd_minus1 <<= 1;
          img->num_ref_pic_active_fwd_minus1 += 1;
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
            input->no_multpred >>= 1;
            MBPairIsField = 0;
            SetStateVariablesForFrameMode();
            img->num_ref_pic_active_fwd_minus1 -= 1;
            img->num_ref_pic_active_fwd_minus1 >>= 1;
            
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
//          img->update_stats = 1;        // Now update the stats                 //! This variable seems never be used
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
//          img->update_stats = 1;        // Now update the stats                 //! This variable seems never be used
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
            input->no_multpred >>= 1;
            img->num_ref_pic_active_fwd_minus1 -= 1;
            img->num_ref_pic_active_fwd_minus1 >>= 1;
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


/*
      if (input->symbol_mode == UVLC)   
        {
          currStream = dataPart->bitstream;
          currStream->bits_to_go = currStream->stored_bits_to_go;
          currStream->byte_pos = currStream->stored_byte_pos;
          currStream->byte_buf = currStream->stored_byte_buf;
//          currStream->tmp_byte_pos = currStream->stored_byte_pos;       // temp store curr byte position
          if (currStream->byte_pos != 0)
            printf ("init_slice: byte_pos %d\n", currStream->byte_pos);
          if (currStream->bits_to_go != 8)
            printf ("init_slice: bits_to_go %d\n", currStream->byte_pos);
        }
      else                      // anything but UVLC and PAR_OF_26L
        {
          if (input->of_mode == PAR_OF_ANNEXB)
            {
              if ((img->current_slice_nr == 0) && (img->pstruct != 2) )
                {
                  currStream = dataPart->bitstream;
                  currStream->bits_to_go = 8;
                  currStream->byte_pos = 0;
                  currStream->byte_buf = 0;
//                  currStream->tmp_byte_pos = 0;
                }
            }
          else
            {
              currStream = dataPart->bitstream;
              currStream->bits_to_go = 8;
              currStream->byte_pos = 0;
              currStream->byte_buf = 0;
            }
        }
*/
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
 *    Memory frees of the Slice structure and of its dependent
 *    data structures
 * \par Input:
 *    Input Parameters struct inp_par *inp,  struct img_par *img
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
    free(img->currentSlice);
  }
}




// JVT-D101: redundant slices
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
// End JVT-D101







/*
//! some leftovers from the old terminate_slice(), no more used

    case PAR_OF_RTP:
      printf ("RTP not implemented, commented code needs more work\n");
      if (input->symbol_mode==CABAC)
      {
        write_terminating_bit (1);
      }

      for (i=0; i<currSlice->max_part_nr; i++)
      {
        currStream = (currSlice->partArr[i]).bitstream;

        if (input->symbol_mode == CABAC)
        {
          eep = &((currSlice->partArr[i]).ee_cabac);
          arienco_done_encoding(eep);
                                        currStream->bits_to_go = eep->Ebits_to_go;
                                        currStream->byte_buf = 0;
                                }
        
        if (currStream->bits_to_go < 8) // trailing bits to process
        {
          currStream->byte_buf <<= currStream->bits_to_go;
          stat->bit_use_stuffingBits[img->type]+=currStream->bits_to_go;
          currStream->streamBuffer[currStream->byte_pos++]=currStream->byte_buf;
          currStream->bits_to_go = 8;
        }
        bytes_written = currStream->byte_pos; // number of written bytes
        
        if(img->type == B_IMG && i > 0)    // to be really sure
          currStream->write_flag = 0;

        if(currStream->write_flag == 1)                                                   
        {
          int Marker;
          int FirstBytePacketType, subFirstBytePacketType;
         
          // Tian Dong (Sept 2002):
          // The SEI message will be encapsulated together with the slice data into a compound
          // packet (the newly name is "aggregation packet", see draft-wenger-avt-rtp-jvt-01.txt)
          FirstBytePacketType = subFirstBytePacketType = 0;
          if (isAggregationPacket())  // has spare reference picture information to be packetized
          {
            FirstBytePacketType = AGGREGATION_PACKET_TYPE; //! compound packet (aggregation packet) See VCEG-N72
            if ( FirstFrameIn2ndIGOP == img->number )
              subFirstBytePacketType = (PAYLOAD_TYPE_IDERP << 4);
            if (input->partition_mode==PAR_DP_1 || img->type == B_IMG)
              subFirstBytePacketType |= 0;
            else
              subFirstBytePacketType |= i+1; //! See VCEG-N72
          }
          else 
          {
            if ( FirstFrameIn2ndIGOP == img->number )
              FirstBytePacketType = (PAYLOAD_TYPE_IDERP << 4);
            if ( input->partition_mode==PAR_DP_1 || img->type == B_IMG)
              FirstBytePacketType |= 0;
            else
              FirstBytePacketType |= i+1; //! See VCEG-N72
          }
          
          // Calculate the Marker bit
          LastPartition=0;
          if (input->partition_mode == PAR_DP_3 && img->type != B_IMG)
          {
            if (currSlice->partArr[1].bitstream->write_flag == 1)  // something in partition 1
              LastPartition=1;
            if (currSlice->partArr[2].bitstream->write_flag == 1)  // something in partition 2
              LastPartition=2;
          }
          Marker = 0;
          if (img->current_mb_nr == img->total_number_mb)   // last MB is in this slice
            if ((input->partition_mode==PAR_DP_1) ||          // Single SLice
                (input->partition_mode==PAR_DP_3 && i == LastPartition ))   // Last partition containing bits
              Marker = 1;

          // and write the RTP packet
          currStream->byte_pos = bytes_written; 
          currStream->bits_to_go = 8;
          currStream->byte_buf = 0;

                                        if (input->symbol_mode==UVLC)  //GB stopbit is inserted in arienco_done_encoding for CABAC
                                                SODBtoRBSP(currStream);

          bytes_written = currStream->byte_pos;
          temp_byte_pos = bytes_written;

          if (input->symbol_mode==CABAC) 
          {
                  eep = &((currSlice->partArr[i]).ee_cabac);
                  bytes_written = RBSPtoEBSP(currStream->streamBuffer, Bytes_After_Header, bytes_written,eep->E);
          }
                else
                  bytes_written = RBSPtoEBSP(currStream->streamBuffer, Bytes_After_Header, bytes_written,0);

          *(stat->em_prev_bits) += (bytes_written - temp_byte_pos) * 8;

          // Tian Dong (Sept 2002):
          if (isAggregationPacket())
            rtp_bytes_written = aggregationRTPWriteBits (Marker, FirstBytePacketType, subFirstBytePacketType, currStream->streamBuffer, bytes_written, out);
          else
            rtp_bytes_written = RTPWriteBits (Marker, FirstBytePacketType, currStream->streamBuffer, bytes_written, out);

          // JVT-D101, write the redundant slice
          if(input->redundant_slice_flag)
          {
            // seems that SODBtoRBSP() and RBSPtoEBSP() did not change bits of slice header, 
            // therefore a bit position in slice header is unchanged.
            img->redundant_pic_cnt = 1;
            modify_redundant_pic_cnt(currStream->streamBuffer);
            rtp_bytes_written = RTPWriteBits (Marker, FirstBytePacketType, currStream->streamBuffer, bytes_written, out);
          }
          // End JVT-D101
        }
        // stat->bit_ctr += 8*bytes_written;
        stat->bit_ctr += 8*bytes_written*(img->redundant_pic_cnt+1); // JVT-D101: modified from the above line

        // Provide the next partition with a 'fresh' buffer
        currStream->stored_bits_to_go = 8;
        currStream->stored_byte_buf   = 0;
        currStream->stored_byte_pos   = 0;
        currStream->bits_to_go = 8;
        currStream->byte_buf   = 0;
        currStream->byte_pos   = 0;
      }
    return 0;
*/



/*!
 *****************************************************************************
 *
 * \brief 
 *    int RTPSliceHeader () write the Slice header for the RTP NAL      
 *
 * \return
 *    Number of bits used by the slice header
 *
 * \para Parameters
 *    none
 *
 * \para Side effects
 *    Slice header as per VCEG-N72r2 is written into the appropriate partition
 *    bit buffer
 *
 * \para Limitations/Shortcomings/Tweaks
 *    The current code does not support the change of picture parameters within
 *    one coded sequence, hence there is only one parameter set necessary.  This
 *    is hard coded to zero.
 *    As per Email discussions on 10/24/01 we decided to include the SliceID
 *    into both the Full Slice pakcet and Partition A packet
 *   
 * \date
 *    October 24, 2001
 *
 * \author
 *    Stephan Wenger   stewe@cs.tu-berlin.de
 *****************************************************************************/

/*
int RTPSliceHeader()
{
  int dP_nr;
  Bitstream *currStream; 
  DataPartition *partition;
  SyntaxElement sym;

  int len = 0;
  int RTPSliceType;

  if(img->type==BS_IMG && img->type!=B_IMG)
    dP_nr= assignSE2partition[input->partition_mode][SE_BFRAME];
  else
    dP_nr= assignSE2partition[input->partition_mode][SE_HEADER];

  currStream = ((img->currentSlice)->partArr[dP_nr]).bitstream;
  partition = &((img->currentSlice)->partArr[dP_nr]);

  assert (input->of_mode == PAR_OF_RTP);
  sym.type = SE_HEADER;         // This will be true for all symbols generated here
  sym.mapping = n_linfo2;       // Mapping rule: Simple code number to len/info


  SYMTRACESTRING("RTP-SH: Parameter Set");
  sym.value1 = 0;
  len += writeSyntaxElement_UVLC (&sym, partition);
  
  // picture structure, (0: progressive, 1: top field, 2: bottom field, 
  // 3: top field first, 4: bottom field first)
  sym.value1 = img->pstruct;
  SYMTRACESTRING("Picture Stucture");
  len += writeSyntaxElement_UVLC (&sym, partition);

  SYMTRACESTRING("RTP-SH: Picture ID");
  sym.value1 = img->currentSlice->picture_id;

  len += writeSyntaxElement_UVLC (&sym, partition);

  switch (img->type)
  {
    case INTER_IMG:
      if (img->types==SP_IMG)
        RTPSliceType = 2;
      else
        RTPSliceType = 0;
      break;
    case INTRA_IMG:
      RTPSliceType = 3;
      break;
    case B_IMG:
    case BS_IMG:
      RTPSliceType = 1;
      break;
    case SP_IMG:
      // That's what I would expect to be the case. But img->type contains INTER_IMG 
      // here (see handling above)
      //! \todo make one variable from img->type and img->types
      RTPSliceType = 2;
      break;
    default:
      printf ("Panic: unknown picture type %d, exit\n", img->type);
      exit (-1);
  }


  SYMTRACESTRING("RTP-SH: Slice Type");
  sym.value1 = RTPSliceType;
  len += writeSyntaxElement_UVLC (&sym, partition);

  // NOTE: following syntax elements have to be moved into nal_unit() and changed from e(v) to u(1)
  if(1)
  {
    sym.value1 = img->disposable_flag;
    SYMTRACESTRING("RTP-SH: NAL disposable_flag");
    len += writeSyntaxElement_UVLC (&sym, partition);
  }
  
  // NOTE: following syntax elements have to be moved into picture_layer_rbsp()
  // weighted_bipred_implicit_flag
  if(img->type==INTER_IMG)
  {
    sym.value1 = input->WeightedPrediction;
    SYMTRACESTRING("PH weighted_pred_flag");
    len += writeSyntaxElement_UVLC (&sym, partition);
  }

  if(img->type==B_IMG || img->type==BS_IMG)
  {
    sym.value1 = (input->WeightedBiprediction == 1)? 1 : 0;
    SYMTRACESTRING("PH weighted_bipred_explicit_flag");
    len += writeSyntaxElement_UVLC (&sym, partition);
  }
  if(img->type==B_IMG || img->type==BS_IMG)
  {
    sym.value1 = (input->WeightedBiprediction == 2)? 1 : 0;
    SYMTRACESTRING("PH weighted_bipred_implicit_flag");
    len += writeSyntaxElement_UVLC (&sym, partition);
  }

  if (img->type==INTER_IMG || img->type==B_IMG || img->type==BS_IMG)
  {
    sym.value1 = img->num_ref_pic_active_fwd_minus1;
    SYMTRACESTRING("RTP-SH: num_ref_pic_active_fwd_minus1");
    len += writeSyntaxElement_UVLC (&sym, partition);
    if (img->type==B_IMG || img->type==BS_IMG)
    {
      sym.value1 = img->num_ref_pic_active_bwd_minus1;
      SYMTRACESTRING("RTP-SH: num_ref_pic_active_bwd_minus1");
      len += writeSyntaxElement_UVLC (&sym, partition);
    }
  }

  SYMTRACESTRING("RTP-SH: FirstMBInSliceX");
  sym.value1 = img->mb_x;
  len += writeSyntaxElement_UVLC (&sym, partition);

  SYMTRACESTRING("RTP-SH: FirstMBinSlinceY");
  sym.value1 = img->mb_y;
  len += writeSyntaxElement_UVLC (&sym, partition);

  // JVT-D101: redundant_pic_cnt should be a ue(v) instead of u(1). However since currently we allow 
  // only 1 redundant slice for each slice, therefore herein encoding it as u(1) is OK.
  if(input->redundant_slice_flag) 
  {
    SYMTRACESTRING("RTP-SH: redundant_pic_cnt");
    rpc_bytes_to_go = currStream->byte_pos;
    rpc_bits_to_go = currStream->bits_to_go;
    sym.bitpattern = img->redundant_pic_cnt = 0; //may be modifed before being written into bitstream file
    sym.len = 1;
    len += writeSyntaxElement_fixed(&sym, partition);
  }
  // End JVT-D101

  if (img->type==B_IMG || img->type==BS_IMG)
  {
    SYMTRACESTRING("RTP-SH DirectSpatialFlag");
    sym.bitpattern = input->direct_type;  // 1 for Spatial Direct, 0 for temporal Direct
    sym.len = 1;
    len += writeSyntaxElement_fixed(&sym, partition);
  }


  SYMTRACESTRING("RTP-SH: InitialQP");
  sym.mapping = dquant_linfo;
  sym.value1 = img->qp - (MAX_QP - MIN_QP +1)/2;
  len += writeSyntaxElement_UVLC (&sym, partition);

  if (img->types==SP_IMG)
  {

    if (img->type!=INTRA_IMG) // Switch Flag only for SP pictures
    {
      SYMTRACESTRING("RTP-SH SWITCH FLAG");
      sym.bitpattern = 0;  // 1 for switching SP, 0 for normal SP
      sym.len = 1;
      len += writeSyntaxElement_fixed(&sym, partition);
    }

    SYMTRACESTRING("RTP-SH: SP InitialQP");
    sym.value1 = img->qpsp - (MAX_QP - MIN_QP +1)/2;
    len += writeSyntaxElement_UVLC (&sym, partition);
  }

  if (input->LFSendParameters)
  {
    SYMTRACESTRING("RTP-SH LF_DISABLE FLAG");
    sym.bitpattern = input->LFDisable;  
    sym.len = 1;
    len += writeSyntaxElement_fixed(&sym, partition);

    if (!input->LFDisable)
    {
      sym.mapping = dquant_linfo;           // Mapping rule: Signed integer
      SYMTRACESTRING("RTP-SH LFAlphaC0OffsetDiv2");
      sym.value1 = input->LFAlphaC0Offset>>1; 
      len += writeSyntaxElement_UVLC (&sym, partition);

      SYMTRACESTRING("RTP-SH LFBetaOffsetDiv2");
      sym.value1 = input->LFBetaOffset>>1; 
      len += writeSyntaxElement_UVLC (&sym, partition);
    }
  }

  sym.mapping = n_linfo2;


  if (input->partition_mode == PAR_DP_3 && img->type != B_IMG && img->type != BS_IMG)
  {
    SYMTRACESTRING("RTP-SH: Slice ID");
    sym.value1 = img->current_slice_nr;
    len += writeSyntaxElement_UVLC (&sym, partition);
  }

  len += writeERPS(&sym, partition);

  // JVT-D097
  if( input->num_slice_groups_minus1 > 0  &&  
      input->mb_allocation_map_type >= 4  &&  
                  input->mb_allocation_map_type <= 6)
  {
    SYMTRACESTRING("RTP-SH: slice_group_change_cycle");
    sym.value1 = slice_group_change_cycle;
    len += writeSyntaxElement_UVLC (&sym, partition);
  }
  // End JVT-D097

  // After this, and when in CABAC mode, terminateSlice() may insert one more
  // UVLC codeword for the number of coded MBs

  return len;
}
*/
