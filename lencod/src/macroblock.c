
/*!
 *************************************************************************************
 * \file macroblock.c
 *
 * \brief
 *    Process one macroblock
 *
 * \author
 *    Main contributors (see contributors.h for copyright, address and affiliation details)
 *    - Inge Lille-Langøy               <inge.lille-langoy@telenor.com>
 *    - Rickard Sjoberg                 <rickard.sjoberg@era.ericsson.se>
 *    - Jani Lainema                    <jani.lainema@nokia.com>
 *    - Sebastian Purreiter             <sebastian.purreiter@mch.siemens.de>
 *    - Detlev Marpe                    <marpe@hhi.de>
 *    - Thomas Wedi                     <wedi@tnt.uni-hannover.de>
 *    - Ragip Kurceren                  <ragip.kurceren@nokia.com>
 *************************************************************************************
 */
#include "contributors.h"

#include <math.h>
#include <stdlib.h>
#include <assert.h>

#include "elements.h"
#include "macroblock.h"
#include "refbuf.h"
#include "fmo.h"
#include "vlc.h"
#include "image.h"
#include "mb_access.h"

 /*!
 ************************************************************************
 * \brief
 *    updates the coordinates for the next macroblock to be processed
 *
 * \input
 *    mb: MB address in scan order
 ************************************************************************
 */

void set_MB_parameters (int mb)
{
  const int number_mb_per_row = img->width / MB_BLOCK_SIZE ;

  img->current_mb_nr = mb;
  img->mb_x = mb % number_mb_per_row;
  img->mb_y = mb / number_mb_per_row;

// printf ("Set_MB_Parameters: mb %d,  mb_x %d,  mb_y %d\n", mb, img->mb_x, img->mb_y);

  // Define vertical positions
  img->block_y = img->mb_y * BLOCK_SIZE;      // vertical luma block position
  img->pix_y   = img->mb_y * MB_BLOCK_SIZE;   // vertical luma macroblock position
  img->pix_c_y = img->mb_y * MB_BLOCK_SIZE/2; // vertical chroma macroblock position

  // Define horizontal positions
  img->block_x = img->mb_x * BLOCK_SIZE;        // luma block
  img->pix_x   = img->mb_x * MB_BLOCK_SIZE;     // luma pixel
  img->block_c_x = img->mb_x * BLOCK_SIZE/2;    // chroma block
  img->pix_c_x   = img->mb_x * MB_BLOCK_SIZE/2; // chroma pixel

  if(input->InterlaceCodingOption >= MB_CODING && mb_adaptive && img->field_mode)
  {
    if(img->top_field)
    {
      img->field_mb_y = img->mb_y/2;
      img->field_pix_y  = img->field_mb_y*MB_BLOCK_SIZE;
      img->field_block_y = img->field_mb_y*BLOCK_SIZE;
      img->field_pix_c_y = img->field_mb_y*MB_BLOCK_SIZE/2;
    }
    else
    {
      img->field_mb_y = (img->mb_y-1)/2;
      img->field_pix_y  = img->field_mb_y*MB_BLOCK_SIZE;
      img->field_block_y = img->field_mb_y*BLOCK_SIZE;
      img->field_pix_c_y = img->field_mb_y*MB_BLOCK_SIZE/2;
    } 
  }
}


int clip1a(int a)
{
  return ((a)>255?255:((a)<0?0:(a)));
}

/*!
 ************************************************************************
 * \brief
 *    updates the coordinates and statistics parameter for the
 *    next macroblock
 ************************************************************************
 */
void proceed2nextMacroblock()
{
#if TRACE
  int use_bitstream_backing = (input->slice_mode == FIXED_RATE || input->slice_mode == CALLBACK);
#endif
  Macroblock *currMB = &img->mb_data[img->current_mb_nr];
  int*        bitCount = currMB->bitcounter;

#if TRACE
  int i;
  if (p_trace)
  {
    fprintf(p_trace, "\n*********** Pic: %i (I/P) MB: %i Slice: %i **********\n\n", frame_no, img->current_mb_nr, img->current_slice_nr);
    if(use_bitstream_backing)
      fprintf(p_trace, "\n*********** Pic: %i (I/P) MB: %i Slice: %i **********\n\n", frame_no, img->current_mb_nr, img->current_slice_nr);
   // Write out the tracestring for each symbol
    for (i=0; i<currMB->currSEnr; i++)
      trace2out(&(img->MB_SyntaxElements[i]));
  }
#endif

  // Update the statistics
  stat->bit_use_mb_type[img->type]      += bitCount[BITS_MB_MODE];
  stat->bit_use_coeffY[img->type]       += bitCount[BITS_COEFF_Y_MB] ;
  stat->tmp_bit_use_cbp[img->type]      += bitCount[BITS_CBP_MB];
  stat->bit_use_coeffC[img->type]       += bitCount[BITS_COEFF_UV_MB];
  stat->bit_use_delta_quant[img->type]  += bitCount[BITS_DELTA_QUANT_MB];

  if (img->type==I_SLICE)
    ++stat->mode_use_intra[currMB->mb_type];
  else
    if (img->type != B_SLICE)
    {
      ++stat->mode_use_inter[0][currMB->mb_type];
      stat->bit_use_mode_inter[0][currMB->mb_type]+= bitCount[BITS_INTER_MB];

    }
    else
    {
      stat->bit_use_mode_inter[1][currMB->mb_type]+= bitCount[BITS_INTER_MB];
      ++stat->mode_use_inter[1][currMB->mb_type];
    }

  // Statistics
  if ((img->type == P_SLICE)||(img->type==SP_SLICE) )
  {
    ++stat->quant0;
    stat->quant1 += img->qp;      // to find average quant for inter frames
  }
}

/*!
 ************************************************************************
 * \brief
 *    initializes the current macroblock
 ************************************************************************
 */
void start_macroblock()
{
  int i,j,k,l;
  int use_bitstream_backing = (input->slice_mode == FIXED_RATE || input->slice_mode == CALLBACK);
  Macroblock *currMB = &img->mb_data[img->current_mb_nr];
  Slice *curr_slice = img->currentSlice;
  DataPartition *dataPart;
  Bitstream *currStream;
  EncodingEnvironmentPtr eep;

  if(use_bitstream_backing)
  {
    // Keep the current state of the bitstreams
    if(!img->cod_counter)
      for (i=0; i<curr_slice->max_part_nr; i++)
      {
        dataPart = &(curr_slice->partArr[i]);
        currStream = dataPart->bitstream;
        currStream->stored_bits_to_go   = currStream->bits_to_go;
        currStream->stored_byte_pos   = currStream->byte_pos;
        currStream->stored_byte_buf   = currStream->byte_buf;

        if (input->symbol_mode ==CABAC)
        {
          eep = &(dataPart->ee_cabac);
          eep->ElowS            = eep->Elow;
          eep->ErangeS           = eep->Erange;
          eep->EbufferS         = eep->Ebuffer;
          eep->Ebits_to_goS     = eep->Ebits_to_go;
          eep->Ebits_to_followS = eep->Ebits_to_follow;
          eep->EcodestrmS       = eep->Ecodestrm;
          eep->Ecodestrm_lenS   = eep->Ecodestrm_len;
          eep->CS               = eep->C;
          eep->BS               = eep->B;
          eep->ES               = eep->E;
        }
      }
  }

  // Save the slice number of this macroblock. When the macroblock below
  // is coded it will use this to decide if prediction for above is possible
  currMB->slice_nr = img->current_slice_nr;

    // Initialize delta qp change from last macroblock. Feature may be used for future rate control
  currMB->delta_qp = 0;
  currMB->qp       = img->qp;       // needed in loop filter (even if constant QP is used)

  // Initialize counter for MB symbols
  currMB->currSEnr=0;

  // If MB is next to a slice boundary, mark neighboring blocks unavailable for prediction
  CheckAvailabilityOfNeighbors(img);

  // Reset vectors before doing motion search in motion_search().
  if (img->type != B_SLICE)
  {
    for (k=0; k < 2; k++)
    {
      for (j=0; j < BLOCK_MULTIPLE; j++)
        for (i=0; i < BLOCK_MULTIPLE; i++)
          tmp_mv[k][img->block_y+j][img->block_x+i+4]=0;
    }
  }

  for (l=0; l<2; l++)
  {
    for (j=0; j < BLOCK_MULTIPLE; j++)
      for (i=0; i < BLOCK_MULTIPLE; i++)
        for (k=0; k < 2; k++)
          enc_picture->mv[l][img->block_x+i][img->block_y+j][k]=0;
  }

  //initialize reference index 
	for (j=0; j < BLOCK_MULTIPLE; j++)
  {
    for (i=0; i < BLOCK_MULTIPLE; i++)
      for (l=0; l<2; l++)
			{
        enc_picture->ref_idx[l][img->block_x+i][img->block_y + j] =-1;
        enc_picture->ref_pic_id[l][img->block_x+i][img->block_y+j] = -1;
      }
    }
  
  // Reset syntax element entries in MB struct
  currMB->mb_type   = 0;
  currMB->cbp_blk   = 0;
  currMB->cbp       = 0;
  currMB->mb_field  = 0;

  for (l=0; l < 2; l++)
    for (j=0; j < BLOCK_MULTIPLE; j++)
      for (i=0; i < BLOCK_MULTIPLE; i++)
        for (k=0; k < 2; k++)
          currMB->mvd[l][j][i][k] = 0;
 
  currMB->cbp_bits   = 0;
  currMB->c_ipred_mode = DC_PRED_8; //GB

  for (i=0; i < (BLOCK_MULTIPLE*BLOCK_MULTIPLE); i++)
    currMB->intra_pred_modes[i] = DC_PRED;

  //initialize the whole MB as INTRA coded
  //Blocks ar set to notINTRA in write_one_macroblock
  if (input->UseConstrainedIntraPred)
  {
    img->intra_block[img->current_mb_nr] = 1;
  }

  // store filtering parameters for this MB; For now, we are using the
  // same offset throughout the sequence
  currMB->lf_disable = input->LFDisableIdc;
  currMB->lf_alpha_c0_offset = input->LFAlphaC0Offset;
  currMB->lf_beta_offset = input->LFBetaOffset;


  // Initialize bitcounters for this macroblock
  if(img->current_mb_nr == 0) // No slice header to account for
  {
    currMB->bitcounter[BITS_HEADER] = 0;
  }
  else if (currMB->slice_nr == img->mb_data[img->current_mb_nr-1].slice_nr) // current MB belongs to the
  // same slice as the last MB
  {
    currMB->bitcounter[BITS_HEADER] = 0;
  }

  currMB->bitcounter[BITS_MB_MODE] = 0;
  currMB->bitcounter[BITS_COEFF_Y_MB] = 0;
  currMB->bitcounter[BITS_INTER_MB] = 0;
  currMB->bitcounter[BITS_CBP_MB] = 0;
  currMB->bitcounter[BITS_DELTA_QUANT_MB] = 0;
  currMB->bitcounter[BITS_COEFF_UV_MB] = 0;

#ifdef _FAST_FULL_ME_
  ResetFastFullIntegerSearch ();
#endif
}

/*!
 ************************************************************************
 * \brief
 *    terminates processing of the current macroblock depending
 *    on the chosen slice mode
 ************************************************************************
 */
void terminate_macroblock(Boolean *end_of_slice, Boolean *recode_macroblock)
{
  int i;
  Slice *currSlice = img->currentSlice;
  Macroblock    *currMB    = &img->mb_data[img->current_mb_nr];
  SyntaxElement *currSE    = &img->MB_SyntaxElements[currMB->currSEnr];
  int *partMap = assignSE2partition[input->partition_mode];
  DataPartition *dataPart;
  Bitstream *currStream;
  int rlc_bits=0;
  EncodingEnvironmentPtr eep;
  int use_bitstream_backing = (input->slice_mode == FIXED_RATE || input->slice_mode == CALLBACK);
  int new_slice = ( (img->current_mb_nr == 0) || img->mb_data[img->current_mb_nr-1].slice_nr != img->current_slice_nr);
  static int skip = FALSE;

  switch(input->slice_mode)
  {
  case NO_SLICES:
    currSlice->num_mb++;
    *recode_macroblock = FALSE;
    if ((currSlice->num_mb) == img->total_number_mb) // maximum number of MBs reached
      *end_of_slice = TRUE;
    break;
  case FIXED_MB:
    // For slice mode one, check if a new slice boundary follows
    currSlice->num_mb++;
    *recode_macroblock = FALSE;
    *end_of_slice = (currSlice->num_mb >= input->slice_argument);  // End of Slice
    *end_of_slice |= (img->current_mb_nr+1) == img->total_number_mb;    // End of Picture condition
    break;

    // For slice modes two and three, check if coding of this macroblock
    // resulted in too many bits for this slice. If so, indicate slice
    // boundary before this macroblock and code the macroblock again
  case FIXED_RATE:
     // in case of skip MBs check if there is a slice boundary
     // only for UVLC (img->cod_counter is always 0 in case of CABAC)
     if(img->cod_counter)
     {
       // write out the skip MBs to know how many bits we need for the RLC
       currSE->value1 = img->cod_counter;
       currSE->mapping = ue_linfo;
       currSE->type = SE_MBTYPE;
	
			 if (img->type != B_SLICE)   dataPart = &(currSlice->partArr[partMap[currSE->type]]);
       else                    dataPart = &(currSlice->partArr[partMap[SE_BFRAME]]);
       dataPart->writeSyntaxElement(  currSE, dataPart);
       rlc_bits=currSE->len;

       currStream = dataPart->bitstream;
       // save the bitstream as it would be if we write the skip MBs
       currStream->bits_to_go_skip  = currStream->bits_to_go;
       currStream->byte_pos_skip    = currStream->byte_pos;
       currStream->byte_buf_skip    = currStream->byte_buf;
       // restore the bitstream
       currStream->bits_to_go = currStream->stored_bits_to_go;
       currStream->byte_pos = currStream->stored_byte_pos;
       currStream->byte_buf = currStream->stored_byte_buf;
       skip = TRUE;
     }
     //! Check if the last coded macroblock fits into the size of the slice
     //! But only if this is not the first macroblock of this slice
     if (!new_slice)
     {
       if(slice_too_big(rlc_bits))
       {
         *recode_macroblock = TRUE;
         *end_of_slice = TRUE;
       }
       else if(!img->cod_counter)
         skip = FALSE;
     }
     // maximum number of MBs
     if ((*recode_macroblock == FALSE) && ((img->current_mb_nr+1) == img->total_number_mb)) 
     {
       *end_of_slice = TRUE;
       if(!img->cod_counter)
         skip = FALSE;
     }
   
     //! (first MB OR first MB in a slice) AND bigger that maximum size of slice
     if (new_slice && slice_too_big(rlc_bits))
     {
       *end_of_slice = TRUE;
       if(!img->cod_counter)
         skip = FALSE;
     }
     if (!*recode_macroblock)
       currSlice->num_mb++;
     break;

  case  CALLBACK:
    if (img->current_mb_nr > 0 && !new_slice)
    {
      if (currSlice->slice_too_big(rlc_bits))
      {
        *recode_macroblock = TRUE;
        *end_of_slice = TRUE;
      }
    }
    if ( (*recode_macroblock == FALSE) && ((img->current_mb_nr+1) == img->total_number_mb) ) // maximum number of MBs
      *end_of_slice = TRUE;
    break;

  case FMO:
    // The FMO slice mode acts like slice mode 1 (fixed maximum #of mbs per slice, in slice_argument)

    currSlice->num_mb++;
    *recode_macroblock = FALSE;
    // Check end-of-slice group condition first
    *end_of_slice = (img->current_mb_nr == FmoGetLastCodedMBOfSliceGroup (FmoMB2SliceGroup (img->current_mb_nr)));
    // Now check maximum # of MBs in slice
    *end_of_slice |= (currSlice->num_mb >= input->slice_argument);
    break;

  default:
    snprintf(errortext, ET_SIZE, "Slice Mode %d not supported", input->slice_mode);
    error(errortext, 600);
  }

  if(*recode_macroblock == TRUE)
  {
    // Restore everything
    for (i=0; i<currSlice->max_part_nr; i++)
    {
      dataPart = &(currSlice->partArr[i]);
      currStream = dataPart->bitstream;
      currStream->bits_to_go = currStream->stored_bits_to_go;
      currStream->byte_pos  = currStream->stored_byte_pos;
      currStream->byte_buf  = currStream->stored_byte_buf;
      if (input->symbol_mode == CABAC)
      {
        eep = &(dataPart->ee_cabac);
        eep->Elow            = eep->ElowS;
        eep->Erange           = eep->ErangeS;
        eep->Ebuffer         = eep->EbufferS;
        eep->Ebits_to_go     = eep->Ebits_to_goS;
        eep->Ebits_to_follow = eep->Ebits_to_followS;
        eep->Ecodestrm       = eep->EcodestrmS;
        eep->Ecodestrm_len   = eep->Ecodestrm_lenS;
        eep->C               = eep->CS;
        eep->B               = eep->BS;
        eep->E               = eep->ES;       
      }
    }
  }

  if(*end_of_slice == TRUE  && skip == TRUE) //! TO 4.11.2001 Skip MBs at the end of this slice
  { 
    //! only for Slice Mode 2 or 3
    // If we still have to write the skip, let's do it!
    if(img->cod_counter && *recode_macroblock == TRUE) //! MB that did not fit in this slice
    { 
      // If recoding is true and we have had skip, 
      // we have to reduce the counter in case of recoding
      img->cod_counter--;
      if(img->cod_counter)
      {
        currSE->value1 = img->cod_counter;
        currSE->mapping = ue_linfo;
        currSE->type = SE_MBTYPE;

				if (img->type != B_SLICE)   dataPart = &(currSlice->partArr[partMap[currSE->type]]);
        else                                             dataPart = &(currSlice->partArr[partMap[SE_BFRAME]]);
        dataPart->writeSyntaxElement(  currSE, dataPart);
        rlc_bits=currSE->len;
        currMB->bitcounter[BITS_MB_MODE]+=rlc_bits;
        img->cod_counter = 0;
      }
    }
    else //! MB that did not fit in this slice anymore is not a Skip MB
    {
			if (img->type != B_SLICE)   dataPart = &(currSlice->partArr[partMap[SE_MBTYPE]]);
      else                                             dataPart = &(currSlice->partArr[partMap[SE_BFRAME]]);
       
      currStream = dataPart->bitstream;
        // update the bitstream
      currStream->bits_to_go = currStream->bits_to_go_skip;
      currStream->byte_pos  = currStream->byte_pos_skip;
      currStream->byte_buf  = currStream->byte_buf_skip;

      // update the statistics
      img->cod_counter = 0;
      skip = FALSE;
    }
  }
  
  //! TO 4.11.2001 Skip MBs at the end of this slice for Slice Mode 0 or 1
  if(*end_of_slice == TRUE && img->cod_counter && !use_bitstream_backing)
  {
    currSE->value1 = img->cod_counter;
    currSE->mapping = ue_linfo;
    currSE->type = SE_MBTYPE;
		if (img->type != B_SLICE)   dataPart = &(currSlice->partArr[partMap[currSE->type]]);
    else                                             dataPart = &(currSlice->partArr[partMap[SE_BFRAME]]);
    dataPart->writeSyntaxElement(  currSE, dataPart);
    rlc_bits=currSE->len;
    currMB->bitcounter[BITS_MB_MODE]+=rlc_bits;
    img->cod_counter = 0;
  }
}

/*!
 *****************************************************************************
 *
 * \brief 
 *    For Slice Mode 2: Checks if one partition of one slice exceeds the 
 *    allowed size
 * 
 * \return
 *    FALSE if all Partitions of this slice are smaller than the allowed size
 *    TRUE is at least one Partition exceeds the limit
 *
 * \para Parameters
 *    
 *    
 *
 * \para Side effects
 *    none
 *
 * \para Other Notes
 *    
 *    
 *
 * \date
 *    4 November 2001
 *
 * \author
 *    Tobias Oelbaum      drehvial@gmx.net
 *****************************************************************************/
 
 int slice_too_big(int rlc_bits)
 {
   Slice *currSlice = img->currentSlice;
   DataPartition *dataPart;
   Bitstream *currStream;
   EncodingEnvironmentPtr eep;
   int i;
   int size_in_bytes;
  
   //! UVLC
   if (input->symbol_mode == UVLC)
   {
     for (i=0; i<currSlice->max_part_nr; i++)
     {
       dataPart = &(currSlice->partArr[i]);
       currStream = dataPart->bitstream;
       size_in_bytes = currStream->byte_pos /*- currStream->tmp_byte_pos*/;

       if (currStream->bits_to_go < 8)
         size_in_bytes++;
       if (currStream->bits_to_go < rlc_bits)
         size_in_bytes++;
       if(size_in_bytes > input->slice_argument)
         return TRUE;
     }
   }
    
   //! CABAC
   if (input->symbol_mode ==CABAC)
   {
     for (i=0; i<currSlice->max_part_nr; i++)
     {
        dataPart= &(currSlice->partArr[i]);
        eep = &(dataPart->ee_cabac);
      
       if( arienco_bits_written(eep) > (input->slice_argument*8))
          return TRUE;
     }
   }
   return FALSE;
 }
/*!
 ************************************************************************
 * \brief
 *    Checks the availability of neighboring macroblocks of
 *    the current macroblock for prediction and context determination;
 *    marks the unavailable MBs for intra prediction in the
 *    ipredmode-array by -1. Only neighboring MBs in the causal
 *    past of the current MB are checked.
 ************************************************************************
 */
 /*
void CheckAvailabilityOfNeighbors()
{
  int i,j;
  const int mb_width = img->width/MB_BLOCK_SIZE;
  const int mb_nr = img->current_mb_nr;
  Macroblock *currMB = &img->mb_data[mb_nr];
  int   pix_y   = img->pix_y;   // For MB level Frame/field coding
  int   block_y = img->block_y;   // For MB level Frame/field coding
  
  if(input->InterlaceCodingOption >= MB_CODING && mb_adaptive && img->field_mode)
  {
    pix_y = img->field_pix_y;
    block_y = img->field_block_y;
  }

  // mark all neighbors as unavailable
  for (i=0; i<3; i++)
    for (j=0; j<3; j++)
    {
      img->mb_data[mb_nr].mb_available[i][j]=NULL;
    }

  img->mb_data[mb_nr].mb_available[1][1]=currMB; // current MB
  

  // Check MB to the left
  if(img->pix_x >= MB_BLOCK_SIZE)
  {
    int remove_prediction = currMB->slice_nr != img->mb_data[mb_nr-1].slice_nr;
    // upper blocks
    if (remove_prediction || (input->UseConstrainedIntraPred && img->intra_block[mb_nr-1][1]==0))
    {
      img->ipredmode[img->block_x][img->block_y+1] = -1;
      img->ipredmode[img->block_x][img->block_y+2] = -1;
      if(input->InterlaceCodingOption >= MB_CODING && mb_adaptive && img->field_mode)  // GB
      {
        if(img->top_field)
        {
          img->ipredmode_top[img->block_x][block_y+1] = -1;
          img->ipredmode_top[img->block_x][block_y+2] = -1;
        }
        else
        {
          img->ipredmode_bot[img->block_x][block_y+1] = -1;
          img->ipredmode_bot[img->block_x][block_y+2] = -1;
        }
      }
    }
    // lower blocks
    if (remove_prediction || (input->UseConstrainedIntraPred && img->intra_block[mb_nr-1][3]==0))
    {
      img->ipredmode[img->block_x][img->block_y+3] = -1;
      img->ipredmode[img->block_x][img->block_y+4] = -1;
      if(input->InterlaceCodingOption >= MB_CODING && mb_adaptive && img->field_mode) //GB
      {
        if(img->top_field)
        {
          img->ipredmode_top[img->block_x][block_y+3] = -1;
          img->ipredmode_top[img->block_x][block_y+4] = -1;
        }
        else
        {
          img->ipredmode_bot[img->block_x][block_y+3] = -1;
          img->ipredmode_bot[img->block_x][block_y+4] = -1;
        }
      }
    }
    if (!remove_prediction)
    {
      currMB->mb_available[1][0]=&(img->mb_data[mb_nr-1]);
    }
  }

  // Check MB above
  if(pix_y >= MB_BLOCK_SIZE) // wrong for MBAFF
  //if(img->pix_y >= MB_BLOCK_SIZE)
  {
    int remove_prediction = currMB->slice_nr != img->mb_data[mb_nr-mb_width].slice_nr;
    // upper blocks
    if (remove_prediction || (input->UseConstrainedIntraPred && img->intra_block[mb_nr-mb_width][2]==0))
    {
      img->ipredmode[img->block_x+1][img->block_y] = -1;
      img->ipredmode[img->block_x+2][img->block_y] = -1;
      if(input->InterlaceCodingOption >= MB_CODING && mb_adaptive && img->field_mode) //GB 
      {
        if(img->top_field)
        {
          img->ipredmode_top[img->block_x+1][block_y] = -1;
          img->ipredmode_top[img->block_x+2][block_y] = -1;
        }
        else
        {
          img->ipredmode_bot[img->block_x+1][block_y] = -1;
          img->ipredmode_bot[img->block_x+2][block_y] = -1;
        }
      }
    }
    // lower blocks
    if (remove_prediction || (input->UseConstrainedIntraPred && img->intra_block[mb_nr-mb_width][3]==0))
    {
      img->ipredmode[img->block_x+3][img->block_y] = -1;
      img->ipredmode[img->block_x+4][img->block_y] = -1;
      if(input->InterlaceCodingOption >= MB_CODING && mb_adaptive && img->field_mode) //GB
      {
        if(img->top_field)
        {
          img->ipredmode_top[img->block_x+3][block_y] = -1;
          img->ipredmode_top[img->block_x+4][block_y] = -1;
        }
        else
        {
          img->ipredmode_bot[img->block_x+3][block_y] = -1;
          img->ipredmode_bot[img->block_x+4][block_y] = -1;
        } 
      }
    }
    if (!remove_prediction)
    {
      currMB->mb_available[0][1]=&(img->mb_data[mb_nr-mb_width]);
    }
  }

  // Check MB left above
  if(img->pix_x >= MB_BLOCK_SIZE && img->pix_y >= MB_BLOCK_SIZE )
  {
    int remove_prediction = currMB->slice_nr != img->mb_data[mb_nr-mb_width-1].slice_nr;

    if (remove_prediction || (input->UseConstrainedIntraPred && img->intra_block[mb_nr-mb_width-1][3]==0))
    {
      img->ipredmode[img->block_x][img->block_y] = -1;
      if(input->InterlaceCodingOption >= MB_CODING && mb_adaptive && img->field_mode) //GB 
      {
        if(img->top_field)
        {
          img->ipredmode_top[img->block_x][block_y] = -1;
        }
        else
        {
          img->ipredmode_bot[img->block_x][block_y] = -1;
        }
      }
    }
    if (!remove_prediction)
    {
      currMB->mb_available[0][0]=&(img->mb_data[mb_nr-mb_width-1]);
    }
  }

  // Check MB right above
  if(pix_y >= MB_BLOCK_SIZE && img->pix_x < (img->width-MB_BLOCK_SIZE ))
  {
    if(currMB->slice_nr == img->mb_data[mb_nr-mb_width+1].slice_nr)
      currMB->mb_available[0][2]=&(img->mb_data[mb_nr-mb_width+1]);
  }

  if(input->InterlaceCodingOption >= MB_CODING && mb_adaptive)
    currMB->mb_available[0][2] = NULL;  // set the prediction from top right MB to zero 

}

*/
/*!
 ************************************************************************
 * \brief
 *    Predict one component of a 4x4 Luma block
 ************************************************************************
 */
void
OneComponentLumaPrediction4x4 (int*   mpred,      //  --> array of prediction values (row by row)
                               int    pic_pix_x,  // <--  absolute horizontal coordinate of 4x4 block
                               int    pic_pix_y,  // <--  absolute vertical   coordinate of 4x4 block
                               int*   mv,         // <--  motion vector
                               int    ref)        // <--  reference frame (0.. / -1:backward)
{
  int incr;
	int block[BLOCK_SIZE][BLOCK_SIZE];
  int     pix_add = 4;
  int     j0      = (pic_pix_y << 2) + mv[1], j1=j0+pix_add, j2=j1+pix_add, j3=j2+pix_add;
  int     i0      = (pic_pix_x << 2) + mv[0], i1=i0+pix_add, i2=i1+pix_add, i3=i2+pix_add;

  pel_t (*get_pel) (pel_t**, int, int) = UMVPelY_14;

  // Tian Dong: PLUS1, June 06, 2002
  incr      = (ref==-1 ? (!img->fld_type&&enc_picture!=enc_frame_picture): direct_mode ? (!img->fld_type&&enc_picture!=enc_frame_picture) : (enc_picture!=enc_frame_picture)) ;
  
  if(input->InterlaceCodingOption >= MB_CODING && mb_adaptive)
    incr    = (ref==-1 ? (img->top_field&&img->field_mode): direct_mode ? (img->top_field&&img->field_mode):img->field_mode);

	get_block(ref<0 ? 0:ref,ref<0 ? listX[1]:listX[0], i0,j0,block);
	//get_block(ref<0 ? 0:ref,ref<0 ? listX[1]:listX[0], i0,j0,block);
 
  *mpred++ = block[0][0];

  *mpred++ = block[1][0];
  *mpred++ = block[2][0];
  *mpred++ = block[3][0];
  *mpred++ = block[0][1];
  *mpred++ = block[1][1];
  *mpred++ = block[2][1];
  *mpred++ = block[3][1];
  *mpred++ = block[0][2];
  *mpred++ = block[1][2];
  *mpred++ = block[2][2];
  *mpred++ = block[3][2];
  *mpred++ = block[0][3];
  *mpred++ = block[1][3];
  *mpred++ = block[2][3];
  *mpred++ = block[3][3];
}

/*!
 ************************************************************************
 * \brief
 *    copy foward/backward prediction values of one component of a 4x4 Luma block
 ************************************************************************
 */

void
copyblock4x4 (int*   mpred,      //  --> array of prediction values (row by row)
              int block[BLOCK_SIZE][BLOCK_SIZE])        
{
  *mpred++ = block[0][0];
  *mpred++ = block[1][0];
  *mpred++ = block[2][0];
  *mpred++ = block[3][0];
  *mpred++ = block[0][1];
  *mpred++ = block[1][1];
  *mpred++ = block[2][1];
  *mpred++ = block[3][1];
  *mpred++ = block[0][2];
  *mpred++ = block[1][2];
  *mpred++ = block[2][2];
  *mpred++ = block[3][2];
  *mpred++ = block[0][3];
  *mpred++ = block[1][3];
  *mpred++ = block[2][3];
  *mpred++ = block[3][3];
}

/*!
 ************************************************************************
 * \brief
 *    Predict one 4x4 Luma block
 ************************************************************************
 */
void
LumaPrediction4x4 (int  block_x,    // <--  relative horizontal block coordinate of 4x4 block
                   int  block_y,    // <--  relative vertical   block coordinate of 4x4 block
                   int  fw_mode,    // <--  forward  prediction mode (1-7, 0=DIRECT if bw_mode=0)
                   int  bw_mode,    // <--  backward prediction mode (1-7, 0=DIRECT if fw_mode=0)
                   int  fw_ref,      // <--  reference frame for forward prediction (-1: Intra4x4 pred. with fw_mode)
                                   int  bw_ref  )    
{
  static int fw_pred[16];
  static int bw_pred[16];

  int  i, j;
  int  block_x4  = block_x+4;
  int  block_y4  = block_y+4;
  int  pic_pix_x = img->pix_x + block_x;
  int  pic_pix_y = img->pix_y + block_y;
  int  bx        = block_x >> 2;
  int  by        = block_y >> 2;
  int* fpred     = fw_pred;
  int* bpred     = bw_pred;
	int  direct    = (fw_mode == 0 && bw_mode == 0 && (img->type == B_SLICE));
  //int  skipped   = (fw_mode == 0 && bw_mode == 0 && (img->type != B_SLICE && img->type != BS_IMG));
	int  skipped   = (fw_mode == 0 && bw_mode == 0 && (img->type != B_SLICE));
  int  *****fmv_array = img->all_mv;    // For MB level frame/field coding
  int  *****bmv_array = img->all_bmv;   // For MB level frame/field coding
  int apply_weights = ( (input->WeightedPrediction && (img->type == P_SLICE || img->type == SP_SLICE)) ||
                       (input->WeightedBiprediction && (img->type == B_SLICE)));                    
		//(input->WeightedBiprediction && (img->type == B_SLICE || img->type == BS_IMG)));
  int fw_ref_idx, bw_ref_idx;
	int block[BLOCK_SIZE][BLOCK_SIZE];

  if(input->InterlaceCodingOption >= MB_CODING && mb_adaptive && img->field_mode)
  {
    if(img->top_field)
    {
      pic_pix_y = img->field_pix_y + block_y;
      fmv_array = img->all_mv_top;
      bmv_array = img->all_bmv_top;
    }
    else
    {
      pic_pix_y = img->field_pix_y + block_y;
      fmv_array = img->all_mv_bot;
      bmv_array = img->all_bmv_bot;
    }

  }

  if (input->direct_type && direct)    
  {
    fw_ref= fwdir_refFrArr[pic_pix_y>>2][pic_pix_x>>2];
    bw_ref= bwdir_refFrArr[pic_pix_y>>2][pic_pix_x>>2];
  }

  if (img->type == B_SLICE && img->nal_reference_idc>0)
  {
    fw_ref_idx = fw_ref;
		bw_ref_idx = (bw_ref < 2) ? 1-bw_ref : bw_ref;
  }
  else
  {
    fw_ref_idx = fw_ref;
    bw_ref_idx = bw_ref;
  }

  direct_mode = direct && input->direct_type==0;


  if (fw_mode ||(direct && (!input->direct_type || fw_ref !=-1) )|| skipped)
  {
		get_block(fw_ref,listX[0], (pic_pix_x << 2)+fmv_array [bx][by][fw_ref][fw_mode][0],(pic_pix_y << 2) + fmv_array [bx][by][fw_ref][fw_mode][1],block);
    copyblock4x4(fw_pred,block);
  }

  if (bw_mode || (direct && (!input->direct_type || bw_ref !=-1) ))
  { 
	  if (input->InterlaceCodingOption == 0)
      get_block(bw_ref,listX[1], (pic_pix_x << 2)+bmv_array [bx][by][bw_ref][bw_mode][0],(pic_pix_y << 2) + bmv_array [bx][by][bw_ref][bw_mode][1],block);
    else
    {	
      if (bw_ref<0) bw_ref = 0;
      get_block(bw_ref,listX[1], (pic_pix_x << 2)+bmv_array [bx][by][bw_ref][bw_mode][0],(pic_pix_y << 2) + bmv_array [bx][by][bw_ref][bw_mode][1],block);
    }
    copyblock4x4(bw_pred,block);
  }

  if (apply_weights)
  {
    if (direct || (fw_mode && bw_mode))
    {
      if (input->direct_type && direct)
      {
        for   (j=block_y; j<block_y4; j++)
          for (i=block_x; i<block_x4; i++)  
            if (fw_ref ==-1)
              img->mpr[i][j] = clip1a(((wp_weight[1][bw_ref_idx][0] * *bpred++ + wp_luma_round) >> luma_log_weight_denom) + wp_offset[1][bw_ref_idx][0]);
            else if (bw_ref ==-1 )
              img->mpr[i][j] = clip1a(((wp_weight[0][fw_ref_idx][0] * *fpred++ + wp_luma_round) >> luma_log_weight_denom) + wp_offset[0][fw_ref_idx][0] );
            else 
              img->mpr[i][j] = clip1a(((wbp_weight[0][fw_ref_idx][bw_ref_idx][0] * *fpred++ + wbp_weight[1][fw_ref_idx][bw_ref_idx][0] * *bpred++ + 2*wp_luma_round) >> (luma_log_weight_denom + 1)) + ((wp_offset[0][fw_ref_idx][0] + wp_offset[1][bw_ref_idx][0] + 1)>>1)); 
      }
      else
        for   (j=block_y; j<block_y4; j++)
          for (i=block_x; i<block_x4; i++)  
              img->mpr[i][j] = clip1a(((wbp_weight[0][fw_ref_idx][bw_ref_idx][0] * *fpred++ + wbp_weight[1][fw_ref_idx][bw_ref_idx][0] * *bpred++ + 2*wp_luma_round) >> (luma_log_weight_denom + 1)) + ((wp_offset[0][fw_ref_idx][0] + wp_offset[1][bw_ref_idx][0] + 1)>>1)); 
    }
		else if (img->type == B_SLICE && img->nal_reference_idc>0)
    {
      for   (j=block_y; j<block_y4; j++)
        for (i=block_x; i<block_x4; i++)  
          img->mpr[i][j] = clip1a((wbp_weight[0][fw_ref_idx][bw_ref_idx][0] * *fpred++ + wbp_weight[1][fw_ref_idx][bw_ref_idx][0] * *bpred++ + wp_offset[0][fw_ref_idx][0] + wp_offset[1][bw_ref_idx][0] + 2*wp_luma_round) >> (luma_log_weight_denom + 1));
    }
    else if (fw_mode || skipped)
    {
      for   (j=block_y; j<block_y4; j++)
        for (i=block_x; i<block_x4; i++)  
          img->mpr[i][j] = clip1a((wp_weight[0][fw_ref_idx][0] * *fpred++ + wp_offset[0][fw_ref_idx][0] + wp_luma_round) >> luma_log_weight_denom);
    }
    else
    {
      for   (j=block_y; j<block_y4; j++)
        for (i=block_x; i<block_x4; i++)  
          img->mpr[i][j] = clip1a((wp_weight[1][bw_ref_idx][0] * *bpred++ + wp_offset[1][bw_ref_idx][0] + wp_luma_round) >> luma_log_weight_denom);
    }
  }
  else
  {
    if (direct || (fw_mode && bw_mode))
    {
      if (input->direct_type && direct)
      {
        for   (j=block_y; j<block_y4; j++)
          for (i=block_x; i<block_x4; i++)  
            if (fw_ref ==-1)
              img->mpr[i][j] = *bpred++;
            else if (bw_ref ==-1 )
              img->mpr[i][j] = *fpred++;
            else 
              img->mpr[i][j] = (*fpred++ + *bpred++ + 1) / 2; 
      }
      else
        for   (j=block_y; j<block_y4; j++)
          for (i=block_x; i<block_x4; i++)  
            img->mpr[i][j] = (*fpred++ + *bpred++ + 1) / 2; 
    }
		else if (img->type == B_SLICE && img->nal_reference_idc>0)
    {
      for   (j=block_y; j<block_y4; j++)
        for (i=block_x; i<block_x4; i++)  
          img->mpr[i][j] = (*fpred++ + *bpred++ + 1) / 2; 
    }
    else if (fw_mode || skipped)
    {
      for   (j=block_y; j<block_y4; j++)
        for (i=block_x; i<block_x4; i++)  img->mpr[i][j] = *fpred++;
    }
    else
    {
      for   (j=block_y; j<block_y4; j++)
        for (i=block_x; i<block_x4; i++)  img->mpr[i][j] = *bpred++;
    }
  }
}

/*!
 ************************************************************************
 * \brief
 *    Residual Coding of an 8x8 Luma block (not for intra)
 ************************************************************************
 */
int                                       //  ==> coefficient cost
LumaResidualCoding8x8 (int  *cbp,         //  --> cbp (updated according to processed 8x8 luminance block)
                       int  *cbp_blk,     //  --> block cbp (updated according to processed 8x8 luminance block)
                       int  block8x8,     // <--  block number of 8x8 block
                       int  fw_mode,      // <--  forward  prediction mode (1-7, 0=DIRECT)
                       int  bw_mode,      // <--  backward prediction mode (1-7, 0=DIRECT)
                       int  fw_refframe,  // <--  reference frame for forward prediction
                       int  bw_refframe   // <--  reference frame for backward prediction
                       )
{
  int    block_y, block_x, pic_pix_y, pic_pix_x, i, j, nonzero, cbp_blk_mask;
  int    coeff_cost = 0;
  int    mb_y       = (block8x8 / 2) << 3;
  int    mb_x       = (block8x8 % 2) << 3;
  int    cbp_mask   = 1 << block8x8;
  int    bxx, byy;                   // indexing curr_blk
  int    scrFlag = 0;                // 0=noSCR, 1=strongSCR, 2=jmSCR
  byte** imgY_original = imgY_org;
  int  pix_y    = img->pix_y;
  int    skipped    = (fw_mode == 0 && bw_mode == 0 && (img->type != B_SLICE));

  if (img->type==B_SLICE)
    scrFlag = 1;
  
  if(input->InterlaceCodingOption >= MB_CODING && mb_adaptive && img->field_mode)
  {
    pix_y     = img->field_pix_y;
    imgY_original = img->top_field ? imgY_org_top:imgY_org_bot;
  }


  //===== loop over 4x4 blocks =====
  for (byy=0, block_y=mb_y; block_y<mb_y+8; byy+=4, block_y+=4)
  {
    pic_pix_y = pix_y + block_y;

    for (bxx=0, block_x=mb_x; block_x<mb_x+8; bxx+=4, block_x+=4)
    {
      pic_pix_x = img->pix_x + block_x;

      cbp_blk_mask = (block_x>>2) + block_y;

      //===== prediction of 4x4 block =====
      LumaPrediction4x4 (block_x, block_y, fw_mode, bw_mode, fw_refframe, bw_refframe);

      //===== get displaced frame difference ======                
      for (j=0; j<4; j++)
      for (i=0; i<4; i++)
      {
        img->m7[i][j] = imgY_original[pic_pix_y+j][pic_pix_x+i] - img->mpr[i+block_x][j+block_y];
      }

      //===== DCT, Quantization, inverse Quantization, IDCT, Reconstruction =====      
      if (img->NoResidueDirect != 1 && !skipped  )
      {
        //===== DCT, Quantization, inverse Quantization, IDCT, Reconstruction =====
        if (img->type!=SP_SLICE)  nonzero = dct_luma   (block_x, block_y, &coeff_cost, 0);
        else                      nonzero = dct_luma_sp(block_x, block_y, &coeff_cost);
        if (nonzero)
        {
          (*cbp_blk) |= 1 << cbp_blk_mask;  // one bit for every 4x4 block
          (*cbp)     |= cbp_mask;           // one bit for the 4x4 blocks of an 8x8 block
        }
      }
    }
  }

  /*
  The purpose of the action below is to prevent that single or 'expensive' coefficients are coded.
  With 4x4 transform there is larger chance that a single coefficient in a 8x8 or 16x16 block may be nonzero.
  A single small (level=1) coefficient in a 8x8 block will cost: 3 or more bits for the coefficient,
  4 bits for EOBs for the 4x4 blocks,possibly also more bits for CBP.  Hence the total 'cost' of that single
  coefficient will typically be 10-12 bits which in a RD consideration is too much to justify the distortion improvement.
  The action below is to watch such 'single' coefficients and set the reconstructed block equal to the prediction according
  to a given criterium.  The action is taken only for inter luma blocks.

  Notice that this is a pure encoder issue and hence does not have any implication on the standard.
  coeff_cost is a parameter set in dct_luma() and accumulated for each 8x8 block.  If level=1 for a coefficient,
  coeff_cost is increased by a number depending on RUN for that coefficient.The numbers are (see also dct_luma()): 3,2,2,1,1,1,0,0,...
  when RUN equals 0,1,2,3,4,5,6, etc.
  If level >1 coeff_cost is increased by 9 (or any number above 3). The threshold is set to 3. This means for example:
  1: If there is one coefficient with (RUN,level)=(0,1) in a 8x8 block this coefficient is discarded.
  2: If there are two coefficients with (RUN,level)=(1,1) and (4,1) the coefficients are also discarded
  sum_cnt_nonz is the accumulation of coeff_cost over a whole macro block.  If sum_cnt_nonz is 5 or less for the whole MB,
  all nonzero coefficients are discarded for the MB and the reconstructed block is set equal to the prediction.
  */

  if (img->NoResidueDirect != 1 && !skipped && coeff_cost <= _LUMA_COEFF_COST_)
  {
    coeff_cost  = 0;
    (*cbp)     &=  (63 - cbp_mask);
    (*cbp_blk) &= ~(51 << (4*block8x8-2*(block8x8%2)));

    for (i=mb_x; i<mb_x+8; i++)
    for (j=mb_y; j<mb_y+8; j++)
    {
      enc_picture->imgY[img->pix_y+j][img->pix_x+i] = img->mpr[i][j];
    }
    if (img->type==SP_SLICE)
    {
      for (i=mb_x; i < mb_x+BLOCK_SIZE*2; i+=BLOCK_SIZE)
        for (j=mb_y; j < mb_y+BLOCK_SIZE*2; j+=BLOCK_SIZE)
          copyblock_sp(i,j);
    }
  }

  return coeff_cost;
}


/*!
 ************************************************************************
 * \brief
 *    Set mode parameters and reference frames for an 8x8 block
 ************************************************************************
 */
void
SetModesAndRefframe (int b8, int* fw_mode, int* bw_mode, int* fw_ref, int* bw_ref)
{
  Macroblock* currMB = &img->mb_data[img->current_mb_nr];
  int         j      = 2*(b8/2);
  int         i      = 2*(b8%2);
  int**     frefarr = refFrArr;   // For MB level field/frame coding
  int**     fw_refarr = fw_refFrArr;  // For MB level field/frame coding
  int**     bw_refarr = bw_refFrArr;  // For MB level field/frame coding
  int     block_y = img->block_y; // For MB level field/frame coding


  if(input->InterlaceCodingOption >= MB_CODING && mb_adaptive && img->field_mode)
  {
    block_y = img->field_block_y;
    if(img->top_field)
    {
      frefarr = refFrArr_top;
      fw_refarr = fw_refFrArr_top;
      bw_refarr = bw_refFrArr_top;
    }
    else
    {
      frefarr = refFrArr_bot;
      fw_refarr = fw_refFrArr_bot;
      bw_refarr = bw_refFrArr_bot;
    }
  }

  *fw_mode = *bw_mode = *fw_ref = *bw_ref = -1;

	if (img->type!=B_SLICE)
  {
    *fw_ref = frefarr[block_y+j][img->block_x+i];
    *bw_ref = 0;
    *bw_mode  = 0;
    *fw_mode  = currMB->b8mode[b8];
  }
  else
  {
    if (currMB->b8pdir[b8]==-1)
    {
      *fw_ref   = -1;
      *bw_ref   = -1;
      *fw_mode  =  0;
      *bw_mode  =  0;
    }
    else if (currMB->b8pdir[b8]==0)
    {
      *fw_ref   = fw_refarr[block_y+j][img->block_x+i];
      *bw_ref   = 0;
      *fw_mode  = currMB->b8mode[b8];
      *bw_mode  = 0;
    }
    else if (currMB->b8pdir[b8]==1)
    {
      *fw_ref   = 0;
      *bw_ref   = bw_refarr[block_y+j][img->block_x+i];
      *fw_mode  = 0;
      *bw_mode  = currMB->b8mode[b8];
    }
    else
    {
      *fw_ref   = fw_refarr[block_y+j][img->block_x+i];
      *bw_ref   = bw_refarr[block_y+j][img->block_x+i];
      *fw_mode  = currMB->b8mode[b8];
      *bw_mode  = currMB->b8mode[b8];

      if (currMB->b8mode[b8]==0) // direct
      {
				if (img->type==B_SLICE && img->nal_reference_idc>0)
        {
            *fw_ref = 0;
            *bw_ref = 0;
        }
        else if (img->type==B_SLICE)
        {
          *fw_ref = max(0,frefarr[block_y+j][img->block_x+i]);
          *bw_ref = 0;
        }
        else
        {
          *fw_ref = max(0,frefarr[block_y+j][img->block_x+i]);
          *bw_ref = 0;
        }
      }
    }
  }
}


/*!
 ************************************************************************
 * \brief
 *    Residual Coding of a Luma macroblock (not for intra)
 ************************************************************************
 */
void
LumaResidualCoding ()
{
  int i,j,block8x8,b8_x,b8_y;
  int fw_mode, bw_mode, refframe;
  int sum_cnt_nonz;
  Macroblock *currMB = &img->mb_data[img->current_mb_nr];

  currMB->cbp     = 0 ;
  currMB->cbp_blk = 0 ;
  sum_cnt_nonz    = 0 ;

  for (block8x8=0; block8x8<4; block8x8++)
  {
    int bw_ref;
    SetModesAndRefframe (block8x8, &fw_mode, &bw_mode, &refframe, &bw_ref);

    sum_cnt_nonz += LumaResidualCoding8x8 (&(currMB->cbp), &(currMB->cbp_blk), block8x8,
                                           fw_mode, bw_mode, refframe, bw_ref);
  }

  if (sum_cnt_nonz <= 5 )
  {
     currMB->cbp     &= 0xfffff0 ;
     currMB->cbp_blk &= 0xff0000 ;
     for (i=0; i < MB_BLOCK_SIZE; i++)
     {
       for (j=0; j < MB_BLOCK_SIZE; j++)
       {
         enc_picture->imgY[img->pix_y+j][img->pix_x+i]=img->mpr[i][j];
       }
     }
     if (img->type==SP_SLICE)
     {
       for(block8x8=0;block8x8<4;block8x8++)
       {
         b8_x=(block8x8&1)<<3;
         b8_y=(block8x8&2)<<2;
         for (i=0;i<8;i+=4)
           for (j=0;j<8;j+=4)
             copyblock_sp(b8_x+i,b8_y+j);
       }
     }
   }
}



/*!
 ************************************************************************
 * \brief
 *    Predict one component of a chroma 4x4 block
 ************************************************************************
 */
void
OneComponentChromaPrediction4x4 (int*     mpred,      //  --> array to store prediction values
                                 int      pix_c_x,    // <--  horizontal pixel coordinate of 4x4 block
                                 int      pix_c_y,    // <--  vertical   pixel coordinate of 4x4 block
                                 int***** mv,         // <--  motion vector array
                                 int      ref,        // <--  reference frame parameter (0.../ -1: backward)
                                 int      blocktype,  // <--  block type
                                 int      uv)         // <--  chroma component
{
  int     i, j, ii, jj, ii0, jj0, ii1, jj1, if0, if1, jf0, jf1;
  int     incr;
  int*    mvb;
  int     refframe  = (ref<0 ?      0 :    ref);
  pel_t** refimage;
  int     je        = pix_c_y + 4;
  int     ie        = pix_c_x + 4;
  int     f1        = 8 , f2=f1-1, f3=f1*f1, f4=f3>>1;
  int     s1        = 3;
  int     img_pic_c_y = img->pix_c_y;
  int   scale   = 1;
  int   pred_dir = (ref<0 ? 1 : 0);
  int   ref_idx;
  int  list_offset;
  StorablePicture **list;

  int curr_mb_field = ((img->MbaffFrameFlag)&&(img->mb_data[img->current_mb_nr].mb_field));

  // find out the correct list offsets
  if (curr_mb_field)
  {
    if(img->current_mb_nr%2)
      list_offset = 4; // top field mb
    else
      list_offset = 2; // bottom field mb
    //max_y_cr = img->height_cr/2-1;
  }
  else
  {
    list_offset = 0;  // no mb aff or frame mb
    //max_y_cr = img->height_cr-1;
  }


  list     = listX[0+list_offset+ pred_dir];

  incr      = (ref==-1 ? (!img->fld_type&&enc_picture!=enc_frame_picture): direct_mode ? (!img->fld_type&&enc_picture!=enc_frame_picture) : (enc_picture!=enc_frame_picture)) ;
  
  if(input->InterlaceCodingOption >= MB_CODING && mb_adaptive && img->field_mode)
  {
    incr  = (ref == -1 ? (img->top_field&&img->field_mode):(direct_mode ? (img->top_field&&img->field_mode):(img->field_mode)));
    scale = 2;
    img_pic_c_y = img->field_pix_c_y;
  }
  
  //refimage  = img->type==B_SLICE? mcef [ref+1+incr][uv] : mcef [ref][uv];
	refimage  = (ref<0 ? listX[LIST_1][0]->imgUV[uv]: listX[LIST_0][ref]->imgUV[uv]);

  for (j=pix_c_y; j<je; j++)
  for (i=pix_c_x; i<ie; i++)
  {
    ref_idx   = enc_picture->ref_idx[LIST_0+pred_dir][pix_c_x>>1][pix_c_y>>1];  
    
    mvb  = mv [(i-img->pix_c_x)>>1][(j-img_pic_c_y)>>1][refframe][blocktype];
    ii   = (i<<s1) + mvb[0];
    jj   = (j<<s1) + mvb[1];
    
    //SW
    /*
    if (!((img->MbaffFrameFlag)&&(img->mb_data[img->current_mb_nr].mb_field)))
      j1=(img->pix_c_y+jj+joff)*f1+mv_array[if1][jf][1];
    else
    {
      if (img->current_mb_nr%2 == 0) 
        j1=(img->pix_c_y)/2*f1 + (jj+joff)*f1+mv_array[if1][jf][1];
      else
        j1=(img->pix_c_y-8)/2*f1 + (jj+joff)*f1 +mv_array[if1][jf][1];
    }
    */
		if (ref_idx!=-1)
      jj += list[ref_idx]->chroma_vector_adjustment;
		else
		{
			if (img->type!= I_SLICE)
			jj += list[0]->chroma_vector_adjustment;
		}


    ii0  = max (0, min (img->width_cr -1, ii>>s1     ));
    jj0  = max (0, min (img->height_cr/scale-1, jj>>s1     ));    // For MB level field/frame -- scale chroma height by 2
    ii1  = max (0, min (img->width_cr -1, (ii+f2)>>s1));
    jj1  = max (0, min (img->height_cr/scale-1, (jj+f2)>>s1));

    if1  = (ii&f2);  if0 = f1-if1;
    jf1  = (jj&f2);  jf0 = f1-jf1;

    *mpred++ = (if0 * jf0 * refimage[jj0][ii0] +
                if1 * jf0 * refimage[jj0][ii1] +
                if0 * jf1 * refimage[jj1][ii0] +
                if1 * jf1 * refimage[jj1][ii1] + f4) / f3;
  }
}


/*!
 ************************************************************************
 * \brief
 *    Predict an intra chroma 4x4 block
 ************************************************************************
 */
void IntraChromaPrediction4x4 (int  uv,       // <-- colour component
                               int  block_x,  // <-- relative horizontal block coordinate of 4x4 block
                               int  block_y)  // <-- relative vertical   block coordinate of 4x4 block
{
  int mode = img->mb_data[img->current_mb_nr].c_ipred_mode;
  int i, j;

  //===== prediction =====
  for (j=block_y; j<block_y+4; j++)
  for (i=block_x; i<block_x+4; i++)
  {
    img->mpr[i][j] = img->mprr_c[uv][mode][i][j];
  }
}



/*!
 ************************************************************************
 * \brief
 *    Predict one chroma 4x4 block
 ************************************************************************
 */
void
ChromaPrediction4x4 (int  uv,           // <-- colour component
                     int  block_x,      // <-- relative horizontal block coordinate of 4x4 block
                     int  block_y,      // <-- relative vertical   block coordinate of 4x4 block
                     int  fw_mode,      // <-- forward  prediction mode (1-7, 0=DIRECT if bw_mode=0)
                     int  bw_mode,      // <-- backward prediction mode (1-7, 0=DIRECT if fw_mode=0)
                     int  fw_ref_frame, // <-- reference frame for forward prediction (if (<0) -> intra prediction)
                     int  bw_ref_frame) // <-- reference frame for backward prediction 
{
  static int fw_pred[16];
  static int bw_pred[16];

  int  i, j;
  int  block_x4  = block_x+4;
  int  block_y4  = block_y+4;
  int  pic_pix_x = img->pix_c_x + block_x;
  int  pic_pix_y = img->pix_c_y + block_y;
  int* fpred     = fw_pred;
  int* bpred     = bw_pred;
 	int  direct    = (fw_mode == 0 && bw_mode == 0 && (img->type == B_SLICE));

  //int  skipped   = (fw_mode == 0 && bw_mode == 0 && (img->type != B_SLICE && img->type != BS_IMG));
	int  skipped   = (fw_mode == 0 && bw_mode == 0 && (img->type != B_SLICE));
  int***** fmv_array = img->all_mv;
  int***** bmv_array = img->all_bmv;
  int fw_ref_idx, bw_ref_idx;
	int apply_weights = ( (input->WeightedPrediction && (img->type == P_SLICE||img->type == SP_SLICE)) ||
		                 (input->WeightedBiprediction && (img->type == B_SLICE)));
//                 (input->WeightedBiprediction && (img->type == B_SLICE || img->type == BS_IMG)));


  if(input->InterlaceCodingOption >= MB_CODING && mb_adaptive && img->field_mode)
  {
    pic_pix_y   = img->field_pix_c_y + block_y;
    if(img->top_field)
    {
      fmv_array = img->all_mv_top;
      bmv_array = img->all_bmv_top;
    }
    else
    {
      fmv_array = img->all_mv_bot;
      bmv_array = img->all_bmv_bot;
    }
  }

  direct_mode = direct && input->direct_type==0;

  //===== INTRA PREDICTION =====
  if (fw_ref_frame < 0)
  {
    IntraChromaPrediction4x4 (uv, block_x, block_y);
    return;
  }
  
  if (input->direct_type && direct)
  {    
    fw_ref_frame = fwdir_refFrArr[pic_pix_y/2][pic_pix_x/2];
    bw_ref_frame = bwdir_refFrArr[pic_pix_y/2][pic_pix_x/2];
  }
//  if (img->type == BS_IMG)
	if (img->type == B_SLICE && img->nal_reference_idc>0)
  {
    fw_ref_idx = fw_ref_frame;
    bw_ref_idx = (bw_ref_frame < 2) ? 1-bw_ref_frame : bw_ref_frame;
  }
  else
  {
    fw_ref_idx = fw_ref_frame;
    bw_ref_idx = bw_ref_frame;
  }


  //===== INTER PREDICTION =====
  if (fw_mode || (direct && (!input->direct_type || fw_ref_frame!=-1)) || skipped)
  {
    OneComponentChromaPrediction4x4 (fw_pred, pic_pix_x, pic_pix_y, fmv_array , fw_ref_frame, fw_mode, uv);
  }
  if (bw_mode || (direct && (!input->direct_type || bw_ref_frame!=-1)))
  {
    //if (img->type == BS_IMG)
		if (img->type == B_SLICE && img->nal_reference_idc>0)
      OneComponentChromaPrediction4x4 (bw_pred, pic_pix_x, pic_pix_y, bmv_array, bw_ref_frame, bw_mode, uv);
    else
      OneComponentChromaPrediction4x4 (bw_pred, pic_pix_x, pic_pix_y, bmv_array,           -1, bw_mode, uv);
  }

  if (apply_weights)
  {
    if (direct || (fw_mode && bw_mode))
    {
      if (direct && input->direct_type)
      {
        for (j=block_y; j<block_y4; j++)
        for (i=block_x; i<block_x4; i++)  
          if (fw_ref_frame==-1)
            img->mpr[i][j] = clip1a(((wp_weight[1][bw_ref_idx][uv+1] * *bpred++  + wp_chroma_round) >> chroma_log_weight_denom) + wp_offset[1][bw_ref_idx][uv+1]);
          else if (bw_ref_frame==-1)
            img->mpr[i][j] =  clip1a(((wp_weight[0][fw_ref_idx][uv+1] * *fpred++  + wp_chroma_round) >> chroma_log_weight_denom) + wp_offset[0][fw_ref_idx][uv+1]);
          else
            img->mpr[i][j] =  clip1a(((wbp_weight[0][fw_ref_idx][bw_ref_idx][uv+1] * *fpred++ + wbp_weight[1][fw_ref_idx][bw_ref_idx][uv+1] * *bpred++ 
                  + 2*wp_chroma_round) >> (chroma_log_weight_denom + 1)) + ((wp_offset[0][fw_ref_idx][uv+1] + wp_offset[1][bw_ref_idx][uv+1] + 1)>>1) );
      }
      else
        for (j=block_y; j<block_y4; j++)
          for (i=block_x; i<block_x4; i++)  
            img->mpr[i][j] = clip1a(((wbp_weight[0][fw_ref_idx][bw_ref_idx][uv+1] * *fpred++ + wbp_weight[1][fw_ref_idx][bw_ref_idx][uv+1] * *bpred++ 
                    + 2*wp_chroma_round) >> (chroma_log_weight_denom + 1)) + ((wp_offset[0][fw_ref_idx][uv+1] + wp_offset[1][bw_ref_idx][uv+1] + 1)>>1));
 
    }
		else if (img->type == B_SLICE && img->nal_reference_idc>0)
    {
      for (j=block_y; j<block_y4; j++)
        for (i=block_x; i<block_x4; i++)  
          img->mpr[i][j] = clip1a(((wbp_weight[0][fw_ref_idx][bw_ref_idx][uv+1] * *fpred++ + wbp_weight[1][fw_ref_idx][bw_ref_idx][uv+1] * *bpred++ 
                     + 2*wp_chroma_round) >> (chroma_log_weight_denom + 1)) + ((wp_offset[0][fw_ref_idx][uv+1] + wp_offset[1][bw_ref_idx][uv+1] + 1)>>1));
    }
    else if (fw_mode || skipped)
    {
      for (j=block_y; j<block_y4; j++)
      for (i=block_x; i<block_x4; i++)  
           img->mpr[i][j] = clip1a(((wp_weight[0][fw_ref_idx][uv+1] * *fpred++ + wp_chroma_round) >> chroma_log_weight_denom) +  wp_offset[0][fw_ref_idx][uv+1]);
    }
    else
    {
      for (j=block_y; j<block_y4; j++)
      for (i=block_x; i<block_x4; i++)  
            img->mpr[i][j] = clip1a(((wp_weight[1][bw_ref_idx][uv+1] * *bpred++ + wp_chroma_round) >> chroma_log_weight_denom) + wp_offset[1][bw_ref_idx][uv+1]);

    }       
  }
  else
  {
    if (direct || (fw_mode && bw_mode))
    {
      if (direct && input->direct_type)
      {
        for (j=block_y; j<block_y4; j++)
          for (i=block_x; i<block_x4; i++)  
            if (fw_ref_frame==-1)
              img->mpr[i][j] = *bpred++;
            else if (bw_ref_frame==-1)
              img->mpr[i][j] = *fpred++;
            else
              img->mpr[i][j] = (*fpred++ + *bpred++ + 1) / 2; 
      }
      else
        for (j=block_y; j<block_y4; j++)
          for (i=block_x; i<block_x4; i++)  
            img->mpr[i][j] = (*fpred++ + *bpred++ + 1) / 2; 
     }
		 else if (img->type == B_SLICE && img->nal_reference_idc>0)
     {
       for (j=block_y; j<block_y4; j++)
         for (i=block_x; i<block_x4; i++)  
           img->mpr[i][j] = (*fpred++ + *bpred++ + 1) / 2; 
     }
     else if (fw_mode || skipped)
     {
       for (j=block_y; j<block_y4; j++)
         for (i=block_x; i<block_x4; i++)  img->mpr[i][j] = *fpred++;
     }
     else
     {
       for (j=block_y; j<block_y4; j++)
         for (i=block_x; i<block_x4; i++)  img->mpr[i][j] = *bpred++;
     }
  }
}




/*!
 ************************************************************************
 * \brief
 *    Chroma residual coding for an macroblock
 ************************************************************************
 */
void ChromaResidualCoding (int* cr_cbp)
{
  int   uv, block8, block_y, block_x, j, i;
  int   fw_mode, bw_mode, refframe;
  int   skipped = (img->mb_data[img->current_mb_nr].mb_type == 0 && (img->type == P_SLICE || img->type == SP_SLICE));
  int   incr = 1, offset = 0; // For MB level field/frame coding 
  int   bw_ref;


  if(input->InterlaceCodingOption >= MB_CODING && mb_adaptive && img->field_mode)
  {
    incr = 2;   // increment every other field  
    if(!img->top_field)
      offset = -7;
  }

  for (*cr_cbp=0, uv=0; uv<2; uv++)
  {
    //===== prediction of chrominance blocks ===d==
    block8 = 0;
    for (block_y=0; block_y<8; block_y+=4)
    for (block_x=0; block_x<8; block_x+=4, block8++)
    {
      SetModesAndRefframe (block8, &fw_mode, &bw_mode, &refframe, &bw_ref);

      ChromaPrediction4x4 (uv, block_x, block_y, fw_mode, bw_mode, refframe, bw_ref);
    }

        // ==== set chroma residue to zero for skip Mode in SP frames 
    if (img->NoResidueDirect)
      for (j=0; j<8; j++)
        for (i=0; i<8; i++)
        {
          enc_picture->imgUV[uv][img->pix_c_y+j][img->pix_c_x+i] = img->mpr[i][j];
        }
        else
    if (skipped && img->type==SP_SLICE)
      for (j=0; j<8; j++)
        for (i=0; i<8; i++)
        {
          img->m7[i][j] = 0;
        }
    else
    if (skipped)
    {
      for (j=0; j<8; j++)
        for (i=0; i<8; i++)
        {
          enc_picture->imgUV[uv][img->pix_c_y+j][img->pix_c_x+i] = img->mpr[i][j];
        }
    }
    else
      for (j=0; j<8; j++)
        for (i=0; i<8; i++)
        {
          img->m7[i][j] = imgUV_org[uv][img->pix_c_y+(j*incr)+offset][img->pix_c_x+i] - img->mpr[i][j];
        }

    //===== DCT, Quantization, inverse Quantization, IDCT, and Reconstruction =====
    //===== Call function for skip mode in SP frames to properly process frame ====
    
    if (skipped && img->type==SP_SLICE)
    {
        *cr_cbp=dct_chroma_sp(uv,*cr_cbp);
    }
    else
    if (!img->NoResidueDirect && !skipped)
    {
      if (img->type!=SP_SLICE || IS_INTRA (&img->mb_data[img->current_mb_nr]))
        *cr_cbp=dct_chroma   (uv,*cr_cbp);
      else
        *cr_cbp=dct_chroma_sp(uv,*cr_cbp);
    }
  }

  //===== update currMB->cbp =====
  img->mb_data[img->current_mb_nr].cbp += ((*cr_cbp)<<4);  
}


/*!
 ************************************************************************
 * \brief
 *    Predict an intra chroma 8x8 block
 ************************************************************************
 */
void IntraChromaPrediction8x8 (int *mb_up, int *mb_left, int*mb_up_left)
{
  Macroblock *currMB = &img->mb_data[img->current_mb_nr];
  int     s, s0, s1, s2, s3, i, j, k;
  pel_t** image;
  int     block_x, block_y, b4;
  int     img_cx            = img->pix_c_x;
  int     img_cy            = img->pix_c_y;
  int     img_cx_1          = img->pix_c_x-1;
  int     img_cx_4          = img->pix_c_x+4;
  int     img_cy_1          = img->pix_c_y-1;
  int     img_cy_4          = img->pix_c_y+4;
  int     mb_nr             = img->current_mb_nr;
  int     mb_width          = img->width/16;
  int     mb_available_up   = (img_cy/BLOCK_SIZE == 0) ? 0 : (img->mb_data[mb_nr].slice_nr==img->mb_data[mb_nr-mb_width].slice_nr);
  int     mb_available_left = (img_cx/BLOCK_SIZE == 0) ? 0 : (img->mb_data[mb_nr].slice_nr==img->mb_data[mb_nr-1]       .slice_nr);
  int     mb_available_up_left = (img_cx/BLOCK_SIZE == 0 || img_cy/BLOCK_SIZE == 0) ? 0 : (img->mb_data[mb_nr].slice_nr==img->mb_data[mb_nr-mb_width-1].slice_nr);
  int     ih,iv;
  int     ib,ic,iaa;
  int     uv;
  int     hline[8], vline[8];
  int     mode;
  int     best_mode = DC_PRED_8;         //just an initilaization here, should always be overwritten
  int     cost;
  int     min_cost;
  int     diff[16];
	int     left_avail;
	PixelPos up;       //!< pixel position p(0,-1)
  PixelPos left[9];  //!< pixel positions p(-1, -1..8)


	for (i=0;i<9;i++)
  {
    getNeighbour(mb_nr, -1 ,  i-1 , 0, &left[i]);
  }
  
  getNeighbour(mb_nr, 0     ,  -1 , 0, &up);

  if(input->InterlaceCodingOption >= MB_CODING && mb_adaptive && img->field_mode)
  {
    img_cy   = img->field_pix_c_y;
    img_cy_1 = img->field_pix_c_y-1;
    img_cy_4 = img->field_pix_c_y+4;
    mb_available_up = (img_cy/BLOCK_SIZE == 0) ? 0 : (img->mb_data[mb_nr].slice_nr==img->mb_data[mb_nr-mb_width].slice_nr);
    mb_available_up_left = (img_cx/BLOCK_SIZE == 0 || img_cy/BLOCK_SIZE == 0) ? 0 : (img->mb_data[mb_nr].slice_nr==img->mb_data[mb_nr-mb_width-1].slice_nr);
  }

  if(input->UseConstrainedIntraPred)
  {
		mb_available_up      = up.available ? img->intra_block[up.mb_addr] : 0;
    for (i=1, left_avail=1; i<9;i++)
      mb_available_left  &= left[i].available ? img->intra_block[left[i].mb_addr]: 0;
    mb_available_up_left = left[0].available ? img->intra_block[left[0].mb_addr]: 0;
  }

  if (mb_up)
    *mb_up = mb_available_up;
  if (mb_left)
    *mb_left = mb_available_left;
  if( mb_up_left )
    *mb_up_left = mb_available_up_left;

  // compute all chroma intra prediction modes for both U and V
  for (uv=0; uv<2; uv++)
  {
/*    if(input->InterlaceCodingOption >= MB_CODING && mb_adaptive && img->field_mode)
    {
      if(img->top_field)
        image = imgUV_top[uv];
      else
        image = imgUV_bot[uv];
    }
    else
*/
    image = enc_picture->imgUV[uv];

    // DC prediction
    for (block_y=0; block_y<8; block_y+=4)
    for (block_x=0; block_x<8; block_x+=4)
    {
      s=128;
      s0=s1=s2=s3=0;
      //===== get prediction value =====
      switch ((block_y>>1) + (block_x>>2))
      {
      case 0:  //===== TOP LEFT =====
        if      (mb_available_up)    for (i=0;i<4;i++)  s0 += image[img_cy_1  ][img_cx  +i];
        if      (mb_available_left)  for (i=0;i<4;i++)  s2 += image[img_cy  +i][img_cx_1  ];
        if      (mb_available_up && mb_available_left)  s  = (s0+s2+4) >> 3;
        else if (mb_available_up)                       s  = (s0   +2) >> 2;
        else if (mb_available_left)                     s  = (s2   +2) >> 2;
        break;
      case 1: //===== TOP RIGHT =====
        if      (mb_available_up)    for (i=0;i<4;i++)  s1 += image[img_cy_1  ][img_cx_4+i];
        else if (mb_available_left)  for (i=0;i<4;i++)  s2 += image[img_cy  +i][img_cx_1  ];
        if      (mb_available_up)                       s  = (s1   +2) >> 2;
        else if (mb_available_left)                     s  = (s2   +2) >> 2;
        break;
      case 2: //===== BOTTOM LEFT =====
        if      (mb_available_left)  for (i=0;i<4;i++)  s3 += image[img_cy_4+i][img_cx_1  ];
        else if (mb_available_up)    for (i=0;i<4;i++)  s0 += image[img_cy_1  ][img_cx  +i];
        if      (mb_available_left)                     s  = (s3   +2) >> 2;
        else if (mb_available_up)                       s  = (s0   +2) >> 2;
        break;
      case 3: //===== BOTTOM RIGHT =====
        if      (mb_available_up)    for (i=0;i<4;i++)  s1 += image[img_cy_1  ][img_cx_4+i];
        if      (mb_available_left)  for (i=0;i<4;i++)  s3 += image[img_cy_4+i][img_cx_1  ];
        if      (mb_available_up && mb_available_left)  s  = (s1+s3+4) >> 3;
        else if (mb_available_up)                       s  = (s1   +2) >> 2;
        else if (mb_available_left)                     s  = (s3   +2) >> 2;
        break;
      }

      //===== prediction =====
      for (j=block_y; j<block_y+4; j++)
      for (i=block_x; i<block_x+4; i++)
      {
        img->mprr_c[uv][DC_PRED_8][i][j] = s;
      }
    }

    // vertical prediction
    if (mb_available_up)
    {
      for (i=0; i<8; i++)
        hline[i] = image[img_cy_1][img_cx+i];
      for (i=0; i<8; i++)
      for (j=0; j<8; j++)
        img->mprr_c[uv][VERT_PRED_8][i][j] = hline[i];
    }

    // horizontal prediction
    if (mb_available_left)
    {
      for (i=0; i<8; i++)
        vline[i] = image[img_cy+i][img_cx_1];
      for (i=0; i<8; i++)
      for (j=0; j<8; j++)
        img->mprr_c[uv][HOR_PRED_8][i][j] = vline[j];
    }

    // plane prediction
    if (mb_available_up_left)
    {
      ih = 4*(hline[7] - image[img_cy_1][img_cx_1]);
      iv = 4*(vline[7] - image[img_cy_1][img_cx_1]);
      for (i=1;i<4;i++)
      {
        ih += i*(hline[3+i] - hline[3-i]);
        iv += i*(vline[3+i] - vline[3-i]);
      }
      ib=(17*ih+16)>>5;
      ic=(17*iv+16)>>5;

      iaa=16*(hline[7]+vline[7]);
      for (j=0; j<8; j++)
      for (i=0; i<8; i++)
        img->mprr_c[uv][PLANE_8][i][j]=max(0,min(255,(iaa+(i-3)*ib +(j-3)*ic + 16)/32));// store plane prediction
    }
  }

  if (!input->rdopt) // the rd-opt part does not work correctly (see encode_one_macroblock)
  {                       // since ipredmodes could be overwritten => encoder-decoder-mismatches
    // pick lowest cost prediction mode
    min_cost = 1<<20;
    for (mode=DC_PRED_8; mode<=PLANE_8; mode++)
    {
      if ((mode==VERT_PRED_8 && !mb_available_up) ||
          (mode==HOR_PRED_8 && !mb_available_left) ||
          (mode==PLANE_8 && (!mb_available_left || !mb_available_up || !mb_available_up_left)))
        continue;

      cost = 0;
      for (uv=0; uv<2; uv++)
      {
        if(input->InterlaceCodingOption >= MB_CODING && mb_adaptive && img->field_mode)
        {
          if(img->top_field)
            image = imgUV_org_top[uv];
          else
            image = imgUV_org_bot[uv];
        }
        else
          image = imgUV_org[uv];
        for (b4=0,block_y=0; block_y<8; block_y+=4)
        for (block_x=0; block_x<8; block_x+=4,b4++)
        {
          for (k=0,j=block_y; j<block_y+4; j++)
          for (i=block_x; i<block_x+4; i++,k++)
          {
            diff[k] = image[img_cy+j][img_cx+i] - img->mprr_c[uv][mode][i][j];
          }
          cost += SATD(diff, input->hadamard);
        }
      }
      if (cost < min_cost)
      {
        best_mode = mode;
        min_cost = cost;
      }
    }

    currMB->c_ipred_mode = best_mode;
  }
}
/*
void IntraChromaPrediction8x8 (int *mb_up, int *mb_left, int*mb_up_left)
{

  Macroblock *currMB = &img->mb_data[img->current_mb_nr];
  int     s, s0, s1, s2, s3, i, j, k;
  pel_t** image;
  int     block_x, block_y;
  int     mb_nr             = img->current_mb_nr;
  int     mb_available_up;
  int     mb_available_left;
  int     mb_available_up_left;
  int     ih,iv;
  int     ib,ic,iaa;
  int     uv;
  int     hline[8], vline[8];
  int     mode;
  int     best_mode = DC_PRED_8;         //just an initilaization here, should always be overwritten
  int     cost;
  int     min_cost;
  int     diff[16];
	PixelPos up;       //!< pixel position p(0,-1)
  PixelPos left[9];  //!< pixel positions p(-1, -1..8)


	for (i=0;i<9;i++)
  {
    getNeighbour(mb_nr, -1 ,  i-1 , 0, &left[i]);
  }
  
  getNeighbour(mb_nr, 0     ,  -1 , 0, &up);

	
  mb_available_up       = up.available;
  mb_available_up_left  = left[0].available;
 	mb_available_left  = left[1].available;

  if(input->UseConstrainedIntraPred)
  {
		mb_available_up      = up.available ? img->intra_block[up.mb_addr] : 0;
    for (i=1, mb_available_left=1; i<9;i++)
      mb_available_left  &= left[i].available ? img->intra_block[left[i].mb_addr]: 0;
    mb_available_up_left = left[0].available ? img->intra_block[left[0].mb_addr]: 0;
  }

  if (mb_up)
    *mb_up = mb_available_up;
  if (mb_left)
    *mb_left = mb_available_left;
  if( mb_up_left )
    *mb_up_left = mb_available_up_left;

  // compute all chroma intra prediction modes for both U and V
  for (uv=0; uv<2; uv++)
  {
    image = enc_picture->imgUV[uv];

    // DC prediction
    for (block_y=0; block_y<8; block_y+=4)
    for (block_x=0; block_x<8; block_x+=4)
    {
      s=128;
      s0=s1=s2=s3=0;
      //===== get prediction value =====
			switch ((block_y>>1) + (block_x>>2))
      {
      case 0:  //===== TOP LEFT =====
        if      (mb_available_up)    for (i=0;i<4;i++)  s0 += image[up.pos_y][up.pos_x + i];
        if      (mb_available_left)  for (i=1;i<5;i++)  s2 += image[left[i].pos_y][left[i].pos_x];
        if      (mb_available_up && mb_available_left)  s  = (s0+s2+4) >> 3;
        else if (mb_available_up)                       s  = (s0   +2) >> 2;
        else if (mb_available_left)                     s  = (s2   +2) >> 2;
        break;
      case 1: //===== TOP RIGHT =====
        if      (mb_available_up)    for (i=4;i<8;i++)  s1 += image[up.pos_y][up.pos_x + i];
        else if (mb_available_left)  for (i=1;i<5;i++)  s2 += image[left[i].pos_y][left[i].pos_x];
        if      (mb_available_up)                       s  = (s1   +2) >> 2;
        else if (mb_available_left)                     s  = (s2   +2) >> 2;
        break;
      case 2: //===== BOTTOM LEFT =====
        if      (mb_available_left)  for (i=5;i<9;i++)  s3 += image[left[i].pos_y][left[i].pos_x];
        else if (mb_available_up)    for (i=0;i<4;i++)  s0 += image[up.pos_y][up.pos_x + i];
        if      (mb_available_left)                     s  = (s3   +2) >> 2;
        else if (mb_available_up)                       s  = (s0   +2) >> 2;
        break;
      case 3: //===== BOTTOM RIGHT =====
        if      (mb_available_up)    for (i=4;i<8;i++)  s1 += image[up.pos_y][up.pos_x + i];
        if      (mb_available_left)  for (i=5;i<9;i++)  s3 += image[left[i].pos_y][left[i].pos_x];
        if      (mb_available_up && mb_available_left)  s  = (s1+s3+4) >> 3;
        else if (mb_available_up)                       s  = (s1   +2) >> 2;
        else if (mb_available_left)                     s  = (s3   +2) >> 2;
        break;
      }


      //===== prediction =====
      for (j=block_y; j<block_y+4; j++)
      for (i=block_x; i<block_x+4; i++)
      {
        img->mprr_c[uv][DC_PRED_8][i][j] = s;
      }
    }

    // vertical prediction
    if (mb_available_up)
    {
      for (i=0; i<8; i++)
        hline[i] = image[up.pos_y][up.pos_x + i];
      for (i=0; i<8; i++)
      for (j=0; j<8; j++)
        img->mprr_c[uv][VERT_PRED_8][i][j] = hline[i];
    }

    // horizontal prediction 
    if (mb_available_left)
    {
      for (i=1; i<9; i++)
        vline[i] = image[left[i].pos_y][left[i].pos_x];
      for (i=0; i<8; i++)
      for (j=0; j<8; j++)
        img->mprr_c[uv][HOR_PRED_8][i][j] = vline[j+1]; 
    }

    // plane prediction 
    if (mb_available_up_left)
    {
      ih = 4*(hline[7] - image[left[0].pos_y][left[0].pos_x]);
      iv = 4*(vline[7+1] - image[left[0].pos_y][left[0].pos_x]);
      for (i=1;i<4;i++)
      {
        ih += i*(hline[3+i] - hline[3-i]);
        iv += i*(vline[3+i+1] - vline[3-i+1]);
      }
      ib=(17*ih+16)>>5;
      ic=(17*iv+16)>>5;

      iaa=16*(hline[7]+vline[7+1]);
      for (j=0; j<8; j++)
      for (i=0; i<8; i++)
        img->mprr_c[uv][PLANE_8][i][j]=max(0,min(255,(iaa+(i-3)*ib +(j-3)*ic + 16)/32));// store plane prediction
    }
  }

  if (!input->rdopt) // the rd-opt part does not work correctly (see encode_one_macroblock)
  {                       // since ipredmodes could be overwritten => encoder-decoder-mismatches
    // pick lowest cost prediction mode
    min_cost = 1<<20;
		for (i=0;i<8;i++)
    {
      getNeighbour(mb_nr, 0 ,  i , 0, &left[i]);
    }
    for (mode=DC_PRED_8; mode<=PLANE_8; mode++)
    {
      if ((mode==VERT_PRED_8 && !mb_available_up) ||
          (mode==HOR_PRED_8 && !mb_available_left) ||
          (mode==PLANE_8 && (!mb_available_left || !mb_available_up || !mb_available_up_left)))
        continue;

      cost = 0;
      for (uv=0; uv<2; uv++)
      {
        image = imgUV_org[uv];
        for (block_y=0; block_y<8; block_y+=4)
        for (block_x=0; block_x<8; block_x+=4)
        {
          for (k=0,j=block_y; j<block_y+4; j++)
          for (i=block_x; i<block_x+4; i++,k++)
          {
            diff[k] = image[left[j].pos_y][left[j].pos_x+i] - img->mprr_c[uv][mode][i][j];
          }
          cost += SATD(diff, input->hadamard);
        }
      }
      if (cost < min_cost)
      {
        best_mode = mode;
        min_cost = cost;
      }
    }

    currMB->c_ipred_mode = best_mode;
  }
 
}

*/
/*!
 ************************************************************************
 * \brief
 *    Set reference frame information in global arrays
 *    depending on mode decision. Used for motion vector prediction.
 ************************************************************************
 */
void SetRefFrameInfo(int refframe, int bwrefframe)
{
  int i,j;

  if (img->type!=B_SLICE)
  {
      for (j=0; j<4; j++)
      for (i=0; i<4; i++)
      {
        refFrArr[img->block_y+j][img->block_x+i] = refframe;
      }
  }
  else
  {
    for (j=0; j<4; j++)
    for (i=0; i<4; i++)
    {
      fw_refFrArr[img->block_y+j][img->block_x+i] = refframe;
      bw_refFrArr[img->block_y+j][img->block_x+i] = bwrefframe;
    }
  }
}




/*!
 ************************************************************************
 * \brief
 *    Check if all reference frames for a macroblock are zero
 ************************************************************************
 */
int
ZeroRef (Macroblock* currMB)
{
  int i,j;
  int block_y = img->block_y;
  int **frefarr = refFrArr;
  
  if(input->InterlaceCodingOption >= MB_CODING && mb_adaptive && img->field_mode)
  {
    block_y = img->field_block_y;
    frefarr = (img->top_field ? refFrArr_top:refFrArr_bot);
  }

  for (j=0; j<4; j++)
  for (i=0; i<4; i++)
  {
    if (frefarr[block_y+j][img->block_x+i]!=0)
    {
        return 0;
    }
  }
  return 1;
}


/*!
 ************************************************************************
 * \brief
 *    Converts macroblock type to coding value
 ************************************************************************
 */
int
MBType2Value (Macroblock* currMB)
{
  static const int dir1offset[3]    =  { 1,  2, 3};
  static const int dir2offset[3][3] = {{ 0,  4,  8},   // 1. block forward
                                       { 6,  2, 10},   // 1. block backward
                                       {12, 14, 16}};  // 1. block bi-directional

  int mbtype, pdir0, pdir1;

	if (img->type!=B_SLICE)
  {
    if      (currMB->mb_type==I4MB)     return (img->type==I_SLICE ? 0 : 6);
    else if (currMB->mb_type==I16MB)    return (img->type==I_SLICE ? 0 : 6) + img->i16offset;
    else if (currMB->mb_type==P8x8)
    {
      if (input->symbol_mode==UVLC && ZeroRef (currMB))  return 5;
      else                                               return 4;
    }
    else                                return currMB->mb_type;
  }
  else
  {
    mbtype = currMB->mb_type;
    pdir0  = currMB->b8pdir[0];
    pdir1  = currMB->b8pdir[3];

    if      (mbtype==0)       return 0;
    else if (mbtype==I4MB)    return 23;
    else if (mbtype==I16MB)   return 23 + img->i16offset;
    else if (mbtype==P8x8)    return 22;
    else if (mbtype==1)       return dir1offset[pdir0];
    else if (mbtype==2)       return 4 + dir2offset[pdir0][pdir1];
    else                      return 5 + dir2offset[pdir0][pdir1];
  }
}



/*!
 ************************************************************************
 * \brief
 *    Writes intra prediction modes for an 8x8 block
 ************************************************************************
 */
int writeIntra4x4Modes(int only_this_block)
{
  int i,j,bs_x,bs_y,ii,jj;
  int block8x8;
  int rate;
  int ipred_array[16],cont_array[16],ipred_number;
  Macroblock    *currMB     = &img->mb_data[img->current_mb_nr];
  SyntaxElement *currSE     = &img->MB_SyntaxElements[currMB->currSEnr];
  int           *bitCount   = currMB->bitcounter;
  Slice         *currSlice  = img->currentSlice;
  DataPartition *dataPart;
  const int     *partMap    = assignSE2partition[input->partition_mode];

  ipred_number=0;
  for(block8x8=0;block8x8<4;block8x8++)
  {
    if( currMB->b8mode[block8x8]==IBLOCK && (only_this_block<0||only_this_block==block8x8) )
    {
      bs_x=bs_y=4;
      ii=(bs_x>>2); // bug fix for solaris. mwi 
      jj=(bs_y>>2); // bug fix for solaris. mwi
      
      for(j=0;j<2;j+=jj)
      {
        for(i=0;i<2;i+=ii)
        {
          ipred_array[ipred_number]=currMB->intra_pred_modes[(block8x8<<2)|(j<<1)|i];
          cont_array[ipred_number]=(block8x8<<2)+(j<<1)+i;
          ipred_number++;
        }
      }
    }
  }
  rate=0;

  for(i=0;i<ipred_number;i++)
  {
    currMB->IntraChromaPredModeFlag = 1;
    currSE->context = cont_array[i];
    currSE->value1  = ipred_array[i];

#if TRACE
    snprintf(currSE->tracestring, TRACESTRING_SIZE, "Intra mode     = %3d %d",currSE->value1,currSE->context);
#endif

    /*--- set symbol type and function pointers ---*/
    if (input->symbol_mode != UVLC)    currSE->writing = writeIntraPredMode_CABAC;
    currSE->type = SE_INTRAPREDMODE;

    /*--- choose data partition ---*/
		if (img->type != B_SLICE) dataPart = &(currSlice->partArr[partMap[SE_INTRAPREDMODE]]);
    else                                           dataPart = &(currSlice->partArr[partMap[SE_BFRAME]]);
    
    /*--- encode and update rate ---*/
    if (input->symbol_mode == UVLC)    writeSyntaxElement_Intra4x4PredictionMode(currSE, dataPart);
    else                               dataPart->writeSyntaxElement (currSE, dataPart);
    bitCount[BITS_COEFF_Y_MB]+=currSE->len;
    rate += currSE->len;
    currSE++;
    currMB->currSEnr++;
  }

  return rate;
}



/*!
 ************************************************************************
 * \brief
 *    Converts 8x8 block tyoe to coding value
 ************************************************************************
 */
int
B8Mode2Value (int b8mode, int b8pdir)
{
  static const int b8start[8] = {0,0,0,0, 1, 4, 5, 10};
  static const int b8inc  [8] = {0,0,0,0, 1, 2, 2, 1};
	
	if (img->type!=B_SLICE)
  {
    return (b8mode-4);
  }
  else
  {
    return b8start[b8mode] + b8inc[b8mode] * b8pdir;
  }
}



/*!
 ************************************************************************
 * \brief
 *    Codes macroblock header
 ************************************************************************
 */
int writeMBHeader (int rdopt)  // GB CHROMA !!!!!!!!
{
  int             i,j;
  int             mb_nr     = img->current_mb_nr;
  Macroblock*     currMB    = &img->mb_data[mb_nr];
  SyntaxElement *currSE    = &img->MB_SyntaxElements[currMB->currSEnr];
  int*            bitCount  = currMB->bitcounter;
  Slice*          currSlice = img->currentSlice;
  DataPartition*  dataPart;
  const int*      partMap   = assignSE2partition[input->partition_mode];
  int             no_bits   = 0;
  int             mb_y      = img->mb_y;
	int             skip      = currMB->mb_type ? 0:((img->type == B_SLICE) ? !currMB->cbp:1);
  int             mb_type;

  if(input->InterlaceCodingOption >= MB_CODING && mb_adaptive && img->field_mode)
    mb_y = img->top_field ? 2*img->field_mb_y:2*img->field_mb_y+1;

  currMB->IntraChromaPredModeFlag = IS_INTRA(currMB);

  currMB->mb_field = img->field_mode;

  // choose the appropriate data partition
 	if (img->type != B_SLICE)   dataPart = &(currSlice->partArr[partMap[SE_MBTYPE]]);
  else                                             dataPart = &(currSlice->partArr[partMap[SE_BFRAME]]);
  
    //=====  BITS FOR MACROBLOCK MODE =====
    if(img->type == I_SLICE)//GB
    {
      // write mb_aff 
      if(input->InterlaceCodingOption >= MB_CODING && (mb_adaptive) && !skip) // check for copy mode, Krit
      {
        if(WriteFrameFieldMBInHeader)
        {
          currSE->value1 = img->field_mode;
          currSE->type   =  SE_MBTYPE;
          
          if (input->symbol_mode==UVLC)   currSE->mapping = ue_linfo;
          else                            currSE->writing = writeFieldModeInfo_CABAC;
          
#if TRACE
          snprintf(currSE->tracestring, TRACESTRING_SIZE, "Field mode = %3d",img->field_mode);
#endif
          if( input->symbol_mode==UVLC)
          {
            currSE->bitpattern = (img->field_mode ? 1 : 0);
            currSE->len = 1;
            writeSyntaxElement2Buf_Fixed(currSE, dataPart->bitstream);
          }
          else
          {
            dataPart->writeSyntaxElement(currSE, dataPart);
          }

          bitCount[BITS_MB_MODE] += currSE->len;
          no_bits                += currSE->len;
          currSE++;
          currMB->currSEnr++;
        }
      }
      
                        // write mb_type
      currSE->value1  = MBType2Value (currMB);
      currSE->type    = SE_MBTYPE;

      if (input->symbol_mode == UVLC)  currSE->mapping = ue_linfo;
      else                             currSE->writing = writeMB_typeInfo_CABAC;

      dataPart->writeSyntaxElement( currSE, dataPart);
#if TRACE
      snprintf(currSE->tracestring, TRACESTRING_SIZE,   "MB mode(%2d,%2d) = %3d",img->mb_x, img->mb_y, currMB->mb_type);
#endif
      bitCount[BITS_MB_MODE] += currSE->len;
      no_bits                += currSE->len;
      currSE++;
      currMB->currSEnr++;
    }
    else if (input->symbol_mode == CABAC)//GB
    {
      // write mb_skip_flag
      mb_type         = MBType2Value (currMB);
      currSE->value1  = mb_type;
      currSE->value2  = currMB->cbp;
      currSE->type    = SE_MBTYPE;
      currSE->writing = writeMB_skip_flagInfo_CABAC;
      dataPart->writeSyntaxElement( currSE, dataPart);
#if TRACE
			if (img->type == B_SLICE)  snprintf(currSE->tracestring, TRACESTRING_SIZE, "B_MB skipflag(%2d,%2d) = %3d",img->mb_x, img->mb_y, (mb_type!=0 ||currMB->cbp!=0));
      else                     snprintf(currSE->tracestring, TRACESTRING_SIZE,   "MB skipflag(%2d,%2d,%d) = %3d",img->mb_x, img->mb_y, currSE->context,(mb_type!=0));
#endif
      bitCount[BITS_MB_MODE] += currSE->len;
      no_bits                += currSE->len;
      currSE++;
      currMB->currSEnr++;

      // write mb_aff
      if(input->InterlaceCodingOption >= MB_CODING && (mb_adaptive) && !skip) // check for copy mode, Krit
      {
        if(WriteFrameFieldMBInHeader)
        {
          currSE->value1 = img->field_mode;
          currSE->type   =  SE_MBTYPE;

          if (input->symbol_mode==UVLC)   currSE->mapping = ue_linfo;
          else                            currSE->writing = writeFieldModeInfo_CABAC;

          if( input->symbol_mode==UVLC)
          {
            currSE->bitpattern = (img->field_mode ? 1 : 0);
            currSE->len = 1;
            writeSyntaxElement2Buf_Fixed(currSE, dataPart->bitstream);
          }
          else
          {
            dataPart->writeSyntaxElement(currSE, dataPart);
          }
#if TRACE
          snprintf(currSE->tracestring, TRACESTRING_SIZE, "Field mode = %3d",img->field_mode);
#endif
          bitCount[BITS_MB_MODE] += currSE->len;
          no_bits                += currSE->len;
          currSE++;
          currMB->currSEnr++;
        }
      }
      
      // write mb_type
			if (currMB->mb_type != 0 || ((img->type == B_SLICE) && currMB->cbp != 0))
      {
        currSE->value1  = mb_type;
        currSE->type    = SE_MBTYPE;
        currSE->writing = writeMB_typeInfo_CABAC;
        dataPart->writeSyntaxElement( currSE, dataPart);
#if TRACE
        //if (img->type == B_SLICE || img->type == BS_IMG)  snprintf(currSE->tracestring, TRACESTRING_SIZE, "B_MB mode(%2d,%2d) = %3d",img->mb_x, img->mb_y, currMB->mb_type);
        if (img->type == B_SLICE)  snprintf(currSE->tracestring, TRACESTRING_SIZE, "B_MB mode(%2d,%2d) = %3d",img->mb_x, img->mb_y, currMB->mb_type);
        else                     snprintf(currSE->tracestring, TRACESTRING_SIZE,   "MB mode(%2d,%2d) = %3d",img->mb_x, img->mb_y,currMB->mb_type);
#endif
        bitCount[BITS_MB_MODE] += currSE->len;
        no_bits                += currSE->len;
        currSE++;
        currMB->currSEnr++;
      }
    }
    //GB
		else if (currMB->mb_type != 0 || ((img->type == B_SLICE) && currMB->cbp != 0))
    {
      //===== Run Length Coding: Non-Skipped macorblock =====
      currSE->value1  = img->cod_counter;
      currSE->mapping = ue_linfo;
      currSE->type    = SE_MBTYPE;
      dataPart->writeSyntaxElement( currSE, dataPart);
#if TRACE
      snprintf(currSE->tracestring, TRACESTRING_SIZE, "MB runlength = %3d",img->cod_counter);
#endif
      bitCount[BITS_MB_MODE] += currSE->len;
      no_bits                += currSE->len;
      currSE++;
      currMB->currSEnr++;
      
      // Reset cod counter
      img->cod_counter = 0;
      
      // write mb_aff
      if(input->InterlaceCodingOption >= MB_CODING && (mb_adaptive) && !skip) // check for copy mode, Krit
      {
        if(WriteFrameFieldMBInHeader)
        {
          currSE->value1 = img->field_mode;
          currSE->type   =  SE_MBTYPE;
          currSE->mapping = ue_linfo;
          
          //dataPart->writeSyntaxElement(currSE, dataPart);
          currSE->bitpattern = (img->field_mode ? 1 : 0);
          currSE->len = 1;
          writeSyntaxElement2Buf_Fixed(currSE, dataPart->bitstream);

#if TRACE
          snprintf(currSE->tracestring, TRACESTRING_SIZE, "Field mode = %3d",img->field_mode);
#endif
          bitCount[BITS_MB_MODE] += currSE->len;
          no_bits                += currSE->len;
          currSE++;
          currMB->currSEnr++;
        }
      }
      // Put out mb mode
      currSE->value1  = MBType2Value (currMB);
      //if (img->type != B_SLICE && img->type != BS_IMG)
			if (img->type != B_SLICE)
      {
        currSE->value1--;
      }
      currSE->mapping = ue_linfo;
      currSE->type    = SE_MBTYPE;

      dataPart->writeSyntaxElement( currSE, dataPart);
#if TRACE
      //if (img->type == B_SLICE || img->type == BS_IMG)   snprintf(currSE->tracestring, TRACESTRING_SIZE, "B_MB mode(%2d,%2d) = %3d",img->mb_x, img->mb_y, currMB->mb_type);
			if (img->type == B_SLICE)   snprintf(currSE->tracestring, TRACESTRING_SIZE, "B_MB mode(%2d,%2d) = %3d",img->mb_x, img->mb_y, currMB->mb_type);
      else                      snprintf(currSE->tracestring, TRACESTRING_SIZE,   "MB mode(%2d,%2d) = %3d",img->mb_x, img->mb_y,currMB->mb_type);
#endif
      bitCount[BITS_MB_MODE] += currSE->len;
      no_bits                += currSE->len;
      currSE++;
      currMB->currSEnr++;
    }
    else
    {
      //Run Length Coding: Skipped macroblock
      img->cod_counter++;

      // CAVLC
      for (j=0; j < 6; j++)
        for (i=0; i < 4; i++)
          img->nz_coeff [img->mb_x ][mb_y ][i][j]=0;


      if(img->current_mb_nr == img->total_number_mb)
      {
        // Put out run
        currSE->value1  = img->cod_counter;
        currSE->mapping = ue_linfo;
        currSE->type    = SE_MBTYPE;

        dataPart->writeSyntaxElement( currSE, dataPart);
#if TRACE
        snprintf(currSE->tracestring, TRACESTRING_SIZE, "MB runlength = %3d",img->cod_counter);
#endif
        bitCount[BITS_MB_MODE] += currSE->len;
        no_bits                += currSE->len;
        currSE++;
        currMB->currSEnr++;

        // Reset cod counter
        img->cod_counter = 0;
      }
    }

  //===== BITS FOR 8x8 SUB-PARTITION MODES =====
  if (IS_P8x8 (currMB))
  {
    //if (img->type != B_SLICE && img->type != BS_IMG) dataPart = &(currSlice->partArr[partMap[SE_MBTYPE]]);
		if (img->type != B_SLICE) dataPart = &(currSlice->partArr[partMap[SE_MBTYPE]]);
    else                                           dataPart = &(currSlice->partArr[partMap[SE_BFRAME]]);
    
    for (i=0; i<4; i++)
    {
      if (input->symbol_mode==UVLC)   currSE->mapping = ue_linfo;
      else                            currSE->writing = writeB8_typeInfo_CABAC;

      currSE->value1  = B8Mode2Value (currMB->b8mode[i], currMB->b8pdir[i]);
      currSE->type    = SE_MBTYPE;
      dataPart->writeSyntaxElement (currSE, dataPart);
#if TRACE
      snprintf(currSE->tracestring, TRACESTRING_SIZE, "8x8 mode/pdir(%2d) = %3d/%d",
        i,currMB->b8mode[i],currMB->b8pdir[i]);
#endif
      bitCount[BITS_MB_MODE]+= currSE->len;
      no_bits               += currSE->len;
      currSE++;
      currMB->currSEnr++;
    }
  }

 //===== BITS FOR INTRA PREDICTION MODES ====
  no_bits += writeIntra4x4Modes(-1);
  //===== BITS FOR CHROMA INTRA PREDICTION MODE ====
  if (currMB->IntraChromaPredModeFlag)
    no_bits += writeChromaIntraPredMode();
  else if(!rdopt) //GB CHROMA !!!!!
    currMB->c_ipred_mode = DC_PRED_8; //setting c_ipred_mode to default is not the right place here
                                      //resetting in rdopt.c (but where ??)
                                      //with cabac and bframes maybe it could crash without this default
                                      //since cabac needs the right neighborhood for the later MBs

  return no_bits;
}

void write_terminating_bit (short bit)
{
  DataPartition*          dataPart;
  const int*              partMap   = assignSE2partition[input->partition_mode];
  EncodingEnvironmentPtr  eep_dp;

  //--- write non-slice termination symbol if the macroblock is not the first one in its slice ---
  //if (img->type != B_SLICE && img->type != BS_IMG) dataPart = &(img->currentSlice->partArr[partMap[SE_MBTYPE]]);
	if (img->type != B_SLICE) dataPart = &(img->currentSlice->partArr[partMap[SE_MBTYPE]]);
  else                                           dataPart = &(img->currentSlice->partArr[partMap[SE_BFRAME]]);
  dataPart->bitstream->write_flag = 1;
  eep_dp                          = &(dataPart->ee_cabac);
  
  biari_encode_symbol_final(eep_dp, bit); 
#if TRACE
  fprintf (p_trace, "      CABAC terminating bit = %d\n",bit);
#endif

}


/*!
 ************************************************************************
 * \brief
 *    Write chroma intra prediction mode.
 ************************************************************************
 */
int writeChromaIntraPredMode()
{
  Macroblock*     currMB    = &img->mb_data[img->current_mb_nr];
  SyntaxElement*  currSE    = &img->MB_SyntaxElements[currMB->currSEnr];
  Slice*          currSlice = img->currentSlice;
  int*            bitCount  = currMB->bitcounter;
  const int*      partMap   = assignSE2partition[input->partition_mode];
  int             rate      = 0;
  DataPartition*  dataPart;

  //===== BITS FOR CHROMA INTRA PREDICTION MODES
  if (input->symbol_mode==UVLC)   currSE->mapping = ue_linfo;
  else                            currSE->writing = writeCIPredMode_CABAC;
  currSE->value1 = currMB->c_ipred_mode;
  currSE->type = SE_INTRAPREDMODE;
  //if (img->type != B_SLICE && img->type != BS_IMG) dataPart = &(currSlice->partArr[partMap[SE_INTRAPREDMODE]]);
	if (img->type != B_SLICE) dataPart = &(currSlice->partArr[partMap[SE_INTRAPREDMODE]]);
  else                                           dataPart = &(currSlice->partArr[partMap[SE_BFRAME]]);
  dataPart->writeSyntaxElement (currSE, dataPart);
  bitCount[BITS_COEFF_UV_MB] += currSE->len;
  rate                    += currSE->len;
#if TRACE
  snprintf(currSE->tracestring, TRACESTRING_SIZE, "Chroma intra pred mode");
#endif
  currSE++;
  currMB->currSEnr++;

  return rate;
}


/*!
 ************************************************************************
 * \brief
 *    Passes the chosen syntax elements to the NAL
 ************************************************************************
 */
void write_one_macroblock (int eos_bit)
{
  Macroblock* currMB   = &img->mb_data[img->current_mb_nr];
  int*        bitCount = currMB->bitcounter;
  int i,j;
  int mb_y = img->mb_y;

  extern int cabac_encoding;

  if(input->InterlaceCodingOption >= MB_CODING && mb_adaptive && img->field_mode)
    mb_y = img->top_field ? 2*img->field_mb_y:2*img->field_mb_y+1;
  
  //--- write non-slice termination symbol if the macroblock is not the first one in its slice ---
  if (input->symbol_mode==CABAC && img->current_mb_nr!=img->currentSlice->start_mb_nr && eos_bit)
  {
    write_terminating_bit (0);
  }

  cabac_encoding = 1;

  //--- write header ---
  writeMBHeader (0); 

  //  Do nothing more if copy and inter mode
  if ((IS_INTERMV (currMB)  || IS_INTRA (currMB)  ) ||
      //((img->type==B_SLICE || img->type==BS_IMG)     && currMB->cbp != 0)  )
			((img->type==B_SLICE)     && currMB->cbp != 0)  )
  {
    writeMotionInfo2NAL  ();
    writeCBPandLumaCoeff ();
    writeChromaCoeff     ();
  }
  else
  { 
    for (j=0; j < 6; j++)
      for (i=0; i < 4; i++)
        img->nz_coeff [img->mb_x ][mb_y ][i][j]=0;  // CAVLC
  }


  //--- constrain intra prediction ---
  if(input->UseConstrainedIntraPred && (img->type==P_SLICE || img->type==B_SLICE))
  {
    if( !IS_NEWINTRA( currMB ) && currMB->mb_type!=I4MB )
    {
      img->intra_block[img->current_mb_nr] = 0;
    }
  }

  //--- set total bit-counter ---
  bitCount[BITS_TOTAL_MB] = bitCount[BITS_MB_MODE] + bitCount[BITS_COEFF_Y_MB]     + bitCount[BITS_INTER_MB]
                          + bitCount[BITS_CBP_MB]  + bitCount[BITS_DELTA_QUANT_MB] + bitCount[BITS_COEFF_UV_MB];
  stat->bit_slice += bitCount[BITS_TOTAL_MB];

  cabac_encoding = 0;
}


/*!
 ************************************************************************
 * \brief
 *    Sets context for reference frame parameter
 ************************************************************************
 */
int BType2CtxRef (int btype)
{
  if (btype<4)   return 0;
  else           return 1;
}


/*!
 ************************************************************************
 * \brief
 *    Codes the reference frame
 ************************************************************************
 */
int writeReferenceFrame (int mode, int i, int j, int fwd_flag, int  ref)
{
  Macroblock*     currMB    = &img->mb_data[img->current_mb_nr];
  SyntaxElement*  currSE    = &img->MB_SyntaxElements[currMB->currSEnr];
  Slice*          currSlice = img->currentSlice;
  int*            bitCount  = currMB->bitcounter;
  const int*      partMap   = assignSE2partition[input->partition_mode];
  int             rate      = 0;
  DataPartition*  dataPart;
  int             num_ref   = ( fwd_flag ? img->num_ref_idx_l0_active : img->num_ref_idx_l1_active);
  int             flag_mode = 0;

  if( num_ref == 1 )
  {
    return 0;
  }
  if ( num_ref == 2 )
  {
    flag_mode = 1;
  }
  /*
  if ((fwd_flag && img->num_ref_pic_active_fwd_minus1==0) || (!fwd_flag &&img->num_ref_pic_active_bwd_minus1==0))
  {
    return 0;
  }
  */

  currSE->value1 = ref;

	currSE->type   = ((img->type==B_SLICE) ? SE_BFRAME : SE_REFFRAME);
  dataPart = &(currSlice->partArr[partMap[currSE->type]]);
  if (input->symbol_mode == UVLC)
  {
    if( flag_mode )
    {
      currSE->bitpattern = 1 - currSE->value1;
      currSE->len = 1;
      writeSyntaxElement2Buf_Fixed(currSE, dataPart->bitstream);
    }
    else
    {
      currSE->mapping = ue_linfo;
      dataPart->writeSyntaxElement (currSE, dataPart);
    }
  }
  else
  {
    currSE->context = BType2CtxRef (mode);
    img->subblock_x = i; // position used for context determination
    img->subblock_y = j; // position used for context determination
    currSE->writing = writeRefFrame_CABAC;
    currSE->value2 = (fwd_flag)? LIST_0:LIST_1;
    dataPart->writeSyntaxElement (currSE, dataPart);
  }

  bitCount[BITS_INTER_MB] += currSE->len;
  rate                    += currSE->len;
#if TRACE
  if (fwd_flag)
  {
    snprintf(currSE->tracestring, TRACESTRING_SIZE, "Fwd Ref frame no %d", currSE->value1);
  }
  else
  {
    snprintf(currSE->tracestring, TRACESTRING_SIZE, "Bwd Ref frame no %d", currSE->value1);
  }
#endif
  currSE++;
  currMB->currSEnr++;

  return rate;
}


/*!
 ************************************************************************
 * \brief
 *    Writes motion vectors of an 8x8 block
 ************************************************************************
 */
int writeMotionVector8x8 (int  i0,
                          int  j0,
                          int  i1,
                          int  j1,
                          int  refframe,
                          int  dmv_flag,
                          int  fwd_flag,
                          int  mv_mode)
{
  int            i, j, k, l, m;
  int            curr_mvd;
  DataPartition* dataPart;
  int            bwflag     = ((refframe<0 || (!fwd_flag))?1:0);
  int            rate       = 0;
  int            step_h     = input->blc_size[mv_mode][0] >> 2;
  int            step_v     = input->blc_size[mv_mode][1] >> 2;
  Macroblock*    currMB     = &img->mb_data[img->current_mb_nr];
  SyntaxElement* currSE     = &img->MB_SyntaxElements[currMB->currSEnr];
  Slice*         currSlice  = img->currentSlice;
  int*           bitCount   = currMB->bitcounter;
  const int*     partMap    = assignSE2partition[input->partition_mode];
  int            refindex   = (refframe<0 ? 0 : refframe);
  int*****       all_mv     = (fwd_flag ? img->all_mv : img->all_bmv);
  //int*****       pred_mv    = ((img->type!=B_SLICE && img->type!=BS_IMG) ? img->mv : (fwd_flag ? img->p_fwMV : img->p_bwMV));
	int*****       pred_mv    = ((img->type!=B_SLICE) ? img->mv : (fwd_flag ? img->p_fwMV : img->p_bwMV));
  int*****       abp_all_dmv = img->abp_all_dmv;
  if (!fwd_flag) bwflag = 1;
  if(input->InterlaceCodingOption >= MB_CODING && mb_adaptive && img->field_mode)
  {
    if(img->top_field)
    {
      all_mv     = (fwd_flag ? img->all_mv_top : img->all_bmv_top);
      //pred_mv    = ((img->type!=B_SLICE && img->type!=BS_IMG) ? img->mv_top : fwd_flag ? img->p_fwMV_top : img->p_bwMV_top);
			pred_mv    = ((img->type!=B_SLICE) ? img->mv_top : fwd_flag ? img->p_fwMV_top : img->p_bwMV_top);
      abp_all_dmv = img->abp_all_dmv_top;
    }
    else
    {
      all_mv     = (fwd_flag ? img->all_mv_bot : img->all_bmv_bot);
//      pred_mv    = ((img->type!=B_SLICE && img->type!=BS_IMG) ? img->mv_bot : fwd_flag ? img->p_fwMV_bot : img->p_bwMV_bot);
			pred_mv    = ((img->type!=B_SLICE) ? img->mv_bot : fwd_flag ? img->p_fwMV_bot : img->p_bwMV_bot);
      abp_all_dmv = img->abp_all_dmv_bot;
    }
  }

  for (j=j0; j<j1; j+=step_v)
  for (i=i0; i<i1; i+=step_h)
  {
    for (k=0; k<2; k++) 
    {
      // (img->type==BS_IMG && dmv_flag)
			if (img->type==B_SLICE && img->nal_reference_idc>0 && dmv_flag)
      {
        curr_mvd = abp_all_dmv[i][j][refindex][mv_mode][k];
      }
      //else if (img->type==BS_IMG && !dmv_flag)
			else if (img->type==B_SLICE && img->nal_reference_idc>0 && !dmv_flag)
      {
        curr_mvd = all_mv[i][j][refindex][mv_mode][k] - pred_mv[i][j][refindex][mv_mode][k];
      }
      else
      {
        curr_mvd = all_mv[i][j][refindex][mv_mode][k] - pred_mv[i][j][refindex][mv_mode][k];
      }

      //--- store (oversampled) mvd ---
      for (l=0; l < step_v; l++) 
      for (m=0; m < step_h; m++)    currMB->mvd[bwflag][j+l][i+m][k] = curr_mvd;

      currSE->value1 = curr_mvd;
      //currSE->type   = ((img->type==B_SLICE || img->type==BS_IMG) ? SE_BFRAME : SE_MVD);
			currSE->type   = ((img->type==B_SLICE) ? SE_BFRAME : SE_MVD);
      if (input->symbol_mode == UVLC)
      {
        currSE->mapping = se_linfo;
      }
      else
      {
        img->subblock_x = i; // position used for context determination
        img->subblock_y = j; // position used for context determination
        currSE->value2  = 2*k+bwflag; // identifies the component and the direction; only used for context determination
        currSE->writing = writeMVD_CABAC;
      }  
      //dataPart = &(currSlice->partArr[partMap[(img->type==B_SLICE || img->type==BS_IMG)? SE_BFRAME : SE_MVD]]);
			dataPart = &(currSlice->partArr[partMap[(img->type==B_SLICE)? SE_BFRAME : SE_MVD]]);
      dataPart->writeSyntaxElement (currSE, dataPart);
#if TRACE
      //if (img->type==BS_IMG && dmv_flag)
			if (img->type==B_SLICE && img->nal_reference_idc>0 && dmv_flag)
      {
        if (!fwd_flag)
        {
          snprintf(currSE->tracestring, TRACESTRING_SIZE, "BDMV(%d) = %3d  (org_mv %3d)",k, curr_mvd, all_mv[i][j][refindex][mv_mode][k]);
        }
      }
      else
      {
        if (fwd_flag)
        {
          snprintf(currSE->tracestring, TRACESTRING_SIZE, "FMVD(%d) = %3d  (org_mv %3d pred_mv %3d) %d",k, curr_mvd, all_mv[i][j][refindex][mv_mode][k], pred_mv[i][j][refindex][mv_mode][k],currSE->value2);
        }
        else
        {
          snprintf(currSE->tracestring, TRACESTRING_SIZE, "BMVD(%d) = %3d  (org_mv %3d pred_mv %3d)",k, curr_mvd, all_mv[i][j][refindex][mv_mode][k], pred_mv[i][j][refindex][mv_mode][k]);
        }
      }
#endif
      bitCount[BITS_INTER_MB] += currSE->len;
      rate                    += currSE->len;
      currSE++;  
      currMB->currSEnr++;
    }
  }

  return rate;
}


/*!
 ************************************************************************
 * \brief
 *    Writes motion info
 ************************************************************************
 */
int writeMotionInfo2NAL ()
{
  int k, j0, i0, refframe;

  Macroblock*     currMB    = &img->mb_data[img->current_mb_nr];
  int             no_bits   = 0;

  //int   bframe          = (img->type==B_SLICE || img->type==BS_IMG);
	int   bframe          = (img->type==B_SLICE);
  //int** refframe_array  = ((img->type==B_SLICE || img->type==BS_IMG) ? fw_refFrArr : refFrArr);
	int** refframe_array  = ((img->type==B_SLICE) ? fw_refFrArr : refFrArr);
  int** bw_refframe_array  = bw_refFrArr;
  int** abp_type_array  = abp_type_FrArr;
#ifdef _ADDITIONAL_REFERENCE_FRAME_
  int   multframe       = (input->num_reference_frames>1 || input->add_ref_frame>0); 
#else
  int   multframe       = (input->num_reference_frames>1); 
#endif
  int   step_h0         = (input->blc_size[IS_P8x8(currMB) ? 4 : currMB->mb_type][0] >> 2);
  int   step_v0         = (input->blc_size[IS_P8x8(currMB) ? 4 : currMB->mb_type][1] >> 2);
  int   block_y     = img->block_y;

  if(input->InterlaceCodingOption==2)
  {
#ifdef _ADDITIONAL_REFERENCE_FRAME_
  multframe       = (2*input->num_reference_frames>1 || input->add_ref_frame>0); 
#else
  multframe       = (2*input->num_reference_frames>1); 
#endif
  }


  if(input->InterlaceCodingOption >= MB_CODING && mb_adaptive && img->field_mode)
  {
    block_y = img->field_block_y;

    if(img->top_field)
    {
      //refframe_array  = ((img->type==B_SLICE || img->type==BS_IMG) ? fw_refFrArr_top : refFrArr_top);
			refframe_array  = ((img->type==B_SLICE) ? fw_refFrArr_top : refFrArr_top);
      bw_refframe_array  = bw_refFrArr_top;
      abp_type_array  = abp_type_FrArr_top;
    }
    else
    {
      //refframe_array  = ((img->type==B_SLICE || img->type==BS_IMG) ? fw_refFrArr_bot : refFrArr_bot);
			refframe_array  = ((img->type==B_SLICE) ? fw_refFrArr_bot : refFrArr_bot);
      bw_refframe_array  = bw_refFrArr_bot;
      abp_type_array  = abp_type_FrArr_bot;
    }
  }




  //=== If multiple ref. frames, write reference frame for the MB ===
  if (IS_INTERMV (currMB) && multframe)
  {
    // if UVLC is turned on, a 8x8 macroblock with all ref=0 in a P-frame is signalled in macroblock mode
    if (!IS_P8x8 (currMB) || !ZeroRef (currMB) || input->symbol_mode==CABAC || bframe)
    {
      for (j0=0; j0<4; j0+=step_v0)
      for (i0=0; i0<4; i0+=step_h0)
      {
        k=j0+(i0/2);

        if ((currMB->b8pdir[k]==0 || currMB->b8pdir[k]==2) && currMB->b8mode[k]!=0)//has forward vector
        {
          no_bits += writeReferenceFrame (currMB->b8mode[k], i0, j0, 1, refframe_array[block_y+j0][img->block_x+i0]);
        }
      }
        for (j0=0; j0<4; j0+=step_v0)
        for (i0=0; i0<4; i0+=step_h0)
        {
          k=j0+(i0/2);
          if ((currMB->b8pdir[k]==1 || currMB->b8pdir[k]==2) && currMB->b8mode[k]!=0)//has backward vector
          {
            no_bits += writeReferenceFrame (currMB->b8mode[k], i0, j0, 0, bw_refFrArr[block_y+j0][img->block_x+i0]);
          }
        }
    }
  }

  //===== write forward motion vectors =====
  if (IS_INTERMV (currMB))
  {
    for (j0=0; j0<4; j0+=step_v0)
    for (i0=0; i0<4; i0+=step_h0)
    {
      k=j0+(i0/2);
      if ((currMB->b8pdir[k]==0 || currMB->b8pdir[k]==2) && currMB->b8mode[k]!=0)//has forward vector
      {
        refframe  = refframe_array[block_y+j0][img->block_x+i0];
        //if(img->type==BS_IMG && currMB->b8pdir[k]==2  && (abp_type_array[block_y+j0][img->block_x+i0]==1 || abp_type_array[block_y+j0][img->block_x+i0]==2))
				if(img->type==B_SLICE && img->nal_reference_idc>0 && currMB->b8pdir[k]==2  && (abp_type_array[block_y+j0][img->block_x+i0]==1 || abp_type_array[block_y+j0][img->block_x+i0]==2))
        {
          no_bits  += writeMotionVector8x8 (i0, j0, i0+step_h0, j0+step_v0, refframe, 0/*MV*/, 1, currMB->b8mode[k]);
        }
        else
        {
          no_bits  += writeMotionVector8x8 (i0, j0, i0+step_h0, j0+step_v0, refframe, 0, 1, currMB->b8mode[k]);
        }
      }
    }
  }


  //===== write backward motion vectors =====
  if (IS_INTERMV (currMB) && bframe)
  {
    for (j0=0; j0<4; j0+=step_v0)
    for (i0=0; i0<4; i0+=step_h0)
    {
      k=j0+(i0/2);
      if ((currMB->b8pdir[k]==1 || currMB->b8pdir[k]==2) && currMB->b8mode[k]!=0)//has backward vector
      {
        refframe  = bw_refframe_array[block_y+j0][img->block_x+i0];
//        if(img->type==BS_IMG && currMB->b8pdir[k]==2 && (abp_type_array[block_y+j0][img->block_x+i0]==1 || abp_type_array[block_y+j0][img->block_x+i0]==2))
				if(img->type==B_SLICE && img->nal_reference_idc>0 && currMB->b8pdir[k]==2 && (abp_type_array[block_y+j0][img->block_x+i0]==1 || abp_type_array[block_y+j0][img->block_x+i0]==2))
        {
          no_bits  += writeMotionVector8x8 (i0, j0, i0+step_h0, j0+step_v0, refframe, 1/*DMV*/, 0, currMB->b8mode[k]);
        }
        else
        {
          no_bits  += writeMotionVector8x8 (i0, j0, i0+step_h0, j0+step_v0, refframe, 0, 0, currMB->b8mode[k]);
        }
      }
    }
  }
  return no_bits;
}



/*!
 ************************************************************************
 * \brief
 *    Writes chrominance coefficients
 ************************************************************************
 */
int writeChromaCoeff ()
{
  int             rate      = 0;
  Macroblock*     currMB    = &img->mb_data[img->current_mb_nr];
  SyntaxElement*  currSE    = &img->MB_SyntaxElements[currMB->currSEnr];
  int*            bitCount  = currMB->bitcounter;
  Slice*          currSlice = img->currentSlice;
  const int*      partMap   = assignSE2partition[input->partition_mode];
  int             cbp       = currMB->cbp;
  DataPartition*  dataPart;

  int   level, run;
  int   i, j, k, uv, mb_x, mb_y, i1, ii, j1, jj;
  int   b8, b4, param;
  int*  ACLevel;
  int*  ACRun;
  int*  DCLevel;
  int*  DCRun;


  //=====
  //=====   D C - C O E F F I C I E N T S
  //=====
  if (cbp > 15)  // check if any chroma bits in coded block pattern is set
  {
    for (uv=0; uv < 2; uv++)
    {

      if (input->symbol_mode == UVLC)
      {
        param = uv;
        rate += writeCoeff4x4_CAVLC (CHROMA_DC, 0, 0, param);
          // CAVLC
      }
      else

      {

        DCLevel = img->cofDC[uv+1][0];
        DCRun   = img->cofDC[uv+1][1];

        level=1;
        for (k=0; k < 5 && level != 0; ++k)
        {
          level = currSE->value1 = DCLevel[k]; // level
          run   = currSE->value2 = DCRun  [k]; // run

          if (input->symbol_mode == UVLC)   currSE->mapping = levrun_linfo_c2x2;
          else                              currSE->writing = writeRunLevel_CABAC;

          currSE->context     = CHROMA_DC;
          currSE->type        = (IS_INTRA(currMB) ? SE_CHR_DC_INTRA : SE_CHR_DC_INTER);
          img->is_intra_block =  IS_INTRA(currMB);
          img->is_v_block     = uv;
    
          // choose the appropriate data partition
          //if (img->type!=B_SLICE && img->type!=BS_IMG)   dataPart = &(currSlice->partArr[partMap[currSE->type]]);
					if (img->type!=B_SLICE)   dataPart = &(currSlice->partArr[partMap[currSE->type]]);
          else                                         dataPart = &(currSlice->partArr[partMap[SE_BFRAME]]);
    
          dataPart->writeSyntaxElement (currSE, dataPart);
          bitCount[BITS_COEFF_UV_MB] += currSE->len;
          rate                       += currSE->len;
#if TRACE
          snprintf(currSE->tracestring, TRACESTRING_SIZE, "2x2 DC Chroma %2d: level =%3d run =%2d",k, level, run);
#endif
          // proceed to next SE 
          currSE++;  
          currMB->currSEnr++;
      }
      }
    }
  }


  //=====
  //=====   A C - C O E F F I C I E N T S
  //=====
  uv=-1;   
  if (cbp >> 4 == 2) // check if chroma bits in coded block pattern = 10b
    {  
    for (mb_y=4; mb_y < 6; mb_y += 2)
    for (mb_x=0; mb_x < 4; mb_x += 2)
      {
      for (j=mb_y; j < mb_y+2; j++)
      {
        jj=j/2;
        j1=j-4;
        for (i=mb_x; i < mb_x+2; i++)
        {
          b8      = 4 + i/2;
          b4      = 2*(j/5)+ (i%2);

          if (input->symbol_mode == UVLC)
          {
            param = i << 4 | j;
            rate += writeCoeff4x4_CAVLC (CHROMA_AC, b8, b4, param);
            // CAVLC
          }
          else

          {

            ACLevel = img->cofAC[b8][b4][0];
            ACRun   = img->cofAC[b8][b4][1];

            ii=i/2;
            i1=i%2;
            level=1;
            uv++;

            img->subblock_y = b4/2;
            img->subblock_x = b4%2;

            for (k=0; k < 16 && level != 0; k++)
            {
              level = currSE->value1 = ACLevel[k]; // level
              run   = currSE->value2 = ACRun  [k]; // run

              if (input->symbol_mode == UVLC)   currSE->mapping = levrun_linfo_inter;
              else                              currSE->writing = writeRunLevel_CABAC;
            
              currSE->context     = CHROMA_AC;
              currSE->type        = (IS_INTRA(currMB) ? SE_CHR_AC_INTRA : SE_CHR_AC_INTER);
              img->is_intra_block =  IS_INTRA(currMB);
              img->is_v_block     = (uv>=4);

              // choose the appropriate data partition
              //if (img->type!=B_SLICE && img->type!=BS_IMG)   dataPart = &(currSlice->partArr[partMap[currSE->type]]);
							if (img->type!=B_SLICE)   dataPart = &(currSlice->partArr[partMap[currSE->type]]);
              else                                         dataPart = &(currSlice->partArr[partMap[SE_BFRAME]]);
        
              dataPart->writeSyntaxElement (currSE, dataPart);
              bitCount[BITS_COEFF_UV_MB] += currSE->len;
              rate                       += currSE->len;
#if TRACE
              snprintf(currSE->tracestring, TRACESTRING_SIZE, "AC Chroma %2d: level =%3d run =%2d",k, level, run);
#endif

              // proceed to next SE 
              currSE++;  
              currMB->currSEnr++;
            }
          }
        }
      }
    }
  }

  return rate;
}



/*!
 ************************************************************************
 * \brief
 *    Writes Luma coeff of an 4x4 block
 ************************************************************************
 */
int writeLumaCoeff4x4_CABAC (int b8, int b4, int intra4x4mode)
{
  int             rate      = 0;
  Macroblock*     currMB    = &img->mb_data[img->current_mb_nr];
  SyntaxElement*  currSE    = &img->MB_SyntaxElements[currMB->currSEnr];
  Slice*          currSlice = img->currentSlice;
  const int*      partMap   = assignSE2partition[input->partition_mode];
  int*            bitCount  = currMB->bitcounter;
  DataPartition*  dataPart;

  int   level, run;
  int   k;
  int*  ACLevel = img->cofAC[b8][b4][0];
  int*  ACRun   = img->cofAC[b8][b4][1];

  img->subblock_x = ((b8&0x1)==0)?(((b4&0x1)==0)?0:1):(((b4&0x1)==0)?2:3); // horiz. position for coeff_count context
  img->subblock_y = (b8<2)?((b4<2)?0:1):((b4<2)?2:3); // vert.  position for coeff_count context

  level=1; // get inside loop
  for(k=0; k<=16 && level !=0; k++)
  {
    level = currSE->value1 = ACLevel[k]; // level
    run   = currSE->value2 = ACRun  [k]; // run
      
    currSE->writing = writeRunLevel_CABAC;

    currSE->context     = LUMA_4x4;
    currSE->type        = (k==0 ? (intra4x4mode?SE_LUM_DC_INTRA:SE_LUM_DC_INTER) : (intra4x4mode?SE_LUM_AC_INTRA:SE_LUM_AC_INTER));
    img->is_intra_block = intra4x4mode;

    // choose the appropriate data partition
    //if (img->type != B_SLICE && img->type != BS_IMG)    dataPart = &(currSlice->partArr[partMap[currSE->type]]);
		if (img->type != B_SLICE)    dataPart = &(currSlice->partArr[partMap[currSE->type]]);
    else                                              dataPart = &(currSlice->partArr[partMap[SE_BFRAME]]);
          
    dataPart->writeSyntaxElement (currSE, dataPart);
    bitCount[BITS_COEFF_Y_MB] += currSE->len;
    rate                      += currSE->len;
#if TRACE
    snprintf(currSE->tracestring, TRACESTRING_SIZE, "Luma sng(%2d) level =%3d run =%2d", k, level,run);
#endif
    /* proceed to next SE */
    currSE++;  
    currMB->currSEnr++;
  }

  return rate;
}


/*!
 ************************************************************************
 * \brief
 *    Writes Luma Coeff of an 8x8 block
 ************************************************************************
 */
int writeLumaCoeff8x8 (int block8x8, int intra4x4mode)
{
  int  block4x4, rate = 0;

  for (block4x4=0; block4x4<4; block4x4++)
  {
    if (input->symbol_mode == UVLC )
      rate += writeCoeff4x4_CAVLC (LUMA, block8x8, block4x4, 0);// CAVLC
    else
      rate += writeLumaCoeff4x4_CABAC (block8x8, block4x4, intra4x4mode);
  }

  return rate;
}



/*!
 ************************************************************************
 * \brief
 *    Writes CBP, DQUANT, and Luma Coefficients of an macroblock
 ************************************************************************
 */
int writeCBPandLumaCoeff ()
{
  int             mb_x, mb_y, i, j, k;
  int             level, run;
  int             rate      = 0;
  Macroblock*     currMB    = &img->mb_data[img->current_mb_nr];
  int*            bitCount  = currMB->bitcounter;
  SyntaxElement*  currSE    = &img->MB_SyntaxElements[currMB->currSEnr];
  Slice*          currSlice = img->currentSlice;
  const int*      partMap   = assignSE2partition[input->partition_mode];
  int             cbp       = currMB->cbp;
  DataPartition*  dataPart;

  int   b8, b4;
  int*  DCLevel = img->cofDC[0][0];
  int*  DCRun   = img->cofDC[0][1];
  int*  ACLevel;
  int*  ACRun;
  int   mb_ypos = img->mb_y;
  if(input->InterlaceCodingOption >= MB_CODING && mb_adaptive && img->field_mode)
    mb_ypos = img->top_field ? 2*img->field_mb_y:2*img->field_mb_y+1;


  if (!IS_NEWINTRA (currMB))
  {
    //=====   C B P   =====
    //---------------------
    currSE->value1 = cbp;
    
    if (IS_OLDINTRA (currMB) || currMB->mb_type == SI4MB)
    {
      if (input->symbol_mode == UVLC)  currSE->mapping = cbp_linfo_intra;
      currSE->type = SE_CBP_INTRA;
    }
    else
    {
      if (input->symbol_mode == UVLC)  currSE->mapping = cbp_linfo_inter;
      currSE->type = SE_CBP_INTER;
    }
    if (input->symbol_mode == CABAC)   currSE->writing = writeCBP_CABAC;
                      
    // choose the appropriate data partition
    //if (img->type != B_SLICE && img->type != BS_IMG) dataPart = &(currSlice->partArr[partMap[currSE->type]]);
		if (img->type != B_SLICE) dataPart = &(currSlice->partArr[partMap[currSE->type]]);
    else                                           dataPart = &(currSlice->partArr[partMap[SE_BFRAME]]);
    
    dataPart->writeSyntaxElement(currSE, dataPart);
    bitCount[BITS_CBP_MB] += currSE->len;
    rate                  += currSE->len;
#if TRACE
    snprintf(currSE->tracestring, TRACESTRING_SIZE, "CBP (%2d,%2d) = %3d",img->mb_x, img->mb_y, cbp);
#endif
    // proceed to next SE
    currSE++;  
    currMB->currSEnr++;
  }

  //=====   DQUANT   =====
  //----------------------
  if (cbp!=0 || IS_NEWINTRA (currMB))
  {
    currSE->value1 = currMB->delta_qp;

    if (input->symbol_mode==UVLC)   currSE->mapping = se_linfo;
    else                            currSE->writing = writeDquant_CABAC;

    if (IS_INTER (currMB))  currSE->type = SE_DELTA_QUANT_INTER;
    else                    currSE->type = SE_DELTA_QUANT_INTRA;


    // choose the appropriate data partition
//    if (img->type != B_SLICE && img->type != BS_IMG)   dataPart = &(img->currentSlice->partArr[partMap[currSE->type]]);
		if (img->type != B_SLICE)   dataPart = &(img->currentSlice->partArr[partMap[currSE->type]]);
    else                                             dataPart = &(img->currentSlice->partArr[partMap[SE_BFRAME]]);

    dataPart->writeSyntaxElement(  currSE, dataPart);
    bitCount[BITS_DELTA_QUANT_MB] += currSE->len;
    rate                          += currSE->len;
#if TRACE
    snprintf(currSE->tracestring, TRACESTRING_SIZE, "Delta QP (%2d,%2d) = %3d",img->mb_x, img->mb_y, currMB->delta_qp);
#endif
    // proceed to next SE
    currSE++;
    currMB->currSEnr++;
  }

  {
    for (j=0; j < 6; j++)
      for (i=0; i < 4; i++)
        img->nz_coeff [img->mb_x ][mb_ypos ][i][j]=0;  // CAVLC
  }

  if (!IS_NEWINTRA (currMB))
  {
    //=====  L U M I N A N C E   =====
    //--------------------------------
    for (i=0; i<4; i++)  if (cbp & (1<<i))
    {
        rate += writeLumaCoeff8x8 (i, (currMB->b8mode[i]==IBLOCK));
    }
  }
  else
  {
    //=====  L U M I N A N C E   f o r   1 6 x 1 6   =====
    //----------------------------------------------------
    // DC coeffs
    if (input->symbol_mode == UVLC)
    {
      rate += writeCoeff4x4_CAVLC (LUMA_INTRA16x16DC, 0, 0, 0);  // CAVLC
    }
    else
    {
      level=1; // get inside loop
      for (k=0; k<=16 && level!=0; k++)
      {
        level = currSE->value1 = DCLevel[k]; // level
        run   = currSE->value2 = DCRun  [k]; // run

        if (input->symbol_mode == UVLC)
        {
          currSE->mapping = levrun_linfo_inter;
        }else{
          currSE->writing = writeRunLevel_CABAC;
        }

        currSE->context     = LUMA_16DC;
        currSE->type        = SE_LUM_DC_INTRA;   // element is of type DC
        img->is_intra_block = 1;

        // choose the appropriate data partition
      //  if (img->type != B_SLICE && img->type != BS_IMG)   dataPart = &(currSlice->partArr[partMap[currSE->type]]);
				if (img->type != B_SLICE)   dataPart = &(currSlice->partArr[partMap[currSE->type]]);
        else                                             dataPart = &(currSlice->partArr[partMap[SE_BFRAME]]);
    
        dataPart->writeSyntaxElement (currSE, dataPart);
        bitCount[BITS_COEFF_Y_MB] += currSE->len;
        rate                      += currSE->len;
#if TRACE
        snprintf(currSE->tracestring, TRACESTRING_SIZE, "DC luma 16x16 sng(%2d) level =%3d run =%2d", k, level, run);
#endif
        // proceed to next SE
        currSE++;
        currMB->currSEnr++;
      }
    }

    // AC coeffs
    if (cbp & 15)
    {
      for (mb_y=0; mb_y < 4; mb_y += 2)
      for (mb_x=0; mb_x < 4; mb_x += 2)
      for (j=mb_y; j < mb_y+2; j++)
      for (i=mb_x; i < mb_x+2; i++)
      {
        b8      = 2*(j/2) + (i/2);
        b4      = 2*(j%2) + (i%2);
        if (input->symbol_mode == UVLC)
        {
          rate += writeCoeff4x4_CAVLC (LUMA_INTRA16x16AC, b8, b4, 0);  // CAVLC
        }
        else
        {
          ACLevel = img->cofAC[b8][b4][0];
          ACRun   = img->cofAC[b8][b4][1];

          img->subblock_y = j;
          img->subblock_x = i;

          level=1; // get inside loop
          for (k=0;k<16 && level !=0;k++)
          {
            level = currSE->value1 = ACLevel[k]; // level
            run   = currSE->value2 = ACRun  [k]; // run

            if (input->symbol_mode == UVLC)
            {
              currSE->mapping = levrun_linfo_inter;
            }else{
              currSE->writing = writeRunLevel_CABAC;
            }
            currSE->context     = LUMA_16AC;
            currSE->type        = SE_LUM_AC_INTRA;   // element is of type AC
            img->is_intra_block = 1;

            // choose the appropriate data partition
//            if (img->type != B_SLICE && img->type != BS_IMG)   dataPart = &(currSlice->partArr[partMap[currSE->type]]);
            if (img->type != B_SLICE)   dataPart = &(currSlice->partArr[partMap[currSE->type]]);
            else                                             dataPart = &(currSlice->partArr[partMap[SE_BFRAME]]);

            dataPart->writeSyntaxElement (currSE, dataPart);
            bitCount[BITS_COEFF_Y_MB] += currSE->len;
            rate                      += currSE->len;
#if TRACE
            snprintf(currSE->tracestring, TRACESTRING_SIZE, "AC luma 16x16 sng(%2d) level =%3d run =%2d", k, level, run);
#endif
            // proceed to next SE
            currSE++;
            currMB->currSEnr++;
          }
        }
      }
    }
  }

  return rate;
}


/*!
 ************************************************************************
 * \brief
 *    Get the Prediction from the Neighboring BLocks for Number of Nonzero Coefficients 
 *    
 *    Luma Blocks
 ************************************************************************
 */
int predict_nnz(int i,int j)
{
  int Left_block,Top_block, pred_nnz;
  int cnt=0;
//  int mb_ypos = img->mb_y;
  int decr    = 1;

  int mb_nr    = img->current_mb_nr;
  int mb_width = img->width/16;

  int mb_available_up   = (img->mb_y == 0 ) ? 0 : (img->mb_data[mb_nr].slice_nr == img->mb_data[mb_nr-mb_width  ].slice_nr);
  int mb_available_left = (img->mb_x == 0 ) ? 0 : (img->mb_data[mb_nr].slice_nr == img->mb_data[mb_nr-1         ].slice_nr);

  if (i)
    Left_block= img->nz_coeff [img->mb_x ][img->mb_y ][i-1][j];
  else
    Left_block= mb_available_left ? img->nz_coeff [img->mb_x-1 ][img->mb_y ][3][j] : -1;
//    Left_block= img->mb_x > 0 ? img->nz_coeff [img->mb_x-1 ][img->mb_y ][3][j] : -1;

  if (j)
    Top_block=  img->nz_coeff [img->mb_x ][img->mb_y ][i][j-1];
  else
  Top_block=  mb_available_up ? img->nz_coeff [img->mb_x ][img->mb_y-decr ][i][3] : -1;
//  Top_block=  mb_ypos > 0 ? img->nz_coeff [img->mb_x ][img->mb_y-decr ][i][3] : -1;
  
  
  pred_nnz=0;
  if (Left_block>-1)
  {
    pred_nnz=Left_block;
    cnt++;
  }
  if (Top_block>-1)
  {
    pred_nnz+=Top_block;
    cnt++;
  }

  if (cnt==2)
    pred_nnz++;

  if (cnt)
    pred_nnz/=cnt; 
  return pred_nnz;
}
/*!
 ************************************************************************
 * \brief
 *    Get the Prediction from the Neighboring BLocks for Number of Nonzero Coefficients 
 *    
 *    Chroma Blocks   
 ************************************************************************
 */
int predict_nnz_chroma(int i,int j)
{
  int Left_block,Top_block, pred_nnz;
  int cnt=0;
  int mb_y = img->mb_y;
//  int mb_ypos = img->mb_y;
  int decr    = 1;

/*  if(input->InterlaceCodingOption >= MB_CODING && mb_adaptive && img->field_mode)
  {
    mb_y = img->top_field ? 2*img->field_mb_y:2*img->field_mb_y+1;
    mb_ypos = img->field_mb_y;
    decr    = 2;
  }*/

  int mb_nr    = img->current_mb_nr;
  int mb_width = img->width/16;

  int mb_available_up   = (img->mb_y == 0 ) ? 0 : (img->mb_data[mb_nr].slice_nr == img->mb_data[mb_nr-mb_width  ].slice_nr);
  int mb_available_left = (img->mb_x == 0 ) ? 0 : (img->mb_data[mb_nr].slice_nr == img->mb_data[mb_nr-1         ].slice_nr);

  if (i==1 || i==3)
    Left_block= img->nz_coeff [img->mb_x ][mb_y ][i-1][j];
  else
    Left_block= mb_available_left ? img->nz_coeff [img->mb_x-1 ][img->mb_y ][i+1][j] : -1;
//    Left_block= img->mb_x > 0 ? img->nz_coeff [img->mb_x-1 ][img->mb_y ][i+1][j] : -1;

  if (j==5)
    Top_block=  img->nz_coeff [img->mb_x ][img->mb_y ][i][j-1];
  else
    Top_block=  mb_available_up ? img->nz_coeff [img->mb_x ][img->mb_y-decr ][i][5] : -1;
//    Top_block=  mb_ypos > 0 ? img->nz_coeff [img->mb_x ][img->mb_y-decr ][i][5] : -1;

//  if(input->InterlaceCodingOption >= MB_CODING && mb_adaptive)
//    Top_block=0;

  pred_nnz=0;
  if (Left_block>-1)
  {
    pred_nnz=Left_block;
    cnt++;
  }
  if (Top_block>-1)
  {
    pred_nnz+=Top_block;
    cnt++;
  }

  if (cnt==2)
    pred_nnz++;

  if (cnt)
    pred_nnz/=cnt; 
  return pred_nnz;
}




/*!
 ************************************************************************
 * \brief
 *    Writes coeff of an 4x4 block (CAVLC)
 *
 * \author
 *    Karl Lillevold <karll@real.com>
 *    contributions by James Au <james@ubvideo.com>
 ************************************************************************
 */

int writeCoeff4x4_CAVLC (int block_type, int b8, int b4, int param)
{
  int           no_bits    = 0;
  Macroblock    *currMB    = &img->mb_data[img->current_mb_nr];
  SyntaxElement *currSE    = &img->MB_SyntaxElements[currMB->currSEnr];
  int           *bitCount  = currMB->bitcounter;
  Slice         *currSlice = img->currentSlice;
  DataPartition *dataPart;
  int           *partMap   = assignSE2partition[input->partition_mode];

  int k,level,run,vlcnum;
  int numcoeff, lastcoeff, numtrailingones; 
  int numones, totzeros, zerosleft, numcoef;
  int numcoeff_vlc;
  int code, level_two_or_higher;
  int dptype = 0, bitcounttype = 0;
  int nnz, max_coeff_num = 0, cdc=0, cac=0;
  int subblock_x, subblock_y;
  char type[15];
  int  mb_y = img->mb_y;

  int incVlc[] = {0,3,6,12,24,48,32768};  // maximum vlc = 6


  int*  pLevel = NULL;
  int*  pRun = NULL;

  if(input->InterlaceCodingOption >= MB_CODING && mb_adaptive && img->field_mode)
    mb_y = img->top_field ? 2*img->field_mb_y:2*img->field_mb_y+1;

  switch (block_type)
  {
  case LUMA:
    max_coeff_num = 16;
    bitcounttype = BITS_COEFF_Y_MB;

    pLevel = img->cofAC[b8][b4][0];
    pRun   = img->cofAC[b8][b4][1];

    sprintf(type, "%s", "Luma");
    if (IS_INTRA (currMB))
    {
      dptype = SE_LUM_AC_INTRA;
    }
    else
    {
      dptype = SE_LUM_AC_INTER;
    }
    break;
  case LUMA_INTRA16x16DC:
    max_coeff_num = 16;
    bitcounttype = BITS_COEFF_Y_MB;

    pLevel = img->cofDC[0][0];
    pRun   = img->cofDC[0][1];

    sprintf(type, "%s", "Lum16DC");
    dptype = SE_LUM_DC_INTRA;
    break;
  case LUMA_INTRA16x16AC:
    max_coeff_num = 15;
    bitcounttype = BITS_COEFF_Y_MB;

    pLevel = img->cofAC[b8][b4][0];
    pRun   = img->cofAC[b8][b4][1];

    sprintf(type, "%s", "Lum16AC");
    dptype = SE_LUM_AC_INTRA;
    break;

  case CHROMA_DC:
    max_coeff_num = 4;
    bitcounttype = BITS_COEFF_UV_MB;
    cdc = 1;

    pLevel = img->cofDC[param+1][0];
    pRun   = img->cofDC[param+1][1];

    sprintf(type, "%s", "ChrDC");
    if (IS_INTRA (currMB))
    {
      dptype = SE_CHR_DC_INTRA;
    }
    else
    {
      dptype = SE_CHR_DC_INTER;
    }
    break;
  case CHROMA_AC:
    max_coeff_num = 15;
    bitcounttype = BITS_COEFF_UV_MB;
    cac = 1;

    pLevel = img->cofAC[b8][b4][0];
    pRun   = img->cofAC[b8][b4][1];

    sprintf(type, "%s", "ChrAC");
    if (IS_INTRA (currMB))
    {
      dptype = SE_CHR_AC_INTRA;
    }
    else
    {
      dptype = SE_CHR_AC_INTER;
    }
    break;
  default:
    error("writeCoeff4x4_CAVLC: Invalid block type", 600);
    break;
  }

  if (img->type == B_SLICE)
  {
    dptype = SE_BFRAME;
  }

  dataPart = &(currSlice->partArr[partMap[dptype]]);

  numcoeff = 0;
  numtrailingones = 0;
  numones = 0;
  lastcoeff = 0;
  totzeros = 0;
  level = 1;

  for(k = 0; (k <= cdc?4:16)&& level !=0; k++)
  {
    level = pLevel[k]; // level
    run   = pRun[k];   // run

    if (level)
    {
      if (run)
        totzeros += run;
      if (abs(level) == 1)
      {
        numtrailingones ++;
        numones ++;
        if (numtrailingones > 3)
        {
          numtrailingones = 3; /* clip to 3 */
        }
      }
      else
      {
        numtrailingones = 0;
      }
      numcoeff ++;
      lastcoeff = k;
    }
  }

  if (!cdc)
  {
    if (!cac)
    {
      // luma
      subblock_x = ((b8&0x1)==0)?(((b4&0x1)==0)?0:1):(((b4&0x1)==0)?2:3); 
        // horiz. position for coeff_count context
      subblock_y = (b8<2)?((b4<2)?0:1):((b4<2)?2:3); 
        // vert.  position for coeff_count context
      nnz = predict_nnz(subblock_x,subblock_y);
    }
    else
    {
      // chroma AC
      subblock_x = param >> 4;
      subblock_y = param & 15;
      nnz = predict_nnz_chroma(subblock_x,subblock_y);
    }

    img->nz_coeff [img->mb_x ][mb_y ][subblock_x][subblock_y] = numcoeff;


    if (nnz < 2)
    {
      numcoeff_vlc = 0;
    }
    else if (nnz < 4)
    {
      numcoeff_vlc = 1;
    }
    else if (nnz < 8)
    {
      numcoeff_vlc = 2;
    }
    else 
    {
      numcoeff_vlc = 3;
    }


  }
  else
  {
    // chroma DC (has its own VLC)
    // numcoeff_vlc not relevant
    numcoeff_vlc = 0;

    subblock_x = param;
    subblock_y = param;
  }

  currSE->type  = dptype;   

  currSE->value1 = numcoeff;
  currSE->value2 = numtrailingones;
  currSE->len = numcoeff_vlc; /* use len to pass vlcnum */

#if TRACE
  snprintf(currSE->tracestring, 
    TRACESTRING_SIZE, "%s # c & tr.1s(%d,%d) vlc=%d #c=%d #t1=%d",
    type, subblock_x, subblock_y, numcoeff_vlc, numcoeff, numtrailingones);
#endif

  if (!cdc)
    writeSyntaxElement_NumCoeffTrailingOnes(currSE, dataPart);
  else
    writeSyntaxElement_NumCoeffTrailingOnesChromaDC(currSE, dataPart);

  bitCount[bitcounttype]+=currSE->len;
  no_bits               +=currSE->len;

  // proceed to next SE
  currSE++;
  currMB->currSEnr++;


  if (!numcoeff)
    return no_bits;

  if (numcoeff)
  {
    code = 0;
    for (k = lastcoeff; k > lastcoeff-numtrailingones; k--)
    {
      level = pLevel[k]; // level
      if (abs(level) > 1)
      {
        printf("ERROR: level > 1\n");
        exit(-1);
      }
      code <<= 1;
      if (level < 0)
      {
        code |= 0x1;
      }
    }

    if (numtrailingones)
    {
      currSE->type  = dptype;   

      currSE->value2 = numtrailingones;
      currSE->value1 = code;

#if TRACE
      snprintf(currSE->tracestring, 
        TRACESTRING_SIZE, "%s trailing ones sign (%d,%d)", 
        type, subblock_x, subblock_y);
#endif

      writeSyntaxElement_VLC (currSE, dataPart);
      bitCount[bitcounttype]+=currSE->len;
      no_bits               +=currSE->len;

      // proceed to next SE
      currSE++;
      currMB->currSEnr++;
    }

    // encode levels
    level_two_or_higher = 1;
    if (numcoeff > 3 && numtrailingones == 3)
      level_two_or_higher = 0;

    if (numcoeff > 10 && numtrailingones < 3)
      vlcnum = 1;
    else
      vlcnum = 0;

    for (k = lastcoeff - numtrailingones; k >= 0; k--)
    {
      level = pLevel[k]; // level

      currSE->value1 = level;
      currSE->type  = dptype;   

  #if TRACE
        snprintf(currSE->tracestring, 
          TRACESTRING_SIZE, "%s lev (%d,%d) k=%d vlc=%d lev=%3d",
            type, subblock_x, subblock_y, k, vlcnum, level);
  #endif

          if (level_two_or_higher)
          {
            if (currSE->value1 > 0)
              currSE->value1 --;
            else
              currSE->value1 ++;
            level_two_or_higher = 0;
          }

      //    encode level
      if (vlcnum == 0)
        writeSyntaxElement_Level_VLC1(currSE, dataPart);
      else
        writeSyntaxElement_Level_VLCN(currSE, vlcnum, dataPart);

      // update VLC table
      if (abs(level)>incVlc[vlcnum])
        vlcnum++;

      if (k == lastcoeff - numtrailingones && abs(level)>3)
        vlcnum = 2;

      bitCount[bitcounttype]+=currSE->len;
      no_bits               +=currSE->len;

      // proceed to next SE
      currSE++;
      currMB->currSEnr++;
    }


    // encode total zeroes
    if (numcoeff < max_coeff_num)
    {

      currSE->type  = dptype;   
      currSE->value1 = totzeros;

      vlcnum = numcoeff-1;

      currSE->len = vlcnum;

#if TRACE
      snprintf(currSE->tracestring, 
        TRACESTRING_SIZE, "%s totalrun (%d,%d) vlc=%d totzeros=%3d",
          type, subblock_x, subblock_y, vlcnum, totzeros);
#endif
      if (!cdc)
        writeSyntaxElement_TotalZeros(currSE, dataPart);
      else
        writeSyntaxElement_TotalZerosChromaDC(currSE, dataPart);

      bitCount[bitcounttype]+=currSE->len;
      no_bits               +=currSE->len;

      // proceed to next SE
      currSE++;
      currMB->currSEnr++;
    }

    // encode run before each coefficient
    zerosleft = totzeros;
    numcoef = numcoeff;
    for (k = lastcoeff; k >= 0; k--)
    {
      run = pRun[k]; // run

      currSE->value1 = run;
      currSE->type  = dptype;   

      // for last coeff, run is remaining totzeros
      // when zerosleft is zero, remaining coeffs have 0 run
      if (numcoeff <= 1 || !zerosleft)
        break;

      if (numcoef > 1 && zerosleft) 
      {

        vlcnum = zerosleft - 1;
        if (vlcnum > RUNBEFORE_NUM-1)
          vlcnum = RUNBEFORE_NUM-1;

        currSE->len = vlcnum;

#if TRACE
        snprintf(currSE->tracestring, 
          TRACESTRING_SIZE, "%s run (%d,%d) k=%d vlc=%d run=%2d",
            type, subblock_x, subblock_y, k, vlcnum, run);
#endif

        writeSyntaxElement_Run(currSE, dataPart);

        bitCount[bitcounttype]+=currSE->len;
        no_bits               +=currSE->len;

        zerosleft -= run;
        numcoef --;

        // proceed to next SE
        currSE++;
        currMB->currSEnr++;
      }
    }
  }

  return no_bits;
}




/*!
 ************************************************************************
 * \brief
 *    Find best 16x16 based intra mode
 *
 * \par Input:
 *    Image parameters, pointer to best 16x16 intra mode
 *
 * \par Output:
 *    best 16x16 based SAD
 ************************************************************************/
int find_sad_16x16(int *intra_mode)
{
  int current_intra_sad_2,best_intra_sad2;
  int M1[16][16],M0[4][4][4][4],M3[4],M4[4][4];

  int i,j,k;
  int ii,jj;
  int incr = 1;
  int offset = 0; // For MB level field/frame coding  
  int mb_nr = img->current_mb_nr;
  
  PixelPos up;          //!< pixel position p(0,-1)
  PixelPos left[17];    //!< pixel positions p(-1, -1..15)

  int up_avail, left_avail, left_up_avail;
	
  for (i=0;i<17;i++)
  {
    getNeighbour(mb_nr, -1 ,  i-1 , 1, &left[i]);
  }
  
  getNeighbour(mb_nr, 0     ,  -1 , 1, &up);
	
  if (!(input->UseConstrainedIntraPred))
	{
    up_avail   = up.available;
    left_avail = left[1].available;
    left_up_avail = left[0].available;
	}
  else
  {
    up_avail      = up.available ? img->intra_block[up.mb_addr] : 0;
    for (i=1, left_avail=1; i<17;i++)
      left_avail  &= left[i].available ? img->intra_block[left[i].mb_addr]: 0;
    left_up_avail = left[0].available ? img->intra_block[left[0].mb_addr]: 0;
  }
	
  if (img->MbaffFrameFlag && img->field_mode)
  {
    incr   = 2;
    if(mb_nr%2)
    {
      offset = -15;
    }
  }

  best_intra_sad2=MAX_VALUE;

  for (k=0;k<4;k++)
  {
    //check if there are neighbours to predict from
    if ((k==0 && !up_avail) || (k==1 && !left_avail) || (k==3 && (!left_avail || !up_avail || !left_up_avail)))
    {
      ; // edge, do nothing
    }
    else
    {
      for (j=0;j<16;j++)
      {
        for (i=0;i<16;i++)
        {
          M1[i][j]=imgY_org[img->pix_y+(j*incr)+offset][img->pix_x+i]-img->mprr_2[k][j][i];
          M0[i%4][i/4][j%4][j/4]=M1[i][j];
        }
      }
      current_intra_sad_2=0;              // no SAD start handicap here
      for (jj=0;jj<4;jj++)
      {
        for (ii=0;ii<4;ii++)
        {
          for (j=0;j<4;j++)
          {
            M3[0]=M0[0][ii][j][jj]+M0[3][ii][j][jj];
            M3[1]=M0[1][ii][j][jj]+M0[2][ii][j][jj];
            M3[2]=M0[1][ii][j][jj]-M0[2][ii][j][jj];
            M3[3]=M0[0][ii][j][jj]-M0[3][ii][j][jj];

            M0[0][ii][j][jj]=M3[0]+M3[1];
            M0[2][ii][j][jj]=M3[0]-M3[1];
            M0[1][ii][j][jj]=M3[2]+M3[3];
            M0[3][ii][j][jj]=M3[3]-M3[2];
          }

          for (i=0;i<4;i++)
          {
            M3[0]=M0[i][ii][0][jj]+M0[i][ii][3][jj];
            M3[1]=M0[i][ii][1][jj]+M0[i][ii][2][jj];
            M3[2]=M0[i][ii][1][jj]-M0[i][ii][2][jj];
            M3[3]=M0[i][ii][0][jj]-M0[i][ii][3][jj];

            M0[i][ii][0][jj]=M3[0]+M3[1];
            M0[i][ii][2][jj]=M3[0]-M3[1];
            M0[i][ii][1][jj]=M3[2]+M3[3];
            M0[i][ii][3][jj]=M3[3]-M3[2];
            for (j=0;j<4;j++)
              if ((i+j)!=0)
                current_intra_sad_2 += abs(M0[i][ii][j][jj]);
          }
        }
      }

      for (j=0;j<4;j++)
        for (i=0;i<4;i++)
          M4[i][j]=M0[0][i][0][j]/4;

        // Hadamard of DC koeff
        for (j=0;j<4;j++)
        {
          M3[0]=M4[0][j]+M4[3][j];
          M3[1]=M4[1][j]+M4[2][j];
          M3[2]=M4[1][j]-M4[2][j];
          M3[3]=M4[0][j]-M4[3][j];

          M4[0][j]=M3[0]+M3[1];
          M4[2][j]=M3[0]-M3[1];
          M4[1][j]=M3[2]+M3[3];
          M4[3][j]=M3[3]-M3[2];
        }

        for (i=0;i<4;i++)
        {
          M3[0]=M4[i][0]+M4[i][3];
          M3[1]=M4[i][1]+M4[i][2];
          M3[2]=M4[i][1]-M4[i][2];
          M3[3]=M4[i][0]-M4[i][3];

          M4[i][0]=M3[0]+M3[1];
          M4[i][2]=M3[0]-M3[1];
          M4[i][1]=M3[2]+M3[3];
          M4[i][3]=M3[3]-M3[2];

          for (j=0;j<4;j++)
            current_intra_sad_2 += abs(M4[i][j]);
        }
        if(current_intra_sad_2 < best_intra_sad2)
        {
          best_intra_sad2=current_intra_sad_2;
          *intra_mode = k; // update best intra mode

        }
    }
  }
  best_intra_sad2 = best_intra_sad2/2;

  return best_intra_sad2;

}
