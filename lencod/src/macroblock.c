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

#include "elements.h"
#include "macroblock.h"
#include "refbuf.h"


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
  const int number_mb_per_row = img->width / MB_BLOCK_SIZE ;
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

/*  if (input->symbol_mode == UVLC)
    stat->bit_ctr += currMB->bitcounter[BITS_TOTAL_MB]; */
  if (img->type==INTRA_IMG)
    ++stat->mode_use_intra[currMB->mb_type];
  else
    if (img->type != B_IMG)
    {
      ++stat->mode_use_inter[0][currMB->mb_type];
      stat->bit_use_mode_inter[0][currMB->mb_type]+= bitCount[BITS_INTER_MB];

    }
    else
    {
      stat->bit_use_mode_inter[1][currMB->mb_type]+= bitCount[BITS_INTER_MB];
      ++stat->mode_use_inter[1][currMB->mb_type];
    }


  // Update coordinates of macroblock
  img->mb_x++;
  if (img->mb_x == number_mb_per_row) // next row of MBs
  {
    img->mb_x = 0; // start processing of next row
    img->mb_y++;
  }
  img->current_mb_nr++;


  // Define vertical positions
  img->block_y = img->mb_y * BLOCK_SIZE;      // vertical luma block position
  img->pix_y   = img->mb_y * MB_BLOCK_SIZE;   // vertical luma macroblock position
  img->pix_c_y = img->mb_y * MB_BLOCK_SIZE/2; // vertical chroma macroblock position

  // Define horizontal positions
  img->block_x = img->mb_x * BLOCK_SIZE;        // luma block
  img->pix_x   = img->mb_x * MB_BLOCK_SIZE;     // luma pixel
  img->block_c_x = img->mb_x * BLOCK_SIZE/2;    // chroma block
  img->pix_c_x   = img->mb_x * MB_BLOCK_SIZE/2; // chroma pixel

  // Statistics
  if ((img->type == INTER_IMG)||(img->types==SP_IMG) )
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
  int x=0, y=0 ;
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
  if (img->type != B_IMG)
  {
    for (k=0; k < 2; k++)
    {
      for (j=0; j < BLOCK_MULTIPLE; j++)
        for (i=0; i < BLOCK_MULTIPLE; i++)
          tmp_mv[k][img->block_y+j][img->block_x+i+4]=0;
    }
  }

  // Reset syntax element entries in MB struct
  currMB->mb_type   = 0;
  currMB->cbp_blk   = 0;
  currMB->cbp       = 0;

  for (l=0; l < 2; l++)
    for (j=0; j < BLOCK_MULTIPLE; j++)
      for (i=0; i < BLOCK_MULTIPLE; i++)
        for (k=0; k < 2; k++)
          currMB->mvd[l][j][i][k] = 0;
 
  for (j=0; j < BLOCK_MULTIPLE; j++)
    for (i=0; i < BLOCK_MULTIPLE; i++)
      currMB->coeffs_count[j][i] = 0;

  for (i=0; i < (BLOCK_MULTIPLE*BLOCK_MULTIPLE); i++)
    currMB->intra_pred_modes[i] = 0;

  if (input->UseConstrainedIntraPred)
  {
    i = img->current_mb_nr;
    img->intra_block[i][0] =img->intra_block[i][1] = img->intra_block[i][2] = img->intra_block[i][3] = 1;
 }

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
  int i,x=0, y=0 ;
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
    *recode_macroblock = FALSE;
    if ((img->current_mb_nr+1) == img->total_number_mb) // maximum number of MBs
      *end_of_slice = TRUE;
    break;
  case FIXED_MB:
    // For slice mode one, check if a new slice boundary follows
    *recode_macroblock = FALSE;
    if ( ((img->current_mb_nr+1) % input->slice_argument == 0) || ((img->current_mb_nr+1) == img->total_number_mb) )
    {
      *end_of_slice = TRUE;
    }
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
       if(img->type == B_IMG)
         dataPart = &(currSlice->partArr[partMap[SE_BFRAME]]);
       else
         dataPart = &(currSlice->partArr[partMap[SE_MBTYPE]]);
       currSE->value1 = img->cod_counter;
       currSE->mapping = n_linfo2;
       currSE->type = SE_MBTYPE;
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
        if(img->type == B_IMG)
          dataPart = &(currSlice->partArr[partMap[SE_BFRAME]]);
        else
          dataPart = &(currSlice->partArr[partMap[SE_MBTYPE]]);
        currSE->value1 = img->cod_counter;
        currSE->mapping = n_linfo2;
        currSE->type = SE_MBTYPE;
        dataPart->writeSyntaxElement(  currSE, dataPart);
        rlc_bits=currSE->len;
        currMB->bitcounter[BITS_MB_MODE]+=rlc_bits;
        img->cod_counter = 0;
      }
    }
    else //! MB that did not fit in this slice anymore is not a Skip MB
    {
      for (i=0; i<currSlice->max_part_nr; i++)
      {
        dataPart = &(currSlice->partArr[i]);
        currStream = dataPart->bitstream;
        // update the bitstream
        currStream->bits_to_go = currStream->bits_to_go_skip;
        currStream->byte_pos  = currStream->byte_pos_skip;
        currStream->byte_buf  = currStream->byte_buf_skip;
      }
      // update the statistics
      img->cod_counter = 0;
      skip = FALSE;
    }
  }
  
  //! TO 4.11.2001 Skip MBs at the end of this slice for Slice Mode 0 or 1
  if(*end_of_slice == TRUE && img->cod_counter && !use_bitstream_backing)
  {
    if(img->type == B_IMG)
      dataPart = &(currSlice->partArr[partMap[SE_BFRAME]]);
    else
      dataPart = &(currSlice->partArr[partMap[SE_MBTYPE]]);
    currSE->value1 = img->cod_counter;
    currSE->mapping = n_linfo2;
    currSE->type = SE_MBTYPE;
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
       size_in_bytes = currStream->byte_pos;
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
void CheckAvailabilityOfNeighbors()
{
  int i,j;
  const int mb_width = img->width/MB_BLOCK_SIZE;
  const int mb_nr = img->current_mb_nr;
  Macroblock *currMB = &img->mb_data[mb_nr];

  // mark all neighbors as unavailable
  for (i=0; i<3; i++)
    for (j=0; j<3; j++)
      img->mb_data[mb_nr].mb_available[i][j]=NULL;
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
    }
    // lower blocks
    if (remove_prediction || (input->UseConstrainedIntraPred && img->intra_block[mb_nr-1][3]==0))
    {
      img->ipredmode[img->block_x][img->block_y+3] = -1;
      img->ipredmode[img->block_x][img->block_y+4] = -1;
    }
    if (!remove_prediction)
    {
      currMB->mb_available[1][0]=&(img->mb_data[mb_nr-1]);
    }
  }


  // Check MB above
  if(img->pix_y >= MB_BLOCK_SIZE)
  {
    int remove_prediction = currMB->slice_nr != img->mb_data[mb_nr-mb_width].slice_nr;
    // upper blocks
    if (remove_prediction || (input->UseConstrainedIntraPred && img->intra_block[mb_nr-mb_width][2]==0))
    {
      img->ipredmode[img->block_x+1][img->block_y] = -1;
      img->ipredmode[img->block_x+2][img->block_y] = -1;
    }
    // lower blocks
    if (remove_prediction || (input->UseConstrainedIntraPred && img->intra_block[mb_nr-mb_width][3]==0))
    {
      img->ipredmode[img->block_x+3][img->block_y] = -1;
      img->ipredmode[img->block_x+4][img->block_y] = -1;
    }
    if (!remove_prediction)
    {
      currMB->mb_available[0][1]=&(img->mb_data[mb_nr-mb_width]);
    }
  }

  // Check MB left above
  if(img->pix_x >= MB_BLOCK_SIZE && img->pix_y  >= MB_BLOCK_SIZE )
  {
    if(currMB->slice_nr == img->mb_data[mb_nr-mb_width-1].slice_nr)
      img->mb_data[mb_nr].mb_available[0][0]=&(img->mb_data[mb_nr-mb_width-1]);
  }

  // Check MB right above
  if(img->pix_y >= MB_BLOCK_SIZE && img->pix_x < (img->width-MB_BLOCK_SIZE ))
  {
    if(currMB->slice_nr == img->mb_data[mb_nr-mb_width+1].slice_nr)
      // currMB->mb_available[0][1]=&(img->mb_data[mb_nr-mb_width+1]);
      currMB->mb_available[0][2]=&(img->mb_data[mb_nr-mb_width+1]);
  }
}





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
  pel_t** ref_pic = (ref==-1 ? mref_P : mref [ref]);
  int     mvshift = (input->mv_res ? 3 : 2);
  int     pix_add = (1 << mvshift);
  int     j0      = (pic_pix_y << mvshift) + mv[1], j1=j0+pix_add, j2=j1+pix_add, j3=j2+pix_add;
  int     i0      = (pic_pix_x << mvshift) + mv[0], i1=i0+pix_add, i2=i1+pix_add, i3=i2+pix_add;

  pel_t (*get_pel) (pel_t**, int, int) = (input->mv_res ? UMVPelY_18 : UMVPelY_14);

  *mpred++ = get_pel (ref_pic, j0, i0);
  *mpred++ = get_pel (ref_pic, j0, i1);
  *mpred++ = get_pel (ref_pic, j0, i2);
  *mpred++ = get_pel (ref_pic, j0, i3);
  *mpred++ = get_pel (ref_pic, j1, i0);
  *mpred++ = get_pel (ref_pic, j1, i1);
  *mpred++ = get_pel (ref_pic, j1, i2);
  *mpred++ = get_pel (ref_pic, j1, i3);
  *mpred++ = get_pel (ref_pic, j2, i0);
  *mpred++ = get_pel (ref_pic, j2, i1);
  *mpred++ = get_pel (ref_pic, j2, i2);
  *mpred++ = get_pel (ref_pic, j2, i3);
  *mpred++ = get_pel (ref_pic, j3, i0);
  *mpred++ = get_pel (ref_pic, j3, i1);
  *mpred++ = get_pel (ref_pic, j3, i2);
  *mpred++ = get_pel (ref_pic, j3, i3);
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
                   int  fw_ref)     // <--  reference frame for forward prediction (-1: Intra4x4 pred. with fw_mode)
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
  int  direct    = (fw_mode == 0 && bw_mode == 0);


  if (fw_mode || direct)
  {
    OneComponentLumaPrediction4x4 (fw_pred, pic_pix_x, pic_pix_y, img->all_mv [bx][by][fw_ref][fw_mode], fw_ref);

  }
  if (bw_mode || direct)
  {
    OneComponentLumaPrediction4x4 (bw_pred, pic_pix_x, pic_pix_y, img->all_bmv[bx][by][     0][bw_mode],     -1);
  }

  if (direct || (fw_mode && bw_mode))
  {
    for   (j=block_y; j<block_y4; j++)
      for (i=block_x; i<block_x4; i++)  img->mpr[i][j] = (int)((*fpred++ + *bpred++) / 2.0 + 0.5);
  }
  else if (fw_mode)
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
                       int  fw_refframe)  // <--  reference frame for forward prediction
{
  int    block_y, block_x, pic_pix_y, pic_pix_x, i, j, nonzero, cbp_blk_mask;
  int    coeff_cost = 0;
  int    mb_y       = (block8x8 / 2) << 3;
  int    mb_x       = (block8x8 % 2) << 3;
  int    cbp_mask   = 1 << block8x8;


  //===== loop over 4x4 blocks =====
  for (block_y=mb_y; block_y<mb_y+8; block_y+=4)
  {
    pic_pix_y = img->pix_y + block_y;

    for (block_x=mb_x; block_x<mb_x+8; block_x+=4)
    {
      pic_pix_x = img->pix_x + block_x;

      cbp_blk_mask = (block_x>>2) + block_y;

      //===== prediction of 4x4 block =====
      LumaPrediction4x4 (block_x, block_y, fw_mode, bw_mode, fw_refframe);

      //===== get displaced frame difference ======                
      for (j=0; j<4; j++)
      for (i=0; i<4; i++)
      {
        img->m7[i][j] = imgY_org[pic_pix_y+j][pic_pix_x+i] - img->mpr[i+block_x][j+block_y];
      }

      //===== DCT, Quantization, inverse Quantization, IDCT, Reconstruction =====
      if (img->types!=SP_IMG)  nonzero = dct_luma   (block_x, block_y, &coeff_cost, 0);
      else                     nonzero = dct_luma_sp(block_x, block_y, &coeff_cost);
      if (nonzero)
      {
        (*cbp_blk) |= 1 << cbp_blk_mask;  // one bit for every 4x4 block
        (*cbp)     |= cbp_mask;           // one bit for the 4x4 blocks of an 8x8 block
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

  if (coeff_cost <= 3)
  {
    coeff_cost  = 0;
    (*cbp)     &=  (63 - cbp_mask);
    (*cbp_blk) &= ~(51 << (4*block8x8-2*(block8x8%2)));

    for (i=mb_x; i<mb_x+8; i++)
    for (j=mb_y; j<mb_y+8; j++)
    {
      imgY[img->pix_y+j][img->pix_x+i] = img->mpr[i][j];
    }
    if (img->types==SP_IMG)
    {
      for (i=mb_x; i < mb_x+BLOCK_SIZE*2; i+=BLOCK_SIZE)
      for (j=mb_y; j < mb_y+BLOCK_SIZE*2; j+=BLOCK_SIZE)
      {
        copyblock_sp(i,j);
      }
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
SetModesAndRefframe (int b8, int* fw_mode, int* bw_mode, int* refframe)
{
  Macroblock* currMB = &img->mb_data[img->current_mb_nr];
  int         j      = 2*(b8/2);
  int         i      = 2*(b8%2);

  if (img->type!=B_IMG)
  {
    *refframe = refFrArr[img->block_y+j][img->block_x+i];
    *bw_mode  = 0;
    *fw_mode  = currMB->b8mode[b8];
  }
  else
  {
    if (currMB->b8pdir[b8]==-1)
    {
      *refframe = -1;
      *fw_mode  =  0;
      *bw_mode  =  0;
    }
    else if (currMB->b8pdir[b8]==0)
    {
      *refframe = fw_refFrArr[img->block_y+j][img->block_x+i];
      *fw_mode  = currMB->b8mode[b8];
      *bw_mode  = 0;
    }
    else if (currMB->b8pdir[b8]==1)
    {
      *refframe = 0;
      *fw_mode  = 0;
      *bw_mode  = currMB->b8mode[b8];
    }
    else
    {
      *refframe = fw_refFrArr[img->block_y+j][img->block_x+i];
      *fw_mode  = currMB->b8mode[b8];
      *bw_mode  = currMB->b8mode[b8];
      if (currMB->b8mode[b8]==0) // direct
      {
        *refframe = max(0,refFrArr[img->block_y+j][img->block_x+i]);
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
  int i,j,block8x8;
  int fw_mode, bw_mode, refframe;
  int sum_cnt_nonz;
  Macroblock *currMB = &img->mb_data[img->current_mb_nr];

  currMB->cbp     = 0 ;
  currMB->cbp_blk = 0 ;
  sum_cnt_nonz    = 0 ;

  for (block8x8=0; block8x8<4; block8x8++)
  {
    SetModesAndRefframe (block8x8, &fw_mode, &bw_mode, &refframe);

    sum_cnt_nonz += LumaResidualCoding8x8 (&(currMB->cbp), &(currMB->cbp_blk), block8x8,
                                           fw_mode, bw_mode, refframe);
  }
  
  if (sum_cnt_nonz <= 5 )
  {
     currMB->cbp     &= 0xfffff0 ;
     currMB->cbp_blk &= 0xff0000 ;
     for (i=0; i < MB_BLOCK_SIZE; i++)
     {
       for (j=0; j < MB_BLOCK_SIZE; j++)
       {
         imgY[img->pix_y+j][img->pix_x+i]=img->mpr[i][j];
       }
     }
     if (img->types==SP_IMG)
     {
       for (i=0; i < MB_BLOCK_SIZE; i+=BLOCK_SIZE)
       for (j=0; j < MB_BLOCK_SIZE; j+=BLOCK_SIZE)
       {
         copyblock_sp(i,j);
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
  int*    mvb;
  int     refframe  = (ref<0 ?      0 :      ref);
  pel_t** refimage  = (ref<0 ? mcef_P : mcef[ref])[uv];
  int     je        = pix_c_y + 4;
  int     ie        = pix_c_x + 4;
  int     f1        =(input->mv_res?16:8), f2=f1-1, f3=f1*f1, f4=f3>>1;
  int     s1        =(input->mv_res? 4:3);


  for (j=pix_c_y; j<je; j++)
  for (i=pix_c_x; i<ie; i++)
  {
    mvb  = mv [(i-img->pix_c_x)>>1][(j-img->pix_c_y)>>1][refframe][blocktype];
    ii   = (i<<s1) + mvb[0];
    jj   = (j<<s1) + mvb[1];

    ii0  = max (0, min (img->width_cr -1, ii>>s1     ));
    jj0  = max (0, min (img->height_cr-1, jj>>s1     ));
    ii1  = max (0, min (img->width_cr -1, (ii+f2)>>s1));
    jj1  = max (0, min (img->height_cr-1, (jj+f2)>>s1));

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
void
IntraChromaPrediction4x4 (int  uv,       // <-- colour component
                          int  block_x,  // <-- relative horizontal block coordinate of 4x4 block
                          int  block_y)  // <-- relative vertical   block coordinate of 4x4 block
{
  int     s=128, s0=0, s1=0, s2=0, s3=0, i, j;
  pel_t** image             = imgUV[uv];
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

  if(input->UseConstrainedIntraPred)
  {
    if (mb_available_up   && (img->intra_block[mb_nr-mb_width][2]==0 || img->intra_block[mb_nr-mb_width][3]==0))
      mb_available_up   = 0;
    if (mb_available_left && (img->intra_block[mb_nr-       1][1]==0 || img->intra_block[mb_nr       -1][3]==0))
      mb_available_left = 0;
  }

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
    img->mpr[i][j] = s;
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
                     int  fw_ref_frame) // <-- reference frame for forward prediction (if (<0) -> intra prediction)
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
  int  direct    = (fw_mode == 0 && bw_mode == 0);

  //===== INTRA PREDICTION =====
  if (fw_ref_frame < 0)
  {
    IntraChromaPrediction4x4 (uv, block_x, block_y);
    return;
  }

  //===== INTER PREDICTION =====
  if (fw_mode || direct)
  {
    OneComponentChromaPrediction4x4 (fw_pred, pic_pix_x, pic_pix_y, img->all_mv , fw_ref_frame, fw_mode, uv);
  }
  if (bw_mode || direct)
  {
    OneComponentChromaPrediction4x4 (bw_pred, pic_pix_x, pic_pix_y, img->all_bmv,           -1, bw_mode, uv);
  }

  if (direct || (fw_mode && bw_mode))
  {
    for (j=block_y; j<block_y4; j++)
    for (i=block_x; i<block_x4; i++)  img->mpr[i][j] = (int)((*fpred++ + *bpred++) / 2.0 + 0.5);
  }
  else if (fw_mode)
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

  for (*cr_cbp=0, uv=0; uv<2; uv++)
  {
    //===== prediction of chrominance blocks ===d==
    block8 = 0;
    for (block_y=0; block_y<8; block_y+=4)
    for (block_x=0; block_x<8; block_x+=4, block8++)
    {
      SetModesAndRefframe (block8, &fw_mode, &bw_mode, &refframe);

      ChromaPrediction4x4 (uv, block_x, block_y, fw_mode, bw_mode, refframe);
    }

    //===== calculation of displaced frame difference =====
    for (j=0; j<8; j++)
    for (i=0; i<8; i++)
    {
      img->m7[i][j] = imgUV_org[uv][img->pix_c_y+j][img->pix_c_x+i] - img->mpr[i][j];
    }

    //===== DCT, Quantization, inverse Quantization, IDCT, and Reconstruction =====
    if (img->types!=SP_IMG || IS_INTRA (&img->mb_data[img->current_mb_nr]))
      *cr_cbp=dct_chroma   (uv,*cr_cbp);
    else
      *cr_cbp=dct_chroma_sp(uv,*cr_cbp);
  }

  //===== update currMB->cbp =====
  img->mb_data[img->current_mb_nr].cbp += ((*cr_cbp)<<4);  
}



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

  if (img->type!=B_IMG)
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
  for (j=0; j<4; j++)
  for (i=0; i<4; i++)
  {
    if (currMB->b8mode[2*(j/2)+(i/2)]!=IBLOCK && refFrArr[img->block_y+j][img->block_x+i]!=0)
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

  if (img->type!=B_IMG)
  {
    if      (currMB->mb_type==I4MB)     return (img->type==INTRA_IMG ? 0 : 6);
    else if (currMB->mb_type==I16MB)    return (img->type==INTRA_IMG ? 0 : 6) + img->i16offset;
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
int
writeIntra4x4Modes (int block8x8)
{
  int i, i0=4*block8x8, i1=i0+4;
  int           rate        = 0;
  Macroblock    *currMB     = &img->mb_data[img->current_mb_nr];
  SyntaxElement *currSE    = &img->MB_SyntaxElements[currMB->currSEnr];
  int           *bitCount   = currMB->bitcounter;
  Slice         *currSlice  = img->currentSlice;
  DataPartition *dataPart;
  const int     *partMap    = assignSE2partition[input->partition_mode];

  for (i=i0; i<i1; )
  {
    currSE->context = i;
    currSE->value1  = currMB->intra_pred_modes[i++];
    currSE->value2  = currMB->intra_pred_modes[i++];

#if TRACE
    snprintf(currSE->tracestring, TRACESTRING_SIZE, "Intra mode     = %3d",IPRED_ORDER[currSE->value1][currSE->value2]);
#endif

    /*--- set symbol type and function pointers ---*/
    if (input->symbol_mode == UVLC)    currSE->mapping = intrapred_linfo;
    else                               currSE->writing = writeIntraPredMode2Buffer_CABAC;
    currSE->type = SE_INTRAPREDMODE;
      
    /*--- choose data partition ---*/
    if (img->type==B_IMG)   dataPart = &(currSlice->partArr[partMap[SE_BFRAME]]);
    else                    dataPart = &(currSlice->partArr[partMap[SE_INTRAPREDMODE]]);
    
    /*--- encode and update rate ---*/
    dataPart->writeSyntaxElement (currSE, dataPart);
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

  if (img->type!=B_IMG)
  {
    return (b8mode==IBLOCK?4:b8mode-4);
  }
  else
  {
    if (b8mode==IBLOCK)  return 13;
    else                 return b8start[b8mode] + b8inc[b8mode] * b8pdir;
  }
}



/*!
 ************************************************************************
 * \brief
 *    Codes macroblock header
 ************************************************************************
 */
int
writeMBHeader ()
{
  int             i;
  int             mb_nr     = img->current_mb_nr;
  Macroblock*     currMB    = &img->mb_data[mb_nr];
  SyntaxElement *currSE    = &img->MB_SyntaxElements[currMB->currSEnr];
  int*            bitCount  = currMB->bitcounter;
  Slice*          currSlice = img->currentSlice;
  DataPartition*  dataPart;
  const int*      partMap   = assignSE2partition[input->partition_mode];
  int             no_bits   = 0;


  // choose the appropriate data partition
  if (img->type == B_IMG)   dataPart = &(currSlice->partArr[partMap[SE_BFRAME]]);
  else                      dataPart = &(currSlice->partArr[partMap[SE_MBTYPE]]);


  //=====  BITS FOR MACROBLOCK MODE =====
  if(img->type == INTRA_IMG || img->types == SP_IMG || input->symbol_mode == CABAC)
  {
    //===== No Run Length Coding of SKIP modes =====
    currSE->value1  = MBType2Value (currMB);
    currSE->type    = SE_MBTYPE;

    if (input->symbol_mode == UVLC)  currSE->mapping = n_linfo2;
    else                             currSE->writing = writeMB_typeInfo2Buffer_CABAC;

    dataPart->writeSyntaxElement( currSE, dataPart);
#if TRACE
    if (img->type == B_IMG)  snprintf(currSE->tracestring, TRACESTRING_SIZE, "B_MB mode(%2d,%2d) = %3d",img->mb_x, img->mb_y, currMB->mb_type);
    else                     snprintf(currSE->tracestring, TRACESTRING_SIZE,   "MB mode(%2d,%2d) = %3d",img->mb_x, img->mb_y,currMB->mb_type);
#endif
    bitCount[BITS_MB_MODE] += currSE->len;
    no_bits                += currSE->len;
    currSE++;
    currMB->currSEnr++;
  }
  else if (currMB->mb_type != 0 || (img->type == B_IMG && currMB->cbp != 0))
  {
    //===== Run Length Coding: Non-Skipped macorblock =====
    currSE->value1  = img->cod_counter;
    currSE->mapping = n_linfo2;
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

    // Put out mb mode
    currSE->value1  = MBType2Value (currMB);
    if(img->type != B_IMG)
    {
      currSE->value1--;
    }
    currSE->mapping = n_linfo2;
    currSE->type    = SE_MBTYPE;

    dataPart->writeSyntaxElement( currSE, dataPart);
#if TRACE
    if (img->type == B_IMG)   snprintf(currSE->tracestring, TRACESTRING_SIZE, "B_MB mode(%2d,%2d) = %3d",img->mb_x, img->mb_y, currMB->mb_type);
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

    if(img->current_mb_nr == img->total_number_mb)
    {
      // Put out run
      currSE->value1  = img->cod_counter;
      currSE->mapping = n_linfo2;
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


  //===== BITS FOR 8x8 SUB-PARTITION MODI =====
  if (IS_P8x8 (currMB))
  {
    if (img->type==B_IMG) dataPart = &(currSlice->partArr[partMap[SE_BFRAME]]);
    else                  dataPart = &(currSlice->partArr[partMap[SE_MBTYPE]]);

    for (i=0; i<4; i++)
    {
      if (input->symbol_mode==UVLC)   currSE->mapping = n_linfo2;
      else                            currSE->writing = writeB8_typeInfo2Buffer_CABAC;

      currSE->value1  = B8Mode2Value (currMB->b8mode[i], currMB->b8pdir[i]);
      currSE->type    = SE_MBTYPE;
      dataPart->writeSyntaxElement (currSE, dataPart);
      bitCount[BITS_MB_MODE]+= currSE->len;
      no_bits               += currSE->len;
      currSE++;
      currMB->currSEnr++;
    }
  }


  //===== BITS FOR INTRA PREDICTION MODI ====
  for (i=0; i<4; i++)
  {
    if (currMB->b8mode[i]==IBLOCK)
    {
      no_bits += writeIntra4x4Modes (i);
    }
  }

  return no_bits;
}



/*!
 ************************************************************************
 * \brief
 *    Passes the chosen syntax elements to the NAL
 ************************************************************************
 */
void write_one_macroblock ()
{
  Macroblock* currMB   = &img->mb_data[img->current_mb_nr];
  int*        bitCount = currMB->bitcounter;

  //--- write header ---
  writeMBHeader ();

  //  Do nothing more if copy and inter mode
  if ((IS_INTERMV (currMB)  || IS_INTRA (currMB)                             ) ||
      (img->type==B_IMG     && input->symbol_mode== CABAC                    ) ||
      (img->type==B_IMG     && input->symbol_mode== UVLC &&  currMB->cbp != 0)  )
  {
    writeMotionInfo2NAL  ();
    writeCBPandLumaCoeff ();
    writeChromaCoeff     ();
  }

  //--- constrain intra prediction ---
  if(input->UseConstrainedIntraPred && img->type==INTER_IMG && img->types != SP_IMG)
  {
    if (!IS_NEWINTRA (currMB) && currMB->b8mode[0]!=IBLOCK) img->intra_block[img->current_mb_nr][0] = 0;
    if (!IS_NEWINTRA (currMB) && currMB->b8mode[1]!=IBLOCK) img->intra_block[img->current_mb_nr][1] = 0;
    if (!IS_NEWINTRA (currMB) && currMB->b8mode[2]!=IBLOCK) img->intra_block[img->current_mb_nr][2] = 0;
    if (!IS_NEWINTRA (currMB) && currMB->b8mode[3]!=IBLOCK) img->intra_block[img->current_mb_nr][3] = 0;
  }

  //--- set total bit-counter ---
  bitCount[BITS_TOTAL_MB] = bitCount[BITS_MB_MODE] + bitCount[BITS_COEFF_Y_MB]     + bitCount[BITS_INTER_MB]
                          + bitCount[BITS_CBP_MB]  + bitCount[BITS_DELTA_QUANT_MB] + bitCount[BITS_COEFF_UV_MB];
  stat->bit_slice += bitCount[BITS_TOTAL_MB];
}




/*!
 ************************************************************************
 * \brief
 *    Sets context for reference frame parameter
 ************************************************************************
 */
int
BType2CtxRef (int btype)
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
int
writeReferenceFrame (int  mode,
                     int  i,
                     int  j,
                     int  ref)
{
  Macroblock*     currMB    = &img->mb_data[img->current_mb_nr];
  SyntaxElement*  currSE    = &img->MB_SyntaxElements[currMB->currSEnr];
  Slice*          currSlice = img->currentSlice;
  int*            bitCount  = currMB->bitcounter;
  const int*      partMap   = assignSE2partition[input->partition_mode];
  int             rate      = 0;
  DataPartition*  dataPart;


  currSE->value1 = ref;
  currSE->type   = (img->type==B_IMG ? SE_BFRAME : SE_REFFRAME);
  if (input->symbol_mode == UVLC)
  {
    currSE->mapping = n_linfo2;
  }
  else
  {
    currSE->context = BType2CtxRef (mode);
    img->subblock_x = i; // position used for context determination
    img->subblock_y = j; // position used for context determination
    currSE->writing = writeRefFrame2Buffer_CABAC;
  }

  dataPart = &(currSlice->partArr[partMap[currSE->type]]);
  dataPart->writeSyntaxElement (currSE, dataPart);
  bitCount[BITS_INTER_MB] += currSE->len;
  rate                    += currSE->len;
#if TRACE
  snprintf(currSE->tracestring, TRACESTRING_SIZE, "Reference frame no %d", currSE->value1);
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
int
writeMotionVector8x8 (int  i0,
                      int  j0,
                      int  i1,
                      int  j1,
                      int  refframe,
                      int  mv_mode)
{
  int            i, j, k, l, m;
  int            curr_mvd;
  DataPartition* dataPart;
  int            bwflag     = (refframe<0?1:0);
  int            rate       = 0;
  int            step_h     = input->blc_size[mv_mode][0] >> 2;
  int            step_v     = input->blc_size[mv_mode][1] >> 2;
  Macroblock*    currMB     = &img->mb_data[img->current_mb_nr];
  SyntaxElement* currSE     = &img->MB_SyntaxElements[currMB->currSEnr];
  Slice*         currSlice  = img->currentSlice;
  int*           bitCount   = currMB->bitcounter;
  const int*     partMap    = assignSE2partition[input->partition_mode];
  int            refindex   = (refframe<0 ? 0 : refframe);
  int*****       all_mv     = (refframe<0 ? img->all_bmv : img->all_mv);
  int*****       pred_mv    = (img->type!=B_IMG ? img->mv : refframe<0 ? img->p_bwMV : img->p_fwMV);

  for (j=j0; j<j1; j+=step_v)
  for (i=i0; i<i1; i+=step_h)
  {
    for (k=0; k<2; k++) 
    {
      curr_mvd = all_mv[i][j][refindex][mv_mode][k] - pred_mv[i][j][refindex][mv_mode][k];

      //--- store (oversampled) mvd ---
      for (l=0; l < step_v; l++) 
      for (m=0; m < step_h; m++)    currMB->mvd[bwflag][j+l][i+m][k] = curr_mvd;

      currSE->value1 = curr_mvd;
      currSE->type   = (img->type==B_IMG ? SE_BFRAME : SE_MVD);
      if (input->symbol_mode == UVLC)
      {
        currSE->mapping = mvd_linfo2;
      }
      else
      {
        img->subblock_x = i; // position used for context determination
        img->subblock_y = j; // position used for context determination
        if (img->type!=B_IMG)
        {
          currSE->value2  = k; // identifies the component and the direction; only used for context determination
          currSE->writing = writeMVD2Buffer_CABAC;
        }
        else
        {
          currSE->value2  = 2*k+bwflag; // identifies the component and the direction; only used for context determination
          currSE->writing = writeBiMVD2Buffer_CABAC;
        }
      }  
      dataPart = &(currSlice->partArr[partMap[img->type==B_IMG ? SE_BFRAME : SE_MVD]]);
      dataPart->writeSyntaxElement (currSE, dataPart);
#if TRACE
      snprintf(currSE->tracestring, TRACESTRING_SIZE, " MVD(%d) = %3d",k, curr_mvd);
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

  int   bframe          = (img->type==B_IMG);
  int** refframe_array  = (img->type==B_IMG ? fw_refFrArr : refFrArr);
#ifdef _ADDITIONAL_REFERENCE_FRAME_
  int   multframe       = (input->no_multpred>1 || input->add_ref_frame>0); 
#else
  int   multframe       = (input->no_multpred>1); 
#endif
  int   step_h0         = (input->blc_size[IS_P8x8(currMB) ? 4 : currMB->mb_type][0] >> 2);
  int   step_v0         = (input->blc_size[IS_P8x8(currMB) ? 4 : currMB->mb_type][1] >> 2);


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
          no_bits += writeReferenceFrame (currMB->b8mode[k], i0, j0, refframe_array[img->block_y+j0][img->block_x+i0]);
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
        refframe  = refframe_array[img->block_y+j0][img->block_x+i0];
        no_bits  += writeMotionVector8x8 (i0, j0, i0+step_h0, j0+step_v0, refframe, currMB->b8mode[k]);
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
        refframe  = refframe_array[img->block_y+j0][img->block_x+i0];
        no_bits  += writeMotionVector8x8 (i0, j0, i0+step_h0, j0+step_v0,       -1, currMB->b8mode[k]);
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
int
writeChromaCoeff ()
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
  int   b8, b4;
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
      DCLevel = img->cofDC[uv+1][0];
      DCRun   = img->cofDC[uv+1][1];

      level=1;
      for (k=0; k < 5 && level != 0; ++k)
      {
        level = currSE->value1 = DCLevel[k]; // level
        run   = currSE->value2 = DCRun  [k]; // run

        if (input->symbol_mode == UVLC)   currSE->mapping = levrun_linfo_c2x2;
        else                              currSE->writing = writeRunLevel2Buffer_CABAC;

        currSE->k = uv;        //ctx for coeff_count

        if (IS_INTRA (currMB))
        {
          currSE->context = 6; // for choosing context model
          currSE->type    = SE_CHR_DC_INTRA;
        }
        else
        {
          currSE->context = 5; // for choosing context model
          currSE->type    = SE_CHR_DC_INTER; 
        }
    
        // choose the appropriate data partition
        if (img->type != B_IMG)   dataPart = &(currSlice->partArr[partMap[currSE->type]]);
        else                      dataPart = &(currSlice->partArr[partMap[SE_BFRAME]]);
    
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
          ACLevel = img->cofAC[b8][b4][0];
          ACRun   = img->cofAC[b8][b4][1];

          ii=i/2;
          i1=i%2;
          level=1;
          uv++;
          for (k=0; k < 16 && level != 0; k++)
          {
            level = currSE->value1 = ACLevel[k]; // level
            run   = currSE->value2 = ACRun  [k]; // run

            if (input->symbol_mode == UVLC)   currSE->mapping = levrun_linfo_inter;
            else                              currSE->writing = writeRunLevel2Buffer_CABAC;
            
            currSE->k=uv;  //ctx for coeff_count

            if (IS_INTRA (currMB))
            {
              currSE->context = 8; // for choosing context model  
              currSE->type  = SE_CHR_AC_INTRA;
            }
            else
            {
              currSE->context = 7; // for choosing context model  
              currSE->type  = SE_CHR_AC_INTER; 
            }

            // choose the appropriate data partition
            if (img->type != B_IMG)   dataPart = &(currSlice->partArr[partMap[currSE->type]]);
            else                      dataPart = &(currSlice->partArr[partMap[SE_BFRAME]]);
        
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

  return rate;
}



/*!
 ************************************************************************
 * \brief
 *    Writes Luma coeff of an 4x4 block
 ************************************************************************
 */
int
writeLumaCoeff4x4 (int b8, int b4, int intra4x4mode)
{
  int             rate      = 0;
  Macroblock*     currMB    = &img->mb_data[img->current_mb_nr];
  SyntaxElement*  currSE    = &img->MB_SyntaxElements[currMB->currSEnr];
  Slice*          currSlice = img->currentSlice;
  const int*      partMap   = assignSE2partition[input->partition_mode];
  int*            bitCount  = currMB->bitcounter;
  DataPartition*  dataPart;

  int   kk,kbeg,kend;
  int   level, run;
  int   k;
  int*  ACLevel = img->cofAC[b8][b4][0];
  int*  ACRun   = img->cofAC[b8][b4][1];

  img->subblock_x = ((b8&0x1)==0)?(((b4&0x1)==0)?0:1):(((b4&0x1)==0)?2:3); // horiz. position for coeff_count context
  img->subblock_y = (b8<2)?((b4<2)?0:1):((b4<2)?2:3); // vert.  position for coeff_count context

  if (intra4x4mode && (img->qp<24 || input->symbol_mode == CABAC))  // double scan, for CABAC always
  {
    for(kk=0;kk<2;kk++)
    {
      kbeg  = kk*9;
      kend  = kbeg+8;
      level = 1; // get inside loop

      for (k=kbeg; k<=kend && level!=0; k++)
      {
        level = currSE->value1 = ACLevel[k];
        run   = currSE->value2 = ACRun  [k];

        if (input->symbol_mode == UVLC)
        {
          currSE->mapping = levrun_linfo_intra;  
        }
        else
        {
          currSE->context = 0; // for choosing context model
          currSE->writing = writeRunLevel2Buffer_CABAC;
        }
              
        if (k == kbeg)
        {
          currSE->type  = SE_LUM_DC_INTRA; // element is of type DC

          // choose the appropriate data partition
          if (img->type!=B_IMG)   dataPart = &(currSlice->partArr[partMap[SE_LUM_DC_INTRA]]);
          else                    dataPart = &(currSlice->partArr[partMap[SE_BFRAME]]);
        }
        else
        {
          currSE->type  = SE_LUM_AC_INTRA;   // element is of type AC
            
          // choose the appropriate data partition
          if (img->type!=B_IMG)   dataPart = &(currSlice->partArr[partMap[SE_LUM_AC_INTRA]]);
          else                    dataPart = &(currSlice->partArr[partMap[SE_BFRAME]]);
        }

        dataPart->writeSyntaxElement (currSE, dataPart);
        bitCount[BITS_COEFF_Y_MB] += currSE->len;
        rate                      += currSE->len;
#if TRACE
        snprintf(currSE->tracestring, TRACESTRING_SIZE, "Luma dbl(%2d,%2d)  level=%3d Run=%2d",kk,k,level,run);
#endif
        // proceed to next SE
        currSE++;  
        currMB->currSEnr++;
      }
    }
  }
  else     // single scan
  {
    level=1; // get inside loop
    for(k=0; k<=16 && level !=0; k++)
    {
      level = currSE->value1 = ACLevel[k]; // level
      run   = currSE->value2 = ACRun  [k]; // run
      
      if (input->symbol_mode == UVLC)  currSE->mapping = levrun_linfo_inter;    
      else                             currSE->writing = writeRunLevel2Buffer_CABAC;

      if (k == 0)
      { 
        if (intra4x4mode)
        {
          currSE->context = 2; // for choosing context model
          currSE->type  = SE_LUM_DC_INTRA;
        }
        else
        {
          currSE->context = 1; // for choosing context model
          currSE->type  = SE_LUM_DC_INTER;
        }
      }
      else
      {
        if (intra4x4mode)
        {
          currSE->context = 2; // for choosing context model
          currSE->type  = SE_LUM_AC_INTRA;
        }
        else
        {
          currSE->context = 1; // for choosing context model
          currSE->type  = SE_LUM_AC_INTER;
        }
      }

      // choose the appropriate data partition
      if (img->type != B_IMG)    dataPart = &(currSlice->partArr[partMap[currSE->type]]);
      else                       dataPart = &(currSlice->partArr[partMap[SE_BFRAME]]);
          
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
    rate += writeLumaCoeff4x4 (block8x8, block4x4, intra4x4mode);
  }

  return rate;
}



/*!
 ************************************************************************
 * \brief
 *    Writes CBP, DQUANT, and Luma Coefficients of an macroblock
 ************************************************************************
 */
int
writeCBPandLumaCoeff ()
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


  if (!IS_NEWINTRA (currMB))
  {
    //=====   C B P   =====
    //---------------------
    currSE->value1 = cbp;
    
    if (IS_OLDINTRA (currMB))
    {
      if (input->symbol_mode == UVLC)  currSE->mapping = cbp_linfo_intra;
      currSE->type = SE_CBP_INTRA;
    }
    else
    {
      if (input->symbol_mode == UVLC)  currSE->mapping = cbp_linfo_inter;
      currSE->type = SE_CBP_INTER;
    }
    if (input->symbol_mode == CABAC)   currSE->writing = writeCBP2Buffer_CABAC;
                      
    // choose the appropriate data partition
    if (img->type==B_IMG)  dataPart = &(img->currentSlice->partArr[partMap[SE_BFRAME]]);
    else                   dataPart = &(img->currentSlice->partArr[partMap[currSE->type]]);

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

    if (IS_INTRA (currMB))
    {
      if (input->symbol_mode==UVLC)   currSE->mapping = dquant_linfo;
      else                            currSE->writing = writeDquant_intra_CABAC;
      currSE->type = SE_DELTA_QUANT_INTRA;
    }
    else
    {
      if (input->symbol_mode==UVLC)   currSE->mapping = dquant_linfo;
      else                            currSE->writing = writeDquant_inter_CABAC;
      currSE->type = SE_DELTA_QUANT_INTER;
    }

    // choose the appropriate data partition
    if (img->type != B_IMG)   dataPart = &(img->currentSlice->partArr[partMap[currSE->type]]);
    else                      dataPart = &(img->currentSlice->partArr[partMap[SE_BFRAME]]);

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
    level=1; // get inside loop
    for (k=0; k<=16 && level!=0; k++)
    {
      level = currSE->value1 = DCLevel[k]; // level
      run   = currSE->value2 = DCRun  [k]; // run

      if (input->symbol_mode == UVLC)
      {
        currSE->mapping = levrun_linfo_inter;
      }
      else
      {
        currSE->context = 3; // for choosing context model
        currSE->writing = writeRunLevel2Buffer_CABAC;
      }
      currSE->type  = SE_LUM_DC_INTRA;   // element is of type DC

      // choose the appropriate data partition
      if (img->type != B_IMG)   dataPart = &(currSlice->partArr[partMap[currSE->type]]);
      else                      dataPart = &(currSlice->partArr[partMap[SE_BFRAME]]);
    
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
        ACLevel = img->cofAC[b8][b4][0];
        ACRun   = img->cofAC[b8][b4][1];

        level=1; // get inside loop
        for (k=0;k<16 && level !=0;k++)
        {
          level = currSE->value1 = ACLevel[k]; // level
          run   = currSE->value2 = ACRun  [k]; // run

          if (input->symbol_mode == UVLC)
          {
            currSE->mapping = levrun_linfo_inter;
          }
          else
          {
            currSE->context = 4; // for choosing context model
            currSE->writing = writeRunLevel2Buffer_CABAC;
          }
          currSE->type  = SE_LUM_AC_INTRA;   // element is of type AC

          // choose the appropriate data partition
          if (img->type != B_IMG)   dataPart = &(currSlice->partArr[partMap[currSE->type]]);
          else                      dataPart = &(currSlice->partArr[partMap[SE_BFRAME]]);

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

  return rate;
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
int find_sad2(int *intra_mode)
{
  int current_intra_sad_2,best_intra_sad2;
  int M1[16][16],M0[4][4][4][4],M3[4],M4[4][4];

  int i,j,k;
  int ii,jj;

  best_intra_sad2=MAX_VALUE;

  for (k=0;k<4;k++)
  {
    int mb_nr = img->current_mb_nr;
    int mb_width = img->width/16;
    int mb_available_up = (img->mb_y == 0) ? 0 : (img->mb_data[mb_nr].slice_nr == img->mb_data[mb_nr-mb_width].slice_nr);
    int mb_available_left = (img->mb_x == 0) ? 0 : (img->mb_data[mb_nr].slice_nr == img->mb_data[mb_nr-1].slice_nr);
    if(input->UseConstrainedIntraPred)
    {
    if (mb_available_up   && (img->intra_block[mb_nr-mb_width][2]==0 || img->intra_block[mb_nr-mb_width][3]==0))
      mb_available_up   = 0;
    if (mb_available_left && (img->intra_block[mb_nr-       1][1]==0 || img->intra_block[mb_nr       -1][3]==0))
      mb_available_left = 0;
    }
    //check if there are neighbours to predict from
    if ((k==0 && !mb_available_up) || (k==1 && !mb_available_left) || (k==3 && (!mb_available_left || !mb_available_up)))
    {
      ; // edge, do nothing
    }
    else
    {
      for (j=0;j<16;j++)
      {
        for (i=0;i<16;i++)
        {
          M1[i][j]=imgY_org[img->pix_y+j][img->pix_x+i]-img->mprr_2[k][j][i];
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
