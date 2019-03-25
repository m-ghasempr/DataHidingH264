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
 *      - Alexis Michael Tourapis       <alexismt@ieee.org>
 ***************************************************************************************
 */

#include "contributors.h"

#include <math.h>
#include <float.h>

#include "global.h"
#include "header.h"
#include "nal.h"
#include "rtp.h"
#include "fmo.h"
#include "vlc.h"
#include "image.h"
#include "cabac.h"
#include "biariencode.h"
#include "elements.h"
#include "macroblock.h"
#include "memalloc.h"
#include "symbol.h"
#include "context_ini.h"
#include "enc_statistics.h"
#include "ratectl.h"
#include "me_epzs.h"
#include "me_epzs_int.h"
#include "wp.h"
#include "slice.h"
#include "rdoq.h"
#include "wp_mcprec.h"
#include "q_offsets.h"
#include "conformance.h"
#include "list_reorder.h"
#include "md_common.h"
#include "mmco.h"
#include "mv_search.h"
#include "quant4x4.h"
#include "quant8x8.h"
#include "quantChroma.h"
#include "rdopt.h"
#include "rdopt_coding_state.h"

// Local declarations
static Slice *malloc_slice(ImageParameters *p_Img, InputParameters *p_Inp);
static void  free_slice   (Slice *currSlice);
static void  set_ref_pic_num(Slice *currSlice);

//! convert from H.263 QP to H.264 quant given by: quant=pow(2,QP/6)

int allocate_block_mem(Slice *currSlice)
{
  int alloc_size = 0;
  alloc_size += get_mem2Dint(&currSlice->tblk4x4, BLOCK_SIZE, BLOCK_SIZE);
  alloc_size += get_mem2Dint(&currSlice->tblk16x16, MB_BLOCK_SIZE, MB_BLOCK_SIZE);
  alloc_size += get_mem4Dint(&currSlice->i16blk4x4, BLOCK_SIZE, BLOCK_SIZE, BLOCK_SIZE, BLOCK_SIZE);
  
  return (alloc_size);
}


void free_block_mem(Slice *currSlice)
{
  free_mem4Dint(currSlice->i16blk4x4);
  free_mem2Dint(currSlice->tblk16x16);
  free_mem2Dint(currSlice->tblk4x4);
}

/*!
 ***********************************************************************
 * \brief
 *    Initializes the p_Img->nz_coeff
 * \par Input:
 *    none
 * \par  Output:
 *    none
 * \ side effects
 *    sets omg->nz_coef[][][][] to -1
 ***********************************************************************
 */
static void CAVLC_init(Slice *currSlice)
{
  memset(&currSlice->p_Img->nz_coeff[0][0][0], 0, currSlice->PicSizeInMbs * 4 * (4 + currSlice->num_blk8x8_uv)* sizeof(int));
}

/*!
 ************************************************************************
 * \brief
 *    Allocate memory for mv
 * \par Input:
 *    Image Parameters ImageParameters *p_Img                             \n
 *    int****** mv
 * \return memory size in bytes
 ************************************************************************
 */
static int get_mem_mv (Slice *currSlice, short ******* mv)
{
  // LIST, reference, block_type, block_y, block_x, component
  get_mem6Dshort(mv, 2, currSlice->max_num_references, 9, 4, 4, 2);

  return 576 * currSlice->max_num_references * sizeof(short); // 2 * ref * 9 * 4 * 4 * 2
}

/*!
 ************************************************************************
 * \brief
 *    Allocate memory for bipredictive mv 
 * \par Input:
 *    Image Parameters ImageParameters *p_Img                             \n
 *    int****** mv
 * \return memory size in bytes
 ************************************************************************
 */
static int get_mem_bipred_mv (Slice *currSlice, short******** bipred_mv) 
{
  get_mem7Dshort(bipred_mv, 2, 2, currSlice->max_num_references, 9, 4, 4, 2);
  
  return 1152 * currSlice->max_num_references * sizeof(short);
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


static int alloc_rddata(Slice *currSlice, RD_DATA *rd_data)
{
  int alloc_size = 0;

  alloc_size += get_mem3Dpel(&(rd_data->rec_mb), 3, MB_BLOCK_SIZE, MB_BLOCK_SIZE);

  alloc_size += get_mem_ACcoeff (currSlice->p_Img, &(rd_data->cofAC));
  alloc_size += get_mem_DCcoeff (&(rd_data->cofDC));  

  if ((currSlice->slice_type != I_SLICE) && (currSlice->slice_type != SI_SLICE))
  {          
    alloc_size += get_mem_mv (currSlice, &(rd_data->all_mv));
  }
  
  // Why is this stored as height_blk * width_blk?
  alloc_size += get_mem2D((byte***)&(rd_data->ipredmode), currSlice->height_blk, currSlice->width_blk);
  alloc_size += get_mem3D((byte****)&(rd_data->refar), 2, 4, 4);

  return alloc_size;
}

static void free_rddata(Slice *currSlice, RD_DATA *rd_data)
{
  free_mem3D((byte***) rd_data->refar);
  free_mem2D((byte**)  rd_data->ipredmode);

  if ((currSlice->slice_type != I_SLICE) && (currSlice->slice_type != SI_SLICE))
  {          
    free_mem_mv (rd_data->all_mv);
  }

  free_mem_DCcoeff (rd_data->cofDC);
  free_mem_ACcoeff (rd_data->cofAC);  

  free_mem3Dpel(rd_data->rec_mb);
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
 *      slice, sets p_Img->currSlice
 ************************************************************************
*/
static int start_slice(Slice *currSlice, StatParameters *cur_stats)
{
  ImageParameters *p_Img = currSlice->p_Img;
  EncodingEnvironmentPtr eep;
  Bitstream *currStream;
  int header_len = 0;
  int i;
  int NumberOfPartitions = (currSlice->partition_mode == PAR_DP_1?1:3);

  //one  partition for an IDR image
  if(currSlice->idr_flag)
  {
    NumberOfPartitions = 1;
  }

  RTPUpdateTimestamp (p_Img, currSlice->frame_no);   // this has no side effects, just leave it for all NALs

  for (i = 0; i < NumberOfPartitions; i++)
  {
    currStream = (currSlice->partArr[i]).bitstream;

    currStream->write_flag = 0;
    if (i==0)     // First partition
      header_len += SliceHeader (currSlice);
    else          // Second/Third partition
      header_len += Partition_BC_Header(currSlice, i);

    //! Initialize CABAC
    if (currSlice->symbol_mode == CABAC)
    {
      eep = &((currSlice->partArr[i]).ee_cabac);
      if (currStream->bits_to_go != 8)
        header_len += currStream->bits_to_go;
      writeVlcByteAlign(p_Img, currStream, cur_stats);
      eep->p_Img = p_Img;
      arienco_start_encoding(eep, currStream->streamBuffer, &(currStream->byte_pos));

      arienco_reset_EC(eep);
    }
    else
    {
      // Initialize CA-VLC
      CAVLC_init(currSlice);
    }
  }

  if(currSlice->symbol_mode == CABAC)
  {
    init_contexts(currSlice);
  }

  return header_len;
}

/*!
************************************************************************
* \brief
*    This creates a NAL unit structures for all data partition of the slice
*
************************************************************************
*/
void create_slice_nalus(Slice *currSlice)
{
  // KS: this is approx. max. allowed code picture size
  //const int buffer_size = 500 + p_Img->FrameSizeInMbs * (128 + 256 * p_Img->bitdepth_luma + 512 * p_Img->bitdepth_chroma);
  int buffer_size = currSlice->partArr[0].bitstream->buffer_size;
  
  int part;
  NALU_t *nalu;

  for (part=0; part< currSlice->max_part_nr; part++)
  {
    if (currSlice->partArr[part].bitstream->write_flag)
    {
      nalu = AllocNALU(buffer_size);
      currSlice->partArr[part].nal_unit = nalu;
      nalu->startcodeprefix_len = 1+ (currSlice->start_mb_nr == 0 && part == 0 ?ZEROBYTES_SHORTSTARTCODE+1:ZEROBYTES_SHORTSTARTCODE);
      nalu->forbidden_bit = 0;

      if (currSlice->idr_flag)
      {
        nalu->nal_unit_type = NALU_TYPE_IDR;
        nalu->nal_reference_idc = NALU_PRIORITY_HIGHEST;
      }
      else
      {
        //different nal header for different partitions
        if(currSlice->partition_mode == 0)
        {
          nalu->nal_unit_type = NALU_TYPE_SLICE;
        }
        else
        {
          nalu->nal_unit_type = (NaluType) (NALU_TYPE_DPA +  part);
        }

        if (currSlice->nal_reference_idc !=0)
        {
          nalu->nal_reference_idc = NALU_PRIORITY_HIGH;
        }
        else
        {
          nalu->nal_reference_idc = NALU_PRIORITY_DISPOSABLE;
        }
      }
    }
    else
    {
      currSlice->partArr[part].nal_unit = NULL;
    }
  }
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
static int terminate_slice(Macroblock *currMB, int lastslice, StatParameters *cur_stats )
{
  Slice *currSlice = currMB->p_slice;
  ImageParameters *p_Img = currMB->p_Img;
  InputParameters *p_Inp = currMB->p_Inp;

  Bitstream *currStream;
  NALU_t    *currNalu;
  EncodingEnvironmentPtr eep;
  int part;
  int tmp_stuffingbits = currMB->bits.mb_stuffing;

  if (currSlice->symbol_mode == CABAC)
    write_terminating_bit (currSlice, 1);      // only once, not for all partitions

  create_slice_nalus(currSlice);

  if (currSlice->slice_type != I_SLICE && currSlice->slice_type != SI_SLICE)
  {
    if (p_Inp->SearchMode == EPZS)
    {
      EPZSStructDelete (currSlice);    
    }
  }

  for (part = 0; part < currSlice->max_part_nr; part++)
  {
    currStream = (currSlice->partArr[part]).bitstream;
    currNalu   = (currSlice->partArr[part]).nal_unit;
    if (currStream->write_flag)
    {
      if (currSlice->symbol_mode == CAVLC)
      {
        SODBtoRBSP(currStream);
        currNalu->len = RBSPtoEBSP(currNalu->buf, currStream->streamBuffer, currStream->byte_pos);
      }
      else     // CABAC
      {
        eep = &((currSlice->partArr[part]).ee_cabac);
        // terminate the arithmetic code
        arienco_done_encoding(currMB, eep);
        set_pic_bin_count(p_Img, eep);

        currStream->bits_to_go = eep->Ebits_to_go;
        currStream->byte_buf = 0;

        currNalu->len = RBSPtoEBSP(currNalu->buf, currStream->streamBuffer, currStream->byte_pos);

        // NumBytesInNALunit is: payload length + 1 byte header
        p_Img->bytes_in_picture += currNalu->len + 1;

        if (lastslice && (part==(currSlice->max_part_nr - 1)))
        {
          addCabacZeroWords(p_Img, currNalu, cur_stats);
        }
      }           // CABAC
    }
  }           // partition loop

  if( currSlice->symbol_mode == CABAC )
  {
    store_contexts(currSlice);
  }

  cur_stats->bit_use_stuffingBits[currSlice->slice_type] += currMB->bits.mb_stuffing - tmp_stuffingbits;

  if (currSlice->slice_type != I_SLICE && currSlice->slice_type != SI_SLICE)
    free_ref_pic_list_reordering_buffer (currSlice);

  return 0;
}

/*!
************************************************************************
* \brief
*    Encodes one slice
* \par
*   returns the number of coded MBs in the SLice
************************************************************************
*/
int encode_one_slice (ImageParameters *p_Img, InputParameters *p_Inp, int SliceGroupId, int TotalCodedMBs)
{
  Boolean end_of_slice = FALSE;
  Boolean recode_macroblock;
  Boolean prev_recode_mb = FALSE;
  int len;
  int NumberOfCodedMBs = 0;
  Macroblock* currMB   = NULL;
  int CurrentMbAddr;
  StatParameters *cur_stats = &p_Img->enc_picture->stats;
  Slice *currSlice = NULL;

  p_Img->Motion_Selected = 0;

  if( IS_INDEPENDENT(p_Inp) )
  {
    change_plane_JV( p_Img, p_Img->colour_plane_id );
  }

  p_Img->cod_counter = 0;

  CurrentMbAddr = FmoGetFirstMacroblockInSlice (p_Img, SliceGroupId);
  // printf ("\n\nEncode_one_slice: PictureID %d SliceGroupId %d  SliceID %d  FirstMB %d \n", p_Img->frame_no, SliceGroupId, p_Img->current_slice_nr, CurrentMbInScanOrder);

  init_slice (p_Img, p_Inp, &currSlice, CurrentMbAddr);
  // Initialize quantization functions based on rounding/quantization method
  // Done here since we may wish to disable adaptive rounding on occasional intervals (even at a frame or gop level).
  init_quant_4x4   (currSlice);
  init_quant_8x8   (currSlice);
  init_quant_Chroma(currSlice);

  currSlice->SetLagrangianMultipliers(p_Img, p_Inp);

  if (currSlice->symbol_mode == CABAC)
  {
    SetCtxModelNumber (currSlice);
  }

  p_Img->checkref = (short) (p_Inp->rdopt && p_Inp->RestrictRef && (p_Img->type==P_SLICE || p_Img->type==SP_SLICE));

  len = start_slice (currSlice, cur_stats);

  // Rate control
  if (p_Inp->RCEnable)
    rc_store_slice_header_bits( p_Img, p_Inp, len );

  // Update statistics
  p_Img->p_Stats->bit_slice += len;
  cur_stats->bit_use_header[currSlice->slice_type] += len;

  if(currSlice->UseRDOQuant == 1 && currSlice->RDOQ_QP_Num > 1)
    get_dQP_table(currSlice);

  while (end_of_slice == FALSE) // loop over macroblocks
  {
    if (p_Img->AdaptiveRounding && p_Inp->AdaptRndPeriod && (p_Img->current_mb_nr % p_Inp->AdaptRndPeriod == 0))
    {
      CalculateOffset4x4Param(p_Img, p_Inp);
      if(p_Inp->Transform8x8Mode)
        CalculateOffset8x8Param(p_Img, p_Inp);
    }

    recode_macroblock = FALSE;
    if(currSlice->UseRDOQuant) // This needs revisit
      currSlice->rddata = &currSlice->rddata_trellis_curr;
    else
      currSlice->rddata = &currSlice->rddata_top_frame_mb;   // store data in top frame MB

    start_macroblock (currSlice,  &currMB, CurrentMbAddr, FALSE);


    if(currSlice->UseRDOQuant)
    {
      trellis_coding(currMB, prev_recode_mb);   
    }
    else
    {
      p_Img->masterQP = p_Img->qp;

      currSlice->encode_one_macroblock (currMB);
      end_encode_one_macroblock(currMB);
      write_macroblock (currMB, 1, prev_recode_mb);    
    }

    end_macroblock (currMB, &end_of_slice, &recode_macroblock);
    prev_recode_mb = recode_macroblock;
    //       printf ("encode_one_slice: mb %d,  slice %d,   bitbuf bytepos %d EOS %d\n",
    //       p_Img->current_mb_nr, p_Img->current_slice_nr,
    //       currSlice->partArr[0].bitstream->byte_pos, end_of_slice);

    if (recode_macroblock == FALSE)       // The final processing of the macroblock has been done
    {
      p_Img->SumFrameQP += currMB->qp;
      CurrentMbAddr = FmoGetNextMBNr (p_Img, CurrentMbAddr);
      if (CurrentMbAddr == -1)   // end of slice
      {
        // printf ("FMO End of Slice Group detected, current MBs %d, force end of slice\n", NumberOfCodedMBs+1);
        end_of_slice = TRUE;
      }
      NumberOfCodedMBs++;       // only here we are sure that the coded MB is actually included in the slice
      next_macroblock (currMB);
    }
    else
    {
      //!Go back to the previous MB to recode it
      p_Img->current_mb_nr = FmoGetPreviousMBNr(p_Img, p_Img->current_mb_nr);
      p_Img->NumberofCodedMacroBlocks--;
      if(p_Img->current_mb_nr == -1 )   // The first MB of the slice group  is too big,
        // which means it's impossible to encode picture using current slice bits restriction
      {
        snprintf (errortext, ET_SIZE, "Error encoding first MB with specified parameter, bits of current MB may be too big");
        error (errortext, 300);
      }
    }
  }


  if ((p_Inp->WPIterMC) && (p_Img->frameOffsetAvail == 0) && p_Img->nal_reference_idc)
  {
    compute_offset(p_Img);
  }
  p_Img->num_ref_idx_l0_active = currSlice->num_ref_idx_active[LIST_0];
  p_Img->num_ref_idx_l1_active = currSlice->num_ref_idx_active[LIST_1];

  terminate_slice (currMB, (NumberOfCodedMBs + TotalCodedMBs >= (int)p_Img->PicSizeInMbs), cur_stats );
  return NumberOfCodedMBs;
}


/*!
************************************************************************
* \brief
*    Encodes one slice (MBAFF Frame)
* \par
*   returns the number of coded MBs in the SLice
************************************************************************
*/
int encode_one_slice_MBAFF (ImageParameters *p_Img, InputParameters *p_Inp, int SliceGroupId, int TotalCodedMBs)
{
  Boolean end_of_slice = FALSE;
  Boolean recode_macroblock;
  Boolean prev_recode_mb = FALSE;
  int len;
  int NumberOfCodedMBs = 0;
  Macroblock* currMB      = NULL;
  int CurrentMbAddr;
  double FrameRDCost = DBL_MAX, FieldRDCost = DBL_MAX;
  StatParameters *cur_stats = &p_Img->enc_picture->stats;
  Slice *currSlice = NULL;

  p_Img->Motion_Selected = 0;

  if( IS_INDEPENDENT(p_Inp) )
  {
    change_plane_JV( p_Img, p_Img->colour_plane_id );
  }

  p_Img->cod_counter = 0;

  CurrentMbAddr = FmoGetFirstMacroblockInSlice (p_Img, SliceGroupId);
  // printf ("\n\nEncode_one_slice: PictureID %d SliceGroupId %d  SliceID %d  FirstMB %d \n", p_Img->frame_no, SliceGroupId, p_Img->current_slice_nr, CurrentMbInScanOrder);

  init_slice (p_Img, p_Inp, &currSlice, CurrentMbAddr);
  // Initialize quantization functions based on rounding/quantization method
  // Done here since we may wish to disable adaptive rounding on occasional intervals (even at a frame or gop level).
  init_quant_4x4   (currSlice);
  init_quant_8x8   (currSlice);
  init_quant_Chroma(currSlice);

  currSlice->SetLagrangianMultipliers(p_Img, p_Inp);

  if (currSlice->symbol_mode == CABAC)
  {
    SetCtxModelNumber (currSlice);
  }

  p_Img->checkref = (short) (p_Inp->rdopt && p_Inp->RestrictRef && (currSlice->slice_type == P_SLICE || currSlice->slice_type == SP_SLICE));

  len = start_slice (currSlice, cur_stats);

  // Rate control
  if (p_Inp->RCEnable)
    rc_store_slice_header_bits( p_Img, p_Inp, len );

  // Update statistics
  p_Img->p_Stats->bit_slice += len;
  cur_stats->bit_use_header[p_Img->type] += len;

  while (end_of_slice == FALSE) // loop over macroblocks
  {

    if (p_Img->AdaptiveRounding && p_Inp->AdaptRndPeriod && (p_Img->current_mb_nr % p_Inp->AdaptRndPeriod == 0))
    {
      CalculateOffset4x4Param(p_Img, p_Inp);
      if(p_Inp->Transform8x8Mode)
        CalculateOffset8x8Param(p_Img, p_Inp);
    }


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
    //!      in macroblock.c, start_macroblock() and end_macroblock() (this code needs
    //!      to be double checked that it works with CA-VLC as well
    //!   2. would it be an option to allocate Bitstreams with zero data in them (or copy the
    //!      already generated bitstream) for the "test coding"?

    p_Img->write_macroblock = FALSE;
    if (p_Inp->MbInterlace == ADAPTIVE_CODING || p_Inp->MbInterlace == FRAME_MB_PAIR_CODING)
    {
      //================ code MB pair as frame MB ================
      //----------------------------------------------------------
      recode_macroblock = FALSE;
      //set mv limits to frame type
      update_mv_limits(p_Img, p_Inp, FALSE);

      p_Img->field_mode = FALSE;  // MB coded as frame
      p_Img->top_field  = FALSE;   // Set top field to 0

      //Rate control
      p_Img->write_macroblock = FALSE;
      p_Img->bot_MB = FALSE;

      // save RC state only when it is going to change
      if ( p_Inp->RCEnable && p_Inp->RCUpdateMode <= MAX_RC_MODE )
      {
        if ( p_Inp->MbInterlace == ADAPTIVE_CODING
          && p_Img->NumberofCodedMacroBlocks > 0 && (p_Img->NumberofCodedMacroBlocks % p_Img->BasicUnit) == 0 )
          rc_copy_quadratic( p_Img, p_Inp, p_Img->p_rc_quad_init, p_Img->p_rc_quad ); // save initial RC status
        if ( p_Inp->MbInterlace == ADAPTIVE_CODING )
          rc_copy_generic( p_Img, p_Img->p_rc_gen_init, p_Img->p_rc_gen ); // save initial RC status
      }

      start_macroblock (currSlice, &currMB, CurrentMbAddr, FALSE);

      currSlice->rddata = &currSlice->rddata_top_frame_mb; // store data in top frame MB
      p_Img->masterQP = p_Img->qp;
      currSlice->encode_one_macroblock (currMB);   // code the MB as frame
      end_encode_one_macroblock(currMB);


      FrameRDCost = currSlice->rddata->min_rdcost;
      //***   Top MB coded as frame MB ***//

      //Rate control
      p_Img->bot_MB = TRUE; //for Rate control

      // go to the bottom MB in the MB pair
      p_Img->field_mode = FALSE;  // MB coded as frame  //GB

      start_macroblock (currSlice, &currMB, CurrentMbAddr + 1, FALSE);
      currSlice->rddata = &currSlice->rddata_bot_frame_mb; // store data in top frame MB
      p_Img->masterQP = p_Img->qp;
      currSlice->encode_one_macroblock (currMB);         // code the MB as frame
      end_encode_one_macroblock(currMB);

      if ( p_Inp->RCEnable && p_Inp->RCUpdateMode <= MAX_RC_MODE )
      {
        if ( p_Inp->MbInterlace == ADAPTIVE_CODING
          && p_Img->NumberofCodedMacroBlocks > 0 && (p_Img->NumberofCodedMacroBlocks % p_Img->BasicUnit) == 0 )
          rc_copy_quadratic( p_Img, p_Inp, p_Img->p_rc_quad_best, p_Img->p_rc_quad ); // restore initial RC status

        if ( p_Inp->MbInterlace == ADAPTIVE_CODING )
          rc_copy_generic( p_Img, p_Img->p_rc_gen_best, p_Img->p_rc_gen ); // save frame RC stats
      }

      FrameRDCost += currSlice->rddata->min_rdcost;
      //***   Bottom MB coded as frame MB ***//
    }

    if ((p_Inp->MbInterlace == ADAPTIVE_CODING) || (p_Inp->MbInterlace == FIELD_CODING))
    {
      //Rate control
      p_Img->bot_MB = FALSE;
      //set mv limits to field type
      update_mv_limits(p_Img, p_Inp, TRUE);

      //=========== start coding the MB pair as a field MB pair =============
      //---------------------------------------------------------------------
      p_Img->field_mode = TRUE;  // MB coded as field
      p_Img->top_field = TRUE;   // Set top field to 1
      p_Img->buf_cycle <<= 1;
      p_Inp->num_ref_frames <<= 1;
      currSlice->num_ref_idx_active[LIST_0] <<= 1;
      currSlice->num_ref_idx_active[LIST_0] += 1;

      if ( p_Inp->RCEnable && p_Inp->RCUpdateMode <= MAX_RC_MODE )
      {
        if ( p_Inp->MbInterlace == ADAPTIVE_CODING
          && p_Img->NumberofCodedMacroBlocks > 0 && (p_Img->NumberofCodedMacroBlocks % p_Img->BasicUnit) == 0 )
          rc_copy_quadratic( p_Img, p_Inp, p_Img->p_rc_quad, p_Img->p_rc_quad_init ); // restore initial RC status

        if ( p_Inp->MbInterlace == ADAPTIVE_CODING )
          rc_copy_generic( p_Img, p_Img->p_rc_gen, p_Img->p_rc_gen_init ); // reset RC stats
      }

      start_macroblock (currSlice, &currMB, CurrentMbAddr, TRUE);

      currSlice->rddata = &currSlice->rddata_top_field_mb; // store data in top frame MB
      //        TopFieldIsSkipped = 0;        // set the top field MB skipped flag to 0
      p_Img->masterQP = p_Img->qp;
      currSlice->encode_one_macroblock (currMB);         // code the MB as field
      end_encode_one_macroblock(currMB);

      FieldRDCost = currSlice->rddata->min_rdcost;
      //***   Top MB coded as field MB ***//
      //Rate control
      p_Img->bot_MB = TRUE;//for Rate control

      p_Img->top_field = FALSE;   // Set top field to 0
      start_macroblock (currSlice, &currMB, CurrentMbAddr+1, TRUE);
      currSlice->rddata = &currSlice->rddata_bot_field_mb; // store data in top frame MB
      p_Img->masterQP = p_Img->qp;
      currSlice->encode_one_macroblock (currMB);         // code the MB as field
      end_encode_one_macroblock(currMB);

      FieldRDCost += currSlice->rddata->min_rdcost;
      //***   Bottom MB coded as field MB ***//
    }

    //Rate control
    p_Img->write_mbaff_frame = 0;  //Rate control

    //=========== decide between frame/field MB pair ============
    //-----------------------------------------------------------
    if ( ((p_Inp->MbInterlace == ADAPTIVE_CODING) && (FrameRDCost < FieldRDCost)) || p_Inp->MbInterlace == FRAME_MB_PAIR_CODING )
    {
      p_Img->field_mode = FALSE;
      p_Img->MBPairIsField = FALSE;
      if ( p_Inp->MbInterlace != FRAME_MB_PAIR_CODING )
      {
        p_Img->buf_cycle >>= 1;
        p_Inp->num_ref_frames >>= 1;
        currSlice->num_ref_idx_active[LIST_0] -= 1;
        currSlice->num_ref_idx_active[LIST_0] >>= 1;
      }

      if ( p_Inp->RCEnable && p_Inp->RCUpdateMode <= MAX_RC_MODE )
      {
        if ( p_Inp->MbInterlace == ADAPTIVE_CODING
          && p_Img->NumberofCodedMacroBlocks > 0 && (p_Img->NumberofCodedMacroBlocks % p_Img->BasicUnit) == 0 )
          rc_copy_quadratic( p_Img, p_Inp, p_Img->p_rc_quad, p_Img->p_rc_quad_best ); // restore initial RC status

        if ( p_Inp->MbInterlace == ADAPTIVE_CODING )
          rc_copy_generic( p_Img, p_Img->p_rc_gen, p_Img->p_rc_gen_best ); // restore frame RC stats
      }

      //Rate control
      p_Img->write_mbaff_frame = 1;  //for Rate control
    }
    else
    {
      p_Img->field_mode = TRUE;
      p_Img->MBPairIsField = TRUE;
    }

    //Rate control
    p_Img->write_macroblock = TRUE;//Rate control

    if (p_Img->MBPairIsField)
      p_Img->top_field = TRUE;
    else
      p_Img->top_field = FALSE;

    //Rate control
    p_Img->bot_MB = FALSE;// for Rate control

    // go back to the Top MB in the MB pair
    start_macroblock (currSlice, &currMB, CurrentMbAddr, p_Img->field_mode);

    currSlice->rddata =  p_Img->field_mode ? &currSlice->rddata_top_field_mb : &currSlice->rddata_top_frame_mb;
    copy_rdopt_data  (currMB);  // copy the MB data for Top MB from the temp buffers
    write_macroblock (currMB, 1, prev_recode_mb);     // write the Top MB data to the bitstream
    end_macroblock   (currMB, &end_of_slice, &recode_macroblock);     // done coding the Top MB
    prev_recode_mb = recode_macroblock;

    if (recode_macroblock == FALSE)       // The final processing of the macroblock has been done
    {
      p_Img->SumFrameQP += currMB->qp;

      CurrentMbAddr = FmoGetNextMBNr (p_Img, CurrentMbAddr);
      if (CurrentMbAddr == -1)   // end of slice
      {
        end_of_slice = TRUE;
      }
      NumberOfCodedMBs++;       // only here we are sure that the coded MB is actually included in the slice
      next_macroblock (currMB);

      //Rate control
      p_Img->bot_MB = TRUE;//for Rate control
      // go to the Bottom MB in the MB pair
      p_Img->top_field = FALSE;
      start_macroblock (currSlice, &currMB, CurrentMbAddr, p_Img->field_mode);

      currSlice->rddata = p_Img->field_mode ? &currSlice->rddata_bot_field_mb : &currSlice->rddata_bot_frame_mb;
      copy_rdopt_data  (currMB);  // copy the MB data for Bottom MB from the temp buffers

      write_macroblock (currMB, 0, prev_recode_mb);     // write the Bottom MB data to the bitstream
      end_macroblock   (currMB, &end_of_slice, &recode_macroblock);     // done coding the Top MB
      prev_recode_mb = recode_macroblock;
      if (recode_macroblock == FALSE)       // The final processing of the macroblock has been done
      {
        p_Img->SumFrameQP += currMB->qp;

        CurrentMbAddr = FmoGetNextMBNr (p_Img, CurrentMbAddr);
        if (CurrentMbAddr == -1)   // end of slice
        {
          end_of_slice = TRUE;
        }
        NumberOfCodedMBs++;       // only here we are sure that the coded MB is actually included in the slice
        next_macroblock (currMB);
      }
      else
      {
        //Go back to the beginning of the macroblock pair to recode it
        p_Img->current_mb_nr = FmoGetPreviousMBNr(p_Img, p_Img->current_mb_nr);
        p_Img->NumberofCodedMacroBlocks -= 2;
        if(p_Img->current_mb_nr == -1 )   // The first MB of the slice group  is too big,
          // which means it's impossible to encode picture using current slice bits restriction
        {
          snprintf (errortext, ET_SIZE, "Error encoding first MB with specified parameter, bits of current MB may be too big");
          error (errortext, 300);
        }
      }
    }
    else
    {
      //!Go back to the previous MB to recode it
      p_Img->current_mb_nr = FmoGetPreviousMBNr(p_Img, p_Img->current_mb_nr);
      p_Img->NumberofCodedMacroBlocks--;
      if(p_Img->current_mb_nr == -1 )   // The first MB of the slice group  is too big,
        // which means it's impossible to encode picture using current slice bits restriction
      {
        snprintf (errortext, ET_SIZE, "Error encoding first MB with specified parameter, bits of current MB may be too big");
        error (errortext, 300);
      }
    }

    if (p_Img->MBPairIsField)    // if MB Pair was coded as field the buffer size variables back to frame mode
    {
      p_Img->buf_cycle >>= 1;
      p_Inp->num_ref_frames >>= 1;
      currSlice->num_ref_idx_active[LIST_0] -= 1;
      currSlice->num_ref_idx_active[LIST_0] >>= 1;
    }

    p_Img->field_mode = p_Img->top_field = FALSE; // reset to frame mode

    if ( !end_of_slice )
    {
      assert( CurrentMbAddr < (int)p_Img->PicSizeInMbs );
      assert( CurrentMbAddr >= 0 );
      if (CurrentMbAddr == FmoGetLastCodedMBOfSliceGroup (p_Img, FmoMB2SliceGroup (p_Img, CurrentMbAddr)))
        end_of_slice = TRUE;        // just in case it doesn't get set in end_macroblock
    }
  }

  p_Img->num_ref_idx_l0_active = currSlice->num_ref_idx_active[LIST_0];
  p_Img->num_ref_idx_l1_active = currSlice->num_ref_idx_active[LIST_1];

  terminate_slice (currMB, (NumberOfCodedMBs + TotalCodedMBs >= (int)p_Img->PicSizeInMbs), cur_stats );
  return NumberOfCodedMBs;
}

static void setup_cabac(Slice *currSlice, char *listXsize)
{
  seq_parameter_set_rbsp_t *active_sps = currSlice->active_sps;
  int i; 

  if (currSlice->slice_type == I_SLICE || currSlice->slice_type == SI_SLICE)
  {
    currSlice->set_modes_and_reframe = set_modes_and_reframe;
    currSlice->writeMB_Skip         = NULL;
    currSlice->writeMB_typeInfo     = writeMB_I_typeInfo_CABAC;
    currSlice->writeB8_typeInfo     = writeB8_typeInfo_CABAC;
    currSlice->writeMotionInfo2NAL  = NULL;
    currSlice->write_MB_layer       = writeMBLayerISlice;
    for (i=0; i<6; i++)
    {
      currSlice->writeRefFrame[i]   = NULL;
    }
  }
  else if (currSlice->slice_type == B_SLICE)
  {
    currSlice->set_modes_and_reframe = set_modes_and_reframe_b_slice;
    currSlice->writeMB_Skip         = writeMB_Bskip_flagInfo_CABAC;
    currSlice->writeMB_typeInfo     = writeMB_B_typeInfo_CABAC;
    currSlice->writeB8_typeInfo     = writeB8_B_typeInfo_CABAC;
    currSlice->writeMotionInfo2NAL  = write_bslice_motion_info_to_NAL;
    currSlice->write_MB_layer       = writeMBLayerBSlice;
    for (i=0; i<6; i++)
    {
      switch (listXsize[i])
      {
      case 0:
        currSlice->writeRefFrame[i]   = NULL;
        break;
      case 1:
        currSlice->writeRefFrame[i]   = writeSE_Dummy;
        break;
      default:
        currSlice->writeRefFrame[i]   = writeRefPic_B_CABAC;
      }
    }
  }
  else
  {
    currSlice->set_modes_and_reframe = set_modes_and_reframe_p_slice;
    currSlice->writeMB_Skip         = writeMB_Pskip_flagInfo_CABAC;
    currSlice->writeMB_typeInfo     = writeMB_P_typeInfo_CABAC;
    currSlice->writeB8_typeInfo     = writeB8_typeInfo_CABAC;
    currSlice->writeMotionInfo2NAL  = write_pslice_motion_info_to_NAL;
    currSlice->write_MB_layer       = writeMBLayerPSlice;

    for (i=0; i<6; i++)
    {
      switch (listXsize[i])
      {
      case 0:
        currSlice->writeRefFrame[i]   = NULL;
        break;
      case 1:
        currSlice->writeRefFrame[i]   = writeSE_Dummy;
        break;
      default:
        currSlice->writeRefFrame[i]   = writeRefPic_P_CABAC;
      }
    }
  }

  currSlice->writeIntraPredMode     = writeIntraPredMode_CABAC;  
  currSlice->writeCoeff16x16        = writeCoeff16x16_CABAC;
  currSlice->writeMVD               = writeMVD_CABAC;
  currSlice->writeCBP               = writeCBP_CABAC;
  currSlice->writeDquant            = writeDquant_CABAC;
  currSlice->writeCIPredMode        = writeCIPredMode_CABAC;
  currSlice->writeFieldModeInfo     = writeFieldModeInfo_CABAC;
  currSlice->writeMB_transform_size = writeMB_transform_size_CABAC;

  if (active_sps->chroma_format_idc == YUV444)
    currSlice->write_and_store_CBP_block_bit = write_and_store_CBP_block_bit_444;
  else
    currSlice->write_and_store_CBP_block_bit = write_and_store_CBP_block_bit;

  memset(currSlice->coeff, 0 , 64 * sizeof(int));
  currSlice->coeff_ctr = 0;
  currSlice->pos       = 0;
}

static void setup_cavlc(Slice *currSlice, char *listXsize)
{
  seq_parameter_set_rbsp_t *active_sps = currSlice->active_sps;

  int i;
  currSlice->writeMB_typeInfo   = writeUVLC_CAVLC;
  currSlice->writeIntraPredMode = writeIntraPredMode_CAVLC;
  currSlice->writeB8_typeInfo   = writeSE_UVLC;
  currSlice->writeCoeff16x16    = writeCoeff16x16_CAVLC;
  for (i=0; i<6; i++)
  {
    switch (listXsize[i])
    {
    case 0:
      currSlice->writeRefFrame[i]   = NULL;
      break;
    case 1:
      currSlice->writeRefFrame[i]   = writeSE_Dummy;
      break;
    case 2:
      currSlice->writeRefFrame[i]   = writeSE_invFlag;
      break;
    default:
      currSlice->writeRefFrame[i]   = writeSE_UVLC;
      break;
    }
  }

  if (currSlice->slice_type == I_SLICE || currSlice->slice_type == SI_SLICE)
  {
    currSlice->set_modes_and_reframe = set_modes_and_reframe;
    currSlice->writeMotionInfo2NAL  = NULL;
    currSlice->write_MB_layer       = writeMBLayerISlice;
  }
  else if (currSlice->slice_type == B_SLICE)
  {
    currSlice->set_modes_and_reframe = set_modes_and_reframe_b_slice;
    currSlice->writeMotionInfo2NAL  = write_bslice_motion_info_to_NAL;
    currSlice->write_MB_layer       = writeMBLayerBSlice;
  }
  else
  {
    currSlice->set_modes_and_reframe = set_modes_and_reframe_p_slice;
    currSlice->writeMotionInfo2NAL  = write_pslice_motion_info_to_NAL;
    currSlice->write_MB_layer       = writeMBLayerPSlice;
  }

  currSlice->writeMVD               = writeSVLC_CAVLC;
  currSlice->writeCBP               = writeCBP_VLC;
  currSlice->writeDquant            = writeSVLC_CAVLC;
  currSlice->writeCIPredMode        = writeUVLC_CAVLC;
  currSlice->writeFieldModeInfo     = writeFlag_CAVLC;
  currSlice->writeMB_transform_size = writeFlag_CAVLC;

  // We should move this to the sequence level and create a new function
  if (active_sps->chroma_format_idc == YUV444)
    currSlice->writeCoeff4x4_CAVLC = writeCoeff4x4_CAVLC_444;
  else
    currSlice->writeCoeff4x4_CAVLC = writeCoeff4x4_CAVLC_normal;
}

/*!
 ************************************************************************
 * \brief
 *    Initializes the parameters for a new slice and
 *     allocates the memory for the coded slice in the Picture structure
 *  \par Side effects:
 *      Adds slice/partition header symbols to the symbol buffer
 *      increments Picture->no_slices, allocates memory for the
 *      slice, sets p_Img->currSlice
 ************************************************************************
 */
void init_slice (ImageParameters *p_Img, InputParameters *p_Inp, Slice **currSlice, int start_mb_addr)
{
  int i,j;
  seq_parameter_set_rbsp_t *active_sps = p_Img->active_sps;
  Picture *currPic = p_Img->currentPicture;
  DataPartition *dataPart;
  Bitstream *currStream;
  int active_ref_lists = (p_Img->MbaffFrameFlag) ? 6 : 2;

  p_Img->current_mb_nr = start_mb_addr;

  // Allocate new Slice in the current Picture, and set p_Img->currentSlice
  assert (currPic != NULL);
  currPic->no_slices++;

  if (currPic->no_slices >= MAXSLICEPERPICTURE)
    error ("Too many slices per picture, increase MAXSLICEPERPICTURE in global.h.", -1);

  currPic->slices[currPic->no_slices - 1] = malloc_slice(p_Img, p_Inp);
  *currSlice = currPic->slices[currPic->no_slices-1];

  p_Img->currentSlice = *currSlice;
  // Using this trick we can basically reference back to the image and input parameter structures
  (*currSlice)->p_Img             = p_Img;
  (*currSlice)->p_Inp             = p_Inp;
  (*currSlice)->active_sps        = active_sps;
  (*currSlice)->active_pps        = p_Img->active_pps;

  //printf("Value %d\n", ((ImageParameters *)(*currSlice)->p_Img)->frame_no);

  (*currSlice)->picture_id        = (p_Img->frame_no & 0xFF); // % 255
  (*currSlice)->slice_nr          = p_Img->current_slice_nr;
  (*currSlice)->idr_flag          = p_Img->currentPicture->idr_flag;
  (*currSlice)->slice_type        = p_Img->type;
  (*currSlice)->frame_no          = p_Img->frame_no;
  (*currSlice)->frame_num         = p_Img->frame_num;
  (*currSlice)->max_frame_num     = p_Img->max_frame_num;
  (*currSlice)->framepoc          = p_Img->framepoc;
  (*currSlice)->ThisPOC           = p_Img->ThisPOC;      
  (*currSlice)->qp                = p_Img->qp;
  (*currSlice)->start_mb_nr       = start_mb_addr;
  (*currSlice)->colour_plane_id   = p_Img->colour_plane_id;

  (*currSlice)->si_frame_indicator =  p_Img->si_frame_indicator;
  (*currSlice)->sp2_frame_indicator = p_Img->sp2_frame_indicator;

  (*currSlice)->P444_joined       = p_Img->P444_joined;
  (*currSlice)->disthres          = p_Inp->disthres;
  (*currSlice)->UseRDOQuant       = p_Inp->UseRDOQuant;
  (*currSlice)->RDOQ_QP_Num       = p_Inp->RDOQ_QP_Num;
  (*currSlice)->Transform8x8Mode  = p_Inp->Transform8x8Mode;

  (*currSlice)->slice_too_big     = dummy_slice_too_big;  
  (*currSlice)->width_blk         = p_Img->width_blk;
  (*currSlice)->height_blk        = p_Img->height_blk;
  (*currSlice)->partition_mode    = (short) p_Inp->partition_mode;
  (*currSlice)->PicSizeInMbs      = p_Img->PicSizeInMbs;
  (*currSlice)->num_blk8x8_uv     = p_Img->num_blk8x8_uv;
  (*currSlice)->nal_reference_idc = p_Img->nal_reference_idc;
  (*currSlice)->bitdepth_luma     = p_Img->bitdepth_luma;
  (*currSlice)->bitdepth_chroma   = p_Img->bitdepth_chroma;  

  for (i = 0; i < (*currSlice)->max_part_nr; i++)
  {
    dataPart = &(*currSlice)->partArr[i];

    currStream = dataPart->bitstream;
    currStream->bits_to_go = 8;
    currStream->byte_pos = 0;
    currStream->byte_buf = 0;
  }

  (*currSlice)->direct_spatial_mv_pred_flag  = p_Img->direct_spatial_mv_pred_flag;
  (*currSlice)->MbaffFrameFlag = p_Img->MbaffFrameFlag;
  (*currSlice)->structure      = p_Img->structure;

  (*currSlice)->num_ref_idx_active[LIST_0] = p_Img->active_pps->num_ref_idx_l0_active_minus1 + 1;
  (*currSlice)->num_ref_idx_active[LIST_1] = p_Img->active_pps->num_ref_idx_l1_active_minus1 + 1;

  (*currSlice)->DFDisableIdc    = p_Img->DFDisableIdc;
  (*currSlice)->DFAlphaC0Offset = p_Img->DFAlphaC0Offset;
  (*currSlice)->DFBetaOffset    = p_Img->DFBetaOffset;

  // primary and redundant slices: number of references overriding.
  if(p_Inp->redundant_pic_flag)
  {
    if(!p_Img->redundant_coding)
    {
      (*currSlice)->num_ref_idx_active[LIST_0] = (char) imin(p_Img->number,p_Inp->NumRefPrimary);
    }
    else
    {
      // 1 reference picture for redundant slices
      (*currSlice)->num_ref_idx_active[LIST_0] = 1;
    }
  }

  // code now also considers fields. Issue whether we should account this within the appropriate input p_Inp directly
  if ((p_Img->type == P_SLICE || p_Img->type == SP_SLICE) && p_Inp->P_List0_refs)
  {
    (*currSlice)->num_ref_idx_active[LIST_0] = (char) imin((*currSlice)->num_ref_idx_active[LIST_0], p_Inp->P_List0_refs * ((p_Img->structure !=0) + 1));
  }

  if (p_Img->type == B_SLICE )
  {
    if (p_Inp->B_List0_refs)
    {
      (*currSlice)->num_ref_idx_active[LIST_0] = (char) imin((*currSlice)->num_ref_idx_active[LIST_0], p_Inp->B_List0_refs * ((p_Img->structure !=0) + 1));
    }
    if (p_Inp->B_List1_refs)
    {
      (*currSlice)->num_ref_idx_active[LIST_1] = (char) imin((*currSlice)->num_ref_idx_active[LIST_1], p_Inp->B_List1_refs * ((p_Img->structure !=0) + 1));
    }
    get_mem3D((byte ****)(void*)&(*currSlice)->direct_ref_idx, (*currSlice)->height_blk, (*currSlice)->width_blk, 2);
    get_mem2D((byte ***) (void*)&(*currSlice)->direct_pdir,    (*currSlice)->height_blk, (*currSlice)->width_blk);
  }

  // generate reference picture lists
  init_lists(*currSlice);

  // assign list 0 size from list size
  (*currSlice)->num_ref_idx_active[LIST_0] = p_Img->listXsize[0];
  (*currSlice)->num_ref_idx_active[LIST_1] = p_Img->listXsize[1];

  if ( p_Inp->WPMCPrecision && p_Inp->WPMCPrecFullRef )
    wpxAdaptRefNum(*currSlice);

  //Perform memory management based on poc distances  
  if (p_Inp->SetFirstAsLongTerm && p_Img->number == 0)
  {
    mmco_long_term(p_Img, p_Img->number);
  }
  else if (p_Img->nal_reference_idc && p_Inp->PocMemoryManagement)
  {
    if (p_Img->structure == FRAME && p_Img->p_Dpb->ref_frames_in_buffer == active_sps->num_ref_frames)
      poc_based_ref_management_frame_pic(p_Img, p_Img->frame_num);
    else if (p_Img->structure == TOP_FIELD && p_Img->p_Dpb->ref_frames_in_buffer== active_sps->num_ref_frames)
      poc_based_ref_management_field_pic(p_Img, (p_Img->frame_num << 1) + 1);      
    else if (p_Img->structure == BOTTOM_FIELD)
      poc_based_ref_management_field_pic(p_Img, (p_Img->frame_num << 1) + 1);
  }

  if (p_Inp->EnableOpenGOP)
  {
    for (i = 0; i<p_Img->listXsize[0]; i++)
    {
      if (p_Img->listX[0][i]->poc < p_Img->last_valid_reference && p_Img->ThisPOC > p_Img->last_valid_reference)
      {
        p_Img->listXsize[0] = (*currSlice)->num_ref_idx_active[LIST_0] = (char) imax(1, i);
        break;
      }
    }

    for (i = 0; i<p_Img->listXsize[1]; i++)
    {
      if (p_Img->listX[1][i]->poc < p_Img->last_valid_reference && p_Img->ThisPOC > p_Img->last_valid_reference)
      {
        p_Img->listXsize[1] = (*currSlice)->num_ref_idx_active[LIST_1] = (char) imax(1,i);
        break;
      }
    }
  }

  init_ref_pic_list_reordering(*currSlice);

    // reference list reordering 
    // RPLR for redundant pictures
    // !KS: that should actually be moved somewhere else
    if(p_Inp->redundant_pic_flag && p_Img->redundant_coding)
    {
      (*currSlice)->ref_pic_list_reordering_flag[LIST_0] = 1;
      (*currSlice)->reordering_of_pic_nums_idc[LIST_0][0] = 0;
      (*currSlice)->reordering_of_pic_nums_idc[LIST_0][1] = 3;
      (*currSlice)->abs_diff_pic_num_minus1[LIST_0][0] = p_Img->redundant_ref_idx - 1;
      (*currSlice)->long_term_pic_idx[LIST_0][0] = 0;
      reorder_ref_pic_list ( *currSlice, p_Img->listX, p_Img->listXsize, LIST_0);
    }
    else if ( (p_Img->type == P_SLICE || p_Img->type == B_SLICE) && p_Inp->WPMCPrecision && p_Img->pWPX->curr_wp_rd_pass->algorithm != WP_REGULAR )
      wp_mcprec_reorder_lists( *currSlice );
    else
      reorder_lists( *currSlice );

  (*currSlice)->max_num_references = (short) p_Img->max_num_references;

  if (((*currSlice)->slice_type != I_SLICE) && ((*currSlice)->slice_type != SI_SLICE))
  {
    get_mem_mv(*currSlice, &(*currSlice)->all_mv);  

    if (p_Inp->BiPredMotionEstimation && ((*currSlice)->slice_type == B_SLICE))
    {
      get_mem_bipred_mv(*currSlice, &(*currSlice)->bipred_mv);
    }

    if (p_Inp->UseRDOQuant && p_Inp->RDOQ_QP_Num > 1)
    {
      if (p_Inp->Transform8x8Mode && p_Inp->RDOQ_CP_MV)
      {
        get_mem4Dmv (&(*currSlice)->tmp_mv8, 2, (*currSlice)->max_num_references, 4, 4);
        get_mem3Dint(&(*currSlice)->motion_cost8, 2, (*currSlice)->max_num_references, 4);
      }
    }
  }

  if (p_Img->MbaffFrameFlag)
    init_mbaff_lists(*currSlice);

  InitWP(p_Img, p_Inp);

  if ((*currSlice)->slice_type == B_SLICE)
    (*currSlice)->weighted_prediction = p_Img->active_pps->weighted_bipred_idc;
  else if (((*currSlice)->slice_type != I_SLICE) || ((*currSlice)->slice_type != SI_SLICE))
    (*currSlice)->weighted_prediction = p_Img->active_pps->weighted_pred_flag;

  if (p_Img->type != I_SLICE && (p_Img->active_pps->weighted_pred_flag == 1 || (p_Img->active_pps->weighted_bipred_idc > 0 && (p_Img->type == B_SLICE))))
  {
    if (p_Img->type == P_SLICE || p_Img->type == SP_SLICE)
    {
      int wp_type = (p_Inp->GenerateMultiplePPS && p_Inp->RDPictureDecision) && (p_Img->enc_picture != p_Img->enc_frame_picture[1]);
      p_Img->EstimateWPPSlice (*currSlice, wp_type);
    }
    else
      p_Img->EstimateWPBSlice (*currSlice);
  }

  set_ref_pic_num(*currSlice);


  if (p_Img->type == B_SLICE)
  {
    int PicWidthInMbs, PicHeightInMapUnits, FrameHeightInMbs;
    int width, height;
    PicWidthInMbs = (active_sps->pic_width_in_mbs_minus1 +1);
    PicHeightInMapUnits = (active_sps->pic_height_in_map_units_minus1 +1);
    FrameHeightInMbs = ( 2 - active_sps->frame_mbs_only_flag ) * PicHeightInMapUnits;

    width = PicWidthInMbs * MB_BLOCK_SIZE;
    height = FrameHeightInMbs * MB_BLOCK_SIZE;

    if( IS_INDEPENDENT(p_Inp) )
    {
       (*currSlice)->p_colocated = alloc_colocated (p_Img->width, p_Img->height, active_sps->mb_adaptive_frame_field_flag);           
      compute_colocated_JV(*currSlice, (*currSlice)->p_colocated, p_Img->listX);
    }
    else
    {
      // Allocation should be based on slice size, not image
      (*currSlice)->p_colocated = alloc_colocated (p_Img->width, p_Img->height, active_sps->mb_adaptive_frame_field_flag);
      compute_colocated(*currSlice, (*currSlice)->p_colocated, p_Img->listX);
    }
  }

  if (p_Img->type != I_SLICE && p_Img->type != SI_SLICE)
  {
    if (p_Inp->SearchMode == EPZS)
    {
      if (((*currSlice)->p_EPZS =  (EPZSParameters*) calloc(1, sizeof(EPZSParameters)))==NULL) 
        no_mem_exit("alloc_img: p_EPZS");  

      EPZSStructInit (*currSlice);
      EPZSSliceInit  (*currSlice);
    }
  }


  if ((*currSlice)->symbol_mode == CAVLC)
  {
    setup_cavlc(*currSlice, p_Img->listXsize);
  }
  else
  {
    setup_cabac(*currSlice, p_Img->listXsize);
  }

  // assign luma common reference picture pointers to be used for ME/sub-pel interpolation

  for(i = 0; i < active_ref_lists; i++)
  {
    for(j = 0; j < p_Img->listXsize[i]; j++)
    {
      if( p_Img->listX[i][j] )
      {
        p_Img->listX[i][j]->p_curr_img     = p_Img->listX[i][j]->p_img    [(short) p_Img->colour_plane_id];
        p_Img->listX[i][j]->p_curr_img_sub = p_Img->listX[i][j]->p_img_sub[(short) p_Img->colour_plane_id];
      }
    }
  }

  if (p_Inp->UseRDOQuant)
  {
    if (((*currSlice)->estBitsCabac = (estBitsCabacStruct*) calloc(NUM_BLOCK_TYPES, sizeof(estBitsCabacStruct)))==NULL) 
      no_mem_exit("init_slice: (*currSlice)->estBitsCabac"); 

    init_rdoq_slice(*currSlice);

    alloc_rddata(*currSlice, &(*currSlice)->rddata_trellis_curr);
    if (p_Inp->RDOQ_QP_Num > 1)
    {
      alloc_rddata(*currSlice, &(*currSlice)->rddata_trellis_best);
    }
  }

  if(p_Img->MbaffFrameFlag)
  {
    alloc_rddata(*currSlice, &(*currSlice)->rddata_top_frame_mb);
    alloc_rddata(*currSlice, &(*currSlice)->rddata_bot_frame_mb);
    if ( p_Inp->MbInterlace != FRAME_MB_PAIR_CODING )
    {
      alloc_rddata(*currSlice, &(*currSlice)->rddata_top_field_mb);
      alloc_rddata(*currSlice, &(*currSlice)->rddata_bot_field_mb);
    }
  }

  if ((*currSlice)->slice_type == B_SLICE)
  {
    (*currSlice)->SetMotionVectorsMB = SetMotionVectorsMBBSlice;
  }
  else if ((*currSlice)->slice_type == P_SLICE || (*currSlice)->slice_type == SP_SLICE)
  {
    (*currSlice)->SetMotionVectorsMB = SetMotionVectorsMBPSlice;
  }
  else
  {
    (*currSlice)->SetMotionVectorsMB = SetMotionVectorsMBISlice;
  }


  if (p_Img->MbaffFrameFlag)
  {
    if (p_Img->direct_spatial_mv_pred_flag)  //spatial direct 
      (*currSlice)->Get_Direct_Motion_Vectors = Get_Direct_MV_Spatial_MBAFF;
    else
      (*currSlice)->Get_Direct_Motion_Vectors = Get_Direct_MV_Temporal;
  }
  else
  {
    if (p_Img->direct_spatial_mv_pred_flag)  //spatial direct 
      (*currSlice)->Get_Direct_Motion_Vectors = Get_Direct_MV_Spatial_Normal;
    else
      (*currSlice)->Get_Direct_Motion_Vectors = Get_Direct_MV_Temporal;
  }

  get_mem3Dpel(&((*currSlice)->mb_pred),   MAX_PLANE, MB_BLOCK_SIZE, MB_BLOCK_SIZE);
  get_mem3Dint(&((*currSlice)->mb_rres),   MAX_PLANE, MB_BLOCK_SIZE, MB_BLOCK_SIZE);
  get_mem3Dint(&((*currSlice)->mb_ores),   MAX_PLANE, MB_BLOCK_SIZE, MB_BLOCK_SIZE);
  get_mem4Dpel(&((*currSlice)->mpr_4x4),   MAX_PLANE, 9, MB_BLOCK_SIZE, MB_BLOCK_SIZE);
  get_mem4Dpel(&((*currSlice)->mpr_8x8),   MAX_PLANE, 9, MB_BLOCK_SIZE, MB_BLOCK_SIZE);
  get_mem4Dpel(&((*currSlice)->mpr_16x16), MAX_PLANE, 5, MB_BLOCK_SIZE, MB_BLOCK_SIZE);

  get_mem_ACcoeff (p_Img, &((*currSlice)->cofAC));
  get_mem_DCcoeff (&((*currSlice)->cofDC));


  allocate_block_mem(*currSlice);
  init_coding_state_methods(*currSlice);
  init_rdopt(*currSlice);

}


/*!
 ************************************************************************
 * \brief
 *    Allocates a slice structure along with its dependent data structures
 * \return
 *    Pointer to a Slice
 ************************************************************************
 */
static Slice *malloc_slice(ImageParameters *p_Img, InputParameters *p_Inp)
{
  int i;
  DataPartition *dataPart;
  Slice *currSlice;
  int cr_size = IS_INDEPENDENT( p_Inp ) ? 0 : 512;

  int buffer_size;

  switch (p_Inp->slice_mode)
  {
  case 2:
    //buffer_size = imax(2 * p_Inp->slice_argument, 500 + (128 + 256 * p_Img->bitdepth_luma + 512 * p_Img->bitdepth_chroma));
    buffer_size = imax(2 * p_Inp->slice_argument, 764);
    break;
  case 1:
    buffer_size = 500 + p_Inp->slice_argument * (128 + 256 * p_Img->bitdepth_luma + cr_size * p_Img->bitdepth_chroma);
    break;
  default:
    buffer_size = 500 + p_Img->FrameSizeInMbs * (128 + 256 * p_Img->bitdepth_luma + cr_size * p_Img->bitdepth_chroma);
    break;
  }

  // KS: this is approx. max. allowed code picture size
  if ((currSlice = (Slice *) calloc(1, sizeof(Slice))) == NULL) no_mem_exit ("malloc_slice: currSlice structure");

  currSlice->p_Img             = p_Img;
  currSlice->p_Inp             = p_Inp;

  if (((currSlice->p_RDO)  = (RDOPTStructure *) calloc(1, sizeof(RDOPTStructure)))==NULL) 
    no_mem_exit("malloc_slice: p_RDO");

  currSlice->symbol_mode  = (char) p_Inp->symbol_mode;

  if (currSlice->symbol_mode == CABAC)
  {
    // create all context models
    currSlice->mot_ctx = create_contexts_MotionInfo ();
    currSlice->tex_ctx = create_contexts_TextureInfo();
  }

  currSlice->max_part_nr = p_Inp->partition_mode==0?1:3;

  //for IDR p_Img there should be only one partition
  if(p_Img->currentPicture->idr_flag)
    currSlice->max_part_nr = 1;

  assignSE2partition[0] = assignSE2partition_NoDP;
  //ZL
  //for IDR p_Img all the syntax element should be mapped to one partition
  if(!p_Img->currentPicture->idr_flag && p_Inp->partition_mode == 1)
    assignSE2partition[1] =  assignSE2partition_DP;
  else
    assignSE2partition[1] =  assignSE2partition_NoDP;

  currSlice->num_mb = 0;          // no coded MBs so far

  if ((currSlice->partArr = (DataPartition *) calloc(currSlice->max_part_nr, sizeof(DataPartition))) == NULL) 
    no_mem_exit ("malloc_slice: partArr");
  for (i=0; i<currSlice->max_part_nr; i++) // loop over all data partitions
  {
    dataPart = &(currSlice->partArr[i]);
    if ((dataPart->bitstream = (Bitstream *) calloc(1, sizeof(Bitstream))) == NULL) 
      no_mem_exit ("malloc_slice: Bitstream");
    if ((dataPart->bitstream->streamBuffer = (byte *) calloc(buffer_size, sizeof(byte))) == NULL) 
      no_mem_exit ("malloc_slice: StreamBuffer");
    dataPart->bitstream->buffer_size = buffer_size;
    // Initialize storage of bitstream parameters
    // Set pointers
    dataPart->p_Slice = currSlice;
    dataPart->p_Img   = p_Img;
    dataPart->p_Inp   = p_Inp;
  }

  if (p_Inp->WeightedPrediction || p_Inp->WeightedBiprediction || p_Inp->GenerateMultiplePPS)
  {
    // Currently only use up to 32 references. Need to use different indicator such as maximum num of references in list
    get_mem3Dshort(&currSlice->wp_weight, 6, MAX_REFERENCE_PICTURES, 3);
    get_mem3Dshort(&currSlice->wp_offset, 6, MAX_REFERENCE_PICTURES, 3);
    get_mem4Dshort(&currSlice->wbp_weight, 6, MAX_REFERENCE_PICTURES, MAX_REFERENCE_PICTURES, 3);
  }

  return currSlice;
}


/*!
 ************************************************************************
 * \brief
 *    This function frees nal units
 *
 ************************************************************************
 */
static void free_nal_unit(Picture *pic)
{
  int partition, slice;
  Slice  *currSlice;

  // loop over all slices of the picture
  for (slice=0; slice < pic->no_slices; slice++)
  {
    currSlice = pic->slices[slice];

    // loop over the partitions
    if (currSlice != NULL)
    {
      for (partition=0; partition < currSlice->max_part_nr; partition++)
      {
        // free only if the partition has content
        if (currSlice->partArr[partition].bitstream->write_flag )
        {
          if (currSlice->partArr[partition].nal_unit != NULL)
          {
            FreeNALU(currSlice->partArr[partition].nal_unit);
            currSlice->partArr[partition].nal_unit = NULL;
          }
        }
      }
    }
  }
}

/*!
 ************************************************************************
 * \brief
 *    Memory frees of all Slice structures and of its dependent
 *    data structures
 * \par Input:
 *    Picture *currPic
 ************************************************************************
 */
void free_slice_list(Picture *currPic)
{
  int i;

  if (currPic !=  NULL)
  {
    free_nal_unit(currPic);

    for (i = 0; i < currPic->no_slices; i++)
    {
      free_slice (currPic->slices[i]);
      currPic->slices[i] = NULL;
    }
  }
}


/*!
 ************************************************************************
 * \brief
 *    Memory frees of the Slice structure and of its dependent
 *    data structures
 * \param currSlice
 *    Slice to be freed
 ************************************************************************
 */
static void free_slice(Slice *currSlice)
{
  if (currSlice != NULL)
  {
    ImageParameters *p_Img = currSlice->p_Img;
    InputParameters *p_Inp = currSlice->p_Inp;

    int i;
    DataPartition *dataPart;

    for (i=0; i<currSlice->max_part_nr; i++) // loop over all data partitions
    {
      dataPart = &(currSlice->partArr[i]);

      if (dataPart != NULL)
      {
        if (dataPart->bitstream != NULL)
        {
          if (dataPart->bitstream->streamBuffer != NULL)
          {
            free(dataPart->bitstream->streamBuffer);       
            dataPart->bitstream->streamBuffer = NULL;
          }
          free(dataPart->bitstream);
          dataPart->bitstream = NULL;
        }
      }
    }

    // free structure for rd-opt. mode decision
    clear_rdopt (currSlice);
    free (currSlice->p_RDO);

    free_mem_ACcoeff (currSlice->cofAC);
    free_mem_DCcoeff (currSlice->cofDC);

    free_mem3Dint(currSlice->mb_rres  );
    free_mem3Dint(currSlice->mb_ores  );
    free_mem3Dpel(currSlice->mb_pred  );
    free_mem4Dpel(currSlice->mpr_16x16);
    free_mem4Dpel(currSlice->mpr_8x8  );
    free_mem4Dpel(currSlice->mpr_4x4  );

    if (currSlice->slice_type == B_SLICE)
    {
      free_colocated(currSlice->p_colocated);
    }


    if (currSlice->partArr != NULL)
    {
      free(currSlice->partArr);
      currSlice->partArr = NULL;
    }

    if (currSlice->symbol_mode == CABAC)
    {
      delete_contexts_MotionInfo(currSlice->mot_ctx);
      delete_contexts_TextureInfo(currSlice->tex_ctx);
    }

    if (p_Inp->WeightedPrediction || p_Inp->WeightedBiprediction || p_Inp->GenerateMultiplePPS)
    {
      free_mem3Dshort(currSlice->wp_weight );
      free_mem3Dshort(currSlice->wp_offset );
      free_mem4Dshort(currSlice->wbp_weight);
    }

    if (currSlice->UseRDOQuant)
    {
      free(currSlice->estBitsCabac);

      free_rddata(currSlice, &currSlice->rddata_trellis_curr);
      if(currSlice->RDOQ_QP_Num > 1)
      {
        free_rddata(currSlice, &currSlice->rddata_trellis_best);
      }
    }

    if(p_Img->MbaffFrameFlag)
    {
      free_rddata(currSlice, &currSlice->rddata_top_frame_mb);
      free_rddata(currSlice, &currSlice->rddata_bot_frame_mb);

      if ( p_Inp->MbInterlace != FRAME_MB_PAIR_CODING )
      {
        free_rddata(currSlice, &currSlice->rddata_top_field_mb);
        free_rddata(currSlice, &currSlice->rddata_bot_field_mb);
      }
    }

    if ((currSlice->slice_type == P_SLICE) || (currSlice->slice_type == SP_SLICE) || (currSlice->slice_type == B_SLICE))
    {
      free_mem_mv (currSlice->all_mv);
      if (p_Inp->BiPredMotionEstimation && (currSlice->slice_type == B_SLICE))
        free_mem_bipred_mv(currSlice->bipred_mv);

      if (currSlice->UseRDOQuant && currSlice->RDOQ_QP_Num > 1)
      {
        if (p_Inp->Transform8x8Mode && p_Inp->RDOQ_CP_MV)
        {
          free_mem4Dmv (currSlice->tmp_mv8);
          free_mem3Dint(currSlice->motion_cost8);
        }
      }
    }

    if (currSlice->slice_type == B_SLICE)
    {
      free_mem3D((byte ***)currSlice->direct_ref_idx);
      free_mem2D((byte **) currSlice->direct_pdir);
    }

    free_block_mem(currSlice);
    free(currSlice);
  }
}

static void set_ref_pic_num(Slice *currSlice)
{
  int i,j;
  StorablePicture *this_ref;
  ImageParameters *p_Img = currSlice->p_Img;

  //! need to add field ref_pic_num that handles field pair.

  for (i=0;i<p_Img->listXsize[LIST_0];i++)
  {
    this_ref = p_Img->listX[LIST_0][i];
    p_Img->enc_picture->ref_pic_num        [LIST_0][i] = this_ref->poc * 2 + ((this_ref->structure==BOTTOM_FIELD)?1:0) ;
    p_Img->enc_picture->frm_ref_pic_num    [LIST_0][i] = this_ref->frame_poc * 2;
    p_Img->enc_picture->top_ref_pic_num    [LIST_0][i] = this_ref->top_poc * 2;
    p_Img->enc_picture->bottom_ref_pic_num [LIST_0][i] = this_ref->bottom_poc * 2 + 1;
  }

  for (i=0;i<p_Img->listXsize[LIST_1];i++)
  {
    this_ref = p_Img->listX[LIST_1][i];
    p_Img->enc_picture->ref_pic_num        [LIST_1][i] = this_ref->poc * 2 + ((this_ref->structure==BOTTOM_FIELD)?1:0);
    p_Img->enc_picture->frm_ref_pic_num    [LIST_1][i] = this_ref->frame_poc * 2;
    p_Img->enc_picture->top_ref_pic_num    [LIST_1][i] = this_ref->top_poc * 2;
    p_Img->enc_picture->bottom_ref_pic_num [LIST_1][i] = this_ref->bottom_poc * 2 + 1;
  }

  if (!currSlice->active_sps->frame_mbs_only_flag && currSlice->structure==FRAME)
  {
    for (j=2;j<6;j++)
    {
      for (i=0;i<p_Img->listXsize[j];i++)
      {
        this_ref = p_Img->listX[j][i];
        p_Img->enc_picture->ref_pic_num[j][i] = this_ref->poc * 2 + ((this_ref->structure==BOTTOM_FIELD)?1:0);
        p_Img->enc_picture->frm_ref_pic_num[j][i] = this_ref->frame_poc * 2 ;
        p_Img->enc_picture->top_ref_pic_num[j][i] = this_ref->top_poc * 2 ;
        p_Img->enc_picture->bottom_ref_pic_num[j][i] = this_ref->bottom_poc * 2 + 1;
      }
    }
  }
}

void UpdateMELambda(ImageParameters *p_Img, InputParameters *p_Inp)
{
  int j, k, qp;
  if (p_Inp->UpdateLambdaChromaME)
  {
    for (j = 0; j < 6; j++)
    {
      for (qp = -p_Img->bitdepth_luma_qp_scale; qp < 52; qp++)
      { 
        for (k = 0; k < 3; k++)
        {
          if ((p_Inp->MEErrorMetric[k] == ERROR_SAD) && (p_Inp->ChromaMEEnable))
          {
            switch(p_Img->yuv_format)
            {
            case YUV420:
              p_Img->lambda_mf[j][qp][k] = (3 * p_Img->lambda_mf[j][qp][k] + 1) >> 1;
              p_Img->lambda_me[j][qp][k] *= 1.5;
              break;
            case YUV422:
              p_Img->lambda_mf[j][qp][k] *= 2;
              p_Img->lambda_me[j][qp][k] *= 2.0;
              break;
            case YUV444:
              p_Img->lambda_mf[j][qp][k] *= 3;
              p_Img->lambda_me[j][qp][k] *= 3.0;
              break;
            default:
              break;
            }
          }
        }
      }
    }
  }
}

void SetLambda(ImageParameters *p_Img, InputParameters *p_Inp, int j, int qp, double lambda_scale)
{
  int k;
  p_Img->lambda_md[j][qp] *= lambda_scale;

  for (k = F_PEL; k <= Q_PEL; k++)
  {
    p_Img->lambda_me[j][qp][k] =  (p_Inp->MEErrorMetric[k] == ERROR_SSE) ? p_Img->lambda_md[j][qp] : sqrt(p_Img->lambda_md[j][qp]);
    p_Img->lambda_mf[j][qp][k] = LAMBDA_FACTOR (p_Img->lambda_me[j][qp][k]);
  }
}

void SetLagrangianMultipliersOn(ImageParameters *p_Img, InputParameters *p_Inp)
{
  int qp, j;
  double qp_temp;
  double lambda_scale = 1.0 - dClip3(0.0,0.5,0.05 * (double) p_Inp->jumpd);

  if (p_Inp->UseExplicitLambdaParams == 1) // consideration of explicit lambda weights.
  {
    for (j = 0; j < 6; j++)
    {
      for (qp = -p_Img->bitdepth_luma_qp_scale; qp < 52; qp++)
      {
        qp_temp = (double)qp + p_Img->bitdepth_luma_qp_scale - SHIFT_QP;

        p_Img->lambda_md[j][qp] = p_Inp->LambdaWeight[j] * pow (2, qp_temp/3.0);
        SetLambda(p_Img, p_Inp, j, qp, ((p_Inp->MEErrorMetric[H_PEL] == ERROR_SATD && p_Inp->MEErrorMetric[Q_PEL] == ERROR_SATD) ? 1.00 : 0.95));

        if (j != B_SLICE)
          p_Img->lambda_md[j][qp] *= lambda_scale;

      }
    }
  }
  else if (p_Inp->UseExplicitLambdaParams == 2) // consideration of fixed lambda values.
  {
    for (j = 0; j < 6; j++)
    {
      for (qp = -p_Img->bitdepth_luma_qp_scale; qp < 52; qp++)
      {
        qp_temp = (double)qp + p_Img->bitdepth_luma_qp_scale - SHIFT_QP;

        p_Img->lambda_md[j][qp] = p_Inp->FixedLambda[j];
        SetLambda(p_Img, p_Inp, j, qp, ((p_Inp->MEErrorMetric[H_PEL] == ERROR_SATD && p_Inp->MEErrorMetric[Q_PEL] == ERROR_SATD) ? 1.00 : 0.95));
      }
    }
  }
  else
  {
    for (j = 0; j < 5; j++)
    {
      for (qp = -p_Img->bitdepth_luma_qp_scale; qp < 52; qp++)
      {
        qp_temp = (double)qp + p_Img->bitdepth_luma_qp_scale - SHIFT_QP;

        if(p_Inp->UseRDOQuant && p_Img->type == I_SLICE && p_Img->qp==qp)
          p_Img->lambda_md[j][qp] = 0.57 * pow (2.0, qp_temp/3.0); 
        else  if (p_Inp->NumberBFrames > 0)
        {
          /*
          if(j != B_SLICE)
               p_Img->lambda_md[j][qp] = 0.85 * pow (2.0, qp_temp/3.0) * ((j == SP_SLICE) ? dClip3(1.4, 3.0,(qp_temp / 12.0)) : 1.0);
          else
          */
          p_Img->lambda_md[j][qp] = 0.68 * pow (2.0, qp_temp/3.0)
          * (j == B_SLICE && p_Img->b_frame_to_code != 0 ? dClip3(2.00, 4.00, (qp_temp / 6.0)) : (j == SP_SLICE) ? dClip3(1.4,3.0,(qp_temp / 12.0)) : 1.0);
        }
        else
          p_Img->lambda_md[j][qp] = 0.85 * pow (2.0, qp_temp/3.0) * ((j == SP_SLICE) ? dClip3(1.4, 3.0,(qp_temp / 12.0)) : 1.0);
        
        // Scale lambda due to hadamard qpel only consideration
        p_Img->lambda_md[j][qp] = ((p_Inp->MEErrorMetric[H_PEL] == ERROR_SATD && p_Inp->MEErrorMetric[Q_PEL] == ERROR_SATD) ? 1.00 : 0.95) * p_Img->lambda_md[j][qp];

        if (j == B_SLICE)
        {
          p_Img->lambda_md[5][qp] = p_Img->lambda_md[j][qp];

          if (p_Img->b_frame_to_code != 0)
          {
            if (p_Inp->HierarchicalCoding == 2)
              p_Img->lambda_md[5][qp] *= (1.0 - dmin(0.4, 0.2 * (double) p_Img->gop_structure[p_Img->b_frame_to_code-1].hierarchy_layer));
            else
              p_Img->lambda_md[5][qp] *= 0.80;

          }
          SetLambda(p_Img, p_Inp, 5, qp, lambda_scale);
        }
        else
          p_Img->lambda_md[j][qp] *= lambda_scale;

        SetLambda(p_Img, p_Inp, j, qp, 1.0);

        if (p_Inp->CtxAdptLagrangeMult == 1)
        {
          int lambda_qp = (qp >= 32 && !p_Inp->RCEnable) ? imax(0, qp - 4) : imax(0, qp - 6);
          p_Img->lambda_mf_factor[j][qp] = log (p_Img->lambda_me[j][lambda_qp][Q_PEL] + 1.0) / log (2.0);
        }
      }
    }
  }

  UpdateMELambda(p_Img, p_Inp);
}


void SetLagrangianMultipliersOff(ImageParameters *p_Img, InputParameters *p_Inp)
{
  int qp, j, k;
  double qp_temp;

  for (j = 0; j < 6; j++)
  {
    for (qp = -p_Img->bitdepth_luma_qp_scale; qp < 52; qp++)
    {
      qp_temp = (double)qp + p_Img->bitdepth_luma_qp_scale - SHIFT_QP;

      switch (p_Inp->UseExplicitLambdaParams)
      {
      case 1:  // explicit lambda weights
        p_Img->lambda_md[j][qp] = sqrt(p_Inp->LambdaWeight[j] * pow (2, qp_temp/3.0));
        break;
      case 2: // explicit lambda
        p_Img->lambda_md[j][qp] = sqrt(p_Inp->FixedLambda[j]);
        break;
      default:
        p_Img->lambda_md[j][qp] = QP2QUANT[imax(0,qp - SHIFT_QP)];
        break;
      }

      for (k = F_PEL; k <= Q_PEL; k++)
      {
        p_Img->lambda_me[j][qp][k]  = (p_Inp->MEErrorMetric[k] == ERROR_SSE) ? (p_Img->lambda_md[j][qp] * p_Img->lambda_md[j][qp]) : p_Img->lambda_md[j][qp];
        p_Img->lambda_mf[j][qp][k]  = LAMBDA_FACTOR (p_Img->lambda_me[j][qp][k]);
      }

      if (p_Inp->CtxAdptLagrangeMult == 1)
      {
        int lambda_qp = (qp >= 32 && !p_Inp->RCEnable) ? imax(0, qp-4) : imax(0, qp-6);
        p_Img->lambda_mf_factor[j][qp] = log (p_Img->lambda_me[j][lambda_qp][Q_PEL] + 1.0) / log (2.0);
      }
    }
  }
  UpdateMELambda(p_Img, p_Inp);
}



