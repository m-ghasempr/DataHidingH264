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
 * \file header.c
 *
 * \brief
 *    H.26L Slice and Sequence headers
 *
 * \author
 *    Main contributors (see contributors.h for copyright, address and affiliation details)
 *      - Stephan Wenger                  <stewe@cs.tu-berlin.de>
 *************************************************************************************
 */

#include <math.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>

#include "elements.h"
#include "header.h"
#include "rtp.h"
#include "mbuffer.h"
#include "encodeiff.h"
#include "defines.h"

// A little trick to avoid those horrible #if TRACE all over the source code
#if TRACE
#define SYMTRACESTRING(s) strncpy(sym->tracestring,s,TRACESTRING_SIZE)
#else
#define SYMTRACESTRING(s) // to nothing
#endif

// local function
static int PutStartCode (Bitstream *s, int zeros_in_startcode);

/*!
 ********************************************************************************************
 * \brief 
 *    Write a slice header
 * \note
 *    We do not have a picture header. All neccessary information is coded at
 *    the slice level.
 ********************************************************************************************
*/
int SliceHeader()
{
  int dP_nr = assignSE2partition[input->partition_mode][SE_HEADER];
  Bitstream *currStream = ((img->currentSlice)->partArr[dP_nr]).bitstream;
  DataPartition *partition = &((img->currentSlice)->partArr[dP_nr]);
  SyntaxElement *sym;
  int len = 0;
  int Quant = img->qp;                  // Hack.  This should be a parameter as soon as Gisle is done
  int MBPosition = img->current_mb_nr;  // Hack.  This should be a paremeter as well.
  int zeros_in_startcode;

  if ((sym=(SyntaxElement*)calloc(1,sizeof(SyntaxElement)))==NULL) no_mem_exit("SliceHeader:sym");

  sym->type = SE_HEADER;                 // This will be true for all symbols generated here
  sym->mapping = n_linfo2;               // Mapping rule: Simple code number to len/info


  // Write a slice start code
  if(MBPosition == 0) //beginning of picture; so use longer start code
    zeros_in_startcode = ZEROBYTES_SHORTSTARTCODE + 1;
  else           // not beginning of picture; so use short start code
    zeros_in_startcode = ZEROBYTES_SHORTSTARTCODE;
  len += PutStartCode (currStream, zeros_in_startcode);
  startcodeprefix_len = (zeros_in_startcode + 1);


  // Now, take care of te UVLC coded data structure described in VCEG-M79
  sym->value1 = 0;               // TRType = 0
  SYMTRACESTRING("SH TemporalReferenceType");
  len += writeSyntaxElement_UVLC (sym, partition);

  sym->value1 = img->tr%256;         // TR, variable length

  SYMTRACESTRING("SH TemporalReference");
  len += writeSyntaxElement_UVLC (sym, partition);

  // Size information.  If the picture is Intra, then transmit the size in MBs,
  // For all other picture type just indicate that it didn't change.
  // Note: this is currently the prudent way to do things, because we do not
  // have anything similar to Annex P of H.263.  However, we should anticipate
  // that one day we may want to have an Annex P type functionality.  Hence, it is
  // unwise to assume that P pictures will never have a size indiciation.  So the
  // one bit for an "unchanged" inidication is well spent.

  if (img->type == INTRA_IMG||img->current_mb_nr==0)
 // if (img->type == INTRA_IMG)
  {
    // Double check that width and height are divisible by 16
    assert (img->width  % 16 == 0);
    assert (img->height % 16 == 0);

    sym->value1 = 1;             // SizeType = Width/Height in MBs
    SYMTRACESTRING("SH FullSizeInformation");
    len += writeSyntaxElement_UVLC (sym, partition);
    sym->value1 = img->width / 16;
    SYMTRACESTRING("SH FullSize-X");
    len += writeSyntaxElement_UVLC (sym, partition);
    sym->value1 = img->height / 16;
    SYMTRACESTRING("SH FullSize-Y");
    len += writeSyntaxElement_UVLC (sym, partition);

  } else
  {    // Not an intra frame -> write "unchanged"
    sym->value1 = 0;             // SizeType = Unchanged
    SYMTRACESTRING("SH SizeUnchanged");
    len += writeSyntaxElement_UVLC (sym, partition);
  }
  
  // AMT - This should follow SH FirstMBInSlice
  select_picture_type (sym);
  SYMTRACESTRING("SH Picture Type Symbol");
  len+= writeSyntaxElement_UVLC (sym, partition);

  // picture structure, (0: progressive, 1: top field, 2: bottom field, 
  // 3: top field first, 4: bottom field first)
  sym->value1 = img->pstruct;
  SYMTRACESTRING("SH Picture Stucture");
  len += writeSyntaxElement_UVLC (sym, partition);

  // Finally, write Reference Picture ID 
  // WYK: Oct. 16, 2001. Now I use this for the reference frame ID (non-B frame ID). 
  // If current frame is a B-frame, it will not change; otherwise it is increased by 1 based
  // on refPicID of the previous frame. Thus, the decoder can know how many  non-B 
  // frames are lost, and then can adjust the reference frame buffers correctly.
  if (1)
  {
    sym->value1 = img->refPicID%16;         // refPicID, variable length
 
    SYMTRACESTRING("SH RefPicID");
    len += writeSyntaxElement_UVLC (sym, partition);
  }

  // NOTE: following syntax elements have to be moved into nal_unit() and changed from e(v) to u(1)
  if(1)
  {
    sym->value1 = img->disposable_flag;
    SYMTRACESTRING("NAL disposable_flag");
    len += writeSyntaxElement_UVLC (sym, partition);
  }
  
  // NOTE: following syntax elements have to be moved into picture_layer_rbsp()
  // explicit_B_prediction_block_weight_indication
  if(img->type==B_IMG || img->type==BS_IMG)
  {
    sym->value1 = (input->BipredictiveWeighting > 0)? (input->BipredictiveWeighting-1) : 0;
    SYMTRACESTRING("PH weighted_bipred_implicit_flag");
    len += writeSyntaxElement_UVLC (sym, partition);
  }

  if (img->type==INTER_IMG || img->type==B_IMG || img->type==BS_IMG)
  {
    sym->value1 = img->num_ref_pic_active_fwd_minus1;
    SYMTRACESTRING("num_ref_pic_active_fwd_minus1");
    len += writeSyntaxElement_UVLC (sym, partition);
    if (img->type==B_IMG || img->type==BS_IMG)
    {
      sym->value1 = img->num_ref_pic_active_bwd_minus1;
      SYMTRACESTRING("num_ref_pic_active_bwd_minus1");
      len += writeSyntaxElement_UVLC (sym, partition);
    }
  }

  // For the GOB address, use straigtforward procedure.  Note that this allows slices 
  // to start at MB addresses up to 2^15 MBs due ot start code emulation problems.  
  // This should be enough for most practical applications,. but probably not for cinema. 
  // Has to be discussed. For the MPEG tests this is irrelevant, because there the rule 
  // will be one silce--one picture.

  // Put MB-Adresse
  assert (MBPosition < (1<<15));
  SYMTRACESTRING("SH FirstMBInSlice");
  sym->value1 = MBPosition;
  len += writeSyntaxElement_UVLC (sym, partition);

  // AMT - Direct Mode Type selection for B pictures
  if (img->type==B_IMG || img->type==BS_IMG)
  {
    SYMTRACESTRING("SH DirectSpatialFlag");
    sym->bitpattern = input->direct_type;  // 1 for Spatial Direct, 0 for temporal Direct
    sym->len = 1;
    len += writeSyntaxElement_fixed(sym, partition);
  }


#ifdef _ABT_FLAG_IN_SLICE_HEADER_
  // Indicate ABT mode to be applied.
  // 0 == ABT off
  // 1 == ABT inter
  // 2 == ABT inter and intra
  SYMTRACESTRING("SH ABTMode");
  sym->value1 = input->abt;
  len += writeSyntaxElement_UVLC (sym, partition);
#endif


  // Put Quant.  It's a bit irrationale that we still put the same quant here, but it's
  // a good provision for the future.  In real-world applications slices typically
  // start with Intra information, and Intra MBs will likely use a different quant
  // than Inter

  // Note:  Traditionally, the QP is a bit mask.  However, at numerically large QPs
  // we usually have high compression and don't want to waste bits, whereas
  // at low QPs this is not as much an issue.  Hence, the QUANT parameter
  // is coded as a UVLC calculated as 31 - QUANT.  That is, the UVLC representation
  // of 31-QUANT is coded, instead of QUANT.
  // Note 3:  In addition to the fields in VCEG-M79 there is one additional header field
  // with the MV resolution.  Currently defined values:
  // 0 == 1/4th pel resolution (old default)
  // 1 == 1/8th pel resolution
  // ... could be enhanced straightforward, it may make sense to define
  // 2 == 1/2 pel resolution
  // 3 == full pel resolution

  SYMTRACESTRING("SH SliceQuant");
  sym->mapping = dquant_linfo;  
  sym->value1 = Quant - (MAX_QP - MIN_QP +1)/2;
  len += writeSyntaxElement_UVLC (sym, partition);
  if (img->types==SP_IMG )  
  {
    if (img->types==SP_IMG && img->type!=INTRA_IMG) // Switch Flag only for SP pictures
    {
      SYMTRACESTRING("SH SWITCH FLAG");
      sym->bitpattern = 0;  // 1 for switching SP, 0 for normal SP
      sym->len = 1;
      len += writeSyntaxElement_fixed(sym, partition);
    }
    
    SYMTRACESTRING("SH SP SliceQuant");
    sym->value1 = img->qpsp - (MAX_QP - MIN_QP +1)/2;
    len += writeSyntaxElement_UVLC (sym, partition);
  }

  if (input->LFSendParameters)
  {
    SYMTRACESTRING("SH LF_DISABLE FLAG");
    sym->bitpattern = input->LFDisable;  /* Turn loop filter on/off on slice basis */
    sym->len = 1;
    len += writeSyntaxElement_fixed(sym, partition);

    if (!input->LFDisable)
    {
      sym->mapping = dquant_linfo;           // Mapping rule: Signed integer
      SYMTRACESTRING("SH LFAlphaC0OffsetDiv2");
      sym->value1 = input->LFAlphaC0Offset>>1; /* Convert from offset to code */
      len += writeSyntaxElement_UVLC (sym, partition);

      SYMTRACESTRING("SH LFBetaOffsetDiv2");
      sym->value1 = input->LFBetaOffset>>1; /* Convert from offset to code */
      len += writeSyntaxElement_UVLC (sym, partition);
    }
  }


  sym->mapping = n_linfo2;
  // Put the Motion Vector resolution as per reflector consensus
  SYMTRACESTRING("SH MVResolution");
  sym->value1 = input->mv_res;
  len += writeSyntaxElement_UVLC (sym, partition);
  
  len+=writeERPS(sym, partition);
  
  free(sym);

  return len;
}

/*!
 ********************************************************************************************
 * \brief 
 *    writes the ERPS syntax elements
 *
 * \return
 *    number of bits used for the ERPS
 ********************************************************************************************
*/
int writeERPS(SyntaxElement *sym, DataPartition *partition)
{
  int len=0;

#ifdef _CHECK_MULTI_BUFFER_1_
  RMPNIbuffer_t *r;
#endif

  /* RPSF: Reference Picture Selection Flags */
  SYMTRACESTRING("SH: Reference Picture Selection Flags");
  sym->value1 = 0;
  len += writeSyntaxElement_UVLC (sym, partition);

  /* PN: Picture Number */
  SYMTRACESTRING("SH: Picture Number");
  sym->value1 = img->pn;
#ifdef FLD0
  sym->value1 =0;
#endif
  len += writeSyntaxElement_UVLC (sym, partition);

#ifdef _CHECK_MULTI_BUFFER_1_

  /* RPSL: Reference Picture Selection Layer */
  SYMTRACESTRING("SH: Reference Picture Selection Layer");
  sym->value1 = 1;
  len += writeSyntaxElement_UVLC (sym, partition);

  if(img->type!=INTRA_IMG)
  {
    // let's mix some reference frames
    if ((img->pn==5)&&(img->type==INTER_IMG))
    {
      r = (RMPNIbuffer_t*)calloc (1,sizeof(RMPNIbuffer_t));
      r->RMPNI=0;
      r->Data=2;
      r->Next=NULL;
      img->currentSlice->rmpni_buffer=r;

    
      // negative ADPN follows
      SYMTRACESTRING("SH: RMPNI");
      sym->value1 = 0;
      len += writeSyntaxElement_UVLC (sym, partition);

      // ADPN
      SYMTRACESTRING("SH: ADPN");
      sym->value1 = 2;
      len += writeSyntaxElement_UVLC (sym, partition);

    }

    // End loop
    SYMTRACESTRING("SH: RMPNI");
    sym->value1 = 3;
    len += writeSyntaxElement_UVLC (sym, partition);
  }
  reorder_mref(img);
  
#else
  /* RPSL: Reference Picture Selection Layer */
  SYMTRACESTRING("SH: Reference Picture Selection Layer");
  sym->value1 = 0;
  len += writeSyntaxElement_UVLC (sym, partition);

#endif

#ifdef _CHECK_MULTI_BUFFER_2_

  SYMTRACESTRING("SH: Reference Picture Bufering Type");
  sym->value1 = 1;
  len += writeSyntaxElement_UVLC (sym, partition);


  // some code to check operation
  if ((img->pn==3) && (img->type==INTER_IMG))
  {

    // check in this frame as long term picture
//    if (img->max_lindex==0)
    {
      // set long term buffer size = 2
      SYMTRACESTRING("SH: MMCO Specify Max Long Term Index");
      // command
      sym->value1 = 4;
      len += writeSyntaxElement_UVLC (sym, partition);
      // size = 2+1 (MLP1)
      sym->value1 = 2+1;
      len += writeSyntaxElement_UVLC (sym, partition);

      img->max_lindex=2;
    }

    // assign a long term index to actual frame
    SYMTRACESTRING("SH: MMCO Assign Long Term Index to a Picture");
    // command
    sym->value1 = 3;
    len += writeSyntaxElement_UVLC (sym, partition);
    // DPN=0 for actual frame 
    sym->value1 = 0;
    len += writeSyntaxElement_UVLC (sym, partition);
    //long term ID
    sym->value1 = img->lindex;
    len += writeSyntaxElement_UVLC (sym, partition);

    // assign local long term
    init_long_term_buffer(2,img);
    init_mref(img);
    init_Refbuf(img);

    assign_long_term_id(3,img->lindex,img);
    
    img->lindex=(img->lindex+1)%img->max_lindex;


  } 
  if ((img->pn==4) && (img->type==INTER_IMG))
  {
    // delete long term picture again
    SYMTRACESTRING("SH: MMCO Mark a Long-Term Picture as Unused");
    // command
    sym->value1 = 2;
    len += writeSyntaxElement_UVLC (sym, partition);
    SYMTRACESTRING("SH: MMCO LPIN");
    // command
    sym->value1 = (img->max_lindex+img->lindex-1)%img->max_lindex;
    len += writeSyntaxElement_UVLC (sym, partition);
  } 

  // end MMCO loop
  SYMTRACESTRING("SH: end loop");
  sym->value1 = 0;
  len += writeSyntaxElement_UVLC (sym, partition);
#else
    /* RPBT: Reference Picture Bufering Type */
    SYMTRACESTRING("SH: Reference Picture Bufering Type");
    sym->value1 = 0;
    len += writeSyntaxElement_UVLC (sym, partition);
#endif 

  return len;
}


int SequenceHeader (FILE *outf)
{
  int LenInBytes = 0;
  int HeaderInfo;         // Binary coded Headerinfo to be written in file
  int ProfileLevelVersionHash;

  // Handle the RTP case
  if (input->of_mode == PAR_OF_RTP)
  {
    // Tian Dong: The following lines will never be reached. June 10, 2002
    if ((LenInBytes = RTPSequenceHeader (outf)) < 0)
    {
      snprintf (errortext, ET_SIZE, "SequenceHeaqder(): Problems writing the RTP Parameter Packet");
      error (errortext, 600);
      return -1;
    }
    else
      return LenInBytes*8;
  }
  // Non RTP-type file formats follow


  switch (input->SequenceHeaderType)
  {
  case 0:
    // No SequenceHeader, do nothing
    return 0;
  case 1:
    // A binary mini Sequence header.  Fixed length 32 bits.  Defined as follows
    // +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
    // |        0          |        1          |        2          | 3 |
    // +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
    // |0|1|2|3|4|5|6|7|8|9|0|1|2|3|4|5|6|7|8|9|0|1|2|3|4|5|6|7|8|9|0|1|
    // +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
    // |     TR Modulus        |     PicIDModulus      |OfMode |part |S|
    // +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
    //
    // S is symbol_mode (UVLC = 0, CABAC = 1), the other are as in the input file

    assert (input->TRModulus < 4096);
    assert (sizeof (int) == 4);
    assert (input->PicIdModulus <4096);
    assert (input->of_mode <16);
    assert (input->partition_mode <8);
    assert (input->symbol_mode < 2);

    ProfileLevelVersionHash = input->of_mode<<0 | input->partition_mode<<4 | input->symbol_mode<<7;
    HeaderInfo = input->TRModulus | input->PicIdModulus<<12 | ProfileLevelVersionHash<<24;

    if (1 != fwrite (&HeaderInfo, 4, 1, outf))
    {
      snprintf (errortext, ET_SIZE, "Error while writing Mini Sequence Header");
      error(errortext, 500);
    }
#if TRACE
    fprintf(p_trace, "Binary Mini Sequence Header 0x%x\n\n", HeaderInfo);
#endif

    return 32;

  case 2:
    // An ASCII representation of many input file parameters, useful for debug purposes
    //
    //
    assert ("To be hacked");
    return -1;
  case 3:
    // An UVLC coded representatioj of the Sequence Header.  To be hacked.
    // Note that all UVLC initialization has already taken place, so just use
    // write_synyaxelement() to put the sequence header symbols.  Has to be double
    // checked whether this is true for cabac as well.  Do we need a new context
    // for that?  Anyway:
    //
    assert ("to be hacked");
    return -1;
  default:
    snprintf (errortext, ET_SIZE, "Unspported Sequence Header Type (should not happen since checked in input module, exiting");
    error (errortext, 600);
    return -1;
  }
}

/********************************************************************************************
 ********************************************************************************************
 *
 * Local Support Functions
 *
 ********************************************************************************************
 ********************************************************************************************/


/*!
 ********************************************************************************************
 * \brief Puts the new Start Code into the Bitstream
 *    
 *
 * \return
 *    number of bits used for the Startcode
 *
 * \note Start-code is of the form N 0x00 bytes, followed by one 0x01 byte.
 *
 *  \param zeros_in_startcode indicates number of zero bytes in start-code
 *
 *  \note Start-code must be put in byte-aligned position
 ********************************************************************************************
*/
static int PutStartCode (Bitstream *s, int zeros_in_startcode)
{
  int i;
  if(s->bits_to_go != 8)
    printf(" Panic: Not byte aligned for putting new startcode - bits_to_go: %d\n", s->bits_to_go);  
  assert(s->bits_to_go==8);
  
  s->byte_buf = 0;
  for(i = 0; i < zeros_in_startcode; i++)
    s->streamBuffer[s->byte_pos++]=s->byte_buf;

  s->byte_buf = 1;
  s->streamBuffer[s->byte_pos++]=s->byte_buf;
  s->byte_buf = 0;
  return (8*zeros_in_startcode+8);
}


// StW Note: This function is a hack.  It would be cleaner if the encoder maintains
// the picture type in the given format.  Note further that I have yet to understand
// why the encoder needs to know whether a picture is predicted from one or more
// reference pictures.

/*!
 ************************************************************************
 * \brief
 *    Selects picture type and codes it to symbol
 ************************************************************************
 */
void select_picture_type(SyntaxElement *symbol)
{
  int multpred;

#ifdef _ADDITIONAL_REFERENCE_FRAME_
  if (input->no_multpred <= 1 && input->add_ref_frame == 0)
#else
  if (input->no_multpred <= 1)
#endif
    multpred=FALSE;
  else
    multpred=TRUE;               // multiple reference frames in motion search

  if (img->type == INTRA_IMG)
  {
    symbol->len=3;
    symbol->inf=1;
    symbol->value1 = 2;
  }
  else if((img->type == INTER_IMG) && (multpred == FALSE) ) // inter single reference frame
  {
    symbol->len=1;
    symbol->inf=0;
    symbol->value1 = 0;
  }
  else if((img->type == INTER_IMG) && (multpred == TRUE)) // inter multiple reference frames
  {
    symbol->len=3;
    symbol->inf=0;
    symbol->value1 = 1;
  }
  else if((img->type == B_IMG) && (multpred == FALSE))
  {
    symbol->len=5;
    symbol->inf=0;
    symbol->value1 = 3;
  }
  else if((img->type == B_IMG) && (multpred == TRUE))
  {
    symbol->len=5;
    symbol->inf=1;
    symbol->value1 = 4;
  }
  else if((img->type == BS_IMG) && (multpred == TRUE))
  {
    symbol->len=5;
    symbol->inf=1;
    symbol->value1 = 4;
  }
  else
  {
    error("Picture Type not supported!",1);
  }

  if((img->types == SP_IMG ) && (multpred == FALSE))
  {
    symbol->len=5;
    symbol->inf=0;
    symbol->value1 = 5;
  }
  else if((img->types == SP_IMG) && (multpred == TRUE))
  {
    symbol->len=5;
    symbol->inf=0;
    symbol->value1 = 6;
  }


#if TRACE
  snprintf(symbol->tracestring, TRACESTRING_SIZE, "Image type = %3d ", img->type);
#endif
  symbol->type = SE_PTYPE;
}



