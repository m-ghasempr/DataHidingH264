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
 *    H.26L Slice and Picture headers
 *
 *************************************************************************************
 */

#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include "global.h"
#include "elements.h"
#include "defines.h"
#include "fmo.h"


#if TRACE
#define SYMTRACESTRING(s) strncpy(sym.tracestring,s,TRACESTRING_SIZE)
#else
#define SYMTRACESTRING(s) // to nothing
#endif


/*!
 ************************************************************************
 * \brief
 *    read the whole Slice header without the Startcode UVLC
 ************************************************************************
 */
int SliceHeader(struct img_par *img, struct inp_par *inp)
{
  Slice *currSlice = img->currentSlice;
  int dP_nr = assignSE2partition[currSlice->dp_mode][SE_HEADER];
  DataPartition *partition = &(currSlice->partArr[dP_nr]);
  SyntaxElement sym;
  int UsedBits= partition->bitstream->frame_bitoffset; // was hardcoded to 31 for previous start-code. This is better.
  static int last_imgtr_frm=0,modulo_ctr_frm=0,last_imgtr_fld=0,modulo_ctr_fld=0;
  static int last_imgtr_frm_b=0,modulo_ctr_frm_b=0,last_imgtr_fld_b=0,modulo_ctr_fld_b=0;
  static int FirstCall = 1;   // used for FMO initialization

  int tmp1;
  RMPNIbuffer_t *tmp_rmpni,*tmp_rmpni2;
  MMCObuffer_t *tmp_mmco,*tmp_mmco2;
  int done;

  sym.type = SE_HEADER;
  sym.mapping = linfo;

  // 1. TRType = 0
  SYMTRACESTRING("SH TemporalReferenceType");
  readSyntaxElement_UVLC (&sym,img,inp,partition);
  // Currently, only the value 0 is supported, hence a simple assert to
  // catch evenntual  encoder problems
  if (sym.value1 != 0)
  {
    snprintf (errortext, ET_SIZE, "Unsupported value %d in Picture Header TemporalReferenceType, len %d info %d", sym.value1, sym.len, sym.inf);
    error(errortext, 600);
  }
  UsedBits += sym.len;

  // 2. TR, variable length
  SYMTRACESTRING("SH TemporalReference");
  readSyntaxElement_UVLC (&sym,img,inp,partition);
  currSlice->picture_id = img->tr = sym.value1;
  UsedBits += sym.len;

  // 3. Size Type
  SYMTRACESTRING ("SH WhichSize...");
  readSyntaxElement_UVLC (&sym,img,inp,partition);
  UsedBits += sym.len;
  if (sym.value1 == 0)
    SYMTRACESTRING("SH:    Size Information Unchanged");
  else
  {
    // width * 16
    SYMTRACESTRING("SH FullSize-X");
    readSyntaxElement_UVLC (&sym,img,inp,partition);
    img->width = sym.value1 * MB_BLOCK_SIZE;
    img->width_cr = img->width/2;
    UsedBits+=sym.len;
    SYMTRACESTRING("SH FullSize-Y");
    readSyntaxElement_UVLC (&sym,img,inp,partition);
    img->height = sym.value1 * MB_BLOCK_SIZE;
    img->height_cr = img->height/2;
    UsedBits+= sym.len;
  }

  // 4. Picture Type indication (I, P, Mult_P, B , Mult_B, SP, Mult_SP)
  SYMTRACESTRING("SH SliceType");
  readSyntaxElement_UVLC (&sym,img,inp,partition);
  currSlice->picture_type = img->type = sym.value1;
  UsedBits += sym.len;

  // Picture Structure indication (0 for frame, 1 for top field and 2 for bottom field)
  SYMTRACESTRING("SH PictureStructure");
  readSyntaxElement_UVLC (&sym,img,inp,partition);
  currSlice->structure = img->structure = sym.value1;
  UsedBits += sym.len;

  img->mb_frame_field_flag = 0;

  if (img->structure == 3)
    img->mb_frame_field_flag = 1;

  if (img->structure == 3)
    img->structure = 0;         // Frame coding with new scan for super MB pair

  // 5. Finally, read Reference Picture ID (same as TR here).  Note that this is an
  // optional field that is not present if the input parameters do not indicate
  // multiframe prediction ??

  // Ok, the above comment is nonsense.  There is no way how a decoder could
  // know that we use multiple reference frames (except probably through a
  // sequence header).  Hence, it's now an if (1) -- PHRefPicID is always present

  // Of course, the decoder can know. It is indicated by the previously decoded
  // parameter "PHPictureType". So, I changed the if-statement again to be
  // compatible with the encoder.

  // WYK: Oct. 16, 2001. Now I use this for the reference frame ID (non-B frame ID). 
  // Thus, we can know how many  non-B frames are lost, and then we can adjust 
  // the reference frame buffers correctly.
  if (1)
  {
    // refPicID, variable length
    SYMTRACESTRING("SH RefPicID");
    readSyntaxElement_UVLC (&sym,img,inp,partition);
    if (img->refPicID != sym.value1)
    {
      img->refPicID_old = img->refPicID;
      img->refPicID = sym.value1;
    }
    UsedBits += sym.len;
  }

  SYMTRACESTRING("disposable_flag");
  readSyntaxElement_UVLC (&sym,img,inp,partition);
  img->disposable_flag = sym.value1;
  UsedBits += sym.len;

  img->num_ref_pic_active_fwd = 0;
  img->num_ref_pic_active_bwd = 0;

  if(img->type==INTER_IMG_1 || img->type==INTER_IMG_MULT || img->type == SP_IMG_1 || img->type == SP_IMG_MULT)
  {
    SYMTRACESTRING("num_ref_pic_active_fwd_minus1");
    readSyntaxElement_UVLC (&sym,img,inp,partition);
    img->num_ref_pic_active_fwd = sym.value1+1;
    UsedBits += sym.len;    
  }
  else if(img->type==B_IMG_1 || img->type==B_IMG_MULT)
  {

    SYMTRACESTRING("explicit_bipred_weigt_indicator");
    readSyntaxElement_UVLC (&sym,img,inp,partition);
    img->explicit_B_prediction = sym.value1;
    UsedBits += sym.len;

    SYMTRACESTRING("num_ref_pic_active_fwd_minus1");
    readSyntaxElement_UVLC (&sym,img,inp,partition);
    img->num_ref_pic_active_fwd = sym.value1+1;
    UsedBits += sym.len;

    SYMTRACESTRING("num_ref_pic_active_bwd_minus1");
    readSyntaxElement_UVLC (&sym,img,inp,partition);
    img->num_ref_pic_active_bwd = sym.value1+1;
    UsedBits += sym.len;
  }

  if(!img->current_slice_nr)
  { 
    if (img->type <= INTRA_IMG || img->type >= SP_IMG_1 || !img->disposable_flag) 
    {
      if (img->structure == FRAME)
      {     
        if(img->tr <last_imgtr_frm) 
          modulo_ctr_frm++;
      
        last_imgtr_frm = img->tr;
        img->tr_frm = img->tr + (256*modulo_ctr_frm);
      }
      else
      {
        if(img->tr <last_imgtr_fld) 
          modulo_ctr_fld++;
      
        last_imgtr_fld = img->tr;
        img->tr_fld = img->tr + (256*modulo_ctr_fld);
      }
    }
    else
    {
      if (img->structure == FRAME)
      {     
        if(img->tr <last_imgtr_frm_b) 
          modulo_ctr_frm_b++;
      
        last_imgtr_frm_b = img->tr;
        img->tr_frm = img->tr + (256*modulo_ctr_frm_b);
      }
      else
      {
        if(img->tr <last_imgtr_fld_b) 
          modulo_ctr_fld_b++;
      
        last_imgtr_fld_b = img->tr;
        img->tr_fld = img->tr + (256*modulo_ctr_fld_b);
      }
    }
    if((img->type != B_IMG_MULT && img->type != B_IMG_1) || !img->disposable_flag) 
    {
      img->pstruct_next_P = img->structure;
      if(img->structure == TOP_FIELD)
      {
        img->imgtr_last_P = img->imgtr_next_P;
        img->imgtr_next_P = img->tr_fld;
      }
      else if(img->structure == FRAME)
      {
        img->imgtr_last_P = img->imgtr_next_P;
        img->imgtr_next_P = 2*img->tr_frm;
      }
    }
  }

  // 6. Get MB-Adresse
  SYMTRACESTRING("SH FirstMBInSlice");
  readSyntaxElement_UVLC (&sym,img,inp,partition);
  currSlice->start_mb_nr = sym.value1;
  UsedBits += sym.len;


  if(img->type==B_IMG_1 || img->type==B_IMG_MULT)
    {
      SYMTRACESTRING("SH DirectSpatialFlag");
      sym.len = 1;
      readSyntaxElement_fixed (&sym,img,inp,partition);
      img->direct_type = sym.inf;
    }

#ifdef _ABT_FLAG_IN_SLICE_HEADER_
  SYMTRACESTRING("SH ABTMode");
  readSyntaxElement_UVLC (&sym,img,inp,partition);
  currSlice->abt = sym.value1;
  UsedBits += sym.len;
#endif

  // 7. Get Quant.
  SYMTRACESTRING("SH SliceQuant");
  sym.mapping = linfo_dquant;
  readSyntaxElement_UVLC (&sym,img,inp,partition);
  currSlice->qp = img->qp = sym.value1 + (MAX_QP - MIN_QP + 1)/2;
  UsedBits += sym.len;
  if(img->type==SP_IMG_1 || img->type==SP_IMG_MULT || img->type == SI_IMG) 
  {
    if(img->type==SP_IMG_1 || img->type==SP_IMG_MULT)
    {
      SYMTRACESTRING("SH SP SWITCH");
      sym.len = 1;
      readSyntaxElement_fixed (&sym,img,inp,partition);
      img->sp_switch = sym.inf;      
    }
    if(img->type==SI_IMG)
      SYMTRACESTRING("SH SI SliceQuant");
    else
      SYMTRACESTRING("SH SP SliceQuant");
    readSyntaxElement_UVLC (&sym,img,inp,partition);
    img->qpsp = sym.value1 + (MAX_QP - MIN_QP + 1)/2;
  }

  if (inp->LFParametersFlag)
  {
    SYMTRACESTRING("SH LF_DISABLE FLAG");
    sym.len = 1;
    readSyntaxElement_fixed (&sym,img,inp,partition);
    currSlice->LFDisable = sym.inf;

    if (!currSlice->LFDisable)
    {
      sym.mapping = linfo_dquant;           // Mapping rule: Signed integer
      SYMTRACESTRING("SH LFAlphaC0OffsetDiv2");
      readSyntaxElement_UVLC (&sym,img,inp,partition);
      currSlice->LFAlphaC0Offset = sym.value1 << 1;
      UsedBits += sym.len;

      SYMTRACESTRING("SH LFBetaOffsetDiv2");
      readSyntaxElement_UVLC (&sym,img,inp,partition);
      currSlice->LFBetaOffset = sym.value1 << 1;
      UsedBits += sym.len;
    }
  }

  sym.mapping = linfo;
  // 8. Get MVResolution
  SYMTRACESTRING("SH MVResolution");
  readSyntaxElement_UVLC (&sym,img,inp,partition);
  img->mv_res = sym.value1;
  UsedBits += sym.len;

  /* KS: Multi-Picture Buffering Syntax */

  // Reference Picture Selection Flags 
  SYMTRACESTRING("SH Reference Picture Selection Flags");
  readSyntaxElement_UVLC (&sym,img,inp,partition);
  UsedBits += sym.len;

  // Picture Number 
  SYMTRACESTRING("SH Picture Number");
  readSyntaxElement_UVLC (&sym,img,inp,partition);
  img->pn=sym.value1;
  UsedBits += sym.len;

  // Reference picture selection layer 
  SYMTRACESTRING("SH Reference picture selection layer");
  readSyntaxElement_UVLC (&sym,img,inp,partition);
  UsedBits += sym.len;
  
  if (sym.value1)
  {
    // read Reference Picture Selection Layer 
    // free old buffer content
    while (img->currentSlice->rmpni_buffer)
    { 
      tmp_rmpni=img->currentSlice->rmpni_buffer;
 
      img->currentSlice->rmpni_buffer=tmp_rmpni->Next;
      free (tmp_rmpni);
    } 
    done=0;
    // if P or B frame RMPNI 

    if ((img->type==INTER_IMG_1)||(img->type==INTER_IMG_MULT)||(img->type==B_IMG_1)||
        (img->type==B_IMG_MULT)||(img->type==SP_IMG_1)||(img->type==SP_IMG_MULT))
    {
      do
      {
    
        SYMTRACESTRING("SH RMPNI");
        readSyntaxElement_UVLC (&sym,img,inp,partition);
        tmp1=sym.value1;
        UsedBits += sym.len;


        // check for illegal values
        if ((tmp1<0)||(tmp1>3))
          error ("Invalid RMPNI operation specified",400);

        if (tmp1!=3)
        {
          printf ("got RMPNI = %d\n",tmp1);
          tmp_rmpni=(RMPNIbuffer_t*)calloc (1,sizeof (RMPNIbuffer_t));
          tmp_rmpni->Next=NULL;
          tmp_rmpni->RMPNI=tmp1;

          // get the additional parameter
          SYMTRACESTRING("SH RMPNI Parameter");
          readSyntaxElement_UVLC (&sym,img,inp,partition);
          tmp_rmpni->Data=sym.value1;
          UsedBits += sym.len;

          // add RMPNI to list
          if (img->currentSlice->rmpni_buffer==NULL) 
          {
            img->currentSlice->rmpni_buffer=tmp_rmpni;
          }
          else
          {
            tmp_rmpni2=img->currentSlice->rmpni_buffer;
            while (tmp_rmpni2->Next!=NULL) 
              tmp_rmpni2=tmp_rmpni2->Next;
            tmp_rmpni2->Next=tmp_rmpni;
          }
        } else
        {
          // free temporary memory (no need to save end loop operation)
          done=1;
        }
      } while (!done);
    }
  }

  // RBPT 
  SYMTRACESTRING("SH RBPT");
  readSyntaxElement_UVLC (&sym,img,inp,partition);
  UsedBits += sym.len;

  // free old buffer content
  while (img->mmco_buffer)
  { 
    tmp_mmco=img->mmco_buffer;

    img->mmco_buffer=tmp_mmco->Next;
    free (tmp_mmco);
  } 

  if (sym.value1)
  {
    // read Memory Management Control Operation 
    do
    {

      tmp_mmco=(MMCObuffer_t*)calloc (1,sizeof (MMCObuffer_t));
      tmp_mmco->Next=NULL;
    
      SYMTRACESTRING("SH MMCO");
      readSyntaxElement_UVLC (&sym,img,inp,partition);
      tmp_mmco->MMCO=sym.value1;
      UsedBits += sym.len;

      switch (tmp_mmco->MMCO)
      {
      case 0:
      case 5:
        break;
      case 1:
        SYMTRACESTRING("SH DPN");
        readSyntaxElement_UVLC (&sym,img,inp,partition);
        tmp_mmco->DPN=sym.value1;
        UsedBits += sym.len;
        break;
      case 2:
        SYMTRACESTRING("SH LPIN");
        readSyntaxElement_UVLC (&sym,img,inp,partition);
        tmp_mmco->LPIN=sym.value1;
        UsedBits += sym.len;
        break;
      case 3:
        SYMTRACESTRING("SH DPN");
        readSyntaxElement_UVLC (&sym,img,inp,partition);
        tmp_mmco->DPN=sym.value1;
        UsedBits += sym.len;
        SYMTRACESTRING("SH LPIN");
        readSyntaxElement_UVLC (&sym,img,inp,partition);
        tmp_mmco->LPIN=sym.value1;
        UsedBits += sym.len;
        break;
      case 4:
        SYMTRACESTRING("SH MLIP1");
        readSyntaxElement_UVLC (&sym,img,inp,partition);
        tmp_mmco->MLIP1=sym.value1;
        UsedBits += sym.len;
        break;
      default:
        error ("Invalid MMCO operation specified",400);
        break;
      }

      // add MMCO to list
      if (img->mmco_buffer==NULL) 
      {
        img->mmco_buffer=tmp_mmco;
      }
      else
      {
        tmp_mmco2=img->mmco_buffer;
        while (tmp_mmco2->Next!=NULL) tmp_mmco2=tmp_mmco2->Next;
        tmp_mmco2->Next=tmp_mmco;
      }
      
    }while (tmp_mmco->MMCO!=0);
  }

  /* end KS */
  img->buf_cycle = inp->buf_cycle+1;
  img->pn=(((img->structure==BOTTOM_FIELD) ? (img->number/2):img->number)%img->buf_cycle);

  img->max_mb_nr = (img->width * img->height) / (MB_BLOCK_SIZE * MB_BLOCK_SIZE);


  //! Note! This software reqiures a valid MBAmap in any case.  When FMO is not
  //! used, the MBAmap consists of all zeros, yielding scan order slices. 
  //! 
  //! We cannot initialize the FMO module with the default (scan order MBAmap)
  //! before we know the picture size.  Hence, the FMO Module is initialized here,
  //! unless we have RTP file format.  The RTP NAL is currently the only one which
  //! implements the ParameterSets.  There, we can initialize the FMO module as
  //! soon as we got the Parameter Set packet.
  
  if (inp->of_mode != PAR_OF_RTP && FirstCall)
  {
    if (currSlice->structure)
    {
      // if first call is field: init with double height
      FmoInit (img->width/16, img->height/8, NULL, 0);   // force a default MBAmap
    }
    else
    {
      FmoInit (img->width/16, img->height/16, NULL, 0);   // force a default MBAmap
    }
    FirstCall = 0;
  }

  return UsedBits;
}

