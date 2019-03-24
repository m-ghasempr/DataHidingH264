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
 ************************************************************************
 * \file  rtp.c
 *
 * \brief
 *    Network Adaptation layer for RTP packets
 *
 * \author
 *    Main contributors (see contributors.h for copyright, address and affiliation details)
 *    - Stephan Wenger   <stewe@cs.tu-berlin.de>
 ************************************************************************
 */


/*!

  A quick guide to the basics of the RTP decoder implementation

  This module contains the RTP packetization, de-packetization, and the
  handling of Parameter Sets, see VCEG-N52 and accompanying documents.
  Note: Compound packets are not yet implemented!

  The interface between every NAL (including the RTP NAL) and the VCL is
  based on Slices.  The slice data structure on which the VCL is working
  is defined in the type Slice (in defines.h).  This type contains the
  various fields of the slice header and a partition array, which itself
  contains the data partitions the slice consists of.  When data
  partitioning is not used, then the whole slice bit string is stored
  in partition #0.  When individual partitions are missing, this is
  indicated by the size of the bit strings in the partition array.
  A complete missing slice (e.g. if a Full Slice packet was lost) is
  indicated in a similar way.  
  
  part of the slice structure is the error indication (ei-flag).  The
  Ei-flag is set in such cases in which at least one partition of a slice
  is damaged or missing.When data partitioning is used, it can happen that
  one partition does not contain any symbols but the ei_flag is cleared,
  which indicates the intentional missing of symbols of that partition.
  A typical example for this behaviour is the Intra Slice, which does not
  have symnbols in its type C partition.

  The VCL requests new data to work on through the call of readSliceRTP().
  This function calls the main state machine of this module in ReadRTPpaacket().

  ReadRTPpacket assumes, when called, that in an error free environment
  a complete slice, consisting of one Full Slice RTP packet, or three Partition
  packets of types A, B, C with consecutive sequence numbers, can be read.
  It first interprets any trailing SUPP and Parameter Update (Header) packets.
  Then it reads one video data packet.  Two cases have to be distinguished:

  1. Type A, or Full Slice packet
  In this case, the PictureID and the macroblock mumbers are used to
  identify the potential loss of a slice.  A slice is lost, when the
  StartMB of the newly read slice header is not equal to the current
  state of the decoder
    1.1 Loss detected
      In this case the last packet is unread (fseek back), and a dummy slice
      containing the missing macroblocks is conveyed to the VCL.  At the next 
      call of the NAL, the same packet is read again, but this time no packet 
      loss is detected by the above algorithm,
    1.2. No loss
      In this case it is checked whether a Full Slice packet or a type A data
      partition was read
        1.2.1 Full Slice
          The Full Slice packet is conveyed to the NAL
        1.2.2 Type A Partition
          The function RTPReadDataPartitionedSlice() is called, which collects
          the remaining type B, C partitions and handles them appropriately.

  Paraneter Update Packets (aka Header packets) are in an SDP-like syntax
  and are interpreted by a simple parser in the function 
  RTPInterpretParameterSetPacket() 

  Each Slice header contaions the information on which parameter set to be used.
  The function RTPSetImgInp() copies the information of the relevant parameter
  set in the VCL's global variables img-> and inp->  IMPORTANT: any changes
  in the semantics of the img-> and inp-> structure members must be represented
  in this function as well!

  A note to the stream-buffer data structure: The stream buffer always contains
  only the contents of the partition in question, and not the slice/partition
  header.  Decoding has to start at bitoffset 0 (UVLC) or bytreoffset 0 (CABAC).

  The remaining functions should be self-explanatory.
  
*/

#include "contributors.h"

#include <assert.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>
#include <string.h>
#include <ctype.h>

#include "global.h"
#include "errorconcealment.h"
#include "elements.h"
#include "bitsbuf.h"
#include "rtp.h"

extern void tracebits(const char *trace_str,  int len,  int info,int value1,
    int value2) ;

extern FILE *bits;
#define MAX_PARAMETER_STRINGLEN 1000

static int CurrentParameterSet = -1;

typedef struct
{
  int FramesToBeEncoded;
  int FrameSkip;
  char SequenceFileName[MAX_PARAMETER_STRINGLEN];
  int NumberBFrames;
} InfoSet_t;

static InfoSet_t InfoSet;

ParameterSet_t ParSet[RTP_MAX_PARAMETER_SET];

//! The following two variables are used to calculate the size of a lost slice.
//! The NAL conveyes this "lost" slice to the VCL with the ei_flag set in order
//! to trigger concealment

static int LastPicID;       //! PicID of the last packet (current state of the decoder 
                            //! before reading the next packet      
static int ExpectedMBNr;    //! MB Nr of the last decoded MB


/*!
 ************************************************************************
 * \brief
 *    read all Partitions of one Slice from RTP packets
 * \return
 *    SOS if this is not a new frame                                    \n
 *    SOP if this is a new frame                                        \n
 *    EOS if end of sequence reached
 ************************************************************************
 */
int readSliceRTP (struct img_par *img, struct inp_par *inp)
{
  Slice *currSlice = img->currentSlice;
  int PartitionMask;

  assert (currSlice != NULL);

  PartitionMask = ReadRTPPacket (img, inp, bits);
/*
{
int i;
for (i=0; i<25; i++)
printf ("%02x ", currSlice->partArr[0].bitstream->streamBuffer[i]);
printf ("\n");
for (i=0; i<25; i++)
printf ("%02x ", currSlice->partArr[1].bitstream->streamBuffer[i]);
printf ("\n");
for (i=0; i<25; i++)
printf ("%02x ", currSlice->partArr[2].bitstream->streamBuffer[i]);
printf ("\n");
}
*/
  if(PartitionMask == -4711)
    return EOS;

  if(currSlice->start_mb_nr != 0)
    return SOS;
  else
    return SOP;
}

  
/*!
 ************************************************************************
 * \brief
 *    read all partitions of one slice from RTP packet stream, also handle
 *    any Parameter Update packets and SuUPP-Packets
 * \return
 *    -1 for EOF                                              \n
 *    Partition-Bitmask otherwise:
 *    Partition Bitmask: 0x1: Full Slice, 0x2: type A, 0x4: Type B 
 *                       0x8: Type C
 ************************************************************************
 */
int ReadRTPPacket (struct img_par *img, struct inp_par *inp, FILE *bits)
{
  Slice *currSlice = img->currentSlice;
  DecodingEnvironmentPtr dep;
  byte *buf;
  static int first=1;
  RTPpacket_t *p, *nextp;
  RTPSliceHeader_t *sh, *nextsh;
  int MBDataIndex;
  int done = 0;
  int err, back=0;
  static int last_pframe=0, bframe_to_code=0;
  int b_interval;
  int FirstMacroblockInSlice;
  static int b_frame = FALSE;
  static int first_slice = FALSE;

  assert (currSlice != NULL);
  assert (bits != 0);

  // Tenporal storage for this function only

  p=alloca (sizeof (RTPpacket_t));      // the RTP packet
  p->packet=alloca (MAXRTPPACKETSIZE);
  p->payload=alloca (MAXRTPPAYLOADLEN);
  nextp=alloca (sizeof (RTPpacket_t));  
  sh=alloca(sizeof (RTPSliceHeader_t));
  sh->RMPNIbuffer=NULL;
  sh->MMCObuffer=NULL;

  nextsh=alloca(sizeof (RTPSliceHeader_t));
  nextsh->RMPNIbuffer=NULL;
  nextsh->MMCObuffer=NULL;

  ExpectedMBNr = img->current_mb_nr;
  LastPicID = img->tr;

  done = 0;
  do  
  {
//    Filepos = ftell (bits);             // to be able to go back one packet
      
    if (RTPReadPacket (p, bits) < 0)    // Read and decompose
      return -4711;

    switch (p->payload[0] & 0xf)
    {
    case 0:       // Full Slice packet
    case 1:       // Partition type A packet
      done = 1;
      break;
 
    case 2:       // Partition B
    case 3:       // Partition C
      // Do nothing.  this results in discarding the unexpected Partition B, C
      // packet, which will later be concealed "automatically", because when
      // interpreting the next Partition A or full slice packet a range of
      // lost blocks is found which will be concealed in the usual manner.
      //
      // If anyonme comes up with an idea how to use coefficients without the
      // header information then this code has to be changed
      printf ("found unexpected Partition %c packet, skipping\n", (p->payload[0]&0xf)==2?'B':'C');
      break;
    case 4:
      //! Compound packets may be handled here, but I (StW) would personally prefer and
      //! recommend to handle them externally, by a pre-processor tool.  For now,
      //! compounds lead to exit()
      printf ("Compound packets not yet implemented, exit\n");
      exit (-700);
      break;
    case 5:
      //! Add implementation of SUPP packets here
      printf ("SUPP packets not yet implemented, skipped\n");
      break;
    case 6:
      printf ("Found header packet\n");

      if ((err = RTPInterpretParameterSetPacket (&p->payload[1], p->paylen-1)) < 0)
      {
        printf ("RTPInterpretParameterSetPacket returns error %d\n", err);
      }
      break;
    default:
      printf ("Undefined packet type %d found, skipped\n", p->payload[0] & 0xf);
      assert (0==1);
      break;
    }
  } while (!done);

  // Here, all the non-video data packets and lonely type B, C partitions
  // are handled.  Now work on the expected type A and full slice packets

  assert ((p->payload[0] & 0xf) < 2);

  if ((p->payload[0] & 0xf) == 0)       // Full Slice packet
  {
    currSlice->ei_flag = 0;
    MBDataIndex = 1;                    // Skip First Byte
    MBDataIndex += RTPInterpretSliceHeader (&p->payload[1], p->paylen-1, 0, sh);
  }
  else                                  // Partition A packet
  {
    currSlice->ei_flag = 0;
    MBDataIndex = 1;                    // Skip First Byte
    MBDataIndex += RTPInterpretSliceHeader (&p->payload[1], p->paylen-1, 1, sh);
  }


  FirstMacroblockInSlice = sh->FirstMBInSliceY * (img->width/16) + 
                               sh->FirstMBInSliceX;   //! assumes picture sizes divisble by 16
  // The purpose of the following if cascade is to check for a lost
  // segment  of macroblocks.  if such a segment is found, a dummy slice
  // without content, but with ei_flag set is generated in order to trigger
  // concealment.
  if(first)
  {
    first = FALSE;
    bframe_to_code = InfoSet.NumberBFrames+1;
    currSlice->max_part_nr = MAX_PART_NR;   // Just to get a start
    ExpectedMBNr = 0;
    LastPicID = -1;
  }
  // slightly different handling for first slice that is received
  // necessary for the case we loose this first slice or even the whole first frame
  if(first_slice == FALSE)
  {
    first_slice = TRUE;
    b_frame = FALSE;
    img->pn = 0;
    if(sh->PictureID == 0)
    {
      if(FirstMacroblockInSlice == 0)
        currSlice->partArr[0].bitstream->ei_flag = 0; // everything seems to be ok.
      else
      {
        back=p->packlen+8;                // Unread the packet
        fseek (bits, -back, SEEK_CUR);    
        currSlice->ei_flag = 1;
        currSlice->dp_mode = PAR_DP_1;
        currSlice->start_mb_nr = ExpectedMBNr;
        currSlice->next_header = RTPGetFollowingSliceHeader (img, nextp, nextsh);
          
        img->tr = currSlice->picture_id = last_pframe = sh->PictureID;
        img->type = INTRA_IMG;
        img->pn = 0;
        bframe_to_code = 1;
 
        currSlice->start_mb_nr = 0;
#if _ERROR_CONCEALMENT_
        currSlice->last_mb_nr = FirstMacroblockInSlice-1;
#else
        currSlice->last_mb_nr = FirstMacroblockInSlice;
#endif
        currSlice->partArr[0].bitstream->bitstream_length=0;
        currSlice->partArr[0].bitstream->read_len=0;
        currSlice->partArr[0].bitstream->code_len=0;
        currSlice->partArr[0].readSyntaxElement = readSyntaxElement_RTP;

        return 0;
      }
    }
    else
    {
      back=p->packlen+8;                // Unread the packet
      fseek (bits, -back, SEEK_CUR);    
      currSlice->ei_flag = 1;
      currSlice->dp_mode = PAR_DP_1;
      currSlice->start_mb_nr = ExpectedMBNr;
      currSlice->next_header = RTPGetFollowingSliceHeader (img, nextp, nextsh);
          
      img->tr = currSlice->picture_id = last_pframe = 0;
      img->type = INTRA_IMG;
      img->pn = 0;
      bframe_to_code = 1;

      img->max_mb_nr = (img->width * img->height) / (MB_BLOCK_SIZE * MB_BLOCK_SIZE);
#if _ERROR_CONCEALMENT_
      currSlice->last_mb_nr = img->max_mb_nr-1;
#else
      currSlice->last_mb_nr = img->max_mb_nr;
#endif
      currSlice->partArr[0].bitstream->bitstream_length=0;
      currSlice->partArr[0].bitstream->read_len=0;
      currSlice->partArr[0].bitstream->code_len=0;
      currSlice->partArr[0].readSyntaxElement = readSyntaxElement_RTP;

      return 0;
    }
  }
  else
  {
    if (LastPicID == sh->PictureID)   // we are in the same picture
    {
      if (ExpectedMBNr == FirstMacroblockInSlice)
      {
        currSlice->partArr[0].bitstream->ei_flag = 0;    // everything seems to be ok.
      }
      else
      {
        if (FirstMacroblockInSlice == 0)
        {
          assert ("weird! PicID wrap around?  Should not happen\n");
        }
        else
        {
          //printf ("SLICE LOSS 1: Slice loss at PicID %d, macoblocks %d to %d\n",LastPicID, ExpectedMBNr, FirstMacroblockInSlice-1);
          back=p->packlen+8;
          fseek (bits, -back, SEEK_CUR);    

          //FirstMacroblockInSlice = (img->width*img->height)/(16*16);
          currSlice->ei_flag = 1;
          currSlice->dp_mode = PAR_DP_1;
          currSlice->start_mb_nr = ExpectedMBNr;
          //! just the same slice we just read!
          currSlice->next_header = RTPGetFollowingSliceHeader (img, nextp, nextsh); 
          assert (currSlice->start_mb_nr == img->current_mb_nr); 
  #if _ERROR_CONCEALMENT_
          currSlice->last_mb_nr = FirstMacroblockInSlice-1;
  #else
          currSlice->last_mb_nr = FirstMacroblockInSlice;
  #endif
          currSlice->partArr[0].bitstream->bitstream_length=0;
          currSlice->partArr[0].bitstream->read_len=0;
          currSlice->partArr[0].bitstream->code_len=0;
          currSlice->partArr[0].bitstream->ei_flag=1;
          
          img->tr = currSlice->picture_id = LastPicID;
         
          return 0;
        }
      }
    }
    else      // we are in a different picture
    {
      b_interval = (int)((float)(InfoSet.FrameSkip +1)/(float)(InfoSet.NumberBFrames +1) + 0.49999);
      if (ExpectedMBNr == 0)    // the old picture was finished
      {
        if (((last_pframe + InfoSet.FrameSkip +1)%256) == sh->PictureID && 
             (bframe_to_code > InfoSet.NumberBFrames)) //! we received a new P-Frame and coded all B-Frames
        {
          if(InfoSet.NumberBFrames)
          {
            bframe_to_code = 1;
            last_pframe = sh->PictureID-(InfoSet.FrameSkip +1);
          }
          else
            last_pframe = sh->PictureID;

          if(last_pframe < 0)
            last_pframe += 256;
          b_frame = FALSE;
        }
        else if(sh->PictureID == ((last_pframe + b_interval*bframe_to_code)%256) && InfoSet.NumberBFrames) //! we received a B-Frame
        {
          bframe_to_code ++;
          b_frame = TRUE;
          if(bframe_to_code > InfoSet.NumberBFrames)
            last_pframe += (InfoSet.FrameSkip +1);
        }
        else //! we lost at least one whole frame
        {
          //printf ("SLICE LOSS 4: Slice loss at PicID %d containing the whole frame\n", LastPicID + InfoSet.FrameSkip +1); 
          back=p->packlen+8;                // Unread the packet
          fseek (bits, -back, SEEK_CUR);    
          currSlice->ei_flag = 1;
          currSlice->dp_mode = PAR_DP_1;
          currSlice->start_mb_nr = ExpectedMBNr;
          currSlice->next_header = RTPGetFollowingSliceHeader (img, nextp, nextsh); 
          currSlice->partArr[0].readSyntaxElement = readSyntaxElement_RTP;
          assert (currSlice->start_mb_nr == img->current_mb_nr); 

  #if _ERROR_CONCEALMENT_
          currSlice->last_mb_nr = img->max_mb_nr-1;
  #else
          currSlice->last_mb_nr = img->max_mb_nr;
  #endif
          currSlice->partArr[0].bitstream->bitstream_length=0;
          currSlice->partArr[0].bitstream->read_len=0;
          currSlice->partArr[0].bitstream->code_len=0;

          if(!InfoSet.NumberBFrames || (bframe_to_code > InfoSet.NumberBFrames)) //! we expect a P-Frame
          {
            img->tr = currSlice->picture_id = (last_pframe + InfoSet.FrameSkip +1)%256;
            img->type = INTER_IMG_1;
            if(InfoSet.NumberBFrames)
            {
              last_pframe = img->tr-(InfoSet.FrameSkip +1);
              bframe_to_code = 1;
            }
            else
              last_pframe = img->tr;
          
            if(last_pframe < 0)
              last_pframe += 256;
            b_frame = FALSE;
            img->pn++;
          }
          else //! we expect a B-Frame
          {
            img->tr = currSlice->picture_id = (last_pframe + b_interval*bframe_to_code)%256;
            img->type = B_IMG_1;
            bframe_to_code ++;
            if(bframe_to_code > InfoSet.NumberBFrames)
              last_pframe += (InfoSet.FrameSkip +1);
          }
        
          return 0;
        }
        if(FirstMacroblockInSlice == 0)
        {
          currSlice->partArr[0].bitstream->ei_flag = 0; // everything seems to be ok.
        }
        else //! slice loss from the begining of the frame
        {
          //printf ("SLICE LOSS 2: Slice loss at PicID %d, beginning of picture to macroblock %d\n", LastPicID, FirstMacroblockInSlice); 
          back=p->packlen+8;                // Unread the packet
          fseek (bits, -back, SEEK_CUR);    
          currSlice->ei_flag = 1;
          currSlice->dp_mode = PAR_DP_1;
          currSlice->start_mb_nr = ExpectedMBNr;
          currSlice->next_header = RTPGetFollowingSliceHeader (img, nextp, nextsh); 
          assert (currSlice->start_mb_nr == img->current_mb_nr); 

  #if _ERROR_CONCEALMENT_
          currSlice->last_mb_nr = FirstMacroblockInSlice-1;
  #else
          currSlice->last_mb_nr = FirstMacroblockInSlice;
  #endif
          currSlice->partArr[0].bitstream->bitstream_length=0;
          currSlice->partArr[0].bitstream->read_len=0;
          currSlice->partArr[0].bitstream->code_len=0;
          currSlice->partArr[0].readSyntaxElement = readSyntaxElement_RTP;
        
          img->tr = currSlice->picture_id = sh->PictureID;
          if(b_frame == FALSE) //! we expect a P-Frame
          {
            img->type = INTER_IMG_1;
            img->pn++;
          }  
          else
            img->type = B_IMG_1;
        
          return 0;
        }
      }
      else //we did not finish the old frame
      {
        //upprintf ("SLICE LOSS 3: Slice loss at PicID %d, macroblocks %d to end of picture\n", LastPicID, ExpectedMBNr); 
        back=p->packlen+8;                // Unread the packet
        fseek (bits, -back, SEEK_CUR);    
        currSlice->ei_flag = 1;
        currSlice->dp_mode = PAR_DP_1;
        currSlice->start_mb_nr = ExpectedMBNr;
        currSlice->next_header = RTPGetFollowingSliceHeader (img, nextp, nextsh); 
        assert (currSlice->start_mb_nr == img->current_mb_nr); 
  #if _ERROR_CONCEALMENT_
        currSlice->last_mb_nr = img->max_mb_nr-1;
  #else
        currSlice->last_mb_nr = img->max_mb_nr;
  #endif
        currSlice->partArr[0].bitstream->bitstream_length=0;
        currSlice->partArr[0].bitstream->read_len=0;
        currSlice->partArr[0].bitstream->code_len=0;
        currSlice->partArr[0].readSyntaxElement = readSyntaxElement_RTP;
      
        img->tr = currSlice->picture_id = LastPicID;
      
        return 0;
      }
    }
  }

  // Here, all concealment is done and we have either a type A partition 
  // packet or a full slice packet, which need to be worked on
    
  RTPUseParameterSet (sh->ParameterSet, img, inp);
  RTPSetImgInp(img, inp, sh);

  free_Partition (currSlice->partArr[0].bitstream);

  assert (p->paylen-MBDataIndex >= 0);
      
  currSlice->partArr[0].bitstream->read_len = 0;
  currSlice->partArr[0].bitstream->code_len = p->paylen-MBDataIndex;        // neu
  currSlice->partArr[0].bitstream->bitstream_length = p->paylen-MBDataIndex;

  memcpy (currSlice->partArr[0].bitstream->streamBuffer, &p->payload[MBDataIndex],p->paylen-MBDataIndex);
  buf = currSlice->partArr[0].bitstream->streamBuffer;

  if(inp->symbol_mode == CABAC)
  {
    dep = &((currSlice->partArr[0]).de_cabac);
    arideco_start_decoding(dep, buf, 0, &currSlice->partArr[0].bitstream->read_len);
  }
      
  currSlice->next_header = RTPGetFollowingSliceHeader (img, nextp, nextsh); // no use for the info in nextp, nextsh yet. 
  
  if ((p->payload[0]&0xf) == 0  || b_frame)         // Full Slice Packet or B-Frame
  {
    currSlice->dp_mode = PAR_DP_1;
    currSlice->max_part_nr=1;
    return 1;
  }
  else
  {
    currSlice->dp_mode = PAR_DP_3;
    currSlice->max_part_nr = 3;
    RTPProcessDataPartitionedSlice (img, inp, bits, p, sh->SliceID);
    return 3;

  }

  return FALSE;
}


/*!
 *****************************************************************************
 *
 * \brief 
 *    DecomposeRTPpacket interprets the RTP packet and writes the various
 *    structure members of the RTPpacket_t structure
 *
 * \return
 *    0 in case of success
 *    negative error code in case of failure
 *
 * \para Parameters
 *    Caller is responsible to allocate enough memory for the generated payload
 *    in parameter->payload. Typically a malloc of paclen-12 bytes is sufficient
 *
 * \para Side effects
 *    none
 *
 * \para Other Notes
 *    Function contains assert() tests for debug purposes (consistency checks
 *    for RTP header fields)
 *
 * \date
 *    30 Spetember 2001
 *
 * \author
 *    Stephan Wenger   stewe@cs.tu-berlin.de
 *****************************************************************************/

int DecomposeRTPpacket (RTPpacket_t *p)

{
  // consistency check 
  assert (p->packlen < 65536 - 28);  // IP, UDP headers
  assert (p->packlen >= 12);         // at least a complete RTP header
  assert (p->payload != NULL);
  assert (p->packet != NULL);

  // Extract header information

  p->v  = p->packet[0] & 0x3;
  p->p  = (p->packet[0] & 0x4) >> 2;
  p->x  = (p->packet[0] & 0x8) >> 3;
  p->cc = (p->packet[0] & 0xf0) >> 4;

  p->m  = p->packet[1] & 0x1;
  p->pt = (p->packet[1] & 0xfe) >> 1;

  p->seq = p->packet[2] | (p->packet[3] << 8);

  memcpy (&p->timestamp, &p->packet[4], 4);// change to shifts for unified byte sex
  memcpy (&p->ssrc, &p->packet[8], 4);// change to shifts for unified byte sex

  // header consistency checks
  if (     (p->v != 2)
        || (p->p != 0)
        || (p->x != 0)
        || (p->cc != 0) )
  {
    printf ("DecomposeRTPpacket, RTP header consistency problem, header follows\n");
    DumpRTPHeader (p);
    return -1;
  }
  p->paylen = p->packlen-12;
  memcpy (p->payload, &p->packet[12], p->paylen);
  return 0;
}

/*!
 *****************************************************************************
 *
 * \brief 
 *    DumpRTPHeader is a debug tool that dumps a human-readable interpretation
 *    of the RTP header
 *
 * \return
 *    n.a.
 * \para Parameters
 *    the RTP packet to be dumped, after DecompositeRTPpacket()
 *
 * \para Side effects
 *    Debug output to stdout
 *
 * \date
 *    30 Spetember 2001
 *
 * \author
 *    Stephan Wenger   stewe@cs.tu-berlin.de
 *****************************************************************************/

void DumpRTPHeader (RTPpacket_t *p)

{
  int i;
  for (i=0; i< 30; i++)
    printf ("%02x ", p->packet[i]);
  printf ("Version (V): %d\n", p->v);
  printf ("Padding (P): %d\n", p->p);
  printf ("Extension (X): %d\n", p->x);
  printf ("CSRC count (CC): %d\n", p->cc);
  printf ("Marker bit (M): %d\n", p->m);
  printf ("Payload Type (PT): %d\n", p->pt);
  printf ("Sequence Number: %d\n", p->seq);
  printf ("Timestamp: %d\n", p->timestamp);
  printf ("SSRC: %d\n", p->ssrc);
}

/*!
 *****************************************************************************
 *
 * \brief 
 *    Parses and interprets the UVLC-coded slice header
 *
 * \return
 *    negative in case of errors, the byte-index where the UVLC/CABAC MB data
 *    starts otherwise
 *
 * \date
 *    27 October, 2001
 *
 * \author
 *    Stephan Wenger   stewe@cs.tu-berlin.de
 *****************************************************************************/

int RTPInterpretSliceHeader (byte *buf, int bufsize, int ReadSliceId, RTPSliceHeader_t *sh)
{
  int len, info, bytes, dummy, bitptr=0;
  int temp, tmp1;
  RMPNIbuffer_t *tmp_rmpni,*tmp_rmpni2;
  MMCObuffer_t *tmp_mmco,*tmp_mmco2;
  int done;
  
  len = GetVLCSymbol(buf, bitptr, &info, bufsize);
  linfo (len, info, &sh->ParameterSet, &dummy);
  bitptr+=len;

  len = GetVLCSymbol(buf, bitptr, &info, bufsize);
  linfo (len, info, &sh->structure, &dummy);
  bitptr+=len;

  len = GetVLCSymbol(buf, bitptr, &info, bufsize);
  linfo (len, info, &sh->PictureID, &dummy);
  bitptr+=len;

  len = GetVLCSymbol(buf, bitptr, &info, bufsize);
  linfo (len, info, &sh->SliceType, &dummy);
  bitptr+=len;

  len = GetVLCSymbol(buf, bitptr, &info, bufsize);
  linfo (len, info, &sh->FirstMBInSliceX, &dummy);
  bitptr+=len;

  len = GetVLCSymbol(buf, bitptr, &info, bufsize);
  linfo (len, info, &sh->FirstMBInSliceY, &dummy);
  bitptr+=len;

  len = GetVLCSymbol(buf, bitptr, &info, bufsize);
  linfo (len, info, &sh->InitialQP, &dummy);
  bitptr+=len;
  sh->InitialQP = MAX_QP-sh->InitialQP;

  if (sh->SliceType==2) // SP Picture
  {
    len = GetVLCSymbol(buf, bitptr, &info, bufsize);
    linfo (len, info, &sh->InitialSPQP, &dummy);
    bitptr+=len;
    sh->InitialSPQP = MAX_QP-sh->InitialSPQP;
    assert (sh->InitialSPQP >=MIN_QP && sh->InitialSPQP <= MAX_QP);
  }

  assert (sh->ParameterSet == 0);     // only for testing, should be deleted as soon as more than one parameter set is generated by trhe encoder
  assert (sh->SliceType > 0 || sh->SliceType < 5);
  assert (sh->InitialQP >=MIN_QP && sh->InitialQP <= MAX_QP);


  if (ReadSliceId)
  {
    len = GetVLCSymbol(buf, bitptr, &info, bufsize);
    linfo (len, info, &sh->SliceID, &dummy);
    bitptr+=len;
  }

  /* KS: Multi-Picture Buffering Syntax */

  /* Reference Picture Selection Flags */
  len = GetVLCSymbol(buf, bitptr, &info, bufsize);
  linfo (len, info, &temp, &dummy);
  bitptr+=len;

  /* Picture Number */
  len = GetVLCSymbol(buf, bitptr, &info, bufsize);
  linfo (len, info, &sh->PictureNum, &dummy);
  bitptr+=len;

  /* Reference picture selection layer */
  len = GetVLCSymbol(buf, bitptr, &info, bufsize);
  linfo (len, info, &temp, &dummy);
  bitptr+=len;

  if (temp)
  {
    /* read Reference Picture Selection Layer */
    // free old buffer content
    while (sh->RMPNIbuffer)
    { 
      tmp_rmpni=sh->RMPNIbuffer;
 
      sh->RMPNIbuffer=tmp_rmpni->Next;
      free (tmp_rmpni);
    } 
    done=0;
    /* if P or B frame RMPNI */

    if ((sh->SliceType>=0)&&(sh->SliceType<=2))
    {
      do
      {
    
        len = GetVLCSymbol(buf, bitptr, &info, bufsize);
        linfo (len, info, &tmp1, &dummy);
        bitptr+=len;


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
          len = GetVLCSymbol(buf, bitptr, &info, bufsize);
          linfo (len, info, &tmp_rmpni->Data, &dummy);
          bitptr+=len;

          // add RMPNI to list
          if (sh->RMPNIbuffer==NULL) 
          {
            sh->RMPNIbuffer=tmp_rmpni;
          }
          else
          {
            tmp_rmpni2=sh->RMPNIbuffer;
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

  len = GetVLCSymbol(buf, bitptr, &info, bufsize);
  linfo (len, info, &sh->RPBT, &dummy);
  bitptr+=len;

  // free old buffer content
  while (sh->MMCObuffer)
  { 
    tmp_mmco=sh->MMCObuffer;
    sh->MMCObuffer=tmp_mmco->Next;
    free (tmp_mmco);
  } 
  
  /* read Memory Management Control Operation */
  if (sh->RPBT)
  {
    do
    {

      tmp_mmco=(MMCObuffer_t*)calloc (1,sizeof (MMCObuffer_t));
      tmp_mmco->Next=NULL;
    
      len = GetVLCSymbol(buf, bitptr, &info, bufsize);
      linfo (len, info, &tmp_mmco->MMCO, &dummy);
      bitptr+=len;

      switch (tmp_mmco->MMCO)
      {
      case 0:
      case 5:
        break;
      case 1:
        len = GetVLCSymbol(buf, bitptr, &info, bufsize);
        linfo (len, info, &tmp_mmco->DPN, &dummy);
        bitptr+=len;
        break;
      case 2:
        len = GetVLCSymbol(buf, bitptr, &info, bufsize);
        linfo (len, info, &tmp_mmco->LPIN, &dummy);
        bitptr+=len;
        break;
      case 3:
        len = GetVLCSymbol(buf, bitptr, &info, bufsize);
        linfo (len, info, &tmp_mmco->DPN, &dummy);
        bitptr+=len;
        len = GetVLCSymbol(buf, bitptr, &info, bufsize);
        linfo (len, info, &tmp_mmco->LPIN, &dummy);
        bitptr+=len;
        break;
      case 4:
        len = GetVLCSymbol(buf, bitptr, &info, bufsize);
        linfo (len, info, &tmp_mmco->MLIP1, &dummy);
        bitptr+=len;
        break;
      default:
        error ("Invalid MMCO operation specified",400);
        break;
      }

      // add MMCO to list
      if (sh->MMCObuffer==NULL) 
      {
        sh->MMCObuffer=tmp_mmco;
      }
      else
      {
        tmp_mmco2=sh->MMCObuffer;
        while (tmp_mmco2->Next!=NULL) tmp_mmco2=tmp_mmco2->Next;
        tmp_mmco2->Next=tmp_mmco;
      }
      
    }while (tmp_mmco->MMCO!=0);
  }
  /* end KS */

  if (ParSet[sh->ParameterSet].EntropyCoding == 1)   // CABAC in use, need to get LastMB
  {
    len = GetVLCSymbol(buf, bitptr, &info, bufsize);
    linfo (len, info, &sh->CABAC_LastMB, &dummy);
    bitptr+=len;
  }

  bytes = bitptr/8;
  if (bitptr%8)
    bytes++;

  return bytes;

}



/*!
 *****************************************************************************
 *
 * \brief 
 *    Parses and interprets the UVLC-coded partition header (Type B and C packets only)
 *
 * \return
 *    negative in case of errors, the byte-index where the UVLC/CABAC MB data
 *    starts otherwise
 * Side effects:
 *    sh->PictureID und sh->SliceID set, all other values unchanged
 *
 * \date
 *    27 October, 2001
 *
 * \author
 *    Stephan Wenger   stewe@cs.tu-berlin.de
 *****************************************************************************/

int RTPInterpretPartitionHeader (byte *buf, int bufsize, RTPSliceHeader_t *sh)
{
  int len, info, bytes, dummy, bitptr=0;
  
  len = GetVLCSymbol(buf, bitptr, &info, bufsize);
  linfo (len, info, &sh->PictureID, &dummy);
  bitptr+=len;

  len = GetVLCSymbol(buf, bitptr, &info, bufsize);
  linfo (len, info, &sh->SliceID, &dummy);
  bitptr+=len;

  bytes = bitptr/8;
  if (bitptr%8)
    bytes++;

  return bytes;

}

/*!
 *****************************************************************************
 *
 * \brief 
 *    Reads and interprets the RTP sequence header, expects a type 6 packet
 *
 * \return
 *
 * Side effects:
 *   sets several fields in the img-> and inp-> structure, see RTPUseParameterSet
 *
 * \date
 *    27 October, 2001
 *
 * \author
 *    Stephan Wenger   stewe@cs.tu-berlin.de
 *****************************************************************************/

// Each RTP File is supposed to start with a type 6 (Header) packet.  It is necessary
// to read this early on in order to allocate the memory for the decoder.  This should
// be fixed some day in such a way that the decoder allocates memory as needed, and
// not statically at the first frame.

int RTPSequenceHeader (struct img_par *img, struct inp_par *inp, FILE *bits)
{
  int TotalPackLen;
  RTPpacket_t *p;
  RTPSliceHeader_t *sh;
  int err;
  int intime=0;

  assert (bits != NULL);

  p=alloca (sizeof (RTPpacket_t));
  sh=alloca(sizeof (RTPSliceHeader_t));

  if (4 != fread (&TotalPackLen,1, 4, bits))
    return -4711;    // EOF inidication
  if (4 != fread (&intime, 1, 4, bits))
    return -4712;

  p->packlen = TotalPackLen;
  p->packet = alloca (p->packlen);
  if (p->packlen != fread (p->packet, 1, p->packlen, bits))
    {
      // The corruption of a packet file is not a case we should handle.
      // In a real-world system, RTP packets may get lost, but they will
      // never get shortened.  Hence, the error checked here cannot occur.
      printf ("RTP File corruption, unexpected end of file, tried to read %d bytes\n", p->packlen);
      return -4713;    // EOF
    }

  p->paylen = p->packlen - 12;          // 12 bytes RTP header
  p->payload = alloca (p->paylen);   

  if (DecomposeRTPpacket (p) < 0)
    {
      // this should never happen, hence exit() is ok.  We probably do not want to attempt
      // to decode a packet that obviously wasn't generated by RTP
      printf ("Errors reported by DecomposePacket(), exit\n");
      exit (-700);
    }

    // Here the packet is ready for interpretation

  assert (p->pt == H26LPAYLOADTYPE);
  assert (p->ssrc == 0x12345678);

  if (p->payload[0] != 6)
  {
    printf ("RTPSequenceHeader: Expect Header Packet (FirstByet = 6), found packet type %d\n", p->payload[0]);
    exit (-1);
  }

  if ((err = RTPInterpretParameterSetPacket (&p->payload[1], p->paylen-1)) < 0)
    {
      printf ("RTPInterpretParameterSetPacket returns error %d\n", err);
    }

  RTPUseParameterSet (0, img, inp);
  img->number = 0;
  
  return 0;
}


/*!
 *****************************************************************************
 *
 * \brief 
 *    Sets various img->, inp-> and currSlice->struct members according to 
 *    the contents of the sh-> slice header structure
 *
 * \return
 *
 * Side effects:
 *    Set img->       qp, current_slice_nr, type, tr
 *    Set inp->
 *    Set currSlice-> qp, start_mb_nr, slice_nr, picture_type, picture_id (CABAC only: last_mb_nr)
 *
 * \date
 *    27 October, 2001
 *
 * \author
 *    Stephan Wenger   stewe@cs.tu-berlin.de
 *****************************************************************************/

void RTPSetImgInp (struct img_par *img, struct inp_par *inp, RTPSliceHeader_t *sh)
{
  static int ActualPictureType;
  Slice *currSlice = img->currentSlice;
  static int last_imgtr_frm=0,modulo_ctr_frm=0,last_imgtr_fld=0,modulo_ctr_fld=0;
  static int last_imgtr_frm_b=0,modulo_ctr_frm_b=0,last_imgtr_fld_b=0,modulo_ctr_fld_b=0;
 
  RMPNIbuffer_t *tmp_rmpni;
  MMCObuffer_t *tmp_mmco;

  img->qp = currSlice->qp = sh->InitialQP;

  if (sh->SliceType==2)
    img->qpsp = sh->InitialSPQP;

  currSlice->start_mb_nr = (img->width/16)*sh->FirstMBInSliceY+sh->FirstMBInSliceX;

  switch (sh->SliceType)
  {
  //! Potential BUG: do we need to distinguish between INTER_IMG_MULT and INTER_IMG?
  //!    similar with B_IMG_! and B_IMG_MULT
  //! also: need to define Slice types for SP images
  //! see VCEG-N72r1 for the Slice types, which are mapped here to img->type
  case 0:
    img->type = currSlice->picture_type = ParSet[CurrentParameterSet].UseMultpred?INTER_IMG_MULT:INTER_IMG_1;
    break;
  case 1:
    img->type = currSlice->picture_type = ParSet[CurrentParameterSet].UseMultpred?B_IMG_MULT:B_IMG_1;
    break;
  case 2:
    img->type = currSlice->picture_type = ParSet[CurrentParameterSet].UseMultpred?SP_IMG_MULT:SP_IMG_1;
    break;
  case 3:
    img->type = currSlice->picture_type = INTRA_IMG;
    break;
  default:
    printf ("Panic: unknown Slice type %d, conceal by loosing slice\n", sh->SliceType);
    currSlice->ei_flag = 1;
  }  
  

  //! The purpose of the following is to check for mixed Slices in one picture.
  //! According to VCEG-N72r1 and common sense this is allowed.  However, the
  //! current software seems to have a problem of some kind, to be checked.  Hence,
  //! printf a warning

  if (currSlice->start_mb_nr == 0)
    ActualPictureType = img->type;
  else
    if (ActualPictureType != img->type)
    {
      printf ("WARNING: mixed Slice types in a single picture -- interesting things may happen :-(\n");
    }

  img->tr = currSlice->picture_id = sh->PictureID;

#if 1
  img->structure = currSlice->structure = sh->structure; //picture structure: 
  
  if (img->type <= INTRA_IMG || img->type >= SP_IMG_1) 
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
  
  if(img->type != B_IMG_MULT) {
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
#endif
  
  currSlice->last_mb_nr = currSlice->start_mb_nr + sh->CABAC_LastMB;

  if (currSlice->last_mb_nr == currSlice->start_mb_nr)
    currSlice->last_mb_nr = img->max_mb_nr;

  /* KS: Multi Frame Buffering Syntax */
  img->pn=sh->PictureNum;

  // clear old slice RMPNI command buffer
  while (img->currentSlice->rmpni_buffer)
  {
    tmp_rmpni = img->currentSlice->rmpni_buffer;
  
    img->currentSlice->rmpni_buffer=tmp_rmpni->Next;
    free (tmp_rmpni);
  }

  img->currentSlice->rmpni_buffer=sh->RMPNIbuffer;

  sh->RMPNIbuffer=NULL;

  // free image MMCO buffer
  while (img->mmco_buffer)
  {
    tmp_mmco=img->mmco_buffer;

    img->mmco_buffer=tmp_mmco->Next;
    free (tmp_mmco);
  }

  // set image mmco bufer to actual MMCO buffer
  img->mmco_buffer=sh->MMCObuffer;
  sh->MMCObuffer=NULL;
}


/*!
 *****************************************************************************
 *
 * \brief 
 *    Function returns the next header type (SOS SOP, EOS)
 *    p-> and sh-> are filled.  p-> does not need memory for payload and packet
 *
 * \return
 *
 * Side effects:
 *
 * \date
 *    27 October, 2001
 *
 * \author
 *    Stephan Wenger   stewe@cs.tu-berlin.de
 *****************************************************************************/



int RTPGetFollowingSliceHeader (struct img_par *img, RTPpacket_t *p, RTPSliceHeader_t *sh)
{
  long int Filepos;
  int TotalPackLen;
  int done=0;
  int intime=0;
  Slice *currSlice = img->currentSlice;
  static unsigned int old_seq=0;
  static int first=0;
  static int i=2;
  
  RTPpacket_t *newp, *nextp;
  RTPSliceHeader_t *nextsh;
  
  assert (p != NULL);
  assert (sh != NULL);

  newp = alloca (sizeof (RTPpacket_t));
  newp->packet = alloca (MAXRTPPACKETSIZE);
  newp->payload = alloca(MAXRTPPAYLOADLEN);
  nextp=alloca (sizeof (RTPpacket_t));
  nextsh=alloca(sizeof (RTPSliceHeader_t));

  Filepos = ftell (bits);

  while (!done)
  {
    if (4 != fread (&TotalPackLen,1, 4, bits))
    {
      fseek (bits, Filepos, SEEK_SET);
      return EOS;    // EOF inidication
    }
    
    if (4 != fread (&intime, 1, 4, bits))
    {
      fseek (bits, Filepos, SEEK_SET);
      printf ("RTPGetFollowingSliceHeader: File corruption, could not read Timestamp\n");
      return EOS;
    }

    newp->packlen = TotalPackLen;
    assert (newp->packlen < MAXRTPPACKETSIZE);
    if (newp->packlen != fread (newp->packet, 1, newp->packlen, bits))
    {
      fseek (bits, Filepos, SEEK_SET);
      return EOS;    // EOF inidication
    }
    DecomposeRTPpacket (newp);
    if (newp->payload [0] == 0 || newp->payload[0] == 1)   // Full Slice or Partition A
      done = 1;
  }
  fseek (bits, Filepos, SEEK_SET);
  p->cc = newp->cc;
  p->m = newp->m;
  p->p = newp->p;
  p->pt = newp->pt;
  p->seq = newp->seq;
  p->ssrc = newp->ssrc;
  p->timestamp = newp->timestamp;
  p->v = newp->v;
  p->x = newp->x;

  if (p->seq != old_seq+i)
    currSlice->next_eiflag =1;
  else
  {
    currSlice->next_eiflag =0;
    old_seq=p->seq;
  }

  if(!first)
  {
    first=1;
    i=1;
  }

  RTPInterpretSliceHeader (&newp->payload[1], newp->packlen, newp->payload[0]==0?0:1, sh);
  if(currSlice->picture_id != sh->PictureID) 
    return (SOP);
  else
    return (SOS);
}
  
/*!
 ************************************************************************
 * \brief
 *    Input file containing Partition structures
 * \return
 *    TRUE if a startcode follows (meaning that all symbols are used).  \n
 *    FALSE otherwise.
 ************************************************************************/
int RTP_startcode_follows(struct img_par *img, struct inp_par *inp)
{
  Slice *currSlice = img->currentSlice;
  int dp_Nr = assignSE2partition[currSlice->dp_mode][SE_MBTYPE];
  DataPartition *dP = &(currSlice->partArr[dp_Nr]);
  Bitstream   *currStream = dP->bitstream;
  byte *buf = currStream->streamBuffer;
  int frame_bitoffset = currStream->frame_bitoffset;
  int info;

  if (currStream->ei_flag)
  {
    //printf ("ei_flag set, img->current_mb_nr %d, currSlice->last_mb_nr %d\n", img->current_mb_nr, currSlice->last_mb_nr);
    return (img->current_mb_nr == currSlice->last_mb_nr);
  }
  else
  {
    if (-1 == GetVLCSymbol (buf, frame_bitoffset, &info, currStream->bitstream_length))
      return TRUE;
    else
      return FALSE;
  }
}


/*!
 ************************************************************************
 * \brief
 *    read next UVLC codeword from SLICE-partition and
 *    map the corresponding syntax element
 *     Add Errorconceilment if necessary
 ************************************************************************
 */
int readSyntaxElement_RTP(SyntaxElement *sym, struct img_par *img, struct inp_par *inp, struct datapartition *dP)
{
  Bitstream   *currStream = dP->bitstream;

  if (RTP_symbols_available(currStream))    // check on existing elements in partition
  {
    RTP_get_symbol(sym, currStream);
  }
  else
  {
    set_ec_flag(sym->type);           // otherwise set error concealment flag
  }
  get_concealed_element(sym);

  sym->mapping(sym->len,sym->inf,&(sym->value1),&(sym->value2));

#if TRACE
  tracebits(sym->tracestring, sym->len, sym->inf, sym->value1, sym->value2);
#endif

  return 1;
}

/*!
 ************************************************************************
 * \brief
 *    checks if thererare symbols to read in the
 *    appropriate partition
 ************************************************************************
 */
int RTP_symbols_available (Bitstream *currStream)
{
  byte *buf = currStream->streamBuffer;
  int frame_bitoffset = currStream->frame_bitoffset;
  int info;

  if (currStream->ei_flag) {
    //printf ("RTP_symbols_available returns FALSE: ei_flag set\n");
    return FALSE;
  }
  if (-1 == GetVLCSymbol (buf, frame_bitoffset, &info, currStream->bitstream_length))
  {
    //printf ("RTP_symbols_available returns FALSE, no more symbols\nframe_bitoffset %d, bitstream_length %d read_len %d\n", currStream->frame_bitoffset, currStream->bitstream_length, currStream->read_len);
    return FALSE;
  }
  else
    return TRUE;
};

/*!
 ************************************************************************
 * \brief
 *    gets info and len of symbol
 ************************************************************************
 */
void RTP_get_symbol(SyntaxElement *sym, Bitstream *currStream)
{
  int frame_bitoffset = currStream->frame_bitoffset;
  byte *buf = currStream->streamBuffer;
  int BitstreamLengthInBytes = currStream->bitstream_length;

  sym->len =  GetVLCSymbol (buf, frame_bitoffset, &(sym->inf), BitstreamLengthInBytes);
  currStream->frame_bitoffset += sym->len;
}



/*!
 ************************************************************************
 * \brief
 *    Interrepts a parameter set packet
 ************************************************************************
 */

#define EXPECT_ATTR 0
#define EXPECT_PARSET_LIST 1
#define EXPECT_PARSET_NO 2
#define EXPECT_STRUCTNAME 3
#define EXPECT_STRUCTVAL_INT 4
#define EXPECT_STRUCTVAL_STRING 5

#define INTERPRET_COPY 100
#define INTERPRET_ENTROPY_CODING 101
#define INTERPRET_MOTION_RESOLUTION 102
#define INTERPRET_INTRA_PREDICTION 103
#define INTERPRET_PARTITIONING_TYPE 104


int RTPInterpretParameterSetPacket (char *buf, int buflen)

{
  // The dumbest possible parser that updates the parameter set variables 
  // static to this module

  int bufp = 0;
  int state = EXPECT_ATTR;
  int interpreter;
  int minus, number;
  int ps;
  char s[MAX_PARAMETER_STRINGLEN];
  void *destin;

  while (bufp < buflen)
  {
// printf ("%d\t%d\t %c%c%c%c%c%c%c%c  ", bufp, state, buf[bufp], buf[bufp+1], buf[bufp+2], buf[bufp+3], buf[bufp+4], buf[bufp+5], buf[bufp+6], buf[bufp+7], buf[bufp+8], buf[bufp+9]);

    switch (state)
    {
    case EXPECT_ATTR:
      if (buf[bufp] == '\004')      // Found UNIX EOF, this is the end marker
        return 0;
      if (strncmp ("a=H26L ", &buf[bufp], 7))
      {
        printf ("Parsing error EXPECT_ATTR in Header Packet: position %d, packet %s\n",
                 bufp, &buf[bufp]);
        return -1;
      }
      bufp += 7;
      state = EXPECT_PARSET_LIST;
      break;

    case EXPECT_PARSET_LIST:
      if (buf[bufp] != '(')
      {
        printf ("Parsing error EXPECT_PARSET in Header Packet: position %d, packet %s\n",
                bufp, &buf[bufp]);
        return -2;
      }
      bufp++;
      state = EXPECT_PARSET_NO;
      break;

    case EXPECT_PARSET_NO:
      number = 0;
      minus = 1;
      if (buf[bufp] == '-')
      {
        minus = -1;
        bufp++;
      }
      while (isdigit (buf[bufp]))
        number = number * 10 + ( (int)buf[bufp++] - (int) '0');
      if (buf[bufp] == ',')
      {
        printf ("Update of more than one prameter set not yet supported\n");
        return -1;
      }
      if (buf[bufp] != ')')
      {
        printf ("Parsing error NO PARSET LISTEND in Header Packet: position %d, packet %s\n",
          bufp, &buf[bufp]);
        return -1;
      }
      if (minus > 0)      // not negative
        ps = number;

      bufp+= 2;     // skip ) and blank
      state = EXPECT_STRUCTNAME;
      break;

    case EXPECT_STRUCTNAME:
      if (1 != sscanf (&buf[bufp], "%100s", s))   // Note the 100, which si MAX_PARAMETER_STRLEN
      {
        printf ("Parsing error EXPECT STRUCTNAME STRING in Header packet: position %d, packet %s\n",
               bufp, &buf[bufp]);
        return -1;
      }
      bufp += strlen (s);
      bufp++;       // Skip the blank


      if (!strncmp (s, "MaxPicID", MAX_PARAMETER_STRINGLEN))
      {
        state = EXPECT_STRUCTVAL_INT;
        interpreter = INTERPRET_COPY;
        destin = &ParSet[ps].MaxPicID;
        break;
      }
      if (!strncmp (s, "BufCycle", MAX_PARAMETER_STRINGLEN))
      {
        state = EXPECT_STRUCTVAL_INT;
        interpreter = INTERPRET_COPY;
        destin = &ParSet[ps].BufCycle;
        break;
      }
      if (!strncmp (s, "MaxPn", MAX_PARAMETER_STRINGLEN))
      {
        state = EXPECT_STRUCTVAL_INT;
        interpreter = INTERPRET_COPY;
        destin = &ParSet[ps].BufCycle;
        break;
      }

      if (!strncmp (s, "UseMultpred", MAX_PARAMETER_STRINGLEN))
      {
        state = EXPECT_STRUCTVAL_INT;
        interpreter = INTERPRET_COPY;
        destin = &ParSet[ps].UseMultpred;
        break;
      }

      if (!strncmp (s, "PixAspectRatioX", MAX_PARAMETER_STRINGLEN))
      {
        state = EXPECT_STRUCTVAL_INT;
        interpreter = INTERPRET_COPY;
        destin = &ParSet[ps].PixAspectRatioX;
        break;
      }
      if (!strncmp (s, "PixAspectRatioY", MAX_PARAMETER_STRINGLEN))
      {
        state = EXPECT_STRUCTVAL_INT;
        interpreter = INTERPRET_COPY;
        destin = &ParSet[ps].PixAspectRatioY;
        break;
      }
      if (!strncmp (s, "DisplayWindowOffsetTop", MAX_PARAMETER_STRINGLEN))
      {
        state = EXPECT_STRUCTVAL_INT;
        interpreter = INTERPRET_COPY;
        destin = &ParSet[ps].DisplayWindowOffsetTop;
        break;
      }
      if (!strncmp (s, "DisplayWindowOffsetBottom", MAX_PARAMETER_STRINGLEN))
      {
        state = EXPECT_STRUCTVAL_INT;
        interpreter = INTERPRET_COPY;
        destin = &ParSet[ps].DisplayWindowOffsetBottom;
        break;
      }
      if (!strncmp (s, "DisplayWindowOffsetRight", MAX_PARAMETER_STRINGLEN))
      {
        state = EXPECT_STRUCTVAL_INT;
        interpreter = INTERPRET_COPY;
        destin = &ParSet[ps].DisplayWindowOffsetRight;
        break;
      }
      if (!strncmp (s, "DisplayWindowOffsetLeft", MAX_PARAMETER_STRINGLEN))
      {
        state = EXPECT_STRUCTVAL_INT;
        interpreter = INTERPRET_COPY;
        destin = &ParSet[ps].DisplayWindowOffsetLeft;
        break;
      }
      if (!strncmp (s, "XSizeMB", MAX_PARAMETER_STRINGLEN))
      {
        state = EXPECT_STRUCTVAL_INT;
        interpreter = INTERPRET_COPY;
        destin = &ParSet[ps].XSizeMB;
        break;
      }
      if (!strncmp (s, "YSizeMB", MAX_PARAMETER_STRINGLEN))
      {
        state = EXPECT_STRUCTVAL_INT;
        interpreter = INTERPRET_COPY;
        destin = &ParSet[ps].YSizeMB;
        break;
      }
      if (!strncmp (s, "EntropyCoding", MAX_PARAMETER_STRINGLEN))
      {
        state = EXPECT_STRUCTVAL_STRING;
        interpreter = INTERPRET_ENTROPY_CODING;
        destin = &ParSet[ps].EntropyCoding;
        break;
      }
      if (!strncmp (s, "MotionResolution", MAX_PARAMETER_STRINGLEN))
      {
        state = EXPECT_STRUCTVAL_STRING;
        interpreter = INTERPRET_MOTION_RESOLUTION;
        destin = &ParSet[ps].MotionResolution;
        break;
      }
      if (!strncmp (s, "PartitioningType", MAX_PARAMETER_STRINGLEN))
      {
        state = EXPECT_STRUCTVAL_STRING;
        interpreter = INTERPRET_PARTITIONING_TYPE;
        destin = &ParSet[ps].PartitioningType;
        break;
      }
      if (!strncmp (s, "IntraPredictionType", MAX_PARAMETER_STRINGLEN))
      {
        state = EXPECT_STRUCTVAL_STRING;
        interpreter = INTERPRET_INTRA_PREDICTION;
        destin = &ParSet[ps].IntraPredictionType;
        break;
      }
      if (!strncmp (s, "HRCParameters", MAX_PARAMETER_STRINGLEN))
      {
        state = EXPECT_STRUCTVAL_INT;
        interpreter = INTERPRET_COPY;
        destin = &ParSet[ps].HRCParameters;
        break;
      }
      if (!strncmp (s, "FramesToBeEncoded", MAX_PARAMETER_STRINGLEN))
      {
        state = EXPECT_STRUCTVAL_INT;
        interpreter = INTERPRET_COPY;
        destin = &InfoSet.FramesToBeEncoded;
        break;
      }
      if (!strncmp (s, "FrameSkip", MAX_PARAMETER_STRINGLEN))
      {
        state = EXPECT_STRUCTVAL_INT;
        interpreter = INTERPRET_COPY;
        destin = &InfoSet.FrameSkip;
        break;
      }
      if (!strncmp (s, "SequenceFileName", MAX_PARAMETER_STRINGLEN))
      {
        state = EXPECT_STRUCTVAL_STRING;
        interpreter = INTERPRET_COPY;
        destin = &InfoSet.SequenceFileName;
        break;
      }
      if (!strncmp (s, "NumberBFrames", MAX_PARAMETER_STRINGLEN))
      {
        state = EXPECT_STRUCTVAL_INT;
        interpreter = INTERPRET_COPY;
        destin = &InfoSet.NumberBFrames;
        break;
      }
     
      // Here, all defined Parameter names are checked.  Anything else is a syntax error
      printf ("Syntax Error: unknown Parameter %s\n", s);
      printf ("Parsing error in Header Packet: position %d, packet %s\n",
          bufp, buf);
      return -3;
      break;        // to make lint happy
    
    case EXPECT_STRUCTVAL_INT:
      if (1!=sscanf (&buf[bufp], "%d", (int *)destin))
      {
        printf ("Parsing error EXPECT STRUCTVAL INT in Header Packet: position %d, packet %s\n",
          bufp, &buf[bufp]);
        return -4;
      }
// printf ("EXPECT_STRCUTVAL_INT: write %d\n", * (int *)destin);
      while (bufp < buflen && buf[bufp] != '\n')    // Skip any trailing whitespace and \n
        bufp++;
      bufp++;
      state=EXPECT_ATTR;
      break;
      
    case EXPECT_STRUCTVAL_STRING:
      if (1 != sscanf (&buf[bufp], "%100s", s))
      {
        printf ("Parsing error EXPECT STRUCTVAL STRING in Header Packet: position %d, packet %s\n",
          bufp, &buf[bufp]);
        return -5;
      }
      while (bufp < buflen && buf[bufp] != '\n')   // Skip any trailing whitespace and \n
        bufp++;
      bufp++;
      state=EXPECT_ATTR;

      switch (interpreter)
      {
      case INTERPRET_COPY:
        // nothing -- handled where it occurs
        break;
      case INTERPRET_ENTROPY_CODING:
        if (!strncmp (s, "UVLC", 4))
          * (int *)destin = 0;
        else
          * (int *)destin = 1;
//        printf ("in INterpret, Entropy COding :%s: results in %d\n", s, *(int *)destin);
        break;
      case INTERPRET_MOTION_RESOLUTION:
        if (!strncmp (s, "quater", 6))
          * (int *)destin = 0;
        else
          * (int *)destin = 1;
        break;
      case INTERPRET_INTRA_PREDICTION:
        if (!strncmp (s, "InterPredicted", 14))
          * (int *)destin = 0;
        else
          * (int *)destin = 1;
        break;
      case INTERPRET_PARTITIONING_TYPE:
        if (!strncmp (s, "one", 3))
          * (int *)destin = 0;
        else
          * (int *)destin = 1;
        break;
      default:
        assert (0==1);
      }
    break;

    default:
      printf ("Parsing error UNDEFINED SYNTAX in Header Packet: position %d, packet %s\n",
        bufp, &buf[bufp]);
      return -1;
    }
//  printf ("\t\t%d\n", bufp);
  }

//  printf ("CurrentParameterSet %d, ps %d\n", CurrentParameterSet, ps);
//  printf ("RTPInterpretParameterPacket: xsize x Ysize, %d x %d, Entropy %d, Motion %d  MaxPicId %d\n",
//    ParSet[ps].XSizeMB, ParSet[ps].YSizeMB, ParSet[ps].EntropyCoding, ParSet[ps].MotionResolution, ParSet[ps].MaxPicID);
  ParSet[ps].Valid = 1;
  return 0;
}


/*!
 ************************************************************************
 * \brief
 *    Update img->xxx with the content of a parameter set, called for 
 *    every slice 
 ************************************************************************
 */



void RTPUseParameterSet (int n, struct img_par *img, struct inp_par *inp)
{
  int status;
  
  if (n == CurrentParameterSet)
    return;   // no change

//  printf ("Use a new parameter set: old %d, new %d\n", CurrentParameterSet, n);
  CurrentParameterSet = n;

  status = RTP_PARAMETER_SET_OK;

  if (CurrentParameterSet < 0 || (CurrentParameterSet > RTP_MAX_PARAMETER_SET))
  {
    printf ("Parameter Set %d out of range, conceal to 0\n", CurrentParameterSet);
    CurrentParameterSet = 0;    // and pray that it works...
    status = RTP_PARAMETER_SET_OUT_OF_RANGE;
  }

  if (!ParSet[CurrentParameterSet].Valid)
  {
    printf ("Try to use uninitialized Parameter Set %d, conceal to 0\n", CurrentParameterSet);
    CurrentParameterSet = 0;
    status = RTP_PARAMETER_SET_INVALID;
  }

  // Now updates global decoder variables with the appropriate parameter set.

  // A full-fledged decoder would make some consistency checks.  For example,
  // it makes sense to change the MotionResolution or the EntropyCode within
  // a picture (img->current_mb_nr != 0).  It doesn't make sense to change
  // the pixel aspect ratio.
  // There is no need to do those checks in an error free environment -- an
  // encoder is simply not supposed to do so.  In error prone environments,
  // however, this is an additional means for error detection.

  // Note: Many parameters are available in both the input-> and img-> structures.
  // Some people seem to use the input-> parameters, others the img-> parameters.
  // This should be cleaned up one day.
  // For now simply copy any updated variables into both structures.

  // MaxPicID: doesn't exist in pinput-> or img->
  inp->buf_cycle = ParSet[CurrentParameterSet].BufCycle;
  img->buf_cycle = inp->buf_cycle+1;      // see init_global_buffers()

  // PixAspectRatioX: doesn't exist
  // PixAspectRatioY: doesn't exist
  // DisplayWindowOffset*: doesn't exist

  // XSizeMB
  img->width = ParSet[CurrentParameterSet].XSizeMB*16;
  img->width_cr = ParSet[CurrentParameterSet].XSizeMB*8;

  //YSizeMB
  img->height = ParSet[CurrentParameterSet].YSizeMB*16;
  img->height_cr = ParSet[CurrentParameterSet].YSizeMB*8;

  // EntropyCoding
  if (ParSet[CurrentParameterSet].EntropyCoding == 0)
    inp->symbol_mode = UVLC;
  else
    inp->symbol_mode = CABAC;

  // MotionResolution
  img->mv_res = ParSet[CurrentParameterSet].MotionResolution;

  // PartitioningType
  inp->partition_mode = ParSet[CurrentParameterSet].PartitioningType;

  // IntraPredictionType
  inp->UseConstrainedIntraPred = ParSet[CurrentParameterSet].IntraPredictionType;
  
  //! img->type: This is calculated by using ParSet[CurrentParameterSet].UseMultpred
  //! and the slice type from the slice header.  It is set in RTPSetImgInp()

  // HRCParameters: Doesn't exist
}


/*!
 *****************************************************************************
 *
 * \brief 
 *    RTPReadPacket reads one packet from file
 *
 * \return
 *    0 in case of success, -1 in case of error
 *
 * \para Paremeters
 *    p: packet data structure, with memory for p->packet allocated
 *
 * Side effects:
 *   - File pointer in bits moved
 *   - p->xxx filled by reading and Decomposepacket()
 *
 * \date
 *    04 November, 2001
 *
 * \author
 *    Stephan Wenger, stewe@cs.tu-berlin.de
 *****************************************************************************/

int RTPReadPacket (RTPpacket_t *p, FILE *bits)
{
  int Filepos, intime;

  assert (p != NULL);
  assert (p->packet != NULL);
  assert (p->payload != NULL);

  Filepos = ftell (bits);
  if (4 != fread (&p->packlen,1, 4, bits))
    {
      printf ("Unable to read 4 bytes for the RTP packet size\n");
      fseek (bits, Filepos, SEEK_SET);
      return -1;
    }
    
  if (4 != fread (&intime, 1, 4, bits))
    {
      fseek (bits, Filepos, SEEK_SET);
      printf ("RTPReadPacket: File corruption, could not read Timestamp, exit\n");
      exit (-1);
    }

  assert (p->packlen < MAXRTPPACKETSIZE);
  if (p->packlen != fread (p->packet, 1, p->packlen, bits))
    {
      printf ("RTPReadPacket: File corruption, could not read %d bytes\n", p->packlen);
      exit (-1);    // EOF inidication
    }
  if (DecomposeRTPpacket (p) < 0)
    {
      // this should never happen, hence exit() is ok.  We probably do not want to attempt
      // to decode a packet that obviously wasn't generated by RTP
      printf ("Errors reported by DecomposePacket(), exit\n");
      exit (-700);
    }
    assert (p->pt == H26LPAYLOADTYPE);
    assert (p->ssrc == 0x12345678);
  return 0;
}


/*!
 *****************************************************************************
 *
 * \brief 
 *    RTPReadDataPartitionedSlice collects all partitiobnss of the slice
 *
 *
 * \return
 *  sequence number of the last packed that has been concealed
 *
 * VCEG-N72r1 may not state it explicitely, but we decided that partitions
 * have to be sent immediately after each other.  Meaning, if A, B, and C are
 * all present, the rtp-timestamp of A is 2 bigger than the one of A.
 * Hence, this function reads all those partitions into their respective 
 * buffers
 *
 * Side effects: many!
 *   - File pointer in bits moved
 *   - currSlice->partArr[0, 1, 2] updated
 *   - img-> and inp-> updated, see RTPSetImgInp() 
 *
 * \date
 *    04 November, 2001
 *
 * \author
 *    Stephan Wenger, stewe@cs.tu-berlin.de
 *****************************************************************************/

#define SEQ_PLUS_1  (a->seq+1 == b->seq)
#define SEQ_PLUS_2  (a->seq+2 == b->seq)
#define SAME_TIME (a->timestamp == b->timestamp)
#define SLICE_NO_OK (b_SliceID == a_SliceID)
#define SAME_SLICE (SAME_TIME && SLICE_NO_OK)
#define TYPE_B (b->payload[0] == 2)
#define TYPE_C (b->payload[0] == 3)

void RTPProcessDataPartitionedSlice (struct img_par *img, struct inp_par *inp, FILE *bits, 
                                     RTPpacket_t *a, int a_SliceID)


{

//!  BUG: need to fix wrap-around-problem for the sequence number
//!       needs to be fixed only for HUGE files (more than 2^^16 packets)
//!       Note: in contrast to RTP spec, encoder starts sequence no at 0 to ease debugging

  RTPpacket_t *b, *c;
  RTPSliceHeader_t *sh;
  long Filepos;
  int StartMBData;
  int b_SliceID, b_PicId, c_SliceID, c_PicID;
  Slice *currSlice = img->currentSlice;

  Filepos = ftell (bits);

  b=alloca(sizeof(RTPpacket_t));
  b->packet=alloca (MAXRTPPACKETSIZE);
  b->payload=alloca (MAXRTPPAYLOADLEN);
  sh=alloca (sizeof(RTPSliceHeader_t));

  free_Partition (currSlice->partArr[1].bitstream);
  free_Partition (currSlice->partArr[2].bitstream);


  if (RTPReadPacket (b, bits) != 0)
  {
    printf ("Error while reading RTP packet, assuming intentional EOF and the last slice contains no type B, C\n");
    img->currentSlice->partArr[1].bitstream->bitstream_length=0;
    img->currentSlice->partArr[1].bitstream->code_len=0;
    img->currentSlice->partArr[1].bitstream->ei_flag=0;
    img->currentSlice->partArr[2].bitstream->bitstream_length=0;
    img->currentSlice->partArr[2].bitstream->code_len=0;
    img->currentSlice->partArr[2].bitstream->ei_flag=0;
    return;
  }

  StartMBData = RTPInterpretPartitionHeader (&b->payload[1], b->paylen, sh);
  //StartMBData++;   // Skip First Byte
  b_SliceID = sh->SliceID;
  b_PicId = sh->PictureID;

  // General: 
  //   Lost partition, identified by missing packet due to sequence number problems, or
  //   Empty partition, identified by not present partition
  // The two cases are signalled to higher layers as follows
  //   Lost Partition: currStream->ei_flag is set, length inidicators don't care
  //   Empty Partition: only length indicator is zero, ei_flag is cleared


  if (b->seq > a->seq+2)    // two or more consecutive packet losses, none, one, or both 
  {                         // partitions lost and no way to figure which ones
    fseek (bits, Filepos, SEEK_SET);    // go back 
    printf ("lost at least two partitions\n");
    img->currentSlice->partArr[1].bitstream->bitstream_length=0;
    img->currentSlice->partArr[1].bitstream->code_len=0;
    img->currentSlice->partArr[1].bitstream->ei_flag=1;
    img->currentSlice->partArr[1].readSyntaxElement = readSyntaxElement_RTP;
    img->currentSlice->partArr[2].bitstream->bitstream_length=0;
    img->currentSlice->partArr[2].bitstream->code_len=0;
    img->currentSlice->partArr[2].bitstream->ei_flag=1;
    img->currentSlice->partArr[2].readSyntaxElement = readSyntaxElement_RTP;
    return;
  }

  if (SEQ_PLUS_2 && !SAME_SLICE) // one packet was lost could be type B partition or type C partition
  {
    fseek (bits, Filepos, SEEK_SET);    // go back 
    printf ("lost one Partition of Type B or Type C\n");
    img->currentSlice->partArr[1].bitstream->bitstream_length=0;
    img->currentSlice->partArr[1].bitstream->code_len=0;
    img->currentSlice->partArr[1].bitstream->ei_flag=1;
    img->currentSlice->partArr[1].readSyntaxElement = readSyntaxElement_RTP;
    img->currentSlice->partArr[2].bitstream->bitstream_length=0;
    img->currentSlice->partArr[2].bitstream->code_len=0;
    img->currentSlice->partArr[2].bitstream->ei_flag=1;
    img->currentSlice->partArr[2].readSyntaxElement = readSyntaxElement_RTP;
    return;
  }

  if (SEQ_PLUS_2 && SAME_SLICE && TYPE_B)
  {
    printf ("Panic: this should never happen: SEQ_PLUS_2 && SAME_SLICE && TYPE_B\n");
    exit (-1);
  }

  if (SEQ_PLUS_2 && SAME_SLICE && TYPE_C) //packet containing a type B partition was lost
  {
    c = b;
    printf ("lost one partition of Type B\n");
    img->currentSlice->partArr[1].bitstream->ei_flag=1;
    img->currentSlice->partArr[1].bitstream->bitstream_length=0;
    img->currentSlice->partArr[1].bitstream->code_len=0;
    img->currentSlice->partArr[1].readSyntaxElement = readSyntaxElement_RTP;
    img->currentSlice->partArr[2].bitstream->ei_flag=0;
    CopyPartitionBitstring (img, c, img->currentSlice->partArr[2].bitstream, 2);    // copy to C
    return;
  }

  if (SEQ_PLUS_1 && !SAME_SLICE)    // Type B and Type C intentionally not coded
  {
    fseek (bits, Filepos, SEEK_SET);    // go back 
    img->currentSlice->partArr[1].bitstream->bitstream_length=0;
    img->currentSlice->partArr[1].bitstream->code_len=0;
    img->currentSlice->partArr[1].bitstream->ei_flag=0;
    img->currentSlice->partArr[2].bitstream->bitstream_length=0;
    img->currentSlice->partArr[2].bitstream->code_len=0;
    img->currentSlice->partArr[2].bitstream->ei_flag=0;
    return;
  }

  if (SEQ_PLUS_1 && SAME_SLICE && TYPE_C) // Type B intentionally not coded
  {
    c = b;
    img->currentSlice->partArr[1].bitstream->bitstream_length=0;
    img->currentSlice->partArr[1].bitstream->code_len=0;
    img->currentSlice->partArr[1].bitstream->ei_flag=0;
    img->currentSlice->partArr[2].bitstream->ei_flag=0;
    img->currentSlice->ei_flag=0;
    CopyPartitionBitstring (img, c, img->currentSlice->partArr[2].bitstream, 2);    // copy to C
    return;
  }

  if (SEQ_PLUS_1 && SAME_SLICE && TYPE_B)
  {
    // the normal case, found type b partition in the right sequence 
    // copy type B and then check for type C
    
    img->currentSlice->partArr[1].bitstream->ei_flag=0;
    CopyPartitionBitstring (img, b, img->currentSlice->partArr[1].bitstream, 1);    // copy to B
    c = b;    // re-use buffer, could as well use b variable
    Filepos = ftell (bits);

    if (RTPReadPacket (c, bits) != 0)
    {
      printf ("Error while reading RTP packet for Partition C, assuming intentional EOF and the last slice contains no type B, C -- B-frame slice?\n");
      img->currentSlice->partArr[2].bitstream->bitstream_length=0;
      img->currentSlice->partArr[2].bitstream->code_len=0;
      img->currentSlice->partArr[2].bitstream->ei_flag=0;
      img->currentSlice->ei_flag=0;
      img->currentSlice->next_header=EOS;
      img->currentSlice->eos_flag=1;
      return;
    }

    StartMBData = RTPInterpretPartitionHeader (&b->payload[1], b->paylen, sh);
    c_SliceID = sh->SliceID;
    c_PicID = sh->PictureID;

    if (a->seq+2 != c->seq)     // packet loss
    {
      fseek (bits, Filepos, SEEK_SET);    // go back 
      printf ("lost one partition of Type C\n");
      img->currentSlice->partArr[2].bitstream->bitstream_length=0;
      img->currentSlice->partArr[2].bitstream->code_len=0;
      img->currentSlice->partArr[2].bitstream->ei_flag=1;
      img->currentSlice->partArr[2].readSyntaxElement = readSyntaxElement_RTP;
      return;
    }
    
    // If we reached this point the sequence number is ok
    // Now check if this is a type C partition or the next type A partition
    if (!TYPE_C                         ||    // not type C packet      or
        c->timestamp != a->timestamp    ||    // incorrect timestamp    or
        c_SliceID != a_SliceID)               // wrong SliceID
    { // Partition C intentrionally not coded, 
      fseek (bits, Filepos, SEEK_SET);    // go back 
      return;
    }
    else    // correct type C
    {
      img->currentSlice->partArr[2].bitstream->ei_flag=0;
      img->currentSlice->ei_flag=0;
      CopyPartitionBitstring (img, c, img->currentSlice->partArr[2].bitstream, 2);    // copy to C
      return;
    }
    assert (1==2);
  }
  printf ("This should never happen\n");
  assert (0==1);
}



void CopyPartitionBitstring (struct img_par *img, RTPpacket_t *p, Bitstream *b, int dP)
{
  int header_bytes;
  RTPSliceHeader_t *sh = alloca (sizeof(RTPSliceHeader_t));

  header_bytes = RTPInterpretPartitionHeader (&p->payload[1], p->paylen-1, sh);
  header_bytes++;     // for the First Byte

  b->bitstream_length = b->code_len = (p->paylen - header_bytes);
  
  //! TO 15.01.2002 for Debug only
  if(b->bitstream_length < 1)
    printf ("Empty Partition\n"); //assert (0==1);
  //! End TO
  
  b->frame_bitoffset = b->read_len = 0;

  memcpy (b->streamBuffer, &p->payload[header_bytes], b->bitstream_length);

  if (ParSet[CurrentParameterSet].EntropyCoding == CABAC) 
  {
    byte *buf;
    DecodingEnvironment *dep;

    buf = b->streamBuffer;
    dep = &((img->currentSlice->partArr[dP]).de_cabac);
    arideco_start_decoding(dep, buf, 0, &b->read_len);
  }
}
