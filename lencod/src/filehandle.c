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
 *    filehandle.c
 * \brief
 *    Handles the operations how to write
 *    the generated symbols on the interim file format or some
 *    other specified output formats
 * \author
 *    Main contributors (see contributors.h for copyright, address and affiliation details)
 *      - Thomas Stockhammer            <stockhammer@ei.tum.de>
 *      - Detlev Marpe                  <marpe@hhi.de>
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

// Global Varibales

static FILE *out;   //!< output file

/*
   The implemented solution for a unified picture header
*/

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
  exit(code);
}

/*!
 ************************************************************************
 * \brief
 *    Enforces Byte Alihnment of the Bitstream b
 * \param b
 *    Bitstream to be byte aligned
 ************************************************************************
 */

static void ByteAlign (Bitstream *b)
{
  if (b->bits_to_go != 8)
  {
    b->byte_buf <<= b->bits_to_go;
    b->bits_to_go=8;
    b->streamBuffer[b->byte_pos++]=b->byte_buf;
    b->byte_buf=0;
  }
}

/*!
 ************************************************************************
 * \brief
 *    This function opens the output files and generates the
 *    appropriate sequence header
 ************************************************************************
 */
int start_sequence()
{
  int len;

  if ((out=fopen(input->outfile,"wb"))==NULL)
  {
    snprintf(errortext, ET_SIZE, "Error open file %s  \n",input->outfile);
    error(errortext,1);
  }

  switch(input->of_mode)
  {
    case PAR_OF_26L:
      len = SequenceHeader(out);
      return 0;
    case PAR_OF_RTP:
      len = RTPSequenceHeader(out);
      return 0;
    case PAR_OF_IFF:
      // Tian Dong: June 10, 2002
      len = initInterimFile();
//      len = SequenceHeader(out);
      return len;
      
    default:
      snprintf(errortext, ET_SIZE, "Output File Mode %d not supported", input->of_mode);
      error(errortext,1);
      return 1;
  }
}

/*!
 ************************************************************************
 * \brief
 *     This function terminates the sequence and closes the
 *     output files
 ************************************************************************
 */
int terminate_sequence()
{
  Bitstream *currStream;

  // Mainly flushing of everything
  // Add termination symbol, etc.

  switch(input->of_mode)
  {
    case PAR_OF_IFF:
      terminateInterimFile(out);
      fclose (out);
      return 0;
    case PAR_OF_26L:
      currStream = ((img->currentSlice)->partArr[0]).bitstream;
      if (input->symbol_mode == UVLC)
      {

        // Current TML File Format
        // Get the trailing bits of the last slice
        currStream->bits_to_go  = currStream->stored_bits_to_go;
        currStream->byte_pos    = currStream->stored_byte_pos;
        currStream->byte_buf    = currStream->stored_byte_buf;

        if (currStream->bits_to_go < 8)   // there are bits left in the last byte
          currStream->streamBuffer[currStream->byte_pos++] = currStream->byte_buf;
        // Write all remaining bits to output bitstream file
        fwrite (currStream->streamBuffer, 1, currStream->byte_pos, out);
        fclose(out);
      }
      else
      {
        // CABAC File Format
        fclose(out);
      }
      return 0;
    case PAR_OF_RTP:
      fclose (out);
      return 0;
    default:
      snprintf(errortext, ET_SIZE, "Output File Mode %d not supported", input->of_mode);
      error(errortext,1);
      return 1;
  }
}

/*!
 ************************************************************************
 *  \brief
 *     This function generates the appropriate slice
 *     header
 *
 *  \return number of bits used for the picture header, including the PSC.
 *
 *  \par Side effects:
 *      Adds picture header symbols to the symbol buffer
 *
 *  \par Remarks:
 *      THIS IS AN INTERIM SOLUTION FOR A PICTURE HEADER, see VCEG-M79
 *                                                                               \par
 *      Generates the Picture Header out of the information in img_par and inp-par,
 *      and writes it to the symbol buffer.  The structure of the Picture Header
 *      is discussed in VCEG-M79.  It is implemented as discussed there.  In addition,
 *      it is preceeded by a picture start code, a UVLC-codeword of LEN=31 and INFO = 0.
 *      It would have been possible to put information into this codeword, similar
 *      to designs earlier than TML 5.9.  But it is deemed that the waste of 15
 *      data bits is acceptable considering the advantages of being able to search
 *      for a picture header quickly and easily, by looking for 30 consecutive, byte-
 *      aligned zero bits.
 *                                                                               \par
 *      The accounting of the header length (variable len) relies on the side effect
 *      of writeUVLCSymbol() that sets len and info in the symbol variable parameter
 ************************************************************************
*/

int start_slice()
{
  EncodingEnvironmentPtr eep;
  Slice *currSlice = img->currentSlice;
  Bitstream *currStream;
  int header_len;
  int i;

  currSlice->num_mb = 0;  // no coded MBs so far

  switch(input->of_mode)
  {
    case PAR_OF_IFF:    // for Interim File Format
      header_len = 0;
      if (input->symbol_mode == UVLC)
      {
        currStream = (currSlice->partArr[0]).bitstream;

        assert (currStream->bits_to_go == 8);
        assert (currStream->byte_buf == 0);
        assert (currStream->byte_pos == 0);

        pCurrPayloadInfo = newPayloadInfo();
        addOnePayloadInfo( pCurrPayloadInfo );
      }
      else
      {                   // CABAC 
        eep = &((currSlice->partArr[0]).ee_cabac);
        currStream = (currSlice->partArr[0]).bitstream;

        assert (currStream->bits_to_go == 8);
        assert (currStream->byte_buf == 0);
        assert (currStream->byte_pos == 0);
        memset(currStream->streamBuffer, 0, 12);

        pCurrPayloadInfo = newPayloadInfo();
        addOnePayloadInfo( pCurrPayloadInfo );
        

        arienco_start_encoding(eep, currStream->streamBuffer, &(currStream->byte_pos),&(currStream->last_startcode),img->type);
        // initialize context models
        init_contexts_MotionInfo (currSlice->mot_ctx);
        init_contexts_TextureInfo(currSlice->tex_ctx);
      }
      return header_len;
      break;

    case PAR_OF_26L:
      if (input->symbol_mode == UVLC)
      {
        currStream = (currSlice->partArr[0]).bitstream;
        currStream->last_startcode=currStream->byte_pos;
        header_len = SliceHeader (0);  // Slice Header without Start Code
       
        return header_len;
      }
      else
      {                   // H.26: CABAC File Format
        eep = &((currSlice->partArr[0]).ee_cabac);
        currStream = (currSlice->partArr[0]).bitstream;

        if((img->current_slice_nr==0)&&(img->pstruct!=2))
        {
          assert (currStream->bits_to_go == 8);
          assert (currStream->byte_buf == 0);
          assert (currStream->byte_pos == 0);
          memset(currStream->streamBuffer, 0, 12);    // fill first 12 bytes with zeros (debug only)
        }

        currStream->last_startcode=currStream->byte_pos;
        header_len = SliceHeader (0);  // Slice Header without Start Code
        
        // Note that SliceHeader() sets the buffer pointers as a side effect
        // Hence no need for adjusting it manually (and preserving space to be patched later

        if (currStream->bits_to_go != 8)
          header_len+=currStream->bits_to_go;
        ByteAlign (currStream);

        // If there is an absolute need to communicate the partition size, this would be
        // the space to insert it
        arienco_start_encoding(eep, currStream->streamBuffer, &(currStream->byte_pos),&(currStream->last_startcode),img->type);
        // initialize context models


        init_contexts_MotionInfo (currSlice->mot_ctx);
        init_contexts_TextureInfo(currSlice->tex_ctx);

        return header_len;

      }
    case PAR_OF_RTP:
      if (img->current_mb_nr == 0)                       // new picture
        RTPUpdateTimestamp (img->tr);

      if (input->symbol_mode == UVLC)
      {
        currStream = currSlice->partArr[0].bitstream;

        assert (currStream != NULL);
        assert (currStream->bits_to_go == 8);
        if(mref!=mref_fld || !img->fld_type) assert (currStream->byte_pos == 0);
        if(mref!=mref_fld || !img->fld_type) assert (currStream->byte_buf == 0);

        header_len=RTPSliceHeader();                      // generate the slice header
        
        if (currStream->bits_to_go != 8)
          header_len+=currStream->bits_to_go;
        ByteAlign (currStream);

        assert (currStream->bits_to_go == 8);

        if (input->partition_mode != PAR_DP_1)
        {
          for (i=1; i<currSlice->max_part_nr; i++)
          {
            // generate the Partition header for all partitions.  terminate_slice()
            // later writes only those partitions which contain more data than
            // header_len (and worries about byte alignment itself!)

            currStream = (currSlice->partArr[i]).bitstream;
            assert (currStream != NULL);
            assert (currStream->bits_to_go == 8);
            assert (currStream->byte_pos == 0);
            assert (currStream->byte_buf == 0);
            RTPPartition_BC_Header(i);
            ByteAlign (currStream);
            currStream->write_flag = 0;
          }
        }
        return header_len;
      }      

      else
      {                   // RTP: CABAC 
        eep = &((currSlice->partArr[0]).ee_cabac);
        currStream = (currSlice->partArr[0]).bitstream;

        assert (currStream->bits_to_go == 8);
        assert (currStream->byte_buf == 0);

        header_len=RTPSliceHeader();                      // generate the slice header

        if (currStream->bits_to_go != 8)
          header_len+=currStream->bits_to_go;
        ByteAlign (currStream);

        arienco_start_encoding(eep, currStream->streamBuffer, &(currStream->byte_pos),&(currStream->last_startcode),img->type);
        currStream->write_flag = 0;

        if(input->partition_mode != PAR_DP_1)
        {
          for (i=1; i<currSlice->max_part_nr; i++)
          {
            eep = &((currSlice->partArr[i]).ee_cabac);
            currStream = (currSlice->partArr[i]).bitstream;

            RTPPartition_BC_Header(i);
            ByteAlign (currStream);

            assert (currStream->bits_to_go == 8);
            assert (currStream->byte_buf == 0);

            arienco_start_encoding(eep, currStream->streamBuffer, &(currStream->byte_pos),&(currStream->last_startcode),img->type);
            currStream->write_flag = 0;
          }
        }
        // initialize context models
        init_contexts_MotionInfo (currSlice->mot_ctx);
        init_contexts_TextureInfo(currSlice->tex_ctx);

        return header_len;
      }
    default:
      snprintf(errortext, ET_SIZE, "Output File Mode %d not supported", input->of_mode);
      error(errortext,1);
      return 1;
  }
}


/*!
 ************************************************************************
 * \brief
 *    This function terminates a slice
 * \return
 *    0 if OK,                                                         \n
 *    1 in case of error
 ************************************************************************
 */
int terminate_slice(int write_out)
{
  int bytes_written;
  Bitstream *currStream;
  Slice *currSlice = img->currentSlice;
  EncodingEnvironmentPtr eep;
  int i;
  int LastPartition;
  int rtp_bytes_written;
  int temp_byte_pos;

  // Mainly flushing of everything
  // Add termination symbol, etc.
  switch(input->of_mode)
  {
     case PAR_OF_IFF:
      assert( box_atm.fpMedia != NULL );

      if (input->symbol_mode == UVLC)
      {
        // Enforce byte alignment of next header: zero bit stuffing
        currStream = (currSlice->partArr[0]).bitstream;

        if (currStream->bits_to_go < 8)
        { // trailing bits to process
          currStream->byte_buf <<= currStream->bits_to_go;
          stat->bit_use_stuffingBits[img->type]+=currStream->bits_to_go;
          currStream->streamBuffer[currStream->byte_pos++]=currStream->byte_buf;
          currStream->bits_to_go = 8;
        }

        if(input->Encapsulated_NAL_Payload)
        {
          SODBtoRBSP(currStream);
        }
        bytes_written = currStream->byte_pos;
        if(input->Encapsulated_NAL_Payload)
        {
          temp_byte_pos = bytes_written;
          bytes_written = RBSPtoEBSP(currStream->streamBuffer, 0, bytes_written);
          *(stat->em_prev_bits) += (bytes_written - temp_byte_pos) * 8;
        }
        stat->bit_ctr += 8*bytes_written;     // actually written bits
        fwrite (currStream->streamBuffer, 1, bytes_written, box_atm.fpMedia );

        currStream->stored_bits_to_go = 8; // store bits_to_go
        currStream->stored_byte_buf   = currStream->byte_buf;   // store current byte
        currStream->stored_byte_pos   = 0; // reset byte position
      }
      else
      {
        // CABAC File Format
        write_terminating_bit (1);

        eep = &((currSlice->partArr[0]).ee_cabac);
        currStream = (currSlice->partArr[0]).bitstream;
        // terminate the arithmetic code
        stat->bit_use_stuffingBits[img->type]+=get_trailing_bits(eep);
        arienco_done_encoding(eep);
        if (eep->Ebits_to_go != 8)
          stat->bit_use_stuffingBits[img->type]+=eep->Ebits_to_go;

        if(input->Encapsulated_NAL_Payload)
        {
          SODBtoRBSP(currStream);
        }
        bytes_written = currStream->byte_pos; // number of written bytes
        if(input->Encapsulated_NAL_Payload)
        {
          temp_byte_pos = bytes_written;
          bytes_written = RBSPtoEBSP(currStream->streamBuffer, 0, bytes_written);
          *(stat->em_prev_bits) += (bytes_written - temp_byte_pos) * 8;
        }
        stat->bit_ctr += 8*bytes_written;     // actually written bits
        fwrite (currStream->streamBuffer, 1, bytes_written, box_atm.fpMedia );

        // Provide the next partition with a 'fresh' buffer
        currStream->stored_bits_to_go = 8;
        currStream->stored_byte_buf   = 0;
        currStream->stored_byte_pos   = 0;
        currStream->bits_to_go = 8;
        currStream->byte_buf   = 0;
        currStream->byte_pos   = 0;
      }

      currPictureInfo.currPictureSize += bytes_written;  // caculate the size of currentPicture
      pCurrPayloadInfo->payloadSize = bytes_written;
      return 0;

   case PAR_OF_26L:


      if (input->symbol_mode == UVLC)
      {
        // Current TML File Format
        // Enforce byte alignment of next header: zero bit stuffing
        currStream = (currSlice->partArr[0]).bitstream;

        if (currStream->bits_to_go < 8)
        { // trailing bits to process
          currStream->byte_buf <<= currStream->bits_to_go;
          stat->bit_use_stuffingBits[img->type]+=currStream->bits_to_go;
          currStream->streamBuffer[currStream->byte_pos++]=currStream->byte_buf;
          currStream->bits_to_go = 8;
        }


        if(!write_out)
        {
          if(input->Encapsulated_NAL_Payload) 
          {
            SODBtoRBSP(currStream);
            temp_byte_pos = currStream->byte_pos;
            currStream->byte_pos = RBSPtoEBSP(currStream->streamBuffer, currStream->last_startcode+startcodeprefix_len, currStream->byte_pos);
            *(stat->em_prev_bits) += (currStream->byte_pos - temp_byte_pos) * 8;
          }
        } 
        else
        {
          bytes_written = currStream->byte_pos;
          stat->bit_ctr += 8*bytes_written;     // actually written bits
          fwrite (currStream->streamBuffer, 1, bytes_written, out);
        }

        currStream->stored_bits_to_go = 8; // store bits_to_go
        currStream->stored_byte_buf   = currStream->byte_buf;   // store current byte
        if(write_out)
          currStream->stored_byte_pos   = 0; // reset byte position
        else
          currStream->stored_byte_pos   = currStream->byte_pos;   // reset byte position

      }
      else
      {
        if(!write_out)
        {
          // CABAC File Format
          write_terminating_bit (1);
          
          eep = &((currSlice->partArr[0]).ee_cabac);
          
          currStream = (currSlice->partArr[0]).bitstream;
          
          // terminate the arithmetic code
          stat->bit_use_stuffingBits[img->type]+=get_trailing_bits(eep);
          arienco_done_encoding(eep);
          
          if (eep->Ebits_to_go != 8)
            stat->bit_use_stuffingBits[img->type]+=eep->Ebits_to_go;
          
          bytes_written = currStream->byte_pos-currStream->tmp_byte_pos;
          if (currStream->bits_to_go < 8) // trailing bits to process
          {
            currStream->byte_buf <<= currStream->bits_to_go;
            stat->bit_use_header[img->type]+=currStream->bits_to_go;
            
            currStream->streamBuffer[bytes_written]= currStream->byte_buf;  // Yue
            bytes_written++;
            currStream->bits_to_go = 8;
          }

          if(input->Encapsulated_NAL_Payload) 
          {
            SODBtoRBSP(currStream);
            temp_byte_pos = currStream->byte_pos;
            currStream->byte_pos = RBSPtoEBSP(currStream->streamBuffer, currStream->last_startcode+startcodeprefix_len, currStream->byte_pos);
            *(stat->em_prev_bits) += (currStream->byte_pos - temp_byte_pos) * 8;
          }
          currStream->tmp_byte_pos = currStream->byte_pos;
        }
        else
        {
          currStream = (currSlice->partArr[0]).bitstream;
          fwrite (currStream->streamBuffer, 1, currStream->byte_pos, out);
          stat->bit_ctr += 8*currStream->byte_pos;
        }
      }
      return 0;

    case PAR_OF_RTP:

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
          stat->bit_use_stuffingBits[img->type]+=get_trailing_bits(eep);
          arienco_done_encoding(eep);
          if (eep->Ebits_to_go != 8)
            stat->bit_use_stuffingBits[img->type]+=eep->Ebits_to_go;
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
          if(input->Encapsulated_NAL_Payload) 
          {              
              currStream->byte_pos = bytes_written; 
              currStream->bits_to_go = 8;
              currStream->byte_buf = 0;
              SODBtoRBSP(currStream);
              bytes_written = currStream->byte_pos;
              temp_byte_pos = bytes_written;
              bytes_written = RBSPtoEBSP(currStream->streamBuffer, Bytes_After_Header, bytes_written);
              *(stat->em_prev_bits) += (bytes_written - temp_byte_pos) * 8;
          }

          // Tian Dong (Sept 2002):
          if (isAggregationPacket())
            rtp_bytes_written = aggregationRTPWriteBits (Marker, FirstBytePacketType, subFirstBytePacketType, currStream->streamBuffer, bytes_written, out);
          else
            rtp_bytes_written = RTPWriteBits (Marker, FirstBytePacketType, currStream->streamBuffer, bytes_written, out);
        }
        stat->bit_ctr += 8*bytes_written;

        // Provide the next partition with a 'fresh' buffer
        currStream->stored_bits_to_go = 8;
        currStream->stored_byte_buf   = 0;
        currStream->stored_byte_pos   = 0;
        currStream->bits_to_go = 8;
        currStream->byte_buf   = 0;
        currStream->byte_pos   = 0;
      }
    return 0;


    default:
      snprintf(errortext, ET_SIZE, "Output File Mode %d not supported", input->of_mode);
      error(errortext,1);
      return 1;
  }
}
