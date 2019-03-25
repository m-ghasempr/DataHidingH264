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

#include "global.h"
#if TRACE
#include <string.h>    // strncpy
#endif
#include "rtp.h"


/*!
 ************************************************************************
 * \brief
 *    Error handling procedure. Print error message to stderr and exit
 *    with supplied code.
 * \param text
 *    Error message
 ************************************************************************
 */
void error(char *text, int code)
{
  fprintf(stderr, "%s\n", text);
  exit(code);
}


/*!
 ************************************************************************
 *  \brief
 *     This function generates the appropriate slice
 *     header
 ************************************************************************
 */
void start_slice(struct img_par *img, struct inp_par *inp)
{
  Slice *currSlice = img->currentSlice;
/*
  switch(inp->FileFormat)
  {
    case PAR_OF_ANNEXB:
      currSlice->dp_mode = PAR_DP_1; //other modes not supported

      if (inp->symbol_mode == UVLC)
      {
        // Current TML File Format
        nal_startcode_follows = uvlc_startcode_follows;
//!        currSlice->readSlice = readSliceUVLC;
        currSlice->partArr[0].readSyntaxElement = readSyntaxElement_UVLC;
      }
      else
      {
        // CABAC File Format
        nal_startcode_follows = cabac_startcode_follows;
        currSlice->readSlice = readSliceCABAC;
        currSlice->partArr[0].readSyntaxElement = readSyntaxElement_CABAC;
      }
      break;
    case PAR_OF_RTP:
      if (inp->symbol_mode == UVLC)
      {
        nal_startcode_follows = RTP_startcode_follows;
        currSlice->readSlice = readSliceRTP;
        for (i=0; i<3; i++)       // always up to three partitions in RTP
        {
          currSlice->partArr[i].readSyntaxElement = readSyntaxElement_RTP;
          currSlice->partArr[i].bitstream->ei_flag = 1;
        }
      }
      else
      {
        // CABAC File Format
        nal_startcode_follows = cabac_startcode_follows;
        currSlice->readSlice = readSliceRTP;
        for (i=0; i<currSlice->max_part_nr; i++)
        {
          currSlice->partArr[i].readSyntaxElement = readSyntaxElement_CABAC;
          currSlice->partArr[i].bitstream->ei_flag = 1;
        }
      }
      break;
    default:
      snprintf(errortext, ET_SIZE, "Input File Mode %d not supported", inp->FileFormat);
      error(errortext,1);
      break;
  }
  */
}


