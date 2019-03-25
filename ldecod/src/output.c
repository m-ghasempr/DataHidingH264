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
 * \file output.c
 *
 * \brief
 *    Output an image and Trance support
 *
 * \author
 *    Main contributors (see contributors.h for copyright, address and affiliation details)
 *    - Inge Lille-Langøy               <inge.lille-langoy@telenor.com>
 *    - Rickard Sjoberg                 <rickard.sjoberg@era.ericsson.se>
 ************************************************************************
 */

#include "contributors.h"

#include <stdlib.h>

#include "global.h"

/*!
 ************************************************************************
 * \brief
 *    Write decoded frame to output file
 ************************************************************************
 */
void write_frame(
  struct img_par *img,  //!< image parameter
  FILE *p_out)          //!< filestream to output file
{
  int i,j;

  /*
   * the mmin, mmax macros are taken out, because it makes no sense due to limited range of data type
   */

  for(i=0;i<img->height;i++)
    for(j=0;j<img->width;j++)
    {
      fputc(imgY[i][j],p_out);
    }
  for(i=0;i<img->height_cr;i++)
    for(j=0;j<img->width_cr;j++)
    {
      fputc(imgUV[0][i][j],p_out);
    }
  for(i=0;i<img->height_cr;i++)
    for(j=0;j<img->width_cr;j++)
    {
      fputc(imgUV[1][i][j],p_out);
    }
  fflush(p_out);
}

/*!
 ************************************************************************
 * \brief
 *    Write previous decoded P frame to output file
 ************************************************************************
 */
void write_prev_Pframe(struct img_par *img, FILE *p_out)
{
  int i,j;

  if (img->structure != FRAME)
  {
    img->height *= 2;
    img->height_cr *= 2;
  }
 
  for(i=0;i<img->height;i++)
    for(j=0;j<img->width;j++)
      fputc(imgY_prev[i][j],p_out);

  for(i=0;i<img->height_cr;i++)
    for(j=0;j<img->width_cr;j++)
      fputc(imgUV_prev[0][i][j],p_out);

  for(i=0;i<img->height_cr;i++)
    for(j=0;j<img->width_cr;j++)
      fputc(imgUV_prev[1][i][j],p_out);
  fflush( p_out  );

  if (img->structure != FRAME)
  {
    img->height /= 2;
    img->height_cr /= 2;
  }
  
}



#if TRACE

static int bitcounter = 0;

/*!
 ************************************************************************
 * \brief
 *    Tracing bitpatterns for symbols
 *    A code word has the following format: 0 Xn...0 X2 0 X1 0 X0 1
 ************************************************************************
 */
void tracebits(
    const char *trace_str,  //!< tracing information, char array describing the symbol
    int len,                //!< length of syntax element in bits
    int info,               //!< infoword of syntax element
    int value1)
{

  int i, chars;
  // int outint = 1;

  if(len>=34)
  {
    snprintf(errortext, ET_SIZE, "Length argument to put too long for trace to work");
    error (errortext, 600);
  }


  putc('@', p_trace);
  chars = fprintf(p_trace, "%i", bitcounter);
  while(chars++ < 6)
    putc(' ',p_trace);

  chars += fprintf(p_trace, "%s", trace_str);
  while(chars++ < 55)
    putc(' ',p_trace);

  // Align bitpattern
  if(len<15)
  {
    for(i=0 ; i<15-len ; i++)
      fputc(' ', p_trace);
  }

  // Print bitpattern
  for(i=0 ; i<len/2 ; i++)
  {
    fputc('0', p_trace);
  }
  // put 1
  fprintf(p_trace, "1");

  // Print bitpattern
  for(i=0 ; i<len/2 ; i++)
  {
      if (0x01 & ( info >> ((len/2-i)-1)))
        fputc('1', p_trace);
      else
        fputc('0', p_trace);
  }

  fprintf(p_trace, "  (%3d)\n", value1);
  bitcounter += len;

  fflush (p_trace);

}

/*!
 ************************************************************************
 * \brief
 *    Tracing bitpatterns 
 ************************************************************************
 */
void tracebits2(
    const char *trace_str,  //!< tracing information, char array describing the symbol
    int len,                //!< length of syntax element in bits
    int info)
{

  int i, chars;
  // int outint = 1;

  if(len>=45)
  {
    snprintf(errortext, ET_SIZE, "Length argument to put too long for trace to work");
    error (errortext, 600);
  }


  putc('@', p_trace);
  chars = fprintf(p_trace, "%i", bitcounter);
  while(chars++ < 6)
    putc(' ',p_trace);
  chars += fprintf(p_trace, "%s", trace_str);
  while(chars++ < 55)
    putc(' ',p_trace);

  // Align bitpattern
  if(len<15)
    for(i=0 ; i<15-len ; i++)
      fputc(' ', p_trace);


  bitcounter += len;
  while (len >= 32)
  {
    for(i=0 ; i<8 ; i++)
    {
      fputc('0', p_trace);
    }
    len -= 8;

  }
  // Print bitpattern
  for(i=0 ; i<len ; i++)
  {
    if (0x01 & ( info >> (len-i-1)))
      fputc('1', p_trace);
    else
      fputc('0', p_trace);
  }

  fprintf(p_trace, "  (%3d)\n", info);

  fflush (p_trace);

}


#endif
