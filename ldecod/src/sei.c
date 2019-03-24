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

#include "contributors.h"

#include <assert.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>
#include <string.h>
#include <ctype.h>

#include "global.h"
#include "memalloc.h"
#include "sei.h"
/*!
 ************************************************************************
 * \file  sei.c
 *
 * \brief
 *    Functions to implement SEI messages
 *
 * \author
 *    Main contributors (see contributors.h for copyright, address and affiliation details)
 *    - Dong Tian   <tian@cs.tut.fi>
 ************************************************************************
 */


 /*!
 ************************************************************************
 *  \basic functions on supplemental enhancement information
 *  \brief
 *     The implementations are based on FCD
 ************************************************************************
 */

/*!
 ************************************************************************
 *  \InterpretSEIMessage
 *  \brief
 *     Interpret the SEI rbsp
 *  \input
 *    msg:      a pointer that point to the sei message.
 *    size:     the size of the sei message
 *    img:      the image pointer
 *  \output
 *    
 ************************************************************************
 */
void InterpretSEIMessage(byte* msg, int size, struct img_par *img)
{
  int payload_type = 0;
  int payload_size = 0;
  int offset = 0;
  byte tmp_byte;
  do
  {
    // sei_message();
    payload_type = 0;
    tmp_byte = msg[offset++];
    while (tmp_byte == 0xFF)
    {
      payload_type += 255;
      tmp_byte = msg[offset++];
    }
    payload_type += tmp_byte;   // this is the last byte

    payload_size = 0;
    tmp_byte = msg[offset++];
    while (tmp_byte == 0xFF)
    {
      payload_size += 255;
      tmp_byte = msg[offset++];
    }
    payload_size += tmp_byte;   // this is the last byte
    
    switch ( payload_type )     // sei_payload( type, size );
    {
    case SEI_SPARE_PICTURE:
      interpret_spare_picture( msg+offset, payload_size, img );
      break;
    case SEI_SUBSEQ_INFORMATION:
      interpret_subsequence_info( msg+offset, payload_size, img );
      break;
    case SEI_SUBSEQ_LAYER_CHARACTERISTICS:
      interpret_subsequence_layer_info( msg+offset, payload_size, img );
      break;
    case SEI_SUBSEQ_CHARACTERISTICS:
      interpret_subsequence_characteristics_info( msg+offset, payload_size, img );
      break;
    default:
      printf("The SEI type %d is not supported in JM yet.\n", payload_type);
      exit(0);
      break;
    }
    offset += payload_size;
    
  } while( msg[offset] != 0x80 );    // more_rbsp_data()  msg[offset] != 0x80
  // ignore the trailing bits rbsp_trailing_bits();
  assert(msg[offset] == 0x80);      // this is the trailing bits
  assert( offset+1 == size );
}

void interpret_spare_picture( byte* payload, int size, struct img_par *img )
{
  int i,x,y;
  SyntaxElement sym;
  Bitstream* buf;
  int bit0, bit1, bitc, no_bit0;
  int delta_frame_num, CurrFrameNum, TargetFrameNum;
  int num_spare_pics;
  int delta_spare_frame_num, CandidateSpareFrameNum, SpareFrameNum;
  int ref_area_indicator;

  int m, n, left, right, top, bottom,directx, directy;
  byte ***map;

// #define WRITE_MAP_IMAGE

#ifdef WRITE_MAP_IMAGE
  int  j, k, i0, j0, tmp, kk;
  char filename[20] = "map_dec.yuv";
  FILE *fp;
  byte** Y;
  static int old_pn=-1;
  static int first = 1;
#endif

  assert( payload!=NULL);
  assert( img!=NULL);

  CurrFrameNum = img->number % MAX_FN;
  // UVLC coded map:
  sym.type = SE_HEADER;       // This will be true for all symbols generated here
  sym.mapping = linfo;        // Mapping rule: Simple code number to len/info

  buf = malloc(sizeof(Bitstream));
  buf->bitstream_length = size;
  buf->streamBuffer = payload;
  buf->frame_bitoffset = 0;

  sym.len = GetVLCSymbol(buf->streamBuffer, buf->frame_bitoffset, &sym.inf, buf->bitstream_length );
  sym.mapping(sym.len, sym.inf, &(sym.value1), &(sym.value2));
  buf->frame_bitoffset += sym.len;
  delta_frame_num = sym.value1;
  TargetFrameNum = CurrFrameNum - delta_frame_num;
  if( TargetFrameNum < 0 )
    TargetFrameNum = MAX_FN + TargetFrameNum;
#ifdef WRITE_MAP_IMAGE
  printf( "TargetFrameNum is %d\n", TargetFrameNum );
#endif

  sym.len = GetVLCSymbol(buf->streamBuffer, buf->frame_bitoffset, &sym.inf, buf->bitstream_length );
  sym.mapping(sym.len, sym.inf, &(sym.value1), &(sym.value2));
  buf->frame_bitoffset += sym.len;
  num_spare_pics = sym.value1 + 1;
#ifdef WRITE_MAP_IMAGE
  printf( "num_spare_pics is %d\n", num_spare_pics );
#endif

  get_mem3D(&map, num_spare_pics, img->height/16, img->width/16);

  for (i=0; i<num_spare_pics; i++)
  {
    if (i==0) 
    {
      CandidateSpareFrameNum = TargetFrameNum - 1;
      if ( CandidateSpareFrameNum < 0 ) CandidateSpareFrameNum = MAX_FN - 1;
    }
    else
      CandidateSpareFrameNum = SpareFrameNum;

    sym.len = GetVLCSymbol(buf->streamBuffer, buf->frame_bitoffset, &sym.inf, buf->bitstream_length );
    sym.mapping(sym.len, sym.inf, &(sym.value1), &(sym.value2));
    buf->frame_bitoffset += sym.len;
    delta_spare_frame_num = sym.value1;

    SpareFrameNum = CandidateSpareFrameNum - delta_spare_frame_num;
    if( SpareFrameNum < 0 )
      SpareFrameNum = MAX_FN + SpareFrameNum;

    sym.len = GetVLCSymbol(buf->streamBuffer, buf->frame_bitoffset, &sym.inf, buf->bitstream_length );
    sym.mapping(sym.len, sym.inf, &(sym.value1), &(sym.value2));
    buf->frame_bitoffset += sym.len;
    ref_area_indicator = sym.value1;

    switch ( ref_area_indicator )
    {
    case 0:   // The whole frame can serve as spare picture
      for (y=0; y<img->height/16; y++)
        for (x=0; x<img->width/16; x++)
          map[i][y][x] = 0;
      break;
    case 1:   // The map is not compressed
      for (y=0; y<img->height/16; y++)
        for (x=0; x<img->width/16; x++)
        {
          sym.len = GetfixedSymbol(buf->streamBuffer, buf->frame_bitoffset, &sym.inf, buf->bitstream_length, 1);
          assert( sym.len == 1 );
          map[i][y][x] = sym.inf;
          buf->frame_bitoffset += sym.len;
        }
      break;
    case 2:   // The map is compressed
      bit0 = 0;
      bit1 = 1;
      bitc = bit0;
      no_bit0 = -1;

      x = ( img->width/16 - 1 ) / 2;
      y = ( img->height/16 - 1 ) / 2;
      left = right = x;
      top = bottom = y;
      directx = 0;
      directy = 1;

      for (m=0; m<img->height/16; m++)
        for (n=0; n<img->width/16; n++)
        {

          if (no_bit0<0)
          {
            sym.len = GetVLCSymbol( buf->streamBuffer, buf->frame_bitoffset, &sym.inf, buf->bitstream_length );
            sym.mapping(sym.len, sym.inf, &(sym.value1), &(sym.value2));
            buf->frame_bitoffset += sym.len;
            no_bit0 = sym.value1;
          }
          if (no_bit0>0) map[i][y][x] = bit0;
          else map[i][y][x] = bit1;
          no_bit0--;

          // go to the next mb:
          if ( directx == -1 && directy == 0 )
          {
            if (x > left) x--;
            else if (x == 0)
            {
              y = bottom + 1;
              bottom++;
              directx = 1;
              directy = 0;
            }
            else if (x == left)
            {
              x--;
              left--;
              directx = 0;
              directy = 1;
            }
          }
          else if ( directx == 1 && directy == 0 )
          {
            if (x < right) x++;
            else if (x == img->width/16 - 1)
            {
              y = top - 1;
              top--;
              directx = -1;
              directy = 0;
            }
            else if (x == right)
            {
              x++;
              right++;
              directx = 0;
              directy = -1;
            }
          }
          else if ( directx == 0 && directy == -1 )
          {
            if ( y > top) y--;
            else if (y == 0)
            {
              x = left - 1;
              left--;
              directx = 0;
              directy = 1;
            }
            else if (y == top)
            {
              y--;
              top--;
              directx = -1;
              directy = 0;
            }
          }
          else if ( directx == 0 && directy == 1 )
          {
            if (y < bottom) y++;
            else if (y == img->height/16 - 1)
            {
              x = right+1;
              right++;
              directx = 0;
              directy = -1;
            }
            else if (y == bottom)
            {
              y++;
              bottom++;
              directx = 1;
              directy = 0;
            }
          }


        }
      break;
    default:
      printf( "Wrong ref_area_indicator %d!\n", ref_area_indicator );
      exit(0);
      break;
    }

  } // end of num_spare_pics

#ifdef WRITE_MAP_IMAGE
  // begin to write map seq
  if ( old_pn != img->number )
  {
    old_pn = img->number;
    get_mem2D(&Y, img->height, img->width);
    if (first)
    {
      fp = fopen( filename, "wb" );
      first = 0;
    }
    else
      fp = fopen( filename, "ab" );
    assert( fp != NULL );
    for (kk=0; kk<num_spare_pics; kk++)
    {
      for (i=0; i < img->height/16; i++)
        for (j=0; j < img->width/16; j++)
        {
          tmp=map[kk][i][j]==0? 255 : 0;
          for (i0=0; i0<16; i0++)
            for (j0=0; j0<16; j0++)
              Y[i*16+i0][j*16+j0]=tmp;
        }

      // write the map image
      for (i=0; i < img->height; i++)
        for (j=0; j < img->width; j++)
          fputc(Y[i][j], fp);

      for (k=0; k < 2; k++)
        for (i=0; i < img->height/2; i++)
          for (j=0; j < img->width/2; j++)
            fputc(128, fp);
    }
    fclose( fp );
    free_mem2D( Y );
  }
  // end of writing map image
#undef WRITE_MAP_IMAGE
#endif

  free_mem3D( map, num_spare_pics );

  free(buf);
}

void interpret_subsequence_info( byte* payload, int size, struct img_par *img )
{
  SyntaxElement sym;
  Bitstream* buf;
  int subseq_layer_num, subseq_id, last_picture_flag, stored_frame_cnt;
// #define PRINT_SUBSEQUENCE_INFO

  sym.type = SE_HEADER;
  sym.mapping = linfo;

  buf = malloc(sizeof(Bitstream));
  buf->bitstream_length = size;
  buf->streamBuffer = payload;
  buf->frame_bitoffset = 0;

  sym.len = GetVLCSymbol(buf->streamBuffer, buf->frame_bitoffset, &sym.inf, buf->bitstream_length );
  sym.mapping(sym.len, sym.inf, &(sym.value1), &(sym.value2));
  buf->frame_bitoffset += sym.len;
  subseq_layer_num = sym.value1;
#ifdef PRINT_SUBSEQUENCE_INFO
  printf("subseq_layer_num = %d\n", subseq_layer_num );
#endif

  sym.len = GetVLCSymbol(buf->streamBuffer, buf->frame_bitoffset, &sym.inf, buf->bitstream_length );
  sym.mapping(sym.len, sym.inf, &(sym.value1), &(sym.value2));
  buf->frame_bitoffset += sym.len;
  subseq_id = sym.value1;
#ifdef PRINT_SUBSEQUENCE_INFO
  printf("subseq_id = %d\n", subseq_id);
#endif

  sym.len = GetfixedSymbol(buf->streamBuffer, buf->frame_bitoffset, &sym.inf, buf->bitstream_length, 1);
  assert( sym.len == 1 );
  buf->frame_bitoffset += sym.len;
  last_picture_flag = sym.inf;
#ifdef PRINT_SUBSEQUENCE_INFO
  printf("last_picture_flag = %d\n", last_picture_flag);
#endif

  if (buf->frame_bitoffset/8 < size)    // more_rbsp_data()
  {
    sym.len = GetVLCSymbol(buf->streamBuffer, buf->frame_bitoffset, &sym.inf, buf->bitstream_length );
    sym.mapping(sym.len, sym.inf, &(sym.value1), &(sym.value2));
    buf->frame_bitoffset += sym.len;
    stored_frame_cnt = sym.value1;
#ifdef PRINT_SUBSEQUENCE_INFO
    printf("stored_frame_cnt = %d\n", stored_frame_cnt);
#endif
  }

  free(buf);

#ifdef PRINT_SUBSEQUENCE_INFO
#undef PRINT_SUBSEQUENCE_INFO
#endif
}

void interpret_subsequence_layer_info( byte* payload, int size, struct img_par *img )
{
  int offset = 0;
  unsigned short bitrate, framerate;
  int i = 0;
// #define PRINT_SUBSEQUENCE_CHAR

  while (offset < size)
  {
    bitrate = *(unsigned short*)&payload[offset];
    offset += 2;
    framerate = *(unsigned short*)&payload[offset];
    offset += 2;
#ifdef PRINT_SUBSEQUENCE_CHAR
    printf("layer %d, bitrate = %d, framerate = %d\n", i++, bitrate, framerate);
#endif
  }
}

void interpret_subsequence_characteristics_info( byte* payload, int size, struct img_par *img )
{
  SyntaxElement sym;
  Bitstream* buf;
  int i;
  int subseq_layer_num, subseq_id, duration_flag, subseq_duration, average_rate_flag, average_bit_rate, average_frame_rate;
  int num_referenced_subseqs, ref_subseq_layer_num, ref_subseq_id;
// #define PRINT_SUBSEQUENCE_CHAR

  sym.type = SE_HEADER;
  sym.mapping = linfo;

  buf = malloc(sizeof(Bitstream));
  buf->bitstream_length = size;
  buf->streamBuffer = payload;
  buf->frame_bitoffset = 0;

  sym.len = GetVLCSymbol(buf->streamBuffer, buf->frame_bitoffset, &sym.inf, buf->bitstream_length );
  sym.mapping(sym.len, sym.inf, &(sym.value1), &(sym.value2));
  buf->frame_bitoffset += sym.len;
  subseq_layer_num = sym.value1;
#ifdef PRINT_SUBSEQUENCE_CHAR
  printf("subseq_layer_num = %d\n", subseq_layer_num );
#endif

  sym.len = GetVLCSymbol(buf->streamBuffer, buf->frame_bitoffset, &sym.inf, buf->bitstream_length );
  sym.mapping(sym.len, sym.inf, &(sym.value1), &(sym.value2));
  buf->frame_bitoffset += sym.len;
  subseq_id = sym.value1;
#ifdef PRINT_SUBSEQUENCE_CHAR
  printf("subseq_id = %d\n", subseq_id);
#endif

  sym.len = GetfixedSymbol(buf->streamBuffer, buf->frame_bitoffset, &sym.inf, buf->bitstream_length, 1);
  assert( sym.len == 1 );
  buf->frame_bitoffset += sym.len;
  duration_flag = sym.inf;
#ifdef PRINT_SUBSEQUENCE_CHAR
  printf("duration_flag = %d\n", duration_flag);
#endif

  if ( duration_flag )
  {
    sym.len = GetfixedSymbol(buf->streamBuffer, buf->frame_bitoffset, &sym.inf, buf->bitstream_length, 32);
    assert( sym.len == 32 );
    buf->frame_bitoffset += sym.len;
    subseq_duration = sym.inf;
#ifdef PRINT_SUBSEQUENCE_CHAR
    printf("subseq_duration = %d\n", subseq_duration);
#endif
  }

  sym.len = GetfixedSymbol(buf->streamBuffer, buf->frame_bitoffset, &sym.inf, buf->bitstream_length, 1);
  assert( sym.len == 1 );
  buf->frame_bitoffset += sym.len;
  average_rate_flag = sym.inf;
#ifdef PRINT_SUBSEQUENCE_CHAR
  printf("average_rate_flag = %d\n", average_rate_flag);
#endif

  if ( average_rate_flag )
  {
    sym.len = GetfixedSymbol(buf->streamBuffer, buf->frame_bitoffset, &sym.inf, buf->bitstream_length, 16);
    assert( sym.len == 16 );
    buf->frame_bitoffset += sym.len;
    average_bit_rate = sym.inf;
#ifdef PRINT_SUBSEQUENCE_CHAR
    printf("average_bit_rate = %d\n", average_bit_rate);
#endif

    sym.len = GetfixedSymbol(buf->streamBuffer, buf->frame_bitoffset, &sym.inf, buf->bitstream_length, 16);
    assert( sym.len == 16 );
    buf->frame_bitoffset += sym.len;
    average_frame_rate = sym.inf;
#ifdef PRINT_SUBSEQUENCE_CHAR
    printf("average_frame_rate = %d\n", average_frame_rate);
#endif
  }

  sym.len = GetVLCSymbol(buf->streamBuffer, buf->frame_bitoffset, &sym.inf, buf->bitstream_length );
  sym.mapping(sym.len, sym.inf, &(sym.value1), &(sym.value2));
  buf->frame_bitoffset += sym.len;
  num_referenced_subseqs = sym.value1;
#ifdef PRINT_SUBSEQUENCE_CHAR
  printf("num_referenced_subseqs = %d\n", num_referenced_subseqs);
#endif

  for (i=0; i<num_referenced_subseqs; i++)
  {
    sym.len = GetVLCSymbol(buf->streamBuffer, buf->frame_bitoffset, &sym.inf, buf->bitstream_length );
    sym.mapping(sym.len, sym.inf, &(sym.value1), &(sym.value2));
    buf->frame_bitoffset += sym.len;
    ref_subseq_layer_num = sym.value1;
#ifdef PRINT_SUBSEQUENCE_CHAR
    printf("ref_subseq_layer_num = %d\n", ref_subseq_layer_num);
#endif

    sym.len = GetVLCSymbol(buf->streamBuffer, buf->frame_bitoffset, &sym.inf, buf->bitstream_length );
    sym.mapping(sym.len, sym.inf, &(sym.value1), &(sym.value2));
    buf->frame_bitoffset += sym.len;
    ref_subseq_id = sym.value1;
#ifdef PRINT_SUBSEQUENCE_CHAR
    printf("ref_subseq_id = %d\n", ref_subseq_id);
#endif
  }

  free( buf );
#ifdef PRINT_SUBSEQUENCE_CHAR
#undef PRINT_SUBSEQUENCE_CHAR
#endif
}