/*
***********************************************************************
*  COPYRIGHT  AND  WARRANTY INFORMATION
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

/*
 *************************************************************************************
 * \file
 *    golomb_dec.c
 *
 * \brief
 *    Description
 *
 * \author
 *    Main contributors (see contributors.h for copyright, address and affiliation details)
 *     -  Achim Dahlhoff      <dahlhoff@ient.rwth-aachen.de>
 *
 * \date
 *    Fri Mar 8 2002
 *
 *  copyright : (C) 2002      Institut und Lehrstuhl für Nachrichtentechnik
 *                            RWTH Aachen University
 *                            52072 Aachen
 *                            Germany
 *************************************************************************************
 */

#include <assert.h>
#include "golomb_dec.h"

unsigned int decode_golomb_word(const unsigned char **buffer,unsigned int *bitoff,unsigned int grad0,unsigned int max_levels)
{
  const unsigned char *rd;
  unsigned int bit,byte,level,databits,t,testbit;
  rd=*buffer;
  bit=*bitoff;
  byte=*rd;
  level=0UL;
  while( level+1UL<max_levels )
  {
    testbit=byte&(1UL<<bit);
    bit = (bit-1UL) & 7UL ;
    if(bit==7UL)byte=*(++rd);
    if( testbit )break;
    level++;
  }
  databits=0UL;
  for( t=0UL ; t<(grad0+level) ; t++ )
  {
    databits = (databits<<1UL) | ((byte>>bit)&1UL) ;
    bit = (bit-1UL) & 7UL ;
    if(bit==7UL)byte=*(++rd);
  }
  *buffer=rd;
  *bitoff=bit;
  return (((1UL<<level)-1UL)<<grad0)+databits;
}

unsigned int decode_multilayer_golomb_word(const unsigned char **buffer,unsigned int *bitoff,const unsigned int *grad0,const unsigned int *max_levels)
{
  unsigned int symbol,partsymbol,tmp;
  symbol=0UL;
  while(1)
  {
    partsymbol=decode_golomb_word(buffer,bitoff,*grad0,*max_levels);
    symbol+=partsymbol;
    tmp=*max_levels;
    if( partsymbol < (((1UL<<tmp)-1UL)<<(*grad0))-1UL )  //not escape symbol?
      break;
    grad0++;
    max_levels++;
  }
  return symbol;
}








int  readSyntaxElement_GOLOMB(SyntaxElement *se, struct img_par *img, struct inp_par *inp, struct datapartition *dp)
{
 Bitstream *currStream;
 int frame_bitoffset;
 unsigned char *buf,*read;
 int BitstreamLengthInBytes;
 unsigned int bit,i;
 unsigned int grad[4],max_lev[4];

  currStream = dp->bitstream;
  frame_bitoffset = currStream->frame_bitoffset;
  buf = (unsigned char*)currStream->streamBuffer;
  BitstreamLengthInBytes = currStream->bitstream_length;

  bit=7UL-(frame_bitoffset&7);
  read=buf+(frame_bitoffset>>3);

  if(!( se->golomb_maxlevels&~0xFF ))
  {
    se->value1=decode_golomb_word(&read,&bit,se->golomb_grad,se->golomb_maxlevels);
  }else{
    for(i=0UL;i<4UL;i++)
    {
      grad[i]=(se->golomb_grad>>(i<<3))&0xFFUL;
      max_lev[i]=(se->golomb_maxlevels>>(i<<3))&0xFFUL;
    }
    se->value1=decode_multilayer_golomb_word(&read,&bit,grad,max_lev);
  }
  se->len=(((read-buf)<<3)+(7-bit))-frame_bitoffset;

  se->value2=0;

  currStream->frame_bitoffset += se->len;

  return 1;
}
