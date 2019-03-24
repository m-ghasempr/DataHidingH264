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
 ***********************************************************************
 *  \file
 *      mbuffer.h
 *
 *  \brief
 *      Frame buffer functions
 *
 *  \author
 *      Main contributors (see contributors.h for copyright, address and affiliation details)
 *      - Karsten Sühring          <suehring@hhi.de>
 ***********************************************************************
 */
#ifndef _BUFFER_H_
#define _BUFFER_H_

#include "global.h"

/*! one single frame */
typedef struct 
{
  int  used;
  int  picID;
  int  lt_picID;
  byte **mref;
  byte ***mcef;
  int  islong;            //<! field is needed for reordering
} Frame;

/* the whole frame buffer structure containing both long and short term frames */
typedef struct
{
  Frame **picbuf_short;
  Frame **picbuf_long;
  int   short_size;
  int   long_size;
  int   short_used;
  int   long_used;
} FrameBuffer;


FrameBuffer *fb;

byte ***mref;                                   //<! these are pointer arrays to the actual structures
byte ****mcef;                                  //<! these are pointer arrays to the actual structures

void init_frame_buffers(struct inp_par *inp, ImageParameters *img);
void free_frame_buffers(struct inp_par *inp, ImageParameters *img);

void init_long_term_buffer(int size, ImageParameters *img);

void assign_long_term_id(int shortID, int longID, ImageParameters *img);

void remove_long_term(int longID);
void remove_short_term(int shortID);

void reset_buffers();

int get_short_ID(int DPN, ImageParameters *img);

void alloc_mref(ImageParameters *img);
void init_mref(ImageParameters *img);
void reorder_mref(ImageParameters *img);

void copy2fb(ImageParameters *img);

#endif