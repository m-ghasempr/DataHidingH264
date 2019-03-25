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
#ifndef _MBUFFER_H_
#define _MBUFFER_H_

#include "global.h"

//! definition a picture (field or frame)
typedef struct storable_picture
{
  PictureStructure structure;

  int         poc;

  int         pic_num;
  int         long_term_pic_num;
  int         long_term_frame_idx;

  int         is_long_term;
  int         used_for_reference;

  int         size_x, size_y, size_x_cr, size_y_cr;
  int         chroma_vector_adjustment;
  int         coded_frame;
  int         mb_adaptive_frame_field_flag;

  byte **     imgY;          //!< Y picture component
  byte **     imgY_qpel;       //!< Y picture component upsampled (Quarter pel)
  byte ***    imgUV;         //!< U and V picture components

  byte *      mb_field;      //!< field macroblock indicator

  int  ***    ref_idx;       //<! reference picture   [list][mb_nr][subblock_x][subblock_y]
                             //   [list][mb_nr][subblock_x][subblock_y]
  int  ****   mv;            //<! motion vector       [list][mb_nr][subblock_x][subblock_y]
  
  struct storable_picture *top_field;     // for mb aff, if frame for referencing the top field
  struct storable_picture *bottom_field;  // for mb aff, if frame for referencing the bottom field
  struct storable_picture *frame;         // for mb aff, if field for referencing the combined frame

} StorablePicture;

//! Frame Stores for Decoded Picture Buffer
typedef struct frame_store
{
  int       is_used;                //<! 0=empty; 1=top; 2=bottom; 3=both fields (or frame)
  int       is_reference;           //<! 0=not used for ref; 1=top used; 2=bottom used; 3=both fields (or frame) used
  int       is_long_term;           //<! 0=not used for ref; 1=top used; 2=bottom used; 3=both fields (or frame) used

  int       is_non_existent;

  unsigned  frame_num;
  unsigned  frame_num_wrap;
  int       long_term_frame_idx;
  int       is_output;
  int       poc;

  StorablePicture *frame;
  StorablePicture *top_field;
  StorablePicture *bottom_field;

} FrameStore;


//! Decoded Picture Buffer
typedef struct decoded_picture_buffer
{
  FrameStore  **fs;
  FrameStore  **fs_ref;
  FrameStore  **fs_ltref;
  unsigned      size;
  unsigned      used_size;
  unsigned      ref_frames_in_buffer;
  unsigned      ltref_frames_in_buffer;
  int           last_output_poc;
  int           max_long_term_pic_idx;
} DecodedPictureBuffer;


extern DecodedPictureBuffer dpb;
extern StorablePicture **listX[6];
extern int listXsize[6];

void             init_dpb(InputParameters *inp);
void             free_dpb();
FrameStore*      alloc_frame_store();
void             free_frame_store(FrameStore* f);
StorablePicture* alloc_storable_picture(PictureStructure type, int size_x, int size_y, int size_x_cr, int size_y_cr);
void             free_storable_picture(StorablePicture* p);
void             store_picture_in_dpb(StorablePicture* p);
void             flush_dpb();

void             dpb_split_field(FrameStore *fs);
void             dpb_combine_field(FrameStore *fs);

void             init_lists(int currSliceType, PictureStructure currPicStructure);
void             reorder_ref_pic_list(StorablePicture **list, int *list_size, 
                                      int num_ref_idx_lX_active_minus1, int *remapping_of_pic_nums_idc, 
                                      int *abs_diff_pic_num_minus1, int *long_term_pic_idx);

void             init_mbaff_lists();
void             alloc_ref_pic_list_reordering_buffer(Slice *currSlice);
void             free_ref_pic_list_reordering_buffer(Slice *currSlice);

#endif

