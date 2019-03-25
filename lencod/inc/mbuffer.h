
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
 *      - Alexis Michael Tourapis  <alexismt@ieee.org>
 ***********************************************************************
 */
#ifndef _MBUFFER_H_
#define _MBUFFER_H_

#include "global.h"
#include "enc_statistics.h"

#define MAX_LIST_SIZE 33

typedef struct picture_stats
{
  double dsum[3];
  double dvar[3];
} PictureStats;

//! definition a picture (field or frame)
struct storable_picture
{
  PictureStructure structure;

  int         poc;
  int         top_poc;
  int         bottom_poc;
  int         frame_poc;
  int         order_num;
  int64       ref_pic_num[6][MAX_LIST_SIZE];
  int64       frm_ref_pic_num[6][MAX_LIST_SIZE];
  int64       top_ref_pic_num[6][MAX_LIST_SIZE];
  int64       bottom_ref_pic_num[6][MAX_LIST_SIZE];
  unsigned    frame_num;
  int         pic_num;
  int         long_term_pic_num;
  int         long_term_frame_idx;

  byte        is_long_term;
  int         used_for_reference;
  int         is_output;
  int         non_existing;

  int         size_x, size_y, size_x_cr, size_y_cr;
  int         size_x_padded, size_y_padded;
  int         size_x_pad, size_y_pad;
  int         size_x_cr_pad, size_y_cr_pad;
  int         chroma_vector_adjustment;
  int         coded_frame;
  int         mb_aff_frame_flag;

  imgpel **   imgY;          //!< Y picture component
  imgpel **** imgY_sub;      //!< Y picture component upsampled (Quarter pel)
  imgpel ***  imgUV;         //!< U and V picture components
  imgpel *****imgUV_sub;     //!< UV picture component upsampled (Quarter/One-Eighth pel)

  //Multiple Decoder Buffers (used if rdopt==3)
  imgpel ***  dec_imgY;       //!< Decoded Y component in multiple hypothetical decoders
  imgpel **** dec_imgUV;      //!< Decoded U and V components in multiple hypothetical decoders
  imgpel ***  p_dec_img[MAX_PLANE];      //!< pointer array for accessing decoded pictures in hypothetical decoders

  byte   ***  mb_error_map;    //!< Map of macroblock errors in hypothetical decoders.

  imgpel **   p_img[MAX_PLANE];          //!< pointer array for accessing imgY/imgUV[]
  imgpel **** p_img_sub[MAX_PLANE];      //!< pointer array for storing top address of imgY_sub/imgUV_sub[]
  imgpel **   p_curr_img;                //!< current int-pel ref. picture area to be used for motion estimation
  imgpel **** p_curr_img_sub;            //!< current sub-pel ref. picture area to be used for motion estimation

  //PicMotionParams2 ***mv_info;    //!< Motion info
  PicMotionParams  motion;    //!< Motion info
  PicMotionParams JVmotion[MAX_PLANE];    //!< Motion info for 4:4:4 independent coding

  int colour_plane_id;                     //!< colour_plane_id to be used for 4:4:4 independent mode encoding

  struct storable_picture *top_field;     // for mb aff, if frame for referencing the top field
  struct storable_picture *bottom_field;  // for mb aff, if frame for referencing the bottom field
  struct storable_picture *frame;         // for mb aff, if field for referencing the combined frame

  int         chroma_format_idc;
  int         chroma_mask_mv_x;
  int         chroma_mask_mv_y;
  int         chroma_shift_y;
  int         chroma_shift_x;
  int         frame_mbs_only_flag;
  int         frame_cropping_flag;
  int         frame_cropping_rect_left_offset;
  int         frame_cropping_rect_right_offset;
  int         frame_cropping_rect_top_offset;
  int         frame_cropping_rect_bottom_offset;

  PictureStats p_stats;
  StatParameters stats;

  int         type;
};

typedef struct storable_picture StorablePicture;

//! definition of motion parameters
struct motion_params
{
  int64 ***   ref_pic_id;    //!< reference picture identifier [list][subblock_y][subblock_x]
  short ****  mv;            //!< motion vector       [list][subblock_x][subblock_y][component]
  char  ***   ref_idx;       //!< reference picture   [list][subblock_y][subblock_x]
  byte **     moving_block;
};

//! definition a picture (field or frame)
struct colocated_params
{
  int         mb_adaptive_frame_field_flag;
  int         size_x, size_y;
  byte        is_long_term;

  struct motion_params frame;
  struct motion_params top;
  struct motion_params bottom;

};

typedef struct motion_params MotionParams;
typedef struct colocated_params ColocatedParams;

//! Frame Stores for Decoded Picture Buffer
struct frame_store
{
  int       is_used;                //!< 0=empty; 1=top; 2=bottom; 3=both fields (or frame)
  int       is_reference;           //!< 0=not used for ref; 1=top used; 2=bottom used; 3=both fields (or frame) used
  int       is_long_term;           //!< 0=not used for ref; 1=top used; 2=bottom used; 3=both fields (or frame) used
  int       is_orig_reference;      //!< original marking by nal_ref_idc: 0=not used for ref; 1=top used; 2=bottom used; 3=both fields (or frame) used

  int       is_non_existent;

  unsigned  frame_num;
  int       frame_num_wrap;
  int       long_term_frame_idx;
  int       is_output;
  int       poc;

  StorablePicture *frame;
  StorablePicture *top_field;
  StorablePicture *bottom_field;
};

typedef struct frame_store FrameStore;


//! Decoded Picture Buffer
struct decoded_picture_buffer
{
  VideoParameters *p_Vid;
  InputParameters *p_Inp;
  FrameStore  **fs;
  FrameStore  **fs_ref;
  FrameStore  **fs_ltref;
  unsigned      size;
  unsigned      used_size;
  unsigned      ref_frames_in_buffer;
  unsigned      ltref_frames_in_buffer;
  int           last_output_poc;
  int           max_long_term_pic_idx;

  int           init_done;

  FrameStore   *last_picture;
};

typedef struct decoded_picture_buffer DecodedPictureBuffer;

extern void             init_dpb                  (VideoParameters *p_Vid, DecodedPictureBuffer *dpb);
extern void             free_dpb                  (VideoParameters *p_Vid, DecodedPictureBuffer *dpb);
extern FrameStore*      alloc_frame_store(void);
extern void             free_frame_store          (VideoParameters *p_Vid, FrameStore* f);
extern StorablePicture* alloc_storable_picture    (VideoParameters *p_Vid, PictureStructure type, int size_x, int size_y, int size_x_cr, int size_y_cr);
extern void             free_storable_picture     (VideoParameters *p_Vid, StorablePicture* p);
extern void             store_picture_in_dpb      (VideoParameters *p_Vid, StorablePicture* p, FrameFormat *output);
extern void             replace_top_pic_with_frame(VideoParameters *p_Vid, StorablePicture* p, FrameFormat *output);
extern void             flush_dpb                 (VideoParameters *p_Vid, FrameFormat *output);
extern void             dpb_split_field           (VideoParameters *p_Vid, FrameStore *fs);
extern void             dpb_combine_field         (VideoParameters *p_Vid, FrameStore *fs);
extern void             dpb_combine_field_yuv     (VideoParameters *p_Vid, FrameStore *fs);
extern void             init_lists                (Slice *currSlice);
extern void             init_lists_single_dir     (Slice *currSlice);
extern void             reorder_ref_pic_list      (Slice *currSlice, StorablePicture **list[6], char list_size[6], int cur_list);
extern void             init_mbaff_lists          (Slice *currSlice);
extern void             alloc_ref_pic_list_reordering_buffer (Slice *currSlice);
extern void             free_ref_pic_list_reordering_buffer  (Slice *currSlice);
extern void             fill_frame_num_gap        (VideoParameters *p_Vid, FrameFormat *output);
extern ColocatedParams* alloc_colocated           (int size_x, int size_y,int mb_adaptive_frame_field_flag);
extern void             free_colocated            (ColocatedParams* p);
extern void             compute_colocated         (Slice *currSlice, ColocatedParams* p, StorablePicture **listX[6]);

// For 4:4:4 independent mode
extern void             compute_colocated_JV      ( Slice *currSlice, ColocatedParams* p, StorablePicture **listX[6]);
extern void             copy_storable_param_JV    ( VideoParameters *p_Vid, int nplane, StorablePicture *d, StorablePicture *s );

#endif

