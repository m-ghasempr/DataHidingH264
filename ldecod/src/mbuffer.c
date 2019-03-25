
/*!
 ***********************************************************************
 *  \file
 *      mbuffer.c
 *
 *  \brief
 *      Frame buffer functions
 *
 *  \author
 *      Main contributors (see contributors.h for copyright, address and affiliation details)
 *      - Karsten Sühring                 <suehring@hhi.de>
 *      - Alexis Tourapis                 <alexismt@ieee.org>
 *      - Jill Boyce                      <jill.boyce@thomson.net>
 *      - Saurav K Bandyopadhyay          <saurav@ieee.org>
 *      - Zhenyu Wu                       <Zhenyu.Wu@thomson.net
 *      - Purvin Pandit                   <Purvin.Pandit@thomson.net>
 *
 ***********************************************************************
 */

#include <limits.h>

#include "global.h"
#include "erc_api.h"
#include "header.h"
#include "image.h"
#include "mbuffer.h"
#include "memalloc.h"
#include "output.h"



static void insert_picture_in_dpb    (ImageParameters *p_Img, FrameStore* fs, StorablePicture* p);
static void output_one_frame_from_dpb(ImageParameters *p_Img);
static void get_smallest_poc         (DecodedPictureBuffer *p_Dpb, int *poc,int * pos);
static void gen_field_ref_ids        (StorablePicture *p);
static int  remove_unused_frame_from_dpb (ImageParameters *p_Img, DecodedPictureBuffer *p_Dpb);
static int  is_used_for_reference    (FrameStore* fs);
static int  is_short_term_reference  (FrameStore* fs);
static int  is_long_term_reference   (FrameStore* fs);

#define MAX_LIST_SIZE 33

/*!
 ************************************************************************
 * \brief
 *    Print out list of pictures in DPB. Used for debug purposes.
 ************************************************************************
 */
void dump_dpb(DecodedPictureBuffer *p_Dpb)
{
#if DUMP_DPB
  unsigned i;

  for (i=0; i<p_Dpb->used_size;i++)
  {
    printf("(");
    printf("fn=%d  ", p_Dpb->fs[i]->frame_num);
    if (p_Dpb->fs[i]->is_used & 1)
    {
      if (p_Dpb->fs[i]->top_field)
        printf("T: poc=%d  ", p_Dpb->fs[i]->top_field->poc);
      else
        printf("T: poc=%d  ", p_Dpb->fs[i]->frame->top_poc);
    }
    if (p_Dpb->fs[i]->is_used & 2)
    {
      if (p_Dpb->fs[i]->bottom_field)
        printf("B: poc=%d  ", p_Dpb->fs[i]->bottom_field->poc);
      else
        printf("B: poc=%d  ", p_Dpb->fs[i]->frame->bottom_poc);
    }
    if (p_Dpb->fs[i]->is_used == 3)
      printf("F: poc=%d  ", p_Dpb->fs[i]->frame->poc);
    printf("G: poc=%d)  ", p_Dpb->fs[i]->poc);
    if (p_Dpb->fs[i]->is_reference) printf ("ref (%d) ", p_Dpb->fs[i]->is_reference);
    if (p_Dpb->fs[i]->is_long_term) printf ("lt_ref (%d) ", p_Dpb->fs[i]->is_reference);
    if (p_Dpb->fs[i]->is_output) printf ("out  ");
    if (p_Dpb->fs[i]->is_used == 3)
    {
      if (p_Dpb->fs[i]->frame->non_existing) printf ("ne  ");
    }
    printf ("\n");
  }
#endif
}

/*!
 ************************************************************************
 * \brief
 *    Returns the size of the dpb depending on level and picture size
 *
 *
 ************************************************************************
 */
int getDpbSize(seq_parameter_set_rbsp_t *active_sps)
{
  int pic_size = (active_sps->pic_width_in_mbs_minus1 + 1) * (active_sps->pic_height_in_map_units_minus1 + 1) * (active_sps->frame_mbs_only_flag?1:2) * 384;

  int size = 0;

  switch (active_sps->level_idc)
  {
  case 9:
    size = 152064;
    break;
  case 10:
    size = 152064;
    break;
  case 11:
    if (!IS_FREXT_PROFILE(active_sps->profile_idc) && (active_sps->constrained_set3_flag == 1))
      size = 152064;
    else
      size = 345600;
    break;
  case 12:
    size = 912384;
    break;
  case 13:
    size = 912384;
    break;
  case 20:
    size = 912384;
    break;
  case 21:
    size = 1824768;
    break;
  case 22:
    size = 3110400;
    break;
  case 30:
    size = 3110400;
    break;
  case 31:
    size = 6912000;
    break;
  case 32:
    size = 7864320;
    break;
  case 40:
    size = 12582912;
    break;
  case 41:
    size = 12582912;
    break;
  case 42:
    size = 13369344;
    break;
  case 50:
    size = 42393600;
    break;
  case 51:
    size = 70778880;
    break;
  default:
    error ("undefined level", 500);
    break;
  }

  size /= pic_size;
  size = imin( size, 16);

  if (active_sps->vui_parameters_present_flag && active_sps->vui_seq_parameters.bitstream_restriction_flag)
  {
    if ((int)active_sps->vui_seq_parameters.max_dec_frame_buffering > size)
    {
      error ("max_dec_frame_buffering larger than MaxDpbSize", 500);
    }
    size = imax (1, active_sps->vui_seq_parameters.max_dec_frame_buffering);
  }

  return size;
}

/*!
 ************************************************************************
 * \brief
 *    Check then number of frames marked "used for reference" and break
 *    if maximum is exceeded
 *
 ************************************************************************
 */
void check_num_ref(DecodedPictureBuffer *p_Dpb)
{
  if ((int)(p_Dpb->ltref_frames_in_buffer +  p_Dpb->ref_frames_in_buffer ) > (imax(1, p_Dpb->num_ref_frames)))
  {
    error ("Max. number of reference frames exceeded. Invalid stream.", 500);
  }
}


/*!
 ************************************************************************
 * \brief
 *    Allocate memory for decoded picture buffer and initialize with sane values.
 *
 ************************************************************************
 */
void init_dpb(ImageParameters *p_Img)
{
  unsigned i,j;
  DecodedPictureBuffer *p_Dpb = p_Img->p_Dpb;
  seq_parameter_set_rbsp_t *active_sps = p_Img->active_sps;

  if (p_Dpb->init_done)
  {
    free_dpb(p_Img);
  }

  p_Dpb->p_Img = p_Img;
  p_Dpb->size  = getDpbSize(active_sps);

  p_Dpb->num_ref_frames = active_sps->num_ref_frames;

  if (p_Dpb->size < active_sps->num_ref_frames)
  {
    error ("DPB size at specified level is smaller than the specified number of reference frames. This is not allowed.\n", 1000);
  }

  p_Dpb->used_size = 0;
  p_Dpb->last_picture = NULL;

  p_Dpb->ref_frames_in_buffer = 0;
  p_Dpb->ltref_frames_in_buffer = 0;

  p_Dpb->fs = calloc(p_Dpb->size, sizeof (FrameStore*));
  if (NULL==p_Dpb->fs)
    no_mem_exit("init_dpb: dpb->fs");

  p_Dpb->fs_ref = calloc(p_Dpb->size, sizeof (FrameStore*));
  if (NULL==p_Dpb->fs_ref)
    no_mem_exit("init_dpb: dpb->fs_ref");

  p_Dpb->fs_ltref = calloc(p_Dpb->size, sizeof (FrameStore*));
  if (NULL==p_Dpb->fs_ltref)
    no_mem_exit("init_dpb: dpb->fs_ltref");

  for (i=0; i<p_Dpb->size; i++)
  {
    p_Dpb->fs[i]       = alloc_frame_store();
    p_Dpb->fs_ref[i]   = NULL;
    p_Dpb->fs_ltref[i] = NULL;
  }

  for (i=0; i<6; i++)
  {
    p_Img->listX[i] = calloc(MAX_LIST_SIZE, sizeof (StorablePicture*)); // +1 for reordering
    if (NULL==p_Img->listX[i])
      no_mem_exit("init_dpb: p_Img->listX[i]");
  }

  /* allocate a dummy storable picture */
  p_Img->no_reference_picture = alloc_storable_picture (p_Img, FRAME, p_Img->width, p_Img->height, p_Img->width_cr, p_Img->height_cr);
  p_Img->no_reference_picture->top_field    = p_Img->no_reference_picture;
  p_Img->no_reference_picture->bottom_field = p_Img->no_reference_picture;
  p_Img->no_reference_picture->frame        = p_Img->no_reference_picture;


  for (j=0;j<6;j++)
  {
    for (i=0; i<MAX_LIST_SIZE; i++)
    {
      p_Img->listX[j][i] = NULL;
    }
    p_Img->listXsize[j]=0;
  }

  p_Dpb->last_output_poc = INT_MIN;

  p_Img->last_has_mmco_5 = 0;

  p_Dpb->init_done = 1;

  // picture error concealment
  if(p_Img->conceal_mode !=0)
      p_Img->last_out_fs = alloc_frame_store();
}
/*!
 ************************************************************************
 * \brief
 *    Free memory for decoded picture buffer.
 ************************************************************************
 */
void free_dpb(ImageParameters *p_Img)
{
  DecodedPictureBuffer *p_Dpb = p_Img->p_Dpb;
  unsigned i;
  if (p_Dpb->fs)
  {
    for (i=0; i<p_Dpb->size; i++)
    {
      free_frame_store(p_Img, p_Dpb->fs[i]);
    }
    free (p_Dpb->fs);
    p_Dpb->fs=NULL;
  }
  if (p_Dpb->fs_ref)
  {
    free (p_Dpb->fs_ref);
  }
  if (p_Dpb->fs_ltref)
  {
    free (p_Dpb->fs_ltref);
  }
  p_Dpb->last_output_poc = INT_MIN;

  for (i=0; i<6; i++)
    if (p_Img->listX[i])
    {
      free (p_Img->listX[i]);
      p_Img->listX[i] = NULL;
    }

  p_Dpb->init_done = 0;

  // picture error concealment
  if(p_Img->conceal_mode != 0)
      free_frame_store(p_Img, p_Img->last_out_fs);

  free_storable_picture(p_Img, p_Img->no_reference_picture);
}


/*!
 ************************************************************************
 * \brief
 *    Allocate memory for decoded picture buffer frame stores an initialize with sane values.
 *
 * \return
 *    the allocated FrameStore structure
 ************************************************************************
 */
FrameStore* alloc_frame_store(void)
{
  FrameStore *f;

  f = calloc (1, sizeof(FrameStore));
  if (NULL==f)
    no_mem_exit("alloc_frame_store: f");

  f->is_used      = 0;
  f->is_reference = 0;
  f->is_long_term = 0;
  f->is_orig_reference = 0;

  f->is_output = 0;

  f->frame        = NULL;;
  f->top_field    = NULL;
  f->bottom_field = NULL;

  return f;
}

void alloc_pic_motion(ImageParameters *p_Img, PicMotionParams *motion, int size_y, int size_x)
{
  seq_parameter_set_rbsp_t *active_sps = p_Img->active_sps;  

  if (active_sps->frame_mbs_only_flag)
  {
    get_mem3Dint64 (&(motion->ref_pic_id), 2, size_y, size_x);
    get_mem3Dint64 (&(motion->ref_id),     2, size_y, size_x);
  }
  else
  {
    get_mem3Dint64 (&(motion->ref_pic_id), 6, size_y, size_x);
    get_mem3Dint64 (&(motion->ref_id),     6, size_y, size_x);
  }

  get_mem4Dshort (&(motion->mv),         2, size_y, size_x, 2);
  get_mem3D      ((byte****)(&(motion->ref_idx)),    2, size_y , size_x);

  motion->mb_field = calloc (size_y * size_x, sizeof(byte));
  if (motion->mb_field == NULL)
    no_mem_exit("alloc_storable_picture: motion->mb_field");

  get_mem2D (&(motion->field_frame), size_y, size_x);
}

/*!
 ************************************************************************
 * \brief
 *    Allocate memory for a stored picture.
 *
 * \param p_Img
 *      image decoding parameters for current picture
 * \param structure
 *    picture structure
 * \param size_x
 *    horizontal luma size
 * \param size_y
 *    vertical luma size
 * \param size_x_cr
 *    horizontal chroma size
 * \param size_y_cr
 *    vertical chroma size
 *
 * \return
 *    the allocated StorablePicture structure
 ************************************************************************
 */
StorablePicture* alloc_storable_picture(ImageParameters *p_Img, PictureStructure structure, int size_x, int size_y, int size_x_cr, int size_y_cr)
{
  seq_parameter_set_rbsp_t *active_sps = p_Img->active_sps;  

  StorablePicture *s;
  int   nplane;

  //printf ("Allocating (%s) picture (x=%d, y=%d, x_cr=%d, y_cr=%d)\n", (type == FRAME)?"FRAME":(type == TOP_FIELD)?"TOP_FIELD":"BOTTOM_FIELD", size_x, size_y, size_x_cr, size_y_cr);

  s = calloc (1, sizeof(StorablePicture));
  if (NULL==s)
    no_mem_exit("alloc_storable_picture: s");

  if (structure!=FRAME)
  {
    size_y    /= 2;
    size_y_cr /= 2;
  }

  s->PicSizeInMbs = (size_x*size_y)/256;
  s->imgUV = NULL;

  get_mem2Dpel (&(s->imgY), size_y, size_x);

  if (active_sps->chroma_format_idc != YUV400)
    get_mem3Dpel (&(s->imgUV), 2, size_y_cr, size_x_cr);
  
  get_mem2Dshort (&(s->slice_id), size_y / MB_BLOCK_SIZE, size_x / MB_BLOCK_SIZE);

  alloc_pic_motion(p_Img, &s->motion, size_y / BLOCK_SIZE, size_x / BLOCK_SIZE);

  if( IS_INDEPENDENT(p_Img) )
  {
    for( nplane=0; nplane<MAX_PLANE; nplane++ )
    {
      alloc_pic_motion(p_Img, &s->JVmotion[nplane], size_y / BLOCK_SIZE, size_x / BLOCK_SIZE);
    }
  }

  s->pic_num=0;
  s->frame_num=0;
  s->long_term_frame_idx=0;
  s->long_term_pic_num=0;
  s->used_for_reference=0;
  s->is_long_term=0;
  s->non_existing=0;
  s->is_output = 0;
  s->max_slice_id = 0;

  s->structure=structure;

  s->size_x = size_x;
  s->size_y = size_y;
  s->size_x_cr = size_x_cr;
  s->size_y_cr = size_y_cr;
  s->size_x_m1 = size_x - 1;
  s->size_y_m1 = size_y - 1;
  s->size_x_cr_m1 = size_x_cr - 1;
  s->size_y_cr_m1 = size_y_cr - 1;

  s->top_field    = p_Img->no_reference_picture;
  s->bottom_field = p_Img->no_reference_picture;
  s->frame        = p_Img->no_reference_picture;

  s->dec_ref_pic_marking_buffer = NULL;

  s->coded_frame                   = 0;
  s->MbaffFrameFlag                = 0;

  s->top_poc = s->bottom_poc = s->poc = 0;
  s->seiHasTone_mapping = 0;

  return s;
}

/*!
 ************************************************************************
 * \brief
 *    Free frame store memory.
 *
 * \param p_Img
 *      image decoding parameters for current picture
 * \param f
 *    FrameStore to be freed
 *
 ************************************************************************
 */
void free_frame_store(ImageParameters *p_Img, FrameStore* f)
{
  if (f)
  {
    if (f->frame)
    {
      free_storable_picture(p_Img, f->frame);
      f->frame=NULL;
    }
    if (f->top_field)
    {
      free_storable_picture(p_Img, f->top_field);
      f->top_field=NULL;
    }
    if (f->bottom_field)
    {
      free_storable_picture(p_Img, f->bottom_field);
      f->bottom_field=NULL;
    }
    free(f);
  }
}

void free_pic_motion(PicMotionParams *motion)
{
  if (motion->ref_pic_id)
  {
    free_mem3Dint64 (motion->ref_pic_id);
    motion->ref_pic_id = NULL;
  }
  if (motion->ref_id)
  {
    free_mem3Dint64 (motion->ref_id);
    motion->ref_id = NULL;
  }

  if (motion->mv)
  {
    free_mem4Dshort (motion->mv);
    motion->mv = NULL;
  }

  if (motion->ref_idx)
  {
    free_mem3D ((byte***)motion->ref_idx);
    motion->ref_idx = NULL;
  }

  if (motion->mb_field)
  {
    free(motion->mb_field);
    motion->mb_field = NULL;
  }

  if (motion->field_frame)
  {
    free_mem2D (motion->field_frame);
    motion->field_frame=NULL;
  }
}


/*!
 ************************************************************************
 * \brief
 *    Free picture memory.
 *
 * \param p_Img
 *      image decoding parameters for current picture
 * \param p
 *    Picture to be freed
 *
 ************************************************************************
 */
void free_storable_picture(ImageParameters *p_Img, StorablePicture* p)
{
  int nplane;
  if (p)
  {
    free_pic_motion(&p->motion);

    if( IS_INDEPENDENT(p_Img) )
    {
      for( nplane=0; nplane<MAX_PLANE; nplane++ )
      {
        free_pic_motion(&p->JVmotion[nplane]);
      }
    }

    if (p->imgY)
    {
      free_mem2Dpel (p->imgY);
      p->imgY=NULL;
    }

    if (p->imgUV)
    {
      free_mem3Dpel (p->imgUV);
      p->imgUV=NULL;
    }

    if (p->slice_id)
    {
      free_mem2Dshort(p->slice_id);
      p->slice_id=NULL;
    }

    if (p->seiHasTone_mapping)
      free(p->tone_mapping_lut);

    free(p);
    p = NULL;
  }
}

/*!
 ************************************************************************
 * \brief
 *    mark FrameStore unused for reference
 *
 ************************************************************************
 */
static void unmark_for_reference(FrameStore* fs)
{

  if (fs->is_used & 1)
  {
    if (fs->top_field)
    {
      fs->top_field->used_for_reference = 0;
    }
  }
  if (fs->is_used & 2)
  {
    if (fs->bottom_field)
    {
      fs->bottom_field->used_for_reference = 0;
    }
  }
  if (fs->is_used == 3)
  {
    if (fs->top_field && fs->bottom_field)
    {
      fs->top_field->used_for_reference = 0;
      fs->bottom_field->used_for_reference = 0;
    }
    fs->frame->used_for_reference = 0;
  }

  fs->is_reference = 0;

  if(fs->frame)
  {
    free_pic_motion(&fs->frame->motion);
  }

  if (fs->top_field)
  {
    free_pic_motion(&fs->top_field->motion);
  }

  if (fs->bottom_field)
  {
    free_pic_motion(&fs->bottom_field->motion);
  }
}


/*!
 ************************************************************************
 * \brief
 *    mark FrameStore unused for reference and reset long term flags
 *
 ************************************************************************
 */
static void unmark_for_long_term_reference(FrameStore* fs)
{

  if (fs->is_used & 1)
  {
    if (fs->top_field)
    {
      fs->top_field->used_for_reference = 0;
      fs->top_field->is_long_term = 0;
    }
  }
  if (fs->is_used & 2)
  {
    if (fs->bottom_field)
    {
      fs->bottom_field->used_for_reference = 0;
      fs->bottom_field->is_long_term = 0;
    }
  }
  if (fs->is_used == 3)
  {
    if (fs->top_field && fs->bottom_field)
    {
      fs->top_field->used_for_reference = 0;
      fs->top_field->is_long_term = 0;
      fs->bottom_field->used_for_reference = 0;
      fs->bottom_field->is_long_term = 0;
    }
    fs->frame->used_for_reference = 0;
    fs->frame->is_long_term = 0;
  }

  fs->is_reference = 0;
  fs->is_long_term = 0;
}


/*!
 ************************************************************************
 * \brief
 *    compares two stored pictures by picture number for qsort in descending order
 *
 ************************************************************************
 */
static inline int compare_pic_by_pic_num_desc( const void *arg1, const void *arg2 )
{
  int pic_num1 = (*(StorablePicture**)arg1)->pic_num;
  int pic_num2 = (*(StorablePicture**)arg2)->pic_num;

  if (pic_num1 < pic_num2)
    return 1;
  if (pic_num1 > pic_num2)
    return -1;
  else
    return 0;
}

/*!
 ************************************************************************
 * \brief
 *    compares two stored pictures by picture number for qsort in descending order
 *
 ************************************************************************
 */
static inline int compare_pic_by_lt_pic_num_asc( const void *arg1, const void *arg2 )
{
  int long_term_pic_num1 = (*(StorablePicture**)arg1)->long_term_pic_num;
  int long_term_pic_num2 = (*(StorablePicture**)arg2)->long_term_pic_num;

  if ( long_term_pic_num1 < long_term_pic_num2)
    return -1;
  if ( long_term_pic_num1 > long_term_pic_num2)
    return 1;
  else
    return 0;
}

/*!
 ************************************************************************
 * \brief
 *    compares two frame stores by pic_num for qsort in descending order
 *
 ************************************************************************
 */
static inline int compare_fs_by_frame_num_desc( const void *arg1, const void *arg2 )
{
  int frame_num_wrap1 = (*(FrameStore**)arg1)->frame_num_wrap;
  int frame_num_wrap2 = (*(FrameStore**)arg2)->frame_num_wrap;
  if ( frame_num_wrap1 < frame_num_wrap2)
    return 1;
  if ( frame_num_wrap1 > frame_num_wrap2)
    return -1;
  else
    return 0;
}


/*!
 ************************************************************************
 * \brief
 *    compares two frame stores by lt_pic_num for qsort in descending order
 *
 ************************************************************************
 */
static inline int compare_fs_by_lt_pic_idx_asc( const void *arg1, const void *arg2 )
{
  int long_term_frame_idx1 = (*(FrameStore**)arg1)->long_term_frame_idx;
  int long_term_frame_idx2 = (*(FrameStore**)arg2)->long_term_frame_idx;

  if ( long_term_frame_idx1 < long_term_frame_idx2)
    return -1;
  if ( long_term_frame_idx1 > long_term_frame_idx2)
    return 1;
  else
    return 0;
}


/*!
 ************************************************************************
 * \brief
 *    compares two stored pictures by poc for qsort in ascending order
 *
 ************************************************************************
 */
static inline int compare_pic_by_poc_asc( const void *arg1, const void *arg2 )
{
  int poc1 = (*(StorablePicture**)arg1)->poc;
  int poc2 = (*(StorablePicture**)arg2)->poc;

  if ( poc1 < poc2)
    return -1;  
  if ( poc1 > poc2)
    return 1;
  else
    return 0;
}


/*!
 ************************************************************************
 * \brief
 *    compares two stored pictures by poc for qsort in descending order
 *
 ************************************************************************
 */
static inline int compare_pic_by_poc_desc( const void *arg1, const void *arg2 )
{
  int poc1 = (*(StorablePicture**)arg1)->poc;
  int poc2 = (*(StorablePicture**)arg2)->poc;

  if (poc1 < poc2)
    return 1;
  if (poc1 > poc2)
    return -1;
  else
    return 0;
}


/*!
 ************************************************************************
 * \brief
 *    compares two frame stores by poc for qsort in ascending order
 *
 ************************************************************************
 */
static inline int compare_fs_by_poc_asc( const void *arg1, const void *arg2 )
{
  int poc1 = (*(FrameStore**)arg1)->poc;
  int poc2 = (*(FrameStore**)arg2)->poc;

  if (poc1 < poc2)
    return -1;  
  if (poc1 > poc2)
    return 1;
  else
    return 0;
}


/*!
 ************************************************************************
 * \brief
 *    compares two frame stores by poc for qsort in descending order
 *
 ************************************************************************
 */
static inline int compare_fs_by_poc_desc( const void *arg1, const void *arg2 )
{
  int poc1 = (*(FrameStore**)arg1)->poc;
  int poc2 = (*(FrameStore**)arg2)->poc;

  if (poc1 < poc2)
    return 1;
  if (poc1 > poc2)
    return -1;
  else
    return 0;
}


/*!
 ************************************************************************
 * \brief
 *    returns true, if picture is short term reference picture
 *
 ************************************************************************
 */
int is_short_ref(StorablePicture *s)
{
  return ((s->used_for_reference) && (!(s->is_long_term)));
}


/*!
 ************************************************************************
 * \brief
 *    returns true, if picture is long term reference picture
 *
 ************************************************************************
 */
int is_long_ref(StorablePicture *s)
{
  return ((s->used_for_reference) && (s->is_long_term));
}


/*!
 ************************************************************************
 * \brief
 *    Generates a alternating field list from a given FrameStore list
 *
 ************************************************************************
 */
static void gen_pic_list_from_frame_list(PictureStructure currStructure, FrameStore **fs_list, int list_idx, StorablePicture **list, char *list_size, int long_term)
{
  int top_idx = 0;
  int bot_idx = 0;

  int (*is_ref)(StorablePicture *s);

  if (long_term)
    is_ref=is_long_ref;
  else
    is_ref=is_short_ref;

  if (currStructure == TOP_FIELD)
  {
    while ((top_idx<list_idx)||(bot_idx<list_idx))
    {
      for ( ; top_idx<list_idx; top_idx++)
      {
        if(fs_list[top_idx]->is_used & 1)
        {
          if(is_ref(fs_list[top_idx]->top_field))
          {
            // short term ref pic
            list[(short) *list_size] = fs_list[top_idx]->top_field;
            (*list_size)++;
            top_idx++;
            break;
          }
        }
      }
      for ( ; bot_idx<list_idx; bot_idx++)
      {
        if(fs_list[bot_idx]->is_used & 2)
        {
          if(is_ref(fs_list[bot_idx]->bottom_field))
          {
            // short term ref pic
            list[(short) *list_size] = fs_list[bot_idx]->bottom_field;
            (*list_size)++;
            bot_idx++;
            break;
          }
        }
      }
    }
  }
  if (currStructure == BOTTOM_FIELD)
  {
    while ((top_idx<list_idx)||(bot_idx<list_idx))
    {
      for ( ; bot_idx<list_idx; bot_idx++)
      {
        if(fs_list[bot_idx]->is_used & 2)
        {
          if(is_ref(fs_list[bot_idx]->bottom_field))
          {
            // short term ref pic
            list[(short) *list_size] = fs_list[bot_idx]->bottom_field;
            (*list_size)++;
            bot_idx++;
            break;
          }
        }
      }
      for ( ; top_idx<list_idx; top_idx++)
      {
        if(fs_list[top_idx]->is_used & 1)
        {
          if(is_ref(fs_list[top_idx]->top_field))
          {
            // short term ref pic
            list[(short) *list_size] = fs_list[top_idx]->top_field;
            (*list_size)++;
            top_idx++;
            break;
          }
        }
      }
    }
  }
}


/*!
 ************************************************************************
 * \brief
 *    Initialize p_Img->listX[0] and list 1 depending on current slice type
 *
 ************************************************************************
 */
void init_lists(Slice *currSlice)
{
  ImageParameters *p_Img = currSlice->p_Img;
  DecodedPictureBuffer *p_Dpb = p_Img->p_Dpb;
  seq_parameter_set_rbsp_t *active_sps = p_Img->active_sps;

  int add_top = 0, add_bottom = 0;
  unsigned i;
  int j;
  int MaxFrameNum = 1 << (active_sps->log2_max_frame_num_minus4 + 4);
  int diff;

  int list0idx = 0;
  int list0idx_1 = 0;
  int listltidx = 0;

  FrameStore **fs_list0;
  FrameStore **fs_list1;
  FrameStore **fs_listlt;

  StorablePicture *tmp_s;

  if (currSlice->structure == FRAME)
  {
    for (i=0; i<p_Dpb->ref_frames_in_buffer; i++)
    {
      if (p_Dpb->fs_ref[i]->is_used==3)
      {
        if ((p_Dpb->fs_ref[i]->frame->used_for_reference)&&(!p_Dpb->fs_ref[i]->frame->is_long_term))
        {
          if( p_Dpb->fs_ref[i]->frame_num > p_Img->frame_num )
          {
            p_Dpb->fs_ref[i]->frame_num_wrap = p_Dpb->fs_ref[i]->frame_num - MaxFrameNum;
          }
          else
          {
            p_Dpb->fs_ref[i]->frame_num_wrap = p_Dpb->fs_ref[i]->frame_num;
          }
          p_Dpb->fs_ref[i]->frame->pic_num = p_Dpb->fs_ref[i]->frame_num_wrap;
        }
      }
    }
    // update long_term_pic_num
    for (i = 0; i < p_Dpb->ltref_frames_in_buffer; i++)
    {
      if (p_Dpb->fs_ltref[i]->is_used==3)
      {
        if (p_Dpb->fs_ltref[i]->frame->is_long_term)
        {
          p_Dpb->fs_ltref[i]->frame->long_term_pic_num = p_Dpb->fs_ltref[i]->frame->long_term_frame_idx;
        }
      }
    }
  }
  else
  {
    if (currSlice->structure == TOP_FIELD)
    {
      add_top    = 1;
      add_bottom = 0;
    }
    else
    {
      add_top    = 0;
      add_bottom = 1;
    }

    for (i=0; i<p_Dpb->ref_frames_in_buffer; i++)
    {
      if (p_Dpb->fs_ref[i]->is_reference)
      {
        if( p_Dpb->fs_ref[i]->frame_num > p_Img->frame_num )
        {
          p_Dpb->fs_ref[i]->frame_num_wrap = p_Dpb->fs_ref[i]->frame_num - MaxFrameNum;
        }
        else
        {
          p_Dpb->fs_ref[i]->frame_num_wrap = p_Dpb->fs_ref[i]->frame_num;
        }
        if (p_Dpb->fs_ref[i]->is_reference & 1)
        {
          p_Dpb->fs_ref[i]->top_field->pic_num = (2 * p_Dpb->fs_ref[i]->frame_num_wrap) + add_top;
        }
        if (p_Dpb->fs_ref[i]->is_reference & 2)
        {
          p_Dpb->fs_ref[i]->bottom_field->pic_num = (2 * p_Dpb->fs_ref[i]->frame_num_wrap) + add_bottom;
        }
      }
    }
    // update long_term_pic_num
    for (i=0; i<p_Dpb->ltref_frames_in_buffer; i++)
    {
      if (p_Dpb->fs_ltref[i]->is_long_term & 1)
      {
        p_Dpb->fs_ltref[i]->top_field->long_term_pic_num = 2 * p_Dpb->fs_ltref[i]->top_field->long_term_frame_idx + add_top;
      }
      if (p_Dpb->fs_ltref[i]->is_long_term & 2)
      {
        p_Dpb->fs_ltref[i]->bottom_field->long_term_pic_num = 2 * p_Dpb->fs_ltref[i]->bottom_field->long_term_frame_idx + add_bottom;
      }
    }
  }

  if ((currSlice->slice_type == I_SLICE)||(currSlice->slice_type == SI_SLICE))
  {
    p_Img->listXsize[0] = 0;
    p_Img->listXsize[1] = 0;
    return;
  }

  if ((currSlice->slice_type == P_SLICE)||(currSlice->slice_type == SP_SLICE))
  {
    // Calculate FrameNumWrap and PicNum
    if (currSlice->structure == FRAME)
    {
      for (i=0; i<p_Dpb->ref_frames_in_buffer; i++)
      {
        if (p_Dpb->fs_ref[i]->is_used==3)
        {
          if ((p_Dpb->fs_ref[i]->frame->used_for_reference)&&(!p_Dpb->fs_ref[i]->frame->is_long_term))
          {
            p_Img->listX[0][list0idx++] = p_Dpb->fs_ref[i]->frame;
          }
        }
      }
      // order list 0 by PicNum
      qsort((void *)p_Img->listX[0], list0idx, sizeof(StorablePicture*), compare_pic_by_pic_num_desc);
      p_Img->listXsize[0] = list0idx;
//      printf("listX[0] (PicNum): "); for (i=0; i<list0idx; i++){printf ("%d  ", p_Img->listX[0][i]->pic_num);} printf("\n");

      // long term handling
      for (i=0; i<p_Dpb->ltref_frames_in_buffer; i++)
      {
        if (p_Dpb->fs_ltref[i]->is_used==3)
        {
          if (p_Dpb->fs_ltref[i]->frame->is_long_term)
          {
            p_Img->listX[0][list0idx++]=p_Dpb->fs_ltref[i]->frame;
          }
        }
      }
      qsort((void *)&p_Img->listX[0][(short) p_Img->listXsize[0]], list0idx - p_Img->listXsize[0], sizeof(StorablePicture*), compare_pic_by_lt_pic_num_asc);
      p_Img->listXsize[0] = list0idx;
    }
    else
    {
      fs_list0 = calloc(p_Dpb->size, sizeof (FrameStore*));
      if (NULL==fs_list0)
         no_mem_exit("init_lists: fs_list0");
      fs_listlt = calloc(p_Dpb->size, sizeof (FrameStore*));
      if (NULL==fs_listlt)
         no_mem_exit("init_lists: fs_listlt");

      for (i=0; i<p_Dpb->ref_frames_in_buffer; i++)
      {
        if (p_Dpb->fs_ref[i]->is_reference)
        {
          fs_list0[list0idx++] = p_Dpb->fs_ref[i];
        }
      }

      qsort((void *)fs_list0, list0idx, sizeof(FrameStore*), compare_fs_by_frame_num_desc);

//      printf("fs_list0 (FrameNum): "); for (i=0; i<list0idx; i++){printf ("%d  ", fs_list0[i]->frame_num_wrap);} printf("\n");

      p_Img->listXsize[0] = 0;
      gen_pic_list_from_frame_list(currSlice->structure, fs_list0, list0idx, p_Img->listX[0], &p_Img->listXsize[0], 0);

//      printf("p_Img->listX[0] (PicNum): "); for (i=0; i<p_Img->listXsize[0]; i++){printf ("%d  ", p_Img->listX[0][i]->pic_num);} printf("\n");

      // long term handling
      for (i=0; i<p_Dpb->ltref_frames_in_buffer; i++)
      {
        fs_listlt[listltidx++]=p_Dpb->fs_ltref[i];
      }

      qsort((void *)fs_listlt, listltidx, sizeof(FrameStore*), compare_fs_by_lt_pic_idx_asc);

      gen_pic_list_from_frame_list(currSlice->structure, fs_listlt, listltidx, p_Img->listX[0], &p_Img->listXsize[0], 1);

      free(fs_list0);
      free(fs_listlt);
    }
    p_Img->listXsize[1] = 0;
  }
  else
  {
    // B-Slice
    if (currSlice->structure == FRAME)
    {
      for (i=0; i<p_Dpb->ref_frames_in_buffer; i++)
      {
        if (p_Dpb->fs_ref[i]->is_used==3)
        {
          if ((p_Dpb->fs_ref[i]->frame->used_for_reference)&&(!p_Dpb->fs_ref[i]->frame->is_long_term))
          {
            if (p_Img->framepoc >= p_Dpb->fs_ref[i]->frame->poc) //!KS use >= for error concealment
//            if (p_Img->framepoc > p_Dpb->fs_ref[i]->frame->poc)
            {
              p_Img->listX[0][list0idx++] = p_Dpb->fs_ref[i]->frame;
            }
          }
        }
      }
      qsort((void *)p_Img->listX[0], list0idx, sizeof(StorablePicture*), compare_pic_by_poc_desc);
      list0idx_1 = list0idx;
      for (i=0; i<p_Dpb->ref_frames_in_buffer; i++)
      {
        if (p_Dpb->fs_ref[i]->is_used==3)
        {
          if ((p_Dpb->fs_ref[i]->frame->used_for_reference)&&(!p_Dpb->fs_ref[i]->frame->is_long_term))
          {
            if (p_Img->framepoc < p_Dpb->fs_ref[i]->frame->poc)
            {
              p_Img->listX[0][list0idx++] = p_Dpb->fs_ref[i]->frame;
            }
          }
        }
      }
      qsort((void *)&p_Img->listX[0][list0idx_1], list0idx-list0idx_1, sizeof(StorablePicture*), compare_pic_by_poc_asc);

      for (j=0; j<list0idx_1; j++)
      {
        p_Img->listX[1][list0idx-list0idx_1+j]=p_Img->listX[0][j];
      }
      for (j=list0idx_1; j<list0idx; j++)
      {
        p_Img->listX[1][j-list0idx_1]=p_Img->listX[0][j];
      }

      p_Img->listXsize[0] = p_Img->listXsize[1] = list0idx;

//      printf("p_Img->listX[0] currPoc=%d (Poc): ", p_Img->framepoc); for (i=0; i<p_Img->listXsize[0]; i++){printf ("%d  ", p_Img->listX[0][i]->poc);} printf("\n");
//      printf("p_Img->listX[1] currPoc=%d (Poc): ", p_Img->framepoc); for (i=0; i<p_Img->listXsize[1]; i++){printf ("%d  ", p_Img->listX[1][i]->poc);} printf("\n");

      // long term handling
      for (i=0; i<p_Dpb->ltref_frames_in_buffer; i++)
      {
        if (p_Dpb->fs_ltref[i]->is_used==3)
        {
          if (p_Dpb->fs_ltref[i]->frame->is_long_term)
          {
            p_Img->listX[0][list0idx]  =p_Dpb->fs_ltref[i]->frame;
            p_Img->listX[1][list0idx++]=p_Dpb->fs_ltref[i]->frame;
          }
        }
      }
      qsort((void *)&p_Img->listX[0][(short) p_Img->listXsize[0]], list0idx-p_Img->listXsize[0], sizeof(StorablePicture*), compare_pic_by_lt_pic_num_asc);
      qsort((void *)&p_Img->listX[1][(short) p_Img->listXsize[0]], list0idx-p_Img->listXsize[0], sizeof(StorablePicture*), compare_pic_by_lt_pic_num_asc);
      p_Img->listXsize[0] = p_Img->listXsize[1] = list0idx;
    }
    else
    {
      fs_list0 = calloc(p_Dpb->size, sizeof (FrameStore*));
      if (NULL==fs_list0)
         no_mem_exit("init_lists: fs_list0");
      fs_list1 = calloc(p_Dpb->size, sizeof (FrameStore*));
      if (NULL==fs_list1)
         no_mem_exit("init_lists: fs_list1");
      fs_listlt = calloc(p_Dpb->size, sizeof (FrameStore*));
      if (NULL==fs_listlt)
         no_mem_exit("init_lists: fs_listlt");

      p_Img->listXsize[0] = 0;
      p_Img->listXsize[1] = 1;

      for (i=0; i<p_Dpb->ref_frames_in_buffer; i++)
      {
        if (p_Dpb->fs_ref[i]->is_used)
        {
          if (p_Img->ThisPOC >= p_Dpb->fs_ref[i]->poc)
          {
            fs_list0[list0idx++] = p_Dpb->fs_ref[i];
          }
        }
      }
      qsort((void *)fs_list0, list0idx, sizeof(FrameStore*), compare_fs_by_poc_desc);
      list0idx_1 = list0idx;
      for (i=0; i<p_Dpb->ref_frames_in_buffer; i++)
      {
        if (p_Dpb->fs_ref[i]->is_used)
        {
          if (p_Img->ThisPOC < p_Dpb->fs_ref[i]->poc)
          {
            fs_list0[list0idx++] = p_Dpb->fs_ref[i];
          }
        }
      }
      qsort((void *)&fs_list0[list0idx_1], list0idx-list0idx_1, sizeof(FrameStore*), compare_fs_by_poc_asc);

      for (j=0; j<list0idx_1; j++)
      {
        fs_list1[list0idx-list0idx_1+j]=fs_list0[j];
      }
      for (j=list0idx_1; j<list0idx; j++)
      {
        fs_list1[j-list0idx_1]=fs_list0[j];
      }

//      printf("fs_list0 currPoc=%d (Poc): ", p_Img->ThisPOC); for (i=0; i<list0idx; i++){printf ("%d  ", fs_list0[i]->poc);} printf("\n");
//      printf("fs_list1 currPoc=%d (Poc): ", p_Img->ThisPOC); for (i=0; i<list0idx; i++){printf ("%d  ", fs_list1[i]->poc);} printf("\n");

      p_Img->listXsize[0] = 0;
      p_Img->listXsize[1] = 0;
      gen_pic_list_from_frame_list(currSlice->structure, fs_list0, list0idx, p_Img->listX[0], &p_Img->listXsize[0], 0);
      gen_pic_list_from_frame_list(currSlice->structure, fs_list1, list0idx, p_Img->listX[1], &p_Img->listXsize[1], 0);

//      printf("p_Img->listX[0] currPoc=%d (Poc): ", p_Img->framepoc); for (i=0; i<p_Img->listXsize[0]; i++){printf ("%d  ", p_Img->listX[0][i]->poc);} printf("\n");
//      printf("p_Img->listX[1] currPoc=%d (Poc): ", p_Img->framepoc); for (i=0; i<p_Img->listXsize[1]; i++){printf ("%d  ", p_Img->listX[1][i]->poc);} printf("\n");

      // long term handling
      for (i=0; i<p_Dpb->ltref_frames_in_buffer; i++)
      {
        fs_listlt[listltidx++]=p_Dpb->fs_ltref[i];
      }

      qsort((void *)fs_listlt, listltidx, sizeof(FrameStore*), compare_fs_by_lt_pic_idx_asc);

      gen_pic_list_from_frame_list(currSlice->structure, fs_listlt, listltidx, p_Img->listX[0], &p_Img->listXsize[0], 1);
      gen_pic_list_from_frame_list(currSlice->structure, fs_listlt, listltidx, p_Img->listX[1], &p_Img->listXsize[1], 1);

      free(fs_list0);
      free(fs_list1);
      free(fs_listlt);
    }
  }

  if ((p_Img->listXsize[0] == p_Img->listXsize[1]) && (p_Img->listXsize[0] > 1))
  {
    // check if lists are identical, if yes swap first two elements of p_Img->listX[1]
    diff=0;
    for (j = 0; j< p_Img->listXsize[0]; j++)
    {
      if (p_Img->listX[0][j]!=p_Img->listX[1][j])
        diff=1;
    }
    if (!diff)
    {
      tmp_s = p_Img->listX[1][0];
      p_Img->listX[1][0]=p_Img->listX[1][1];
      p_Img->listX[1][1]=tmp_s;
    }
  }
  // set max size
  p_Img->listXsize[0] = imin (p_Img->listXsize[0], currSlice->num_ref_idx_l0_active);
  p_Img->listXsize[1] = imin (p_Img->listXsize[1], currSlice->num_ref_idx_l1_active);

  // set the unused list entries to NULL
  for (i=p_Img->listXsize[0]; i< (MAX_LIST_SIZE) ; i++)
  {
      p_Img->listX[0][i] = p_Img->no_reference_picture;

  }
  for (i=p_Img->listXsize[1]; i< (MAX_LIST_SIZE) ; i++)
  {
      p_Img->listX[1][i] = p_Img->no_reference_picture;
  }
}

/*!
 ************************************************************************
 * \brief
 *    Initialize listX[2..5] from lists 0 and 1
 *    listX[2]: list0 for current_field==top
 *    listX[3]: list1 for current_field==top
 *    listX[4]: list0 for current_field==bottom
 *    listX[5]: list1 for current_field==bottom
 *
 ************************************************************************
 */
void init_mbaff_lists(ImageParameters *p_Img)
{
  unsigned j;
  int i;

  for (i=2;i<6;i++)
  {
    for (j=0; j<MAX_LIST_SIZE; j++)
    {
      p_Img->listX[i][j] = p_Img->no_reference_picture;
    }
    p_Img->listXsize[i]=0;
  }

  for (i=0; i<p_Img->listXsize[0]; i++)
  {
    p_Img->listX[2][2*i  ] = p_Img->listX[0][i]->top_field;
    p_Img->listX[2][2*i+1] = p_Img->listX[0][i]->bottom_field;
    p_Img->listX[4][2*i  ] = p_Img->listX[0][i]->bottom_field;
    p_Img->listX[4][2*i+1] = p_Img->listX[0][i]->top_field;
  }
  p_Img->listXsize[2]=p_Img->listXsize[4]=p_Img->listXsize[0] * 2;

  for (i=0; i<p_Img->listXsize[1]; i++)
  {
    p_Img->listX[3][2*i  ] = p_Img->listX[1][i]->top_field;
    p_Img->listX[3][2*i+1] = p_Img->listX[1][i]->bottom_field;
    p_Img->listX[5][2*i  ] = p_Img->listX[1][i]->bottom_field;
    p_Img->listX[5][2*i+1] = p_Img->listX[1][i]->top_field;
  }
  p_Img->listXsize[3]=p_Img->listXsize[5]=p_Img->listXsize[1] * 2;
}

 /*!
 ************************************************************************
 * \brief
 *    Returns short term pic with given picNum
 *
 ************************************************************************
 */
static StorablePicture*  get_short_term_pic(ImageParameters *p_Img, int picNum)
{
  DecodedPictureBuffer *p_Dpb = p_Img->p_Dpb;
  unsigned i;

  for (i=0; i<p_Dpb->ref_frames_in_buffer; i++)
  {
    if (p_Img->structure==FRAME)
    {
      if (p_Dpb->fs_ref[i]->is_reference == 3)
        if ((!p_Dpb->fs_ref[i]->frame->is_long_term)&&(p_Dpb->fs_ref[i]->frame->pic_num == picNum))
          return p_Dpb->fs_ref[i]->frame;
    }
    else
    {
      if (p_Dpb->fs_ref[i]->is_reference & 1)
        if ((!p_Dpb->fs_ref[i]->top_field->is_long_term)&&(p_Dpb->fs_ref[i]->top_field->pic_num == picNum))
          return p_Dpb->fs_ref[i]->top_field;
      if (p_Dpb->fs_ref[i]->is_reference & 2)
        if ((!p_Dpb->fs_ref[i]->bottom_field->is_long_term)&&(p_Dpb->fs_ref[i]->bottom_field->pic_num == picNum))
          return p_Dpb->fs_ref[i]->bottom_field;
    }
  }

  return p_Img->no_reference_picture;
}

/*!
 ************************************************************************
 * \brief
 *    Returns short term pic with given LongtermPicNum
 *
 ************************************************************************
 */
static StorablePicture*  get_long_term_pic(ImageParameters *p_Img, int LongtermPicNum)
{
  DecodedPictureBuffer *p_Dpb = p_Img->p_Dpb;
  unsigned i;

  for (i=0; i<p_Dpb->ltref_frames_in_buffer; i++)
  {
    if (p_Img->structure==FRAME)
    {
      if (p_Dpb->fs_ltref[i]->is_reference == 3)
        if ((p_Dpb->fs_ltref[i]->frame->is_long_term)&&(p_Dpb->fs_ltref[i]->frame->long_term_pic_num == LongtermPicNum))
          return p_Dpb->fs_ltref[i]->frame;
    }
    else
    {
      if (p_Dpb->fs_ltref[i]->is_reference & 1)
        if ((p_Dpb->fs_ltref[i]->top_field->is_long_term)&&(p_Dpb->fs_ltref[i]->top_field->long_term_pic_num == LongtermPicNum))
          return p_Dpb->fs_ltref[i]->top_field;
      if (p_Dpb->fs_ltref[i]->is_reference & 2)
        if ((p_Dpb->fs_ltref[i]->bottom_field->is_long_term)&&(p_Dpb->fs_ltref[i]->bottom_field->long_term_pic_num == LongtermPicNum))
          return p_Dpb->fs_ltref[i]->bottom_field;
    }
  }
  return NULL;
}

/*!
 ************************************************************************
 * \brief
 *    Reordering process for short-term reference pictures
 *
 ************************************************************************
 */
static void reorder_short_term(ImageParameters *p_Img, StorablePicture **RefPicListX, int num_ref_idx_lX_active_minus1, int picNumLX, int *refIdxLX)
{
  int cIdx, nIdx;

  StorablePicture *picLX;

  picLX = get_short_term_pic(p_Img, picNumLX);

  for( cIdx = num_ref_idx_lX_active_minus1+1; cIdx > *refIdxLX; cIdx-- )
    RefPicListX[ cIdx ] = RefPicListX[ cIdx - 1];

  RefPicListX[ (*refIdxLX)++ ] = picLX;

  nIdx = *refIdxLX;

  for( cIdx = *refIdxLX; cIdx <= num_ref_idx_lX_active_minus1+1; cIdx++ )
    if (RefPicListX[ cIdx ])
      if( (RefPicListX[ cIdx ]->is_long_term ) ||  (RefPicListX[ cIdx ]->pic_num != picNumLX ))
        RefPicListX[ nIdx++ ] = RefPicListX[ cIdx ];

}


/*!
 ************************************************************************
 * \brief
 *    Reordering process for long-term reference pictures
 *
 ************************************************************************
 */
static void reorder_long_term(ImageParameters *p_Img, StorablePicture **RefPicListX, int num_ref_idx_lX_active_minus1, int LongTermPicNum, int *refIdxLX)
{
  int cIdx, nIdx;

  StorablePicture *picLX;

  picLX = get_long_term_pic(p_Img, LongTermPicNum);

  for( cIdx = num_ref_idx_lX_active_minus1+1; cIdx > *refIdxLX; cIdx-- )
    RefPicListX[ cIdx ] = RefPicListX[ cIdx - 1];

  RefPicListX[ (*refIdxLX)++ ] = picLX;

  nIdx = *refIdxLX;

  for( cIdx = *refIdxLX; cIdx <= num_ref_idx_lX_active_minus1+1; cIdx++ )
    if (RefPicListX[ cIdx ])
      if( (!RefPicListX[ cIdx ]->is_long_term ) ||  (RefPicListX[ cIdx ]->long_term_pic_num != LongTermPicNum ))
        RefPicListX[ nIdx++ ] = RefPicListX[ cIdx ];
}


/*!
 ************************************************************************
 * \brief
 *    Reordering process for reference picture lists
 *
 ************************************************************************
 */
void reorder_ref_pic_list(ImageParameters *p_Img, StorablePicture **list, char *list_size, int num_ref_idx_lX_active_minus1, int *reordering_of_pic_nums_idc, int *abs_diff_pic_num_minus1, int *long_term_pic_idx)
{
  int i;

  int maxPicNum, currPicNum, picNumLXNoWrap, picNumLXPred, picNumLX;
  int refIdxLX = 0;

  if (p_Img->structure==FRAME)
  {
    maxPicNum  = p_Img->MaxFrameNum;
    currPicNum = p_Img->frame_num;
  }
  else
  {
    maxPicNum  = 2 * p_Img->MaxFrameNum;
    currPicNum = 2 * p_Img->frame_num + 1;
  }

  picNumLXPred = currPicNum;

  for (i=0; reordering_of_pic_nums_idc[i]!=3; i++)
  {
    if (reordering_of_pic_nums_idc[i]>3)
      error ("Invalid remapping_of_pic_nums_idc command", 500);

    if (reordering_of_pic_nums_idc[i] < 2)
    {
      if (reordering_of_pic_nums_idc[i] == 0)
      {
        if( picNumLXPred - ( abs_diff_pic_num_minus1[i] + 1 ) < 0 )
          picNumLXNoWrap = picNumLXPred - ( abs_diff_pic_num_minus1[i] + 1 ) + maxPicNum;
        else
          picNumLXNoWrap = picNumLXPred - ( abs_diff_pic_num_minus1[i] + 1 );
      }
      else // (remapping_of_pic_nums_idc[i] == 1)
      {
        if( picNumLXPred + ( abs_diff_pic_num_minus1[i] + 1 )  >=  maxPicNum )
          picNumLXNoWrap = picNumLXPred + ( abs_diff_pic_num_minus1[i] + 1 ) - maxPicNum;
        else
          picNumLXNoWrap = picNumLXPred + ( abs_diff_pic_num_minus1[i] + 1 );
      }
      picNumLXPred = picNumLXNoWrap;

      if( picNumLXNoWrap > currPicNum )
        picNumLX = picNumLXNoWrap - maxPicNum;
      else
        picNumLX = picNumLXNoWrap;

      reorder_short_term(p_Img, list, num_ref_idx_lX_active_minus1, picNumLX, &refIdxLX);
    }
    else //(remapping_of_pic_nums_idc[i] == 2)
    {
      reorder_long_term(p_Img, list, num_ref_idx_lX_active_minus1, long_term_pic_idx[i], &refIdxLX);
    }

  }
  // that's a definition
  *list_size = num_ref_idx_lX_active_minus1 + 1;
}



/*!
 ************************************************************************
 * \brief
 *    Update the list of frame stores that contain reference frames/fields
 *
 ************************************************************************
 */
void update_ref_list(DecodedPictureBuffer *p_Dpb)
{
  unsigned i, j;
  for (i=0, j=0; i<p_Dpb->used_size; i++)
  {
    if (is_short_term_reference(p_Dpb->fs[i]))
    {
      p_Dpb->fs_ref[j++]=p_Dpb->fs[i];
    }
  }

  p_Dpb->ref_frames_in_buffer = j;

  while (j<p_Dpb->size)
  {
    p_Dpb->fs_ref[j++]=NULL;
  }
}


/*!
 ************************************************************************
 * \brief
 *    Update the list of frame stores that contain long-term reference
 *    frames/fields
 *
 ************************************************************************
 */
void update_ltref_list(DecodedPictureBuffer *p_Dpb)
{
  unsigned i, j;
  for (i=0, j=0; i<p_Dpb->used_size; i++)
  {
    if (is_long_term_reference(p_Dpb->fs[i]))
    {
      p_Dpb->fs_ltref[j++]=p_Dpb->fs[i];
    }
  }

  p_Dpb->ltref_frames_in_buffer=j;

  while (j<p_Dpb->size)
  {
    p_Dpb->fs_ltref[j++]=NULL;
  }
}

/*!
 ************************************************************************
 * \brief
 *    Perform Memory management for idr pictures
 *
 ************************************************************************
 */
static void idr_memory_management(ImageParameters *p_Img, StorablePicture* p)
{
  DecodedPictureBuffer *p_Dpb = p_Img->p_Dpb;
  unsigned i;

  assert (p->idr_flag);

  if (p->no_output_of_prior_pics_flag)
  {
    // free all stored pictures
    for (i=0; i<p_Dpb->used_size; i++)
    {
      // reset all reference settings
      free_frame_store(p_Img, p_Dpb->fs[i]);
      p_Dpb->fs[i] = alloc_frame_store();
    }
    for (i=0; i<p_Dpb->ref_frames_in_buffer; i++)
    {
      p_Dpb->fs_ref[i]=NULL;
    }
    for (i=0; i<p_Dpb->ltref_frames_in_buffer; i++)
    {
      p_Dpb->fs_ltref[i]=NULL;
    }
    p_Dpb->used_size=0;
  }
  else
  {
    flush_dpb(p_Img);
  }
  p_Dpb->last_picture = NULL;

  update_ref_list(p_Dpb);
  update_ltref_list(p_Dpb);
  p_Dpb->last_output_poc = INT_MIN;

  if (p->long_term_reference_flag)
  {
    p_Dpb->max_long_term_pic_idx = 0;
    p->is_long_term           = 1;
    p->long_term_frame_idx    = 0;
  }
  else
  {
    p_Dpb->max_long_term_pic_idx = -1;
    p->is_long_term           = 0;
  }
}

/*!
 ************************************************************************
 * \brief
 *    Perform Sliding window decoded reference picture marking process
 *
 ************************************************************************
 */
static void sliding_window_memory_management(DecodedPictureBuffer *p_Dpb, StorablePicture* p)
{
  unsigned i;

  assert (!p->idr_flag);
  // if this is a reference pic with sliding sliding window, unmark first ref frame
  if (p_Dpb->ref_frames_in_buffer==p_Dpb->num_ref_frames - p_Dpb->ltref_frames_in_buffer)
  {
    for (i=0; i<p_Dpb->used_size;i++)
    {
      if (p_Dpb->fs[i]->is_reference  && (!(p_Dpb->fs[i]->is_long_term)))
      {
        unmark_for_reference(p_Dpb->fs[i]);
        update_ref_list(p_Dpb);
        break;
      }
    }
  }

  p->is_long_term = 0;
}

/*!
 ************************************************************************
 * \brief
 *    Calculate picNumX
 ************************************************************************
 */
static int get_pic_num_x (StorablePicture *p, int difference_of_pic_nums_minus1)
{
  int currPicNum;

  if (p->structure == FRAME)
    currPicNum = p->frame_num;
  else
    currPicNum = 2 * p->frame_num + 1;

  return currPicNum - (difference_of_pic_nums_minus1 + 1);
}


/*!
 ************************************************************************
 * \brief
 *    Adaptive Memory Management: Mark short term picture unused
 ************************************************************************
 */
static void mm_unmark_short_term_for_reference(DecodedPictureBuffer *p_Dpb, StorablePicture *p, int difference_of_pic_nums_minus1)
{
  int picNumX;

  unsigned i;

  picNumX = get_pic_num_x(p, difference_of_pic_nums_minus1);

  for (i=0; i<p_Dpb->ref_frames_in_buffer; i++)
  {
    if (p->structure == FRAME)
    {
      if ((p_Dpb->fs_ref[i]->is_reference==3) && (p_Dpb->fs_ref[i]->is_long_term==0))
      {
        if (p_Dpb->fs_ref[i]->frame->pic_num == picNumX)
        {
          unmark_for_reference(p_Dpb->fs_ref[i]);
          return;
        }
      }
    }
    else
    {
      if ((p_Dpb->fs_ref[i]->is_reference & 1) && (!(p_Dpb->fs_ref[i]->is_long_term & 1)))
      {
        if (p_Dpb->fs_ref[i]->top_field->pic_num == picNumX)
        {
          p_Dpb->fs_ref[i]->top_field->used_for_reference = 0;
          p_Dpb->fs_ref[i]->is_reference &= 2;
          if (p_Dpb->fs_ref[i]->is_used == 3)
          {
            p_Dpb->fs_ref[i]->frame->used_for_reference = 0;
          }
          return;
        }
      }
      if ((p_Dpb->fs_ref[i]->is_reference & 2) && (!(p_Dpb->fs_ref[i]->is_long_term & 2)))
      {
        if (p_Dpb->fs_ref[i]->bottom_field->pic_num == picNumX)
        {
          p_Dpb->fs_ref[i]->bottom_field->used_for_reference = 0;
          p_Dpb->fs_ref[i]->is_reference &= 1;
          if (p_Dpb->fs_ref[i]->is_used == 3)
          {
            p_Dpb->fs_ref[i]->frame->used_for_reference = 0;
          }
          return;
        }
      }
    }
  }
}


/*!
 ************************************************************************
 * \brief
 *    Adaptive Memory Management: Mark long term picture unused
 ************************************************************************
 */
static void mm_unmark_long_term_for_reference(DecodedPictureBuffer *p_Dpb, StorablePicture *p, int long_term_pic_num)
{
  unsigned i;
  for (i=0; i<p_Dpb->ltref_frames_in_buffer; i++)
  {
    if (p->structure == FRAME)
    {
      if ((p_Dpb->fs_ltref[i]->is_reference==3) && (p_Dpb->fs_ltref[i]->is_long_term==3))
      {
        if (p_Dpb->fs_ltref[i]->frame->long_term_pic_num == long_term_pic_num)
        {
          unmark_for_long_term_reference(p_Dpb->fs_ltref[i]);
        }
      }
    }
    else
    {
      if ((p_Dpb->fs_ltref[i]->is_reference & 1) && ((p_Dpb->fs_ltref[i]->is_long_term & 1)))
      {
        if (p_Dpb->fs_ltref[i]->top_field->long_term_pic_num == long_term_pic_num)
        {
          p_Dpb->fs_ltref[i]->top_field->used_for_reference = 0;
          p_Dpb->fs_ltref[i]->top_field->is_long_term = 0;
          p_Dpb->fs_ltref[i]->is_reference &= 2;
          p_Dpb->fs_ltref[i]->is_long_term &= 2;
          if (p_Dpb->fs_ltref[i]->is_used == 3)
          {
            p_Dpb->fs_ltref[i]->frame->used_for_reference = 0;
            p_Dpb->fs_ltref[i]->frame->is_long_term = 0;
          }
          return;
        }
      }
      if ((p_Dpb->fs_ltref[i]->is_reference & 2) && ((p_Dpb->fs_ltref[i]->is_long_term & 2)))
      {
        if (p_Dpb->fs_ltref[i]->bottom_field->long_term_pic_num == long_term_pic_num)
        {
          p_Dpb->fs_ltref[i]->bottom_field->used_for_reference = 0;
          p_Dpb->fs_ltref[i]->bottom_field->is_long_term = 0;
          p_Dpb->fs_ltref[i]->is_reference &= 1;
          p_Dpb->fs_ltref[i]->is_long_term &= 1;
          if (p_Dpb->fs_ltref[i]->is_used == 3)
          {
            p_Dpb->fs_ltref[i]->frame->used_for_reference = 0;
            p_Dpb->fs_ltref[i]->frame->is_long_term = 0;
          }
          return;
        }
      }
    }
  }
}


/*!
 ************************************************************************
 * \brief
 *    Mark a long-term reference frame or complementary field pair unused for referemce
 ************************************************************************
 */
static void unmark_long_term_frame_for_reference_by_frame_idx(DecodedPictureBuffer *p_Dpb, int long_term_frame_idx)
{
  unsigned i;
  for(i=0; i<p_Dpb->ltref_frames_in_buffer; i++)
  {
    if (p_Dpb->fs_ltref[i]->long_term_frame_idx == long_term_frame_idx)
      unmark_for_long_term_reference(p_Dpb->fs_ltref[i]);
  }
}

/*!
 ************************************************************************
 * \brief
 *    Mark a long-term reference field unused for reference only if it's not
 *    the complementary field of the picture indicated by picNumX
 ************************************************************************
 */
static void unmark_long_term_field_for_reference_by_frame_idx(ImageParameters *p_Img, PictureStructure structure, int long_term_frame_idx, int mark_current, unsigned curr_frame_num, int curr_pic_num)
{
  unsigned i;
  DecodedPictureBuffer *p_Dpb = p_Img->p_Dpb;

  assert(structure!=FRAME);
  if (curr_pic_num<0)
    curr_pic_num+=(2*p_Img->MaxFrameNum);

  for(i=0; i<p_Dpb->ltref_frames_in_buffer; i++)
  {
    if (p_Dpb->fs_ltref[i]->long_term_frame_idx == long_term_frame_idx)
    {
      if (structure == TOP_FIELD)
      {
        if ((p_Dpb->fs_ltref[i]->is_long_term == 3))
        {
          unmark_for_long_term_reference(p_Dpb->fs_ltref[i]);
        }
        else
        {
          if ((p_Dpb->fs_ltref[i]->is_long_term == 1))
          {
            unmark_for_long_term_reference(p_Dpb->fs_ltref[i]);
          }
          else
          {
            if (mark_current)
            {
              if (p_Dpb->last_picture)
              {
                if ( ( p_Dpb->last_picture != p_Dpb->fs_ltref[i] )|| p_Dpb->last_picture->frame_num != curr_frame_num)
                  unmark_for_long_term_reference(p_Dpb->fs_ltref[i]);
              }
              else
              {
                unmark_for_long_term_reference(p_Dpb->fs_ltref[i]);
              }
            }
            else
            {
              if ((p_Dpb->fs_ltref[i]->frame_num) != (unsigned)(curr_pic_num/2))
              {
                unmark_for_long_term_reference(p_Dpb->fs_ltref[i]);
              }
            }
          }
        }
      }
      if (structure == BOTTOM_FIELD)
      {
        if ((p_Dpb->fs_ltref[i]->is_long_term == 3))
        {
          unmark_for_long_term_reference(p_Dpb->fs_ltref[i]);
        }
        else
        {
          if ((p_Dpb->fs_ltref[i]->is_long_term == 2))
          {
            unmark_for_long_term_reference(p_Dpb->fs_ltref[i]);
          }
          else
          {
            if (mark_current)
            {
              if (p_Dpb->last_picture)
              {
                if ( ( p_Dpb->last_picture != p_Dpb->fs_ltref[i] )|| p_Dpb->last_picture->frame_num != curr_frame_num)
                  unmark_for_long_term_reference(p_Dpb->fs_ltref[i]);
              }
              else
              {
                unmark_for_long_term_reference(p_Dpb->fs_ltref[i]);
              }
            }
            else
            {
              if ((p_Dpb->fs_ltref[i]->frame_num) != (unsigned)(curr_pic_num/2))
              {
                unmark_for_long_term_reference(p_Dpb->fs_ltref[i]);
              }
            }
          }
        }
      }
    }
  }
}


/*!
 ************************************************************************
 * \brief
 *    mark a picture as long-term reference
 ************************************************************************
 */
static void mark_pic_long_term(DecodedPictureBuffer *p_Dpb, StorablePicture* p, int long_term_frame_idx, int picNumX)
{
  unsigned i;
  int add_top, add_bottom;

  if (p->structure == FRAME)
  {
    for (i=0; i<p_Dpb->ref_frames_in_buffer; i++)
    {
      if (p_Dpb->fs_ref[i]->is_reference == 3)
      {
        if ((!p_Dpb->fs_ref[i]->frame->is_long_term)&&(p_Dpb->fs_ref[i]->frame->pic_num == picNumX))
        {
          p_Dpb->fs_ref[i]->long_term_frame_idx = p_Dpb->fs_ref[i]->frame->long_term_frame_idx
                                             = long_term_frame_idx;
          p_Dpb->fs_ref[i]->frame->long_term_pic_num = long_term_frame_idx;
          p_Dpb->fs_ref[i]->frame->is_long_term = 1;

          if (p_Dpb->fs_ref[i]->top_field && p_Dpb->fs_ref[i]->bottom_field)
          {
            p_Dpb->fs_ref[i]->top_field->long_term_frame_idx = p_Dpb->fs_ref[i]->bottom_field->long_term_frame_idx
                                                          = long_term_frame_idx;
            p_Dpb->fs_ref[i]->top_field->long_term_pic_num = long_term_frame_idx;
            p_Dpb->fs_ref[i]->bottom_field->long_term_pic_num = long_term_frame_idx;

            p_Dpb->fs_ref[i]->top_field->is_long_term = p_Dpb->fs_ref[i]->bottom_field->is_long_term
                                                   = 1;

          }
          p_Dpb->fs_ref[i]->is_long_term = 3;
          return;
        }
      }
    }
    printf ("Warning: reference frame for long term marking not found\n");
  }
  else
  {
    if (p->structure == TOP_FIELD)
    {
      add_top    = 1;
      add_bottom = 0;
    }
    else
    {
      add_top    = 0;
      add_bottom = 1;
    }
    for (i=0; i<p_Dpb->ref_frames_in_buffer; i++)
    {
      if (p_Dpb->fs_ref[i]->is_reference & 1)
      {
        if ((!p_Dpb->fs_ref[i]->top_field->is_long_term)&&(p_Dpb->fs_ref[i]->top_field->pic_num == picNumX))
        {
          if ((p_Dpb->fs_ref[i]->is_long_term) && (p_Dpb->fs_ref[i]->long_term_frame_idx != long_term_frame_idx))
          {
              printf ("Warning: assigning long_term_frame_idx different from other field\n");
          }

          p_Dpb->fs_ref[i]->long_term_frame_idx = p_Dpb->fs_ref[i]->top_field->long_term_frame_idx
                                             = long_term_frame_idx;
          p_Dpb->fs_ref[i]->top_field->long_term_pic_num = 2 * long_term_frame_idx + add_top;
          p_Dpb->fs_ref[i]->top_field->is_long_term = 1;
          p_Dpb->fs_ref[i]->is_long_term |= 1;
          if (p_Dpb->fs_ref[i]->is_long_term == 3)
          {
            p_Dpb->fs_ref[i]->frame->is_long_term = 1;
            p_Dpb->fs_ref[i]->frame->long_term_frame_idx = p_Dpb->fs_ref[i]->frame->long_term_pic_num = long_term_frame_idx;
          }
          return;
        }
      }
      if (p_Dpb->fs_ref[i]->is_reference & 2)
      {
        if ((!p_Dpb->fs_ref[i]->bottom_field->is_long_term)&&(p_Dpb->fs_ref[i]->bottom_field->pic_num == picNumX))
        {
          if ((p_Dpb->fs_ref[i]->is_long_term) && (p_Dpb->fs_ref[i]->long_term_frame_idx != long_term_frame_idx))
          {
              printf ("Warning: assigning long_term_frame_idx different from other field\n");
          }

          p_Dpb->fs_ref[i]->long_term_frame_idx = p_Dpb->fs_ref[i]->bottom_field->long_term_frame_idx
                                             = long_term_frame_idx;
          p_Dpb->fs_ref[i]->bottom_field->long_term_pic_num = 2 * long_term_frame_idx + add_bottom;
          p_Dpb->fs_ref[i]->bottom_field->is_long_term = 1;
          p_Dpb->fs_ref[i]->is_long_term |= 2;
          if (p_Dpb->fs_ref[i]->is_long_term == 3)
          {
            p_Dpb->fs_ref[i]->frame->is_long_term = 1;
            p_Dpb->fs_ref[i]->frame->long_term_frame_idx = p_Dpb->fs_ref[i]->frame->long_term_pic_num = long_term_frame_idx;
          }
          return;
        }
      }
    }
    printf ("Warning: reference field for long term marking not found\n");
  }
}


/*!
 ************************************************************************
 * \brief
 *    Assign a long term frame index to a short term picture
 ************************************************************************
 */
static void mm_assign_long_term_frame_idx(ImageParameters *p_Img, StorablePicture* p, int difference_of_pic_nums_minus1, int long_term_frame_idx)
{
  int picNumX;
  DecodedPictureBuffer *p_Dpb = p_Img->p_Dpb;

  picNumX = get_pic_num_x(p, difference_of_pic_nums_minus1);

  // remove frames/fields with same long_term_frame_idx
  if (p->structure == FRAME)
  {
    unmark_long_term_frame_for_reference_by_frame_idx(p_Dpb, long_term_frame_idx);
  }
  else
  {
    unsigned i;
    PictureStructure structure = FRAME;

    for (i=0; i<p_Dpb->ref_frames_in_buffer; i++)
    {
      if (p_Dpb->fs_ref[i]->is_reference & 1)
      {
        if (p_Dpb->fs_ref[i]->top_field->pic_num == picNumX)
        {
          structure = TOP_FIELD;
          break;
        }
      }
      if (p_Dpb->fs_ref[i]->is_reference & 2)
      {
        if (p_Dpb->fs_ref[i]->bottom_field->pic_num == picNumX)
        {
          structure = BOTTOM_FIELD;
          break;
        }
      }
    }
    if (structure==FRAME)
    {
      error ("field for long term marking not found",200);
    }

    unmark_long_term_field_for_reference_by_frame_idx(p_Img, structure, long_term_frame_idx, 0, 0, picNumX);
  }

  mark_pic_long_term(p_Dpb, p, long_term_frame_idx, picNumX);
}

/*!
 ************************************************************************
 * \brief
 *    Set new max long_term_frame_idx
 ************************************************************************
 */
void mm_update_max_long_term_frame_idx(DecodedPictureBuffer *p_Dpb, int max_long_term_frame_idx_plus1)
{
  unsigned i;

  p_Dpb->max_long_term_pic_idx = max_long_term_frame_idx_plus1 - 1;

  // check for invalid frames
  for (i=0; i<p_Dpb->ltref_frames_in_buffer; i++)
  {
    if (p_Dpb->fs_ltref[i]->long_term_frame_idx > p_Dpb->max_long_term_pic_idx)
    {
      unmark_for_long_term_reference(p_Dpb->fs_ltref[i]);
    }
  }
}


/*!
 ************************************************************************
 * \brief
 *    Mark all long term reference pictures unused for reference
 ************************************************************************
 */
static void mm_unmark_all_long_term_for_reference (DecodedPictureBuffer *p_Dpb)
{
  mm_update_max_long_term_frame_idx(p_Dpb, 0);
}

/*!
 ************************************************************************
 * \brief
 *    Mark all short term reference pictures unused for reference
 ************************************************************************
 */
static void mm_unmark_all_short_term_for_reference (DecodedPictureBuffer *p_Dpb)
{
  unsigned int i;
  for (i=0; i<p_Dpb->ref_frames_in_buffer; i++)
  {
    unmark_for_reference(p_Dpb->fs_ref[i]);
  }
  update_ref_list(p_Dpb);
}


/*!
 ************************************************************************
 * \brief
 *    Mark the current picture used for long term reference
 ************************************************************************
 */
static void mm_mark_current_picture_long_term(ImageParameters *p_Img, DecodedPictureBuffer *p_Dpb, StorablePicture *p, int long_term_frame_idx)
{
  // remove long term pictures with same long_term_frame_idx
  if (p->structure == FRAME)
  {
    unmark_long_term_frame_for_reference_by_frame_idx(p_Dpb, long_term_frame_idx);
  }
  else
  {
    unmark_long_term_field_for_reference_by_frame_idx(p_Img, p->structure, long_term_frame_idx, 1, p->pic_num, 0);
  }

  p->is_long_term = 1;
  p->long_term_frame_idx = long_term_frame_idx;
}


/*!
 ************************************************************************
 * \brief
 *    Perform Adaptive memory control decoded reference picture marking process
 ************************************************************************
 */
static void adaptive_memory_management(ImageParameters *p_Img, StorablePicture* p)
{
  DecRefPicMarking_t *tmp_drpm;
  DecodedPictureBuffer *p_Dpb = p_Img->p_Dpb;

  p_Img->last_has_mmco_5 = 0;

  assert (!p->idr_flag);
  assert (p->adaptive_ref_pic_buffering_flag);

  while (p->dec_ref_pic_marking_buffer)
  {
    tmp_drpm = p->dec_ref_pic_marking_buffer;
    switch (tmp_drpm->memory_management_control_operation)
    {
      case 0:
        if (tmp_drpm->Next != NULL)
        {
          error ("memory_management_control_operation = 0 not last operation in buffer", 500);
        }
        break;
      case 1:
        mm_unmark_short_term_for_reference(p_Dpb, p, tmp_drpm->difference_of_pic_nums_minus1);
        update_ref_list(p_Dpb);
        break;
      case 2:
        mm_unmark_long_term_for_reference(p_Dpb, p, tmp_drpm->long_term_pic_num);
        update_ltref_list(p_Dpb);
        break;
      case 3:
        mm_assign_long_term_frame_idx(p_Img, p, tmp_drpm->difference_of_pic_nums_minus1, tmp_drpm->long_term_frame_idx);
        update_ref_list(p_Dpb);
        update_ltref_list(p_Dpb);
        break;
      case 4:
        mm_update_max_long_term_frame_idx (p_Dpb, tmp_drpm->max_long_term_frame_idx_plus1);
        update_ltref_list(p_Dpb);
        break;
      case 5:
        mm_unmark_all_short_term_for_reference(p_Dpb);
        mm_unmark_all_long_term_for_reference(p_Dpb);
       p_Img->last_has_mmco_5 = 1;
        break;
      case 6:
        mm_mark_current_picture_long_term(p_Img, p_Dpb, p, tmp_drpm->long_term_frame_idx);
        check_num_ref(p_Dpb);
        break;
      default:
        error ("invalid memory_management_control_operation in buffer", 500);
    }
    p->dec_ref_pic_marking_buffer = tmp_drpm->Next;
    free (tmp_drpm);
  }
  if ( p_Img->last_has_mmco_5 )
  {
    p->pic_num = p->frame_num = 0;

    switch (p->structure)
    {
    case TOP_FIELD:
      {
        p->poc = p->top_poc = p_Img->toppoc =0;
        break;
      }
    case BOTTOM_FIELD:
      {
        p->poc = p->bottom_poc = p_Img->bottompoc = 0;
        break;
      }
    case FRAME:
      {
        p->top_poc    -= p->poc;
        p->bottom_poc -= p->poc;

        p_Img->toppoc = p->top_poc;
        p_Img->bottompoc = p->bottom_poc;

        p->poc = imin (p->top_poc, p->bottom_poc);
        p_Img->framepoc = p->poc;
        break;
      }
    }
    p_Img->ThisPOC = p->poc;
    flush_dpb(p_Img);
  }
}


/*!
 ************************************************************************
 * \brief
 *    Store a picture in DPB. This includes cheking for space in DPB and
 *    flushing frames.
 *    If we received a frame, we need to check for a new store, if we
 *    got a field, check if it's the second field of an already allocated
 *    store.
 *
 * \param p_Img
 *      image decoding parameters for current picture
 * \param p
 *    Picture to be stored
 *
 ************************************************************************
 */
void store_picture_in_dpb(ImageParameters *p_Img, StorablePicture* p)
{
  DecodedPictureBuffer *p_Dpb = p_Img->p_Dpb;
  unsigned i;
  int poc, pos;
  // picture error concealment
  
  // diagnostics
  //printf ("Storing (%s) non-ref pic with frame_num #%d\n", (p->type == FRAME)?"FRAME":(p->type == TOP_FIELD)?"TOP_FIELD":"BOTTOM_FIELD", p->pic_num);
  // if frame, check for new store,
  assert (p!=NULL);

  p_Img->last_has_mmco_5=0;
  p_Img->last_pic_bottom_field = (p->structure == BOTTOM_FIELD);

  if (p->idr_flag)
  {
    idr_memory_management(p_Img, p);
  // picture error concealment
    memset(p_Img->pocs_in_dpb, 0, sizeof(int)*100);
  }
  else
  {
    // adaptive memory management
    if (p->used_for_reference && (p->adaptive_ref_pic_buffering_flag))
      adaptive_memory_management(p_Img, p);
  }

  if ((p->structure==TOP_FIELD)||(p->structure==BOTTOM_FIELD))
  {
    // check for frame store with same pic_number
    if (p_Dpb->last_picture)
    {
      if ((int)p_Dpb->last_picture->frame_num == p->pic_num)
      {
        if (((p->structure==TOP_FIELD)&&(p_Dpb->last_picture->is_used==2))||((p->structure==BOTTOM_FIELD)&&(p_Dpb->last_picture->is_used==1)))
        {
          if ((p->used_for_reference && (p_Dpb->last_picture->is_orig_reference!=0))||
              (!p->used_for_reference && (p_Dpb->last_picture->is_orig_reference==0)))
          {
            insert_picture_in_dpb(p_Img, p_Dpb->last_picture, p);            
            update_ref_list(p_Dpb);
            update_ltref_list(p_Dpb);
            dump_dpb(p_Dpb);
            p_Dpb->last_picture = NULL;
            return;
          }
        }
      }
    }
  }

  // this is a frame or a field which has no stored complementary field

  // sliding window, if necessary
  if ((!p->idr_flag)&&(p->used_for_reference && (!p->adaptive_ref_pic_buffering_flag)))
  {
    sliding_window_memory_management(p_Dpb, p);
  }

  // picture error concealment
  if(p_Img->conceal_mode != 0)
    for(i=0;i<p_Dpb->size;i++)
      if(p_Dpb->fs[i]->is_reference)
        p_Dpb->fs[i]->concealment_reference = 1;

  // first try to remove unused frames
  if (p_Dpb->used_size==p_Dpb->size)
  {
    // picture error concealment
    if (p_Img->conceal_mode != 0)
      conceal_non_ref_pics(p_Img, 2);

    remove_unused_frame_from_dpb(p_Img, p_Dpb);

    if(p_Img->conceal_mode != 0)
      sliding_window_poc_management(p_Dpb, p);
  }

  // then output frames until one can be removed
  while (p_Dpb->used_size == p_Dpb->size)
  {
    // non-reference frames may be output directly
    if (!p->used_for_reference)
    {
      get_smallest_poc(p_Dpb, &poc, &pos);
      if ((-1==pos) || (p->poc < poc))
      {
        direct_output(p_Img, p, p_Img->p_out);
        return;
      }
    }
    // flush a frame
    output_one_frame_from_dpb(p_Img);
  }

  // check for duplicate frame number in short term reference buffer
  if ((p->used_for_reference)&&(!p->is_long_term))
  {
    for (i=0; i<p_Dpb->ref_frames_in_buffer; i++)
    {
      if (p_Dpb->fs_ref[i]->frame_num == p->frame_num)
      {
        error("duplicate frame_num in short-term reference picture buffer", 500);
      }
    }

  }
  // store at end of buffer
  insert_picture_in_dpb(p_Img, p_Dpb->fs[p_Dpb->used_size],p);

  // picture error concealment
  if (p->idr_flag)
  {
      p_Img->earlier_missing_poc = 0;
  }

  if (p->structure != FRAME)
  {
    p_Dpb->last_picture = p_Dpb->fs[p_Dpb->used_size];
  }
  else
  {
    p_Dpb->last_picture = NULL;
  }

  p_Dpb->used_size++;

  if(p_Img->conceal_mode != 0)
      p_Img->pocs_in_dpb[p_Dpb->used_size-1] = p->poc;

  update_ref_list(p_Dpb);
  update_ltref_list(p_Dpb);

  check_num_ref(p_Dpb);

  dump_dpb(p_Dpb);
}

/*!
 ************************************************************************
 * \brief
 *    Insert the picture into the DPB. A free DPB position is necessary
 *    for frames, .
 *
 * \param p_Img
 *      image decoding parameters for current picture
 * \param fs
 *    FrameStore into which the picture will be inserted
 * \param p
 *    StorablePicture to be inserted
 *
 ************************************************************************
 */
static void insert_picture_in_dpb(ImageParameters *p_Img, FrameStore* fs, StorablePicture* p)
{
  InputParameters *p_Inp = p_Img->p_Inp;
//  printf ("insert (%s) pic with frame_num #%d, poc %d\n", (p->structure == FRAME)?"FRAME":(p->structure == TOP_FIELD)?"TOP_FIELD":"BOTTOM_FIELD", p->pic_num, p->poc);
  assert (p!=NULL);
  assert (fs!=NULL);
  switch (p->structure)
  {
  case FRAME:
    fs->frame = p;
    fs->is_used = 3;
    if (p->used_for_reference)
    {
      fs->is_reference = 3;
      fs->is_orig_reference = 3;
      if (p->is_long_term)
      {
        fs->is_long_term = 3;
        fs->long_term_frame_idx = p->long_term_frame_idx;
      }
    }
    // generate field views
    dpb_split_field(p_Img, fs);
    break;
  case TOP_FIELD:
    fs->top_field = p;
    fs->is_used |= 1;
    if (p->used_for_reference)
    {
      fs->is_reference |= 1;
      fs->is_orig_reference |= 1;
      if (p->is_long_term)
      {
        fs->is_long_term |= 1;
        fs->long_term_frame_idx = p->long_term_frame_idx;
      }
    }
    if (fs->is_used == 3)
    {
      // generate frame view
      dpb_combine_field(p_Img, fs);
    } else
    {
      fs->poc = p->poc;
      gen_field_ref_ids(p);
    }
    break;
  case BOTTOM_FIELD:
    fs->bottom_field = p;
    fs->is_used |= 2;
    if (p->used_for_reference)
    {
      fs->is_reference |= 2;
      fs->is_orig_reference |= 2;
      if (p->is_long_term)
      {
        fs->is_long_term |= 2;
        fs->long_term_frame_idx = p->long_term_frame_idx;
      }
    }
    if (fs->is_used == 3)
    {
      // generate frame view
      dpb_combine_field(p_Img, fs);
    } 
    else
    {
      fs->poc = p->poc;
      gen_field_ref_ids(p);
    }
    break;
  }
  fs->frame_num = p->pic_num;
  fs->recovery_frame = p->recovery_frame;

  fs->is_output = p->is_output;

  if (fs->is_used==3)
  {
    calculate_frame_no(p_Img, p);
    if (-1 != p_Img->p_ref && !p_Inp->silent)
      find_snr(p_Img, fs->frame, &p_Img->p_ref);
  }
}

/*!
 ************************************************************************
 * \brief
 *    Check if one of the frames/fields in frame store is used for reference
 ************************************************************************
 */
static int is_used_for_reference(FrameStore* fs)
{
  if (fs->is_reference)
  {
    return 1;
  }

  if (fs->is_used == 3) // frame
  {
    if (fs->frame->used_for_reference)
    {
      return 1;
    }
  }

  if (fs->is_used & 1) // top field
  {
    if (fs->top_field)
    {
      if (fs->top_field->used_for_reference)
      {
        return 1;
      }
    }
  }

  if (fs->is_used & 2) // bottom field
  {
    if (fs->bottom_field)
    {
      if (fs->bottom_field->used_for_reference)
      {
        return 1;
      }
    }
  }
  return 0;
}


/*!
 ************************************************************************
 * \brief
 *    Check if one of the frames/fields in frame store is used for short-term reference
 ************************************************************************
 */
static int is_short_term_reference(FrameStore* fs)
{

  if (fs->is_used==3) // frame
  {
    if ((fs->frame->used_for_reference)&&(!fs->frame->is_long_term))
    {
      return 1;
    }
  }

  if (fs->is_used & 1) // top field
  {
    if (fs->top_field)
    {
      if ((fs->top_field->used_for_reference)&&(!fs->top_field->is_long_term))
      {
        return 1;
      }
    }
  }

  if (fs->is_used & 2) // bottom field
  {
    if (fs->bottom_field)
    {
      if ((fs->bottom_field->used_for_reference)&&(!fs->bottom_field->is_long_term))
      {
        return 1;
      }
    }
  }
  return 0;
}


/*!
 ************************************************************************
 * \brief
 *    Check if one of the frames/fields in frame store is used for short-term reference
 ************************************************************************
 */
static int is_long_term_reference(FrameStore* fs)
{

  if (fs->is_used==3) // frame
  {
    if ((fs->frame->used_for_reference)&&(fs->frame->is_long_term))
    {
      return 1;
    }
  }

  if (fs->is_used & 1) // top field
  {
    if (fs->top_field)
    {
      if ((fs->top_field->used_for_reference)&&(fs->top_field->is_long_term))
      {
        return 1;
      }
    }
  }

  if (fs->is_used & 2) // bottom field
  {
    if (fs->bottom_field)
    {
      if ((fs->bottom_field->used_for_reference)&&(fs->bottom_field->is_long_term))
      {
        return 1;
      }
    }
  }
  return 0;
}


/*!
 ************************************************************************
 * \brief
 *    remove one frame from DPB
 ************************************************************************
 */
static void remove_frame_from_dpb(ImageParameters *p_Img, int pos)
{
  DecodedPictureBuffer *p_Dpb = p_Img->p_Dpb;
  FrameStore* fs = p_Dpb->fs[pos];
  FrameStore* tmp;
  unsigned i;

//  printf ("remove frame with frame_num #%d\n", fs->frame_num);
  switch (fs->is_used)
  {
  case 3:
    free_storable_picture(p_Img, fs->frame);
    free_storable_picture(p_Img, fs->top_field);
    free_storable_picture(p_Img, fs->bottom_field);
    fs->frame=NULL;
    fs->top_field=NULL;
    fs->bottom_field=NULL;
    break;
  case 2:
    free_storable_picture(p_Img, fs->bottom_field);
    fs->bottom_field=NULL;
    break;
  case 1:
    free_storable_picture(p_Img, fs->top_field);
    fs->top_field=NULL;
    break;
  case 0:
    break;
  default:
    error("invalid frame store type",500);
  }
  fs->is_used = 0;
  fs->is_long_term = 0;
  fs->is_reference = 0;
  fs->is_orig_reference = 0;

  // move empty framestore to end of buffer
  tmp = p_Dpb->fs[pos];

  for (i=pos; i<p_Dpb->used_size-1;i++)
  {
    p_Dpb->fs[i] = p_Dpb->fs[i+1];
  }
  p_Dpb->fs[p_Dpb->used_size-1] = tmp;
  p_Dpb->used_size--;
}

/*!
 ************************************************************************
 * \brief
 *    find smallest POC in the DPB.
 ************************************************************************
 */
static void get_smallest_poc(DecodedPictureBuffer *p_Dpb, int *poc,int * pos)
{
  unsigned i;

  if (p_Dpb->used_size<1)
  {
    error("Cannot determine smallest POC, DPB empty.",150);
  }

  *pos=-1;
  *poc = INT_MAX;
  for (i=0; i<p_Dpb->used_size; i++)
  {
    if ((*poc > p_Dpb->fs[i]->poc)&&(!p_Dpb->fs[i]->is_output))
    {
      *poc = p_Dpb->fs[i]->poc;
      *pos=i;
    }
  }
}

/*!
 ************************************************************************
 * \brief
 *    Remove a picture from DPB which is no longer needed.
 ************************************************************************
 */
static int remove_unused_frame_from_dpb(ImageParameters *p_Img, DecodedPictureBuffer *p_Dpb)
{
  unsigned i;

  // check for frames that were already output and no longer used for reference
  for (i = 0; i < p_Dpb->used_size; i++)
  {
    if (p_Dpb->fs[i]->is_output && (!is_used_for_reference(p_Dpb->fs[i])))
    {
      remove_frame_from_dpb(p_Img, i);
      return 1;
    }
  }
  return 0;
}

/*!
 ************************************************************************
 * \brief
 *    Output one picture stored in the DPB.
 ************************************************************************
 */
static void output_one_frame_from_dpb(ImageParameters *p_Img)
{
  DecodedPictureBuffer *p_Dpb = p_Img->p_Dpb;
  int poc, pos;
  //diagnostics
  if (p_Dpb->used_size<1)
  {
    error("Cannot output frame, DPB empty.",150);
  }

  // find smallest POC
  get_smallest_poc(p_Dpb, &poc, &pos);

  if(pos==-1)
  {
    error("no frames for output available", 150);
  }

  // call the output function
//  printf ("output frame with frame_num #%d, poc %d (dpb. p_Dpb->size=%d, p_Dpb->used_size=%d)\n", p_Dpb->fs[pos]->frame_num, p_Dpb->fs[pos]->frame->poc, p_Dpb->size, p_Dpb->used_size);

  // picture error concealment
  if(p_Img->conceal_mode != 0)
  {
    if(p_Dpb->last_output_poc == 0)
    {
      write_lost_ref_after_idr(p_Img, pos);
    }
    write_lost_non_ref_pic(p_Img, poc, p_Img->p_out);
  }

// JVT-P072 ends

  write_stored_frame(p_Img, p_Dpb->fs[pos], p_Img->p_out);

  // picture error concealment
  if(p_Img->conceal_mode == 0)
    if (p_Dpb->last_output_poc >= poc)
    {
      error ("output POC must be in ascending order", 150);
    }
  p_Dpb->last_output_poc = poc;
  // free frame store and move empty store to end of buffer
  if (!is_used_for_reference(p_Dpb->fs[pos]))
  {
    remove_frame_from_dpb(p_Img, pos);
  }
}



/*!
 ************************************************************************
 * \brief
 *    All stored picture are output. Should be called to empty the buffer
 ************************************************************************
 */
void flush_dpb(ImageParameters *p_Img)
{
  unsigned i;
  DecodedPictureBuffer *p_Dpb = p_Img->p_Dpb;

  //diagnostics
//  printf("Flush remaining frames from dpb. p_Dpb->size=%d, p_Dpb->used_size=%d\n",p_Dpb->size,p_Dpb->used_size);

//  if(p_Img->conceal_mode == 0)
  if (p_Img->conceal_mode != 0)
    conceal_non_ref_pics(p_Img, 0);

  // mark all frames unused
  for (i=0; i<p_Dpb->used_size; i++)
  {
    unmark_for_reference (p_Dpb->fs[i]);
  }

  while (remove_unused_frame_from_dpb(p_Img, p_Dpb)) ;

  // output frames in POC order
  while (p_Dpb->used_size)
  {
    output_one_frame_from_dpb(p_Img);
  }

  p_Dpb->last_output_poc = INT_MIN;
}


static void gen_field_ref_ids(StorablePicture *p)
{
  int i,j, dummylist0, dummylist1;
   //! Generate Frame parameters from field information.
  for (i=0 ; i<p->size_x/4 ; i++)
  {
    for (j=0 ; j<p->size_y/4 ; j++)
    {
        dummylist0= p->motion.ref_idx[LIST_0][j][i];
        dummylist1= p->motion.ref_idx[LIST_1][j][i];
        //! association with id already known for fields.
        p->motion.ref_id[LIST_0][j][i] = (dummylist0>=0)? p->ref_pic_num[p->slice_id[j>>2][i>>2]][LIST_0][dummylist0] : 0;
        p->motion.ref_id[LIST_1][j][i] = (dummylist1>=0)? p->ref_pic_num[p->slice_id[j>>2][i>>2]][LIST_1][dummylist1] : 0;
        p->motion.field_frame[j][i]=1;
    }
  }
}

/*!
 ************************************************************************
 * \brief
 *    Extract top field from a frame
 ************************************************************************
 */
void dpb_split_field(ImageParameters *p_Img, FrameStore *fs)
{
  int i, j, ii, jj, jj4;
  int idiv,jdiv;
  int currentmb;
  int dummylist0, dummylist1;
  int twosz16 = 2 * (fs->frame->size_x >> 4);
  StorablePicture *fs_top, *fs_btm; 
  StorablePicture *frame = fs->frame;


  fs->poc = frame->poc;

  if (!frame->frame_mbs_only_flag)
  {
    fs_top = fs->top_field    = alloc_storable_picture(p_Img, TOP_FIELD,    frame->size_x, frame->size_y, frame->size_x_cr, frame->size_y_cr);
    fs_btm = fs->bottom_field = alloc_storable_picture(p_Img, BOTTOM_FIELD, frame->size_x, frame->size_y, frame->size_x_cr, frame->size_y_cr);

    for (i = 0; i < (frame->size_y>>1); i++)
    {
      memcpy(fs_top->imgY[i], frame->imgY[i*2], frame->size_x*sizeof(imgpel));
    }

    for (i = 0; i< (frame->size_y_cr>>1); i++)
    {
      memcpy(fs_top->imgUV[0][i], frame->imgUV[0][i*2], frame->size_x_cr*sizeof(imgpel));
      memcpy(fs_top->imgUV[1][i], frame->imgUV[1][i*2], frame->size_x_cr*sizeof(imgpel));
    }

    for (i = 0; i < (frame->size_y>>1); i++)
    {
      memcpy(fs_btm->imgY[i], frame->imgY[i*2 + 1], frame->size_x*sizeof(imgpel));
    }

    for (i = 0; i < (frame->size_y_cr>>1); i++)
    {
      memcpy(fs_btm->imgUV[0][i], frame->imgUV[0][i*2 + 1], frame->size_x_cr*sizeof(imgpel));
      memcpy(fs_btm->imgUV[1][i], frame->imgUV[1][i*2 + 1], frame->size_x_cr*sizeof(imgpel));
    }

    fs_top->poc = frame->top_poc;
    fs_btm->poc = frame->bottom_poc;

    fs_top->frame_poc =  frame->frame_poc;

    fs_top->bottom_poc = fs_btm->bottom_poc =  frame->bottom_poc;
    fs_top->top_poc    = fs_btm->top_poc    =  frame->top_poc;
    fs_btm->frame_poc  = frame->frame_poc;

    fs_top->used_for_reference = fs_btm->used_for_reference
                                      = frame->used_for_reference;
    fs_top->is_long_term = fs_btm->is_long_term
                                = frame->is_long_term;
    fs->long_term_frame_idx = fs_top->long_term_frame_idx
                            = fs_btm->long_term_frame_idx
                            = frame->long_term_frame_idx;

    fs_top->coded_frame = fs_btm->coded_frame = 1;
    fs_top->MbaffFrameFlag = fs_btm->MbaffFrameFlag
                                  = frame->MbaffFrameFlag;

    frame->top_field    = fs_top;
    frame->bottom_field = fs_btm;

    fs_top->bottom_field = fs_btm;
    fs_top->frame        = frame;
    fs_btm->top_field = fs_top;
    fs_btm->frame     = frame;

    fs_top->chroma_format_idc = fs_btm->chroma_format_idc = frame->chroma_format_idc;

    //store reference picture index
    for (j=0; j<=frame->max_slice_id; j++)
    {
      memcpy(&fs_top->ref_pic_num[j][LIST_0][0], &frame->ref_pic_num[j][2 + LIST_0][0], 66 * sizeof(int64));
      //memcpy(&fs_top->ref_pic_num[j][LIST_1][0], &frame->ref_pic_num[j][2 + LIST_1][0], 33 * sizeof(int64));            
      memcpy(&fs_btm->ref_pic_num[j][LIST_0][0], &frame->ref_pic_num[j][4 + LIST_0][0], 66 * sizeof(int64));
      //memcpy(&fs_btm->ref_pic_num[j][LIST_1][0], &frame->ref_pic_num[j][4 + LIST_1][0], 33 * sizeof(int64));
    }
  }
  else
  {
    fs_top=NULL;
    fs_btm=NULL;
    frame->top_field=NULL;
    frame->bottom_field=NULL;
  }

  if (!frame->MbaffFrameFlag)
  {
    for (j = 0; (j < frame->size_y >> 2) ; j++)
    {
      jdiv = j >> 2;
      for (i = 0 ; i < (frame->size_x >> 2) ; i++)
      {
        idiv = (i >> 2);

        dummylist0 = frame->motion.ref_idx[LIST_0][j][i];
        dummylist1 = frame->motion.ref_idx[LIST_1][j][i];
        frame->motion.ref_id[LIST_0][j][i] = (dummylist0>=0)? frame->ref_pic_num[frame->slice_id[jdiv][idiv]][LIST_0][dummylist0] : -1;
        frame->motion.ref_id[LIST_1][j][i] = (dummylist1>=0)? frame->ref_pic_num[frame->slice_id[jdiv][idiv]][LIST_1][dummylist1] : -1;
      }
    }
  }
  else
  {
    for (j = 0; (j < frame->size_y >> 2) ; j++)
    {
      jdiv = j >> 2;
      for (i = 0 ; i < (frame->size_x >> 2) ; i++)
      {
        idiv = (i >> 2);
        currentmb = twosz16*(jdiv >> 1)+ (idiv)*2 + (jdiv & 0x01);

        if (frame->motion.mb_field[currentmb])
        {
          int list_offset = currentmb%2? 4: 2;
          dummylist0 = frame->motion.ref_idx[LIST_0][j][i];
          dummylist1 = frame->motion.ref_idx[LIST_1][j][i];
          //! association with id already known for fields.
          frame->motion.ref_id[LIST_0 + list_offset][j][i] = (dummylist0>=0)? frame->ref_pic_num[frame->slice_id[jdiv][idiv]][LIST_0 + list_offset][dummylist0] : 0;
          frame->motion.ref_id[LIST_1 + list_offset][j][i] = (dummylist1>=0)? frame->ref_pic_num[frame->slice_id[jdiv][idiv]][LIST_1 + list_offset][dummylist1] : 0;
          //! need to make association with frames
          frame->motion.ref_id[LIST_0][j][i] = (dummylist0>=0)? frame->frm_ref_pic_num[frame->slice_id[jdiv][idiv]][LIST_0 + list_offset][dummylist0] : 0;
          frame->motion.ref_id[LIST_1][j][i] = (dummylist1>=0)? frame->frm_ref_pic_num[frame->slice_id[jdiv][idiv]][LIST_1 + list_offset][dummylist1] : 0;

        }
        else
        {
          dummylist0 = frame->motion.ref_idx[LIST_0][j][i];
          dummylist1 = frame->motion.ref_idx[LIST_1][j][i];
          frame->motion.ref_id[LIST_0][j][i] = (dummylist0>=0)? frame->ref_pic_num[frame->slice_id[jdiv][idiv]][LIST_0][dummylist0] : -1;
          frame->motion.ref_id[LIST_1][j][i] = (dummylist1>=0)? frame->ref_pic_num[frame->slice_id[jdiv][idiv]][LIST_1][dummylist1] : -1;
        }
      }
    }
  }

  if (!frame->frame_mbs_only_flag && frame->MbaffFrameFlag)
  {
    PicMotionParams *frm_motion = &frame->motion;
    PicMotionParams *top_motion = &fs_top->motion;
    PicMotionParams *btm_motion = &fs_btm->motion;
    for (j=0 ; j< (frame->size_y >> 3); j++)
    {
      jj = (j >> 2)*8 + (j & 0x03);
      jj4 = jj + 4;
      jdiv = (j >> 1);
      for (i=0 ; i < (frame->size_x>>2); i++)
      {
        idiv = (i >> 2);

        currentmb = twosz16*(jdiv >> 1)+ (idiv)*2 + (jdiv & 0x01);
        // Assign field mvs attached to MB-Frame buffer to the proper buffer
        if (frm_motion->mb_field[currentmb])
        {
          btm_motion->field_frame[j][i]  = top_motion->field_frame[j][i]=1;
          frm_motion->field_frame[2*j][i] = frm_motion->field_frame[2*j+1][i]=1;

          btm_motion->mv[LIST_0][j][i][0] = frm_motion->mv[LIST_0][jj4][i][0];
          btm_motion->mv[LIST_0][j][i][1] = frm_motion->mv[LIST_0][jj4][i][1];
          btm_motion->mv[LIST_1][j][i][0] = frm_motion->mv[LIST_1][jj4][i][0];
          btm_motion->mv[LIST_1][j][i][1] = frm_motion->mv[LIST_1][jj4][i][1];
          btm_motion->ref_idx[LIST_0][j][i] = frm_motion->ref_idx[LIST_0][jj4][i];
          btm_motion->ref_idx[LIST_1][j][i] = frm_motion->ref_idx[LIST_1][jj4][i];
          btm_motion->ref_id[LIST_0][j][i] = frm_motion->ref_id[LIST_0+4][jj4][i];
          btm_motion->ref_id[LIST_1][j][i] = frm_motion->ref_id[LIST_1+4][jj4][i];


          top_motion->mv[LIST_0][j][i][0] = frm_motion->mv[LIST_0][jj][i][0];
          top_motion->mv[LIST_0][j][i][1] = frm_motion->mv[LIST_0][jj][i][1];
          top_motion->mv[LIST_1][j][i][0] = frm_motion->mv[LIST_1][jj][i][0];
          top_motion->mv[LIST_1][j][i][1] = frm_motion->mv[LIST_1][jj][i][1];
          top_motion->ref_idx[LIST_0][j][i] = frm_motion->ref_idx[LIST_0][jj][i];
          top_motion->ref_idx[LIST_1][j][i] = frm_motion->ref_idx[LIST_1][jj][i];
          top_motion->ref_id[LIST_0][j][i] = frm_motion->ref_id[LIST_0+2][jj][i];
          top_motion->ref_id[LIST_1][j][i] = frm_motion->ref_id[LIST_1+2][jj][i];
        }
      }
    }
  }

  //! Generate field MVs from Frame MVs
  if (!frame->frame_mbs_only_flag)
  {
    for (j=0 ; j < (frame->size_y >> 3) ; j++)
    {
      jj = 2* RSD(j);
      jdiv = (j >> 1);
      for (i=0 ; i < (frame->size_x >> 2) ; i++)
      {
        ii = RSD(i);
        idiv = (i >> 2);

        currentmb = twosz16 * (jdiv >> 1)+ (idiv)*2 + (jdiv & 0x01);

        if (!frame->MbaffFrameFlag  || !frame->motion.mb_field[currentmb])
        {
          frame->motion.field_frame[2*j+1][i] = frame->motion.field_frame[2*j][i]=0;

          fs_top->motion.field_frame[j][i] = fs_btm->motion.field_frame[j][i] = 0;

          fs_top->motion.mv[LIST_0][j][i][0] = fs_btm->motion.mv[LIST_0][j][i][0] = frame->motion.mv[LIST_0][jj][ii][0];
          fs_top->motion.mv[LIST_0][j][i][1] = fs_btm->motion.mv[LIST_0][j][i][1] = frame->motion.mv[LIST_0][jj][ii][1];
          fs_top->motion.mv[LIST_1][j][i][0] = fs_btm->motion.mv[LIST_1][j][i][0] = frame->motion.mv[LIST_1][jj][ii][0];
          fs_top->motion.mv[LIST_1][j][i][1] = fs_btm->motion.mv[LIST_1][j][i][1] = frame->motion.mv[LIST_1][jj][ii][1];

          // Scaling of references is done here since it will not affect spatial direct (2*0 =0)
          if (frame->motion.ref_idx[LIST_0][jj][ii] == -1)
            fs_top->motion.ref_idx[LIST_0][j][i] = fs_btm->motion.ref_idx[LIST_0][j][i] = - 1;
          else
          {
            dummylist0=fs_top->motion.ref_idx[LIST_0][j][i] = fs_btm->motion.ref_idx[LIST_0][j][i] = frame->motion.ref_idx[LIST_0][jj][ii] ;
            fs_top->motion.ref_id[LIST_0][j][i] = (dummylist0>=0)? frame->top_ref_pic_num[frame->slice_id[jj>>2][ii>>2]][LIST_0][dummylist0] : 0;
            fs_btm->motion.ref_id[LIST_0][j][i] = (dummylist0>=0)? frame->bottom_ref_pic_num[frame->slice_id[jj>>2][ii>>2]][LIST_0][dummylist0] : 0;
          }

          if (frame->motion.ref_idx[LIST_1][jj][ii] == -1)
            fs_top->motion.ref_idx[LIST_1][j][i] = fs_btm->motion.ref_idx[LIST_1][j][i] = - 1;
          else
          {
            dummylist1=fs_top->motion.ref_idx[LIST_1][j][i] = fs_btm->motion.ref_idx[LIST_1][j][i] = frame->motion.ref_idx[LIST_1][jj][ii];

            fs_top->motion.ref_id[LIST_1][j][i] = (dummylist1>=0)? frame->top_ref_pic_num[frame->slice_id[jj>>2][ii>>2]][LIST_1][dummylist1] : 0;
            fs_btm->motion.ref_id[LIST_1][j][i] = (dummylist1>=0)? frame->bottom_ref_pic_num[frame->slice_id[jj>>2][ii>>2]][LIST_1][dummylist1] : 0;
          }
        }
        else
        {
          frame->motion.field_frame[2*j+1][i] = frame->motion.field_frame[2*j][i]= frame->motion.mb_field[currentmb];
        }
      }
    }
  }
  else
  {
    memset( &(frame->motion.field_frame[0][0]), 0, (frame->size_y * frame->size_x >> 4) * sizeof(byte));
  }
}


/*!
 ************************************************************************
 * \brief
 *    Generate a frame from top and bottom fields,
 *    YUV components and display information only
 ************************************************************************
 */
void dpb_combine_field_yuv(ImageParameters *p_Img, FrameStore *fs)
{
  int i, j;

  fs->frame = alloc_storable_picture(p_Img, FRAME, fs->top_field->size_x, fs->top_field->size_y*2, fs->top_field->size_x_cr, fs->top_field->size_y_cr*2);

  for (i=0; i<fs->top_field->size_y; i++)
  {
    memcpy(fs->frame->imgY[i*2],     fs->top_field->imgY[i]   , fs->top_field->size_x * sizeof(imgpel));     // top field
    memcpy(fs->frame->imgY[i*2 + 1], fs->bottom_field->imgY[i], fs->bottom_field->size_x * sizeof(imgpel)); // bottom field
  }

  for (j = 0; j < 2; j++)
  {
    for (i=0; i<fs->top_field->size_y_cr; i++)
    {
      memcpy(fs->frame->imgUV[j][i*2],     fs->top_field->imgUV[j][i],    fs->top_field->size_x_cr*sizeof(imgpel));
      memcpy(fs->frame->imgUV[j][i*2 + 1], fs->bottom_field->imgUV[j][i], fs->bottom_field->size_x_cr*sizeof(imgpel));
    }
  }

  fs->poc=fs->frame->poc =fs->frame->frame_poc = imin (fs->top_field->poc, fs->bottom_field->poc);

  fs->bottom_field->frame_poc=fs->top_field->frame_poc=fs->frame->poc;

  fs->bottom_field->top_poc=fs->frame->top_poc=fs->top_field->poc;
  fs->top_field->bottom_poc=fs->frame->bottom_poc=fs->bottom_field->poc;

  fs->frame->used_for_reference = (fs->top_field->used_for_reference && fs->bottom_field->used_for_reference );
  fs->frame->is_long_term = (fs->top_field->is_long_term && fs->bottom_field->is_long_term );

  if (fs->frame->is_long_term)
    fs->frame->long_term_frame_idx = fs->long_term_frame_idx;

  fs->frame->top_field    = fs->top_field;
  fs->frame->bottom_field = fs->bottom_field;

  fs->frame->coded_frame = 0;

  fs->frame->chroma_format_idc = fs->top_field->chroma_format_idc;
  fs->frame->frame_cropping_flag = fs->top_field->frame_cropping_flag;
  if (fs->frame->frame_cropping_flag)
  {
    fs->frame->frame_cropping_rect_top_offset = fs->top_field->frame_cropping_rect_top_offset;
    fs->frame->frame_cropping_rect_bottom_offset = fs->top_field->frame_cropping_rect_bottom_offset;
    fs->frame->frame_cropping_rect_left_offset = fs->top_field->frame_cropping_rect_left_offset;
    fs->frame->frame_cropping_rect_right_offset = fs->top_field->frame_cropping_rect_right_offset;
  }

  fs->top_field->frame = fs->bottom_field->frame = fs->frame;
}


/*!
 ************************************************************************
 * \brief
 *    Generate a frame from top and bottom fields
 ************************************************************************
 */
void dpb_combine_field(ImageParameters *p_Img, FrameStore *fs)
{
  int i,j, k, jj, jj4;
  int dummylist0, dummylist1;

  dpb_combine_field_yuv(p_Img, fs);


  //combine field for frame
  for (j=0; j<=(imax(fs->top_field->max_slice_id, fs->bottom_field->max_slice_id)); j++)
  {
    for (k = LIST_0; k <= LIST_1; k++)
    {
      for (i=0;i<16;i++)
      {
        fs->frame->ref_pic_num[j][k][i]=  i64min ((fs->top_field->ref_pic_num[j][k][2*i]/2)*2, (fs->bottom_field->ref_pic_num[j][k][2*i]/2)*2);
      }
    }
  }

   //! Use inference flag to remap mvs/references

  //! Generate Frame parameters from field information.
  for (j=0 ; j < (fs->top_field->size_y >> 2) ; j++)
  {
    jj = 8*(j >> 2) + (j & 0x03);
    jj4 = jj + 4;
    for (i=0 ; i< (fs->top_field->size_x >> 2) ; i++)
    {
      fs->frame->motion.field_frame[jj][i]= fs->frame->motion.field_frame[jj4][i]=1;

      fs->frame->motion.mv[LIST_0][jj][i][0] = fs->top_field->motion.mv[LIST_0][j][i][0];
      fs->frame->motion.mv[LIST_0][jj][i][1] = fs->top_field->motion.mv[LIST_0][j][i][1];
      fs->frame->motion.mv[LIST_1][jj][i][0] = fs->top_field->motion.mv[LIST_1][j][i][0];
      fs->frame->motion.mv[LIST_1][jj][i][1] = fs->top_field->motion.mv[LIST_1][j][i][1];

      dummylist0=fs->frame->motion.ref_idx[LIST_0][jj][i]  = fs->top_field->motion.ref_idx[LIST_0][j][i];
      dummylist1=fs->frame->motion.ref_idx[LIST_1][jj][i]  = fs->top_field->motion.ref_idx[LIST_1][j][i];

      //! association with id already known for fields.
      fs->top_field->motion.ref_id[LIST_0][j][i] = (dummylist0>=0)? fs->top_field->ref_pic_num[fs->top_field->slice_id[j>>2][i>>2]][LIST_0][dummylist0] : 0;
      fs->top_field->motion.ref_id[LIST_1][j][i] = (dummylist1>=0)? fs->top_field->ref_pic_num[fs->top_field->slice_id[j>>2][i>>2]][LIST_1][dummylist1] : 0;

      //! need to make association with frames
      fs->frame->motion.ref_id[LIST_0][jj][i] = (dummylist0>=0)? fs->top_field->frm_ref_pic_num[fs->top_field->slice_id[j>>2][i>>2]][LIST_0][dummylist0] : 0;
      fs->frame->motion.ref_id[LIST_1][jj][i] = (dummylist1>=0)? fs->top_field->frm_ref_pic_num[fs->top_field->slice_id[j>>2][i>>2]][LIST_1][dummylist1] : 0;

      fs->frame->motion.mv[LIST_0][jj4][i][0] = fs->bottom_field->motion.mv[LIST_0][j][i][0];
      fs->frame->motion.mv[LIST_0][jj4][i][1] = fs->bottom_field->motion.mv[LIST_0][j][i][1] ;
      fs->frame->motion.mv[LIST_1][jj4][i][0] = fs->bottom_field->motion.mv[LIST_1][j][i][0];
      fs->frame->motion.mv[LIST_1][jj4][i][1] = fs->bottom_field->motion.mv[LIST_1][j][i][1] ;

      dummylist0=fs->frame->motion.ref_idx[LIST_0][jj4][i]  = fs->bottom_field->motion.ref_idx[LIST_0][j][i];
      dummylist1=fs->frame->motion.ref_idx[LIST_1][jj4][i]  = fs->bottom_field->motion.ref_idx[LIST_1][j][i];

      fs->bottom_field->motion.ref_id[LIST_0][j][i] = (dummylist0>=0)? fs->bottom_field->ref_pic_num[fs->bottom_field->slice_id[j>>2][i>>2]][LIST_0][dummylist0] : 0;
      fs->bottom_field->motion.ref_id[LIST_1][j][i] = (dummylist1>=0)? fs->bottom_field->ref_pic_num[fs->bottom_field->slice_id[j>>2][i>>2]][LIST_1][dummylist1] : 0;

      //! need to make association with frames
      fs->frame->motion.ref_id[LIST_0][jj4][i] = (dummylist0>=0)? fs->bottom_field->frm_ref_pic_num[fs->bottom_field->slice_id[j>>2][i>>2]][LIST_0][dummylist0] : -1;
      fs->frame->motion.ref_id[LIST_1][jj4][i] = (dummylist1>=0)? fs->bottom_field->frm_ref_pic_num[fs->bottom_field->slice_id[j>>2][i>>2]][LIST_1][dummylist1] : -1;

      fs->top_field->motion.field_frame[j][i]=1;
      fs->bottom_field->motion.field_frame[j][i]=1;
    }
  }
}


/*!
 ************************************************************************
 * \brief
 *    Allocate memory for buffering of reference picture reordering commands
 ************************************************************************
 */
void alloc_ref_pic_list_reordering_buffer(Slice *currSlice)
{
  ImageParameters *p_Img = currSlice->p_Img;
  int size = currSlice->num_ref_idx_l0_active + 1;

  if (p_Img->type!=I_SLICE && p_Img->type!=SI_SLICE)
  {
    if ((currSlice->reordering_of_pic_nums_idc_l0 = calloc(size,sizeof(int)))==NULL) no_mem_exit("alloc_ref_pic_list_reordering_buffer: reordering_of_pic_nums_idc_l0");
    if ((currSlice->abs_diff_pic_num_minus1_l0 = calloc(size,sizeof(int)))==NULL) no_mem_exit("alloc_ref_pic_list_reordering_buffer: abs_diff_pic_num_minus1_l0");
    if ((currSlice->long_term_pic_idx_l0 = calloc(size,sizeof(int)))==NULL) no_mem_exit("alloc_ref_pic_list_reordering_buffer: long_term_pic_idx_l0");
  }
  else
  {
    currSlice->reordering_of_pic_nums_idc_l0 = NULL;
    currSlice->abs_diff_pic_num_minus1_l0 = NULL;
    currSlice->long_term_pic_idx_l0 = NULL;
  }

  size = currSlice->num_ref_idx_l1_active+1;

  if (p_Img->type==B_SLICE)
  {
    if ((currSlice->reordering_of_pic_nums_idc_l1 = calloc(size,sizeof(int)))==NULL) no_mem_exit("alloc_ref_pic_list_reordering_buffer: reordering_of_pic_nums_idc_l1");
    if ((currSlice->abs_diff_pic_num_minus1_l1 = calloc(size,sizeof(int)))==NULL) no_mem_exit("alloc_ref_pic_list_reordering_buffer: abs_diff_pic_num_minus1_l1");
    if ((currSlice->long_term_pic_idx_l1 = calloc(size,sizeof(int)))==NULL) no_mem_exit("alloc_ref_pic_list_reordering_buffer: long_term_pic_idx_l1");
  }
  else
  {
    currSlice->reordering_of_pic_nums_idc_l1 = NULL;
    currSlice->abs_diff_pic_num_minus1_l1 = NULL;
    currSlice->long_term_pic_idx_l1 = NULL;
  }
}


/*!
 ************************************************************************
 * \brief
 *    Free memory for buffering of reference picture reordering commands
 ************************************************************************
 */
void free_ref_pic_list_reordering_buffer(Slice *currSlice)
{

  if (currSlice->reordering_of_pic_nums_idc_l0)
    free(currSlice->reordering_of_pic_nums_idc_l0);
  if (currSlice->abs_diff_pic_num_minus1_l0)
    free(currSlice->abs_diff_pic_num_minus1_l0);
  if (currSlice->long_term_pic_idx_l0)
    free(currSlice->long_term_pic_idx_l0);

  currSlice->reordering_of_pic_nums_idc_l0 = NULL;
  currSlice->abs_diff_pic_num_minus1_l0 = NULL;
  currSlice->long_term_pic_idx_l0 = NULL;

  if (currSlice->reordering_of_pic_nums_idc_l1)
    free(currSlice->reordering_of_pic_nums_idc_l1);
  if (currSlice->abs_diff_pic_num_minus1_l1)
    free(currSlice->abs_diff_pic_num_minus1_l1);
  if (currSlice->long_term_pic_idx_l1)
    free(currSlice->long_term_pic_idx_l1);

  currSlice->reordering_of_pic_nums_idc_l1 = NULL;
  currSlice->abs_diff_pic_num_minus1_l1 = NULL;
  currSlice->long_term_pic_idx_l1 = NULL;
}

/*!
 ************************************************************************
 * \brief
 *      Tian Dong
 *          June 13, 2002, Modifed on July 30, 2003
 *
 *      If a gap in frame_num is found, try to fill the gap
 * \param p_Img
 *
 ************************************************************************
 */
void fill_frame_num_gap(ImageParameters *p_Img)
{
  seq_parameter_set_rbsp_t *active_sps = p_Img->active_sps;
  
  int CurrFrameNum;
  int UnusedShortTermFrameNum;
  StorablePicture *picture = NULL;
  int tmp1 = p_Img->delta_pic_order_cnt[0];
  int tmp2 = p_Img->delta_pic_order_cnt[1];
  p_Img->delta_pic_order_cnt[0] = p_Img->delta_pic_order_cnt[1] = 0;

//  printf("A gap in frame number is found, try to fill it.\n");

  UnusedShortTermFrameNum = (p_Img->pre_frame_num + 1) % p_Img->MaxFrameNum;
  CurrFrameNum = p_Img->frame_num;

  while (CurrFrameNum != UnusedShortTermFrameNum)
  {
    picture = alloc_storable_picture (p_Img, FRAME, p_Img->width, p_Img->height, p_Img->width_cr, p_Img->height_cr);
    picture->coded_frame = 1;
    picture->pic_num = UnusedShortTermFrameNum;
    picture->frame_num = UnusedShortTermFrameNum;
    picture->non_existing = 1;
    picture->is_output = 1;
    picture->used_for_reference = 1;

    picture->adaptive_ref_pic_buffering_flag = 0;

    p_Img->frame_num = UnusedShortTermFrameNum;
    if (active_sps->pic_order_cnt_type!=0)
    {
      decode_poc(p_Img);
    }
    picture->top_poc=p_Img->toppoc;
    picture->bottom_poc=p_Img->bottompoc;
    picture->frame_poc=p_Img->framepoc;
    picture->poc=p_Img->framepoc;

    store_picture_in_dpb(p_Img, picture);

    picture=NULL;
    p_Img->pre_frame_num = UnusedShortTermFrameNum;
    UnusedShortTermFrameNum = (UnusedShortTermFrameNum + 1) % p_Img->MaxFrameNum;
  }
  p_Img->delta_pic_order_cnt[0] = tmp1;
  p_Img->delta_pic_order_cnt[1] = tmp2;
  p_Img->frame_num = CurrFrameNum;

}

/*!
 ************************************************************************
 * \brief
 *    Allocate motion parameter memory for colocated structure
 *
 ************************************************************************
 */
void alloc_motion_params(MotionParams *ftype, int size_y, int size_x)
{
  get_mem3Dint64 (&(ftype->ref_pic_id), 2, size_y, size_x);
  get_mem4Dshort (&(ftype->mv)        , 2, size_y, size_x, 2);
  get_mem3D      ((byte****)(&(ftype->ref_idx)) , 2, size_y, size_x);
  get_mem2D      (&(ftype->moving_block) , size_y, size_x);
}

/*!
 ************************************************************************
 * \brief
 *    Allocate co-located memory
 *
 * \param size_x
 *    horizontal luma size
 * \param size_y
 *    vertical luma size
 * \param mb_adaptive_frame_field_flag
 *    flag that indicates macroblock adaptive frame/field coding
 *
 * \return
 *    the allocated StorablePicture structure
 ************************************************************************
 */
ColocatedParams* alloc_colocated(int size_x, int size_y, int mb_adaptive_frame_field_flag)
{
  ColocatedParams *s;

  s = calloc(1, sizeof(ColocatedParams));
  if (NULL == s)
    no_mem_exit("alloc_colocated: s");

  s->size_x = size_x;
  s->size_y = size_y;

  alloc_motion_params(&s->frame, size_y / BLOCK_SIZE, size_x / BLOCK_SIZE);

  if (mb_adaptive_frame_field_flag)
  {
    alloc_motion_params(&s->top   , size_y / (BLOCK_SIZE * 2), size_x / BLOCK_SIZE);
    alloc_motion_params(&s->bottom, size_y / (BLOCK_SIZE * 2), size_x / BLOCK_SIZE);
  }

  s->mb_adaptive_frame_field_flag  = mb_adaptive_frame_field_flag;

  return s;
}

/*!
 ************************************************************************
 * \brief
 *    Free co-located memory.
 *
 * \param p
 *    Picture to be freed
 *
 ************************************************************************
 */
void free_colocated(ColocatedParams* p)
{
  if (p)
  {
    free_mem3D      ((byte***)p->frame.ref_idx);
    free_mem3Dint64 (p->frame.ref_pic_id);
    free_mem4Dshort (p->frame.mv);

    if (p->frame.moving_block)
    {
      free_mem2D (p->frame.moving_block);
      p->frame.moving_block=NULL;
    }

    if (p->mb_adaptive_frame_field_flag)
    {
      free_mem3D      ((byte***)p->top.ref_idx);
      free_mem3Dint64 (p->top.ref_pic_id);
      free_mem4Dshort (p->top.mv);

      if (p->top.moving_block)
      {
        free_mem2D (p->top.moving_block);
        p->top.moving_block=NULL;
      }

      free_mem3D      ((byte***)p->bottom.ref_idx);
      free_mem3Dint64 (p->bottom.ref_pic_id);
      free_mem4Dshort (p->bottom.mv);

      if (p->bottom.moving_block)
      {
        free_mem2D (p->bottom.moving_block);
        p->bottom.moving_block=NULL;
      }
    }

    free(p);

    p = NULL;
  }
}

/*!
 ************************************************************************
 * \brief
 *    Compute co-located motion info
 *
 ************************************************************************
 */

void compute_colocated (Slice *currSlice, ColocatedParams* p, StorablePicture **listX[6])
{
  StorablePicture *fs = listX[LIST_1 ][0];
  StorablePicture *fs_top = fs, *fs_bottom = fs;
  int i,j, ii, jj, jdiv;
  int fs_size_x4 = (fs->size_x >> 2);
  int fs_size_y4 = (fs->size_y >> 2);
  MotionParams *p_motion = &p->frame;
  PicMotionParams *p_frm_motion = &fs->motion;
  ImageParameters *p_Img = currSlice->p_Img;
  seq_parameter_set_rbsp_t *active_sps = p_Img->active_sps;


  if (currSlice->MbaffFrameFlag)
  {
    fs_top = listX[LIST_1 + 2][0];
    fs_bottom = listX[LIST_1 + 4][0];
  }
  else
  {
    if (p_Img->field_pic_flag)
    {
      if ((p_Img->structure != fs->structure) && (fs->coded_frame))
      {
        if (p_Img->structure==TOP_FIELD)
        {
          fs_top = fs_bottom = fs = listX[LIST_1 ][0]->top_field;
        }
        else
        {
          fs_top = fs_bottom = fs = listX[LIST_1 ][0]->bottom_field;
        }
      }
      p_frm_motion = &fs->motion;
    }
  }  

  if (!active_sps->frame_mbs_only_flag || active_sps->direct_8x8_inference_flag)
  {
    if (!currSlice->MbaffFrameFlag)
    { 
      int k;

      for (k = LIST_0; k<=LIST_1; k++)
      {
        for (j = 0; j < (fs->size_y>>2); j++)
        {
          for (i = 0 ; i < fs_size_x4 ; i++)
          {
            p_motion->mv[k][j][i][0]      = p_frm_motion->mv[k][j][i][0];
            p_motion->mv[k][j][i][1]      = p_frm_motion->mv[k][j][i][1];
            p_motion->ref_idx[k][j][i]    = p_frm_motion->ref_idx[k][j][i];
            p_motion->ref_pic_id[k][j][i] = p_frm_motion->ref_id[k][j][i];
          }
        }
      }
      p->is_long_term = fs->is_long_term;
    }
    else
    {
      for (j=0 ; j < (fs->size_y>>2); j++)
      {
        jdiv = (j>>1);
        jj = jdiv + ((j>>3)<<2);
        for (i=0 ; i < fs_size_x4 ; i++)
        {
          if (p_frm_motion->field_frame[j][i])
          {
            //! Assign frame buffers for field MBs
            //! Check whether we should use top or bottom field mvs.
            //! Depending on the assigned poc values.

            if (iabs(p_Img->dec_picture->poc - fs_bottom->poc)> iabs(p_Img->dec_picture->poc -fs_top->poc) )
            {
              p_motion->mv[LIST_0][j][i][0]      = fs_top->motion.mv[LIST_0][jdiv][i][0];
              p_motion->mv[LIST_0][j][i][1]      = fs_top->motion.mv[LIST_0][jdiv][i][1] ;
              p_motion->mv[LIST_1][j][i][0]      = fs_top->motion.mv[LIST_1][jdiv][i][0];
              p_motion->mv[LIST_1][j][i][1]      = fs_top->motion.mv[LIST_1][jdiv][i][1] ;
              p_motion->ref_idx[LIST_0][j][i]    = fs_top->motion.ref_idx[LIST_0][jdiv][i];
              p_motion->ref_idx[LIST_1][j][i]    = fs_top->motion.ref_idx[LIST_1][jdiv][i];
              p_motion->ref_pic_id[LIST_0][j][i] = p_frm_motion->ref_id[LIST_0][jj][i];
              p_motion->ref_pic_id[LIST_1][j][i] = p_frm_motion->ref_id[LIST_1][jj][i];

              p->is_long_term             = fs_top->is_long_term;
            }
            else
            {
              p_motion->mv[LIST_0][j][i][0]      = fs_bottom->motion.mv[LIST_0][jdiv][i][0];
              p_motion->mv[LIST_0][j][i][1]      = fs_bottom->motion.mv[LIST_0][jdiv][i][1] ;
              p_motion->mv[LIST_1][j][i][0]      = fs_bottom->motion.mv[LIST_1][jdiv][i][0];
              p_motion->mv[LIST_1][j][i][1]      = fs_bottom->motion.mv[LIST_1][jdiv][i][1] ;
              p_motion->ref_idx[LIST_0][j][i]    = fs_bottom->motion.ref_idx[LIST_0][jdiv][i];
              p_motion->ref_idx[LIST_1][j][i]    = fs_bottom->motion.ref_idx[LIST_1][jdiv][i];
              p_motion->ref_pic_id[LIST_0][j][i] = p_frm_motion->ref_id[LIST_0][jj + 4][i];
              p_motion->ref_pic_id[LIST_1][j][i] = p_frm_motion->ref_id[LIST_1][jj + 4][i];

              p->is_long_term             = fs_bottom->is_long_term;
            }
          }
          else
          {
            p_motion->mv[LIST_0][j][i][0]      = p_frm_motion->mv[LIST_0][j][i][0];
            p_motion->mv[LIST_0][j][i][1]      = p_frm_motion->mv[LIST_0][j][i][1] ;
            p_motion->mv[LIST_1][j][i][0]      = p_frm_motion->mv[LIST_1][j][i][0];
            p_motion->mv[LIST_1][j][i][1]      = p_frm_motion->mv[LIST_1][j][i][1] ;
            p_motion->ref_idx[LIST_0][j][i]    = p_frm_motion->ref_idx[LIST_0][j][i];
            p_motion->ref_idx[LIST_1][j][i]    = p_frm_motion->ref_idx[LIST_1][j][i];
            p_motion->ref_pic_id[LIST_0][j][i] = p_frm_motion->ref_id[LIST_0][j][i];
            p_motion->ref_pic_id[LIST_1][j][i] = p_frm_motion->ref_id[LIST_1][j][i];

            p->is_long_term             = fs->is_long_term;
          }
        }
      }
    }
  }

  //! Generate field MVs from Frame MVs
  if (p_Img->structure || currSlice->MbaffFrameFlag)
  {
    for (j = 0; j < fs->size_y >> 3; j++)
    {
      jj = RSD(j);
      for (i = 0 ; i < fs->size_x >> 2; i++)
      {
        ii = RSD(i);
        //! Do nothing if macroblock as field coded in MB-AFF
        if (!currSlice->MbaffFrameFlag )
        {
          p_motion->mv[LIST_0][j][i][0] = p_frm_motion->mv[LIST_0][jj][ii][0];
          p_motion->mv[LIST_0][j][i][1] = p_frm_motion->mv[LIST_0][jj][ii][1];
          p_motion->mv[LIST_1][j][i][0] = p_frm_motion->mv[LIST_1][jj][ii][0];
          p_motion->mv[LIST_1][j][i][1] = p_frm_motion->mv[LIST_1][jj][ii][1];

          // Scaling of references is done here since it will not affect spatial direct (2*0 =0)

          if (p_frm_motion->ref_idx[LIST_0][jj][ii] == -1)
          {
            p_motion->ref_idx   [LIST_0][j][i] = -1;
            p_motion->ref_pic_id[LIST_0][j][i] = -1;
          }
          else
          {
            p_motion->ref_idx   [LIST_0][j][i] = p_frm_motion->ref_idx[LIST_0][jj][ii] ;
            p_motion->ref_pic_id[LIST_0][j][i] = p_frm_motion->ref_id [LIST_0][jj][ii];
          }

          if (p_frm_motion->ref_idx[LIST_1][jj][ii] == -1)
          {
            p_motion->ref_idx   [LIST_1][j][i] = -1;
            p_motion->ref_pic_id[LIST_1][j][i] = -1;
          }
          else
          {
            p_motion->ref_idx   [LIST_1][j][i] = p_frm_motion->ref_idx[LIST_1][jj][ii];
            p_motion->ref_pic_id[LIST_1][j][i] = p_frm_motion->ref_id [LIST_1][jj][ii];
          }

          p->is_long_term = fs->is_long_term;

          if (currSlice->direct_spatial_mv_pred_flag == 1)
          {
            p_motion->moving_block[j][i] =
              !((!p->is_long_term
              && ((p_motion->ref_idx[LIST_0][j][i] == 0)
              &&  (iabs(p_motion->mv[LIST_0][j][i][0])>>1 == 0)
              &&  (iabs(p_motion->mv[LIST_0][j][i][1])>>1 == 0)))
              || ((p_motion->ref_idx[LIST_0][j][i] == -1)
              &&  (p_motion->ref_idx[LIST_1][j][i] == 0)
              &&  (iabs(p_motion->mv[LIST_1][j][i][0])>>1 == 0)
              &&  (iabs(p_motion->mv[LIST_1][j][i][1])>>1 == 0)));
          }
        }
        else
        {
          p->bottom.mv[LIST_0][j][i][0] = fs_bottom->motion.mv[LIST_0][jj][ii][0];
          p->bottom.mv[LIST_0][j][i][1] = fs_bottom->motion.mv[LIST_0][jj][ii][1];
          p->bottom.mv[LIST_1][j][i][0] = fs_bottom->motion.mv[LIST_1][jj][ii][0];
          p->bottom.mv[LIST_1][j][i][1] = fs_bottom->motion.mv[LIST_1][jj][ii][1];
          p->bottom.ref_idx[LIST_0][j][i] = fs_bottom->motion.ref_idx[LIST_0][jj][ii];
          p->bottom.ref_idx[LIST_1][j][i] = fs_bottom->motion.ref_idx[LIST_1][jj][ii];
          p->bottom.ref_pic_id[LIST_0][j][i] = fs_bottom->motion.ref_id[LIST_0][jj][ii];
          p->bottom.ref_pic_id[LIST_1][j][i] = fs_bottom->motion.ref_id[LIST_1][jj][ii];

          if (currSlice->direct_spatial_mv_pred_flag == 1)
          {
            p->bottom.moving_block[j][i] =
              !((!fs_bottom->is_long_term
              && ((p->bottom.ref_idx[LIST_0][j][i] == 0)
              &&  (iabs(p->bottom.mv[LIST_0][j][i][0])>>1 == 0)
              &&  (iabs(p->bottom.mv[LIST_0][j][i][1])>>1 == 0)))
              || ((p->bottom.ref_idx[LIST_0][j][i] == -1)
              &&  (p->bottom.ref_idx[LIST_1][j][i] == 0)
              &&  (iabs(p->bottom.mv[LIST_1][j][i][0])>>1 == 0)
              &&  (iabs(p->bottom.mv[LIST_1][j][i][1])>>1 == 0)));
          }

          p->top.mv[LIST_0][j][i][0] = fs_top->motion.mv[LIST_0][jj][ii][0];
          p->top.mv[LIST_0][j][i][1] = fs_top->motion.mv[LIST_0][jj][ii][1];
          p->top.mv[LIST_1][j][i][0] = fs_top->motion.mv[LIST_1][jj][ii][0];
          p->top.mv[LIST_1][j][i][1] = fs_top->motion.mv[LIST_1][jj][ii][1];
          p->top.ref_idx[LIST_0][j][i] = fs_top->motion.ref_idx[LIST_0][jj][ii];
          p->top.ref_idx[LIST_1][j][i] = fs_top->motion.ref_idx[LIST_1][jj][ii];
          p->top.ref_pic_id[LIST_0][j][i] = fs_top->motion.ref_id[LIST_0][jj][ii];
          p->top.ref_pic_id[LIST_1][j][i] = fs_top->motion.ref_id[LIST_1][jj][ii];

          if (currSlice->direct_spatial_mv_pred_flag == 1)
          {
            p->top.moving_block[j][i] =
              !((!fs_top->is_long_term
              && ((p->top.ref_idx[LIST_0][j][i] == 0)
              &&  (iabs(p->top.mv[LIST_0][j][i][0])>>1 == 0)
              &&  (iabs(p->top.mv[LIST_0][j][i][1])>>1 == 0)))
              || ((p->top.ref_idx[LIST_0][j][i] == -1)
              &&  (p->top.ref_idx[LIST_1][j][i] == 0)
              &&  (iabs(p->top.mv[LIST_1][j][i][0])>>1 == 0)
              &&  (iabs(p->top.mv[LIST_1][j][i][1])>>1 == 0)));
          }

          if ((currSlice->direct_spatial_mv_pred_flag == 0 ) && !p_frm_motion->field_frame[2*j][i])
          {
            p->top.mv[LIST_0][j][i][1] /= 2;
            p->top.mv[LIST_1][j][i][1] /= 2;
            p->bottom.mv[LIST_0][j][i][1] /= 2;
            p->bottom.mv[LIST_1][j][i][1] /= 2;
          }
        }
      }
    }
  }

  //if (!active_sps->frame_mbs_only_flag || active_sps->direct_8x8_inference_flag)
  if (!active_sps->frame_mbs_only_flag)
  {
    //! Use inference flag to remap mvs/references
    //! Frame with field co-located
    if (!p_Img->structure)
    {
      for (j=0 ; j < fs_size_y4; j++)
      {
        jdiv = (j >> 1);
        jj   = jdiv + ((j >> 3) << 2);
        for (i = 0 ; i < fs_size_x4; i++)
        {
          if (p_frm_motion->field_frame[j][i])
          {
            if (iabs(p_Img->dec_picture->poc - fs->bottom_field->poc) > iabs(p_Img->dec_picture->poc - fs->top_field->poc))
            {
              p_motion->mv[LIST_0][j][i][0] = fs->top_field->motion.mv[LIST_0][jdiv][i][0];
              p_motion->mv[LIST_0][j][i][1] = fs->top_field->motion.mv[LIST_0][jdiv][i][1] ;
              p_motion->mv[LIST_1][j][i][0] = fs->top_field->motion.mv[LIST_1][jdiv][i][0];
              p_motion->mv[LIST_1][j][i][1] = fs->top_field->motion.mv[LIST_1][jdiv][i][1] ;

              p_motion->ref_idx[LIST_0][j][i]    = fs->top_field->motion.ref_idx[LIST_0][jdiv][i];
              p_motion->ref_idx[LIST_1][j][i]    = fs->top_field->motion.ref_idx[LIST_1][jdiv][i];
              p_motion->ref_pic_id[LIST_0][j][i] = p_frm_motion->ref_id[LIST_0][jj][i];
              p_motion->ref_pic_id[LIST_1][j][i] = p_frm_motion->ref_id[LIST_1][jj][i];
              p->is_long_term             = fs->top_field->is_long_term;
            }
            else
            {
              p_motion->mv[LIST_0][j][i][0] = fs->bottom_field->motion.mv[LIST_0][jdiv][i][0];
              p_motion->mv[LIST_0][j][i][1] = fs->bottom_field->motion.mv[LIST_0][jdiv][i][1] ;
              p_motion->mv[LIST_1][j][i][0] = fs->bottom_field->motion.mv[LIST_1][jdiv][i][0];
              p_motion->mv[LIST_1][j][i][1] = fs->bottom_field->motion.mv[LIST_1][jdiv][i][1] ;

              p_motion->ref_idx[LIST_0][j][i]  = fs->bottom_field->motion.ref_idx[LIST_0][jdiv][i];
              p_motion->ref_idx[LIST_1][j][i]  = fs->bottom_field->motion.ref_idx[LIST_1][jdiv][i];
              p_motion->ref_pic_id[LIST_0][j][i] = p_frm_motion->ref_id[LIST_0][jj + 4][i];
              p_motion->ref_pic_id[LIST_1][j][i] = p_frm_motion->ref_id[LIST_1][jj + 4][i];
              p->is_long_term             = fs->bottom_field->is_long_term;
            }
          }
        }
      }
    }
  }

  p->is_long_term = fs->is_long_term;

  if (!active_sps->frame_mbs_only_flag || active_sps->direct_8x8_inference_flag)
  {
    if (currSlice->direct_spatial_mv_pred_flag == 1)
    {
      for (j=0 ; j < (fs->size_y>>2); j++)
      {
        jj = RSD(j);
        for (i=0 ; i < (fs->size_x>>2); i++)
        {
          ii = RSD(i);

          p_motion->mv[LIST_0][j][i][0]=p_motion->mv[LIST_0][jj][ii][0];
          p_motion->mv[LIST_0][j][i][1]=p_motion->mv[LIST_0][jj][ii][1];
          p_motion->mv[LIST_1][j][i][0]=p_motion->mv[LIST_1][jj][ii][0];
          p_motion->mv[LIST_1][j][i][1]=p_motion->mv[LIST_1][jj][ii][1];

          p_motion->ref_idx[LIST_0][j][i]=p_motion->ref_idx[LIST_0][jj][ii];
          p_motion->ref_idx[LIST_1][j][i]=p_motion->ref_idx[LIST_1][jj][ii];
          p_motion->ref_pic_id[LIST_0][j][i] = p_motion->ref_pic_id[LIST_0][jj][ii];
          p_motion->ref_pic_id[LIST_1][j][i] = p_motion->ref_pic_id[LIST_1][jj][ii];

          p_motion->moving_block[j][i]= (byte) (
            !((!p->is_long_term
            && ((p_motion->ref_idx[LIST_0][j][i] == 0)
            &&  (iabs(p_motion->mv[LIST_0][j][i][0])>>1 == 0)
            &&  (iabs(p_motion->mv[LIST_0][j][i][1])>>1 == 0)))
            || ((p_motion->ref_idx[LIST_0][j][i] == -1)
            &&  (p_motion->ref_idx[LIST_1][j][i] == 0)
            &&  (iabs(p_motion->mv[LIST_1][j][i][0])>>1 == 0)
            &&  (iabs(p_motion->mv[LIST_1][j][i][1])>>1 == 0))));
        }
      }
    }
    else
    {
      for (j=0 ; j < (fs->size_y>>2); j++)
      {
        jj = RSD(j);
        for (i=0 ; i < (fs->size_x>>2); i++)
        {
          ii = RSD(i);

          p_motion->mv[LIST_0][j][i][0]=p_motion->mv[LIST_0][jj][ii][0];
          p_motion->mv[LIST_0][j][i][1]=p_motion->mv[LIST_0][jj][ii][1];
          p_motion->mv[LIST_1][j][i][0]=p_motion->mv[LIST_1][jj][ii][0];
          p_motion->mv[LIST_1][j][i][1]=p_motion->mv[LIST_1][jj][ii][1];

          p_motion->ref_idx[LIST_0][j][i]=p_motion->ref_idx[LIST_0][jj][ii];
          p_motion->ref_idx[LIST_1][j][i]=p_motion->ref_idx[LIST_1][jj][ii];
          p_motion->ref_pic_id[LIST_0][j][i] = p_motion->ref_pic_id[LIST_0][jj][ii];
          p_motion->ref_pic_id[LIST_1][j][i] = p_motion->ref_pic_id[LIST_1][jj][ii];
        }
      }
    }
  }
  else
  {
    memcpy(&p_motion->mv[LIST_0][0][0][0], &p_frm_motion->mv[LIST_0][0][0][0], 4 * fs_size_y4 * fs_size_x4 * sizeof(short));
    memcpy(p_motion->ref_idx[LIST_0][0],    p_frm_motion->ref_idx[LIST_0][0],  2 * fs_size_y4 * fs_size_x4 * sizeof(char));
    memcpy(p_motion->ref_pic_id[LIST_0][0], p_frm_motion->ref_id [LIST_0][0],  2 * fs_size_y4 * fs_size_x4 * sizeof(int64));

    if (currSlice->direct_spatial_mv_pred_flag == 1)
    {
      for (j=0 ; j < fs_size_y4; j++)
      {
        for (i=0 ; i < fs_size_x4; i++)
        {
          p_motion->moving_block[j][i]=
            !((!p->is_long_term
            && ((p_motion->ref_idx[LIST_0][j][i] == 0)
            &&  (iabs(p_motion->mv[LIST_0][j][i][0])>>1 == 0)
            &&  (iabs(p_motion->mv[LIST_0][j][i][1])>>1 == 0)))
            || ((p_motion->ref_idx[LIST_0][j][i] == -1)
            &&  (p_motion->ref_idx[LIST_1][j][i] == 0)
            &&  (iabs(p_motion->mv[LIST_1][j][i][0])>>1 == 0)
            &&  (iabs(p_motion->mv[LIST_1][j][i][1])>>1 == 0)));
        }
      }
    }
  }

  if (currSlice->direct_spatial_mv_pred_flag == 0)
  {    
    if (currSlice->MbaffFrameFlag || !p_Img->structure)
    {
      for (j=0 ; j < fs_size_y4; j++)
      {
        for (i=0 ; i < fs_size_x4 ; i++)
        {
          if (p_frm_motion->field_frame[j][i])
          {
            p_motion->mv[LIST_0][j][i][1] *= 2;
            p_motion->mv[LIST_1][j][i][1] *= 2;
          }
        }
      }
    }
    else  if (p_Img->structure)
    {
      for (j=0 ; j < fs_size_y4; j++)
      {
        for (i=0 ; i < fs_size_x4 ; i++)
        {
          if (!p_frm_motion->field_frame[j][i])
          {
            p_motion->mv[LIST_0][j][i][1] /= 2;
            p_motion->mv[LIST_1][j][i][1] /= 2;
          }
        }
      }
    }

    for (j=0; j<2 + (currSlice->MbaffFrameFlag * 4);j+=2)
    {
      for (i=0; i<p_Img->listXsize[j];i++)
      {
        int prescale, iTRb, iTRp;

        if (j==0)
        {
          iTRb = iClip3( -128, 127, p_Img->dec_picture->poc - listX[LIST_0 + j][i]->poc );
        }
        else if (j == 2)
        {
          iTRb = iClip3( -128, 127, p_Img->dec_picture->top_poc - listX[LIST_0 + j][i]->poc );
        }
        else
        {
          iTRb = iClip3( -128, 127, p_Img->dec_picture->bottom_poc - listX[LIST_0 + j][i]->poc );
        }

        iTRp = iClip3( -128, 127,  listX[LIST_1 + j][0]->poc - listX[LIST_0 + j][i]->poc);

        if (iTRp!=0)
        {
          prescale = ( 16384 + iabs( iTRp / 2 ) ) / iTRp;
          currSlice->mvscale[j][i] = iClip3( -1024, 1023, ( iTRb * prescale + 32 ) >> 6 ) ;
        }
        else
        {
          currSlice->mvscale[j][i] = 9999;
        }
      }
    }
  }
}

/*!
 ************************************************************************
 * \brief
 *    Compute co-located motion info
 *    for 4:4:4 Independent mode
 *
 ************************************************************************
 */

void compute_colocated_JV(Slice *currSlice, ColocatedParams* p, StorablePicture **listX[6])
{
  ImageParameters *p_Img = currSlice->p_Img;
  seq_parameter_set_rbsp_t *active_sps = p_Img->active_sps;

  StorablePicture *fs, *fs_top, *fs_bottom;
  int i,j, ii, jj, jdiv;
  int np = p_Img->colour_plane_id;

  fs_top = fs_bottom = fs = listX[LIST_1 ][0];

  if (currSlice->MbaffFrameFlag)
  {
    fs_top= listX[LIST_1 + 2][0];
    fs_bottom= listX[LIST_1 + 4][0];
  }
  else
  {
    if (p_Img->field_pic_flag)
    {
      if ((p_Img->structure != fs->structure) && (fs->coded_frame))
      {
        if (p_Img->structure==TOP_FIELD)
        {
          fs_top=fs_bottom=fs = listX[LIST_1 ][0]->top_field;
        }
        else
        {
          fs_top=fs_bottom=fs = listX[LIST_1 ][0]->bottom_field;
        }
      }
    }
  }

  if (!active_sps->frame_mbs_only_flag || active_sps->direct_8x8_inference_flag)
  {
    for (j=0 ; j<fs->size_y/4 ; j++)
    {      
      jdiv = j/2;
      jj = j/2 + 4 * (j/8);
      for (i=0 ; i<fs->size_x/4 ; i++)
      {

        if (currSlice->MbaffFrameFlag && fs->motion.field_frame[j][i])
        {
          //! Assign frame buffers for field MBs
          //! Check whether we should use top or bottom field mvs.
          //! Depending on the assigned poc values.

          if (iabs(p_Img->dec_picture->poc - fs_bottom->poc)> iabs(p_Img->dec_picture->poc -fs_top->poc) )
          {
            p->frame.mv[LIST_0][j][i][0]      = fs_top->JVmotion[np].mv[LIST_0][jdiv][i][0];
            p->frame.mv[LIST_0][j][i][1]      = fs_top->JVmotion[np].mv[LIST_0][jdiv][i][1] ;
            p->frame.mv[LIST_1][j][i][0]      = fs_top->JVmotion[np].mv[LIST_1][jdiv][i][0];
            p->frame.mv[LIST_1][j][i][1]      = fs_top->JVmotion[np].mv[LIST_1][jdiv][i][1] ;
            p->frame.ref_idx[LIST_0][j][i]    = fs_top->JVmotion[np].ref_idx[LIST_0][jdiv][i];
            p->frame.ref_idx[LIST_1][j][i]    = fs_top->JVmotion[np].ref_idx[LIST_1][jdiv][i];
            p->frame.ref_pic_id[LIST_0][j][i] = fs->JVmotion[np].ref_id[LIST_0][jj][i];
            p->frame.ref_pic_id[LIST_1][j][i] = fs->JVmotion[np].ref_id[LIST_1][jj][i];

            p->is_long_term             = fs_top->is_long_term;
          }
          else
          {
            p->frame.mv[LIST_0][j][i][0]      = fs_bottom->JVmotion[np].mv[LIST_0][jdiv][i][0];
            p->frame.mv[LIST_0][j][i][1]      = fs_bottom->JVmotion[np].mv[LIST_0][jdiv][i][1] ;
            p->frame.mv[LIST_1][j][i][0]      = fs_bottom->JVmotion[np].mv[LIST_1][jdiv][i][0];
            p->frame.mv[LIST_1][j][i][1]      = fs_bottom->JVmotion[np].mv[LIST_1][jdiv][i][1] ;
            p->frame.ref_idx[LIST_0][j][i]    = fs_bottom->JVmotion[np].ref_idx[LIST_0][jdiv][i];
            p->frame.ref_idx[LIST_1][j][i]    = fs_bottom->JVmotion[np].ref_idx[LIST_1][jdiv][i];
            p->frame.ref_pic_id[LIST_0][j][i] = fs->JVmotion[np].ref_id[LIST_0][jj + 4][i];
            p->frame.ref_pic_id[LIST_1][j][i] = fs->JVmotion[np].ref_id[LIST_1][jj + 4][i];

            p->is_long_term             = fs_bottom->is_long_term;
          }
        }
        else
        {
          p->frame.mv[LIST_0][j][i][0]      = fs->JVmotion[np].mv[LIST_0][j][i][0];
          p->frame.mv[LIST_0][j][i][1]      = fs->JVmotion[np].mv[LIST_0][j][i][1] ;
          p->frame.mv[LIST_1][j][i][0]      = fs->JVmotion[np].mv[LIST_1][j][i][0];
          p->frame.mv[LIST_1][j][i][1]      = fs->JVmotion[np].mv[LIST_1][j][i][1] ;
          p->frame.ref_idx[LIST_0][j][i]    = fs->JVmotion[np].ref_idx[LIST_0][j][i];
          p->frame.ref_idx[LIST_1][j][i]    = fs->JVmotion[np].ref_idx[LIST_1][j][i];
          p->frame.ref_pic_id[LIST_0][j][i] = fs->JVmotion[np].ref_id[LIST_0][j][i];
          p->frame.ref_pic_id[LIST_1][j][i] = fs->JVmotion[np].ref_id[LIST_1][j][i];

          p->is_long_term             = fs->is_long_term;
        }
      }
    }
  }


  //! Generate field MVs from Frame MVs
  if (p_Img->structure || currSlice->MbaffFrameFlag)
  {
    for (j=0 ; j<fs->size_y/8 ; j++)
    {
      jj = RSD(j);
      for (i=0 ; i<fs->size_x/4 ; i++)
      {
        ii = RSD(i);
        //! Do nothing if macroblock as field coded in MB-AFF
        if (!currSlice->MbaffFrameFlag )
        {
          p->frame.mv[LIST_0][j][i][0] = fs->JVmotion[np].mv[LIST_0][jj][ii][0];
          p->frame.mv[LIST_0][j][i][1] = fs->JVmotion[np].mv[LIST_0][jj][ii][1];
          p->frame.mv[LIST_1][j][i][0] = fs->JVmotion[np].mv[LIST_1][jj][ii][0];
          p->frame.mv[LIST_1][j][i][1] = fs->JVmotion[np].mv[LIST_1][jj][ii][1];

          // Scaling of references is done here since it will not affect spatial direct (2*0 =0)

          if (fs->JVmotion[np].ref_idx[LIST_0][jj][ii] == -1)
          {
            p->frame.ref_idx   [LIST_0][j][i] = -1;
            p->frame.ref_pic_id[LIST_0][j][i] = -1;
          }
          else
          {
            p->frame.ref_idx   [LIST_0][j][i] = fs->JVmotion[np].ref_idx[LIST_0][jj][ii] ;
            p->frame.ref_pic_id[LIST_0][j][i] = fs->JVmotion[np].ref_id [LIST_0][jj][ii];
          }

          if (fs->JVmotion[np].ref_idx[LIST_1][jj][ii] == -1)
          {
            p->frame.ref_idx   [LIST_1][j][i] = -1;
            p->frame.ref_pic_id[LIST_1][j][i] = -1;
          }
          else
          {
            p->frame.ref_idx   [LIST_1][j][i] = fs->JVmotion[np].ref_idx[LIST_1][jj][ii];
            p->frame.ref_pic_id[LIST_1][j][i] = fs->JVmotion[np].ref_id [LIST_1][jj][ii];
          }

          p->is_long_term = fs->is_long_term;

          if (currSlice->direct_spatial_mv_pred_flag == 1)
          {
            p->frame.moving_block[j][i] =
              !((!p->is_long_term
              && ((p->frame.ref_idx[LIST_0][j][i] == 0)
              &&  (iabs(p->frame.mv[LIST_0][j][i][0])>>1 == 0)
              &&  (iabs(p->frame.mv[LIST_0][j][i][1])>>1 == 0)))
              || ((p->frame.ref_idx[LIST_0][j][i] == -1)
              &&  (p->frame.ref_idx[LIST_1][j][i] == 0)
              &&  (iabs(p->frame.mv[LIST_1][j][i][0])>>1 == 0)
              &&  (iabs(p->frame.mv[LIST_1][j][i][1])>>1 == 0)));
          }
        }
        else
        {
          p->bottom.mv[LIST_0][j][i][0] = fs_bottom->JVmotion[np].mv[LIST_0][jj][ii][0];
          p->bottom.mv[LIST_0][j][i][1] = fs_bottom->JVmotion[np].mv[LIST_0][jj][ii][1];
          p->bottom.mv[LIST_1][j][i][0] = fs_bottom->JVmotion[np].mv[LIST_1][jj][ii][0];
          p->bottom.mv[LIST_1][j][i][1] = fs_bottom->JVmotion[np].mv[LIST_1][jj][ii][1];
          p->bottom.ref_idx[LIST_0][j][i] = fs_bottom->JVmotion[np].ref_idx[LIST_0][jj][ii];
          p->bottom.ref_idx[LIST_1][j][i] = fs_bottom->JVmotion[np].ref_idx[LIST_1][jj][ii];
          p->bottom.ref_pic_id[LIST_0][j][i] = fs_bottom->JVmotion[np].ref_id[LIST_0][jj][ii];
          p->bottom.ref_pic_id[LIST_1][j][i] = fs_bottom->JVmotion[np].ref_id[LIST_1][jj][ii];

          if (currSlice->direct_spatial_mv_pred_flag == 1)
          {
            p->bottom.moving_block[j][i] =
              !((!fs_bottom->is_long_term
              && ((p->bottom.ref_idx[LIST_0][j][i] == 0)
              &&  (iabs(p->bottom.mv[LIST_0][j][i][0])>>1 == 0)
              &&  (iabs(p->bottom.mv[LIST_0][j][i][1])>>1 == 0)))
              || ((p->bottom.ref_idx[LIST_0][j][i] == -1)
              &&  (p->bottom.ref_idx[LIST_1][j][i] == 0)
              &&  (iabs(p->bottom.mv[LIST_1][j][i][0])>>1 == 0)
              &&  (iabs(p->bottom.mv[LIST_1][j][i][1])>>1 == 0)));
          }

          p->top.mv[LIST_0][j][i][0] = fs_top->JVmotion[np].mv[LIST_0][jj][ii][0];
          p->top.mv[LIST_0][j][i][1] = fs_top->JVmotion[np].mv[LIST_0][jj][ii][1];
          p->top.mv[LIST_1][j][i][0] = fs_top->JVmotion[np].mv[LIST_1][jj][ii][0];
          p->top.mv[LIST_1][j][i][1] = fs_top->JVmotion[np].mv[LIST_1][jj][ii][1];
          p->top.ref_idx[LIST_0][j][i] = fs_top->JVmotion[np].ref_idx[LIST_0][jj][ii];
          p->top.ref_idx[LIST_1][j][i] = fs_top->JVmotion[np].ref_idx[LIST_1][jj][ii];
          p->top.ref_pic_id[LIST_0][j][i] = fs_top->JVmotion[np].ref_id[LIST_0][jj][ii];
          p->top.ref_pic_id[LIST_1][j][i] = fs_top->JVmotion[np].ref_id[LIST_1][jj][ii];

          if (currSlice->direct_spatial_mv_pred_flag == 1)
          {
            p->top.moving_block[j][i] =
              !((!fs_top->is_long_term
              && ((p->top.ref_idx[LIST_0][j][i] == 0)
              &&  (iabs(p->top.mv[LIST_0][j][i][0])>>1 == 0)
              &&  (iabs(p->top.mv[LIST_0][j][i][1])>>1 == 0)))
              || ((p->top.ref_idx[LIST_0][j][i] == -1)
              &&  (p->top.ref_idx[LIST_1][j][i] == 0)
              &&  (iabs(p->top.mv[LIST_1][j][i][0])>>1 == 0)
              &&  (iabs(p->top.mv[LIST_1][j][i][1])>>1 == 0)));
          }

          if ((currSlice->direct_spatial_mv_pred_flag == 0 ) && !fs->motion.field_frame[2*j][i])
          {
            p->top.mv[LIST_0][j][i][1] /= 2;
            p->top.mv[LIST_1][j][i][1] /= 2;
            p->bottom.mv[LIST_0][j][i][1] /= 2;
            p->bottom.mv[LIST_1][j][i][1] /= 2;
          }

        }
      }
    }
  }


  if (!active_sps->frame_mbs_only_flag || active_sps->direct_8x8_inference_flag)
  {
    //! Use inference flag to remap mvs/references
    //! Frame with field co-located

    if (!p_Img->structure)
    {
      for (j=0 ; j<fs->size_y/4 ; j++)
      {
        jdiv = j/2;
        jj = j/2 + 4*(j/8);
        for (i=0 ; i<fs->size_x/4 ; i++)
        {

          if (fs->motion.field_frame[j][i])
          {
            if (iabs(p_Img->dec_picture->poc - fs->bottom_field->poc) > iabs(p_Img->dec_picture->poc - fs->top_field->poc))
            {
              p->frame.mv[LIST_0][j][i][0] = fs->top_field->JVmotion[np].mv[LIST_0][jdiv][i][0];
              p->frame.mv[LIST_0][j][i][1] = fs->top_field->JVmotion[np].mv[LIST_0][jdiv][i][1] ;
              p->frame.mv[LIST_1][j][i][0] = fs->top_field->JVmotion[np].mv[LIST_1][jdiv][i][0];
              p->frame.mv[LIST_1][j][i][1] = fs->top_field->JVmotion[np].mv[LIST_1][jdiv][i][1] ;

              p->frame.ref_idx[LIST_0][j][i]  = fs->top_field->JVmotion[np].ref_idx[LIST_0][jdiv][i];
              p->frame.ref_idx[LIST_1][j][i]  = fs->top_field->JVmotion[np].ref_idx[LIST_1][jdiv][i];
              p->frame.ref_pic_id[LIST_0][j][i]   = fs->JVmotion[np].ref_id[LIST_0][jj][i];
              p->frame.ref_pic_id[LIST_1][j][i]   = fs->JVmotion[np].ref_id[LIST_1][jj][i];
              p->is_long_term               = fs->top_field->is_long_term;
            }
            else
            {
              p->frame.mv[LIST_0][j][i][0] = fs->bottom_field->JVmotion[np].mv[LIST_0][jdiv][i][0];
              p->frame.mv[LIST_0][j][i][1] = fs->bottom_field->JVmotion[np].mv[LIST_0][jdiv][i][1] ;
              p->frame.mv[LIST_1][j][i][0] = fs->bottom_field->JVmotion[np].mv[LIST_1][jdiv][i][0];
              p->frame.mv[LIST_1][j][i][1] = fs->bottom_field->JVmotion[np].mv[LIST_1][jdiv][i][1] ;

              p->frame.ref_idx[LIST_0][j][i]  = fs->bottom_field->JVmotion[np].ref_idx[LIST_0][jdiv][i];
              p->frame.ref_idx[LIST_1][j][i]  = fs->bottom_field->JVmotion[np].ref_idx[LIST_1][jdiv][i];
              p->frame.ref_pic_id[LIST_0][j][i] = fs->JVmotion[np].ref_id[LIST_0][jj + 4][i];
              p->frame.ref_pic_id[LIST_1][j][i] = fs->JVmotion[np].ref_id[LIST_1][jj + 4][i];
              p->is_long_term             = fs->bottom_field->is_long_term;
            }
          }
        }
      }
    }
  }

  p->is_long_term = fs->is_long_term;

  if (!active_sps->frame_mbs_only_flag || active_sps->direct_8x8_inference_flag)
  {
    for (j=0 ; j<fs->size_y/4 ; j++)
    {
      jj = RSD(j);
      for (i=0 ; i<fs->size_x/4 ; i++)
      {
        ii = RSD(i);

        p->frame.mv[LIST_0][j][i][0] = p->frame.mv[LIST_0][jj][ii][0];
        p->frame.mv[LIST_0][j][i][1] = p->frame.mv[LIST_0][jj][ii][1];
        p->frame.mv[LIST_1][j][i][0] = p->frame.mv[LIST_1][jj][ii][0];
        p->frame.mv[LIST_1][j][i][1] = p->frame.mv[LIST_1][jj][ii][1];

        p->frame.ref_idx[LIST_0][j][i]=p->frame.ref_idx[LIST_0][jj][ii];
        p->frame.ref_idx[LIST_1][j][i]=p->frame.ref_idx[LIST_1][jj][ii];
        p->frame.ref_pic_id[LIST_0][j][i] = p->frame.ref_pic_id[LIST_0][jj][ii];
        p->frame.ref_pic_id[LIST_1][j][i] = p->frame.ref_pic_id[LIST_1][jj][ii];

        if (currSlice->direct_spatial_mv_pred_flag == 1)
        {
          p->frame.moving_block[j][i]= (byte) (
            !((!p->is_long_term
            && ((p->frame.ref_idx[LIST_0][j][i] == 0)
            &&  (iabs(p->frame.mv[LIST_0][j][i][0])>>1 == 0)
            &&  (iabs(p->frame.mv[LIST_0][j][i][1])>>1 == 0)))
            || ((p->frame.ref_idx[LIST_0][j][i] == -1)
            &&  (p->frame.ref_idx[LIST_1][j][i] == 0)
            &&  (iabs(p->frame.mv[LIST_1][j][i][0])>>1 == 0)
            &&  (iabs(p->frame.mv[LIST_1][j][i][1])>>1 == 0))));
        }
      }
    }
  }
  else
  {
    for (j=0 ; j<fs->size_y/4 ; j++)
    {
      jj = RSD(j);
      for (i=0 ; i<fs->size_x/4 ; i++)
      {
        ii = RSD(i);
        //! Use inference flag to remap mvs/references
        p->frame.mv[LIST_0][j][i][0] = fs->JVmotion[np].mv[LIST_0][j][i][0];
        p->frame.mv[LIST_0][j][i][1] = fs->JVmotion[np].mv[LIST_0][j][i][1];
        p->frame.mv[LIST_1][j][i][0] = fs->JVmotion[np].mv[LIST_1][j][i][0];
        p->frame.mv[LIST_1][j][i][1] = fs->JVmotion[np].mv[LIST_1][j][i][1];

        p->frame.ref_idx[LIST_0][j][i] = fs->JVmotion[np].ref_idx[LIST_0][j][i];
        p->frame.ref_idx[LIST_1][j][i] = fs->JVmotion[np].ref_idx[LIST_1][j][i];
        p->frame.ref_pic_id[LIST_0][j][i] = fs->JVmotion[np].ref_id[LIST_0][j][i];
        p->frame.ref_pic_id[LIST_1][j][i] = fs->JVmotion[np].ref_id[LIST_1][j][i];

        if (currSlice->direct_spatial_mv_pred_flag == 1)
        {
          p->frame.moving_block[j][i]= (byte) (
            !((!p->is_long_term
            && ((p->frame.ref_idx[LIST_0][j][i] == 0)
            &&  (iabs(p->frame.mv[LIST_0][j][i][0])>>1 == 0)
            &&  (iabs(p->frame.mv[LIST_0][j][i][1])>>1 == 0)))
            || ((p->frame.ref_idx[LIST_0][j][i] == -1)
            &&  (p->frame.ref_idx[LIST_1][j][i] == 0)
            &&  (iabs(p->frame.mv[LIST_1][j][i][0])>>1 == 0)
            &&  (iabs(p->frame.mv[LIST_1][j][i][1])>>1 == 0))));
        }
      }
    }
  }


  if (currSlice->direct_spatial_mv_pred_flag == 0)
  {
    for (j=0 ; j<fs->size_y/4 ; j++)
    {
      for (i=0 ; i<fs->size_x/4 ; i++)
      {
        if ((!currSlice->MbaffFrameFlag &&!p_Img->structure && fs->motion.field_frame[j][i]) || (currSlice->MbaffFrameFlag && fs->motion.field_frame[j][i]))
        {
          p->frame.mv[LIST_0][j][i][1] *= 2;
          p->frame.mv[LIST_1][j][i][1] *= 2;
        }
        else  if (p_Img->structure && !fs->motion.field_frame[j][i])
        {
          p->frame.mv[LIST_0][j][i][1] /= 2;
          p->frame.mv[LIST_1][j][i][1] /= 2;
        }

      }
    }

    for (j=0; j<2 + (currSlice->MbaffFrameFlag * 4);j+=2)
    {
      for (i=0; i<p_Img->listXsize[j];i++)
      {
        int prescale, iTRb, iTRp;

        if (j==0)
        {
          iTRb = iClip3( -128, 127, p_Img->dec_picture->poc - listX[LIST_0 + j][i]->poc );
        }
        else if (j == 2)
        {
          iTRb = iClip3( -128, 127, p_Img->dec_picture->top_poc - listX[LIST_0 + j][i]->poc );
        }
        else
        {
          iTRb = iClip3( -128, 127, p_Img->dec_picture->bottom_poc - listX[LIST_0 + j][i]->poc );
        }

        iTRp = iClip3( -128, 127,  listX[LIST_1 + j][0]->poc - listX[LIST_0 + j][i]->poc);

        if (iTRp!=0)
        {
          prescale = ( 16384 + iabs( iTRp / 2 ) ) / iTRp;
          currSlice->mvscale[j][i] = iClip3( -1024, 1023, ( iTRb * prescale + 32 ) >> 6 ) ;
        }
        else
        {
          currSlice->mvscale[j][i] = 9999;
        }
      }
    }
  }
}

void copy_storable_param_JV( ImageParameters *p_Img, PicMotionParams *JVplane, PicMotionParams *motion )
{
  seq_parameter_set_rbsp_t *active_sps = p_Img->active_sps;

  int md_size = (p_Img->height / BLOCK_SIZE) * (p_Img->width / BLOCK_SIZE);
  int ref_size = active_sps->frame_mbs_only_flag ? 2 * md_size : 6 * md_size;

  // copy ref_idx
  memcpy( JVplane->ref_idx[0][0], motion->ref_idx[0][0], 2 * md_size * sizeof(byte));


  // copy ref_pic_id
  memcpy( JVplane->ref_pic_id[0][0], motion->ref_pic_id[0][0], ref_size * sizeof(int64));

  // copy motion.ref_id
  memcpy( JVplane->ref_id[0][0], motion->ref_id[0][0], ref_size * sizeof(int64));

  // copy mv
  memcpy( JVplane->mv[0][0][0], motion->mv[0][0][0], 2 * md_size * 2 * sizeof(short) );
}
