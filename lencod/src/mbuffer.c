
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
 ***********************************************************************
 */

#include <limits.h>

#include "global.h"
#include "mbuffer.h"
#include "memalloc.h"
#include "output.h"
#include "image.h"
#include "nalucommon.h"
#include "img_luma.h"
#include "img_chroma.h"

extern void init_stats                   (InputParameters *p_Inp, StatParameters *stats);
static void insert_picture_in_dpb        (VideoParameters *p_Vid, FrameStore* fs, StorablePicture* p);
static void output_one_frame_from_dpb    (DecodedPictureBuffer *p_Dpb, FrameFormat *output);
static void get_smallest_poc             (DecodedPictureBuffer *p_Dpb, int *poc,int * pos);
static void gen_field_ref_ids            (StorablePicture *p);
static int  is_used_for_reference        (FrameStore* fs);
static int  remove_unused_frame_from_dpb (DecodedPictureBuffer *p_Dpb);
static int  is_short_term_reference      (FrameStore* fs);
static int  is_long_term_reference       (FrameStore* fs);


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
    if (!IS_FREXT_PROFILE(active_sps->profile_idc)&&(active_sps->constrained_set3_flag == 1))
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
  return imin( size, 16);
}

/*!
 ************************************************************************
 * \brief
 *    Check then number of frames marked "used for reference" and break
 *    if maximum is exceeded
 *
 ************************************************************************
 */
void check_num_ref(VideoParameters *p_Vid, DecodedPictureBuffer *p_Dpb)
{
  if ((int)(p_Dpb->ltref_frames_in_buffer +  p_Dpb->ref_frames_in_buffer ) > (imax(1,p_Vid->num_ref_frames)))
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
void init_dpb(VideoParameters *p_Vid, DecodedPictureBuffer *p_Dpb)
{
  unsigned i,j;

  p_Dpb->p_Vid = p_Vid;
  p_Dpb->p_Inp = p_Vid->p_Inp;

  if (p_Dpb->init_done)
  {
    free_dpb(p_Vid, p_Dpb);
  }

  p_Dpb->size = getDpbSize(p_Vid->active_sps);

  if (p_Dpb->size < (unsigned int) (p_Vid->p_Inp)->num_ref_frames)
  {
    error ("DPB size at specified level is smaller than the specified number of reference frames. This is not allowed.\n", 1000);
  }

  p_Dpb->used_size = 0;
  p_Dpb->last_picture = NULL;

  p_Dpb->ref_frames_in_buffer = 0;
  p_Dpb->ltref_frames_in_buffer = 0;

  p_Dpb->fs = calloc(p_Dpb->size, sizeof (FrameStore*));
  if (NULL==p_Dpb->fs)
    no_mem_exit("init_dpb: p_Dpb->fs");

  p_Dpb->fs_ref = calloc(p_Dpb->size, sizeof (FrameStore*));
  if (NULL==p_Dpb->fs_ref)
    no_mem_exit("init_dpb: p_Dpb->fs_ref");

  p_Dpb->fs_ltref = calloc(p_Dpb->size, sizeof (FrameStore*));
  if (NULL==p_Dpb->fs_ltref)
    no_mem_exit("init_dpb: p_Dpb->fs_ltref");

  for (i = 0; i < p_Dpb->size; i++)
  {
    p_Dpb->fs[i]       = alloc_frame_store();
    p_Dpb->fs_ref[i]   = NULL;
    p_Dpb->fs_ltref[i] = NULL;
  }

  for (i = 0; i < 6; i++)
  {
    p_Vid->listX[i] = calloc(MAX_LIST_SIZE, sizeof (StorablePicture*)); // +1 for reordering
    if (NULL==p_Vid->listX[i])
      no_mem_exit("init_dpb: p_Vid->listX[i]");
  }

  for (j = 0; j < 6; j++)
  {
    for (i = 0; i < MAX_LIST_SIZE; i++)
    {
      p_Vid->listX[j][i] = NULL;
    }
    p_Vid->listXsize[j]=0;
  }

  p_Dpb->last_output_poc = INT_MIN;

  p_Vid->last_has_mmco_5 = 0;

  p_Dpb->init_done = 1;
}


/*!
 ************************************************************************
 * \brief
 *    Free memory for decoded picture buffer.
 ************************************************************************
 */
void free_dpb(VideoParameters *p_Vid, DecodedPictureBuffer *p_Dpb)
{
  unsigned i;
  if (p_Dpb->fs)
  {
    for (i=0; i<p_Dpb->size; i++)
    {
      free_frame_store(p_Vid, p_Dpb->fs[i]);
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
  {
    if (p_Vid->listX[i])
    {
      free (p_Vid->listX[i]);
      p_Vid->listX[i] = NULL;
    }
  }

  p_Dpb->init_done = 0;
}


/*!
 ************************************************************************
 * \brief
 *    Allocate memory for decoded picture buffer frame stores and initialize with sane values.
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

void alloc_pic_motion(VideoParameters *p_Vid, PicMotionParams *motion, int size_y, int size_x)
{
  if (p_Vid->active_sps->frame_mbs_only_flag)
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
  get_mem3D((byte****)(&(motion->ref_idx)), 2, size_y, size_x);

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
 * \param p_Vid
 *    VideoParameters
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
StorablePicture* alloc_storable_picture(VideoParameters *p_Vid, PictureStructure structure, int size_x, int size_y, int size_x_cr, int size_y_cr)
{
  StorablePicture *s;
  int   nplane;
  int   dec;
  InputParameters *p_Inp = p_Vid->p_Inp;
  int   ndec = p_Inp->NoOfDecoders;

  //printf ("Allocating (%s) picture (x=%d, y=%d, x_cr=%d, y_cr=%d)\n", (type == FRAME)?"FRAME":(type == TOP_FIELD)?"TOP_FIELD":"BOTTOM_FIELD", size_x, size_y, size_x_cr, size_y_cr);

  s = calloc (1, sizeof(StorablePicture));
  if (NULL==s)
    no_mem_exit("alloc_storable_picture: s");

  s->imgUV      = NULL;
  s->imgY_sub   = NULL;
  s->imgUV_sub  = NULL;
  
  s->p_img_sub[0] = NULL;
  s->p_img_sub[1] = NULL;
  s->p_img_sub[2] = NULL;

  s->dec_imgY     = NULL;
  s->dec_imgUV    = NULL;
  s->mb_error_map = NULL;
  s->p_dec_img[0] = NULL;
  s->p_dec_img[1] = NULL;
  s->p_dec_img[2] = NULL;

  if ((p_Inp->rdopt == 3) && (ndec > 0))  //PPAHA_TODO: Need to verify whether both checks are necessary.
  {
    get_mem3D(&(s->mb_error_map), ndec, size_y/MB_BLOCK_SIZE, size_x/MB_BLOCK_SIZE);
    get_mem3Dpel(&(s->dec_imgY), ndec, size_y, size_x);

    // This seems somewhat inefficient. Why not allocate array as [ndec][x] where x goes from 0 to 2?
    if ((s->p_dec_img[0] = (imgpel***)calloc(ndec,sizeof(imgpel**))) == NULL)
    {
      no_mem_exit("mbuffer.c: p_dec_img[0]");
    }

    if (p_Vid->yuv_format != YUV400)
    {
      get_mem4Dpel(&(s->dec_imgUV), ndec, 2, size_y_cr, size_x_cr);
      if ((s->p_dec_img[1] = (imgpel***)calloc(ndec,sizeof(imgpel**))) == NULL)
      {
        no_mem_exit("mbuffer.c: p_dec_img[1]");
      }
      if ((s->p_dec_img[2] = (imgpel***)calloc(ndec,sizeof(imgpel**))) == NULL)
      {
        no_mem_exit("mbuffer.c: p_dec_img[2]");
      }
    }
  
    for (dec = 0; dec < ndec; dec++)
    {
      s->p_dec_img[0][dec] = s->dec_imgY[dec];
    }

    if (p_Vid->yuv_format != YUV400)
    {
      for (dec = 0; dec < ndec; dec++)
      {
        s->p_dec_img[1][dec] = s->dec_imgUV[dec][0];
        s->p_dec_img[2][dec] = s->dec_imgUV[dec][1];
      }
    }
  }

  get_mem2Dpel (&(s->imgY), size_y, size_x);
    
  s->p_img[0] = s->imgY;
  s->p_curr_img = s->p_img[0];

  if (p_Vid->yuv_format != YUV400)
  {
    get_mem3Dpel (&(s->imgUV), 2, size_y_cr, size_x_cr);
    s->p_img[1] = s->imgUV[0];
    s->p_img[2] = s->imgUV[1];
  }
    
  s->p_curr_img_sub = s->p_img_sub[0];

  /*
  if (p_Inp->MbInterlace)
  get_mem3Dmp    (&s->mv_info, size_y, size_x, 6);
  else
  get_mem3Dmp    (&s->mv_info, size_y, size_x, 2);
  */

  alloc_pic_motion(p_Vid, &s->motion, size_y / BLOCK_SIZE, size_x / BLOCK_SIZE);

  if( IS_INDEPENDENT(p_Inp) )
  {
    for( nplane=0; nplane<MAX_PLANE; nplane++ )
    {
      alloc_pic_motion(p_Vid, &s->JVmotion[nplane], size_y / BLOCK_SIZE, size_x / BLOCK_SIZE);
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

  s->structure=structure;

  s->size_x = size_x;
  s->size_y = size_y;
  s->size_x_padded = size_x + 2 * IMG_PAD_SIZE;
  s->size_y_padded = size_y + 2 * IMG_PAD_SIZE;
  s->size_x_pad = size_x + 2 * IMG_PAD_SIZE - 1 - MB_BLOCK_SIZE;
  s->size_y_pad = size_y + 2 * IMG_PAD_SIZE - 1 - MB_BLOCK_SIZE;
  s->size_x_cr = size_x_cr;
  s->size_y_cr = size_y_cr;
  s->size_x_cr_pad = (int) (size_x_cr - 1) + (p_Vid->pad_size_uv_x << 1) - (p_Vid->mb_cr_size_x);
  s->size_y_cr_pad = (int) (size_y_cr - 1) + (p_Vid->pad_size_uv_y << 1) - (p_Vid->mb_cr_size_y);

  s->top_field    = NULL;
  s->bottom_field = NULL;
  s->frame        = NULL;

  s->coded_frame    = 0;
  s->mb_aff_frame_flag = 0;

  init_stats(p_Inp, &s->stats);
  return s;
}

/*!
 ************************************************************************
 * \brief
 *    Free frame store memory.
 *
 * \param p_Vid
 *    VideoParameters
 * \param f
 *    FrameStore to be freed
 *
 ************************************************************************
 */
void free_frame_store(VideoParameters *p_Vid, FrameStore* f)
{
  if (f)
  {
    if (f->frame)
    {
      free_storable_picture(p_Vid, f->frame);
      f->frame=NULL;
    }
    if (f->top_field)
    {
      free_storable_picture(p_Vid, f->top_field);
      f->top_field=NULL;
    }
    if (f->bottom_field)
    {
      free_storable_picture(p_Vid, f->bottom_field);
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
    motion->mb_field=NULL;
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
 * \param p_Vid
 *    VideoParameters
 * \param p
 *    Picture to be freed
 *
 ************************************************************************
 */
void free_storable_picture(VideoParameters *p_Vid, StorablePicture* p)
{
  int nplane;
  if (p)
  {
    InputParameters *p_Inp = p_Vid->p_Inp;
    //free_mem3Dmp (p->mv_info);
    free_pic_motion(&p->motion);

    if( IS_INDEPENDENT(p_Inp) )
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

    if( IS_INDEPENDENT(p_Inp) )
    {
      if (p->imgY_sub)
      {
        free_mem4Dpel (p->imgY_sub);
        p->imgY_sub=NULL;
      }
      if (p->imgUV_sub)
      {
        free_mem5Dpel (p->imgUV_sub);
        p->imgUV_sub = NULL;
      }

      p->p_curr_img     = NULL;
      p->p_curr_img_sub = NULL;
      p->p_img[0]       = NULL;
      p->p_img[1]       = NULL;
      p->p_img[2]       = NULL;
      p->p_img_sub[0]   = NULL;
      p->p_img_sub[1]   = NULL;
      p->p_img_sub[2]   = NULL;
    }
    else
    {
      if (p->imgY_sub)
      {
        free_mem4Dpel (p->imgY_sub);
        p->imgY_sub=NULL;
      }

      if ( p->imgUV_sub && p_Vid->yuv_format != YUV400 && p_Inp->ChromaMCBuffer )
      {
        free_mem5Dpel (p->imgUV_sub);
        p->imgUV_sub = NULL;
      }
    }

    if (p->imgUV)
    {
      free_mem3Dpel (p->imgUV);
      p->imgUV=NULL;
    }

    if (p->dec_imgY)
    {
      free_mem3Dpel(p->dec_imgY);
    }
    if (p->dec_imgUV)
    {
      free_mem4Dpel(p->dec_imgUV);
    }
    for (nplane = 0; nplane < 3; nplane++)
    {
      if (p->p_dec_img[nplane])
      {  
        free(p->p_dec_img[nplane]);
        p->p_dec_img[nplane] = NULL;
      }
    }
    if (p->mb_error_map)
    {
      free_mem3D(p->mb_error_map);
    }

    p->mb_error_map = NULL;
    p->dec_imgY   = NULL;
    p->dec_imgUV  = NULL;

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
    if (fs->frame->imgY_sub)
    {
      free_mem4Dpel (fs->frame->imgY_sub);
      fs->frame->imgY_sub=NULL;
    }

    if (fs->frame->imgUV_sub)
    {
      free_mem5Dpel (fs->frame->imgUV_sub);
      fs->frame->imgUV_sub = NULL;
    }

    free_pic_motion(&fs->frame->motion);
  }

  if (fs->top_field)
  {
    if (fs->top_field->imgY_sub)
    {
      free_mem4Dpel (fs->top_field->imgY_sub);
      fs->top_field->imgY_sub=NULL;
    }

    if (fs->top_field->imgUV_sub)
    {
      free_mem5Dpel (fs->top_field->imgUV_sub);
      fs->top_field->imgUV_sub = NULL;
    }
    
    free_pic_motion(&fs->top_field->motion);
  }
  if (fs->bottom_field)
  {
    if (fs->bottom_field->imgY_sub)
    {
      free_mem4Dpel (fs->bottom_field->imgY_sub);
      fs->bottom_field->imgY_sub=NULL;
    }
    if (fs->bottom_field->imgUV_sub)
    {
      free_mem5Dpel (fs->bottom_field->imgUV_sub);
      fs->bottom_field->imgUV_sub = NULL;
    }

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
  if ( pic_num1 < pic_num2)
    return 1;
  if ( pic_num1 > pic_num2)
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

  if ( poc1 < poc2)
    return 1;
  if ( poc1 > poc2)
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
 *    Initialize p_Vid->listX[0] and list 1 depending on current picture type
 *
 ************************************************************************
 */
void init_lists(Slice *currSlice)
{
  VideoParameters *p_Vid = currSlice->p_Vid;
  DecodedPictureBuffer *p_Dpb = p_Vid->p_Dpb;

  int add_top = 0, add_bottom = 0;
  unsigned int i;
  int j, diff;

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
          if( p_Dpb->fs_ref[i]->frame_num > currSlice->frame_num )
          {
            p_Dpb->fs_ref[i]->frame_num_wrap = p_Dpb->fs_ref[i]->frame_num - currSlice->max_frame_num;
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
    for (i=0; i<p_Dpb->ltref_frames_in_buffer; i++)
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
        if( p_Dpb->fs_ref[i]->frame_num > currSlice->frame_num )
        {
          p_Dpb->fs_ref[i]->frame_num_wrap = p_Dpb->fs_ref[i]->frame_num - currSlice->max_frame_num;
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

  if ((currSlice->slice_type == I_SLICE)||(p_Vid->type == SI_SLICE))
  {
    p_Vid->listXsize[0] = 0;
    p_Vid->listXsize[1] = 0;
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
            p_Vid->listX[0][list0idx++] = p_Dpb->fs_ref[i]->frame;
          }
        }
      }
      // order list 0 by PicNum
      qsort((void *)p_Vid->listX[0], list0idx, sizeof(StorablePicture*), compare_pic_by_pic_num_desc);
      p_Vid->listXsize[0] = (char) list0idx;
      //printf("p_Vid->listX[0] (PicNum): "); for (i=0; i<list0idx; i++){printf ("%d  ", p_Vid->listX[0][i]->pic_num);} printf("\n");

      // long term handling
      for (i=0; i<p_Dpb->ltref_frames_in_buffer; i++)
      {
        if (p_Dpb->fs_ltref[i]->is_used==3)
        {
          if (p_Dpb->fs_ltref[i]->frame->is_long_term)
          {
            p_Vid->listX[0][list0idx++]=p_Dpb->fs_ltref[i]->frame;            
          }
        }
      }
      qsort((void *)&p_Vid->listX[0][(short) p_Vid->listXsize[0]], list0idx - p_Vid->listXsize[0], sizeof(StorablePicture*), compare_pic_by_lt_pic_num_asc);
      p_Vid->listXsize[0] = (char) list0idx;

      //printf("p_Vid->listX[0] currPoc=%d (Poc): ", currSlice->framepoc); for (i=0; i<(unsigned int) p_Vid->listXsize[0]; i++){printf ("%d  ", p_Vid->listX[0][i]->poc);} printf("\n");
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

      //printf("fs_list0 (FrameNum): "); for (i=0; i<list0idx; i++){printf ("%d  ", fs_list0[i]->frame_num_wrap);} printf("\n");

      p_Vid->listXsize[0] = 0;
      gen_pic_list_from_frame_list(currSlice->structure, fs_list0, list0idx, p_Vid->listX[0], &p_Vid->listXsize[0], 0);

      //printf("p_Vid->listX[0] (PicNum): "); for (i=0; i < p_Vid->listXsize[0]; i++){printf ("%d  ", p_Vid->listX[0][i]->pic_num);} printf("\n");

      // long term handling
      for (i=0; i<p_Dpb->ltref_frames_in_buffer; i++)
      {
        fs_listlt[listltidx++]=p_Dpb->fs_ltref[i];
      }

      qsort((void *)fs_listlt, listltidx, sizeof(FrameStore*), compare_fs_by_lt_pic_idx_asc);

      gen_pic_list_from_frame_list(currSlice->structure, fs_listlt, listltidx, p_Vid->listX[0], &p_Vid->listXsize[0], 1);

      free(fs_list0);
      free(fs_listlt);
    }
    p_Vid->listXsize[1] = 0;
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
            if (currSlice->framepoc > p_Dpb->fs_ref[i]->frame->poc)
            {
              p_Vid->listX[0][list0idx++] = p_Dpb->fs_ref[i]->frame;
            }
          }
        }
      }
      qsort((void *)p_Vid->listX[0], list0idx, sizeof(StorablePicture*), compare_pic_by_poc_desc);

      //get the backward reference picture (POC>current POC) in list0;
      list0idx_1 = list0idx;
      for (i=0; i<p_Dpb->ref_frames_in_buffer; i++)
      {
        if (p_Dpb->fs_ref[i]->is_used==3)
        {
          if ((p_Dpb->fs_ref[i]->frame->used_for_reference)&&(!p_Dpb->fs_ref[i]->frame->is_long_term))
          {
            if (currSlice->framepoc < p_Dpb->fs_ref[i]->frame->poc)
            {
              p_Vid->listX[0][list0idx++] = p_Dpb->fs_ref[i]->frame;
            }
          }
        }
      }
      qsort((void *)&p_Vid->listX[0][list0idx_1], list0idx-list0idx_1, sizeof(StorablePicture*), compare_pic_by_poc_asc);

      for (j=0; j<list0idx_1; j++)
      {
        p_Vid->listX[1][list0idx-list0idx_1+j]=p_Vid->listX[0][j];
      }
      for (j=list0idx_1; j<list0idx; j++)
      {
        p_Vid->listX[1][j-list0idx_1]=p_Vid->listX[0][j];
      }

      p_Vid->listXsize[0] = p_Vid->listXsize[1] = (char) list0idx;

//      printf("p_Vid->listX[0] currPoc=%d (Poc): ", currSlice->framepoc); for (i=0; i<p_Vid->listXsize[0]; i++){printf ("%d  ", p_Vid->listX[0][i]->poc);} printf("\n");
//      printf("p_Vid->listX[1] currPoc=%d (Poc): ", currSlice->framepoc); for (i=0; i<p_Vid->listXsize[1]; i++){printf ("%d  ", p_Vid->listX[1][i]->poc);} printf("\n");

      // long term handling
      for (i=0; i<p_Dpb->ltref_frames_in_buffer; i++)
      {
        if (p_Dpb->fs_ltref[i]->is_used==3)
        {
          if (p_Dpb->fs_ltref[i]->frame->is_long_term)
          {
            p_Vid->listX[0][list0idx]   = p_Dpb->fs_ltref[i]->frame;
            p_Vid->listX[1][list0idx++] = p_Dpb->fs_ltref[i]->frame;
          }
        }
      }
      qsort((void *)&p_Vid->listX[0][(short) p_Vid->listXsize[0]], list0idx - p_Vid->listXsize[0], sizeof(StorablePicture*), compare_pic_by_lt_pic_num_asc);
      qsort((void *)&p_Vid->listX[1][(short) p_Vid->listXsize[0]], list0idx - p_Vid->listXsize[0], sizeof(StorablePicture*), compare_pic_by_lt_pic_num_asc);
      p_Vid->listXsize[0] = p_Vid->listXsize[1] = (char) list0idx;
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

      p_Vid->listXsize[0] = 0;
      p_Vid->listXsize[1] = 1;

      for (i=0; i<p_Dpb->ref_frames_in_buffer; i++)
      {
        if (p_Dpb->fs_ref[i]->is_used)
        {
          if (currSlice->ThisPOC >= p_Dpb->fs_ref[i]->poc)
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
          if (currSlice->ThisPOC < p_Dpb->fs_ref[i]->poc)
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

//      printf("fs_list0 currPoc=%d (Poc): ", currSlice->ThisPOC); for (i=0; i<list0idx; i++){printf ("%d  ", fs_list0[i]->poc);} printf("\n");
//      printf("fs_list1 currPoc=%d (Poc): ", currSlice->ThisPOC); for (i=0; i<list0idx; i++){printf ("%d  ", fs_list1[i]->poc);} printf("\n");

      p_Vid->listXsize[0] = 0;
      p_Vid->listXsize[1] = 0;
      gen_pic_list_from_frame_list(currSlice->structure, fs_list0, list0idx, p_Vid->listX[0], &p_Vid->listXsize[0], 0);
      gen_pic_list_from_frame_list(currSlice->structure, fs_list1, list0idx, p_Vid->listX[1], &p_Vid->listXsize[1], 0);

//      printf("p_Vid->listX[0] currPoc=%d (Poc): ", currSlice->framepoc); for (i=0; i<p_Vid->listXsize[0]; i++){printf ("%d  ", p_Vid->listX[0][i]->poc);} printf("\n");
//      printf("p_Vid->listX[1] currPoc=%d (Poc): ", currSlice->framepoc); for (i=0; i<p_Vid->listXsize[1]; i++){printf ("%d  ", p_Vid->listX[1][i]->poc);} printf("\n");

      // long term handling
      for (i=0; i<p_Dpb->ltref_frames_in_buffer; i++)
      {
        fs_listlt[listltidx++]=p_Dpb->fs_ltref[i];
      }

      qsort((void *)fs_listlt, listltidx, sizeof(FrameStore*), compare_fs_by_lt_pic_idx_asc);

      gen_pic_list_from_frame_list(currSlice->structure, fs_listlt, listltidx, p_Vid->listX[0], &p_Vid->listXsize[0], 1);
      gen_pic_list_from_frame_list(currSlice->structure, fs_listlt, listltidx, p_Vid->listX[1], &p_Vid->listXsize[1], 1);

      free(fs_list0);
      free(fs_list1);
      free(fs_listlt);
    }
  }

  if ((p_Vid->listXsize[0] == p_Vid->listXsize[1]) && (p_Vid->listXsize[0] > 1))
  {
    // check if lists are identical, if yes swap first two elements of p_Vid->listX[1]
    diff=0;
    for (j = 0; j< p_Vid->listXsize[0]; j++)
    {
      if (p_Vid->listX[0][j]!=p_Vid->listX[1][j])
        diff=1;
    }
    if (!diff)
    {
      tmp_s = p_Vid->listX[1][0];
      p_Vid->listX[1][0]=p_Vid->listX[1][1];
      p_Vid->listX[1][1]=tmp_s;
    }
  }

  // set max size
  p_Vid->listXsize[0] = (char) imin (p_Vid->listXsize[0], currSlice->num_ref_idx_active[LIST_0]);
  p_Vid->listXsize[1] = (char) imin (p_Vid->listXsize[1], currSlice->num_ref_idx_active[LIST_1]);

  // set the unused list entries to NULL
  for (i=p_Vid->listXsize[0]; i< (MAX_LIST_SIZE) ; i++)
  {
    p_Vid->listX[0][i] = NULL;
  }
  for (i=p_Vid->listXsize[1]; i< (MAX_LIST_SIZE) ; i++)
  {
    p_Vid->listX[1][i] = NULL;
  }
}

/*!
 ************************************************************************
 * \brief
 *    Initialize p_Vid->listX[2..5] from lists 0 and 1
 *    p_Vid->listX[2]: list0 for current_field==top
 *    p_Vid->listX[3]: list1 for current_field==top
 *    p_Vid->listX[4]: list0 for current_field==bottom
 *    p_Vid->listX[5]: list1 for current_field==bottom
 *
 ************************************************************************
 */
void init_mbaff_lists(Slice *currSlice)
{
  // for the time being listX is part of p_Vid
  VideoParameters *p_Vid = currSlice->p_Vid;
  unsigned j;
  int i;

  for (i=2;i<6;i++)
  {
    for (j=0; j<MAX_LIST_SIZE; j++)
    {
      p_Vid->listX[i][j] = NULL;
    }
    p_Vid->listXsize[i]=0;
  }

  for (i=0; i<p_Vid->listXsize[0]; i++)
  {
    p_Vid->listX[2][2*i]  =p_Vid->listX[0][i]->top_field;
    p_Vid->listX[2][2*i+1]=p_Vid->listX[0][i]->bottom_field;
    p_Vid->listX[4][2*i]  =p_Vid->listX[0][i]->bottom_field;
    p_Vid->listX[4][2*i+1]=p_Vid->listX[0][i]->top_field;
  }
  p_Vid->listXsize[2]=p_Vid->listXsize[4]=p_Vid->listXsize[0] * 2;

  for (i=0; i<p_Vid->listXsize[1]; i++)
  {
    p_Vid->listX[3][2*i]  =p_Vid->listX[1][i]->top_field;
    p_Vid->listX[3][2*i+1]=p_Vid->listX[1][i]->bottom_field;
    p_Vid->listX[5][2*i]  =p_Vid->listX[1][i]->bottom_field;
    p_Vid->listX[5][2*i+1]=p_Vid->listX[1][i]->top_field;
  }
  p_Vid->listXsize[3]=p_Vid->listXsize[5]=p_Vid->listXsize[1] * 2;
}

 /*!
 ************************************************************************
 * \brief
 *    Returns short term pic with given picNum
 *
 ************************************************************************
 */
static StorablePicture*  get_short_term_pic(Slice *currSlice, DecodedPictureBuffer *p_Dpb, int picNum)
{
  unsigned i;

  for (i=0; i < p_Dpb->ref_frames_in_buffer; i++)
  {
    if (currSlice->structure == FRAME)
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
  return NULL;
}

/*!
 ************************************************************************
 * \brief
 *    Returns short term pic with given LongtermPicNum
 *
 ************************************************************************
 */
static StorablePicture*  get_long_term_pic(Slice *currSlice, DecodedPictureBuffer *p_Dpb, int LongtermPicNum)
{
  unsigned i;

  for (i=0; i<p_Dpb->ltref_frames_in_buffer; i++)
  {
    if (currSlice->structure==FRAME)
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
static void reorder_short_term(Slice *currSlice, DecodedPictureBuffer *p_Dpb, StorablePicture **RefPicListX, int cur_list, int picNumLX, int *refIdxLX)
{
  int cIdx, nIdx;

  StorablePicture *picLX;

  picLX = get_short_term_pic(currSlice, p_Dpb, picNumLX);

  for( cIdx = currSlice->num_ref_idx_active[cur_list]; cIdx > *refIdxLX; cIdx-- )
    RefPicListX[ cIdx ] = RefPicListX[ cIdx - 1];

  RefPicListX[ (*refIdxLX)++ ] = picLX;

  nIdx = *refIdxLX;

  for( cIdx = *refIdxLX; cIdx <= currSlice->num_ref_idx_active[cur_list]; cIdx++ )
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
static void reorder_long_term(Slice *currSlice, DecodedPictureBuffer *p_Dpb, StorablePicture **RefPicListX, int cur_list, int frame_no, int *refIdxLX)
{
  int cIdx, nIdx;
  int LongTermPicNum = currSlice->long_term_pic_idx[cur_list][frame_no];

  StorablePicture *picLX;

  picLX = get_long_term_pic(currSlice, p_Dpb, LongTermPicNum);

  for( cIdx = currSlice->num_ref_idx_active[cur_list]; cIdx > *refIdxLX; cIdx-- )
    RefPicListX[ cIdx ] = RefPicListX[ cIdx - 1];

  RefPicListX[ (*refIdxLX)++ ] = picLX;

  nIdx = *refIdxLX;

  for( cIdx = *refIdxLX; cIdx <= currSlice->num_ref_idx_active[cur_list]; cIdx++ )
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
void reorder_ref_pic_list(Slice *currSlice, StorablePicture **list[6], char list_size[6], int cur_list)
{
  int i;

  int maxPicNum, currPicNum, picNumLXNoWrap, picNumLXPred, picNumLX;
  int refIdxLX = 0;
  int *reordering_of_pic_nums_idc  = currSlice->reordering_of_pic_nums_idc[cur_list];
  int *abs_diff_pic_num_minus1 = currSlice->abs_diff_pic_num_minus1[cur_list];
  VideoParameters *p_Vid = currSlice->p_Vid;
  DecodedPictureBuffer *p_Dpb = p_Vid->p_Dpb;

  if (currSlice->structure==FRAME)
  {
    maxPicNum  = currSlice->max_frame_num;
    currPicNum = currSlice->frame_num;
  }
  else
  {
    maxPicNum  = 2 * currSlice->max_frame_num;
    currPicNum = 2 * currSlice->frame_num + 1;
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
      else // (reordering_of_pic_nums_idc[i] == 1)
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

      reorder_short_term(currSlice, p_Dpb, list[cur_list], cur_list, picNumLX, &refIdxLX);
    }
    else //(reordering_of_pic_nums_idc[i] == 2)
    {
      reorder_long_term (currSlice, p_Dpb, list[cur_list], cur_list, i, &refIdxLX);
    }
  }

  // that's a definition
  list_size[cur_list] = currSlice->num_ref_idx_active[cur_list];
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
      p_Dpb->fs_ltref[j++] = p_Dpb->fs[i];
    }
  }

  p_Dpb->ltref_frames_in_buffer = j;

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
static void idr_memory_management(DecodedPictureBuffer *p_Dpb, StorablePicture* p, FrameFormat *output)
{
  unsigned i;
  VideoParameters *p_Vid = p_Dpb->p_Vid;

  assert (p_Vid->currentPicture->idr_flag);

  if (p_Vid->no_output_of_prior_pics_flag)
  {
    // free all stored pictures
    for (i=0; i<p_Dpb->used_size; i++)
    {
      // reset all reference settings
      free_frame_store(p_Vid, p_Dpb->fs[i]);
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
    flush_dpb(p_Vid, output);
  }
  p_Dpb->last_picture = NULL;

  update_ref_list(p_Dpb);
  update_ltref_list(p_Dpb);
  p_Dpb->last_output_poc = INT_MIN;

  if (p_Vid->long_term_reference_flag)
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
static void sliding_window_memory_management(VideoParameters *p_Vid, DecodedPictureBuffer *p_Dpb, StorablePicture* p)
{
  unsigned i;

  // if this is a reference pic with sliding sliding window, unmark first ref frame
  if (p_Dpb->ref_frames_in_buffer==p_Vid->active_sps->num_ref_frames - p_Dpb->ltref_frames_in_buffer)
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
static void unmark_long_term_field_for_reference_by_frame_idx(VideoParameters *p_Vid, DecodedPictureBuffer *p_Dpb, PictureStructure structure, int long_term_frame_idx, int mark_current, unsigned curr_frame_num, int curr_pic_num)
{
  unsigned i;

  assert(structure!=FRAME);
  if (curr_pic_num<0)
    curr_pic_num += (2 * p_Vid->max_frame_num);

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
              if ((p_Dpb->fs_ltref[i]->frame_num) != (unsigned)(curr_pic_num >> 1))
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
              if ((p_Dpb->fs_ltref[i]->frame_num) != (unsigned)(curr_pic_num >> 1))
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
static void mm_assign_long_term_frame_idx(VideoParameters *p_Vid, StorablePicture* p, int difference_of_pic_nums_minus1, int long_term_frame_idx)
{
  DecodedPictureBuffer *p_Dpb = p_Vid->p_Dpb;
  int picNumX = get_pic_num_x(p, difference_of_pic_nums_minus1);

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

    unmark_long_term_field_for_reference_by_frame_idx(p_Vid, p_Dpb, structure, long_term_frame_idx, 0, 0, picNumX);
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
static void mm_mark_current_picture_long_term(VideoParameters *p_Vid, DecodedPictureBuffer *p_Dpb, StorablePicture *p, int long_term_frame_idx)
{
  // remove long term pictures with same long_term_frame_idx
  if (p->structure == FRAME)
  {
    unmark_long_term_frame_for_reference_by_frame_idx(p_Dpb, long_term_frame_idx);
  }
  else
  {
    unmark_long_term_field_for_reference_by_frame_idx(p_Vid, p_Dpb, p->structure, long_term_frame_idx, 1, p->pic_num, 0);
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
static void adaptive_memory_management(DecodedPictureBuffer *p_Dpb, StorablePicture* p, FrameFormat *output)
{
  DecRefPicMarking_t *tmp_drpm;
  VideoParameters *p_Vid = p_Dpb->p_Vid;

  p_Vid->last_has_mmco_5 = 0;

  assert (!p_Vid->currentPicture->idr_flag);
  assert (p_Vid->adaptive_ref_pic_buffering_flag);

  while (p_Vid->dec_ref_pic_marking_buffer)
  {
    tmp_drpm = p_Vid->dec_ref_pic_marking_buffer;
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
        mm_assign_long_term_frame_idx(p_Vid, p, tmp_drpm->difference_of_pic_nums_minus1, tmp_drpm->long_term_frame_idx);
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
        p_Vid->last_has_mmco_5 = 1;
        break;
      case 6:
        mm_mark_current_picture_long_term(p_Vid, p_Dpb, p, tmp_drpm->long_term_frame_idx);
        check_num_ref(p_Vid, p_Dpb);
        break;
      default:
        error ("invalid memory_management_control_operation in buffer", 500);
    }
    p_Vid->dec_ref_pic_marking_buffer = tmp_drpm->Next;
    free (tmp_drpm);
  }
  if ( p_Vid->last_has_mmco_5 )
  {
    p->pic_num = p->frame_num = 0;

    switch (p->structure)
    {
    case TOP_FIELD:
      {
        p->poc = p->top_poc = p_Vid->toppoc =0;
        break;
      }
    case BOTTOM_FIELD:
      {
        p->poc = p->bottom_poc = p_Vid->bottompoc = 0;
        break;
      }
    case FRAME:
      {
        p->top_poc    -= p->poc;
        p->bottom_poc -= p->poc;

        p_Vid->toppoc = p->top_poc;
        p_Vid->bottompoc = p->bottom_poc;

        p->poc = imin (p->top_poc, p->bottom_poc);
        p_Vid->framepoc = p->poc;
        break;
      }
    }
    p_Vid->ThisPOC = p->poc;
    flush_dpb(p_Vid, output);
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
 * \param p_Vid
 *    VideoParameters
 * \param p
 *    Picture to be stored
 * \param output
 *    FrameFormat for output
 *
 ************************************************************************
 */
void store_picture_in_dpb(VideoParameters *p_Vid, StorablePicture* p, FrameFormat *output)
{
  unsigned i;
  int poc, pos;
  DecodedPictureBuffer *p_Dpb = p_Vid->p_Dpb;
  // diagnostics
  //printf ("Storing (%s) non-ref pic with frame_num #%d\n", (p->type == FRAME)?"FRAME":(p->type == TOP_FIELD)?"TOP_FIELD":"BOTTOM_FIELD", p->pic_num);
  // if frame, check for new store,
  assert (p!=NULL);

  p->used_for_reference = (p_Vid->nal_reference_idc != NALU_PRIORITY_DISPOSABLE);
  p->type = p_Vid->type;

  p_Vid->last_has_mmco_5=0;
  p_Vid->last_pic_bottom_field = (p_Vid->structure == BOTTOM_FIELD);

  if (p_Vid->currentPicture->idr_flag)
    idr_memory_management(p_Dpb, p, output);
  else
  {
    // adaptive memory management
    if (p->used_for_reference && (p_Vid->adaptive_ref_pic_buffering_flag))
      adaptive_memory_management(p_Dpb, p, output);
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
            insert_picture_in_dpb(p_Vid, p_Dpb->last_picture, p);
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
  if ((!p_Vid->currentPicture->idr_flag)&&(p->used_for_reference && (!p_Vid->adaptive_ref_pic_buffering_flag)))
  {
    sliding_window_memory_management(p_Vid, p_Dpb, p);
  }

  // first try to remove unused frames
  if (p_Dpb->used_size==p_Dpb->size)
  {
    remove_unused_frame_from_dpb(p_Dpb);
  }

  // then output frames until one can be removed
  while (p_Dpb->used_size==p_Dpb->size)
  {
    // non-reference frames may be output directly
    if (!p->used_for_reference)
    {
      get_smallest_poc(p_Dpb, &poc, &pos);
      if ((-1==pos) || (p->poc < poc))
      {
        direct_output(p_Vid, p, output, p_Vid->p_dec);
        return;
      }
    }
    // flush a frame
    output_one_frame_from_dpb(p_Dpb, output);
  }

  // check for duplicate frame number in short term reference buffer
  if ((p->used_for_reference)&&(!p->is_long_term))
  {
    for (i=0; i<p_Dpb->ref_frames_in_buffer; i++)
    {
      if (p_Dpb->fs_ref[i]->frame_num == p->frame_num)
      {
        error("duplicate frame_num im short-term reference picture buffer", 500);
      }
    }

  }
  // store at end of buffer
  insert_picture_in_dpb(p_Vid, p_Dpb->fs[p_Dpb->used_size],p);

  if (p->structure != FRAME)
  {
    p_Dpb->last_picture = p_Dpb->fs[p_Dpb->used_size];
  }
  else
  {
    p_Dpb->last_picture = NULL;
  }

  p_Dpb->used_size++;

  update_ref_list(p_Dpb);
  update_ltref_list(p_Dpb);

  check_num_ref(p_Vid, p_Dpb);

  dump_dpb(p_Dpb);
}


/*!
 ************************************************************************
 * \brief
 *    Insert the frame picture into the if the top field has already
 *    been stored for the coding decision
 *
 * \param p_Vid
 *    VideoParameters
 * \param p
 *    StorablePicture to be inserted
 * \param output
 *    FrameFormat for output
 *
 ************************************************************************
 */
void replace_top_pic_with_frame(VideoParameters *p_Vid, StorablePicture* p, FrameFormat *output)
{
  DecodedPictureBuffer *p_Dpb = p_Vid->p_Dpb;
  FrameStore* fs = NULL;
  unsigned i, found;

  assert (p!=NULL);
  assert (p->structure==FRAME);

  p->used_for_reference = (p_Vid->nal_reference_idc != NALU_PRIORITY_DISPOSABLE);
  p->type = p_Vid->type;
  // upsample a reference picture
  if (p->used_for_reference)
  {
    if( IS_INDEPENDENT(p_Vid->p_Inp) )
    {
      UnifiedOneForthPix_JV(p_Vid, 0, p);
      UnifiedOneForthPix_JV(p_Vid, 1, p);
      UnifiedOneForthPix_JV(p_Vid, 2, p);
    }
    else
    {
      UnifiedOneForthPix(p_Vid, p);
    }
  }

  found=0;

  for (i = 0; i < p_Dpb->used_size; i++)
  {
    if((p_Dpb->fs[i]->frame_num == p_Vid->frame_num)&&(p_Dpb->fs[i]->is_used==1))
    {
      found=1;
      fs = p_Dpb->fs[i];
      break;
    }
  }

  if (!found)
  {
    // this should only happen for non-reference pictures when the dpb is full of reference pics
    direct_output_paff(p_Vid, p, output, p_Vid->p_dec);
  }
  else
  {
    free_storable_picture(p_Vid, fs->top_field);
    fs->top_field=NULL;
    fs->frame=p;
    fs->is_used = 3;
    if (p->used_for_reference)
    {
      fs->is_reference = 3;
      if (p->is_long_term)
      {
        fs->is_long_term = 3;
      }
    }
    // generate field views
    dpb_split_field(p_Vid, fs);
    update_ref_list(p_Dpb);
    update_ltref_list(p_Dpb);
  }
}


/*!
 ************************************************************************
 * \brief
 *    Insert the picture into the DPB. A free DPB position is necessary
 *    for frames, .
 *
 * \param p_Vid
 *    VideoParameters
 * \param fs
 *    FrameStore into which the picture will be inserted
 * \param p
 *    StorablePicture to be inserted
 *
 ************************************************************************
 */
static void insert_picture_in_dpb(VideoParameters *p_Vid, FrameStore* fs, StorablePicture* p)
{
  //  printf ("insert (%s) pic with frame_num #%d, poc %d\n", (p->structure == FRAME)?"FRAME":(p->structure == TOP_FIELD)?"TOP_FIELD":"BOTTOM_FIELD", p->pic_num, p->poc);
  assert (p!=NULL);
  assert (fs!=NULL);

  // upsample a reference picture
  if (p->used_for_reference)
  {    
    if( IS_INDEPENDENT(p_Vid->p_Inp) )
    {
      UnifiedOneForthPix_JV(p_Vid, 0, p);
      UnifiedOneForthPix_JV(p_Vid, 1, p);
      UnifiedOneForthPix_JV(p_Vid, 2, p);
    }
    else
    {
      UnifiedOneForthPix(p_Vid, p);
    }
  }

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
    dpb_split_field(p_Vid, fs);
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
      dpb_combine_field(p_Vid, fs);
    }
    else
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
      dpb_combine_field(p_Vid, fs);
    }
    else
    {
      fs->poc = p->poc;
      gen_field_ref_ids(p);
    }
    break;
  }
  fs->frame_num = p->pic_num;
  fs->is_output = p->is_output;

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
static void remove_frame_from_dpb(DecodedPictureBuffer *p_Dpb, int pos)
{  
  VideoParameters *p_Vid = p_Dpb->p_Vid;
  FrameStore* fs = p_Dpb->fs[pos];
  FrameStore* tmp;
  unsigned i;

//  printf ("remove frame with frame_num #%d\n", fs->frame_num);
  switch (fs->is_used)
  {
  case 3:
    free_storable_picture(p_Vid, fs->frame);
    free_storable_picture(p_Vid, fs->top_field);
    free_storable_picture(p_Vid, fs->bottom_field);
    fs->frame=NULL;
    fs->top_field=NULL;
    fs->bottom_field=NULL;
    break;
  case 2:
    free_storable_picture(p_Vid, fs->bottom_field);
    fs->bottom_field=NULL;
    break;
  case 1:
    free_storable_picture(p_Vid, fs->top_field);
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
  for (i=0; i < p_Dpb->used_size; i++)
  {
    if ((*poc>p_Dpb->fs[i]->poc)&&(!p_Dpb->fs[i]->is_output))
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
static int remove_unused_frame_from_dpb(DecodedPictureBuffer *p_Dpb)
{
  unsigned i;

  // check for frames that were already output and no longer used for reference
  for (i=0; i<p_Dpb->used_size; i++)
  {
    if (p_Dpb->fs[i]->is_output && (!is_used_for_reference(p_Dpb->fs[i])))
    {
      remove_frame_from_dpb(p_Dpb, i);
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
static void output_one_frame_from_dpb(DecodedPictureBuffer *p_Dpb, FrameFormat *output)
{
  VideoParameters *p_Vid = p_Dpb->p_Vid;
  InputParameters *p_Inp = p_Dpb->p_Inp;
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
  //  printf ("output frame with frame_num #%d, poc %d (p_Dpb-> p_Dpb->size=%d, p_Dpb->used_size=%d)\n", p_Dpb->fs[pos]->frame_num, p_Dpb->fs[pos]->frame->poc, p_Dpb->size, p_Dpb->used_size);
  write_stored_frame(p_Vid, p_Dpb->fs[pos], output, p_Vid->p_dec);

  // if redundant picture in use, output POC may be not in ascending order
  if(p_Inp->redundant_pic_flag == 0)
  {
    if (p_Dpb->last_output_poc >= poc)
    {
      error ("output POC must be in ascending order", 150);
    }
  }
  p_Dpb->last_output_poc = poc;

  // free frame store and move empty store to end of buffer
  if (!is_used_for_reference(p_Dpb->fs[pos]))
  {
    remove_frame_from_dpb(p_Dpb, pos);
  }
}



/*!
 ************************************************************************
 * \brief
 *    All stored picture are output. Should be called to empty the buffer
 ************************************************************************
 */
void flush_dpb(VideoParameters *p_Vid, FrameFormat *output)
{
  DecodedPictureBuffer *p_Dpb = p_Vid->p_Dpb;
  unsigned i;

  //diagnostics
//  printf("Flush remaining frames from p_Dpb. p_Dpb->size=%d, p_Dpb->used_size=%d\n",p_Dpb->size,p_Dpb->used_size);

  // mark all frames unused
  for (i=0; i<p_Dpb->used_size; i++)
  {
    unmark_for_reference (p_Dpb->fs[i]);
  }

  while (remove_unused_frame_from_dpb(p_Dpb)) ;

  // output frames in POC order
  while (p_Dpb->used_size)
  {
    output_one_frame_from_dpb(p_Dpb, output);
  }

  p_Dpb->last_output_poc = INT_MIN;
}


static void gen_field_ref_ids(StorablePicture *p)
{
  int i,j, dummylist0, dummylist1;
   //! Generate Frame parameters from field information.
  for (i=0 ; i<p->size_x >> 2 ; i++)
  {
    for (j=0 ; j<p->size_y >> 2 ; j++)
    {
        dummylist0= p->motion.ref_idx[LIST_0][j][i];
        dummylist1= p->motion.ref_idx[LIST_1][j][i];
        //! association with id already known for fields.
        p->motion.ref_id[LIST_0][j][i] = (dummylist0>=0)? p->ref_pic_num[LIST_0][dummylist0] : 0;
        p->motion.ref_id[LIST_1][j][i] = (dummylist1>=0)? p->ref_pic_num[LIST_1][dummylist1] : 0;
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
void dpb_split_field(VideoParameters *p_Vid, FrameStore *fs)
{
  int i, j, k, ii, jj, jj4;
  int idiv,jdiv;
  int currentmb;
  int dummylist0,dummylist1;
  int twosz16 = 2*(fs->frame->size_x>>4);

  fs->poc = fs->frame->poc;

  if (!p_Vid->active_sps->frame_mbs_only_flag)
  {
    fs->top_field    = alloc_storable_picture(p_Vid, TOP_FIELD,    fs->frame->size_x, fs->frame->size_y >> 1, fs->frame->size_x_cr, fs->frame->size_y_cr >> 1);
    fs->bottom_field = alloc_storable_picture(p_Vid, BOTTOM_FIELD, fs->frame->size_x, fs->frame->size_y >> 1, fs->frame->size_x_cr, fs->frame->size_y_cr >> 1);

    for (i=0; i<fs->frame->size_y >> 1; i++)
    {
      memcpy(fs->top_field->imgY[i], fs->frame->imgY[i*2], fs->frame->size_x*sizeof(imgpel));
      memcpy(fs->bottom_field->imgY[i], fs->frame->imgY[i*2 + 1], fs->frame->size_x*sizeof(imgpel));
    }

    for (k = 0; k < 2; k++)
    {
      for (i=0; i<fs->frame->size_y_cr >> 1; i++)
      {
        memcpy(fs->top_field->imgUV[k][i], fs->frame->imgUV[k][i*2], fs->frame->size_x_cr*sizeof(imgpel));
        memcpy(fs->bottom_field->imgUV[k][i], fs->frame->imgUV[k][i*2 + 1], fs->frame->size_x_cr*sizeof(imgpel));
      }
    }

    if( IS_INDEPENDENT(p_Vid->p_Inp) )
    {
      UnifiedOneForthPix_JV(p_Vid, 0, fs->top_field);
      UnifiedOneForthPix_JV(p_Vid, 1, fs->top_field);
      UnifiedOneForthPix_JV(p_Vid, 2, fs->top_field);
      UnifiedOneForthPix_JV(p_Vid, 0, fs->bottom_field);
      UnifiedOneForthPix_JV(p_Vid, 1, fs->bottom_field);
      UnifiedOneForthPix_JV(p_Vid, 2, fs->bottom_field);
    }
    else
    {
      UnifiedOneForthPix(p_Vid, fs->top_field);
      UnifiedOneForthPix(p_Vid, fs->bottom_field);
    }

    fs->top_field->poc = fs->frame->top_poc;
    fs->bottom_field->poc =  fs->frame->bottom_poc;

    fs->top_field->frame_poc =  fs->frame->frame_poc;

    fs->top_field->bottom_poc =fs->bottom_field->bottom_poc =  fs->frame->bottom_poc;
    fs->top_field->top_poc =fs->bottom_field->top_poc =  fs->frame->top_poc;
    fs->bottom_field->frame_poc =  fs->frame->frame_poc;

    fs->top_field->used_for_reference = fs->bottom_field->used_for_reference
                                      = fs->frame->used_for_reference;
    fs->top_field->is_long_term = fs->bottom_field->is_long_term
                                = fs->frame->is_long_term;
    fs->long_term_frame_idx = fs->top_field->long_term_frame_idx
                            = fs->bottom_field->long_term_frame_idx
                            = fs->frame->long_term_frame_idx;

    fs->top_field->coded_frame = fs->bottom_field->coded_frame = 1;
    fs->top_field->mb_aff_frame_flag = fs->bottom_field->mb_aff_frame_flag
                                  = fs->frame->mb_aff_frame_flag;

    fs->frame->top_field    = fs->top_field;
    fs->frame->bottom_field = fs->bottom_field;

    fs->top_field->bottom_field = fs->bottom_field;
    fs->top_field->frame        = fs->frame;
    fs->bottom_field->top_field = fs->top_field;
    fs->bottom_field->frame     = fs->frame;

    fs->top_field->chroma_format_idc = fs->bottom_field->chroma_format_idc = fs->frame->chroma_format_idc;
    fs->top_field->chroma_mask_mv_x  = fs->bottom_field->chroma_mask_mv_x  = fs->frame->chroma_mask_mv_x;
    fs->top_field->chroma_mask_mv_y  = fs->bottom_field->chroma_mask_mv_y  = fs->frame->chroma_mask_mv_y;
    fs->top_field->chroma_shift_x    = fs->bottom_field->chroma_shift_x    = fs->frame->chroma_shift_x;
    fs->top_field->chroma_shift_y    = fs->bottom_field->chroma_shift_y    = fs->frame->chroma_shift_y;

    //store reference picture index
    memcpy(fs->top_field->ref_pic_num[LIST_1]   , fs->frame->ref_pic_num[2 + LIST_1], 2 * p_Vid->listXsize[LIST_1] * sizeof(int64));
    memcpy(fs->bottom_field->ref_pic_num[LIST_1], fs->frame->ref_pic_num[4 + LIST_1], 2 * p_Vid->listXsize[LIST_1] * sizeof(int64));
    memcpy(fs->top_field->ref_pic_num[LIST_0]   , fs->frame->ref_pic_num[2 + LIST_0], 2 * p_Vid->listXsize[LIST_0] * sizeof(int64));
    memcpy(fs->bottom_field->ref_pic_num[LIST_0], fs->frame->ref_pic_num[4 + LIST_0], 2 * p_Vid->listXsize[LIST_0] * sizeof(int64));

  }
  else
  {
    fs->top_field=NULL;
    fs->bottom_field=NULL;
    fs->frame->top_field=NULL;
    fs->frame->bottom_field=NULL;
  }

  if (!fs->frame->mb_aff_frame_flag)
  {
    for (j=0 ; j<fs->frame->size_y >> 2 ; j++)
    {
      for (i=0 ; i<fs->frame->size_x >> 2 ; i++)
      {
        dummylist0 = fs->frame->motion.ref_idx[LIST_0][j][i];
        dummylist1 = fs->frame->motion.ref_idx[LIST_1][j][i];
        fs->frame->motion.ref_id[LIST_0][j][i] = (dummylist0>=0)? fs->frame->ref_pic_num[LIST_0][dummylist0] : -1;
        fs->frame->motion.ref_id[LIST_1][j][i] = (dummylist1>=0)? fs->frame->ref_pic_num[LIST_1][dummylist1] : -1;
      }
    }
  }
  else
  {
    for (j = 0 ; j < (fs->frame->size_y >> 2); j++)
    {
      jdiv = (j >> 2);
      for (i=0 ; i < (fs->frame->size_x >> 2); i++)
      {
        idiv = (i >> 2);
        currentmb = twosz16*(jdiv >> 1)+ (idiv)*2 + (jdiv & 0x01);

        if (fs->frame->motion.mb_field[currentmb])
        {
          int list_offset = (currentmb & 0x01)? 4: 2;
          dummylist0 = fs->frame->motion.ref_idx[LIST_0][j][i];
          dummylist1 = fs->frame->motion.ref_idx[LIST_1][j][i];
          //! association with id already known for fields.
          fs->frame->motion.ref_id[LIST_0 + list_offset][j][i] = (dummylist0>=0)? fs->frame->ref_pic_num[LIST_0 + list_offset][dummylist0] : 0;
          fs->frame->motion.ref_id[LIST_1 + list_offset][j][i] = (dummylist1>=0)? fs->frame->ref_pic_num[LIST_1 + list_offset][dummylist1] : 0;
          //! need to make association with frames
          fs->frame->motion.ref_id[LIST_0][j][i] = (dummylist0>=0)? fs->frame->frm_ref_pic_num[LIST_0 + list_offset][dummylist0] : 0;
          fs->frame->motion.ref_id[LIST_1][j][i] = (dummylist1>=0)? fs->frame->frm_ref_pic_num[LIST_1 + list_offset][dummylist1] : 0;

        }
        else
        {
          dummylist0 = fs->frame->motion.ref_idx[LIST_0][j][i];
          dummylist1 = fs->frame->motion.ref_idx[LIST_1][j][i];
          fs->frame->motion.ref_id[LIST_0][j][i] = (dummylist0>=0)? fs->frame->ref_pic_num[LIST_0][dummylist0] : -1;
          fs->frame->motion.ref_id[LIST_1][j][i] = (dummylist1>=0)? fs->frame->ref_pic_num[LIST_1][dummylist1] : -1;
        }
      }
    }
  }

  if (!p_Vid->active_sps->frame_mbs_only_flag && fs->frame->mb_aff_frame_flag)
  {
    for (j=0 ; j<fs->frame->size_y >> 3; j++)
    {
      jj = (j >> 2)*8 + (j & 0x03);
      jj4 = jj + 4;
      jdiv = (j >> 1);
      for (i=0 ; i < (fs->frame->size_x >> 2); i++)
      {
        idiv=i >> 2;

        currentmb = twosz16*(jdiv >> 1)+ (idiv)*2 + (jdiv & 0x01);
        // Assign field mvs attached to MB-Frame buffer to the proper buffer
        if (fs->frame->motion.mb_field[currentmb])
        {
          fs->bottom_field->motion.field_frame[j][i] = fs->top_field->motion.field_frame[j][i]=1;
          fs->frame->motion.field_frame[2*j][i] = fs->frame->motion.field_frame[2*j+1][i]=1;

          fs->bottom_field->motion.mv[LIST_0][j][i][0] = fs->frame->motion.mv[LIST_0][jj4][i][0];
          fs->bottom_field->motion.mv[LIST_0][j][i][1] = fs->frame->motion.mv[LIST_0][jj4][i][1];
          fs->bottom_field->motion.mv[LIST_1][j][i][0] = fs->frame->motion.mv[LIST_1][jj4][i][0];
          fs->bottom_field->motion.mv[LIST_1][j][i][1] = fs->frame->motion.mv[LIST_1][jj4][i][1];
          fs->bottom_field->motion.ref_idx[LIST_0][j][i] = fs->frame->motion.ref_idx[LIST_0][jj4][i];
          fs->bottom_field->motion.ref_idx[LIST_1][j][i] = fs->frame->motion.ref_idx[LIST_1][jj4][i];
          fs->bottom_field->motion.ref_id[LIST_0][j][i] = fs->frame->motion.ref_id[LIST_0+4][jj4][i];
          fs->bottom_field->motion.ref_id[LIST_1][j][i] = fs->frame->motion.ref_id[LIST_1+4][jj4][i];


          fs->top_field->motion.mv[LIST_0][j][i][0] = fs->frame->motion.mv[LIST_0][jj][i][0];
          fs->top_field->motion.mv[LIST_0][j][i][1] = fs->frame->motion.mv[LIST_0][jj][i][1];
          fs->top_field->motion.mv[LIST_1][j][i][0] = fs->frame->motion.mv[LIST_1][jj][i][0];
          fs->top_field->motion.mv[LIST_1][j][i][1] = fs->frame->motion.mv[LIST_1][jj][i][1];
          fs->top_field->motion.ref_idx[LIST_0][j][i] = fs->frame->motion.ref_idx[LIST_0][jj][i];
          fs->top_field->motion.ref_idx[LIST_1][j][i] = fs->frame->motion.ref_idx[LIST_1][jj][i];
          fs->top_field->motion.ref_id[LIST_0][j][i] = fs->frame->motion.ref_id[LIST_0+2][jj][i];
          fs->top_field->motion.ref_id[LIST_1][j][i] = fs->frame->motion.ref_id[LIST_1+2][jj][i];
        }
      }
    }
  }

  //! Generate field MVs from Frame MVs
  if (!p_Vid->active_sps->frame_mbs_only_flag)
  {
    for (j=0 ; j<fs->frame->size_y >> 3 ; j++)
    {
      jj = 2* RSD(j);
      jdiv = j >> 1;
      for (i=0 ; i<fs->frame->size_x >> 2 ; i++)
      {
        ii = RSD(i);
        idiv = i >> 2;

        currentmb = twosz16*(jdiv >> 1)+ (idiv)*2 + (jdiv & 0x01);

        if (!fs->frame->mb_aff_frame_flag  || !fs->frame->motion.mb_field[currentmb])
        {
          fs->frame->motion.field_frame[2*j+1][i] = fs->frame->motion.field_frame[2*j][i]=0;

          fs->top_field->motion.field_frame[j][i] = fs->bottom_field->motion.field_frame[j][i] = 0;

          fs->top_field->motion.mv[LIST_0][j][i][0] = fs->bottom_field->motion.mv[LIST_0][j][i][0] = fs->frame->motion.mv[LIST_0][jj][ii][0];
          fs->top_field->motion.mv[LIST_0][j][i][1] = fs->bottom_field->motion.mv[LIST_0][j][i][1] = fs->frame->motion.mv[LIST_0][jj][ii][1];
          fs->top_field->motion.mv[LIST_1][j][i][0] = fs->bottom_field->motion.mv[LIST_1][j][i][0] = fs->frame->motion.mv[LIST_1][jj][ii][0];
          fs->top_field->motion.mv[LIST_1][j][i][1] = fs->bottom_field->motion.mv[LIST_1][j][i][1] = fs->frame->motion.mv[LIST_1][jj][ii][1];

          // Scaling of references is done here since it will not affect spatial direct (2*0 =0)
          if (fs->frame->motion.ref_idx[LIST_0][jj][ii] == -1)
            fs->top_field->motion.ref_idx[LIST_0][j][i] = fs->bottom_field->motion.ref_idx[LIST_0][j][i] = - 1;
          else
          {
            dummylist0=fs->top_field->motion.ref_idx[LIST_0][j][i] = fs->bottom_field->motion.ref_idx[LIST_0][j][i] = fs->frame->motion.ref_idx[LIST_0][jj][ii];
            fs->top_field   ->motion.ref_id[LIST_0][j][i] = (dummylist0>=0)? fs->frame->top_ref_pic_num[LIST_0][dummylist0] : 0;
            fs->bottom_field->motion.ref_id[LIST_0][j][i] = (dummylist0>=0)? fs->frame->bottom_ref_pic_num[LIST_0][dummylist0] : 0;
          }

          if (fs->frame->motion.ref_idx[LIST_1][jj][ii] == -1)
            fs->top_field->motion.ref_idx[LIST_1][j][i] = fs->bottom_field->motion.ref_idx[LIST_1][j][i] = - 1;
          else
          {
            dummylist1=fs->top_field->motion.ref_idx[LIST_1][j][i] = fs->bottom_field->motion.ref_idx[LIST_1][j][i] = fs->frame->motion.ref_idx[LIST_1][jj][ii];

            fs->top_field   ->motion.ref_id[LIST_1][j][i] = (dummylist1>=0)? fs->frame->top_ref_pic_num[LIST_1][dummylist1] : 0;
            fs->bottom_field->motion.ref_id[LIST_1][j][i] = (dummylist1>=0)? fs->frame->bottom_ref_pic_num[LIST_1][dummylist1] : 0;
          }
        }
        else
        {
          fs->frame->motion.field_frame[2*j+1][i] = fs->frame->motion.field_frame[2*j][i]= fs->frame->motion.mb_field[currentmb];
        }
      }
    }
  }
  else
  {
    memset( &(fs->frame->motion.field_frame[0][0]), 0, fs->frame->size_y * (fs->frame->size_x >>4) * sizeof(byte));
  }
}


/*!
 ************************************************************************
 * \brief
 *    Generate a frame from top and bottom fields,
 *    YUV components and display information only
 ************************************************************************
 */
void dpb_combine_field_yuv(VideoParameters *p_Vid, FrameStore *fs)
{
  int i, j;

  fs->frame = alloc_storable_picture(p_Vid, FRAME, fs->top_field->size_x, fs->top_field->size_y*2, fs->top_field->size_x_cr, fs->top_field->size_y_cr*2);

  for (i=0; i<fs->top_field->size_y; i++)
  {
    memcpy(fs->frame->imgY[i*2],     fs->top_field->imgY[i]   , fs->top_field->size_x*sizeof(imgpel));     // top field
    memcpy(fs->frame->imgY[i*2 + 1], fs->bottom_field->imgY[i], fs->bottom_field->size_x*sizeof(imgpel)); // bottom field
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
  fs->frame->chroma_mask_mv_x  = fs->top_field->chroma_mask_mv_x;
  fs->frame->chroma_mask_mv_y  = fs->top_field->chroma_mask_mv_y;
  fs->frame->chroma_shift_x    = fs->top_field->chroma_shift_x;
  fs->frame->chroma_shift_y    = fs->top_field->chroma_shift_y;

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
void dpb_combine_field(VideoParameters *p_Vid, FrameStore *fs)
{
  int i,j, jj, jj4;
  int dummylist0, dummylist1;

  dpb_combine_field_yuv(p_Vid, fs);

  if( IS_INDEPENDENT(p_Vid->p_Inp) )
  {
    UnifiedOneForthPix_JV(p_Vid, 0, fs->frame);
    UnifiedOneForthPix_JV(p_Vid, 1, fs->frame);
    UnifiedOneForthPix_JV(p_Vid, 2, fs->frame);
  }
  else
  {
    UnifiedOneForthPix(p_Vid, fs->frame);
  }

  //combine field for frame
  for (i=0;i<(p_Vid->listXsize[LIST_1]+1) >> 1;i++)
  {
    fs->frame->ref_pic_num[LIST_1][i]= i64min ((fs->top_field->ref_pic_num[LIST_1][2*i] >> 1)*2, (fs->bottom_field->ref_pic_num[LIST_1][2*i] >> 1)*2);
  }

  for (i=0;i<(p_Vid->listXsize[LIST_0]+1) >> 1;i++)
  {
    fs->frame->ref_pic_num[LIST_0][i]= i64min ((fs->top_field->ref_pic_num[LIST_0][2*i] >> 1)*2, (fs->bottom_field->ref_pic_num[LIST_0][2*i] >> 1)*2);
  }

   //! Use inference flag to remap mvs/references

  //! Generate Frame parameters from field information.
  for (j=0 ; j<fs->top_field->size_y >> 2 ; j++)
  {
    jj = 8*(j >> 2) + (j & 0x03);
    jj4 = jj + 4;
    for (i=0 ; i<fs->top_field->size_x >> 2 ; i++)
    {
      fs->frame->motion.field_frame[jj][i]= fs->frame->motion.field_frame[jj4][i]=1;

      fs->frame->motion.mv[LIST_0][jj][i][0] = fs->top_field->motion.mv[LIST_0][j][i][0];
      fs->frame->motion.mv[LIST_0][jj][i][1] = fs->top_field->motion.mv[LIST_0][j][i][1] ;
      fs->frame->motion.mv[LIST_1][jj][i][0] = fs->top_field->motion.mv[LIST_1][j][i][0];
      fs->frame->motion.mv[LIST_1][jj][i][1] = fs->top_field->motion.mv[LIST_1][j][i][1] ;

      dummylist0=fs->frame->motion.ref_idx[LIST_0][jj][i]  = fs->top_field->motion.ref_idx[LIST_0][j][i];
      dummylist1=fs->frame->motion.ref_idx[LIST_1][jj][i]  = fs->top_field->motion.ref_idx[LIST_1][j][i];

      //! association with id already known for fields.
      fs->top_field->motion.ref_id[LIST_0][j][i] = (dummylist0>=0)? fs->top_field->ref_pic_num[LIST_0][dummylist0] : 0;
      fs->top_field->motion.ref_id[LIST_1][j][i] = (dummylist1>=0)? fs->top_field->ref_pic_num[LIST_1][dummylist1] : 0;

      //! need to make association with frames
      fs->frame->motion.ref_id[LIST_0][jj][i] = (dummylist0>=0)? fs->top_field->frm_ref_pic_num[LIST_0][dummylist0] : 0;
      fs->frame->motion.ref_id[LIST_1][jj][i] = (dummylist1>=0)? fs->top_field->frm_ref_pic_num[LIST_1][dummylist1] : 0;

      fs->frame->motion.mv[LIST_0][jj4][i][0] = fs->bottom_field->motion.mv[LIST_0][j][i][0];
      fs->frame->motion.mv[LIST_0][jj4][i][1] = fs->bottom_field->motion.mv[LIST_0][j][i][1] ;
      fs->frame->motion.mv[LIST_1][jj4][i][0] = fs->bottom_field->motion.mv[LIST_1][j][i][0];
      fs->frame->motion.mv[LIST_1][jj4][i][1] = fs->bottom_field->motion.mv[LIST_1][j][i][1] ;

      dummylist0=fs->frame->motion.ref_idx[LIST_0][jj4][i]  = fs->bottom_field->motion.ref_idx[LIST_0][j][i];
      dummylist1=fs->frame->motion.ref_idx[LIST_1][jj4][i]  = fs->bottom_field->motion.ref_idx[LIST_1][j][i];

      fs->bottom_field->motion.ref_id[LIST_0][j][i] = (dummylist0>=0)? fs->bottom_field->ref_pic_num[LIST_0][dummylist0] : 0;
      fs->bottom_field->motion.ref_id[LIST_1][j][i] = (dummylist1>=0)? fs->bottom_field->ref_pic_num[LIST_1][dummylist1] : 0;

      //! need to make association with frames
      fs->frame->motion.ref_id[LIST_0][jj4][i] = (dummylist0>=0)? fs->bottom_field->frm_ref_pic_num[LIST_0][dummylist0] : -1;
      fs->frame->motion.ref_id[LIST_1][jj4][i] = (dummylist1>=0)? fs->bottom_field->frm_ref_pic_num[LIST_1][dummylist1] : -1;

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
  int size = currSlice->num_ref_idx_active[LIST_0] + 1;

  if (currSlice->slice_type!=I_SLICE && currSlice->slice_type != SI_SLICE)
  {
    if ((currSlice->reordering_of_pic_nums_idc[LIST_0] = calloc(size, sizeof(int)))==NULL) no_mem_exit("alloc_ref_pic_list_reordering_buffer: remapping_of_pic_nums_idc_l0");
    if ((currSlice->abs_diff_pic_num_minus1[LIST_0] = calloc(size, sizeof(int)))==NULL) no_mem_exit("alloc_ref_pic_list_reordering_buffer: abs_diff_pic_num_minus1_l0");
    if ((currSlice->long_term_pic_idx[LIST_0] = calloc(size, sizeof(int)))==NULL) no_mem_exit("alloc_ref_pic_list_reordering_buffer: long_term_pic_idx_l0");
  }
  else
  {
    currSlice->reordering_of_pic_nums_idc[LIST_0] = NULL;
    currSlice->abs_diff_pic_num_minus1[LIST_0] = NULL;
    currSlice->long_term_pic_idx[LIST_0] = NULL;
  }

  size = currSlice->num_ref_idx_active[LIST_1] + 1;

  if (currSlice->slice_type == B_SLICE)
  {
    if ((currSlice->reordering_of_pic_nums_idc[LIST_1] = calloc(size, sizeof(int)))==NULL) no_mem_exit("alloc_ref_pic_list_reordering_buffer: remapping_of_pic_nums_idc_l1");
    if ((currSlice->abs_diff_pic_num_minus1[LIST_1] = calloc(size, sizeof(int)))==NULL) no_mem_exit("alloc_ref_pic_list_reordering_buffer: abs_diff_pic_num_minus1_l1");
    if ((currSlice->long_term_pic_idx[LIST_1] = calloc(size, sizeof(int)))==NULL) no_mem_exit("alloc_ref_pic_list_reordering_buffer: long_term_pic_idx_l1");
  }
  else
  {
    currSlice->reordering_of_pic_nums_idc[LIST_1] = NULL;
    currSlice->abs_diff_pic_num_minus1[LIST_1] = NULL;
    currSlice->long_term_pic_idx[LIST_1] = NULL;
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

  if (currSlice->reordering_of_pic_nums_idc[LIST_0])
    free(currSlice->reordering_of_pic_nums_idc[LIST_0]);
  if (currSlice->abs_diff_pic_num_minus1[LIST_0])
    free(currSlice->abs_diff_pic_num_minus1[LIST_0]);
  if (currSlice->long_term_pic_idx[LIST_0])
    free(currSlice->long_term_pic_idx[LIST_0]);

  currSlice->reordering_of_pic_nums_idc[LIST_0] = NULL;
  currSlice->abs_diff_pic_num_minus1[LIST_0] = NULL;
  currSlice->long_term_pic_idx[LIST_0] = NULL;

  if (currSlice->reordering_of_pic_nums_idc[LIST_1])
    free(currSlice->reordering_of_pic_nums_idc[LIST_1]);
  if (currSlice->abs_diff_pic_num_minus1[LIST_1])
    free(currSlice->abs_diff_pic_num_minus1[LIST_1]);
  if (currSlice->long_term_pic_idx[LIST_1])
    free(currSlice->long_term_pic_idx[LIST_1]);

  currSlice->reordering_of_pic_nums_idc[LIST_1] = NULL;
  currSlice->abs_diff_pic_num_minus1[LIST_1] = NULL;
  currSlice->long_term_pic_idx[LIST_1] = NULL;
}

/*!
 ************************************************************************
 * \brief
 *      Tian Dong
 *          June 13, 2002, Modifed on July 30, 2003
 *
 *      If a gap in frame_num is found, try to fill the gap
 * \param p_Vid
 *    VideoParameters structure
 * \param output
 *    FrameFormat for ouput
 *
 ************************************************************************
 */
void fill_frame_num_gap(VideoParameters *p_Vid, FrameFormat *output)
{
  int CurrFrameNum;
  int UnusedShortTermFrameNum;
  StorablePicture *picture = NULL;
  int nal_ref_idc_bak;

//  printf("A gap in frame number is found, try to fill it.\n");

  nal_ref_idc_bak = p_Vid->nal_reference_idc;
  p_Vid->nal_reference_idc = NALU_PRIORITY_LOW;

  UnusedShortTermFrameNum = (p_Vid->pre_frame_num + 1) % p_Vid->max_frame_num;
  CurrFrameNum = p_Vid->frame_num;

  while (CurrFrameNum != UnusedShortTermFrameNum)
  {
    picture = alloc_storable_picture (p_Vid, FRAME, p_Vid->width, p_Vid->height, p_Vid->width_cr, p_Vid->height_cr);
    picture->coded_frame = 1;
    picture->pic_num = UnusedShortTermFrameNum;
    picture->non_existing = 1;
    picture->is_output = 1;

    p_Vid->adaptive_ref_pic_buffering_flag = 0;

    store_picture_in_dpb(p_Vid, picture, output);

    picture=NULL;
    UnusedShortTermFrameNum = (UnusedShortTermFrameNum + 1) % p_Vid->max_frame_num;
  }

  p_Vid->nal_reference_idc = nal_ref_idc_bak;
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
 *    Free motion parameter memory for colocated structure
 *
 ************************************************************************
 */
void free_motion_params(MotionParams *ftype)
{
  free_mem3Dint64 (ftype->ref_pic_id);
  free_mem3D      ((byte***)ftype->ref_idx);
  free_mem4Dshort (ftype->mv);

  if (ftype->moving_block)
  {
    free_mem2D (ftype->moving_block);
    ftype->moving_block=NULL;
  }
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
    free_motion_params(&p->frame);

    if (p->mb_adaptive_frame_field_flag)
    {
      free_motion_params(&p->top   );
      free_motion_params(&p->bottom);
    }

    free(p);

    p=NULL;
  }
}

/*!
 ************************************************************************
 * \brief
 *    Compute co-located motion info
 *
 ************************************************************************
 */

void compute_colocated(Slice *currSlice, ColocatedParams* p, StorablePicture **listX[6])
{
  VideoParameters *p_Vid = currSlice->p_Vid;
  StorablePicture *fs, *fs_top, *fs_bottom;
  MotionParams *ftype = NULL;
  int i,j, ii, jj, jdiv;

  fs_top = fs_bottom = fs = listX[LIST_1 ][0];

  if (currSlice->mb_aff_frame_flag)
  {
    fs_top= listX[LIST_1 + 2][0];
    fs_bottom= listX[LIST_1 + 4][0];
  }
  else
  {
    if (currSlice->structure!=FRAME)
    {
      if ((currSlice->structure != fs->structure) && (fs->coded_frame))
      {
        if (currSlice->structure==TOP_FIELD)
        {
          fs_top = fs_bottom = fs = listX[LIST_1 ][0]->top_field;
        }
        else
        {
          fs_top = fs_bottom = fs = listX[LIST_1 ][0]->bottom_field;
        }
      }
    }
  }

  if (!p_Vid->active_sps->frame_mbs_only_flag || p_Vid->active_sps->direct_8x8_inference_flag)
  {
    ftype = &p->frame;
    for (j=0 ; j<fs->size_y >> 2 ; j++)
    {
      jdiv = j>>1;
      jj = (j >> 1) + 4 * (j >> 3);
      for (i=0 ; i<fs->size_x >> 2 ; i++)
      {

        if (currSlice->mb_aff_frame_flag && fs->motion.field_frame[j][i])
        {
          //! Assign frame buffers for field MBs
          //! Check whether we should use top or bottom field mvs.
          //! Depending on the assigned poc values.

          if (iabs(p_Vid->enc_picture->poc - fs_bottom->poc) > iabs(p_Vid->enc_picture->poc - fs_top->poc) )
          {
            ftype->mv[LIST_0][j][i][0]    = fs_top->motion.mv[LIST_0][jdiv][i][0];
            ftype->mv[LIST_0][j][i][1]    = fs_top->motion.mv[LIST_0][jdiv][i][1] ;
            ftype->mv[LIST_1][j][i][0]    = fs_top->motion.mv[LIST_1][jdiv][i][0];
            ftype->mv[LIST_1][j][i][1]    = fs_top->motion.mv[LIST_1][jdiv][i][1] ;
            ftype->ref_idx[LIST_0][j][i]  = fs_top->motion.ref_idx[LIST_0][jdiv][i];
            ftype->ref_idx[LIST_1][j][i]  = fs_top->motion.ref_idx[LIST_1][jdiv][i];
            ftype->ref_pic_id[LIST_0][j][i]   = fs->motion.ref_id[LIST_0][jj][i];
            ftype->ref_pic_id[LIST_1][j][i]   = fs->motion.ref_id[LIST_1][jj][i];

            p->is_long_term             = fs_top->is_long_term;
          }
          else
          {
            ftype->mv[LIST_0][j][i][0]      = fs_bottom->motion.mv[LIST_0][jdiv][i][0];
            ftype->mv[LIST_0][j][i][1]      = fs_bottom->motion.mv[LIST_0][jdiv][i][1] ;
            ftype->mv[LIST_1][j][i][0]      = fs_bottom->motion.mv[LIST_1][jdiv][i][0];
            ftype->mv[LIST_1][j][i][1]      = fs_bottom->motion.mv[LIST_1][jdiv][i][1] ;
            ftype->ref_idx[LIST_0][j][i]    = fs_bottom->motion.ref_idx[LIST_0][jdiv][i];
            ftype->ref_idx[LIST_1][j][i]    = fs_bottom->motion.ref_idx[LIST_1][jdiv][i];
            ftype->ref_pic_id[LIST_0][j][i] = fs->motion.ref_id[LIST_0][jj + 4][i];
            ftype->ref_pic_id[LIST_1][j][i] = fs->motion.ref_id[LIST_1][jj + 4][i];

            p->is_long_term             = fs_bottom->is_long_term;
          }
        }
        else
        {
          ftype->mv[LIST_0][j][i][0]      = fs->motion.mv[LIST_0][j][i][0];
          ftype->mv[LIST_0][j][i][1]      = fs->motion.mv[LIST_0][j][i][1] ;
          ftype->mv[LIST_1][j][i][0]      = fs->motion.mv[LIST_1][j][i][0];
          ftype->mv[LIST_1][j][i][1]      = fs->motion.mv[LIST_1][j][i][1] ;
          ftype->ref_idx[LIST_0][j][i]    = fs->motion.ref_idx[LIST_0][j][i];
          ftype->ref_idx[LIST_1][j][i]    = fs->motion.ref_idx[LIST_1][j][i];
          ftype->ref_pic_id[LIST_0][j][i] = fs->motion.ref_id[LIST_0][j][i];
          ftype->ref_pic_id[LIST_1][j][i] = fs->motion.ref_id[LIST_1][j][i];

          p->is_long_term             = fs->is_long_term;
        }
      }
    }
  }


  //! Generate field MVs from Frame MVs
  if (currSlice->structure || currSlice->mb_aff_frame_flag)
  {
    ftype = &p->frame;
    for (j=0 ; j<fs->size_y >> 3 ; j++)
    {
      jj = RSD(j);
      for (i=0 ; i<fs->size_x >> 2 ; i++)
      {
        ii = RSD(i);
        //! Do nothing if macroblock as field coded in MB-AFF
        if (!currSlice->mb_aff_frame_flag )
        {
          ftype->mv[LIST_0][j][i][0] = fs->motion.mv[LIST_0][jj][ii][0];
          ftype->mv[LIST_0][j][i][1] = fs->motion.mv[LIST_0][jj][ii][1];
          ftype->mv[LIST_1][j][i][0] = fs->motion.mv[LIST_1][jj][ii][0];
          ftype->mv[LIST_1][j][i][1] = fs->motion.mv[LIST_1][jj][ii][1];

          // Scaling of references is done here since it will not affect spatial direct (2*0 =0)

          if (fs->motion.ref_idx[LIST_0][jj][ii] == -1)
          {
            ftype->ref_idx   [LIST_0][j][i] = -1;
            ftype->ref_pic_id[LIST_0][j][i] = -1;
          }
          else
          {
            ftype->ref_idx   [LIST_0][j][i] = fs->motion.ref_idx[LIST_0][jj][ii] ;
            ftype->ref_pic_id[LIST_0][j][i] = fs->motion.ref_id [LIST_0][jj][ii];
          }

          if (fs->motion.ref_idx[LIST_1][jj][ii] == -1)
          {
            ftype->ref_idx   [LIST_1][j][i] = -1;
            ftype->ref_pic_id[LIST_1][j][i] = -1;
          }
          else
          {
            ftype->ref_idx   [LIST_1][j][i] = fs->motion.ref_idx[LIST_1][jj][ii];
            ftype->ref_pic_id[LIST_1][j][i] = fs->motion.ref_id [LIST_1][jj][ii];
          }

          p->is_long_term = fs->is_long_term;

          if (currSlice->direct_spatial_mv_pred_flag == 1)
          {
            ftype->moving_block[j][i] = (byte)
              !((!p->is_long_term
              && ((ftype->ref_idx[LIST_0][j][i] == 0)
              &&  (iabs(ftype->mv[LIST_0][j][i][0])>>1 == 0)
              &&  (iabs(ftype->mv[LIST_0][j][i][1])>>1 == 0)))
              || ((ftype->ref_idx[LIST_0][j][i] == -1)
              &&  (ftype->ref_idx[LIST_1][j][i] == 0)
              &&  (iabs(ftype->mv[LIST_1][j][i][0])>>1 == 0)
              &&  (iabs(ftype->mv[LIST_1][j][i][1])>>1 == 0)));
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
            p->bottom.moving_block[j][i] = (byte)
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
            p->top.moving_block[j][i] = (byte)
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


  if (!p_Vid->active_sps->frame_mbs_only_flag || p_Vid->active_sps->direct_8x8_inference_flag)
  {
    //! Use inference flag to remap mvs/references
    //! Frame with field co-located

    if (!currSlice->structure)
    {
      ftype = &p->frame;
      for (j=0 ; j < (fs->size_y>>2) ; j++)
      {
        jdiv = j>>1;
        jj = (j>>1) + 4*(j>>3);
        for (i=0 ; i < (fs->size_x>>2) ; i++)
        {

          if (fs->motion.field_frame[j][i])
          {
            if (iabs(p_Vid->enc_picture->poc - fs->bottom_field->poc) > iabs(p_Vid->enc_picture->poc - fs->top_field->poc))
            {
              ftype->mv[LIST_0][j][i][0] = fs->top_field->motion.mv[LIST_0][jdiv][i][0];
              ftype->mv[LIST_0][j][i][1] = fs->top_field->motion.mv[LIST_0][jdiv][i][1] ;
              ftype->mv[LIST_1][j][i][0] = fs->top_field->motion.mv[LIST_1][jdiv][i][0];
              ftype->mv[LIST_1][j][i][1] = fs->top_field->motion.mv[LIST_1][jdiv][i][1] ;

              ftype->ref_idx[LIST_0][j][i]  = fs->top_field->motion.ref_idx[LIST_0][jdiv][i];
              ftype->ref_idx[LIST_1][j][i]  = fs->top_field->motion.ref_idx[LIST_1][jdiv][i];
              ftype->ref_pic_id[LIST_0][j][i]   = fs->motion.ref_id[LIST_0][jj][i];
              ftype->ref_pic_id[LIST_1][j][i]   = fs->motion.ref_id[LIST_1][jj][i];
              p->is_long_term               = fs->top_field->is_long_term;
            }
            else
            {
              ftype->mv[LIST_0][j][i][0] = fs->bottom_field->motion.mv[LIST_0][jdiv][i][0];
              ftype->mv[LIST_0][j][i][1] = fs->bottom_field->motion.mv[LIST_0][jdiv][i][1] ;
              ftype->mv[LIST_1][j][i][0] = fs->bottom_field->motion.mv[LIST_1][jdiv][i][0];
              ftype->mv[LIST_1][j][i][1] = fs->bottom_field->motion.mv[LIST_1][jdiv][i][1] ;

              ftype->ref_idx[LIST_0][j][i]  = fs->bottom_field->motion.ref_idx[LIST_0][jdiv][i];
              ftype->ref_idx[LIST_1][j][i]  = fs->bottom_field->motion.ref_idx[LIST_1][jdiv][i];
              ftype->ref_pic_id[LIST_0][j][i] = fs->motion.ref_id[LIST_0][jj + 4][i];
              ftype->ref_pic_id[LIST_1][j][i] = fs->motion.ref_id[LIST_1][jj + 4][i];
              p->is_long_term             = fs->bottom_field->is_long_term;
            }
          }
        }
      }
    }
  }


  p->is_long_term = fs->is_long_term;

  if (!p_Vid->active_sps->frame_mbs_only_flag || p_Vid->active_sps->direct_8x8_inference_flag)
  {
    ftype = &p->frame;
    for (j=0 ; j < (fs->size_y>>2) ; j++)
    {
      jj = RSD(j);
      for (i=0 ; i < (fs->size_x>>2) ; i++)
      {
        ii = RSD(i);

        ftype->mv[LIST_0][j][i][0]=ftype->mv[LIST_0][jj][ii][0];
        ftype->mv[LIST_0][j][i][1]=ftype->mv[LIST_0][jj][ii][1];
        ftype->mv[LIST_1][j][i][0]=ftype->mv[LIST_1][jj][ii][0];
        ftype->mv[LIST_1][j][i][1]=ftype->mv[LIST_1][jj][ii][1];

        ftype->ref_idx[LIST_0][j][i]=ftype->ref_idx[LIST_0][jj][ii];
        ftype->ref_idx[LIST_1][j][i]=ftype->ref_idx[LIST_1][jj][ii];
        ftype->ref_pic_id[LIST_0][j][i] = ftype->ref_pic_id[LIST_0][jj][ii];
        ftype->ref_pic_id[LIST_1][j][i] = ftype->ref_pic_id[LIST_1][jj][ii];

        if (currSlice->direct_spatial_mv_pred_flag == 1)
        {
          ftype->moving_block[j][i]= (byte)
            !((!p->is_long_term
            && ((ftype->ref_idx[LIST_0][j][i] == 0)
            &&  (iabs(ftype->mv[LIST_0][j][i][0])>>1 == 0)
            &&  (iabs(ftype->mv[LIST_0][j][i][1])>>1 == 0)))
            || ((ftype->ref_idx[LIST_0][j][i] == -1)
            &&  (ftype->ref_idx[LIST_1][j][i] == 0)
            &&  (iabs(ftype->mv[LIST_1][j][i][0])>>1 == 0)
            &&  (iabs(ftype->mv[LIST_1][j][i][1])>>1 == 0)));
        }
      }
    }
  }
  else
  {
    ftype = &p->frame;
    for (j=0 ; j<fs->size_y >> 2 ; j++)
    {
      jj = RSD(j);
      for (i=0 ; i<fs->size_x >> 2 ; i++)
      {
        ii = RSD(i);
        //! Use inference flag to remap mvs/references
        ftype->mv[LIST_0][j][i][0]=fs->motion.mv[LIST_0][j][i][0];
        ftype->mv[LIST_0][j][i][1]=fs->motion.mv[LIST_0][j][i][1];
        ftype->mv[LIST_1][j][i][0]=fs->motion.mv[LIST_1][j][i][0];
        ftype->mv[LIST_1][j][i][1]=fs->motion.mv[LIST_1][j][i][1];

        ftype->ref_idx[LIST_0][j][i]=fs->motion.ref_idx[LIST_0][j][i];
        ftype->ref_idx[LIST_1][j][i]=fs->motion.ref_idx[LIST_1][j][i];
        ftype->ref_pic_id[LIST_0][j][i] = fs->motion.ref_id[LIST_0][j][i];
        ftype->ref_pic_id[LIST_1][j][i] = fs->motion.ref_id[LIST_1][j][i];

        if (currSlice->direct_spatial_mv_pred_flag == 1)
        {
          ftype->moving_block[j][i]= (byte)
            !((!p->is_long_term
            && ((ftype->ref_idx[LIST_0][j][i] == 0)
            &&  (iabs(ftype->mv[LIST_0][j][i][0])>>1 == 0)
            &&  (iabs(ftype->mv[LIST_0][j][i][1])>>1 == 0)))
            || ((ftype->ref_idx[LIST_0][j][i] == -1)
            &&  (ftype->ref_idx[LIST_1][j][i] == 0)
            &&  (iabs(ftype->mv[LIST_1][j][i][0])>>1 == 0)
            &&  (iabs(ftype->mv[LIST_1][j][i][1])>>1 == 0)));
        }
      }
    }
  }


  if (currSlice->direct_spatial_mv_pred_flag ==0)
  {
    ftype = &p->frame;
    for (j=0 ; j<fs->size_y >> 2 ; j++)
    {
      for (i=0 ; i<fs->size_x >> 2 ; i++)
      {
        if ((!currSlice->mb_aff_frame_flag &&!currSlice->structure && fs->motion.field_frame[j][i]) || (currSlice->mb_aff_frame_flag && fs->motion.field_frame[j][i]))
        {
          ftype->mv[LIST_0][j][i][1] *= 2;
          ftype->mv[LIST_1][j][i][1] *= 2;
        }
        else  if (currSlice->structure && !fs->motion.field_frame[j][i])
        {
          ftype->mv[LIST_0][j][i][1] /= 2;
          ftype->mv[LIST_1][j][i][1] /= 2;
        }

      }
    }

    for (j=0; j<2 + (currSlice->mb_aff_frame_flag * 4);j+=2)
    {
      for (i=0; i<p_Vid->listXsize[j];i++)
      {
        int prescale, iTRb, iTRp;

        if (j==0)
        {
          iTRb = iClip3( -128, 127, p_Vid->enc_picture->poc - listX[LIST_0 + j][i]->poc );
        }
        else if (j == 2)
        {
          iTRb = iClip3( -128, 127, p_Vid->enc_picture->top_poc - listX[LIST_0 + j][i]->poc );
        }
        else
        {
          iTRb = iClip3( -128, 127, p_Vid->enc_picture->bottom_poc - listX[LIST_0 + j][i]->poc );
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
  VideoParameters *p_Vid = currSlice->p_Vid;
  StorablePicture *fs, *fs_top, *fs_bottom;
  MotionParams *ftype;
  int i,j, ii, jj, jdiv;
  int np = currSlice->p_Vid->colour_plane_id;

  fs_top = fs_bottom = fs = listX[LIST_1 ][0];

  if (currSlice->mb_aff_frame_flag)
  {
    fs_top= listX[LIST_1 + 2][0];
    fs_bottom= listX[LIST_1 + 4][0];
  }
  else
  {
    if (currSlice->structure!=FRAME)
    {
      if ((currSlice->structure != fs->structure) && (fs->coded_frame))
      {
        if (currSlice->structure==TOP_FIELD)
        {
          fs_top = fs_bottom = fs = listX[LIST_1 ][0]->top_field;
        }
        else
        {
          fs_top = fs_bottom = fs = listX[LIST_1 ][0]->bottom_field;
        }
      }
    }
  }

  if (!p_Vid->active_sps->frame_mbs_only_flag || p_Vid->active_sps->direct_8x8_inference_flag)
  {
    ftype = &p->frame;
    for (j=0 ; j<fs->size_y >> 2 ; j++)
    {
      jdiv = j >> 1;
      jj = (j >> 1) + 4 * (j >> 3);
      for (i=0 ; i<fs->size_x >> 2 ; i++)
      {

        if (currSlice->mb_aff_frame_flag && fs->motion.field_frame[j][i])
        {
          //! Assign frame buffers for field MBs
          //! Check whether we should use top or bottom field mvs.
          //! Depending on the assigned poc values.

          if (iabs(p_Vid->enc_picture->poc - fs_bottom->poc) > iabs(p_Vid->enc_picture->poc - fs_top->poc) )
          {
            ftype->mv[LIST_0][j][i][0]    = fs_top->JVmotion[np].mv[LIST_0][jdiv][i][0];
            ftype->mv[LIST_0][j][i][1]    = fs_top->JVmotion[np].mv[LIST_0][jdiv][i][1] ;
            ftype->mv[LIST_1][j][i][0]    = fs_top->JVmotion[np].mv[LIST_1][jdiv][i][0];
            ftype->mv[LIST_1][j][i][1]    = fs_top->JVmotion[np].mv[LIST_1][jdiv][i][1] ;
            ftype->ref_idx[LIST_0][j][i]  = fs_top->JVmotion[np].ref_idx[LIST_0][jdiv][i];
            ftype->ref_idx[LIST_1][j][i]  = fs_top->JVmotion[np].ref_idx[LIST_1][jdiv][i];
            ftype->ref_pic_id[LIST_0][j][i]   = fs->JVmotion[np].ref_id[LIST_0][jj][i];
            ftype->ref_pic_id[LIST_1][j][i]   = fs->JVmotion[np].ref_id[LIST_1][jj][i];

            p->is_long_term             = fs_top->is_long_term;
          }
          else
          {
            ftype->mv[LIST_0][j][i][0]      = fs_bottom->JVmotion[np].mv[LIST_0][jdiv][i][0];
            ftype->mv[LIST_0][j][i][1]      = fs_bottom->JVmotion[np].mv[LIST_0][jdiv][i][1] ;
            ftype->mv[LIST_1][j][i][0]      = fs_bottom->JVmotion[np].mv[LIST_1][jdiv][i][0];
            ftype->mv[LIST_1][j][i][1]      = fs_bottom->JVmotion[np].mv[LIST_1][jdiv][i][1] ;
            ftype->ref_idx[LIST_0][j][i]    = fs_bottom->JVmotion[np].ref_idx[LIST_0][jdiv][i];
            ftype->ref_idx[LIST_1][j][i]    = fs_bottom->JVmotion[np].ref_idx[LIST_1][jdiv][i];
            ftype->ref_pic_id[LIST_0][j][i] = fs->JVmotion[np].ref_id[LIST_0][jj + 4][i];
            ftype->ref_pic_id[LIST_1][j][i] = fs->JVmotion[np].ref_id[LIST_1][jj + 4][i];

            p->is_long_term             = fs_bottom->is_long_term;
          }
        }
        else
        {
          ftype->mv[LIST_0][j][i][0]      = fs->JVmotion[np].mv[LIST_0][j][i][0];
          ftype->mv[LIST_0][j][i][1]      = fs->JVmotion[np].mv[LIST_0][j][i][1] ;
          ftype->mv[LIST_1][j][i][0]      = fs->JVmotion[np].mv[LIST_1][j][i][0];
          ftype->mv[LIST_1][j][i][1]      = fs->JVmotion[np].mv[LIST_1][j][i][1] ;
          ftype->ref_idx[LIST_0][j][i]    = fs->JVmotion[np].ref_idx[LIST_0][j][i];
          ftype->ref_idx[LIST_1][j][i]    = fs->JVmotion[np].ref_idx[LIST_1][j][i];
          ftype->ref_pic_id[LIST_0][j][i] = fs->JVmotion[np].ref_id[LIST_0][j][i];
          ftype->ref_pic_id[LIST_1][j][i] = fs->JVmotion[np].ref_id[LIST_1][j][i];

          p->is_long_term             = fs->is_long_term;
        }
      }
    }
  }


  //! Generate field MVs from Frame MVs
  if (currSlice->structure || currSlice->mb_aff_frame_flag)
  {
    ftype = &p->frame;
    for (j=0 ; j<fs->size_y >> 3 ; j++)
    {
      jj = RSD(j);
      for (i=0 ; i<fs->size_x >> 2 ; i++)
      {
        ii = RSD(i);
        //! Do nothing if macroblock as field coded in MB-AFF
        if (!currSlice->mb_aff_frame_flag )
        {
          ftype->mv[LIST_0][j][i][0] = fs->JVmotion[np].mv[LIST_0][jj][ii][0];
          ftype->mv[LIST_0][j][i][1] = fs->JVmotion[np].mv[LIST_0][jj][ii][1];
          ftype->mv[LIST_1][j][i][0] = fs->JVmotion[np].mv[LIST_1][jj][ii][0];
          ftype->mv[LIST_1][j][i][1] = fs->JVmotion[np].mv[LIST_1][jj][ii][1];

          // Scaling of references is done here since it will not affect spatial direct (2*0 =0)

          if (fs->motion.ref_idx[LIST_0][jj][ii] == -1)
          {
            ftype->ref_idx   [LIST_0][j][i] = -1;
            ftype->ref_pic_id[LIST_0][j][i] = -1;
          }
          else
          {
            ftype->ref_idx   [LIST_0][j][i] = fs->JVmotion[np].ref_idx[LIST_0][jj][ii] ;
            ftype->ref_pic_id[LIST_0][j][i] = fs->JVmotion[np].ref_id[LIST_0][jj][ii];
          }

          if (fs->motion.ref_idx[LIST_1][jj][ii] == -1)
          {
            ftype->ref_idx   [LIST_1][j][i] = -1;
            ftype->ref_pic_id[LIST_1][j][i] = -1;
          }
          else
          {
            ftype->ref_idx   [LIST_1][j][i] = fs->JVmotion[np].ref_idx[LIST_1][jj][ii];
            ftype->ref_pic_id[LIST_1][j][i] = fs->JVmotion[np].ref_id[LIST_1][jj][ii];
          }

          p->is_long_term = fs->is_long_term;

          if (currSlice->direct_spatial_mv_pred_flag == 1)
          {
            ftype->moving_block[j][i] = (byte)
              !((!p->is_long_term
              && ((ftype->ref_idx[LIST_0][j][i] == 0)
              &&  (iabs(ftype->mv[LIST_0][j][i][0])>>1 == 0)
              &&  (iabs(ftype->mv[LIST_0][j][i][1])>>1 == 0)))
              || ((ftype->ref_idx[LIST_0][j][i] == -1)
              &&  (ftype->ref_idx[LIST_1][j][i] == 0)
              &&  (iabs(ftype->mv[LIST_1][j][i][0])>>1 == 0)
              &&  (iabs(ftype->mv[LIST_1][j][i][1])>>1 == 0)));
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
            p->bottom.moving_block[j][i] = (byte)
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
            p->top.moving_block[j][i] = (byte)
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


  if (!p_Vid->active_sps->frame_mbs_only_flag || p_Vid->active_sps->direct_8x8_inference_flag)
  {
    //! Use inference flag to remap mvs/references
    //! Frame with field co-located

    if (!currSlice->structure)
    {
      ftype = &p->frame;
      for (j=0 ; j < (fs->size_y>>2) ; j++)
      {
        jdiv = j>>1;
        jj = (j>>1) + 4*(j>>3);
        for (i=0 ; i < (fs->size_x>>2) ; i++)
        {

          if (fs->motion.field_frame[j][i])
          {
            if (iabs(p_Vid->enc_picture->poc - fs->bottom_field->poc) > iabs(p_Vid->enc_picture->poc - fs->top_field->poc))
            {
              ftype->mv[LIST_0][j][i][0] = fs->top_field->JVmotion[np].mv[LIST_0][jdiv][i][0];
              ftype->mv[LIST_0][j][i][1] = fs->top_field->JVmotion[np].mv[LIST_0][jdiv][i][1] ;
              ftype->mv[LIST_1][j][i][0] = fs->top_field->JVmotion[np].mv[LIST_1][jdiv][i][0];
              ftype->mv[LIST_1][j][i][1] = fs->top_field->JVmotion[np].mv[LIST_1][jdiv][i][1] ;

              ftype->ref_idx[LIST_0][j][i]  = fs->top_field->JVmotion[np].ref_idx[LIST_0][jdiv][i];
              ftype->ref_idx[LIST_1][j][i]  = fs->top_field->JVmotion[np].ref_idx[LIST_1][jdiv][i];
              ftype->ref_pic_id[LIST_0][j][i]   = fs->JVmotion[np].ref_id[LIST_0][jj][i];
              ftype->ref_pic_id[LIST_1][j][i]   = fs->JVmotion[np].ref_id[LIST_1][jj][i];
              p->is_long_term               = fs->top_field->is_long_term;
            }
            else
            {
              ftype->mv[LIST_0][j][i][0] = fs->bottom_field->JVmotion[np].mv[LIST_0][jdiv][i][0];
              ftype->mv[LIST_0][j][i][1] = fs->bottom_field->JVmotion[np].mv[LIST_0][jdiv][i][1] ;
              ftype->mv[LIST_1][j][i][0] = fs->bottom_field->JVmotion[np].mv[LIST_1][jdiv][i][0];
              ftype->mv[LIST_1][j][i][1] = fs->bottom_field->JVmotion[np].mv[LIST_1][jdiv][i][1] ;

              ftype->ref_idx[LIST_0][j][i]  = fs->bottom_field->JVmotion[np].ref_idx[LIST_0][jdiv][i];
              ftype->ref_idx[LIST_1][j][i]  = fs->bottom_field->JVmotion[np].ref_idx[LIST_1][jdiv][i];
              ftype->ref_pic_id[LIST_0][j][i] = fs->JVmotion[np].ref_id[LIST_0][jj + 4][i];
              ftype->ref_pic_id[LIST_1][j][i] = fs->JVmotion[np].ref_id[LIST_1][jj + 4][i];
              p->is_long_term             = fs->bottom_field->is_long_term;
            }
          }
        }
      }
    }
  }


  p->is_long_term = fs->is_long_term;

  if (!p_Vid->active_sps->frame_mbs_only_flag || p_Vid->active_sps->direct_8x8_inference_flag)
  {
    ftype = &p->frame;
    for (j=0 ; j < (fs->size_y>>2) ; j++)
    {
      jj = RSD(j);
      for (i=0 ; i < (fs->size_x>>2) ; i++)
      {
        ii = RSD(i);

        ftype->mv[LIST_0][j][i][0]=ftype->mv[LIST_0][jj][ii][0];
        ftype->mv[LIST_0][j][i][1]=ftype->mv[LIST_0][jj][ii][1];
        ftype->mv[LIST_1][j][i][0]=ftype->mv[LIST_1][jj][ii][0];
        ftype->mv[LIST_1][j][i][1]=ftype->mv[LIST_1][jj][ii][1];

        ftype->ref_idx[LIST_0][j][i]=ftype->ref_idx[LIST_0][jj][ii];
        ftype->ref_idx[LIST_1][j][i]=ftype->ref_idx[LIST_1][jj][ii];
        ftype->ref_pic_id[LIST_0][j][i] = ftype->ref_pic_id[LIST_0][jj][ii];
        ftype->ref_pic_id[LIST_1][j][i] = ftype->ref_pic_id[LIST_1][jj][ii];

        if (currSlice->direct_spatial_mv_pred_flag == 1)
        {
          ftype->moving_block[j][i]= (byte)
            !((!p->is_long_term
            && ((ftype->ref_idx[LIST_0][j][i] == 0)
            &&  (iabs(ftype->mv[LIST_0][j][i][0])>>1 == 0)
            &&  (iabs(ftype->mv[LIST_0][j][i][1])>>1 == 0)))
            || ((ftype->ref_idx[LIST_0][j][i] == -1)
            &&  (ftype->ref_idx[LIST_1][j][i] == 0)
            &&  (iabs(ftype->mv[LIST_1][j][i][0])>>1 == 0)
            &&  (iabs(ftype->mv[LIST_1][j][i][1])>>1 == 0)));
        }
      }
    }
  }
  else
  {
    ftype = &p->frame;
    for (j=0 ; j<fs->size_y >> 2 ; j++)
    {
      jj = RSD(j);
      for (i=0 ; i<fs->size_x >> 2 ; i++)
      {
        ii = RSD(i);
        //! Use inference flag to remap mvs/references
        ftype->mv[LIST_0][j][i][0]=fs->JVmotion[np].mv[LIST_0][j][i][0];
        ftype->mv[LIST_0][j][i][1]=fs->JVmotion[np].mv[LIST_0][j][i][1];
        ftype->mv[LIST_1][j][i][0]=fs->JVmotion[np].mv[LIST_1][j][i][0];
        ftype->mv[LIST_1][j][i][1]=fs->JVmotion[np].mv[LIST_1][j][i][1];

        ftype->ref_idx[LIST_0][j][i]=fs->JVmotion[np].ref_idx[LIST_0][j][i];
        ftype->ref_idx[LIST_1][j][i]=fs->JVmotion[np].ref_idx[LIST_1][j][i];
        ftype->ref_pic_id[LIST_0][j][i] = fs->JVmotion[np].ref_id[LIST_0][j][i];
        ftype->ref_pic_id[LIST_1][j][i] = fs->JVmotion[np].ref_id[LIST_1][j][i];

        if (currSlice->direct_spatial_mv_pred_flag == 1)
        {
          ftype->moving_block[j][i]= (byte)
            !((!p->is_long_term
            && ((ftype->ref_idx[LIST_0][j][i] == 0)
            &&  (iabs(ftype->mv[LIST_0][j][i][0])>>1 == 0)
            &&  (iabs(ftype->mv[LIST_0][j][i][1])>>1 == 0)))
            || ((ftype->ref_idx[LIST_0][j][i] == -1)
            &&  (ftype->ref_idx[LIST_1][j][i] == 0)
            &&  (iabs(ftype->mv[LIST_1][j][i][0])>>1 == 0)
            &&  (iabs(ftype->mv[LIST_1][j][i][1])>>1 == 0)));
        }
      }
    }
  }


  if (currSlice->direct_spatial_mv_pred_flag ==0)
  {
    ftype = &p->frame;
    for (j=0 ; j<fs->size_y >> 2 ; j++)
    {
      for (i=0 ; i<fs->size_x >> 2 ; i++)
      {
        if ((!currSlice->mb_aff_frame_flag &&!currSlice->structure && fs->motion.field_frame[j][i]) || (currSlice->mb_aff_frame_flag && fs->motion.field_frame[j][i]))
        {
          ftype->mv[LIST_0][j][i][1] *= 2;
          ftype->mv[LIST_1][j][i][1] *= 2;
        }
        else  if (currSlice->structure && !fs->motion.field_frame[j][i])
        {
          ftype->mv[LIST_0][j][i][1] /= 2;
          ftype->mv[LIST_1][j][i][1] /= 2;
        }

      }
    }

    for (j=0; j<2 + (currSlice->mb_aff_frame_flag * 4);j+=2)
    {
      for (i=0; i<p_Vid->listXsize[j];i++)
      {
        int prescale, iTRb, iTRp;

        if (j==0)
        {
          iTRb = iClip3( -128, 127, p_Vid->enc_picture->poc - listX[LIST_0 + j][i]->poc );
        }
        else if (j == 2)
        {
          iTRb = iClip3( -128, 127, p_Vid->enc_picture->top_poc - listX[LIST_0 + j][i]->poc );
        }
        else
        {
          iTRb = iClip3( -128, 127, p_Vid->enc_picture->bottom_poc - listX[LIST_0 + j][i]->poc );
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
 *    Copy StorablePicture parameters
 *    for 4:4:4 Independent mode
 *
 ************************************************************************
 */
void copy_storable_param_JV( VideoParameters *p_Vid, int nplane, StorablePicture *d, StorablePicture *s )
{
  int md_size = (p_Vid->height / BLOCK_SIZE) * (p_Vid->width / BLOCK_SIZE);
  int ref_size = p_Vid->active_sps->frame_mbs_only_flag ? 2 * md_size : 6 * md_size;

  // copy ref_idx
  memcpy( d->JVmotion[nplane].ref_idx[0][0], s->motion.ref_idx[0][0], 2 * md_size * sizeof(byte) );

  // copy ref_pic_id
  memcpy( d->JVmotion[nplane].ref_pic_id[0][0], s->motion.ref_pic_id[0][0], ref_size * sizeof(int64) );

  // copy motion.ref_id
  memcpy( d->JVmotion[nplane].ref_id[0][0], s->motion.ref_id[0][0], ref_size * sizeof(int64));

  // copy mv
  memcpy( d->JVmotion[nplane].mv[0][0][0], s->motion.mv[0][0][0], 2 * md_size * 2 * sizeof(short) );
}
