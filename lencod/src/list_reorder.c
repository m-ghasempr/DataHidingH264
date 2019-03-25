/*!
 **************************************************************************************
 * \file
 *    slice.c
 * \brief
 *    generate the slice header, setup the bit buffer for slices,
 *    and generates the slice NALU(s)

 * \author
 *    Main contributors (see contributors.h for copyright, address and affiliation details)
 *     - Athanasios Leontaris            <aleon@dolby.com>
 *     - Karsten Sühring                 <suehring@hhi.de> 
 *     - Alexis Michael Tourapis         <alexismt@ieee.org> 
 ***************************************************************************************
 */

#include "contributors.h"

#include <math.h>
#include <float.h>

#include "global.h"
#include "image.h"
#include "wp.h"
#include "list_reorder.h"

/*!
 ************************************************************************
 * \brief
 *    init_ref_pic_list_reordering initializations should go here
 ************************************************************************
 */
void init_ref_pic_list_reordering(Slice* currSlice)
{
  currSlice->ref_pic_list_reordering_flag[LIST_0] = 0;
  currSlice->ref_pic_list_reordering_flag[LIST_1] = 0;

  currSlice->poc_ref_pic_reorder_frame = poc_ref_pic_reorder_frame_default;

}

/*!
************************************************************************
* \brief
*    decide reference picture reordering, Frame only
************************************************************************
*/
void poc_ref_pic_reorder_frame_default(Slice *currSlice, StorablePicture **list, 
                               unsigned num_ref_idx_lX_active, int *reordering_of_pic_nums_idc, int *abs_diff_pic_num_minus1, int *long_term_pic_idx, int list_no)
{
  VideoParameters *p_Vid = currSlice->p_Vid;
  InputParameters *p_Inp = currSlice->p_Inp;
  StorablePicture *p_Enc_Pic = p_Vid->enc_picture;

  unsigned int i,j,k;

  int currPicNum, picNumLXPred;

  int default_order[32];
  int re_order[32];
  int tmp_reorder[32];
  int list_sign[32];
  int reorder_stop, no_reorder;
  int poc_diff[32];
  int tmp_value, diff;

  int abs_poc_dist;
  int maxPicNum;
  unsigned int num_refs;

  maxPicNum  = p_Vid->max_frame_num;
  currPicNum = p_Vid->frame_num;

  picNumLXPred = currPicNum;

  // First assign default list order.
  for (i=0; i<num_ref_idx_lX_active; i++)
  {
    default_order[i] = list[i]->pic_num;
  }

  // Now access all references in buffer and assign them
  // to a potential reordering list. For each one of these
  // references compute the poc distance compared to current
  // frame.
  num_refs = p_Vid->p_Dpb->ref_frames_in_buffer;
  for (i=0; i<p_Vid->p_Dpb->ref_frames_in_buffer; i++)
  {
    poc_diff[i] = 0xFFFF;
    re_order[i] = p_Vid->p_Dpb->fs_ref[i]->frame->pic_num;

    if (p_Vid->p_Dpb->fs_ref[i]->is_used==3 && (p_Vid->p_Dpb->fs_ref[i]->frame->used_for_reference)&&(!p_Vid->p_Dpb->fs_ref[i]->frame->is_long_term))
    {
      abs_poc_dist = iabs(p_Vid->p_Dpb->fs_ref[i]->frame->poc - p_Enc_Pic->poc) ;
      poc_diff[i] = abs_poc_dist;
      if (list_no == LIST_0)
      {
        list_sign[i] = (p_Enc_Pic->poc < p_Vid->p_Dpb->fs_ref[i]->frame->poc) ? +1 : -1;
      }
      else
      {
        list_sign[i] = (p_Enc_Pic->poc > p_Vid->p_Dpb->fs_ref[i]->frame->poc) ? +1 : -1;
      }
    }
  }

  // now sort these references based on poc (temporal) distance
  for (i = 0; i < (num_refs - 1); i++)
  {
    for (j = i + 1; j < num_refs; j++)
    {
      if (poc_diff[i]>poc_diff[j] || (poc_diff[i] == poc_diff[j] && list_sign[j] > list_sign[i]))
      {

        tmp_value = poc_diff[i];
        poc_diff[i] = poc_diff[j];
        poc_diff[j] = tmp_value;
        tmp_value  = re_order[i];
        re_order[i] = re_order[j];
        re_order[j] = tmp_value ;
        tmp_value  = list_sign[i];
        list_sign[i] = list_sign[j];
        list_sign[j] = tmp_value;
      }
    }
  }

  // populate list with selections from the pre-analysis stage
  if ( p_Inp->WPMCPrecision 
    && p_Vid->pWPX->curr_wp_rd_pass->algorithm != WP_REGULAR
    && p_Vid->pWPX->num_wp_ref_list[list_no] )
  {
    for (i=0; i<num_ref_idx_lX_active; i++)
    {
      re_order[i] = p_Vid->pWPX->wp_ref_list[list_no][i].PicNum;
    }
  }

  // Check versus default list to see if any
  // change has happened
  no_reorder = 1;
  for (i=0; i<num_ref_idx_lX_active; i++)
  {
    if (default_order[i] != re_order[i])
    {
      no_reorder = 0;
    }
  }

  // If different, then signal reordering
  if (no_reorder==0)
  {
    for (i=0; i<num_ref_idx_lX_active; i++)
    {
      diff = re_order[i]-picNumLXPred;
      if (diff <= 0)
      {
        reordering_of_pic_nums_idc[i] = 0;
        abs_diff_pic_num_minus1[i] = iabs(diff)-1;
        if (abs_diff_pic_num_minus1[i] < 0)
          abs_diff_pic_num_minus1[i] = maxPicNum -1;
      }
      else
      {
        reordering_of_pic_nums_idc[i] = 1;
        abs_diff_pic_num_minus1[i] = iabs(diff)-1;
      }
      picNumLXPred = re_order[i];

      tmp_reorder[i] = re_order[i];

      k = i;
      for (j=i; j<num_ref_idx_lX_active; j++)
      {
        if (default_order[j] != re_order[i])
        {
          ++k;
          tmp_reorder[k] = default_order[j];
        }
      }
      reorder_stop = 1;
      for(j=i+1; j<num_ref_idx_lX_active; j++)
      {
        if (tmp_reorder[j] != re_order[j])
        {
          reorder_stop = 0;
          break;
        }
      }

      if (reorder_stop==1)
      {
        ++i;
        break;
      }
      memcpy ( default_order, tmp_reorder, num_ref_idx_lX_active * sizeof(int));

    }
    reordering_of_pic_nums_idc[i] = 3;

    memcpy ( default_order, tmp_reorder, num_ref_idx_lX_active * sizeof(int));

    if (list_no==0)
    {
      currSlice->ref_pic_list_reordering_flag[LIST_0] = 1;
    }
    else
    {
      currSlice->ref_pic_list_reordering_flag[LIST_1] = 1;
    }
  }
}

/*!
************************************************************************
* \brief
*    decide reference picture reordering, Field only
************************************************************************
*/
void poc_ref_pic_reorder_field(Slice *currSlice, StorablePicture **list, unsigned num_ref_idx_lX_active, int *reordering_of_pic_nums_idc, int *abs_diff_pic_num_minus1, int *long_term_pic_idx, int list_no)
{
  VideoParameters *p_Vid = currSlice->p_Vid;
  //InputParameters *p_Inp = currSlice->p_Inp;
  StorablePicture *p_Enc_Pic = p_Vid->enc_picture;

  unsigned int i,j,k;

  int currPicNum, picNumLXPred;

  int default_order[32];
  int re_order[32];
  int tmp_reorder[32];
  int list_sign[32];  
  int poc_diff[32];
  int fld_type[32];

  int reorder_stop, no_reorder;
  int tmp_value, diff;

  int abs_poc_dist;
  int maxPicNum;
  unsigned int num_refs;

  int field_used[2] = {1, 2};
  int fld, idx, num_flds;

  unsigned int top_idx = 0;
  unsigned int bot_idx = 0;
  unsigned int list_size = 0;

  StorablePicture *pField[2]; // 0: TOP_FIELD, 1: BOTTOM_FIELD
  FrameStore      *pFrameStore; 

  maxPicNum  = 2 * p_Vid->max_frame_num;
  currPicNum = 2 * p_Vid->frame_num + 1;

  picNumLXPred = currPicNum;

  // First assign default list order.
  for (i=0; i<num_ref_idx_lX_active; i++)
  {
    default_order[i] = list[i]->pic_num;
  }

  // Now access all references in buffer and assign them
  // to a potential reordering list. For each one of these
  // references compute the poc distance compared to current
  // frame.  
  // look for eligible fields
  idx = 0;

  for (i=0; i<p_Vid->p_Dpb->ref_frames_in_buffer; i++)
  {
    pFrameStore = p_Vid->p_Dpb->fs_ref[i];
    pField[0]   = pFrameStore->top_field;
    pField[1]   = pFrameStore->bottom_field;
    num_flds    = (currSlice->structure == BOTTOM_FIELD && (p_Enc_Pic->poc == (pField[0]->poc + 1) ) ) ? 1 : 2;

    poc_diff[2*i    ] = 0xFFFF;
    poc_diff[2*i + 1] = 0xFFFF;

    for ( fld = 0; fld < num_flds; fld++ )
    {
      if ( (pFrameStore->is_used & field_used[fld]) && pField[fld]->used_for_reference && !(pField[fld]->is_long_term) )
      {
        abs_poc_dist = iabs(pField[fld]->poc - p_Enc_Pic->poc) ;
        poc_diff[idx] = abs_poc_dist;
        re_order[idx] = pField[fld]->pic_num;
        fld_type[idx] = fld + 1;

        if (list_no == LIST_0)
        {
          list_sign[idx] = (p_Enc_Pic->poc < pField[fld]->poc) ? +1 : -1;
        }
        else
        {
          list_sign[idx] = (p_Enc_Pic->poc > pField[fld]->poc) ? +1 : -1;
        }
        idx++;
      }
    }
  }
  num_refs = idx;

  // now sort these references based on poc (temporal) distance
  for (i=0; i < num_refs-1; i++)
  {
    for (j = (i + 1); j < num_refs; j++)
    {
      if (poc_diff[i] > poc_diff[j] || (poc_diff[i] == poc_diff[j] && list_sign[j] > list_sign[i]))
      {
        // poc_diff
        tmp_value   = poc_diff[i];
        poc_diff[i] = poc_diff[j];
        poc_diff[j] = tmp_value;
        // re_order (PicNum)
        tmp_value   = re_order[i];
        re_order[i] = re_order[j];
        re_order[j] = tmp_value;
        // list_sign
        tmp_value    = list_sign[i];
        list_sign[i] = list_sign[j];
        list_sign[j] = tmp_value;
        // fld_type
        tmp_value   = fld_type[i];
        fld_type[i] = fld_type[j];
        fld_type[j] = tmp_value ;
      }
    }
  }

  if (currSlice->structure == TOP_FIELD)
  {
    while ((top_idx < num_refs)||(bot_idx < num_refs))
    {
      for ( ; top_idx < num_refs; top_idx++)
      {
        if ( fld_type[top_idx] == TOP_FIELD )
        {
          tmp_reorder[list_size] = re_order[top_idx];
          list_size++;
          top_idx++;
          break;
        }
      }
      for ( ; bot_idx < num_refs; bot_idx++)
      {
        if ( fld_type[bot_idx] == BOTTOM_FIELD )
        {
          tmp_reorder[list_size] = re_order[bot_idx];
          list_size++;
          bot_idx++;
          break;
        }
      }
    }
  }
  if (currSlice->structure == BOTTOM_FIELD)
  {
    while ((top_idx < num_refs)||(bot_idx < num_refs))
    {
      for ( ; bot_idx < num_refs; bot_idx++)
      {
        if ( fld_type[bot_idx] == BOTTOM_FIELD )
        {
          tmp_reorder[list_size] = re_order[bot_idx];
          list_size++;
          bot_idx++;
          break;
        }
      }
      for ( ; top_idx < num_refs; top_idx++)
      {
        if ( fld_type[top_idx] == TOP_FIELD )
        {
          tmp_reorder[list_size] = re_order[top_idx];
          list_size++;
          top_idx++;
          break;
        }
      }
    }
  }

  // copy to final matrix
  list_size = imin( list_size, 32 );
  for ( i = 0; i < list_size; i++ )
  {
    re_order[i] = tmp_reorder[i];
  }

  // Check versus default list to see if any
  // change has happened
  no_reorder = 1;
  for (i=0; i<num_ref_idx_lX_active; i++)
  {
    if (default_order[i] != re_order[i])
    {
      no_reorder = 0;
    }
  }

  // If different, then signal reordering
  if (no_reorder == 0)
  {
    for (i=0; i<num_ref_idx_lX_active; i++)
    {
      diff = re_order[i] - picNumLXPred;
      if (diff <= 0)
      {
        reordering_of_pic_nums_idc[i] = 0;
        abs_diff_pic_num_minus1[i] = iabs(diff)-1;
        if (abs_diff_pic_num_minus1[i] < 0)
          abs_diff_pic_num_minus1[i] = maxPicNum -1;
      }
      else
      {
        reordering_of_pic_nums_idc[i] = 1;
        abs_diff_pic_num_minus1[i] = iabs(diff)-1;
      }
      picNumLXPred = re_order[i];

      tmp_reorder[i] = re_order[i];

      k = i;
      for (j = i; j < num_ref_idx_lX_active; j++)
      {
        if (default_order[j] != re_order[i])
        {
          ++k;
          tmp_reorder[k] = default_order[j];
        }
      }
      reorder_stop = 1;
      for(j=i+1; j<num_ref_idx_lX_active; j++)
      {
        if (tmp_reorder[j] != re_order[j])
        {
          reorder_stop = 0;
          break;
        }
      }

      if (reorder_stop==1)
      {
        ++i;
        break;
      }

      memcpy ( default_order, tmp_reorder, num_ref_idx_lX_active * sizeof(int));
    }
    reordering_of_pic_nums_idc[i] = 3;

    memcpy ( default_order, tmp_reorder, num_ref_idx_lX_active * sizeof(int));

    if (list_no==0)
    {
      currSlice->ref_pic_list_reordering_flag[LIST_0] = 1;
    }
    else
    {
      currSlice->ref_pic_list_reordering_flag[LIST_1] = 1;
    }
  }
}

/*!
************************************************************************
* \brief
*    Reorder lists
************************************************************************
*/

void reorder_lists( Slice *currSlice )
{
  VideoParameters *p_Vid = currSlice->p_Vid;
  InputParameters *p_Inp = currSlice->p_Inp;

  // Perform reordering based on poc distances for HierarchicalCoding
  if ( currSlice->slice_type == P_SLICE && p_Inp->ReferenceReorder)
  {
    int i, num_ref;

    alloc_ref_pic_list_reordering_buffer(currSlice);

    for (i = 0; i < currSlice->num_ref_idx_active[LIST_0] + 1; i++)
    {
      currSlice->reordering_of_pic_nums_idc[LIST_0][i] = 3;
      currSlice->abs_diff_pic_num_minus1[LIST_0][i] = 0;
      currSlice->long_term_pic_idx[LIST_0][i] = 0;
    }

    num_ref = currSlice->num_ref_idx_active[LIST_0];
    if ( currSlice->structure == FRAME )
      currSlice->poc_ref_pic_reorder_frame(currSlice, p_Vid->listX[LIST_0], num_ref,
      currSlice->reordering_of_pic_nums_idc[LIST_0],
      currSlice->abs_diff_pic_num_minus1[LIST_0],
      currSlice->long_term_pic_idx[LIST_0], LIST_0);
    else
    {
      poc_ref_pic_reorder_field(currSlice, p_Vid->listX[LIST_0], num_ref,
        currSlice->reordering_of_pic_nums_idc[LIST_0],
        currSlice->abs_diff_pic_num_minus1[LIST_0],
        currSlice->long_term_pic_idx[LIST_0], LIST_0);
    }
    //reference picture reordering
    reorder_ref_pic_list ( currSlice, p_Vid->listX, p_Vid->listXsize, LIST_0);
  }
}

/*!
************************************************************************
* \brief
*    Reorder lists
************************************************************************
*/

void wp_mcprec_reorder_lists( Slice *currSlice )
{

  VideoParameters *p_Vid = currSlice->p_Vid;
  int i, num_ref;

  wpxModifyRefPicList( currSlice );

  alloc_ref_pic_list_reordering_buffer(currSlice);

  for (i = 0; i < currSlice->num_ref_idx_active[LIST_0] + 1; i++)
  {
    currSlice->reordering_of_pic_nums_idc[LIST_0][i] = 3;
    currSlice->abs_diff_pic_num_minus1[LIST_0][i] = 0;
    currSlice->long_term_pic_idx[LIST_0][i] = 0;
  }

  if (currSlice->slice_type == B_SLICE) // type should be part of currSlice not p_Vid
  {
    for (i = 0; i < currSlice->num_ref_idx_active[LIST_1] + 1; i++)
    {
      currSlice->reordering_of_pic_nums_idc[LIST_1][i] = 3;
      currSlice->abs_diff_pic_num_minus1[LIST_1][i] = 0;
      currSlice->long_term_pic_idx[LIST_1][i] = 0;
    }
  }

  // LIST_0
  num_ref = currSlice->num_ref_idx_active[LIST_0];

  currSlice->poc_ref_pic_reorder_frame(currSlice, p_Vid->listX[LIST_0], num_ref,
    currSlice->reordering_of_pic_nums_idc[LIST_0],
    currSlice->abs_diff_pic_num_minus1[LIST_0],
    currSlice->long_term_pic_idx[LIST_0], LIST_0);
  // reference picture reordering
  reorder_ref_pic_list ( currSlice, p_Vid->listX, p_Vid->listXsize, LIST_0);

  if ( currSlice->slice_type == B_SLICE )
  {
    // LIST_1
    num_ref = currSlice->num_ref_idx_active[LIST_1];

    currSlice->poc_ref_pic_reorder_frame(currSlice, p_Vid->listX[LIST_1], num_ref,
      currSlice->reordering_of_pic_nums_idc[LIST_1],
      currSlice->abs_diff_pic_num_minus1[LIST_1],
      currSlice->long_term_pic_idx[LIST_1], LIST_1);
    // reference picture reordering
    reorder_ref_pic_list ( currSlice, p_Vid->listX, p_Vid->listXsize, LIST_1);
  }
}

