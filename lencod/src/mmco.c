/*!
 *************************************************************************************
 * \file mmco.c
 *
 * \brief
 *    MMCO example operations.
 *
 * \author
 *    Main contributors (see contributors.h for copyright, address and affiliation details)
 *     - Alexis Michael Tourapis                     <alexismt@ieee.org>
 *     - Athanasios Leontaris                        <aleon@dolby.com>
 *************************************************************************************
 */

#include "contributors.h"

#include <ctype.h>
#include <limits.h>
#include "global.h"

#include "image.h"
#include "nalucommon.h"
#include "report.h"


void mmco_long_term(VideoParameters *p_Vid, int current_pic_num)
{
  DecRefPicMarking_t *tmp_drpm,*tmp_drpm2;

  if ( !p_Vid->currentPicture->idr_flag )
  {
    if (p_Vid->dec_ref_pic_marking_buffer!=NULL)
      return;

    if (NULL==(tmp_drpm=(DecRefPicMarking_t*)calloc (1,sizeof (DecRefPicMarking_t)))) 
      no_mem_exit("poc_based_ref_management: tmp_drpm");

    tmp_drpm->Next=NULL;

    tmp_drpm->memory_management_control_operation = 0;

    if (NULL==(tmp_drpm2=(DecRefPicMarking_t*)calloc (1,sizeof (DecRefPicMarking_t)))) 
      no_mem_exit("poc_based_ref_management: tmp_drpm2");
    tmp_drpm2->Next=tmp_drpm;

    tmp_drpm2->memory_management_control_operation = 3;
    tmp_drpm2->long_term_frame_idx = current_pic_num;
    p_Vid->dec_ref_pic_marking_buffer = tmp_drpm2;
  }
  else
  {
    p_Vid->long_term_reference_flag = TRUE;
  }
}

/*!
************************************************************************
* \brief
*    POC-based reference management (FRAME)
************************************************************************
*/

void poc_based_ref_management_frame_pic(DecodedPictureBuffer *p_Dpb, int current_pic_num)
{
  VideoParameters *p_Vid = p_Dpb->p_Vid;
  unsigned i;
  int pic_num = 0;

  int min_poc=INT_MAX;
  DecRefPicMarking_t *tmp_drpm,*tmp_drpm2;

  if (p_Vid->dec_ref_pic_marking_buffer!=NULL)
    return;

  if ( p_Vid->currentPicture->idr_flag )
    return;

  if ((p_Dpb->ref_frames_in_buffer + p_Dpb->ltref_frames_in_buffer)==0)
    return;

  for (i = 0; i < p_Dpb->used_size; i++)
  {
    if (p_Dpb->fs[i]->is_reference  && (!(p_Dpb->fs[i]->is_long_term)) && p_Dpb->fs[i]->poc < min_poc)
    {
      min_poc = p_Dpb->fs[i]->frame->poc ;
      pic_num = p_Dpb->fs[i]->frame->pic_num;
    }
  }

  if (NULL==(tmp_drpm=(DecRefPicMarking_t*)calloc (1,sizeof (DecRefPicMarking_t)))) 
    no_mem_exit("poc_based_ref_management: tmp_drpm");
  tmp_drpm->Next=NULL;

  tmp_drpm->memory_management_control_operation = 0;

  if (NULL==(tmp_drpm2=(DecRefPicMarking_t*)calloc (1,sizeof (DecRefPicMarking_t)))) 
    no_mem_exit("poc_based_ref_management: tmp_drpm2");
  tmp_drpm2->Next=tmp_drpm;

  tmp_drpm2->memory_management_control_operation = 1;

  if ( p_Vid->num_of_layers == 1 )
    tmp_drpm2->difference_of_pic_nums_minus1 = current_pic_num - pic_num - 1;
  else
    tmp_drpm2->difference_of_pic_nums_minus1 = (current_pic_num - pic_num)/2 - 1;

  p_Vid->dec_ref_pic_marking_buffer = tmp_drpm2;
}

/*!
************************************************************************
* \brief
*    POC-based reference management (FIELD)
************************************************************************
*/

void poc_based_ref_management_field_pic(DecodedPictureBuffer *p_Dpb, int current_pic_num)
{
  VideoParameters *p_Vid = p_Dpb->p_Vid;
  unsigned int i; 
  int pic_num1 = 0, pic_num2 = 0;

  int min_poc=INT_MAX;
  DecRefPicMarking_t *tmp_drpm,*tmp_drpm2, *tmp_drpm3;

  if (p_Vid->dec_ref_pic_marking_buffer!=NULL)
    return;

  if ( p_Vid->currentPicture->idr_flag )
    return;

  if ((p_Dpb->ref_frames_in_buffer + p_Dpb->ltref_frames_in_buffer)==0)
    return;

  if ( p_Vid->structure == TOP_FIELD )
  {
    for (i=0; i<p_Dpb->used_size;i++)
    {
      if (p_Dpb->fs[i]->is_reference && (!(p_Dpb->fs[i]->is_long_term)) && p_Dpb->fs[i]->poc < min_poc)
      {      
        min_poc  = p_Dpb->fs[i]->poc;
        pic_num1 = p_Dpb->fs[i]->top_field->pic_num;
        pic_num2 = p_Dpb->fs[i]->bottom_field->pic_num;
      }
    }
  }

  if (NULL==(tmp_drpm=(DecRefPicMarking_t*)calloc (1,sizeof (DecRefPicMarking_t)))) no_mem_exit("poc_based_ref_management_field_pic: tmp_drpm");
  tmp_drpm->Next=NULL;
  tmp_drpm->memory_management_control_operation = 0;

  if ( p_Vid->structure == BOTTOM_FIELD )
  {
    p_Vid->dec_ref_pic_marking_buffer = tmp_drpm;
    return;
  }

  if (NULL==(tmp_drpm2=(DecRefPicMarking_t*)calloc (1,sizeof (DecRefPicMarking_t)))) no_mem_exit("poc_based_ref_management_field_pic: tmp_drpm2");
  tmp_drpm2->Next=tmp_drpm;
  tmp_drpm2->memory_management_control_operation = 1;
  tmp_drpm2->difference_of_pic_nums_minus1 = current_pic_num - pic_num1 - 1;

  if (NULL==(tmp_drpm3=(DecRefPicMarking_t*)calloc (1,sizeof (DecRefPicMarking_t)))) no_mem_exit("poc_based_ref_management_field_pic: tmp_drpm3");
  tmp_drpm3->Next=tmp_drpm2;
  tmp_drpm3->memory_management_control_operation = 1;
  tmp_drpm3->difference_of_pic_nums_minus1 = current_pic_num - pic_num2 - 1;

  p_Vid->dec_ref_pic_marking_buffer = tmp_drpm3;
}


/*!
************************************************************************
* \brief
*    Temporal layer-based reference management (FRAME)
************************************************************************
*/

void tlyr_based_ref_management_frame_pic(VideoParameters *p_Vid, int current_pic_num)
{
  unsigned i, first = 1;
  DecRefPicMarking_t *drpm = NULL, *current_drpm = NULL, *tmp_drpm = NULL;
  struct decoded_picture_buffer *p_Dpb = p_Vid->p_Dpb_layer[0];

  if (p_Vid->dec_ref_pic_marking_buffer!=NULL)
    return;

  if ( p_Vid->currentPicture->idr_flag )
    return;

  if ((p_Dpb->ref_frames_in_buffer + p_Dpb->ltref_frames_in_buffer)==0)
    return;

  for (i = 0; i < p_Dpb->used_size; i++)
  {
    if (p_Dpb->fs[i]->is_reference && (!(p_Dpb->fs[i]->is_long_term)) && p_Dpb->fs[i]->frame->temporal_layer > p_Vid->enc_picture->temporal_layer)
    {
      if (NULL == (tmp_drpm=(DecRefPicMarking_t*)calloc (1,sizeof (DecRefPicMarking_t)))) 
        no_mem_exit("poc_based_ref_management: tmp_drpm2");
      tmp_drpm->memory_management_control_operation = 1;
      tmp_drpm->difference_of_pic_nums_minus1 = current_pic_num - p_Dpb->fs[i]->frame->pic_num - 1;
      
      if (first) 
      {
        drpm = current_drpm = tmp_drpm;
        first = 0;
      }
      else 
      {
        current_drpm->Next = tmp_drpm;
        current_drpm = current_drpm->Next;
      }
    }
  }

  if (first)
    return;

  if (NULL==(tmp_drpm=(DecRefPicMarking_t*)calloc (1,sizeof (DecRefPicMarking_t)))) 
    no_mem_exit("poc_based_ref_management: tmp_drpm");
  tmp_drpm->Next=NULL;

  tmp_drpm->memory_management_control_operation = 0;

  current_drpm->Next=tmp_drpm;

  p_Vid->dec_ref_pic_marking_buffer = drpm;
}

