
/*!
 *************************************************************************************
 * \file mc_prediction.c
 *
 * \brief
 *    Functions for motion compensated prediction
 *
 * \author
 *      Main contributors (see contributors.h for copyright, 
 *                         address and affiliation details)
 *      - Alexis Michael Tourapis  <alexismt@ieee.org>
 *
 *************************************************************************************
 */
#include "global.h"
#include "block.h"
#include "mc_prediction.h"
#include "mbuffer.h"
#include "mb_access.h"
#include "macroblock.h"
#include "memalloc.h"

int allocate_pred_mem(Slice *currSlice)
{
  int alloc_size = 0;
  alloc_size += get_mem2Dpel(&currSlice->tmp_block_l0, MB_BLOCK_SIZE, MB_BLOCK_SIZE);
  alloc_size += get_mem2Dpel(&currSlice->tmp_block_l1, MB_BLOCK_SIZE, MB_BLOCK_SIZE);
  alloc_size += get_mem2Dint(&currSlice->tmp_res, MB_BLOCK_SIZE + 5, MB_BLOCK_SIZE + 5);

  return (alloc_size);
}

void free_pred_mem(Slice *currSlice)
{
  free_mem2Dint(currSlice->tmp_res);
  free_mem2Dpel(currSlice->tmp_block_l0);
  free_mem2Dpel(currSlice->tmp_block_l1);
}

static const int COEF[6] = { 1, -5, 20, 20, -5, 1 };
/*!
 ************************************************************************
 * \brief
 *    block single list prediction
 ************************************************************************
 */
static inline void mc_prediction(imgpel **mb_pred,
                    int ver_block_size, 
                    int hor_block_size,
                    int ioff,
                    imgpel **block)
{
  int jj;

  if (hor_block_size == MB_BLOCK_SIZE)
  {
    memcpy(&(mb_pred[0][ioff]), &(block[0][0]), hor_block_size * ver_block_size * sizeof(imgpel));
  }
  else
  {
    for(jj = 0; jj < ver_block_size; jj++)
    {
      memcpy(&(mb_pred[jj][ioff]), &(block[jj][0]), hor_block_size * sizeof(imgpel));
    }
  }
}

/*!
 ************************************************************************
 * \brief
 *    block single list weighted prediction
 ************************************************************************
 */
static inline void weighted_mc_prediction(imgpel **mb_pred,
                            int ver_block_size, 
                            int hor_block_size,
                            int ioff,
                            imgpel **block, 
                            int wp_scale,
                            int wp_offset,
                            int weight_denom,
                            int color_clip)
{
  int ii, jj;
  imgpel *mpr, *b0;
  
  for(jj=0;jj<ver_block_size;jj++)
  {
    mpr = &mb_pred[jj][ioff];
    b0 = block[jj];
    for(ii=0;ii<hor_block_size;ii++)
      *(mpr++) = (imgpel) iClip1(color_clip, (rshift_rnd((wp_scale * *(b0++)), weight_denom)  + wp_offset ));
  }
}


/*!
 ************************************************************************
 * \brief
 *    block biprediction
 ************************************************************************
 */
static inline void bi_prediction(imgpel **mb_pred,  
                                 imgpel **block_l0, 
                                 imgpel **block_l1,
                                 int ver_block_size, 
                                 int hor_block_size,
                                 int ioff)
{
  int ii, jj;
  imgpel *mpr, *b0, *b1;

  for(jj = 0;jj < ver_block_size;jj++)
  {
    mpr = &mb_pred[jj][ioff];
    b0 = block_l0[jj];
    b1 = block_l1[jj];
    for(ii = 0; ii < hor_block_size;ii++)
      *(mpr++) = (imgpel) rshift_rnd_sf((*(b0++) + *(b1++)), 1);
  }
}

/*!
 ************************************************************************
 * \brief
 *    block weighted biprediction
 ************************************************************************
 */
static inline void weighted_bi_prediction(imgpel **mb_pred, 
                                          imgpel **block_l0, 
                                          imgpel **block_l1,
                                          int ver_block_size, 
                                          int hor_block_size,
                                          int ioff,
                                          int wp_scale_l0,
                                          int wp_scale_l1,
                                          int wp_offset,
                                          int weight_denom,
                                          int color_clip)
{
  int ii, jj;
  imgpel *mpr, *b0, *b1;
  
  for(jj = 0; jj < ver_block_size; jj++)
  {
    mpr = &mb_pred[jj][ioff];    
    b0  = block_l0[jj];
    b1  = block_l1[jj];

    for(ii=0;ii<hor_block_size;ii++)
      *(mpr++) = (imgpel) iClip1(color_clip, (rshift_rnd((wp_scale_l0 * *(b0++) + wp_scale_l1 * *(b1++)), weight_denom) + wp_offset));
  }
}

/*!
 ************************************************************************
 * \brief
 *    Interpolation of 1/4 subpixel
 ************************************************************************
 */ 
void get_block_luma(Macroblock *currMB, ColorPlane pl, StorablePicture *curr_ref, int x_pos, int y_pos, int hor_block_size, int ver_block_size, imgpel **block)
{  
  
  imgpel **cur_imgY = curr_ref->imgY;
  
  int dx = (x_pos & 3), dy = (y_pos & 3);
  int i, j, jj;

  ImageParameters *p_Img = currMB->p_Img;
  StorablePicture *dec_picture = p_Img->dec_picture;
  int shift_x  = dec_picture->size_x;
  int maxold_x = dec_picture->size_x_m1;
  int maxold_y = (dec_picture->motion.mb_field[currMB->mbAddrX]) ? (dec_picture->size_y >> 1) - 1 : dec_picture->size_y_m1;   
  
  if (curr_ref == p_Img->no_reference_picture && p_Img->framepoc < p_Img->recovery_poc)
  {
    printf("list[ref_frame] is equal to 'no reference picture' before RAP\n");

    /* fill the block with sample value 128 */
    for (j = 0; j < ver_block_size; j++)
      for (i = 0; i < hor_block_size; i++)
        block[j][i] = 128;

    return;
  }

  if( IS_INDEPENDENT(p_Img) )
  {
    switch( p_Img->colour_plane_id )
    {
    case    0:
      cur_imgY = curr_ref->imgY;
      break;
    case    1:
      cur_imgY = curr_ref->imgUV[0];
      break;
    case    2:
      cur_imgY = curr_ref->imgUV[1];
      break;
    }
  }
  else if (pl==PLANE_Y)
  {
    cur_imgY = curr_ref->imgY;
  }
  else
  {
    cur_imgY = curr_ref->imgUV[pl-1]; 
  }

  x_pos = x_pos >> 2;
  y_pos = y_pos >> 2;

  if ( (y_pos > 1) && (y_pos < maxold_y - 2 - ver_block_size) && (x_pos > 1) && (x_pos < maxold_x - 2 - hor_block_size))
  {
    cur_imgY = &cur_imgY[ y_pos];
    if (dx == 0 && dy == 0)
    {  /* fullpel position */
      for (j = 0; j < ver_block_size; j++)
      {        
        memcpy(&(block[j][0]), &(cur_imgY[j][x_pos]), hor_block_size * sizeof(imgpel));
      }
    }
    else
    { /* other positions */

      int *tmp_line;
      imgpel *p0, *p1, *p2, *p3, *p4, *p5;
      int    *x0, *x1, *x2, *x3, *x4, *x5;  
      int max_imgpel_value = p_Img->max_imgpel_value_comp[pl];
      imgpel *orig_line, *cur_lineY;

      int result;
      if (dy == 0)
      { /* No vertical interpolation */
        for (j = 0; j < ver_block_size; j++)
        {
          p0 = &cur_imgY[j][x_pos - 2];
          p1 = p0 + 1;
          p2 = p1 + 1;
          p3 = p2 + 1;
          p4 = p3 + 1;
          p5 = p4 + 1;
          orig_line = block[j];

          for (i = 0; i < hor_block_size; i++)
          {        
            result  = (*(p0++) + *(p5++)) * COEF[0]
            + (*(p1++) + *(p4++)) * COEF[1]
            + (*(p2++) + *(p3++)) * COEF[2];

            *orig_line++ = (imgpel) iClip1(max_imgpel_value, ((result + 16)>>5));
          }
        }

        if ((dx & 0x01) == 1)
        {          
          for (j = 0; j < ver_block_size; j++)
          {
            cur_lineY = &(cur_imgY[j][x_pos + (dx >> 1)]);
            orig_line = block[j];
            for (i = 0; i < hor_block_size; i++)
            {
              *orig_line = (imgpel) ((*orig_line + *(cur_lineY++) + 1 ) >> 1);
              orig_line++;
            }
          }
        }
      }
      else if (dx == 0)
      {  /* No horizontal interpolation */        
        p0 = &(cur_imgY[ - 2][x_pos]);
        for (j = 0; j < ver_block_size; j++)
        {                  
          p1 = p0 + shift_x;          
          p2 = p1 + shift_x;
          p3 = p2 + shift_x;
          p4 = p3 + shift_x;
          p5 = p4 + shift_x;
          orig_line = block[j];

          for (i = 0; i < hor_block_size; i++)
          {
            result  = (*(p0++) + *(p5++)) * COEF[0]
            + (*(p1++) + *(p4++)) * COEF[1]
            + (*(p2++) + *(p3++)) * COEF[2];

            *orig_line++ = (imgpel) iClip1(max_imgpel_value, ((result + 16)>>5));
          }
          p0 = p1 - hor_block_size;
        }

        if ((dy&1) == 1)
        {
          jj = (dy >> 1);
          for (j = 0; j < ver_block_size; j++)
          {
            cur_lineY = &(cur_imgY[jj++][x_pos]);
            orig_line = block[j];
            for (i = 0; i < hor_block_size; i++)
            {
              *orig_line = (imgpel) ((*orig_line + *(cur_lineY++) + 1 ) >> 1);
              orig_line++;
            }
          }
        }
      }
      else if (dx == 2)
      {  /* Vertical & horizontal interpolation */
        int **tmp_res = currMB->p_Slice->tmp_res;

        jj = - 2;
        for (j = 0; j < ver_block_size + 5; j++)
        {
          p0 = &cur_imgY[jj++][x_pos - 2];
          p1 = p0 + 1;
          p2 = p1 + 1;
          p3 = p2 + 1;
          p4 = p3 + 1;
          p5 = p4 + 1;          
          tmp_line  = tmp_res[j];

          for (i = 0; i < hor_block_size; i++)
          {        
            *(tmp_line++) = (*(p0++) + *(p5++)) * COEF[0]
            + (*(p1++) + *(p4++)) * COEF[1]
            + (*(p2++) + *(p3++)) * COEF[2];
          }
        }
        
        for (j = 0; j < ver_block_size; j++)
        {
          x0 = tmp_res[j    ];
          x1 = tmp_res[j + 1];
          x2 = tmp_res[j + 2];
          x3 = tmp_res[j + 3];
          x4 = tmp_res[j + 4];
          x5 = tmp_res[j + 5];
          orig_line = block[j];

          for (i = 0; i < hor_block_size; i++)
          {
            result  = (*x0++ + *x5++) * COEF[0]
            + (*x1++ + *x4++) * COEF[1]
            + (*x2++ + *x3++) * COEF[2];

            *(orig_line++) = (imgpel) iClip1(max_imgpel_value, ((result+512)>>10));
          }
        }

        if ((dy&1) == 1)
        {
          jj = 2 + (dy>>1);
          for (j = 0; j < ver_block_size; j++)
          {            
            tmp_line  = tmp_res[jj++];
            orig_line = block[j];
            for (i = 0; i < hor_block_size; i++)
            {
              *orig_line = (imgpel) ((*orig_line + iClip1(max_imgpel_value, ((*(tmp_line++) + 16) >> 5)) + 1 )>> 1);
              orig_line++;
            }
          }
        }
      }
      else if (dy == 2)
      {  /* Horizontal & vertical interpolation */
        int **tmp_res = currMB->p_Slice->tmp_res;
        p0 = &(cur_imgY[ -2][x_pos - 2]);
        for (j = 0; j < ver_block_size; j++)
        {                    
          p1 = p0 + shift_x;
          p2 = p1 + shift_x;
          p3 = p2 + shift_x;
          p4 = p3 + shift_x;
          p5 = p4 + shift_x;
          tmp_line  = tmp_res[j];

          for (i = 0; i < hor_block_size + 5; i++)
          {
            *(tmp_line++)  = (*(p0++) + *(p5++)) * COEF[0]
            + (*(p1++) + *(p4++)) * COEF[1]
            + (*(p2++) + *(p3++)) * COEF[2];
          }
          p0 = p1 - (hor_block_size + 5);
        }

        for (j = 0; j < ver_block_size; j++)
        {
          orig_line = block[j];
          x0 = tmp_res[j];
          x1 = x0 + 1;
          x2 = x1 + 1;
          x3 = x2 + 1;
          x4 = x3 + 1;
          x5 = x4 + 1;

          for (i = 0; i < hor_block_size; i++)
          {
            result  = (*(x0++) + *(x5++)) * COEF[0]
            + (*(x1++) + *(x4++)) * COEF[1]
            + (*(x2++) + *(x3++)) * COEF[2];

            *(orig_line++) = (imgpel) iClip1(max_imgpel_value, ((result + 512)>>10));
          }
        }

        if ((dx&1) == 1)
        {
          for (j = 0; j < ver_block_size; j++)
          {
            tmp_line  = &tmp_res[j][2 + (dx>>1)];
            orig_line = block[j];
            for (i = 0; i < hor_block_size; i++)
            {
              *orig_line = (imgpel) ((*orig_line + iClip1(max_imgpel_value, ((*(tmp_line++) + 16)>>5))+1)>>1);
              orig_line ++;
            }
          }
        }
      }
      else
      {  /* Diagonal interpolation */
        jj = (dy == 1 ? 0 : 1);

        for (j = 0; j < ver_block_size; j++)
        {
          p0 = &cur_imgY[jj++][x_pos - 2];
          p1 = p0 + 1;
          p2 = p1 + 1;
          p3 = p2 + 1;
          p4 = p3 + 1;
          p5 = p4 + 1;

          orig_line = block[j];

          for (i = 0; i < hor_block_size; i++)
          {        
            result  = (*(p0++) + *(p5++)) * COEF[0]
            + (*(p1++) + *(p4++)) * COEF[1]
            + (*(p2++) + *(p3++)) * COEF[2];

            *(orig_line++) = (imgpel) iClip1(max_imgpel_value, ((result + 16)>>5));
          }
        }

        p0 = &(cur_imgY[-2][(dx == 1 ? x_pos : x_pos + 1)]);
        for (j = 0; j < ver_block_size; j++)
        {        
          p1 = p0 + shift_x;
          p2 = p1 + shift_x;
          p3 = p2 + shift_x;
          p4 = p3 + shift_x;
          p5 = p4 + shift_x;
          orig_line = block[j];

          for (i = 0; i < hor_block_size; i++)
          {
            result  = (*(p0++) + *(p5++)) * COEF[0]
            + (*(p1++) + *(p4++)) * COEF[1]
            + (*(p2++) + *(p3++)) * COEF[2];

            *orig_line = (imgpel) ((*orig_line + iClip1(max_imgpel_value, ((result + 16) >> 5)) + 1) >> 1);
            orig_line++;
          }
          p0 = p1 - hor_block_size ;
        }      
      }
    }
  }
  else // unsafe positions
  {
    int result;
    imgpel *orig_line, *cur_lineY;
    if (dx == 0 && dy == 0)
    {  /* fullpel position */
      for (j = 0; j < ver_block_size; j++)
      {
        cur_lineY = cur_imgY[iClip3(0, maxold_y, y_pos + j)];
        orig_line = block[j];
        for (i = 0; i < hor_block_size; i++)
        {
          *(orig_line++) = cur_lineY[iClip3(0, maxold_x, x_pos + i )];
        }
      }
    }
    else
    { /* other positions */

      int max_imgpel_value = p_Img->max_imgpel_value_comp[pl];
      int *tmp_line;            
      int ipos_m2, ipos_m1, ipos, ipos_p1, ipos_p2, ipos_p3;

      if (dy == 0)
      { /* No vertical interpolation */
        int tmp_pos = x_pos - 2;
        for (i = 0; i < hor_block_size; i++)
        {        
          ipos_m2 = iClip3(0, maxold_x, tmp_pos++);
          ipos_m1 = iClip3(0, maxold_x, tmp_pos++);
          ipos    = iClip3(0, maxold_x, tmp_pos++);
          ipos_p1 = iClip3(0, maxold_x, tmp_pos++);
          ipos_p2 = iClip3(0, maxold_x, tmp_pos++);
          ipos_p3 = iClip3(0, maxold_x, tmp_pos  );
          tmp_pos -= 4;          

          for (j = 0; j < ver_block_size; j++)
          {
            cur_lineY = cur_imgY[iClip3(0,maxold_y,y_pos+j)];

            result  = (cur_lineY[ipos_m2] + cur_lineY[ipos_p3]) * COEF[0];
            result += (cur_lineY[ipos_m1] + cur_lineY[ipos_p2]) * COEF[1];
            result += (cur_lineY[ipos   ] + cur_lineY[ipos_p1]) * COEF[2];

            block[j][i] = (imgpel) iClip1(max_imgpel_value, ((result + 16)>>5));
          }
        }

        if ((dx&1) == 1)
        {
          for (j = 0; j < ver_block_size; j++)
          {
            cur_lineY = cur_imgY[iClip3(0,maxold_y,y_pos+j)];
            orig_line = block[j];
            for (i = 0; i < hor_block_size; i++)
            {
              *orig_line = (imgpel) ((*orig_line + cur_lineY[iClip3(0, maxold_x, x_pos + i + (dx >> 1))] + 1 ) >> 1);
              orig_line++;
            }
          }
        }
      }
      else if (dx == 0)
      {  /* No horizontal interpolation */
        int pres_x;
        imgpel *p0, *p1, *p2, *p3, *p4, *p5;
        int tmp_pos = y_pos - 2;

        for (j = 0; j < ver_block_size; j++)
        { 
          p0 = cur_imgY[iClip3(0, maxold_y, tmp_pos++)];
          p1 = cur_imgY[iClip3(0, maxold_y, tmp_pos++)];
          p2 = cur_imgY[iClip3(0, maxold_y, tmp_pos++)];
          p3 = cur_imgY[iClip3(0, maxold_y, tmp_pos++)];
          p4 = cur_imgY[iClip3(0, maxold_y, tmp_pos++)];
          p5 = cur_imgY[iClip3(0, maxold_y, tmp_pos  )];

          tmp_pos -= 4;
          orig_line = block[j];

          for (i = 0; i < hor_block_size; i++)
          {
            pres_x = iClip3(0,maxold_x,x_pos+i);

            result  = (p0[pres_x] + p5[pres_x]) * COEF[0];
            result += (p1[pres_x] + p4[pres_x]) * COEF[1];
            result += (p2[pres_x] + p3[pres_x]) * COEF[2];
            *(orig_line++) = (imgpel) iClip1(max_imgpel_value, ((result+16)>>5));
          }
        }

        if ((dy&1) == 1)
        {
          for (j = 0; j < ver_block_size; j++)
          {
            cur_lineY = cur_imgY[iClip3(0,maxold_y,y_pos+j+(dy>>1))];
            orig_line = block[j];
            for (i = 0; i < hor_block_size; i++)
            {
              *orig_line = (imgpel) ((*orig_line + cur_lineY[iClip3(0, maxold_x, x_pos + i)] + 1 ) >> 1);
              orig_line++;
            }
          }
        }
      }
      else if (dx == 2)
      {  /* Vertical & horizontal interpolation */
        int    *x0, *x1, *x2, *x3, *x4, *x5;  
        int **tmp_res = currMB->p_Slice->tmp_res;
        int tmp_pos = x_pos - 2;
        for (i = 0; i < hor_block_size; i++)
        {        
          ipos_m2 = iClip3(0, maxold_x, tmp_pos++);
          ipos_m1 = iClip3(0, maxold_x, tmp_pos++);
          ipos    = iClip3(0, maxold_x, tmp_pos++);
          ipos_p1 = iClip3(0, maxold_x, tmp_pos++);
          ipos_p2 = iClip3(0, maxold_x, tmp_pos++);
          ipos_p3 = iClip3(0, maxold_x, tmp_pos  );
          tmp_pos -= 4;

          for (j = 0; j < ver_block_size + 5; j++)
          {
            cur_lineY = cur_imgY[iClip3(0,maxold_y,y_pos + j - 2)];

            tmp_res[j][i]  = (cur_lineY[ipos_m2] + cur_lineY[ipos_p3]) * COEF[0];
            tmp_res[j][i] += (cur_lineY[ipos_m1] + cur_lineY[ipos_p2]) * COEF[1];
            tmp_res[j][i] += (cur_lineY[ipos   ] + cur_lineY[ipos_p1]) * COEF[2];
          }
        }

        for (j = 0; j < ver_block_size; j++)
        {
          x0 = tmp_res[j    ];
          x1 = tmp_res[j + 1];
          x2 = tmp_res[j + 2];
          x3 = tmp_res[j + 3];
          x4 = tmp_res[j + 4];
          x5 = tmp_res[j + 5];
          orig_line = block[j];

          for (i = 0; i < hor_block_size; i++)
          {
            result  = (*x0++ + *x5++) * COEF[0]
            + (*x1++ + *x4++) * COEF[1]
            + (*x2++ + *x3++) * COEF[2];

            *(orig_line++) = (imgpel) iClip1(max_imgpel_value, ((result+512)>>10));
          }
        }

        if ((dy&1) == 1)
        {
          for (j = 0; j < ver_block_size; j++)
          {
            tmp_line  = tmp_res[j + 2 + (dy>>1)];            
            orig_line = block[j];
            for (i = 0; i < hor_block_size; i++)
            {
              *orig_line = (imgpel) ((*orig_line + iClip1(max_imgpel_value, ((*(tmp_line++) + 16) >> 5)) + 1 )>>1);
              orig_line++;
            }
          }
        }
      }
      else if (dy == 2)
      {  /* Horizontal & vertical interpolation */
        int    *x0, *x1, *x2, *x3, *x4, *x5;  
        imgpel *p0, *p1, *p2, *p3, *p4, *p5;
        int pres_x;
        int **tmp_res = currMB->p_Slice->tmp_res;

        int tmp_pos = y_pos - 2;

        for (j = 0; j < ver_block_size; j++)
        {
          p0 = cur_imgY[iClip3(0, maxold_y, tmp_pos++)];
          p1 = cur_imgY[iClip3(0, maxold_y, tmp_pos++)];
          p2 = cur_imgY[iClip3(0, maxold_y, tmp_pos++)];
          p3 = cur_imgY[iClip3(0, maxold_y, tmp_pos++)];
          p4 = cur_imgY[iClip3(0, maxold_y, tmp_pos++)];
          p5 = cur_imgY[iClip3(0, maxold_y, tmp_pos  )];

          tmp_pos -= 4;

          for (i = 0; i < hor_block_size + 5; i++)
          {
            pres_x = iClip3(0,maxold_x, x_pos + i - 2);
            result  = (p0[pres_x] + p5[pres_x])*COEF[0]
            + (p1[pres_x] + p4[pres_x])*COEF[1]
            + (p2[pres_x] + p3[pres_x])*COEF[2];
            tmp_res[j][i] = result;
          }
        }

        for (j = 0; j < ver_block_size; j++)
        {
          orig_line = block[j];
          tmp_line  = tmp_res[j];
          x0 = tmp_res[j];
          x1 = x0 + 1;
          x2 = x1 + 1;
          x3 = x2 + 1;
          x4 = x3 + 1;
          x5 = x4 + 1;

          for (i = 0; i < hor_block_size; i++)
          {
            result  = (*(x0++) + *(x5++)) * COEF[0]
            + (*(x1++) + *(x4++)) * COEF[1]
            + (*(x2++) + *(x3++)) * COEF[2];

            *(orig_line++) = (imgpel) iClip1(max_imgpel_value, ((result + 512)>>10));
          }
        }

        if ((dx&1) == 1)
        {
          for (j = 0; j < ver_block_size; j++)
          {
            tmp_line  = &tmp_res[j][2 + (dx>>1)];
            orig_line = block[j];            
            for (i = 0; i < hor_block_size; i++)
            {
              *orig_line = (imgpel) ((*orig_line + iClip1(max_imgpel_value, ((*(tmp_line++) + 16)>>5)) + 1)>> 1);
              orig_line++;
            }
          }
        }
      }
      else
      {  /* Diagonal interpolation */
        int pres_x;
        imgpel *p0, *p1, *p2, *p3, *p4, *p5;

        int tmp_pos = x_pos - 2;

        for (i = 0; i < hor_block_size; i++)
        {
          ipos_m2 = iClip3(0, maxold_x, tmp_pos++);
          ipos_m1 = iClip3(0, maxold_x, tmp_pos++);
          ipos    = iClip3(0, maxold_x, tmp_pos++);
          ipos_p1 = iClip3(0, maxold_x, tmp_pos++);
          ipos_p2 = iClip3(0, maxold_x, tmp_pos++);
          ipos_p3 = iClip3(0, maxold_x, tmp_pos  );
          tmp_pos -= 4;

          for (j = 0; j < ver_block_size; j++)
          {
            cur_lineY = cur_imgY[iClip3(0,maxold_y,(dy == 1 ? y_pos+j : y_pos+j+1))];

            result  = (cur_lineY[ipos_m2] + cur_lineY[ipos_p3]) * COEF[0];
            result += (cur_lineY[ipos_m1] + cur_lineY[ipos_p2]) * COEF[1];
            result += (cur_lineY[ipos   ] + cur_lineY[ipos_p1]) * COEF[2];

            block[j][i] = (imgpel) iClip1(max_imgpel_value, ((result+16)>>5));
          }
        }

        tmp_pos = y_pos - 2;
        for (j = 0; j < ver_block_size; j++)
        {      
          p0 = cur_imgY[iClip3(0, maxold_y, tmp_pos++)];
          p1 = cur_imgY[iClip3(0, maxold_y, tmp_pos++)];
          p2 = cur_imgY[iClip3(0, maxold_y, tmp_pos++)];
          p3 = cur_imgY[iClip3(0, maxold_y, tmp_pos++)];
          p4 = cur_imgY[iClip3(0, maxold_y, tmp_pos++)];
          p5 = cur_imgY[iClip3(0, maxold_y, tmp_pos  )];

          tmp_pos -= 4;
          orig_line = block[j];

          for (i = 0; i < hor_block_size; i++)
          {
            pres_x = dx == 1 ? x_pos+i : x_pos+i+1;
            pres_x = iClip3(0, maxold_x, pres_x);

            result  = (p0[pres_x] + p5[pres_x]) * COEF[0];
            result += (p1[pres_x] + p4[pres_x]) * COEF[1];
            result += (p2[pres_x] + p3[pres_x]) * COEF[2];

            *orig_line = (imgpel) ((*orig_line + iClip1(max_imgpel_value, ((result+16)>>5)) + 1 ) >> 1);
            orig_line++;
          }
        }      
      }
    }
  }
}

void get_block_chroma(Macroblock *currMB, int uv, StorablePicture *curr_ref, int x_pos, int y_pos, int hor_block_size, int ver_block_size, imgpel **block)
{
  ImageParameters *p_Img = currMB->p_Img;
  int subpel_x    = p_Img->mb_cr_size_x == 8 ? 7 : 3;
  int subpel_y    = p_Img->mb_cr_size_y == 8 ? 7 : 3;
  int shiftpel_x  = p_Img->mb_cr_size_x == 8 ? 3 : 2;
  int shiftpel_y  = p_Img->mb_cr_size_y == 8 ? 3 : 2;

  int dx = (x_pos & subpel_x);
  int dy = (y_pos & subpel_y);


  int i, j;
  StorablePicture *dec_picture = p_Img->dec_picture;
  int maxold_x = dec_picture->size_x_cr_m1;
  int maxold_y = (dec_picture->motion.mb_field[p_Img->current_mb_nr]) ? (dec_picture->size_y_cr >> 1) - 1 : dec_picture->size_y_cr_m1;
  
  imgpel **cur_img, *blk_line;


  if (curr_ref == p_Img->no_reference_picture && p_Img->framepoc < p_Img->recovery_poc)
  {
    printf("list[ref_frame] is equal to 'no reference picture' before RAP\n");

    /* fill the block with sample value 128 */
    for (j = 0; j < ver_block_size; j++)
      for (i = 0; i < hor_block_size; i++)
        block[j][i] = 128;
    return;
  }

  cur_img = curr_ref->imgUV[uv];

  x_pos = x_pos >> shiftpel_x;
  y_pos = y_pos >> shiftpel_y;

  if ((y_pos >= 0) && (y_pos < maxold_y - ver_block_size) && (x_pos >= 0) && (x_pos < maxold_x - hor_block_size))
  {
    if (dx == 0 && dy == 0)
    {  /* fullpel position */
      for (j = 0; j < ver_block_size; j++)
      {        
        memcpy(&(block[j][0]), &(cur_img[ y_pos + j ][x_pos]), hor_block_size * sizeof(imgpel));
      }
    }
    else if (dx == 0)
    { 
      int    max_imgpel_value = p_Img->max_imgpel_value_comp[uv + 1];

      int total_scale = shiftpel_x + shiftpel_y;
      int dxcur = (subpel_x + 1 - dx);
      int dycur = (subpel_y + 1 - dy);

      int w00 = dxcur * dycur;
      int w01 = dxcur * dy;
      imgpel *cur_line, *cur_line_p1;
      int result;

      for (j = 0; j < ver_block_size; j++)
      {
        cur_line    = &cur_img[y_pos + j    ][x_pos];
        cur_line_p1 = &cur_img[y_pos + j + 1][x_pos];
        blk_line = block[j];

        for (i = 0; i < hor_block_size; i++)
        {
          result = (w00 * *cur_line++ + w01 * *cur_line_p1++);
          *(blk_line++) = (imgpel) iClip1(max_imgpel_value, rshift_rnd_sf(result, total_scale));
        }
      }
    }
    else if (dy == 0)
    { 
      int    max_imgpel_value = p_Img->max_imgpel_value_comp[uv + 1];
      int total_scale = shiftpel_x + shiftpel_y;
      int dxcur = (subpel_x + 1 - dx);
      int dycur = (subpel_y + 1 - dy);

      int w00 = dxcur * dycur;
      int w10 = dx * dycur;
      imgpel *cur_line, *cur_line_p1;
      int result;

      for (j = 0; j < ver_block_size; j++)
      {
        cur_line    = &cur_img[y_pos + j][x_pos];
        cur_line_p1 = cur_line + 1;
        blk_line = block[j];

        for (i = 0; i < hor_block_size; i++)
        {
          result = (w00 * *cur_line++ + w10 * *cur_line_p1++);
          *(blk_line++) = (imgpel) iClip1(max_imgpel_value, rshift_rnd_sf(result, total_scale));
        }
      }
    }
    else
    { /* other positions */
      int    max_imgpel_value = p_Img->max_imgpel_value_comp[uv + 1];
      int total_scale = shiftpel_x + shiftpel_y;
      int dxcur = (subpel_x + 1 - dx);
      int dycur = (subpel_y + 1 - dy);

      int w00 = dxcur * dycur;
      int w01 = dxcur * dy;
      int w10 = dx * dycur;
      int w11 = dx * dy;
      imgpel *cur_line, *cur_line_p1;
      int result;

      for (j = 0; j < ver_block_size; j++)
      {
        cur_line    = &cur_img[y_pos + j    ][x_pos];
        cur_line_p1 = &cur_img[y_pos + j + 1][x_pos];
        blk_line = block[j];

        for (i = 0; i < hor_block_size; i++)
        {
          result  = (w00 * *(cur_line++) + w01 * *(cur_line_p1++));
          result += (w10 * *(cur_line  ) + w11 * *(cur_line_p1  ));
          *(blk_line++) = (imgpel) iClip1(max_imgpel_value, rshift_rnd_sf(result, total_scale));
        }
      }
    }
  }
  else // unsafe positions
  {
    int    max_imgpel_value = p_Img->max_imgpel_value_comp[uv + 1];
    int total_scale = shiftpel_x + shiftpel_y;
    int dxcur = (subpel_x + 1 - dx);
    int dycur = (subpel_y + 1 - dy);

    int w00 = dxcur * dycur;
    int w01 = dxcur * dy;
    int w10 = dx * dycur;
    int w11 = dx * dy;
    imgpel *cur_line, *cur_line_p1;
    int result;

    if (dx == 0 && dy == 0)
    {  /* fullpel position */
      for (j = 0; j < ver_block_size; j++)
      {
        cur_line = cur_img[iClip3(0, maxold_y, y_pos + j)];
        blk_line = block[j];
        for (i = 0; i < hor_block_size; i++)
        {
          *(blk_line++) = cur_line[iClip3(0, maxold_x, x_pos + i )];
        }
      }
    }
    else if (dx == 0)
    { 
      int tmp_pos, ipos;
      for (j = 0; j < ver_block_size; j++)
      {
        tmp_pos = x_pos;
        cur_line    = cur_img[iClip3(0, maxold_y, y_pos + j)];
        cur_line_p1 = cur_img[iClip3(0, maxold_y, y_pos + j + 1)];

        blk_line = block[j];

        for (i = 0; i < hor_block_size; i++)
        {
          ipos    = iClip3(0, maxold_x, tmp_pos++);

          result = (w00 * cur_line[ipos] + w01 * cur_line_p1[ipos]);
          *(blk_line++) = (imgpel) iClip1(max_imgpel_value, rshift_rnd_sf(result, total_scale));
        }
      }      
    }
    else if (dy == 0)
    { 
      int tmp_pos, ipos, ipos_p1;
      for (j = 0; j < ver_block_size; j++)
      {
        tmp_pos = x_pos;
        cur_line    = cur_img[iClip3(0, maxold_y, y_pos + j)];        
        blk_line = block[j];

        for (i = 0; i < hor_block_size; i++)
        {
          ipos    = iClip3(0, maxold_x, tmp_pos++);
          ipos_p1 = iClip3(0, maxold_x, tmp_pos  );

          result = (w00 * cur_line[ipos   ] + w10 * cur_line[ipos_p1]);
          *(blk_line++) = (imgpel)iClip1(max_imgpel_value, rshift_rnd_sf(result, total_scale));
        }
      }      
    }
    else
    { /* other positions */ 
      int tmp_pos, ipos, ipos_p1;
      for (j = 0; j < ver_block_size; j++)
      {
        cur_line    = cur_img[iClip3(0, maxold_y, y_pos + j)];
        cur_line_p1 = cur_img[iClip3(0, maxold_y, y_pos + j + 1)];
        tmp_pos = x_pos;
        blk_line = block[j];

        for (i = 0; i < hor_block_size; i++)
        {
          ipos    = iClip3(0, maxold_x, tmp_pos++);
          ipos_p1 = iClip3(0, maxold_x, tmp_pos  );

          result = (
            w00 * cur_line   [ipos   ] + 
            w10 * cur_line   [ipos_p1] +
            w01 * cur_line_p1[ipos   ] +
            w11 * cur_line_p1[ipos_p1]);
          *(blk_line++) = (imgpel) iClip1(max_imgpel_value, rshift_rnd_sf(result, total_scale));
        }
      }      
    }
  }
}


void intra_cr_decoding(Macroblock *currMB, int yuv, int smb)
{
  ImageParameters *p_Img = currMB->p_Img;
  Slice *currSlice = currMB->p_Slice;
  StorablePicture *dec_picture = p_Img->dec_picture;
  imgpel **curUV;
  int uv;
  int b8,b4;
  int ioff, joff, ii, jj;

  for(uv = 0; uv < 2; uv++)
  {
    currMB->itrans_4x4 = (currMB->is_lossless == FALSE) ? itrans4x4 : itrans4x4_ls;

    curUV = dec_picture->imgUV[uv];
    intrapred_chroma(currMB, uv);

    if (!smb && (currMB->cbp >> 4))
    {
      for (b8 = 0; b8 < (p_Img->num_uv_blocks); b8++)
      {
        for(b4 = 0; b4 < 4; b4++)
        {
          joff = subblk_offset_y[yuv][b8][b4];          
          ioff = subblk_offset_x[yuv][b8][b4];          

          currMB->itrans_4x4(currMB, (ColorPlane) (uv + 1), ioff, joff);

          copy_image_data_4x4(&curUV[currMB->pix_c_y + joff], &(currSlice->mb_rec[uv + 1][joff]), currMB->pix_c_x + ioff, ioff);
        }
      }
    }
    else if ((currMB->cbp >> 4) == 0)
    {      
      for (b8 = 0; b8 < (p_Img->num_uv_blocks); b8++)
      {
        for(b4 = 0; b4 < 4; b4++)
        {
          joff = subblk_offset_y[yuv][b8][b4];
          ioff = subblk_offset_x[yuv][b8][b4];          

          copy_image_data_4x4(&curUV[currMB->pix_c_y + joff], &(currSlice->mb_pred[uv + 1][joff]), currMB->pix_c_x + ioff, ioff);
        }
      }
    }
    else
    {
      itrans_sp_cr(currMB, uv);

      for (joff  = 0; joff < 8; joff += 4)
      {
        for(ioff = 0; ioff < 8;ioff+=4)
        {          
          currMB->itrans_4x4(currMB, (ColorPlane) (uv + 1), ioff, joff);

          for(jj = joff; jj < joff + 4; jj++)
          {
            for(ii = ioff; ii < ioff + 4; ii++)
            {
              curUV[currMB->pix_c_y+jj][ii + currMB->pix_c_x] = (imgpel) currSlice->mb_rres[uv+1][jj][ii];
            }
          }
        }
      }
    }
  }
}

void prepare_direct_params(Macroblock *currMB, StorablePicture *dec_picture, short pmvl0[2], short pmvl1[2],char *l0_rFrame, char *l1_rFrame)
{
  ImageParameters *p_Img = currMB->p_Img;
  Slice *currSlice = currMB->p_Slice;
  char l0_rFrameL, l0_rFrameU, l0_rFrameUR;
  char l1_rFrameL, l1_rFrameU, l1_rFrameUR;
  PicMotionParams *motion = &dec_picture->motion;
  
  PixelPos mb[4];

  get_neighbors(currMB, mb, 0, 0, 16);

  if (!currSlice->MbaffFrameFlag)
  {
    l0_rFrameL  = (char) (mb[0].available ? motion->ref_idx[LIST_0][mb[0].pos_y][mb[0].pos_x] : -1);
    l0_rFrameU  = (char) (mb[1].available ? motion->ref_idx[LIST_0][mb[1].pos_y][mb[1].pos_x] : -1);
    l0_rFrameUR = (char) (mb[2].available ? motion->ref_idx[LIST_0][mb[2].pos_y][mb[2].pos_x] : -1);

    l1_rFrameL  = (char) (mb[0].available ? motion->ref_idx[LIST_1][mb[0].pos_y][mb[0].pos_x] : -1);
    l1_rFrameU  = (char) (mb[1].available ? motion->ref_idx[LIST_1][mb[1].pos_y][mb[1].pos_x] : -1);
    l1_rFrameUR = (char) (mb[2].available ? motion->ref_idx[LIST_1][mb[2].pos_y][mb[2].pos_x] : -1);
  }
  else
  {
    if (currMB->mb_field)
    {
      l0_rFrameL = (char) (mb[0].available 
        ? p_Img->mb_data[mb[0].mb_addr].mb_field  || motion->ref_idx[LIST_0][mb[0].pos_y][mb[0].pos_x] < 0
        ? motion->ref_idx[LIST_0][mb[0].pos_y][mb[0].pos_x] 
        : motion->ref_idx[LIST_0][mb[0].pos_y][mb[0].pos_x] * 2: -1);

      l0_rFrameU = (char) (mb[1].available 
        ? p_Img->mb_data[mb[1].mb_addr].mb_field || motion->ref_idx[LIST_0][mb[1].pos_y][mb[1].pos_x] < 0
        ? motion->ref_idx[LIST_0][mb[1].pos_y][mb[1].pos_x] 
        : motion->ref_idx[LIST_0][mb[1].pos_y][mb[1].pos_x] * 2: -1);

       l0_rFrameUR = (char) (mb[2].available 
         ? p_Img->mb_data[mb[2].mb_addr].mb_field || motion->ref_idx[LIST_0][mb[2].pos_y][mb[2].pos_x] < 0 
         ? motion->ref_idx[LIST_0][mb[2].pos_y][mb[2].pos_x] 
         : motion->ref_idx[LIST_0][mb[2].pos_y][mb[2].pos_x] * 2: -1);

       l1_rFrameL = (char) (mb[0].available 
         ? p_Img->mb_data[mb[0].mb_addr].mb_field || motion->ref_idx[LIST_1][mb[0].pos_y][mb[0].pos_x]  < 0 
         ? motion->ref_idx[LIST_1][mb[0].pos_y][mb[0].pos_x] 
         : motion->ref_idx[LIST_1][mb[0].pos_y][mb[0].pos_x] * 2: -1);

       l1_rFrameU = (char) (mb[1].available 
         ? p_Img->mb_data[mb[1].mb_addr].mb_field || motion->ref_idx[LIST_1][mb[1].pos_y][mb[1].pos_x]  < 0 
         ? motion->ref_idx[LIST_1][mb[1].pos_y][mb[1].pos_x] 
         : motion->ref_idx[LIST_1][mb[1].pos_y][mb[1].pos_x] * 2: -1);

       l1_rFrameUR = (char) (mb[2].available 
         ? p_Img->mb_data[mb[2].mb_addr].mb_field || motion->ref_idx[LIST_1][mb[2].pos_y][mb[2].pos_x] < 0
         ? motion->ref_idx[LIST_1][mb[2].pos_y][mb[2].pos_x] 
         : motion->ref_idx[LIST_1][mb[2].pos_y][mb[2].pos_x] * 2: -1);
    }
    else
    {
      l0_rFrameL = (char) (mb[0].available 
        ? p_Img->mb_data[mb[0].mb_addr].mb_field || motion->ref_idx[LIST_0][mb[0].pos_y][mb[0].pos_x]  < 0 
        ? motion->ref_idx[LIST_0][mb[0].pos_y][mb[0].pos_x] >> 1 
        : motion->ref_idx[LIST_0][mb[0].pos_y][mb[0].pos_x]: -1);

      l0_rFrameU = (char) (mb[1].available 
        ? p_Img->mb_data[mb[1].mb_addr].mb_field || motion->ref_idx[LIST_0][mb[1].pos_y][mb[1].pos_x] < 0 
        ? motion->ref_idx[LIST_0][mb[1].pos_y][mb[1].pos_x] >> 1 
        : motion->ref_idx[LIST_0][mb[1].pos_y][mb[1].pos_x] : -1);

      l0_rFrameUR = (char) (mb[2].available 
        ? p_Img->mb_data[mb[2].mb_addr].mb_field || motion->ref_idx[LIST_0][mb[2].pos_y][mb[2].pos_x] < 0 
        ? motion->ref_idx[LIST_0][mb[2].pos_y][mb[2].pos_x] >> 1 
        : motion->ref_idx[LIST_0][mb[2].pos_y][mb[2].pos_x] : -1);

      l1_rFrameL = (char) (mb[0].available 
        ? p_Img->mb_data[mb[0].mb_addr].mb_field || motion->ref_idx[LIST_1][mb[0].pos_y][mb[0].pos_x] < 0 
        ? motion->ref_idx[LIST_1][mb[0].pos_y][mb[0].pos_x] >> 1 
        : motion->ref_idx[LIST_1][mb[0].pos_y][mb[0].pos_x] : -1);

      l1_rFrameU = (char) (mb[1].available 
        ? p_Img->mb_data[mb[1].mb_addr].mb_field || motion->ref_idx[LIST_1][mb[1].pos_y][mb[1].pos_x] < 0 
        ? motion->ref_idx[LIST_1][mb[1].pos_y][mb[1].pos_x] >> 1 
        : motion->ref_idx[LIST_1][mb[1].pos_y][mb[1].pos_x] : -1);

      l1_rFrameUR = (char) (mb[2].available 
        ? p_Img->mb_data[mb[2].mb_addr].mb_field || motion->ref_idx[LIST_1][mb[2].pos_y][mb[2].pos_x] < 0 
        ? motion->ref_idx[LIST_1][mb[2].pos_y][mb[2].pos_x] >> 1
        : motion->ref_idx[LIST_1][mb[2].pos_y][mb[2].pos_x] : -1);
    }
  }

  *l0_rFrame = (char) ((l0_rFrameL >= 0 && l0_rFrameU >= 0)  ? imin(l0_rFrameL,l0_rFrameU) : imax(l0_rFrameL,l0_rFrameU));
  *l0_rFrame = (char) ((*l0_rFrame >= 0 && l0_rFrameUR >= 0) ? imin(*l0_rFrame,l0_rFrameUR): imax(*l0_rFrame,l0_rFrameUR));

  *l1_rFrame = (char) ((l1_rFrameL >= 0 && l1_rFrameU >= 0)  ? imin(l1_rFrameL,l1_rFrameU) : imax(l1_rFrameL,l1_rFrameU));
  *l1_rFrame = (char) ((*l1_rFrame >= 0 && l1_rFrameUR >= 0) ? imin(*l1_rFrame,l1_rFrameUR): imax(*l1_rFrame,l1_rFrameUR));

  if (*l0_rFrame >=0)
    currMB->GetMVPredictor (currMB, mb, pmvl0, *l0_rFrame, motion->ref_idx[LIST_0], motion->mv[LIST_0], 0, 0, 16, 16);

  if (*l1_rFrame >=0)
    currMB->GetMVPredictor (currMB, mb, pmvl1, *l1_rFrame, motion->ref_idx[LIST_1], motion->mv[LIST_1], 0, 0, 16, 16);
}

void check_motion_vector_range(ImageParameters *p_Img, short mv_x, short mv_y)
{
  if (mv_x > 8191 || mv_x < -8192)
  {
    fprintf(stderr,"WARNING! Horizontal motion vector %d is out of allowed range {-8192, 8191} in picture %d, macroblock %d\n", mv_x, p_Img->number, p_Img->current_mb_nr);
    //error("invalid stream: too big horizontal motion vector", 500);
  }

  if (mv_y > (p_Img->max_mb_vmv_r - 1) || mv_y < (-p_Img->max_mb_vmv_r))
  {
    fprintf(stderr,"WARNING! Vertical motion vector %d is out of allowed range {%d, %d} in picture %d, macroblock %d\n", mv_y, (-p_Img->max_mb_vmv_r), (p_Img->max_mb_vmv_r - 1), p_Img->number, p_Img->current_mb_nr);
    //error("invalid stream: too big vertical motion vector", 500);
  }
}

void perform_mc(Macroblock *currMB, ColorPlane pl, StorablePicture *dec_picture, int pred_dir, int i, int j, int list_offset, int block_size_x, int block_size_y, int curr_mb_field)
{
  ImageParameters *p_Img = currMB->p_Img;  
    seq_parameter_set_rbsp_t *active_sps = p_Img->active_sps;

  Slice *currSlice = currMB->p_Slice;
  int vec1_x=0, vec1_y=0, vec2_x=0, vec2_y=0;
  int vec1_y_cr = 0, vec2_y_cr = 0;
  int alpha_l0, alpha_l1, wp_offset;
  static const int mv_mul = 16; // 4 * 4
  int        max_imgpel_value = p_Img->max_imgpel_value_comp[pl];
  
  int i4   = currMB->block_x + i;
  int j4   = currMB->block_y + j;
  int ioff = (i << 2);
  int joff = (j << 2);         
  
  assert (pred_dir<=2);

  if (pred_dir != 2)
  {
    //===== Single List Prediction =====
    short       ref_idx = dec_picture->motion.ref_idx[pred_dir][j4][i4];
    short       ref_idx_wp = ref_idx;
    short      *mv_array = dec_picture->motion.mv[pred_dir][j4][i4];
    StorablePicture *list = p_Img->listX[list_offset + pred_dir][ref_idx];

    check_motion_vector_range(p_Img, mv_array[0], mv_array[1]);

    vec1_x = i4 * mv_mul + mv_array[0];
    vec1_y = (currMB->block_y_aff + j) * mv_mul + mv_array[1];

    get_block_luma (currMB, pl, list, vec1_x, vec1_y, block_size_x, block_size_y, currSlice->tmp_block_l0); 

    if (currSlice->apply_weights)
    {
      if (curr_mb_field && ((p_Img->active_pps->weighted_pred_flag&&(p_Img->type==P_SLICE|| p_Img->type == SP_SLICE))||
         (p_Img->active_pps->weighted_bipred_idc==1 && (p_Img->type==B_SLICE))))
      {
        ref_idx_wp >>=1;
      }

      alpha_l0  = currSlice->wp_weight[pred_dir][ref_idx_wp][0];
      wp_offset = currSlice->wp_offset[pred_dir][ref_idx_wp][0];

      weighted_mc_prediction(&currSlice->mb_pred[pl][joff], block_size_y, block_size_x, ioff, currSlice->tmp_block_l0, alpha_l0, wp_offset, currSlice->luma_log2_weight_denom, max_imgpel_value);
    }
    else
    {
      mc_prediction(&currSlice->mb_pred[pl][joff], block_size_y, block_size_x, ioff, currSlice->tmp_block_l0); 
    }

    if ((dec_picture->chroma_format_idc != YUV400) && (dec_picture->chroma_format_idc != YUV444) ) 
    {
      int uv;

      int ioff_cr = (p_Img->mb_cr_size_x == MB_BLOCK_SIZE) ? ioff : ioff >> 1;
      int joff_cr = (p_Img->mb_cr_size_y == MB_BLOCK_SIZE) ? joff : joff >> 1;
      int block_size_x_cr = p_Img->mb_cr_size_x == MB_BLOCK_SIZE ? block_size_x : block_size_x >> 1;
      int block_size_y_cr = p_Img->mb_cr_size_y == MB_BLOCK_SIZE ? block_size_y : block_size_y >> 1;

      vec1_y_cr = vec1_y + ((active_sps->chroma_format_idc == 1)? list->chroma_vector_adjustment : 0);

      for(uv=0;uv<2;uv++)
      {
        get_block_chroma (currMB, uv, list, vec1_x, vec1_y_cr, block_size_x_cr, block_size_y_cr, currSlice->tmp_block_l0);

        if (currSlice->apply_weights)
        {
          alpha_l0  = currSlice->wp_weight[pred_dir][ref_idx_wp][uv + 1];
          wp_offset = currSlice->wp_offset[pred_dir][ref_idx_wp][uv + 1];

          weighted_mc_prediction(&currSlice->mb_pred[uv + 1][joff_cr], block_size_y_cr, block_size_x_cr, ioff_cr, currSlice->tmp_block_l0, alpha_l0, wp_offset, currSlice->chroma_log2_weight_denom, p_Img->max_imgpel_value_comp[uv + 1]);
        }
        else
        {
          mc_prediction(&currSlice->mb_pred[uv + 1][joff_cr], block_size_y_cr, block_size_x_cr, ioff_cr, currSlice->tmp_block_l0);
        }
      }
    }
  }
  else
  {
    //===== BI-PREDICTION =====
    short *l0_mv_array = dec_picture->motion.mv[LIST_0][j4][i4];
    short *l1_mv_array = dec_picture->motion.mv[LIST_1][j4][i4];

    short l0_refframe = dec_picture->motion.ref_idx[LIST_0][j4][i4];
    short l0_ref_idx  = l0_refframe;
    short l1_refframe = dec_picture->motion.ref_idx[LIST_1][j4][i4];
    short l1_ref_idx  = l1_refframe;

    check_motion_vector_range(p_Img, l0_mv_array[0], l0_mv_array[1]);
    check_motion_vector_range(p_Img, l1_mv_array[0], l1_mv_array[1]);
    vec1_x = i4 * mv_mul + l0_mv_array[0];
    vec2_x = i4 * mv_mul + l1_mv_array[0];

    vec1_y = (currMB->block_y_aff + j) * mv_mul + l0_mv_array[1];
    vec2_y = (currMB->block_y_aff + j) * mv_mul + l1_mv_array[1];

    get_block_luma(currMB, pl, p_Img->listX[LIST_0 + list_offset][l0_refframe], vec1_x, vec1_y, block_size_x, block_size_y, currSlice->tmp_block_l0);  
    get_block_luma(currMB, pl, p_Img->listX[LIST_1 + list_offset][l1_refframe], vec2_x, vec2_y, block_size_x, block_size_y, currSlice->tmp_block_l1);  

    if(currSlice->apply_weights)
    {
      int wt_list_offset = (p_Img->active_pps->weighted_bipred_idc==2)? list_offset : 0;

      // This code existed in the original. Seems pointless but copying it here for reference and in case temporal direct breaks.
      // if (mv_mode==0 && currSlice->direct_spatial_mv_pred_flag==0 ) l1_ref_idx=0;    
      if (((p_Img->active_pps->weighted_pred_flag&&(p_Img->type==P_SLICE|| p_Img->type == SP_SLICE))||
        (p_Img->active_pps->weighted_bipred_idc==1 && (p_Img->type==B_SLICE))) && curr_mb_field)
      {
        l0_ref_idx >>=1;
        l1_ref_idx >>=1;
      }

      alpha_l0  =   currSlice->wbp_weight[LIST_0 + wt_list_offset][l0_ref_idx][l1_ref_idx][0];
      alpha_l1  =   currSlice->wbp_weight[LIST_1 + wt_list_offset][l0_ref_idx][l1_ref_idx][0];
      wp_offset = ((currSlice->wp_offset [LIST_0 + wt_list_offset][l0_ref_idx][0] + currSlice->wp_offset[LIST_1 + wt_list_offset][l1_ref_idx][0] + 1) >>1);

      weighted_bi_prediction(&currSlice->mb_pred[pl][joff], currSlice->tmp_block_l0, currSlice->tmp_block_l1, block_size_y, block_size_x, ioff, alpha_l0, alpha_l1, wp_offset, (currSlice->luma_log2_weight_denom + 1), max_imgpel_value);
    }
    else
    { 
      bi_prediction(&currSlice->mb_pred[pl][joff], currSlice->tmp_block_l0, currSlice->tmp_block_l1, block_size_y, block_size_x, ioff); 
    }

    if ((dec_picture->chroma_format_idc != YUV400) && (dec_picture->chroma_format_idc != YUV444) ) 
    {
      int uv;

      int ioff_cr = p_Img->mb_cr_size_x == MB_BLOCK_SIZE ? ioff : ioff >> 1;
      int joff_cr = p_Img->mb_cr_size_y == MB_BLOCK_SIZE ? joff : joff >> 1;
      int block_size_x_cr = p_Img->mb_cr_size_x == MB_BLOCK_SIZE ? block_size_x : block_size_x >> 1;
      int block_size_y_cr = p_Img->mb_cr_size_y == MB_BLOCK_SIZE ? block_size_y : block_size_y >> 1;

      vec1_y_cr = vec1_y + ((active_sps->chroma_format_idc == 1)? p_Img->listX[LIST_0 + list_offset][l0_refframe]->chroma_vector_adjustment : 0);
      vec2_y_cr = vec2_y + ((active_sps->chroma_format_idc == 1)? p_Img->listX[LIST_1 + list_offset][l1_refframe]->chroma_vector_adjustment : 0);

      for(uv=0;uv<2;uv++)
      {
        get_block_chroma (currMB, uv, p_Img->listX[LIST_0 + list_offset][l0_refframe], vec1_x, vec1_y_cr, block_size_x_cr, block_size_y_cr, currSlice->tmp_block_l0);
        get_block_chroma (currMB, uv, p_Img->listX[LIST_1 + list_offset][l1_refframe], vec2_x, vec2_y_cr, block_size_x_cr, block_size_y_cr, currSlice->tmp_block_l1);

        if(currSlice->apply_weights)
        {
          int wt_list_offset = (p_Img->active_pps->weighted_bipred_idc==2)? list_offset : 0;

          alpha_l0  =   currSlice->wbp_weight[LIST_0 + wt_list_offset][l0_ref_idx][l1_ref_idx][uv + 1];
          alpha_l1  =   currSlice->wbp_weight[LIST_1 + wt_list_offset][l0_ref_idx][l1_ref_idx][uv + 1];
          wp_offset = ((currSlice->wp_offset [LIST_0 + wt_list_offset][l0_ref_idx][uv + 1] + currSlice->wp_offset[LIST_1 + wt_list_offset][l1_ref_idx][uv + 1] + 1) >>1);

          weighted_bi_prediction(&currSlice->mb_pred[uv+1][joff_cr], currSlice->tmp_block_l0, currSlice->tmp_block_l1, block_size_y_cr, block_size_x_cr, ioff_cr, alpha_l0, alpha_l1, wp_offset, (currSlice->chroma_log2_weight_denom + 1), p_Img->max_imgpel_value_comp[uv + 1]);
        }
        else
        {
          bi_prediction(&currSlice->mb_pred[uv + 1][joff_cr], currSlice->tmp_block_l0, currSlice->tmp_block_l1, block_size_y_cr, block_size_x_cr, ioff_cr);
        }
      }
    }      
  }
}

