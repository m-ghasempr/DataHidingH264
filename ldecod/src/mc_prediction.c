
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
#include <assert.h>
#include <string.h>

#include "global.h"
#include "block.h"
#include "mc_prediction.h"
#include "mbuffer.h"
#include "mb_access.h"

extern StorablePicture *no_reference_picture;
extern const unsigned char subblk_offset_y[3][8][4];
extern const unsigned char subblk_offset_x[3][8][4];
extern const unsigned char cofuv_blk_y[3][8][4];
extern const unsigned char cofuv_blk_x[3][8][4];
/*!
 ************************************************************************
 * \brief
 *    block single list prediction
 ************************************************************************
 */
static void mc_prediction(struct img_par *img, 
                    int yuv,
                    int ver_block_size, 
                    int hor_block_size,
                    int joff,
                    int ioff,
                    imgpel block[MB_BLOCK_SIZE][MB_BLOCK_SIZE])
{
  static int jj;
  imgpel (*mpr) [16] = img->mpr[yuv];

  for(jj = 0; jj < ver_block_size; jj++)
  {
    memcpy(&(mpr[jj + joff][ioff]), &(block[jj][0]), hor_block_size * sizeof(imgpel));
  }
}

/*!
 ************************************************************************
 * \brief
 *    block single list weighted prediction
 ************************************************************************
 */
static void weighted_mc_prediction(struct img_par *img,
                            int yuv, 
                            int ver_block_size, 
                            int hor_block_size,
                            int joff,
                            int ioff,
                            imgpel block[MB_BLOCK_SIZE][MB_BLOCK_SIZE], 
                            int wp_scale,
                            int wp_offset,
                            int weight_denom,
                            int color_clip)
{
  static int ii, jj;
  imgpel *mpr;
  
  for(jj=0;jj<ver_block_size;jj++)
  {
    mpr = img->mpr[yuv][jj + joff];
    for(ii=0;ii<hor_block_size;ii++)
      mpr[ii+ioff] = iClip1(color_clip, 
      (rshift_rnd((wp_scale *  block[jj][ii]), weight_denom)  + wp_offset ));
  }
}


/*!
 ************************************************************************
 * \brief
 *    block biprediction
 ************************************************************************
 */
static void bi_prediction(struct img_par *img, 
                    int yuv,
                    int ver_block_size, 
                    int hor_block_size,
                    int joff,
                    int ioff,
                    imgpel block_l0[MB_BLOCK_SIZE][MB_BLOCK_SIZE], 
                    imgpel block_l1[MB_BLOCK_SIZE][MB_BLOCK_SIZE])
{
  static int ii, jj;
  imgpel *mpr;

  for(jj = 0;jj < ver_block_size;jj++)
  {
    mpr = img->mpr[yuv][jj + joff];    
    for(ii = 0; ii < hor_block_size;ii++)
      mpr[ii + ioff] = rshift_rnd_sf(block_l0[jj][ii] + block_l1[jj][ii], 1);
  }
}

/*!
 ************************************************************************
 * \brief
 *    block weighted biprediction
 ************************************************************************
 */
static void weighted_bi_prediction(struct img_par *img, 
                            int yuv,
                            int ver_block_size, 
                            int hor_block_size,
                            int joff,
                            int ioff,
                            imgpel block_l0[MB_BLOCK_SIZE][MB_BLOCK_SIZE], 
                            imgpel block_l1[MB_BLOCK_SIZE][MB_BLOCK_SIZE],
                            int wp_scale_l0,
                            int wp_scale_l1,
                            int wp_offset,
                            int weight_denom,
                            int color_clip)
{
  static int ii, jj;
  imgpel *mpr;
  
  for(jj=0;jj<ver_block_size;jj++)
  {
    mpr = img->mpr[yuv][jj + joff];    
    for(ii=0;ii<hor_block_size;ii++)
      mpr[ii+ioff] = (int)iClip1(color_clip, 
      (rshift_rnd((wp_scale_l0 * block_l0[jj][ii] + wp_scale_l1 * block_l1[jj][ii]), 
      weight_denom) + wp_offset));
  }
}

/*!
 ************************************************************************
 * \brief
 *    Interpolation of 1/4 subpixel
 ************************************************************************
 */
void get_block_luma(int ref_frame, StorablePicture **list, int x_pos, int y_pos, int hor_block_size, int ver_block_size, struct img_par *img, imgpel block[MB_BLOCK_SIZE][MB_BLOCK_SIZE])
{
  int dx = (x_pos & 3), dy = (y_pos & 3);
  int i, j;
  int shift_x = dec_picture->size_x;
  int maxold_x = dec_picture->size_x_m1;
  int maxold_y = (dec_picture->mb_field[img->current_mb_nr]) ? (dec_picture->size_y >> 1) - 1 : dec_picture->size_y_m1;
  int result;
  int pres_x, pres_y;

  static int tmp_res[21][21];
  static const int COEF[6] = { 1, -5, 20, 20, -5, 1  };
  StorablePicture *curr_ref = list[ref_frame];
  static imgpel **cur_imgY, *cur_lineY;
  static imgpel *cur_lineY_m2, *cur_lineY_m1, *cur_lineY_p1, *cur_lineY_p2, *cur_lineY_p3;
  int tmp_pos;
  static int jpos_m2, jpos_m1, jpos, jpos_p1, jpos_p2, jpos_p3;
  static int ipos_m2, ipos_m1, ipos, ipos_p1, ipos_p2, ipos_p3;

  if (curr_ref == no_reference_picture && img->framepoc < img->recovery_poc)
  {
    printf("list[ref_frame] is equal to 'no reference picture' before RAP\n");

    /* fill the block with sample value 128 */
    for (j = 0; j < ver_block_size; j++)
      for (i = 0; i < hor_block_size; i++)
        block[j][i] = 128;
    return;
  }

  if( IS_INDEPENDENT(img) )
  {
    switch( img->colour_plane_id )
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
  else
  {
    cur_imgY = curr_ref->imgY;
  }

  x_pos = x_pos >> 2;
  y_pos = y_pos >> 2;

  if ( (y_pos > 1) && (y_pos < maxold_y - 2 - ver_block_size) && (x_pos > 1) && (x_pos < maxold_x - 2 - hor_block_size))
  {
    if (dx == 0 && dy == 0)
    {  /* fullpel position */
      for (j = 0; j < ver_block_size; j++)
      {        
        memcpy(&(block[j][0]), &(cur_imgY[ y_pos + j ][x_pos]), hor_block_size * sizeof(imgpel));
      }
    }
    else
    { /* other positions */

      if (dy == 0)
      { /* No vertical interpolation */
        tmp_pos = x_pos - 2;
        for (i = 0; i < hor_block_size; i++)
        {        
          ipos_m2 = tmp_pos++;
          ipos_m1 = tmp_pos++;
          ipos    = tmp_pos++;
          ipos_p1 = tmp_pos++;
          ipos_p2 = tmp_pos++;
          ipos_p3 = tmp_pos  ;
          tmp_pos -= 4;

          for (j = 0; j < ver_block_size; j++)
          {
            cur_lineY = cur_imgY[y_pos + j];

            result  = (cur_lineY[ipos_m2] + cur_lineY[ipos_p3]) * COEF[0];
            result += (cur_lineY[ipos_m1] + cur_lineY[ipos_p2]) * COEF[1];
            result += (cur_lineY[ipos   ] + cur_lineY[ipos_p1]) * COEF[2];

            block[j][i] = iClip1(img->max_imgpel_value, ((result + 16)>>5));
          }
        }

        if ((dx&1) == 1)
        {
          for (j = 0; j < ver_block_size; j++)
          {
            cur_lineY = &(cur_imgY[ y_pos + j ][x_pos + (dx >> 1)]);
            for (i = 0; i < hor_block_size; i++)
            {
              block[j][i] = (block[j][i] + *(cur_lineY++) + 1 )>>1;
            }
          }
        }
      }
      else if (dx == 0)
      {  /* No horizontal interpolation */        
        cur_lineY_m2 = &(cur_imgY[y_pos - 2][x_pos]);
        for (j = 0; j < ver_block_size; j++)
        {                  
          cur_lineY_m1 = cur_lineY_m2 + shift_x;          
          cur_lineY    = cur_lineY_m1 + shift_x;
          cur_lineY_p1 = cur_lineY    + shift_x;
          cur_lineY_p2 = cur_lineY_p1 + shift_x;
          cur_lineY_p3 = cur_lineY_p2 + shift_x;

          for (i = 0; i < hor_block_size; i++)
          {
            result  = (*(cur_lineY_m2++) + *(cur_lineY_p3++)) * COEF[0];
            result += (*(cur_lineY_m1++) + *(cur_lineY_p2++)) * COEF[1];
            result += (*(cur_lineY   ++) + *(cur_lineY_p1++)) * COEF[2];

            block[j][i] = iClip1(img->max_imgpel_value, ((result+16)>>5));
          }
          cur_lineY_m2 = cur_lineY_m1 - hor_block_size;
        }

        if ((dy&1) == 1)
        {
          for (j = 0; j < ver_block_size; j++)
          {
            cur_lineY = &(cur_imgY[ y_pos + j + (dy >> 1)][x_pos]);
            for (i = 0; i < hor_block_size; i++)
            {
              block[j][i] = (block[j][i] + *(cur_lineY++) + 1 )>>1;
            }
          }
        }
      }
      else if (dx == 2)
      {  /* Vertical & horizontal interpolation */
        tmp_pos = x_pos - 2;
        for (i = 0; i < hor_block_size; i++)
        {        
          ipos_m2 = tmp_pos++;
          ipos_m1 = tmp_pos++;
          ipos    = tmp_pos++;
          ipos_p1 = tmp_pos++;
          ipos_p2 = tmp_pos++;
          ipos_p3 = tmp_pos  ;
          tmp_pos -= 4;

          for (j = 0; j < ver_block_size + 5; j++)
          {
            cur_lineY = cur_imgY[y_pos + j - 2];

            tmp_res[j][i]  = (cur_lineY[ipos_m2] + cur_lineY[ipos_p3]) * COEF[0];
            tmp_res[j][i] += (cur_lineY[ipos_m1] + cur_lineY[ipos_p2]) * COEF[1];
            tmp_res[j][i] += (cur_lineY[ipos   ] + cur_lineY[ipos_p1]) * COEF[2];
          }
        }

        for (j = 0; j < ver_block_size; j++)
        {
          jpos_m2 = j    ;
          jpos_m1 = j + 1;
          jpos    = j + 2;
          jpos_p1 = j + 3;
          jpos_p2 = j + 4;
          jpos_p3 = j + 5;

          for (i = 0; i < hor_block_size; i++)
          {
            result  = (tmp_res[jpos_m2][i] + tmp_res[jpos_p3][i]) * COEF[0];
            result += (tmp_res[jpos_m1][i] + tmp_res[jpos_p2][i]) * COEF[1];
            result += (tmp_res[jpos   ][i] + tmp_res[jpos_p1][i]) * COEF[2];

            block[j][i] = iClip1(img->max_imgpel_value, ((result+512)>>10));
          }
        }

        if ((dy&1) == 1)
        {
          for (j = 0; j < ver_block_size; j++)
          {
            pres_y = j+2+(dy>>1);
            for (i = 0; i < hor_block_size; i++)
            {
              block[j][i] = (block[j][i] + iClip1(img->max_imgpel_value, ((tmp_res[pres_y][i] + 16) >> 5)) + 1 )>>1;
            }
          }
        }
      }
      else if (dy == 2)
      {  /* Horizontal & vertical interpolation */
        cur_lineY_m2 = &(cur_imgY[y_pos - 2][x_pos - 2]);
        for (j = 0; j < ver_block_size; j++)
        {                    
          cur_lineY_m1 = cur_lineY_m2 + shift_x;
          cur_lineY    = cur_lineY_m1 + shift_x;
          cur_lineY_p1 = cur_lineY    + shift_x;
          cur_lineY_p2 = cur_lineY_p1 + shift_x;
          cur_lineY_p3 = cur_lineY_p2 + shift_x;

          for (i = 0; i < hor_block_size + 5; i++)
          {
            tmp_res[j][i]  = (*(cur_lineY_m2++) + *(cur_lineY_p3++)) * COEF[0];
            tmp_res[j][i] += (*(cur_lineY_m1++) + *(cur_lineY_p2++)) * COEF[1];
            tmp_res[j][i] += (*(cur_lineY   ++) + *(cur_lineY_p1++)) * COEF[2];
          }
          cur_lineY_m2 = cur_lineY_m1 - (hor_block_size + 5);
        }

        for (j = 0; j < ver_block_size; j++)
        {
          for (i = 0; i < hor_block_size; i++)
          {
            result  = (tmp_res[j][i    ] + tmp_res[j][i + 5]) * COEF[0];
            result += (tmp_res[j][i + 1] + tmp_res[j][i + 4]) * COEF[1];
            result += (tmp_res[j][i + 2] + tmp_res[j][i + 3]) * COEF[2];
            block[j][i] = iClip1(img->max_imgpel_value, ((result+512)>>10));
          }
        }

        if ((dx&1) == 1)
        {
          for (j = 0; j < ver_block_size; j++)
          {
            for (i = 0; i < hor_block_size; i++)
            {
              block[j][i] = (block[j][i] + iClip1(img->max_imgpel_value, ((tmp_res[j][i + 2 + (dx>>1)] + 16)>>5))+1)>>1;
            }
          }
        }
      }
      else
      {  /* Diagonal interpolation */
        tmp_pos = x_pos - 2;
        for (i = 0; i < hor_block_size; i++)
        {
          ipos_m2 = tmp_pos++;
          ipos_m1 = tmp_pos++;
          ipos    = tmp_pos++;
          ipos_p1 = tmp_pos++;
          ipos_p2 = tmp_pos++;
          ipos_p3 = tmp_pos  ;
          tmp_pos -= 4;

          for (j = 0; j < ver_block_size; j++)
          {
            cur_lineY = cur_imgY[(dy == 1 ? y_pos + j : y_pos + j + 1)];

            result  = (cur_lineY[ipos_m2] + cur_lineY[ipos_p3]) * COEF[0];
            result += (cur_lineY[ipos_m1] + cur_lineY[ipos_p2]) * COEF[1];
            result += (cur_lineY[ipos   ] + cur_lineY[ipos_p1]) * COEF[2];

            block[j][i] = iClip1(img->max_imgpel_value, ((result+16)>>5));
          }
        }

        cur_lineY_m2 = &(cur_imgY[y_pos - 2][(dx == 1 ? x_pos : x_pos + 1)]);
        for (j = 0; j < ver_block_size; j++)
        {        
          cur_lineY_m1 = cur_lineY_m2 + shift_x;
          cur_lineY    = cur_lineY_m1 + shift_x;
          cur_lineY_p1 = cur_lineY    + shift_x;
          cur_lineY_p2 = cur_lineY_p1 + shift_x;
          cur_lineY_p3 = cur_lineY_p2 + shift_x;

          for (i = 0; i < hor_block_size; i++)
          {
            result  = (*(cur_lineY_m2++) + *(cur_lineY_p3++)) * COEF[0];
            result += (*(cur_lineY_m1++) + *(cur_lineY_p2++)) * COEF[1];
            result += (*(cur_lineY   ++) + *(cur_lineY_p1++)) * COEF[2];

            block[j][i] = (block[j][i] + iClip1(img->max_imgpel_value, ((result + 16) >> 5)) +1 ) >>1;
          }
          cur_lineY_m2 = cur_lineY_m1 - hor_block_size ;
        }      
      }
    }
  }
  else // unsafe positions
  {
    if (dx == 0 && dy == 0)
    {  /* fullpel position */
      for (j = 0; j < ver_block_size; j++)
      {
        cur_lineY = cur_imgY[iClip3(0, maxold_y, y_pos + j)];
        for (i = 0; i < hor_block_size; i++)
        {
          block[j][i] = cur_lineY[iClip3(0, maxold_x, x_pos + i )];
        }
      }
    }
    else
    { /* other positions */

      if (dy == 0)
      { /* No vertical interpolation */
        tmp_pos = x_pos - 2;
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

            block[j][i] = iClip1(img->max_imgpel_value, ((result + 16)>>5));
          }
        }

        if ((dx&1) == 1)
        {
          for (j = 0; j < ver_block_size; j++)
          {
            cur_lineY = cur_imgY[iClip3(0,maxold_y,y_pos+j)];
            for (i = 0; i < hor_block_size; i++)
            {
              block[j][i] = (block[j][i] + cur_lineY[iClip3(0, maxold_x, x_pos + i + (dx >> 1))] + 1 )>>1;
            }
          }
        }
      }
      else if (dx == 0)
      {  /* No horizontal interpolation */
        tmp_pos = y_pos - 2;
        for (j = 0; j < ver_block_size; j++)
        {        
          jpos_m2 = iClip3(0, maxold_y, tmp_pos++);
          jpos_m1 = iClip3(0, maxold_y, tmp_pos++);
          jpos    = iClip3(0, maxold_y, tmp_pos++);
          jpos_p1 = iClip3(0, maxold_y, tmp_pos++);
          jpos_p2 = iClip3(0, maxold_y, tmp_pos++);
          jpos_p3 = iClip3(0, maxold_y, tmp_pos  );
          tmp_pos -= 4;

          for (i = 0; i < hor_block_size; i++)
          {
            pres_x = iClip3(0,maxold_x,x_pos+i);

            result  = (cur_imgY[jpos_m2][pres_x] + cur_imgY[jpos_p3][pres_x]) * COEF[0];
            result += (cur_imgY[jpos_m1][pres_x] + cur_imgY[jpos_p2][pres_x]) * COEF[1];
            result += (cur_imgY[jpos   ][pres_x] + cur_imgY[jpos_p1][pres_x]) * COEF[2];
            block[j][i] = iClip1(img->max_imgpel_value, ((result+16)>>5));
          }
        }

        if ((dy&1) == 1)
        {
          for (j = 0; j < ver_block_size; j++)
          {
            cur_lineY = cur_imgY[iClip3(0,maxold_y,y_pos+j+(dy>>1))];
            for (i = 0; i < hor_block_size; i++)
            {
              block[j][i] = (block[j][i] + cur_lineY[iClip3(0, maxold_x, x_pos + i)] + 1 )>>1;
            }
          }
        }
      }
      else if (dx == 2)
      {  /* Vertical & horizontal interpolation */
        tmp_pos = x_pos - 2;
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
          jpos_m2 = j    ;
          jpos_m1 = j + 1;
          jpos    = j + 2;
          jpos_p1 = j + 3;
          jpos_p2 = j + 4;
          jpos_p3 = j + 5;

          for (i = 0; i < hor_block_size; i++)
          {
            result  = (tmp_res[jpos_m2][i] + tmp_res[jpos_p3][i]) * COEF[0];
            result += (tmp_res[jpos_m1][i] + tmp_res[jpos_p2][i]) * COEF[1];
            result += (tmp_res[jpos   ][i] + tmp_res[jpos_p1][i]) * COEF[2];

            block[j][i] = iClip1(img->max_imgpel_value, ((result+512)>>10));
          }
        }

        if ((dy&1) == 1)
        {
          for (j = 0; j < ver_block_size; j++)
          {
            pres_y = j+2+(dy>>1);
            for (i = 0; i < hor_block_size; i++)
            {
              block[j][i] = (block[j][i] + iClip1(img->max_imgpel_value, ((tmp_res[pres_y][i] + 16) >> 5)) + 1 )>>1;
            }
          }
        }
      }
      else if (dy == 2)
      {  /* Horizontal & vertical interpolation */

        tmp_pos = y_pos - 2;
        for (j = 0; j < ver_block_size; j++)
        {

          jpos_m2 = iClip3(0, maxold_y, tmp_pos++);
          jpos_m1 = iClip3(0, maxold_y, tmp_pos++);
          jpos    = iClip3(0, maxold_y, tmp_pos++);
          jpos_p1 = iClip3(0, maxold_y, tmp_pos++);
          jpos_p2 = iClip3(0, maxold_y, tmp_pos++);
          jpos_p3 = iClip3(0, maxold_y, tmp_pos  );
          tmp_pos -= 4;

          for (i = 0; i < hor_block_size + 5; i++)
          {
            pres_x = iClip3(0,maxold_x, x_pos + i - 2);
            tmp_res[j][i]  = (cur_imgY[jpos_m2][pres_x] + cur_imgY[jpos_p3][pres_x])*COEF[0];
            tmp_res[j][i] += (cur_imgY[jpos_m1][pres_x] + cur_imgY[jpos_p2][pres_x])*COEF[1];
            tmp_res[j][i] += (cur_imgY[jpos   ][pres_x] + cur_imgY[jpos_p1][pres_x])*COEF[2];
          }
        }

        for (j = 0; j < ver_block_size; j++)
        {
          for (i = 0; i < hor_block_size; i++)
          {
            result  = (tmp_res[j][i    ] + tmp_res[j][i + 5]) * COEF[0];
            result += (tmp_res[j][i + 1] + tmp_res[j][i + 4]) * COEF[1];
            result += (tmp_res[j][i + 2] + tmp_res[j][i + 3]) * COEF[2];
            block[j][i] = iClip1(img->max_imgpel_value, ((result+512)>>10));
          }
        }

        if ((dx&1) == 1)
        {
          for (j = 0; j < ver_block_size; j++)
          {
            for (i = 0; i < hor_block_size; i++)
            {
              block[j][i] = (block[j][i] + iClip1(img->max_imgpel_value, ((tmp_res[j][i + 2 + (dx>>1)] + 16)>>5))+1)>>1;
            }
          }
        }
      }
      else
      {  /* Diagonal interpolation */
        tmp_pos = x_pos - 2;
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

            block[j][i] = iClip1(img->max_imgpel_value, ((result+16)>>5));
          }
        }

        tmp_pos = y_pos - 2;
        for (j = 0; j < ver_block_size; j++)
        {        
          jpos_m2 = iClip3(0, maxold_y, tmp_pos++);
          jpos_m1 = iClip3(0, maxold_y, tmp_pos++);
          jpos    = iClip3(0, maxold_y, tmp_pos++);
          jpos_p1 = iClip3(0, maxold_y, tmp_pos++);
          jpos_p2 = iClip3(0, maxold_y, tmp_pos++);
          jpos_p3 = iClip3(0, maxold_y, tmp_pos  );
          tmp_pos -= 4;

          for (i = 0; i < hor_block_size; i++)
          {
            pres_x = dx == 1 ? x_pos+i : x_pos+i+1;
            pres_x = iClip3(0, maxold_x, pres_x);

            result  = (cur_imgY[jpos_m2][pres_x] + cur_imgY[jpos_p3][pres_x]) * COEF[0];
            result += (cur_imgY[jpos_m1][pres_x] + cur_imgY[jpos_p2][pres_x]) * COEF[1];
            result += (cur_imgY[jpos   ][pres_x] + cur_imgY[jpos_p1][pres_x]) * COEF[2];

            block[j][i] = (block[j][i] + iClip1(img->max_imgpel_value, ((result+16)>>5)) +1 ) >>1;
          }
        }      
      }
    }
  }
}

void get_block_chroma(int uv, int ref_frame, StorablePicture **list, int x_pos, int y_pos, int hor_block_size, int ver_block_size, struct img_par *img, imgpel block[MB_BLOCK_SIZE][MB_BLOCK_SIZE])
{
  int subpel_x = img->mb_cr_size_x == 8 ? 7 : 3;
  int subpel_y = img->mb_cr_size_y == 8 ? 7 : 3;
  int shiftpel_x  = img->mb_cr_size_x == 8 ? 3 : 2;
  int shiftpel_y  = img->mb_cr_size_y == 8 ? 3 : 2;
  int total_scale = shiftpel_x + shiftpel_y;

  int dx = (x_pos & subpel_x);
  int dy = (y_pos & subpel_y);
  int dxcur = (subpel_x + 1 - dx);
  int dycur = (subpel_y + 1 - dy);

  int i, j;
  int maxold_x = dec_picture->size_x_cr_m1;
  int maxold_y = (dec_picture->mb_field[img->current_mb_nr]) ? (dec_picture->size_y_cr >> 1) - 1 : dec_picture->size_y_cr_m1;
  int result;

  StorablePicture *curr_ref = list[ref_frame];
  static imgpel **cur_img;
  static imgpel *cur_line, *cur_line_p1;
  int tmp_pos;
  static int ipos, ipos_p1;

  if (curr_ref == no_reference_picture && img->framepoc < img->recovery_poc)
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
      for (j = 0; j < ver_block_size; j++)
      {
        cur_line    = cur_img[y_pos + j];
        cur_line_p1 = cur_img[y_pos + j + 1];
        tmp_pos = x_pos;

        for (i = 0; i < hor_block_size; i++)
        {
          ipos    = tmp_pos++;

          result = (
            dxcur * dycur * cur_line   [ipos   ] + 
            dxcur * dy    * cur_line_p1[ipos   ]);
          block[j][i] = iClip1(img->max_imgpel_value_uv, rshift_rnd_sf(result, total_scale));
        }
      }
    }
    else if (dy == 0)
    { 
      for (j = 0; j < ver_block_size; j++)
      {
        cur_line    = cur_img[y_pos + j];
        tmp_pos = x_pos;

        for (i = 0; i < hor_block_size; i++)
        {
          ipos    = tmp_pos++;
          ipos_p1 = tmp_pos;

          result = (
            dxcur * dycur * cur_line   [ipos   ] + 
            dx    * dycur * cur_line   [ipos_p1]);
          block[j][i] = iClip1(img->max_imgpel_value_uv, rshift_rnd_sf(result, total_scale));
        }
      }
    }
    else
    { /* other positions */
      for (j = 0; j < ver_block_size; j++)
      {
        cur_line    = cur_img[y_pos + j];
        cur_line_p1 = cur_img[y_pos + j + 1];
        tmp_pos = x_pos;

        for (i = 0; i < hor_block_size; i++)
        {
          ipos    = tmp_pos++;
          ipos_p1 = tmp_pos;

          result = (
            dxcur * dycur * cur_line   [ipos   ] + 
            dx    * dycur * cur_line   [ipos_p1] +
            dxcur * dy    * cur_line_p1[ipos   ] +
            dx    * dy    * cur_line_p1[ipos_p1]);
          block[j][i] = iClip1(img->max_imgpel_value_uv, rshift_rnd_sf(result, total_scale));
        }
      }
    }
  }
  else // unsafe positions
  {
    if (dx == 0 && dy == 0)
    {  /* fullpel position */
      for (j = 0; j < ver_block_size; j++)
      {
        cur_line = cur_img[iClip3(0, maxold_y, y_pos + j)];
        for (i = 0; i < hor_block_size; i++)
        {
          block[j][i] = cur_line[iClip3(0, maxold_x, x_pos + i )];
        }
      }
    }
    else if (dx == 0)
    { 
      for (j = 0; j < ver_block_size; j++)
      {
        cur_line    = cur_img[iClip3(0, maxold_y, y_pos + j)];
        cur_line_p1 = cur_img[iClip3(0, maxold_y, y_pos + j + 1)];
        tmp_pos = x_pos;

        for (i = 0; i < hor_block_size; i++)
        {
          ipos    = iClip3(0, maxold_x, tmp_pos++);

          result = (
            dxcur * dycur * cur_line   [ipos   ] + 
            dxcur * dy    * cur_line_p1[ipos   ]);
          block[j][i] = iClip1(img->max_imgpel_value_uv, rshift_rnd_sf(result, total_scale));
        }
      }      
    }
    else if (dy == 0)
    { 
      for (j = 0; j < ver_block_size; j++)
      {
        cur_line    = cur_img[iClip3(0, maxold_y, y_pos + j)];
        tmp_pos = x_pos;

        for (i = 0; i < hor_block_size; i++)
        {
          ipos    = iClip3(0, maxold_x, tmp_pos++);
          ipos_p1 = iClip3(0, maxold_x, tmp_pos  );


          result = (
            dxcur * dycur * cur_line   [ipos   ] + 
            dx    * dycur * cur_line   [ipos_p1]);
          block[j][i] = iClip1(img->max_imgpel_value_uv, rshift_rnd_sf(result, total_scale));
        }
      }      
    }
    else
    { /* other positions */ 
      for (j = 0; j < ver_block_size; j++)
      {
        cur_line    = cur_img[iClip3(0, maxold_y, y_pos + j)];
        cur_line_p1 = cur_img[iClip3(0, maxold_y, y_pos + j + 1)];
        tmp_pos = x_pos;

        for (i = 0; i < hor_block_size; i++)
        {
          ipos    = iClip3(0, maxold_x, tmp_pos++);
          ipos_p1 = iClip3(0, maxold_x, tmp_pos  );


          result = (
            dxcur * dycur * cur_line   [ipos   ] + 
            dx    * dycur * cur_line   [ipos_p1] +
            dxcur * dy    * cur_line_p1[ipos   ] +
            dx    * dy    * cur_line_p1[ipos_p1]);
          block[j][i] = iClip1(img->max_imgpel_value_uv, rshift_rnd_sf(result, total_scale));
        }
      }      
    }
  }
}


void intra_cr_decoding(int yuv, struct img_par *img, int smb, Macroblock *currMB)
{
  imgpel **curUV;
  int uv_shift, uv;
  int b8,b4;
  int i4, j4, ioff, joff, ii, jj, i, j;

  for(uv = 0; uv < 2; uv++)
  {
    curUV = dec_picture->imgUV[uv];
    uv_shift = uv * (img->num_uv_blocks);
    intrapred_chroma(img, uv);

    if (!smb && (currMB->cbp >> 4))
    {
      for (b8 = 0; b8 < (img->num_uv_blocks); b8++)
      {
        for(b4 = 0; b4 < 4; b4++)
        {
          joff = subblk_offset_y[yuv][b8][b4];
          j4   = img->pix_c_y + joff;
          ioff = subblk_offset_x[yuv][b8][b4];
          i4   = img->pix_c_x + ioff;

          itrans(img, ioff, joff, cofuv_blk_x[yuv][b8+uv_shift][b4], cofuv_blk_y[yuv][b8+uv_shift][b4], 1, uv + 1);

          for(jj=0;jj<4;jj++)
            for(ii=0;ii<4;ii++)
            {
              curUV[j4+jj][i4+ii]=img->m7[jj][ii];
            }
        }
      }
    }
    else if ((currMB->cbp >> 4) == 0)
    {
      for (b8 = 0; b8 < (img->num_uv_blocks); b8++)
      {
        for(b4 = 0; b4 < 4; b4++)
        {
          joff = subblk_offset_y[yuv][b8][b4];
          j4   = img->pix_c_y + joff;
          ioff = subblk_offset_x[yuv][b8][b4];
          i4   = img->pix_c_x + ioff;

          for(jj = 0; jj < 4;jj++)
            memcpy(&(curUV[j4+jj][i4]), &(img->mpr[uv + 1][joff + jj][ioff]), BLOCK_SIZE * sizeof(imgpel));
        }
      }
    }
    else
    {
      itrans_sp_cr(img, uv);

      for (j = 4; j < 6; j++)
      {
        joff = (j - 4) * 4;
        j4   = img->pix_c_y+joff;
        for(i = 0; i < 2;i++)
        {
          ioff = i * 4;
          i4   = img->pix_c_x + ioff;
          itrans(img, ioff, joff, 2 * uv + i,j, 1, uv + 1);

          for(jj = 0; jj < 4; jj++)
            for(ii = 0; ii < 4; ii++)
            {
              curUV[j4+jj][i4+ii]=img->m7[jj][ii];
            }
        }
      }
    }
  }
}

void prepare_direct_params(StorablePicture *dec_picture, struct img_par *img, int mb_nr, short pmvl0[2], short pmvl1[2],char *l0_rFrame, char *l1_rFrame)
{
  char l0_rFrameL, l0_rFrameU, l0_rFrameUL, l0_rFrameUR;
  char l1_rFrameL, l1_rFrameU, l1_rFrameUL, l1_rFrameUR;

  PixelPos mb_left, mb_up, mb_upleft, mb_upright;

  getLuma4x4Neighbour(mb_nr, -1,  0, &mb_left);
  getLuma4x4Neighbour(mb_nr,  0, -1, &mb_up);
  getLuma4x4Neighbour(mb_nr, 16, -1, &mb_upright);
  getLuma4x4Neighbour(mb_nr, -1, -1, &mb_upleft);

  if (!img->MbaffFrameFlag)
  {
    l0_rFrameL  = mb_left.available    ? dec_picture->ref_idx[LIST_0][mb_left.pos_y][mb_left.pos_x]       : -1;
    l0_rFrameU  = mb_up.available      ? dec_picture->ref_idx[LIST_0][mb_up.pos_y][mb_up.pos_x]           : -1;
    l0_rFrameUL = mb_upleft.available  ? dec_picture->ref_idx[LIST_0][mb_upleft.pos_y][mb_upleft.pos_x]   : -1;
    l0_rFrameUR = mb_upright.available ? dec_picture->ref_idx[LIST_0][mb_upright.pos_y][mb_upright.pos_x] : l0_rFrameUL;

    l1_rFrameL = mb_left.available     ? dec_picture->ref_idx[LIST_1][mb_left.pos_y][mb_left.pos_x]       : -1;
    l1_rFrameU = mb_up.available       ? dec_picture->ref_idx[LIST_1][mb_up.pos_y][mb_up.pos_x]           : -1;
    l1_rFrameUL = mb_upleft.available  ? dec_picture->ref_idx[LIST_1][mb_upleft.pos_y][mb_upleft.pos_x]   : -1;
    l1_rFrameUR = mb_upright.available ? dec_picture->ref_idx[LIST_1][mb_upright.pos_y][mb_upright.pos_x] : l1_rFrameUL;
  }
  else
  {
    if (img->mb_data[mb_nr].mb_field)
    {
      l0_rFrameL = mb_left.available 
        ? img->mb_data[mb_left.mb_addr].mb_field  || dec_picture->ref_idx[LIST_0][mb_left.pos_y][mb_left.pos_x] < 0
        ? dec_picture->ref_idx[LIST_0][mb_left.pos_y][mb_left.pos_x] 
        : dec_picture->ref_idx[LIST_0][mb_left.pos_y][mb_left.pos_x] * 2: -1;

      l0_rFrameU = mb_up.available 
        ? img->mb_data[mb_up.mb_addr].mb_field || dec_picture->ref_idx[LIST_0][mb_up.pos_y][mb_up.pos_x] < 0
        ? dec_picture->ref_idx[LIST_0][mb_up.pos_y][mb_up.pos_x] 
        : dec_picture->ref_idx[LIST_0][mb_up.pos_y][mb_up.pos_x] * 2: -1;

       l0_rFrameUL = mb_upleft.available 
         ? img->mb_data[mb_upleft.mb_addr].mb_field || dec_picture->ref_idx[LIST_0][mb_upleft.pos_y][mb_upleft.pos_x] < 0
         ? dec_picture->ref_idx[LIST_0][mb_upleft.pos_y][mb_upleft.pos_x] 
         : dec_picture->ref_idx[LIST_0][mb_upleft.pos_y][mb_upleft.pos_x] *2: -1;

       l0_rFrameUR = mb_upright.available 
         ? img->mb_data[mb_upright.mb_addr].mb_field || dec_picture->ref_idx[LIST_0][mb_upright.pos_y][mb_upright.pos_x] < 0 
         ? dec_picture->ref_idx[LIST_0][mb_upright.pos_y][mb_upright.pos_x] 
         : dec_picture->ref_idx[LIST_0][mb_upright.pos_y][mb_upright.pos_x] * 2: l0_rFrameUL;

       l1_rFrameL = mb_left.available 
         ? img->mb_data[mb_left.mb_addr].mb_field || dec_picture->ref_idx[LIST_1][mb_left.pos_y][mb_left.pos_x]  < 0 
         ? dec_picture->ref_idx[LIST_1][mb_left.pos_y][mb_left.pos_x] 
         : dec_picture->ref_idx[LIST_1][mb_left.pos_y][mb_left.pos_x] * 2: -1;

       l1_rFrameU = mb_up.available 
         ? img->mb_data[mb_up.mb_addr].mb_field || dec_picture->ref_idx[LIST_1][mb_up.pos_y][mb_up.pos_x]  < 0 
         ? dec_picture->ref_idx[LIST_1][mb_up.pos_y][mb_up.pos_x] 
         : dec_picture->ref_idx[LIST_1][mb_up.pos_y][mb_up.pos_x] * 2: -1;

       l1_rFrameUL = mb_upleft.available 
         ? img->mb_data[mb_upleft.mb_addr].mb_field || dec_picture->ref_idx[LIST_1][mb_upleft.pos_y][mb_upleft.pos_x]  < 0 
         ? dec_picture->ref_idx[LIST_1][mb_upleft.pos_y][mb_upleft.pos_x] 
         : dec_picture->ref_idx[LIST_1][mb_upleft.pos_y][mb_upleft.pos_x] *2 : -1;

       l1_rFrameUR = mb_upright.available 
         ? img->mb_data[mb_upright.mb_addr].mb_field || dec_picture->ref_idx[LIST_1][mb_upright.pos_y][mb_upright.pos_x] < 0
         ? dec_picture->ref_idx[LIST_1][mb_upright.pos_y][mb_upright.pos_x] 
         : dec_picture->ref_idx[LIST_1][mb_upright.pos_y][mb_upright.pos_x] * 2: l1_rFrameUL;
    }
    else
    {
      l0_rFrameL = mb_left.available 
        ? img->mb_data[mb_left.mb_addr].mb_field || dec_picture->ref_idx[LIST_0][mb_left.pos_y][mb_left.pos_x]  < 0 
        ? dec_picture->ref_idx[LIST_0][mb_left.pos_y][mb_left.pos_x] >> 1 
        : dec_picture->ref_idx[LIST_0][mb_left.pos_y][mb_left.pos_x]: -1;

      l0_rFrameU = mb_up.available 
        ? img->mb_data[mb_up.mb_addr].mb_field || dec_picture->ref_idx[LIST_0][mb_up.pos_y][mb_up.pos_x] < 0 
        ? dec_picture->ref_idx[LIST_0][mb_up.pos_y][mb_up.pos_x] >> 1 
        : dec_picture->ref_idx[LIST_0][mb_up.pos_y][mb_up.pos_x] : -1;

      l0_rFrameUL = mb_upleft.available 
        ? img->mb_data[mb_upleft.mb_addr].mb_field || dec_picture->ref_idx[LIST_0][mb_upleft.pos_y][mb_upleft.pos_x] < 0 
        ? dec_picture->ref_idx[LIST_0][mb_upleft.pos_y][mb_upleft.pos_x]>> 1 
        : dec_picture->ref_idx[LIST_0][mb_upleft.pos_y][mb_upleft.pos_x] : -1;

      l0_rFrameUR = mb_upright.available 
        ? img->mb_data[mb_upright.mb_addr].mb_field || dec_picture->ref_idx[LIST_0][mb_upright.pos_y][mb_upright.pos_x] < 0 
        ? dec_picture->ref_idx[LIST_0][mb_upright.pos_y][mb_upright.pos_x] >> 1 
        : dec_picture->ref_idx[LIST_0][mb_upright.pos_y][mb_upright.pos_x] : l0_rFrameUL;

      l1_rFrameL = mb_left.available 
        ? img->mb_data[mb_left.mb_addr].mb_field || dec_picture->ref_idx[LIST_1][mb_left.pos_y][mb_left.pos_x] < 0 
        ? dec_picture->ref_idx[LIST_1][mb_left.pos_y][mb_left.pos_x] >> 1 
        : dec_picture->ref_idx[LIST_1][mb_left.pos_y][mb_left.pos_x] : -1;

      l1_rFrameU = mb_up.available 
        ? img->mb_data[mb_up.mb_addr].mb_field || dec_picture->ref_idx[LIST_1][mb_up.pos_y][mb_up.pos_x] < 0 
        ? dec_picture->ref_idx[LIST_1][mb_up.pos_y][mb_up.pos_x] >> 1 
        : dec_picture->ref_idx[LIST_1][mb_up.pos_y][mb_up.pos_x] : -1;

      l1_rFrameUL = mb_upleft.available 
        ? img->mb_data[mb_upleft.mb_addr].mb_field || dec_picture->ref_idx[LIST_1][mb_upleft.pos_y][mb_upleft.pos_x] < 0 
        ? dec_picture->ref_idx[LIST_1][mb_upleft.pos_y][mb_upleft.pos_x] >> 1 
        : dec_picture->ref_idx[LIST_1][mb_upleft.pos_y][mb_upleft.pos_x] : -1;

      l1_rFrameUR = mb_upright.available 
        ? img->mb_data[mb_upright.mb_addr].mb_field || dec_picture->ref_idx[LIST_1][mb_upright.pos_y][mb_upright.pos_x] < 0 
        ? dec_picture->ref_idx[LIST_1][mb_upright.pos_y][mb_upright.pos_x] >> 1
        : dec_picture->ref_idx[LIST_1][mb_upright.pos_y][mb_upright.pos_x] : l1_rFrameUL;
    }
  }

  *l0_rFrame = (l0_rFrameL >= 0 && l0_rFrameU >= 0) ? imin(l0_rFrameL,l0_rFrameU): imax(l0_rFrameL,l0_rFrameU);
  *l0_rFrame = (*l0_rFrame >= 0 && l0_rFrameUR >= 0) ? imin(*l0_rFrame,l0_rFrameUR): imax(*l0_rFrame,l0_rFrameUR);

  *l1_rFrame = (l1_rFrameL >= 0 && l1_rFrameU >= 0) ? imin(l1_rFrameL,l1_rFrameU): imax(l1_rFrameL,l1_rFrameU);
  *l1_rFrame = (*l1_rFrame >= 0 && l1_rFrameUR >= 0) ? imin(*l1_rFrame,l1_rFrameUR): imax(*l1_rFrame,l1_rFrameUR);

  if (*l0_rFrame >=0)
    SetMotionVectorPredictor (img, pmvl0, *l0_rFrame, LIST_0, dec_picture->ref_idx, dec_picture->mv, 0, 0, 16, 16);

  if (*l1_rFrame >=0)
    SetMotionVectorPredictor (img, pmvl1, *l1_rFrame, LIST_1, dec_picture->ref_idx, dec_picture->mv, 0, 0, 16, 16);
}


void perform_mc(StorablePicture *dec_picture, struct img_par *img, Macroblock *currMB, int pred_dir, int i, int j, int list_offset, int block_size_x, int block_size_y, int curr_mb_field)
{
  static imgpel tmp_block_l0[MB_BLOCK_SIZE][MB_BLOCK_SIZE];
  static imgpel tmp_block_l1[MB_BLOCK_SIZE][MB_BLOCK_SIZE];

  static int vec1_x=0, vec1_y=0, vec2_x=0, vec2_y=0;
  static int vec1_y_cr = 0, vec2_y_cr = 0;
  static int alpha_l0, alpha_l1, wp_offset;
  static const int mv_mul = 4;
  
  int i4   = img->block_x + i;
  int j4   = img->block_y + j;
  int ioff = (i << 2);
  int joff = (j << 2);         
  
  assert (pred_dir<=2);

  if (pred_dir != 2)
  {
    //===== Single List Prediction =====
    short          ref_idx = dec_picture->ref_idx[LIST_0 + pred_dir][j4][i4];
    short          ref_idx_wp = ref_idx;
    short      ***mv_array = dec_picture->mv[LIST_0 + pred_dir];
    StorablePicture **list = listX[LIST_0 + list_offset + pred_dir];

    vec1_x = i4 * 4 * mv_mul + mv_array[j4][i4][0];
    vec1_y = (img->block_y_aff + j) * 4 * mv_mul + mv_array[j4][i4][1];

    get_block_luma (ref_idx, list, vec1_x, vec1_y, block_size_x, block_size_y, img, tmp_block_l0);

    if (img->apply_weights)
    {
      if (curr_mb_field && ((active_pps->weighted_pred_flag&&(img->type==P_SLICE|| img->type == SP_SLICE))||
         (active_pps->weighted_bipred_idc==1 && (img->type==B_SLICE))))
      {
        ref_idx_wp >>=1;
      }

      alpha_l0  = img->wp_weight[pred_dir][ref_idx_wp][0];
      wp_offset = img->wp_offset[pred_dir][ref_idx_wp][0];

      weighted_mc_prediction(img, LumaComp, block_size_y, block_size_x, joff, ioff, tmp_block_l0, alpha_l0, wp_offset, img->luma_log2_weight_denom, img->max_imgpel_value);
    }
    else
    {
      mc_prediction(img, LumaComp, block_size_y, block_size_x, joff, ioff, tmp_block_l0);
    }

    if ((dec_picture->chroma_format_idc != YUV400) && !IS_INDEPENDENT(img))
    {
      imgpel **curUV;
      int uv;

      int ioff_cr = img->mb_cr_size_x == MB_BLOCK_SIZE ? ioff : ioff >> 1;
      int joff_cr = img->mb_cr_size_y == MB_BLOCK_SIZE ? joff : joff >> 1;
      int block_size_x_cr = img->mb_cr_size_x == MB_BLOCK_SIZE ? block_size_x : block_size_x >> 1;
      int block_size_y_cr = img->mb_cr_size_y == MB_BLOCK_SIZE ? block_size_y : block_size_y >> 1;

      vec1_y_cr = vec1_y + ((active_sps->chroma_format_idc == 1)? list[ref_idx]->chroma_vector_adjustment : 0);

      for(uv=0;uv<2;uv++)
      {
        curUV = dec_picture->imgUV[uv]; 

        get_block_chroma (uv, ref_idx, list, vec1_x, vec1_y_cr, block_size_x_cr, block_size_y_cr, img, tmp_block_l0);

        if (img->apply_weights)
        {
          alpha_l0  = img->wp_weight[pred_dir][ref_idx_wp][uv + 1];
          wp_offset = img->wp_offset[pred_dir][ref_idx_wp][uv + 1];

          weighted_mc_prediction(img, uv + 1, block_size_y_cr, block_size_x_cr, joff_cr, ioff_cr, tmp_block_l0, alpha_l0, wp_offset, img->chroma_log2_weight_denom, img->max_imgpel_value_uv);
        }
        else
        {
          mc_prediction(img, uv + 1, block_size_y_cr, block_size_x_cr, joff_cr, ioff_cr, tmp_block_l0);
        }
      }
    }
  }
  else
  {
    //===== BI-PREDICTION =====
    short ***l0_mv_array = dec_picture->mv[LIST_0];
    short ***l1_mv_array = dec_picture->mv[LIST_1];

    short l0_refframe = dec_picture->ref_idx[LIST_0][j4][i4];
    short l1_refframe = dec_picture->ref_idx[LIST_1][j4][i4];
    short l0_ref_idx  = l0_refframe;
    short l1_ref_idx  = l1_refframe;

    vec1_x = i4 * 4 * mv_mul + l0_mv_array[j4][i4][0];
    vec2_x = i4 * 4 * mv_mul + l1_mv_array[j4][i4][0];

    vec1_y = (img->block_y_aff + j) * 4 * mv_mul + l0_mv_array[j4][i4][1];
    vec2_y = (img->block_y_aff + j) * 4 * mv_mul + l1_mv_array[j4][i4][1];

    get_block_luma(l0_refframe, listX[LIST_0 + list_offset], vec1_x, vec1_y, block_size_x, block_size_y, img, tmp_block_l0);
    get_block_luma(l1_refframe, listX[LIST_1 + list_offset], vec2_x, vec2_y, block_size_x, block_size_y, img, tmp_block_l1);

    if(img->apply_weights)
    {
      int wt_list_offset = (active_pps->weighted_bipred_idc==2)? list_offset : 0;

      // This code existed in the original. Seems pointless but copying it here for reference and in case temporal direct breaks.
      // if (mv_mode==0 && img->direct_spatial_mv_pred_flag==0 ) l1_ref_idx=0;    
      if (((active_pps->weighted_pred_flag&&(img->type==P_SLICE|| img->type == SP_SLICE))||
        (active_pps->weighted_bipred_idc==1 && (img->type==B_SLICE))) && curr_mb_field)
      {
        l0_ref_idx >>=1;
        l1_ref_idx >>=1;
      }

      alpha_l0  =   img->wbp_weight[LIST_0 + wt_list_offset][l0_ref_idx][l1_ref_idx][0];
      alpha_l1  =   img->wbp_weight[LIST_1 + wt_list_offset][l0_ref_idx][l1_ref_idx][0];
      wp_offset = ((img->wp_offset [LIST_0 + wt_list_offset][l0_ref_idx][0] + img->wp_offset[LIST_1 + wt_list_offset][l1_ref_idx][0] + 1) >>1);

      weighted_bi_prediction(img, LumaComp, block_size_y, block_size_x, joff, ioff, tmp_block_l0, tmp_block_l1, alpha_l0, alpha_l1, wp_offset, (img->luma_log2_weight_denom + 1), img->max_imgpel_value);
    }
    else
    {
      bi_prediction(img, LumaComp, block_size_y, block_size_x, joff, ioff, tmp_block_l0, tmp_block_l1);
    }

    if ((dec_picture->chroma_format_idc != YUV400) && !IS_INDEPENDENT(img))
    {
      imgpel **curUV;
      int uv;

      int ioff_cr = img->mb_cr_size_x == MB_BLOCK_SIZE ? ioff : ioff >> 1;
      int joff_cr = img->mb_cr_size_y == MB_BLOCK_SIZE ? joff : joff >> 1;
      int block_size_x_cr = img->mb_cr_size_x == MB_BLOCK_SIZE ? block_size_x : block_size_x >> 1;
      int block_size_y_cr = img->mb_cr_size_y == MB_BLOCK_SIZE ? block_size_y : block_size_y >> 1;

      vec1_y_cr = vec1_y + ((active_sps->chroma_format_idc == 1)? listX[LIST_0 + list_offset][l0_refframe]->chroma_vector_adjustment : 0);
      vec2_y_cr = vec2_y + ((active_sps->chroma_format_idc == 1)? listX[LIST_1 + list_offset][l1_refframe]->chroma_vector_adjustment : 0);

      for(uv=0;uv<2;uv++)
      {
        curUV = dec_picture->imgUV[uv]; 
        get_block_chroma (uv, l0_refframe, listX[LIST_0 + list_offset], vec1_x, vec1_y_cr, block_size_x_cr, block_size_y_cr, img, tmp_block_l0);
        get_block_chroma (uv, l1_refframe, listX[LIST_1 + list_offset], vec2_x, vec2_y_cr, block_size_x_cr, block_size_y_cr, img, tmp_block_l1);

        if(img->apply_weights)
        {
          int wt_list_offset = (active_pps->weighted_bipred_idc==2)? list_offset : 0;

          alpha_l0  =   img->wbp_weight[LIST_0 + wt_list_offset][l0_ref_idx][l1_ref_idx][uv + 1];
          alpha_l1  =   img->wbp_weight[LIST_1 + wt_list_offset][l0_ref_idx][l1_ref_idx][uv + 1];
          wp_offset = ((img->wp_offset [LIST_0 + wt_list_offset][l0_ref_idx][uv + 1] + img->wp_offset[LIST_1 + wt_list_offset][l1_ref_idx][uv + 1] + 1) >>1);

          weighted_bi_prediction(img, uv + 1, block_size_y_cr, block_size_x_cr, joff_cr, ioff_cr, tmp_block_l0, tmp_block_l1, alpha_l0, alpha_l1, wp_offset, (img->chroma_log2_weight_denom + 1), img->max_imgpel_value_uv);
        }
        else
        {
          bi_prediction(img, uv + 1, block_size_y_cr, block_size_x_cr, joff_cr, ioff_cr, tmp_block_l0, tmp_block_l1);
        }
      }
    }      
  }
}

