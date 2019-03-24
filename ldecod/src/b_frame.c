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
 *************************************************************************************
 * \file b_frame.c
 *
 * \brief
 *    B picture decoding
 *
 * \author
 *    Main contributors (see contributors.h for copyright, address and affiliation details)
 *    - Byeong-Moon Jeon                <jeonbm@lge.com>
 *    - Yoon-Seong Soh                  <yunsung@lge.com>
 *    - Thomas Wedi                     <wedi@tnt.uni-hannover.de>
 *************************************************************************************
 */

#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "global.h"
#include "mbuffer.h"
#include "b_frame.h"

#define POS 0

/*!
 ************************************************************************
 * \brief
 *    Copy decoded P frame to temporary image array
 ************************************************************************
 */
void copy_Pframe(struct img_par *img, int postfilter)
{
  int i,j;

  /*
   * the mmin, mmax macros are taken out, because it makes no sense due to limited range of data type
   */

  if(postfilter)
  {
    for(i=0;i<img->height;i++)
      for(j=0;j<img->width;j++)
      {
        imgY_prev[i][j] = imgY_pf[i][j];
      }
    for(i=0;i<img->height_cr;i++)
      for(j=0;j<img->width_cr;j++)
      {
        imgUV_prev[0][i][j] = imgUV_pf[0][i][j];
      }
    for(i=0;i<img->height_cr;i++)
      for(j=0;j<img->width_cr;j++)
      {
        imgUV_prev[1][i][j] = imgUV_pf[1][i][j];
      }
  }
  else
  {
    for(i=0;i<img->height;i++)
      for(j=0;j<img->width;j++)
      {
        imgY_prev[i][j] = imgY[i][j];
      }
    for(i=0;i<img->height_cr;i++)
      for(j=0;j<img->width_cr;j++)
      {
        imgUV_prev[0][i][j] = imgUV[0][i][j];
      }
    for(i=0;i<img->height_cr;i++)
      for(j=0;j<img->width_cr;j++)
      {
        imgUV_prev[1][i][j] = imgUV[1][i][j];
      }
  }
}




/*!
 ************************************************************************
 * \brief
 *    init macroblock B frames
 ************************************************************************
 */
void init_macroblock_Bframe(struct img_par *img)
{
  int i,j,k;
  int img_block_y;
  int j2 = img->block_y; 
  Macroblock *currMB = &img->mb_data[img->current_mb_nr];

  if (img->mb_frame_field_flag)
    j2 = img->block_y/2 - 2*(img->current_mb_nr%2);

  // reset vectors and pred. modes
  for (i=0;i<BLOCK_SIZE;i++)
  {
    for(j=0;j<BLOCK_SIZE;j++)
    {
      img->fw_mv_frm[img->block_x+i+4][img->block_y+j][0]=img->fw_mv_frm[img->block_x+i+4][img->block_y+j][1]=0;
      img->bw_mv_frm[img->block_x+i+4][img->block_y+j][0]=img->bw_mv_frm[img->block_x+i+4][img->block_y+j][1]=0;
      img->dfMV [img->block_x+i+4][img->block_y+j][0]=img->dfMV[img->block_x+i+4][img->block_y+j][1]=0;
      img->dbMV [img->block_x+i+4][img->block_y+j][0]=img->dbMV[img->block_x+i+4][img->block_y+j][1]=0;
      img->ipredmode[img->block_x+i+1][img->block_y+j+1] = 0;
      if (img->mb_frame_field_flag)
      {
        img->dfMV_top[img->block_x+i+4][img->block_y+j][0]=img->dfMV_top[img->block_x+i+4][img->block_y+j][1]=0;
        img->dbMV_top[img->block_x+i+4][img->block_y+j][0]=img->dbMV_top[img->block_x+i+4][img->block_y+j][1]=0;
        
        img->dfMV_bot[img->block_x+i+4][img->block_y+j][0]=img->dfMV_bot[img->block_x+i+4][img->block_y+j][1]=0;
        img->dbMV_bot[img->block_x+i+4][img->block_y+j][0]=img->dbMV_bot[img->block_x+i+4][img->block_y+j][1]=0;
      }
      if (img->current_mb_nr%2==0 && img->mb_frame_field_flag)
      {
        img->ipredmode_top[img->block_x+i+1][j2+j+1] = 0;
        img->fw_mv_top[img->block_x+i+4][j2+j][0]=img->fw_mv_top[img->block_x+i+4][j2+j][1]=0;
        img->bw_mv_top[img->block_x+i+4][j2+j][0]=img->bw_mv_top[img->block_x+i+4][j2+j][1]=0;
      }
      else if (img->mb_frame_field_flag)
      {
        img->ipredmode_bot[img->block_x+i+1][j2+j+1] = 0;
        img->fw_mv_bot[img->block_x+i+4][j2+j][0]=img->fw_mv_bot[img->block_x+i+4][j2+j][1]=0;
        img->bw_mv_bot[img->block_x+i+4][j2+j][0]=img->bw_mv_bot[img->block_x+i+4][j2+j][1]=0;
      }
      else
      {
        img->fw_mv[img->block_x+i+4][img->block_y+j][0]=img->fw_mv[img->block_x+i+4][img->block_y+j][1]=0;
        img->bw_mv[img->block_x+i+4][img->block_y+j][0]=img->bw_mv[img->block_x+i+4][img->block_y+j][1]=0;
      }
    }
  }

  // Set the reference frame information for motion vector prediction
  if (IS_INTRA (currMB) || IS_DIRECT (currMB))
  {
    if (img->structure != FRAME || img->mb_frame_field_flag == 0)
    {
      for (j=0; j<4; j++)
        for (i=0; i<4; i++)
        {
          img->fw_refFrArr[img->block_y+j][img->block_x+i] = -1;
          img->bw_refFrArr[img->block_y+j][img->block_x+i] = -1;
        }
    }
    else 
    {
      if (img->current_mb_nr%2==0)
      {
        img_block_y   = img->block_y/2;
        for (j=0; j<4; j++)
          for (i=0; i<4; i++)
          {
            img->fw_refFrArr_top[img_block_y+j][img->block_x+i] = -1;
            img->fw_refFrArr_frm[img->block_y+j][img->block_x+i] = -1;
            img->bw_refFrArr_top[img_block_y+j][img->block_x+i] = -1;
            img->bw_refFrArr_frm[img->block_y+j][img->block_x+i] = -1;
          }
      }
      else
      {
        img_block_y   = (img->block_y-4)/2;
        for (j=0; j<4; j++)
          for (i=0; i<4; i++)
          {
            img->fw_refFrArr_bot[img_block_y+j][img->block_x+i] = -1;
            img->fw_refFrArr_frm[img->block_y+j][img->block_x+i] = -1;
            img->bw_refFrArr_bot[img_block_y+j][img->block_x+i] = -1;
            img->bw_refFrArr_frm[img->block_y+j][img->block_x+i] = -1;
          }
      }
    }
  }
  else
  {
    for(j=0;j<4;j++)
      for(i=0;i<4;i++)
      {
        k=2*(j/2)+(i/2);
        if(img->structure != FRAME)
        {
          img->fw_refFrArr[img->block_y+j][img->block_x+i] = ((currMB->b8pdir[k]==0||currMB->b8pdir[k]==2)&&currMB->b8mode[k]!=0?2:-1);
          img->bw_refFrArr[img->block_y+j][img->block_x+i] = ((currMB->b8pdir[k]==1||currMB->b8pdir[k]==2)&&currMB->b8mode[k]!=0?0:-1);
        }
        else if (img->mb_frame_field_flag)
        {
          if (img->current_mb_nr%2==0)
          {
            img_block_y   = img->block_y/2;
            img->fw_refFrArr_top[img_block_y+j][img->block_x+i] = ((currMB->b8pdir[k]==0||currMB->b8pdir[k]==2)&&currMB->b8mode[k]!=0?0:-1);
            img->bw_refFrArr_top[img_block_y+j][img->block_x+i] = ((currMB->b8pdir[k]==1||currMB->b8pdir[k]==2)&&currMB->b8mode[k]!=0?0:-1);
            img->fw_refFrArr_frm[img->block_y+j][img->block_x+i] = ((currMB->b8pdir[k]==0||currMB->b8pdir[k]==2)&&currMB->b8mode[k]!=0?0:-1);
            img->bw_refFrArr_frm[img->block_y+j][img->block_x+i] = ((currMB->b8pdir[k]==1||currMB->b8pdir[k]==2)&&currMB->b8mode[k]!=0?0:-1);
          }
          else
          {
            img_block_y   = (img->block_y-4)/2;
            img->fw_refFrArr_bot[img_block_y+j][img->block_x+i] = ((currMB->b8pdir[k]==0||currMB->b8pdir[k]==2)&&currMB->b8mode[k]!=0?0:-1);
            img->bw_refFrArr_bot[img_block_y+j][img->block_x+i] = ((currMB->b8pdir[k]==1||currMB->b8pdir[k]==2)&&currMB->b8mode[k]!=0?0:-1);
            img->fw_refFrArr_frm[img->block_y+j][img->block_x+i] = ((currMB->b8pdir[k]==0||currMB->b8pdir[k]==2)&&currMB->b8mode[k]!=0?0:-1);
            img->bw_refFrArr_frm[img->block_y+j][img->block_x+i] = ((currMB->b8pdir[k]==1||currMB->b8pdir[k]==2)&&currMB->b8mode[k]!=0?0:-1);
          }
        }
        else
        {
          img->fw_refFrArr[img->block_y+j][img->block_x+i] = ((currMB->b8pdir[k]==0||currMB->b8pdir[k]==2)&&currMB->b8mode[k]!=0?0:-1);
          img->bw_refFrArr[img->block_y+j][img->block_x+i] = ((currMB->b8pdir[k]==1||currMB->b8pdir[k]==2)&&currMB->b8mode[k]!=0?0:-1);
        }
      }
  }
}
