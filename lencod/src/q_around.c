
/*!
 *************************************************************************************
 * \file q_around.c
 *
 * \brief
 *    Quantization Adaptive Rounding (JVT-N011)
 * \author
 *    - Alexis Michael Tourapis    <alexismt@ieee.org>
 *    - Ehsan Maani                <emaan@dolby.com>
 *
 *************************************************************************************
 */

#include "global.h"
#include "memalloc.h"
#include "q_offsets.h"
#include "q_around.h"

static const int AdaptRndCrPos[2][5] =
{
  //  P,   B,   I,  SP,  SI
  {   4,   7,   1,   4,   1}, // Intra MB
  {  10,  13,  10,  10,  10}  // Inter MB
};

static const int AdaptRndPos[4][5] =
{
  //  P,   B,   I,  SP,  SI
  {   3,   6,   0,   3,   0}, // 4x4 Intra MB
  {   1,   2,   0,   1,   2}, // 8x8 Intra MB
  {   9,  12,   9,   9,   9}, // 4x4 Inter MB
  {   3,   4,   3,   3,   3}, // 8x8 Inter MB
};

int   **fadjust8x8 = NULL, **fadjust4x4 = NULL, ***fadjust4x4Cr = NULL, ***fadjust8x8Cr = NULL;


/*!
************************************************************************
* \brief
*    copy data from to a temporary location
************************************************************************
*/
void store_adaptive_rounding_4x4 (ImageParameters *img, int****ARCofAdj, int mode, int block_y, int block_x)
{
  int j;

  for (j = block_y; j < block_y + 4; j++)
    memcpy(&ARCofAdj[0][DUMMY][j][block_x], &ARCofAdj[0][mode][j][block_x], BLOCK_SIZE * sizeof(int));

  if (img->P444_joined)
  {
    for (j = block_y; j < block_y + 4; j++)
    {
      memcpy(&ARCofAdj[1][DUMMY][j][block_x],&ARCofAdj[1][mode][j][block_x], BLOCK_SIZE * sizeof(int));
    }
    for (j = block_y; j < block_y + 4; j++)
    {
      memcpy(&ARCofAdj[2][DUMMY][j][block_x],&ARCofAdj[2][mode][j][block_x], BLOCK_SIZE * sizeof(int));
    }
  }
}


/*!
************************************************************************
* \brief
*    copy data from a temporary location to the actual location
************************************************************************
*/
void update_adaptive_rounding_4x4 (ImageParameters *img, int****ARCofAdj , int mode, int block_y, int block_x)
{
  int j;

  for (j = block_y; j < block_y + BLOCK_SIZE; j++)
    memcpy (&ARCofAdj[0][mode][j][block_x],&ARCofAdj[0][DUMMY][j][block_x], BLOCK_SIZE * sizeof(int));

  if (img->P444_joined)
  {
    for (j = 0; j < block_y + BLOCK_SIZE; j++)
    {
      memcpy (&ARCofAdj[1][mode][j][block_x], &ARCofAdj[1][DUMMY][j][block_x], BLOCK_SIZE * sizeof(int));
    }
    for (j = 0; j < block_y + BLOCK_SIZE; j++)
    {
      memcpy (&ARCofAdj[2][mode][j][block_x], &ARCofAdj[2][DUMMY][j][block_x], BLOCK_SIZE * sizeof(int));
    }
  }
}

/*!
************************************************************************
* \brief
*    copy data from to a temporary location
************************************************************************
*/
void store_adaptive_rounding_16x16 (ImageParameters *img, int****ARCofAdj, int mode)
{
    memcpy(&ARCofAdj[0][DUMMY][0][0], &ARCofAdj[0][mode][0][0], MB_PIXELS * sizeof(int));
    if (img->P444_joined)
    {
      memcpy(&ARCofAdj[1][DUMMY][0][0], &ARCofAdj[1][mode][0][0], MB_PIXELS * sizeof(int));
      memcpy(&ARCofAdj[2][DUMMY][0][0], &ARCofAdj[2][mode][0][0], MB_PIXELS * sizeof(int));
    }
}


/*!
************************************************************************
* \brief
*    copy data from a temporary location to the actual location
************************************************************************
*/
void update_adaptive_rounding_16x16(ImageParameters *img, int****ARCofAdj , int mode)
{  
  memcpy (&ARCofAdj[0][mode][0][0],&ARCofAdj[0][DUMMY][0][0], MB_PIXELS * sizeof(int));
  if (img->P444_joined)
  {
    memcpy (&ARCofAdj[1][mode][0][0], &ARCofAdj[1][DUMMY][0][0], MB_PIXELS * sizeof(int));
    memcpy (&ARCofAdj[2][mode][0][0], &ARCofAdj[2][DUMMY][0][0], MB_PIXELS * sizeof(int));
  }
}




/*!
************************************************************************
* \brief
*    update rounding offsets based on JVT-N011
************************************************************************
*/
void update_offset_params(Macroblock *currMB, int mode, int luma_transform_size_8x8_flag)
{
  int is_inter = (mode != I4MB)&&(mode != I16MB) && (mode != I8MB);
  int luma_pos = AdaptRndPos[(is_inter<<1) + luma_transform_size_8x8_flag][img->type];
  int i,j;
  int qp = currMB->qp + img->bitdepth_luma_qp_scale;
  int cur_qp = params->AdaptRoundingFixed ? 0 : qp;
  int temp = 0;
  int offsetRange = 1 << (OffsetBits - 1);
  int blk_mask = 0x03 + (luma_transform_size_8x8_flag<<2);
  int blk_shift = 2 + luma_transform_size_8x8_flag;
  short **offsetList = luma_transform_size_8x8_flag ? OffsetList8x8[cur_qp] : OffsetList4x4[cur_qp];

  int **fAdjust = luma_transform_size_8x8_flag ? img->ARCofAdj8x8[0][mode] : img->ARCofAdj4x4[0][mode];
  
  if (mode == IPCM) return;

  if( (active_sps->chroma_format_idc == YUV444) && IS_INDEPENDENT(params) )
  {
    if( luma_transform_size_8x8_flag )  // 8x8
      luma_pos += 5 * img->colour_plane_id;
    else  // 4x4
      luma_pos += img->colour_plane_id;
  }
 
  for (j=0; j < MB_BLOCK_SIZE; j++)
  {
    int j_pos = ((j & blk_mask)<<blk_shift);
    for (i=0; i < MB_BLOCK_SIZE; i++)
    {
      temp = j_pos + (i & blk_mask);
      offsetList[luma_pos][temp] += fAdjust[j][i];
      offsetList[luma_pos][temp]  = iClip3(0, offsetRange, offsetList[luma_pos][temp]);
    }
  }

  if(img->P444_joined)
  { 
    int **fAdjustCbCr;
    int uv;

    for(uv = 0; uv < 2; uv++)
    {
      luma_pos = AdaptRndPos[(is_inter<<1) + luma_transform_size_8x8_flag][img->type];
      fAdjustCbCr = luma_transform_size_8x8_flag ? img->ARCofAdj8x8[uv + 1][mode] : img->ARCofAdj4x4[uv + 1][mode];
      if(luma_transform_size_8x8_flag )  // 8x8
        luma_pos += 5 * (uv+1);
      else  // 4x4
        luma_pos += (uv+1);
      for (j=0; j < MB_BLOCK_SIZE; j++)
      {
        int j_pos = ((j & blk_mask)<<blk_shift);
        for (i=0; i < MB_BLOCK_SIZE; i++)
        {
          temp = j_pos + (i & blk_mask);
          offsetList[luma_pos][temp] += fAdjustCbCr[j][i];
          offsetList[luma_pos][temp] = iClip3(0,offsetRange,offsetList[luma_pos][temp]);
        }
      }
    }
  }

  if ((params->AdaptRndChroma) && (img->yuv_format == YUV420 || img->yuv_format == YUV422 ))
  {  
    int u_pos = AdaptRndCrPos[is_inter][img->type];
    int v_pos = u_pos + 1;
    int jpos;

    //int ***fAdjustCr = is_inter ? bestOffset.InterFAdjust4x4Cr : bestOffset.IntraFAdjust4x4Cr;

    int **fAdjustCb = (luma_transform_size_8x8_flag && mode == P8x8 )? img->ARCofAdj4x4[1][4] : img->ARCofAdj4x4[1][mode];
    int **fAdjustCr = (luma_transform_size_8x8_flag && mode == P8x8 )? img->ARCofAdj4x4[2][4] : img->ARCofAdj4x4[2][mode];

  
    for (j = 0; j < img->mb_cr_size_y; j++)
    {
      jpos = ((j & 0x03)<<2);
      for (i = 0; i < img->mb_cr_size_x; i++)
      {
        temp = jpos + (i & 0x03);
        OffsetList4x4[cur_qp][u_pos][temp] += fAdjustCb[j][i];
        OffsetList4x4[cur_qp][u_pos][temp]  = iClip3(0,offsetRange,OffsetList4x4[cur_qp][u_pos][temp]);
        OffsetList4x4[cur_qp][v_pos][temp] += fAdjustCr[j][i];
        OffsetList4x4[cur_qp][v_pos][temp]  = iClip3(0,offsetRange,OffsetList4x4[cur_qp][v_pos][temp]);
      }
    }
  }
}


/*!
************************************************************************
* \brief
*    resets adaptive rounding params to zero
************************************************************************
*/
void reset_adaptive_rounding()
{
  int maxplane;
  maxplane = (img->yuv_format == YUV400)?  1 : 3;
  memset(&img->ARCofAdj4x4[0][0][0][0], 0, maxplane * MAXMODE * MB_PIXELS * sizeof (int));
  maxplane = (img->yuv_format == YUV400)?  1 : (img->P444_joined ? 3 : 1);
  memset(&img->ARCofAdj8x8[0][0][0][0], 0, maxplane * MAXMODE * MB_PIXELS * sizeof (int));
}

/*!
************************************************************************
* \brief
*    resets adaptive rounding params for the direct mode to zero
************************************************************************
*/
void reset_adaptive_rounding_direct()
{
  int maxplane, pl;
  maxplane = (img->yuv_format == YUV400)?  1 : 3;
  for (pl = 0; pl < maxplane; pl++)
    memset(&img->ARCofAdj4x4[pl][DUMMY][0][0], 0, MB_PIXELS * sizeof (int));

  maxplane = (img->yuv_format == YUV400)?  1 : (img->P444_joined ? 3 : 1);
  for (pl = 0; pl < maxplane; pl++)
    memset(&img->ARCofAdj8x8[pl][DUMMY][0][0], 0, MB_PIXELS * sizeof (int));
 
}

/*!
************************************************************************
* \brief
*    copy data from a temporary location to the actual location for mode P8x8
************************************************************************
*/
void update_adaptive_rounding_8x8(RD_8x8DATA* dataTr, int**** ARCofAdj)
{
  int b8, pl, plane, j0, i0, j;

  if (params->AdaptiveRounding)
  {
    plane = (img->P444_joined) ?  3 : 1;
    for (pl = 0; pl < plane; pl++)
    {
      for (b8 = 0; b8 < 4; b8++)
      {
        j0   = 8*(b8 >> 1);
        i0   = 8*(b8 & 0x01);
        if (dataTr->part8x8mode[b8])
        {
          for (j = j0; j < j0 + BLOCK_SIZE_8x8; j++)
          {
            memcpy(&ARCofAdj[pl][P8x8][j][i0], &ARCofAdj[pl][dataTr->part8x8mode[b8]][j][i0], BLOCK_SIZE_8x8 * sizeof(int));
          }
        }
      }
    }
  }
}


