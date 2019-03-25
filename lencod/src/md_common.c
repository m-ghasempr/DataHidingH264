/*!
 ***************************************************************************
 * \file md_common.c
 *
 * \brief
 *    Mode Decision common functions
 *
 * \author
 *    - Alexis Michael Tourapis    <alexismt@ieee.org>
 * \date
 *    04. October 2008
 **************************************************************************
 */

#include "global.h"
#include "image.h"

/*!
 *************************************************************************************
 * \brief
 *    Copy rdo mv info to 16xN picture mv buffer
 *************************************************************************************
 */
static inline void CopyMVBlock16(short ***enc_mv, short ***rdo_mv, int start, int end)
{
  int j;
  for (j = start; j < end; j++)
  {
    memcpy(enc_mv[img->block_y + j][img->block_x], rdo_mv[j][0], BLOCK_MULTIPLE * 2 * sizeof(short));
  }
}

/*!
 *************************************************************************************
 * \brief
 *    Reset corresponding 16xN block picture mv buffer
 *************************************************************************************
 */
static inline void ResetMVBlock16(short ***enc_mv, int start, int end)
{
  int j;
  for (j = img->block_y + start; j < img->block_y + end; j++)
  {
    memset(enc_mv[j][img->block_x], 0, BLOCK_MULTIPLE * 2 * sizeof(short));
  }
}

/*!
 *************************************************************************************
 * \brief
 *    Reset corresponding 16xN block picture ref buffer (intra)
 *************************************************************************************
 */
static inline void ResetRefBlock16(char **enc_ref, int start, int end)
{
  int j;
  for (j = img->block_y + start; j < img->block_y + end; j++)
  {
    memset(&enc_ref[j][img->block_x], -1, BLOCK_MULTIPLE * sizeof(char));
  }
}

/*!
 *************************************************************************************
 * \brief
 *    Copy rdo mv info to 8xN picture mv buffer
 *************************************************************************************
 */
static inline void CopyMVBlock8(short ***enc_mv, short ***rdo_mv, int start, int end, int offset)
{
  int j;
  int block_x = img->block_x + offset;
  for (j = start; j < end; j++)
  {
    memcpy(enc_mv[img->block_y + j][block_x], rdo_mv[j][offset], 2 * 2 * sizeof(short));
  }
}

/*!
 *************************************************************************************
 * \brief
 *    Reset corresponding 8xN block picture mv buffer
 *************************************************************************************
 */
static inline void ResetMVBlock8(short ***enc_mv, int start, int end, int offset)
{
  int j;
  int block_x = img->block_x + offset;
  for (j = img->block_y + start; j < img->block_y + end; j++)
  {
    memset(enc_mv[j][block_x], 0, 2 * 2 * sizeof(short));
  }
}

/*!
 *************************************************************************************
 * \brief
 *    Reset corresponding 8xN block picture ref buffer (intra)
 *************************************************************************************
 */
static inline void ResetRefBlock8(char **enc_ref, int start, int end, int offset)
{
  int j;
  int block_x = img->block_x + offset;
  for (j = img->block_y + start; j < img->block_y + end; j++)
  {
    memset(&enc_ref[j][block_x], -1, 2 * sizeof(char));
  }
}

/*!
 *************************************************************************************
 * \brief
 *    Sets motion vectors for a macroblock in a P Slice
 *************************************************************************************
 */
static void SetMotionVectorsMBPSlice (ImageParameters *img, PicMotionParams *motion, Macroblock* currMB)
{
  short *****all_mv  = img->all_mv[LIST_0];
  int  l0_ref, mode8;
  short ***rdo_mv = motion->mv[LIST_0];

  if (img->MbaffFrameFlag || (params->UseRDOQuant && params->RDOQ_QP_Num > 1))
  {
    memcpy(&rdopt->all_mv[LIST_0][0][0][0][0][0], &all_mv[0][0][0][0][0], img->max_num_references * 9 * MB_BLOCK_PARTITIONS * 2 * sizeof(short));
  }

  if (currMB->mb_type == PSKIP) // Skip mode
  {
    CopyMVBlock16(rdo_mv, all_mv[0][PSKIP], 0, 4);
  }
  else if (currMB->mb_type == P16x16) // 16x16
  {
    l0_ref = motion->ref_idx[LIST_0][img->block_y][img->block_x];
    CopyMVBlock16(rdo_mv, all_mv[l0_ref][P16x16], 0, 4);
  }
  else if (currMB->mb_type == P16x8) // 16x8
  {
    l0_ref = motion->ref_idx[LIST_0][img->block_y][img->block_x];
    CopyMVBlock16(rdo_mv, all_mv[l0_ref][P16x8], 0, 2);

    l0_ref = motion->ref_idx[LIST_0][img->block_y + 2][img->block_x];
    CopyMVBlock16(rdo_mv, all_mv[l0_ref][P16x8], 2, 4);
  }
  else if (currMB->mb_type == P8x16) // 8x16
  {
    l0_ref = motion->ref_idx[LIST_0][img->block_y][img->block_x];
    CopyMVBlock8(rdo_mv, all_mv[l0_ref][P8x16], 0, 4, 0);
    
    l0_ref = motion->ref_idx[LIST_0][img->block_y][img->block_x + 2];
    CopyMVBlock8(rdo_mv, all_mv[l0_ref][P8x16], 0, 4, 2);
  } 
  else if (currMB->mb_type == P8x8) // 8x8
  {
    mode8 = currMB->b8mode[0];
    l0_ref = motion->ref_idx[LIST_0][img->block_y    ][img->block_x    ];     
    CopyMVBlock8(rdo_mv, all_mv[l0_ref][mode8], 0, 2, 0);
    
    mode8 = currMB->b8mode[1];
    l0_ref = motion->ref_idx[LIST_0][img->block_y    ][img->block_x + 2];
    CopyMVBlock8(rdo_mv, all_mv[l0_ref][mode8], 0, 2, 2);
    
    mode8 = currMB->b8mode[2];
    l0_ref = motion->ref_idx[LIST_0][img->block_y + 2][img->block_x    ];
    CopyMVBlock8(rdo_mv, all_mv[l0_ref][mode8], 2, 4, 0);

    mode8 = currMB->b8mode[3];
    l0_ref = motion->ref_idx[LIST_0][img->block_y + 2][img->block_x + 2];
    CopyMVBlock8(rdo_mv, all_mv[l0_ref][mode8], 2, 4, 2);
  }
  else // Intra modes
  {
    ResetMVBlock16 (rdo_mv, 0, 4);
    ResetRefBlock16(motion->ref_idx[LIST_0], 0, 4);
  }
}

/*!
 *************************************************************************************
 * \brief
 *    Sets motion vectors for a 16x8 partition in a B Slice
 *************************************************************************************
 */
static void SetMVBSlice16x8(PicMotionParams *motion, Macroblock* currMB, int pos)
{
  int l0_ref, l1_ref;
  int pdir = currMB->b8pdir[pos];
  if (pdir == LIST_0) 
  {
    l0_ref = motion->ref_idx[LIST_0][img->block_y + pos][img->block_x];
    CopyMVBlock16(motion->mv[LIST_0], img->all_mv [LIST_0][l0_ref][P16x8], pos, pos + 2);
    ResetMVBlock16 (motion->mv[LIST_1], pos, pos + 2);
    ResetRefBlock16(motion->ref_idx[LIST_1], pos, pos + 2);
  }
  else if (pdir == LIST_1)
  {
    ResetMVBlock16 (motion->mv[LIST_0], pos, pos + 2);
    ResetRefBlock16(motion->ref_idx[LIST_0], pos, pos + 2);
    l1_ref = motion->ref_idx[LIST_1][img->block_y + pos][img->block_x];
    CopyMVBlock16 (motion->mv[LIST_1], img->all_mv [LIST_1][l1_ref][P16x8], pos, pos + 2);
  }
  else
  {
    int bipred_me = currMB->bipred_me[pos];
    short ******all_mv = bipred_me ? img->bipred_mv[bipred_me - 1]: img->all_mv;
    l0_ref = motion->ref_idx[LIST_0][img->block_y + pos][img->block_x];
    CopyMVBlock16 (motion->mv[LIST_0], all_mv [LIST_0][l0_ref][P16x8], pos, pos + 2);
    
    l1_ref = motion->ref_idx[LIST_1][img->block_y + pos][img->block_x];
    CopyMVBlock16 (motion->mv[LIST_1], all_mv [LIST_1][l1_ref][P16x8], pos, pos + 2);

    if (bipred_me && (img->MbaffFrameFlag || (params->UseRDOQuant && params->RDOQ_QP_Num > 1)))
    {
      memcpy(rdopt->all_mv [LIST_0][l0_ref][P16x8][pos][0], all_mv [LIST_0][l0_ref][P16x8][pos][0], 2 * BLOCK_MULTIPLE * 2 * sizeof(short));
      memcpy(rdopt->all_mv [LIST_1][l1_ref][P16x8][pos][0], all_mv [LIST_1][l1_ref][P16x8][pos][0], 2 * BLOCK_MULTIPLE * 2 * sizeof(short));
    }
  }
}

/*!
 *************************************************************************************
 * \brief
 *    Sets motion vectors for a 8x16 partition in a B Slice
 *************************************************************************************
 */
static void SetMVBSlice8x16(PicMotionParams *motion, Macroblock* currMB, int pos)
{
  int l0_ref, l1_ref;
  int pdir = currMB->b8pdir[pos >> 1];
  if (pdir == LIST_0) 
  {
    l0_ref = motion->ref_idx[LIST_0][img->block_y][img->block_x + pos];
    CopyMVBlock8 (motion->mv[LIST_0], img->all_mv [LIST_0][l0_ref][P8x16], 0, 4, pos);
    ResetMVBlock8(motion->mv[LIST_1], 0, 4, pos);
    ResetRefBlock8(motion->ref_idx[LIST_1], 0, 4, pos);
  }
  else if (pdir == LIST_1)
  {
    ResetMVBlock8(motion->mv[LIST_0], 0, 4, pos);
    ResetRefBlock8(motion->ref_idx[LIST_0], 0, 4, pos);
    l1_ref = motion->ref_idx[LIST_1][img->block_y][img->block_x + pos];    
    CopyMVBlock8 (motion->mv[LIST_1], img->all_mv [LIST_1][l1_ref][P8x16], 0, 4, pos);
  }
  else
  {
    int bipred_me = currMB->bipred_me[pos >> 1];
    short ******all_mv = bipred_me ? img->bipred_mv[bipred_me - 1]: img->all_mv;
    l0_ref = motion->ref_idx[LIST_0][img->block_y][img->block_x + pos];
    CopyMVBlock8(motion->mv [LIST_0], all_mv [LIST_0][l0_ref][P8x16], 0, 4, pos);
    l1_ref = motion->ref_idx[LIST_1][img->block_y][img->block_x + pos];    
    CopyMVBlock8(motion->mv [LIST_1], all_mv [LIST_1][l1_ref][P8x16], 0, 4, pos);

    if (bipred_me && (img->MbaffFrameFlag || (params->UseRDOQuant && params->RDOQ_QP_Num > 1)))
    {
      int j;
      for (j = 0; j < BLOCK_MULTIPLE; j++)
      {
        memcpy(rdopt->all_mv [LIST_0][l0_ref][P8x16][j][pos], all_mv [LIST_0][l0_ref][P8x16][j][pos], 2 * 2 * sizeof(short));
      }
      for (j = 0; j < BLOCK_MULTIPLE; j++)
      {
        memcpy(rdopt->all_mv [LIST_1][l1_ref][P8x16][j][pos], all_mv [LIST_1][l1_ref][P8x16][j][pos], 2 * 2 * sizeof(short));
      }
    }
  }
}

/*!
 *************************************************************************************
 * \brief
 *    Sets motion vectors for a 8x8 partition in a B Slice
 *************************************************************************************
 */
static void SetMVBSlice8x8(PicMotionParams *motion, Macroblock* currMB, int pos_y, int pos_x)
{
  int l0_ref, l1_ref;
  int pos = pos_y + (pos_x >> 1);
  int pdir = currMB->b8pdir[pos];
  int mode = currMB->b8mode[pos];
  int block_y = img->block_y + pos_y;
  int block_x = img->block_x + pos_x;

  if (pdir == LIST_0) 
  {
    l0_ref = motion->ref_idx[LIST_0][block_y][block_x];
    CopyMVBlock8 (motion->mv[LIST_0], img->all_mv [LIST_0][l0_ref][mode], pos_y, pos_y + 2, pos_x);
    ResetMVBlock8(motion->mv[LIST_1], pos_y, pos_y + 2, pos_x);
    ResetRefBlock8(motion->ref_idx[LIST_1], pos_y, pos_y + 2, pos_x);
  }
  else if (pdir == LIST_1)
  {
    ResetMVBlock8(motion->mv[LIST_0], pos_y, pos_y + 2, pos_x);
    ResetRefBlock8(motion->ref_idx[LIST_0], pos_y, pos_y + 2, pos_x);
    l1_ref = motion->ref_idx[LIST_1][block_y][block_x];    
    CopyMVBlock8 (motion->mv[LIST_1], img->all_mv [LIST_1][l1_ref][mode], pos_y, pos_y + 2, pos_x);
  }
  else if (pdir == BI_PRED)
  {
    int bipred_me = currMB->bipred_me[pos];
    short ******all_mv = bipred_me ? img->bipred_mv[bipred_me - 1]: img->all_mv;
    l0_ref = motion->ref_idx[LIST_0][block_y][block_x];
    CopyMVBlock8(motion->mv [LIST_0], all_mv [LIST_0][l0_ref][mode], pos_y, pos_y + 2, pos_x);
    
    l1_ref = motion->ref_idx[LIST_1][block_y][block_x];    
    CopyMVBlock8(motion->mv [LIST_1], all_mv [LIST_1][l1_ref][mode], pos_y, pos_y + 2, pos_x);

    if (bipred_me && (img->MbaffFrameFlag || (params->UseRDOQuant && params->RDOQ_QP_Num > 1)))
    {
      memcpy(rdopt->all_mv [LIST_0][l0_ref][mode][pos_y    ][pos_x], all_mv [LIST_0][l0_ref][mode][pos_y    ][pos_x], 2 * 2 * sizeof(short));
      memcpy(rdopt->all_mv [LIST_0][l0_ref][mode][pos_y + 1][pos_x], all_mv [LIST_0][l0_ref][mode][pos_y + 1][pos_x], 2 * 2 * sizeof(short));
      memcpy(rdopt->all_mv [LIST_1][l1_ref][mode][pos_y    ][pos_x], all_mv [LIST_1][l1_ref][mode][pos_y    ][pos_x], 2 * 2 * sizeof(short));
      memcpy(rdopt->all_mv [LIST_1][l1_ref][mode][pos_y + 1][pos_x], all_mv [LIST_1][l1_ref][mode][pos_y + 1][pos_x], 2 * 2 * sizeof(short));
    }
  }
  else // Invalid direct modes (out of range mvs). Adding this here for precaution purposes.
  {
    ResetMVBlock8(motion->mv[LIST_0], pos_y, pos_y + 2, pos_x);
    ResetRefBlock8(motion->ref_idx[LIST_0], pos_y, pos_y + 2, pos_x);
    ResetMVBlock8(motion->mv[LIST_1], pos_y, pos_y + 2, pos_x);
    ResetRefBlock8(motion->ref_idx[LIST_1], pos_y, pos_y + 2, pos_x);
  }
}

/*!
 *************************************************************************************
 * \brief
 *    Sets motion vectors for a macroblock in a P Slice
 *************************************************************************************
 */
static void SetMotionVectorsMBBSlice (ImageParameters *img, PicMotionParams *motion, Macroblock* currMB)
{
  int  l0_ref, l1_ref;

  // copy all the motion vectors into rdopt structure
  // Can simplify this by copying the MV's of the best mode (TBD)
  // Should maybe add code to check for Intra only profiles

  if (img->MbaffFrameFlag || (params->UseRDOQuant && params->RDOQ_QP_Num > 1))
  {
    memcpy(&rdopt->all_mv[LIST_0][0][0][0][0][0], &img->all_mv[LIST_0][0][0][0][0][0], 2 * img->max_num_references * 9 * MB_BLOCK_PARTITIONS * 2 * sizeof(short));
  }

  if (currMB->mb_type == P16x16) // 16x16
  {
    int pdir = currMB->b8pdir[0];      
    if (pdir == LIST_0) 
    {
      l0_ref = motion->ref_idx[LIST_0][img->block_y][img->block_x];
      CopyMVBlock16 (motion->mv[LIST_0], img->all_mv [LIST_0][l0_ref][P16x16], 0, 4);
      ResetMVBlock16(motion->mv[LIST_1], 0, 4);
    }
    else if (pdir == LIST_1)
    {
      l1_ref = motion->ref_idx[LIST_1][img->block_y][img->block_x];
      ResetMVBlock16(motion->mv[LIST_0], 0, 4);
      CopyMVBlock16 (motion->mv[LIST_1], img->all_mv [LIST_1][l1_ref][P16x16], 0, 4);
    }
    else
    {
      int bipred_me = currMB->bipred_me[0];
      short ******all_mv  = bipred_me ? img->bipred_mv[bipred_me - 1]: img->all_mv;
      l0_ref = motion->ref_idx[LIST_0][img->block_y][img->block_x];
      CopyMVBlock16 (motion->mv[LIST_0], all_mv [LIST_0][l0_ref][P16x16], 0, 4);
      l1_ref = motion->ref_idx[LIST_1][img->block_y][img->block_x];
      CopyMVBlock16 (motion->mv[LIST_1], all_mv [LIST_1][l1_ref][P16x16], 0, 4);

      // Is this necessary here? Can this be moved somewhere else?
      if (bipred_me && (img->MbaffFrameFlag || (params->UseRDOQuant && params->RDOQ_QP_Num > 1)))
      {
        memcpy(rdopt->all_mv [LIST_0][l0_ref][P16x16][0][0], all_mv [LIST_0][l0_ref][P16x16][0][0], MB_BLOCK_PARTITIONS * 2 * sizeof(short));
        memcpy(rdopt->all_mv [LIST_1][l1_ref][P16x16][0][0], all_mv [LIST_1][l1_ref][P16x16][0][0], MB_BLOCK_PARTITIONS * 2 * sizeof(short));
      }
    }
  }
  else if (currMB->mb_type == P16x8) // 16x8
  {        
    SetMVBSlice16x8(motion, currMB, 0);
    SetMVBSlice16x8(motion, currMB, 2);
  }
  else if (currMB->mb_type == P8x16) // 16x8
  {        
    SetMVBSlice8x16(motion, currMB, 0);
    SetMVBSlice8x16(motion, currMB, 2);
  }
  else if (currMB->mb_type == P8x8 || currMB->mb_type == BSKIP_DIRECT) // 8x8 & Direct/SKIP
  { 
    SetMVBSlice8x8(motion, currMB, 0, 0);
    SetMVBSlice8x8(motion, currMB, 0, 2);
    SetMVBSlice8x8(motion, currMB, 2, 0);
    SetMVBSlice8x8(motion, currMB, 2, 2);
  }
  else // Intra modes
  {
    ResetMVBlock16(motion->mv[LIST_0], 0, 4);
    ResetRefBlock16(motion->ref_idx[LIST_0], 0, 4);
    ResetMVBlock16(motion->mv[LIST_1], 0, 4);
    ResetRefBlock16(motion->ref_idx[LIST_1], 0, 4);
  }
}

/*!
 *************************************************************************************
 * \brief
 *    Sets motion vectors for a macroblock
 *************************************************************************************
 */
void SetMotionVectorsMB (ImageParameters *img, PicMotionParams *motion, Macroblock *currMB)
{
  if (img->type == B_SLICE)
  {
    SetMotionVectorsMBBSlice (img, motion, currMB);
  }
  else if (img->type == P_SLICE || img->type == SP_SLICE)
  {
    SetMotionVectorsMBPSlice (img, motion, currMB);
  }
}


