
/*!
 *************************************************************************************
 * \file loopFilter.c
 *
 * \brief
 *    Filter to reduce blocking artifacts on a macroblock level.
 *    The filter strength is QP dependent.
 *
 * \author
 *    Contributors:
 *    - Peter List       Peter.List@t-systems.de:  Original code                                 (13-Aug-2001)
 *    - Jani Lainema     Jani.Lainema@nokia.com:   Some bug fixing, removal of recursiveness     (16-Aug-2001)
 *    - Peter List       Peter.List@t-systems.de:  inplace filtering and various simplifications (10-Jan-2002)
 *    - Anthony Joch     anthony@ubvideo.com:      Simplified switching between filters and
 *                                                 non-recursive default filter.                 (08-Jul-2002)
 *    - Cristina Gomila  cristina.gomila@thomson.net: Simplification of the chroma deblocking
 *                                                    from JVT-E089                              (21-Nov-2002)
 *    - Alexis Michael Tourapis atour@dolby.com:   Speed/Architecture improvements               (08-Feb-2007)
 *************************************************************************************
 */

#include "global.h"
#include "image.h"
#include "mb_access.h"
//#include "loopfilter.h"

/*********************************************************************************************************/

// NOTE: In principle, the alpha and beta tables are calculated with the formulas below
//       Alpha( qp ) = 0.8 * (2^(qp/6)  -  1)
//       Beta ( qp ) = 0.5 * qp  -  7

// The tables actually used have been "hand optimized" though (by Anthony Joch). So, the
// table values might be a little different to formula-generated values. Also, the first
// few values of both tables is set to zero to force the filter off at low qp’s

static const byte ALPHA_TABLE[52]  = {0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,4,4,5,6,  7,8,9,10,12,13,15,17,  20,22,25,28,32,36,40,45,  50,56,63,71,80,90,101,113,  127,144,162,182,203,226,255,255} ;
static const byte  BETA_TABLE[52]  = {0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,2,2,2,3,  3,3,3, 4, 4, 4, 6, 6,   7, 7, 8, 8, 9, 9,10,10,  11,11,12,12,13,13, 14, 14,   15, 15, 16, 16, 17, 17, 18, 18} ;
static const byte CLIP_TAB[52][5]  =
{
  { 0, 0, 0, 0, 0},{ 0, 0, 0, 0, 0},{ 0, 0, 0, 0, 0},{ 0, 0, 0, 0, 0},{ 0, 0, 0, 0, 0},{ 0, 0, 0, 0, 0},{ 0, 0, 0, 0, 0},{ 0, 0, 0, 0, 0},
  { 0, 0, 0, 0, 0},{ 0, 0, 0, 0, 0},{ 0, 0, 0, 0, 0},{ 0, 0, 0, 0, 0},{ 0, 0, 0, 0, 0},{ 0, 0, 0, 0, 0},{ 0, 0, 0, 0, 0},{ 0, 0, 0, 0, 0},
  { 0, 0, 0, 0, 0},{ 0, 0, 0, 1, 1},{ 0, 0, 0, 1, 1},{ 0, 0, 0, 1, 1},{ 0, 0, 0, 1, 1},{ 0, 0, 1, 1, 1},{ 0, 0, 1, 1, 1},{ 0, 1, 1, 1, 1},
  { 0, 1, 1, 1, 1},{ 0, 1, 1, 1, 1},{ 0, 1, 1, 1, 1},{ 0, 1, 1, 2, 2},{ 0, 1, 1, 2, 2},{ 0, 1, 1, 2, 2},{ 0, 1, 1, 2, 2},{ 0, 1, 2, 3, 3},
  { 0, 1, 2, 3, 3},{ 0, 2, 2, 3, 3},{ 0, 2, 2, 4, 4},{ 0, 2, 3, 4, 4},{ 0, 2, 3, 4, 4},{ 0, 3, 3, 5, 5},{ 0, 3, 4, 6, 6},{ 0, 3, 4, 6, 6},
  { 0, 4, 5, 7, 7},{ 0, 4, 5, 8, 8},{ 0, 4, 6, 9, 9},{ 0, 5, 7,10,10},{ 0, 6, 8,11,11},{ 0, 6, 8,13,13},{ 0, 7,10,14,14},{ 0, 8,11,16,16},
  { 0, 9,12,18,18},{ 0,10,13,20,20},{ 0,11,15,23,23},{ 0,13,17,25,25}
} ;

static const char chroma_edge[2][4][4] = //[dir][edge][yuv_format]
{ { {-4, 0, 0, 0},
    {-4,-4,-4, 4},
    {-4, 4, 4, 8},
    {-4,-4,-4, 12}},

  { {-4, 0,  0,  0},
    {-4,-4,  4,  4},
    {-4, 4,  8,  8},
    {-4,-4, 12, 12}}};

static const int pelnum_cr[2][4] =  {{0,8,16,16}, {0,8, 8,16}};  //[dir:0=vert, 1=hor.][yuv_format]

static void GetStrengthVer      (byte Strength[MB_BLOCK_SIZE], Macroblock *MbQ, int dir,int edge, int mvlimit);
static void GetStrengthHor      (byte Strength[MB_BLOCK_SIZE], Macroblock *MbQ, int dir,int edge, int mvlimit);
static void EdgeLoopLumaVer     (ColorPlane pl, imgpel** Img, byte Strength[MB_BLOCK_SIZE], Macroblock *MbQ, int dir, int edge, int width);
static void EdgeLoopLumaHor     (ColorPlane pl, imgpel** Img, byte Strength[MB_BLOCK_SIZE], Macroblock *MbQ, int dir, int edge, int width);
static void EdgeLoopChromaVer   (imgpel** Img, byte Strength[MB_BLOCK_SIZE], Macroblock *MbQ, int dir, int edge, int width, int uv);
static void EdgeLoopChromaHor   (imgpel** Img, byte Strength[MB_BLOCK_SIZE], Macroblock *MbQ, int dir, int edge, int width, int uv);
static void GetStrengthMBAff    (byte Strength[MB_BLOCK_SIZE], Macroblock *MbQ, int dir,int edge, int mvlimit);
static void EdgeLoopLumaMBAff   (ColorPlane pl, imgpel** Img, byte Strength[MB_BLOCK_SIZE], Macroblock *MbQ, int dir, int edge, int width);
static void EdgeLoopChromaMBAff (imgpel** Img, byte Strength[MB_BLOCK_SIZE], Macroblock *MbQ, int dir, int edge, int width, int uv);

static void DeblockMb(VideoParameters *p_Vid, imgpel **imgY, imgpel ***imgUV, int MbQAddr);


static inline int compare_mvs(const MotionVector *mv0, const MotionVector *mv1, int mvlimit)
{
  return ((iabs( mv0->mv_x - mv1->mv_x) >= 4) | (iabs( mv0->mv_y - mv1->mv_y) >= mvlimit));
}

static void  init_Deblock(VideoParameters *p_Vid)
{
  unsigned int i;
  if (p_Vid->mb_aff_frame_flag == 1) 
  {
    p_Vid->GetStrengthVer    = GetStrengthMBAff;
    p_Vid->GetStrengthHor    = GetStrengthMBAff;
    p_Vid->EdgeLoopLumaVer   = EdgeLoopLumaMBAff;
    p_Vid->EdgeLoopLumaHor   = EdgeLoopLumaMBAff;
    p_Vid->EdgeLoopChromaVer = EdgeLoopChromaMBAff;
    p_Vid->EdgeLoopChromaHor = EdgeLoopChromaMBAff;
  }
  else
  {
    p_Vid->GetStrengthVer    = GetStrengthVer;
    p_Vid->GetStrengthHor    = GetStrengthHor;
    p_Vid->EdgeLoopLumaVer   = EdgeLoopLumaVer;
    p_Vid->EdgeLoopLumaHor   = EdgeLoopLumaHor;
    p_Vid->EdgeLoopChromaVer = EdgeLoopChromaVer;
    p_Vid->EdgeLoopChromaHor = EdgeLoopChromaHor;
  }

  for (i=0; i < p_Vid->PicSizeInMbs; i++)
  {
    if (p_Vid->mb_data[i].mb_type==IPCM)
    {
      p_Vid->mb_data[i].qp = 0;
      p_Vid->mb_data[i].qpc[0] = 0;
      p_Vid->mb_data[i].qpc[1] = 0;
    }
  }
}

#define GROUP_SIZE  1

#if (JM_PARALLEL_DEBLOCK == 0)
/*!
 *****************************************************************************************
 * \brief
 *    Filter all macroblocks in order of increasing macroblock address.
 *****************************************************************************************
 */
void DeblockFrame(VideoParameters *p_Vid, imgpel **imgY, imgpel ***imgUV)
{
  unsigned int i;
  init_Deblock(p_Vid);
  for (i=0; i < p_Vid->PicSizeInMbs; i++)
  {
    DeblockMb( p_Vid, imgY, imgUV, i ) ;
  }
}
#else
static void DeblockParallel(VideoParameters *p_Vid, imgpel **imgY, imgpel ***imgUV, unsigned int column, int block, int n_last)
{
  int i, j;
  
  for (j = 0; j < GROUP_SIZE; j++)
  {
    i = block++ * (p_Vid->PicWidthInMbs - 2) + column;

    DeblockMb( p_Vid, imgY, imgUV, i ) ;
    if (block == n_last) break;
  }
}

/*!
 *****************************************************************************************
 * \brief
 *    Filter all macroblocks in a diagonal manner to enable parallelization.
 *****************************************************************************************
 */
void DeblockFrame(VideoParameters *p_Vid, imgpel **imgY, imgpel ***imgUV)
{
  int iheightMBs =(p_Vid->PicSizeInMbs/p_Vid->PicWidthInMbs);
  unsigned int i, k = p_Vid->PicWidthInMbs + 2 * (iheightMBs - 1);
  init_Deblock(p_Vid);

  for (i = 0; i < k; i++)
  {
    int nn;    
    int n_last = imin(iheightMBs, (i >> 1) + 1);
    int n_start = (i < p_Vid->PicWidthInMbs) ? 0 : ((i - p_Vid->PicWidthInMbs) >> 1) + 1;

#if defined(OPENMP)
#pragma omp parallel for
#endif
    for (nn = n_start; nn < n_last; nn += GROUP_SIZE)
      DeblockParallel(p_Vid, imgY, imgUV, i, nn, n_last);
  }
}
#endif

/*!
 *****************************************************************************************
 * \brief
 *    Deblocking filter for one macroblock.
 *****************************************************************************************
 */

static void DeblockMb(VideoParameters *p_Vid, imgpel **imgY, imgpel ***imgUV, int MbQAddr)
{
  int           edge;
  byte          Strength[16];
  short         mb_x, mb_y;

  int           filterNon8x8LumaEdgesFlag[4] = {1,1,1,1};
  int           filterLeftMbEdgeFlag;
  int           filterTopMbEdgeFlag;
  int           i, StrengthSum;
  int           edge_cr;
  Macroblock    *MbQ = &(p_Vid->mb_data[MbQAddr]) ; // current Mb
  Slice  *currSlice = MbQ->p_Slice;
  int           mvlimit = (p_Vid->structure!=FRAME) || (p_Vid->mb_aff_frame_flag && MbQ->mb_field) ? 2 : 4;
  seq_parameter_set_rbsp_t *active_sps = p_Vid->active_sps;
  p_Vid->mixedModeEdgeFlag = 0;

  // return, if filter is disabled
  if (MbQ->DFDisableIdc == 1) 
  {
    MbQ->DeblockCall = 0;
    return;
  }

  MbQ->DeblockCall = 1;
  get_mb_pos (p_Vid, MbQAddr, p_Vid->mb_size[IS_LUMA], &mb_x, &mb_y);

  if (MbQ->mb_type == I8MB)
    assert(MbQ->luma_transform_size_8x8_flag);

  filterNon8x8LumaEdgesFlag[1] =
    filterNon8x8LumaEdgesFlag[3] = !(MbQ->luma_transform_size_8x8_flag);

  filterLeftMbEdgeFlag = (mb_x != 0);
  filterTopMbEdgeFlag  = (mb_y != 0);

  if (p_Vid->mb_aff_frame_flag && mb_y == MB_BLOCK_SIZE && MbQ->mb_field)
    filterTopMbEdgeFlag = 0;

  if (MbQ->DFDisableIdc==2)
  {
    // don't filter at slice boundaries
    filterLeftMbEdgeFlag = MbQ->mbAvailA;
    // if this the bottom of a frame macroblock pair then always filter the top edge
    filterTopMbEdgeFlag  = (p_Vid->mb_aff_frame_flag && !MbQ->mb_field && (MbQAddr & 0x01)) ? 1 : MbQ->mbAvailB;
  }

  CheckAvailabilityOfNeighbors(MbQ);

  // Vertical deblocking
  for (edge = 0; edge < 4 ; ++edge )
  {
    // If cbp == 0 then deblocking for some macroblock types could be skipped                                                                  
    if (MbQ->cbp == 0)
    {
      if (filterNon8x8LumaEdgesFlag[edge] == 0 && active_sps->chroma_format_idc!=YUV444)
        continue;
      else if (edge > 0)
      {
        if ((MbQ->mb_type == 0 && currSlice->slice_type == P_SLICE) || (MbQ->mb_type == 1) || (MbQ->mb_type == 2))
          continue;
        else if ((edge & 0x01) && ((MbQ->mb_type == 3) || ((edge & 0x01) && MbQ->mb_type == 0 && currSlice->slice_type == B_SLICE && active_sps->direct_8x8_inference_flag)))
          continue;
      }
    }


    if( edge || filterLeftMbEdgeFlag )
    {
      // Strength for 4 blks in 1 stripe
      p_Vid->GetStrengthVer(Strength, MbQ, 0, edge << 2, mvlimit); // Strength for 4 blks in 1 stripe
      StrengthSum = Strength[0];

      for (i = 1; i < MB_BLOCK_SIZE && StrengthSum == 0 ; ++i)
      {
        StrengthSum += (int) Strength[i];
      }

      if( StrengthSum )                      // only if one of the 16 Strength bytes is != 0
      {
        if (filterNon8x8LumaEdgesFlag[edge])
        {
          p_Vid->EdgeLoopLumaVer( PLANE_Y, imgY, Strength, MbQ, 0, edge << 2, p_Vid->width_padded) ;
          if (p_Vid->P444_joined)
          {
            p_Vid->EdgeLoopLumaVer(PLANE_U, imgUV[0], Strength, MbQ, 0, edge << 2, p_Vid->width_padded);
            p_Vid->EdgeLoopLumaVer(PLANE_V, imgUV[1], Strength, MbQ, 0, edge << 2, p_Vid->width_padded);
          }
        }
        if(p_Vid->yuv_format==YUV420 || p_Vid->yuv_format==YUV422 )
        {
          edge_cr = chroma_edge[0][edge][p_Vid->yuv_format];
          if( (imgUV != NULL) && (edge_cr >= 0))
          {
            p_Vid->EdgeLoopChromaVer( imgUV[0], Strength, MbQ, 0, edge_cr, p_Vid->width_cr+p_Vid->pad_size_uv_x*2, 0);
            p_Vid->EdgeLoopChromaVer( imgUV[1], Strength, MbQ, 0, edge_cr, p_Vid->width_cr+p_Vid->pad_size_uv_x*2, 1);
          }
        }
      }        
    }
  }//end edge

  // horizontal deblocking  
  for( edge = 0; edge < 4 ; ++edge )
  {
    // If cbp == 0 then deblocking for some macroblock types could be skipped
    if (MbQ->cbp == 0)
    {
      if (filterNon8x8LumaEdgesFlag[edge] == 0 && active_sps->chroma_format_idc==YUV420)
        continue;
      else if (edge > 0)
      {
        if (((MbQ->mb_type == PSKIP && currSlice->slice_type == P_SLICE) || (MbQ->mb_type == P16x16) || (MbQ->mb_type == P8x16)))
          continue;
        else if ((edge & 0x01) && (( MbQ->mb_type == P16x8) || (currSlice->slice_type == B_SLICE && MbQ->mb_type == BSKIP_DIRECT && active_sps->direct_8x8_inference_flag)))
          continue;
      }
    }

    if( edge || filterTopMbEdgeFlag )
    {
      // Strength for 4 blks in 1 stripe
      p_Vid->GetStrengthHor(Strength, MbQ, 1, edge << 2, mvlimit);

      StrengthSum = Strength[0];

      for (i = 1; i < MB_BLOCK_SIZE && StrengthSum == 0 ; ++i)
      {
        StrengthSum += (int) Strength[i];
      }

      if( StrengthSum )                      // only if one of the 16 Strength bytes is != 0
      {
        if (filterNon8x8LumaEdgesFlag[edge])
        {
          p_Vid->EdgeLoopLumaHor( PLANE_Y, imgY, Strength, MbQ, 1, edge << 2, p_Vid->width_padded) ;
          if (p_Vid->P444_joined)
          {
            p_Vid->EdgeLoopLumaHor(PLANE_U, imgUV[0], Strength, MbQ, 1, edge << 2, p_Vid->width_padded);
            p_Vid->EdgeLoopLumaHor(PLANE_V, imgUV[1], Strength, MbQ, 1, edge << 2, p_Vid->width_padded);
          }
        }
        if(p_Vid->yuv_format==YUV420 || p_Vid->yuv_format==YUV422 )
        {
          edge_cr = chroma_edge[1][edge][p_Vid->yuv_format];
          if( (imgUV != NULL) && (edge_cr >= 0))
          {
            p_Vid->EdgeLoopChromaHor( imgUV[0], Strength, MbQ, 1, edge_cr, p_Vid->width_cr+p_Vid->pad_size_uv_x*2, 0);
            p_Vid->EdgeLoopChromaHor( imgUV[1], Strength, MbQ, 1, edge_cr, p_Vid->width_cr+p_Vid->pad_size_uv_x*2, 1);
          }
        }
      }

      if (!edge && !MbQ->mb_field && p_Vid->mixedModeEdgeFlag) 
      {
        // this is the extra horizontal edge between a frame macroblock pair and a field above it
        MbQ->DeblockCall = 2;
        p_Vid->GetStrengthHor(Strength, MbQ, 1, MB_BLOCK_SIZE, mvlimit); // Strength for 4 blks in 1 stripe
        //if( *((int*)Strength) )                      // only if one of the 4 Strength bytes is != 0
        {
          if (filterNon8x8LumaEdgesFlag[edge])
          {
            p_Vid->EdgeLoopLumaHor( PLANE_Y, imgY, Strength, MbQ, 1, MB_BLOCK_SIZE, p_Vid->width_padded) ;
            if (p_Vid->P444_joined)
            {
              p_Vid->EdgeLoopLumaHor(PLANE_U, imgUV[0], Strength, MbQ, 1, MB_BLOCK_SIZE, p_Vid->width_padded) ;
              p_Vid->EdgeLoopLumaHor(PLANE_V, imgUV[1], Strength, MbQ, 1, MB_BLOCK_SIZE, p_Vid->width_padded) ;
            }
          }
          if( p_Vid->yuv_format == YUV420 || p_Vid->yuv_format==YUV422 )
          {
            edge_cr = chroma_edge[1][edge][p_Vid->yuv_format];
            if( (imgUV != NULL) && (edge_cr >= 0))
            {
              p_Vid->EdgeLoopChromaHor( imgUV[0], Strength, MbQ, 1, MB_BLOCK_SIZE, p_Vid->width_cr+p_Vid->pad_size_uv_x*2, 0) ;
              p_Vid->EdgeLoopChromaHor( imgUV[1], Strength, MbQ, 1, MB_BLOCK_SIZE, p_Vid->width_cr+p_Vid->pad_size_uv_x*2, 1) ;
            }
          }
        }
        MbQ->DeblockCall = 1;
      }
    }
  }//end edge

  MbQ->DeblockCall = 0;
}

  /*!
 *********************************************************************************************
 * \brief
 *    returns a buffer of 16 Strength values for one stripe in a mb (for different Frame or Field types)
 *********************************************************************************************
 */
static void GetStrengthVer(byte Strength[MB_BLOCK_SIZE], Macroblock *MbQ, int dir, int edge, int mvlimit)
{
  VideoParameters *p_Vid = MbQ->p_Vid;
  Slice *currSlice = MbQ->p_Slice;
  int     StrValue;

  if ((currSlice->slice_type==SP_SLICE)||(currSlice->slice_type==SI_SLICE) )
  {
    // Set strength to either 3 or 4 regardless of pixel position
    StrValue = (edge == 0 && (((p_Vid->structure==FRAME)) || ((p_Vid->structure != FRAME)))) ? 4 : 3;
    memset(Strength, (byte) StrValue, MB_BLOCK_SIZE * sizeof(byte));
  }
  else
  {        
    if (!(MbQ->mb_type==I4MB||MbQ->mb_type==I8MB||MbQ->mb_type==I16MB||MbQ->mb_type==IPCM))
    {
      PixelPos pixMB;
      Macroblock *MbP;
      int xQ = edge - 1;

      getNonAffNeighbour(MbQ, xQ, 0, p_Vid->mb_size[IS_LUMA], &pixMB);
      MbP = (edge) ? MbQ : &(p_Vid->mb_data[pixMB.mb_addr]);

      if (!(MbP->mb_type==I4MB||MbP->mb_type==I8MB||MbP->mb_type==I16MB||MbP->mb_type==IPCM))
      {
        PixelPos pixP = pixMB;
        PicMotionParams **mv_info = p_Vid->enc_picture->mv_info;

        int      blkP, blkQ, idx;

        short    mb_x, mb_y;

        get_mb_block_pos_normal (MbQ->mbAddrX, &mb_x, &mb_y);
        mb_x <<= 2;
        mb_y <<= 2;

        xQ ++;

        for( idx = 0 ; idx < MB_BLOCK_SIZE ; idx += BLOCK_SIZE )
        {
          pixP.y = (short) (pixMB.y + idx);
          pixP.pos_y = (short) (pixMB.pos_y + idx);

          blkQ = (idx & 0xFFFC) + (xQ >> 2);
          blkP = (pixP.y & 0xFFFC) + (pixP.x >> 2);

          if (((MbQ->cbp_blk & i64_power2(blkQ)) != 0) || ((MbP->cbp_blk & i64_power2(blkP)) != 0))
            StrValue = 2;
          else if (edge && ((MbQ->mb_type == 1)  || (MbQ->mb_type == 2)))
            StrValue = 0; // if internal edge of certain types, we already know StrValue should be 0
          else // for everything else, if no coefs, but vector difference >= 1 set Strength=1
          {
            int blk_y  = mb_y + (blkQ >> 2);
            int blk_x  = mb_x + (blkQ  & 3);
            int blk_y2 = pixP.pos_y >> 2;
            int blk_x2 = pixP.pos_x >> 2;

            PicMotionParams *mv_info_p = &mv_info[blk_y ][blk_x ];
            PicMotionParams *mv_info_q = &mv_info[blk_y2][blk_x2];

            StorablePicturePtr ref_p0 = mv_info_p->ref_idx[LIST_0] == -1 ? NULL : mv_info_p->ref_pic[LIST_0];
            StorablePicturePtr ref_q0 = mv_info_q->ref_idx[LIST_0] == -1 ? NULL : mv_info_q->ref_pic[LIST_0];
            StorablePicturePtr ref_p1 = mv_info_p->ref_idx[LIST_1] == -1 ? NULL : mv_info_p->ref_pic[LIST_1];
            StorablePicturePtr ref_q1 = mv_info_q->ref_idx[LIST_1] == -1 ? NULL : mv_info_q->ref_pic[LIST_1];

            if ( ((ref_p0==ref_q0) && (ref_p1==ref_q1)) || ((ref_p0==ref_q1) && (ref_p1==ref_q0)))
            {
              // L0 and L1 reference pictures of p0 are different; q0 as well
              if (ref_p0 != ref_p1)
              {
                // compare MV for the same reference picture
                if (ref_p0 == ref_q0)
                {
                  StrValue = 
                    compare_mvs(&mv_info_p->mv[LIST_0], &mv_info_q->mv[LIST_0], mvlimit) |
                    compare_mvs(&mv_info_p->mv[LIST_1], &mv_info_q->mv[LIST_1], mvlimit);
                }
                else
                {
                  StrValue = 
                    compare_mvs(&mv_info_p->mv[LIST_0], &mv_info_q->mv[LIST_1], mvlimit) |
                    compare_mvs(&mv_info_p->mv[LIST_1], &mv_info_q->mv[LIST_0], mvlimit);
                }
              }
              else
              { // L0 and L1 reference pictures of p0 are the same; q0 as well
                StrValue = ((
                  compare_mvs(&mv_info_p->mv[LIST_0], &mv_info_q->mv[LIST_0], mvlimit) |
                  compare_mvs(&mv_info_p->mv[LIST_1], &mv_info_q->mv[LIST_1], mvlimit))
                  && (
                  compare_mvs(&mv_info_p->mv[LIST_0], &mv_info_q->mv[LIST_1], mvlimit) |
                  compare_mvs(&mv_info_p->mv[LIST_1], &mv_info_q->mv[LIST_0], mvlimit)
                  ));
              }
            }
            else
              StrValue = 1;
          }
          *(int*)(Strength+idx) = (StrValue<<24)|(StrValue<<16)|(StrValue<<8)|StrValue;
        }
      }
      else
      {
        // Start with Strength=3. or Strength=4 for Mb-edge
        StrValue = (edge == 0 && ((((p_Vid->structure==FRAME))) || ((p_Vid->structure != FRAME)))) ? 4 : 3;
        memset(Strength, (byte) StrValue, MB_BLOCK_SIZE * sizeof(byte));
      }      
    }
    else
    {
      // Start with Strength=3. or Strength=4 for Mb-edge
      StrValue = (edge == 0 && ((((p_Vid->structure==FRAME))) || ((p_Vid->structure != FRAME)))) ? 4 : 3;
      memset(Strength, (byte) StrValue, MB_BLOCK_SIZE * sizeof(byte));
    }      
  }
}

  /*!
 *********************************************************************************************
 * \brief
 *    returns a buffer of 16 Strength values for one stripe in a mb (for different Frame or Field types)
 *********************************************************************************************
 */
static void GetStrengthHor(byte Strength[MB_BLOCK_SIZE], Macroblock *MbQ, int dir, int edge, int mvlimit)
{
  int     StrValue;
  Slice *currSlice = MbQ->p_Slice;
  VideoParameters *p_Vid = MbQ->p_Vid;

  if ((currSlice->slice_type==SP_SLICE)||(currSlice->slice_type==SI_SLICE) )
  {
    // Set strength to either 3 or 4 regardless of pixel position
    StrValue = (edge == 0 && (((p_Vid->structure==FRAME)))) ? 4 : 3;
    memset(Strength, (byte) StrValue, MB_BLOCK_SIZE * sizeof(byte));
  }
  else
  {        
    if (!(MbQ->mb_type==I4MB||MbQ->mb_type==I8MB||MbQ->mb_type==I16MB||MbQ->mb_type==IPCM))
    {
      PixelPos pixMB;
      Macroblock *MbP;
      int yQ = (edge < 16 ? edge - 1: 0);

      getNonAffNeighbour(MbQ, 0, yQ, p_Vid->mb_size[IS_LUMA], &pixMB);
      MbP = (edge) ? MbQ : &(p_Vid->mb_data[pixMB.mb_addr]);

      if (!(MbP->mb_type==I4MB||MbP->mb_type==I8MB||MbP->mb_type==I16MB||MbP->mb_type==IPCM))
      {
        PixelPos pixP = pixMB;
        PicMotionParams **mv_info = p_Vid->enc_picture->mv_info;

        int      blkP, blkQ, idx;

        short    mb_x, mb_y;

        get_mb_block_pos_normal (MbQ->mbAddrX, &mb_x, &mb_y);
        mb_x <<= 2;
        mb_y <<= 2;
        yQ ++;

        for( idx = 0 ; idx < MB_BLOCK_SIZE ; idx += BLOCK_SIZE )
        {
          pixP.x = (short) (pixMB.x + idx);
          pixP.pos_x =  (short) (pixMB.pos_x + idx);

          blkQ = (yQ & 0xFFFC) + (idx >> 2);
          blkP = (pixP.y & 0xFFFC) + (pixP.x >> 2);

          if( ((MbQ->cbp_blk & i64_power2(blkQ)) != 0) || ((MbP->cbp_blk & i64_power2(blkP)) != 0) )
            StrValue = 2;
          else if (edge && ((MbQ->mb_type == 1)  || (MbQ->mb_type == 3)))
            StrValue = 0; // if internal edge of certain types, we already know StrValue should be 0
          else // for everything else, if no coefs, but vector difference >= 1 set Strength=1
          {
            int blk_y  = mb_y + (blkQ >> 2);
            int blk_x  = mb_x + (blkQ  & 3);
            int blk_y2 = pixP.pos_y >> 2;
            int blk_x2 = pixP.pos_x >> 2;

            PicMotionParams *mv_info_p = &mv_info[blk_y ][blk_x ];
            PicMotionParams *mv_info_q = &mv_info[blk_y2][blk_x2];

            StorablePicturePtr ref_p0 = mv_info_p->ref_idx[LIST_0] == -1 ? NULL : mv_info_p->ref_pic[LIST_0];
            StorablePicturePtr ref_q0 = mv_info_q->ref_idx[LIST_0] == -1 ? NULL : mv_info_q->ref_pic[LIST_0];
            StorablePicturePtr ref_p1 = mv_info_p->ref_idx[LIST_1] == -1 ? NULL : mv_info_p->ref_pic[LIST_1];
            StorablePicturePtr ref_q1 = mv_info_q->ref_idx[LIST_1] == -1 ? NULL : mv_info_q->ref_pic[LIST_1];

            if ( ((ref_p0==ref_q0) && (ref_p1==ref_q1)) || ((ref_p0==ref_q1) && (ref_p1==ref_q0)))
            {
              // L0 and L1 reference pictures of p0 are different; q0 as well
              if (ref_p0 != ref_p1)
              {
                // compare MV for the same reference picture
                if (ref_p0 == ref_q0)
                {
                  StrValue = 
                    compare_mvs(&mv_info_p->mv[LIST_0], &mv_info_q->mv[LIST_0], mvlimit) |
                    compare_mvs(&mv_info_p->mv[LIST_1], &mv_info_q->mv[LIST_1], mvlimit);
                }
                else
                {
                  StrValue = 
                    compare_mvs(&mv_info_p->mv[LIST_0], &mv_info_q->mv[LIST_1], mvlimit) |
                    compare_mvs(&mv_info_p->mv[LIST_1], &mv_info_q->mv[LIST_0], mvlimit);
                }
              }
              else
              { // L0 and L1 reference pictures of p0 are the same; q0 as well
                StrValue = ((
                  compare_mvs(&mv_info_p->mv[LIST_0], &mv_info_q->mv[LIST_0], mvlimit) |
                  compare_mvs(&mv_info_p->mv[LIST_1], &mv_info_q->mv[LIST_1], mvlimit))
                  && (
                  compare_mvs(&mv_info_p->mv[LIST_0], &mv_info_q->mv[LIST_1], mvlimit) |
                  compare_mvs(&mv_info_p->mv[LIST_1], &mv_info_q->mv[LIST_0], mvlimit)
                  ));
              }
            }
            else
              StrValue = 1;
          }
          *(int*)(Strength+idx) = (StrValue<<24)|(StrValue<<16)|(StrValue<<8)|StrValue;
        }
      }
      else
      {
        // Start with Strength=3. or Strength=4 for Mb-edge
        StrValue = (edge == 0 && (p_Vid->structure == FRAME)) ? 4 : 3;
        memset(Strength, (byte) StrValue, MB_BLOCK_SIZE * sizeof(byte));
      }      
    }
    else
    {
      // Start with Strength=3. or Strength=4 for Mb-edge
      StrValue = (edge == 0 && (p_Vid->structure == FRAME)) ? 4 : 3;
      memset(Strength, (byte) StrValue, MB_BLOCK_SIZE * sizeof(byte));
    }      
  }
}

/*!
 *********************************************************************************************
 * \brief
 *    returns a buffer of 16 Strength values for one stripe in a mb (for MBAFF)
 *********************************************************************************************
 */
static void GetStrengthMBAff(byte Strength[MB_BLOCK_SIZE], Macroblock *MbQ, int dir, int edge, int mvlimit)
{
  short  blkP, blkQ, idx;
  short  blk_x, blk_x2, blk_y, blk_y2 ;

  int    xQ, yQ;
  short  mb_x, mb_y;

  Macroblock *MbP;

  PixelPos pixP;
  int dir_m1 = (1 - dir);
  VideoParameters *p_Vid = MbQ->p_Vid;
  PicMotionParams **mv_info = p_Vid->enc_picture->mv_info;


  for( idx = 0; idx < 16; ++idx )
  {    
    xQ = dir ? idx : edge;
    yQ = dir ? (edge < MB_BLOCK_SIZE ? edge : 1) : idx;
    p_Vid->getNeighbour(MbQ, xQ - dir_m1, yQ - dir, p_Vid->mb_size[IS_LUMA], &pixP);
    blkQ = (short) ((yQ & 0xFFFC) + (xQ >> 2));
    blkP = (short) ((pixP.y & 0xFFFC) + (pixP.x >> 2));

    MbP = &(p_Vid->mb_data[pixP.mb_addr]);
    p_Vid->mixedModeEdgeFlag = (byte) (MbQ->mb_field != MbP->mb_field);   

    if (p_Vid->type==SP_SLICE || p_Vid->type==SI_SLICE)
    {
      Strength[idx] = (edge == 0 && (((!p_Vid->mb_aff_frame_flag && (p_Vid->structure==FRAME)) ||
        (p_Vid->mb_aff_frame_flag && !MbP->mb_field && !MbQ->mb_field)) ||
        ((p_Vid->mb_aff_frame_flag || (p_Vid->structure != FRAME)) && !dir))) ? 4 : 3;
    }
    else
    {
      // Start with Strength=3. or Strength=4 for Mb-edge
      Strength[idx] = (edge == 0 && (((!p_Vid->mb_aff_frame_flag && (p_Vid->structure==FRAME)) ||
        (p_Vid->mb_aff_frame_flag && !MbP->mb_field && !MbQ->mb_field)) ||
        ((p_Vid->mb_aff_frame_flag || (p_Vid->structure!=FRAME)) && !dir))) ? 4 : 3;

      if(  !(MbP->mb_type==I4MB || MbP->mb_type==I16MB || MbP->mb_type==I8MB || MbP->mb_type==IPCM)
        && !(MbQ->mb_type==I4MB || MbQ->mb_type==I16MB || MbQ->mb_type==I8MB || MbQ->mb_type==IPCM) )
      {        
        if( ((MbQ->cbp_blk &  ((int64)1 << blkQ )) != 0) || ((MbP->cbp_blk &  ((int64)1 << blkP)) != 0) )
          Strength[idx] = 2 ;
        else
        {
          // if no coefs, but vector difference >= 1 set Strength=1
          // if this is a mixed mode edge then one set of reference pictures will be frame and the
          // other will be field
          if (p_Vid->mixedModeEdgeFlag)
          {
            (Strength[idx] = 1);
          }
          else
          {
            p_Vid->get_mb_block_pos (MbQ->mbAddrX, &mb_x, &mb_y);
            blk_y  = (short) ((mb_y<<2) + (blkQ >> 2));
            blk_x  = (short) ((mb_x<<2) + (blkQ  & 3));
            blk_y2 = (short) (pixP.pos_y >> 2);
            blk_x2 = (short) (pixP.pos_x >> 2);
            {
              PicMotionParams *mv_info_p = &mv_info[blk_y ][blk_x ];
              PicMotionParams *mv_info_q = &mv_info[blk_y2][blk_x2];

            StorablePicturePtr ref_p0 = mv_info_p->ref_idx[LIST_0] == -1 ? NULL : mv_info_p->ref_pic[LIST_0];
            StorablePicturePtr ref_q0 = mv_info_q->ref_idx[LIST_0] == -1 ? NULL : mv_info_q->ref_pic[LIST_0];
            StorablePicturePtr ref_p1 = mv_info_p->ref_idx[LIST_1] == -1 ? NULL : mv_info_p->ref_pic[LIST_1];
            StorablePicturePtr ref_q1 = mv_info_q->ref_idx[LIST_1] == -1 ? NULL : mv_info_q->ref_pic[LIST_1];

              if ( ((ref_p0==ref_q0) && (ref_p1==ref_q1)) ||
                ((ref_p0==ref_q1) && (ref_p1==ref_q0)))
              {
                Strength[idx]=0;
                // L0 and L1 reference pictures of p0 are different; q0 as well
                if (ref_p0 != ref_p1)
                {
                  // compare MV for the same reference picture
                  if (ref_p0==ref_q0)
                  {
                    Strength[idx] =  (byte) (
                      compare_mvs(&mv_info_p->mv[LIST_0], &mv_info_q->mv[LIST_0], mvlimit) ||
                      compare_mvs(&mv_info_p->mv[LIST_1], &mv_info_q->mv[LIST_1], mvlimit));
                  }
                  else
                  {
                    Strength[idx] =  (byte) (
                      compare_mvs(&mv_info_p->mv[LIST_0], &mv_info_q->mv[LIST_1], mvlimit) ||
                      compare_mvs(&mv_info_p->mv[LIST_1], &mv_info_q->mv[LIST_0], mvlimit));
                  }
                }
                else
                { // L0 and L1 reference pictures of p0 are the same; q0 as well

                  Strength[idx] = (byte) ((
                    compare_mvs(&mv_info_p->mv[LIST_0], &mv_info_q->mv[LIST_0], mvlimit) ||
                    compare_mvs(&mv_info_p->mv[LIST_1], &mv_info_q->mv[LIST_1], mvlimit))
                    &&(
                    compare_mvs(&mv_info_p->mv[LIST_0], &mv_info_q->mv[LIST_1], mvlimit) ||
                    compare_mvs(&mv_info_p->mv[LIST_1], &mv_info_q->mv[LIST_0], mvlimit)));
                }
              }
              else
              {
                Strength[idx] = 1;
              }
            }
          }
        }
      }
    }
  }
}

/*!
 *****************************************************************************************
 * \brief
 *    Filters 16 pel block edge of Frame or Field coded MBs 
 *****************************************************************************************
 */
static void EdgeLoopLumaVer(ColorPlane pl, imgpel** Img, byte Strength[16], Macroblock *MbQ, int dir, int edge, int width)
{
  VideoParameters *p_Vid = MbQ->p_Vid;

  PixelPos pixMB1;
  getNonAffNeighbour(MbQ, edge - 1, 0, p_Vid->mb_size[IS_LUMA], &pixMB1); 

  if (pixMB1.available || (MbQ->DFDisableIdc== 0))
  {   
    int bitdepth_scale   = pl ? p_Vid->bitdepth_scale[IS_CHROMA] : p_Vid->bitdepth_scale[IS_LUMA];

    Macroblock *MbP  = &(p_Vid->mb_data[pixMB1.mb_addr]);

    // Average QP of the two blocks
    int QP = pl? ((MbP->qpc[pl-1] + MbQ->qpc[pl-1] + 1) >> 1) : (MbP->qp + MbQ->qp + 1) >> 1;

    int indexA = iClip3(0, MAX_QP, QP + MbQ->DFAlphaC0Offset);
    int indexB = iClip3(0, MAX_QP, QP + MbQ->DFBetaOffset);

    int Alpha  = ALPHA_TABLE[indexA] * bitdepth_scale;
    int Beta   = BETA_TABLE [indexB] * bitdepth_scale;

    if ((Alpha | Beta )!= 0)
    {
      const byte *ClipTab = CLIP_TAB[indexA];
      int max_imgpel_value = p_Vid->max_pel_value_comp[pl];      

      int pos_x1 = pixMB1.pos_x;
      imgpel **cur_img = &Img[pixMB1.pos_y];
      int pel;

      for( pel = 0 ; pel < MB_BLOCK_SIZE ; pel += 4 )
      {
        if(*Strength == 4 )    // INTRA strong filtering
        {
          int i;
          for( i = 0 ; i < BLOCK_SIZE ; ++i )
          {
            imgpel *SrcPtrP = *(cur_img++) + pos_x1;
            imgpel *SrcPtrQ = SrcPtrP + 1;
            imgpel  L0 = *SrcPtrP;
            imgpel  R0 = *SrcPtrQ;

            if( iabs( R0 - L0 ) < Alpha )
            {          
              imgpel  R1 = *(SrcPtrQ + 1);
              imgpel  L1 = *(SrcPtrP - 1);
              if ((iabs( R0 - R1) < Beta)  && (iabs(L0 - L1) < Beta))
              {        
                imgpel  R2 = *(SrcPtrQ + 2);
                imgpel  L2 = *(SrcPtrP - 2);
                int RL0 = L0 + R0;
                int small_gap = (iabs( R0 - L0 ) < ((Alpha >> 2) + 2));
                int aq  = ( iabs( R0 - R2) < Beta ) & small_gap;
                int ap  = ( iabs( L0 - L2) < Beta ) & small_gap;

                if (ap)
                {
                  imgpel  L3 = *(SrcPtrP - 3);
                  *(SrcPtrP--) = (imgpel)  (( R1 + ((L1 + RL0) << 1) +  L2 + 4) >> 3);
                  *(SrcPtrP--) = (imgpel)  (( L2 + L1 + RL0 + 2) >> 2);
                  *(SrcPtrP  ) = (imgpel) ((((L3 + L2) <<1) + L2 + L1 + RL0 + 4) >> 3);                
                }
                else
                {
                  *SrcPtrP = (imgpel) (((L1 << 1) + L0 + R1 + 2) >> 2);
                }

                if (aq)
                {
                  imgpel  R3 = *(SrcPtrQ + 3);
                  *(SrcPtrQ++) = (imgpel) (( L1 + ((R1 + RL0) << 1) +  R2 + 4) >> 3);
                  *(SrcPtrQ++) = (imgpel) (( R2 + R0 + L0 + R1 + 2) >> 2);
                  *(SrcPtrQ  ) = (imgpel) ((((R3 + R2) <<1) + R2 + R1 + RL0 + 4) >> 3);
                }
                else
                {
                  *SrcPtrQ = (imgpel) (((R1 << 1) + R0 + L1 + 2) >> 2);
                }
              }
            }              
          }
        }
        else if( *Strength != 0) // normal filtering
        {              
          int C0  = ClipTab[ *Strength ] * bitdepth_scale;
          int i;
          imgpel *SrcPtrP, *SrcPtrQ;
          int edge_diff;
          for( i = 0 ; i < BLOCK_SIZE ; ++i )
          {             
            SrcPtrP = *(cur_img++) + pos_x1;
            SrcPtrQ = SrcPtrP + 1;
            edge_diff = *SrcPtrQ - *SrcPtrP;

            if( iabs( edge_diff ) < Alpha )
            {          
              imgpel  *SrcPtrQ1 = SrcPtrQ + 1;
              imgpel  *SrcPtrP1 = SrcPtrP - 1;

              if ((iabs( *SrcPtrQ - *SrcPtrQ1) < Beta)  && (iabs(*SrcPtrP - *SrcPtrP1) < Beta))
              {                          
                int RL0 = (*SrcPtrP + *SrcPtrQ + 1) >> 1;
                imgpel  R2 = *(SrcPtrQ1 + 1);
                imgpel  L2 = *(SrcPtrP1 - 1);

                int aq  = (iabs(*SrcPtrQ - R2) < Beta);
                int ap  = (iabs(*SrcPtrP - L2) < Beta);

                int tc0  = (C0 + ap + aq) ;
                int dif = iClip3( -tc0, tc0, (((edge_diff) << 2) + (*SrcPtrP1 - *SrcPtrQ1) + 4) >> 3 );

                if( ap )
                  *SrcPtrP1 += iClip3( -C0,  C0, (L2 + RL0 - (*SrcPtrP1<<1)) >> 1 );

                if (dif != 0)
                {
                  *SrcPtrP = (imgpel) iClip1(max_imgpel_value, *SrcPtrP + dif);
                  *SrcPtrQ = (imgpel) iClip1(max_imgpel_value, *SrcPtrQ - dif);
                }
                if( aq  )
                  *SrcPtrQ1 += iClip3( -C0,  C0, (R2 + RL0 - (*SrcPtrQ1<<1)) >> 1 );          
              }
            }
          }
        }
        else
        {
          cur_img += 4;
        }
        Strength += 4;
      }
    }
  }  
}


/*!
 *****************************************************************************************
 * \brief
 *    Filters 16 pel block edge of Frame or Field coded MBs 
 *****************************************************************************************
 */
static void EdgeLoopLumaHor(ColorPlane pl, imgpel** Img, byte Strength[16], Macroblock *MbQ, int dir, int edge, int width)
{
  VideoParameters *p_Vid = MbQ->p_Vid;

  PixelPos pixMB1;
  getNonAffNeighbour(MbQ, 0, (edge < MB_BLOCK_SIZE ? edge - 1: 0), p_Vid->mb_size[IS_LUMA], &pixMB1); 

  if (pixMB1.available || (MbQ->DFDisableIdc== 0))
  {   
    int bitdepth_scale   = pl ? p_Vid->bitdepth_scale[IS_CHROMA] : p_Vid->bitdepth_scale[IS_LUMA];

    Macroblock *MbP  = &(p_Vid->mb_data[pixMB1.mb_addr]);

    // Average QP of the two blocks
    int QP = pl? ((MbP->qpc[pl-1] + MbQ->qpc[pl-1] + 1) >> 1) : (MbP->qp + MbQ->qp + 1) >> 1;

    int indexA = iClip3(0, MAX_QP, QP + MbQ->DFAlphaC0Offset);
    int indexB = iClip3(0, MAX_QP, QP + MbQ->DFBetaOffset);

    int Alpha  = ALPHA_TABLE[indexA] * bitdepth_scale;
    int Beta   = BETA_TABLE [indexB] * bitdepth_scale;

    if ((Alpha | Beta )!= 0)
    {
      const byte *ClipTab = CLIP_TAB[indexA];
      int max_imgpel_value = p_Vid->max_pel_value_comp[pl];
      imgpel *imgP = &Img[pixMB1.pos_y    ][pixMB1.pos_x];
      imgpel *imgQ = imgP + width;
      int pel;

      for( pel = 0 ; pel < MB_BLOCK_SIZE ; pel += 4 )
      {
        if(*Strength == 4 )    // INTRA strong filtering
        {
          int pixel;
          int inc_dim2 = width * 2;
          int inc_dim3 = width * 3;
          for( pixel = 0 ; pixel < BLOCK_SIZE ; ++pixel )
          {
            imgpel *SrcPtrP = imgP++;
            imgpel *SrcPtrQ = imgQ++;
            imgpel  L0 = *SrcPtrP;
            imgpel  R0 = *SrcPtrQ;

            if( iabs( R0 - L0 ) < Alpha )
            {          
              imgpel  L1 = *(SrcPtrP - width);
              imgpel  R1 = *(SrcPtrQ + width);

              if ((iabs( R0 - R1) < Beta)  && (iabs(L0 - L1) < Beta))
              {        
                imgpel  L2 = *(SrcPtrP - inc_dim2);
                imgpel  R2 = *(SrcPtrQ + inc_dim2);
                int RL0 = L0 + R0;
                int small_gap = (iabs( R0 - L0 ) < ((Alpha >> 2) + 2));
                int aq  = ( iabs( R0 - R2) < Beta ) & small_gap;
                int ap  = ( iabs( L0 - L2) < Beta ) & small_gap;

                if (ap)
                {
                  imgpel  L3 = *(SrcPtrP - inc_dim3);
                  *SrcPtrP   = (imgpel)  (( R1 + ((L1 + RL0) << 1) +  L2 + 4) >> 3);
                  *(SrcPtrP -= width) = (imgpel)  (( L2 + L1 + RL0 + 2) >> 2);
                  *(SrcPtrP -  width) = (imgpel) ((((L3 + L2) <<1) + L2 + L1 + RL0 + 4) >> 3);                
                }
                else
                {
                  *SrcPtrP = (imgpel) (((L1 << 1) + L0 + R1 + 2) >> 2);
                }

                if (aq)
                {
                  imgpel  R3 = *(SrcPtrQ + inc_dim3);
                  *(SrcPtrQ            ) = (imgpel) (( L1 + ((R1 + RL0) << 1) +  R2 + 4) >> 3);
                  *(SrcPtrQ += width ) = (imgpel) (( R2 + R0 + L0 + R1 + 2) >> 2);
                  *(SrcPtrQ +  width ) = (imgpel) ((((R3 + R2) <<1) + R2 + R1 + RL0 + 4) >> 3);
                }
                else
                {
                  *SrcPtrQ = (imgpel) (((R1 << 1) + R0 + L1 + 2) >> 2);
                }
              }
            }              
          }
        }
        else if( *Strength != 0) // normal filtering
        {              
          int C0  = ClipTab[ *Strength ] * bitdepth_scale;
          int i;
          imgpel *SrcPtrP, *SrcPtrQ;
          int edge_diff;
          for( i= 0 ; i < BLOCK_SIZE ; ++i )
          {             
            SrcPtrP = imgP++;
            SrcPtrQ = imgQ++;
            edge_diff = *SrcPtrQ - *SrcPtrP;

            if( iabs( edge_diff ) < Alpha )
            {          
              imgpel  *SrcPtrQ1 = SrcPtrQ + width;
              imgpel  *SrcPtrP1 = SrcPtrP - width;

              if ((iabs( *SrcPtrQ - *SrcPtrQ1) < Beta)  && (iabs(*SrcPtrP - *SrcPtrP1) < Beta))
              {                          
                int RL0 = (*SrcPtrP + *SrcPtrQ + 1) >> 1;
                imgpel  R2 = *(SrcPtrQ1 + width);
                imgpel  L2 = *(SrcPtrP1 - width);

                int aq  = (iabs(*SrcPtrQ - R2) < Beta);
                int ap  = (iabs(*SrcPtrP - L2) < Beta);

                int tc0  = (C0 + ap + aq) ;
                int dif = iClip3( -tc0, tc0, (((edge_diff) << 2) + (*SrcPtrP1 - *SrcPtrQ1) + 4) >> 3 );

                if( ap )
                  *SrcPtrP1 += iClip3( -C0,  C0, (L2 + RL0 - (*SrcPtrP1<<1)) >> 1 );

                if (dif != 0)
                {
                  *SrcPtrP = (imgpel) iClip1(max_imgpel_value, *SrcPtrP + dif);
                  *SrcPtrQ = (imgpel) iClip1(max_imgpel_value, *SrcPtrQ - dif);
                }

                if( aq  )
                  *SrcPtrQ1 += iClip3( -C0,  C0, (R2 + RL0 - (*SrcPtrQ1<<1)) >> 1 );          
              }
            }
          }
        }
        else
        {
          imgP += 4;
          imgQ += 4;
        }
        Strength += 4;
      }
    }
  }
}

/*!
 *****************************************************************************************
 * \brief
 *    Filters 16 pel block edge of Super MB Frame coded MBs
 *****************************************************************************************
 */
static void EdgeLoopLumaMBAff(ColorPlane pl, imgpel** Img, byte Strength[16], Macroblock *MbQ, int dir, int edge, int width)
{
  int      pel, ap = 0, aq = 0, Strng ;
  int      incP, incQ;
  int      C0, tc0, dif;
  imgpel   L2 = 0, L1, L0, R0, R1, R2 = 0, L3, R3;
  int      RL0;
  int      Alpha = 0, Beta = 0 ;
  const byte* ClipTab = NULL;
  int      small_gap;
  int      indexA, indexB;
  int      PelNum = MB_BLOCK_SIZE;

  int      QP;
  int      xQ, yQ;

  PixelPos pixP, pixQ;
  int      dir_m1 = (1 - dir);
  VideoParameters *p_Vid = MbQ->p_Vid;
  int      bitdepth_scale = pl? p_Vid->bitdepth_scale[IS_CHROMA] : p_Vid->bitdepth_scale[IS_LUMA];
  int      max_imgpel_value = p_Vid->max_pel_value_comp[pl];

  int AlphaC0Offset = MbQ->DFAlphaC0Offset;
  int BetaOffset = MbQ->DFBetaOffset;
  byte fieldModeFilteringFlag;

  Macroblock *MbP;
  imgpel   *SrcPtrP, *SrcPtrQ;

  for( pel = 0 ; pel < PelNum ; ++pel )
  {
    xQ = dir ? pel : edge;
    yQ = dir ? (edge < 16 ? edge : 1) : pel;
    p_Vid->getNeighbour(MbQ, xQ - (dir_m1), yQ - dir, p_Vid->mb_size[IS_LUMA], &pixP);     

    if (pixP.available || (MbQ->DFDisableIdc== 0))
    {
      if( (Strng = Strength[pel]) != 0)
      {
        p_Vid->getNeighbour(MbQ, xQ, yQ, p_Vid->mb_size[IS_LUMA], &pixQ);

        MbP = &(p_Vid->mb_data[pixP.mb_addr]);
        fieldModeFilteringFlag = (byte) (MbQ->mb_field || MbP->mb_field);

        incQ    = dir ? ((fieldModeFilteringFlag && !MbQ->mb_field) ? 2 * width : width) : 1;
        incP    = dir ? ((fieldModeFilteringFlag && !MbP->mb_field) ? 2 * width : width) : 1;
        SrcPtrQ = &(Img[pixQ.pos_y][pixQ.pos_x]);
        SrcPtrP = &(Img[pixP.pos_y][pixP.pos_x]);

        // Average QP of the two blocks
        QP = pl? ((MbP->qpc[pl-1] + MbQ->qpc[pl-1] + 1) >> 1) : (MbP->qp + MbQ->qp + 1) >> 1;

        indexA = iClip3(0, MAX_QP, QP + AlphaC0Offset);
        indexB = iClip3(0, MAX_QP, QP + BetaOffset);

        Alpha   = ALPHA_TABLE[indexA] * bitdepth_scale;
        Beta    = BETA_TABLE [indexB] * bitdepth_scale;
        ClipTab = CLIP_TAB[indexA];

        L3  = SrcPtrP[-incP*3];
        L2  = SrcPtrP[-incP*2];
        L1  = SrcPtrP[-incP];
        L0  = SrcPtrP[0] ;
        R0  = SrcPtrQ[0] ;      
        R1  = SrcPtrQ[ incQ];      
        R2  = SrcPtrQ[ incQ*2];
        R3  = SrcPtrQ[ incQ*3];

        if( iabs( R0 - L0 ) < Alpha )
        {          
          if ((iabs( R0 - R1) < Beta )   && (iabs(L0 - L1) < Beta ))
          {
            if(Strng == 4 )    // INTRA strong filtering
            {
              RL0 = L0 + R0;
              small_gap = (iabs( R0 - L0 ) < ((Alpha >> 2) + 2));
              aq  = ( iabs( R0 - R2) < Beta ) & small_gap;               
              ap  = ( iabs( L0 - L2) < Beta ) & small_gap;

              if (ap)
              {
                SrcPtrP[-incP * 2] = (imgpel) ((((L3 + L2) << 1) + L2 + L1 + RL0 + 4) >> 3);
                SrcPtrP[-incP    ] = (imgpel) (( L2 + L1 + L0 + R0 + 2) >> 2);
                SrcPtrP[    0    ] = (imgpel) (( R1 + ((L1 + RL0) << 1) +  L2 + 4) >> 3);
              }
              else
              {
                SrcPtrP[     0     ] = (imgpel) (((L1 << 1) + L0 + R1 + 2) >> 2) ;
              }

              if (aq)
              {
                SrcPtrQ[    0     ] = (imgpel) (( L1 + ((R1 + RL0) << 1) +  R2 + 4) >> 3);
                SrcPtrQ[ incQ     ] = (imgpel) (( R2 + R0 + R1 + L0 + 2) >> 2);
                SrcPtrQ[ incQ * 2 ] = (imgpel) ((((R3 + R2) << 1) + R2 + R1 + RL0 + 4) >> 3);
              }
              else
              {
                SrcPtrQ[    0     ] = (imgpel) (((R1 << 1) + R0 + L1 + 2) >> 2);
              }
            }
            else   // normal filtering
            {              
              RL0 = (L0 + R0 + 1) >> 1;
              aq  = (iabs( R0 - R2) < Beta);
              ap  = (iabs( L0 - L2) < Beta);

              C0  = ClipTab[ Strng ] * bitdepth_scale;
              tc0  = (C0 + ap + aq) ;
              dif = iClip3( -tc0, tc0, (((R0 - L0) << 2) + (L1 - R1) + 4) >> 3) ;

              if( ap )
                *(SrcPtrP - incP) += iClip3( -C0,  C0, ( L2 + RL0 - (L1 << 1)) >> 1 ) ;

              *SrcPtrP  = (imgpel) iClip1 (max_imgpel_value, L0 + dif) ;
              *SrcPtrQ  = (imgpel) iClip1 (max_imgpel_value, R0 - dif) ;

              if( aq  )
                *(SrcPtrQ + incQ) += iClip3( -C0,  C0, ( R2 + RL0 - (R1 << 1)) >> 1 ) ;
            }            
          }
        }
      }
    }
  }
}

/*!
 *****************************************************************************************
 * \brief
 *    Filters chroma block edge for Frame or Field coded pictures
 *****************************************************************************************
 */
static void EdgeLoopChromaVer(imgpel** Img, byte Strength[16], Macroblock *MbQ, int dir, int edge, int width, int uv)
{
  VideoParameters *p_Vid = MbQ->p_Vid;  

  int xQ = edge - 1;
  int yQ = 0;  
  PixelPos pixMB1;

  getNonAffNeighbour(MbQ, xQ, yQ, p_Vid->mb_size[IS_CHROMA], &pixMB1);

  if (pixMB1.available || (MbQ->DFDisableIdc == 0))
  {
    int      bitdepth_scale   = p_Vid->bitdepth_scale[IS_CHROMA];
    int      max_imgpel_value = p_Vid->max_pel_value_comp[uv + 1];

    int AlphaC0Offset = MbQ->DFAlphaC0Offset;
    int BetaOffset = MbQ->DFBetaOffset;
    Macroblock *MbP = &(p_Vid->mb_data[pixMB1.mb_addr]);

    // Average QP of the two blocks
    int QP = (MbP->qpc[uv] + MbQ->qpc[uv] + 1) >> 1;

    int indexA = iClip3(0, MAX_QP, QP + AlphaC0Offset);
    int indexB = iClip3(0, MAX_QP, QP + BetaOffset);

    int Alpha   = ALPHA_TABLE[indexA] * bitdepth_scale;
    int Beta    = BETA_TABLE [indexB] * bitdepth_scale;

    if ((Alpha | Beta) != 0)
    {
      const int PelNum = pelnum_cr[0][p_Vid->yuv_format];
      const     byte *ClipTab = CLIP_TAB[indexA];

      int pel;
      int pos_x1 = pixMB1.pos_x;
      imgpel **cur_img = &Img[pixMB1.pos_y];

      for( pel = 0 ; pel < PelNum ; ++pel )
      {
        int Strng = Strength[(PelNum == 8) ? (((pel >> 1) << 2) + (pel & 0x01)) : pel];

        if( Strng != 0)
        {
          imgpel *SrcPtrP = *cur_img + pos_x1;
          imgpel *SrcPtrQ = SrcPtrP + 1;
          int edge_diff = *SrcPtrQ - *SrcPtrP;

          if ( iabs( edge_diff ) < Alpha ) 
          {
            imgpel R1  = *(SrcPtrQ + 1);
            if ( iabs(*SrcPtrQ - R1) < Beta )  
            {
              imgpel L1  = *(SrcPtrP - 1);
              if ( iabs(*SrcPtrP - L1) < Beta )
              {
                if( Strng == 4 )    // INTRA strong filtering
                {
                  *SrcPtrP = (imgpel) ( ((L1 << 1) + *SrcPtrP + R1 + 2) >> 2 );
                  *SrcPtrQ = (imgpel) ( ((R1 << 1) + *SrcPtrQ + L1 + 2) >> 2 );
                }
                else
                {
                  int tc0  = ClipTab[ Strng ] * bitdepth_scale + 1;
                  int dif = iClip3( -tc0, tc0, ( ((edge_diff) << 2) + (L1 - R1) + 4) >> 3 );

                  if (dif != 0)
                  {
                    *SrcPtrP = (imgpel) iClip1 ( max_imgpel_value, *SrcPtrP + dif) ;
                    *SrcPtrQ = (imgpel) iClip1 ( max_imgpel_value, *SrcPtrQ - dif) ;
                  }
                }
              }
            }
          }
        }

        cur_img++;
      }
    }
  }
}


/*!
 *****************************************************************************************
 * \brief
 *    Filters chroma block edge for Frame or Field coded pictures
 *****************************************************************************************
 */
static void EdgeLoopChromaHor(imgpel** Img, byte Strength[16], Macroblock *MbQ, int dir, int edge, int width, int uv)
{
  VideoParameters *p_Vid = MbQ->p_Vid;  

  int xQ = 0;
  int yQ = (edge < 16 ? edge - 1: 0);
  PixelPos pixMB1;

  getNonAffNeighbour(MbQ, xQ, yQ, p_Vid->mb_size[IS_CHROMA], &pixMB1);

  if (pixMB1.available || (MbQ->DFDisableIdc == 0))
  {
    int      bitdepth_scale   = p_Vid->bitdepth_scale[IS_CHROMA];
    int      max_imgpel_value = p_Vid->max_pel_value_comp[uv + 1];

    int AlphaC0Offset = MbQ->DFAlphaC0Offset;
    int BetaOffset = MbQ->DFBetaOffset;
    Macroblock *MbP = &(p_Vid->mb_data[pixMB1.mb_addr]);

    // Average QP of the two blocks
    int QP = (MbP->qpc[uv] + MbQ->qpc[uv] + 1) >> 1;

    int indexA = iClip3(0, MAX_QP, QP + AlphaC0Offset);
    int indexB = iClip3(0, MAX_QP, QP + BetaOffset);

    int Alpha   = ALPHA_TABLE[indexA] * bitdepth_scale;
    int Beta    = BETA_TABLE [indexB] * bitdepth_scale;

    if ((Alpha | Beta) != 0)
    {
      const int PelNum = pelnum_cr[1][p_Vid->yuv_format];
      const     byte *ClipTab = CLIP_TAB[indexA];

      int pel;

      imgpel *imgP = &Img[pixMB1.pos_y    ][pixMB1.pos_x];
      imgpel *imgQ = imgP + width ;

      for( pel = 0 ; pel < PelNum ; ++pel )
      {
        int Strng = Strength[(PelNum == 8) ? (((pel >> 1) << 2) + (pel & 0x01)) : pel];

        if( Strng != 0)
        {
          imgpel *SrcPtrP = imgP;
          imgpel *SrcPtrQ = imgQ;
          int edge_diff = *imgQ - *imgP;

          if ( iabs( edge_diff ) < Alpha ) 
          {
            imgpel R1  = *(SrcPtrQ + width);
            if ( iabs(*SrcPtrQ - R1) < Beta )  
            {
              imgpel L1  = *(SrcPtrP - width);
              if ( iabs(*SrcPtrP - L1) < Beta )
              {
                if( Strng == 4 )    // INTRA strong filtering
                {
                  *SrcPtrP = (imgpel) ( ((L1 << 1) + *SrcPtrP + R1 + 2) >> 2 );
                  *SrcPtrQ = (imgpel) ( ((R1 << 1) + *SrcPtrQ + L1 + 2) >> 2 );
                }
                else
                {
                  int tc0  = ClipTab[ Strng ] * bitdepth_scale + 1;
                  int dif = iClip3( -tc0, tc0, ( ((edge_diff) << 2) + (L1 - R1) + 4) >> 3 );

                  if (dif != 0)
                  {
                    *SrcPtrP = (imgpel) iClip1 ( max_imgpel_value, *SrcPtrP + dif) ;
                    *SrcPtrQ = (imgpel) iClip1 ( max_imgpel_value, *SrcPtrQ - dif) ;
                  }
                }
              }
            }
          }
        }
        imgP++;
        imgQ++;
      }
    }
  }
}

/*!
*****************************************************************************************
* \brief
*    Filters chroma block edge for MBAFF types
*****************************************************************************************
 */
static void EdgeLoopChromaMBAff(imgpel** Img, byte Strength[16], Macroblock *MbQ, int dir, int edge, int width, int uv)
{
  int      pel, Strng ;
  int      incP, incQ;
  int      C0, tc0, dif;
  imgpel   L1, L0, R0, R1;
  int      Alpha = 0, Beta = 0;
  const byte* ClipTab = NULL;
  int      indexA, indexB;
  VideoParameters *p_Vid = MbQ->p_Vid;
  int      PelNum = pelnum_cr[dir][p_Vid->yuv_format];
  int      StrengthIdx;
  int      QP;
  int      xQ, yQ;
  PixelPos pixP, pixQ;
  int      dir_m1 = 1 - dir;
  int      bitdepth_scale = p_Vid->bitdepth_scale[IS_CHROMA];
  int      max_imgpel_value = p_Vid->max_pel_value_comp[uv + 1];
  
  int      AlphaC0Offset = MbQ->DFAlphaC0Offset;
  int      BetaOffset    = MbQ->DFBetaOffset;
  byte fieldModeFilteringFlag;
  Macroblock *MbP;
  imgpel   *SrcPtrP, *SrcPtrQ;

  for( pel = 0 ; pel < PelNum ; ++pel )
  {
    xQ = dir ? pel : edge;
    yQ = dir ? (edge < 16? edge : 1) : pel;
    p_Vid->getNeighbour(MbQ, xQ, yQ, p_Vid->mb_size[IS_CHROMA], &pixQ);
    p_Vid->getNeighbour(MbQ, xQ - (dir_m1), yQ - dir, p_Vid->mb_size[IS_CHROMA], &pixP);    
    MbP = &(p_Vid->mb_data[pixP.mb_addr]);    
    StrengthIdx = (PelNum == 8) ? ((MbQ->mb_field && !MbP->mb_field) ? pel << 1 :((pel >> 1) << 2) + (pel & 0x01)) : pel;

    if (pixP.available || (MbQ->DFDisableIdc == 0))
    {
      if( (Strng = Strength[StrengthIdx]) != 0)
      {
        fieldModeFilteringFlag = (byte) (MbQ->mb_field || MbP->mb_field);
        incQ = dir ? ((fieldModeFilteringFlag && !MbQ->mb_field) ? 2 * width : width) : 1;
        incP = dir ? ((fieldModeFilteringFlag && !MbP->mb_field) ? 2 * width : width) : 1;
        SrcPtrQ = &(Img[pixQ.pos_y][pixQ.pos_x]);
        SrcPtrP = &(Img[pixP.pos_y][pixP.pos_x]);

        // Average QP of the two blocks
        QP = (MbP->qpc[uv] + MbQ->qpc[uv] + 1) >> 1;

        indexA = iClip3(0, MAX_QP, QP + AlphaC0Offset);
        indexB = iClip3(0, MAX_QP, QP + BetaOffset);

        Alpha   = ALPHA_TABLE[indexA] * bitdepth_scale;
        Beta    = BETA_TABLE [indexB] * bitdepth_scale;
        ClipTab = CLIP_TAB[indexA];

        L1  = SrcPtrP[-incP];
        L0  = SrcPtrP[0] ;
        R0  = SrcPtrQ[0] ;      
        R1  = SrcPtrQ[ incQ];      

        if( iabs( R0 - L0 ) < Alpha )
        {          
          if( ((iabs( R0 - R1) - Beta < 0)  && (iabs(L0 - L1) - Beta < 0 ))  )
          {
            if( Strng == 4 )    // INTRA strong filtering
            {
              SrcPtrQ[0] = (imgpel) ( ((R1 << 1) + R0 + L1 + 2) >> 2 );
              SrcPtrP[0] = (imgpel) ( ((L1 << 1) + L0 + R1 + 2) >> 2 );
            }
            else
            {
              C0  = ClipTab[ Strng ] * bitdepth_scale;
              tc0  = (C0 + 1);
              dif = iClip3( -tc0, tc0, ( ((R0 - L0) << 2) + (L1 - R1) + 4) >> 3 );

              SrcPtrP[0] = (imgpel) iClip1 ( max_imgpel_value, L0 + dif );
              SrcPtrQ[0] = (imgpel) iClip1 ( max_imgpel_value, R0 - dif );
            }
          }
        }
      }
    }
  }
}


