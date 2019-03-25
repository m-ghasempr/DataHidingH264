
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
#include "loopfilter.h"

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


static void GetStrengthNormal (byte Strength[MB_BLOCK_SIZE], Macroblock *MbQ, int dir,int edge, int mvlimit,StorablePicture *p);
static void GetStrengthMBAff  (byte Strength[MB_BLOCK_SIZE], Macroblock *MbQ, int dir,int edge, int mvlimit,StorablePicture *p);
static void EdgeLoopLumaNormal(ColorPlane pl, imgpel** Img, byte Strength[16], Macroblock *MbQ, int dir, int edge, StorablePicture *p);
static void EdgeLoopLumaMBAff (ColorPlane pl, imgpel** Img, byte Strength[16], Macroblock *MbQ, int dir, int edge, StorablePicture *p);
static void EdgeLoopChromaNormal(imgpel** Img, byte Strength[16], Macroblock *MbQ, int dir, int edge, int uv, StorablePicture *p);
static void EdgeLoopChromaMBAff (imgpel** Img, byte Strength[16], Macroblock *MbQ, int dir, int edge, int uv, StorablePicture *p);
static void DeblockMb(VideoParameters *p_Vid, StorablePicture *p, int MbQAddr);

/*!
 *****************************************************************************************
 * \brief
 *    Filter all macroblocks in order of increasing macroblock address.
 *****************************************************************************************
 */
void DeblockPicture(VideoParameters *p_Vid, StorablePicture *p)
{
  unsigned i;

  if (p->mb_aff_frame_flag == 1) 
  {
    p_Vid->GetStrength    = GetStrengthMBAff;
    p_Vid->EdgeLoopLuma   = EdgeLoopLumaMBAff;
    p_Vid->EdgeLoopChroma = EdgeLoopChromaMBAff;
  }
  else
  {
    p_Vid->GetStrength    = GetStrengthNormal;
    p_Vid->EdgeLoopLuma   = EdgeLoopLumaNormal;
    p_Vid->EdgeLoopChroma = EdgeLoopChromaNormal;
  }

  for (i = 0; i < p->PicSizeInMbs; ++i)
  {
    DeblockMb( p_Vid, p, i ) ;
  }
}


/*!
 *****************************************************************************************
 * \brief
 *    Deblocking filter for one macroblock.
 *****************************************************************************************
 */

static void DeblockMb(VideoParameters *p_Vid, StorablePicture *p, int MbQAddr)
{
  int           EdgeCondition;
  int           dir, edge;
  byte          Strength[16];
  short         mb_x, mb_y;
  
  int           filterNon8x8LumaEdgesFlag[4] = {1,1,1,1};
  int           filterLeftMbEdgeFlag;
  int           filterTopMbEdgeFlag;
  int           fieldModeMbFlag;
  int           mvlimit = 4;
  int           i, StrengthSum;
  Macroblock    *MbQ = &(p_Vid->mb_data[MbQAddr]) ; // current Mb
  imgpel **imgY   = p->imgY;
  imgpel ***imgUV = p->imgUV;
  
  int           edge_cr;

  // return, if filter is disabled
  if (MbQ->DFDisableIdc==1) 
  {
    p_Vid->DeblockCall = 0;
    return;
  }

  p_Vid->DeblockCall = 1;
  get_mb_pos (p_Vid, MbQAddr, p_Vid->mb_size[IS_LUMA], &mb_x, &mb_y);
  
  filterLeftMbEdgeFlag = (mb_x != 0);
  filterTopMbEdgeFlag  = (mb_y != 0);
  
  if (MbQ->mb_type == I8MB)
    assert(MbQ->luma_transform_size_8x8_flag);
  
  filterNon8x8LumaEdgesFlag[1] =
  filterNon8x8LumaEdgesFlag[3] = !(MbQ->luma_transform_size_8x8_flag);
  
  if (p->mb_aff_frame_flag && mb_y == MB_BLOCK_SIZE && MbQ->mb_field)
    filterTopMbEdgeFlag = 0;
  
  fieldModeMbFlag = (p->structure!=FRAME) || (p->mb_aff_frame_flag && MbQ->mb_field);
  if (fieldModeMbFlag)
    mvlimit = 2;
  
  if (MbQ->DFDisableIdc==2)
  {
    // don't filter at slice boundaries
    filterLeftMbEdgeFlag = MbQ->mbAvailA;
    // if this the bottom of a frame macroblock pair then always filter the top edge
    filterTopMbEdgeFlag  = (p->mb_aff_frame_flag && !MbQ->mb_field && (MbQAddr & 0x01)) ? 1 : MbQ->mbAvailB;
  }

  CheckAvailabilityOfNeighbors(MbQ);

  for( dir = 0 ; dir < 2 ; ++dir )                                                      // filter first vertical edges, followed by horizontal 
  {
    EdgeCondition = (dir && filterTopMbEdgeFlag) || (!dir && filterLeftMbEdgeFlag); // can not filter beyond picture boundaries
    for( edge=0; edge<4 ; ++edge )                                                  // first 4 vertical strips of 16 pel
    {                                                                               // then  4 horizontal
      if( edge || EdgeCondition )
      {
        edge_cr = chroma_edge[dir][edge][p->chroma_format_idc];

        p_Vid->GetStrength(Strength, MbQ, dir, edge << 2, mvlimit, p); // Strength for 4 blks in 1 stripe
        StrengthSum = Strength[0];
        for (i = 1; i < MB_BLOCK_SIZE && StrengthSum == 0 ; ++i)
        {
          StrengthSum += (int) Strength[i];
        }
        
        if( StrengthSum )                      // only if one of the 16 Strength bytes is != 0
        {
          if (filterNon8x8LumaEdgesFlag[edge])
          {
            p_Vid->EdgeLoopLuma( PLANE_Y, imgY, Strength, MbQ, dir, edge << 2, p) ;
            if( p_Vid->active_sps->chroma_format_idc==YUV444 && !IS_INDEPENDENT(p_Vid) )
            {
              p_Vid->EdgeLoopLuma(PLANE_U, imgUV[0], Strength, MbQ, dir, edge << 2, p);
              p_Vid->EdgeLoopLuma(PLANE_V, imgUV[1], Strength, MbQ, dir, edge << 2, p);
            }
          }
          if (p_Vid->active_sps->chroma_format_idc==YUV420 || p_Vid->active_sps->chroma_format_idc==YUV422)
          {
            if( (imgUV != NULL) && (edge_cr >= 0))
            {
              p_Vid->EdgeLoopChroma( imgUV[0], Strength, MbQ, dir, edge_cr, 0, p);
              p_Vid->EdgeLoopChroma( imgUV[1], Strength, MbQ, dir, edge_cr, 1, p);
            }
          }
        }
                
        if (dir && !edge && !MbQ->mb_field && p_Vid->mixedModeEdgeFlag) 
        {
          // this is the extra horizontal edge between a frame macroblock pair and a field above it
          p_Vid->DeblockCall = 2;
          p_Vid->GetStrength(Strength, MbQ, dir, MB_BLOCK_SIZE, mvlimit, p); // Strength for 4 blks in 1 stripe
          //if( *((int*)Strength) )                      // only if one of the 4 Strength bytes is != 0
          {            
            if (filterNon8x8LumaEdgesFlag[edge])
            {             
              p_Vid->EdgeLoopLuma(PLANE_Y, imgY, Strength, MbQ, dir, MB_BLOCK_SIZE, p) ;
              if( p_Vid->active_sps->chroma_format_idc==YUV444 && !IS_INDEPENDENT(p_Vid) )
              {
                p_Vid->EdgeLoopLuma(PLANE_U, imgUV[0], Strength, MbQ, dir, MB_BLOCK_SIZE, p) ;
                p_Vid->EdgeLoopLuma(PLANE_V, imgUV[1], Strength, MbQ, dir, MB_BLOCK_SIZE, p) ;
              }
            }
            if (p_Vid->active_sps->chroma_format_idc==YUV420 || p_Vid->active_sps->chroma_format_idc==YUV422) 
            {
              if( (imgUV != NULL) && (edge_cr >= 0))
              {
                p_Vid->EdgeLoopChroma( imgUV[0], Strength, MbQ, dir, MB_BLOCK_SIZE, 0, p) ;
                p_Vid->EdgeLoopChroma( imgUV[1], Strength, MbQ, dir, MB_BLOCK_SIZE, 1, p) ;
              }
            }
          }
          p_Vid->DeblockCall = 1;
        }
      }
    }//end edge
  }//end loop dir
  
  p_Vid->DeblockCall = 0;
}

#define ANY_INTRA (MbP->mb_type==I4MB||MbP->mb_type==I8MB||MbP->mb_type==I16MB||MbP->mb_type==IPCM||MbQ->mb_type==I4MB||MbQ->mb_type==I8MB||MbQ->mb_type==I16MB||MbQ->mb_type==IPCM)
 
/*!
 *********************************************************************************************
 * \brief
 *    returns a buffer of 16 Strength values for one stripe in a mb (for different Frame or Field types)
 *********************************************************************************************
 */
static void GetStrengthNormal(byte Strength[MB_BLOCK_SIZE], Macroblock *MbQ, int dir, int edge, int mvlimit, StorablePicture *p)
{
  PixelPos pixP, pixMB;
  byte     StrValue;
  Macroblock *MbP;

  if ((p->slice_type==SP_SLICE)||(p->slice_type==SI_SLICE) )
  { 
    // Set strength to either 3 or 4 regardless of pixel position
    StrValue = (edge == 0 && (((p->structure==FRAME)) || ((p->structure != FRAME) && !dir))) ? 4 : 3;
    memset(&Strength[0], (byte) StrValue, MB_BLOCK_SIZE * sizeof(byte));
  }
  else
  {    
    VideoParameters *p_Vid = MbQ->p_Vid;
    int xQ = dir ? 0 : edge - 1;
    int yQ = dir ? (edge < 16 ? edge - 1: 0) : 0;

    p_Vid->getNeighbour(MbQ, xQ, yQ, p_Vid->mb_size[IS_LUMA], &pixMB);
    pixP = pixMB;
    MbP = &(p_Vid->mb_data[pixP.mb_addr]);

    if (!(MbP->mb_type==I4MB||MbP->mb_type==I8MB||MbP->mb_type==I16MB||MbP->mb_type==IPCM||MbQ->mb_type==I4MB||MbQ->mb_type==I8MB||MbQ->mb_type==I16MB||MbQ->mb_type==IPCM))
    {
      PicMotionParams *motion = &p->motion;
      int64    ref_p0,ref_p1,ref_q0,ref_q1;
      int      blkP, blkQ, idx;
      int      blk_x, blk_x2, blk_y, blk_y2 ;

      int64  **list0_refPicIdArr = motion->ref_pic_id[LIST_0];
      int64  **list1_refPicIdArr = motion->ref_pic_id[LIST_1];
      short  ***list0_mv = motion->mv[LIST_0];
      short  ***list1_mv = motion->mv[LIST_1];
      char   **list0_refIdxArr = motion->ref_idx[LIST_0];
      char   **list1_refIdxArr = motion->ref_idx[LIST_1];
      short    mb_x, mb_y;

      p_Vid->get_mb_block_pos (MbQ->mbAddrX, &mb_x, &mb_y);
      mb_x <<= 2;
      mb_y <<= 2;

      yQ += dir;
      xQ += (1 - dir);

      for( idx = 0 ; idx < MB_BLOCK_SIZE ; idx += BLOCK_SIZE )
      {
        if (dir)
        {
          xQ = idx;
          pixP.x = (short) (pixMB.x + idx);
          pixP.pos_x = (short) (pixMB.pos_x + idx);
        }
        else
        {
          yQ = idx;
          pixP.y = (short) (pixMB.y + idx);
          pixP.pos_y = (short) (pixMB.pos_y + idx);
        }

        blkQ = (yQ & 0xFFFC) + (xQ >> 2);
        blkP = (pixP.y & 0xFFFC) + (pixP.x >> 2);

        if( ((MbQ->cbp_blk[0] & ((int64)1 << blkQ )) != 0) || ((MbP->cbp_blk[0] & ((int64)1 << blkP)) != 0) )
          StrValue = 2;
        else
        {
          // if no coefs, but vector difference >= 1 set Strength=1
          // if this is a mixed mode edge then one set of reference pictures will be frame and the
          // other will be field          
          blk_y  = mb_y + (blkQ >> 2);
          blk_x  = mb_x + (blkQ  & 3);
          blk_y2 = pixP.pos_y >> 2;
          blk_x2 = pixP.pos_x >> 2;

          ref_p0 = list0_refIdxArr[blk_y ][blk_x ] < 0 ? INT64_MIN : list0_refPicIdArr[blk_y ][blk_x ];
          ref_q0 = list0_refIdxArr[blk_y2][blk_x2] < 0 ? INT64_MIN : list0_refPicIdArr[blk_y2][blk_x2];
          ref_p1 = list1_refIdxArr[blk_y ][blk_x ] < 0 ? INT64_MIN : list1_refPicIdArr[blk_y ][blk_x ];
          ref_q1 = list1_refIdxArr[blk_y2][blk_x2] < 0 ? INT64_MIN : list1_refPicIdArr[blk_y2][blk_x2];
          if ( ((ref_p0==ref_q0) && (ref_p1==ref_q1)) || ((ref_p0==ref_q1) && (ref_p1==ref_q0)))
          {
            // L0 and L1 reference pictures of p0 are different; q0 as well
            if (ref_p0 != ref_p1)
            {
              // compare MV for the same reference picture
              if (ref_p0 == ref_q0)
              {
                if (ref_p0 == INT64_MIN)
                {
                  StrValue =  (byte) (
                    (iabs( list1_mv[blk_y][blk_x][0] - list1_mv[blk_y2][blk_x2][0]) >= 4) |
                    (iabs( list1_mv[blk_y][blk_x][1] - list1_mv[blk_y2][blk_x2][1]) >= mvlimit));
                }
                else if (ref_p1 == INT64_MIN)
                {
                  StrValue =  (byte) (
                    (iabs( list0_mv[blk_y][blk_x][0] - list0_mv[blk_y2][blk_x2][0]) >= 4) |
                    (iabs( list0_mv[blk_y][blk_x][1] - list0_mv[blk_y2][blk_x2][1]) >= mvlimit));
                }
                else
                {
                  StrValue =  (byte) (
                    (iabs( list0_mv[blk_y][blk_x][0] - list0_mv[blk_y2][blk_x2][0]) >= 4) |
                    (iabs( list0_mv[blk_y][blk_x][1] - list0_mv[blk_y2][blk_x2][1]) >= mvlimit) |
                    (iabs( list1_mv[blk_y][blk_x][0] - list1_mv[blk_y2][blk_x2][0]) >= 4) |
                    (iabs( list1_mv[blk_y][blk_x][1] - list1_mv[blk_y2][blk_x2][1]) >= mvlimit));
                }
              }
              else
              {
                StrValue =  (byte) (
                  (iabs( list0_mv[blk_y][blk_x][0] - list1_mv[blk_y2][blk_x2][0]) >= 4) |
                  (iabs( list0_mv[blk_y][blk_x][1] - list1_mv[blk_y2][blk_x2][1]) >= mvlimit) |
                  (iabs( list1_mv[blk_y][blk_x][0] - list0_mv[blk_y2][blk_x2][0]) >= 4) |
                  (iabs( list1_mv[blk_y][blk_x][1] - list0_mv[blk_y2][blk_x2][1]) >= mvlimit));
              }
            }
            else
            { // L0 and L1 reference pictures of p0 are the same; q0 as well
              StrValue = (byte) (
                ((iabs( list0_mv[blk_y][blk_x][0] - list0_mv[blk_y2][blk_x2][0]) >= 4) |
                (iabs( list0_mv[blk_y][blk_x][1] - list0_mv[blk_y2][blk_x2][1]) >= mvlimit ) |
                (iabs( list1_mv[blk_y][blk_x][0] - list1_mv[blk_y2][blk_x2][0]) >= 4) |
                (iabs( list1_mv[blk_y][blk_x][1] - list1_mv[blk_y2][blk_x2][1]) >= mvlimit))
                &&
                ((iabs( list0_mv[blk_y][blk_x][0] - list1_mv[blk_y2][blk_x2][0]) >= 4) |
                (iabs( list0_mv[blk_y][blk_x][1] - list1_mv[blk_y2][blk_x2][1]) >= mvlimit) |
                (iabs( list1_mv[blk_y][blk_x][0] - list0_mv[blk_y2][blk_x2][0]) >= 4) |
                (iabs( list1_mv[blk_y][blk_x][1] - list0_mv[blk_y2][blk_x2][1]) >= mvlimit)));
            }
          }
          else
          {
            StrValue = 1;
          }
        }
        memset(&Strength[idx], (byte) StrValue, BLOCK_SIZE * sizeof(byte));
      }
    }
    else
    {
      // Start with Strength=3. or Strength=4 for Mb-edge
      StrValue = (edge == 0 && ((((p->structure==FRAME))) || ((p->structure != FRAME) && !dir))) ? 4 : 3;
      memset(&Strength[0], (byte) StrValue, MB_BLOCK_SIZE * sizeof(byte));
    }      
  }
}

/*!
 *********************************************************************************************
 * \brief
 *    returns a buffer of 16 Strength values for one stripe in a mb (for MBAFF)
 *********************************************************************************************
 */
static void GetStrengthMBAff(byte Strength[MB_BLOCK_SIZE], Macroblock *MbQ, int dir, int edge, int mvlimit, StorablePicture *p)
{
  short  blkP, blkQ, idx;
  short  blk_x, blk_x2, blk_y, blk_y2 ;
  int64  ref_p0,ref_p1,ref_q0,ref_q1;
  int    xQ, yQ;
  short  mb_x, mb_y;
  int64  **list0_refPicIdArr, **list1_refPicIdArr;
  short  ***list0_mv, ***list1_mv;
  char   **list0_refIdxArr, **list1_refIdxArr;
  Macroblock *MbP;

  PixelPos pixP;
  int dir_m1 = (1 - dir);

  PicMotionParams *motion = &p->motion;
  list0_mv = motion->mv[LIST_0];
  list1_mv = motion->mv[LIST_1];
  list0_refIdxArr = motion->ref_idx[LIST_0];
  list1_refIdxArr = motion->ref_idx[LIST_1];
  list0_refPicIdArr = motion->ref_pic_id[LIST_0];
  list1_refPicIdArr = motion->ref_pic_id[LIST_1];

  for( idx = 0; idx < 16; ++idx )
  {
    VideoParameters *p_Vid = MbQ->p_Vid;
    xQ = dir ? idx : edge;
    yQ = dir ? (edge < MB_BLOCK_SIZE ? edge : 1) : idx;
    p_Vid->getNeighbour(MbQ, xQ - dir_m1, yQ - dir, p_Vid->mb_size[IS_LUMA], &pixP);
    blkQ = (short) ((yQ & 0xFFFC) + (xQ >> 2));
    blkP = (short) ((pixP.y & 0xFFFC) + (pixP.x >> 2));

    MbP = &(p_Vid->mb_data[pixP.mb_addr]);
    p_Vid->mixedModeEdgeFlag = (byte) (MbQ->mb_field != MbP->mb_field);   

    if ((p->slice_type==SP_SLICE)||(p->slice_type==SI_SLICE) )
    {
      Strength[idx] = (edge == 0 && (((!p->mb_aff_frame_flag && (p->structure==FRAME)) ||
      (p->mb_aff_frame_flag && !MbP->mb_field && !MbQ->mb_field)) ||
      ((p->mb_aff_frame_flag || (p->structure != FRAME)) && !dir))) ? 4 : 3;
    }
    else
    {
      // Start with Strength=3. or Strength=4 for Mb-edge
      Strength[idx] = (edge == 0 && (((!p->mb_aff_frame_flag && (p->structure==FRAME)) ||
        (p->mb_aff_frame_flag && !MbP->mb_field && !MbQ->mb_field)) ||
        ((p->mb_aff_frame_flag || (p->structure!=FRAME)) && !dir))) ? 4 : 3;

      if(  !(MbP->mb_type==I4MB || MbP->mb_type==I16MB || MbP->mb_type==I8MB || MbP->mb_type==IPCM)
        && !(MbQ->mb_type==I4MB || MbQ->mb_type==I16MB || MbQ->mb_type==I8MB || MbQ->mb_type==IPCM) )
      {
        if( ((MbQ->cbp_blk[0] &  ((int64)1 << blkQ )) != 0) || ((MbP->cbp_blk[0] &  ((int64)1 << blkP)) != 0) )
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
              ref_p0 = list0_refIdxArr[blk_y ][blk_x ] < 0 ? INT64_MIN : list0_refPicIdArr[blk_y ][blk_x];
              ref_q0 = list0_refIdxArr[blk_y2][blk_x2] < 0 ? INT64_MIN : list0_refPicIdArr[blk_y2][blk_x2];
              ref_p1 = list1_refIdxArr[blk_y ][blk_x ] < 0 ? INT64_MIN : list1_refPicIdArr[blk_y ][blk_x];
              ref_q1 = list1_refIdxArr[blk_y2][blk_x2]<0 ? INT64_MIN : list1_refPicIdArr[blk_y2][blk_x2];
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
                      (iabs( list0_mv[blk_y][blk_x][0] - list0_mv[blk_y2][blk_x2][0]) >= 4) |
                      (iabs( list0_mv[blk_y][blk_x][1] - list0_mv[blk_y2][blk_x2][1]) >= mvlimit) |
                      (iabs( list1_mv[blk_y][blk_x][0] - list1_mv[blk_y2][blk_x2][0]) >= 4) |
                      (iabs( list1_mv[blk_y][blk_x][1] - list1_mv[blk_y2][blk_x2][1]) >= mvlimit));
                  }
                  else
                  {
                    Strength[idx] =  (byte) (
                      (iabs( list0_mv[blk_y][blk_x][0] - list1_mv[blk_y2][blk_x2][0]) >= 4) |
                      (iabs( list0_mv[blk_y][blk_x][1] - list1_mv[blk_y2][blk_x2][1]) >= mvlimit) |
                      (iabs( list1_mv[blk_y][blk_x][0] - list0_mv[blk_y2][blk_x2][0]) >= 4) |
                      (iabs( list1_mv[blk_y][blk_x][1] - list0_mv[blk_y2][blk_x2][1]) >= mvlimit));
                  }
                }
                else
                { // L0 and L1 reference pictures of p0 are the same; q0 as well

                  Strength[idx] = (byte) (
                    ((iabs( list0_mv[blk_y][blk_x][0] - list0_mv[blk_y2][blk_x2][0]) >= 4) |
                    (iabs( list0_mv[blk_y][blk_x][1] - list0_mv[blk_y2][blk_x2][1]) >= mvlimit ) |
                    (iabs( list1_mv[blk_y][blk_x][0] - list1_mv[blk_y2][blk_x2][0]) >= 4) |
                    (iabs( list1_mv[blk_y][blk_x][1] - list1_mv[blk_y2][blk_x2][1]) >= mvlimit))
                    &&
                    ((iabs( list0_mv[blk_y][blk_x][0] - list1_mv[blk_y2][blk_x2][0]) >= 4) |
                    (iabs( list0_mv[blk_y][blk_x][1] - list1_mv[blk_y2][blk_x2][1]) >= mvlimit) |
                    (iabs( list1_mv[blk_y][blk_x][0] - list0_mv[blk_y2][blk_x2][0]) >= 4) |
                    (iabs( list1_mv[blk_y][blk_x][1] - list0_mv[blk_y2][blk_x2][1]) >= mvlimit)));
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
static void EdgeLoopLumaNormal(ColorPlane pl, imgpel** Img, byte Strength[16], Macroblock *MbQ, 
              int dir, int edge, StorablePicture *p)
{
  VideoParameters *p_Vid = MbQ->p_Vid;
  int      xQ = dir ? 0 : edge - 1;
  int      yQ = dir ? (edge < MB_BLOCK_SIZE ? edge - 1: 0) : 0;

  PixelPos pixMB1;
  p_Vid->getNeighbour(MbQ, xQ, yQ, p_Vid->mb_size[IS_LUMA], &pixMB1); 

  if (pixMB1.available || (MbQ->DFDisableIdc== 0))
  {   
    int dirM1 = (dir - 1);
    int bitdepth_scale   = pl ? p_Vid->bitdepth_scale[IS_CHROMA] : p_Vid->bitdepth_scale[IS_LUMA];
    
    PixelPos   pixP = pixMB1;
    Macroblock *MbP  = &(p_Vid->mb_data[pixP.mb_addr]);

    // Average QP of the two blocks
    int QP = pl? ((MbP->qpc[pl-1] + MbQ->qpc[pl-1] + 1) >> 1) : (MbP->qp + MbQ->qp + 1) >> 1;

    int indexA = iClip3(0, MAX_QP, QP + MbQ->DFAlphaC0Offset);
    int indexB = iClip3(0, MAX_QP, QP + MbQ->DFBetaOffset);

    int Alpha  = ALPHA_TABLE[indexA] * bitdepth_scale;
    int Beta   = BETA_TABLE [indexB] * bitdepth_scale;

    PixelPos pixMB2;
    p_Vid->getNeighbour(MbQ, (xQ -= dirM1), (yQ += dir), p_Vid->mb_size[IS_LUMA], &pixMB2);

    if (Alpha != 0 && Beta !=0)
    {
      const byte *ClipTab = CLIP_TAB   [indexA];
      PixelPos pixQ = pixMB2;
      int max_imgpel_value = p_Vid->max_pel_value_comp[pl];
      int width = p->size_x;
      int inc_dim = dir ? width : 1;
      int pel;

      for( pel = 0 ; pel < MB_BLOCK_SIZE ; ++pel )
      {
        if( *Strength != 0)
        {
          imgpel *SrcPtrP = &Img[pixP.pos_y][pixP.pos_x];
          imgpel  L0 = *SrcPtrP;
          
          imgpel *SrcPtrQ = &(Img[pixQ.pos_y][pixQ.pos_x]);
          imgpel  R0 = *SrcPtrQ;

          if( iabs( R0 - L0 ) < Alpha )
          {          
            imgpel  R1 = *(SrcPtrQ + inc_dim);
            imgpel  L1 = *(SrcPtrP - inc_dim);
            if ((iabs( R0 - R1) < Beta)  && (iabs(L0 - L1) < Beta))
            {        
              int inc_dim2 = inc_dim * 2;
              imgpel  R2 = *(SrcPtrQ + inc_dim2);
              imgpel  L2 = *(SrcPtrP - inc_dim2);
              
              if(*Strength == 4 )    // INTRA strong filtering
              {
                int inc_dim3 = inc_dim * 3;
                int RL0 = L0 + R0;
                int small_gap = (iabs( R0 - L0 ) < ((Alpha >> 2) + 2));
                int aq  = ( iabs( R0 - R2) < Beta ) & small_gap;
                int ap  = ( iabs( L0 - L2) < Beta ) & small_gap;

                if (ap)
                {
                  imgpel  L3 = *(SrcPtrP - inc_dim3);
                  *SrcPtrP              = (imgpel)  (( R1 + ((L1 + RL0) << 1) +  L2 + 4) >> 3);
                  *(SrcPtrP -= inc_dim) = (imgpel)  (( L2 + L1 + RL0 + 2) >> 2);
                  *(SrcPtrP -  inc_dim) = (imgpel) ((((L3 + L2) <<1) + L2 + L1 + RL0 + 4) >> 3);                
                }
                else
                {
                  *SrcPtrP = (imgpel) (((L1 << 1) + L0 + R1 + 2) >> 2) ;                
                }

                if (aq)
                {
                   imgpel  R3 = *(SrcPtrQ + inc_dim3);
                  *(SrcPtrQ            ) = (imgpel) (( L1 + ((R1 + RL0) << 1) +  R2 + 4) >> 3);
                  *(SrcPtrQ += inc_dim ) = (imgpel) (( R2 + R0 + L0 + R1 + 2) >> 2);
                  *(SrcPtrQ +  inc_dim ) = (imgpel) ((((R3 + R2) <<1) + R2 + R1 + RL0 + 4) >> 3);
                }
                else
                {
                  *SrcPtrQ = (imgpel) (((R1 << 1) + R0 + L1 + 2) >> 2);
                }
              }
              else   // normal filtering
              {              
                int RL0 = (L0 + R0 + 1) >> 1;
                int aq  = (iabs(R0 - R2) < Beta);
                int ap  = (iabs(L0 - L2) < Beta);

                int C0  = ClipTab[ *Strength ] * bitdepth_scale;
                int tc0  = (C0 + ap + aq) ;
                int dif = iClip3( -tc0, tc0, (((R0 - L0) << 2) + (L1 - R1) + 4) >> 3 );

                if( ap )
                  *(SrcPtrP - inc_dim) += iClip3( -C0,  C0, (L2 + RL0 - (L1<<1)) >> 1 );
                *SrcPtrP = (imgpel) iClip1(max_imgpel_value, L0 + dif);
                 
                *SrcPtrQ = (imgpel) iClip1(max_imgpel_value, R0 - dif);
                if( aq  )
                  *(SrcPtrQ + inc_dim) += iClip3( -C0,  C0, (R2 + RL0 - (R1<<1)) >> 1 );
              }            
            }
          }
        }
        Strength++;
        pixP.pos_x = (short) (pixP.pos_x + dir);
        pixP.pos_y = (short) (pixP.pos_y - dirM1);
        pixQ.pos_x = (short) (pixQ.pos_x + dir);      
        pixQ.pos_y = (short) (pixQ.pos_y - dirM1);
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
static void EdgeLoopLumaMBAff(ColorPlane pl, imgpel** Img, byte Strength[16], Macroblock *MbQ, 
              int dir, int edge, StorablePicture *p)
{
  int      width = p->size_x;
  int      pel, ap = 0, aq = 0, Strng ;
  int      incP, incQ;
  int      C0, tc0, dif;
  imgpel   L2 = 0, L1, L0, R0, R1, R2 = 0, L3, R3;
  int      RL0;
  int      Alpha = 0, Beta = 0 ;
  const byte* ClipTab = NULL;
  int      small_gap;
  int      indexA, indexB;
  int      PelNum = pl? pelnum_cr[dir][p->chroma_format_idc] : MB_BLOCK_SIZE;

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
static void EdgeLoopChromaNormal(imgpel** Img, byte Strength[16], Macroblock *MbQ, int dir, int edge, int uv, StorablePicture *p)
{ 
  VideoParameters *p_Vid = MbQ->p_Vid;  

  int xQ = dir ? 0 : edge - 1;
  int yQ = dir ? (edge < 16 ? edge - 1: 0) : 0;  
  PixelPos pixMB1;
  
  p_Vid->getNeighbour(MbQ, xQ, yQ, p_Vid->mb_size[IS_CHROMA], &pixMB1);

  if (pixMB1.available || (MbQ->DFDisableIdc == 0))
  {
    int      bitdepth_scale   = p_Vid->bitdepth_scale[IS_CHROMA];
    int      max_imgpel_value = p_Vid->max_pel_value_comp[uv + 1];

    int AlphaC0Offset = MbQ->DFAlphaC0Offset;
    int BetaOffset = MbQ->DFBetaOffset;
    int width = p->size_x_cr;
    int dirM1 = dir - 1;
    PixelPos pixP = pixMB1;
    Macroblock *MbP = &(p_Vid->mb_data[pixP.mb_addr]);

    // Average QP of the two blocks
    int QP = (MbP->qpc[uv] + MbQ->qpc[uv] + 1) >> 1;

    int indexA = iClip3(0, MAX_QP, QP + AlphaC0Offset);
    int indexB = iClip3(0, MAX_QP, QP + BetaOffset);

    int Alpha   = ALPHA_TABLE[indexA] * bitdepth_scale;
    int Beta    = BETA_TABLE [indexB] * bitdepth_scale;

    if (Alpha !=0 && Beta != 0)
    {
      const int PelNum = pelnum_cr[dir][p->chroma_format_idc];
      const     byte *ClipTab = CLIP_TAB[indexA];
      int       inc_dim = dir ? width : 1;
      int pel;
      PixelPos pixQ, pixMB2;    

      p_Vid->getNeighbour(MbQ, (xQ -= dirM1), (yQ += dir), p_Vid->mb_size[IS_CHROMA], &pixMB2);
      pixQ = pixMB2;

      for( pel = 0 ; pel < PelNum ; ++pel )
      {
        int Strng = Strength[(PelNum == 8) ? (((pel >> 1) << 2) + (pel & 0x01)) : pel];

        if( Strng != 0)
        {
          imgpel *SrcPtrP = &(Img[pixP.pos_y][pixP.pos_x]);
          imgpel  L0  = *SrcPtrP;
          imgpel *SrcPtrQ = &(Img[pixQ.pos_y][pixQ.pos_x]);
          imgpel  R0  = *SrcPtrQ;

          if ( iabs( R0 - L0 ) < Alpha ) 
          {
            imgpel R1  = *(SrcPtrQ + inc_dim);
            if ( iabs(R0 - R1) < Beta )  
            {
              imgpel L1  = *(SrcPtrP - inc_dim);
              if ( iabs(L0 - L1) < Beta )
              {
                if( Strng == 4 )    // INTRA strong filtering
                {
                  *SrcPtrP = (imgpel) ( ((L1 << 1) + L0 + R1 + 2) >> 2 );
                  *SrcPtrQ = (imgpel) ( ((R1 << 1) + R0 + L1 + 2) >> 2 );
                }
                else
                {
                  int tc0  = ClipTab[ Strng ] * bitdepth_scale + 1;
                  int dif = iClip3( -tc0, tc0, ( ((R0 - L0) << 2) + (L1 - R1) + 4) >> 3 );

                  *SrcPtrP = (imgpel) iClip1 ( max_imgpel_value, L0 + dif) ;
                  *SrcPtrQ = (imgpel) iClip1 ( max_imgpel_value, R0 - dif) ;
                }
              }
            }
          }
        }
      pixP.pos_x = (short) (pixP.pos_x + dir);
      pixP.pos_y = (short) (pixP.pos_y - dirM1);
      pixQ.pos_x = (short) (pixQ.pos_x + dir);
      pixQ.pos_y = (short) (pixQ.pos_y - dirM1);
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
static void EdgeLoopChromaMBAff(imgpel** Img, byte Strength[16], Macroblock *MbQ, int dir, int edge, int uv, StorablePicture *p)
{
  int      pel, Strng ;
  int      incP, incQ;
  int      C0, tc0, dif;
  imgpel   L1, L0, R0, R1;
  int      Alpha = 0, Beta = 0;
  const byte* ClipTab = NULL;
  int      indexA, indexB;
  VideoParameters *p_Vid = MbQ->p_Vid;
  int      PelNum = pelnum_cr[dir][p->chroma_format_idc];
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
  int      width = p->size_x_cr;

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
          //if( ((iabs( R0 - R1) - Beta )  & (iabs(L0 - L1) - Beta )) < 0  )
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


