
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


void GetStrengthNormal (byte Strength[16], Macroblock *MbQ, int dir,int edge, int mvlimit,StorablePicture *p);
void GetStrengthMBAff  (byte Strength[16], Macroblock *MbQ, int dir,int edge, int mvlimit,StorablePicture *p);
void EdgeLoopLumaNormal(ColorPlane pl, imgpel** Img, byte Strength[16],Macroblock *MbQ, int dir, int edge, StorablePicture *p);
void EdgeLoopLumaMBAff (ColorPlane pl, imgpel** Img, byte Strength[16],Macroblock *MbQ, int dir, int edge, StorablePicture *p);
void EdgeLoopChromaNormal(imgpel** Img, byte Strength[16],Macroblock *MbQ, int dir, int edge, int uv, StorablePicture *p);
void EdgeLoopChromaMBAff(imgpel** Img, byte Strength[16],Macroblock *MbQ, int dir, int edge, int uv, StorablePicture *p);
void DeblockMb(ImageParameters *p_Img, StorablePicture *p, int MbQAddr);
int compute_deblock_strength(char **list0_refIdxArr, char **list1_refIdxArr, int64 **list0_refPicIdArr, int64 **list1_refPicIdArr, 
                             short  ***list0_mv, short  ***list1_mv, int blk_y, int blk_x, int blk_y2, int blk_x2, int mvlimit);

/*!
 *****************************************************************************************
 * \brief
 *    Filter all macroblocks in order of increasing macroblock address.
 *****************************************************************************************
 */
void DeblockPicture(ImageParameters *p_Img, StorablePicture *p)
{
  unsigned i;

  if (p->MbaffFrameFlag == 1) 
  {
    p_Img->GetStrength    = GetStrengthMBAff;
    p_Img->EdgeLoopLuma   = EdgeLoopLumaMBAff;
    p_Img->EdgeLoopChroma = EdgeLoopChromaMBAff;
  }
  else
  {
    p_Img->GetStrength    = GetStrengthNormal;
    p_Img->EdgeLoopLuma   = EdgeLoopLumaNormal;
    p_Img->EdgeLoopChroma = EdgeLoopChromaNormal;
  }

  for (i = 0; i < p->PicSizeInMbs; ++i)
  {
    DeblockMb( p_Img, p, i ) ;
  }
}


/*!
 *****************************************************************************************
 * \brief
 *    Deblocking filter for one macroblock.
 *****************************************************************************************
 */

void DeblockMb(ImageParameters *p_Img, StorablePicture *p, int MbQAddr)
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
  Macroblock    *MbQ = &(p_Img->mb_data[MbQAddr]) ; // current Mb
  imgpel **imgY   = p->imgY;
  imgpel ***imgUV = p->imgUV;
  
  int           edge_cr;

  // return, if filter is disabled
  if (MbQ->DFDisableIdc==1) 
  {
    p_Img->DeblockCall = 0;
    return;
  }

  p_Img->DeblockCall = 1;
  get_mb_pos (p_Img, MbQAddr, p_Img->mb_size[IS_LUMA], &mb_x, &mb_y);
  
  filterLeftMbEdgeFlag = (mb_x != 0);
  filterTopMbEdgeFlag  = (mb_y != 0);
  
  if (MbQ->mb_type == I8MB)
    assert(MbQ->luma_transform_size_8x8_flag);
  
  filterNon8x8LumaEdgesFlag[1] =
  filterNon8x8LumaEdgesFlag[3] = !(MbQ->luma_transform_size_8x8_flag);
  
  if (p->MbaffFrameFlag && mb_y == MB_BLOCK_SIZE && MbQ->mb_field)
    filterTopMbEdgeFlag = 0;
  
  fieldModeMbFlag = (p->structure!=FRAME) || (p->MbaffFrameFlag && MbQ->mb_field);
  if (fieldModeMbFlag)
    mvlimit = 2;
  
  if (MbQ->DFDisableIdc==2)
  {
    // don't filter at slice boundaries
    filterLeftMbEdgeFlag = MbQ->mbAvailA;
    // if this the bottom of a frame macroblock pair then always filter the top edge
    filterTopMbEdgeFlag  = (p->MbaffFrameFlag && !MbQ->mb_field && (MbQAddr & 0x01)) ? 1 : MbQ->mbAvailB;
  }

  CheckAvailabilityOfNeighbors(MbQ);

  for( dir = 0 ; dir < 2 ; ++dir )                                                      // filter first vertical edges, followed by horizontal 
  {
    EdgeCondition = (dir && filterTopMbEdgeFlag) || (!dir && filterLeftMbEdgeFlag); // can not filter beyond picture boundaries
    for( edge=0; edge<4 ; ++edge )                                            // first 4 vertical strips of 16 pel
    {                                                                               // then  4 horizontal
      if( edge || EdgeCondition )
      {
        edge_cr = chroma_edge[dir][edge][p->chroma_format_idc];

        p_Img->GetStrength(Strength, MbQ, dir, edge << 2, mvlimit, p); // Strength for 4 blks in 1 stripe
        StrengthSum = Strength[0];
        for (i = 1; i < MB_BLOCK_SIZE && StrengthSum == 0 ; ++i)
        {
          StrengthSum += (int) Strength[i];
        }
        
        if( StrengthSum )                      // only if one of the 16 Strength bytes is != 0
        {
          if (filterNon8x8LumaEdgesFlag[edge])
          {
            p_Img->EdgeLoopLuma( PLANE_Y, imgY, Strength, MbQ, dir, edge << 2, p) ;
            if( p_Img->active_sps->chroma_format_idc==YUV444 && !IS_INDEPENDENT(p_Img) )
            {
              p_Img->EdgeLoopLuma(PLANE_U, imgUV[0], Strength, MbQ, dir, edge << 2, p);
              p_Img->EdgeLoopLuma(PLANE_V, imgUV[1], Strength, MbQ, dir, edge << 2, p);
            }
          }
          if (p_Img->active_sps->chroma_format_idc==YUV420 || p_Img->active_sps->chroma_format_idc==YUV422)
          {
            if( (imgUV != NULL) && (edge_cr >= 0))
            {
              p_Img->EdgeLoopChroma( imgUV[0], Strength, MbQ, dir, edge_cr, 0, p);
              p_Img->EdgeLoopChroma( imgUV[1], Strength, MbQ, dir, edge_cr, 1, p);
            }
          }
        }
                
        if (dir && !edge && !MbQ->mb_field && p_Img->mixedModeEdgeFlag) 
        {
          // this is the extra horizontal edge between a frame macroblock pair and a field above it
          p_Img->DeblockCall = 2;
          p_Img->GetStrength(Strength, MbQ, dir, MB_BLOCK_SIZE, mvlimit, p); // Strength for 4 blks in 1 stripe
          //if( *((int*)Strength) )                      // only if one of the 4 Strength bytes is != 0
          {            
            if (filterNon8x8LumaEdgesFlag[edge])
            {             
              p_Img->EdgeLoopLuma(PLANE_Y, imgY, Strength, MbQ, dir, MB_BLOCK_SIZE, p) ;
              if( p_Img->active_sps->chroma_format_idc==YUV444 && !IS_INDEPENDENT(p_Img) )
              {
                p_Img->EdgeLoopLuma(PLANE_U, imgUV[0], Strength, MbQ, dir, MB_BLOCK_SIZE, p) ;
                p_Img->EdgeLoopLuma(PLANE_V, imgUV[1], Strength, MbQ, dir, MB_BLOCK_SIZE, p) ;
              }
            }
            if (p_Img->active_sps->chroma_format_idc==YUV420 || p_Img->active_sps->chroma_format_idc==YUV422) 
            {
              if( (imgUV != NULL) && (edge_cr >= 0))
              {
                p_Img->EdgeLoopChroma( imgUV[0], Strength, MbQ, dir, MB_BLOCK_SIZE, 0, p) ;
                p_Img->EdgeLoopChroma( imgUV[1], Strength, MbQ, dir, MB_BLOCK_SIZE, 1, p) ;
              }
            }
          }
          p_Img->DeblockCall = 1;
        }
      }
    }//end edge
  }//end loop dir
  
  p_Img->DeblockCall = 0;
}

  /*!
 *********************************************************************************************
 * \brief
 *    returns a buffer of 16 Strength values for one stripe in a mb (for different Frame or Field types)
 *********************************************************************************************
 */

#define ANY_INTRA (MbP->mb_type==I4MB||MbP->mb_type==I8MB||MbP->mb_type==I16MB||MbP->mb_type==IPCM||MbQ->mb_type==I4MB||MbQ->mb_type==I8MB||MbQ->mb_type==I16MB||MbQ->mb_type==IPCM)

void GetStrengthNormal(byte Strength[MB_BLOCK_SIZE], Macroblock *MbQ, int dir, int edge, int mvlimit, StorablePicture *p)
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
    ImageParameters *p_Img = MbQ->p_Img;
    int xQ = dir ? 0 : edge - 1;
    int yQ = dir ? (edge < 16 ? edge - 1: 0) : 0;

    p_Img->getNeighbour(MbQ, xQ, yQ, p_Img->mb_size[IS_LUMA], &pixMB);
    pixP = pixMB;
    MbP = &(p_Img->mb_data[pixP.mb_addr]);

    if (!ANY_INTRA)
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

      p_Img->get_mb_block_pos (MbQ->mbAddrX, &mb_x, &mb_y);
      mb_x <<= 2;
      mb_y <<= 2;

      yQ += dir;
      xQ += (1 - dir);

      for( idx = 0 ; idx < MB_BLOCK_SIZE ; idx += BLOCK_SIZE )
      {
        if (dir)
        {
          xQ = idx;
          pixP.x = pixMB.x + idx;
          pixP.pos_x = pixMB.pos_x + idx;
        }
        else
        {
          yQ = idx;
          pixP.y = pixMB.y + idx;
          pixP.pos_y = pixMB.pos_y + idx;
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
void GetStrengthMBAff(byte Strength[16], Macroblock *MbQ, int dir, int edge, int mvlimit, StorablePicture *p)
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
    ImageParameters *p_Img = MbQ->p_Img;
    xQ = dir ? idx : edge;
    yQ = dir ? (edge < MB_BLOCK_SIZE ? edge : 1) : idx;
    p_Img->getNeighbour(MbQ, xQ - dir_m1, yQ - dir, p_Img->mb_size[IS_LUMA], &pixP);
    blkQ = (yQ & 0xFFFC) + (xQ >> 2);
    blkP = (pixP.y & 0xFFFC) + (pixP.x >> 2);

    MbP = &(p_Img->mb_data[pixP.mb_addr]);
    p_Img->mixedModeEdgeFlag = (byte) (MbQ->mb_field != MbP->mb_field);   

    if ((p->slice_type==SP_SLICE)||(p->slice_type==SI_SLICE) )
    {
      Strength[idx] = (edge == 0 && (((!p->MbaffFrameFlag && (p->structure==FRAME)) ||
      (p->MbaffFrameFlag && !MbP->mb_field && !MbQ->mb_field)) ||
      ((p->MbaffFrameFlag || (p->structure != FRAME)) && !dir))) ? 4 : 3;
    }
    else
    {
      // Start with Strength=3. or Strength=4 for Mb-edge
      Strength[idx] = (edge == 0 && (((!p->MbaffFrameFlag && (p->structure==FRAME)) ||
        (p->MbaffFrameFlag && !MbP->mb_field && !MbQ->mb_field)) ||
        ((p->MbaffFrameFlag || (p->structure!=FRAME)) && !dir))) ? 4 : 3;

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
          if (p_Img->mixedModeEdgeFlag)
          {
            (Strength[idx] = 1);
          }
          else
          {
            p_Img->get_mb_block_pos (MbQ->mbAddrX, &mb_x, &mb_y);
            blk_y  = (mb_y<<2) + (blkQ >> 2) ;
            blk_x  = (mb_x<<2) + (blkQ  & 3) ;
            blk_y2 = pixP.pos_y >> 2;
            blk_x2 = pixP.pos_x >> 2;
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
void EdgeLoopLumaNormal(ColorPlane pl, imgpel** Img, byte Strength[16], Macroblock *MbQ, 
              int dir, int edge, StorablePicture *p)
{
  imgpel   L2, L1, L0, R0, R1, R2, L3, R3;
  PixelPos pixP, pixQ, pixMB1, pixMB2;
  int      C0, tc0, dif, RL0;
  int      pel, ap, aq, Strng;
  int      Alpha, Beta, small_gap;
  int      indexA, indexB;
  int      QP;
  const    byte* ClipTab;
  int      inc_dim, inc_dim3;
  int      xQ = dir ? 0 : edge - 1;
  int      yQ = dir ? (edge < MB_BLOCK_SIZE ? edge - 1: 0) : 0;

  Macroblock *MbP;
  imgpel   *SrcPtrP, *SrcPtrQ;

  ImageParameters *p_Img = MbQ->p_Img;

  p_Img->getNeighbour(MbQ, xQ, yQ, p_Img->mb_size[IS_LUMA], &pixMB1); 

  if (pixMB1.available || (MbQ->DFDisableIdc== 0))
  {
    int AlphaC0Offset = MbQ->DFAlphaC0Offset;
    int BetaOffset = MbQ->DFBetaOffset;
    int  dirM1 = (dir - 1);
    int  bitdepth_scale   = pl ? p_Img->bitdepth_scale[IS_CHROMA] : p_Img->bitdepth_scale[IS_LUMA];
    int  max_imgpel_value = p_Img->max_imgpel_value_comp[pl];
    int  width = p->size_x;
    pixP = pixMB1;
    MbP  = &(p_Img->mb_data[pixP.mb_addr]);
    inc_dim = dir ? width : 1;
    inc_dim3 = inc_dim * 3;

    yQ += dir;
    xQ -= dirM1;
    p_Img->getNeighbour(MbQ, xQ, yQ, p_Img->mb_size[IS_LUMA], &pixMB2);
    pixQ = pixMB2;

    // Average QP of the two blocks
    QP = pl? ((MbP->qpc[pl-1] + MbQ->qpc[pl-1] + 1) >> 1) : (MbP->qp + MbQ->qp + 1) >> 1;

    indexA = iClip3(0, MAX_QP, QP + AlphaC0Offset);
    indexB = iClip3(0, MAX_QP, QP + BetaOffset);

    Alpha   = ALPHA_TABLE[indexA] * bitdepth_scale;
    Beta    = BETA_TABLE [indexB] * bitdepth_scale;
    ClipTab = CLIP_TAB[indexA];

    for( pel = 0 ; pel < MB_BLOCK_SIZE ; ++pel )
    {
      if( (Strng = *(Strength++)) != 0)
      {
        SrcPtrQ = &(Img[pixQ.pos_y][pixQ.pos_x]);
        SrcPtrP = &(Img[pixP.pos_y][pixP.pos_x]);

        L3 = *(SrcPtrP -= inc_dim3);
        L2 = *(SrcPtrP += inc_dim);
        L1 = *(SrcPtrP += inc_dim);
        L0 = *(SrcPtrP += inc_dim);
        R0 = *SrcPtrQ;
        R1 = *(SrcPtrQ += inc_dim);
        R2 = *(SrcPtrQ += inc_dim);
        R3 = *(SrcPtrQ += inc_dim);        

        if( iabs( R0 - L0 ) < Alpha )
        {          
          if ((iabs( R0 - R1) < Beta)  && (iabs(L0 - L1) < Beta))
          {            
            if(Strng == 4 )    // INTRA strong filtering
            {
              RL0 = L0 + R0;
              small_gap = (iabs( R0 - L0 ) < ((Alpha >> 2) + 2));
              aq  = ( iabs( R0 - R2) < Beta ) & small_gap;
              ap  = ( iabs( L0 - L2) < Beta ) & small_gap;

              if (ap)
              {
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
                *(SrcPtrQ -= inc_dim3) = (imgpel) (( L1 + ((R1 + RL0) << 1) +  R2 + 4) >> 3);
                *(SrcPtrQ += inc_dim ) = (imgpel) (( R2 + R0 + L0 + R1 + 2) >> 2);
                *(SrcPtrQ +  inc_dim ) = (imgpel) ((((R3 + R2) <<1) + R2 + R1 + RL0 + 4) >> 3);
              }
              else
              {
                *(SrcPtrQ -  inc_dim3) = (imgpel) (((R1 << 1) + R0 + L1 + 2) >> 2);
              }
            }
            else   // normal filtering
            {              
              RL0 = (L0 + R0 + 1) >> 1;
              aq  = (iabs(R0 - R2) < Beta);
              ap  = (iabs(L0 - L2) < Beta);

              C0  = ClipTab[ Strng ] * bitdepth_scale;
              tc0  = (C0 + ap + aq) ;
              dif = iClip3( -tc0, tc0, (((R0 - L0) << 2) + (L1 - R1) + 4) >> 3 );

              if( ap )
                *(SrcPtrP - inc_dim) += iClip3( -C0,  C0, (L2 + RL0 - (L1<<1)) >> 1 );

              *SrcPtrP               = (imgpel) iClip1(max_imgpel_value, L0 + dif);
              *(SrcPtrQ -= inc_dim3) = (imgpel) iClip1(max_imgpel_value, R0 - dif);

              if( aq  )
                *(SrcPtrQ + inc_dim) += iClip3( -C0,  C0, (R2 + RL0 - (R1<<1)) >> 1 );
            }            
          }
        }
      }
      pixP.pos_x += dir;
      pixQ.pos_x += dir;
      pixP.pos_y -= dirM1;
      pixQ.pos_y -= dirM1;
    }
  }
}

/*!
 *****************************************************************************************
 * \brief
 *    Filters 16 pel block edge of Super MB Frame coded MBs
 *****************************************************************************************
 */
void EdgeLoopLumaMBAff(ColorPlane pl, imgpel** Img, byte Strength[16], Macroblock *MbQ, 
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
  ImageParameters *p_Img = MbQ->p_Img;
  int      bitdepth_scale = pl? p_Img->bitdepth_scale[IS_CHROMA] : p_Img->bitdepth_scale[IS_LUMA];
  int      max_imgpel_value = p_Img->max_imgpel_value_comp[pl];

  int AlphaC0Offset = MbQ->DFAlphaC0Offset;
  int BetaOffset = MbQ->DFBetaOffset;
  byte fieldModeFilteringFlag;

  Macroblock *MbP;
  imgpel   *SrcPtrP, *SrcPtrQ;

  for( pel = 0 ; pel < PelNum ; ++pel )
  {
    xQ = dir ? pel : edge;
    yQ = dir ? (edge < 16 ? edge : 1) : pel;
    p_Img->getNeighbour(MbQ, xQ - (dir_m1), yQ - dir, p_Img->mb_size[IS_LUMA], &pixP);     

    if (pixP.available || (MbQ->DFDisableIdc== 0))
    {
      if( (Strng = Strength[pel]) != 0)
      {
        p_Img->getNeighbour(MbQ, xQ, yQ, p_Img->mb_size[IS_LUMA], &pixQ);

        MbP = &(p_Img->mb_data[pixP.mb_addr]);
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
void EdgeLoopChromaNormal(imgpel** Img, byte Strength[16], Macroblock *MbQ, int dir, int edge, int uv, StorablePicture *p)
{
  int      pel, Strng ;
  int      inc_dim;
  int      tc0, dif;
  imgpel   L1, L0, R0, R1;
  int      Alpha, Beta;
  const    byte* ClipTab;
  int      indexA, indexB;
  int      StrengthIdx;
  PixelPos pixP, pixQ, pixMB1, pixMB2;
  int      QP;
  int      PelNum = pelnum_cr[dir][p->chroma_format_idc];
  ImageParameters *p_Img = MbQ->p_Img;
  
  int      bitdepth_scale   = p_Img->bitdepth_scale[IS_CHROMA];
  int      max_imgpel_value = p_Img->max_imgpel_value_comp[uv + 1];

  int xQ = dir ? 0 : edge - 1;
  int yQ = dir ? (edge < 16 ? edge - 1: 0) : 0;
  Macroblock *MbP;
  imgpel   *SrcPtrP, *SrcPtrQ;
  
  p_Img->getNeighbour(MbQ, xQ, yQ, p_Img->mb_size[IS_CHROMA], &pixMB1);

  if (pixMB1.available || (MbQ->DFDisableIdc == 0))
  {
    int AlphaC0Offset = MbQ->DFAlphaC0Offset;
    int BetaOffset = MbQ->DFBetaOffset;
    int width = p->size_x_cr;
    int dirM1 = dir - 1;
    pixP = pixMB1;
    MbP = &(p_Img->mb_data[pixP.mb_addr]);
    yQ += dir;
    xQ -= dirM1;
    inc_dim = dir ? width : 1;

    p_Img->getNeighbour(MbQ, xQ, yQ, p_Img->mb_size[IS_CHROMA], &pixMB2);
    pixQ = pixMB2;

    // Average QP of the two blocks
    QP = (MbP->qpc[uv] + MbQ->qpc[uv] + 1) >> 1;

    indexA = iClip3(0, MAX_QP, QP + AlphaC0Offset);
    indexB = iClip3(0, MAX_QP, QP + BetaOffset);

    Alpha   = ALPHA_TABLE[indexA] * bitdepth_scale;
    Beta    = BETA_TABLE [indexB] * bitdepth_scale;
    ClipTab = CLIP_TAB[indexA];

    for( pel = 0 ; pel < PelNum ; ++pel )
    {
      StrengthIdx = (PelNum == 8) ? (((pel >> 1) << 2) + (pel & 0x01)) : pel;

      if( (Strng = Strength[StrengthIdx]) != 0)
      {
        SrcPtrP = &(Img[pixP.pos_y][pixP.pos_x]);
        L1  = *(SrcPtrP - inc_dim);
        L0  = *SrcPtrP;
        SrcPtrQ = &(Img[pixQ.pos_y][pixQ.pos_x]);
        R0  = *SrcPtrQ;
        R1  = *(SrcPtrQ + inc_dim);

        if (( iabs( R0 - L0 ) < Alpha ) && ( iabs(R0 - R1) < Beta )  && ( iabs(L0 - L1) < Beta )  )
        {
          if( Strng == 4 )    // INTRA strong filtering
          {
            *SrcPtrP = (imgpel) ( ((L1 << 1) + L0 + R1 + 2) >> 2 );
            *SrcPtrQ = (imgpel) ( ((R1 << 1) + R0 + L1 + 2) >> 2 );
          }
          else
          {
            tc0  = ClipTab[ Strng ] * bitdepth_scale + 1;
            dif = iClip3( -tc0, tc0, ( ((R0 - L0) << 2) + (L1 - R1) + 4) >> 3 );

            *SrcPtrP = (imgpel) iClip1 ( max_imgpel_value, L0 + dif) ;
            *SrcPtrQ = (imgpel) iClip1 ( max_imgpel_value, R0 - dif) ;
          }
        }
      }
      pixP.pos_x += dir;
      pixQ.pos_x += dir;
      pixP.pos_y -= dirM1;
      pixQ.pos_y -= dirM1;
    }
  }
}

/*!
*****************************************************************************************
* \brief
*    Filters chroma block edge for MBAFF types
*****************************************************************************************
 */
void EdgeLoopChromaMBAff(imgpel** Img, byte Strength[16], Macroblock *MbQ, int dir, int edge, int uv, StorablePicture *p)
{
  int      pel, Strng ;
  int      incP, incQ;
  int      C0, tc0, dif;
  imgpel   L1, L0, R0, R1;
  int      Alpha = 0, Beta = 0;
  const byte* ClipTab = NULL;
  int      indexA, indexB;
  int      PelNum = pelnum_cr[dir][p->chroma_format_idc];
  int      StrengthIdx;
  int      QP;
  int      xQ, yQ;
  PixelPos pixP, pixQ;
  int      dir_m1 = 1 - dir;
  ImageParameters *p_Img = MbQ->p_Img;
  int      bitdepth_scale = p_Img->bitdepth_scale[IS_CHROMA];
  int      max_imgpel_value = p_Img->max_imgpel_value_comp[uv + 1];
  
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
    p_Img->getNeighbour(MbQ, xQ, yQ, p_Img->mb_size[IS_CHROMA], &pixQ);
    p_Img->getNeighbour(MbQ, xQ - (dir_m1), yQ - dir, p_Img->mb_size[IS_CHROMA], &pixP);    
    MbP = &(p_Img->mb_data[pixP.mb_addr]);    
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
          if( ((iabs( R0 - R1) - Beta )  & (iabs(L0 - L1) - Beta )) < 0  )
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

int compute_deblock_strength(char **list0_refIdxArr, 
                             char **list1_refIdxArr, 
                             int64 **list0_refPicIdArr, 
                             int64 **list1_refPicIdArr, 
                             short  ***list0_mv,
                             short  ***list1_mv,
                             int blk_y, int blk_x, 
                             int blk_y2, int blk_x2, 
                             int mvlimit)
{
  int64 ref_p0, ref_q0, ref_p1, ref_q1;
  byte StrValue;

  ref_p0 = list0_refIdxArr[blk_y] [blk_x] <0 ? INT64_MIN : list0_refPicIdArr[blk_y] [blk_x];
  ref_q0 = list0_refIdxArr[blk_y2][blk_x2]<0 ? INT64_MIN : list0_refPicIdArr[blk_y2][blk_x2];
  ref_p1 = list1_refIdxArr[blk_y] [blk_x] <0 ? INT64_MIN : list1_refPicIdArr[blk_y] [blk_x];
  ref_q1 = list1_refIdxArr[blk_y2][blk_x2]<0 ? INT64_MIN : list1_refPicIdArr[blk_y2][blk_x2];
  if ( ((ref_p0==ref_q0) && (ref_p1==ref_q1)) || ((ref_p0==ref_q1) && (ref_p1==ref_q0)))
  {
    // L0 and L1 reference pictures of p0 are different; q0 as well
    if (ref_p0 != ref_p1)
    {
      // compare MV for the same reference picture
      if (ref_p0 == ref_q0)
      {
        StrValue =  (byte) (
          (iabs( list0_mv[blk_y][blk_x][0] - list0_mv[blk_y2][blk_x2][0]) >= 4) |
          (iabs( list0_mv[blk_y][blk_x][1] - list0_mv[blk_y2][blk_x2][1]) >= mvlimit) |
          (iabs( list1_mv[blk_y][blk_x][0] - list1_mv[blk_y2][blk_x2][0]) >= 4) |
          (iabs( list1_mv[blk_y][blk_x][1] - list1_mv[blk_y2][blk_x2][1]) >= mvlimit));                  
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
 return StrValue;
}

