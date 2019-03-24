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
 * \file loopFilter.c
 *
 * \brief
 *    Filter to reduce blocking artifacts on a macroblock level.
 *    The filter strengh is QP dependent.
 *
 * \author
 *    Contributors:
 *    - Peter List      Peter.List@t-systems.de:  Original code                                 (13-Aug-2001)
 *    - Jani Lainema    Jani.Lainema@nokia.com:   Some bug fixing, removal of recusiveness      (16-Aug-2001)
 *    - Peter List      Peter.List@t-systems.de:  inplace filtering and various simplifications (10-Jan-2002)
 *    - Anthony Joch    anthony@ubvideo.com:      Simplified switching between filters and 
 *                                                non-recursive default filter.                 (08-Jul-2002)
 *
 *************************************************************************************
 */

#include <stdlib.h>
#include <string.h>
#include "global.h"
extern const byte QP_SCALE_CR[52] ;

/*********************************************************************************************************/

#define  IClip( Min, Max, Val) (((Val)<(Min))? (Min):(((Val)>(Max))? (Max):(Val)))

// NOTE: to change the tables below for instance when the QP doubling is changed from 6 to 8 values 
//       send an e-mail to Peter.List@t-systems.com to get a little programm that calculates them automatically 


byte ALPHA_TABLE[40]  = {0,0,0,0,0,0,0,0,  4,4,5,6,7,9,10,12,  14,17,20,24,28,33,39,46,  55,65,76,90,106,126,148,175,  207,245,255,255,255,255,255,255} ;
byte  BETA_TABLE[40]  = {0,0,0,0,0,0,0,0,  3,3,3,4,4,4, 6, 6,   7, 7, 8, 8, 9, 9,10,10,  11,11,12,12, 13, 13, 14, 14,   15, 15, 16, 16, 17, 17, 18, 18} ;
byte CLIP_TAB[40][5]  =
 {{ 0, 0, 0, 0, 0},{ 0, 0, 0, 0, 0},{ 0, 0, 0, 0, 0},{ 0, 0, 0, 0, 0},{ 0, 0, 0, 0, 0},{ 0, 0, 0, 1, 1},{ 0, 0, 0, 1, 1},{ 0, 0, 0, 1, 1},
  { 0, 0, 0, 1, 1},{ 0, 0, 1, 1, 1},{ 0, 0, 1, 1, 1},{ 0, 1, 1, 1, 1},{ 0, 1, 1, 1, 1},{ 0, 1, 1, 1, 1},{ 0, 1, 1, 1, 1},{ 0, 1, 1, 2, 2},
  { 0, 1, 1, 2, 2},{ 0, 1, 1, 2, 2},{ 0, 1, 1, 2, 2},{ 0, 1, 2, 3, 3},{ 0, 1, 2, 3, 3},{ 0, 2, 2, 3, 3},{ 0, 2, 2, 4, 4},{ 0, 2, 3, 4, 4},
  { 0, 2, 3, 4, 4},{ 0, 3, 3, 5, 5},{ 0, 3, 4, 6, 6},{ 0, 3, 4, 6, 6},{ 0, 4, 5, 7, 7},{ 0, 4, 5, 8, 8},{ 0, 4, 6, 9, 9},{ 0, 5, 7,10,10},
  { 0, 6, 8,11,11},{ 0, 6, 8,13,13},{ 0, 7,10,14,14},{ 0, 8,11,16,16},{ 0, 9,12,18,18},{ 0,10,13,20,20},{ 0,11,15,23,23},{ 0,13,17,25,25}} ;


void GetStrength(byte Strength[4],byte LargeBlockEdge[4],struct img_par *img,Macroblock* MbP,Macroblock* MbQ,int dir,int edge,int mb_y,int mb_x,byte blkmode[2][2]);
void EdgeLoop(byte* SrcPtr,byte Strength[4],byte LargeBlockEdge[4],int QP,int dir,int width,int yuv);
void DeblockMb(ImageParameters *img, byte **imgY, byte ***imgUV, int mb_y, int mb_x) ;
void get_deblock_modes(ImageParameters *img,Macroblock *MB,byte blkmode[2][2]);

/*!
 *****************************************************************************************
 * \brief
 *    The main MB-filtering function
 *****************************************************************************************
 */
void DeblockFrame(ImageParameters *img, byte **imgY, byte ***imgUV)
  {
  int       mb_x, mb_y ;

  for( mb_y=0 ; mb_y<(img->height>>4) ; mb_y++ )
    for( mb_x=0 ; mb_x<(img->width>>4) ; mb_x++ )
      DeblockMb( img, imgY, imgUV, mb_y, mb_x ) ;
  } 


  /*!
 *****************************************************************************************
 * \brief
 *    Deblocks one macroblock
 *****************************************************************************************
 */

void DeblockMb(ImageParameters *img, byte **imgY, byte ***imgUV, int mb_y, int mb_x)
{
  int           EdgeCondition;
  int           dir,edge,QP;
  byte          Strength[4], *SrcY, *SrcU, *SrcV ;
  byte          LargeBlockEdge[4];
  byte          blkmode[2][2];
  Macroblock    *MbP, *MbQ ; 
  int           sizey;
  
  SrcY = imgY    [mb_y<<4] + (mb_x<<4) ;                                                      // pointers to source
  SrcU = imgUV[0][mb_y<<3] + (mb_x<<3) ;
  SrcV = imgUV[1][mb_y<<3] + (mb_x<<3) ;

  if (img->mb_frame_field_flag)
    MbQ  = &img->mb_data[((mb_y/2)*(img->width>>3))+(mb_y%2)+mb_x*2];                            // current Mb
  else
    MbQ  = &img->mb_data[mb_y*(img->width>>4) + mb_x] ;                                                 // current Mb

  get_deblock_modes(img,MbQ,blkmode);

  for( dir=0 ; dir<2 ; dir++ )                                             // vertical edges, than horicontal edges
  {
    EdgeCondition = (dir && mb_y) || (!dir && mb_x)  ;                    // can not filter beyond frame boundaries
    for( edge=0 ; edge<4 ; edge++ )                                            // first 4 vertical strips of 16 pel
    {                                                                                         // then  4 horicontal
      if( edge || EdgeCondition )
      {
        
        sizey = mb_y%2 ? 1:2*img->width/MB_BLOCK_SIZE-1;
        if (img->mb_frame_field_flag)
          MbP = (edge)? MbQ : ((dir)? (MbQ-sizey) : (MbQ-2) ) ;       // MbP = Mb of the remote 4x4 block
        else
          MbP = (edge)? MbQ : ((dir)? (MbQ -(img->width>>4))  : (MbQ-1) ) ;       // MbP = Mb of the remote 4x4 block

        QP  = max( 0, ((MbP->qp-SHIFT_QP) + (MbQ->qp-SHIFT_QP) ) >> 1) ;                   // Average QP of the two blocks
        GetStrength(Strength,LargeBlockEdge,img,MbP,MbQ,dir,edge,mb_y<<2,mb_x<<2,blkmode); // Strength for 4 blks in 1 stripe
        if( *((int*)Strength) )  // && (QP>= 8) )                    // only if one of the 4 Strength bytes is != 0
        {
          EdgeLoop( SrcY + (edge<<2)* ((dir)? img->width:1 ), Strength, LargeBlockEdge, QP, dir, img->width, 0) ; 
          if( (imgUV != NULL) && !(edge & 1) )
          {
            EdgeLoop( SrcU +  (edge<<1) * ((dir)? img->width_cr:1 ), Strength, LargeBlockEdge, QP_SCALE_CR[QP+SHIFT_QP]-SHIFT_QP, dir, img->width_cr, 1 ) ; 
            EdgeLoop( SrcV +  (edge<<1) * ((dir)? img->width_cr:1 ), Strength, LargeBlockEdge, QP_SCALE_CR[QP+SHIFT_QP]-SHIFT_QP, dir, img->width_cr, 1 ) ; 
          }
        }
      }
    }//end edge
  }//end loop dir
}

void get_deblock_modes(ImageParameters *img,Macroblock *MB,byte blkmode[2][2])
{
 int i,j,mode;

  if( MB->mb_type || !(img->type==INTER_IMG_1||img->type==INTER_IMG_MULT||img->type==SP_IMG_1||img->type==SP_IMG_MULT) ) //if not mode 0 or not a (S)P-frame
  {
    //not a copy block in a P-frame.
    for(j=0;j<2;j++)
      for(i=0;i<2;i++)
      {
        mode=3; //no ABT, always deblock all (depending on strength)
        if(MB->useABT[(j<<1)|i])
        {
          mode=MB->abt_mode[(j<<1)|i]; //ABT
          if(mode<0)mode=0;
        }
        blkmode[j][i]=mode;
      }
  }
  else
  {
    //is a copy block in a P-frame
    //unlike in the encoder, the useABT[] is valid here. (and we do not have the input params here)
    if(MB->useABT[0])
      for(i=0;i<4;i++)
        blkmode[i>>1][i&1]=4+i;
    else
      for(i=0;i<4;i++)
        blkmode[i>>1][i&1]=3;   //no ABT, always deblock all (depending on strength)
  }
}


  /*!
 *********************************************************************************************
 * \brief
 *    returns a buffer of 4 Strength values for one stripe in a mb (for different Frame types)
 *********************************************************************************************
 */

int  ININT_STRENGTH[4] = {0x04040404, 0x03030303, 0x03030303, 0x03030303} ; 
byte BLK_NUM[2][4][4]  = {{{0,4,8,12},{1,5,9,13},{2,6,10,14},{3,7,11,15}},{{0,1,2,3},{4,5,6,7},{8,9,10,11},{12,13,14,15}}} ;
byte BLK_4_TO_8[16]    = {0,0,1,1,0,0,1,1,2,2,3,3,2,2,3,3} ;
//   indices to ABT_has_edge[][][]: ABTmode, subblockX, subblockY.
//   content of ABT_has_edge[][][]: bit #0: left edge present, bit #1: top edge present, bit #2: right edge present, bit#3: bottom edge present.
//   In this array, the values 4 to 7 for ABTmode are used to indicate 8x8 blocks of which two outer edges do not need to be deblocked.
static const char ABT_has_edge[8][2][2]={
  {{0x03,0x06},{0x09,0x0C}} , //8*8
  {{0x0B,0x0E},{0x0B,0x0E}} , //8*4
  {{0x07,0x07},{0x0D,0x0D}} , //4*8
  {{0x0F,0x0F},{0x0F,0x0F}} , //4*4
  {{0x03,0x02},{0x01,0x00}} , //UL
  {{0x02,0x06},{0x00,0x04}} , //UR
  {{0x01,0x00},{0x09,0x08}} , //BL
  {{0x00,0x04},{0x08,0x0C}}   //BR
};

void GetStrength(byte Strength[4],byte LargeBlockEdge[4],struct img_par *img,Macroblock* MbP,Macroblock* MbQ,int dir,int edge,int block_y,int block_x,byte blkmode[2][2])
{
  int    blkP, blkQ, idx ;
  int    blk_x, blk_x2, blk_y, blk_y2 ;
  int    LBcount;
  int    blk_mode1, blk_mode2;
  int    ***fw_mv = img->fw_mv;
  int    ***bw_mv = img->bw_mv;
  int    **fw_refFrArr = img->fw_refFrArr;
  int    **bw_refFrArr = img->bw_refFrArr;


  int    mvshift = (img->mv_res ? 3 : 2);                      // Consideration of mv resolution for filter strength

  *((int*)Strength) = ININT_STRENGTH[edge] ;                     // Start with Strength=3. or Strength=4 for Mb-edge

  if(img->mb_frame_field_flag)
  {
        fw_mv = img->fw_mv_frm;
        bw_mv = img->bw_mv_frm;
        fw_refFrArr = img->fw_refFrArr_frm;
        bw_refFrArr = img->bw_refFrArr_frm;
  }


  for( idx=0 ; idx<4 ; idx++ )
  {                                                                                       // if not intra or SP-frame
    LargeBlockEdge[idx]=0;
    blkQ = BLK_NUM[dir][ edge       ][idx] ;                 // if one of the 4x4 blocks has coefs.    set Strength=2
    blkP = BLK_NUM[dir][(edge-1) & 3][idx] ; 
    blk_y  = block_y + (blkQ >> 2) ;                         // moved here for ABT
    blk_x  = block_x + (blkQ  & 3)+4 ;                       // moved here for ABT
    blk_mode1=blkmode[(blk_y&2)>>1][(blk_x&2)>>1];
    
    if(  ABT_has_edge[ blk_mode1 ][ blk_y&1 ][ blk_x&1 ]  &  (1<<dir)  )
    {
      //is outer edge of a block
      LargeBlockEdge[idx]=0;
      blk_y2 = blk_y -  dir ;                       // moved here for ABT
      blk_x2 = blk_x - !dir ;                       // moved here for ABT
      blk_mode2=blkmode[(blk_y2&2)>>1][(blk_x2&2)>>1];
      //if neighboring blocks are larger than 4x4 pixels, increase filter strength.
      LBcount=0;
      if(  !( ABT_has_edge[blk_mode1][blk_y&1][blk_x&1] & (1<<(dir+2)) )  )//if current block is 8 pix deep
        LBcount++;
      if(  !( ABT_has_edge[blk_mode2][blk_y2&1][blk_x2&1] & (1<<dir) )  )//if other block is 8 pix deep
        LBcount++;
      if(  (10>>dir)  !=  ( ABT_has_edge[blk_mode1][blk_y&1][blk_x&1] & (10>>dir) )  )//if current block a wide block
        LBcount++;
      if(  (10>>dir)  !=  ( ABT_has_edge[blk_mode2][blk_y2&1][blk_x2&1] & (10>>dir) )  )//if other block a wide block
        LBcount++;
      if(LBcount>3)LBcount=3;
      LargeBlockEdge[idx]=LBcount;

      if( (   img->type != SP_IMG_1) && (img->type != SP_IMG_MULT)  && (img->type != SI_IMG)
          && !(MbP->b8mode[ BLK_4_TO_8[blkP] ]==IBLOCK || MbP->mb_type==I16MB)
          && !(MbQ->b8mode[ BLK_4_TO_8[blkQ] ]==IBLOCK || MbQ->mb_type==I16MB) )
      {
        if( ((MbQ->cbp_blk &  (1 << blkQ )) != 0) || ((MbP->cbp_blk &  (1 << blkP)) != 0) )
          Strength[idx] = 2 ;
        else
        {                                                     // if no coefs, but vector difference >= 1 set Strength=1 
          if( (img->type == B_IMG_1)  || (img->type == B_IMG_MULT) )
          {
            Strength[idx] =  (abs( fw_mv[blk_x][blk_y][0] - fw_mv[blk_x2][blk_y2][0]) >= (1 << mvshift)) |
                             (abs( fw_mv[blk_x][blk_y][1] - fw_mv[blk_x2][blk_y2][1]) >= (1 << mvshift)) |
                             (abs( bw_mv[blk_x][blk_y][0] - bw_mv[blk_x2][blk_y2][0]) >= (1 << mvshift)) |
                             (abs( bw_mv[blk_x][blk_y][1] - bw_mv[blk_x2][blk_y2][1]) >= (1 << mvshift)) |
                                  (fw_refFrArr[blk_y][blk_x-4]   !=   fw_refFrArr[blk_y2][blk_x2-4] )|
                                  (bw_refFrArr[blk_y][blk_x-4]   !=   bw_refFrArr[blk_y2][blk_x2-4] );
          }
          else
            Strength[idx] =    (abs( img->mv[blk_x][blk_y][0] - img->mv[blk_x2][blk_y2][0]) >= (1 << mvshift) ) |
                               (abs( img->mv[blk_x][blk_y][1] - img->mv[blk_x2][blk_y2][1]) >= (1 << mvshift) ) |
                                    (refFrArr[blk_y][blk_x-4]   !=   refFrArr[blk_y2][blk_x2-4] );
        }
      }
    }
    else
    {
      //is inside a larger transform block (ABT). Do not filter here.
      Strength[idx]=0;
    }
  }
}


/*!
 *****************************************************************************************
 * \brief
 *    Filters one edge of 16 (luma) or 8 (chroma) pel
 *****************************************************************************************
 */
void EdgeLoop(byte* SrcPtr,byte Strength[4],byte LargeBlockEdge[4],int QP,int dir,int width,int yuv)
{
  int      pel, ap, aq, PtrInc, Strng ;
  int      inc, inc2, inc3, inc4 ;
  int      C0, c0, Delta, dif, AbsDelta ;
  int      L2, L1, L0, R0, R1, R2, RL0 ;
  int      Alpha = 0, Beta = 0 ;
  byte*    ClipTab = NULL;   

  PtrInc  = dir?      1 : width ;
  inc     = dir?  width : 1 ;                     // vertical filtering increment to next pixel is 1 else width
  inc2    = inc<<1 ;    
  inc3    = inc + inc2 ;    
  inc4    = inc<<2 ;

  for( pel=0 ; pel<16 ; pel++ )
  {
    if(!(pel&3))
    {
      c0=QP;                    // QP is passed to EdgeLoop() in 'old range', no SHIFT_QP correction
      if( LargeBlockEdge[pel>>2] )
      {
        c0+=LargeBlockEdge[pel>>2];
        if(c0>MAX_QP-SHIFT_QP)c0=MAX_QP-SHIFT_QP;
      }
      Alpha=ALPHA_TABLE[c0];
      Beta=BETA_TABLE[c0];  
      ClipTab=CLIP_TAB[c0];
    }
    if( (Strng = Strength[pel >> 2]) )
      {
      L0  = SrcPtr [-inc ] ;
      R0  = SrcPtr [    0] ;
      AbsDelta  = abs( Delta = R0 - L0 )  ;

      if( AbsDelta < Alpha )
        {
        C0  = ClipTab[ Strng ] ;
        L1  = SrcPtr[-inc2] ;
        R1  = SrcPtr[ inc ] ;
        if( ((abs( R0 - R1) - Beta )  & (abs(L0 - L1) - Beta )) < 0  ) 
          {
          L2  = SrcPtr[-inc3] ;
          R2  = SrcPtr[ inc2] ;
          aq  = (abs( R0 - R2) - Beta ) < 0  ;
          ap  = (abs( L0 - L2) - Beta ) < 0  ;

          RL0             = L0 + R0 ;

          if(Strng == 4 )    // INTRA strong filtering
            {
            SrcPtr[   0 ]   = aq ? ( L1 + ((R1 + RL0) << 1) +  SrcPtr[ inc2] + 4) >> 3 : ((R1 << 1) + R0 + L1 + 2) >> 2 ;
            SrcPtr[-inc ]   = ap ? ( R1 + ((L1 + RL0) << 1) +  SrcPtr[-inc3] + 4) >> 3 : ((L1 << 1) + L0 + R1 + 2) >> 2 ;

            SrcPtr[ inc ]   = aq ? ( SrcPtr[ inc3] + ((SrcPtr[ inc2] + R0 + R1) << 1) + L0 + 4) >> 3 : R1;
            SrcPtr[-inc2]   = ap ? ( SrcPtr[-inc4] + ((SrcPtr[-inc3] + L1 + L0) << 1) + R0 + 4) >> 3 : L1;

            SrcPtr[ inc2] = (aq && !yuv) ? (((SrcPtr[ inc3] + SrcPtr[ inc2]) <<1) + SrcPtr[ inc2] + R1 + RL0 + 4) >> 3 : R2;
            SrcPtr[-inc3] = (ap && !yuv) ? (((SrcPtr[-inc4] + SrcPtr[-inc3]) <<1) + SrcPtr[-inc3] + L1 + RL0 + 4) >> 3 : L2;
            }
          else                                                                                   // normal filtering
            {
            c0               = C0 + ap + aq ;
            dif              = IClip( -c0, c0, ( (Delta << 2) + (L1 - R1) + 4) >> 3 ) ;
            SrcPtr[  -inc ]  = IClip(0, 255, L0 + dif) ;
            SrcPtr[     0 ]  = IClip(0, 255, R0 - dif) ;

            if( !yuv )
              {
              if( ap )
                SrcPtr[-inc2] += IClip( -C0,  C0, ( L2 + (RL0 >> 1) - (L1<<1)) >> 1 ) ;
              if( aq  )
                SrcPtr[  inc] += IClip( -C0,  C0, ( R2 + (RL0 >> 1) - (R1<<1)) >> 1 ) ;
              } ;
            } ;
          } ; 
        } ;
      SrcPtr += PtrInc ;      // Increment to next set of pixel
      pel    += yuv ;
      } 
    else
      {
      SrcPtr += PtrInc << (2 - yuv) ;
      pel    += 3 ;
      }  ;
  }
}
