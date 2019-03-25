
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
 *    - Peter List       Peter.List@t-systems.de:  Original code                                 (13-Aug-2001)
 *    - Jani Lainema     Jani.Lainema@nokia.com:   Some bug fixing, removal of recusiveness      (16-Aug-2001)
 *    - Peter List       Peter.List@t-systems.de:  inplace filtering and various simplifications (10-Jan-2002)
 *    - Anthony Joch     anthony@ubvideo.com:      Simplified switching between filters and 
 *                                                 non-recursive default filter.                 (08-Jul-2002)
 *    - Cristina Gomila  cristina.gomila@thomson.net: Simplification of the chroma deblocking
 *                                                    from JVT-E089                              (21-Nov-2002)
 *************************************************************************************
 */
#include <stdlib.h>
#include <string.h>
#include "global.h"
#include "image.h"
#include "mb_access.h"

extern const byte QP_SCALE_CR[52] ;

/*********************************************************************************************************/

#define  IClip( Min, Max, Val) (((Val)<(Min))? (Min):(((Val)>(Max))? (Max):(Val)))

// NOTE: to change the tables below for instance when the QP doubling is changed from 6 to 8 values 
//       send an e-mail to Peter.List@t-systems.com to get a little programm that calculates them automatically 

byte ALPHA_TABLE[52]  = {0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,4,4,5,6,  7,8,9,10,12,13,15,17,  20,22,25,28,32,36,40,45,  50,56,63,71,80,90,101,113,  127,144,162,182,203,226,255,255} ;
byte  BETA_TABLE[52]  = {0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,2,2,2,3,  3,3,3, 4, 4, 4, 6, 6,   7, 7, 8, 8, 9, 9,10,10,  11,11,12,12,13,13, 14, 14,   15, 15, 16, 16, 17, 17, 18, 18} ;
byte CLIP_TAB[52][5]  =
{
  { 0, 0, 0, 0, 0},{ 0, 0, 0, 0, 0},{ 0, 0, 0, 0, 0},{ 0, 0, 0, 0, 0},{ 0, 0, 0, 0, 0},{ 0, 0, 0, 0, 0},{ 0, 0, 0, 0, 0},{ 0, 0, 0, 0, 0},
  { 0, 0, 0, 0, 0},{ 0, 0, 0, 0, 0},{ 0, 0, 0, 0, 0},{ 0, 0, 0, 0, 0},{ 0, 0, 0, 0, 0},{ 0, 0, 0, 0, 0},{ 0, 0, 0, 0, 0},{ 0, 0, 0, 0, 0},
  { 0, 0, 0, 0, 0},{ 0, 0, 0, 1, 1},{ 0, 0, 0, 1, 1},{ 0, 0, 0, 1, 1},{ 0, 0, 0, 1, 1},{ 0, 0, 1, 1, 1},{ 0, 0, 1, 1, 1},{ 0, 1, 1, 1, 1},
  { 0, 1, 1, 1, 1},{ 0, 1, 1, 1, 1},{ 0, 1, 1, 1, 1},{ 0, 1, 1, 2, 2},{ 0, 1, 1, 2, 2},{ 0, 1, 1, 2, 2},{ 0, 1, 1, 2, 2},{ 0, 1, 2, 3, 3},
  { 0, 1, 2, 3, 3},{ 0, 2, 2, 3, 3},{ 0, 2, 2, 4, 4},{ 0, 2, 3, 4, 4},{ 0, 2, 3, 4, 4},{ 0, 3, 3, 5, 5},{ 0, 3, 4, 6, 6},{ 0, 3, 4, 6, 6},
  { 0, 4, 5, 7, 7},{ 0, 4, 5, 8, 8},{ 0, 4, 6, 9, 9},{ 0, 5, 7,10,10},{ 0, 6, 8,11,11},{ 0, 6, 8,13,13},{ 0, 7,10,14,14},{ 0, 8,11,16,16},
  { 0, 9,12,18,18},{ 0,10,13,20,20},{ 0,11,15,23,23},{ 0,13,17,25,25}
} ;


void GetStrength(byte Strength[4],struct img_par *img,Macroblock* MbP,Macroblock* MbQ,int dir,int edge,int mb_y,int mb_x);
void EdgeLoop(byte* SrcPtr,byte Strength[4], int QP, int AlphaC0Offset, int BetaOffset, int dir,int width,int yuv);
void DeblockMb(ImageParameters *img, byte **imgY, byte ***imgUV, int mb_y, int mb_x) ;

/*!
 *****************************************************************************************
 * \brief
 *    Filter all macroblocks in order of increasing macroblock address.
 *****************************************************************************************
 */
void DeblockFrame(ImageParameters *img, byte **imgY, byte ***imgUV)
{
  unsigned i;
  int mb_x, mb_y ;

  for (i=0; i<img->PicSizeInMbs; i++)
  {
    get_mb_block_pos(i, &mb_x, &mb_y);
    DeblockMb( img, imgY, imgUV, mb_y, mb_x ) ;
  }
} 


/*!
 *****************************************************************************************
 * \brief
 *    Deblocking filter for one macroblock.
 *****************************************************************************************
 */

void DeblockMb(ImageParameters *img, byte **imgY, byte ***imgUV, int mb_y, int mb_x)
{
  int           EdgeCondition;
  int           dir,edge,QP;
  byte          Strength[4], *SrcY, *SrcU, *SrcV ;
  Macroblock    *MbP, *MbQ ; 
  int           sizey;
  int           QPC;

  int           filterLeftMbEdgeFlag  = (mb_x != 0);
  int           filterTopMbEdgeFlag   = (mb_y != 0);
  int           fieldModeMbFlag;
  
  SrcY = imgY    [mb_y<<4] + (mb_x<<4) ;                                                      // pointers to source
  SrcU = imgUV[0][mb_y<<3] + (mb_x<<3) ;
  SrcV = imgUV[1][mb_y<<3] + (mb_x<<3) ;

  MbQ  = &img->mb_data[mb_y*img->PicWidthInMbs + mb_x] ;                                                 // current Mb

  fieldModeMbFlag       = img->field_pic_flag || (img->MbaffFrameFlag && MbQ->mb_field);

  // return, if filter is disabled
  if (MbQ->LFDisableIdc==1) return;

  if (MbQ->LFDisableIdc==2)
  {
    // don't filter at slice boundaries
    filterLeftMbEdgeFlag = MbQ->mbAvailA;
    filterTopMbEdgeFlag  = MbQ->mbAvailB;;
  }

  for( dir=0 ; dir<2 ; dir++ )                                             // vertical edges, than horicontal edges
  {
    EdgeCondition = (dir && filterTopMbEdgeFlag) || (!dir && filterLeftMbEdgeFlag); // can not filter beyond picture boundaries
    for( edge=0 ; edge<4 ; edge++ )                                            // first 4 vertical strips of 16 pel
    {                                                                                         // then  4 horicontal
      if( edge || EdgeCondition )
      {
        sizey = mb_y%2 ? 1:2*img->width/MB_BLOCK_SIZE-1;
        MbP = (edge)? MbQ : ((dir)? (MbQ -(img->width>>4))  : (MbQ-1) ) ;       // MbP = Mb of the remote 4x4 block

        QP  = ( MbP->qp + MbQ->qp + 1) >> 1 ;                   // Average QP of the two blocks
        GetStrength(Strength,img,MbP,MbQ,dir,edge,mb_y<<2,mb_x<<2); // Strength for 4 blks in 1 stripe
        if( *((int*)Strength) )                      // only if one of the 4 Strength bytes is != 0
        {
          EdgeLoop( SrcY + (edge<<2)* ((dir)? img->width:1 ), Strength, QP, MbQ->LFAlphaC0Offset, MbQ->LFBetaOffset, dir, img->width, 0) ; 
          if( (imgUV != NULL) && !(edge & 1) )
          {
            QPC  = (QP_SCALE_CR[MbP->qp] + QP_SCALE_CR[MbQ->qp] + 1) >> 1;
            EdgeLoop( SrcU +  (edge<<1) * ((dir)? img->width_cr:1 ), Strength, QPC, MbQ->LFAlphaC0Offset, MbQ->LFBetaOffset, dir, img->width_cr, 1 ) ; 
            EdgeLoop( SrcV +  (edge<<1) * ((dir)? img->width_cr:1 ), Strength, QPC, MbQ->LFAlphaC0Offset, MbQ->LFBetaOffset, dir, img->width_cr, 1 ) ; 
          }
        }
      }
    }//end edge
  }//end loop dir
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

void GetStrength(byte Strength[4],struct img_par *img,Macroblock* MbP,Macroblock* MbQ,int dir,int edge,int block_y,int block_x)
{
  int    blkP, blkQ, idx ;
  int    blk_x, blk_x2, blk_y, blk_y2 ;
  int    ***list0_mv = dec_picture->mv[LIST_0];
  int    ***list1_mv = dec_picture->mv[LIST_1];
  int    **list0_refIdxArr = dec_picture->ref_idx[LIST_0];
  int    **list1_refIdxArr = dec_picture->ref_idx[LIST_1];
	int    **list0_refPicIdArr = dec_picture->ref_pic_id[LIST_0];
	int    **list1_refPicIdArr = dec_picture->ref_pic_id[LIST_1];

  if ((img->structure==FRAME) || !dir)
    *((int*)Strength) = ININT_STRENGTH[edge] ;                     // Start with Strength=3. or Strength=4 for Mb-edge
  else
    *((int*)Strength) = ININT_STRENGTH[3] ;                        //  Strength=3. for fields


  for( idx=0 ; idx<4 ; idx++ )
  {                                                                                       // if not intra or SP-frame
    blkQ = BLK_NUM[dir][ edge       ][idx] ;                 // if one of the 4x4 blocks has coefs.    set Strength=2
    blkP = BLK_NUM[dir][(edge-1) & 3][idx] ; 
    
    if(  !(MbP->mb_type==I4MB || MbP->mb_type==I16MB)
          && !(MbQ->mb_type==I4MB || MbQ->mb_type==I16MB) )
    {
      if( ((MbQ->cbp_blk &  (1 << blkQ )) != 0) || ((MbP->cbp_blk &  (1 << blkP)) != 0) )
        Strength[idx] = 2 ;
      else
      {                                                     // if no coefs, but vector difference >= 1 set Strength=1 
        blk_y  = block_y + (blkQ >> 2) ;   blk_y2 = blk_y -  dir ;
        blk_x  = block_x + (blkQ  & 3) ;   blk_x2 = blk_x - !dir ;
        if( (img->type == B_SLICE) )
        {
          int ref_p0,ref_p1,ref_q0,ref_q1;	    
          ref_p0 = list0_refIdxArr[blk_x][blk_y]<0 ? -1 : list0_refPicIdArr[blk_x][blk_y];
          ref_q0 = list0_refIdxArr[blk_x2][blk_y2]<0 ? -1 : list0_refPicIdArr[blk_x2][blk_y2];
          ref_p1 = list1_refIdxArr[blk_x][blk_y]<0 ? -1 : list1_refPicIdArr[blk_x][blk_y];
          ref_q1 = list1_refIdxArr[blk_x2][blk_y2]<0 ? -1 : list1_refPicIdArr[blk_x2][blk_y2];
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
                Strength[idx] =  (abs( list0_mv[blk_x][blk_y][0] - list0_mv[blk_x2][blk_y2][0]) >= 4) |
                  (abs( list0_mv[blk_x][blk_y][1] - list0_mv[blk_x2][blk_y2][1]) >= 4) |
                  (abs( list1_mv[blk_x][blk_y][0] - list1_mv[blk_x2][blk_y2][0]) >= 4) |
                  (abs( list1_mv[blk_x][blk_y][1] - list1_mv[blk_x2][blk_y2][1]) >= 4);
              }
              else 
              {
                Strength[idx] =  (abs( list0_mv[blk_x][blk_y][0] - list1_mv[blk_x2][blk_y2][0]) >= 4) |
                  (abs( list0_mv[blk_x][blk_y][1] - list1_mv[blk_x2][blk_y2][1]) >= 4) |
                  (abs( list1_mv[blk_x][blk_y][0] - list0_mv[blk_x2][blk_y2][0]) >= 4) |
                  (abs( list1_mv[blk_x][blk_y][1] - list0_mv[blk_x2][blk_y2][1]) >= 4);
              }	
            }
            else 
            {	// L0 and L1 reference pictures of p0 are the same; q0 as well
              
              Strength[idx] =  ((abs( list0_mv[blk_x][blk_y][0] - list0_mv[blk_x2][blk_y2][0]) >= 4) |
                (abs( list0_mv[blk_x][blk_y][1] - list0_mv[blk_x2][blk_y2][1]) >= 4 ) |
                (abs( list1_mv[blk_x][blk_y][0] - list1_mv[blk_x2][blk_y2][0]) >= 4) |
                (abs( list1_mv[blk_x][blk_y][1] - list1_mv[blk_x2][blk_y2][1]) >= 4))
                &&
                ((abs( list0_mv[blk_x][blk_y][0] - list1_mv[blk_x2][blk_y2][0]) >= 4) |
                (abs( list0_mv[blk_x][blk_y][1] - list1_mv[blk_x2][blk_y2][1]) >= 4) |
                (abs( list1_mv[blk_x][blk_y][0] - list0_mv[blk_x2][blk_y2][0]) >= 4) |
                (abs( list1_mv[blk_x][blk_y][1] - list0_mv[blk_x2][blk_y2][1]) >= 4));
            }				
          }
          else 
          {
            Strength[idx] = 1;				
          }		
        }
        else	
        {	// P slice
          int ref_p0,ref_q0;	    
          ref_p0 = list0_refIdxArr[blk_x][blk_y]<0 ? -1 : list0_refPicIdArr[blk_x][blk_y];
          ref_q0 = list0_refIdxArr[blk_x2][blk_y2]<0 ? -1 : list0_refPicIdArr[blk_x2][blk_y2];
          Strength[idx] =  (ref_p0 != ref_q0 ) |
            (abs( list0_mv[blk_x][blk_y][0] - list0_mv[blk_x2][blk_y2][0]) >= 4 ) |
            (abs( list0_mv[blk_x][blk_y][1] - list0_mv[blk_x2][blk_y2][1]) >= 4 );
        }
      }
    }
  }
}


/*!
 *****************************************************************************************
 * \brief
 *    Filters one edge of 16 (luma) or 8 (chroma) pel
 *****************************************************************************************
 */
void EdgeLoop(byte* SrcPtr,byte Strength[4],int QP,
              int AlphaC0Offset, int BetaOffset, int dir,int width,int yuv)
{
  int      pel, ap = 0, aq = 0, PtrInc, Strng ;
  int      inc, inc2, inc3, inc4 ;
  int      C0, c0, Delta, dif, AbsDelta ;
  int      L2 = 0, L1, L0, R0, R1, R2 = 0, RL0 ;
  int      Alpha = 0, Beta = 0 ;
  byte*    ClipTab = NULL;   
  int      small_gap;
  int      indexA, indexB;
  
  
  PtrInc  = dir?      1 : width ;
  inc     = dir?  width : 1 ;                     // vertical filtering increment to next pixel is 1 else width
  inc2    = inc<<1 ;    
  inc3    = inc + inc2 ;    
  inc4    = inc<<2 ;
  
  for( pel=0 ; pel<16 ; pel++ )
  {
    if(!(pel&3))
    {
      indexA = IClip(0, MAX_QP, QP + AlphaC0Offset);
      indexB = IClip(0, MAX_QP, QP + BetaOffset);
      
      Alpha=ALPHA_TABLE[indexA];
      Beta=BETA_TABLE[indexB];  
      ClipTab=CLIP_TAB[indexA];
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
          if( !yuv)
          {
            L2  = SrcPtr[-inc3] ;
            R2  = SrcPtr[ inc2] ;
            aq  = (abs( R0 - R2) - Beta ) < 0  ;
            ap  = (abs( L0 - L2) - Beta ) < 0  ;
          }
          
          RL0             = L0 + R0 ;
          
          if(Strng == 4 )    // INTRA strong filtering
          {
            if( yuv)  // Chroma
            {
              SrcPtr[   0 ] = ((R1 << 1) + R0 + L1 + 2) >> 2; 
              SrcPtr[-inc ] = ((L1 << 1) + L0 + R1 + 2) >> 2;                                           
            }
            else  // Luma
            {
              small_gap = (AbsDelta < ((Alpha >> 2) + 2));
              
              aq &= small_gap;
              ap &= small_gap;
              
              SrcPtr[   0 ]   = aq ? ( L1 + ((R1 + RL0) << 1) +  SrcPtr[ inc2] + 4) >> 3 : ((R1 << 1) + R0 + L1 + 2) >> 2 ;
              SrcPtr[-inc ]   = ap ? ( R1 + ((L1 + RL0) << 1) +  SrcPtr[-inc3] + 4) >> 3 : ((L1 << 1) + L0 + R1 + 2) >> 2 ;
              
              SrcPtr[ inc ] =   aq  ? ( SrcPtr[ inc2] + R0 + R1 + L0 + 2) >> 2 : SrcPtr[ inc ];
              SrcPtr[-inc2] =   ap  ? ( SrcPtr[-inc3] + L1 + L0 + R0 + 2) >> 2 : SrcPtr[-inc2];
              
              SrcPtr[ inc2] = aq ? (((SrcPtr[ inc3] + SrcPtr[ inc2]) <<1) + SrcPtr[ inc2] + R1 + RL0 + 4) >> 3 : R2;
              SrcPtr[-inc3] = ap ? (((SrcPtr[-inc4] + SrcPtr[-inc3]) <<1) + SrcPtr[-inc3] + L1 + RL0 + 4) >> 3 : L2;
            }
          }
          else                                                                                   // normal filtering
          {
            c0               = yuv? (C0+1):(C0 + ap + aq) ;
            dif              = IClip( -c0, c0, ( (Delta << 2) + (L1 - R1) + 4) >> 3 ) ;
            SrcPtr[  -inc ]  = IClip(0, 255, L0 + dif) ;
            SrcPtr[     0 ]  = IClip(0, 255, R0 - dif) ;
            
            if( !yuv )
            {
              if( ap )
                SrcPtr[-inc2] += IClip( -C0,  C0, ( L2 + ((RL0 + 1) >> 1) - (L1<<1)) >> 1 ) ;
              if( aq  )
                SrcPtr[  inc] += IClip( -C0,  C0, ( R2 + ((RL0 + 1) >> 1) - (R1<<1)) >> 1 ) ;
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


