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
 ***********************************************************************
 *  \file
 *      block.c
 *
 *  \brief
 *      Block functions
 *
 *  \author
 *      Main contributors (see contributors.h for copyright, address and affiliation details)
 *      - Inge Lille-Langøy          <inge.lille-langoy@telenor.com>
 *      - Rickard Sjoberg            <rickard.sjoberg@era.ericsson.se>
 ***********************************************************************
 */

#include "contributors.h"

#include <stdlib.h>
#include <math.h>

#include "block.h"

/*!
 ***********************************************************************
 * \brief
 *    makes and returns 4x4 blocks with all 5 intra prediction modes
 *
 * \return
 *    DECODING_OK   decoding of intraprediction mode was sucessfull            \n
 *    SEARCH_SYNC   search next sync element as errors while decoding occured
 ***********************************************************************
 */
int intrapred(struct img_par *img,  //!< image parameters
              int ioff,             //!< ??
              int joff,             //!< ??
              int img_block_x,      //!< ??
              int img_block_y)      //!< ??
{
  int js0=0,js1,js2,i,j;
  int img_y=0,img_x=0;
  int ia[7];

  byte predmode = img->ipredmode[img_block_x+1][img_block_y+1];

  int block_available_up = (img->ipredmode[img_block_x+1][img_block_y] >=0);
  int block_available_left = (img->ipredmode[img_block_x][img_block_y+1] >=0);

  img_x=img_block_x*4;
  img_y=img_block_y*4;

  for (i=0; i<7; i++)
    ia[i]=0;

  switch (predmode)
  {
  case DC_PRED:                         // DC prediction

    js1=0;
    js2=0;
    for (i=0;i<BLOCK_SIZE;i++)
    {
      if (block_available_up)
        js1=js1+imgY[img_y-1][img_x+i];
      if (block_available_left)
        js2=js2+imgY[img_y+i][img_x-1];
    }
    if(block_available_left && block_available_up)
      js0=(js1+js2+4)/8;                // make DC prediction from both block above and to the left
    if(block_available_left && !block_available_up)
      js0=(js2+2)/4;                    // make DC prediction from block to the left
    if(!block_available_left && block_available_up)
      js0=(js1+2)/4;                    // make DC prediction from block above
    if(!block_available_left && !block_available_up)
      js0=128;                          // current block is the upper left
    // no block to make DC pred.from
    for(i=0;i<BLOCK_SIZE;i++)
      for(j=0;j<BLOCK_SIZE;j++)
        img->mpr[i+ioff][j+joff]=js0;             // store predicted 4x4 block
    break;

  case VERT_PRED:                       // vertical prediction from block above
    for(j=0;j<BLOCK_SIZE;j++)
      for(i=0;i<BLOCK_SIZE;i++)
        img->mpr[i+ioff][j+joff]=imgY[img_y-1][img_x+i];// store predicted 4x4 block
    break;

  case HOR_PRED:                        // horisontal prediction from left block
    for(j=0;j<BLOCK_SIZE;j++)
      for(i=0;i<BLOCK_SIZE;i++)
        img->mpr[i+ioff][j+joff]=imgY[img_y+j][img_x-1]; // store predicted 4x4 block
    break;

  case DIAG_PRED_LR_45:       // diagonal prediction from left to right 45 degree
    ia[0]=(imgY[img_y+3][img_x-1]+2*imgY[img_y+2][img_x-1]+imgY[img_y+1][img_x-1]+2)/4;
    ia[1]=(imgY[img_y+2][img_x-1]+2*imgY[img_y+1][img_x-1]+imgY[img_y+0][img_x-1]+2)/4;
    ia[2]=(imgY[img_y+1][img_x-1]+2*imgY[img_y+0][img_x-1]+imgY[img_y-1][img_x-1]+2)/4;
    ia[3]=(imgY[img_y+0][img_x-1]+2*imgY[img_y-1][img_x-1]+imgY[img_y-1][img_x+0]+2)/4;
    ia[4]=(imgY[img_y-1][img_x-1]+2*imgY[img_y-1][img_x+0]+imgY[img_y-1][img_x+1]+2)/4;
    ia[5]=(imgY[img_y-1][img_x+0]+2*imgY[img_y-1][img_x+1]+imgY[img_y-1][img_x+2]+2)/4;
    ia[6]=(imgY[img_y-1][img_x+1]+2*imgY[img_y-1][img_x+2]+imgY[img_y-1][img_x+3]+2)/4;
    for(i=0;i<BLOCK_SIZE;i++)
      for(j=0;j<BLOCK_SIZE;j++)
        img->mpr[i+ioff][j+joff]=ia[i-j+3];
    break;

  case DIAG_PRED_RL: // diagonal prediction 22.5 deg to vertical plane
    ia[0]=(imgY[img_y-1][img_x+0]+imgY[img_y-1][img_x+1])/2;
    ia[1]= imgY[img_y-1][img_x+1];
    ia[2]=(imgY[img_y-1][img_x+1]+imgY[img_y-1][img_x+2])/2;
    ia[3]= imgY[img_y-1][img_x+2];
    ia[4]=(imgY[img_y-1][img_x+2]+imgY[img_y-1][img_x+3])/2;
    ia[5]= imgY[img_y-1][img_x+3];
    for(i=0;i<BLOCK_SIZE;i++)
      for(j=0;j<BLOCK_SIZE;j++)
        img->mpr[i+ioff][j+joff]=ia[MAP[j][i]];
    break;

  case  DIAG_PRED_LR:// diagonal prediction -22.5 deg to horizontal plane
    ia[0]=(imgY[img_y+0][img_x-1]+imgY[img_y+1][img_x-1])/2;
    ia[1]= imgY[img_y+1][img_x-1];
    ia[2]=(imgY[img_y+1][img_x-1]+imgY[img_y+2][img_x-1])/2;
    ia[3]= imgY[img_y+2][img_x-1];
    ia[4]=(imgY[img_y+2][img_x-1]+imgY[img_y+3][img_x-1])/2;
    ia[5]= imgY[img_y+3][img_x-1];
    for(i=0;i<BLOCK_SIZE;i++)
      for(j=0;j<BLOCK_SIZE;j++)
        img->mpr[i+ioff][j+joff]=ia[MAP[i][j]];
    break;

  default:
    printf("Error: illegal prediction mode input: %d\n",predmode);
    return SEARCH_SYNC;
    break;
  }
  return DECODING_OK;
}


/*!
 ***********************************************************************
 * \return
 *    best SAD
 ***********************************************************************
 */
int intrapred_luma_2(struct img_par *img, //!< image parameters
                     int predmode)        //!< prediction mode
{
  int s0=0,s1,s2;

  int i,j;

  int ih,iv;
  int ib,ic,iaa;

  int mb_nr = img->current_mb_nr;
  int mb_width = img->width/16;
  int mb_available_up = (img->mb_y == 0) ? 0 : (img->mb_data[mb_nr].slice_nr == img->mb_data[mb_nr-mb_width].slice_nr);
  int mb_available_left = (img->mb_x == 0) ? 0 : (img->mb_data[mb_nr].slice_nr == img->mb_data[mb_nr-1].slice_nr);

  if(img->UseConstrainedIntraPred)
  {
    // If neighbor MB is Inter, flag as not available
    if (mb_available_up && (img->intra_mb[mb_nr-mb_width] ==0))
      mb_available_up = 0;
    if (mb_available_left && (img->intra_mb[mb_nr-1] ==0))
      mb_available_left = 0;
  }

  s1=s2=0;

  switch (predmode)
  {
  case VERT_PRED_16:                       // vertical prediction from block above
    for(j=0;j<MB_BLOCK_SIZE;j++)
      for(i=0;i<MB_BLOCK_SIZE;i++)
        img->mpr[i][j]=imgY[img->pix_y-1][img->pix_x+i];// store predicted 16x16 block
    break;

  case HOR_PRED_16:                        // horisontal prediction from left block
    for(j=0;j<MB_BLOCK_SIZE;j++)
      for(i=0;i<MB_BLOCK_SIZE;i++)
        img->mpr[i][j]=imgY[img->pix_y+j][img->pix_x-1]; // store predicted 16x16 block
    break;

  case DC_PRED_16:                         // DC prediction
    s1=s2=0;
    for (i=0; i < MB_BLOCK_SIZE; i++)
    {
      if (mb_available_up)
        s1 += imgY[img->pix_y-1][img->pix_x+i];    // sum hor pix
      if (mb_available_left)
        s2 += imgY[img->pix_y+i][img->pix_x-1];    // sum vert pix
    }
    if (mb_available_up && mb_available_left)
      s0=(s1+s2+16)/(2*MB_BLOCK_SIZE);       // no edge
    if (!mb_available_up && mb_available_left)
      s0=(s2+8)/MB_BLOCK_SIZE;              // upper edge
    if (mb_available_up && !mb_available_left)
      s0=(s1+8)/MB_BLOCK_SIZE;              // left edge
    if (!mb_available_up && !mb_available_left)
      s0=128;                            // top left corner, nothing to predict from
    for(i=0;i<MB_BLOCK_SIZE;i++)
      for(j=0;j<MB_BLOCK_SIZE;j++)
      {
        img->mpr[i][j]=s0;
      }
    break;
  case PLANE_16:// 16 bit integer plan pred
    ih=0;
    iv=0;
    for (i=1;i<9;i++)
    {
      ih += i*(imgY[img->pix_y-1][img->pix_x+7+i] - imgY[img->pix_y-1][img->pix_x+7-i]);
      iv += i*(imgY[img->pix_y+7+i][img->pix_x-1] - imgY[img->pix_y+7-i][img->pix_x-1]);
    }
    ib=5*(ih/4)/16;
    ic=5*(iv/4)/16;

    iaa=16*(imgY[img->pix_y-1][img->pix_x+15]+imgY[img->pix_y+15][img->pix_x-1]);
    for (j=0;j< MB_BLOCK_SIZE;j++)
    {
      for (i=0;i< MB_BLOCK_SIZE;i++)
      {
        img->mpr[i][j]=max(0,min((iaa+(i-7)*ib +(j-7)*ic + 16)/32,255));
      }
    }// store plane prediction
    break;

  default:
    {                                    // indication of fault in bitstream,exit
      printf("Error: illegal prediction mode input: %d\n",predmode);
      return SEARCH_SYNC;
    }
  }

  return DECODING_OK;
}

/*!
 ***********************************************************************
 * \brief
 *    Inverse 4x4 transformation, transforms cof to m7
 ***********************************************************************
 */
void itrans(struct img_par *img, //!< image parameters
            int ioff,            //!< index to 4x4 block
            int joff,            //!<
            int i0,              //!<
            int j0)              //!<
{
  int i,j,i1,j1;
  int m5[4];
  int m6[4];

  // horizontal
  for (j=0;j<BLOCK_SIZE;j++)
  {
    for (i=0;i<BLOCK_SIZE;i++)
    {
      m5[i]=img->cof[i0][j0][i][j];
    }
    m6[0]=(m5[0]+m5[2])*13;
    m6[1]=(m5[0]-m5[2])*13;
    m6[2]=m5[1]*7-m5[3]*17;
    m6[3]=m5[1]*17+m5[3]*7;

    for (i=0;i<2;i++)
    {
      i1=3-i;
      img->m7[i][j]=m6[i]+m6[i1];
      img->m7[i1][j]=m6[i]-m6[i1];
    }
  }
  // vertical
  for (i=0;i<BLOCK_SIZE;i++)
  {
    for (j=0;j<BLOCK_SIZE;j++)
      m5[j]=img->m7[i][j];

    m6[0]=(m5[0]+m5[2])*13;
    m6[1]=(m5[0]-m5[2])*13;
    m6[2]=m5[1]*7-m5[3]*17;
    m6[3]=m5[1]*17+m5[3]*7;

    for (j=0;j<2;j++)
    {
      j1=3-j;
      img->m7[i][j]=mmax(0,mmin(255,(m6[j]+m6[j1]+img->mpr[i+ioff][j+joff]*JQQ+JQQ2)/JQQ));
      img->m7[i][j1]=mmax(0,mmin(255,(m6[j]-m6[j1]+img->mpr[i+ioff][j1+joff]*JQQ+JQQ2)/JQQ));
    }
  }

}


/*!
 ***********************************************************************
 * \brief
 *    invers  transform
 ***********************************************************************
 */
void itrans_2(
   struct img_par *img) //!< image parameters
{
  int i,j,i1,j1;
  int M5[4];
  int M6[4];

  // horizontal
  for (j=0;j<4;j++)
  {
    for (i=0;i<4;i++)
      M5[i]=img->cof[i][j][0][0];

    M6[0]=(M5[0]+M5[2])*13;
    M6[1]=(M5[0]-M5[2])*13;
    M6[2]= M5[1]*7 -M5[3]*17;
    M6[3]= M5[1]*17+M5[3]*7;

    for (i=0;i<2;i++)
    {
      i1=3-i;
      img->cof[i ][j][0][0]= M6[i]+M6[i1];
      img->cof[i1][j][0][0]=M6[i]-M6[i1];
    }
  }

  // vertical
  for (i=0;i<4;i++)
  {
    for (j=0;j<4;j++)
      M5[j]=img->cof[i][j][0][0];

    M6[0]=(M5[0]+M5[2])*13;
    M6[1]=(M5[0]-M5[2])*13;
    M6[2]= M5[1]*7 -M5[3]*17;
    M6[3]= M5[1]*17+M5[3]*7;

    for (j=0;j<2;j++)
    {
      j1=3-j;
      img->cof[i][j][0][0] = ((M6[j]+M6[j1])/8) *JQ1[img->qp];
      img->cof[i][j1][0][0]= ((M6[j]-M6[j1])/8) *JQ1[img->qp];
    }
  }
  for (j=0;j<4;j++)
  {
    for (i=0;i<4;i++)
    {
      img->cof[i][j][0][0] = 3 * img->cof[i][j][0][0]/256;
    }
  }
}


void itrans_sp(struct img_par *img,  //!< image parameters
               int ioff,             //!< index to 4x4 block
               int joff,             //!<
               int i0,               //!<
               int j0)               //!<
{
  int i,j,i1,j1;
  int m5[4];
  int m6[4];
  int predicted_block[BLOCK_SIZE][BLOCK_SIZE],Fq1q2,ilev;

  for (j=0; j< BLOCK_SIZE; j++)
    for (i=0; i< BLOCK_SIZE; i++)
    {
      predicted_block[i][j]=img->mpr[i+ioff][j+joff];
    }
  for (j=0; j < BLOCK_SIZE; j++)
  {
    for (i=0; i < 2; i++)
    {
      i1=3-i;
      m5[i]=predicted_block[i][j]+predicted_block[i1][j];
      m5[i1]=predicted_block[i][j]-predicted_block[i1][j];
    }
    predicted_block[0][j]=(m5[0]+m5[1])*13;
    predicted_block[2][j]=(m5[0]-m5[1])*13;
    predicted_block[1][j]=m5[3]*17+m5[2]*7;
    predicted_block[3][j]=m5[3]*7-m5[2]*17;
  }

  //  Vertival transform

  for (i=0; i < BLOCK_SIZE; i++)
  {
    for (j=0; j < 2; j++)
    {
      j1=3-j;
      m5[j]=predicted_block[i][j]+predicted_block[i][j1];
      m5[j1]=predicted_block[i][j]-predicted_block[i][j1];
    }
    predicted_block[i][0]=(m5[0]+m5[1])*13;
    predicted_block[i][2]=(m5[0]-m5[1])*13;
    predicted_block[i][1]=m5[3]*17+m5[2]*7;
    predicted_block[i][3]=m5[3]*7-m5[2]*17;
  }

  Fq1q2=(JQQ*JQ[img->qpsp]+JQ[img->qp]/2)/JQ[img->qp];

  for (j=0;j<BLOCK_SIZE;j++)
    for (i=0;i<BLOCK_SIZE;i++)
    {
      ilev=img->cof[i0][j0][i][j]/JQ1[img->qp]*Fq1q2+predicted_block[i][j]*JQ[img->qpsp];
      img->cof[i0][j0][i][j]=sign((abs(ilev)+JQQ2)/ JQQ,ilev)*JQ1[img->qpsp];
    }
  // horizontal
  for (j=0;j<BLOCK_SIZE;j++)
  {
    for (i=0;i<BLOCK_SIZE;i++)
    {
      m5[i]=img->cof[i0][j0][i][j];
    }
    m6[0]=(m5[0]+m5[2])*13;
    m6[1]=(m5[0]-m5[2])*13;
    m6[2]=m5[1]*7-m5[3]*17;
    m6[3]=m5[1]*17+m5[3]*7;

    for (i=0;i<2;i++)
    {
      i1=3-i;
      img->m7[i][j]=m6[i]+m6[i1];
      img->m7[i1][j]=m6[i]-m6[i1];
    }
  }
  // vertical
  for (i=0;i<BLOCK_SIZE;i++)
  {
    for (j=0;j<BLOCK_SIZE;j++)
      m5[j]=img->m7[i][j];

    m6[0]=(m5[0]+m5[2])*13;
    m6[1]=(m5[0]-m5[2])*13;
    m6[2]=m5[1]*7-m5[3]*17;
    m6[3]=m5[1]*17+m5[3]*7;

    for (j=0;j<2;j++)
    {
      j1=3-j;
      img->m7[i][j]=mmax(0,mmin(255,(m6[j]+m6[j1]+JQQ2)/JQQ));
      img->m7[i][j1]=mmax(0,mmin(255,(m6[j]-m6[j1]+JQQ2)/JQQ));
    }
  }

}

/*!
 ***********************************************************************
 * \brief
 *    The routine performs transform,quantization,inverse transform, adds the diff.
 *    to the prediction and writes the result to the decoded luma frame. Includes the
 *    RD constrained quantization also.
 *
 * \par Input:
 *    block_x,block_y: Block position inside a macro block (0,4,8,12).
 *
 * \par Output:
 *    nonzero: 0 if no levels are nonzero.  1 if there are nonzero levels. \n
 *    coeff_cost: Counter for nonzero coefficients, used to discard expencive levels.
 ************************************************************************
 */
void copyblock_sp(struct img_par *img,int block_x,int block_y)
{
  int sign(int a,int b);

  int i,j,i1,j1,m5[4],m6[4];

  int predicted_block[BLOCK_SIZE][BLOCK_SIZE],quant_set1;

  //  Horizontal transform
  for (j=0; j< BLOCK_SIZE; j++)
    for (i=0; i< BLOCK_SIZE; i++)
    {
      predicted_block[i][j]=img->mpr[i+block_x][j+block_y];
    }

  for (j=0; j < BLOCK_SIZE; j++)
  {
    for (i=0; i < 2; i++)
    {
      i1=3-i;
      m5[i]=predicted_block[i][j]+predicted_block[i1][j];
      m5[i1]=predicted_block[i][j]-predicted_block[i1][j];
    }
    predicted_block[0][j]=(m5[0]+m5[1])*13;
    predicted_block[2][j]=(m5[0]-m5[1])*13;
    predicted_block[1][j]=m5[3]*17+m5[2]*7;
    predicted_block[3][j]=m5[3]*7-m5[2]*17;
  }

  //  Vertical transform

  for (i=0; i < BLOCK_SIZE; i++)
  {
    for (j=0; j < 2; j++)
    {
      j1=3-j;
      m5[j]=predicted_block[i][j]+predicted_block[i][j1];
      m5[j1]=predicted_block[i][j]-predicted_block[i][j1];
    }
    predicted_block[i][0]=(m5[0]+m5[1])*13;
    predicted_block[i][2]=(m5[0]-m5[1])*13;
    predicted_block[i][1]=m5[3]*17+m5[2]*7;
    predicted_block[i][3]=m5[3]*7-m5[2]*17;
  }

  // Quant
  quant_set1=img->qpsp;
  for (j=0;j < BLOCK_SIZE; j++)
    for (i=0; i < BLOCK_SIZE; i++)
      img->m7[i][j]=sign((abs(predicted_block[i][j])*JQ[quant_set1]+JQQ2)/ JQQ,predicted_block[i][j])*JQ1[quant_set1];

  //     IDCT.
  //     horizontal

  for (j=0; j < BLOCK_SIZE; j++)
  {
    for (i=0; i < BLOCK_SIZE; i++)
    {
      m5[i]=img->m7[i][j];
    }
    m6[0]=(m5[0]+m5[2])*13;
    m6[1]=(m5[0]-m5[2])*13;
    m6[2]=m5[1]*7-m5[3]*17;
    m6[3]=m5[1]*17+m5[3]*7;

    for (i=0; i < 2; i++)
    {
      i1=3-i;
      img->m7[i][j]=m6[i]+m6[i1];
      img->m7[i1][j]=m6[i]-m6[i1];
    }
  }

  //  vertical

  for (i=0; i < BLOCK_SIZE; i++)
  {
    for (j=0; j < BLOCK_SIZE; j++)
    {
      m5[j]=img->m7[i][j];
    }
    m6[0]=(m5[0]+m5[2])*13;
    m6[1]=(m5[0]-m5[2])*13;
    m6[2]=m5[1]*7-m5[3]*17;
    m6[3]=m5[1]*17+m5[3]*7;

    for (j=0; j < 2; j++)
    {
      j1=3-j;
      img->m7[i][j] =min(255,max(0,(m6[j]+m6[j1]+JQQ2)/JQQ));
      img->m7[i][j1]=min(255,max(0,(m6[j]-m6[j1]+JQQ2)/JQQ));
    }
  }

  //  Decoded block moved to frame memory

  for (j=0; j < BLOCK_SIZE; j++)
    for (i=0; i < BLOCK_SIZE; i++)
      imgY[img->pix_y+block_y+j][img->pix_x+block_x+i]=img->m7[i][j];

}

void itrans_sp_chroma(struct img_par *img,int ll)
{
  int i,j,i1,j2,ilev,n2,n1,j1,mb_y,qp_const,quant_set;
  int m5[BLOCK_SIZE];
  int predicted_chroma_block[MB_BLOCK_SIZE/2][MB_BLOCK_SIZE/2],Fq1q2,mp1[BLOCK_SIZE],quant_set1;

  qp_const=JQQ4;

  for (j=0; j < MB_BLOCK_SIZE/2; j++)
    for (i=0; i < MB_BLOCK_SIZE/2; i++)
    {
      predicted_chroma_block[i][j]=img->mpr[i][j];
      img->mpr[i][j]=0;
    }

  for (n2=0; n2 <= BLOCK_SIZE; n2 += BLOCK_SIZE)
  {
    for (n1=0; n1 <= BLOCK_SIZE; n1 += BLOCK_SIZE)
    {

      //  Horizontal transform.
      for (j=0; j < BLOCK_SIZE; j++)
      {
        mb_y=n2+j;
        for (i=0; i < 2; i++)
        {
          i1=3-i;
          m5[i]=predicted_chroma_block[i+n1][mb_y]+predicted_chroma_block[i1+n1][mb_y];
          m5[i1]=predicted_chroma_block[i+n1][mb_y]-predicted_chroma_block[i1+n1][mb_y];
        }
        predicted_chroma_block[n1][mb_y]=(m5[0]+m5[1])*13;
        predicted_chroma_block[n1+2][mb_y]=(m5[0]-m5[1])*13;
        predicted_chroma_block[n1+1][mb_y]=m5[3]*17+m5[2]*7;
        predicted_chroma_block[n1+3][mb_y]=m5[3]*7-m5[2]*17;
      }

      //  Vertical transform.

      for (i=0; i < BLOCK_SIZE; i++)
      {
        j1=n1+i;
        for (j=0; j < 2; j++)
        {
          j2=3-j;
          m5[j]=predicted_chroma_block[j1][n2+j]+predicted_chroma_block[j1][n2+j2];
          m5[j2]=predicted_chroma_block[j1][n2+j]-predicted_chroma_block[j1][n2+j2];
        }
        predicted_chroma_block[j1][n2+0]=(m5[0]+m5[1])*13;
        predicted_chroma_block[j1][n2+2]=(m5[0]-m5[1])*13;
        predicted_chroma_block[j1][n2+1]=m5[3]*17+m5[2]*7;
        predicted_chroma_block[j1][n2+3]=m5[3]*7-m5[2]*17;
      }
    }
  }

  //     2X2 transform of DC coeffs.
  mp1[0]=(predicted_chroma_block[0][0]+predicted_chroma_block[4][0]+predicted_chroma_block[0][4]+predicted_chroma_block[4][4])/2;
  mp1[1]=(predicted_chroma_block[0][0]-predicted_chroma_block[4][0]+predicted_chroma_block[0][4]-predicted_chroma_block[4][4])/2;
  mp1[2]=(predicted_chroma_block[0][0]+predicted_chroma_block[4][0]-predicted_chroma_block[0][4]-predicted_chroma_block[4][4])/2;
  mp1[3]=(predicted_chroma_block[0][0]-predicted_chroma_block[4][0]-predicted_chroma_block[0][4]+predicted_chroma_block[4][4])/2;

  quant_set=QP_SCALE_CR[img->qp];
  quant_set1=QP_SCALE_CR[img->qpsp];
  Fq1q2=(JQQ*JQ[quant_set1]+JQ[quant_set]/2)/JQ[quant_set];

  ilev=img->cof[0+ll][4][0][0]*Fq1q2+mp1[0]*JQ[quant_set1];
  mp1[0]=sign((abs(ilev)+JQQ2)/ JQQ,ilev)*JQ1[quant_set1];
  ilev=img->cof[1+ll][4][0][0]*Fq1q2+mp1[1]*JQ[quant_set1];
  mp1[1]=sign((abs(ilev)+JQQ2)/ JQQ,ilev)*JQ1[quant_set1];
  ilev=img->cof[0+ll][5][0][0]*Fq1q2+mp1[2]*JQ[quant_set1];
  mp1[2]=sign((abs(ilev)+JQQ2)/ JQQ,ilev)*JQ1[quant_set1];
  ilev=img->cof[1+ll][5][0][0]*Fq1q2+mp1[3]*JQ[quant_set1];
  mp1[3]=sign((abs(ilev)+JQQ2)/ JQQ,ilev)*JQ1[quant_set1];

  img->cof[0+ll][4][0][0]=(mp1[0]+mp1[1]+mp1[2]+mp1[3])/2;
  img->cof[1+ll][4][0][0]=(mp1[0]-mp1[1]+mp1[2]-mp1[3])/2;
  img->cof[0+ll][5][0][0]=(mp1[0]+mp1[1]-mp1[2]-mp1[3])/2;
  img->cof[1+ll][5][0][0]=(mp1[0]-mp1[1]-mp1[2]+mp1[3])/2;

  for (n2=0; n2 <= BLOCK_SIZE; n2 += BLOCK_SIZE)
    for (n1=0; n1 <= BLOCK_SIZE; n1 += BLOCK_SIZE)
      for (i=0;i< BLOCK_SIZE; i++)
        for (j=0;j< BLOCK_SIZE; j++)
        if ((i!=0) || (j!=0))
        {
          ilev=img->cof[n1/BLOCK_SIZE+ll][4+n2/BLOCK_SIZE][i][j]/JQ1[QP_SCALE_CR[img->qp]]*Fq1q2+predicted_chroma_block[n1+i][n2+j]*JQ[quant_set1];
          img->cof[n1/BLOCK_SIZE+ll][4+n2/BLOCK_SIZE][i][j]=sign((abs(ilev)+JQQ2)/ JQQ,ilev)*JQ1[quant_set1];
        }
}

int sign(int a , int b)
{
  int x;

  x=abs(a);
  if (b>0)
    return(x);
  else return(-x);
}
