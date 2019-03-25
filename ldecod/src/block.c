
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
 *      - Inge Lille-Langoy          <inge.lille-langoy@telenor.com>
 *      - Rickard Sjoberg            <rickard.sjoberg@era.ericsson.se>
 ***********************************************************************
 */

#include "contributors.h"

#include "global.h"
#include "block.h"
#include "blk_prediction.h"
#include "image.h"
#include "mb_access.h"
#include "transform.h"
#include "quant.h"
#include "memalloc.h"

/*!
 ***********************************************************************
 * \brief
 *    Inverse 4x4 transformation, transforms cof to mb_rres
 ***********************************************************************
 */
void itrans4x4(Macroblock *currMB,   //!< current macroblock
               ColorPlane pl,        //!< used color plane
               int ioff,             //!< index to 4x4 block
               int joff)             //!< index to 4x4 block
{
  Slice *currSlice = currMB->p_Slice;
  int    **mb_rres = currSlice->mb_rres[pl];

  inverse4x4(currSlice->cof[pl],mb_rres,joff,ioff);

  sample_reconstruct (&currSlice->mb_rec[pl][joff], &currSlice->mb_pred[pl][joff], &mb_rres[joff], ioff, ioff, BLOCK_SIZE, BLOCK_SIZE, currMB->p_Vid->max_pel_value_comp[pl], DQ_BITS);
}

/*!
 ****************************************************************************
 * \brief
 *    Inverse 4x4 lossless_qpprime transformation, transforms cof to mb_rres
 ****************************************************************************
 */
void itrans4x4_ls(Macroblock *currMB,   //!< current macroblock
                  ColorPlane pl,        //!< Color plane (for 4:4:4)                  
                  int ioff,             //!< index to 4x4 block
                  int joff)             //!< index to 4x4 block
{
  int i,j;

  Slice *currSlice = currMB->p_Slice;
  VideoParameters *p_Vid = currMB->p_Vid;
  int max_imgpel_value = p_Vid->max_pel_value_comp[pl];

  imgpel **mb_pred = currSlice->mb_pred[pl];
  imgpel **mb_rec  = currSlice->mb_rec[pl];
  int    **mb_rres = currSlice->mb_rres [pl];

  inverse4x4(currSlice->cof[pl], mb_rres, joff, ioff);

  for (j = joff; j < joff + BLOCK_SIZE; ++j)
  {
    for (i = ioff; i < ioff + BLOCK_SIZE; ++i)
    {      
      mb_rec[j][i] = (imgpel) iClip1(max_imgpel_value, mb_pred[j][i] + rshift_rnd_sf(mb_rres[j][i], DQ_BITS));
    }
  }
}

/*!
************************************************************************
* \brief
*    Inverse residual DPCM for Intra lossless coding
*
************************************************************************
*/
void Inv_Residual_trans_4x4(Macroblock *currMB,   //!< current macroblock
                            ColorPlane pl,        //!< used color plane
                            int ioff,             //!< index to 4x4 block
                            int joff)             //!< index to 4x4 block
{
  int i,j;
  int temp[4][4];
  Slice *currSlice = currMB->p_Slice;
  imgpel **mb_pred = currSlice->mb_pred[pl];
  imgpel **mb_rec  = currSlice->mb_rec[pl];
  int    **mb_rres = currSlice->mb_rres[pl];
  int    **cof     = currSlice->cof[pl];

  if(currMB->ipmode_DPCM == VERT_PRED)
  {
    for(i=0; i<4; ++i)
    {
      temp[0][i] = cof[joff + 0][ioff + i];
      temp[1][i] = cof[joff + 1][ioff + i] + temp[0][i];
      temp[2][i] = cof[joff + 2][ioff + i] + temp[1][i];
      temp[3][i] = cof[joff + 3][ioff + i] + temp[2][i];
    }

    for(i=0; i<4; ++i)
    {
      mb_rres[joff    ][ioff + i]=temp[0][i];
      mb_rres[joff + 1][ioff + i]=temp[1][i];
      mb_rres[joff + 2][ioff + i]=temp[2][i];
      mb_rres[joff + 3][ioff + i]=temp[3][i];
    }
  }
  else if(currMB->ipmode_DPCM == HOR_PRED)
  {
    for(j=0; j<4; ++j)
    {
      temp[j][0] = cof[joff + j][ioff    ];
      temp[j][1] = cof[joff + j][ioff + 1] + temp[j][0];
      temp[j][2] = cof[joff + j][ioff + 2] + temp[j][1];
      temp[j][3] = cof[joff + j][ioff + 3] + temp[j][2];
    }

    for(j=0; j<4; ++j)
    {
      mb_rres[joff + j][ioff    ]=temp[j][0];
      mb_rres[joff + j][ioff + 1]=temp[j][1];
      mb_rres[joff + j][ioff + 2]=temp[j][2];
      mb_rres[joff + j][ioff + 3]=temp[j][3];
    }
  }
  else
  {
    for (j = joff; j < joff + BLOCK_SIZE; ++j)
      for (i = ioff; i < ioff + BLOCK_SIZE; ++i)
        mb_rres[j][i] = cof[j][i];
  }

  for (j = joff; j < joff + BLOCK_SIZE; ++j)
  {
    for (i = ioff; i < ioff + BLOCK_SIZE; ++i)
    {
      mb_rec[j][i] = (imgpel) (mb_rres[j][i] + mb_pred[j][i]);
    }
  }
}

/*!
************************************************************************
* \brief
*    Inverse residual DPCM for Intra lossless coding
*
* \par Input:
*    ioff_x,joff_y: Block position inside a macro block (0,8).
************************************************************************
*/
//For residual DPCM
void Inv_Residual_trans_8x8(Macroblock *currMB, ColorPlane pl, int ioff,int joff)
{
  Slice *currSlice = currMB->p_Slice;
  int i, j;
  int temp[8][8];
  imgpel **mb_pred = currSlice->mb_pred[pl];
  imgpel **mb_rec  = currSlice->mb_rec[pl];
  int    **mb_rres = currSlice->mb_rres[pl];

  if(currMB->ipmode_DPCM == VERT_PRED)
  {
    for(i=0; i<8; ++i)
    {
      temp[0][i] = mb_rres[joff + 0][ioff + i];
      temp[1][i] = mb_rres[joff + 1][ioff + i] + temp[0][i];
      temp[2][i] = mb_rres[joff + 2][ioff + i] + temp[1][i];
      temp[3][i] = mb_rres[joff + 3][ioff + i] + temp[2][i];
      temp[4][i] = mb_rres[joff + 4][ioff + i] + temp[3][i];
      temp[5][i] = mb_rres[joff + 5][ioff + i] + temp[4][i];
      temp[6][i] = mb_rres[joff + 6][ioff + i] + temp[5][i];
      temp[7][i] = mb_rres[joff + 7][ioff + i] + temp[6][i];
    }
    for(i=0; i<8; ++i)
    {
      mb_rres[joff  ][ioff+i]=temp[0][i];
      mb_rres[joff+1][ioff+i]=temp[1][i];
      mb_rres[joff+2][ioff+i]=temp[2][i];
      mb_rres[joff+3][ioff+i]=temp[3][i];
      mb_rres[joff+4][ioff+i]=temp[4][i];
      mb_rres[joff+5][ioff+i]=temp[5][i];
      mb_rres[joff+6][ioff+i]=temp[6][i];
      mb_rres[joff+7][ioff+i]=temp[7][i];
    }
  }
  else if(currMB->ipmode_DPCM == HOR_PRED)//HOR_PRED
  {
    for(i=0; i<8; ++i)
    {
      temp[i][0] = mb_rres[joff + i][ioff + 0];
      temp[i][1] = mb_rres[joff + i][ioff + 1] + temp[i][0];
      temp[i][2] = mb_rres[joff + i][ioff + 2] + temp[i][1];
      temp[i][3] = mb_rres[joff + i][ioff + 3] + temp[i][2];
      temp[i][4] = mb_rres[joff + i][ioff + 4] + temp[i][3];
      temp[i][5] = mb_rres[joff + i][ioff + 5] + temp[i][4];
      temp[i][6] = mb_rres[joff + i][ioff + 6] + temp[i][5];
      temp[i][7] = mb_rres[joff + i][ioff + 7] + temp[i][6];
    }
    for(i=0; i<8; ++i)
    {
      mb_rres[joff+i][ioff+0]=temp[i][0];
      mb_rres[joff+i][ioff+1]=temp[i][1];
      mb_rres[joff+i][ioff+2]=temp[i][2];
      mb_rres[joff+i][ioff+3]=temp[i][3];
      mb_rres[joff+i][ioff+4]=temp[i][4];
      mb_rres[joff+i][ioff+5]=temp[i][5];
      mb_rres[joff+i][ioff+6]=temp[i][6];
      mb_rres[joff+i][ioff+7]=temp[i][7];
    }
  }

  for (j = joff; j < joff + BLOCK_SIZE*2; ++j)
  {
    for (i = ioff; i < ioff + BLOCK_SIZE*2; ++i)
    {
      mb_rec [j][i]  = (imgpel) (mb_rres[j][i] + mb_pred[j][i]);
    }
  }
}

/*!
 ***********************************************************************
 * \brief
 *    Luma DC inverse transform
 ***********************************************************************
 */ 
void itrans_2(Macroblock *currMB,    //!< current macroblock
              ColorPlane pl)         //!< used color plane
{
  Slice *currSlice = currMB->p_Slice;
  VideoParameters *p_Vid = currMB->p_Vid;
  int i,j;

  int transform_pl = IS_INDEPENDENT(p_Vid) ? PLANE_Y /*p_Vid->colour_plane_id*/ : pl;
  int **cof = currSlice->cof[transform_pl];
  int qp_scaled = currMB->qp_scaled[pl];

  int qp_per = p_Vid->qp_per_matrix[ qp_scaled ];
  int qp_rem = p_Vid->qp_rem_matrix[ qp_scaled ];      

  int invLevelScale = currSlice->InvLevelScale4x4_Intra[pl][qp_rem][0][0];
  int **M4;
  get_mem2Dint(&M4, BLOCK_SIZE, BLOCK_SIZE);
  
  // horizontal
  for (j=0; j < 4;++j) 
  {
    for (i=0; i < 4;++i) 
    {
      M4[j][i]=cof[j<<2][i<<2];
    }
  }

  ihadamard4x4(M4, M4);

  // vertical
  for (j=0; j < 4;++j) 
  {
    for (i=0; i < 4;++i) 
    {
      cof[j<<2][i<<2] = rshift_rnd((( M4[j][i] * invLevelScale) << qp_per), 6);
    }
  }

  free_mem2Dint(M4);
}


void itrans_sp(Macroblock *currMB,   //!< current macroblock
               ColorPlane pl,        //!< used color plane
               int ioff,             //!< index to 4x4 block
               int joff)             //!< index to 4x4 block
{
  VideoParameters *p_Vid = currMB->p_Vid;
  Slice *currSlice = currMB->p_Slice;
  int i,j;  
  int ilev, icof;

  int qp = (currSlice->slice_type == SI_SLICE) ? currSlice->qs : p_Vid->qp;
  int qp_per = p_Vid->qp_per_matrix[ qp ];
  int qp_rem = p_Vid->qp_rem_matrix[ qp ];

  int qp_per_sp = p_Vid->qp_per_matrix[ currSlice->qs ];
  int qp_rem_sp = p_Vid->qp_rem_matrix[ currSlice->qs ];
  int q_bits_sp = Q_BITS + qp_per_sp;

  imgpel **mb_pred = currSlice->mb_pred[pl];
  imgpel **mb_rec  = currSlice->mb_rec[pl];
  int    **mb_rres = currSlice->mb_rres[pl];
  int    **cof     = currSlice->cof[pl];
  int max_imgpel_value = p_Vid->max_pel_value_comp[pl];

  const int (*InvLevelScale4x4)  [4] = dequant_coef[qp_rem];
  const int (*InvLevelScale4x4SP)[4] = dequant_coef[qp_rem_sp];  
  int **PBlock;  

  get_mem2Dint(&PBlock, MB_BLOCK_SIZE, MB_BLOCK_SIZE);

  for (j=0; j< BLOCK_SIZE; ++j)
    for (i=0; i< BLOCK_SIZE; ++i)
      PBlock[j][i] = mb_pred[j+joff][i+ioff];

  forward4x4(PBlock, PBlock, 0, 0);

  if(p_Vid->sp_switch || currSlice->slice_type==SI_SLICE)
  {    
    for (j=0;j<BLOCK_SIZE;++j)
    {
      for (i=0;i<BLOCK_SIZE;++i)
      {
        // recovering coefficient since they are already dequantized earlier
        icof = (cof[joff + j][ioff + i] >> qp_per) / InvLevelScale4x4[j][i];
        //icof = ((cof[joff + j][ioff + i] * quant_coef[qp_rem][j][i])>> (qp_per + 15)) ;
        // icof  = rshift_rnd_sf(cof[joff + j][ioff + i] * quant_coef[qp_rem][j][i], qp_per + 15);
        ilev  = rshift_rnd_sf(iabs(PBlock[j][i]) * quant_coef[qp_rem_sp][j][i], q_bits_sp);
        ilev  = isignab(ilev, PBlock[j][i]) + icof;
        cof[joff + j][ioff + i] = ilev * InvLevelScale4x4SP[j][i] << qp_per_sp;
      }
    }
  }
  else
  {
    for (j=0;j<BLOCK_SIZE;++j)
    {
      for (i=0;i<BLOCK_SIZE;++i)
      {
        // recovering coefficient since they are already dequantized earlier
        icof = (cof[joff + j][ioff + i] >> qp_per) / InvLevelScale4x4[j][i];
        //icof = cof[joff + j][ioff + i];
        //icof  = rshift_rnd_sf(cof[joff + j][ioff + i] * quant_coef[qp_rem][j][i], qp_per + 15);
        ilev = PBlock[j][i] + ((icof * InvLevelScale4x4[j][i] * A[j][i] <<  qp_per) >> 6);
        ilev  = isign(ilev) * rshift_rnd_sf(iabs(ilev) * quant_coef[qp_rem_sp][j][i], q_bits_sp);
        //cof[joff + j][ioff + i] = ilev * InvLevelScale4x4SP[j][i] << qp_per_sp;
        cof[joff + j][ioff + i] = ilev * InvLevelScale4x4SP[j][i] << qp_per_sp;
      }
    }
  }

  inverse4x4(cof, mb_rres, joff, ioff);

  for (j=joff; j<joff +BLOCK_SIZE;++j)
  {
    for (i=ioff; i < ioff + BLOCK_SIZE;++i)
    {
      mb_rec[j][i] = (imgpel) iClip1(max_imgpel_value,rshift_rnd_sf(mb_rres[j][i], DQ_BITS));
    }
  }  

  free_mem2Dint(PBlock);
}


void itrans_sp_cr(Macroblock *currMB, int uv)
{
  Slice *currSlice = currMB->p_Slice;
  VideoParameters *p_Vid = currMB->p_Vid;
  int i,j,ilev, icof, n2,n1;
  int mp1[BLOCK_SIZE];
  int qp_per,qp_rem;
  int qp_per_sp,qp_rem_sp,q_bits_sp;
  imgpel **mb_pred = currSlice->mb_pred[uv + 1];
  int    **cof = currSlice->cof[uv + 1];
  int **PBlock;  

  get_mem2Dint(&PBlock, MB_BLOCK_SIZE, MB_BLOCK_SIZE);


  qp_per    = p_Vid->qp_per_matrix[ ((p_Vid->qp < 0 ? p_Vid->qp : QP_SCALE_CR[p_Vid->qp]))];
  qp_rem    = p_Vid->qp_rem_matrix[ ((p_Vid->qp < 0 ? p_Vid->qp : QP_SCALE_CR[p_Vid->qp]))];

  qp_per_sp = p_Vid->qp_per_matrix[ ((currSlice->qs < 0 ? currSlice->qs : QP_SCALE_CR[currSlice->qs]))];
  qp_rem_sp = p_Vid->qp_rem_matrix[ ((currSlice->qs < 0 ? currSlice->qs : QP_SCALE_CR[currSlice->qs]))];
  q_bits_sp = Q_BITS + qp_per_sp;  

  if (currSlice->slice_type == SI_SLICE)
  {
    qp_per = qp_per_sp;
    qp_rem = qp_rem_sp;
  }

  for (j=0; j < p_Vid->mb_cr_size_y; ++j)
  {
    for (i=0; i < p_Vid->mb_cr_size_x; ++i)
    {
      PBlock[j][i] = mb_pred[j][i];
      mb_pred[j][i] = 0;
    }
  }

  for (n2=0; n2 < p_Vid->mb_cr_size_y; n2 += BLOCK_SIZE)
  {
    for (n1=0; n1 < p_Vid->mb_cr_size_x; n1 += BLOCK_SIZE)
    {
      forward4x4(PBlock, PBlock, n2, n1);
    }
  }

  //     2X2 transform of DC coeffs.
  mp1[0] = (PBlock[0][0] + PBlock[4][0] + PBlock[0][4] + PBlock[4][4]);
  mp1[1] = (PBlock[0][0] - PBlock[4][0] + PBlock[0][4] - PBlock[4][4]);
  mp1[2] = (PBlock[0][0] + PBlock[4][0] - PBlock[0][4] - PBlock[4][4]);
  mp1[3] = (PBlock[0][0] - PBlock[4][0] - PBlock[0][4] + PBlock[4][4]);

  if (p_Vid->sp_switch || currSlice->slice_type == SI_SLICE)  
  {        
    for (n2=0; n2 < 2; ++n2 )
    {
      for (n1=0; n1 < 2; ++n1 )
      {
        //quantization fo predicted block
        ilev = rshift_rnd_sf(iabs (mp1[n1+n2*2]) * quant_coef[qp_rem_sp][0][0], q_bits_sp + 1);
        //addition
        ilev = isignab(ilev, mp1[n1+n2*2]) + cof[n2<<2][n1<<2];
        //dequantization
        mp1[n1+n2*2] =ilev * dequant_coef[qp_rem_sp][0][0] << qp_per_sp;
      }
    }

    for (n2 = 0; n2 < p_Vid->mb_cr_size_y; n2 += BLOCK_SIZE)
    {
      for (n1 = 0; n1 < p_Vid->mb_cr_size_x; n1 += BLOCK_SIZE)
      {
        for (j = 0; j < BLOCK_SIZE; ++j)
        {
          for (i = 0; i < BLOCK_SIZE; ++i)
          {
            // recovering coefficient since they are already dequantized earlier
            cof[n2 + j][n1 + i] = (cof[n2 + j][n1 + i] >> qp_per) / dequant_coef[qp_rem][j][i];

            //quantization of the predicted block
            ilev = rshift_rnd_sf(iabs(PBlock[n2 + j][n1 + i]) * quant_coef[qp_rem_sp][j][i], q_bits_sp);
            //addition of the residual
            ilev = isignab(ilev,PBlock[n2 + j][n1 + i]) + cof[n2 + j][n1 + i];
            // Inverse quantization
            cof[n2 + j][n1 + i] = ilev * dequant_coef[qp_rem_sp][j][i] << qp_per_sp;
          }
        }
      }
    }
  }
  else
  {
    for (n2=0; n2 < 2; ++n2 )
    {
      for (n1=0; n1 < 2; ++n1 )
      {
        ilev = mp1[n1+n2*2] + (((cof[n2<<2][n1<<2] * dequant_coef[qp_rem][0][0] * A[0][0]) << qp_per) >> 5);
        ilev = isign(ilev) * rshift_rnd_sf(iabs(ilev) * quant_coef[qp_rem_sp][0][0], q_bits_sp + 1);
        //ilev = isignab(rshift_rnd_sf(iabs(ilev)* quant_coef[qp_rem_sp][0][0], q_bits_sp + 1), ilev);
        mp1[n1+n2*2] = ilev * dequant_coef[qp_rem_sp][0][0] << qp_per_sp;
      }
    }

    for (n2 = 0; n2 < p_Vid->mb_cr_size_y; n2 += BLOCK_SIZE)
    {
      for (n1 = 0; n1 < p_Vid->mb_cr_size_x; n1 += BLOCK_SIZE)
      {
        for (j = 0; j< BLOCK_SIZE; ++j)
        {
          for (i = 0; i< BLOCK_SIZE; ++i)
          {
            // recovering coefficient since they are already dequantized earlier
            //icof = ((((cof[n2 + j][n1 + i] << 4) + qp_per/2)>> qp_per) + dequant_coef[qp_rem][j][i]/2) / dequant_coef[qp_rem][j][i];
            icof = (cof[n2 + j][n1 + i] >> qp_per) / dequant_coef[qp_rem][j][i];
            //dequantization and addition of the predicted block      
            ilev = PBlock[n2 + j][n1 + i] + ((icof * dequant_coef[qp_rem][j][i] * A[j][i] << qp_per) >> 6);
            //quantization and dequantization
            ilev = isign(ilev) * rshift_rnd_sf(iabs(ilev) * quant_coef[qp_rem_sp][j][i], q_bits_sp);
            cof[n2 + j][n1 + i] = ilev * dequant_coef[qp_rem_sp][j][i] << qp_per_sp;
            //printf( " %d %d %d\n", j, i, quant_coef[qp_rem_sp][j][i]);
          }
        }
      }
    }
  }

  cof[0][0] = (mp1[0] + mp1[1] + mp1[2] + mp1[3]) >> 1;
  cof[0][4] = (mp1[0] + mp1[1] - mp1[2] - mp1[3]) >> 1;
  cof[4][0] = (mp1[0] - mp1[1] + mp1[2] - mp1[3]) >> 1;
  cof[4][4] = (mp1[0] - mp1[1] - mp1[2] + mp1[3]) >> 1;

  free_mem2Dint(PBlock);
}

void iMBtrans4x4(Macroblock *currMB, ColorPlane pl, int smb)
{
  Slice *currSlice = currMB->p_Slice;
  VideoParameters *p_Vid = currMB->p_Vid;

  StorablePicture *dec_picture = p_Vid->dec_picture;
  int jj, ii;
  int block8x8;
  int k;  

  imgpel **curr_img = pl ? dec_picture->imgUV[pl - 1]: dec_picture->imgY;

  // =============== 4x4 itrans ================
  // -------------------------------------------

  currMB->itrans_4x4 = (smb) ? itrans_sp : ((currMB->is_lossless == FALSE) ? itrans4x4 : Inv_Residual_trans_4x4);

  for (block8x8=0; block8x8 < MB_BLOCK_SIZE; block8x8 += 4)
  { 
    for (k = block8x8; k < block8x8 + 4; ++k )
    {
      jj = ((decode_block_scan[k] >> 2) & 3) << BLOCK_SHIFT;
      ii = (decode_block_scan[k] & 3) << BLOCK_SHIFT;

      currMB->itrans_4x4(currMB, pl, ii, jj);   // use DCT transform and make 4x4 block mb_rres from prediction block mb_pred
    }
  }

  // construct picture from 4x4 blocks
  copy_image_data_16x16(&curr_img[currMB->pix_y], currSlice->mb_rec[pl], currMB->pix_x, 0);
}

void iMBtrans8x8(Macroblock *currMB, ColorPlane pl)
{
  Slice *currSlice = currMB->p_Slice;
  VideoParameters *p_Vid = currMB->p_Vid;
  StorablePicture *dec_picture = p_Vid->dec_picture;
  imgpel **curr_img = pl ? dec_picture->imgUV[pl - 1]: dec_picture->imgY;

  int block8x8;
  int ioff, joff;

  for (block8x8=0; block8x8<4; ++block8x8)
  {
    // =============== 8x8 itrans ================
    // -------------------------------------------
    ioff = 8 * (block8x8 & 0x01);
    joff = 8 * (block8x8 >> 1);

    itrans8x8(currMB, pl, ioff, joff);      // use DCT transform and make 8x8 block mb_rres from prediction block mb_pred
  }

  copy_image_data_16x16(&curr_img[currMB->pix_y], currSlice->mb_rec[pl], currMB->pix_x, 0);
}

void iTransform(Macroblock *currMB, ColorPlane pl, int smb)
{
  Slice *currSlice = currMB->p_Slice;
  VideoParameters *p_Vid = currMB->p_Vid;
  StorablePicture *dec_picture = p_Vid->dec_picture;
  imgpel **curr_img;
  int uv = pl-1; 

  if ((currMB->cbp & 15) != 0 || smb)
  {
    if(currMB->luma_transform_size_8x8_flag == 0) // 4x4 inverse transform
    {
      iMBtrans4x4(currMB, pl, smb); 
    }
    else // 8x8 inverse transform
    {  
      iMBtrans8x8(currMB, pl);    
    }
  }
  else
  {
    curr_img = pl ? dec_picture->imgUV[uv] : dec_picture->imgY;
    copy_image_data_16x16(&curr_img[currMB->pix_y], currSlice->mb_pred[pl], currMB->pix_x, 0);
  }

  if ((dec_picture->chroma_format_idc != YUV400) && (dec_picture->chroma_format_idc != YUV444)) 
  {
    imgpel **curUV;
    int b4, b8;
    int ioff, joff;
    imgpel **mb_rec;

    for(uv=0;uv<2;++uv)
    {
      currMB->itrans_4x4 = (currMB->is_lossless == FALSE) ? itrans4x4 : itrans4x4_ls;

      // =============== 4x4 itrans ================
      // -------------------------------------------
      curUV = &dec_picture->imgUV[uv][currMB->pix_c_y]; 
      mb_rec = currSlice->mb_rec[uv+1];

      if (!smb && (currMB->cbp>>4))
      {
        for (b8 = 0; b8 < (p_Vid->num_uv_blocks); ++b8)
        {
          for(b4 = 0; b4 < 4; ++b4)
          {
            joff = subblk_offset_y[1][b8][b4];
            ioff = subblk_offset_x[1][b8][b4];

            currMB->itrans_4x4(currMB, (ColorPlane) (uv + 1), ioff, joff);
          }
        }

        copy_image_data(curUV, mb_rec, currMB->pix_c_x, 0, p_Vid->mb_size[1][0], p_Vid->mb_size[1][1]);
      }
      else if (smb)
      {
        itrans_sp_cr(currMB, uv);

        for (joff = 0; joff < p_Vid->mb_cr_size_y; joff += BLOCK_SIZE)
        {
          for(ioff = 0; ioff < p_Vid->mb_cr_size_x ;ioff += BLOCK_SIZE)
          {
            currMB->itrans_4x4(currMB, (ColorPlane) (uv + 1), ioff, joff);
          }
        }

        copy_image_data(curUV, mb_rec, currMB->pix_c_x, 0, p_Vid->mb_size[1][0], p_Vid->mb_size[1][1]);
      }
      else 
      {
        copy_image_data(curUV, currSlice->mb_pred[uv + 1], currMB->pix_c_x, 0, p_Vid->mb_size[1][0], p_Vid->mb_size[1][1]);
      }
    }
  }
}

/*!
 *************************************************************************************
 * \brief
 *    Copy ImgPel Data from one structure to another (16x16)
 *************************************************************************************
 */
void copy_image_data_16x16(imgpel  **imgBuf1, imgpel  **imgBuf2, int off1, int off2)
{
  int j;
  for(j=0; j<MB_BLOCK_SIZE; ++j)
  {
    memcpy(&imgBuf1[j][off1], &imgBuf2[j][off2], MB_BLOCK_SIZE * sizeof (imgpel));
  }
}

/*!
 *************************************************************************************
 * \brief
 *    Copy ImgPel Data from one structure to another (8x8)
 *************************************************************************************
 */
void copy_image_data_8x8(imgpel  **imgBuf1, imgpel  **imgBuf2, int off1, int off2)
{
  int j;
  for(j = 0; j < BLOCK_SIZE_8x8; ++j)
  {
    memcpy(&imgBuf1[j][off1], &imgBuf2[j][off2], BLOCK_SIZE_8x8 * sizeof (imgpel));
  }
}


/*!
 *************************************************************************************
 * \brief
 *    Copy ImgPel Data from one structure to another (4x4)
 *************************************************************************************
 */
void copy_image_data_4x4(imgpel  **imgBuf1, imgpel  **imgBuf2, int off1, int off2)
{
  int j;
  for(j = 0; j < BLOCK_SIZE; ++j)
  {
    memcpy(&imgBuf1[j][off1], &imgBuf2[j][off2], BLOCK_SIZE * sizeof (imgpel));
  }
}

/*!
 *************************************************************************************
 * \brief
 *    Copy ImgPel Data from one structure to another (8x8)
 *************************************************************************************
 */
void copy_image_data(imgpel  **imgBuf1, imgpel  **imgBuf2, int off1, int off2, int width, int height)
{
  int j;
  for(j = 0; j < height; ++j)
  {
    memcpy(&imgBuf1[j][off1], &imgBuf2[j][off2], width * sizeof (imgpel));
  }
}
