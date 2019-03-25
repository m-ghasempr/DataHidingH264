/*!
 *************************************************************************************
 * \file intra16x16_pred.c
 *
 * \brief
 *    Functions for intra 8x8 prediction
 *
 * \author
 *      Main contributors (see contributors.h for copyright, 
 *                         address and affiliation details)
 *      - Yuri Vatis
 *      - Jan Muenster
 *      - Alexis Michael Tourapis  <alexismt@ieee.org>
 *
 *************************************************************************************
 */
#include "global.h"
#include "intra16x16_pred.h"
#include "mb_access.h"
#include "image.h"

/*!
 ***********************************************************************
 * \brief
 *    makes and returns 16x16 DC prediction mode
 *
 * \return
 *    DECODING_OK   decoding of intraprediction mode was sucessfull            \n
 *
 ***********************************************************************
 */
static inline int intra16x16_dc_pred(ImageParameters *img,  //!< image parameters 
                                   Macroblock *currMB, 
                                    ColorPlane pl)
{
  int s0 = 0, s1 = 0, s2 = 0;

  int i,j;

  imgpel **imgY = (pl) ? dec_picture->imgUV[pl - 1] : dec_picture->imgY;
  imgpel **mb_pred = &(img->mb_pred[pl][0]); 

  PixelPos up;          //!< pixel position p(0,-1)
  PixelPos left[17];    //!< pixel positions p(-1, -1..15)

  int up_avail, left_avail, left_up_avail;

  s1=s2=0;

  for (i=0;i<17;i++)
  {
    getNeighbour(currMB, -1,  i-1, img->mb_size[IS_LUMA], &left[i]);
  }
  getNeighbour(currMB,    0,   -1, img->mb_size[IS_LUMA], &up);

  if (!active_pps->constrained_intra_pred_flag)
  {
    up_avail      = up.available;
    left_avail    = left[1].available;
    left_up_avail = left[0].available;
  }
  else
  {
    up_avail      = up.available ? img->intra_block[up.mb_addr] : 0;
    for (i = 1, left_avail = 1; i < 17; i++)
      left_avail  &= left[i].available ? img->intra_block[left[i].mb_addr]: 0;
    left_up_avail = left[0].available ? img->intra_block[left[0].mb_addr]: 0;
  }

  for (i = 0; i < MB_BLOCK_SIZE; i++)
  {
    if (up_avail)
      s1 += imgY[up.pos_y][up.pos_x+i];    // sum hor pix
    if (left_avail)
      s2 += imgY[left[i + 1].pos_y][left[i + 1].pos_x];    // sum vert pix
  }
  if (up_avail && left_avail)
    s0 = (s1 + s2 + 16)>>5;       // no edge
  else if (!up_avail && left_avail)
    s0 = (s2 + 8)>>4;              // upper edge
  else if (up_avail && !left_avail)
    s0 = (s1 + 8)>>4;              // left edge
  else
    s0 = img->dc_pred_value_comp[pl];                            // top left corner, nothing to predict from

  for(j = 0; j < MB_BLOCK_SIZE; j++)
  {
    for(i = 0; i < MB_BLOCK_SIZE; i++)
    {
      mb_pred[j][i]=(imgpel) s0;
    }
  }

  return DECODING_OK;
}


/*!
 ***********************************************************************
 * \brief
 *    makes and returns 16x16 vertical prediction mode
 *
 * \return
 *    DECODING_OK   decoding of intraprediction mode was sucessfull            \n
 *
 ***********************************************************************
 */
static inline int intra16x16_vert_pred(ImageParameters *img,  //!< image parameters 
                                       Macroblock *currMB, 
                                       ColorPlane pl)
{
  int j;

  imgpel **imgY = (pl) ? dec_picture->imgUV[pl - 1] : dec_picture->imgY;

  PixelPos up;          //!< pixel position p(0,-1)

  int up_avail;

  getNeighbour(currMB,    0,   -1, img->mb_size[IS_LUMA], &up);

  if (!active_pps->constrained_intra_pred_flag)
  {
    up_avail = up.available;
  }
  else
  {
    up_avail = up.available ? img->intra_block[up.mb_addr] : 0;
  }

  if (!up_avail)
    error ("invalid 16x16 intra pred Mode VERT_PRED_16",500);

  for(j=0;j<MB_BLOCK_SIZE;j++)
    memcpy(img->mb_pred[pl][j], &(imgY[up.pos_y][up.pos_x]), MB_BLOCK_SIZE * sizeof(imgpel));

  return DECODING_OK;
}


/*!
 ***********************************************************************
 * \brief
 *    makes and returns 16x16 horizontal prediction mode
 *
 * \return
 *    DECODING_OK   decoding of intraprediction mode was sucessfull            \n
 *
 ***********************************************************************
 */
static inline int intra16x16_hor_pred(ImageParameters *img,  //!< image parameters 
                                    Macroblock *currMB, 
                                    ColorPlane pl)
{
  int i,j;

  imgpel **imgY = (pl) ? dec_picture->imgUV[pl - 1] : dec_picture->imgY;
  imgpel **mb_pred = &(img->mb_pred[pl][0]); 
  imgpel prediction;

  PixelPos left[17];    //!< pixel positions p(-1, -1..15)

  int left_avail, left_up_avail;

  for (i=0;i<17;i++)
  {
    getNeighbour(currMB, -1,  i-1, img->mb_size[IS_LUMA], &left[i]);
  }

  if (!active_pps->constrained_intra_pred_flag)
  {
    left_avail    = left[1].available;
    left_up_avail = left[0].available;
  }
  else
  {
    for (i = 1, left_avail = 1; i < 17; i++)
      left_avail  &= left[i].available ? img->intra_block[left[i].mb_addr]: 0;
    left_up_avail = left[0].available ? img->intra_block[left[0].mb_addr]: 0;
  }

  if (!left_avail)
    error ("invalid 16x16 intra pred Mode HOR_PRED_16",500);

  for(j = 0; j < MB_BLOCK_SIZE; j++)
  {
    prediction = imgY[left[j+1].pos_y][left[j+1].pos_x];
    for(i = 0; i < MB_BLOCK_SIZE; i++)
      mb_pred[j][i]= prediction; // store predicted 16x16 block
  }

  return DECODING_OK;
}

/*!
 ***********************************************************************
 * \brief
 *    makes and returns 16x16 horizontal prediction mode
 *
 * \return
 *    DECODING_OK   decoding of intraprediction mode was sucessfull            \n
 *
 ***********************************************************************
 */
static inline int intra16x16_plane_pred(ImageParameters *img,  //!< image parameters 
                                    Macroblock *currMB, 
                                    ColorPlane pl)
{
  int i,j;

  int ih = 0, iv = 0;
  int ib,ic,iaa;

  imgpel **imgY = (pl) ? dec_picture->imgUV[pl - 1] : dec_picture->imgY;
  imgpel **mb_pred = &(img->mb_pred[pl][0]); 
  imgpel *mpr_line;
  int max_imgpel_value = img->max_imgpel_value_comp[pl];

  PixelPos up;          //!< pixel position p(0,-1)
  PixelPos left[17];    //!< pixel positions p(-1, -1..15)

  int up_avail, left_avail, left_up_avail;

  for (i=0;i<17;i++)
  {
    getNeighbour(currMB, -1,  i-1, img->mb_size[IS_LUMA], &left[i]);
  }
  getNeighbour(currMB,    0,   -1, img->mb_size[IS_LUMA], &up);

  if (!active_pps->constrained_intra_pred_flag)
  {
    up_avail      = up.available;
    left_avail    = left[1].available;
    left_up_avail = left[0].available;
  }
  else
  {
    up_avail      = up.available ? img->intra_block[up.mb_addr] : 0;
    for (i = 1, left_avail = 1; i < 17; i++)
      left_avail  &= left[i].available ? img->intra_block[left[i].mb_addr]: 0;
    left_up_avail = left[0].available ? img->intra_block[left[0].mb_addr]: 0;
  }

  if (!up_avail || !left_up_avail  || !left_avail)
    error ("invalid 16x16 intra pred Mode PLANE_16",500);

  mpr_line = &imgY[up.pos_y][up.pos_x+7];
  for (i = 1; i < 8; i++)
  {
    ih += i*(mpr_line[i] - mpr_line[-i]);
    iv += i*(imgY[left[8+i].pos_y][left[8+i].pos_x] - imgY[left[8-i].pos_y][left[8-i].pos_x]);
  }

  ih += 8*(mpr_line[8] - imgY[left[0].pos_y][left[0].pos_x]);
  iv += 8*(imgY[left[16].pos_y][left[16].pos_x] - imgY[left[0].pos_y][left[0].pos_x]);

  ib=(5 * ih + 32)>>6;
  ic=(5 * iv + 32)>>6;

  iaa=16 * (mpr_line[8] + imgY[left[16].pos_y][left[16].pos_x]);
  for (j = 0;j < MB_BLOCK_SIZE; j++)
  {
    for (i = 0;i < MB_BLOCK_SIZE; i++)
    {
      mb_pred[j][i] = (imgpel) iClip1(max_imgpel_value, ((iaa + (i - 7) * ib + (j - 7) * ic + 16) >> 5));
    }
  }// store plane prediction

  return DECODING_OK;
}

/*!
 ***********************************************************************
 * \brief
 *    makes and returns 16x16 intra prediction blocks 
 *
 * \return
 *    DECODING_OK   decoding of intraprediction mode was sucessfull            \n
 *    SEARCH_SYNC   search next sync element as errors while decoding occured
 ***********************************************************************
 */
int intrapred16x16(ImageParameters *img, //!< image parameters
                         Macroblock *currMB,  //!< Current Macroblock
                         ColorPlane pl,       //!< Current colorplane (for 4:4:4)                         
                         int predmode)        //!< prediction mode
{
  switch (predmode)
  {
  case VERT_PRED_16:                       // vertical prediction from block above
    return (intra16x16_vert_pred(img, currMB, pl));
    break;
  case HOR_PRED_16:                        // horizontal prediction from left block
    return (intra16x16_hor_pred(img, currMB, pl));
    break;
  case DC_PRED_16:                         // DC prediction
    return (intra16x16_dc_pred(img, currMB, pl));
    break;
  case PLANE_16:// 16 bit integer plan pred
    return (intra16x16_plane_pred(img, currMB, pl));
    break;
  default:
    {                                    // indication of fault in bitstream,exit
      printf("illegal 16x16 intra prediction mode input: %d\n",predmode);
      return SEARCH_SYNC;
    }
  }

  return DECODING_OK;
}

