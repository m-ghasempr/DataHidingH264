/*!
*************************************************************************************
* \file img_process.c
*
* \brief
*    Input data Image Processing functions
*
* \author
*    Main contributors (see contributors.h for copyright, address and affiliation details)
*      - Alexis Michael Tourapis <alexis.tourapis@dolby.com>
*
*************************************************************************************
*/

#include "contributors.h"
#include "global.h"
#include "img_process.h"
#include "io_image.h"
#include "memalloc.h"


static inline void ResetImage(ImageData *imgOut)
{
  memset(imgOut->frm_data[0][0], 0, imgOut->format.height * imgOut->format.width * sizeof (imgpel));

  if (imgOut->format.yuv_format != YUV400)
  {
    if (sizeof(imgpel) == sizeof(char))
    {
      memset(imgOut->frm_data[1][0], 128, imgOut->format.height_cr * imgOut->format.width_cr * sizeof (imgpel));
      memset(imgOut->frm_data[2][0], 128, imgOut->format.height_cr * imgOut->format.width_cr * sizeof (imgpel));
    }
    else
    {
      int i, j, k;
      imgpel med_value;
      for (k = 1; k <=2; k++)
      {
        med_value = (imgpel) (imgOut->format.max_value[k] + 1) >> 1;
        for (j = 0; j < imgOut->format.height_cr; j++)
        {
          for (i = 0; i < imgOut->format.width_cr; i++)
          {
            imgOut->frm_data[k][j][i] = med_value;
          }
        }
      }
    }
  }
}


static inline void CPImage(ImageData *imgOut, ImageData *imgIn)
{
  memcpy(imgOut->frm_data[0][0], imgIn->frm_data[0][0], imgIn->format.height * imgIn->format.width * sizeof (imgpel));

  if (imgIn->format.yuv_format != YUV400)
  {
    memcpy(imgOut->frm_data[1][0], imgIn->frm_data[1][0], imgIn->format.height_cr * imgIn->format.width_cr * sizeof (imgpel));
    memcpy(imgOut->frm_data[2][0], imgIn->frm_data[2][0], imgIn->format.height_cr * imgIn->format.width_cr * sizeof (imgpel));
  }
}

// to be modified
static inline void FilterImage(ImageData *imgOut, ImageData *imgIn)
{
  memcpy(imgOut->frm_data[0][0], imgIn->frm_data[0][0], imgIn->format.height * imgIn->format.width * sizeof (imgpel));

  if (imgIn->format.yuv_format != YUV400)
  {
    memcpy(imgOut->frm_data[1][0], imgIn->frm_data[1][0], imgIn->format.height_cr * imgIn->format.width_cr * sizeof (imgpel));
    memcpy(imgOut->frm_data[2][0], imgIn->frm_data[2][0], imgIn->format.height_cr * imgIn->format.width_cr * sizeof (imgpel));
  }
}

// to be modified
static inline void FilterImageSep(ImageData *imgOut, ImageData *imgIn)
{
  int i, j;
  static const int SepFilter[6] = {1, -5, 20, 20, -5, 1};
  int max_width  = imgOut->format.width - 1;
  int max_height = imgOut->format.height - 1;

  int **temp_data; // temp memory for filtering. Could be allocated once to speed up code
  
  get_mem2Dint(&temp_data, imgIn->format.height, imgIn->format.width);

  // implementation was not optimized. only just implemented as proof of concept
  // horizontal filtering
  for (j = 0; j < imgOut->format.height; j++)
  {
    for (i = 0; i < imgOut->format.width; i++)
    {
      temp_data[j][i] = 
        SepFilter[0] * imgIn->frm_data[0][j][iClip3(0, max_width, i - 2)] +
        SepFilter[1] * imgIn->frm_data[0][j][iClip3(0, max_width, i - 1)] +
        SepFilter[2] * imgIn->frm_data[0][j][iClip3(0, max_width, i    )] +
        SepFilter[3] * imgIn->frm_data[0][j][iClip3(0, max_width, i + 1)] +
        SepFilter[4] * imgIn->frm_data[0][j][iClip3(0, max_width, i + 2)] +
        SepFilter[5] * imgIn->frm_data[0][j][iClip3(0, max_width, i + 3)];
    }
  }

  for (j = 0; j < imgOut->format.height; j++)
  {
    for (i = 0; i < imgOut->format.width; i++)
    {
      imgOut->frm_data[0][j][i] = (imgpel) iClip3(0, imgOut->format.max_value[0], rshift_rnd_sign(
        SepFilter[0] * temp_data[iClip3(0, max_height, j - 2)][i] +
        SepFilter[1] * temp_data[iClip3(0, max_height, j - 1)][i] +
        SepFilter[2] * temp_data[iClip3(0, max_height, j    )][i] +
        SepFilter[3] * temp_data[iClip3(0, max_height, j + 1)][i] +
        SepFilter[4] * temp_data[iClip3(0, max_height, j + 2)][i] +
        SepFilter[5] * temp_data[iClip3(0, max_height, j + 3)][i], 10));
    }
  }

  if (imgOut->format.yuv_format != YUV400)
  {
    int k;
    max_width  = imgOut->format.width_cr - 1;
    max_height = imgOut->format.height_cr - 1;

    for (k = 1; k <=2; k++)
    {
      // horizontal filtering
      for (j = 0; j < imgOut->format.height_cr; j++)
      {
        for (i = 0; i < imgOut->format.width_cr; i++)
        {
          temp_data[j][i] = 
            SepFilter[0] * imgIn->frm_data[k][j][iClip3(0, max_width, i - 2)] +
            SepFilter[1] * imgIn->frm_data[k][j][iClip3(0, max_width, i - 1)] +
            SepFilter[2] * imgIn->frm_data[k][j][iClip3(0, max_width, i    )] +
            SepFilter[3] * imgIn->frm_data[k][j][iClip3(0, max_width, i + 1)] +
            SepFilter[4] * imgIn->frm_data[k][j][iClip3(0, max_width, i + 2)] +
            SepFilter[5] * imgIn->frm_data[k][j][iClip3(0, max_width, i + 3)];
        }
      }

      for (j = 0; j < imgOut->format.height_cr; j++)
      {
        for (i = 0; i < imgOut->format.width_cr; i++)
        {
          imgOut->frm_data[k][j][i] = (imgpel) iClip3(0, imgOut->format.max_value[k], rshift_rnd_sign(
            SepFilter[0] * temp_data[iClip3(0, max_height, j - 2)][i] +
            SepFilter[1] * temp_data[iClip3(0, max_height, j - 1)][i] +
            SepFilter[2] * temp_data[iClip3(0, max_height, j    )][i] +
            SepFilter[3] * temp_data[iClip3(0, max_height, j + 1)][i] +
            SepFilter[4] * temp_data[iClip3(0, max_height, j + 2)][i] +
            SepFilter[5] * temp_data[iClip3(0, max_height, j + 3)][i], 10));
        }
      }
    }
  }

  free_mem2Dint(temp_data);
}


// to be modified
static inline void MuxImages(ImageData *imgOut, ImageData *imgIn0, ImageData *imgIn1, ImageData *Map)
{
  int i, j;
  for (j = 0; j < imgOut->format.height; j++)
  {
    for (i = 0; i < imgOut->format.width; i++)
    {
      imgOut->frm_data[0][j][i] = (imgpel) rshift_rnd_sf(imgIn0->frm_data[0][j][i] * (Map->format.max_value[0] - Map->frm_data[0][j][i]) + imgIn1->frm_data[0][j][i] * Map->frm_data[0][j][i], Map->format.bit_depth[0]);
    }
  }
  
  if (imgOut->format.yuv_format != YUV400)
  {
    int k;
    for (k = 1; k <=2; k++)
    {
      for (j = 0; j < imgOut->format.height_cr; j++)
      {
        for (i = 0; i < imgOut->format.width_cr; i++)
        {
          imgOut->frm_data[k][j][i] = (imgpel) rshift_rnd_sf(imgIn0->frm_data[k][j][i] * (Map->format.max_value[k] - Map->frm_data[k][j][i]) + imgIn1->frm_data[k][j][i] * Map->frm_data[k][j][i], Map->format.bit_depth[k]);
        }
      }
    }
  }
}

static inline void YV12toYUV(ImageData *imgOut, ImageData *imgIn)
{
  memcpy(imgOut->frm_data[0][0], imgIn->frm_data[0][0], imgIn->format.height * imgIn->format.width * sizeof (imgpel));

  if (imgIn->format.yuv_format != YUV400)
  {
    memcpy(imgOut->frm_data[1][0], imgIn->frm_data[2][0], imgIn->format.height_cr * imgIn->format.width_cr * sizeof (imgpel));
    memcpy(imgOut->frm_data[2][0], imgIn->frm_data[1][0], imgIn->format.height_cr * imgIn->format.width_cr * sizeof (imgpel));
  }
}

int InitProcessImage( VideoParameters *p_Vid, InputParameters *p_Inp)
{
  int memory_size = 0;
  switch( p_Inp->ProcessInput )
  {
  default:
    break;
  }
  return memory_size;
}

void ClearProcessImage( VideoParameters *p_Vid, InputParameters *p_Inp)
{
  switch( p_Inp->ProcessInput )
  {
  default:
    break;
  }   
}

void ProcessImage( VideoParameters *p_Vid, InputParameters *p_Inp )
{
  switch( p_Inp->ProcessInput )
  {
  default:
  case 0:
    CPImage(&p_Vid->imgData, &p_Vid->imgData0);
    break;
  case 1:
    FilterImage(&p_Vid->imgData, &p_Vid->imgData0);
    break;
  case 2:
    YV12toYUV(&p_Vid->imgData, &p_Vid->imgData0);
    break;
  case 3:
    MuxImages(&p_Vid->imgData, &p_Vid->imgData0, &p_Vid->imgData1, &p_Vid->imgData2);
    break;
  case 4:
    FilterImageSep(&p_Vid->imgData, &p_Vid->imgData0);
    break;
  }
}



