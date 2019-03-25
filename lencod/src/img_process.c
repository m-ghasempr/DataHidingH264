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
#include "image.h"
#include "img_process.h"


static inline void CPImage(ImageData *imgOut, ImageData *imgIn)
{
  memcpy(imgOut->frm_data[0][0], imgIn->frm_data[0][0], img->height * img->width * sizeof (imgpel));

  if (img->yuv_format != YUV400)
  {
    memcpy(imgOut->frm_data[1][0], imgIn->frm_data[1][0], img->height_cr * img->width_cr * sizeof (imgpel));
    memcpy(imgOut->frm_data[2][0], imgIn->frm_data[2][0], img->height_cr * img->width_cr * sizeof (imgpel));
  }
}

// to be modified
static inline void FilterImage(ImageData *imgOut, ImageData *imgIn)
{
  memcpy(imgOut->frm_data[0][0], imgIn->frm_data[0][0], img->height * img->width * sizeof (imgpel));

  if (img->yuv_format != YUV400)
  {
    memcpy(imgOut->frm_data[1][0], imgIn->frm_data[1][0], img->height_cr * img->width_cr * sizeof (imgpel));
    memcpy(imgOut->frm_data[2][0], imgIn->frm_data[2][0], img->height_cr * img->width_cr * sizeof (imgpel));
  }
}

static inline void YV12toYUV(ImageData *imgOut, ImageData *imgIn)
{
  memcpy(imgOut->frm_data[0][0], imgIn->frm_data[0][0], img->height * img->width * sizeof (imgpel));

  if (img->yuv_format != YUV400)
  {
    memcpy(imgOut->frm_data[1][0], imgIn->frm_data[2][0], img->height_cr * img->width_cr * sizeof (imgpel));
    memcpy(imgOut->frm_data[2][0], imgIn->frm_data[1][0], img->height_cr * img->width_cr * sizeof (imgpel));
  }
}

void ProcessImage( InputParameters *params)
{
  switch( params->ProcessInput  )
  {
  default:
  case 0:
    CPImage(&imgData, &imgData1);
    break;
  case 1:
    FilterImage(&imgData, &imgData1);
    break;
  case 2:
    YV12toYUV(&imgData, &imgData1);
    break;
  }
}


