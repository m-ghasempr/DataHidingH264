
/*!
 *************************************************************************************
 * \file img_distortion.c
 *
 * \brief
 *    YUV to RGB color conversion
 *
 * \author
 *    Main contributors (see contributors.h for copyright, address and affiliation details)
 *     - Woo-Shik Kim                    <wooshik.kim@usc.edu>
 *************************************************************************************
 */
#include "contributors.h"
#include "global.h"
#include "memalloc.h"
#include "img_distortion.h"

#define YUV2RGB_YOFFSET

//YUV to RGB conversion
#ifdef YUV2RGB_YOFFSET
#define OFFSET_Y 16
static const float K0 = 1.164f;
static const float K1 = 1.596f;
static const float K2 = 0.391f;
static const float K3 = 0.813f;
static const float K4 = 2.018f;
#else
static const float K0 = 1.000f;
static const float K1 = 1.402f;
static const float K2 = 0.34414f;
static const float K3 = 0.71414f;
static const float K4 = 1.772f;
#endif

int create_RGB_memory(ImageParameters *p_Img)
{
  int memory_size = 0;
  int j;
  for( j = 0; j < 3; j++ )
  {
    memory_size += get_mem2Dpel (&p_Img->imgRGB_src.data[j], p_Img->height, p_Img->width);
  }
  for( j = 0; j < 3; j++ )
  {
    memory_size += get_mem2Dpel (&p_Img->imgRGB_ref.data[j], p_Img->height, p_Img->width);
  }
  
  return memory_size;
}

void delete_RGB_memory(ImageParameters *p_Img)
{
  int i;
  for( i = 0; i < 3; i++ )
  {
    free_mem2Dpel(p_Img->imgRGB_src.data[i]);
  }
  for( i = 0; i < 3; i++ )
  {
    free_mem2Dpel(p_Img->imgRGB_ref.data[i]);
  }
}

void init_YUVtoRGB(ImageParameters *p_Img, InputParameters *p_Inp)
{ 
  float conv_scale = (float) (65536.0f);

  p_Img->wka0 = float2int(  conv_scale * K0);
  p_Img->wka1 = float2int(  conv_scale * K1);
  p_Img->wka2 = float2int( -conv_scale * K2);
  p_Img->wka3 = float2int( -conv_scale * K3);
  p_Img->wka4 = float2int(  conv_scale * K4);

#ifdef YUV2RGB_YOFFSET
  p_Img->offset_y = OFFSET_Y << (p_Inp->output.bit_depth[0] - 8);
  p_Img->offset_cr = 1 << (p_Inp->output.bit_depth[0] - 1);
#endif
}

/*! 
*************************************************************************************
* \brief
*    YUV to RGB conversion
*    ITU 601 with and without offset consideration
*    Upsampling by repetition of chroma samples in case of 4:2:0 and 4:2:2
*    Method not support for 4:0:0 content
*************************************************************************************
*/
void YUVtoRGB(ImageParameters *p_Img, ImageStructure *YUV, ImageStructure *RGB)
{
  int i, j, j_cr, i_cr;
  int sy, su, sv;
  int wbuv, wguv, wruv;
  imgpel *Y, *U, *V, *R, *G, *B;
  FrameFormat format = YUV->format;
  int width = format.width;
  int height = format.height;
  int max_value = format.max_value[0];

  // Color conversion
  for (j = 0; j < height; j++)
  {
    j_cr = j >> p_Img->shift_cr_y;
    Y = YUV->data[0][j];
    U = YUV->data[1][j_cr];
    V = YUV->data[2][j_cr];
    R = RGB->data[0][j];
    G = RGB->data[1][j];
    B = RGB->data[2][j];

    for (i = 0; i < width; i++)
    {
      i_cr = i >> p_Img->shift_cr_x;

      su = U[i_cr] - p_Img->offset_cr;
      sv = V[i_cr] - p_Img->offset_cr;

      wruv = p_Img->wka1 * sv;
      wguv = p_Img->wka2 * su + p_Img->wka3 * sv;
      wbuv = p_Img->wka4 * su;

#ifdef YUV2RGB_YOFFSET // Y offset value of 16 is considered
      sy = p_Img->wka0 * (Y[i] - p_Img->offset_y);
#else
      sy = p_Img->wka0 * Y[i];
#endif

      R[i] = (imgpel) iClip1( max_value, rshift_rnd(sy + wruv, 16));
      G[i] = (imgpel) iClip1( max_value, rshift_rnd(sy + wguv, 16));
      B[i] = (imgpel) iClip1( max_value, rshift_rnd(sy + wbuv, 16));
    }
  }
  // Setting RGB FrameFormat
  RGB->format = format;  // copy format information from YUV to RGB
  RGB->format.yuv_format  = YUV444;
  RGB->format.color_model = CM_RGB;
  RGB->format.height_cr   = format.height;
  RGB->format.width_cr    = format.width;
  for (i = 1; i < 3; i++)
  {
    RGB->format.size_cmp[i]     = format.size_cmp[0];
    RGB->format.bit_depth[i]    = format.bit_depth[0];
    RGB->format.max_value[i]    = max_value;
    RGB->format.max_value_sq[i] = format.max_value_sq[0];
  }
}
