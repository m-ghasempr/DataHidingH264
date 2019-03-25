
/*!
 ************************************************************************
 * \file refbuf.h
 *
 * \brief
 *    Declarations of the reference frame buffer types and functions
 ************************************************************************
 */
#ifndef _REBUF_H_
#define _REBUF_H_

// global sub-pel image access variables
extern int height_pad, width_pad;
extern int height_pad_cr, width_pad_cr;

/*!
 ************************************************************************
 * \brief
 *    Yields a pel line _pointer_ from one of the 16 sub-images
 *    Input does not require subpixel image indices
 ************************************************************************
 */
static inline imgpel *UMVLine4X (imgpel ****Pic, int y, int x, int height, int width)
{
  return &(Pic[(y & 0x03)][(x & 0x03)][iClip3( 0, height, y >> 2)][iClip3( 0, width, x >> 2)]);
}

/*!
 ************************************************************************
 * \brief
 *    Yields a pel line _pointer_ from one of the 16 sub-images
 *    Input does not require subpixel image indices
 ************************************************************************
 */
static inline imgpel *FastLine4X (imgpel ****Pic, int y, int x)
{
  return &(Pic[(y & 0x03)][(x & 0x03)][y >> 2][x >> 2]);
}

/*!
 ************************************************************************
 * \brief
 *    Yields a pel line _pointer_ from one of the 16 (4:4:4), 32 (4:2:2),
 *    or 64 (4:2:0) sub-images
 *    Input does not require subpixel image indices
 ************************************************************************
 */
static inline imgpel *UMVLine8X_chroma (imgpel ****Pic, int y, int x, int height, int width)
{
  return &(Pic[y & chroma_mask_mv_y][x & chroma_mask_mv_x][iClip3 (0, height, y >> chroma_shift_y)][iClip3 (0, width , x >> chroma_shift_x)]);
}

/*!
 ************************************************************************
 * \brief
 *    Yields a pel line _pointer_ from one of the 16 (4:4:4), 32 (4:2:2),
 *    or 64 (4:2:0) sub-images
 *    Input does not require subpixel image indices
 ************************************************************************
 */
static inline imgpel *FastLine8X_chroma (imgpel ****Pic, int y, int x)
{
  return &(Pic[y & chroma_mask_mv_y][x & chroma_mask_mv_x][y >> chroma_shift_y][x >> chroma_shift_x]);
}


#endif

