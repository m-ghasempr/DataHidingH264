
/*!
 ************************************************************************
 * \file refbuf.c
 *
 * \brief
 *    Declarations of teh reference frame buffer types and functions
 ************************************************************************
 */



#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <assert.h>

#include "refbuf.h"

#define CACHELINESIZE 32


/*!
 ************************************************************************
 * \brief
 *    Reference buffer write routines
 ************************************************************************
 */
void PutPel_14 (pel_t **Pic, int y, int x, pel_t val)
{
  Pic [IMG_PAD_SIZE*4+y][IMG_PAD_SIZE*4+x] = val;
}

void PutPel_11 (pel_t *Pic, int y, int x, pel_t val)
{
  Pic [y*img->width+x] = val;
}


/*!
 ************************************************************************
 * \brief
 *    Reference buffer read, Full pel
 ************************************************************************
 */
pel_t FastPelY_11 (pel_t *Pic, int y, int x)
{
  return Pic [y*img->width+x];
}


pel_t *FastLine16Y_11 (pel_t *Pic, int y, int x)
{
  return &Pic [y*img->width+x];
}

pel_t *FastLineX (int dummy, pel_t* Pic, int y, int x)
{
  return Pic + y*img->width + x;
}

pel_t UMVPelY_11 (pel_t *Pic, int y, int x)
{
  if (x < 0)
  {
    if (y < 0)
      return Pic [0];
    if (y >= img->height)
      return Pic [(img->height-1) * img->width];
    return Pic [y*img->width];
  }

  if (x >= img->width)
  {
    if (y < 0)
      return Pic [img->width-1];
    if (y >= img->height)
      return Pic [img->height * img->width -1];
    return Pic [(y+1)*img->width -1 ];
  }

  if (y < 0)    // note: corner pixels were already processed
    return Pic [x];
  if (y >= img->height)
    return Pic [(img->height-1)*img->width+x];

  return Pic [y*img->width+x];
}

/*!
 ************************************************************************
 * \note
 *    The following function is NOT reentrant!  Use a buffer
 *    provided by the caller to change that (but it costs a memcpy()...
 ************************************************************************
 */
static pel_t line[16];

pel_t *UMVLine16Y_11 (pel_t *Pic, int y, int x)
{
  int i, maxx;
  pel_t *Picy;

  Picy = &Pic [max(0,min(img->height-1,y)) * img->width];

  if (x < 0) {                    // Left edge ?

    maxx = min(0,x+16);
    for (i = x; i < maxx; i++)
      line[i-x] = Picy [0];       // Replicate left edge pixel

    maxx = x+16;
    for (i = 0; i < maxx; i++)    // Copy non-edge pixels
      line[i-x] = Picy [i];
  }
  else if (x > img->width-16)  {  // Right edge ?

    maxx = img->width;
    for (i = x; i < maxx; i++)
      line[i-x] = Picy [i];       // Copy non-edge pixels

    maxx = x+16;
    for (i = max(img->width,x); i < maxx; i++)
      line[i-x] = Picy [img->width-1];  // Replicate right edge pixel
  }
  else                            // No edge
    return &Picy [x];

  return line;
}


pel_t *UMVLineX (int size, pel_t* Pic, int y, int x)
{
  int i, maxx;
  pel_t *Picy;

  Picy = Pic + max(0,min(img->height-1,y)) * img->width;

  if (x < 0)                            // Left edge
  {
    maxx = min(0,x+size);
    for (i = x; i < maxx; i++)
    {
      line[i-x] = Picy [0];             // Replicate left edge pixel
    }
    maxx = x+size;
    for (i = 0; i < maxx; i++)          // Copy non-edge pixels
      line[i-x] = Picy [i];
  }
  else if (x > img->width-size)         // Right edge
  {
    maxx = img->width;
    for (i = x; i < maxx; i++)
    {
      line[i-x] = Picy [i];             // Copy non-edge pixels
    }
    maxx = x+size;
    for (i = max(img->width,x); i < maxx; i++)
    {
      line[i-x] = Picy [img->width-1];  // Replicate right edge pixel
    }
  }
  else                                  // No edge
  {
    return Picy + x;
  }

  return line;
}

/*!
 ************************************************************************
 * \brief
 *    Reference buffer read, 1/2 pel
 ************************************************************************
 */
pel_t FastPelY_12 (pel_t **Pic, int y, int x)
{
  return Pic [IMG_PAD_SIZE*4+(y<<1)][IMG_PAD_SIZE*4+(x<<1)];
}


pel_t UMVPelY_12 (pel_t **Pic, int y, int x)
{
  return UMVPelY_14 (Pic, y*2, x*2);
}


/*!
 ************************************************************************
 * \brief
 *    Reference buffer, 1/4 pel
 ************************************************************************
 */
pel_t UMVPelY_14 (pel_t **Pic, int y, int x)
{
  int width4  = ((img->width+2*IMG_PAD_SIZE-1)<<2);
  int height4 = ((img->height+2*IMG_PAD_SIZE-1)<<2);

  x = x + IMG_PAD_SIZE*4;
  y = y + IMG_PAD_SIZE*4;

  if (x < 0)
  {
    if (y < 0)
      return Pic [y&3][x&3];
    if (y > height4)
      return Pic [height4+(y&3)][x&3];
    return Pic [y][x&3];
  }

  if (x > width4)
  {
    if (y < 0)
      return Pic [y&3][width4+(x&3)];
    if (y > height4)
      return Pic [height4+(y&3)][width4+(x&3)];
    return Pic [y][width4+(x&3)];
  }

  if (y < 0)    // note: corner pixels were already processed
    return Pic [y&3][x];
  if (y > height4)
    return Pic [height4+(y&3)][x];

  return Pic [y][x];
}

pel_t FastPelY_14 (pel_t **Pic, int y, int x)
{
  return Pic [IMG_PAD_SIZE*4+y][IMG_PAD_SIZE*4+x];
}


