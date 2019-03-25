
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

pel_t UMVPelY_14 (pel_t **Pic, int y, int x, int height, int width);
pel_t FastPelY_14 (pel_t **Pic, int y, int x, int height, int width);

pel_t* UMVLineX  (int, pel_t*, int, int, int, int);
pel_t* UMVLineX2 (int, pel_t*, int, int, int, int);
pel_t* FastLineX (int, pel_t*, int, int, int, int);

pel_t FastPelY14  (pel_t **Pic, int y, int x, int height, int width);
pel_t *UMVLine4X  (pel_t **Pic, int y, int x, int height4, int width4);
pel_t *UMVLine4X2 (pel_t **Pic, int y, int x, int height4, int width4);

pel_t *FastLine4X (pel_t **Pic, int y, int x, int height, int width);

pel_t *UMVLine16Y_11 (pel_t *Pic, int y, int x, int height, int width);
pel_t *FastLine16Y_11 (pel_t *Pic, int y, int x, int height, int width);

void PutPel_14 (pel_t **Pic, int y, int x, pel_t val);
void PutPel_11 (pel_t *Pic, int pel_pos, pel_t val);

#endif

