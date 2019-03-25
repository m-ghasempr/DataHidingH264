
/*!
 ************************************************************************
 * \file conformance.h
 *
 * \brief
 *   Level & Profile Related definitions  
 *
 * \author
 *    Alexis Michael Tourapis         <alexismt@ieee.org>       \n
 *
 ************************************************************************
 */

#ifndef _CONFORMANCE_H_
#define _CONFORMANCE_H_

void ProfileCheck(void);
void LevelCheck(void);

void update_mv_limits(ImageParameters *img, byte is_field);
void clip_mv_range(ImageParameters *img, int search_range, MotionVector *mv, int res);
int  out_of_bounds_mvs(ImageParameters *img, short mv[2]);
void test_clip_mvs(ImageParameters *img, short mv[2], Boolean write_mb);

int InvalidWeightsForBiPrediction(Block8x8Info* b8x8info, int mode);
int InvalidMotionVectors(Block8x8Info* b8x8info, int mode);

#endif

