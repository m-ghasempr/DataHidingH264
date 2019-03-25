
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

extern void    ProfileCheck         (InputParameters *p_Inp);
extern void    LevelCheck           (ImageParameters *p_Img, InputParameters *p_Inp);
extern void    update_mv_limits     (ImageParameters *p_Img, InputParameters *p_Inp, byte is_field);
extern void    clip_mv_range        (ImageParameters *p_Img, int search_range, MotionVector *mv, int res);
extern int     out_of_bounds_mvs    (ImageParameters *p_Img, short mv[2]);
extern void    test_clip_mvs        (ImageParameters *p_Img, short mv[2], Boolean write_mb);
extern Boolean CheckPredictionParams(Macroblock  *currMB, Block8x8Info *b8x8info, int mode);

extern unsigned int getMaxMBPS(unsigned int levelIdc);
extern unsigned int getMinCR  (unsigned int levelIdc);
extern unsigned int getMaxBR  (unsigned int levelIdc);
extern unsigned int getMaxCPB (unsigned int levelIdc);

#endif

