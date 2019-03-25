
/*!
 ************************************************************************
 * \file
 *     me_fullfast.h
 *
 * \author
 *    Alexis Michael Tourapis        <alexis.tourapis@dolby.com>
 *
 * \date
 *    9 September 2006
 *
 * \brief
 *    Headerfile for Fast Full Search motion estimation
 **************************************************************************
 */


#ifndef _ME_FULLFAST_H_
#define _ME_FULLFAST_H_
int FastFullPelBlockMotionSearch (Macroblock *currMB, imgpel* orig_pic, short ref, int list, 
                                  char ***refPic, short ****tmp_mv, int pic_pix_x, int pic_pix_y, int blocktype, MotionVector *pred_mv, MotionVector *mv,
                              int search_range,  int min_mcost, int lambda_factor, int apply_weights);
void InitializeFastFullIntegerSearch (void);
void ResetFastFullIntegerSearch (void);
void ClearFastFullIntegerSearch (void);

#endif

