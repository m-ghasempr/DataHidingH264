
/*!
 ************************************************************************
 * \file
 *     me_fullsearch.h
 *
 * \author
 *    Alexis Michael Tourapis        <alexis.tourapis@dolby.com>
 *
 * \date
 *    9 September 2006
 *
 * \brief
 *    Headerfile for Full Search motion estimation
 **************************************************************************
 */


#ifndef _ME_FULLSEARCH_H_
#define _ME_FULLSEARCH_H_
extern distblk FullPelBlockMotionSearch (Macroblock *currMB, MotionVector *pred_mv, MEBlock *mv_block,
                                     distblk min_mcost, int lambda_factor);
extern distblk FullPelBlockMotionBiPred (Macroblock *currMB, int list, 
                                     MotionVector *pred_mv1, MotionVector *pred_mv2, MotionVector *mv1, MotionVector *mv2, MEBlock *,
                                     int search_range, distblk min_mcost, int lambda_factor);
extern distblk SubPelBlockMotionSearch  (Macroblock *currMB, MotionVector *pred_mv, MEBlock *mv_block, 
                                     distblk min_mcost, int* lambda_factor);
extern distblk SubPelBlockSearchBiPred  (Macroblock *currMB, MEBlock *mv_block, int list, 
                                     MotionVector *pred_mv1, MotionVector *pred_mv2, MotionVector *mv1, MotionVector *mv2, 
                                     distblk min_mcost, int* lambda_factor);
#endif

