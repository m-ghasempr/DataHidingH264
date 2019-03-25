
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
extern int FullPelBlockMotionSearch (Macroblock *currMB, imgpel* orig_pic, short ref, int list, 
                                     char ***refPic, short ****tmp_mv, int pic_pix_x, int pic_pix_y,
                                     int blocktype, MotionVector *pred_mv, MotionVector *mv, 
                                     int search_range,  int min_mcost, int lambda_factor, int apply_weights);
extern int FullPelBlockMotionBiPred (Macroblock *currMB, imgpel* orig_pic, short ref, int list, 
                                     char ***refPic, short ****tmp_mv, int pic_pix_x, int pic_pix_y,
                                     int blocktype, MotionVector *pred_mv1, MotionVector *pred_mv2, MotionVector *mv1, MotionVector *mv2,
                                     int search_range, int min_mcost, int iteration_no, int lambda_factor, int apply_weights);
extern int SubPelBlockMotionSearch  (imgpel* orig_pic, short ref, int list, int list_offset, int pic_pix_x, int pic_pix_y,
                                     int blocktype, MotionVector *pred_mv, MotionVector *mv,
                                     int search_pos2, int search_pos4, int min_mcost, int* lambda_factor, int apply_weights);
extern int SubPelBlockSearchBiPred  (imgpel* orig_pic, short ref, int list, int pic_pix_x, int pic_pix_y,
                                     int blocktype, MotionVector *pred_mv1, MotionVector *pred_mv2, MotionVector *mv1, MotionVector *mv2, 
                                     int search_pos2, int search_pos4, int min_mcost, int* lambda_factor, int apply_weights);

#endif

