
/*!
 ************************************************************************
 * \file mv-search.h
 *
 * \brief
 *   array definition for motion search
 *
 * \author
 *    Inge Lille-Langoy               <inge.lille-langoy@telenor.com>   \n
 *    Alexis Michael Tourapis         <alexis.tourapis@dolby.com>       \n
 *
 ************************************************************************
 */

#ifndef _MV_SEARCH_H_
#define _MV_SEARCH_H_

typedef struct {
  int    min_x;
  int    max_x;
  int    min_y;
  int    max_y;
} SearchWindow;


extern int* mvbits;

extern int *byte_abs;
extern SearchWindow searchrange;

extern void get_neighbors(Macroblock *currMB, 
                          PixelPos *block_a, PixelPos *block_b, PixelPos *block_c, PixelPos *block_d,
                          int mb_x, int mb_y, int blockshape_x);

extern void set_access_method(int *access_method, MotionVector *blk, int min_x, int min_y, int max_x, int max_y);

extern int (*BiPredME)      (Macroblock *, imgpel *, short, int, char  ***, short  ****,
                       int, int, int, MotionVector *, MotionVector *, MotionVector *, MotionVector *, int, int, int, int, int);

extern int (*SubPelBiPredME)(imgpel* orig_pic, short ref, int list, int pic_pix_x, int pic_pix_y,
                             int blocktype, MotionVector *pred_mv1, MotionVector *pred_mv2, MotionVector *mv1, MotionVector *mv2, 
                             int search_pos2, int search_pos4, int min_mcost, int* lambda_factor, int apply_weights);
extern int (*SubPelME)      (imgpel* orig_pic, short ref, int list, int list_offset, int pic_pix_x, int pic_pix_y, 
                             int blocktype, MotionVector *pred_mv, MotionVector *mv, 
                             int search_pos2, int search_pos4, int min_mcost, int* lambda_factor, int apply_weights);

extern void SetMotionVectorPredictor (Macroblock *currMB, short  pmv[2], char   **refPic, short  ***tmp_mv,
                               short  ref_frame, int list, int block_x, int block_y, int blockshape_x, int blockshape_y);
extern void GetMotionVectorPredictor (Macroblock *currMB, PixelPos *block_a, PixelPos *block_b, PixelPos *block_c, 
                               short  pmv[2], char **refPic, short ***tmp_mv, short ref_frame, int mb_x, int mb_y, int bx, int by);

extern void PrepareMEParams      (int apply_weights, int ChromaMEEnable, int list, int ref);
extern void PrepareBiPredMEParams(int apply_weights, int ChromaMEEnable, int list, int list_offset, int ref);

extern void Init_Motion_Search_Module (void);
extern void Clear_Motion_Search_Module (void);

#endif

