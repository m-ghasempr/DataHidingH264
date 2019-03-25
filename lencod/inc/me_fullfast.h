
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
typedef struct me_full_fast
{
  int          **search_setup_done;     //!< flag if all block SAD's have been calculated yet
  MotionVector **search_center; //!< absolute search center for fast full motion search
  MotionVector **search_center_padded; //!< absolute search center for fast full motion search
  int          **pos_00;             //!< position of (0,0) vector
  distpel   *****BlockSAD;        //!< SAD for all blocksize, ref. frames and motion vectors
  int          **max_search_range;
} MEFullFast;

extern distblk FastFullPelBlockMotionSearch (Macroblock *currMB, MotionVector *pred_mv, MEBlock *mv_block, distblk min_mcost, int lambda_factor);
extern void InitializeFastFullIntegerSearch (VideoParameters *p_Vid, InputParameters *p_Inp);
extern void ResetFastFullIntegerSearch      (VideoParameters *p_Vid);
extern void ClearFastFullIntegerSearch      (VideoParameters *p_Vid);


#endif

