/*!
 ***************************************************************************
 * \file
 *    me_distortion.h
 *
 * \author
 *    Alexis Michael Tourapis        <alexis.tourapis@dolby.com>
 *    Athanasios Leontaris           <aleon@dolby.com>
 *
 * \date
 *    11. August 2006
 *
 * \brief
 *    Headerfile for motion estimation distortion
 **************************************************************************
 */

#ifndef _ME_DISTORTION_H_
#define _ME_DISTORTION_H_

extern int HadamardSAD4x4(int* diff);
extern int HadamardSAD8x8(int* diff);
// SAD functions
extern int computeSAD         (StorablePicture *ref1, MEBlock*, int, MotionVector *);
extern int computeSAD4x4      (StorablePicture *ref1, MEBlock*, int, MotionVector *);
extern int computeSAD4x8      (StorablePicture *ref1, MEBlock*, int, MotionVector *);
extern int computeSAD8x4      (StorablePicture *ref1, MEBlock*, int, MotionVector *);
extern int computeSAD8x8      (StorablePicture *ref1, MEBlock*, int, MotionVector *);
extern int computeSAD8x16     (StorablePicture *ref1, MEBlock*, int, MotionVector *);
extern int computeSAD16x8     (StorablePicture *ref1, MEBlock*, int, MotionVector *);
extern int computeSAD16x16    (StorablePicture *ref1, MEBlock*, int, MotionVector *);
// Weighted Prediction SAD functions
extern int computeSADWP       (StorablePicture *ref1, MEBlock*, int, MotionVector *);
extern int computeSADWP4x4    (StorablePicture *ref1, MEBlock*, int, MotionVector *);
extern int computeSADWP4x8    (StorablePicture *ref1, MEBlock*, int, MotionVector *);
extern int computeSADWP8x4    (StorablePicture *ref1, MEBlock*, int, MotionVector *);
extern int computeSADWP8x8    (StorablePicture *ref1, MEBlock*, int, MotionVector *);
extern int computeSADWP8x16   (StorablePicture *ref1, MEBlock*, int, MotionVector *);
extern int computeSADWP16x8   (StorablePicture *ref1, MEBlock*, int, MotionVector *);
extern int computeSADWP16x16  (StorablePicture *ref1, MEBlock*, int, MotionVector *);
// SATD
extern int computeSATD        (StorablePicture *ref1, MEBlock*, int, MotionVector *);
extern int computeSAT4x4D     (StorablePicture *ref1, MEBlock*, int, MotionVector *);
extern int computeSAT8x8D     (StorablePicture *ref1, MEBlock*, int, MotionVector *);
extern int computeSAT8x8D8x16 (StorablePicture *ref1, MEBlock*, int, MotionVector *);
extern int computeSAT4x4D8x16 (StorablePicture *ref1, MEBlock*, int, MotionVector *);
extern int computeSAT8x8D16x8 (StorablePicture *ref1, MEBlock*, int, MotionVector *);
extern int computeSAT4x4D16x8 (StorablePicture *ref1, MEBlock*, int, MotionVector *);
extern int computeSAT8x8D16x16(StorablePicture *ref1, MEBlock*, int, MotionVector *);
extern int computeSAT4x4D16x16(StorablePicture *ref1, MEBlock*, int, MotionVector *);
extern int computeSATDWP      (StorablePicture *ref1, MEBlock*, int, MotionVector *);
// SSE
extern int computeSSE         (StorablePicture *ref1, MEBlock*, int, MotionVector *);
extern int computeSSEWP       (StorablePicture *ref1, MEBlock*, int, MotionVector *);

// Bipred SAD
extern int computeBiPredSAD1      (StorablePicture *ref1, StorablePicture *ref2, MEBlock*, int, MotionVector *, MotionVector *);
extern int computeBiPred4x4SAD1   (StorablePicture *ref1, StorablePicture *ref2, MEBlock*, int, MotionVector *, MotionVector *);
extern int computeBiPred4x8SAD1   (StorablePicture *ref1, StorablePicture *ref2, MEBlock*, int, MotionVector *, MotionVector *);
extern int computeBiPred8x4SAD1   (StorablePicture *ref1, StorablePicture *ref2, MEBlock*, int, MotionVector *, MotionVector *);
extern int computeBiPred8x8SAD1   (StorablePicture *ref1, StorablePicture *ref2, MEBlock*, int, MotionVector *, MotionVector *);
extern int computeBiPred8x16SAD1  (StorablePicture *ref1, StorablePicture *ref2, MEBlock*, int, MotionVector *, MotionVector *);
extern int computeBiPred16x8SAD1  (StorablePicture *ref1, StorablePicture *ref2, MEBlock*, int, MotionVector *, MotionVector *);
extern int computeBiPred16x16SAD1 (StorablePicture *ref1, StorablePicture *ref2, MEBlock*, int, MotionVector *, MotionVector *);
extern int computeBiPredSAD2      (StorablePicture *ref1, StorablePicture *ref2, MEBlock*, int, MotionVector *, MotionVector *);

// Bipred SATD
extern int computeBiPredSATD1     (StorablePicture *ref1, StorablePicture *ref2, MEBlock*, int, MotionVector *, MotionVector *);
extern int computeBiPred8x8SATD1  (StorablePicture *ref1, StorablePicture *ref2, MEBlock*, int, MotionVector *, MotionVector *);
extern int computeBiPred8x16SATD1 (StorablePicture *ref1, StorablePicture *ref2, MEBlock*, int, MotionVector *, MotionVector *);
extern int computeBiPred16x8SATD1 (StorablePicture *ref1, StorablePicture *ref2, MEBlock*, int, MotionVector *, MotionVector *);
extern int computeBiPred16x16SATD1(StorablePicture *ref1, StorablePicture *ref2, MEBlock*, int, MotionVector *, MotionVector *);
extern int computeBiPredSATD2     (StorablePicture *ref1, StorablePicture *ref2, MEBlock*, int, MotionVector *, MotionVector *);

// Bipred SSE
extern int computeBiPredSSE1      (StorablePicture *ref1, StorablePicture *ref2, MEBlock*, int, MotionVector *, MotionVector *);
extern int computeBiPredSSE2      (StorablePicture *ref1, StorablePicture *ref2, MEBlock*, int, MotionVector *, MotionVector *);


extern void select_distortion   (ImageParameters *p_Img, InputParameters *p_Inp);

#endif
