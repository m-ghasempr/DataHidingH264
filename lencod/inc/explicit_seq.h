
/*!
 *************************************************************************************
 * \file explicit_seq.h
 *
 * \brief
 *    Functions for explicit sequence support
 *
 * \author
 *     Main contributors (see contributors.h for copyright, address and affiliation details)
 *     - Alexis Michael Tourapis          <alexismt@ieee.org>
 *************************************************************************************
 */

#ifndef _EXPLICIT_SEQ_H_
#define _EXPLICIT_SEQ_H_

// Definition of structures used for explicit sequence representation
typedef struct
{
  int seq_number;
  int slice_type;
  int is_idr;
  int reference_idc;
  int frame_qp;
  int deblocking;
  int is_field;
} ExpFrameInfo;

typedef struct
{
  int no_frames;
  ExpFrameInfo *info;
} ExpSeqInfo;

extern FILE       *expSFile;
extern ExpSeqInfo *expSeq;

extern void ReadExplicitSeqFile    (ExpSeqInfo *seq_info, int coding_index);
extern void OpenExplicitSeqFile    (InputParameters *pparams);
extern void CloseExplicitSeqFile   (void);
extern void ExplicitUpdateImgParams(ExpFrameInfo *info, ImageParameters *p_img);
#endif
