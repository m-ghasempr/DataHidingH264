
/*!
 ***************************************************************************
 * \file
 *    cabac.h
 *
 * \brief
 *    Header file for entropy coding routines
 *
 * \author
 *    Detlev Marpe                                                         \n
 *    Copyright (C) 2000 HEINRICH HERTZ INSTITUTE All Rights Reserved.
 *
 * \date
 *    21. Oct 2000 (Changes by Tobias Oelbaum 28.08.2001)
 ***************************************************************************
 */

#ifndef _CABAC_H_
#define _CABAC_H_

#include "global.h"

extern MotionInfoContexts*  create_contexts_MotionInfo(void);
extern TextureInfoContexts* create_contexts_TextureInfo(void);
extern void delete_contexts_MotionInfo(MotionInfoContexts *enco_ctx);
extern void delete_contexts_TextureInfo(TextureInfoContexts *enco_ctx);

extern void cabac_new_slice(void);

extern void readMB_typeInfo_CABAC           (Macroblock *currMB, SyntaxElement *se, ImageParameters *img, DecodingEnvironmentPtr dep_dp);
extern void readB8_typeInfo_CABAC           (Macroblock *currMB, SyntaxElement *se, ImageParameters *img, DecodingEnvironmentPtr dep_dp);
extern void readIntraPredMode_CABAC         (Macroblock *currMB, SyntaxElement *se, ImageParameters *img, DecodingEnvironmentPtr dep_dp);
extern void readRefFrame_CABAC              (Macroblock *currMB, SyntaxElement *se, ImageParameters *img, DecodingEnvironmentPtr dep_dp);
extern void readMVD_CABAC                   (Macroblock *currMB, SyntaxElement *se, ImageParameters *img, DecodingEnvironmentPtr dep_dp);
extern void readCBP_CABAC                   (Macroblock *currMB, SyntaxElement *se, ImageParameters *img, DecodingEnvironmentPtr dep_dp);
extern void readRunLevel_CABAC              (Macroblock *currMB, SyntaxElement *se, ImageParameters *img, DecodingEnvironmentPtr dep_dp);
extern void readDquant_CABAC                (Macroblock *currMB, SyntaxElement *se, ImageParameters *img, DecodingEnvironmentPtr dep_dp);
extern void readCIPredMode_CABAC            (Macroblock *currMB, SyntaxElement *se, ImageParameters *img, DecodingEnvironmentPtr dep_dp);
extern void readMB_skip_flagInfo_CABAC      (Macroblock *currMB, SyntaxElement *se, ImageParameters *img, DecodingEnvironmentPtr dep_dp);
extern void readFieldModeInfo_CABAC         (Macroblock *currMB, SyntaxElement *se, ImageParameters *img, DecodingEnvironmentPtr dep_dp);
extern void readMB_transform_size_flag_CABAC(Macroblock *currMB, SyntaxElement *se, ImageParameters *img, DecodingEnvironmentPtr dep_dp);

extern void readIPCM_CABAC(struct datapartition *dP);

extern int  cabac_startcode_follows(Slice *currSlice, int eos_bit);

extern int  readSyntaxElement_CABAC         (SyntaxElement *se, ImageParameters *img, DataPartition *this_dataPart);

extern int check_next_mb_and_get_field_mode_CABAC( SyntaxElement *se, ImageParameters *img, DataPartition  *act_dp);

extern void CheckAvailabilityOfNeighborsCABAC(Macroblock *currMB);


#endif  // _CABAC_H_

