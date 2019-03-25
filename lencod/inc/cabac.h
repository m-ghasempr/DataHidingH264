
/*!
 ***************************************************************************
 * \file
 *    cabac.h
 *
 * \brief
 *    Headerfile for entropy coding routines
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

extern const byte maxpos       [];
extern const byte c1isdc       [];
extern const byte type2ctx_bcbp[];
extern const byte type2ctx_map [];
extern const byte type2ctx_last[];
extern const byte type2ctx_one [];
extern const byte type2ctx_abs [];
extern const byte max_c2       [];

extern const byte  pos2ctx_map8x8  [];
extern const byte  pos2ctx_map8x4  [];
extern const byte  pos2ctx_map4x4  [];
extern const byte  pos2ctx_map2x4c [];
extern const byte  pos2ctx_map4x4c [];
extern const byte* pos2ctx_map     [];
extern const byte  pos2ctx_map8x8i [];
extern const byte  pos2ctx_map8x4i [];
extern const byte  pos2ctx_map4x8i [];
extern const byte* pos2ctx_map_int [];
extern const byte  pos2ctx_last8x8 [];
extern const byte  pos2ctx_last8x4 [];
extern const byte  pos2ctx_last4x4 [];
extern const byte  pos2ctx_last2x4c[];
extern const byte  pos2ctx_last4x4c[];
extern const byte* pos2ctx_last    [];

// CABAC
MotionInfoContexts* create_contexts_MotionInfo(void);
TextureInfoContexts* create_contexts_TextureInfo(void);
void delete_contexts_MotionInfo(MotionInfoContexts *enco_ctx);
void delete_contexts_TextureInfo(TextureInfoContexts *enco_ctx);
void writeHeaderToBuffer(void);
void writeMB_typeInfo_CABAC(SyntaxElement *se, DataPartition *dp);
void writeIntraPredMode_CABAC(SyntaxElement *se, DataPartition *dp);
void writeB8_typeInfo_CABAC(SyntaxElement *se, DataPartition *dp);
void writeRefFrame_CABAC(SyntaxElement *se, DataPartition *dp);
void writeMVD_CABAC(SyntaxElement *se, DataPartition *dp);
void writeCBP_CABAC(Macroblock *currMB, SyntaxElement *se, DataPartition *dp);
void writeDquant_CABAC(SyntaxElement *se, DataPartition *dp);
void writeRunLevel_CABAC(Macroblock* currMB, SyntaxElement *se, DataPartition *dp);
void writeCIPredMode_CABAC(SyntaxElement *se, DataPartition *dp);
void print_ctx_TextureInfo(TextureInfoContexts *enco_ctx);
void writeMB_skip_flagInfo_CABAC(Macroblock *currMB, SyntaxElement *se, DataPartition *dp);
void writeFieldModeInfo_CABAC(SyntaxElement *se, DataPartition *dp); //GB
void writeCBP_BIT_CABAC (Macroblock* currMB, int b8, int bit, int cbp, int inter, EncodingEnvironmentPtr eep_dp, TextureInfoContexts *ctx);
void cabac_new_slice(void);
void CheckAvailabilityOfNeighborsCABAC(Macroblock* currMB);

void writeMB_transform_size_CABAC(SyntaxElement *se, DataPartition *dp);


#endif  // CABAC_H

