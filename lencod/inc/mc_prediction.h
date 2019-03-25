
/*!
 ************************************************************************
 * \file
 *    mc_prediction.h
 *
 * \brief
 *    motion compensation header
 *
 * \author
 *      Main contributors (see contributors.h for copyright, 
 *                         address and affiliation details)
 *      - Alexis Michael Tourapis  <alexismt@ieee.org>
 *
 *************************************************************************************
 */

#ifndef _MC_PREDICTION_H_
#define _MC_PREDICTION_H_
#include "mbuffer.h"

extern void luma_prediction       ( Macroblock* currMB, int, int, int, int, int, int[2], char *, short );
extern void luma_prediction_bi    ( Macroblock* currMB, int, int, int, int, int, int, short, short, int );
extern void chroma_prediction     ( Macroblock* currMB, int, int, int, int, int, int, int, int, short, short, short );
extern void chroma_prediction_4x4 ( Macroblock* currMB, int, int, int, int, int, int, short, short, short);   

extern void OneComponentChromaPrediction4x4_regenerate (Macroblock *currMB, imgpel* , int , int , short*** , StorablePicture *listX, int );
extern void OneComponentChromaPrediction4x4_retrieve   (Macroblock *currMB, imgpel* , int , int , short*** , StorablePicture *listX, int );

extern void intra_chroma_prediction (Macroblock *currMB, int*, int*, int*);

#endif

