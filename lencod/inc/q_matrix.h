
/*!
 ***************************************************************************
 * \file
 *    q_matrix.h
 *
 * \brief
 *    Headerfile for q_matrix array
 *
 * \date
 *    07. Apr 2004
 ***************************************************************************
 */

#ifndef _Q_MATRIX_H_
#define _Q_MATRIX_H_

static const char MatrixType4x4[6][20] =
{
  "INTRA4X4_LUMA",
  "INTRA4X4_CHROMAU",
  "INTRA4X4_CHROMAV",
  "INTER4X4_LUMA",
  "INTER4X4_CHROMAU",
  "INTER4X4_CHROMAV"
};

static const char MatrixType8x8[2][20] =
{
  "INTRA8X8_LUMA",
  "INTER8X8_LUMA",
};

int LevelScale4x4Luma_Intra[6][4][4];
int LevelScale4x4Chroma_Intra[2][6][4][4];

int LevelScale4x4Luma_Inter[6][4][4];
int LevelScale4x4Chroma_Inter[2][6][4][4];

int LevelScale8x8Luma_Intra[6][8][8];

int LevelScale8x8Luma_Inter[6][8][8];

int InvLevelScale4x4Luma_Intra[6][4][4];
int InvLevelScale4x4Chroma_Intra[2][6][4][4];

int InvLevelScale4x4Luma_Inter[6][4][4];
int InvLevelScale4x4Chroma_Inter[2][6][4][4];

int InvLevelScale8x8Luma_Intra[6][8][8];

int InvLevelScale8x8Luma_Inter[6][8][8];

short ScalingList4x4input[6][16];
short ScalingList8x8input[2][64];
short ScalingList4x4[6][16];
short ScalingList8x8[2][64];

short UseDefaultScalingMatrix4x4Flag[6];
short UseDefaultScalingMatrix8x8Flag[2];

static const byte ZZ_SCAN[16]  =
{  0,  1,  4,  8,  5,  2,  3,  6,  9, 12, 13, 10,  7, 11, 14, 15
};

static const byte ZZ_SCAN8[64] =
{  0,  1,  8, 16,  9,  2,  3, 10, 17, 24, 32, 25, 18, 11,  4,  5,
   12, 19, 26, 33, 40, 48, 41, 34, 27, 20, 13,  6,  7, 14, 21, 28,
   35, 42, 49, 56, 57, 50, 43, 36, 29, 22, 15, 23, 30, 37, 44, 51,
   58, 59, 52, 45, 38, 31, 39, 46, 53, 60, 61, 54, 47, 55, 62, 63
};

static const short Quant_intra_default[16] =
{
 6,13,20,28,
13,20,28,32,
20,28,32,37,
28,32,37,42
};

static const short Quant_inter_default[16] =
{
10,14,20,24,
14,20,24,27,
20,24,27,30,
24,27,30,34
};

static const short Quant8_intra_default[64] =
{
 6,10,13,16,18,23,25,27,
10,11,16,18,23,25,27,29,
13,16,18,23,25,27,29,31,
16,18,23,25,27,29,31,33,
18,23,25,27,29,31,33,36,
23,25,27,29,31,33,36,38,
25,27,29,31,33,36,38,40,
27,29,31,33,36,38,40,42
};

static const short Quant8_inter_default[64] =
{
 9,13,15,17,19,21,22,24,
13,13,17,19,21,22,24,25,
15,17,19,21,22,24,25,27,
17,19,21,22,24,25,27,28,
19,21,22,24,25,27,28,30,
21,22,24,25,27,28,30,32,
22,24,25,27,28,30,32,33,
24,25,27,28,30,32,33,35
};

void Init_QMatrix (void);

#endif
