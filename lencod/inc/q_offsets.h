
/*!
 ***************************************************************************
 * \file
 *    q_offsets.h
 *
 * \brief
 *    Headerfile for q_offsets array
 *
 * \date
 *    18. Nov 2004
 ***************************************************************************
 */

#ifndef _Q_OFFSETS_H_
#define _Q_OFFSETS_H_

static const char OffsetType4x4[9][24] =
{
  "INTRA4X4_LUMA_INTRA",
  "INTRA4X4_CHROMAU_INTRA",
  "INTRA4X4_CHROMAV_INTRA",
  "INTRA4X4_LUMA_INTER",
  "INTRA4X4_CHROMAU_INTER",
  "INTRA4X4_CHROMAV_INTER",
  "INTER4X4_LUMA",
  "INTER4X4_CHROMAU",
  "INTER4X4_CHROMAV"
};

static const char OffsetType8x8[3][24] =
{
  "INTRA8X8_LUMA_INTRA",
  "INTRA8X8_LUMA_INTER",
  "INTER8X8_LUMA",
};


int LevelOffset4x4Luma_Intra[13][4][4];
int LevelOffset4x4Chroma_Intra[2][13][4][4];

int LevelOffset4x4Luma_Inter[13][4][4];
int LevelOffset4x4Chroma_Inter[2][13][4][4];

int LevelOffset8x8Luma_Intra[13][8][8];

int LevelOffset8x8Luma_Inter[13][8][8];

short OffsetList4x4input[9][16];
short OffsetList8x8input[3][64];
short OffsetList4x4[9][16];
short OffsetList8x8[3][64];

//short UseDefaultOffsetMatrix4x4Flag[6];
//short UseDefaultOffsetMatrix8x8Flag[2];

static const short Offset_intra_default_intra[16] =
{
  341,341,341,341,
  341,341,341,341,
  341,341,341,341,
  341,341,341,341
};

static const short Offset_intra_default_inter[16] =
{
  171,171,171,171,
  171,171,171,171,
  171,171,171,171,
  171,171,171,171,
};

static const short Offset_inter_default[16] =
{
  171,171,171,171,
  171,171,171,171,
  171,171,171,171,
  171,171,171,171,
};

static const short Offset8_intra_default_intra[64] =
{
  341,341,341,341,341,341,341,341,
  341,341,341,341,341,341,341,341,
  341,341,341,341,341,341,341,341,
  341,341,341,341,341,341,341,341,
  341,341,341,341,341,341,341,341,
  341,341,341,341,341,341,341,341,
  341,341,341,341,341,341,341,341,
  341,341,341,341,341,341,341,341
};

static const short Offset8_intra_default_inter[64] =
{
  171,171,171,171,171,171,171,171,
  171,171,171,171,171,171,171,171,
  171,171,171,171,171,171,171,171,
  171,171,171,171,171,171,171,171,
  171,171,171,171,171,171,171,171,
  171,171,171,171,171,171,171,171,
  171,171,171,171,171,171,171,171,
  171,171,171,171,171,171,171,171
};

static const short Offset8_inter_default[64] =
{
  171,171,171,171,171,171,171,171,
  171,171,171,171,171,171,171,171,
  171,171,171,171,171,171,171,171,
  171,171,171,171,171,171,171,171,
  171,171,171,171,171,171,171,171,
  171,171,171,171,171,171,171,171,
  171,171,171,171,171,171,171,171,
  171,171,171,171,171,171,171,171
};

void Init_QOffsetMatrix (void);

#endif
