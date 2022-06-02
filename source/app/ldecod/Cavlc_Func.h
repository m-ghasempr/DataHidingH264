#pragma once
#include "Enc_Entropy.h"

int is_intra(Enc_Macroblock *curr_MB);
int symbol2vlc(Enc_SyntaxElement *sym);
void  writeUVLC2buffer(Enc_SyntaxElement *se, Enc_Bitstream *currStream);
int writeSyntaxElement_Run(Enc_SyntaxElement *se, Enc_DataPartition *dp);
int writeSyntaxElement_TotalZeros(Enc_SyntaxElement *se, Enc_DataPartition *dp);
int writeSyntaxElement_TotalZerosChromaDC(Enc_VideoParameters *p_Vid, Enc_SyntaxElement *se, Enc_DataPartition *dp);
int writeSyntaxElement_Level_VLC1(Enc_SyntaxElement *se, Enc_DataPartition *dp, int profile_idc);
int writeSyntaxElement_Level_VLCN(Enc_SyntaxElement *se, int vlc, Enc_DataPartition *dp, int profile_idc);
int writeSyntaxElement_VLC(Enc_SyntaxElement *se, Enc_DataPartition *dp);
int writeSyntaxElement_NumCoeffTrailingOnesChromaDC(Enc_VideoParameters *p_Vid, Enc_SyntaxElement *se, Enc_DataPartition *dp);
int writeSyntaxElement_NumCoeffTrailingOnes(Enc_SyntaxElement *se, Enc_DataPartition *dp);
void get4x4NeighbourV(Enc_Macroblock *currMB, int block_x, int block_y, int mb_size[2], PixelPos *pix);
int predict_nnz_chromaI(Enc_Macroblock *currMB, int i, int j);
int predict_nnzI(Enc_Macroblock *currMB, int block_type, int i, int j, VideoParameters *Vid);
int writeCoeff4x4_CAVLC_normal(Enc_Macroblock* currMB, int block_type, int b8, int b4, int param, int *lv, int *rn, VideoParameters *Vid);
void reset_mb_nz_coeff(Enc_VideoParameters *p_Vid, int mb_number);
void init_Data(char bufB, int off, Macroblock *currMB, Slice *currSlice, VideoParameters *p_Vid, Enc_Macroblock *curr_MBI, Enc_Slice *currSliceI, Enc_VideoParameters *p_VidI);
void getNonAffNeighbourV(Enc_Macroblock *currMB, int xN, int yN, int mb_size[2], PixelPos *pix);
int ReadBit(Slice *currSlice);
void WriteBit(int Bit);
void WriteFrame(int startBit, int endBit, Slice *currSlice);
int ReadPLNZV(int numcoeff, int *Run);
int Embeding(int *point, char *StringI, int *lv, int *rn, int numcoeff, int PLNZ);

int startOff, endOff;