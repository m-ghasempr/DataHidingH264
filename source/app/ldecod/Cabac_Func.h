#pragma once

#include "types.h"
#include "typedefs.h"
#include "Enc_Entropy.h"

#define B_BITS         10    // Number of bits to represent the whole coding interval
#define BITS_TO_LOAD   16
#define MAX_BITS       26          //(B_BITS + BITS_TO_LOAD)
#define MASK_BITS      18          //(MAX_BITS - 8)
#define ONE            0x04000000  //(1 << MAX_BITS)
#define ONE_M1         0x03FFFFFF  //(ONE - 1)
#define HALF           0x01FE      //(1 << (B_BITS-1)) - 2
#define QUARTER        0x0100      //(1 << (B_BITS-2))
#define MIN_BITS_TO_GO 0
#define B_LOAD_MASK    0xFFFF      // ((1<<BITS_TO_LOAD) - 1)
#define B_BITS    10      // Number of bits to represent the whole coding interval
#define HALF      0x01FE  //(1 << (B_BITS-1)) - 2
#define QUARTER   0x0100  //(1 << (B_BITS-2))

void reset_coding_state_cabacE(EncodingEnvironment *currMB, EncodingEnvironment *cs, int *len);
void store_coding_state_cabacE(EncodingEnvironment *currMB, EncodingEnvironment *cs, int *len);
void biari_encode_symbol(EncodingEnvironmentPtr eep, int symbol, BiContextTypePtr bi_ct);
void biari_encode_symbol_eq_prob(EncodingEnvironmentPtr eep, int symbol);
static void exp_golomb_encode_eq_prob(EncodingEnvironmentPtr eep_dp, unsigned int symbol, int k);
void unary_exp_golomb_level_encode(EncodingEnvironmentPtr eep_dp, unsigned int symbol, BiContextTypePtr ctx);
int writeCoeff4x4_CABAC(Enc_Macroblock* currMB, ColorPlane plane, int iX, int iY, int *ACLevel, int *ACRun);
void initData(Macroblock *currMB, Enc_Macroblock *currMBI);
int RunLevelAC(int block_y, int block_x, int *level, int *run, int *cof, int *PLNZ, int Y, int X);
void unary_exp_golomb_mv_encode(EncodingEnvironmentPtr eep_dp, unsigned int symbol, BiContextTypePtr ctx, unsigned int max_bin);
void biari_encode_symbol_final(EncodingEnvironmentPtr eep, int symbol);
void unary_bin_encode(EncodingEnvironmentPtr eep_dp, unsigned int symbol, BiContextTypePtr ctx, int ctx_offset);
void arienco_done_encoding(Enc_Macroblock *currMB, EncodingEnvironmentPtr eep);
int Embedding_CABAC(int *point, char *StringI, int *lv, int *rn, int numcoeff, int PLNZ);