/*!
 *************************************************************************************
 * \file rdo_quant.c
 *
 * \brief
 *    Rate Distortion Optimized Quantization based on VCEG-AH21
 *************************************************************************************
 */

#include "contributors.h"

#include <math.h>
#include <float.h>

#include "global.h"
#include "image.h"
#include "fmo.h"
#include "macroblock.h"
#include "mb_access.h"
#include "ratectl.h"
#include "rdo_quant.h"

const int entropyBits[128]= 
{
     895,    943,    994,   1048,   1105,   1165,   1228,   1294, 
    1364,   1439,   1517,   1599,   1686,   1778,   1875,   1978, 
    2086,   2200,   2321,   2448,   2583,   2725,   2876,   3034, 
    3202,   3380,   3568,   3767,   3977,   4199,   4435,   4684, 
    4948,   5228,   5525,   5840,   6173,   6527,   6903,   7303, 
    7727,   8178,   8658,   9169,   9714,  10294,  10914,  11575, 
   12282,  13038,  13849,  14717,  15650,  16653,  17734,  18899, 
   20159,  21523,  23005,  24617,  26378,  28306,  30426,  32768, 
   32768,  35232,  37696,  40159,  42623,  45087,  47551,  50015, 
   52479,  54942,  57406,  59870,  62334,  64798,  67262,  69725, 
   72189,  74653,  77117,  79581,  82044,  84508,  86972,  89436, 
   91900,  94363,  96827,  99291, 101755, 104219, 106683, 109146, 
  111610, 114074, 116538, 119002, 121465, 123929, 126393, 128857, 
  131321, 133785, 136248, 138712, 141176, 143640, 146104, 148568, 
  151031, 153495, 155959, 158423, 160887, 163351, 165814, 168278, 
  170742, 173207, 175669, 178134, 180598, 183061, 185525, 187989
};

extern const int maxpos       [];
extern const int type2ctx_bcbp[];
extern const int type2ctx_map [];
extern const int type2ctx_last[]; 
extern const int type2ctx_one []; // 7
extern const int type2ctx_abs []; // 7
extern const int max_c2       []; // 9

extern const int  pos2ctx_map8x8 [];
extern const int  pos2ctx_map8x4 [];
extern const int  pos2ctx_map4x4 [];
extern const int  pos2ctx_map2x4c[];
extern const int  pos2ctx_map4x4c[];
extern const int* pos2ctx_map    [];

extern const int  pos2ctx_last8x8 [];
extern const int  pos2ctx_last8x4 [];
extern const int  pos2ctx_last4x4 [];
extern const int  pos2ctx_last2x4c[];
extern const int  pos2ctx_last4x4c[];
extern const int* pos2ctx_last    [];


extern void SetMotionVectorPredictor (Macroblock *currMB, short  pmv[2], char   **refPic,
                                      short  ***tmp_mv, short  ref_frame,
                                      int    list,      int mb_x, int mb_y, 
                                      int    blockshape_x, int blockshape_y);
extern int *mvbits;


void precalculate_unary_exp_golomb_level()
{
  int state, ctx_state0, ctx_state1, estBits0, estBits1, symbol;

  for (state=0; state<=63; state++)
  {
    // symbol 0 is MPS
    ctx_state0=64+state;
    estBits0=entropyBits[127-ctx_state0];
    ctx_state1=63-state;
    estBits1=entropyBits[127-ctx_state1];

    for (symbol=0; symbol<MAX_PREC_COEFF; symbol++)
    {
      precalcUnaryLevelTab[ctx_state0][symbol]=est_unary_exp_golomb_level_bits(symbol, estBits0, estBits1);

      // symbol 0 is LPS
      precalcUnaryLevelTab[ctx_state1][symbol]=est_unary_exp_golomb_level_bits(symbol, estBits1, estBits0);
    }
  }
}

int est_unary_exp_golomb_level_bits(unsigned int symbol, int bits0, int bits1)
{
  unsigned int l,k;
  unsigned int exp_start = 13; // 15-2 : 0,1 level decision always sent
  int estBits;

  if (symbol==0)
  {
    return (bits0);
  }
  else
  {
    estBits=bits1;
    l=symbol;
    k=1;
    while (((--l)>0) && (++k <= exp_start))
    {
      estBits += bits1;
    }
    if (symbol < exp_start)
    {
      estBits += bits0;
    }
    else 
    {
      estBits += est_exp_golomb_encode_eq_prob(symbol-exp_start);
    }
  }
  return(estBits);
}

/*!
****************************************************************************
* \brief
*    estimate exp golomb bit cost 
****************************************************************************
*/
int est_exp_golomb_encode_eq_prob(unsigned int symbol)
{
  int k=0, estBits=0;

  while(1)
  {
    if (symbol >= (unsigned int)(1<<k))   
    {
      estBits++;
      symbol = symbol - (1<<k);
      k++;
    }
    else                  
    {
      estBits++;  
      while (k--)  
      {
        estBits++;
      }
      break;
    }
  }
  return(estBits);
}

/*!
****************************************************************************
* \brief
*   estimate bit cost for CBP, significant map and significant coefficients
****************************************************************************
*/
void estRunLevel_CABAC (Macroblock *currMB, int context) // marta - writes CABAC run/level 
{
  DataPartition*  dataPart = &(img->currentSlice->partArr[0]); // assumed that no DP is used (table assignSE2partition_NoDP)
  EncodingEnvironmentPtr eep_dp = &(dataPart->ee_cabac); 

  est_CBP_block_bit  (currMB, eep_dp, context);      
  //===== encode significance map =====
  est_significance_map         (currMB, eep_dp, context);      
  //===== encode significant coefficients =====
  est_significant_coefficients (currMB, eep_dp, context);
}

/*!
****************************************************************************
* \brief
*    estimate bit cost for each CBP bit
****************************************************************************
*/
void est_CBP_block_bit (Macroblock* currMB, EncodingEnvironmentPtr eep_dp, int type)
{
  int ctx;
  short cbp_bit;

  for (ctx=0; ctx<=3; ctx++)
  {
    cbp_bit=0;
    estBitsCabac[type].blockCbpBits[ctx][cbp_bit]=biari_no_bits(eep_dp, cbp_bit, img->currentSlice->tex_ctx->bcbp_contexts[type2ctx_bcbp[type]]+ctx);

    cbp_bit=1;
    estBitsCabac[type].blockCbpBits[ctx][cbp_bit]=biari_no_bits(eep_dp, cbp_bit, img->currentSlice->tex_ctx->bcbp_contexts[type2ctx_bcbp[type]]+ctx);
  }
}

/*!
****************************************************************************
* \brief
*    estimate CABAC bit cost for significant coefficient map
****************************************************************************
*/
void est_significance_map(Macroblock* currMB, EncodingEnvironmentPtr eep_dp, int type)
{
  int   k;
  unsigned short sig, last;
  int   k1  = maxpos[type]-1;
#if ENABLE_FIELD_CTX
  int   fld = ( img->structure!=FRAME || currMB->mb_field );
#else
  int   fld = 0;
#endif
  BiContextTypePtr  map_ctx   = img->currentSlice->tex_ctx->map_contexts [fld][type2ctx_map [type]];
  BiContextTypePtr  last_ctx  = img->currentSlice->tex_ctx->last_contexts[fld][type2ctx_last[type]];


  for (k=0; k<k1; k++) // if last coeff is reached, it has to be significant
  {
    sig   = 0;     
    estBitsCabac[type].significantBits[pos2ctx_map[type][k]][sig]=biari_no_bits  (eep_dp, sig,  map_ctx+pos2ctx_map     [type][k]);

    sig   = 1;     
    estBitsCabac[type].significantBits[pos2ctx_map[type][k]][sig]=biari_no_bits  (eep_dp, sig,  map_ctx+pos2ctx_map     [type][k]);

    last=0;
    estBitsCabac[type].lastBits[pos2ctx_last[type][k]][last]=biari_no_bits(eep_dp, last, last_ctx+pos2ctx_last[type][k]);

    last=1;
    estBitsCabac[type].lastBits[pos2ctx_last[type][k]][last]=biari_no_bits(eep_dp, last, last_ctx+pos2ctx_last[type][k]);
  }
  // if last coeff is reached, it has to be significant
  estBitsCabac[type].significantBits[pos2ctx_map[type][k1]][0]=0;
  estBitsCabac[type].significantBits[pos2ctx_map[type][k1]][1]=0;
  estBitsCabac[type].lastBits[pos2ctx_last[type][k1]][0]=0;
  estBitsCabac[type].lastBits[pos2ctx_last[type][k1]][1]=0;
}

/*!
****************************************************************************
* \brief
*    estimate bit cost of significant coefficient
****************************************************************************
*/
void est_significant_coefficients (Macroblock* currMB, EncodingEnvironmentPtr eep_dp,  int type)
{
  int   ctx;
  short greater_one;
  int maxCtx=imin(4, max_c2[type]);

  for (ctx=0; ctx<=4; ctx++){    
    greater_one=0;
    estBitsCabac[type].greaterOneBits[0][ctx][greater_one]=
      biari_no_bits (eep_dp, greater_one, img->currentSlice->tex_ctx->one_contexts[type2ctx_one[type]] + ctx);

    greater_one=1;
    estBitsCabac[type].greaterOneBits[0][ctx][greater_one]=
      biari_no_bits (eep_dp, greater_one, img->currentSlice->tex_ctx->one_contexts[type2ctx_one[type]] + ctx);
  }

  for (ctx=0; ctx<=maxCtx; ctx++){
    estBitsCabac[type].greaterOneBits[1][ctx][0]=
      biari_no_bits(eep_dp, 0, img->currentSlice->tex_ctx->abs_contexts[type2ctx_abs[type]] + ctx);

    estBitsCabac[type].greaterOneState[ctx]=biari_state(eep_dp, 0, img->currentSlice->tex_ctx->abs_contexts[type2ctx_abs[type]] + ctx);

    estBitsCabac[type].greaterOneBits[1][ctx][1]=
      biari_no_bits(eep_dp, 1, img->currentSlice->tex_ctx->abs_contexts[type2ctx_abs[type]] + ctx);
  }
}

int biari_no_bits(EncodingEnvironmentPtr eep, signed short symbol, BiContextTypePtr bi_ct )
{
  int ctx_state, estBits;

  symbol = (short) (symbol != 0);

  ctx_state=(symbol==bi_ct->MPS)?64+bi_ct->state:63-bi_ct->state;
  estBits=entropyBits[127-ctx_state];

  return(estBits);
}

int biari_state(EncodingEnvironmentPtr eep, signed short symbol, BiContextTypePtr bi_ct )
{ 
  int ctx_state;

  symbol = (short) (symbol != 0);
  ctx_state=(symbol==bi_ct->MPS)?64+bi_ct->state:63-bi_ct->state;

  return(ctx_state);
}


/*!
****************************************************************************
* \brief
*    estimate CABAC CBP bits
****************************************************************************
*/
int est_write_and_store_CBP_block_bit(Macroblock* currMB, int type, int block_y, int block_x) // marta - CBP
{
#define BIT_SET(x,n)  ((int)(((x)&((int64)1<<(n)))>>(n)))

  int bit, default_bit = (IS_INTRA(currMB) ? 1 : 0);
  int upper_bit   = default_bit;
  int left_bit    = default_bit;
  int ctx, estBits;

  int bit_pos_a   = 0;
  int bit_pos_b   = 0;
  int *mb_size = img->mb_size[IS_LUMA];

  if (type!=LUMA_8x8)
  {
    PixelPos block_a, block_b;

    get4x4Neighbour (currMB, block_x - 1 , block_y    , mb_size, &block_a);    
    get4x4Neighbour (currMB, block_x     , block_y - 1, mb_size, &block_b);    

    if (block_a.available)
      bit_pos_a = 4*block_a.y + block_a.x;
    if (block_b.available)
      bit_pos_b = 4*block_b.y + block_b.x;

    bit = 1; // 4x4: bit=1

    if (block_b.available)
    {
      if(img->mb_data[block_b.mb_addr].mb_type==IPCM)
        upper_bit=1;
      else
        upper_bit = BIT_SET(img->mb_data[block_b.mb_addr].cbp_bits[0],bit+bit_pos_b);
    }


    if (block_a.available)
    {
      if(img->mb_data[block_a.mb_addr].mb_type==IPCM)
        left_bit=1;
      else
        left_bit = BIT_SET(img->mb_data[block_a.mb_addr].cbp_bits[0],bit+bit_pos_a);
    }

    ctx = 2*upper_bit+left_bit;
    //===== encode symbol =====
    estBits=estBitsCabac[type].blockCbpBits[ctx][0]-estBitsCabac[type].blockCbpBits[ctx][1];
  }
  else
  {
    estBits=0;
  }
  return(estBits);
}

/*!
****************************************************************************
* \brief
*    Rate distortion optimized trellis quantization
****************************************************************************
*/
void est_writeRunLevel_CABAC(levelDataStruct levelData[], int levelTabMin[], int type, double lambda, int kInit, int kStop, 
                             int noCoeff, int estCBP)
{
  int   k, i;
  int estBits;
  double lagr, lagrMin=0, lagrTabMin, lagrTab;
  int   c1 = 1, c2 = 0, c1Tab[3], c2Tab[3];
  int   iBest, levelTab[64];
  int   ctx, greater_one, last, maxK;
  double   lagrAcc, lagrLastMin=0, lagrLast;
  int      kBest=0, kStart, first;

  maxK=maxpos[type];
  for (k=0; k<=maxK; k++)
  {
    levelTabMin[k]=0;
  }

  if (noCoeff>0)
  {
    if (noCoeff>1)
    {
      kStart=kInit; kBest=0; first=1; 

      lagrAcc=0; 
      for (k=kStart; k<=kStop; k++)
      {
        lagrAcc+=levelData[k].errLevel[0];
      }

      if (levelData[kStart].noLevels>2)
      { 
        lagrAcc-=levelData[kStart].errLevel[0];
        lagrLastMin=lambda*(estBitsCabac[type].lastBits[pos2ctx_last[type][kStart]][1]-estBitsCabac[type].lastBits[pos2ctx_last[type][kStart]][0])+lagrAcc;

        kBest=kStart;
        kStart=kStart+1;
        first=0;
      }

      for (k=kStart; k<=kStop; k++)
      {
        lagrMin=levelData[k].errLevel[0]+lambda*estBitsCabac[type].significantBits[pos2ctx_map[type][k]][0];

        lagrAcc-=levelData[k].errLevel[0];
        if (levelData[k].noLevels>1)
        { 
          estBits=SIGN_BITS+estBitsCabac[type].significantBits[pos2ctx_map[type][k]][1]+
            estBitsCabac[type].greaterOneBits[0][4][0];

          lagrLast=levelData[k].errLevel[1]+lambda*(estBits+estBitsCabac[type].lastBits[pos2ctx_last[type][k]][1])+lagrAcc;
          lagr=levelData[k].errLevel[1]+lambda*(estBits+estBitsCabac[type].lastBits[pos2ctx_last[type][k]][0]);

          lagrMin=(lagr<lagrMin)?lagr:lagrMin;

          if (lagrLast<lagrLastMin || first==1)
          {
            kBest=k;
            first=0;
            lagrLastMin=lagrLast;
          }

        }
        lagrAcc+=lagrMin;
      }

      kStart=kBest;
    }
    else
    {
      kStart=kStop;
    }

    lagrTabMin=0;
    for (k=0; k<=kStart; k++)
    {
      lagrTabMin+=levelData[k].errLevel[0];
    }
    // Initial Lagrangian calculation
    lagrTab=0;

    //////////////////////////

    lagrTabMin+=(lambda*estCBP);
    iBest=0; first=1;
    for (k=kStart; k>=0; k--)
    {
      last=(k==kStart);
      if (!last)
      {
        lagrMin=levelData[k].errLevel[0]+lambda*estBitsCabac[type].significantBits[pos2ctx_map[type][k]][0];
        iBest=0;
        first=0;
      }

      for (i=1; i < levelData[k].noLevels; i++)
      {
        estBits=SIGN_BITS+estBitsCabac[type].significantBits[pos2ctx_map[type][k]][1];
        estBits+=estBitsCabac[type].lastBits[pos2ctx_last[type][k]][last];

        // greater than 1
        greater_one = (levelData[k].level[i]>1);

        c1Tab[i]=c1;   c2Tab[i]=c2;

        ctx = imin(c1Tab[i],4);  
        estBits+=estBitsCabac[type].greaterOneBits[0][ctx][greater_one];

        // magnitude if greater than 1
        if (greater_one)
        {
          ctx = imin(c2Tab[i], max_c2[type]);
          if ((levelData[k].level[i]-2)<MAX_PREC_COEFF)
          {
            estBits+=precalcUnaryLevelTab[estBitsCabac[type].greaterOneState[ctx]][levelData[k].level[i]-2];
          }
          else
          {
            estBits+=est_unary_exp_golomb_level_encode(levelData[k].level[i]-2, ctx, type);
          }

          c1Tab[i] = 0;
          c2Tab[i]++;
        }
        else if (c1Tab[i])
        {
          c1Tab[i]++;
        }

        lagr=levelData[k].errLevel[i]+lambda*estBits;
        if (lagr<lagrMin || first==1)
        {
          iBest=i;
          lagrMin=lagr;
          first=0;
        }
      }

      if (iBest>0)
      {
        c1=c1Tab[iBest]; c2=c2Tab[iBest];
      }

      levelTab[k]=levelData[k].level[iBest];
      lagrTab+=lagrMin;
    }
    ///////////////////////////////////

    if (lagrTab<lagrTabMin)
    {
      for (k=0; k<=kStart; k++)
      {
        levelTabMin[k]=levelTab[k];
      }
    }
  }
}

/*!
****************************************************************************
* \brief
*    estimate unary exp golomb bit cost
****************************************************************************
*/
int est_unary_exp_golomb_level_encode(unsigned int symbol, int ctx, int type)
{
  unsigned int l,k;
  unsigned int exp_start = 13; // 15-2 : 0,1 level decision always sent
  int estBits;

  if (symbol==0)
  {
    estBits=estBitsCabac[type].greaterOneBits[1][ctx][0];
    return (estBits);
  }
  else
  {
    estBits=estBitsCabac[type].greaterOneBits[1][ctx][1];
    l=symbol;
    k=1;
    while (((--l)>0) && (++k <= exp_start))
    {
      estBits+=estBitsCabac[type].greaterOneBits[1][ctx][1];
    }
    if (symbol < exp_start)
    {
      estBits+=estBitsCabac[type].greaterOneBits[1][ctx][0];
    }
    else 
    {
      estBits+=est_exp_golomb_encode_eq_prob(symbol-exp_start);
    }
  }
  return(estBits);
}

#define RDOQ_BASE 0
void trellis_mp(Macroblock *currMB, int CurrentMbAddr, Boolean prev_recode_mb)
{
  int masterQP = 0, deltaQP;
  int qp_left, qp_up;
#if RDOQ_BASE
  const int deltaQPTabB[] = {0,  1, -1,  2, 3, -2, 4,  5, -3};
  const int deltaQPTabP[] = {0, -1,  1, -2, 2, -3, 3, -4,  4};
#endif
  int   deltaQPCnt; 
  int   qp_anchor; 
  int   prev_mb = FmoGetPreviousMBNr(img->current_mb_nr);
  int   qp_offset = (img->type == B_SLICE) ? (params->RDOQ_QP_Num / 3): (params->RDOQ_QP_Num >> 1);

  masterQP = img->masterQP = img->qp;
  Motion_Selected = 0;
  rddata_trellis_best.min_rdcost = 1e30;

  estRunLevel_CABAC(currMB, LUMA_4x4); 
  estRunLevel_CABAC(currMB, LUMA_16AC);
  if (params->Transform8x8Mode)
    estRunLevel_CABAC(currMB, LUMA_8x8);

  qp_left   = (currMB->mb_available_left) ? currMB->mb_available_left->qp : img->masterQP;
  qp_up     = (currMB->mb_available_up)   ? currMB->mb_available_up->qp   : img->masterQP;
  qp_anchor = (qp_left + qp_up + 1)>>1;

  for (deltaQPCnt=0; deltaQPCnt < params->RDOQ_QP_Num; deltaQPCnt++)
  {
    rdopt = &rddata_trellis_curr;
#if RDOQ_BASE
    if (img->type == B_SLICE)
      deltaQP = deltaQPTabB[deltaQPCnt];      
    else
      deltaQP = deltaQPTabP[deltaQPCnt];
#else

    // It seems that pushing the masterQP as first helps things when fast me is enabled. 
    // Could there be an issue with motion estimation?
    if (deltaQPCnt == 0)
      deltaQP = 0;
    else if (deltaQPCnt <= qp_offset)
      deltaQP = deltaQPCnt - 1 - qp_offset;
    else
      deltaQP = deltaQPCnt - qp_offset;
    //printf("qp %d %d %d\n", deltaQP,  deltaQPCnt, masterQP);
#endif

    img->qp = iClip3(-img->bitdepth_luma_qp_scale, 51, masterQP + deltaQP);

#if 0
    if(deltaQP != 0 && !(img->qp - qp_anchor >= -2 && img->qp - qp_anchor <= 1) && currMB->mb_available_left && currMB->mb_available_up && img->type == P_SLICE)
      continue; 
    if(deltaQP != 0 && !(img->qp - qp_anchor >= -1 && img->qp - qp_anchor <= 2) && currMB->mb_available_left && currMB->mb_available_up && img->type == B_SLICE)
      continue;
#endif
    if (img->current_mb_nr ==0 && deltaQP != 0)
      continue;
    
    reset_macroblock(currMB, prev_mb);
    currMB->qp       = img->qp;
    currMB->delta_qp = currMB->qp - currMB->prev_qp;
    update_qp (currMB);    
    
    delta_qp_mbaff[currMB->mb_field][img->bot_MB] = currMB->delta_qp;
    qp_mbaff      [currMB->mb_field][img->bot_MB] = currMB->qp;

    encode_one_macroblock (currMB);

    if ( rddata_trellis_curr.min_rdcost < rddata_trellis_best.min_rdcost)
      copy_rddata_trellis(&rddata_trellis_best,rdopt);
    
    if (params->RDOQ_CP_MV)
      Motion_Selected = 1;

#if (!RDOQ_BASE)
    if ((params->RDOQ_Fast) && (img->qp - rddata_trellis_best.qp > 1))
      break;
    if ((params->RDOQ_Fast) && (rddata_trellis_curr.cbp == 0) && (rddata_trellis_curr.mb_type != 0))
      break;
#endif
  }

  reset_macroblock(currMB, prev_mb);
  rdopt = &rddata_trellis_best;

  copy_rdopt_data (currMB, FALSE);  // copy the MB data for Top MB from the temp buffers
  write_one_macroblock (currMB, 1, prev_recode_mb);
  img->qp = masterQP;
}

void trellis_sp(Macroblock *currMB, int CurrentMbAddr, Boolean prev_recode_mb)
{
  img->masterQP = img->qp;

  estRunLevel_CABAC(currMB, LUMA_4x4); 
  estRunLevel_CABAC(currMB, LUMA_16AC);
  if (params->Transform8x8Mode)
    estRunLevel_CABAC(currMB, LUMA_8x8);

  encode_one_macroblock (currMB);
  write_one_macroblock (currMB, 1, prev_recode_mb);    
}

void trellis_coding(Macroblock *currMB, int CurrentMbAddr, Boolean prev_recode_mb)
{
  if (params->RDOQ_QP_Num > 1)
  {
    trellis_mp(currMB, CurrentMbAddr, prev_recode_mb);   
  }
  else
  {
    trellis_sp(currMB, CurrentMbAddr, prev_recode_mb);   
  }
}


void RDOQ_update_mode(RD_PARAMS *enc_mb, int bslice)
{
  int i;
  for(i=0; i<MAXMODE; i++)
    enc_mb->valid[i] = 0;

  enc_mb->valid[rdopt->mb_type] = 1;

  if(rdopt->mb_type  == P8x8)
  {            
    enc_mb->valid[4] = (params->InterSearch[bslice][4]);
    enc_mb->valid[5] = (params->InterSearch[bslice][5] && !(params->Transform8x8Mode==2));
    enc_mb->valid[6] = (params->InterSearch[bslice][6] && !(params->Transform8x8Mode==2));
    enc_mb->valid[7] = (params->InterSearch[bslice][7] && !(params->Transform8x8Mode==2));
  }
}

void copy_rddata_trellis (RD_DATA *dest, RD_DATA *src)
{
  int j; 

  dest->min_rdcost = src->min_rdcost;
  dest->min_dcost  = src->min_dcost;

  memcpy(&dest->rec_mbY[0][0],&src->rec_mbY[0][0], MB_PIXELS * sizeof(imgpel));

  if (img->yuv_format != YUV400) 
  {
    // we could allocate these dynamically to improve performance.
    memcpy(&dest->rec_mb_cr[0][0][0],&src->rec_mb_cr[0][0][0], 2 * MB_PIXELS * sizeof(imgpel));
  }

  memcpy(&dest->cofAC[0][0][0][0], &src->cofAC[0][0][0][0], (4 + img->num_blk8x8_uv) * 4 * 2 * 65 * sizeof(int));
  memcpy(&dest->cofDC[0][0][0], &src->cofDC[0][0][0], 3 * 2 * 18 * sizeof(int));

  dest->mb_type = src->mb_type;
  memcpy(dest->b8mode, src->b8mode, BLOCK_MULTIPLE * sizeof(short));
  memcpy(dest->b8pdir, src->b8pdir, BLOCK_MULTIPLE * sizeof(short));
  dest->cbp = src->cbp;
  dest->mode = src->mode;
  dest->i16offset = src->i16offset;
  dest->c_ipred_mode = src->c_ipred_mode;
  dest->luma_transform_size_8x8_flag = src->luma_transform_size_8x8_flag;
  dest->NoMbPartLessThan8x8Flag = src->NoMbPartLessThan8x8Flag;
  dest->qp = src->qp;

  dest->prev_qp = src->prev_qp;
  dest->prev_dqp = src->prev_dqp;
  dest->delta_qp = src->delta_qp;
  dest->prev_cbp = src->prev_cbp;
  dest->cbp_blk = src->cbp_blk;
  dest->bi_pred_me = src->bi_pred_me;

  if (img->type != I_SLICE)
  {

    memcpy(&dest->all_mv [0][0][0][0][0][0], &src->all_mv [0][0][0][0][0][0], 2 * img->max_num_references * 9 * 4 * 4 * 2 * sizeof(short));
    memcpy(&dest->pred_mv[0][0][0][0][0][0], &src->pred_mv[0][0][0][0][0][0], 2 * img->max_num_references * 9 * 4 * 4 * 2 * sizeof(short));
  }

  memcpy(dest->intra_pred_modes,src->intra_pred_modes, MB_BLOCK_PARTITIONS * sizeof(char));
  memcpy(dest->intra_pred_modes8x8,src->intra_pred_modes8x8, MB_BLOCK_PARTITIONS * sizeof(char));
  for(j = img->block_y; j < img->block_y + BLOCK_MULTIPLE; j++)
    memcpy(&dest->ipredmode[j][img->block_x],&src->ipredmode[j][img->block_x], BLOCK_MULTIPLE * sizeof(char));
  memcpy(&dest->refar[LIST_0][0][0], &src->refar[LIST_0][0][0], 2 * BLOCK_MULTIPLE * BLOCK_MULTIPLE * sizeof(char));
}                            

void updateMV_mp(int *m_cost, short ref, int list, int h, int v, int blocktype, int *lambda_factor, int block8x8)
{
  int       i, j;
  int       bsx       = params->blc_size[blocktype][0];
  int       bsy       = params->blc_size[blocktype][1];
  short     tmp_pred_mv[2];
  short*    pred_mv = img->pred_mv[list][ref][blocktype][v][h];
  short     all_mv[2];
  Macroblock *currMB = &img->mb_data[img->current_mb_nr];
  if ( (params->Transform8x8Mode == 1) && (blocktype == 4) && currMB->luma_transform_size_8x8_flag)
  {
    all_mv[0] = tmp_mv8[list][ref][v][h][0];
    all_mv[1] = tmp_mv8[list][ref][v][h][1];
    tmp_pred_mv[0] = tmp_pmv8[list][ref][v][h][0];
    tmp_pred_mv[1] = tmp_pmv8[list][ref][v][h][1];
    *m_cost   = motion_cost8[list][ref][block8x8];
  }
  else
  {
    all_mv[0] = rddata_trellis_best.all_mv[list][ref][blocktype][v][h][0];
    all_mv[1] = rddata_trellis_best.all_mv[list][ref][blocktype][v][h][1];
    tmp_pred_mv[0] = rddata_trellis_best.pred_mv[list][ref][blocktype][v][h][0];
    tmp_pred_mv[1] = rddata_trellis_best.pred_mv[list][ref][blocktype][v][h][1];
  }

  for (j = 0; j < (bsy>>2); j++)
  {
    for (i = 0; i < (bsx>>2); i++) 
      memcpy(img->all_mv[list][ref][blocktype][v+j][h+i], all_mv, 2 * sizeof(short));
  }

  SetMotionVectorPredictor (currMB, pred_mv, enc_picture->ref_idx[list], enc_picture->mv[list], ref, list, h<<2, v<<2, bsx, bsy);
  
  if ( (tmp_pred_mv[0] != pred_mv[0]) || (tmp_pred_mv[1] != pred_mv[1]) )
  {
    *m_cost -= MV_COST_SMP (lambda_factor[H_PEL], all_mv[0], all_mv[1], tmp_pred_mv[0], tmp_pred_mv[1]);
    *m_cost += MV_COST_SMP (lambda_factor[H_PEL], all_mv[0], all_mv[1], pred_mv[0], pred_mv[1]);
  }
}



