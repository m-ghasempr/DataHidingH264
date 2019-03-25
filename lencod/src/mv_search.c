
/*!
 *************************************************************************************
 * \file mv_search.c
 *
 * \brief
 *    Motion Vector Search, unified for B and P Pictures
 *
 * \author
 *    Main contributors (see contributors.h for copyright, address and affiliation details)
 *      - Stephan Wenger                  <stewe@cs.tu-berlin.de>
 *      - Inge Lille-Langoy               <inge.lille-langoy@telenor.com>
 *      - Rickard Sjoberg                 <rickard.sjoberg@era.ericsson.se>
 *      - Stephan Wenger                  <stewe@cs.tu-berlin.de>
 *      - Jani Lainema                    <jani.lainema@nokia.com>
 *      - Detlev Marpe                    <marpe@hhi.de>
 *      - Thomas Wedi                     <wedi@tnt.uni-hannover.de>
 *      - Heiko Schwarz                   <hschwarz@hhi.de>
 *      - Alexis Michael Tourapis         <alexismt@ieee.org>
 *
 *************************************************************************************
*/

#include "contributors.h"

#include <math.h>
#include <limits.h>
#include <time.h>

#include "global.h"

#include "image.h"
#include "mv_search.h"
#include "refbuf.h"
#include "memalloc.h"
#include "mb_access.h"
#include "macroblock.h"
#include "mc_prediction.h"
#include "conformance.h"
#include "mode_decision.h"

// Motion estimation distortion header file
#include "me_distortion.h"

// Motion estimation search algorithms
#include "me_epzs.h"
#include "me_epzs_int.h"
#include "me_fullfast.h"
#include "me_fullsearch.h"
#include "me_umhex.h"
#include "me_umhexsmp.h"
#include "rdoq.h"


static const short bx0[5][4] = {{0,0,0,0}, {0,0,0,0}, {0,0,0,0}, {0,2,0,0}, {0,2,0,2}};
static const short by0[5][4] = {{0,0,0,0}, {0,0,0,0}, {0,2,0,0}, {0,0,0,0}, {0,0,2,2}};



static distblk GetSkipCostMB          (Macroblock *currMB, int lambda);
static distblk BiPredBlockMotionSearch(Macroblock *currMB, MEBlock *, MotionVector*, int, int , int*);

/*!
 ************************************************************************
 * \brief
 *    Set search range. This needs to be changed to provide 2D support
 ************************************************************************
 */
void get_search_range(MEBlock *mv_block, InputParameters *p_Inp, short ref, int blocktype)
{
  mv_block->searchRange = mv_block->p_Vid->searchRange;
  //----- set search range ---
  if (p_Inp->full_search == 1)
  {
    int scale = (imin(ref,1)+1);
    mv_block->searchRange.min_x /= scale;
    mv_block->searchRange.max_x /= scale;
    mv_block->searchRange.min_y /= scale;
    mv_block->searchRange.max_y /= scale;
  }
  else if  (p_Inp->full_search != 2)
  {
    int scale = ((imin(ref,1)+1) * imin(2,blocktype));
    mv_block->searchRange.min_x /= scale;
    mv_block->searchRange.max_x /= scale;
    mv_block->searchRange.min_y /= scale;
    mv_block->searchRange.max_y /= scale;
  }
}

/*!
 ************************************************************************
 * \brief
 *    Set search range. This needs to be changed to provide 2D support
 ************************************************************************
 */
static inline void set_me_parameters( char  **ref_array, short ***mv_array ,short *all_mv, short ref, int step_h, int step_v, int pic_block_x)
{
  int i, j;
  for (j = 0; j < step_v; j++)
  {
    memset(&ref_array [j][pic_block_x], ref, step_h * sizeof(char));
  }

  // Set first line
  for (i=pic_block_x; i<pic_block_x + step_h; i++)
  {
    memcpy(mv_array  [0][i], all_mv, 2* sizeof(short));
  }
  // Set remaining lines 
  for (j = 1; j < step_v; j++)
  {
    memcpy(mv_array  [j][pic_block_x], mv_array  [j - 1][pic_block_x], 2 * step_h * sizeof(short));
  }
}

/*!
 ************************************************************************
 * \brief
 *    Set ME access method
 ************************************************************************
 */
void set_access_method(int *access_method, MotionVector *blk, int min_x, int min_y, int max_x, int max_y)
{
  if ( (blk->mv_x > min_x) && (blk->mv_x < max_x) && (blk->mv_y > min_y) && (blk->mv_y < max_y))
  {
    *access_method = FAST_ACCESS;
  }
  else
  {
    *access_method = UMV_ACCESS;
  }
}

/*!
 ************************************************************************
 * \brief
 *    Initialize ME engine
 ************************************************************************
 */
void init_ME_engine(Macroblock *currMB)
{
  InputParameters *p_Inp = currMB->p_Inp;
  switch (p_Inp->SearchMode)
  {
   case EPZS:
     if (p_Inp->EPZSSubPelGrid)
     {
       currMB->IntPelME       = EPZSIntPelBlockMotionSearch;
       currMB->SubPelME       = (p_Inp->EPZSSubPelME) ? EPZSSubPelBlockMotionSearch : SubPelBlockMotionSearch;
       currMB->BiPredME       = EPZSIntBiPredBlockMotionSearch;
       currMB->SubPelBiPredME = (p_Inp->EPZSSubPelMEBiPred) ? EPZSSubPelBlockSearchBiPred : SubPelBlockSearchBiPred;
       
     }
     else
     {
       currMB->IntPelME       = EPZSPelBlockMotionSearch;
       currMB->BiPredME       = EPZSBiPredBlockMotionSearch;
       currMB->SubPelBiPredME = (p_Inp->EPZSSubPelMEBiPred) ? EPZSSubPelBlockSearchBiPred : SubPelBlockSearchBiPred;
       currMB->SubPelME       = (p_Inp->EPZSSubPelME) ? EPZSSubPelBlockMotionSearch : SubPelBlockMotionSearch;
     }
     break;
   case UM_HEX:
     currMB->IntPelME       = UMHEXIntegerPelBlockMotionSearch;
     currMB->BiPredME       = UMHEXBipredIntegerPelBlockMotionSearch;
     currMB->SubPelBiPredME = SubPelBlockSearchBiPred;
     currMB->SubPelME       = UMHEXSubPelBlockME;
     break;
   case UM_HEX_SIMPLE:
     currMB->IntPelME       = smpUMHEXIntegerPelBlockMotionSearch;
     currMB->BiPredME       = smpUMHEXBipredIntegerPelBlockMotionSearch;
     currMB->SubPelBiPredME = SubPelBlockSearchBiPred;
     currMB->SubPelME       = smpUMHEXSubPelBlockME;
     break;
   case FULL_SEARCH:
     currMB->IntPelME       = FullPelBlockMotionSearch;
     currMB->BiPredME       = FullPelBlockMotionBiPred;
     currMB->SubPelBiPredME = SubPelBlockSearchBiPred;
     currMB->SubPelME       = SubPelBlockMotionSearch;
     break;
   case FAST_FULL_SEARCH:
   default:
     currMB->IntPelME       = FastFullPelBlockMotionSearch;
     currMB->BiPredME       = FullPelBlockMotionBiPred;
     currMB->SubPelBiPredME = SubPelBlockSearchBiPred;
     currMB->SubPelME       = SubPelBlockMotionSearch;
     break;
  }
}

/*!
 ************************************************************************
 * \brief
 *    Prepare Motion Estimation parameters for single list ME
 ************************************************************************
 */
void PrepareMEParams(Slice *currSlice, MEBlock *mv_block, int ChromaMEEnable, int list, int ref)
{
  if (mv_block->apply_weights)
  {
    mv_block->weight_luma = currSlice->wp_weight[list][ref][0];
    mv_block->offset_luma = currSlice->wp_offset[list][ref][0];

    if ( ChromaMEEnable)
    {
      mv_block->weight_cr[0] = currSlice->wp_weight[list][ref][1];
      mv_block->weight_cr[1] = currSlice->wp_weight[list][ref][2];
      mv_block->offset_cr[0] = currSlice->wp_offset[list][ref][1];
      mv_block->offset_cr[1] = currSlice->wp_offset[list][ref][2];
    }
  }
}

/*!
 ************************************************************************
 * \brief
 *    Prepare Motion Estimation parameters for bipred list ME
 ************************************************************************
 */
void PrepareBiPredMEParams(Slice *currSlice, MEBlock *mv_block, int ChromaMEEnable, int list, int list_offset, int ref)
{
  if (mv_block->apply_weights)
  {
    if (list == LIST_0)
    {
      mv_block->weight1  = currSlice->wbp_weight[list_offset         ][ref][0][0];
      mv_block->weight2  = currSlice->wbp_weight[list_offset + LIST_1][ref][0][0];
      mv_block->offsetBi = (currSlice->wp_offset[list_offset         ][ref][0] + currSlice->wp_offset[list_offset + LIST_1][ref][0] + 1)>>1;

      if ( ChromaMEEnable)
      {
        mv_block->weight1_cr[0] = currSlice->wbp_weight[list_offset         ][ref][0][1];
        mv_block->weight1_cr[1] = currSlice->wbp_weight[list_offset         ][ref][0][2];
        mv_block->weight2_cr[0] = currSlice->wbp_weight[list_offset + LIST_1][ref][0][1];
        mv_block->weight2_cr[1] = currSlice->wbp_weight[list_offset + LIST_1][ref][0][2];
   
        mv_block->offsetBi_cr[0] = (currSlice->wp_offset[list_offset        ][ref][1] + currSlice->wp_offset[list_offset + LIST_1][ref][1] + 1) >> 1;
        mv_block->offsetBi_cr[1] = (currSlice->wp_offset[list_offset        ][ref][2] + currSlice->wp_offset[list_offset + LIST_1][ref][2] + 1) >> 1;
      }
    }
    else
    {
      mv_block->weight1  = currSlice->wbp_weight[list_offset + LIST_1][0  ][ref][0];
      mv_block->weight2  = currSlice->wbp_weight[list_offset         ][0  ][ref][0];
      mv_block->offsetBi = (currSlice->wp_offset[list_offset + LIST_1][0][0] + currSlice->wp_offset[list_offset][0][0] + 1)>>1;

      if ( ChromaMEEnable)
      {
        mv_block->weight1_cr[0] = currSlice->wbp_weight[list_offset + LIST_1][0  ][ref][1];
        mv_block->weight1_cr[1] = currSlice->wbp_weight[list_offset + LIST_1][0  ][ref][2];
        mv_block->weight2_cr[0] = currSlice->wbp_weight[list_offset         ][0  ][ref][1];
        mv_block->weight2_cr[1] = currSlice->wbp_weight[list_offset         ][0  ][ref][2];

        mv_block->offsetBi_cr[0] = (currSlice->wp_offset[list_offset + LIST_1][0  ][1] + currSlice->wp_offset[list_offset         ][0  ][1] + 1) >> 1;
        mv_block->offsetBi_cr[1] = (currSlice->wp_offset[list_offset + LIST_1][0  ][2] + currSlice->wp_offset[list_offset         ][0  ][2] + 1) >> 1;
      }
    }
  }
  else
  {
    mv_block->weight1 = (short) (1 << currSlice->luma_log_weight_denom);
    mv_block->weight2 = (short) (1 << currSlice->luma_log_weight_denom);
    mv_block->offsetBi = 0;
    if ( ChromaMEEnable)
    {
      mv_block->weight1_cr[0] = 1<<currSlice->chroma_log_weight_denom;
      mv_block->weight1_cr[1] = 1<<currSlice->chroma_log_weight_denom;
      mv_block->weight2_cr[0] = 1<<currSlice->chroma_log_weight_denom;
      mv_block->weight2_cr[1] = 1<<currSlice->chroma_log_weight_denom;
      mv_block->offsetBi_cr[0] = 0;
      mv_block->offsetBi_cr[1] = 0;
    }
  }
}

/*!
 ************************************************************************
 * \brief
 *    Get current block spatial neighbors
 ************************************************************************
 */
void get_neighbors(Macroblock *currMB,       // <--  current Macroblock
                   PixelPos   *block,        // <--> neighbor blocks
                   int         mb_x,         // <--  block x position
                   int         mb_y,         // <--  block y position
                   int         blockshape_x  // <--  block width
                   )
{
  VideoParameters *p_Vid = currMB->p_Vid;
  int *mb_size = p_Vid->mb_size[IS_LUMA];

  get4x4Neighbour(currMB, mb_x - 1,            mb_y    , mb_size, &block[0]);
  get4x4Neighbour(currMB, mb_x,                mb_y - 1, mb_size, &block[1]);
  get4x4Neighbour(currMB, mb_x + blockshape_x, mb_y - 1, mb_size, &block[2]);
  get4x4Neighbour(currMB, mb_x - 1,            mb_y - 1, mb_size, &block[3]);

  if (mb_y > 0)
  {
    if (mb_x < 8)  // first column of 8x8 blocks
    {
      if (mb_y == 8 )
      {
        if (blockshape_x == MB_BLOCK_SIZE)      
          block[2].available  = 0;
      }
      else if (mb_x + blockshape_x == 8)
      {
        block[2].available = 0;
      }
    }
    else if (mb_x + blockshape_x == MB_BLOCK_SIZE)
    {
      block[2].available = 0;
    }
  }

  if (!block[2].available)
  {
    block[2] = block[3];
  }
}

/*!
************************************************************************
* \brief
*    Initialize the motion search
************************************************************************
*/
void Init_Motion_Search_Module (VideoParameters *p_Vid, InputParameters *p_Inp)
{
  int bits;
  int i_min, i_max,k;
  int i, l;

  int search_range               = p_Inp->search_range;
  int max_search_points          = imax(9, (2 * search_range + 1) * (2 * search_range + 1));
  int max_ref_bits               = 1 + 2 * (int)floor(log(imax(16, p_Vid->max_num_references + 1)) / log(2) + 1e-10);
  int max_ref                    = (1<<((max_ref_bits>>1)+1))-1;
  int number_of_subpel_positions = 4 * (2*search_range+3);
  int max_mv_bits                = 3 + 2 * (int)ceil (log(number_of_subpel_positions + 1) / log(2) + 1e-10);
  int max_mvd                    = (1<<( max_mv_bits >>1)   ) - 1;
  
  p_Vid->max_mvd = max_mvd;
  p_Vid->imgpel_abs_range          = (imax(p_Vid->max_pel_value_comp[0],p_Vid->max_pel_value_comp[1]) + 1) * 64;

  //=====   CREATE ARRAYS   =====
  //-----------------------------
  if ((p_Vid->spiral_search = (MotionVector*)calloc(max_search_points, sizeof(MotionVector))) == NULL)
    no_mem_exit("Init_Motion_Search_Module: p_Vid->spiral_search");
  if ((p_Vid->spiral_hpel_search = (MotionVector*)calloc(max_search_points, sizeof(MotionVector))) == NULL)
    no_mem_exit("Init_Motion_Search_Module: p_Vid->spiral_hpel_search");
  if ((p_Vid->spiral_qpel_search = (MotionVector*)calloc(max_search_points, sizeof(MotionVector))) == NULL)
    no_mem_exit("Init_Motion_Search_Module: p_Vid->spiral_qpel_search");

  if ((p_Vid->mvbits = (int*)calloc(2 * max_mvd + 1, sizeof(int))) == NULL)
    no_mem_exit("Init_Motion_Search_Module: p_Vid->mvbits");

  if ((p_Vid->refbits = (int*)calloc(max_ref, sizeof(int))) == NULL)
    no_mem_exit("Init_Motion_Search_Module: p_Vid->refbits");

#if (JM_MEM_DISTORTION)
  if ((p_Vid->imgpel_abs = (int*)calloc(p_Vid->imgpel_abs_range, sizeof(int))) == NULL)
    no_mem_exit("Init_Motion_Search_Module: p_Vid->imgpel_abs");
  if ((p_Vid->imgpel_quad = (int*)calloc(p_Vid->imgpel_abs_range, sizeof(int))) == NULL)
    no_mem_exit("Init_Motion_Search_Module: p_Vid->imgpel_quad");
  p_Vid->imgpel_abs  += p_Vid->imgpel_abs_range / 2;
  p_Vid->imgpel_quad += p_Vid->imgpel_abs_range / 2;
#endif

  if (p_Vid->max_num_references)
    get_mem4Ddistblk (&p_Vid->motion_cost, 8, 2, p_Vid->max_num_references, 4);

  //--- set array offsets ---
  p_Vid->mvbits      += max_mvd;

  //=====   INIT ARRAYS   =====
  //---------------------------
  //--- init array: motion vector bits ---
  p_Vid->mvbits[0] = 1;
  for (bits = 3; bits <= max_mv_bits; bits += 2)
  {
    i_max = (short) (1 << (bits >> 1));
    i_min = i_max >> 1;

    for (i = i_min; i < i_max; i++)
      p_Vid->mvbits[-i] = p_Vid->mvbits[i] = bits;
  }

  //--- init array: reference frame bits ---
  p_Vid->refbits[0] = 1;
  for (bits=3; bits<=max_ref_bits; bits+=2)
  {
    i_max = (short) (1 << ((bits >> 1) + 1)) - 1;
    i_min = i_max >> 1;

    for (i = i_min; i < i_max; i++)
      p_Vid->refbits[i] = bits;
  }

#if (JM_MEM_DISTORTION)
  //--- init array: absolute value ---
  p_Vid->imgpel_abs[0] = 0;

  for (i=1; i<p_Vid->imgpel_abs_range / 2; i++)
  {
    p_Vid->imgpel_abs[i] = p_Vid->imgpel_abs[-i] = i;
  }

  //--- init array: square value ---
  p_Vid->imgpel_quad[0] = 0;

  for (i=1; i<p_Vid->imgpel_abs_range / 2; i++)
  {
    p_Vid->imgpel_quad[i] = p_Vid->imgpel_quad[-i] = i * i;
  }
#endif

  //--- init array: search pattern ---
  p_Vid->spiral_search[0].mv_x = p_Vid->spiral_search[0].mv_y = 0;
  p_Vid->spiral_hpel_search[0].mv_x = p_Vid->spiral_hpel_search[0].mv_y = 0;
  p_Vid->spiral_qpel_search[0].mv_x = p_Vid->spiral_qpel_search[0].mv_y = 0;

  for (k=1, l=1; l <= imax(1,search_range); l++)
  {
    for (i=-l+1; i< l; i++)
    {
      p_Vid->spiral_search[k].mv_x =     (short)  i;
      p_Vid->spiral_search[k].mv_y =     (short) -l;
      p_Vid->spiral_hpel_search[k].mv_x =   (short) (i<<1);
      p_Vid->spiral_hpel_search[k].mv_y =   (short) -(l<<1);
      p_Vid->spiral_qpel_search[k].mv_x =   (short)  (i<<2);
      p_Vid->spiral_qpel_search[k++].mv_y = (short) -(l<<2);
      p_Vid->spiral_search[k].mv_x =     (short)  i;
      p_Vid->spiral_search[k].mv_y =     (short)  l;
      p_Vid->spiral_hpel_search[k].mv_x =   (short) (i<<1);
      p_Vid->spiral_hpel_search[k].mv_y =   (short) (l<<1);
      p_Vid->spiral_qpel_search[k].mv_x =   (short) (i<<2);
      p_Vid->spiral_qpel_search[k++].mv_y = (short) (l<<2);
    }
    for (i=-l;   i<=l; i++)
    {
      p_Vid->spiral_search[k].mv_x =     (short) -l;
      p_Vid->spiral_search[k].mv_y =     (short)  i;
      p_Vid->spiral_hpel_search[k].mv_x =   (short) -(l<<1);
      p_Vid->spiral_hpel_search[k].mv_y =   (short)  (i<<1);
      p_Vid->spiral_qpel_search[k].mv_x =   (short) -(l<<2);
      p_Vid->spiral_qpel_search[k++].mv_y = (short)  (i<<2);
      p_Vid->spiral_search[k].mv_x =     (short)  l;
      p_Vid->spiral_search[k].mv_y =     (short)  i;
      p_Vid->spiral_hpel_search[k].mv_x =   (short) (l<<1);
      p_Vid->spiral_hpel_search[k].mv_y =   (short) (i<<1);
      p_Vid->spiral_qpel_search[k].mv_x =   (short) (l<<2);
      p_Vid->spiral_qpel_search[k++].mv_y = (short) (i<<2);
    }
  }

  // set global variable prior to ME
  p_Vid->start_me_refinement_hp = (p_Inp->ChromaMEEnable == 1 || p_Inp->MEErrorMetric[F_PEL] != p_Inp->MEErrorMetric[H_PEL] ) ? 0 : 1;
  p_Vid->start_me_refinement_qp = (p_Inp->ChromaMEEnable == 1 || p_Inp->MEErrorMetric[H_PEL] != p_Inp->MEErrorMetric[Q_PEL] ) ? 0 : 1;

  select_distortion(p_Vid, p_Inp);

  // Setup Distortion Metrics depending on refinement level
  for (i=0; i<3; i++)
  {
    switch(p_Inp->MEErrorMetric[i])
    {
    case ERROR_SAD:
      p_Vid->computeUniPred[i] = computeSAD;
      p_Vid->computeUniPred[i + 3] = computeSADWP;
      p_Vid->computeBiPred1[i] = computeBiPredSAD1;
      p_Vid->computeBiPred2[i] = computeBiPredSAD2;
      break;
    case ERROR_SSE:
      p_Vid->computeUniPred[i] = computeSSE;
      p_Vid->computeUniPred[i + 3] = computeSSEWP;
      p_Vid->computeBiPred1[i] = computeBiPredSSE1;
      p_Vid->computeBiPred2[i] = computeBiPredSSE2;
      break;
    case ERROR_SATD :
    default:
      p_Vid->computeUniPred[i] = computeSATD;
      p_Vid->computeUniPred[i + 3] = computeSATDWP;
      p_Vid->computeBiPred1[i] = computeBiPredSATD1;
      p_Vid->computeBiPred2[i] = computeBiPredSATD2;
      break;
    }
  }
  if (!p_Inp->IntraProfile)
  {
    if(p_Inp->SearchMode == FAST_FULL_SEARCH)
      InitializeFastFullIntegerSearch (p_Vid, p_Inp);

    if (p_Inp->SearchMode == UM_HEX)
      UMHEX_DefineThreshold(p_Vid);
  }
}

/*!
 ************************************************************************
 * \brief
 *    Free memory used by motion search
 ************************************************************************
 */
void
Clear_Motion_Search_Module (VideoParameters *p_Vid, InputParameters *p_Inp)
{
  //int search_range               = p_Inp->search_range;
  //int number_of_subpel_positions = 4 * (2*search_range+3);
  //int max_mv_bits                = 3 + 2 * (int)ceil (log(number_of_subpel_positions + 1) / log(2) + 1e-10);
  int max_mvd                    = p_Vid->max_mvd; //(1<<( max_mv_bits >>1)   ) - 1;


  //--- correct array offset ---
  p_Vid->mvbits      -= max_mvd;
#if (JM_MEM_DISTORTION)
  p_Vid->imgpel_abs  -= p_Vid->imgpel_abs_range / 2;
  p_Vid->imgpel_quad -= p_Vid->imgpel_abs_range / 2;
#endif

  //--- delete arrays ---
  free (p_Vid->spiral_search);
  free (p_Vid->spiral_hpel_search);
  free (p_Vid->spiral_qpel_search);
  free (p_Vid->mvbits);
  free (p_Vid->refbits);

#if (JM_MEM_DISTORTION)
  free (p_Vid->imgpel_abs);
  free (p_Vid->imgpel_quad);
#endif

  if (p_Vid->motion_cost)
    free_mem4Ddistblk (p_Vid->motion_cost);

  if ((p_Inp->SearchMode == FAST_FULL_SEARCH) && (!p_Inp->IntraProfile) )
    ClearFastFullIntegerSearch (p_Vid);
}
static inline int mv_bits_cost(VideoParameters *p_Vid, short ***all_mv, short ***p_mv, int by, int bx, int step_v0, int step_v, int step_h0, int step_h, int mvd_bits)
{
  int v, h;
  for (v=by; v<by + step_v0; v+=step_v)
  {
    for (h=bx; h<bx + step_h0; h+=step_h)
    {
      mvd_bits += (int) p_Vid->mvbits[ all_mv[v][h][0] - p_mv[v][h][0] ];
      mvd_bits += (int) p_Vid->mvbits[ all_mv[v][h][1] - p_mv[v][h][1] ];
    }
  }
  return mvd_bits;
}

static inline int mv_bit_cost(Macroblock *currMB, short ***all_mv, int cur_list, short cur_ref, int by, int bx, int step_v0, int step_v, int step_h0, int step_h, int mvd_bits)
{
  int v, h;
  short predMV[2]; 
  PixelPos block[4];  // neighbor blocks
  VideoParameters *p_Vid = currMB->p_Vid;
  PicMotionParams *motion = &p_Vid->enc_picture->motion;

  for (v=by; v<by + step_v0; v+=step_v)
  {
    for (h=bx; h<bx + step_h0; h+=step_h)
    {

      get_neighbors(currMB, block, h, v, step_h);
      // Lets recompute MV predictor. This should avoid any problems with alterations of the motion vectors after ME
      currMB->GetMVPredictor (currMB, block, predMV, cur_ref, motion->ref_idx[cur_list], motion->mv[cur_list], h, v, step_h, step_v);

      mvd_bits += p_Vid->mvbits[ all_mv[v][h][0] - predMV[0] ];
      mvd_bits += p_Vid->mvbits[ all_mv[v][h][1] - predMV[1] ];
    }
  }

  return mvd_bits;
}

/*!
 ***********************************************************************
 * \brief
 *    Motion Cost for Bidirectional modes
 ***********************************************************************
 */
distblk BPredPartitionCost (Macroblock *currMB,
                        int   blocktype,
                        int   block8x8,
                        short ref_l0,
                        short ref_l1,
                        int   lambda_factor,
                        int   list)
{
  VideoParameters *p_Vid = currMB->p_Vid;
  Slice *currSlice = currMB->p_slice;
  imgpel **cur_img = p_Vid->pCurImg;

  int   curr_blk[MB_BLOCK_SIZE][MB_BLOCK_SIZE]; // ABT pred.error buffer
  int bsx = (short) imin(block_size[blocktype][0], 8);
  int bsy = (short) imin(block_size[blocktype][1], 8);

  short pic_pix_x, pic_pix_y;
  short  v, h;
  int i, j, k;
  distblk mcost;

  int   mvd_bits  = 0;
  int diff64[64];


  short parttype  = (short) (blocktype < 4 ? blocktype : 4);
  short step_h0   = (part_size[ parttype][0]);
  short step_v0   = (part_size[ parttype][1]);
  short step_h    = (part_size[blocktype][0]);
  short step_v    = (part_size[blocktype][1]);
  int   bxx, byy;                               // indexing curr_blk
  short by0_part = by0[parttype][block8x8];
  short bx0_part = bx0[parttype][block8x8];

  short   ***all_mv_l0 = currSlice->bipred_mv[list][LIST_0][ref_l0][blocktype]; 
  short   ***all_mv_l1 = currSlice->bipred_mv[list][LIST_1][ref_l1][blocktype]; 
  imgpel  **mb_pred    = currSlice->mb_pred[0];

  // List0 
  mvd_bits = mv_bit_cost(currMB, all_mv_l0, LIST_0, ref_l0, by0_part, bx0_part, step_v0, step_v, step_h0, step_h, mvd_bits);
  // List1
  mvd_bits = mv_bit_cost(currMB, all_mv_l1, LIST_1, ref_l1, by0_part, bx0_part, step_v0, step_v, step_h0, step_h, mvd_bits);

  mcost = weighted_cost (lambda_factor, mvd_bits);

  //----- cost of residual signal -----
  if ((!currSlice->p_Inp->Transform8x8Mode) || (blocktype>4))
  {
    for (byy=0, v = by0_part << 2; v < (by0_part + step_v0) << 2; byy += 4, v += 4)
    {

      pic_pix_y = currMB->opix_y + v;
      for (bxx=0, h = (bx0_part << 2); h < (bx0_part + step_h0) << 2; bxx += 4, h += 4)
      {
        pic_pix_x = currMB->pix_x + h;
        luma_prediction_bi (currMB, h, v, 4, 4, blocktype, blocktype, ref_l0, ref_l1, list);

        for (k = j = 0; j < 4; j++)
        {
          for (i = 0; i < 4; i++)
            diff64[k++] = cur_img[pic_pix_y+j][pic_pix_x+i] - mb_pred[j + v][i + h];
        }
        mcost += p_Vid->distortion4x4 (diff64, DISTBLK_MAX);
      }
    }
  }
  else
  {
    for (byy=0, v = by0_part << 2; v < (by0_part + step_v0) << 2; byy += 4, v += 4)
    {

      pic_pix_y = currMB->opix_y + v;
      for (bxx=0, h = (bx0_part << 2); h < (bx0_part + step_h0) << 2; bxx += 4, h += 4)
      {
        pic_pix_x = currMB->pix_x + h;
        luma_prediction_bi (currMB, h, v, 4, 4, blocktype, blocktype, ref_l0, ref_l1, list);

        for (k = j = 0; j < 4; j++)
        {
          for (i = 0; i < 4; i++)
            curr_blk[byy+j][bxx+i] = cur_img[pic_pix_y+j][pic_pix_x+i] - mb_pred[j + v][i + h];
        }
      }
    }

    for (byy=0; byy < block_size[parttype][1]; byy += bsy)
    {
      for (bxx=0; bxx < block_size[parttype][0]; bxx += bsx)
      {
        for (k=0, j = byy; j < byy + 8; j++, k += 8)
          memcpy(&diff64[k], &(curr_blk[j][bxx]), 8 * sizeof(int));

        mcost += p_Vid->distortion8x8(diff64, DISTBLK_MAX);
      }
    }
  }
  return mcost;
}

void update_mv_block(Macroblock *currMB, MEBlock *mv_block, int h, int v)
{
  mv_block->block_x      = (short) h;
  mv_block->block_y      = (short) v;
  mv_block->pos_x        = (short) (currMB->pix_x  + (h << 2));
  mv_block->pos_y        = (short) (currMB->opix_y + (v << 2));
  mv_block->pos_x2       = (short) (mv_block->pos_x >> 2);
  mv_block->pos_y2       = (short) (mv_block->pos_y >> 2);
#if (PAD_AFTER)
  mv_block->pos_x_padded = (short) (mv_block->pos_x << 2);
  mv_block->pos_y_padded = (short) (mv_block->pos_y << 2);
#else
  mv_block->pos_x_padded = (short) (mv_block->pos_x << 2) + IMG_PAD_SIZE_TIMES4;
  mv_block->pos_y_padded = (short) (mv_block->pos_y << 2) + IMG_PAD_SIZE_TIMES4;
#endif
  mv_block->pos_cr_x     = (short) (mv_block->pos_x >> currMB->p_Vid->shift_cr_x);
  mv_block->pos_cr_y     = (short) (mv_block->pos_y >> currMB->p_Vid->shift_cr_y);
}

/*!
 ***********************************************************************
 * \brief
 *    Init motion vector block
 ***********************************************************************
 */
void init_mv_block(Macroblock *currMB, MEBlock *mv_block, short blocktype, int list, char ref_idx, short mb_x, short mb_y)
{
  InputParameters *p_Inp = currMB->p_Inp;
  VideoParameters *p_Vid = currMB->p_Vid;
  Slice *currSlice = currMB->p_slice;
  mv_block->blocktype         = blocktype;
  mv_block->blocksize_x       = block_size[blocktype][0];  // horizontal block size
  mv_block->blocksize_y       = block_size[blocktype][1];  // vertical block size
  // update position info
  update_mv_block(currMB, mv_block, mb_x, mb_y);

  mv_block->list              = (char) list;
  mv_block->ref_idx           = ref_idx;

  mv_block->mv[LIST_0].mv_x   = 0;
  mv_block->mv[LIST_0].mv_y   = 0;
  mv_block->mv[LIST_1].mv_x   = 0;
  mv_block->mv[LIST_1].mv_y   = 0;
  // Init WP parameters
  mv_block->p_Vid             = p_Vid;
  mv_block->p_slice           = currSlice;
  mv_block->cost              = INT_MAX;
  mv_block->search_pos2       = 9;
  mv_block->search_pos4       = 9;

  if ((mv_block->orig_pic = (imgpel**)calloc(3, sizeof(imgpel *))) == NULL)
    no_mem_exit("init_mv_block: mv_block->orig_pic");

  get_mem1Dpel(&(mv_block->orig_pic[0]), mv_block->blocksize_x * mv_block->blocksize_y);
  
  mv_block->ChromaMEEnable = p_Inp->ChromaMEEnable;

  mv_block->apply_bi_weights = p_Inp->UseWeightedReferenceME && ((currSlice->slice_type == B_SLICE) && p_Vid->active_pps->weighted_bipred_idc != 0);
  mv_block->apply_weights    = p_Inp->UseWeightedReferenceME && ( currSlice->weighted_prediction != 0 );

  if (p_Inp->ChromaMEEnable)
  {    
    mv_block->blocksize_cr_x = (short) (mv_block->blocksize_x >> p_Vid->shift_cr_x);
    mv_block->blocksize_cr_y = (short) (mv_block->blocksize_y >> p_Vid->shift_cr_y);

    mv_block->ChromaMEWeight = p_Inp->ChromaMEWeight;
    get_mem1Dpel(&(mv_block->orig_pic[1]), mv_block->blocksize_cr_x * mv_block->blocksize_cr_y);
    get_mem1Dpel(&(mv_block->orig_pic[2]), mv_block->blocksize_cr_x * mv_block->blocksize_cr_y);
  }

  if (mv_block->apply_weights)
  {
    mv_block->computePredFPel   = p_Vid->computeUniPred[F_PEL + 3];
    mv_block->computePredHPel   = p_Vid->computeUniPred[H_PEL + 3];
    mv_block->computePredQPel   = p_Vid->computeUniPred[Q_PEL + 3];
    mv_block->computeBiPredFPel = p_Vid->computeBiPred2[F_PEL];
    mv_block->computeBiPredHPel = p_Vid->computeBiPred2[H_PEL];
    mv_block->computeBiPredQPel = p_Vid->computeBiPred2[Q_PEL];    
  }
  else
  {
    mv_block->computePredFPel   = p_Vid->computeUniPred[F_PEL];
    mv_block->computePredHPel   = p_Vid->computeUniPred[H_PEL];
    mv_block->computePredQPel   = p_Vid->computeUniPred[Q_PEL];
    mv_block->computeBiPredFPel = p_Vid->computeBiPred1[F_PEL];
    mv_block->computeBiPredHPel = p_Vid->computeBiPred1[H_PEL];
    mv_block->computeBiPredQPel = p_Vid->computeBiPred1[Q_PEL];
  }
}

/*!
 ***********************************************************************
 * \brief
 *    free motion vector block
 ***********************************************************************
 */
void free_mv_block(InputParameters *p_Inp, MEBlock *mv_block)
{
  if (mv_block->orig_pic)
  {
    free_mem1Dpel(mv_block->orig_pic[0]);
    if (p_Inp->ChromaMEEnable)
    {
      free_mem1Dpel(mv_block->orig_pic[1]);
      free_mem1Dpel(mv_block->orig_pic[2]);
    }
    free(mv_block->orig_pic);
  }
}


void get_original_block(VideoParameters *p_Vid, MEBlock *mv_block)
{
  //==================================
  //=====   GET ORIGINAL BLOCK   =====
  //==================================
  imgpel *orig_pic_tmp = mv_block->orig_pic[0];
  int   bsx       = mv_block->blocksize_x;
  int   pic_pix_x = mv_block->pos_x;
  int   i, j;
  imgpel **cur_img = &p_Vid->pCurImg[mv_block->pos_y];

  for (j = 0; j < mv_block->blocksize_y; j++)
  {
    memcpy(orig_pic_tmp,&cur_img[j][pic_pix_x], bsx * sizeof(imgpel));
    orig_pic_tmp += bsx;
  }

  if ( p_Vid->p_Inp->ChromaMEEnable )
  {
    bsx       = mv_block->blocksize_cr_x;
    pic_pix_x = mv_block->pos_cr_x;

    // copy the original cmp1 and cmp2 data to the orig_pic matrix
    for ( i = 1; i<=2; i++)
    {
      cur_img = &p_Vid->pImgOrg[i][mv_block->pos_cr_y];
      orig_pic_tmp = mv_block->orig_pic[i];
      for (j = 0; j < mv_block->blocksize_cr_y; j++)
      {
        memcpy(orig_pic_tmp, &(cur_img[j][pic_pix_x]), bsx * sizeof(imgpel));
        orig_pic_tmp += bsx;
      }
    }
  }
}

void CheckSearchRange(VideoParameters *p_Vid, MotionVector *pPredMV, MotionVector *pSWC, MEBlock *mv_block)
{
   int iMaxMVD = p_Vid->max_mvd-2;
   int left, right, top, down;
   left = pSWC->mv_x+mv_block->searchRange.min_x;
   right = pSWC->mv_x+mv_block->searchRange.max_x;
   top = pSWC->mv_y+mv_block->searchRange.min_y;
   down = pSWC->mv_y+mv_block->searchRange.max_y;
   //left;
   if(left < pPredMV->mv_x-iMaxMVD)
     left =pPredMV->mv_x-iMaxMVD;
   else if(left > pPredMV->mv_x+iMaxMVD)
     left =pPredMV->mv_x+iMaxMVD;
   //right;
   if(right < pPredMV->mv_x-iMaxMVD)
     right =pPredMV->mv_x-iMaxMVD;
   else if(right > pPredMV->mv_x+iMaxMVD)
     right =pPredMV->mv_x+iMaxMVD;
   
   //top;
   if(top < pPredMV->mv_y-iMaxMVD)
     top =pPredMV->mv_y-iMaxMVD;
   else if(top > pPredMV->mv_y+iMaxMVD)
     top =pPredMV->mv_y+iMaxMVD;
   //down;
   if(down < pPredMV->mv_y-iMaxMVD)
     down =pPredMV->mv_y-iMaxMVD;
   else if(down > pPredMV->mv_y+iMaxMVD)
     down =pPredMV->mv_y+iMaxMVD;
   
   if(left<right && top<down)
   {
     pSWC->mv_x = (short) ((left + right)>>1);
     pSWC->mv_y = (short) ((top + down)>>1);
     mv_block->searchRange.min_x = left - pSWC->mv_x;
     mv_block->searchRange.max_x = imin(pSWC->mv_x-left, right-pSWC->mv_x);
     mv_block->searchRange.min_y = top - pSWC->mv_y;
     mv_block->searchRange.max_y = imin(pSWC->mv_y-top, down-pSWC->mv_y);
   }
   else
   {
      *pSWC = *pPredMV;
   }
}
/*!
 ***********************************************************************
 * \brief
 *    Block motion search
 ***********************************************************************
 */
distblk                                         //!< minimum motion cost after search
BlockMotionSearch (Macroblock *currMB,      //!< Current Macroblock
                   MEBlock   *mv_block,     //!< Motion estimation information block
                   int       mb_x,          //!< x-coordinate inside macroblock
                   int       mb_y,          //!< y-coordinate inside macroblock
                   int*      lambda_factor) //!< lagrangian parameter for determining motion cost
{
  // each 48-pel line stores the 16 luma pels (at 0) followed by 8 or 16 crcb[0] (at 16) and crcb[1] (at 32) pels
  // depending on the type of chroma subsampling used: YUV 4:4:4, 4:2:2, and 4:2:0
  Slice *currSlice = currMB->p_slice;
  VideoParameters *p_Vid = currMB->p_Vid;
  InputParameters *p_Inp = currMB->p_Inp;

  short pred_mv[2];
  int   i, j;
  distblk   max_value = DISTBLK_MAX;
  distblk   min_mcost = max_value;
  int   block_x   = (mb_x>>2);
  int   block_y   = (mb_y>>2);

  int   bsx       = mv_block->blocksize_x;
  int   bsy       = mv_block->blocksize_y;

  short pic_pix_x = (short) (currMB->pix_x + mb_x);

  int  blocktype = mv_block->blocktype;
  int  list = mv_block->list;
  short ref = mv_block->ref_idx;
  MotionVector *mv = &mv_block->mv[list], pred; 

  short***   all_mv = &currSlice->all_mv[list][ref][blocktype][block_y];
  PicMotionParams *motion = &p_Vid->enc_picture->motion;
  distblk *prevSad = (p_Inp->SearchMode == EPZS)? currSlice->p_EPZS->distortion[list + currMB->list_offset][blocktype - 1]: NULL;

  get_neighbors(currMB, mv_block->block, mb_x, mb_y, bsx);



  PrepareMEParams(currSlice, mv_block, p_Inp->ChromaMEEnable, list + currMB->list_offset, ref);

  //==================================
  //=====   GET ORIGINAL BLOCK   =====
  //==================================
  if (blocktype > 4)
    get_original_block(p_Vid, mv_block);

  //===========================================
  //=====   GET MOTION VECTOR PREDICTOR   =====
  //===========================================
  if (p_Inp->SearchMode == UM_HEX)
  {
    p_Vid->p_UMHex->UMHEX_blocktype = blocktype;
    p_Vid->p_UMHex->bipred_flag = 0;
    UMHEXSetMotionVectorPredictor(currMB, pred_mv, motion->ref_idx[list], motion->mv[list], ref, list, mb_x, mb_y, bsx, bsy, mv_block);
  }
  else if (p_Inp->SearchMode == UM_HEX_SIMPLE)
  {
    smpUMHEX_setup(currMB, ref, list, block_y, block_x, blocktype, currSlice->all_mv );
    currMB->GetMVPredictor (currMB, mv_block->block, pred_mv, ref, motion->ref_idx[list], motion->mv[list], mb_x, mb_y, bsx, bsy);
  }
  else
  {
    currMB->GetMVPredictor (currMB, mv_block->block, pred_mv, ref, motion->ref_idx[list], motion->mv[list], mb_x, mb_y, bsx, bsy);
  }

  pred.mv_x = pred_mv[0];
  pred.mv_y = pred_mv[1];

  //==================================
  //=====   INTEGER-PEL SEARCH   =====
  //==================================
  if (p_Inp->EPZSSubPelGrid)
  {
    *mv = pred;
  }
  else
  {
#if (JM_INT_DIVIDE)
    mv->mv_x = (short) (((pred.mv_x  + 2) >> 2) * 4);
    mv->mv_y = (short) (((pred.mv_y  + 2) >> 2) * 4);
#else
    mv->mv_x = (short) ((pred.mv_x / 4) * 4);
    mv->mv_y = (short) ((pred.mv_y / 4) * 4);
#endif
  }

  if (!p_Inp->rdopt)
  {
    MotionVector center = *mv;
    //--- adjust search center so that the (0,0)-vector is inside ---
    mv->mv_x = (short) iClip3 (mv_block->searchRange.min_x, mv_block->searchRange.max_x, mv->mv_x);
    mv->mv_y = (short) iClip3 (mv_block->searchRange.min_y, mv_block->searchRange.max_y, mv->mv_y);
    //mvbits overflow checking;
    if((mv->mv_x != center.mv_x) || (mv->mv_y != center.mv_y))
      CheckSearchRange(p_Vid, &center, mv, mv_block);
  }

  // valid search range limits could be precomputed once during the initialization process
  clip_mv_range(p_Vid, 0, mv, Q_PEL);

  //--- perform motion search ---
  min_mcost = currMB->IntPelME (currMB, &pred, mv_block, min_mcost, lambda_factor[F_PEL]);

  //==============================
  //=====   SUB-PEL SEARCH   =====
  //============================== 
  mv_block->ChromaMEEnable = (p_Inp->ChromaMEEnable == ME_YUV_FP_SP ) ? 1 : 0; // set it externally

  if (!p_Inp->DisableSubpelME)
  {
    if (p_Inp->SearchMode != EPZS || (ref == 0 || currSlice->structure != FRAME || (ref > 0 && min_mcost < 3.5 * prevSad[pic_pix_x >> 2])))
    {
      if ( !p_Vid->start_me_refinement_hp )
      {
        min_mcost = max_value;
      }
      min_mcost =  currMB->SubPelME (currMB, &pred, mv_block, min_mcost, lambda_factor);
    }
  }

  // clip mvs after me is performed (is not exactly the best)
  // better solution is to modify search window appropriately
  clip_mv_range(p_Vid, 0, mv, Q_PEL);

  if (!p_Inp->rdopt)
  {
    // Get the skip mode cost
    if (blocktype == 1 && (currSlice->slice_type == P_SLICE|| (currSlice->slice_type == SP_SLICE) ))
    {
      distblk cost;
      FindSkipModeMotionVector (currMB);

      cost  = GetSkipCostMB (currMB, lambda_factor[Q_PEL]);
      if (cost < min_mcost)
      {
        min_mcost = cost;
        mv->mv_x = currSlice->all_mv [0][0][0][0][0][0];
        mv->mv_y = currSlice->all_mv [0][0][0][0][0][1];
      }
    } 
  }

  //===============================================
  //=====   SET MV'S AND RETURN MOTION COST   =====
  //===============================================

  // Set first line
  for (i=block_x; i < block_x + (bsx>>2); i++)
  {
    all_mv[0][i][0] = mv->mv_x;
    all_mv[0][i][1] = mv->mv_y;
  }

  // set all other lines
  for (j=1; j < (bsy>>2); j++)
  {
    memcpy(all_mv[j][block_x], all_mv[0][block_x], (bsx>>2) * 2 * sizeof(short));
  }


  // Bipred ME consideration: returns minimum bipred cost
  if (currSlice->slice_type == B_SLICE && is_bipred_enabled(p_Inp, blocktype) && (ref == 0)) 
  {
    BiPredBlockMotionSearch(currMB, mv_block, &pred, mb_x, mb_y, lambda_factor);
  }

  return min_mcost;
}


/*!
 ***********************************************************************
 * \brief
 *    Bi-predictive motion search
 ***********************************************************************
 */
static distblk BiPredBlockMotionSearch(Macroblock *currMB,      //!< Current Macroblock
                                   MEBlock  *mv_block,
                                   MotionVector *pred_mv,     //!< current list motion vector predictor
                                   int       mb_x,            //!< x-coordinate inside macroblock
                                   int       mb_y,            //!< y-coordinate inside macroblock
                                   int*      lambda_factor)   //!< lagrangian parameter for determining motion cost
{
  VideoParameters *p_Vid     = currMB->p_Vid;
  InputParameters *p_Inp     = currMB->p_Inp;
  Slice           *currSlice = currMB->p_slice;
  int         list = mv_block->list;
  int         i, j;
  short       bipred_type = list ? 0 : 1;
  short****** bipred_mv = currSlice->bipred_mv[bipred_type];
  distblk     min_mcostbi = DISTBLK_MAX;
  MotionVector *mv = &mv_block->mv[list];
  MotionVector bimv, tempmv;
  MotionVector pred_mv1, pred_mv2, pred_bi;
  MotionVector bi_mv1 = { 0, 0}, bi_mv2 = { 0, 0};
  short       iterlist = (short) list;
  short       pred_mv_bi[2];
  int         block_x   = (mb_x>>2);
  int         block_y   = (mb_y>>2);
  int  blocktype = mv_block->blocktype;
  int         bsx       = mv_block->blocksize_x;
  int         bsy       = mv_block->blocksize_y;
  //PixelPos    block[4];  // neighbor blocks
  PicMotionParams *motion = &p_Vid->enc_picture->motion;
  
  //get_neighbors(currMB, mv_block->block, mb_x, mb_y, bsx);

  if (p_Inp->SearchMode == UM_HEX)
  {
    p_Vid->p_UMHex->bipred_flag = 1;
    UMHEXSetMotionVectorPredictor(currMB, pred_mv_bi, motion->ref_idx[list ^ 1], motion->mv[list ^ 1], 0, list ^ 1, mb_x, mb_y, bsx, bsy, mv_block);
  }
  else
    currMB->GetMVPredictor (currMB, mv_block->block, pred_mv_bi, 0, motion->ref_idx[list ^ 1], motion->mv[list ^ 1], mb_x, mb_y, bsx, bsy);

  pred_bi.mv_x = pred_mv_bi[0];
  pred_bi.mv_y = pred_mv_bi[1];

  if ((p_Inp->SearchMode != EPZS) || (p_Inp->EPZSSubPelGrid == 0))
  {
    mv->mv_x = ((mv->mv_x  + 2) >> 2) * 4;
    mv->mv_y = ((mv->mv_y  + 2) >> 2) * 4;
    bimv.mv_x = ((pred_bi.mv_x  + 2) >> 2) * 4;
    bimv.mv_y = ((pred_bi.mv_y  + 2) >> 2) * 4;
  }
  else
  {
    bimv = pred_bi;
  }

  //Bi-predictive motion Refinements
  for (mv_block->iteration_no = 0; mv_block->iteration_no <= p_Inp->BiPredMERefinements; mv_block->iteration_no++)
  {
    if (mv_block->iteration_no & 0x01)
    {
      pred_mv1  = *pred_mv;
      pred_mv2  = pred_bi;
      bi_mv1    = *mv;
      bi_mv2    = bimv;
      iterlist  = (short) list;
    }
    else
    {
      pred_mv1  = pred_bi;
      pred_mv2  = *pred_mv;
      bi_mv1    = bimv;
      bi_mv2    = *mv;
      iterlist = (short) (list ^ 1);
    }

    tempmv = bi_mv1;

    PrepareBiPredMEParams(currSlice, mv_block, mv_block->ChromaMEEnable, iterlist, currMB->list_offset, mv_block->ref_idx);
    // Get bipred mvs for list iterlist given previously computed mvs from other list
    min_mcostbi = currMB->BiPredME (currMB, iterlist, 
      &pred_mv1, &pred_mv2, &bi_mv1, &bi_mv2, mv_block, 
      (p_Inp->BiPredMESearchRange <<2)>>mv_block->iteration_no, min_mcostbi, lambda_factor[F_PEL]);

    if (mv_block->iteration_no > 0 && (tempmv.mv_x == bi_mv1.mv_x) && (tempmv.mv_y == bi_mv1.mv_y))
    {
      break;
    }
  }

  if (!p_Inp->DisableSubpelME)
  {
    if (p_Inp->BiPredMESubPel)
    {
      min_mcostbi = DISTBLK_MAX;
      PrepareBiPredMEParams(currSlice, mv_block, mv_block->ChromaMEEnable, iterlist, currMB->list_offset, mv_block->ref_idx);

      min_mcostbi =  currMB->SubPelBiPredME (currMB, mv_block, iterlist, &pred_mv1, &pred_mv2, &bi_mv1, &bi_mv2, min_mcostbi, lambda_factor);
    }

    if (p_Inp->BiPredMESubPel==2)
    {
      min_mcostbi = DISTBLK_MAX;
      PrepareBiPredMEParams(currSlice, mv_block, mv_block->ChromaMEEnable, iterlist ^ 1, currMB->list_offset, mv_block->ref_idx);

      min_mcostbi =  currMB->SubPelBiPredME (currMB, mv_block, iterlist ^ 1, &pred_mv2, &pred_mv1, &bi_mv2, &bi_mv1, min_mcostbi, lambda_factor);
    }
  }

  clip_mv_range(p_Vid, 0, &bi_mv1, Q_PEL);
  clip_mv_range(p_Vid, 0, &bi_mv2, Q_PEL);

  for (j=block_y; j < block_y + (bsy>>2); j++)
  {
    for (i=block_x ; i < block_x + (bsx>>2); i++)
    {
      bipred_mv[iterlist    ][(short) mv_block->ref_idx][blocktype][j][i][0] = bi_mv1.mv_x;
      bipred_mv[iterlist    ][(short) mv_block->ref_idx][blocktype][j][i][1] = bi_mv1.mv_y;
      bipred_mv[iterlist ^ 1][(short) mv_block->ref_idx][blocktype][j][i][0] = bi_mv2.mv_x;
      bipred_mv[iterlist ^ 1][(short) mv_block->ref_idx][blocktype][j][i][1] = bi_mv2.mv_y;
    }
  }
  return min_mcostbi;
}

/*!
 ***********************************************************************
 * \brief
 *    Motion Cost for Bidirectional modes
 ***********************************************************************
 */
distblk BIDPartitionCost (Macroblock *currMB, 
                      int   blocktype,
                      int   block8x8,
                      char  cur_ref[2],
                      int   lambda_factor)
{
  VideoParameters *p_Vid = currMB->p_Vid;
  Slice *currSlice = currMB->p_slice;
  imgpel **cur_img = p_Vid->pCurImg;

  int   curr_blk[MB_BLOCK_SIZE][MB_BLOCK_SIZE]; // ABT pred.error buffer
  int   bsx       = imin(block_size[blocktype][0],8);
  int   bsy       = imin(block_size[blocktype][1],8);

  short pic_pix_x, pic_pix_y, block_x, block_y;
  int   v, h, i, j, k;
  distblk mcost;

  int   mvd_bits  = 0;

  int   parttype  = (blocktype < 4 ? blocktype : 4);
  int   step_h0   = (part_size[ parttype][0]);
  int   step_v0   = (part_size[ parttype][1]);
  int   step_h    = (part_size[blocktype][0]);
  int   step_v    = (part_size[blocktype][1]);
  int   bxx, byy;                               // indexing curr_blk
  int   bx = bx0[parttype][block8x8];
  int   by = by0[parttype][block8x8];
  short   *** all_mv_l0 = currSlice->all_mv [LIST_0][(int) cur_ref[LIST_0]][blocktype];
  short   *** all_mv_l1 = currSlice->all_mv [LIST_1][(int) cur_ref[LIST_1]][blocktype];
  short bipred_me =  0; //no bipred for this case 
  imgpel  **mb_pred = currSlice->mb_pred[0];
  int diff64[64];

  int   list_mode[2];
  list_mode[0] = blocktype;
  list_mode[1] = blocktype;

  //----- cost for motion vector bits -----
  // Should write a separate, small function to do this processing
  // List0 
  // mvd_bits = mv_bits_cost(p_Vid, all_mv_l0, p_mv_l0, by, bx, step_v0, step_v, step_h0, step_h, mvd_bits);
  mvd_bits = mv_bit_cost(currMB, all_mv_l0, LIST_0, cur_ref[LIST_0], by, bx, step_v0, step_v, step_h0, step_h, mvd_bits);
  // List1
  // mvd_bits = mv_bits_cost(p_Vid, all_mv_l1, p_mv_l1, by, bx, step_v0, step_v, step_h0, step_h, mvd_bits);
  mvd_bits = mv_bit_cost(currMB, all_mv_l1, LIST_1, cur_ref[LIST_1], by, bx, step_v0, step_v, step_h0, step_h, mvd_bits);

  mcost = weighted_cost (lambda_factor, mvd_bits);

  //----- cost of residual signal -----
  if ((!currSlice->p_Inp->Transform8x8Mode) || (blocktype>4))
  {
    for (byy=0, v=by; v<by + step_v0; byy+=4, v++)
    {
      pic_pix_y = (short) (currMB->opix_y + (block_y = (short) (v<<2)));
      for (bxx=0, h=bx; h<bx + step_h0; bxx+=4, h++)
      {
        pic_pix_x = (short) (currMB->pix_x + (block_x = (short) (h<<2)));
        luma_prediction (currMB, block_x, block_y, 4, 4, 2, list_mode, cur_ref, bipred_me);

        for (k=j=0; j<4; j++)
        {
          for (  i=0; i<4; i++)
            diff64[k++] = curr_blk[byy+j][bxx+i] =
            cur_img[pic_pix_y+j][pic_pix_x+i] - mb_pred[j+block_y][i+block_x];
        }

        mcost += p_Vid->distortion4x4 (diff64, DISTBLK_MAX);
      }
    }
  }
  else
  {
    for (byy=0, v= (by << 2); v < (by + step_v0) << 2; byy += 4, v += 4)
    {
      pic_pix_y = (short) (currMB->opix_y + v);
      for (bxx=0, h = (bx << 2); h < (bx + step_h0) << 2; bxx+=4, h += 4)
      {
        pic_pix_x = (short) (currMB->pix_x + h);
        luma_prediction (currMB, h, v, 4, 4, 2, list_mode, cur_ref, bipred_me);

        for (k=j=0; j<4; j++)
        {
          for (  i=0; i<4; i++)
            diff64[k++] = curr_blk[byy+j][bxx+i] =
            cur_img[pic_pix_y+j][pic_pix_x+i] - mb_pred[j + v][i + h];
        }
      }
    }

    for (byy=0; byy < block_size[parttype][1]; byy+=bsy)
    {
      for (bxx=0; bxx<block_size[parttype][0]; bxx+=bsx)
      {
        for (k=0, j=byy;j<byy + 8;j++, k += 8)
          memcpy(&diff64[k], &(curr_blk[j][bxx]), 8 * sizeof(int));

        mcost += p_Vid->distortion8x8(diff64, DISTBLK_MAX);
      }
    }
  }
  return mcost;
}

/*!
 ************************************************************************
 * \brief
 *    Get cost for skip mode for an macroblock
 ************************************************************************
 */
static distblk GetSkipCostMB (Macroblock *currMB, int lambda)
{
  Slice *currSlice = currMB->p_slice;
  VideoParameters *p_Vid = currMB->p_Vid;
  InputParameters *p_Inp = currMB->p_Inp;
  int block_y, block_x, pic_pix_y, pic_pix_x, i, j, k;
  distblk cost = 0;
  int curr_diff[8][8];
  int mb_x, mb_y;
  int block;
  imgpel  **mb_pred = currSlice->mb_pred[0];
  char  cur_ref[2] = {0, 0};
  int   list_mode[2] = {0, 0};
  int diff  [16];
  int diff64[64];

  for(block = 0;block < 4;block++)
  {
    mb_y    = (block >>   1)<<3;
    mb_x    = (block & 0x01)<<3;
    for (block_y = mb_y; block_y < mb_y+8; block_y += 4)
    {
      pic_pix_y = currMB->opix_y + block_y;
      for (block_x = mb_x; block_x < mb_x + 8; block_x += 4)
      {
        pic_pix_x = currMB->pix_x + block_x;

        //===== prediction of 4x4 block =====
        luma_prediction (currMB, block_x, block_y, 4, 4, 0, list_mode, cur_ref, 0);

        //===== get displaced frame difference ======
        for (k = j = 0; j < 4; j++)
        {
          for (i = 0; i < 4; i++, k++)
          {
            diff[k] = curr_diff[block_y-mb_y+j][block_x-mb_x+i] = p_Vid->pCurImg[pic_pix_y+j][pic_pix_x+i] - mb_pred[j+block_y][i+block_x];
          }
        }

        if(!((p_Inp->rdopt == 0) && (p_Inp->Transform8x8Mode)))
          cost += p_Vid->distortion4x4 (diff, DISTBLK_MAX);
      }
    }

    if((p_Inp->rdopt == 0) && (p_Inp->Transform8x8Mode))
    {
      for(k=j=0; j<8; j++, k+=8)
        memcpy(&diff64[k], &(curr_diff[j]), 8 * sizeof(int));
      cost += p_Vid->distortion8x8 (diff64, DISTBLK_MAX);
    }
  }

  //cost -= ((lambda_factor[Q_PEL] + 4096) >> 13);
  cost -= weight_cost(lambda, 8);

  return cost;
}

/*!
 ************************************************************************
 * \brief
 *    Find motion vector for the Skip mode
 ************************************************************************
 */
void FindSkipModeMotionVector (Macroblock *currMB)
{
  Slice *currSlice = currMB->p_slice;
  VideoParameters *p_Vid = currMB->p_Vid;
  PicMotionParams *motion = &p_Vid->enc_picture->motion;
  int   bx, by;
  short ***all_mv = currSlice->all_mv[0][0][0];

  short pmv[2];

  int zeroMotionAbove;
  int zeroMotionLeft;
  PixelPos mb[4];
  int      a_mv_y = 0;
  int      a_ref_idx = 0;
  int      b_mv_y = 0;
  int      b_ref_idx = 0;
  short    ***mv = motion->mv[LIST_0];

  get_neighbors(currMB, mb, 0, 0, 16);

  if (mb[0].available)
  {
    a_mv_y    = mv[mb[0].pos_y][mb[0].pos_x][1];
    a_ref_idx = motion->ref_idx[LIST_0][mb[0].pos_y][mb[0].pos_x];

    if (currMB->mb_field && !p_Vid->mb_data[mb[0].mb_addr].mb_field)
    {
      a_mv_y    /=2;
      a_ref_idx *=2;
    }
    if (!currMB->mb_field && p_Vid->mb_data[mb[0].mb_addr].mb_field)
    {
      a_mv_y    *= 2;
      a_ref_idx >>=1;
    }
  }

  if (mb[1].available)
  {
    b_mv_y    = mv[mb[1].pos_y][mb[1].pos_x][1];
    b_ref_idx = motion->ref_idx[LIST_0][mb[1].pos_y][mb[1].pos_x];

    if (currMB->mb_field && !p_Vid->mb_data[mb[1].mb_addr].mb_field)
    {
      b_mv_y    /=2;
      b_ref_idx *=2;
    }
    if (!currMB->mb_field && p_Vid->mb_data[mb[1].mb_addr].mb_field)
    {
      b_mv_y    *=2;
      b_ref_idx >>=1;
    }
  }

  zeroMotionLeft  = !mb[0].available ? 1 : a_ref_idx==0 && mv[mb[0].pos_y][mb[0].pos_x][0]==0 && a_mv_y==0 ? 1 : 0;
  zeroMotionAbove = !mb[1].available ? 1 : b_ref_idx==0 && mv[mb[1].pos_y][mb[1].pos_x][0]==0 && b_mv_y==0 ? 1 : 0;

  if (zeroMotionAbove || zeroMotionLeft)
  {
    memset(all_mv [0][0], 0, 32 * sizeof(short)); // 4 * 4 * 2
  }
  else
  {
    currMB->GetMVPredictor (currMB, mb, pmv, 0, motion->ref_idx[LIST_0], mv, 0, 0, 16, 16);

    for (bx = 0;bx < 4;bx++)
    {
      memcpy(all_mv [0][bx], pmv, 2* sizeof(short));
    }

    for (by = 1;by < 4;by++)
      memcpy(all_mv [by][0], all_mv [0][0], 4 * 2* sizeof(short));
  }
}

/*!
 ************************************************************************
 * \brief
 *    Get cost for direct mode for an 8x8 block
 ************************************************************************
 */
distblk GetDirectCost8x8 (Macroblock *currMB, int block, distblk *cost8x8)
{
  Slice *currSlice = currMB->p_slice; 
  VideoParameters *p_Vid = currMB->p_Vid;
  InputParameters *p_Inp = currMB->p_Inp;
  int block_y, block_x, pic_pix_y, pic_pix_x, i, j, k;
  int curr_diff[8][8];
  distblk cost  = 0;
  int mb_y  = (block >> 1)<<3;
  int mb_x  = (block & 0x01)<<3;
  short bipred_me  = 0;
  imgpel  **mb_pred = currSlice->mb_pred[0];
  int   list_mode[2] = {0, 0};
  int diff  [16];
  int diff64[64];


  for (block_y=mb_y; block_y < mb_y + 8; block_y += 4)
  {
    pic_pix_y = currMB->opix_y + block_y;

    for (block_x=mb_x; block_x<mb_x+8; block_x+=4)
    {
      pic_pix_x = currMB->pix_x + block_x;

      if (currSlice->direct_pdir[pic_pix_y>>2][pic_pix_x>>2]<0)
      {
        *cost8x8=DISTBLK_MAX;
        return DISTBLK_MAX; //mode not allowed
      }

      //===== prediction of 4x4 block =====

      luma_prediction (currMB, block_x, block_y, 4, 4,
        currSlice->direct_pdir[pic_pix_y>>2][pic_pix_x>>2], list_mode,
        currSlice->direct_ref_idx[pic_pix_y>>2][pic_pix_x>>2], bipred_me);

      //===== get displaced frame difference ======
      for (k=j=0; j<4; j++)
        for (i=0; i<4; i++, k++)
        {
          diff[k] = curr_diff[block_y-mb_y+j][block_x - mb_x+i] =
            p_Vid->pCurImg[pic_pix_y+j][pic_pix_x+i] - mb_pred[j+block_y][i+block_x];
        }
        cost += p_Vid->distortion4x4 (diff, DISTBLK_MAX);
    }
  }

  if((p_Inp->rdopt == 0) && (p_Inp->Transform8x8Mode))
  {
    k=0;
    for(j=0; j<8; j++, k+=8)
      memcpy(&diff64[k], &(curr_diff[j]), 8 * sizeof(int));          

    *cost8x8 += p_Vid->distortion8x8 (diff64, DISTBLK_MAX);
  }

  return cost;
}



/*!
 ************************************************************************
 * \brief
 *    Get cost for direct mode for an macroblock
 ************************************************************************
 */
distblk GetDirectCostMB (Macroblock *currMB)
{
  Slice *currSlice = currMB->p_slice; 
  InputParameters *p_Inp = currSlice->p_Inp;
  int i;
  distblk cost = 0;
  distblk cost8x8 = 0;
  int bslice = currSlice->slice_type == B_SLICE;

  for (i=0; i<4; i++)
  {
    cost += GetDirectCost8x8 (currMB, i, &cost8x8);
    if (cost8x8 == DISTBLK_MAX) return DISTBLK_MAX;
  }

  switch(p_Inp->Transform8x8Mode)
  {
  case 1: // Mixture of 8x8 & 4x4 transform
    if((cost8x8 < cost)||
      !(p_Inp->InterSearch[bslice][5] &&
      p_Inp->InterSearch[bslice][6] &&
      p_Inp->InterSearch[bslice][7])
      )
    {
      cost = cost8x8; //return 8x8 cost
    }
    break;
  case 2: // 8x8 Transform only
    cost = cost8x8;
    break;
  default: // 4x4 Transform only
    break;
  }

  return cost;
}


/*!
 ************************************************************************
 * \brief
 *    Motion search for a macroblock partition
 ************************************************************************
 */
void PartitionMotionSearch (Macroblock *currMB,
                            int    blocktype,
                            int    block8x8,
                            int    *lambda_factor)
{
  VideoParameters *p_Vid = currMB->p_Vid;
  InputParameters *p_Inp = currMB->p_Inp;
  Slice *currSlice = currMB->p_slice;
  PicMotionParams *motion = &p_Vid->enc_picture->motion;

  char  **ref_array;
  short ***mv_array;
  short ref = 0;
  int   step_h    = (part_size[blocktype][0]);
  int   step_v    = (part_size[blocktype][1]);
  int   list = LIST_0;
  int   numlists  = (currSlice->slice_type == B_SLICE) ? 2 : 1;
  int   list_offset = currMB->list_offset;
  distblk *m_cost;
  short   by = by0[blocktype][block8x8];
  short bx = bx0[blocktype][block8x8];
  short pic_block_y = currMB->block_y + by;
  short pic_block_x = currMB->block_x + bx;
  MEBlock  mv_block;
  //int ref_lambda = (p_Inp->rdopt) ? lambda_factor[Q_PEL] :  lambda_factor[Q_PEL] >> 2;

#if GET_METIME
  TIME_T me_time_start;
  TIME_T me_time_end;
  int64 me_tmp_time;
  gettime( &me_time_start );    // start time ms
#endif

  if (p_Vid->Motion_Selected == 1)
  {
    //===== LOOP OVER REFERENCE FRAMES =====
    for (list=0; list<numlists;list++)
    {
      //----- set arrays -----
      ref_array = &motion->ref_idx[list][pic_block_y];
      mv_array  = &motion->mv     [list][pic_block_y];

      for (ref=0; ref < p_Vid->listXsize[list+list_offset]; ref++)
      {
        m_cost = &p_Vid->motion_cost[blocktype][list][ref][block8x8];

        //===== LOOP OVER SUB MACRO BLOCK partitions
        updateMV_mp(currMB, m_cost, ref, list, bx, by, blocktype, block8x8);
        set_me_parameters(ref_array, mv_array, currSlice->all_mv[list][ref][blocktype][by][bx], ref, step_h, step_v, pic_block_x);
      }
    }
  }
  else
  {
    //int ref_pics_valid=-1;

    // Set flag for 8x8 Hadamard consideration for SATD (only used when 8x8 integer DCT is used for encoding)
    mv_block.test8x8 = p_Inp->Transform8x8Mode;

    init_mv_block(currMB, &mv_block, (short) blocktype, list, (char) ref, bx, by);

    if (p_Inp->SearchMode == EPZS)
    {
      if (p_Inp->EPZSSubPelGrid)
        currMB->IntPelME = EPZSIntPelBlockMotionSearch;
      else
        currMB->IntPelME = EPZSPelBlockMotionSearch;
    }

    get_original_block(p_Vid, &mv_block);
    //--- motion search for block ---   
    {
      //===== LOOP OVER REFERENCE FRAMES =====
      for (list = 0; list < numlists; list++)
      {
        //----- set arrays -----
        ref_array = &motion->ref_idx[list][pic_block_y];
        mv_array  = &motion->mv     [list][pic_block_y];
        mv_block.list = (char) list;
        for (ref=0; ref < p_Vid->listXsize[list+list_offset]; ref++) 
        {
            mv_block.ref_idx = (char) ref;
            m_cost = &p_Vid->motion_cost[blocktype][list][ref][block8x8];

            //----- set search range ---
            get_search_range(&mv_block, p_Inp, ref, blocktype);

            //===== LOOP OVER MACROBLOCK partitions        
            *m_cost = BlockMotionSearch (currMB, &mv_block, bx<<2, by<<2, lambda_factor);             
            //--- set motion vectors and reference frame ---
            set_me_parameters(ref_array, mv_array, currSlice->all_mv[list][ref][blocktype][by][bx], ref, step_h, step_v, pic_block_x);        
        }

      }
    }

    free_mv_block(p_Inp, &mv_block);
  }

#if GET_METIME
  gettime(&me_time_end);   // end time ms
  me_tmp_time = timediff (&me_time_start, &me_time_end);
  p_Vid->me_tot_time += me_tmp_time;
  p_Vid->me_time += me_tmp_time;
#endif
}

/*!
 ************************************************************************
 * \brief
 *    Motion search for a submacroblock partition
 ************************************************************************
 */
void SubPartitionMotionSearch (Macroblock *currMB,
                               int    blocktype,
                               int    block8x8,
                               int    *lambda_factor)
{
  Slice *currSlice = currMB->p_slice;
  VideoParameters *p_Vid = currSlice->p_Vid;
  InputParameters *p_Inp = currSlice->p_Inp;
  PicMotionParams *motion = &p_Vid->enc_picture->motion;

  char  **ref_array;
  short ***mv_array;
  short *all_mv;
  short ref = 0;
  int   v, h;
  int   pic_block_y;
  int   parttype  = 4;
  short step_h0   = (part_size[ parttype][0]);
  short step_v0   = (part_size[ parttype][1]);
  short step_h    = (part_size[blocktype][0]);
  short step_v    = (part_size[blocktype][1]);
  short list = LIST_0;
  int   numlists  = (currSlice->slice_type == B_SLICE) ? 2 : 1;
  int   list_offset = currMB->list_offset;
  distblk   mcost;
  distblk   *m_cost;
  short by = by0[parttype][block8x8];
  short bx = bx0[parttype][block8x8];
  MEBlock  mv_block;
  //int ref_lambda = (p_Inp->rdopt) ? lambda_factor[Q_PEL] :  lambda_factor[Q_PEL] >> 2;

#if GET_METIME
  TIME_T me_time_start;
  TIME_T me_time_end;
  int64 me_tmp_time;
  gettime( &me_time_start );    // start time ms
#endif

  if (p_Vid->Motion_Selected == 1)
  {
    //===== LOOP OVER REFERENCE FRAMES =====
    for (list=0; list<numlists;list++)
    {
      ref_array = motion->ref_idx[list];
      mv_array  = motion->mv[list];
      for (ref=0; ref < p_Vid->listXsize[list+list_offset]; ref++)
      {
        m_cost = &p_Vid->motion_cost[blocktype][list][ref][block8x8];

        //===== LOOP OVER SUB MACRO BLOCK partitions
        for (v=by; v<by + step_v0; v += step_v)
        {
          pic_block_y = currMB->block_y + v;
          for (h=bx; h<bx+step_h0; h+=step_h)
          {
            all_mv = currSlice->all_mv[list][ref][blocktype][v][h];

            updateMV_mp(currMB, m_cost, ref, list, h, v, blocktype, block8x8);

            //--- set motion vectors and reference frame (for motion vector prediction) ---
            set_me_parameters(&ref_array [pic_block_y], &mv_array [pic_block_y], all_mv, ref, step_h, step_v, currMB->block_x + h);
          } // h
        } // v
      }
    }
  }
  else
  {
    //int ref_pics_valid = -1;
    // Set if 8x8 transform will be used if SATD is used
    mv_block.test8x8 = p_Inp->Transform8x8Mode && blocktype == 4;

    if (p_Inp->SearchMode == EPZS)
    {
      if (p_Inp->EPZSSubPelGrid)
      {
        if (blocktype > 4)
          currMB->IntPelME = EPZSIntPelBlockMotionSearchSubMB;
        else
          currMB->IntPelME = EPZSIntPelBlockMotionSearch;
      }
      else
      {
        if (blocktype > 4)
          currMB->IntPelME = EPZSPelBlockMotionSearchSubMB;
        else
          currMB->IntPelME = EPZSPelBlockMotionSearch;
      }
    }

    init_mv_block(currMB, &mv_block, (short) blocktype, list, (char) ref, bx, by);
    if (blocktype == 4)
      get_original_block(p_Vid, &mv_block);


    //===== LOOP OVER REFERENCE FRAMES =====
    for (list=0; list<numlists;list++)
    {
      mv_block.list = (char) list;
      //----- set arrays -----
      ref_array = motion->ref_idx[list];
      mv_array  = motion->mv[list];

      for (ref=0; ref < p_Vid->listXsize[list+list_offset]; ref++)
      {
        mv_block.ref_idx = (char) ref;
        m_cost = &p_Vid->motion_cost[blocktype][list][ref][block8x8];
        //----- set search range ---
        get_search_range(&mv_block, p_Inp, ref, blocktype);

        //----- init motion cost -----
        *m_cost = 0;

        //===== LOOP OVER SUB MACRO BLOCK partitions
        for (v=by; v<by + step_v0; v += step_v)
        {
          pic_block_y = currMB->block_y + v;

          for (h=bx; h<bx+step_h0; h+=step_h)
          {
            all_mv = currSlice->all_mv[list][ref][blocktype][v][h];

            //--- motion search for block ---          
            {
              update_mv_block(currMB, &mv_block, h, v);
              //----- set search range ---
              get_search_range(&mv_block, p_Inp, ref, blocktype);

              mcost = BlockMotionSearch (currMB, &mv_block, h<<2, v<<2, lambda_factor);

              *m_cost += mcost;

            }

            //--- set motion vectors and reference frame (for motion vector prediction) ---
            set_me_parameters(&ref_array [pic_block_y], &mv_array [pic_block_y], all_mv, ref, step_h, step_v, currMB->block_x + h);
          }
        }

        if ( (p_Inp->Transform8x8Mode == 1) && p_Inp->RDOQ_CP_MV && (blocktype == 4) && currMB->luma_transform_size_8x8_flag)
        {
          currSlice->tmp_mv8[list][ref][by][bx].mv_x = currSlice->all_mv[list][ref][blocktype][by][bx][0];
          currSlice->tmp_mv8[list][ref][by][bx].mv_y = currSlice->all_mv[list][ref][blocktype][by][bx][1];
          currSlice->motion_cost8[list][ref][block8x8] = *m_cost;
        }
        else if ( (p_Inp->Transform8x8Mode == 1) && p_Inp->RDOQ_CP_MV && (blocktype == 4) && currMB->luma_transform_size_8x8_flag == 0)
        {
          currSlice->tmp_mv4[list][ref][by][bx].mv_x = currSlice->all_mv[list][ref][blocktype][by][bx][0];
          currSlice->tmp_mv4[list][ref][by][bx].mv_y = currSlice->all_mv[list][ref][blocktype][by][bx][1];
          currSlice->motion_cost4[list][ref][block8x8] = *m_cost;
        }
      }
    }

    free_mv_block(p_Inp, &mv_block);
  }

#if GET_METIME
  gettime(&me_time_end);   // end time ms
  me_tmp_time = timediff (&me_time_start, &me_time_end);
  p_Vid->me_tot_time += me_tmp_time;
  p_Vid->me_time += me_tmp_time;
#endif
}

/*!
 ************************************************************************
 * \brief
 *    Calculate Temporal Direct Mode Motion Vectors
 ************************************************************************
 */
void Get_Direct_MV_Temporal (Macroblock *currMB)
{
  Slice *currSlice = currMB->p_slice; 
  int   block_x, block_y, pic_block_x, pic_block_y, opic_block_x, opic_block_y;
  short ******all_mvs;
  int   mv_scale;
  int refList;
  int ref_idx;
  VideoParameters *p_Vid = currMB->p_Vid;
  int64 *refpic = p_Vid->enc_picture->ref_pic_num[LIST_0 +currMB->list_offset];  

  MotionParams *colocated;  

  if (currMB->list_offset)
  {
    if(currMB->mbAddrX & 0x01)
    {
      colocated = &currSlice->p_colocated->bottom;
    }
    else
    {
      colocated = &currSlice->p_colocated->top;
    }
  }
  else
  {
    colocated = &currSlice->p_colocated->frame;
  }


  //temporal direct mode copy from decoder
  for (block_y = 0; block_y < 4; block_y++)
  {
    pic_block_y  = currMB->block_y + block_y;
    opic_block_y = (currMB->opix_y >> 2) + block_y;

    for (block_x = 0; block_x < 4; block_x++)
    {
      pic_block_x  = currMB->block_x + block_x;
      opic_block_x = (currMB->pix_x>>2) + block_x;
      all_mvs = currSlice->all_mv;

      refList = (colocated->ref_idx[LIST_0][opic_block_y][opic_block_x]== -1 ? LIST_1 : LIST_0);
      ref_idx = colocated->ref_idx[refList][opic_block_y][opic_block_x];

      // next P is intra mode
      if (ref_idx==-1)
      {
        memset(all_mvs[LIST_0][0][0][block_y][block_x], 0, 2* sizeof(short));
        memset(all_mvs[LIST_1][0][0][block_y][block_x], 0, 2* sizeof(short));
        currSlice->direct_ref_idx[pic_block_y][pic_block_x][LIST_0] = 0;
        currSlice->direct_ref_idx[pic_block_y][pic_block_x][LIST_1] = 0;
        currSlice->direct_pdir[pic_block_y][pic_block_x] = 2;
      }
      // next P is skip or inter mode
      else
      {
        int mapped_idx=INVALIDINDEX;
        int iref;

        for (iref = 0; iref < imin(currSlice->num_ref_idx_active[LIST_0], p_Vid->listXsize[LIST_0 + currMB->list_offset]); iref++)
        {
          if (refpic[iref]==colocated->ref_pic_id[refList ][opic_block_y][opic_block_x])
          {
            mapped_idx=iref;
            break;
          }
          else //! invalid index. Default to zero even though this case should not happen
          {
            mapped_idx=INVALIDINDEX;
          }
        }

        if (mapped_idx !=INVALIDINDEX)
        {
          mv_scale = currSlice->mvscale[LIST_0+currMB->list_offset][mapped_idx];

          if (mv_scale==9999)
          {
            // forward
            memcpy(all_mvs[LIST_0][0][0][block_y][block_x], colocated->mv[refList][opic_block_y][opic_block_x], 2* sizeof(short));
            // backward
            memset(all_mvs[LIST_1][0][0][block_y][block_x], 0, 2* sizeof(short));
          }
          else
          {
            // forward
            all_mvs[LIST_0][mapped_idx][0][block_y][block_x][0] = (short) ((mv_scale * colocated->mv[refList][opic_block_y][opic_block_x][0] + 128) >> 8);
            all_mvs[LIST_0][mapped_idx][0][block_y][block_x][1] = (short) ((mv_scale * colocated->mv[refList][opic_block_y][opic_block_x][1] + 128) >> 8);
            // backward
            all_mvs[LIST_1][         0][0][block_y][block_x][0] = (short) (((mv_scale - 256)* colocated->mv[refList][opic_block_y][opic_block_x][0] + 128) >> 8);
            all_mvs[LIST_1][         0][0][block_y][block_x][1] = (short) (((mv_scale - 256)* colocated->mv[refList][opic_block_y][opic_block_x][1] + 128) >> 8);
          }

          // Test Level Limits if satisfied.
          if ( out_of_bounds_mvs(p_Vid, all_mvs[LIST_0][mapped_idx][0][block_y][block_x])|| out_of_bounds_mvs(p_Vid, all_mvs[LIST_1][0][0][block_y][block_x]))
          {
            currSlice->direct_ref_idx[pic_block_y][pic_block_x][LIST_0] = -1;
            currSlice->direct_ref_idx[pic_block_y][pic_block_x][LIST_1] = -1;
            currSlice->direct_pdir[pic_block_y][pic_block_x] = -1;
          }
          else
          {
            currSlice->direct_ref_idx[pic_block_y][pic_block_x][LIST_0] = (char) mapped_idx;
            currSlice->direct_ref_idx[pic_block_y][pic_block_x][LIST_1] = 0;
            currSlice->direct_pdir[pic_block_y][pic_block_x] = 2;
          }
        }
        else
        {
          currSlice->direct_ref_idx[pic_block_y][pic_block_x][LIST_0] = -1;
          currSlice->direct_ref_idx[pic_block_y][pic_block_x][LIST_1] = -1;
          currSlice->direct_pdir[pic_block_y][pic_block_x] = -1;
        }
      }
      if (p_Vid->active_pps->weighted_bipred_idc == 1 && currSlice->direct_pdir[pic_block_y][pic_block_x] == 2)
      {
        int weight_sum, i;
        short l0_refX = currSlice->direct_ref_idx[pic_block_y][pic_block_x][LIST_0];
        short l1_refX = currSlice->direct_ref_idx[pic_block_y][pic_block_x][LIST_1];
        for (i=0;i< (p_Vid->active_sps->chroma_format_idc == YUV400 ? 1 : 3); i++)
        {
          weight_sum = currSlice->wbp_weight[0][l0_refX][l1_refX][i] + currSlice->wbp_weight[1][l0_refX][l1_refX][i];
          if (weight_sum < -128 ||  weight_sum > 127)
          {
            currSlice->direct_ref_idx[pic_block_y][pic_block_x][LIST_0] = -1;
            currSlice->direct_ref_idx[pic_block_y][pic_block_x][LIST_1] = -1;
            currSlice->direct_pdir   [pic_block_y][pic_block_x]         = -1;
            break;
          }
        }
      }
    }
  }
}

/*!
************************************************************************
* \brief
*    Calculate Spatial Direct Mode Motion Vectors 
************************************************************************
*/
void Get_Direct_MV_Spatial_Normal (Macroblock *currMB)
{
  Slice *currSlice = currMB->p_slice; 
  VideoParameters *p_Vid = currMB->p_Vid;
  PicMotionParams *motion = &p_Vid->enc_picture->motion;
  short l0_refA, l0_refB, l0_refC;
  short l1_refA, l1_refB, l1_refC;
  short l0_refX,l1_refX;
  short pmvfw[2]={0,0},pmvbw[2]={0,0};

  int   block_x, block_y, pic_block_x, pic_block_y, opic_block_x, opic_block_y;
  short ******all_mvs;
  char  *direct_ref_idx;

  MotionParams *colocated;
  char  **     ref_pic_l0 = motion->ref_idx[LIST_0];
  char  **     ref_pic_l1 = motion->ref_idx[LIST_1];

  PixelPos mb[4];  
  get_neighbors(currMB, mb, 0, 0, 16);

  if (currMB->list_offset)
  {
    if(currMB->mbAddrX & 0x01)
    {
      colocated = &currSlice->p_colocated->bottom;
    }
    else
    {
      colocated = &currSlice->p_colocated->top;
    }
  }
  else
  {
    colocated = &currSlice->p_colocated->frame;
  }

  l0_refA = mb[0].available ? ref_pic_l0[mb[0].pos_y][mb[0].pos_x] : -1;
  l0_refB = mb[1].available ? ref_pic_l0[mb[1].pos_y][mb[1].pos_x] : -1;
  l0_refC = mb[2].available ? ref_pic_l0[mb[2].pos_y][mb[2].pos_x] : -1;

  l1_refA = mb[0].available ? ref_pic_l1[mb[0].pos_y][mb[0].pos_x] : -1;
  l1_refB = mb[1].available ? ref_pic_l1[mb[1].pos_y][mb[1].pos_x] : -1;
  l1_refC = mb[2].available ? ref_pic_l1[mb[2].pos_y][mb[2].pos_x] : -1;

  l0_refX = (short) ((l0_refA >= 0 && l0_refB >= 0) ? imin(l0_refA,l0_refB): imax(l0_refA,l0_refB));
  l0_refX = (short) ((l0_refX >= 0 && l0_refC >= 0) ? imin(l0_refX,l0_refC): imax(l0_refX,l0_refC));

  l1_refX = (short) ((l1_refA >= 0 && l1_refB >= 0) ? imin(l1_refA,l1_refB): imax(l1_refA,l1_refB));
  l1_refX = (short) ((l1_refX >= 0 && l1_refC >= 0) ? imin(l1_refX,l1_refC): imax(l1_refX,l1_refC));

  if (l0_refX >=0)
    currMB->GetMVPredictor (currMB, mb, pmvfw, l0_refX, motion->ref_idx[LIST_0], motion->mv[LIST_0], 0, 0, 16, 16);

  if (l1_refX >=0)
    currMB->GetMVPredictor (currMB, mb, pmvbw, l1_refX, motion->ref_idx[LIST_1], motion->mv[LIST_1], 0, 0, 16, 16);

  for (block_y=0; block_y<4; block_y++)
  {
    pic_block_y  = currMB->block_y + block_y;
    opic_block_y = (currMB->opix_y >> 2) + block_y;

    for (block_x=0; block_x<4; block_x++)
    {
      pic_block_x  = currMB->block_x + block_x;
      direct_ref_idx = currSlice->direct_ref_idx[pic_block_y][pic_block_x];
      opic_block_x = (currMB->pix_x >> 2) + block_x;

      all_mvs = currSlice->all_mv;

      if (l0_refX >=0)
      {
        if (!l0_refX  && !colocated->moving_block[opic_block_y][opic_block_x])
        {
          memset(all_mvs[LIST_0][0][0][block_y][block_x], 0, 2 * sizeof(short));
          direct_ref_idx[LIST_0]=0;
        }
        else
        {
          memcpy(all_mvs[LIST_0][l0_refX][0][block_y][block_x], pmvfw, 2 * sizeof(short));
          direct_ref_idx[LIST_0]= (char)l0_refX;
        }
      }
      else
      {
        memset(all_mvs[LIST_0][0][0][block_y][block_x], 0, 2 * sizeof(short));
        direct_ref_idx[LIST_0]=-1;
      }

      if (l1_refX >=0)
      {
        if(l1_refX==0 && !colocated->moving_block[opic_block_y][opic_block_x])
        {
          memset(all_mvs[LIST_1][0][0][block_y][block_x], 0, 2 * sizeof(short));
          direct_ref_idx[LIST_1]= (char)l1_refX;
        }
        else
        {
          memcpy(all_mvs[LIST_1][l1_refX][0][block_y][block_x], pmvbw, 2 * sizeof(short));
          direct_ref_idx[LIST_1] = (char)l1_refX;
        }
      }
      else
      {
        memset(all_mvs[LIST_1][0][0][block_y][block_x], 0, 2 * sizeof(short));
        direct_ref_idx[LIST_1] = -1;
      }

      // Test Level Limits if satisfied.
      if (l0_refX < 0 && l1_refX < 0)
      {
        direct_ref_idx[LIST_0] = direct_ref_idx[LIST_1] = 0;
        l0_refX = 0;
        l1_refX = 0;
      }

      if      (direct_ref_idx[LIST_1] == -1)
        currSlice->direct_pdir[pic_block_y][pic_block_x] = 0;
      else if (direct_ref_idx[LIST_0] == -1)
        currSlice->direct_pdir[pic_block_y][pic_block_x] = 1;
      else if (p_Vid->active_pps->weighted_bipred_idc == 1)
      {
        int weight_sum, i;
        Boolean invalid_wp = FALSE;
        for (i=0;i< (p_Vid->active_sps->chroma_format_idc == YUV400 ? 1 : 3); i++)
        {
          weight_sum = currSlice->wbp_weight[0][l0_refX][l1_refX][i] + currSlice->wbp_weight[1][l0_refX][l1_refX][i];
          if (weight_sum < -128 ||  weight_sum > 127)
          {
            invalid_wp = TRUE;
            break;
          }
        }
        if (invalid_wp == FALSE)
          currSlice->direct_pdir[pic_block_y][pic_block_x] = 2;
        else
        {
          direct_ref_idx[LIST_0] = -1;
          direct_ref_idx[LIST_1] = -1;
          currSlice->direct_pdir           [pic_block_y][pic_block_x] = -1;
        }
      }
      else
        currSlice->direct_pdir[pic_block_y][pic_block_x] = 2;
    }
  }
}


/*!
************************************************************************
* \brief
*    Calculate Spatial Direct Mode Motion Vectors 
************************************************************************
*/
void Get_Direct_MV_Spatial_MBAFF (Macroblock *currMB)
{
  short l0_refA, l0_refB, l0_refC;
  short l1_refA, l1_refB, l1_refC;
  short l0_refX,l1_refX;
  short pmvfw[2]={0,0},pmvbw[2]={0,0};

  int   block_x, block_y, pic_block_x, pic_block_y, opic_block_x, opic_block_y;
  short ******all_mvs;
  char  *direct_ref_idx;

  MotionParams *colocated;
  Slice *currSlice = currMB->p_slice;
  VideoParameters *p_Vid = currMB->p_Vid;
  PicMotionParams *motion = &p_Vid->enc_picture->motion;

  char  **ref_pic_l0 = motion->ref_idx[LIST_0];
  char  **ref_pic_l1 = motion->ref_idx[LIST_1];

  PixelPos mb[4];  
  get_neighbors(currMB, mb, 0, 0, 16);

  if (currMB->list_offset)
  {
    if(currMB->mbAddrX & 0x01)
    {
      colocated = &currSlice->p_colocated->bottom;
    }
    else
    {
      colocated = &currSlice->p_colocated->top;
    }
  }
  else
  {
    colocated = &currSlice->p_colocated->frame;
  }

  if (currMB->mb_field)
  {
    l0_refA = mb[0].available
      ? (p_Vid->mb_data[mb[0].mb_addr].mb_field  || ref_pic_l0[mb[0].pos_y][mb[0].pos_x] < 0
      ?  ref_pic_l0[mb[0].pos_y][mb[0].pos_x]
    :  ref_pic_l0[mb[0].pos_y][mb[0].pos_x] * 2) : -1;

    l0_refB = mb[1].available
      ? (p_Vid->mb_data[mb[1].mb_addr].mb_field || ref_pic_l0[mb[1].pos_y][mb[1].pos_x] < 0
      ?  ref_pic_l0[mb[1].pos_y][mb[1].pos_x]
    :  ref_pic_l0[mb[1].pos_y][mb[1].pos_x] * 2) : -1;

    l0_refC = mb[2].available
      ? (p_Vid->mb_data[mb[2].mb_addr].mb_field || ref_pic_l0[mb[2].pos_y][mb[2].pos_x] < 0
      ?  ref_pic_l0[mb[2].pos_y][mb[2].pos_x]
    :  ref_pic_l0[mb[2].pos_y][mb[2].pos_x] * 2) : -1;

    l1_refA = mb[0].available
      ? (p_Vid->mb_data[mb[0].mb_addr].mb_field || ref_pic_l1[mb[0].pos_y][mb[0].pos_x] < 0
      ?  ref_pic_l1[mb[0].pos_y][mb[0].pos_x]
    :  ref_pic_l1[mb[0].pos_y][mb[0].pos_x] * 2) : -1;

    l1_refB = mb[1].available
      ? (p_Vid->mb_data[mb[1].mb_addr].mb_field || ref_pic_l1[mb[1].pos_y][mb[1].pos_x] < 0
      ?  ref_pic_l1[mb[1].pos_y][mb[1].pos_x]
    :  ref_pic_l1[mb[1].pos_y][mb[1].pos_x] * 2) : -1;

    l1_refC = mb[2].available
      ? (p_Vid->mb_data[mb[2].mb_addr].mb_field || ref_pic_l1[mb[2].pos_y][mb[2].pos_x] < 0
      ?  ref_pic_l1[mb[2].pos_y][mb[2].pos_x]
    :  ref_pic_l1[mb[2].pos_y][mb[2].pos_x] * 2) : -1;
  }
  else
  {
    l0_refA = mb[0].available
      ? (p_Vid->mb_data[mb[0].mb_addr].mb_field || ref_pic_l0[mb[0].pos_y][mb[0].pos_x]  < 0
      ?  ref_pic_l0[mb[0].pos_y][mb[0].pos_x] >> 1
      :  ref_pic_l0[mb[0].pos_y][mb[0].pos_x]) : -1;

    l0_refB = mb[1].available
      ? (p_Vid->mb_data[mb[1].mb_addr].mb_field || ref_pic_l0[mb[1].pos_y][mb[1].pos_x] < 0
      ?  ref_pic_l0[mb[1].pos_y][mb[1].pos_x] >> 1
      :  ref_pic_l0[mb[1].pos_y][mb[1].pos_x]) : -1;

    l0_refC = mb[2].available
      ? (p_Vid->mb_data[mb[2].mb_addr].mb_field || ref_pic_l0[mb[2].pos_y][mb[2].pos_x] < 0
      ?  ref_pic_l0[mb[2].pos_y][mb[2].pos_x] >> 1
      :  ref_pic_l0[mb[2].pos_y][mb[2].pos_x]) : -1;

    l1_refA = mb[0].available
      ? (p_Vid->mb_data[mb[0].mb_addr].mb_field || ref_pic_l1[mb[0].pos_y][mb[0].pos_x] < 0
      ?  ref_pic_l1[mb[0].pos_y][mb[0].pos_x] >> 1
      :  ref_pic_l1[mb[0].pos_y][mb[0].pos_x]) : -1;

    l1_refB = mb[1].available
      ? (p_Vid->mb_data[mb[1].mb_addr].mb_field || ref_pic_l1[mb[1].pos_y][mb[1].pos_x] < 0
      ?  ref_pic_l1[mb[1].pos_y][mb[1].pos_x] >> 1
      :  ref_pic_l1[mb[1].pos_y][mb[1].pos_x]) : -1;

    l1_refC = mb[2].available
      ? (p_Vid->mb_data[mb[2].mb_addr].mb_field || ref_pic_l1[mb[2].pos_y][mb[2].pos_x] < 0
      ?  ref_pic_l1[mb[2].pos_y][mb[2].pos_x] >> 1
      :  ref_pic_l1[mb[2].pos_y][mb[2].pos_x]) : -1;
  }

  l0_refX = (short) ((l0_refA >= 0 && l0_refB >= 0) ? imin(l0_refA,l0_refB): imax(l0_refA,l0_refB));
  l0_refX = (short) ((l0_refX >= 0 && l0_refC >= 0) ? imin(l0_refX,l0_refC): imax(l0_refX,l0_refC));

  l1_refX = (short) ((l1_refA >= 0 && l1_refB >= 0) ? imin(l1_refA,l1_refB): imax(l1_refA,l1_refB));
  l1_refX = (short) ((l1_refX >= 0 && l1_refC >= 0) ? imin(l1_refX,l1_refC): imax(l1_refX,l1_refC));

  if (l0_refX >=0)
    currMB->GetMVPredictor (currMB, mb, pmvfw, l0_refX, motion->ref_idx[LIST_0], motion->mv[LIST_0], 0, 0, 16, 16);

  if (l1_refX >=0)
    currMB->GetMVPredictor (currMB, mb, pmvbw, l1_refX, motion->ref_idx[LIST_1], motion->mv[LIST_1], 0, 0, 16, 16);

  for (block_y=0; block_y<4; block_y++)
  {
    pic_block_y  = currMB->block_y + block_y;
    opic_block_y = (currMB->opix_y >> 2) + block_y;

    for (block_x=0; block_x<4; block_x++)
    {
      pic_block_x  = currMB->block_x + block_x;
      direct_ref_idx = currSlice->direct_ref_idx[pic_block_y][pic_block_x];
      opic_block_x = (currMB->pix_x >> 2) + block_x;

      all_mvs = currSlice->all_mv;

      if (l0_refX >=0)
      {
        if (!l0_refX  && !colocated->moving_block[opic_block_y][opic_block_x])
        {
          memset(all_mvs[LIST_0][0][0][block_y][block_x], 0, 2 * sizeof(short));
          direct_ref_idx[LIST_0] = 0;
        }
        else
        {
          memcpy(all_mvs[LIST_0][l0_refX][0][block_y][block_x], pmvfw, 2 * sizeof(short));
          direct_ref_idx[LIST_0] = (char)l0_refX;
        }
      }
      else
      {
        memset(all_mvs[LIST_0][0][0][block_y][block_x], 0, 2 * sizeof(short));
        direct_ref_idx[LIST_0] = -1;
      }

      if (l1_refX >=0)
      {
        if(l1_refX==0 && !colocated->moving_block[opic_block_y][opic_block_x])
        {
          memset(all_mvs[LIST_1][0][0][block_y][block_x], 0, 2 * sizeof(short));
          direct_ref_idx[LIST_1] = (char)l1_refX;
        }
        else
        {
          memcpy(all_mvs[LIST_1][l1_refX][0][block_y][block_x], pmvbw, 2 * sizeof(short));
          direct_ref_idx[LIST_1] = (char)l1_refX;
        }
      }
      else
      {
        memset(all_mvs[LIST_1][0][0][block_y][block_x], 0, 2 * sizeof(short));
        direct_ref_idx[LIST_1] = -1;
      }

     // Test Level Limits if satisfied.

      // Test Level Limits if satisfied.
      if ((out_of_bounds_mvs(p_Vid, all_mvs[LIST_0][l0_refX < 0? 0 : l0_refX][0][block_y][block_x])
        ||  out_of_bounds_mvs(p_Vid, all_mvs[LIST_1][l1_refX < 0? 0 : l1_refX][0][block_y][block_x])))
      {
        direct_ref_idx[LIST_0] = -1;
        direct_ref_idx[LIST_1] = -1;
        currSlice->direct_pdir   [pic_block_y][pic_block_x]         = -1;
      }     
      else
      {
        if (l0_refX < 0 && l1_refX < 0)
        {
          direct_ref_idx[LIST_0] = direct_ref_idx[LIST_1] = 0;
          l0_refX = 0;
          l1_refX = 0;
        }

          if      (direct_ref_idx[LIST_1] == -1)
            currSlice->direct_pdir[pic_block_y][pic_block_x] = 0;
          else if (direct_ref_idx[LIST_0] == -1)
            currSlice->direct_pdir[pic_block_y][pic_block_x] = 1;
          else if (p_Vid->active_pps->weighted_bipred_idc == 1)
          {
            int weight_sum, i;
            Boolean invalid_wp = FALSE;
            for (i=0;i< (p_Vid->active_sps->chroma_format_idc == YUV400 ? 1 : 3); i++)
            {
              weight_sum = currSlice->wbp_weight[0][l0_refX][l1_refX][i] + currSlice->wbp_weight[1][l0_refX][l1_refX][i];
              if (weight_sum < -128 ||  weight_sum > 127)
              {
                invalid_wp = TRUE;
                break;
              }
            }
            if (invalid_wp == FALSE)
              currSlice->direct_pdir[pic_block_y][pic_block_x] = 2;
            else
            {
              direct_ref_idx[LIST_0] = -1;
              direct_ref_idx[LIST_1] = -1;
              currSlice->direct_pdir           [pic_block_y][pic_block_x] = -1;
            }
          }
          else
            currSlice->direct_pdir[pic_block_y][pic_block_x] = 2;
        }
    }
  }
}
