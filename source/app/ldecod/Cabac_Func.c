#include "Cabac_Func.h"

static const int assignSE2partition[][SE_MAX_ELEMENTS] =
{
	// 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19  // element number (do not uncomment)
	{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },   //!< all elements in one partition no data partitioning
	{ 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 2, 2, 2, 2, 0, 0, 0, 0 }    //!< three partitions per slice
};

static const short maxpos[] = { 15, 14, 63, 31, 31, 15, 3, 14, 7, 15, 15, 14, 63, 31, 31, 15, 15, 14, 63, 31, 31, 15 };
static const short c1isdc[] = { 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1 };
static const short type2ctx_bcbp[] = { 0, 1, 2, 3, 3, 4, 5, 6, 5, 5, 10, 11, 12, 13, 13, 14, 16, 17, 18, 19, 19, 20 };
static const short type2ctx_map[] = { 0, 1, 2, 3, 4, 5, 6, 7, 6, 6, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21 }; // 8
static const short type2ctx_last[] = { 0, 1, 2, 3, 4, 5, 6, 7, 6, 6, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21 }; // 8
static const short type2ctx_one[] = { 0, 1, 2, 3, 3, 4, 5, 6, 5, 5, 10, 11, 12, 13, 13, 14, 16, 17, 18, 19, 19, 20 }; // 7
static const short type2ctx_abs[] = { 0, 1, 2, 3, 3, 4, 5, 6, 5, 5, 10, 11, 12, 13, 13, 14, 16, 17, 18, 19, 19, 20 }; // 7
static const short max_c2[] = { 4, 4, 4, 4, 4, 4, 3, 4, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4 }; // 9
static const byte renorm_table_32[32] = { 6, 5, 4, 4, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
static const byte  pos2ctx_map8x8[] = { 0, 1, 2, 3, 4, 5, 5, 4, 4, 3, 3, 4, 4, 4, 5, 5,
4, 4, 4, 4, 3, 3, 6, 7, 7, 7, 8, 9, 10, 9, 8, 7,
7, 6, 11, 12, 13, 11, 6, 7, 8, 9, 14, 10, 9, 8, 6, 11,
12, 13, 11, 6, 9, 14, 10, 9, 11, 12, 13, 11, 14, 10, 12, 14 }; // 15 CTX
static const byte  pos2ctx_map8x4[] = { 0, 1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 9, 8, 6, 7, 8,
9, 10, 11, 9, 8, 6, 12, 8, 9, 10, 11, 9, 13, 13, 14, 14 }; // 15 CTX
static const byte  pos2ctx_map4x4[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 14 }; // 15 CTX
static const byte  pos2ctx_map2x4c[] = { 0, 0, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2 }; // 15 CTX
static const byte  pos2ctx_map4x4c[] = { 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2 }; // 15 CTX
static const byte* pos2ctx_map[] = { pos2ctx_map4x4, pos2ctx_map4x4, pos2ctx_map8x8, pos2ctx_map8x4,
pos2ctx_map8x4, pos2ctx_map4x4, pos2ctx_map4x4, pos2ctx_map4x4,
pos2ctx_map2x4c, pos2ctx_map4x4c,
pos2ctx_map4x4, pos2ctx_map4x4, pos2ctx_map8x8, pos2ctx_map8x4,
pos2ctx_map8x4, pos2ctx_map4x4,
pos2ctx_map4x4, pos2ctx_map4x4, pos2ctx_map8x8, pos2ctx_map8x4,
pos2ctx_map8x4, pos2ctx_map4x4 };

static const byte  pos2ctx_last8x8[] = { 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4,
5, 5, 5, 5, 6, 6, 6, 6, 7, 7, 7, 7, 8, 8, 8, 8 }; //  9 CTX
static const byte  pos2ctx_last8x4[] = { 0, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2,
3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8 }; //  9 CTX

static const byte  pos2ctx_last4x4[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 }; // 15 CTX
static const byte  pos2ctx_last2x4c[] = { 0, 0, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2 }; // 15 CTX
static const byte  pos2ctx_last4x4c[] = { 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2 }; // 15 CTX
static const byte* pos2ctx_last[] = { pos2ctx_last4x4, pos2ctx_last4x4, pos2ctx_last8x8, pos2ctx_last8x4,
pos2ctx_last8x4, pos2ctx_last4x4, pos2ctx_last4x4, pos2ctx_last4x4,
pos2ctx_last2x4c, pos2ctx_last4x4c,
pos2ctx_last4x4, pos2ctx_last4x4, pos2ctx_last8x8, pos2ctx_last8x4,
pos2ctx_last8x4, pos2ctx_last4x4,
pos2ctx_last4x4, pos2ctx_last4x4, pos2ctx_last8x8, pos2ctx_last8x4,
pos2ctx_last8x4, pos2ctx_last4x4 };

static const byte  pos2ctx_map8x8i[] = { 0, 1, 1, 2, 2, 3, 3, 4, 5, 6, 7, 7, 7, 8, 4, 5,
6, 9, 10, 10, 8, 11, 12, 11, 9, 9, 10, 10, 8, 11, 12, 11,
9, 9, 10, 10, 8, 11, 12, 11, 9, 9, 10, 10, 8, 13, 13, 9,
9, 10, 10, 8, 13, 13, 9, 9, 10, 10, 14, 14, 14, 14, 14, 14 }; // 15 CTX
static const byte  pos2ctx_map8x4i[] = { 0, 1, 2, 3, 4, 5, 6, 3, 4, 5, 6, 3, 4, 7, 6, 8,
9, 7, 6, 8, 9, 10, 11, 12, 12, 10, 11, 13, 13, 14, 14, 14 }; // 15 CTX
static const byte  pos2ctx_map4x8i[] = { 0, 1, 1, 1, 2, 3, 3, 4, 4, 4, 5, 6, 2, 7, 7, 8,
8, 8, 5, 6, 9, 10, 10, 11, 11, 11, 12, 13, 13, 14, 14, 14 }; // 15 CTX
static const byte* pos2ctx_map_int[] = { pos2ctx_map4x4, pos2ctx_map4x4, pos2ctx_map8x8i, pos2ctx_map8x4i,
pos2ctx_map4x8i, pos2ctx_map4x4, pos2ctx_map4x4, pos2ctx_map4x4,
pos2ctx_map2x4c, pos2ctx_map4x4c,
pos2ctx_map4x4, pos2ctx_map4x4, pos2ctx_map8x8i, pos2ctx_map8x4i,
pos2ctx_map8x4i, pos2ctx_map4x4,
pos2ctx_map4x4, pos2ctx_map4x4, pos2ctx_map8x8i, pos2ctx_map8x4i,
pos2ctx_map8x4i, pos2ctx_map4x4 };

void reset_coding_state_cabac(Macroblock *currMB, Dec_CSobj *cs)
{
	Slice *currSlice = currMB->p_Slice;
	DataPartition *partArr = &currSlice->partArr[0];
	partArr->bitstream->bitstream_length = cs->bitstream[0].bitstream_length;
	partArr->bitstream->code_len = cs->bitstream[0].code_len;
	partArr->bitstream->ei_flag = cs->bitstream[0].ei_flag;
	partArr->bitstream->frame_bitoffset = cs->bitstream[0].frame_bitoffset;
	partArr->bitstream->read_len = cs->bitstream[0].read_len;
	partArr->de_cabac.DbitsLeft = cs->decenv[0].DbitsLeft;
	partArr->de_cabac.Drange = cs->decenv[0].Drange;
	partArr->de_cabac.Dvalue = cs->decenv[0].Dvalue;
	*partArr->de_cabac.Dcodestrm_len = *cs->decenv[0].Dcodestrm_len;

	*currSlice->mot_ctx = *cs->mot_ctx;
	*currSlice->tex_ctx = *cs->tex_ctx;
}

void store_coding_state_cabac(Macroblock *currMB, Dec_CSobj *cs)
{
	Slice *currSlice = currMB->p_Slice;
	DataPartition *partArr = &currSlice->partArr[0];

	cs->bitstream[0].bitstream_length = partArr->bitstream->bitstream_length;
	cs->bitstream[0].code_len = partArr->bitstream->code_len;
	cs->bitstream[0].ei_flag = partArr->bitstream->ei_flag;
	cs->bitstream[0].frame_bitoffset = partArr->bitstream->frame_bitoffset;
	cs->bitstream[0].read_len = partArr->bitstream->read_len;
	cs->decenv[0].DbitsLeft = partArr->de_cabac.DbitsLeft;
	cs->decenv[0].Drange = partArr->de_cabac.Drange;
	cs->decenv[0].Dvalue = partArr->de_cabac.Dvalue;
	cs->decenv[0].Dcodestrm_len = partArr->de_cabac.Dcodestrm_len;

	//=== contexts for binary arithmetic coding ===
	*cs->mot_ctx = *currSlice->mot_ctx;
	*cs->tex_ctx = *currSlice->tex_ctx;
}

void reset_coding_state_cabacE(EncodingEnvironment *ee_cabac, EncodingEnvironment *cs, int *len)
{
	ee_cabac->C = cs->C;
	ee_cabac->E = cs->E;
	ee_cabac->Ebits_to_go = cs->Ebits_to_go;
	ee_cabac->Ebuffer = cs->Ebuffer;
	ee_cabac->Echunks_outstanding = cs->Echunks_outstanding;
	ee_cabac->Elow = cs->Elow;
	ee_cabac->Epbuf = cs->Epbuf;
	ee_cabac->Erange = cs->Erange;
	*ee_cabac->Ecodestrm_len = *len;
}
void store_coding_state_cabacE(EncodingEnvironment *ee_cabac, EncodingEnvironment *cs, int *len)
{
	cs->C = ee_cabac->C;
	cs->E = ee_cabac->E;
	cs->Ebits_to_go = ee_cabac->Ebits_to_go;
	cs->Ebuffer = ee_cabac->Ebuffer;
	cs->Echunks_outstanding = ee_cabac->Echunks_outstanding;
	cs->Elow = ee_cabac->Elow;
	cs->Epbuf = ee_cabac->Epbuf;
	cs->Erange = ee_cabac->Erange;
	*len = *ee_cabac->Ecodestrm_len;
}
/*!
************************************************************************
* \brief
*    get neighbouring positions for non-aff coding
* \param currMB
*   current macroblock
* \param xN
*    input x position
* \param yN
*    input y position
* \param mb_size
*    Macroblock size in pixel (according to luma or chroma MB access)
* \param pix
*    returns position informations
************************************************************************
*/
void getNonAffNeighbourI(Enc_Macroblock *currMB, int xN, int yN, int mb_size[2], PixelPos *pix)
{
	BlockPos *PicPos = currMB->p_Vid->PicPos;
	if (xN < 0)
	{
		if (yN < 0)
		{
			pix->mb_addr = currMB->mbAddrD;
			pix->available = currMB->mbAvailD;
		}
		else if ((yN >= 0) && (yN < mb_size[1]))
		{
			pix->mb_addr = currMB->mbAddrA;
			pix->available = currMB->mbAvailA;
		}
		else
		{
			pix->available = FALSE;
		}
	}
	else if ((xN >= 0) && (xN < mb_size[0]))
	{
		if (yN<0)
		{
			pix->mb_addr = currMB->mbAddrB;
			pix->available = currMB->mbAvailB;
		}
		else if (((yN >= 0) && (yN < mb_size[1])))
		{
			pix->mb_addr = currMB->mbAddrX;
			pix->available = TRUE;
		}
		else
		{
			pix->available = FALSE;
		}
	}
	else if ((xN >= mb_size[0]) && (yN < 0))
	{
		pix->mb_addr = currMB->mbAddrC;
		pix->available = currMB->mbAvailC;
	}
	else
	{
		pix->available = FALSE;
	}

	if (pix->available || currMB->DeblockCall)
	{
		pix->x = (short)(xN & (mb_size[0] - 1));
		pix->y = (short)(yN & (mb_size[1] - 1));
		pix->pos_x = (short)(pix->x + PicPos[pix->mb_addr].x * mb_size[0]);
		pix->pos_y = (short)(pix->y + PicPos[pix->mb_addr].y * mb_size[1]);
	}
}

/*!
************************************************************************
* \brief
*    get neighboring 4x4 block
* \param currMB
*   current macroblock
* \param block_x
*    input x block position
* \param block_y
*    input y block position
* \param mb_size
*    Macroblock size in pixel (according to luma or chroma MB access)
* \param pix
*    returns position informations
************************************************************************
*/
void get4x4NeighbourI(Enc_Macroblock *currMB, int block_x, int block_y, int mb_size[2], PixelPos *pix)
{
	//currMB->p_Vid->getNeighbour(currMB, block_x, block_y, mb_size, pix);
	getNonAffNeighbourI(currMB, block_x, block_y, mb_size, pix);
	//getAffNeighbourI(currMB, block_x, block_y, mb_size, pix);
	if (pix->available)
	{
		pix->x >>= 2;
		pix->y >>= 2;
		pix->pos_x >>= 2;
		pix->pos_y >>= 2;
	}
}

static forceinline void put_buffer(EncodingEnvironmentPtr eep)
{
	while (eep->Epbuf >= 0)
	{
		eep->Ecodestrm[(*eep->Ecodestrm_len)++] = (byte)((eep->Ebuffer >> ((eep->Epbuf--) << 3)) & 0xFF);
		if (eep->Ecodestrm[(*eep->Ecodestrm_len) - 1] <= 3)
			if (eep->Ecodestrm[(*eep->Ecodestrm_len) - 2] == 0)
				if (eep->Ecodestrm[(*eep->Ecodestrm_len) - 3] == 0)
				{
					eep->Ecodestrm[(*eep->Ecodestrm_len)] = eep->Ecodestrm[(*eep->Ecodestrm_len) - 1];
					eep->Ecodestrm[(*eep->Ecodestrm_len) - 1] = 0x03;
					(*eep->Ecodestrm_len)++;
				}
	}
	while (eep->C > 7)
	{
		eep->C -= 8;
		++(eep->E);
	}
	eep->Ebuffer = 0;
}

static inline void put_one_word(EncodingEnvironmentPtr eep, int b)
{
	if (eep->Epbuf >= 3)
	{
		put_buffer(eep);
	}
	eep->Ebuffer <<= 16;
	eep->Ebuffer += (b);

	eep->Epbuf += 2;
}

static forceinline void propagate_carry(EncodingEnvironmentPtr eep)
{
	++(eep->Ebuffer);
	while (eep->Echunks_outstanding > 0)
	{
		put_one_word(eep, 0);
		--(eep->Echunks_outstanding);
	}
}
static inline void put_last_chunk_plus_outstanding(EncodingEnvironmentPtr eep, unsigned int l)
{
	while (eep->Echunks_outstanding > 0)
	{
		put_one_word(eep, 0xFFFF);
		--(eep->Echunks_outstanding);
	}
	put_one_word(eep, l);
}

/*!
************************************************************************
* \brief
*    Actually arithmetic encoding of one binary symbol by using
*    the probability estimate of its associated context model
************************************************************************
*/
void biari_encode_symbol(EncodingEnvironmentPtr eep, int symbol, BiContextTypePtr bi_ct)
{
	static const byte rLPS_table_64x4[64][4] =
	{
		{ 128, 176, 208, 240 },
		{ 128, 167, 197, 227 },
		{ 128, 158, 187, 216 },
		{ 123, 150, 178, 205 },
		{ 116, 142, 169, 195 },
		{ 111, 135, 160, 185 },
		{ 105, 128, 152, 175 },
		{ 100, 122, 144, 166 },
		{ 95, 116, 137, 158 },
		{ 90, 110, 130, 150 },
		{ 85, 104, 123, 142 },
		{ 81, 99, 117, 135 },
		{ 77, 94, 111, 128 },
		{ 73, 89, 105, 122 },
		{ 69, 85, 100, 116 },
		{ 66, 80, 95, 110 },
		{ 62, 76, 90, 104 },
		{ 59, 72, 86, 99 },
		{ 56, 69, 81, 94 },
		{ 53, 65, 77, 89 },
		{ 51, 62, 73, 85 },
		{ 48, 59, 69, 80 },
		{ 46, 56, 66, 76 },
		{ 43, 53, 63, 72 },
		{ 41, 50, 59, 69 },
		{ 39, 48, 56, 65 },
		{ 37, 45, 54, 62 },
		{ 35, 43, 51, 59 },
		{ 33, 41, 48, 56 },
		{ 32, 39, 46, 53 },
		{ 30, 37, 43, 50 },
		{ 29, 35, 41, 48 },
		{ 27, 33, 39, 45 },
		{ 26, 31, 37, 43 },
		{ 24, 30, 35, 41 },
		{ 23, 28, 33, 39 },
		{ 22, 27, 32, 37 },
		{ 21, 26, 30, 35 },
		{ 20, 24, 29, 33 },
		{ 19, 23, 27, 31 },
		{ 18, 22, 26, 30 },
		{ 17, 21, 25, 28 },
		{ 16, 20, 23, 27 },
		{ 15, 19, 22, 25 },
		{ 14, 18, 21, 24 },
		{ 14, 17, 20, 23 },
		{ 13, 16, 19, 22 },
		{ 12, 15, 18, 21 },
		{ 12, 14, 17, 20 },
		{ 11, 14, 16, 19 },
		{ 11, 13, 15, 18 },
		{ 10, 12, 15, 17 },
		{ 10, 12, 14, 16 },
		{ 9, 11, 13, 15 },
		{ 9, 11, 12, 14 },
		{ 8, 10, 12, 14 },
		{ 8, 9, 11, 13 },
		{ 7, 9, 11, 12 },
		{ 7, 9, 10, 12 },
		{ 7, 8, 10, 11 },
		{ 6, 8, 9, 11 },
		{ 6, 7, 9, 10 },
		{ 6, 7, 8, 9 },
		{ 2, 2, 2, 2 }
	};
	static const byte AC_next_state_MPS_64[64] =
	{
		1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
		11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
		21, 22, 23, 24, 25, 26, 27, 28, 29, 30,
		31, 32, 33, 34, 35, 36, 37, 38, 39, 40,
		41, 42, 43, 44, 45, 46, 47, 48, 49, 50,
		51, 52, 53, 54, 55, 56, 57, 58, 59, 60,
		61, 62, 62, 63
	};

	static const byte AC_next_state_LPS_64[64] =
	{
		0, 0, 1, 2, 2, 4, 4, 5, 6, 7,
		8, 9, 9, 11, 11, 12, 13, 13, 15, 15,
		16, 16, 18, 18, 19, 19, 21, 21, 22, 22,
		23, 24, 24, 25, 26, 26, 27, 27, 28, 29,
		29, 30, 30, 30, 31, 32, 32, 33, 33, 33,
		34, 34, 35, 35, 35, 36, 36, 36, 37, 37,
		37, 38, 38, 63
	};

	unsigned int low = eep->Elow;
	unsigned int range = eep->Erange;
	int bl = eep->Ebits_to_go;
	unsigned int rLPS = rLPS_table_64x4[bi_ct->state][(range >> 6) & 3];

	range -= rLPS;

	++(eep->C);
	//bi_ct->count += eep->p_Vid->cabac_encoding;

	/* covers all cases where code does not bother to shift down symbol to be
	* either 0 or 1, e.g. in some cases for cbp, mb_Type etc the code simply
	* masks off the bit position and passes in the resulting value */
	//symbol = (short) (symbol != 0);

	if ((symbol != 0) == bi_ct->MPS)  //MPS
	{
		bi_ct->state = AC_next_state_MPS_64[bi_ct->state]; // next state

		if (range >= 0x100) // no renorm
		{
			eep->Erange = range;
			return;
		}
		else
		{
			range <<= 1;
			if (--bl > 0)  // renorm once, no output
			{
				eep->Erange = range;
				eep->Ebits_to_go = bl;
				return;
			}
		}
	}
	else         //LPS
	{
		unsigned int renorm = renorm_table_32[(rLPS >> 3) & 0x1F];

		low += range << bl;
		range = (rLPS << renorm);
		bl -= renorm;

		if (!bi_ct->state)
			bi_ct->MPS ^= 0x01;               // switch MPS if necessary

		bi_ct->state = AC_next_state_LPS_64[bi_ct->state]; // next state

		if (low >= 0x04000000) // output of carry needed
		{
			low -= 0x04000000;
			propagate_carry(eep);
		}

		if (bl > 0)
		{
			eep->Elow = low;
			eep->Erange = range;
			eep->Ebits_to_go = bl;
			return;
		}
	}

	//renorm needed
	eep->Elow = (low << 16)& (0x03FFFFFF);
	low = (low >> 10) & 0xFFFF; // mask out the 8/16 MSBs for output

	if (low < 0xFFFF) // no carry possible, output now
	{
		put_last_chunk_plus_outstanding(eep, low);
	}
	else          // low == "FF.."; keep it, may affect future carry
	{
		++(eep->Echunks_outstanding);
	}
	eep->Erange = range;
	eep->Ebits_to_go = bl + 16;
}

/*!
****************************************************************************
* \brief
*    Write CBP4-BIT
****************************************************************************
*/
void write_and_store_CBP_block_bit(Enc_Macroblock* currMB, EncodingEnvironmentPtr eep_dp, int type, int cbp_bit, TextureInfoContexts*  tex_ctx)
{
	Enc_VideoParameters *p_Vid = currMB->p_Vid;
	//Enc_Macroblock *mb_data = p_Vid->mb_data;
	BiContextTypePtr bct;
	int y_ac = (type == LUMA_16AC || type == LUMA_8x8 || type == LUMA_8x4 || type == LUMA_4x8 || type == LUMA_4x4);
	int y_dc = (type == LUMA_16DC);
	int u_ac = (type == CHROMA_AC && !p_Vid->is_v_block);
	int v_ac = (type == CHROMA_AC &&  p_Vid->is_v_block);
	int chroma_dc = (type == CHROMA_DC || type == CHROMA_DC_2x4 || type == CHROMA_DC_4x4);
	int u_dc = (chroma_dc && !p_Vid->is_v_block);
	int v_dc = (chroma_dc &&  p_Vid->is_v_block);
	int j = (y_ac || u_ac || v_ac ? currMB->subblock_y : 0);
	int i = (y_ac || u_ac || v_ac ? currMB->subblock_x : 0);
	int bit = (y_dc ? 0 : y_ac ? 1 : u_dc ? 17 : v_dc ? 18 : u_ac ? 19 : 23);
	int default_bit = (currMB->is_intra_block ? 1 : 0);
	int upper_bit = default_bit;
	int left_bit = default_bit;
	int ctx;

	int bit_pos_a = 0;
	int bit_pos_b = 0;

	PixelPos block_a, block_b;

	if (y_ac || y_dc)
	{
		get4x4NeighbourI(currMB, i - 1, j, p_Vid->mb_size[IS_LUMA], &block_a);
		get4x4NeighbourI(currMB, i, j - 1, p_Vid->mb_size[IS_LUMA], &block_b);
		if (y_ac)
		{
			if (block_a.available)
				bit_pos_a = 4 * block_a.y + block_a.x;
			if (block_b.available)
				bit_pos_b = 4 * block_b.y + block_b.x;
		}
	}
	else
	{
		get4x4NeighbourI(currMB, i - 1, j, p_Vid->mb_size[IS_CHROMA], &block_a);
		get4x4NeighbourI(currMB, i, j - 1, p_Vid->mb_size[IS_CHROMA], &block_b);
		if (u_ac || v_ac)
		{
			if (block_a.available)
				bit_pos_a = (block_a.y << 2) + block_a.x;
			if (block_b.available)
				bit_pos_b = (block_b.y << 2) + block_b.x;
		}
	}

	bit = (y_dc ? 0 : y_ac ? 1 + j + (i >> 2) : u_dc ? 17 : v_dc ? 18 : u_ac ? 19 + j + (i >> 2) : 35 + j + (i >> 2));
	//--- set bits for current block ---
	if (cbp_bit)
	{
		if (type == LUMA_8x8)
		{
			Enc_MBs[(p_Vid->width / 16) * currMB->mb_y + currMB->mb_x].cbp_bits[0] |= ((int64)0x33 << bit);
		}
		else if (type == LUMA_8x4)
		{
			Enc_MBs[(p_Vid->width / 16) * currMB->mb_y + currMB->mb_x].cbp_bits[0] |= ((int64)0x03 << bit);
		}
		else if (type == LUMA_4x8)
		{
			Enc_MBs[(p_Vid->width / 16) * currMB->mb_y + currMB->mb_x].cbp_bits[0] |= ((int64)0x11 << bit);
		}
		else
		{
			Enc_MBs[(p_Vid->width / 16) * currMB->mb_y + currMB->mb_x].cbp_bits[0] |= ((int64)0x01 << bit);
		}
	}

	bit = (y_dc ? 0 : y_ac ? 1 : u_dc ? 17 : v_dc ? 18 : u_ac ? 19 : 35);

	if (type != LUMA_8x8)
	{
		if (block_b.available)
		{
			if (Enc_MBs[block_b.mb_addr].mb_type == IPCM)
				upper_bit = 1;
			else
				upper_bit = get_bit(Enc_MBs[block_b.mb_addr].cbp_bits[0], bit + bit_pos_b);
		}


		if (block_a.available)
		{
			if (Enc_MBs[block_a.mb_addr].mb_type == IPCM)
				left_bit = 1;
			else
				left_bit = get_bit(Enc_MBs[block_a.mb_addr].cbp_bits[0], bit + bit_pos_a);
		}

		ctx = (upper_bit << 1) + left_bit;


		//===== encode symbol =====
		bct = tex_ctx->bcbp_contexts[type2ctx_bcbp[type]] + ctx;
		biari_encode_symbol(eep_dp, cbp_bit, bct);
	}
}

/*!
****************************************************************************
* \brief
*    Write Significance MAP
****************************************************************************
*/
void write_significance_map(Enc_Macroblock* currMB, EncodingEnvironmentPtr eep_dp, int type, int coeff[], int coeff_ctr, TextureInfoContexts*  tex_ctx)
{
	int   k;
	int   sig, last;
	int   k0 = 0;
	int   k1 = maxpos[type];

#if ENABLE_FIELD_CTX
	Enc_VideoParameters *p_Vid = currMB->p_Vid;
	int               fld = (p_Vid->structure != FRAME || currMB->mb_field);
#else
	int               fld = 0;
#endif
	BiContextTypePtr  map_ctx = tex_ctx->map_contexts[fld][type2ctx_map[type]];
	BiContextTypePtr  last_ctx = tex_ctx->last_contexts[fld][type2ctx_last[type]];
	const byte *pos2ctxmap = fld ? pos2ctx_map_int[type] : pos2ctx_map[type];
	const byte *pos2ctxlast = pos2ctx_last[type];

	if (!c1isdc[type])
	{
		++k0;
		++k1;
		--coeff;
	}

	for (k = k0; k<k1; ++k) // if last coeff is reached, it has to be significant
	{
		sig = (coeff[k] != 0);
		biari_encode_symbol(eep_dp, sig, map_ctx + pos2ctxmap[k]);
		if (sig)
		{
			last = (--coeff_ctr == 0);

			biari_encode_symbol(eep_dp, last, last_ctx + pos2ctxlast[k]);
			if (last)
				return;
		}
	}
}

/*!
************************************************************************
* \brief
*    Arithmetic encoding of one binary symbol assuming
*    a fixed prob. distribution with p(symbol) = 0.5
************************************************************************
*/
void biari_encode_symbol_eq_prob(EncodingEnvironmentPtr eep, int symbol)
{
	unsigned int low = eep->Elow;
	--(eep->Ebits_to_go);
	++(eep->C);

	if (symbol != 0)
	{
		low += eep->Erange << eep->Ebits_to_go;
		if (low >= 0x04000000) // output of carry needed
		{
			low -= 0x04000000;
			propagate_carry(eep);
		}
	}
	if (eep->Ebits_to_go == 0)  // renorm needed
	{
		eep->Elow = (low << 16)& (0x03FFFFFF);
		low = (low >> 10) & 0xFFFF; // mask out the 8/16 MSBs for output
		if (low < 0xFFFF)      // no carry possible, output now
		{
			put_last_chunk_plus_outstanding(eep, low);
		}
		else          // low == "FF"; keep it, may affect future carry
		{
			++(eep->Echunks_outstanding);
		}

		eep->Ebits_to_go = 16;
		return;
	}
	else         // no renorm needed
	{
		eep->Elow = low;
		return;
	}
}

/*!
************************************************************************
* \brief
*    Exp Golomb binarization and encoding
************************************************************************
*/
static void exp_golomb_encode_eq_prob(EncodingEnvironmentPtr eep_dp, unsigned int symbol, int k)
{
	for (;;)
	{
		if (symbol >= (unsigned int)(1 << k))
		{
			biari_encode_symbol_eq_prob(eep_dp, 1);   //first unary part
			symbol = symbol - (1 << k);
			k++;
		}
		else
		{
			biari_encode_symbol_eq_prob(eep_dp, 0);   //now terminated zero of unary part
			while (k--)                               //next binary part
				biari_encode_symbol_eq_prob(eep_dp, ((symbol >> k) & 1));
			break;
		}
	}
}

/*!
************************************************************************
* \brief
*    Exp-Golomb for Level Encoding
*
************************************************************************/
void unary_exp_golomb_level_encode(EncodingEnvironmentPtr eep_dp,
	unsigned int symbol,
	BiContextTypePtr ctx)
{
	if (symbol == 0)
	{
		biari_encode_symbol(eep_dp, 0, ctx);
		return;
	}
	else
	{
		unsigned int l = symbol;
		unsigned int k = 1;

		biari_encode_symbol(eep_dp, 1, ctx);
		while (((--l)>0) && (++k <= 13))
			biari_encode_symbol(eep_dp, 1, ctx);
		if (symbol < 13)
			biari_encode_symbol(eep_dp, 0, ctx);
		else
			exp_golomb_encode_eq_prob(eep_dp, symbol - 13, 0);
	}
}

/*!
****************************************************************************
* \brief
*    Write Levels
****************************************************************************
*/
void write_significant_coefficients(EncodingEnvironmentPtr eep_dp, int type, int coeff[], TextureInfoContexts*  tex_ctx, int coeff_cnt)
{
	BiContextType *one_contexts = tex_ctx->one_contexts[type2ctx_one[type]];
	BiContextType *abs_contexts = tex_ctx->abs_contexts[type2ctx_abs[type]];
	int   absLevel;
	int   ctx;
	short sign;
	int greater_one;
	int   c1 = 1;
	int   c2 = 0;
	int   i;

	for (i = maxpos[type]; (i >= 0) && (coeff_cnt); i--)
	{
		if (coeff[i] != 0)
		{
			coeff_cnt--;

			if (coeff[i] > 0)
			{
				absLevel = coeff[i];
				sign = 0;
			}
			else
			{
				absLevel = -coeff[i];
				sign = 1;
			}

			greater_one = (absLevel > 1);

			//--- if coefficient is one ---
			ctx = imin(c1, 4);
			biari_encode_symbol(eep_dp, greater_one, one_contexts + ctx);

			if (greater_one)
			{
				ctx = imin(c2++, max_c2[type]);
				unary_exp_golomb_level_encode(eep_dp, absLevel - 2, abs_contexts + ctx);
				c1 = 0;
			}
			else if (c1)
			{
				c1++;
			}
			biari_encode_symbol_eq_prob(eep_dp, sign);
		}
	}
}

/*!
****************************************************************************
* \brief
*    Write Block-Transform Coefficients
****************************************************************************
*/
void writeRunLevel_CABAC(Enc_Macroblock* currMB, Enc_SyntaxElement *se, Enc_DataPartition *dp)
{
	Enc_Slice *currSlice = currMB->p_Slice;
	EncodingEnvironmentPtr eep_dp = &(dp->ee_cabac);

	//--- accumulate run-level information ---
	if (se->value1 != 0)
	{
		currSlice->pos += se->value2;
		currSlice->coeff[currSlice->pos++] = se->value1;
		++currSlice->coeff_ctr;
	}
	else
	{
		TextureInfoContexts *tex_ctx = currSlice->tex_ctx;
		//===== encode CBP-BIT =====
		if (currSlice->coeff_ctr > 0)
		{
			//currSlice->write_and_store_CBP_block_bit  (currMB, eep_dp, se->context, 1, tex_ctx);
			write_and_store_CBP_block_bit(currMB, eep_dp, se->context, 1, tex_ctx);
			//===== encode significance map =====
			write_significance_map(currMB, eep_dp, se->context, currSlice->coeff, currSlice->coeff_ctr, tex_ctx);
			//===== encode significant coefficients =====
			write_significant_coefficients(eep_dp, se->context, currSlice->coeff, tex_ctx, currSlice->coeff_ctr);

			memset(currSlice->coeff, 0, 64 * sizeof(int));
		}
		else
		{
			//currSlice->write_and_store_CBP_block_bit  (currMB, eep_dp, se->context, 0, tex_ctx);
			write_and_store_CBP_block_bit(currMB, eep_dp, se->context, 0, tex_ctx);
		}

		//--- reset counters ---
		currSlice->pos = currSlice->coeff_ctr = 0;
	}

	//dp->bitstream->write_flag = 1;
}

/*!
************************************************************************
* \brief
*    Writes Coeffs of an 4x4 block
************************************************************************
*/
int writeCoeff4x4_CABAC(Enc_Macroblock* currMB, ColorPlane plane, int iX, int iY, int *ACLevel, int *ACRun)
{
	Enc_Slice* currSlice = currMB->p_Slice;

	int             rate = 0;
	Enc_SyntaxElement   se;
	const int*      partMap = assignSE2partition[currSlice->partition_mode];
	Enc_DataPartition*  dataPart;

	int level;
	int k;

	currMB->subblock_x = iX;
	currMB->subblock_y = iY;
	se.context = ((plane == 0) ? (LUMA_4x4) : ((plane == 1) ? CB_4x4 : CR_4x4));

	// DC
	se.type = (currMB->is_intra_block ? SE_LUM_DC_INTRA : SE_LUM_DC_INTER);

	// choose the appropriate data partition
	dataPart = &(currSlice->partArr[partMap[se.type]]);

	level = se.value1 = ACLevel[0]; // level
	se.value2 = ACRun[0]; // run
	writeRunLevel_CABAC(currMB, &se, dataPart);
	rate += se.len;

	// AC Coefficients
	se.type = (currMB->is_intra_block ? SE_LUM_AC_INTRA : SE_LUM_AC_INTER);
	// choose the appropriate data partition
	dataPart = &(currSlice->partArr[partMap[se.type]]);

	for (k = 1; k <= 16 && level != 0; k++)
	{
		level = se.value1 = ACLevel[k]; // level
		se.value2 = ACRun[k]; // run    
		writeRunLevel_CABAC(currMB, &se, dataPart);

		rate += se.len;
	}

	return rate;
}
void initData(Macroblock *currMB, Enc_Macroblock *currMBI)
{
	Enc_MB->p_Slice->partition_mode = 0;
	Enc_MB->p_Slice->pos = 0;
	Enc_MB->p_Slice->coeff_ctr = 0;

	currMBI->p_Vid->is_v_block = 0;
	currMBI->is_intra_block = currMB->is_intra_block;
	currMBI->block_x = currMB->block_x;
	currMBI->block_y = currMB->block_y;
	currMBI->mb_x = currMB->mb.x;
	currMBI->mb_y = currMB->mb.y;
	currMBI->p_Vid->mb_size[0][0] = currMB->p_Vid->mb_size[0][0];
	currMBI->p_Vid->mb_size[0][1] = currMB->p_Vid->mb_size[0][1];
	currMBI->p_Vid->mb_size[1][0] = currMB->p_Vid->mb_size[1][0];
	currMBI->p_Vid->mb_size[1][1] = currMB->p_Vid->mb_size[1][1];
	currMBI->p_Vid->mb_size[2][0] = currMB->p_Vid->mb_size[2][0];
	currMBI->p_Vid->mb_size[2][1] = currMB->p_Vid->mb_size[2][1];
	currMBI->p_Vid->yuv_format = (ColorFormat)currMB->p_Vid->yuv_format;	//YUV420
	Enc_MB->mbAddrA = currMB->mbAddrA;
	Enc_MB->mbAddrB = currMB->mbAddrB;
	Enc_MB->mbAddrC = currMB->mbAddrC;
	Enc_MB->mbAddrD = currMB->mbAddrD;
	Enc_MB->mbAddrX = currMB->mbAddrX;
	Enc_MB->mbAvailA = currMB->mbAvailA;
	Enc_MB->mbAvailB = currMB->mbAvailB;
	Enc_MB->mbAvailC = currMB->mbAvailC;
	Enc_MB->mbAvailD = currMB->mbAvailD;
	Enc_MB->DeblockCall = currMB->DeblockCall;
	Enc_MB->mb_type = currMB->mb_type;
	Enc_MB->p_Vid->PicPos->x = currMB->p_Vid->PicPos->x;
	Enc_MB->p_Vid->PicPos->y = currMB->p_Vid->PicPos->y;
	currMBI->mb_field = currMB->mb_field;
	currMBI->DeblockCall = currMB->DeblockCall;
	currMBI->p_Vid->structure = currMB->p_Vid->structure;
}
int RunLevelAC(int block_y, int block_x, int *level, int *run, int *cof, int *PLNZ, int Y, int X)
{
	int coeffcount = 0;
	int runtemp = 0, i;
	int tempC;
	for (i = 0; i < 16; i++)
	{
		level[i] = 0;
		run[i] = 0;
		tempC = cof[block_y * 32 + block_x * 16 + i];
		if (tempC)
		{
			level[coeffcount] = tempC;
			run[coeffcount++] = runtemp;
			runtemp = 0;
			*PLNZ = i + 1;
		}
		else
		{
			runtemp++;
		}
	}
	return coeffcount;
}/*
 int RunLevelAC(int block_y, int block_x, int *level, int *run, int **cof, int *PLNZ, int Y, int X)
 {
 int coeffcount = 0;
 int runtemp = 0, i;
 int tempC;
 int scanP[16];	//xy
 scanP[0] = 00; scanP[4] = 11; scanP[8] = 12; scanP[12] = 31;
 scanP[1] = 10; scanP[5] = 20; scanP[9] = 03; scanP[13] = 32;
 scanP[2] = 01; scanP[6] = 30; scanP[10] = 13; scanP[14] = 23;
 scanP[3] = 02; scanP[7] = 21; scanP[11] = 22; scanP[15] = 33;
 for (i = 0; i < 16; i++)
 {
 level[i] = 0;
 run[i] = 0;
 tempC = cof[(scanP[i] % 10) + (4 * block_y)][(scanP[i] / 10) + (4 * block_x) + X];
 if (tempC)
 {
 level[coeffcount] = tempC;
 run[coeffcount++] = runtemp;
 runtemp = 0;
 *PLNZ = i + 1;
 }
 else
 {
 runtemp++;
 }
 }
 return coeffcount;
 }*/
 /*
 ************************************************************************
 * \brief
 *    Exp-Golomb for MV Encoding
 *
 ************************************************************************/
void unary_exp_golomb_mv_encode(EncodingEnvironmentPtr eep_dp,
	unsigned int symbol,
	BiContextTypePtr ctx,
	unsigned int max_bin)
{
	if (symbol == 0)
	{
		biari_encode_symbol(eep_dp, 0, ctx);
		return;
	}
	else
	{
		unsigned int bin = 1;
		unsigned int l = symbol, k = 1;
		biari_encode_symbol(eep_dp, 1, ctx++);

		while (((--l)>0) && (++k <= 8))
		{
			biari_encode_symbol(eep_dp, 1, ctx);
			if ((++bin) == 2)
				++ctx;
			if (bin == max_bin)
				++ctx;
		}
		if (symbol < 8)
			biari_encode_symbol(eep_dp, 0, ctx);
		else
			exp_golomb_encode_eq_prob(eep_dp, symbol - 8, 3);
	}
}
/*!
************************************************************************
* \brief
*    Arithmetic encoding for last symbol before termination
************************************************************************
*/
void biari_encode_symbol_final(EncodingEnvironmentPtr eep, int symbol)
{
	unsigned int range = eep->Erange - 2;
	unsigned int low = eep->Elow;
	int bl = eep->Ebits_to_go;

	++(eep->C);

	if (symbol == 0) // MPS
	{
		if (range >= 0x0100) // no renorm
		{
			eep->Erange = range;
			return;
		}
		else
		{
			range <<= 1;
			if (--bl > 0)  // renorm once, no output
			{
				eep->Erange = range;
				eep->Ebits_to_go = bl;
				return;
			}
		}
	}
	else     // LPS
	{
		low += (range << bl);
		range = 2;

		if (low >= 0x04000000) // output of carry needed
		{
			low -= 0x04000000; // remove MSB, i.e., carry bit
			propagate_carry(eep);
		}
		bl -= 7; // 7 left shifts needed to renormalize

		range <<= 7;
		if (bl > 0)
		{
			eep->Erange = range;
			eep->Elow = low;
			eep->Ebits_to_go = bl;
			return;
		}
	}


	//renorm needed
	eep->Elow = (low << 16) & (0x03FFFFFF);
	low = (low >> 10) & 0xFFFF; // mask out the 8/16 MSBs
	if (low < 0xFFFF)
	{  // no carry possible, output now
		put_last_chunk_plus_outstanding(eep, low);
	}
	else
	{  // low == "FF"; keep it, may affect future carry
		++(eep->Echunks_outstanding);
	}

	eep->Erange = range;
	bl += 16;
	eep->Ebits_to_go = bl;
}
void unary_bin_encode(EncodingEnvironmentPtr eep_dp,
	unsigned int symbol,
	BiContextTypePtr ctx,
	int ctx_offset)
{
	if (symbol == 0)
	{
		biari_encode_symbol(eep_dp, 0, ctx);
		return;
	}
	else
	{
		biari_encode_symbol(eep_dp, 1, ctx);
		ctx += ctx_offset;
		while ((--symbol) > 0)
			biari_encode_symbol(eep_dp, 1, ctx);
		biari_encode_symbol(eep_dp, 0, ctx);
	}
}
void put_one_byte_final(EncodingEnvironmentPtr eep, unsigned int b)
{
	eep->Ecodestrm[(*eep->Ecodestrm_len)++] = (byte)b;
}

void put_one_byte(EncodingEnvironmentPtr eep, int b)
{
	if (eep->Epbuf == 3)
	{
		put_buffer(eep);
	}
	eep->Ebuffer <<= 8;
	eep->Ebuffer += (b);

	++(eep->Epbuf);
}

void put_last_chunk_plus_outstanding_final(EncodingEnvironmentPtr eep, unsigned int l)
{
	while (eep->Echunks_outstanding > 0)
	{
		put_one_word(eep, 0xFFFF);
		--(eep->Echunks_outstanding);
	}
	put_one_byte(eep, l);
}

void arienco_done_encoding(Enc_Macroblock *currMB, EncodingEnvironmentPtr eep)
{
	unsigned int low = eep->Elow;
	int remaining_bits = 16 - eep->Ebits_to_go; // output (2 + remaining) bits for terminating the codeword + one stop bit
	unsigned char mask;
	BitCounter *mbBits = &currMB->bits;

	//p_Vid->pic_bin_count += eep->E*8 + eep->C; // no of processed bins

	if (remaining_bits <= 5) // one terminating byte 
	{
		mbBits->mb_stuffing += (5 - remaining_bits);
		mask = (unsigned char)(255 - ((1 << (6 - remaining_bits)) - 1));
		low = (low >> (MASK_BITS)) & mask; // mask out the (2+remaining_bits) MSBs
		low += (32 >> remaining_bits);       // put the terminating stop bit '1'

		put_last_chunk_plus_outstanding_final(eep, low);
		put_buffer(eep);
	}
	else if (remaining_bits <= 13)            // two terminating bytes
	{
		mbBits->mb_stuffing += (13 - remaining_bits);
		put_last_chunk_plus_outstanding_final(eep, ((low >> (MASK_BITS)) & 0xFF)); // mask out the 8 MSBs for output

		put_buffer(eep);
		if (remaining_bits > 6)
		{
			mask = (unsigned char)(255 - ((1 << (14 - remaining_bits)) - 1));
			low = (low >> (B_BITS)) & mask;
			low += (0x2000 >> remaining_bits);     // put the terminating stop bit '1'
			put_one_byte_final(eep, low);
		}
		else
		{
			put_one_byte_final(eep, 128); // second byte contains terminating stop bit '1' only
		}
	}
	else             // three terminating bytes
	{
		put_last_chunk_plus_outstanding(eep, ((low >> B_BITS) & B_LOAD_MASK)); // mask out the 16 MSBs for output
		put_buffer(eep);
		mbBits->mb_stuffing += (21 - remaining_bits);

		if (remaining_bits > 14)
		{
			mask = (unsigned char)(255 - ((1 << (22 - remaining_bits)) - 1));
			low = (low >> (MAX_BITS - 24)) & mask;
			low += (0x200000 >> remaining_bits);       // put the terminating stop bit '1'
			put_one_byte_final(eep, low);
		}
		else
		{
			put_one_byte_final(eep, 128); // third byte contains terminating stop bit '1' only
		}
	}
	eep->Ebits_to_go = 8;
}
int Embedding_CABAC(int *point, char *StringI, int *lv, int *rn, int numcoeff, int PLNZ)
{
	int ret = 0;
	unsigned char Si = (StringI[(*point) >> 3] >> (7 - (*point) & 7)) & 0x1;

	if (((PLNZ & 0x1) == 0) && Si)		//even
	{
		ret = 1;
		if (PLNZ == 16)				//Pim - 1
		{
			if (rn[numcoeff - 1] != 0)
				rn[numcoeff - 1]--;
			else
				lv[numcoeff - 1] = 0;
		}
		else						//Pi + 1
			rn[numcoeff - 1] += 1;
	}
	else if (((PLNZ & 0x1) != 0) && (Si == 0))				//odd		//Pi + 1
	{
		ret = 1;
		rn[numcoeff - 1] += 1;
	}
	*point += 1;
	return ret;
}