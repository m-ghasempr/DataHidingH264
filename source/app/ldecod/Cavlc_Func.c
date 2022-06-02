#include "Cavlc_Func.h"

int is_intra(Enc_Macroblock *curr_MB)
{
	return ((curr_MB)->mb_type == SI4MB || (curr_MB)->mb_type == I4MB || (curr_MB)->mb_type == I16MB || (curr_MB)->mb_type == I8MB || (curr_MB)->mb_type == IPCM);
}

/*!
************************************************************************
* \brief
*    Makes code word and passes it back
*
* \par Input:
*    Info   : Xn..X2 X1 X0                                             \n
*    Length : Total number of bits in the codeword
************************************************************************
*/

int symbol2vlc(Enc_SyntaxElement *sym)
{
	int info_len = sym->len;

	// Convert info into a bitpattern int
	sym->bitpattern = 0;

	// vlc coding
	while (--info_len >= 0)
	{
		sym->bitpattern <<= 1;
		sym->bitpattern |= (0x01 & (sym->inf >> info_len));
	}
	return 0;
}

/*!
************************************************************************
* \brief
*    writes UVLC code to the appropriate buffer
************************************************************************
*/
void writeUVLC2buffer(Enc_SyntaxElement *se, Enc_Bitstream *currStream)
{
	unsigned int mask = 1 << (se->len - 1);
	byte *byte_buf = &currStream->byte_buf;
	int *bits_to_go = &currStream->bits_to_go;
	int i;

	// Add the new bits to the bitstream.
	// Write out a byte if it is full
	if (se->len < 33)
	{
		for (i = 0; i < se->len; i++)
		{
			*byte_buf <<= 1;

			if (se->bitpattern & mask)
				*byte_buf |= 1;

			mask >>= 1;

			if ((--(*bits_to_go)) == 0)
			{
				*bits_to_go = 8;
				currStream->streamBuffer[currStream->byte_pos++] = *byte_buf;
				*byte_buf = 0;
			}
		}
	}
	else
	{
		// zeros
		for (i = 0; i < (se->len - 32); i++)
		{
			*byte_buf <<= 1;

			if ((--(*bits_to_go)) == 0)
			{
				*bits_to_go = 8;
				currStream->streamBuffer[currStream->byte_pos++] = *byte_buf;
				*byte_buf = 0;
			}
		}
		// actual info
		mask = (unsigned int)1 << 31;
		for (i = 0; i < 32; i++)
		{
			*byte_buf <<= 1;

			if (se->bitpattern & mask)
				*byte_buf |= 1;

			mask >>= 1;

			if ((--(*bits_to_go)) == 0)
			{
				*bits_to_go = 8;
				currStream->streamBuffer[currStream->byte_pos++] = *byte_buf;
				*byte_buf = 0;
			}
		}
	}
}

/*!
************************************************************************
* \brief
*    write VLC for Run Before Next Coefficient, VLC0
************************************************************************
*/
int writeSyntaxElement_Run(Enc_SyntaxElement *se, Enc_DataPartition *dp)
{
	static const byte lentab[TOTRUN_NUM][16] =
	{
		{ 1, 1 },
		{ 1, 2, 2 },
		{ 2, 2, 2, 2 },
		{ 2, 2, 2, 3, 3 },
		{ 2, 2, 3, 3, 3, 3 },
		{ 2, 3, 3, 3, 3, 3, 3 },
		{ 3, 3, 3, 3, 3, 3, 3, 4, 5, 6, 7, 8, 9, 10, 11 },
	};

	static const byte codtab[TOTRUN_NUM][16] =
	{
		{ 1, 0 },
		{ 1, 1, 0 },
		{ 3, 2, 1, 0 },
		{ 3, 2, 1, 1, 0 },
		{ 3, 2, 3, 2, 1, 0 },
		{ 3, 0, 1, 3, 2, 5, 4 },
		{ 7, 6, 5, 4, 3, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
	};
	int vlcnum = se->len;

	// se->value1 : run
	se->len = lentab[vlcnum][se->value1];
	se->inf = codtab[vlcnum][se->value1];

	if (se->len == 0)
	{
		printf("ERROR: (run) not valid: (%d)\n", se->value1);
		exit(-1);
	}

	symbol2vlc(se);

	writeUVLC2buffer(se, dp->bitstream);

	if (se->type != SE_HEADER)
		dp->bitstream->write_flag = 1;

#if TRACE
	if (dp->bitstream->trace_enabled)
		trace2out(se);
#endif

	return (se->len);
}

/*!
************************************************************************
* \brief
*    write VLC for TotalZeros
************************************************************************
*/
int writeSyntaxElement_TotalZeros(Enc_SyntaxElement *se, Enc_DataPartition *dp)
{
	static const byte lentab[TOTRUN_NUM][16] =
	{
		{ 1, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 9 },
		{ 3, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 6, 6, 6, 6 },
		{ 4, 3, 3, 3, 4, 4, 3, 3, 4, 5, 5, 6, 5, 6 },
		{ 5, 3, 4, 4, 3, 3, 3, 4, 3, 4, 5, 5, 5 },
		{ 4, 4, 4, 3, 3, 3, 3, 3, 4, 5, 4, 5 },
		{ 6, 5, 3, 3, 3, 3, 3, 3, 4, 3, 6 },
		{ 6, 5, 3, 3, 3, 2, 3, 4, 3, 6 },
		{ 6, 4, 5, 3, 2, 2, 3, 3, 6 },
		{ 6, 6, 4, 2, 2, 3, 2, 5 },
		{ 5, 5, 3, 2, 2, 2, 4 },
		{ 4, 4, 3, 3, 1, 3 },
		{ 4, 4, 2, 1, 3 },
		{ 3, 3, 1, 2 },
		{ 2, 2, 1 },
		{ 1, 1 },
	};

	static const byte codtab[TOTRUN_NUM][16] =
	{
		{ 1, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2, 1 },
		{ 7, 6, 5, 4, 3, 5, 4, 3, 2, 3, 2, 3, 2, 1, 0 },
		{ 5, 7, 6, 5, 4, 3, 4, 3, 2, 3, 2, 1, 1, 0 },
		{ 3, 7, 5, 4, 6, 5, 4, 3, 3, 2, 2, 1, 0 },
		{ 5, 4, 3, 7, 6, 5, 4, 3, 2, 1, 1, 0 },
		{ 1, 1, 7, 6, 5, 4, 3, 2, 1, 1, 0 },
		{ 1, 1, 5, 4, 3, 3, 2, 1, 1, 0 },
		{ 1, 1, 1, 3, 3, 2, 2, 1, 0 },
		{ 1, 0, 1, 3, 2, 1, 1, 1, },
		{ 1, 0, 1, 3, 2, 1, 1, },
		{ 0, 1, 1, 2, 1, 3 },
		{ 0, 1, 1, 1, 1 },
		{ 0, 1, 1, 1 },
		{ 0, 1, 1 },
		{ 0, 1 },
	};
	int vlcnum = se->len;

	// se->value1 : TotalZeros
	se->len = lentab[vlcnum][se->value1];
	se->inf = codtab[vlcnum][se->value1];

	if (se->len == 0)
	{
		printf("ERROR: (TotalZeros) not valid: (%d)\n", se->value1);
		exit(-1);
	}

	symbol2vlc(se);

	writeUVLC2buffer(se, dp->bitstream);

	if (se->type != SE_HEADER)
		dp->bitstream->write_flag = 1;

#if TRACE
	if (dp->bitstream->trace_enabled)
		trace2out(se);
#endif

	return (se->len);
}


/*!
************************************************************************
* \brief
*    write VLC for TotalZeros for Chroma DC
************************************************************************
*/
int writeSyntaxElement_TotalZerosChromaDC(Enc_VideoParameters *p_Vid, Enc_SyntaxElement *se, Enc_DataPartition *dp)
{
	static const byte lentab[3][TOTRUN_NUM][16] =
	{
		//YUV420
		{ { 1, 2, 3, 3 },
		{ 1, 2, 2 },
		{ 1, 1 } },
		//YUV422
		{ { 1, 3, 3, 4, 4, 4, 5, 5 },
		{ 3, 2, 3, 3, 3, 3, 3 },
		{ 3, 3, 2, 2, 3, 3 },
		{ 3, 2, 2, 2, 3 },
		{ 2, 2, 2, 2 },
		{ 2, 2, 1 },
		{ 1, 1 } },
		//YUV444
		{ { 1, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 9 },
		{ 3, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 6, 6, 6, 6 },
		{ 4, 3, 3, 3, 4, 4, 3, 3, 4, 5, 5, 6, 5, 6 },
		{ 5, 3, 4, 4, 3, 3, 3, 4, 3, 4, 5, 5, 5 },
		{ 4, 4, 4, 3, 3, 3, 3, 3, 4, 5, 4, 5 },
		{ 6, 5, 3, 3, 3, 3, 3, 3, 4, 3, 6 },
		{ 6, 5, 3, 3, 3, 2, 3, 4, 3, 6 },
		{ 6, 4, 5, 3, 2, 2, 3, 3, 6 },
		{ 6, 6, 4, 2, 2, 3, 2, 5 },
		{ 5, 5, 3, 2, 2, 2, 4 },
		{ 4, 4, 3, 3, 1, 3 },
		{ 4, 4, 2, 1, 3 },
		{ 3, 3, 1, 2 },
		{ 2, 2, 1 },
		{ 1, 1 } }
	};

	static const byte codtab[3][TOTRUN_NUM][16] =
	{
		//YUV420
		{ { 1, 1, 1, 0 },
		{ 1, 1, 0 },
		{ 1, 0 } },
		//YUV422
		{ { 1, 2, 3, 2, 3, 1, 1, 0 },
		{ 0, 1, 1, 4, 5, 6, 7 },
		{ 0, 1, 1, 2, 6, 7 },
		{ 6, 0, 1, 2, 7 },
		{ 0, 1, 2, 3 },
		{ 0, 1, 1 },
		{ 0, 1 } },
		//YUV444
		{ { 1, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2, 1 },
		{ 7, 6, 5, 4, 3, 5, 4, 3, 2, 3, 2, 3, 2, 1, 0 },
		{ 5, 7, 6, 5, 4, 3, 4, 3, 2, 3, 2, 1, 1, 0 },
		{ 3, 7, 5, 4, 6, 5, 4, 3, 3, 2, 2, 1, 0 },
		{ 5, 4, 3, 7, 6, 5, 4, 3, 2, 1, 1, 0 },
		{ 1, 1, 7, 6, 5, 4, 3, 2, 1, 1, 0 },
		{ 1, 1, 5, 4, 3, 3, 2, 1, 1, 0 },
		{ 1, 1, 1, 3, 3, 2, 2, 1, 0 },
		{ 1, 0, 1, 3, 2, 1, 1, 1, },
		{ 1, 0, 1, 3, 2, 1, 1, },
		{ 0, 1, 1, 2, 1, 3 },
		{ 0, 1, 1, 1, 1 },
		{ 0, 1, 1, 1 },
		{ 0, 1, 1 },
		{ 0, 1 } }
	};
	int vlcnum = se->len;
	int yuv = p_Vid->yuv_format - 1;

	// se->value1 : TotalZeros
	se->len = lentab[yuv][vlcnum][se->value1];
	se->inf = codtab[yuv][vlcnum][se->value1];

	if (se->len == 0)
	{
		printf("ERROR: (TotalZeros) not valid: (%d)\n", se->value1);
		exit(-1);
	}

	symbol2vlc(se);

	writeUVLC2buffer(se, dp->bitstream);

	if (se->type != SE_HEADER)
		dp->bitstream->write_flag = 1;

#if TRACE
	if (dp->bitstream->trace_enabled)
		trace2out(se);
#endif

	return (se->len);
}

/*!
************************************************************************
* \brief
*    write VLC for Coeff Level (VLC1)
************************************************************************
*/
int writeSyntaxElement_Level_VLC1(Enc_SyntaxElement *se, Enc_DataPartition *dp, int profile_idc)
{
	int level = se->value1;
	int sign = (level < 0 ? 1 : 0);
	int levabs = iabs(level);

	if (levabs < 8)
	{
		se->len = levabs * 2 + sign - 1;
		se->inf = 1;
	}
	else if (levabs < 16)
	{
		// escape code1
		se->len = 19;
		se->inf = 16 | ((levabs << 1) - 16) | sign;
	}
	else
	{
		int iMask = 4096, numPrefix = 0;
		int levabsm16 = levabs + 2032;

		// escape code2
		if ((levabsm16) >= 4096)
		{
			numPrefix++;
			while ((levabsm16) >= (4096 << numPrefix))
			{
				numPrefix++;
			}
		}

		iMask <<= numPrefix;
		se->inf = iMask | ((levabsm16 << 1) - iMask) | sign;

		/* Assert to make sure that the code fits in the VLC */
		/* make sure that we are in High Profile to represent level_prefix > 15 */
		if (numPrefix > 0 && !is_FREXT_profile(profile_idc))
		{
			////error( "level_prefix must be <= 15 except in High Profile\n",  1000 );
			se->len = 0x0000FFFF; // This can be some other big number
			return (se->len);
		}

		se->len = 28 + (numPrefix << 1);
	}

	symbol2vlc(se);

	writeUVLC2buffer(se, dp->bitstream);

	if (se->type != SE_HEADER)
		dp->bitstream->write_flag = 1;

#if TRACE
	if (dp->bitstream->trace_enabled)
		trace2out(se);
#endif

	return (se->len);
}


/*!
************************************************************************
* \brief
*    write VLC for Coeff Level
************************************************************************
*/
int writeSyntaxElement_Level_VLCN(Enc_SyntaxElement *se, int vlc, Enc_DataPartition *dp, int profile_idc)
{
	int level = se->value1;
	int sign = (level < 0 ? 1 : 0);
	int levabs = iabs(level) - 1;

	int shift = vlc - 1;
	int escape = (15 << shift);

	if (levabs < escape)
	{
		int sufmask = ~((0xffffffff) << shift);
		int suffix = (levabs)& sufmask;

		se->len = ((levabs) >> shift) + 1 + vlc;
		se->inf = (2 << shift) | (suffix << 1) | sign;
	}
	else
	{
		int iMask = 4096;
		int levabsesc = levabs - escape + 2048;
		int numPrefix = 0;

		if ((levabsesc) >= 4096)
		{
			numPrefix++;
			while ((levabsesc) >= (4096 << numPrefix))
			{
				numPrefix++;
			}
		}

		iMask <<= numPrefix;
		se->inf = iMask | ((levabsesc << 1) - iMask) | sign;

		/* Assert to make sure that the code fits in the VLC */
		/* make sure that we are in High Profile to represent level_prefix > 15 */
		if (numPrefix > 0 && !is_FREXT_profile(profile_idc))
		{
			////error( "level_prefix must be <= 15 except in High Profile\n",  1000 );
			se->len = 0x0000FFFF; // This can be some other big number
			return (se->len);
		}
		se->len = 28 + (numPrefix << 1);
	}

	symbol2vlc(se);

	writeUVLC2buffer(se, dp->bitstream);

	if (se->type != SE_HEADER)
		dp->bitstream->write_flag = 1;

#if TRACE
	if (dp->bitstream->trace_enabled)
		trace2out(se);
#endif

	return (se->len);
}

/*!
************************************************************************
* \brief
*    generates VLC code and passes the codeword to the buffer
************************************************************************
*/
int writeSyntaxElement_VLC(Enc_SyntaxElement *se, Enc_DataPartition *dp)
{
	se->inf = se->value1;
	se->len = se->value2;
	symbol2vlc(se);

	writeUVLC2buffer(se, dp->bitstream);

	if (se->type != SE_HEADER)
		dp->bitstream->write_flag = 1;

#if TRACE
	if (dp->bitstream->trace_enabled)
		trace2out(se);
#endif

	return (se->len);
}

/*!
************************************************************************
* \brief
*    write VLC for NumCoeff and TrailingOnes for Chroma DC
************************************************************************
*/
int writeSyntaxElement_NumCoeffTrailingOnesChromaDC(Enc_VideoParameters *p_Vid, Enc_SyntaxElement *se, Enc_DataPartition *dp)
{
	static const byte lentab[3][4][17] =
	{
		//YUV420
		{ { 2, 6, 6, 6, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
		{ 0, 1, 6, 7, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
		{ 0, 0, 3, 7, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
		{ 0, 0, 0, 6, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 } },
		//YUV422
		{ { 1, 7, 7, 9, 9, 10, 11, 12, 13, 0, 0, 0, 0, 0, 0, 0, 0 },
		{ 0, 2, 7, 7, 9, 10, 11, 12, 12, 0, 0, 0, 0, 0, 0, 0, 0 },
		{ 0, 0, 3, 7, 7, 9, 10, 11, 12, 0, 0, 0, 0, 0, 0, 0, 0 },
		{ 0, 0, 0, 5, 6, 7, 7, 10, 11, 0, 0, 0, 0, 0, 0, 0, 0 } },
		//YUV444
		{ { 1, 6, 8, 9, 10, 11, 13, 13, 13, 14, 14, 15, 15, 16, 16, 16, 16 },
		{ 0, 2, 6, 8, 9, 10, 11, 13, 13, 14, 14, 15, 15, 15, 16, 16, 16 },
		{ 0, 0, 3, 7, 8, 9, 10, 11, 13, 13, 14, 14, 15, 15, 16, 16, 16 },
		{ 0, 0, 0, 5, 6, 7, 8, 9, 10, 11, 13, 14, 14, 15, 15, 16, 16 } }
	};

	static const byte codtab[3][4][17] =
	{
		//YUV420
		{ { 1, 7, 4, 3, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
		{ 0, 1, 6, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
		{ 0, 0, 1, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
		{ 0, 0, 0, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 } },
		//YUV422
		{ { 1, 15, 14, 7, 6, 7, 7, 7, 7, 0, 0, 0, 0, 0, 0, 0, 0 },
		{ 0, 1, 13, 12, 5, 6, 6, 6, 5, 0, 0, 0, 0, 0, 0, 0, 0 },
		{ 0, 0, 1, 11, 10, 4, 5, 5, 4, 0, 0, 0, 0, 0, 0, 0, 0 },
		{ 0, 0, 0, 1, 1, 9, 8, 4, 4, 0, 0, 0, 0, 0, 0, 0, 0 } },
		//YUV444
		{ { 1, 5, 7, 7, 7, 7, 15, 11, 8, 15, 11, 15, 11, 15, 11, 7, 4 },
		{ 0, 1, 4, 6, 6, 6, 6, 14, 10, 14, 10, 14, 10, 1, 14, 10, 6 },
		{ 0, 0, 1, 5, 5, 5, 5, 5, 13, 9, 13, 9, 13, 9, 13, 9, 5 },
		{ 0, 0, 0, 3, 3, 4, 4, 4, 4, 4, 12, 12, 8, 12, 8, 12, 8 } }

	};
	int yuv = p_Vid->yuv_format - 1;

	// se->value1 : numcoeff
	// se->value2 : numtrailingones
	se->len = lentab[yuv][se->value2][se->value1];
	se->inf = codtab[yuv][se->value2][se->value1];

	if (se->len == 0)
	{
		printf("ERROR: (numcoeff,trailingones) not valid: (%d, %d)\n",
			se->value1, se->value2);
		exit(-1);
	}

	symbol2vlc(se);

	writeUVLC2buffer(se, dp->bitstream);

	if (se->type != SE_HEADER)
		dp->bitstream->write_flag = 1;

#if TRACE
	if (dp->bitstream->trace_enabled)
		trace2out(se);
#endif

	return (se->len);
}

/*!
************************************************************************
* \brief
*    write VLC for NumCoeff and TrailingOnes
************************************************************************
*/

int writeSyntaxElement_NumCoeffTrailingOnes(Enc_SyntaxElement *se, Enc_DataPartition *dp)
{
	static const byte lentab[3][4][17] =
	{
		{   // 0702
			{ 1, 6, 8, 9, 10, 11, 13, 13, 13, 14, 14, 15, 15, 16, 16, 16, 16 },
			{ 0, 2, 6, 8, 9, 10, 11, 13, 13, 14, 14, 15, 15, 15, 16, 16, 16 },
			{ 0, 0, 3, 7, 8, 9, 10, 11, 13, 13, 14, 14, 15, 15, 16, 16, 16 },
			{ 0, 0, 0, 5, 6, 7, 8, 9, 10, 11, 13, 14, 14, 15, 15, 16, 16 },
		},
		{
			{ 2, 6, 6, 7, 8, 8, 9, 11, 11, 12, 12, 12, 13, 13, 13, 14, 14 },
			{ 0, 2, 5, 6, 6, 7, 8, 9, 11, 11, 12, 12, 13, 13, 14, 14, 14 },
			{ 0, 0, 3, 6, 6, 7, 8, 9, 11, 11, 12, 12, 13, 13, 13, 14, 14 },
			{ 0, 0, 0, 4, 4, 5, 6, 6, 7, 9, 11, 11, 12, 13, 13, 13, 14 },
		},
		{
			{ 4, 6, 6, 6, 7, 7, 7, 7, 8, 8, 9, 9, 9, 10, 10, 10, 10 },
			{ 0, 4, 5, 5, 5, 5, 6, 6, 7, 8, 8, 9, 9, 9, 10, 10, 10 },
			{ 0, 0, 4, 5, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 10 },
			{ 0, 0, 0, 4, 4, 4, 4, 4, 5, 6, 7, 8, 8, 9, 10, 10, 10 },
		},

	};

	static const byte codtab[3][4][17] =
	{
		{
			{ 1, 5, 7, 7, 7, 7, 15, 11, 8, 15, 11, 15, 11, 15, 11, 7, 4 },
			{ 0, 1, 4, 6, 6, 6, 6, 14, 10, 14, 10, 14, 10, 1, 14, 10, 6 },
			{ 0, 0, 1, 5, 5, 5, 5, 5, 13, 9, 13, 9, 13, 9, 13, 9, 5 },
			{ 0, 0, 0, 3, 3, 4, 4, 4, 4, 4, 12, 12, 8, 12, 8, 12, 8 },
		},
		{
			{ 3, 11, 7, 7, 7, 4, 7, 15, 11, 15, 11, 8, 15, 11, 7, 9, 7 },
			{ 0, 2, 7, 10, 6, 6, 6, 6, 14, 10, 14, 10, 14, 10, 11, 8, 6 },
			{ 0, 0, 3, 9, 5, 5, 5, 5, 13, 9, 13, 9, 13, 9, 6, 10, 5 },
			{ 0, 0, 0, 5, 4, 6, 8, 4, 4, 4, 12, 8, 12, 12, 8, 1, 4 },
		},
		{
			{ 15, 15, 11, 8, 15, 11, 9, 8, 15, 11, 15, 11, 8, 13, 9, 5, 1 },
			{ 0, 14, 15, 12, 10, 8, 14, 10, 14, 14, 10, 14, 10, 7, 12, 8, 4 },
			{ 0, 0, 13, 14, 11, 9, 13, 9, 13, 10, 13, 9, 13, 9, 11, 7, 3 },
			{ 0, 0, 0, 12, 11, 10, 9, 8, 13, 12, 12, 12, 8, 12, 10, 6, 2 },
		},
	};
	int vlcnum = se->len;

	// se->value1 : numcoeff
	// se->value2 : numtrailingones

	if (vlcnum == 3)
	{
		se->len = 6;  // 4 + 2 bit FLC
		if (se->value1 > 0)
		{
			se->inf = ((se->value1 - 1) << 2) | se->value2;
		}
		else
		{
			se->inf = 3;
		}
	}
	else
	{
		se->len = lentab[vlcnum][se->value2][se->value1];
		se->inf = codtab[vlcnum][se->value2][se->value1];
	}

	if (se->len == 0)
	{
		printf("ERROR: (numcoeff,trailingones) not valid: vlc=%d (%d, %d)\n",
			vlcnum, se->value1, se->value2);
		exit(-1);
	}

	symbol2vlc(se);

	writeUVLC2buffer(se, dp->bitstream);

	if (se->type != SE_HEADER)
		dp->bitstream->write_flag = 1;

#if TRACE
	if (dp->bitstream->trace_enabled)
		trace2out(se);
#endif

	return (se->len);
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
void getNonAffNeighbourV(Enc_Macroblock *currMB, int xN, int yN, int mb_size[2], PixelPos *pix)
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
void get4x4NeighbourV(Enc_Macroblock *currMB, int block_x, int block_y, int mb_size[2], PixelPos *pix)
{
	//currMB->p_Vid->getNeighbourI(currMB, block_x, block_y, mb_size, pix);
	getNonAffNeighbourV(currMB, block_x, block_y, mb_size, pix);
	if (pix->available)
	{
		pix->x >>= 2;
		pix->y >>= 2;
		pix->pos_x >>= 2;
		pix->pos_y >>= 2;
	}
}

/*!
************************************************************************
* \brief
*    Get the Prediction from the Neighboring Blocks for Number of Nonzero Coefficients
*
*    Chroma Blocks
************************************************************************
*/
int predict_nnz_chromaI(Enc_Macroblock *currMB, int i, int j)
{
	Enc_Slice *currSlice = currMB->p_Slice;
	Enc_VideoParameters *p_Vid = currMB->p_Vid;
	PixelPos pix;

	int pred_nnz = 0;
	int cnt = 0;

	if (p_Vid->yuv_format != YUV444)
	{
		//YUV420 and YUV422
		// left block
		get4x4NeighbourV(currMB, ((i & 0x01) << 2) - 1, ((j - 4) << 2), p_Vid->mb_size[IS_CHROMA], &pix);

		if (is_intra(currMB) && pix.available && p_Vid->active_pps->constrained_intra_pred_flag && ((currSlice->partition_mode != 0) && !currSlice->idr_flag))
		{
			pix.available &= p_Vid->intra_block[pix.mb_addr];
			if (!pix.available)
				cnt++;
		}

		if (pix.available)
		{
			pred_nnz = p_Vid->nz_coeff[pix.mb_addr][2 * (i >> 1) + pix.x][4 + pix.y];
			cnt++;
		}

		// top block
		get4x4NeighbourV(currMB, ((i & 0x01) << 2), ((j - 4) << 2) - 1, p_Vid->mb_size[IS_CHROMA], &pix);

		if (is_intra(currMB) && pix.available && p_Vid->active_pps->constrained_intra_pred_flag && ((currSlice->partition_mode != 0) && !currSlice->idr_flag))
		{
			pix.available &= p_Vid->intra_block[pix.mb_addr];
			if (!pix.available)
				cnt++;
		}

		if (pix.available)
		{
			pred_nnz += p_Vid->nz_coeff[pix.mb_addr][2 * (i >> 1) + pix.x][4 + pix.y];
			cnt++;
		}
	}


	if (cnt == 2)
	{
		pred_nnz++;
		pred_nnz >>= 1;
	}

	return pred_nnz;
}

/*!
************************************************************************
* \brief
*    Get the Prediction from the Neighboring Blocks for Number of Nonzero Coefficients
*
*    Luma Blocks
************************************************************************
*/
int predict_nnzI(Enc_Macroblock *currMB, int block_type, int i, int j, VideoParameters *Vid)
{
	Enc_Slice *currSlice = currMB->p_Slice;
	Enc_VideoParameters *p_Vid = currMB->p_Vid;
	PixelPos pix;

	int pred_nnz = 0;
	int cnt = 0;

	// left block
	get4x4NeighbourV(currMB, (i << 2) - 1, (j << 2), p_Vid->mb_size[IS_LUMA], &pix);

	if (is_intra(currMB) && pix.available && p_Vid->active_pps->constrained_intra_pred_flag && ((currSlice->partition_mode != 0) && !currSlice->idr_flag))
	{
		pix.available &= p_Vid->intra_block[pix.mb_addr];
		if (!pix.available)
			cnt++;
	}
	if (pix.available)
	{
		if ((p_Vid->nz_coeff[pix.mb_addr][pix.x][pix.y] != Vid->nz_coeff[pix.mb_addr][0][pix.y][pix.x]))
			cnt = cnt;
		switch (block_type)
		{
		case LUMA:
			pred_nnz = p_Vid->nz_coeff[pix.mb_addr][pix.x][pix.y];
			cnt++;
			break;
		case CB:
			pred_nnz = p_Vid->nz_coeff[pix.mb_addr][pix.x][4 + pix.y];
			cnt++;
			break;
		case CR:
			pred_nnz = p_Vid->nz_coeff[pix.mb_addr][pix.x][8 + pix.y];
			cnt++;
			break;
		default:
			////error("writeCoeff4x4_CAVLC: Invalid block type", 600);
			break;
		}
	}

	// top block
	get4x4NeighbourV(currMB, (i << 2), (j << 2) - 1, p_Vid->mb_size[IS_LUMA], &pix);

	if (is_intra(currMB) && pix.available && p_Vid->active_pps->constrained_intra_pred_flag && ((currSlice->partition_mode != 0) && !currSlice->idr_flag))
	{
		pix.available &= p_Vid->intra_block[pix.mb_addr];
		if (!pix.available)
			cnt++;
	}

	if (pix.available)
	{
		if ((p_Vid->nz_coeff[pix.mb_addr][pix.x][pix.y] != Vid->nz_coeff[pix.mb_addr][0][pix.y][pix.x]))
			cnt = cnt;
		switch (block_type)
		{
		case LUMA:
			pred_nnz += p_Vid->nz_coeff[pix.mb_addr][pix.x][pix.y];
			cnt++;
			break;
		case CB:
			pred_nnz += p_Vid->nz_coeff[pix.mb_addr][pix.x][4 + pix.y];
			cnt++;
			break;
		case CR:
			pred_nnz += p_Vid->nz_coeff[pix.mb_addr][pix.x][8 + pix.y];
			cnt++;
			break;
		default:
			////error("writeCoeff4x4_CAVLC: Invalid block type", 600);
			break;
		}
	}

	if (cnt == 2)
	{
		pred_nnz++;
		pred_nnz >>= 1;
	}

	return pred_nnz;
}

/*!
************************************************************************
* \brief
*    Writes coeff of an 4x4 block (CAVLC)
*
* \author
*    Karl Lillevold <karll@real.com>
*    contributions by James Au <james@ubvideo.com>
************************************************************************
*/
int writeCoeff4x4_CAVLC_normal(Enc_Macroblock* currMB, int block_type, int iX, int iY, int param, int *lv, int *rn, VideoParameters *Vid)
{
	Enc_Slice* currSlice = currMB->p_Slice;
	Enc_VideoParameters *p_Vid = currSlice->p_Vid;
	int           no_bits = 0;
	Enc_SyntaxElement se;
	Enc_DataPartition *dataPart;
	//  const int *partMap   = assignSE2partitionI[currSlice->partition_mode];
	const int assignSE2partition_NoDP_I[SE_MAX_ELEMENTS] =
	{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
	const int *partMap = assignSE2partition_NoDP_I;

	int k, level = 1, run = 0, vlcnum;
	int numcoeff = 0, lastcoeff = 0, numtrailingones = 0;
	int numones = 0, totzeros = 0, zerosleft, numcoef;
	int numcoeff_vlc;
	int code, level_two_or_higher;
	int dptype = 0;
	int nnz, max_coeff_num = 0, cdc = 0, cac = 0;
	int subblock_x, subblock_y;
	int *mb_bits_coeff = &currMB->bits.mb_y_coeff;
#if TRACE
	char type[15];
#endif

	static const int incVlc[] = { 0, 3, 6, 12, 24, 48, 32768 };  // maximum vlc = 6


	int*  pLevel = NULL;
	int*  pRun = NULL;
	switch (block_type)
	{
	case LUMA:
		max_coeff_num = 16;

		pLevel = lv;//currSlice->cofAC[b8][b4][0];
		pRun = rn;//currSlice->cofAC[b8][b4][1];
#if TRACE
		sprintf(type, "%s", "Luma");
#endif
		dptype = (is_intra(currMB)) ? SE_LUM_AC_INTRA : SE_LUM_AC_INTER;
		break;

	case CHROMA_AC:
		max_coeff_num = 15;
		mb_bits_coeff = &currMB->bits.mb_uv_coeff;
		cac = 1;

		//pLevel = currSlice->cofAC[b8][b4][0];
		//pRun   = currSlice->cofAC[b8][b4][1];
#if TRACE
		sprintf(type, "%s", "ChrAC");
#endif
		dptype = (is_intra(currMB)) ? SE_CHR_AC_INTRA : SE_CHR_AC_INTER;
		break;

	case CHROMA_DC:
		max_coeff_num = p_Vid->num_cdc_coeff;
		mb_bits_coeff = &currMB->bits.mb_uv_coeff;
		cdc = 1;

		pLevel = currSlice->cofDC[param + 1][0];
		pRun = currSlice->cofDC[param + 1][1];
#if TRACE
		sprintf(type, "%s", "ChrDC");
#endif
		dptype = (is_intra(currMB)) ? SE_CHR_DC_INTRA : SE_CHR_DC_INTER;
		break;

	case LUMA_INTRA16x16AC:
		max_coeff_num = 15;

		pLevel = lv;//currSlice->cofAC[b8][b4][0];
		pRun = rn;//currSlice->cofAC[b8][b4][1];
#if TRACE
		sprintf(type, "%s", "Lum16AC");
#endif
		dptype = SE_LUM_AC_INTRA;
		break;

	case LUMA_INTRA16x16DC:
		max_coeff_num = 16;

		pLevel = lv;//currSlice->cofDC[0][0];
		pRun = rn;//currSlice->cofDC[0][1];
#if TRACE
		sprintf(type, "%s", "Lum16DC");
#endif
		dptype = SE_LUM_DC_INTRA;
		break;


	default:
		////error("writeCoeff4x4_CAVLC: Invalid block type", 600);
		break;
	}

	dataPart = &(currSlice->partArr[partMap[dptype]]);

	for (k = 0; (k < ((cdc) ? p_Vid->num_cdc_coeff : 16)) && level != 0; k++)
	{
		level = pLevel[k]; // level
		run = pRun[k];   // run

		if (level)
		{

			totzeros += run;
			if (iabs(level) == 1)
			{
				numones++;
				numtrailingones++;
				numtrailingones = imin(numtrailingones, 3); // clip to 3
			}
			else
			{
				numtrailingones = 0;
			}
			numcoeff++;
			lastcoeff = k;
		}
	}
	if (!cdc)
	{
		if (!cac)
		{
			// luma
			subblock_x = iX;//((b8 & 0x1) == 0) ? (((b4 & 0x1) == 0) ? 0 : 1) : (((b4 & 0x1) == 0) ? 2 : 3);
							// horiz. position for coeff_count context
			subblock_y = iY;//(b8 < 2) ? ((b4 < 2) ? 0 : 1) : ((b4 < 2) ? 2 : 3);
							// vert.  position for coeff_count context
			nnz = predict_nnzI(currMB, LUMA, subblock_x, subblock_y, Vid);
		}
		else
		{
			// chroma AC
			subblock_x = param >> 4;
			subblock_y = param & 15;
			nnz = predict_nnz_chromaI(currMB, subblock_x, subblock_y);
		}
		p_Vid->nz_coeff[currMB->mbAddrX][subblock_x][subblock_y] = numcoeff;
		numcoeff_vlc = (nnz < 2) ? 0 : ((nnz < 4) ? 1 : ((nnz < 8) ? 2 : 3));
	}
	else
	{
		// chroma DC (has its own VLC)
		// numcoeff_vlc not relevant
		numcoeff_vlc = 0;

		subblock_x = param;
		subblock_y = param;
	}

	se.type = dptype;

	se.value1 = numcoeff;
	se.value2 = numtrailingones;
	se.len = numcoeff_vlc; /* use len to pass vlcnum */

#if TRACE
	snprintf(se.tracestring,
		TRACESTRING_SIZE, "%s # c & tr.1s(%d,%d) vlc=%d #c=%d #t1=%d",
		type, subblock_x, subblock_y, numcoeff_vlc, numcoeff, numtrailingones);
#endif

	if (!cdc)
		writeSyntaxElement_NumCoeffTrailingOnes(&se, dataPart);
	else
		writeSyntaxElement_NumCoeffTrailingOnesChromaDC(p_Vid, &se, dataPart);

	*mb_bits_coeff += se.len;
	no_bits += se.len;

	if (!numcoeff)
		return no_bits;

	if (numcoeff)
	{
		code = 0;
		for (k = lastcoeff; k > lastcoeff - numtrailingones; k--)
		{
			level = pLevel[k]; // level
#ifdef  _DEBUG
			if (iabs(level) > 1)
			{
				printf("ERROR: level > 1\n");
				exit(-1);
			}
#endif
			code <<= 1;

			code |= (level < 0);
		}

		if (numtrailingones)
		{
			se.type = dptype;

			se.value2 = numtrailingones;
			se.value1 = code;

#if TRACE
			snprintf(se.tracestring,
				TRACESTRING_SIZE, "%s trailing ones sign (%d,%d)",
				type, subblock_x, subblock_y);
#endif

			writeSyntaxElement_VLC(&se, dataPart);
			*mb_bits_coeff += se.len;
			no_bits += se.len;

		}

		// encode levels
		level_two_or_higher = (numcoeff > 3 && numtrailingones == 3) ? 0 : 1;

		vlcnum = (numcoeff > 10 && numtrailingones < 3) ? 1 : 0;

		for (k = lastcoeff - numtrailingones; k >= 0; k--)
		{
			level = pLevel[k]; // level

			se.value1 = level;
			se.type = dptype;

#if TRACE
			snprintf(se.tracestring,
				TRACESTRING_SIZE, "%s lev (%d,%d) k=%d vlc=%d lev=%3d",
				type, subblock_x, subblock_y, k, vlcnum, level);
#endif

			if (level_two_or_higher)
			{
				level_two_or_higher = 0;

				if (se.value1 > 0)
					se.value1--;
				else
					se.value1++;
			}

			//    encode level

			if (vlcnum == 0)
				writeSyntaxElement_Level_VLC1(&se, dataPart, 100);//p_Vid->active_sps->profile_idc);
			else
				writeSyntaxElement_Level_VLCN(&se, vlcnum, dataPart, 100);//p_Vid->active_sps->profile_idc);

																		  // update VLC table
			if (iabs(level) > incVlc[vlcnum])
				vlcnum++;

			if ((k == lastcoeff - numtrailingones) && iabs(level) > 3)
				vlcnum = 2;

			*mb_bits_coeff += se.len;
			no_bits += se.len;
		}

		// encode total zeroes
		if (numcoeff < max_coeff_num)
		{

			se.type = dptype;
			se.value1 = totzeros;

			vlcnum = numcoeff - 1;

			se.len = vlcnum;

#if TRACE
			snprintf(se.tracestring,
				TRACESTRING_SIZE, "%s totalrun (%d,%d) vlc=%d totzeros=%3d",
				type, subblock_x, subblock_y, vlcnum, totzeros);
#endif
			if (!cdc)
				writeSyntaxElement_TotalZeros(&se, dataPart);
			else
				writeSyntaxElement_TotalZerosChromaDC(p_Vid, &se, dataPart);

			*mb_bits_coeff += se.len;
			no_bits += se.len;
		}

		// encode run before each coefficient
		zerosleft = totzeros;
		numcoef = numcoeff;
		for (k = lastcoeff; k >= 0; k--)
		{
			run = pRun[k]; // run

			se.value1 = run;
			se.type = dptype;

			// for last coeff, run is remaining totzeros
			// when zerosleft is zero, remaining coeffs have 0 run
			if ((!zerosleft) || (numcoeff <= 1))
				break;

			if (numcoef > 1 && zerosleft)
			{
				vlcnum = imin(zerosleft - 1, RUNBEFORE_NUM_M1);
				se.len = vlcnum;

#if TRACE
				snprintf(se.tracestring,
					TRACESTRING_SIZE, "%s run (%d,%d) k=%d vlc=%d run=%2d",
					type, subblock_x, subblock_y, k, vlcnum, run);
#endif

				writeSyntaxElement_Run(&se, dataPart);

				*mb_bits_coeff += se.len;
				no_bits += se.len;

				zerosleft -= run;
				numcoef--;
			}
		}
	}
	return no_bits;
}

/*!
************************************************************************
* \brief
*    Resets the nz_coeff parameters for a macroblock
************************************************************************/
void reset_mb_nz_coeff(Enc_VideoParameters *p_Vid, int mb_number)
{
	memset(&p_Vid->nz_coeff[mb_number][0][0], 0, BLOCK_SIZE * (BLOCK_SIZE + p_Vid->num_blk8x8_uv) * sizeof(int));
}

void init_Data(char bufB, int off, Macroblock *currMB, Slice *currSlice, VideoParameters *p_Vid, Enc_Macroblock *curr_MBI, Enc_Slice *currSliceI, Enc_VideoParameters *p_VidI)
{
	currSliceI->partition_mode = 0;
	p_VidI->active_pps->constrained_intra_pred_flag = p_Vid->active_pps->constrained_intra_pred_flag;
	p_VidI->mb_size[0][0] = p_Vid->mb_size[0][0];
	p_VidI->mb_size[0][1] = p_Vid->mb_size[0][1];
	p_VidI->mb_size[1][0] = p_Vid->mb_size[1][0];
	p_VidI->mb_size[1][1] = p_Vid->mb_size[1][1];
	p_VidI->mb_size[2][0] = p_Vid->mb_size[2][0];
	p_VidI->mb_size[2][1] = p_Vid->mb_size[2][1];
	p_VidI->yuv_format = (ColorFormat)p_Vid->yuv_format;	//YUV420
	p_VidI->intra_block = (short *)p_Vid->intra_block;
	curr_MBI->mbAddrA = currMB->mbAddrA;
	curr_MBI->mbAddrB = currMB->mbAddrB;
	curr_MBI->mbAddrC = currMB->mbAddrC;
	curr_MBI->mbAddrD = currMB->mbAddrD;
	curr_MBI->mbAddrX = currMB->mbAddrX;
	curr_MBI->mbAvailA = currMB->mbAvailA;
	curr_MBI->mbAvailB = currMB->mbAvailB;
	curr_MBI->mbAvailC = currMB->mbAvailC;
	curr_MBI->mbAvailD = currMB->mbAvailD;
	curr_MBI->DeblockCall = currMB->DeblockCall;
	curr_MBI->mb_type = currMB->mb_type;
	currSliceI->idr_flag = currSlice->idr_flag;

	Enc_VStream->streamBuffer = outputBuffer.BuffFrame;//currSlice->partArr->bitstream->streamBuffer;//point

	Enc_VStream->byte_pos = outputBuffer.pBuffFrame;//currSlice->partArr->bitstream->frame_bitoffset >> 3;
	Enc_VStream->bits_to_go = off;// - (currSlice->partArr->bitstream->frame_bitoffset & 7);
	Enc_VStream->byte_buf = bufB;//currSlice->partArr->bitstream->streamBuffer[currStreamI->byte_pos] >> currStreamI->bits_to_go;

	curr_MBI->p_Slice = currSliceI;
	curr_MBI->p_Vid = p_VidI;
	curr_MBI->p_Vid->PicPos->x = currMB->p_Vid->PicPos->x;
	curr_MBI->p_Vid->PicPos->y = currMB->p_Vid->PicPos->y;
	currSliceI->p_Vid = p_VidI;
	Enc_dp->bitstream = Enc_VStream;
	currSliceI->partArr = Enc_dp;
}
int ReadPLNZV(int numcoeff, int *Run)
{
	int i;
	int PLNZ = numcoeff;
	int *Rn = Run;
	for (i = 0; i <= numcoeff; i++)
	{
		PLNZ += *Rn++;
	}
	return PLNZ;
}

void WriteFrame(int startBit, int endBit, Slice *currSlice)
{
	while (endBit != startBit)
	{
		if (outputBuffer.pBuffFrame == 0x110e)
			endBit = endBit;
		WriteBit(ReadBit(currSlice));
		startBit++;
	}
}
void WriteBit(int Bit)
{
	outputBuffer.byteBuffFrame |= Bit << (outputBuffer.bits2Go - 1);
	outputBuffer.bits2Go--;
	if (outputBuffer.bits2Go == 0)
	{
		outputBuffer.bits2Go = 8;
		outputBuffer.BuffFrame[outputBuffer.pBuffFrame++] = outputBuffer.byteBuffFrame;
		outputBuffer.byteBuffFrame = 0;
	}
}
int ReadBit(Slice *currSlice)
{
	unsigned char tempByte;
	tempByte = currSlice->partArr->bitstream->streamBuffer[writepoint >> 3];
	return (tempByte >> (7 - (writepoint++ & 7))) & 0x1;
}
int Embeding(int *point, char *StringI, int *lv, int *rn, int numcoeff, int PLNZ)
{
	int i, j;
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
			{
				if (numcoeff != 16)
				{
					lv[numcoeff - 1] = 0;
					for (i = 15; i >= 0; i--)
					{
						if (rn[i] != 0)
						{
							rn[i]--;
							break;
						}
					}
					j = i;
					for (i = (numcoeff - 1); i >= j; i--)
					{
						lv[i + 1] = lv[i];
					}
					lv[j] = 1;
				}
				else
					lv[15] = (lv[15] < 0) ? (lv[15] * -1) : lv[15];
			}
		}
		else						//Pi + 1
			rn[numcoeff - 1] += 1;
	}
	else if (((PLNZ & 0x1) != 0) && (Si == 0))				//odd		//Pi + 1
	{
		ret = 1;
		rn[numcoeff - 1] += 1;
	}
	else if ((numcoeff == 16) && (Si == 0))
	{
		lv[15] = (lv[15] > 0) ? (lv[15] * -1) : lv[15];
	}
	(*point)++;
	return ret;
}
