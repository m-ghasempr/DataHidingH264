/*
***********************************************************************
* COPYRIGHT AND WARRANTY INFORMATION
*
* Copyright 2001, International Telecommunications Union, Geneva
*
* DISCLAIMER OF WARRANTY
*
* These software programs are available to the user without any
* license fee or royalty on an "as is" basis. The ITU disclaims
* any and all warranties, whether express, implied, or
* statutory, including any implied warranties of merchantability
* or of fitness for a particular purpose.  In no event shall the
* contributor or the ITU be liable for any incidental, punitive, or
* consequential damages of any kind whatsoever arising from the
* use of these programs.
*
* This disclaimer of warranty extends to the user of these programs
* and user's customers, employees, agents, transferees, successors,
* and assigns.
*
* The ITU does not represent or warrant that the programs furnished
* hereunder are free of infringement of any third-party patents.
* Commercial implementations of ITU-T Recommendations, including
* shareware, may be subject to royalty fees to patent holders.
* Information regarding the ITU-T patent policy is available from
* the ITU Web site at http://www.itu.int.
*
* THIS IS NOT A GRANT OF PATENT RIGHTS - SEE THE ITU-T PATENT POLICY.
************************************************************************
*/

/*!
 *************************************************************************************
 * \file vlc.h
 *
 * \brief
 *    Prototypes for VLC coding funtions
 * \author
 *     Karsten Suehring
 *************************************************************************************
 */

#ifndef _VLC_H_
#define _VLC_H_

int u_1 (char *tracestring, int value, DataPartition *part);
int se_v (char *tracestring, int value, DataPartition *part);
int ue_v (char *tracestring, int value, DataPartition *part);
int u_v (int n, char *tracestring, int value, DataPartition *part);


void levrun_linfo_c2x2(int level,int run,int *len,int *info);
void levrun_linfo_intra(int level,int run,int *len,int *info);
void levrun_linfo_inter(int level,int run,int *len,int *info);

int   writeSyntaxElement_UVLC(SyntaxElement *se, DataPartition *this_dataPart);
int   writeSyntaxElement_fixed(SyntaxElement *se, DataPartition *this_dataPart);
int   writeSyntaxElement2Buf_UVLC(SyntaxElement *se, Bitstream* this_streamBuffer );
void  writeUVLC2buffer(SyntaxElement *se, Bitstream *currStream);
int   writeSyntaxElement2Buf_Fixed(SyntaxElement *se, Bitstream* this_streamBuffer );
int   symbol2uvlc(SyntaxElement *se);
void  ue_linfo(int n, int dummy, int *len,int *info);
void  se_linfo(int mvd, int dummy, int *len,int *info);
void  cbp_linfo_intra(int cbp, int dummy, int *len,int *info);
void  cbp_linfo_inter(int cbp, int dummy, int *len,int *info);

// CAVLC
void  CAVLC_init();
int writeCoeff4x4_CAVLC (int block_type, int b8, int b4, int param);

int   writeSyntaxElement_VLC(SyntaxElement *se, DataPartition *this_dataPart);
int   writeSyntaxElement_TotalZeros(SyntaxElement *se, DataPartition *this_dataPart);
int   writeSyntaxElement_TotalZerosChromaDC(SyntaxElement *se, DataPartition *this_dataPart);
int   writeSyntaxElement_Run(SyntaxElement *se, DataPartition *this_dataPart);
int   writeSyntaxElement_NumCoeffTrailingOnes(SyntaxElement *se, DataPartition *this_dataPart);
int   writeSyntaxElement_NumCoeffTrailingOnesChromaDC(SyntaxElement *se, DataPartition *this_dataPart);
int   writeSyntaxElement_Level_VLC1(SyntaxElement *se, DataPartition *this_dataPart);
int   writeSyntaxElement_Level_VLCN(SyntaxElement *se, int vlc, DataPartition *this_dataPart);
int   writeSyntaxElement_Intra4x4PredictionMode(SyntaxElement *se, DataPartition *this_dataPart);

#endif

