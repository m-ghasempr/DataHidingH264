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
 * \file rtp.h
 *
 * \brief
 *    Prototypes for rtp.c
 *************************************************************************************
 */

#ifndef _RTP_H_
#define _RTP_H_

#define RTP_MAX_PARAMETER_SET 10    /*!< Maximum number of supported parameter sets */

#define RTP_PARAMETER_SET_OK  0
#define RTP_PARAMETER_SET_OUT_OF_RANGE (-1)
#define RTP_PARAMETER_SET_INVALID (-2)

#define RTP_MAX_STRING_LEN 4096

#include "mbuffer.h"

typedef struct
{
  int Valid;            //!< 0: invalid, else valid, not yet used
  int UseMultpred;
  int MaxPicID;                 
  int BufCycle;
  int PixAspectRatioX;
  int PixAspectRatioY;
  int DisplayWindowOffsetTop;
  int DisplayWindowOffsetBottom;
  int DisplayWindowOffsetRight;
  int DisplayWindowOffsetLeft;
  int XSizeMB;
  int YSizeMB;
  int EntropyCoding;
  int MotionResolution;
  int PartitioningType;
  int IntraPredictionType;
  int HRCParameters;
  int *MBAmap;
} ParameterSet_t;

typedef struct
{
  int FramesToBeEncoded;
  int FrameSkip;
  char SequenceFileName [RTP_MAX_STRING_LEN];
  int NumberBFrames;
} InformationParameterSet_t;

#define MAXRTPPAYLOADLEN  (65536 - 40)    //!< Maximum payload size of an RTP packet */
#define MAXRTPPACKETSIZE  (65536 - 28)    //!< Maximum size of an RTP packet incl. header */
#define H26LPAYLOADTYPE 105               //!< RTP paylaod type fixed here for simplicity*/
#define H26LSSRC 0x12345678               //!< SSRC, chosen to simplify debugging */
#define RTP_TR_TIMESTAMP_MULT 1000        //!< should be something like 27 Mhz / 29.97 Hz */

typedef struct 
{
  unsigned int v;          //!< Version, 2 bits, MUST be 0x2
  unsigned int p;          //!< Padding bit, Padding MUST NOT be used
  unsigned int x;          //!< Extension, MUST be zero
  unsigned int cc;         /*!< CSRC count, normally 0 in the absence 
                                of RTP mixers */
  unsigned int m;          //!< Marker bit
  unsigned int pt;         //!< 7 bits, Payload Type, dynamically established 
  unsigned int seq;        /*!< RTP sequence number, incremented by one for
                                each sent packet */
  unsigned int old_seq;    //!< to detect wether packets were lost
  unsigned int timestamp;  //!< timestamp, 27 MHz for H.26L
  unsigned int ssrc;       //!< Synchronization Source, chosen randomly
  byte *       payload;    //!< the payload including payload headers
  unsigned int paylen;     //!< length of payload in bytes
  byte *       packet;     //!< complete packet including header and payload
  unsigned int packlen;    //!< length of packet, typically paylen+12
} RTPpacket_t;


typedef struct
{
  int ParameterSet;
  int structure;
  int PictureID;
  int PictureNum;
  int SliceType;
  int FirstMBInSliceX;
  int FirstMBInSliceY;
  int Direct_type;
  int InitialQP;
  int SwitchSP;  
  int InitialSPQP;  
  int SliceID;            //!< not used for single Slice packets, see VCEG-N72
  int RPBT;
  RMPNIbuffer_t *RMPNIbuffer;
  MMCObuffer_t  *MMCObuffer;
  int disposable_flag;                          //!< flag for disposable frame, 1:disposable
  int num_ref_pic_active_fwd;                   //!< number of forward reference
  int num_ref_pic_active_bwd;                   //!< number of backward reference
  int explicit_B_prediction;                    //!< type of weight for bi-pred, 0:average 1:implicit
} RTPSliceHeader_t;


extern ParameterSet_t ParSet[];

int  ReadRTPPacket (struct img_par *img, struct inp_par *inp, FILE *bits);
int  readSyntaxElement_RTP(SyntaxElement *sym, struct img_par *img, struct inp_par *inp, struct datapartition *dP);
int  RTP_startcode_follows(struct img_par *img, struct inp_par *inp);
void RTP_get_symbol(struct img_par *img, struct inp_par *inp, struct datapartition *dP, SyntaxElement *sym, Bitstream *currStream);
int  RTP_symbols_available (Bitstream *currStream);
int  RTPInterpretParameterSetPacket (char *buf, int buflen);
void  RTPUseParameterSet (int n, struct img_par *img, struct inp_par *inp);
int  RTPReadPartitions (struct img_par *img, struct inp_par *inp, FILE *bits);
int  DecomposeRTPpacket (RTPpacket_t *p);
void DumpRTPHeader (RTPpacket_t *p);
int  RTPInterpretSliceHeader (byte *buf, int bufsize, int ReadSliceId, RTPSliceHeader_t *sh);
int  RTPInterpretPartitionHeader (byte *buf, int bufsize, RTPSliceHeader_t *sh);
int  readSliceRTP (struct img_par *img, struct inp_par *inp);
int  RTPSequenceHeader (struct img_par *img, struct inp_par *inp, FILE *bits);
void RTPSetImgInp(struct img_par *img, struct inp_par *inp, RTPSliceHeader_t *sh);
int  RTPGetFollowingSliceHeader (struct img_par *img, RTPpacket_t *p, RTPSliceHeader_t *sh);
int  get_lastMB(struct img_par *img, RTPSliceHeader_t *sh, RTPpacket_t *p);
int RTPReadPacket (RTPpacket_t *p, FILE *bits);
void RTPProcessDataPartitionedSlice (struct img_par *img, struct inp_par *inp, FILE *bits, 
                                     RTPpacket_t *a, int a_SliceID);
void CopyPartitionBitstring (struct img_par *img, RTPpacket_t *p, Bitstream *b, int dP, struct inp_par *inp);

#endif
