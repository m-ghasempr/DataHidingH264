
/*!
 ***************************************************************************
 *
 * \file rtp.h
 *
 * \brief
 *    Definition of structures and functions to handle RTP headers.  For a
 *    description of RTP see RFC1889 on http://www.ietf.org
 *
 * \date
 *    30 September 2001
 *
 * \author
 *    Stephan Wenger   stewe@cs.tu-berlin.de
 **************************************************************************/

#ifndef _RTP_H_
#define _RTP_H_

#include "nalu.h"

#define MAXRTPPAYLOADLEN  (65536 - 40)    //!< Maximum payload size of an RTP packet
#define MAXRTPPACKETSIZE  (65536 - 28)    //!< Maximum size of an RTP packet incl. header
#define H264PAYLOADTYPE 105               //!< RTP paylaod type fixed here for simplicity
#define H264SSRC 0x12345678               //!< SSRC, chosen to simplify debugging
#define RTP_TR_TIMESTAMP_MULT 1000        //!< should be something like 27 Mhz / 29.97 Hz

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
  unsigned int timestamp;  //!< timestamp, 27 MHz for H.264
  unsigned int ssrc;       //!< Synchronization Source, chosen randomly
  byte *       payload;    //!< the payload including payload headers
  unsigned int paylen;     //!< length of payload in bytes
  byte *       packet;     //!< complete packet including header and payload
  unsigned int packlen;    //!< length of packet, typically paylen+12
} RTPpacket_t;

#if 0
extern int  ComposeRTPPacket (RTPpacket_t *p);
extern int  DecomposeRTPpacket (RTPpacket_t *p);
extern int  WriteRTPPacket (RTPpacket_t *p, FILE *f);
extern void DumpRTPHeader (RTPpacket_t *p);
extern int  RTPWriteBits (int Marker, int PacketType, void * bitstream,
                   int BitStreamLenInByte, FILE *out);

extern Boolean isAggregationPacket(VideoParameters *p_Vid);
extern int aggregationRTPWriteBits (int Marker, int PacketType, int subPacketType, void * bitstream, int BitStreamLenInByte, FILE *out);

extern void begin_sub_sequence_rtp(VideoParameters *p_Vid, InputParameters *p_Inp);
extern void end_sub_sequence_rtp  (VideoParameters *p_Vid, InputParameters *p_Inp);
#endif

extern void RTPUpdateTimestamp (VideoParameters *p_Vid, int tr);
extern void OpenRTPFile        (VideoParameters *p_Vid, char *Filename);
extern void CloseRTPFile       (VideoParameters *p_Vid);
extern int WriteRTPNALU        (VideoParameters *p_Vid, NALU_t *n);



#endif

