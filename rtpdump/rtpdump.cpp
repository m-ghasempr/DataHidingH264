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

// rtpdump.cpp : Defines the entry point for the console application.
//

#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <malloc.h>

#include "stdafx.h"

int main(int argc, char* argv[])

{
  unsigned int bufsize, pacno=0;;
  unsigned char buf[65000];
  int i, intime;
  FILE *f;


  if (argc != 2)
  {
    printf ("Usage: %s <H.26L RTP packet file\n");
    return -1;
  }

  if (NULL == (f = fopen (argv[1], "rb")))
  {
    printf ("%s: cannot open H.26L packet file %s for reading\n", argv[0], argv[1]);
    return -2;
  }

  while (1)
  {
    if (4 != fread (&bufsize, 1, 4, f))
      return 0;
    if (4 != fread (&intime, 1, 4, f))
    {
      printf ("Panic, cannot read timestamp, old software version file?\n");
      return -1;
    }
    printf ("\n\npacket #%4d containing %5d bytes\n", pacno++, bufsize);
    if (bufsize != fread (buf, 1, bufsize, f))
    {
      printf ("Problems while reading buffer, exit\n");
      return -3;
    }

    for (i=0; i< 25; i++)
      printf ("%02x ", buf[i]);
    printf ("\n");
    printf ("Version (V): %d\n", buf[0] & 0x3);
    printf ("Padding (P): %d\n", (buf[0] & 0x4) >> 2);
    printf ("Extension (X): %d\n", (buf[0] & 0x8) >> 3);
    printf ("CSRC count (CC): %d\n", (buf[0] & 0xf0) >> 4);
    printf ("Marker bit (M): %d\n", buf[1] & 0x1);
    printf ("Payload Type (PT): %d\n", (buf[1] & 0xfe) >> 1);
    printf ("Sequence Number: %d\n", buf[2] | (buf[3] << 8));
    printf ("Timestamp: %d\n", buf[4] | (buf[5] <<8) | (buf[6] <<16) | (buf[7] <<24));
    printf ("SSRC: %d\n", buf[8] | (buf[9] <<8) | (buf[10] <<16) | (buf[11] <<24));

    printf ("First Byte: 0x%x\n", buf[12]);
  }
}

