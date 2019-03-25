/*!
 ***************************************************************************
 * \file transform.c
 *
 * \brief
 *    Transform functions
 *
 * \author
 *    Main contributors (see contributors.h for copyright, address and affiliation details)
 *    - Alexis Michael Tourapis
 * \date
 *    01. July 2007
 **************************************************************************
 */

#include "global.h"
#include "transform.h"

static int tmp[64];

void forward4x4(int (*block)[16], int (*tblock)[16], int pos_y, int pos_x)
{
  int i, ii;  
  int *pTmp = tmp, *pblock;
  static int p0,p1,p2,p3;
  static int t0,t1,t2,t3;

  // Horizontal
  for (i=pos_y; i < pos_y + BLOCK_SIZE; i++)
  {
    pblock = &block[i][pos_x];
    p0 = *(pblock++);
    p1 = *(pblock++);
    p2 = *(pblock++);
    p3 = *(pblock  );

    t0 = p0 + p3;
    t1 = p1 + p2;
    t2 = p1 - p2;
    t3 = p0 - p3;

    *(pTmp++) =  t0 + t1;
    *(pTmp++) = (t3 << 1) + t2;
    *(pTmp++) =  t0 - t1;    
    *(pTmp++) =  t3 - (t2 << 1);
  }

  // Vertical 
  for (i=0; i < BLOCK_SIZE; i++)
  {
    pTmp = tmp + i;
    p0 = *pTmp;
    p1 = *(pTmp += BLOCK_SIZE);
    p2 = *(pTmp += BLOCK_SIZE);
    p3 = *(pTmp += BLOCK_SIZE);

    t0 = p0 + p3;
    t1 = p1 + p2;
    t2 = p1 - p2;
    t3 = p0 - p3;

    ii = pos_x + i;
    tblock[pos_y    ][ii] = t0 +  t1;
    tblock[pos_y + 1][ii] = t2 + (t3 << 1);
    tblock[pos_y + 2][ii] = t0 -  t1;
    tblock[pos_y + 3][ii] = t3 - (t2 << 1);
  }
}


void inverse4x4(int (*tblock)[16], int (*block)[16], int pos_y, int pos_x)
{
  int i, ii;  
  int *pTmp = tmp, *pblock;
  static int p0,p1,p2,p3;
  static int t0,t1,t2,t3;

  // Horizontal
  for (i = pos_y; i < pos_y + BLOCK_SIZE; i++)
  {
    pblock = &tblock[i][pos_x];
    t0 = *(pblock++);
    t1 = *(pblock++);
    t2 = *(pblock++);
    t3 = *(pblock  );

    p0 =  t0 + t2;
    p1 =  t0 - t2;
    p2 = (t1 >> 1) - t3;
    p3 =  t1 + (t3 >> 1);

    *(pTmp++) = p0 + p3;
    *(pTmp++) = p1 + p2;
    *(pTmp++) = p1 - p2;
    *(pTmp++) = p0 - p3;
  }

  //  Vertical 
  for (i = 0; i < BLOCK_SIZE; i++)
  {
    pTmp = tmp + i;
    t0 = *pTmp;
    t1 = *(pTmp += BLOCK_SIZE);
    t2 = *(pTmp += BLOCK_SIZE);
    t3 = *(pTmp += BLOCK_SIZE);

    p0 = t0 + t2;
    p1 = t0 - t2;
    p2 =(t1 >> 1) - t3;
    p3 = t1 + (t3 >> 1);

    ii = i + pos_x;
    block[pos_y    ][ii] = p0 + p3;
    block[pos_y + 1][ii] = p1 + p2;
    block[pos_y + 2][ii] = p1 - p2;
    block[pos_y + 3][ii] = p0 - p3;
  }
}


void hadamard4x4(int (*block)[4], int (*tblock)[4])
{
  int i;
  int *pTmp = tmp, *pblock;
  static int p0,p1,p2,p3;
  static int t0,t1,t2,t3;

  // Horizontal
  for (i = 0; i < BLOCK_SIZE; i++)
  {
    pblock = block[i];
    p0 = *(pblock++);
    p1 = *(pblock++);
    p2 = *(pblock++);
    p3 = *(pblock  );

    t0 = p0 + p3;
    t1 = p1 + p2;
    t2 = p1 - p2;
    t3 = p0 - p3;

    *(pTmp++) = t0 + t1;
    *(pTmp++) = t3 + t2;
    *(pTmp++) = t0 - t1;    
    *(pTmp++) = t3 - t2;
  }

  // Vertical 
  for (i = 0; i < BLOCK_SIZE; i++)
  {
    pTmp = tmp + i;
    p0 = *pTmp;
    p1 = *(pTmp += BLOCK_SIZE);
    p2 = *(pTmp += BLOCK_SIZE);
    p3 = *(pTmp += BLOCK_SIZE);

    t0 = p0 + p3;
    t1 = p1 + p2;
    t2 = p1 - p2;
    t3 = p0 - p3;

    tblock[0][i] = (t0 + t1) >> 1;
    tblock[1][i] = (t2 + t3) >> 1;
    tblock[2][i] = (t0 - t1) >> 1;
    tblock[3][i] = (t3 - t2) >> 1;
  }
}


void ihadamard4x4(int (*tblock)[4], int (*block)[4])
{
  int i;  
  int *pTmp = tmp, *pblock;
  static int p0,p1,p2,p3;
  static int t0,t1,t2,t3;

  // Horizontal
  for (i = 0; i < BLOCK_SIZE; i++)
  {
    pblock = tblock[i];
    t0 = *(pblock++);
    t1 = *(pblock++);
    t2 = *(pblock++);
    t3 = *(pblock  );

    p0 = t0 + t2;
    p1 = t0 - t2;
    p2 = t1 - t3;
    p3 = t1 + t3;

    *(pTmp++) = p0 + p3;
    *(pTmp++) = p1 + p2;
    *(pTmp++) = p1 - p2;
    *(pTmp++) = p0 - p3;
  }

  //  Vertical 
  for (i = 0; i < BLOCK_SIZE; i++)
  {
    pTmp = tmp + i;
    t0 = *pTmp;
    t1 = *(pTmp += BLOCK_SIZE);
    t2 = *(pTmp += BLOCK_SIZE);
    t3 = *(pTmp += BLOCK_SIZE);

    p0 = t0 + t2;
    p1 = t0 - t2;
    p2 = t1 - t3;
    p3 = t1 + t3;
    
    block[0][i] = p0 + p3;
    block[1][i] = p1 + p2;
    block[2][i] = p1 - p2;
    block[3][i] = p0 - p3;
  }
}


void forward8x8(int (*block)[16], int (*tblock)[16], int pos_y, int pos_x)
{
  int i, ii;  
  int *pTmp = tmp, *pblock;
  static int a0, a1, a2, a3;
  static int p0, p1, p2, p3, p4, p5 ,p6, p7;
  static int b0, b1, b2, b3, b4, b5, b6, b7;

  // Horizontal
  for (i=pos_y; i < pos_y + BLOCK_SIZE_8x8; i++)
  {
    pblock = &block[i][pos_x];
    p0 = *(pblock++);
    p1 = *(pblock++);
    p2 = *(pblock++);
    p3 = *(pblock++);
    p4 = *(pblock++);
    p5 = *(pblock++);
    p6 = *(pblock++);
    p7 = *(pblock  );

    a0 = p0 + p7;
    a1 = p1 + p6;
    a2 = p2 + p5;
    a3 = p3 + p4;

    b0 = a0 + a3;
    b1 = a1 + a2;
    b2 = a0 - a3;
    b3 = a1 - a2;

    a0 = p0 - p7;
    a1 = p1 - p6;
    a2 = p2 - p5;
    a3 = p3 - p4;

    b4 = a1 + a2 + ((a0 >> 1) + a0);
    b5 = a0 - a3 - ((a2 >> 1) + a2);
    b6 = a0 + a3 - ((a1 >> 1) + a1);
    b7 = a1 - a2 + ((a3 >> 1) + a3);

    *(pTmp++) =  b0 + b1;
    *(pTmp++) =  b4 + (b7 >> 2);
    *(pTmp++) =  b2 + (b3 >> 1);
    *(pTmp++) =  b5 + (b6 >> 2);
    *(pTmp++) =  b0 - b1;
    *(pTmp++) =  b6 - (b5 >> 2);
    *(pTmp++) = (b2 >> 1) - b3;                 
    *(pTmp++) = (b4 >> 2) - b7;
  }

  // Vertical 
  for (i=0; i < BLOCK_SIZE_8x8; i++)
  {
    pTmp = tmp + i;
    p0 = *pTmp;
    p1 = *(pTmp += BLOCK_SIZE_8x8);
    p2 = *(pTmp += BLOCK_SIZE_8x8);
    p3 = *(pTmp += BLOCK_SIZE_8x8);
    p4 = *(pTmp += BLOCK_SIZE_8x8);
    p5 = *(pTmp += BLOCK_SIZE_8x8);
    p6 = *(pTmp += BLOCK_SIZE_8x8);
    p7 = *(pTmp += BLOCK_SIZE_8x8);

    a0 = p0 + p7;
    a1 = p1 + p6;
    a2 = p2 + p5;
    a3 = p3 + p4;

    b0 = a0 + a3;
    b1 = a1 + a2;
    b2 = a0 - a3;
    b3 = a1 - a2;

    a0 = p0 - p7;
    a1 = p1 - p6;
    a2 = p2 - p5;
    a3 = p3 - p4;

    b4 = a1 + a2 + ((a0 >> 1) + a0);
    b5 = a0 - a3 - ((a2 >> 1) + a2);
    b6 = a0 + a3 - ((a1 >> 1) + a1);
    b7 = a1 - a2 + ((a3 >> 1) + a3);

    ii = pos_x + i;
    tblock[pos_y    ][ii] =  b0 + b1;
    tblock[pos_y + 1][ii] =  b4 + (b7 >> 2);
    tblock[pos_y + 2][ii] =  b2 + (b3 >> 1);
    tblock[pos_y + 3][ii] =  b5 + (b6 >> 2);
    tblock[pos_y + 4][ii] =  b0 - b1;
    tblock[pos_y + 5][ii] =  b6 - (b5 >> 2);
    tblock[pos_y + 6][ii] = (b2 >> 1) - b3;
    tblock[pos_y + 7][ii] = (b4 >> 2) - b7;
  }
}

void inverse8x8(int (*tblock)[16], int (*block)[16], int pos_y, int pos_x)
{
  int i, ii;
  int *pTmp = tmp, *pblock;
  static int a0, a1, a2, a3;
  static int p0, p1, p2, p3, p4, p5 ,p6, p7;  
  static int b0, b1, b2, b3, b4, b5, b6, b7;

  // Horizontal  
  for (i=pos_y; i < pos_y + BLOCK_SIZE_8x8; i++)
  {
    pblock = &tblock[i][pos_x];
    p0 = *(pblock++);
    p1 = *(pblock++);
    p2 = *(pblock++);
    p3 = *(pblock++);
    p4 = *(pblock++);
    p5 = *(pblock++);
    p6 = *(pblock++);
    p7 = *(pblock  );

    a0 = p0 + p4;
    a1 = p0 - p4;
    a2 = p6 - (p2 >> 1);
    a3 = p2 + (p6 >> 1);

    b0 =  a0 + a3;
    b2 =  a1 - a2;
    b4 =  a1 + a2;
    b6 =  a0 - a3;

    a0 = -p3 + p5 - p7 - (p7 >> 1);    
    a1 =  p1 + p7 - p3 - (p3 >> 1);    
    a2 = -p1 + p7 + p5 + (p5 >> 1);    
    a3 =  p3 + p5 + p1 + (p1 >> 1);

    
    b1 =  a0 + (a3>>2);    
    b3 =  a1 + (a2>>2);    
    b5 =  a2 - (a1>>2);
    b7 =  a3 - (a0>>2);                

    *(pTmp++) = b0 + b7;
    *(pTmp++) = b2 - b5;
    *(pTmp++) = b4 + b3;
    *(pTmp++) = b6 + b1;
    *(pTmp++) = b6 - b1;
    *(pTmp++) = b4 - b3;
    *(pTmp++) = b2 + b5;
    *(pTmp++) = b0 - b7;
  }

  //  Vertical 
  for (i=0; i < BLOCK_SIZE_8x8; i++)
  {
    pTmp = tmp + i;
    p0 = *pTmp;
    p1 = *(pTmp += BLOCK_SIZE_8x8);
    p2 = *(pTmp += BLOCK_SIZE_8x8);
    p3 = *(pTmp += BLOCK_SIZE_8x8);
    p4 = *(pTmp += BLOCK_SIZE_8x8);
    p5 = *(pTmp += BLOCK_SIZE_8x8);
    p6 = *(pTmp += BLOCK_SIZE_8x8);
    p7 = *(pTmp += BLOCK_SIZE_8x8);

    a0 =  p0 + p4;
    a1 =  p0 - p4;
    a2 =  p6 - (p2>>1);
    a3 =  p2 + (p6>>1);

    b0 = a0 + a3;
    b2 = a1 - a2;
    b4 = a1 + a2;
    b6 = a0 - a3;

    a0 = -p3 + p5 - p7 - (p7 >> 1);
    a1 =  p1 + p7 - p3 - (p3 >> 1);
    a2 = -p1 + p7 + p5 + (p5 >> 1);
    a3 =  p3 + p5 + p1 + (p1 >> 1);


    b1 =  a0 + (a3 >> 2);
    b7 =  a3 - (a0 >> 2);
    b3 =  a1 + (a2 >> 2);
    b5 =  a2 - (a1 >> 2);

    ii = i + pos_x;
    block[pos_y    ][ii] = b0 + b7;
    block[pos_y + 1][ii] = b2 - b5;
    block[pos_y + 2][ii] = b4 + b3;
    block[pos_y + 3][ii] = b6 + b1;
    block[pos_y + 4][ii] = b6 - b1;
    block[pos_y + 5][ii] = b4 - b3;
    block[pos_y + 6][ii] = b2 + b5;
    block[pos_y + 7][ii] = b0 - b7;
  }
}

