
/*!
 ***************************************************************************
 *
 * \file transform.h
 *
 * \brief
*    prototypes of transform functions
  *
 * \date
 *    10 July 2007
 *
 * \author
 *    Main contributors (see contributors.h for copyright, address and affiliation details)
 *    Alexis Michael Tourapis
 **************************************************************************/

#ifndef _TRANSFORM_H_
#define _TRANSFORM_H_

void forward4x4   (int **block , int **tblock, int pos_y, int pos_x);
void inverse4x4   (int **tblock, int **block, int pos_y, int pos_x);
void forward8x8   (int **block , int **tblock, int pos_y, int pos_x);
void inverse8x8   (int **tblock, int **block, int pos_y, int pos_x);
void hadamard4x4  (int **block , int **tblock);
void ihadamard4x4 (int **tblock, int **block);
void hadamard4x2  (int **block , int **tblock);
void ihadamard4x2 (int **tblock, int **block);
void hadamard2x2  (int **block , int tblock[4]);
void ihadamard2x2 (int block[4], int tblock[4]);

#endif //_TRANSFORM_H_
