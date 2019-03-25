
/*!
 ***************************************************************************
 *
 * \file filehandle.h
 *
 * \brief
 *
 * \date
 *
 * \author
 *    Main contributors (see contributors.h for copyright, address and affiliation details)
 *
 **************************************************************************/

#ifndef _FILEHANDLE_H_
#define _FILEHANDLE_H_

extern int  rewrite_paramsets (ImageParameters *p_Img, InputParameters *p_Inp);
extern int  start_sequence    (ImageParameters *p_Img, InputParameters *p_Inp);
extern int  terminate_sequence(ImageParameters *p_Img, InputParameters *p_Inp);
extern int  write_PPS         (ImageParameters *p_Img, InputParameters *p_Inp, int len, int PPS_id);


#endif
