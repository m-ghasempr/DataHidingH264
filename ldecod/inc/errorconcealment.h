

/*!
 ****************************************************************************
 * \file errorconcealment.h
 *
 * \brief
 *    Header file for errorconcealment.c
 *
 ****************************************************************************
 */

#ifndef _ERRORCONCEALMENT_H_
#define _ERRORCONCEALMENT_H_

extern int set_ec_flag(ImageParameters *p_Img, int se);
extern void reset_ec_flags(ImageParameters *p_Img);
extern int get_concealed_element(ImageParameters *p_Img, SyntaxElement *sym);

#endif

