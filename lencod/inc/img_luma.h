/*!
 ***************************************************************************
 * \file
 *    img_luma.h
 *
 * \author
 *    Athanasios Leontaris           <aleon@dolby.com>
 *    Alexis Michael Tourapis        <alexis.tourapis@dolby.com>
 *
 * \date
 *    4. October 2006
 *
 * \brief
 *    Headerfile for luma interpolation functions
 **************************************************************************
 */

#ifndef _IMG_LUMA_H_
#define _IMG_LUMA_H_

extern void getSubImagesLuma       ( ImageParameters *p_Img, InputParameters *p_Inp, StorablePicture *s );
extern void getSubImageInteger     ( StorablePicture *s, imgpel **dstImg, imgpel **srcImg);
extern void getHorSubImageSixTap   ( ImageParameters *p_Img, StorablePicture *s, imgpel **dst_imgY, imgpel **ref_imgY);
extern void getVerSubImageSixTap   ( ImageParameters *p_Img, StorablePicture *s, imgpel **dst_imgY, imgpel **ref_imgY);
extern void getVerSubImageSixTapTmp( ImageParameters *p_Img, StorablePicture *s, imgpel **dst_imgY);
extern void getSubImageBiLinear    ( StorablePicture *s, imgpel **dstImg, imgpel **srcImgL, imgpel **srcImgR);
extern void getHorSubImageBiLinear ( StorablePicture *s, imgpel **dstImg, imgpel **srcImgL, imgpel **srcImgR);
extern void getVerSubImageBiLinear ( StorablePicture *s, imgpel **dstImg, imgpel **srcImgT, imgpel **srcImgB);
extern void getDiagSubImageBiLinear( StorablePicture *s, imgpel **dstImg, imgpel **srcImgT, imgpel **srcImgB);
#endif // _IMG_LUMA_H_
