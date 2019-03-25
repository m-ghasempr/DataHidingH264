/*!
************************************************************************
* \file io_tiff.h
*
* \brief
*    I/O functions related to tiff images
*    Part of code was based on libtiff (see http://www.libtiff.org/)
*
* \author
*     - Alexis Michael Tourapis         <alexismt@ieee.org>
*
************************************************************************
*/

#ifndef _IO_TIFF_H_
#define _IO_TIFF_H_
// See TIFF 6.0 Specification
// http://partners.adobe.com/public/developer/tiff/index.html

// TIFF header.
typedef struct tiff_header
{
  unsigned char  border[2];   //!< Byte order used (2 bytes). Valid values "II" (little endian) and "MM" (big endian). 
  uint16 version;     //!< TIFF identification number, i.e. 42 (2 bytes)
  unsigned int   IFDoff;      //!< The offset (in bytes) of the first Image File Directory/IFD (4 bytes).
} TIFFHeader;


// TIFF IFD structure
typedef struct tiff_idf_entry
{
  uint16 tIFD_tag;     //!< Identification tag
  uint16 tIFD_type;    //!< field type.
  unsigned int   tIFD_count;   //!< The number of values, Count of the indicated Type.
  unsigned int   tIFD_offset;  //!< The Value Offset, the file offset (in bytes) of the Value for the field.
} TIFFIFDEntry;

// TIFF Types
typedef enum {
  TIFF_BYTE       = 1,  // 8-bit unsigned integer.
  TIFF_ASCII      = 2,  // 8-bit byte that contains a 7-bit ASCII code.
  TIFF_SHORT      = 3,  // 16-bit (2-byte) unsigned integer.
  TIFF_LONG       = 4,  // 32-bit (4-byte) unsigned integer.
  TIFF_RATIONAL   = 5,  // Two LONGs: the first represents the numerator of
                        // a fraction; the second, the denominator
  TIFF_SBYTE      = 6,  // An 8-bit signed (twos-complement) integer.
  TIFF_UNDEFINED  = 7,  // An 8-bit undefined byte.
  TIFF_SSHORT     = 8,  // A 16-bit (2-byte) signed (twos-complement) integer.
  TIFF_SLONG      = 9,  // A 32-bit (4-byte) signed (twos-complement) integer.
  TIFF_SRATIONAL  = 10, // Two SLONG’s: the first represents the numerator of 
                        // a fraction, the second the denominator.
  TIFF_FLOAT      = 11, // Single precision (4-byte) IEEE format
  TIFF_DOUBLE     = 12  // Double precision (8-byte) IEEE format  
} TIFFType;
// TIFF Types

typedef enum {
  TIFFTAG_SUBFILETYPE            = 254, /* subfile data descriptor */
  TIFFTAG_OSUBFILETYPE           = 255,
  TIFFTAG_IMAGEWIDTH             = 256, /* image width in pixels */
  TIFFTAG_IMAGELENGTH            = 257, /* image height in pixels */
  TIFFTAG_BITSPERSAMPLE          = 258, /* bits per channel (sample) */
  TIFFTAG_COMPRESSION            = 259, /* data compression technique */
  TIFFTAG_PHOTOMETRIC            = 262, /* photometric interpretation */
  TIFFTAG_THRESHHOLDING          = 263, /* +thresholding used on data */
  TIFFTAG_CELLWIDTH              = 264, /* +dithering matrix width */
  TIFFTAG_CELLLENGTH             = 265, /* +dithering matrix height */
  TIFFTAG_FILLORDER              = 266, /* data order within a byte */
  TIFFTAG_DOCUMENTNAME           = 269, /* name of doc. image is from */
  TIFFTAG_IMAGEDESCRIPTION       = 270, /* info about image */
  TIFFTAG_MAKE                   = 271, /* scanner manufacturer name */
  TIFFTAG_MODEL                  = 272, /* scanner model name/number */
  TIFFTAG_STRIPOFFSETS           = 273, /* offsets to data strips */
  TIFFTAG_ORIENTATION            = 274, /* +image orientation */
  TIFFTAG_SAMPLESPERPIXEL        = 277, /* samples per pixel */
  TIFFTAG_ROWSPERSTRIP           = 278, /* rows per strip of data */
  TIFFTAG_STRIPBYTECOUNTS        = 279, /* bytes counts for strips */
  TIFFTAG_MINSAMPLEVALUE         = 280, /* +minimum sample value */
  TIFFTAG_MAXSAMPLEVALUE         = 281, /* +maximum sample value */
  TIFFTAG_XRESOLUTION            = 282, /* pixels/resolution in x */
  TIFFTAG_YRESOLUTION            = 283, /* pixels/resolution in y */
  TIFFTAG_PLANARCONFIG           = 284, /* storage organization */
  TIFFTAG_PAGENAME               = 285, /* page name image is from */
  TIFFTAG_XPOSITION              = 286, /* x page offset of image lhs */
  TIFFTAG_YPOSITION              = 287, /* y page offset of image lhs */
  TIFFTAG_FREEOFFSETS            = 288, /* +byte offset to free block */
  TIFFTAG_FREEBYTECOUNTS         = 289, /* +sizes of free blocks */
  TIFFTAG_GRAYRESPONSEUNIT       = 290, /* $gray scale curve accuracy */
  TIFFTAG_GRAYRESPONSECURVE      = 291, /* $gray scale response curve */
  TIFFTAG_GROUP3OPTIONS          = 292, /* 32 flag bits */
  TIFFTAG_T4OPTIONS              = 292, /* TIFF 6.0 proper name alias */
  TIFFTAG_GROUP4OPTIONS          = 293, /* 32 flag bits */
  TIFFTAG_T6OPTIONS              = 293    , /* TIFF 6.0 proper name */
  TIFFTAG_RESOLUTIONUNIT         = 296, /* units of resolutions */
  TIFFTAG_PAGENUMBER             = 297, /* page numbers of multi-page */
  TIFFTAG_COLORRESPONSEUNIT      = 300, /* $color curve accuracy */
  TIFFTAG_TRANSFERFUNCTION       = 301, /* !colorimetry info */
  TIFFTAG_SOFTWARE               = 305, /* name & release */
  TIFFTAG_DATETIME               = 306, /* creation date and time */
  TIFFTAG_ARTIST                 = 315, /* creator of image */
  TIFFTAG_HOSTCOMPUTER           = 316, /* machine where created */
  TIFFTAG_PREDICTOR              = 317, /* prediction scheme w/ LZW */
  TIFFTAG_WHITEPOINT             = 318, /* image white point */
  TIFFTAG_PRIMARYCHROMATICITIES  = 319, /* !primary chromaticities */
  TIFFTAG_COLORMAP               = 320, /* RGB map for pallette image */
  TIFFTAG_HALFTONEHINTS          = 321, /* !highlight+shadow info */
  TIFFTAG_TILEWIDTH              = 322, /* !tile width in pixels */
  TIFFTAG_TILELENGTH             = 323, /* !tile height in pixels */
  TIFFTAG_TILEOFFSETS            = 324, /* !offsets to data tiles */
  TIFFTAG_TILEBYTECOUNTS         = 325, /* !byte counts for tiles */
  TIFFTAG_BADFAXLINES            = 326, /* lines w/ wrong pixel count */
  TIFFTAG_CLEANFAXDATA           = 327, /* regenerated line info */
  TIFFTAG_CONSECUTIVEBADFAXLINES = 328, /* max consecutive bad lines */
  TIFFTAG_SUBIFD                 = 330, /* subimage descriptors */
  TIFFTAG_INKSET                 = 332, /* !inks in separated image */
  TIFFTAG_INKNAMES               = 333, /* !ascii names of inks */
  TIFFTAG_NUMBEROFINKS           = 334, /* !number of inks */
  TIFFTAG_DOTRANGE               = 336, /* !0% and 100% dot codes */
  TIFFTAG_TARGETPRINTER          = 337, /* !separation target */
  TIFFTAG_EXTRASAMPLES           = 338, /* !info about extra samples */
  TIFFTAG_SAMPLEFORMAT           = 339, /* !data sample format */
  TIFFTAG_SMINSAMPLEVALUE        = 340, /* !variable MinSampleValue */
  TIFFTAG_SMAXSAMPLEVALUE        = 341  /* !variable MaxSampleValue */
} TIFFTag;

typedef enum {
  COMPRESSION_NONE      = 1, /* dump mode (only mode that will be initially supported) */
  COMPRESSION_CCITTRLE  = 2, /* CCITT modified Huffman RLE */
  COMPRESSION_CCITTFAX3 = 3, /* CCITT Group 3 fax encoding */
  COMPRESSION_CCITT_T4  = 3, /* CCITT T.4 (TIFF 6 name) */
  COMPRESSION_CCITTFAX4 = 4, /* CCITT Group 4 fax encoding */
  COMPRESSION_CCITT_T6  = 4, /* CCITT T.6 (TIFF 6 name) */
  COMPRESSION_LZW       = 5, /* Lempel-Ziv  & Welch */
  COMPRESSION_OJPEG     = 6, /* !6.0 JPEG */
  COMPRESSION_JPEG      = 7  /* %JPEG DCT compression */
} TIFFCompression;


extern int ReadTIFFImage (InputParameters *p_Inp, VideoDataFile *input_file, int FrameNoInFile, FrameFormat *source, unsigned char *buf);

#endif

