
/*!
 ***********************************************************************
 *  \file
 *     configfile.h
 *  \brief
 *     Prototypes for configfile.c and definitions of used structures.
 ***********************************************************************
 */

#ifndef _CONFIGFILE_H_
#define _CONFIGFILE_H_


#define DEFAULTCONFIGFILENAME "encoder.cfg"

#define PROFILE_IDC     88
#define LEVEL_IDC       21


typedef struct {
  char *TokenName;
  void *Place;
  int Type;
} Mapping;



InputParameters configinput;


#ifdef INCLUDED_BY_CONFIGFILE_C

Mapping Map[] = {
    {"ProfileIDC",               &configinput.ProfileIDC,              0},
    {"LevelIDC",                 &configinput.LevelIDC,                0},
    {"FrameRate",                &configinput.FrameRate,                    2},
    {"IDRIntraEnable",           &configinput.idr_enable,              0},
    {"StartFrame",               &configinput.start_frame,             0},
    {"IntraPeriod",              &configinput.intra_period,            0},
    {"FramesToBeEncoded",        &configinput.no_frames,               0},
    {"QPISlice",                 &configinput.qp0,                          0},
    {"QPPSlice",                 &configinput.qpN,                          0},
    {"FrameSkip",                &configinput.jumpd,                   0},
    {"UseHadamard",              &configinput.hadamard,                0},
    {"SearchRange",              &configinput.search_range,            0},
    {"NumberReferenceFrames",    &configinput.num_ref_frames,               0},
    {"PList0References",         &configinput.P_List0_refs,            0},
    {"BList0References",         &configinput.B_List0_refs,            0},
    {"BList1References",         &configinput.B_List1_refs,            0},
    {"Log2MaxFrameNum",          &configinput.Log2MaxFrameNum,         0},
    {"SourceWidth",              &configinput.img_width,               0},
    {"SourceHeight",             &configinput.img_height,              0},
    {"MbLineIntraUpdate",        &configinput.intra_upd,               0},
    {"SliceMode",                &configinput.slice_mode,              0},
    {"SliceArgument",            &configinput.slice_argument,          0},
    {"UseConstrainedIntraPred",  &configinput.UseConstrainedIntraPred, 0},
    {"InputFile",                &configinput.infile,                  1},
    {"InputHeaderLength",        &configinput.infile_header,           0},
    {"OutputFile",               &configinput.outfile,                 1},
    {"ReconFile",                &configinput.ReconFile,               1},
    {"TraceFile",                &configinput.TraceFile,               1},
    {"NumberBFrames",            &configinput.successive_Bframe,       0},
    {"QPBSlice",                 &configinput.qpB,                          0},
    {"BStoredQPOffset",          &configinput.qpBSoffset,                   0},         
    {"DirectModeType",           &configinput.direct_spatial_mv_pred_flag,  0},
    {"DirectInferenceFlag",      &configinput.directInferenceFlag,     0},
    {"SPPicturePeriodicity",     &configinput.sp_periodicity,          0},
    {"QPSPSlice",                &configinput.qpsp,                         0},
    {"QPSP2Slice",               &configinput.qpsp_pred,                    0},
    {"SymbolMode",               &configinput.symbol_mode,             0},
    {"OutFileMode",              &configinput.of_mode,                 0},
    {"PartitionMode",            &configinput.partition_mode,          0},
    {"PictureTypeSequence",      &configinput.PictureTypeSequence,     1},
    {"InterSearch16x16",         &configinput.InterSearch16x16,        0},
    {"InterSearch16x8",          &configinput.InterSearch16x8 ,        0},
    {"InterSearch8x16",          &configinput.InterSearch8x16,         0},
    {"InterSearch8x8",           &configinput.InterSearch8x8 ,         0},
    {"InterSearch8x4",           &configinput.InterSearch8x4,          0},
    {"InterSearch4x8",           &configinput.InterSearch4x8,          0},
    {"InterSearch4x4",           &configinput.InterSearch4x4,          0},
    {"IntraDisableInterOnly",    &configinput.IntraDisableInterOnly,   0},
    {"Intra4x4ParDisable",       &configinput.Intra4x4ParDisable,      0},
    {"Intra4x4DiagDisable",      &configinput.Intra4x4DiagDisable,     0},
    {"Intra4x4DirDisable",       &configinput.Intra4x4DirDisable,      0},
    {"Intra16x16ParDisable",     &configinput.Intra16x16ParDisable,    0},
    {"Intra16x16PlaneDisable",   &configinput.Intra16x16PlaneDisable,  0},
    {"ChromaIntraDisable",       &configinput.ChromaIntraDisable,      0},

#ifdef _FULL_SEARCH_RANGE_
    {"RestrictSearchRange",      &configinput.full_search,             0},
#endif
#ifdef _ADAPT_LAST_GROUP_
    {"LastFrameNumber",          &configinput.last_frame,              0},
#endif
#ifdef _CHANGE_QP_
    {"ChangeQPI",                &configinput.qp02,                    0},
    {"ChangeQPP",                &configinput.qpN2,                    0},
    {"ChangeQPB",                &configinput.qpB2,                    0},
    {"ChangeQPBs",               &configinput.qpBs2offset,             0},
    {"ChangeQPStart",            &configinput.qp2start,                0},
#endif
    {"RDOptimization",           &configinput.rdopt,                   0},
    {"DisableThresholding",      &configinput.disthres,                0},
    {"DisableBSkipRDO",          &configinput.nobskip,                 0},
    {"LossRateA",                &configinput.LossRateA,               0},
    {"LossRateB",                &configinput.LossRateB,               0},
    {"LossRateC",                &configinput.LossRateC,               0},
    {"NumberOfDecoders",         &configinput.NoOfDecoders,            0},
    {"RestrictRefFrames",        &configinput.RestrictRef ,            0},
#ifdef _LEAKYBUCKET_
    {"NumberofLeakyBuckets",     &configinput.NumberLeakyBuckets,      0},
    {"LeakyBucketRateFile",      &configinput.LeakyBucketRateFile,     1},
    {"LeakyBucketParamFile",     &configinput.LeakyBucketParamFile,    1},
#endif
    {"PicInterlace",             &configinput.PicInterlace,            0},
    {"MbInterlace",              &configinput.MbInterlace,             0},

    {"IntraBottom",              &configinput.IntraBottom,             0},

    {"NumberFramesInEnhancementLayerSubSequence", &configinput.NumFramesInELSubSeq, 0},
    {"NumberOfFrameInSecondIGOP",&configinput.NumFrameIn2ndIGOP,       0},
    {"RandomIntraMBRefresh",     &configinput.RandomIntraMBRefresh,    0},
		
		
    {"WeightedPrediction",       &configinput.WeightedPrediction,      0},
    {"WeightedBiprediction",     &configinput.WeightedBiprediction,    0},
    {"UseWeightedReferenceME",   &configinput.UseWeightedReferenceME,  0},
    {"StoredBPictures",          &configinput.StoredBPictures,         0},
    {"PyramidCoding",            &configinput.PyramidCoding,           0},
    {"ExplicitPyramidFormat",    &configinput.ExplicitPyramidFormat,   1},
    {"PyramidRefReorder",        &configinput.PyramidRefReorder,       0},
    {"PocMemoryManagement",      &configinput.PocMemoryManagement,     0},

    {"LoopFilterParametersFlag", &configinput.LFSendParameters,        0},
    {"LoopFilterDisable",        &configinput.LFDisableIdc,            0},
    {"LoopFilterAlphaC0Offset",  &configinput.LFAlphaC0Offset,         0},
    {"LoopFilterBetaOffset",     &configinput.LFBetaOffset,            0},
    {"SparePictureOption",       &configinput.SparePictureOption,      0},
    {"SparePictureDetectionThr", &configinput.SPDetectionThreshold,    0},
    {"SparePicturePercentageThr",&configinput.SPPercentageThreshold,   0},

    {"num_slice_groups_minus1",           &configinput.num_slice_groups_minus1,           0},
    {"slice_group_map_type",              &configinput.slice_group_map_type,              0},
    {"slice_group_change_direction_flag", &configinput.slice_group_change_direction_flag, 0},
    {"slice_group_change_rate_minus1",    &configinput.slice_group_change_rate_minus1,    0},
    {"SliceGroupConfigFileName",          &configinput.SliceGroupConfigFileName,          1},
		

    {"UseRedundantSlice",        &configinput.redundant_slice_flag,    0},
    {"PicOrderCntType",          &configinput.pic_order_cnt_type,      0},

    {"ContextInitMethod",        &configinput.context_init_method,     0},
    {"FixedModelNumber",         &configinput.model_number,            0},

    {"Transform8x8Mode",         &configinput.AllowTransform8x8,       0},
    {"ReportFrameStats",         &configinput.ReportFrameStats,        0},
    // Rate Control
    {"RateControlEnable",        &configinput.RCEnable,                0},
    {"Bitrate",                  &configinput.bit_rate,                0},
    {"InitialQP",                &configinput.SeinitialQP,             0},
    {"BasicUnit",                &configinput.basicunit,               0},
    {"ChannelType",              &configinput.channel_type,            0},

    // Q_Matrix
    {"QmatrixFile",              &configinput.QmatrixFile,                1},
    {"ScalingMatrixPresentFlag", &configinput.ScalingMatrixPresentFlag,   0},
    {"ScalingListPresentFlag0",  &configinput.ScalingListPresentFlag[0],  0},
    {"ScalingListPresentFlag1",  &configinput.ScalingListPresentFlag[1],  0},
    {"ScalingListPresentFlag2",  &configinput.ScalingListPresentFlag[2],  0},
    {"ScalingListPresentFlag3",  &configinput.ScalingListPresentFlag[3],  0},
    {"ScalingListPresentFlag4",  &configinput.ScalingListPresentFlag[4],  0},
    {"ScalingListPresentFlag5",  &configinput.ScalingListPresentFlag[5],  0},
    {"ScalingListPresentFlag6",  &configinput.ScalingListPresentFlag[6],  0},
    {"ScalingListPresentFlag7",  &configinput.ScalingListPresentFlag[7],  0},

    // Fast ME enable
    {"UseFME",                   &configinput.FMEnable,                0},

    {"ChromaQPOffset",           &configinput.chroma_qp_index_offset,  0},    

    // Fidelity Range Extensions
    {"BitDepthLuma",             &configinput.BitDepthLuma,            0},
    {"BitDepthChroma",           &configinput.BitDepthChroma,          0},
    {"YUVFormat",                &configinput.yuv_format,              0},
    {"RGBInput",                 &configinput.rgb_input_flag,          0},
    {"CbQPOffset",               &configinput.cb_qp_index_offset,      0},
    {"CrQPOffset",               &configinput.cr_qp_index_offset,      0},
   
    // Lossless Coding
    {"QPPrimeYZeroTransformBypassFlag", &configinput.lossless_qpprime_y_zero_flag, 0},

    // Residue Color Transform
    {"ResidueTransformFlag",     &configinput.residue_transform_flag , 0},

    {NULL,                       NULL,                                -1}
};

#endif

#ifndef INCLUDED_BY_CONFIGFILE_C
extern Mapping Map[];
#endif


void Configure (int ac, char *av[]);
void PatchInputNoFrames();

#endif

