H.264 Data Embedder
===============

This software embeds data into H.264 coded videos. The embedded data can be extracted by <a target="_blank" href="https://github.com/mohghasem/DataExtractorH264">Data Extractor</a>. It uses quantized DCT coefficients to embed data and can embeds multiple data into the video. For more information, please read the paper that is cited below.

### Paper presenting video data hiding:
M. Ghasempour and M. Ghanbari, "A Low Complexity System for Multiple Data Embedding Into H.264 Coded Video Bit-Stream," in IEEE Transactions on Circuits and Systems for Video Technology, vol. 30, no. 11, pp. 4009-4019, Nov. 2020, [DOI: 10.1109/TCSVT.2019.2947545](https://doi.org/10.1109/TCSVT.2019.2947545)

- [Build](#build)
- [Usage](#usage)
- [Copyright](#copyright)
- [References](#references)

Build
----------------------------------------------

1. You will need to check [JM repository](https://vcgit.hhi.fraunhofer.de/jvet/JM) to see how to build the project.

Usage
----------------------------------------------
```
./ldecod (or ldecod.exe on Windows) [-d JM_CONFIG_FILE] [-MD INPUT_TEXT_FILE] [-TF FRAME_TYPE] [-TH THRESHOLD_VALUE]
```
- **JM_CONFIG_FILE**: The original config file for H.264/AVC JM reference software
- **INPUT_TEXT_FILE**: The path and name of the data file that should be embedded into the video.
- **FRAME_TYPE**: The frame type that is used for embedding data. The value can be from 0 to 4. 0: **B-frames**, 1: **I-frames**, 2: **P-frames**, 3: **Limited P-frames (only the last P-frame of each GOP to prevent drift distortion)** 4: **Only for multiple data insertion**
- **THRESHOLD_VALUE**: The threshold value specifies which block is selected for embedding data. The blocks, which the position of last non-zero coefficient is higher than the threshold value, are embedding candidate blocks. It can be 1 to 14. The higher threshold value, the lower embedding capacity can be achieved, but the quality of video is preserved.

Example:
```
./ldecod -d decoder.cfg -MD input.txt -TF 0 -TH 7
```
The example command embeds "input.txt" file into the macroblocks of B-frames, which the position of last non-zero coefficient is higher than 7.

Copyright
----------------------------------------------
Permission is hereby granted, without written agreement and without license or royalty fees, to use, copy, modify, and distribute the software provided and its documentation for research purpose only, provided that this copyright notice and the original authors' names appear on all copies and supporting documentation. The software provided may not be commercially distributed.

References
----------------------------------------------
- T. Wiegand, G. J. Sullivan, G. Bjontegaard and A. Luthra, “Overview of the H.264/AVC video coding standard,” IEEE Transactions on Circuits and Systems for Video Technology, 13(7), 560–576, 2003.
- M. Fallahpour, S. Shirmohammadi, M. Ghanbari, “A high capacity data hiding algorithm for H.264/AVC video,” Security and Communication Networks, 8(16), 2947-2955, 2015.
