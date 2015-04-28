# Reference Matlab/Octave implementations of feature extraction algorithms

The scripts provided in this software package were written to perform the
feature extraction in automatic speech recogniton experiments and to
evaluate the obtained recognition performance in [1].

This includes the calculation of:
- Logarithmically scaled Mel-spectrograms (LogMS) in *log_mel_spectrogram.m*
- Mel-frequency cepstral coefficient (MFCC) features in *mfcc_feature_extraction.m*
- Gabor filter bank features (GBFB) features in *gbfb_feature_extraction.m*
- Separable Gabor filter bank (SGBFB) features in *sgbfb_feature_extraction.m*
- Histogram equalization (HEQ) in *heq.m*
- Equal-performance SNR increase (EPSI) in *epsi.m*

Detailed explanations of the corresponding concepts are provided in [1].
The SGBFB feature extraction is closely related to the GBFB feature extraction
which was introduced in [2].
For an overview of current publications and further developments of the SGBFB
front-end, visit http://medi.uni-oldenburg.de/SGBFB

*play_demo.m* is demonstrates the use of the feature extraction scripts.

## References

[1] M.R. Schädler and B. Kollmeier,
    "Separable spectro-temporal Gabor filter bank features: Reducing the
    complexity of robust features for automatic speech recognition",
    J. Acoust. Soc. Am. Volume 137, Issue 4, pp. 2047-2059,
    DOI: 10.1121/1.4916618,
    URL: http://link.aip.org/link/?JAS/137/2047 (2015)

[2] M.R. Schädler, B.T. Meyer, B. Kollmeier
    "Spectro-temporal modulation subspace-spanning filter bank features for
    robust automatic speech recognition",
    J. Acoust. Soc. Am. Volume 131, Issue 5, pp. 4134-4151,
    DOI: 10.1121/1.3699200,
    URL: http://link.aip.org/link/?JAS/131/4134 (2012)
