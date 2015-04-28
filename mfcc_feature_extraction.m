function mfcc_features = mfcc_feature_extraction(signal, fs)
  mfcc_features = heq(mfcc(log_mel_spectrogram(signal, fs)));
end
