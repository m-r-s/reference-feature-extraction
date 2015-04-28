function gbfb_features = gbfb_feature_extraction(signal, fs)
  gbfb_features = heq(gbfb(log_mel_spectrogram(signal, fs)));
end
