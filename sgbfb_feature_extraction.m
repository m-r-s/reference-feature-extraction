function sgbfb_features = sgbfb_feature_extraction(signal, fs)
  sgbfb_features = heq(sgbfb(log_mel_spectrogram(signal, fs)));
end
