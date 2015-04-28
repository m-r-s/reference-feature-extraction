close all;
clear;
clc;

% Sample rate
fs = 16000;

% Lets generate a 1000 Hz windowed tone signal
sinewave = sin(2*pi*linspace(0,2000/2,fs/2)) .* hanning(fs/2).';

% Tone in Noise :)
signal = 0.1 .* randn(1,fs);
signal(1+fs/4:fs/4+fs/2) = signal(1+fs/4:fs/4+fs/2) + sinewave;

% Extract the feature vectors
mfcc_features = mfcc_feature_extraction(signal,fs);
gbfb_features = gbfb_feature_extraction(signal,fs);
sgbfb_features = sgbfb_feature_extraction(signal,fs);

% Get the corresponding log Mel-spectrogram
logms = log_mel_spectrogram(signal,fs);

% Plot the signal representations
figure;
subplot(4,1,1);
imagesc(logms);
axis xy;
ylabel('LogMS');

subplot(4,1,2);
imagesc(mfcc_features);
axis xy;
ylabel('MFCC');

subplot(4,1,3);
imagesc(gbfb_features);
axis xy;
ylabel('GBFB');

subplot(4,1,4);
imagesc(sgbfb_features);
axis xy;
ylabel('SGBFB');

