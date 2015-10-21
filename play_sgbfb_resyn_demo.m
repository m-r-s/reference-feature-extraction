close all;
clear;
clc;

%% Preparation

% Lets generate a 1000 Hz windowed tone signal
fs = 16000;
signal = sin(2*pi*linspace(0,2000/2,fs/2)) .* hanning(fs/2).';

% Or just load some signal
%[signal, fs] = wavread('...');

% Calculate the log Mel-spectrogram
log_melspec = log_mel_spectrogram(signal, fs);

% Use SGBFB features (only for two phase combinations and without histogram equalization)
feature_extraction = @(x) sgbfb(x, [], [], [], [], {[0 0], [pi/2 pi/2]});


%% Learn resynthesis matrix as SGBFB is linear
context = 20; % Maximum temporal context in frames that SGBFB filters have
num_context = 20*2+1;
num_noise_samples = 20000;
noise = randn(size(log_melspec,1), num_noise_samples);
noise_features = feature_extraction(noise);
% Create context
noise_context = zeros(size(noise).*[num_context,1]);
for i=-context:context
  noise_context(1+size(noise,1)*(context+i):size(noise,1)*(context+i+1),max(1,1+i):min(end,end+i)) = noise(:,max(1,1+i):min(end,end+i));
end
sgbfb_resyn_matrix = noise_context(:,1+context:end-context) / noise_features(:,1+context:end-context);


%% Resynthesis
features = feature_extraction(log_melspec);

% Apply pseudo-inverse SGBFB feature extraction
log_melspec_resyn_context = sgbfb_resyn_matrix * features;

% Merge/average context
log_melspec_resyn = zeros(size(log_melspec_resyn_context)./[num_context,1]);
for i=-context:context
  log_melspec_resyn(:,max(1,1+i):min(end,end+i)) = log_melspec_resyn(:,max(1,1+i):min(end,end+i)) + ...
    log_melspec_resyn_context(1+size(log_melspec_resyn,1)*(context+i):size(log_melspec_resyn,1)*(context+i+1),max(1,1+i):min(end,end+i));
end
log_melspec_resyn = log_melspec_resyn./num_context;


%% Compare results
subplot(2,1,1);
imagesc(log_melspec); colorbar;

subplot(2,1,2);
imagesc(log_melspec_resyn); colorbar;

