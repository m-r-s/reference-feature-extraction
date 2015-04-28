function features = mfcc(log_mel_spec, num_coeff, delta, ddelta)
% usage: features = mfcc(log_mel_spec)
%   log_mel_spec   log Mel-spectrogram or similar spectro-temporal representation
%
% usage: features = mfcc(log_mel_spec, num_coeff, delta, ddelta)
%   num_coeff      number of MFCCs to output
%   delta          temporal context to calculate first discrete derivative (0 - off)
%   ddelta         temporal context to calculate second discrete derivative (0 - off)
%
% - Mel frequency cepstral coefficients v1.0 -
%
% This script extracts Mel frequency cepstral coefficient features
% from a spectro-temporal representation, e.g., from a log Mel-spectrogram.
% The sample rate of the spectro-temporal representation is assumed to be
% 100 Hz and the spectral resolution is assumed to be about 1 ERB per channel.
% A detailed explanation is given in [1].
%
% Copyright (C) 2015 Marc René Schädler
% E-mail marc.r.schaedler@uni-oldenburg.de
% Institute Carl-von-Ossietzky University Oldenburg, Germany
%
%-----------------------------------------------------------------------------
%
% Release Notes:
% v1.0 - Inital release
%

%% Default settings and checks

% Get number of bands from input
num_bands = size(log_mel_spec,1);

% Default number of MFCCs covers at least 13/23 cycles per band
if nargin < 2 || isempty(num_coeff)
  min_cycl_band = 13./23;
  num_coeff = ceil(num_bands.*min_cycl_band);
end

% Default context for first discrete derivative
if nargin < 3 || isempty(delta)
  delta = 2;
end

% Default context for second discrete derivative
if nargin < 4 || isempty(ddelta)
  if delta > 0
    ddelta = 2;
  else
    ddelta = 0;
  end
end

% Maximum context
context = delta+ddelta;

%% Calculate Mel frequency cepstral coefficients

% Temporally pad log Mel-spectrogram by repeating first and last frames
if context > 0
  log_mel_spec = [repmat(log_mel_spec(:,1),1,context) log_mel_spec repmat(log_mel_spec(:,end),1,context)];
end

% Perform DCT and select coefficients
mfccs = dct(log_mel_spec);
mfccs = mfccs(1:num_coeff,:);
features = mfccs;

% Temporal processing
if delta > 0
  delta_filter = linspace(-1,1,2.*delta+1);
  deltas = conv2(mfccs, delta_filter, 'same'); ...
  features = [mfccs; deltas];
  if ddelta > 0
    ddelta_filter = linspace(-1,1,2.*ddelta+1);
    ddeltas = conv2(deltas, ddelta_filter, 'same'); ...
    features = [mfccs; deltas; ddeltas];
  end
end

% Remove padded context
if context > 0
  features = features(:,(1+context):(end-context));
end
end
