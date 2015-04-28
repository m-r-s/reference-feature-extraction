function features = sgbfb(log_mel_spec, omega_max, size_max, nu, distance, phases)
% usage: features = sgbfb(log_mel_spec)
%   log_mel_spec    log Mel-spectrogram or similar spectro-temporal representation
%
% usage: features = sgbfb(log_mel_spec, omega_max, size_max, nu, distance, phases)
%   omega_max       upper center modulation frequencies in [radian radian]
%   size_max        maximum filter size (lower center modulation frequency) in [bands, frames]
%   nu              half-waves under the envelope in [spectral temporal] dimension
%   distance        spacing of filters (<1)
%   phases          Phase of [spectral temporal] modulation filters (may be a cell with more than one pair)
%
% - Separated Gabor filter bank v1.0 -
%
% This script extracts separable Gabor filter bank features
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

% Default upper center modulation frequencies
if nargin < 2 || isempty(omega_max)
  omega_max = [pi/2 pi/2];
end

% Get number of bands from input
num_bands = size(log_mel_spec,1);

% Default number of half-waves under the envelope
if nargin < 4 || isempty(nu)
  nu = [3.5 3.5];
end

% Default maximum filter size (lower center modulation frequencies)
if nargin < 3 || isempty(size_max)
  size_max = [3*num_bands 40];
end

% Default spacing of filters
if nargin < 5 || isempty(distance)
  distance = [0.3 0.2];
end

% Default Phases of modulation filters [spectral temporal]
if nargin < 6 || isempty(phases)
  phases = {[0 0], [0 pi/2], [pi/2 0], [pi/2 pi/2]};
end

if ~iscell(phases)
  phases = {phases};
end

% Maximum context
context = floor(size_max(2)/2);


%% Calculate separated Gabor filter bank features

% Temporally pad log Mel-spectrogram by repeating first and last frames
log_mel_spec = [repmat(log_mel_spec(:,1),1,context) log_mel_spec repmat(log_mel_spec(:,end),1,context)];

% Iterate feature extraction for all phase pairs
streams = length(phases);
features = cell(streams,1);
for i=1:streams
  % Get phase pair: [spectral temporal]
  phase = phases{i};
  
  % Note: The order of spectral and temporal modulation processing
  % does not matter
  
  % Spectral 1D GBFB processing
  features_spec = cell2mat(gbfb1d(log_mel_spec, ...
    omega_max(1), size_max(1), nu(1), distance(1), phase(1), 1, 1));
  
  % Temporal 1D GBFB processing
  features_spec_temp = cell2mat(gbfb1d(features_spec.', ...
    omega_max(2), size_max(2), nu(2), distance(2), phase(2), 1, 0).').';
  
  features{i} = features_spec_temp;
end
features = cell2mat(features);

% Remove padded context
features = features(:,(1+context):(end-context));
end
 

%%%%%%%%%% 1D Gabor filter bank (1D-GBFB) functions %%%%%%%%%%
% Most of this code is modified from the reference GBFB implementation [2]

function [out rep_output] = gbfb1d(in, omega_max, size_max, nu, distance, phase, dc, rep)
% One dimensional Gabor filter bank filtering

% Get the center modulation frequencies
omega = gbfb_axis(omega_max,size_max,nu,distance);

% Add the low-pass filter by adding frequency 0
if dc
  omega = [0 omega];
end

% For each center modulation frequency
out = cell(length(omega),1);
rep_output = cell(length(omega),1);
for i=1:length(omega)
  % Generate a 1D Gabor filter
  gfilter = gfilter_gen(omega(i), nu, phase, size_max);
  
  % Filter the input
  in_filtered = conv2(in, gfilter, 'same');
  
  % If only representative channels should be kept, do so
  if rep
    [in_filtered rep_output{i}] = gfilter_rep(gfilter, in_filtered);
  end
  
  out{i} = in_filtered;
end
end


function omega = gbfb_axis(omega_max, size_max, nu, distance)
% Calculates the center modulation frequencies (c.f. [2])
% omega_max   Maximum center modulation frequency in rad
% size_max    Maximum extension of a Gabor filter in samples
% nu          Number of half-waves under the envelope
% distance    Distance between adjacent filters
omega_min = (pi * nu) / size_max;
c = distance * 8 / nu;
space = (1 + c/2) / (1 - c/2);
count = 0;
omega(1) = omega_max;
while omega(end)/space > omega_min
  omega(1+count) = omega_max/space^count;
  count = count + 1;
end
omega = fliplr(omega);
end


function gfilter = gfilter_gen(omega, nu, phi, size_max)
% Generates a 1D Gabor filter function
% omega       Center frequency in rad
% nu          Number of half-waves under the envelope
% phi         Phase of gfilter at the center sample in rad
% size_max    Maximum size to determine when to use the envelope
%             as the filter function

% Caching whitelist (feel free to add Matlab versions)
caching = is_octave();

if caching
  % Build a config id string
  config = strrep(sprintf('c%.0f', [omega nu phi size_max]*10000),'-','_');
  % Load cache
  persistent cache;
end

% Only generate filters which are not cached
if ~caching || isempty(cache) || ~isfield(cache, config)
  w = 2*pi / abs(omega) * nu / 2;
  if w > size_max
    w = size_max;
    omega = 0;
  end
  envelope = hann_win(w); % C.f. Equation 1a in [1]
  win_size = length(envelope);
  x_0 = (win_size+1) / 2;
  x = 1:win_size;
  sinusoid = exp(1i.*(omega*(x - x_0) + phi)); % C.f. Equation 1b in [1]
  gfilter  = real(envelope(:) .* sinusoid(:)); % C.f. Equation 1c in [1]
  envelope_mean = mean(mean(envelope));
  gfilter_mean = mean(mean(gfilter));
  if (omega ~= 0)
    gfilter = gfilter - envelope./envelope_mean .* gfilter_mean;
  end
  gfilter = gfilter ./ max(abs(fft(gfilter)));
  
  if caching
    % Save to cache
    cache.(config) = gfilter;
  end
else
  % Load from cache
  gfilter = cache.(config);
end
end


function [out, k_idx] = gfilter_rep(gfilter, in)
% Sub-sampling of 'in' in the first dimension at a rate equal to a quarter 
% of the size of 'gfilter' in the same dimension
k_factor = floor(1/4 * size(gfilter,1));
if k_factor < 1
  k_factor = 1;
end
k_offset = mod(floor(size(in,1)/2),k_factor);
k_idx = (1+k_offset):k_factor:size(in,1);
out = in(k_idx,:);
end


function window_function = hann_win(width)
% A hanning window of "width" with the maximum centered on the center sample
x_center = 0.5;
step = 1/width;
right = x_center:step:1;
left = x_center:-step:0;
x_values = [left(end:-1:1) right(2:end)].';
valid_values_mask = (x_values > 0) & (x_values < 1);
window_function = 0.5 * (1 - ( cos(2*pi*x_values(valid_values_mask))));
end


function r = is_octave ()
  persistent x;
  if (isempty (x))
    x = exist ('OCTAVE_VERSION', 'builtin');
  end
  r = x;
end
