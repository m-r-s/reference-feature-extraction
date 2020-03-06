function [epsi epsi_std] = epsi(snr1, correct1, total1, snr2, correct2, total2)
% usage: [epsi epsi_std] = epsi(snr1, correct1, total1, snr2, correct2, total2)
%   snr       Vector of signal-to-noise ratios (SNR) in dB
%   correct   Vector of number of correct decisions at the corresponding SNRs
%   total     Total number of decisions (may also be a vector)
%   epsi      Equal-performance SNR increase (EPSI)
%   epsi_std  Estimated standard deviation of the EPSI
%
% - Equal-performance SNR increase v1.0 -
%
% This script calculates the increase in signal-to-noise ratio (SNR) which
% is needed on average in order to make System2 perform as well as System1.
% If the equal-performance SNR increase (EPSI) is positive the SNR needs
% to be increased, which would indicate that System2 is less robust than
% System1. A negative EPSI would inidcate that System2 is more robust than
% System1 because the SNR for System2 needs to be decrased to match the
% average recognition performance of System1 and System2.
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

%%%%%%%%%%%%%%%%%%%%%%%%% Example with data from [1] %%%%%%%%%%%%%%%%%%%%%%%%%
%  SNR           = [   -6    -3    +0    +3    +6    +9]; % SNR in dB
%  HSR           = [90.30 93.00 93.80 95.30 96.80 98.80]; % percent correct
%  MFCC          = [68.66 74.58 82.25 87.50 89.08 92.00]; % percent correct
%  GBFB          = [71.41 77.75 84.25 88.91 92.16 92.66]; % percent correct
%  ASR_decisions = 1200; % Number of independent decisions
%  HSR_decisions = 1200; % Number of independent decisions
%  HSR_correct   = HSR./100.*HSR_decisions;
%  MFCC_correct  = MFCC./100.*ASR_decisions;
%  GBFB_correct  = GBFB./100.*ASR_decisions;
%
%  % Calculate the EPSI of GBFB over HSR
%  [HSR_GBFB_epsi   HSR_GBFB_epsi_std] = epsi( ...
%    SNR,   HSR_correct,    HSR_decisions, ...
%    SNR,   GBFB_correct,   ASR_decisions);
%  [HSR_MFCC_epsi   HSR_MFCC_epsi_std] = epsi( ...
%    SNR,   HSR_correct,    HSR_decisions, ...
%    SNR,   MFCC_correct,   ASR_decisions);
%  [GBFB_MFCC_epsi  GBFB_MFCC_epsi_std] = epsi( ...
%    SNR,   GBFB_correct,   ASR_decisions, ...
%    SNR,   MFCC_correct,   ASR_decisions);
%    
%  disp(['The SNR needs to be increased by ' ...
%    num2str(HSR_GBFB_epsi) '+-' num2str(HSR_GBFB_epsi_std) ' dB' ...
%    ' for GBFB in order to make GBFB and HSR perform equally.']);
%  disp(['The SNR needs to be increased by ' ...
%    num2str(HSR_MFCC_epsi) '+-' num2str(HSR_MFCC_epsi_std) ' dB' ...
%    ' for MFCC in order to make MFCC and HSR perform equally.']);
%  disp(['The SNR needs to be increased by ' ...
%    num2str(GBFB_MFCC_epsi) '+-' num2str(GBFB_MFCC_epsi_std) ' dB' ...
%    ' for MFCC in order to make MFCC and GBFB perform equally.']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EPSI implementation %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of bootstrap samples for error estimation
num_bootstrap = 1000; 

% Make sure all data is stored in column vectors
snr1 = snr1(:);
correct1 = correct1(:);
total1 = total1(:);
snr2 = snr2(:);
correct2 = correct2(:);
total2 = total2(:);

% Expand total numbers if no SNR-dependent vectors are supplied
if numel(total1) == 1
  total1 = total1.*ones(size(snr1));
end
if numel(total2) == 1
  total2 = total2.*ones(size(snr2));
end

% Sort results by SNR and estimate the uncertainties
[snr1 sortidx1] = sort(snr1);
correct1 = correct1(sortidx1);
total1 = total1(sortidx1);
performance1 = correct1./total1;
uncertainty1 = estd(round(total1), round(correct1), num_bootstrap);

[snr2 sortidx2] = sort(snr2);
correct2 = correct2(sortidx2);
total2 = total2(sortidx2);
performance2 = correct2./total2;
uncertainty2 = estd(round(total2), round(correct2), num_bootstrap);

% Calculate the shift in SNR to match the performance of the two systems
epsi = shiftmetric(snr1, performance1, snr2, performance2);

% Estimate the EPSI standard deviation by bootstrapping and assuming normally
% distributed deviations on the performance scale.
% Please note: The shift is independent of the scale of the y-axis but its
% estimated deviations are not. Ideally the performance would be reported on
% a scale on which the deviations could be expected to be normally distributed
% in order to better match this assumtion.
epsi_std = NaN;
if ~isnan(epsi)
  bootstrap = zeros(num_bootstrap,1);
  for k=1:num_bootstrap
    performance1_tmp = performance1 + randn(size(performance1)).*uncertainty1;
    performance2_tmp = performance2 + randn(size(performance2)).*uncertainty2;
    bootstrap(k) = shiftmetric(snr1, performance1_tmp, snr2, performance2_tmp);
  end
  if any(~isfinite(bootstrap))
    warning('Ignoring infinite numbers during error estimation!');
  end
  epsi_std = std(bootstrap(isfinite(bootstrap)));
end
end


function shift = shiftmetric(x1, y1, x2, y2)
% Calulcates the average shift of two graphs, y1(x1) and y2(x2), on the x-axis.
% The average is taken in steps of 0.5 over the x-range in which both graphs
% are defined. Hence, the calculated shift does not depend on the scale of the
% y-axis.

% Make sure that data comes in column vectors
x1 = x1(:);
y1 = y1(:);
x2 = x2(:);
y2 = y2(:);

% Make sure the performance graphs are monotonic (c.f. Equation 5 in [1])
y1 = monotonify(-y1(end:-1:1),0.0001);
y2 = monotonify(-y2(end:-1:1),0.0001);
y1 = -y1(end:-1:1);
y2 = -y2(end:-1:1);

shift = NaN;
assert(all(size(x1) == size(y1)) && all(size(x2) == size(y2)), ...
  'x and y vectors must have same lengths');

% Get the overlapping window on the y-axis.
% This is the region in which the two graphs can be compared.
yrange = [max(min(y1),min(y2)); min(max(y1),max(y2))];

% Interpolate the corrsponding ranges on the x-axis.
xrange1 = interp1(y1, x1, yrange);
xrange2 = interp1(y2, x2, yrange);

% Sample the x-ranges in 0.5 dB steps
xrange1 = [ceil(xrange1(1)*2)/2; floor(xrange1(2)*2)/2];
xpoints1 = (xrange1(1):0.5:xrange1(2)).';
xrange2 = [ceil(xrange2(1)*2)/2; floor(xrange2(2)*2)/2];
xpoints2 = (xrange2(1):0.5:xrange2(2)).';

% Integrate the shifts in both directions if sampling points exist
if ~isempty(xpoints1) && ~isempty(xpoints2)
  % Look up the y-values at the sample points
  ypoints1 = interp1(x1, y1, xpoints1, 'linear','extrap');
  ypoints2 = interp1(x2, y2, xpoints2, 'linear','extrap');
  % Get the x-values of the "other" system at the same y-value
  xpoints12 = interp1(y2, x2, ypoints1, 'linear','extrap');
  xpoints21 = interp1(y1, x1, ypoints2, 'linear','extrap');
  % Average over the difference on the x-values
  shift1 = mean(xpoints12-xpoints1);
  shift2 = mean(xpoints21-xpoints2);
  % Weight the calculated differences from both graphs equally
  shift = 0.5.*(shift1-shift2);
end

end


function y = monotonify(x, epsilon)
% Guarantees monotonically increasing output
y = x;
if nargin < 2
  epsilon = eps;
end
for i=2:length(x)
  y(i) = max(x(i-1)+epsilon, y(i));
end
end


function rate_std = estd(total, correct, n)
% Estimates the standard deviation of the correct-rate of a finite number
% of binary decisions by bootstrapping
if numel(total) > 1
  % Handle the situation when total is a vector
  if numel(n) == 1
    n = n.*ones(size(total));
  end
  rate_std = arrayfun(@estd, total, correct, n);
else
  % Generate a vector of ones and zeros reflecting correct and false decisions
  x = [zeros(1,total-correct) ones(1,correct)];
  % Generate n random selections drawing total samples from x
  selections = ceil(rand(total,n)*total);
  % Get the rates for these selections
  rates = mean(x(selections));
  % Get the standard deviation of the artificially generated samples
  rate_std = std(rates);
end
end
