function out = heq(in, points)
% usage: out = heq(in, points)
%   in        values
%   points    number of fixed points for mapping
%   out       normalized values
%
% - Histogram equalization v1.0 -
%
% This script performs a histogram equalization (HEQ) of each row of 'in'
% using 'points' fixed percentiles/quantiles for the mapping.
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

% Default number of fixed points
if nargin < 2 || isempty(points)
  points = 100;
end

% Calculate the expected minimum and maximum quantiles
% when drawing 'context' samples from the unknown distribution
% C.f. Equation 4 in [1]
context = size(in,2);
exptected_min = 1-context./(context+1);
exptected_max = context./(context+1);

% Define source and target quantiles
source_quantiles = linspace(0, 1, points);
target_quantiles = linspace(exptected_min, exptected_max, points);

% Get source quantiles from data
quantiles = quantile(in, source_quantiles, 2);

% Allocate memory
out = zeros(size(in));

% Loop over all rows
for i = 1:size(in,1)
  % Handle the case if all input values are almost equal
  if (quantiles(i,end) - quantiles(i,1)) < 100.*eps
    out(i,:) = 0.5;
  else
    mask = [true diff(quantiles(i,:))>0];
    out(i,:) = interp1q(quantiles(i,mask).', target_quantiles(mask).', in(i,:).');
  end
end

% Map the quantiles to the Gaussian distribution
% using the inverse error function
out = erfinv(out.*2 - 1);
end
