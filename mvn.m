function out = mvn(in)
% usage out = mean_variance_norm(in)
%
% Perform mean and variance normalization of each row of 'in'
%
% Copyright (C) 2015 Marc René Schädler
%

out = bsxfun(@minus, in, mean(in,2));
out = bsxfun(@times, out, 1./sqrt(mean(out.^2,2)));
end
