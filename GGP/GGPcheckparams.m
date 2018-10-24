function GGPcheckparams(alpha, sigma, tau)

%GGPcheckparams checks the parameters of a GGP.
% GGPcheckparams(alpha, sigma, tau)
%
%   Valid range is 
%       alpha>0, tau>0, sigma in (-infinity, 1) 
%           or
%       alpha>0, tau=0, sigma in (0,1)
% 
%   The function returns an error if the parameters are outside of this
%   range.

% Copyright (C) Francois Caron, University of Oxford
% caron@stats.ox.ac.uk
% April 2015
%--------------------------------------------------------------------------

if tau<0
    error('tau must be nonnegative');
end
if alpha<=0
    error('alpha must be strictly positive')
end
if sigma>=1
    error('sigma must be in (-infinity, 1)');
end
if (tau<1e-8 && sigma<1e-8)
    error('If tau==0, sigma must be in (0,1)')
end

end
