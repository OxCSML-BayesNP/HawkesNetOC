function [out, logout] = GGPkappa(n, z, alpha, sigma, tau)

% GGPkappa returns the nth moment of a tilted GGP.
% [out, logout] = GGPkappa(n, z, alpha, sigma, tau)
%
%   kappa(n,z) = int_0^infty w^n * exp(-zw) * rho(w)dw
%              = alpha / (z+tau)^(n-sigma) * gamma(n-sigma)/gamma(1-sigma) 
%       where rho(w) = alpha/gamma(1-sigma) * w^(-1-sigma) * exp(-tau*w)
%       is the L�vy measure of a generalized gamma process
% -------------------------------------------------------------------------
% INPUTS
%   - n: strictly positive integer
%   - z: positive scalar
%   - alpha: positive scalar
%   - sigma: real in (-Inf, 1)
%   - tau: positive scalar
% 
% OUTPUTS
%   - out: kappa(n,z)
%   - logout: log(kappa(n,z))
% -------------------------------------------------------------------------
% EXAMPLE
% n = 10; z = 1.5; alpha = 100; sigma = 0.5; tau = 1;
% out = GGPkappa(n, z, alpha, sigma, tau);

% Copyright (C) Francois Caron, University of Oxford
% caron@stats.ox.ac.uk
% April 2015
%--------------------------------------------------------------------------

logout = log(alpha) - (n-sigma).*log(z+tau) + gammaln(n-sigma) - gammaln(1-sigma);
out = exp(logout);

end  