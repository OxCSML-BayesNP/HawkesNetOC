function out = GGPpsi(t, alpha, sigma, tau)

%GGPPSI returns the Laplace exponent of a GGP.
% out = GGPPSI(t, alpha, sigma, tau)
%
% The Laplace exponent of a GGP evaluated at t is
% psi(t) = -log ( E[exp(-t * sum_i w_i)] )
%        = alpha/sigma * ( (t+tau)^sigma - tau^sigma))
% where
% (w_i)_{i=1,2,..} are the points of a Poisson process on R_+ of mean measure 
% rho(dw) = alpha/Gamma(1-sigma) * w^{-1-sigma} * exp(-tau*w)dw
% -------------------------------------------------------------------------
% INPUTS
%   - t: vector of positive scalars of length n
%   - alpha: positive scalar
%   - sigma: real in (-inf, 1)
%   - tau: positive scalar
%
% OUTPUT
%   - out: Laplace exponent evaluated at the values t
% -------------------------------------------------------------------------
% EXAMPLE
% t = .1:.1:10;
% alpha = 100; sigma = 0.5; tau = 1;
% out = GGPpsi(t, alpha, sigma, tau);

% Copyright (C) Francois Caron, University of Oxford
% caron@stats.ox.ac.uk
% April 2015
%--------------------------------------------------------------------------

if (sigma==0) % gamma process
    out = alpha * log( 1+ t/tau );
else
    out = alpha/sigma * ((t+tau).^sigma - tau^sigma);
end
    