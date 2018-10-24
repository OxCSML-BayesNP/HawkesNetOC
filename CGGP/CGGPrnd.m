function [W, T, w0, beta] = CGGPrnd(p, alpha, sigma, tau, Fdist, gamma, varargin)

% CGGPrnd samples points of a generalized gamma process.
% [W, T, w0, beta] = CGGPrnd(p, alpha, sigma, tau, Fdist)
% [__] = CGGPrnd(p, alpha, sigma, tau, Fdist, gamma)
% [__] = CGGPrnd(p, alpha, sigma, tau, Fdist, gamma, T)
% [__] = CGGPrnd(p, alpha, sigma, tau, Fdist, gamma, T, maxiter)
%
% Samples the points of the CGGP with Levy measure
%   nu(w(1:p)) = exp(-sum(gamma.*w)) int_0^infty w0^(-p) F(w(1:p)/w0) rho0(w0)dw0
% where rho0(w0) = alpha/Gamma(1-sigma) * w0^(-1-sigma) * exp(-tau*w0)
%
% For sigma>=0, it samples points s.t. w0 is above the threshold T>0 using the adaptive
% thinning strategy described in Favaro and Teh (2013).
% -------------------------------------------------------------------------
% INPUTS
%   - p: positive integer
%   - alpha: positive scalar
%   - sigma: real in (-Inf, 1)
%   - tau: positive scalar
%   - Fdist: distribution F
%
% Optional inputs
%   - gamma: positive tilting parameters (default=0)
%   - T: truncation threshold; positive scalar
%   - maxiter: maximum number of iterations for the adaptive thinning
%     strategy (default=1e8)
%
% OUTPUTS
%   - W: points of the GGP
%   - T: threshold
%
% See also GGPrnd
% -------------------------------------------------------------------------
% EXAMPLE
% p = 3; alpha = 100; sigma = 0.5; tau = 1e-4;
% Fdist.name = 'gamma'; Fdist.param = 0.01;
% W = CGGPrnd(p, alpha, sigma, tau, Fdist);

% -------------------------------------------------------------------------
% Reference:
% S. Favaro and Y.W. Teh. MCMC for normalized random measure mixture
% models. Statistical Science, vol.28(3), pp.335-359, 2013.

% Copyright (c) F. Caron (University of Oxford), A. Todeschini (Inria), and 
% X. Miscouridou (University of Oxford)
% caron@stats.ox.ac.uk
% adrien.todeschini@gmail.com
% xenia.miscouridou@spc.ox.ac.uk
% September 2017
%--------------------------------------------------------------------------

if nargin<6
    gamma = zeros(p,1);
elseif any(size(gamma)~=[p,1])
    error('gamma must be a scalar or vector of size [p,1]')
end


% Check the parameters of the CGGP
CGGPcheckparams(p, alpha, sigma, tau, Fdist, gamma);

switch(Fdist.name)
    case 'gamma'
        if isnumeric(Fdist.param)
            a = Fdist.param';
            b = a;
        else
            a = Fdist.param.a';
            b = Fdist.param.b';
        end
end

% Sample w0
%[w0, T] = tiltGGPrndNEW(log(alpha), sigma, tau, a, b, gamma', varargin{:});
[w0, T] = tiltGGPrnd(log(alpha), sigma, tau, a, b, gamma', varargin{:});

beta = scoreCGGPrnd(p, w0, Fdist, gamma);

W = exp(bsxfun(@plus, log(w0), log(beta)));
end
