function S = GGPsumrnd(alpha, sigma, tau)

%GGPsumrnd samples from the distribution of the total mass of a GGP.
% S = GGPsumrnd(alpha, sigma, tau)
%
%   It generates a realization of the random variable S with Laplace
%   transform
%   E[e^-(t*S)] = exp(-alpha/sigma * [(t+tau)^sigma - tau^sigma])
% -------------------------------------------------------------------------
% INPUTS
%   - alpha: positive scalar
%   - sigma: real in (-Inf, 1)
%   - tau: positive scalar
%
% OUTPUTS
%   - S: positive scalar
% -------------------------------------------------------------------------
% EXAMPLE
% alpha = 100; sigma = 0.5; tau = 1;
% S = GGPsumrnd(alpha, sigma, tau);
% -------------------------------------------------------------------------

% Copyright (C) Francois Caron, University of Oxford
% caron@stats.ox.ac.uk
% April 2015
%--------------------------------------------------------------------------

if sigma<-10^-8
    % Compound Poisson case
    % S is distributed from a Poisson mixture of gamma variables
    K = poissrnd(-alpha/sigma/tau^(-sigma));
    S = gamrnd(-sigma*K, 1/tau);
elseif sigma < 10^-8
    % Gamma process case
    % S is gamma distributed
    S = gamrnd(alpha, 1/tau);
elseif sigma==0.5 && tau==0
    % Inverse Gaussian process case
    % S is distributed from an inverse Gaussian distribution
    lambda = 2*alpha^2;
    mu = alpha/sqrt(tau);
    S = igaussrnd(mu, lambda, 1, 1);
else
    % General case
    % S is distributed from an exponentially tilted stable distribution
    S = etstablernd(alpha/sigma, sigma, tau);
end
