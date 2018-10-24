function [S_w, T, S_w_small] = CGGPsumrnd(p, logalpha, sigma, tau, Fdist, gamma, T, varargin)

% CGGPsumrnd samples from the distribution of the total mass of a CGGP.
% S_w = CGGPsumrnd(p, logalpha, sigma, tau, Fdist, gamma, T)
%
%   It generates a realization of the random variable S with Laplace
%   transform: 
%   E[e^-(sum(t.*S))] = exp(-alpha \psi(t)), t=(t_1,.., t_p)
%   \psi(t) see Eq. (40) in 
%   Reference: 
%   A. Todeschini, X. Miscouridou and F. Caron (2016) <https://arxiv.org/abs/1602.02114 Exchangeable Random Measures for Sparse and Modular Graphs with Overlapping Communities>. arXiv:1602.02114.
%
% -------------------------------------------------------------------------
% INPUTS
%   - p: positive integer
%   - logalpha: scalar
%   - sigma: real in (-Inf, 1)
%   - tau: positive scalar
%   - Fdist: distribution F
%
% Optional inputs
%   - gamma: positive tilting parameters (default=0)
%   - T: truncation threshold; positive scalar
%
% OUTPUTS
%   - S_w: positive vector of size [1,p]
%   - T: threshold
%
% See also CGGPrnd
% -------------------------------------------------------------------------
% EXAMPLE
% p = 3; alpha = 100; sigma = 0.5; tau = 1e-4;
% Fdist.name = 'gamma'; Fdist.param = 0.01;
% S_w = CGGPsumrnd(p, log(alpha), sigma, tau, Fdist);
% -------------------------------------------------------------------------

% Copyright (c) F. Caron (University of Oxford), A. Todeschini (Inria), and 
% X. Miscouridou (University of Oxford)
% caron@stats.ox.ac.uk
% adrien.todeschini@gmail.com
% xenia.miscouridou@spc.ox.ac.uk
% September 2017
%--------------------------------------------------------------------------


if nargin < 6
    gamma = zeros(p,1);
elseif any(size(gamma)~=[p,1])
    error('gamma must be a scalar or vector of size [p,1]')
end

if nargin < 7
    T = 1e-5;
end

switch Fdist.name
    case 'gamma'
        if isnumeric(Fdist.param)
            a = Fdist.param';
            b = a;
        else
            a = Fdist.param.a';
            b = Fdist.param.b';
        end
    otherwise
        error('Unknown distribution F')
end


if sigma < 0
    %% Finite activity
    % Compound Poisson parametric case (when sigma<0)    
    T = 0;
    S_w = zeros(1, p);
    S_w_small = zeros(1,p);
    rate = exp( logalpha - log(-sigma) + sigma*log(tau) ); 
    K = poissrnd(rate);
    if K>1e8
        warning('Generating %d jumps: alpha=%.2f, sigma=%.2e, tau=%.2f', K, exp( logalpha), sigma, tau);
    end
    gamma_b = gamma'./b;
    % avoid creating too big arrays
    while K > 0
        n = min(1e7, K);
        K = K-n;
        w0 = gamrnd(-sigma, 1/tau, n, 1);
        w0 = w0(w0>0);
        accept = rand(size(w0)) < prod(bsxfun(@power, 1+bsxfun(@times, w0, gamma_b), -a), 2);
        w0 = w0(accept);
        
        if ~isempty(w0)
            beta = scoreCGGPrnd(p, w0, Fdist, gamma);
            try
                w = exp(bsxfun(@plus, log(w0), log(beta)));
            catch
                keyboard
            end
            S_w = S_w + sum(w, 1);
        end
    end
    return;
end
    
if sigma == 0 % gamma process
    error('******* Gamma Process (sigma=0) not implemented yet **********')
end

%% Infinite activity

if exist('tCGGPsumrnd', 'file')==3 %%% mexfile
    S_w = tCGGPsumrnd(p, logalpha, sigma, tau, a, b, gamma', T, double(randi(intmax)));
else
    [w, T] = CGGPrnd(p, exp(logalpha), sigma, tau, Fdist, gamma, T);
    S_w = sum(w,1);
end

% Gaussian approximation below T
S_w_small = CGGPsumsmallrnd(p, logalpha, sigma, tau, a, b, gamma', T, varargin{:});

S_w = S_w + S_w_small;
