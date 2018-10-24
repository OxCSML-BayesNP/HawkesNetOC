function [N, T] = tiltGGPrnd(logalpha, sigma, tau, a, b, gamma, T, maxiter)

% tiltGGPrnd samples points of a generalized gamma process with exponentially
% tilted Levy measure and it also handles the case of tau=0
%
%   Levy measure:
%   \tilted{\rho}_0 = M(-w_0 \gamma_1,...,-w_0 \gamma_p)\rho_0 
%   where rho_0 = alpha/Gamma(1-sigma) * w^(-1-sigma) * exp(-tau*w)
%   and M(t_1,...t_k) is the mgf of F from which betas are sampled:
%   beta_1,...,beta_p  ~ F, w_0 ~ \tilted{rho_0}, w_1,...,w_p ~

% For sigma>=0, it samples points above the threshold T>0 using a thinning,
% or adaptive thinning strategy described in Favaro and Teh (2013).
% -------------------------------------------------------------------------
% INPUTS
%   - logalpha: scalar
%   - sigma: real in (-Inf, 1)
%   - tau: nonnegative scalar
%   - a: shape parameters of the gamma score distribution F
%   - b: rate parameters of the gamma score distribution F
%   - gamma: vectors of size [1,p]; positive tilting parameters
%
% Optional inputs
%   - T: truncation threshold; positive scalar
%   - maxiter: maximum number of iterations for the adaptive thinning
%     strategy (default=1e8); if nargin<8, it uses thinning
%
% OUTPUTS
%   - N: points of the GGP
%   - T: threshold
% -------------------------------------------------------------------------
%
% [N, T] = tiltGGPrnd(alpha, sigma, tau, a, b, gamma, T)
% [__] = tiltGGPrnd(alpha, sigma, tau, a, b, gamma, T, maxiter)

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


if nargin < 8
    maxiter = 1e8;
elseif maxiter<=0
    error('maxiter must be a strictly positive integer');
end

if all(gamma==0)
    if nargin < 7
        [N, T] = GGPrnd(exp(logalpha), sigma, tau);
    elseif nargin < 8
        [N, T] = GGPrnd(exp(logalpha), sigma, tau, T);
    else
        [N, T] = GGPrnd(exp(logalpha), sigma, tau, T, maxiter);
    end
    return
end

% Check the parameters of the GGP
GGPcheckparams(exp(logalpha), sigma, tau);

gamma_b = gamma./b;

%% Finite activity
if sigma < -1e-8
    % Compound Poisson parametric case (when sigma<0)  
    T = 0;  
    rate = exp( logalpha - log(-sigma) + sigma*log(tau) ); 
    K = poissrnd(rate);
    if K>1e8
        warning('Generating %d jumps', K)
    end
    N = gamrnd(-sigma, 1/tau, K, 1);
    N = N(N>0);
    accept = rand(size(N)) < prod(bsxfun(@power, 1+bsxfun(@times, N, gamma_b), -a), 2);
    N = N(accept);
    return;
end


%% Infinite activity

%%%%%%%% FIXME
sigma = max(sigma, 0); %%% set sigma=0 if in [-1e-8, 0]

Njumps = [];
if nargin<7
    % set the threshold automatically so that we sample of the order Njumps jumps
    % Number of jumps of order alpha/sigma/Gamma(1-sigma) * T^{-sigma} for sigma>0
    % and alpha*log(T) for sigma=0
    if sigma>.1
        Njumps = 20000; % Expected number of jumps
        T = exp(1/sigma*(logalpha - log(sigma) - gammaln(1-sigma) - log(Njumps)));
    else
        T = 1e-10;
    end
elseif T<=0
    error('Threshold T must be strictly positive');
end
    
if isempty(Njumps)
    if sigma>1e-3
        Njumps = floor(exp(logalpha - log(sigma) -gammaln(1-sigma) - sigma*log(T)));
    else
        Njumps = floor(-exp(logalpha)*log(T));
    end
end

if Njumps > 1e7
%%% TEMP %%%%%%%%%%
    warning('T too small - Its value was increased at %f', T*2);
    [N, T] = tiltGGPrnd(logalpha, sigma, tau, a, b, gamma, T*2, maxiter);
    return
end

%% Adaptive thinning
% with truncation level T
% Method described in Favaro and Teh, 2012
N = zeros(ceil(Njumps+3*sqrt(Njumps)), 1);
t = T;
k = 0;

if tau < 1e-8
    %% case tau==0
    
    log_cst = logalpha - gammaln(1-sigma) - log(sigma);
    msigma = -sigma;
    msigmainv = -1/sigma;
    
    for i=1:maxiter
        log_mgf = -sum( a.*log1p(t.*gamma_b) );
        log_r = log(-log(rand)) - log_mgf - log_cst;
        
        if log_r > msigma*log(t)
            completed = true;
            break;
        end
        
        t_new = exp( msigmainv * log( t^msigma - exp(log_r) ) );
        log_mgf_new = -sum( a.*log1p(t_new.*gamma_b) );
        
        if log(rand) < log_mgf_new - log_mgf
            k = k+1;
            N(k,1) = t_new;
        end
        t = t_new;
    end
    
else
    %% case tau>0
    
    log_cst = logalpha - gammaln(1-sigma) - log(tau);
    sigmap1 = sigma+1;
    tauinv = 1/tau;
    
    for i=1:maxiter
        log_mgf = -sum( a.*log1p(t.*gamma_b) );
        log_r = log(-log(rand));
        log_G = log_mgf + log_cst - sigmap1*log(t) - tau*t;
        
        if log_r > log_G
            completed = true;
            break;
        end
        
        t_new = t - log1p(-exp(log_r-log_G)) * tauinv;
        log_mgf_new = -sum( a.*log1p(t_new.*gamma_b) );
        
        if log(rand) < sigmap1*(log(t)-log(t_new)) + log_mgf_new - log_mgf
            k = k+1;
            N(k,1) = t_new;
        end
        t = t_new;
    end
    
end

N = N(1:k,1);

% If too many computions, we increase the threshold T and rerun
if ~completed
    warning('T too small - Its value was increased at %f', T*10);
    [N, T] = tiltGGPrnd(logalpha, sigma, tau, a, b, gamma, T*10, maxiter);
end

end
