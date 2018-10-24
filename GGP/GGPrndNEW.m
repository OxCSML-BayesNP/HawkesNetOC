function [N, T] = GGPrndNEW(alpha, sigma, tau, T, maxiter)

% GGPrndNEW samples points of a generalized gamma process exactly as
% GGPrnd.m with the addition that:
% 1) it handles the case of tau=0
% 2) for \sigma \in [-1e-8,1e-8] the jumps are sampled from the truncated
% measure to avoid numerical instabilities
% [N, T] = GGPrndNEW(alpha, sigma, tau, T)
%
% Samples the points of the GGP with Levy measure
%   alpha/Gamma(1-sigma) * w^(-1-sigma) * exp(-tau*w)
%
% For sigma>=0, it samples points above the threshold T>0 using a thinning,
% or adaptive thinning strategy described in Favaro and Teh (2013).
% -------------------------------------------------------------------------
% INPUTS
%   - alpha: positive scalar
%   - sigma: real in (-Inf, 1)
%   - tau: non negative scalar
%
% Optional inputs
%   - T: truncation threshold; positive scalar
%   - maxiter: maximum number of iterations for the adaptive thinning
%     strategy (default=1e8); if nargin<5, it uses thinning
%
% OUTPUTS
%   - N: points of the GGP
%   - T: threshold
% -------------------------------------------------------------------------
% EXAMPLE
% alpha = 100; sigma = 0.5; tau = 1e-4;
% N = GGPrnd(alpha, sigma, tau);

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


%%%%%%%%%%%%%%%%%%
%%% TODO: return expected number of connections not generated
%%% TODO: C++ with openmp implementation
%%% FC: ADD TRUNCATION for finite-activity CRM when sigma is close to
%%% zero to avoid numerical problems
%%%%%%%%%%%%%%%%%


% Check the parameters of the GGP
GGPcheckparams(alpha, sigma, tau);

%% Finite activity GGP
if sigma<-1e-8
    rate = exp( log(alpha) - log(-sigma) + sigma*log(tau) );
    K = poissrnd(rate);
    N = gamrnd(-sigma, 1/tau, K, 1);
    N = N(N>0);
    T = 0;
    return;
end

%% Infinite activity GGP


if nargin<5 % Use simple thinning
    sigma = max(sigma, 0); %%% set sigma=0 if in [-1e-8, 0] - This needs to be fixed
    
    if ~(sigma==0) & tau ==0
        lograte = log(alpha)  - log(sigma) - gammaln(1-sigma) + log(T^(-sigma));
        Njumps = poissrnd(exp(lograte));
        N = T*rand(Njumps, 1).^(-1/sigma);
        return;
    end
    a = 5; %truncation level for pareto
    % Use a truncated Pareto on [0,a]
    if ~(sigma==0)
        lograte = log(alpha)  - log(sigma) - tau*T - gammaln(1-sigma) + log(T^(-sigma) - a^(-sigma));
        Njumps = poissrnd(exp(lograte));
        log_N1 = - 1/sigma * log(-(rand(Njumps, 1) * (a^sigma-T^sigma)-a^sigma)./ (a*T)^sigma); % Sample from truncated Pareto
    else
        lograte = log(alpha)  - tau*T - gammaln(1-sigma) + log(log(a) - log(T));
        Njumps = poissrnd(exp(lograte));
        log_N1 = (rand(Njumps, 1) * (log(a)-log(T)) + log(T)); 
    end   
    
    N1 = exp(log_N1);
    ind1  = log(rand(Njumps, 1)) < tau*(T - N1);sum(ind1)/length(ind1);
    
    % Use a truncated exponential on (a,+infty)
    lograte = log(alpha) - tau*a - (1+sigma)*log(a) - log(tau) - gammaln(1-sigma);
    Njumps = poissrnd(exp(lograte));
    log_N2 = texprnd(tau, a, Njumps, 1);    % Sample from truncated exponential        
    ind2  = log(rand(Njumps, 1)) < -(1+sigma)* (log_N2 - log(a) );
    N = [N1(ind1); exp(log_N2(ind2))];
    

%     if T > sigma/tau
%         % Uses truncated exponential
%         lograte = log(alpha) - tau*T - (1+sigma)*log(T) - log(tau) - gammaln(1-sigma);
%         Njumps = poissrnd(exp(lograte))
%         log_N = texprnd(tau, T, Njumps, 1);    % Sample from truncated exponential        
%         ind  = log(rand(Njumps, 1)) < -(1+sigma)* (log_N - log(T) );
%         N = exp(log_N(ind));
%     else
%         % Uses Pareto
%         lograte = log(alpha) - tau*T - log(sigma) - sigma*log(T) - gammaln(1-sigma);
%         Njumps = poissrnd(exp(lograte));
%         log_N = log(T) - 1/sigma * log(rand(Njumps, 1)); % Sample from Pareto
%         N = exp(log_N);
%         ind  = log(rand(Njumps, 1)) < tau*(T - N);
%         N = N(ind);
%         
%     end

else % Use adaptive thinning
    %sigma = max(sigma, 0); %%% set sigma=0 if in [-1e-8, 0]


    Njumps = [];
    if nargin<4
        % set the threshold automatically so that we sample of the order Njumps jumps
        % Number of jumps of order alpha/sigma/Gamma(1-sigma) * T^{-sigma} for sigma>0
        % and alpha*log(T) for sigma=0
        if sigma>.1
            Njumps = 20000; % Expected number of jumps
            T = exp(1/sigma*(log(alpha) - log(sigma) - gammaln(1-sigma) - log(Njumps)));
        else
            T = 1e-10;
        end
    elseif T<=0
        error('Threshold T must be strictly positive');
    end

    if isempty(Njumps)
        if sigma>1e-3
            Njumps = floor(exp(log(alpha) - log(sigma) -gammaln(1-sigma) - sigma*log(T)));
        else
            Njumps = floor(-alpha*log(T));
        end
    end

    if Njumps > 1e8
        warning('Expected number of jumps = %d - press key if you wish to continue', Njumps);
        pause
    end

    if nargin < 5
        maxiter = 1e8;
    elseif maxiter<=0
        error('maxiter must be a strictly positive integer');
    end

    N = zeros(ceil(Njumps+3*sqrt(Njumps)), 1);
    t = T;
    k = 0;

    if tau<1e-8
        %% case tau==0
        log_cst = log(alpha) - gammaln(1-sigma) - log(sigma);
        msigma = -sigma;
        msigmainv = -1/sigma;

        for k=1:maxiter
            log_r = log(-log(rand)) - log_cst;
            if log_r > msigma * log(t)
                completed = true;
                break;
            end
            t = exp( msigmainv * log( t^msigma - exp(log_r) ) );
            N(k) = t;
        end
    else
        %% case tau>0
        % Adaptive thinning strategy
        log_cst = log(alpha)-gammaln(1-sigma)-log(tau);
        sigmap1 = 1+sigma;
        tauinv = 1/tau;

        for i=1:maxiter
            log_r = log(-log(rand)); % Sample exponential random variable e of unit rate
            log_G = log_cst-sigmap1*log(t)-tau*t;
            if log_r > log_G
                completed = true;
                break;
            end
            t_new = t-tauinv*log(1-exp(log_r-log_G));
            if log(rand) < sigmap1*(log(t)-log(t_new))
                k = k+1;
                N(k) = t_new;
            end
            t = t_new;
        end
    end

    N = N(1:k);

    % If too many computions, we increase the threshold T and rerun
    if ~completed
        T = T*10;
        warning('T too small - Its was value increased at %f', T);
        N = GGPrndNEW(alpha, sigma, tau, T);
    end

end


