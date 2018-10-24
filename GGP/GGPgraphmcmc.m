function [samples, stats] = GGPgraphmcmc(G, modelparam, mcmcparam, typegraph, verbose)

%GGPgraphmcmc runs a MCMC sampler for the GGP graph model
% [samples, stats] = GGPgraphmcmc(G, modelparam, mcmcparam, typegraph, verbose)
%
% -------------------------------------------------------------------------
% INPUTS
%   - G: sparse logical adjacency matrix
%   - modelparam: structure of model parameters with the following fields:
%       - alpha: if scalar, the value of alpha. If vector of length
%         2, parameters of the gamma prior over alpha
%       - sigma: if scalar, the value of sigma. If vector of length
%         2, parameters of the gamma prior over (1-sigma)
%       - tau: if scalar, the value of tau. If vector of length
%         2, parameters of the gamma prior over tau
%   - mcmcparam: structure of mcmc parameters with the following fields:
%       - niter: number of MCMC iterations
%       - nburn: number of burn-in iterations
%       - thin: thinning of the MCMC output
%       - leapfrog.L: number of leapfrog steps
%       - leapfrog.epsilon: leapfrog stepsize
%       - latent.MH_nb: number of MH iterations for latent (if 0: Gibbs update)
%       - hyper.MH_nb: number of MH iterations for hyperparameters
%       - hyper.rw_std: standard deviation of the random walk
%       - store_w: logical. If true, returns MCMC draws of w
%   - typegraph: type of graph ('undirected' or 'simple')
% Optional inputs:
%   - verbose: logical. If true (default), print information
%
% OUTPUTS
%   - samples: structure with the MCMC samples for the variables
%       - w
%       - w_rem
%       - alpha
%       - logalpha
%       - sigma
%       - tau
%   - stats: structure with summary stats about the MCMC algorithm
%       - rate: acceptance rate of the HMC step at each iteration
%       - rate2: acceptance rate of the MH for the hyperparameters at
%         each iteration
%
% See also graphmcmc, graphmodel
% -------------------------------------------------------------------------

% Copyright (C) Francois Caron, University of Oxford
% caron@stats.ox.ac.uk
% April 2015
%--------------------------------------------------------------------------

if nargin<5
    verbose = true;
end

% Check G
if ~isequal(G,G') || ~issparse(G) || ~islogical(G)
    error('Adjacency matrix G must be a symmetric sparse logical matrix');
end

if strcmp(typegraph, 'simple')
    issimple = true;
else
    issimple = false;
end

% Hyperparameter initialization
if numel(modelparam.alpha)==2
    alpha = 100*rand;
    estim.alpha = 1;
else
    alpha = modelparam.alpha;
    estim.alpha = 0;
end
logalpha = log(alpha);
if numel(modelparam.sigma)==2
    sigma = 2*rand - 1;
    estim.sigma = 1;
else
    sigma = modelparam.sigma;
    estim.sigma = 0;
end
if numel(modelparam.tau)==2
    tau = 10*rand;
    estim.tau = 1;
else
    tau = modelparam.tau;
    estim.tau = 0;
end

K = size(G, 1);
if issimple % If no self-loops
    G2 = triu((G+G')>0, 1);
else
    G2 = triu((G+G')>0);
end
[ind1, ind2] = find(G2);

% Initialisation
n = randi(10,size(ind1) );
count = sparse(ind1, ind2, n, K, K);
N = sum(count,1)' + sum(count, 2);
w = gamrnd(1,1,K,1);%ones(K, 1);
logw = log(w);
w_rem = gamrnd(1,1);%1;

% Parameters of the MCMC algorithm
niter = mcmcparam.niter;
nburn = mcmcparam.nburn;
thin = mcmcparam.thin;
L = mcmcparam.leapfrog.L; % Number of leapfrog steps
epsilon = mcmcparam.leapfrog.epsilon/K^(1/4); % Leapfrog stepsize
% Choice of update for the latent (Gibbs/MH)
if mcmcparam.latent.MH_nb==0
    % Gibbs update
    update_n = @(logw, d, K, count, ind1, ind2)...
        update_n_Gibbs(logw, K, ind1, ind2);
else
    % Metropolis-Hastings update
    update_n = @(logw, d, K, count, ind1, ind2) ...
        update_n_MH(logw, d, K, count, ind1, ind2, mcmcparam.latent.MH_nb);
end

% To store MCMC samples
n_samples = (niter-nburn)/thin;
if mcmcparam.store_w
    samples.w = zeros(K, n_samples, 'double');
else
    samples.w = [];
end
samples.w_rem = zeros(1, n_samples);
samples.alpha = zeros(1, n_samples);
samples.logalpha = zeros(1, n_samples);
samples.sigma = zeros(1, n_samples);
samples.tau = zeros(1, n_samples);

rate = zeros(niter, 1);
rate2 = zeros(niter, 1);

for i=1:niter
    if verbose && rem(i, 2000)==0
        fprintf('i=%i; ', i);
        fprintf('alpha=%.2f; ', alpha);
        fprintf('sigma=%.3f; ', sigma);
        fprintf('tau=%.2f;\n', tau);
    end

    % Update w using Hamiltonian Monte Carlo
    [w, logw, rate(i)] = update_w(w, logw, w_rem, N, L, epsilon, sigma, tau, issimple);
    if i<mcmcparam.leapfrog.nadapt % Adapt the stepsize
        epsilon = exp(log(epsilon) + .01*(mean(rate(1:i)) - 0.6));
    end

    % Update w_rem and hyperparameters using Metropolis-Hastings
    if rem(i,2)==0 % alternate between random walk and gamma proposals for alpha
        rw_alpha = true;
    else
        rw_alpha = false;
    end
    [w_rem, alpha, logalpha, sigma, tau, rate2(i)] = update_hyper(w, logw, w_rem, alpha, logalpha, sigma,...
        tau, estim, modelparam, mcmcparam.hyper.MH_nb, mcmcparam.hyper.rw_std, rw_alpha);

    % Update the counts using Gibbs or Metropolis-Hastings step
    [N, n, count] = update_n(logw, n, K, count, ind1, ind2);

    if isnan(alpha)
        keyboard
    end

    % Store output
    if (i>nburn && rem((i-nburn),thin)==0)
        ind = ((i-nburn)/thin);
        if mcmcparam.store_w
            samples.w(:, ind) = w;
        end
        samples.w_rem(:, ind) = w_rem;
        samples.logalpha(:, ind) = logalpha;
        samples.alpha(:, ind) = alpha;
        samples.tau(:, ind) = tau;
        samples.sigma(:, ind) = sigma;
    end

end

stats.rate = rate;
stats.rate2 = rate2;

end


%% ------------------------------------------------------------------------
% Subfunctions
% ------------------------------------------------------------------------


%% Gradient
function out = grad_U(N, w,w_rem, sigma, tau, issimple)

if issimple
    out = - (N - sigma) + w*(tau+2*sum(w)+2*w_rem) - 2*w.^2;
else
    out = - (N - sigma) + w*(tau+2*sum(w)+2*w_rem);
end

end

%% Update of the weights w
function [w, logw, rate] = update_w(w, logw, w_rem, N, L, epsilon, sigma, tau, issimple)

% L = 30;
% Update w
sum_w = sum(w);
sumall_w = sum_w + w_rem;

K = numel(w);
% wprop = w;
logwprop = logw;
p = randn(K, 1);
grad1 = grad_U(N, w, w_rem, sigma, tau, issimple);
pprop = p - epsilon * grad1/2;
for lp=1:L
%             wprop = exp(log(wprop) + epsilon*pprop);%wprop = exp(log(wprop) + epsilon*pprop./m_leapfrog);
    logwprop = logwprop + epsilon*pprop;
    if lp~=L
        pprop = pprop  - epsilon * grad_U(N, exp(logwprop), w_rem, sigma, tau, issimple);
    end
end
wprop = exp(logwprop);
pprop = pprop - epsilon/2 * grad_U(N, wprop, w_rem, sigma, tau, issimple);
% min(pprop)
% min(logwprop)

sum_wprop = sum(wprop);
sumall_wprop = sum_wprop + w_rem;

temp1 = - sumall_wprop^2 + sumall_w^2 ...
    + sum((N-sigma-1).*(logwprop - logw) )...
    - tau * (sum_wprop - sum_w);

logaccept = temp1 -.5*sum(pprop.^2-p.^2) -sum(logw) + sum(logwprop);
if issimple % If simple graph, do not take into account self-connections
    logaccept = logaccept + sum(wprop.^2) - sum(w.^2);
end

if isnan(logaccept)
    logaccept = -Inf;
end

if log(rand)<logaccept
    w = wprop;
    logw = logwprop;
end
rate = min(1, exp(logaccept));

end

%% Update of the latent counts (Gibbs)
function [N, d, count] = update_n_Gibbs(logw, K, ind1, ind2)

% Update the latent counts from the full conditional
% rate_poi = 2*w(ind1).*w(ind2) - w(ind1).^2.*(ind1==ind2);
lograte_poi = log(2) + logw(ind1) + logw(ind2);
lograte_poi(ind1==ind2) = 2*logw(ind1(ind1==ind2));
d = tpoissrnd(exp(lograte_poi));
count = sparse(ind1, ind2, d, K, K);
N = sum(count,1)' + sum(count, 2);
end

%% Update of the latent counts (MH)
function [N, d, count] = update_n_MH(logw, d, K, count, ind1, ind2, nbMH)

for i=1:nbMH
    % Metropolis-Hastings update for the latent counts
    lograte_poi = log(2) + logw(ind1) + logw(ind2);
    lograte_poi(ind1==ind2) = 2*logw(ind1(ind1==ind2));
    ind = (d ==1);
    dprop = d;
    dprop(ind) = 2;
    if sum(~ind)>0
        dprop(~ind) = dprop(~ind) + 2*randi(2,[sum(~ind),1]) - 3;
    end

    logqprop = zeros(size(ind));
    logqprop(~ind) = log(.5);

    indbis = (dprop==1);
    logq = zeros(size(indbis));
    if sum(~indbis)>0
        logq(~indbis) = log(.5);
    end

    diff_d = (dprop - d);
    logaccept_d = diff_d.* lograte_poi ...
        - gammaln(dprop +1) + gammaln(d + 1)...
        -logqprop + logq;
    % This is slightly faster than the above, but not much
    %     logaccept_d = zeros(size(d));
    %     indter = (dprop>d);
    %     logaccept_d(indter) = lograte_poi(indter) - log(dprop(indter));
    %     logaccept_d(~indter) = - lograte_poi(~indter) + log(d(~indter));
    %     logaccept_d = logaccept_d -logqprop + logq;

    indaccept = (log(rand(size(logaccept_d))) <logaccept_d);
    d(indaccept) = dprop(indaccept);
    count = count +  sparse(ind1(indaccept),ind2(indaccept),  diff_d(indaccept), K, K);
    N = sum(count,1)' + sum(count, 2);
end

end


%% Update of the hyperparameters
function [w_rem, alpha, logalpha, sigma, tau, rate2] = update_hyper(w, logw, w_rem, alpha, logalpha, sigma,...
    tau, estim, hyper, nbMH, rw_std, rw_alpha)

K = numel(w);
for nn=1:nbMH
    sum_w = sum(w);
    sumall_w = sum_w + w_rem;

    % Sample (alpha,sigma,tau,w*) from the proposal distribution
    if estim.sigma
        sigmaprop = 1-exp(log(1-sigma) + rw_std(1)*randn);
    else
        sigmaprop = sigma;
    end
    if estim.tau
        tauprop = exp(log(tau) + rw_std(2)*randn);
    else
        tauprop = tau;
    end
    if sigmaprop>-1
        if estim.alpha
            if ~rw_alpha % gamma proposal
                alphaprop = gamrnd(K, 1/( GGPpsi(2*sum_w + 2*w_rem, 1, sigmaprop, tauprop) ));
            else
                alphaprop = alpha.*exp(.02*randn);
            end
            logalphaprop = log(alphaprop);
        else
            alphaprop = alpha;
            logalphaprop = logalpha;
        end

            wprop_rem = GGPsumrnd(alphaprop, sigmaprop, tauprop + 2*sum_w + 2*w_rem);
    else % more stable numerically as alpha can take very large values in that case, we sample alpha2=alpha*tau^sigma
        if estim.alpha
            if ~rw_alpha % gamma proposal
                alpha2prop = gamrnd(K, 1/( GGPpsi((2*sum_w + 2*w_rem)/tauprop, 1, sigmaprop, 1) ));
                logalphaprop = log(alpha2prop) - sigmaprop*log(tauprop);
            else % random walk proposal
                logalphaprop = logalpha + .02*randn;
            end
            alphaprop = exp(logalphaprop);
            rate_K = exp( logalphaprop - log(-sigmaprop) + sigmaprop*log(tauprop + 2*sum_w + 2*w_rem ) );
            num_clust = poissrnd(rate_K);
            wprop_rem = gamrnd(-sigmaprop* num_clust, 1/(tauprop+ 2*sum_w + 2*w_rem));
        else
            alphaprop = alpha;
            logalphaprop = logalpha;
            wprop_rem = GGPsumrnd(alphaprop, sigmaprop, tauprop + 2*sum_w + 2*w_rem);
        end
    end

    % Compute the acceptance probability
    sum_wprop = sum(w);
    sumall_wprop = sum_wprop + wprop_rem;

    temp1 = - sumall_wprop^2 + sumall_w^2 ...
        + (sigma - sigmaprop) * sum(logw)...
        - (tauprop - tau - 2*wprop_rem + 2*w_rem) * sum_w;
    temp2 =   K* (gammaln(1-sigma) - gammaln(1-sigmaprop));

     logaccept = temp1 + temp2;
     if estim.alpha
         if ~rw_alpha
             logaccept = logaccept ...
                 + K * (log(GGPpsi((2*sum_wprop + 2*wprop_rem)/tau, 1, sigma, 1) ) + sigma*log(tau)...
                 - log(GGPpsi((2*sum_w + 2*w_rem)/tauprop, 1, sigmaprop, 1) ) - sigmaprop*log(tauprop) );
         else
            logaccept = logaccept ...
                 - exp(logalphaprop + sigmaprop*log(tauprop))* GGPpsi((2*sum_w + 2*w_rem)/tauprop, 1, sigmaprop, 1) ...
            + exp(logalpha + sigma*log(tau)) * GGPpsi((2*sum_wprop + 2*wprop_rem)/tau, 1, sigma, 1)...
            + K*(logalphaprop - logalpha);
         end

         if hyper.alpha(1)>0
             logaccept = logaccept + hyper.alpha(1)*( logalphaprop - logalpha);
         end
         if hyper.alpha(2)>0
              logaccept = logaccept - hyper.alpha(2) * (alphaprop - alpha);
         end
     else
         logaccept = logaccept ...
             - GGPpsi(2*sum_w + 2*w_rem, alphaprop, sigmaprop, tauprop) ...
            + GGPpsi(2*sum_wprop + 2*wprop_rem, alpha, sigma, tau);
     end
     if estim.tau
         logaccept = logaccept ...
             + hyper.tau(1)*( log(tauprop) - log(tau)) - hyper.tau(2) * (tauprop - tau);
     end
     if estim.sigma
         logaccept = logaccept ...
             + hyper.sigma(1)*( log(1 - sigmaprop) - log(1-sigma)) ...
             - hyper.sigma(2) * (1 - sigmaprop - 1 + sigma);
     end

     if isnan(logaccept)
         keyboard
     end
     % Accept step
     if log(rand)<logaccept
        w_rem = wprop_rem;
        alpha = alphaprop;
        logalpha = logalphaprop;
        sigma = sigmaprop;
        tau = tauprop;
    end
end
rate2 = min(1, exp(logaccept));
end
