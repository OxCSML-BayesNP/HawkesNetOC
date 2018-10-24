function [samples, stats] = GGPbipgraphmcmc(G, modelparam, mcmcparam, verbose)

%GGPbipgraphmcmc runs a MCMC sampler for the GGP bipartite graph model
% [samples, stats] = GGPbipgraphmcmc(G, modelparam, mcmcparam, verbose)
%
% -------------------------------------------------------------------------
% INPUTS
%   - G: sparse logical adjacency matrix
%   - modelparam: structure of model parameters with the following fields:
%       - alpha: cell of length 2, corresponding to each type of node.
%         If alpha{i} is scalar, the value of alpha{i}. If vector of length
%         2, parameters of the gamma prior over alpha{i}
%       - sigma: cell of length 2, corresponding to each type of node.
%         if sigma{i} is scalar, the value of sigma. If vector of length
%         2, parameters of the gamma prior over (1-sigma{i})
%       - tau: cell of length 2, corresponding to each type of node.
%         if tau{i} scalar, the value of tau. If vector of length
%         2, parameters of the gamma prior over tau
%   - mcmcparam: structure of mcmc parameters with the following fields:
%       - niter: number of MCMC iterations
%       - nburn: number of burn-in iterations
%       - thin: thinning of the MCMC output
%       - hyper.MH_nb: number of MH iterations
%       - hyper.rw_std: standard deviation of the random walk
%       - store_w: logical. If true, returns MCMC draws of w1 and w2
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

if nargin<4
    verbose = true;
end

% Check G
if ~issparse(G) || ~islogical(G)
    error('Adjacency matrix G must be a sparse logical matrix');
end

G = sparse(G);
[N1, N2] = size(G);
[ind_w1, ind_w2] = find(G);

% struct array of logicals indicating whether to estimate parameters
[alpha1, sigma1, tau1] = GGPgetparams(modelparam(1).alpha,...
    modelparam(1).sigma, modelparam(1).tau, 'empty');
[alpha2, sigma2, tau2] = GGPgetparams(modelparam(2).alpha,...
    modelparam(2).sigma, modelparam(2).tau, 'empty');

estim(1).alpha = isempty(alpha1);
estim(1).sigma = isempty(sigma1);
estim(1).tau = isempty(tau1);
estim(2).alpha = isempty(alpha2);
estim(2).sigma = isempty(sigma2);
estim(2).tau = isempty(tau2);

% Parameters of the MCMC algorithm
niter = mcmcparam.niter;
nburn = mcmcparam.nburn;
thin = mcmcparam.thin;

% To store MCMC samples
n_samples = (niter-nburn)/thin;
if mcmcparam.store_w
    samples.w1 = zeros(N1, n_samples, 'double');
    samples.w2 = zeros(N2, n_samples, 'double');
else
    samples.w1 = [];
    samples.w2 = [];
end
samples.w1_rem = zeros(1, n_samples);
samples.w2_rem = zeros(1, n_samples);
samples.alpha1 = zeros(1, n_samples);
samples.sigma1 = zeros(1, n_samples);
samples.tau1 = zeros(1, n_samples);
samples.alpha2 = zeros(1, n_samples);
samples.sigma2 = zeros(1, n_samples);
samples.tau2 = zeros(1, n_samples);

% Initialize
[alpha1, sigma1, tau1] = GGPgetparams(modelparam(1).alpha,...
    modelparam(1).sigma, modelparam(1).tau, 'init');
[alpha2, sigma2, tau2] = GGPgetparams(modelparam(2).alpha,...
    modelparam(2).sigma, modelparam(2).tau, 'init');

logalpha1 = log(alpha1);
logalpha2 = log(alpha2);

w1 = exp(randn(N1, 1));%ones(N1, 1);
w2 = exp(randn(N2, 1));%ones(N2, 1);
w1_rem = exp(randn);%1;

w1_rep = sparse(ind_w1, ind_w2, w1(ind_w1), N1, N2);
m_w2 = full(sum(G, 1));
m_w1 = full(sum(G, 2)');
tic
for i=1:niter
    if rem(i, 2000)==0 && verbose
        fprintf('i=%i; ', i);
        fprintf('alpha1=%.2f; ', alpha1);
        fprintf('sigma1=%.3f; ', sigma1);
        fprintf('tau1=%.2f; ', tau1);
        fprintf('alpha2=%.2f; ', alpha2);
        fprintf('sigma2=%.3f; ', sigma2);
        fprintf('tau2=%.2f;\n', tau2);
    end
    
    % Sample U
    Umod = sample_U(w1, w2, ind_w1, ind_w2); % Umod = U-1 to have a sparse vector
    
    % Sample hyperparameters of w2
    [alpha2, logalpha2, sigma2, tau2] =...
        update_hyper(m_w2, Umod, w1, w1_rem, w1_rep, alpha2, logalpha2, sigma2, tau2,...
        estim(2), modelparam(2), mcmcparam.hyper.rw_std, mcmcparam.hyper.MH_nb);
    
    % Sample w2
    [w2, w2_rem, w2_rep] = sample_w(m_w2, Umod, w1, w1_rem, w1_rep, alpha2, sigma2, tau2);
    
    % Sample hyperparameters of w1
    [alpha1, logalpha1, sigma1, tau1] =...
        update_hyper(m_w1, Umod', w2, w2_rem, w2_rep, alpha1, logalpha1, sigma1, tau1,...
        estim(1), modelparam(1), mcmcparam.hyper.rw_std, mcmcparam.hyper.MH_nb);
    
    % Sample w1
    [w1, w1_rem, w1_rep] = sample_w(m_w1, Umod', w2, w2_rem, w2_rep, alpha1, sigma1, tau1);
    
    % Print some information
    if i==10
        time = toc * niter/10;
        hours = floor(time/3600);
        minutes = (time - hours*3600)/60;
        fprintf('-----------------------------------\n')
        fprintf('Start MCMC for bipartite GGP graphs\n')
        fprintf('Nb of nodes: (%d,%d) - Nb of edges: %d\n',size(G, 1),size(G, 2), full(sum(G(:))));
        fprintf('Number of iterations: %d\n', niter)
        fprintf('Estimated computation time: %.0f hour(s) %.0f minute(s)\n',hours,minutes);
        fprintf('Estimated end of computation: %s \n', datestr(now + time/3600/24));
        fprintf('-----------------------------------\n')
    end
    
    % Store output
    if (i>nburn && rem((i-nburn),thin)==0)
        ind = ((i-nburn)/thin);
        if mcmcparam.store_w
            samples.w1(:, ind) = w1;
            samples.w2(:, ind) = w2;
        end
        samples.w1_rem(:, ind) = w1_rem;
        samples.w2_rem(:, ind) = w2_rem;
        samples.alpha1(:, ind) = alpha1;
        samples.sigma1(:, ind) = sigma1;
        samples.tau1(:, ind) = tau1;
        samples.alpha2(:, ind) = alpha2;
        samples.sigma2(:, ind) = sigma2;
        samples.tau2(:, ind) = tau2;
    end
    
end

samples.logalpha1 = log(samples.alpha1);
samples.logalpha2 = log(samples.alpha2);

stats = struct();

time = toc;
hours = floor(time/3600);
minutes = (time - hours*3600)/60;
fprintf('-----------------------------------\n')
fprintf('End MCMC for bipartite GGP graphs\n')
fprintf('Computation time: %.0f hour(s) %.0f minute(s)\n',hours,minutes);
fprintf('-----------------------------------\n')
end

%% ------------------------------------------------------------------------
% MAJOR SUBFUNCTIONS
% -------------------------------------------------------------------------


function Umod = sample_U(w1, w2, ind_w1, ind_w2)

N1 = numel(w1);
N2 = numel(w2);

% Sample U conditional on w1 and w2
w1_w2 = w1(ind_w1).*w2(ind_w2);
% U = sparse(ind_w1, ind_w2, rexprnd(w1_w2, 1), N1, N2);
Umod = sparse(ind_w1, ind_w2, rexprnd(w1_w2, 1) - 1, N1, N2); % Umod = U- 1 to have a sparse matrix
end

function [w2, w2_rem, w2_rep] = sample_w(m_w2, Umod, w1, w1_rem, w1_rep, alpha, sigma, tau)

[N1, N2] = size(Umod);
[ind_w1, ind_w2] = find(Umod);

sum_w1 = sum(w1);
% w1_rep = sparse(ind_w1, ind_w2, w1(ind_w1), N1, N2);
w1_U = w1_rem + sum_w1 + sum(w1_rep.* Umod); % sum_i w1_i * u_ij
w2 = gamrnd(m_w2 - sigma, 1./(tau+ w1_U))';

w2_rem = GGPsumrnd(alpha, sigma, tau + sum_w1 + w1_rem);

w2_rep = sparse(ind_w2, ind_w1, w2(ind_w2), N2, N1);
end

function [alpha, logalpha, sigma, tau, rate2] =...
    update_hyper(m_w2, Umod, w1, w1_rem, w1_rep, alpha, logalpha, sigma, tau,...
    estim, hyper, rw_std, MH_nb)

N2 = numel(m_w2);
sum_w1 = sum(w1) + w1_rem;
w1_U = sum_w1 + sum(w1_rep.* Umod);
for nn=1:MH_nb
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
            alphaprop = gamrnd(N2, 1/(GGPpsi(sum_w1, 1, sigmaprop, tauprop) ));
            logalphaprop = log(alphaprop);
        else
            alphaprop = alpha;
            logalphaprop = logalpha;
        end
    else % more stable numerically as alpha can take very large values in that case, we sample alpha2=alpha*tau^sigma
        if estim.alpha
            alpha2prop = gamrnd(N2, 1/( GGPpsi(sum_w1/tauprop, 1, sigmaprop, 1) ));%exp(log(alpha) + rw_std*randn);
            logalphaprop = log(alpha2prop) - sigmaprop*log(tauprop);
            alphaprop = exp(logalphaprop);
        else
            alphaprop = alpha;
            logalphaprop = logalpha;
        end
    end
    
    [~, logkappa] = GGPkappa(m_w2, w1_U, 1, sigma, tau);
    [~, logkappaprop] = GGPkappa(m_w2, w1_U, 1, sigmaprop, tauprop);
    logaccept = sum(logkappaprop - logkappa);
    
    if estim.alpha
        logaccept = logaccept ...
            + N2 * (log(GGPpsi((sum_w1)/tau, 1, sigma, 1) ) + sigma*log(tau)...
            - log(GGPpsi((sum_w1)/tauprop, 1, sigmaprop, 1) ) - sigmaprop*log(tauprop) );
        if hyper.alpha(1)>0
            logaccept = logaccept + hyper.alpha(1)*( logalphaprop - logalpha);
        end
        if hyper.alpha(2)>0
            logaccept = logaccept - hyper.alpha(2) * (alphaprop - alpha);
        end
    else
        logaccept = logaccept ...
            - GGPpsi(sum_w1, alphaprop, sigmaprop, tauprop) ...
            + GGPpsi(sum_w1, alpha, sigma, tau);
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
    
    if log(rand)<logaccept
        alpha = alphaprop;
        logalpha = logalphaprop;
        sigma = sigmaprop;
        tau = tauprop;
    end
end
rate2 = min(1, exp(logaccept));
end
