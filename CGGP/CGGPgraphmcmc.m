function [samples, stats] = CGGPgraphmcmc(G, modelparam, mcmcparam, typegraph, varargin)

% CGGPgraphmcmc runs a MCMC sampler for the CGGP graph model
% [samples, stats] = CGGPgraphmcmc(G, modelparam, mcmcparam, typegraph, verbose)
%
% -------------------------------------------------------------------------
% INPUTS
%   - G: sparse logical adjacency matrix
%   - modelparam: struct array of length 2, corresponding to each type of node
%       constaining model parameters with the following fields:
%       - p: positive integer. number of features
%       - alpha: If scalar, the value of alpha is fixed. If vector of length 2,
%         parameters of the gamma prior over alpha.
%       - sigma: If scalar, the value of sigma is fixed. If vector of length 2,
%         parameters of the gamma prior over (1-sigma).
%       - tau: If scalar, the value of tau is fixed. If vector of length 2,
%         parameters of the gamma prior over tau.
%       - Fdist: struct with fields 'name' and 'param'
%           - Case Fdist.name = 'gamma':
%             If Fdist.param is numeric, Gamma(Fdist.param, Fdist.param)
%             If Fdist.param is struct, Gamma(Fdist.param.a, Fdist.param.b)
%             If the parameter has one column, its value is fixed.
%             If two columns, parameters of the gamma prior over the parameter.
%             If one row, the p components of the parameter are considered equal.
%       - gamma: If column vector of length p, the value of gamma is fixed.
%         If matrix of size [p,2], parameters of the gamma prior over gamma.
%   - mcmcparam: structure of mcmc parameters with the following fields:
%       - niter: number of MCMC iterations
%       - nburn: number of burn-in iterations
%       - thin: thinning of the MCMC output
%       - hyper.MH_nb: number of MH iterations
%       - hyper.rw_std: standard deviation of the random walk
%       - store_w: logical. If true, returns MCMC draws of w1 and w2
%   - typegraph: string with value 'simple' (no self-loops) or undirected
%                (with self-loops)
% Optional inputs:
%   - verbose: logical. If true (default), print information
%
% OUTPUTS
%   - samples: structure array of length 2, corresponding to each type of node
%     with the MCMC samples for the variables
%       - w (if mcmcparam.store_w=true)
%       - w_rem
%       - alpha
%       - logalpha
%       - sigma
%       - tau
%       - gamma
%       - Fparam
%   - stats: structure with summary stats about the MCMC algorithm
%       - rate: acceptance rate of the HMC step at each iteration
%       - rate2: acceptance rate of the MH for the hyperparameters at
%         each iteration
%
% See also graphmcmc, graphmodel
% -------------------------------------------------------------------------

% Copyright (c) F. Caron (University of Oxford), A. Todeschini (Inria), and 
% X. Miscouridou (University of Oxford)
% caron@stats.ox.ac.uk
% adrien.todeschini@gmail.com
% xenia.miscouridou@spc.ox.ac.uk
% September 2017
%--------------------------------------------------------------------------

%% parse input arguments
ip = inputParser;

addRequired(ip, 'G', @(x) islogical(x) || isnumeric(x));
addRequired(ip, 'modelparam', @isstruct);
addRequired(ip, 'mcmcparam', @isstruct);
addOptional(ip, 'verbose', true, @islogical);
addParameter(ip, 'T', 1e-2, @isnumeric); % truncation level for adaptive thinning
addParameter(ip, 'like_model', 'ber', @(x) validatestring(x, {'ber', 'poiss'}));
addParameter(ip, 'sumsmall_arg', {}, @iscell);
addParameter(ip, 'integral_arg', {}, @iscell);
addParameter(ip, 'r_miss', 1, @isnumeric);
addParameter(ip, 'init', struct(), @isstruct);

parse(ip, G, modelparam, mcmcparam, varargin{:});
verbose = ip.Results.verbose;
T = ip.Results.T;
like_model = ip.Results.like_model;
sumsmall_arg = ip.Results.sumsmall_arg;
integral_arg = ip.Results.integral_arg;
r_miss = ip.Results.r_miss;
init = ip.Results.init;

if r_miss~=1
    error('*** FIXME *** not implemented yet')
    %%%%%%%%% TODO: take into account prior ratio different from 1
    %%%%%%%%% r_miss = p(miss|z=1)/p(miss|z=0)
end
if modelparam.observe_all && modelparam.necessarily_sparse
    error('observe_all and necessarily_sparse can''t be both true')
end

if strcmp(typegraph, 'simple')
    issimple = true;
else
    issimple = false;
end

%% Check G
if ~ismatrix(G) || ~issparse(G)
    error('G must be a sparse symmetric matrix');
end

% Handle missing values
G_ismiss = sparse(isnan(G));
if ~isequal(G_ismiss, G_ismiss')
    error('Missing values must be a symmetric sparse matrix');
end
nmiss = nnz(G_ismiss);
if nmiss>0
    error('*** FIXME *** not implemented yet');
    [i_miss, j_miss] = find(G_ismiss);
end

G(G_ismiss) = 0;
if ~isequal(G, G')
    error('G must be a sparse symmetric matrix');
end

if strcmp(like_model, 'ber')
    G = logical(G);
else
    error('*** FIXME *** not implemented yet')
end

if ~modelparam.observe_all && any(sum(G, 2)==0)
    error('G must not have empty rows or columns')
end

K = size(G, 1);

if issimple % If no self-loops
    [ind1, ind2] = find(triu(G, 1));
else
    [ind1, ind2] = find(triu(G));
end

nedges = numel(ind1);

%% Number of features
p = modelparam.p;

%% struct array of logicals indicating whether to estimate parameters
S = struct;
estim = struct;
% replace hyperparameters by empty if unobserved
[S.alpha, S.sigma, S.tau, S.Fdist, S.gamma] = CGGPgetparams(modelparam.p, modelparam.alpha,...
    modelparam.sigma, modelparam.tau, modelparam.Fdist, modelparam.gamma, 'empty', ...
    modelparam.observe_all, modelparam.necessarily_sparse);

% estim field is true if S field is empty
vnames = {'alpha', 'sigma', 'tau', 'gamma'};
for v = 1:numel(vnames)
    estim.(vnames{v}) = isempty(S.(vnames{v}));
end
if isnumeric(S.Fdist.param)
    estim.Fparam = isempty(S.Fdist.param);
elseif isstruct(S.Fdist.param)
    Fparnames = fieldnames(S.Fdist.param);
    for v=1:numel(Fparnames)
        estim.Fparam.(Fparnames{v}) =  isempty(S.Fdist.param.(Fparnames{v}));
    end
else
    error('Fdist must be either numeric or struct')
end

estim.w_rem = ~modelparam.observe_all; % If all nodes are observed (including those with zero connection), estimate it

%% Initial values
if ~isempty(fieldnames(init))
    %%% TODO: check fields
    S = init;
else
    [S.alpha, S.sigma, S.tau, S.Fdist, S.gamma] = CGGPgetparams(modelparam.p, modelparam.alpha,...
        modelparam.sigma, modelparam.tau, modelparam.Fdist, modelparam.gamma, 'init', ...
        modelparam.observe_all, modelparam.necessarily_sparse);
    S.logalpha = log(S.alpha);
    S.w_rem = zeros(1, p);
    if estim.w_rem
        S.w_rem = gamrnd(1, 1, [1, p]); 
    end
    S.w0 = gamrnd(1, 1, [K, 1]); 
    S.beta = gamrnd(1/sqrt(p), 1, [K, p]); 
end
S.logw0 = log(S.w0);
S.logbeta = log(S.beta);
S.logw = bsxfun(@plus, S.logw0, S.logbeta);
S.w = exp(S.logw);

if p==1 && strcmp(like_model, 'poiss') && nmiss==0
    S.m = sum(G, 1)/2;
end

if any(S.gamma ~= 0)
    error('*** FIXME *** case gamma_k>0 not supported yet');
end

%% Parameters of the MCMC algorithm
niter = mcmcparam.niter;
nburn = mcmcparam.nburn;
nlatent = mcmcparam.latent.nlatent;
thin = mcmcparam.thin;
L = mcmcparam.leapfrog.L; % Number of leapfrog steps
epsilon = mcmcparam.leapfrog.epsilon/(K*(p+1))^(1/4); % Leapfrog stepsize
rw_std = mcmcparam.hyper.rw_std; % hyperparameters random walk stepsizes
ntotalmass = mcmcparam.totalmass.ntotalmass;

%% To store MCMC samples
nsamples = floor((niter-nburn)/thin);
samples = struct;
if mcmcparam.store_w
    samples.w = NaN(K, p, nsamples, 'double');
else
    samples.w = [];
end
samples.w_rem = NaN(1, p, nsamples);
samples.alpha = NaN(1, 1, nsamples);
samples.logalpha = NaN(1, 1, nsamples);
samples.sigma = NaN(1, 1, nsamples);
samples.tau = NaN(1, 1, nsamples);
samples.gamma = NaN(1, p, nsamples);
if isnumeric(S.Fdist.param)
    samples.Fparam = NaN([size(S.Fdist.param), nsamples]);
else
    Fparnames = fieldnames(S.Fdist.param);
    for v=1:numel(Fparnames)
        samples.Fparam.(Fparnames{v}) = NaN([size(S.Fdist.param.(Fparnames{v})), nsamples]);
    end
end

rate_hmc = zeros(niter, 1);
rate_hyper = zeros(niter, 1);
epsilon_all = zeros(niter, 1);

%% MCMC iterations
for i = 1:niter
    % Update the counts using Gibbs or Metropolis-Hastings step
    if mod(i-1, nlatent)==0
        S.m = update_m_graph(S.logw, ind1, ind2);
    end
    
    % Update w using Hamiltonian Monte Carlo
    [S.w, S.logw, S.w0, S.logw0, S.beta, S.logbeta, rate_hmc(i)] = update_w_graph(...
        S.w, S.logw, S.w0, S.logw0, S.beta, S.logbeta, S.w_rem, S.m, ...
        L, epsilon, S.sigma, S.tau, S.Fdist, issimple);
    
    if i<mcmcparam.leapfrog.nadapt % Adapt the stepsize
        epsilon = exp(log(epsilon) + mcmcparam.leapfrog.adaptrate* ...
            (mean(rate_hmc(max(1,i-mcmcparam.leapfrog.adaptwidth+1):i)) - 0.65));
    elseif i==mcmcparam.leapfrog.nadapt
        epsilon = min(epsilon_all(ceil(i/2):i-1)); 
    end
    epsilon_all(i) = epsilon;
    
    
    if ~modelparam.observe_all % Update hyperparam every ntotalmass iterations
        estim.w_rem = (mod(i-1, ntotalmass)==0);
        estim.alpha = (mod(i-1, ntotalmass)==0);
        estim.tau = (mod(i-1, ntotalmass)==0);
        estim.sigma = (mod(i-1, ntotalmass)==0);
    end        
        
    % Update hyperparameters using MH
    rw_alpha = rem(i,2)==0; % alternate between random walk and gamma proposals for alpha
    
    [S.w_rem, S.alpha, S.logalpha, S.sigma, S.tau, S.Fdist, S.gamma, rate_hyper(i)] = ...
        update_hyper_graph(sum(S.w), S.beta, S.logbeta, S.w0, S.logw0, S.w_rem,...
        S.alpha, S.logalpha, S.sigma, S.tau, S.Fdist, S.gamma, T, estim, ...
        modelparam, rw_std, rw_alpha,...
        mcmcparam.hyper.MH_nb, modelparam.necessarily_sparse, modelparam.observe_all, sumsmall_arg, integral_arg{:});
    
    if i<mcmcparam.hyper.nadapt % Adapt the stepsize
        rw_std = structfun(@(x) exp(log(x) + mcmcparam.hyper.adaptrate*...
            (mean(rate_hyper(max(1,i-mcmcparam.hyper.adaptwidth):i)) - 0.23)), ...
            rw_std, 'uniformoutput', false);
    end
    
    % Print current state
    printind = min(2000, round(niter/5)); % print output every 2000 or less
    if verbose && rem(i, printind)==0
        print_state(i, S, mean(rate_hyper(i-printind+1:i)), mean(rate_hmc(i-printind+1:i)), epsilon, rw_std);
    end
    
    % Store output
    if (i>nburn && rem((i-nburn),thin)==0)
        ind = ((i-nburn)/thin);
        samples = store_samples(ind, S, samples, mcmcparam);
    end
end

stats.rate = rate_hmc;
stats.rate2 = rate_hyper;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Utility subfunctions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function print_state(i, S, rate_hyper, rate_hmc, epsilon, rw_std)
fprintf('i=%i ', i);
fprintf('alp=%.2f ', exp(S.logalpha));
fprintf('sig=%.3f ', S.sigma);
fprintf('tau=%.2f ', S.tau);
switch S.Fdist.name
    case 'gamma'
        if isnumeric(S.Fdist.param)
            fprintf('a=');
            fprintf('%.2f ', S.Fdist.param);
        else
            fprintf('a=');
            fprintf('%.2f ', S.Fdist.param.a);
            fprintf('b=');
            fprintf('%.2f ', S.Fdist.param.b);
        end
    otherwise
        error('Unknown distribution F')
end
fprintf('w*=');
fprintf('%.2f ', S.w_rem);
fprintf('b2=');fprintf('%.2f ', S.Fdist.param.b*S.tau);
fprintf('alp2=%.2f ', exp(S.logalpha + S.sigma*log(S.tau)));
fprintf('rhmc=%.2f ', rate_hmc);
fprintf('rhyp=%.2f ', rate_hyper);
fprintf('eps=%.2g ', epsilon);
fprintf('rwsd=%.2g\n', rw_std.alpha);
end

function samples = store_samples(ind, S, samples, mcmcparam)
if mcmcparam.store_w
    samples.w(:, :, ind) = S.w;
end
samples.w_rem(:, :, ind) = S.w_rem;
samples.alpha(:, :, ind) = S.alpha;
samples.logalpha(:, :, ind) = S.logalpha;
samples.tau(:, :, ind) = S.tau;
samples.sigma(:, :, ind) = S.sigma;
samples.gamma(:, :, ind) = S.gamma;
if isnumeric(S.Fdist.param)
    samples.Fparam(:, :, ind) = S.Fdist.param;
else
    Fparnames = fieldnames(S.Fdist.param);
    for v=1:numel(Fparnames)
        samples.Fparam.(Fparnames{v})(:, :, ind) = S.Fdist.param.(Fparnames{v});
    end
end
samples.logalpha2(:, :, ind) = S.logalpha + S.sigma*log(S.tau);
samples.Fparam.b2(:, :, ind) = S.Fdist.param.b*S.tau;
end

