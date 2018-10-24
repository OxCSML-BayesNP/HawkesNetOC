function [estimates, C_st] = graphest(objmcmc, nburn, verbose, cost)
 
% graphest returns point estimates of the graph parameters using an absolute 
% cost function invariant to the features order (label switching)
% [samples_all, estimates] = GRAPHEST(objmcmc, nburn)
%
% -------------------------------------------------------------------------
% INPUTS
%   - objmcmc: an object of the class graphmcmc
% Optional input
%   - nburn: number of MCMC iterations to remove to get estimates (default:0)
%
% OUTPUTS
%   - samples_all: structure containing the MCMC samples for all variables,
%       concatenated over all MCMC chains
%   - estimates: structure containing the median estimates for the
%       different parameter (one field per parameter)
%
% See also GRAPHMCMC, GRAPHMCMC.GRAPHMCMC, GRAPHMCMCSAMPLES
% -------------------------------------------------------------------------
 
% Copyright (c) F. Caron (University of Oxford), A. Todeschini (Inria), and 
% X. Miscouridou (University of Oxford)
% caron@stats.ox.ac.uk
% adrien.todeschini@gmail.com
% xenia.miscouridou@spc.ox.ac.uk
%
% September 2015
%--------------------------------------------------------------------------
 
if nargin <2
    nburn = 0;
end
if nargin<3
    verbose = true;
end
if nargin<4
    cost = 'p3';
end
names = fieldnames(objmcmc.samples);
nsamples = size(objmcmc.samples(1).(names{1}), ndims(objmcmc.samples(1).(names{1})));
nchains = size(objmcmc.samples, 2);
nburn = floor(nburn/objmcmc.settings.thin);
 
if nburn>=nsamples
    error('nburn is larger than the chain length')
end
 
%% get all samples (from several chains)
samples_all = combine(discard(objmcmc.samples, nburn));
nsamples_all = size(samples_all(1).(names{1}), ndims(samples_all(1).(names{1})));
 
%% get estimates
if (isfield(objmcmc.samples, 'w') && ndims(objmcmc.samples(1).w)==3 && size(objmcmc.samples(1).w, 2)>1) 
    % case CGGP with p>1 features
    
    %%% use a label switching invariant cost function
    switch cost
        case 'p4'
            cost_fun = @cost_p4;
        case 'p3'
            cost_fun = @cost_p3;
        case 'p2'
            cost_fun = @cost_p2;
        case 'plogp'
            cost_fun = @cost_plogp;
        otherwise
            error('Unknown cost function')
    end
    
    samples_all = arraystruct2structarray(samples_all);
    
    C_st = zeros(nsamples-nburn, nchains);
    
    if verbose
        fprintf('%s\n', repmat('-',1,35))
        fprintf('Start parameters estimation for CGGP graphs: %d samples\n', nsamples_all)
    end
    tstart = tic;
    
    for i=1:nsamples_all
        %%% Monte-Carlo approximation of posterior cost for sample i
        C = zeros(nsamples_all,1);
        sample_est = samples_all(:,1,i);
        % Parallel computations
        parfor (j=1:nsamples_all, getpoolsize())
            C(j) = cost_fun(sample_est, samples_all(:,1,j));
        end
        C_st(i) = mean(C);
        if verbose
            progress(i, tstart, nsamples_all, 35);
        end
    end
    
    if verbose
        fprintf('|\n')
        fprintf('End parameters estimation for CGGP graphs\n')
        fprintf('Computation time: %.1f hours\n', toc(tstart)/3600);
        fprintf('%s\n', repmat('-',1,35))
    end
    
    [~, i_est] = min(C_st(:));
    
    estimates = samples_all(:,1,i_est);
   
elseif (isfield(objmcmc.samples, 'w') && size(objmcmc.samples(1).w, 2)== 1) 
    %%% take the median
    for t = ntypes:-1:1 %% loop over types of nodes
        for j=1:numel(names)
            name = names{j};
            if isempty(objmcmc.samples(t,1).(name))
                estimates(t,1).(name) = [];
            else
                nd = ndims(samples_all(t,1).(name));
                if isnumeric(objmcmc.samples(t,1).(name))
                    estimates(t,1).(name) = median(samples_all(t,1).(name), nd);
                else
                    fn = fieldnames(objmcmc.samples(t,1).(name));
                    for v=1:numel(fn)
                        estimates(t,1).(name).(fn{v}) = median(samples_all(t,1).(name).(fn{v}), nd);
                    end
                end
            end
        end
    end
end
 
if (isfield(objmcmc.samples, 'pi') && ndims(objmcmc.samples(1).pi)==3 )
    
    cost_fun = @cost_mmsb;   
    samples_all = arraystruct2structarray(samples_all);
    
    C_st = zeros(nsamples-nburn, nchains);
    
    if verbose
        fprintf('%s\n', repmat('-',1,35))
        fprintf('Start parameters estimation for MMSB graph: %d samples\n', nsamples_all)
    end
    tstart = tic;
    
    for i=1:nsamples_all
        
        % Monte-Carlo approximation of posterior cost for sample i
        C = zeros(nsamples_all,1);
        sample_est = samples_all(:,1,i);
        
        % Parallel computations
        parfor (j=1:nsamples_all, getpoolsize())
            C(j) = cost_fun(sample_est, samples_all(:,1,j));
        end
        C_st(i) = mean(C);
        if verbose
            progress(i, tstart, nsamples_all, 35);
        end
    end    
    
    [~, i_est] = min(C_st(:));
    estimates = samples_all(:,1,i_est);
       
    if verbose
        fprintf('|\n')
        fprintf('End parameters estimation for MMSB graphs\n')
        fprintf('Computation time: %.1f hours\n', toc(tstart)/3600);
        fprintf('%s\n', repmat('-',1,35))
    end
   
end
end
 
function progress(i, tstart, niter, s)
if i==5
    est_time = toc(tstart) * niter/i/3600;
    fprintf('Estimated end of computation: %s (%.1f hours)\n', datestr(now + toc(tstart) * (niter-i)/i/3600/24), est_time);
    fprintf('|%s|\n|', repmat('-',1,(min(niter,s-2))))
    fprintf(repmat('*', 1, floor(i/ceil(niter/min(niter,s-2)))))
end
if i>5 && mod(i, ceil(niter/min(niter,s)))==0
    fprintf('*')
end
end
 
 
%% cost function
 
function c = cost_p4(param_est, param_true)
 
% NOT GOOD TAKES TOO MUCH TIME
c = sum(sum( (param_est(1).w*param_est(2).w'- param_true(1).w*param_true(2).w').^2 ))...
    + ( (sum(sum(param_est(1).w, 1) + param_est(1).w_rem))*(sum(sum(param_est(2).w, 1) + param_est(2).w_rem))...
    - (sum(sum(param_true(1).w, 1) + param_true(1).w_rem))*(sum(sum(param_true(2).w, 1) + param_true(2).w_rem)) ).^2;
 
end
 
 
% the cost is minimized over permutations of features
function c = cost_p3(param_est, param_true)
 
if numel(param_est) > 1 % RENORMALIZE
    bias1 = sum(param_est(1).w, 1) + param_est(1).w_rem;
    bias2 = sum(param_est(2).w, 1) + param_est(2).w_rem;
    param_est(2).w = bsxfun(@times, param_est(2).w, bias1);
    param_est(2).w_rem = param_est(2).w_rem.*bias1;
    param_est(1).w = bsxfun(@times, param_est(1).w, bias2);
    param_est(1).w_rem = param_est(1).w_rem.*bias2;
    
    bias1 = sum(param_true(1).w, 1) + param_true(1).w_rem;
    bias2 = sum(param_true(2).w, 1) + param_true(2).w_rem;
    param_true(2).w = bsxfun(@times, param_true(2).w, bias1);
    param_true(2).w_rem = param_true(2).w_rem.*bias1;
    param_true(1).w = bsxfun(@times, param_true(1).w, bias2);
    param_true(1).w_rem = param_true(1).w_rem.*bias2;
end
 
p = size(param_true(1).w, 2);
 
%%% matrix of costs for all permutations of features
%%% C(i,j) = C(param_est{i}, param_true{j})
%%% where {k} indicates feature k
C = zeros(p);
for k=1:p
    for t=1:numel(param_est) %% loop over types of nodes
        C(k,:) = C(k,:) + sum(abs(bsxfun(@minus, param_est(t).w(:,k), param_true(t).w)))...
            + abs(param_est(t).w_rem(k)-param_true(t).w_rem);
    end
end
%%% Hungarian algorithm: finds the permutation with minimal cost
%%% c is the minimal cost for each row
[~, c] = munkres(C);
 
end
 
 
%% cost function
% the cost is minimized over permutations of features
function c = cost_p2(param_est, param_true)
 
if numel(param_est) > 1 % RENORMALIZE
    bias1 = sum(param_est(1).w, 1) + param_est(1).w_rem;
    bias2 = sum(param_est(2).w, 1) + param_est(2).w_rem;
    param_est(2).w = bsxfun(@times, param_est(2).w, bias1);
    param_est(2).w_rem = param_est(2).w_rem.*bias1;
    param_est(1).w = bsxfun(@times, param_est(1).w, bias2);
    param_est(1).w_rem = param_est(1).w_rem.*bias2;
    
    bias1 = sum(param_true(1).w, 1) + param_true(1).w_rem;
    bias2 = sum(param_true(2).w, 1) + param_true(2).w_rem;
    param_true(2).w = bsxfun(@times, param_true(2).w, bias1);
    param_true(2).w_rem = param_true(2).w_rem.*bias1;
    param_true(1).w = bsxfun(@times, param_true(1).w, bias2);
    param_true(1).w_rem = param_true(1).w_rem.*bias2;
end
 
p = size(param_true(1).w, 2);
 
I = true(1,p); % index of unassigned features
 
c = 0;
for k=1:p
    %%% for each true feature assign the estimated feature that minimize
    %%%
    C = 0;
    for t=1:numel(param_est) %% loop over types of nodes
        C = C + sum(abs(bsxfun(@minus, param_est(t).w(:,I), param_true(t).w(:,k))))...
            + abs(param_est(t).w_rem(I) - param_true(t).w_rem(k));
    end
    [cmin, imin] = min(C);
    c = c + cmin;
    ok = I(I);
    ok(imin) = false;
    I(I) = ok;
end
 
end
 
%% Proxy of cost function
% simplifies the permutation problem by ordering the features by descending
% overall weight
function c = cost_plogp(param_est, param_true)
 
if numel(param_est) > 1 % RENORMALIZE
    bias1 = sum(param_est(1).w, 1) + param_est(1).w_rem;
    bias2 = sum(param_est(2).w, 1) + param_est(2).w_rem;
    param_est(2).w = bsxfun(@times, param_est(2).w, bias1);
    param_est(2).w_rem = param_est(2).w_rem.*bias1;
    param_est(1).w = bsxfun(@times, param_est(1).w, bias2);
    param_est(1).w_rem = param_est(1).w_rem.*bias2;
    
    bias1 = sum(param_true(1).w, 1) + param_true(1).w_rem;
    bias2 = sum(param_true(2).w, 1) + param_true(2).w_rem;
    param_true(2).w = bsxfun(@times, param_true(2).w, bias1);
    param_true(2).w_rem = param_true(2).w_rem.*bias1;
    param_true(1).w = bsxfun(@times, param_true(1).w, bias2);
    param_true(1).w_rem = param_true(1).w_rem.*bias2;
end
 
%%% overall weight of each feature
S_est = 1;
S_true = 1;
for t=1:numel(param_est) %% loop over types of nodes
    S_est = S_est .* (sum(param_est(t).w, 1) + param_est(1).w_rem);
    S_true = S_est .* (sum(param_true(t).w, 1) + param_true(1).w_rem);
end
 
%%% sort features by descending overall weight
[~,ind_est] = sort(S_est, 'descend');
[~,ind_true] = sort(S_true, 'descend');
 
%%% compute cost on permuted features
c = 0;
for t=1:numel(param_est) %% loop over types of nodes
    c = c + sum(sum(abs(param_est(t).w(:,ind_est) - param_true(t).w(:,ind_true))))...
        + sum(abs(param_est(t).w_rem(ind_est) - param_true(t).w_rem(ind_true)));
end
end
 
function c = cost_mmsb(param_est, param_true)
if numel(param_est) > 1 % RENORMALIZE
    bias1 = sum(param_est(1).w, 1) + param_est(1).w_rem;
    bias2 = sum(param_est(2).w, 1) + param_est(2).w_rem;
    param_est(2).w = bsxfun(@times, param_est(2).w, bias1);
    param_est(2).w_rem = param_est(2).w_rem.*bias1;
    param_est(1).w = bsxfun(@times, param_est(1).w, bias2);
    param_est(1).w_rem = param_est(1).w_rem.*bias2;
    
    bias1 = sum(param_true(1).w, 1) + param_true(1).w_rem;
    bias2 = sum(param_true(2).w, 1) + param_true(2).w_rem;
    param_true(2).w = bsxfun(@times, param_true(2).w, bias1);
    param_true(2).w_rem = param_true(2).w_rem.*bias1;
    param_true(1).w = bsxfun(@times, param_true(1).w, bias2);
    param_true(1).w_rem = param_true(1).w_rem.*bias2;
end
 
p = size(param_true(1).pi, 2);
 
%%% matrix of costs for all permutations of features
%%% C(i,j) = C(param_est{i}, param_true{j})
%%% where {k} indicates feature k
C = zeros(p);
for k=1:p
    for t=1:numel(param_est) %% loop over types of nodes
        C(k,:) = C(k,:) + sum(abs(bsxfun(@minus, param_est(t).pi(:,k), param_true(t).pi)));
    end
end
%%% Hungarian algorithm: finds the permutation with minimal cost
%%% c is the minimal cost for each row
[~, c] = munkres(C);
end