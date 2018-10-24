function [logpost, loglik] = logpost(objmcmc, G)

% logpost computes the log posterior (up to a constant) 
% for all samples of parameters from the MCMC posterior 
% given the observed adjacency matrix G
% ---------------------------------------------
% INPUTS:
%   - objmcmc: object of class graphmcmc
%   - G: observed binary adjacency matrix
% ---------------------------------------------
% OUTPUT:
%   - logpost: logposterior
%   - loglik: loglikelihood
% -------------------------------------------------------------------------

% Copyright (c) F. Caron (University of Oxford), A. Todeschini (Inria), and 
% X. Miscouridou (University of Oxford)
% caron@stats.ox.ac.uk
% adrien.todeschini@gmail.com
% xenia.miscouridou@spc.ox.ac.uk
% September 2017
%--------------------------------------------------------------------------


% Computes the logposterior (up to a constant)

if ~objmcmc.prior.param(1).observe_all || (numel(objmcmc.prior.param)>1 && ~objmcmc.prior.param(2).observe_all)
    error('Not implemented');
end

samples = arraystruct2structarray(objmcmc.samples);

[~, nchains, nsamples] = size(samples);
logpost = zeros(nsamples*nchains, 1);
loglik = zeros(nsamples*nchains, 1);

for ch=1:nchains
    samples_all = arraystruct2structarray(combine(objmcmc.samples));
    parfor (i=1:nsamples*nchains, getpoolsize())
%     for i=1:nsamples*nchains
        [logpost(i), loglik(i)] = logpostcggp(G, samples_all(:,1,i), objmcmc.prior);
    end
end

logpost = reshape(logpost, [nsamples nchains]);
loglik = reshape(loglik, [nsamples nchains]);