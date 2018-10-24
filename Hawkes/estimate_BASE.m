function estimates = estimate_BASE(G, p, nchains, istest)
% estimate_BASE estimates the base intensities for the Hawkes Proceses
% -------------------------------------------------------------------------
% INPUTS
%   - G: sparse logical matrix, the corresponding graph
%   - p: positive integer. number of features
%   - nchains: num of chains
%
% OUTPUTS
%   - estimates: the estimated base intensities 
% Copyright (c) F. Caron (University of Oxford), A. Todeschini (Inria), and 
% X. Miscouridou (University of Oxford)
% caron@stats.ox.ac.uk
% adrien.todeschini@gmail.com
% xenia.miscouridou@spc.ox.ac.uk
% September 2017
%--------------------------------------------------------------------------
%
% Reference
% A. Todeschini, X. Miscouridou, F. Caron,
% Exchangeable Random Measures for Sparse and Modular Graphs with Overlapping Communities
% arXiv:1602.02114, 2016.

    addpath ./CGGP ./GGP .//utils
    
    %% Prior distribution 
    objprior =  graphmodel('CGGP', p);

    %% Posterior inference
    %
    % Parameters of the MCMC algorithm
    if istest
        niterinit = 1000;
        niter = 500;
        nsamples = 100; % Nb of Monte Carlo samples to return
    else
        niterinit = 10000;
        niter = 10000;
        nsamples = 500;
    end   
    nburn = floor(3*niter/4); 
    thin = ceil((niter-nburn)/nsamples);
    verbose = true;

    % Create the graphMCMC object
    objmcmc = graphmcmc(objprior, niter, 0, thin, nchains); 
    % Note: nburn is set to zero here in order to store samples in the transient regime of the MCMC

    % Run initialisation
    init = graphinit(objmcmc, G, niterinit);

    % Run MCMC sampler
    objmcmc = graphmcmcsamples(objmcmc, G, verbose, init);

    % discard burnin to compute estimates
    objmcmc_all = objmcmc;
    objmcmc.samples = discard(objmcmc_all.samples, floor(nburn/objmcmc_all.settings.thin));
    objmcmc.settings.nburn = nburn;

    % Get estimates and cost
    [estimates, ~] = graphest(objmcmc);
end
