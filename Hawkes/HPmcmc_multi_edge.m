function [samples] = HPmcmc_multi_edge(train_data_forw,  train_data_backw, precomputed_diff,  objmodel, objmcmc)
% HPmcmc_multi_edge runs an MCMC sampler for the mutually exciting Hawkes processes
% conditionally on the value of the base intensity
% It performs Gibbs sampling with MH steps for the two parameters
% eta,delta.
%
% [ samples ] = HPmcmc_multi_edge(train_data_forw,  train_data_backw, precomputed_diff,  objmodel, objmcmc)
%
% -------------------------------------------------------------------------
% INPUTS
%   
%   - train_data_forw:  the forward event times of the processes
%   - train_data_backw: the backward event times of the processes
%   - precomputed_diff: a precomputed term to speed up inference
%   - objmodel: struct with the model parameters and hyperparameters
%   - objmcmc: struct with the mcmc hyperparameters and parameters 
% 
% OUTPUTS
%   - samples: structure array of length 2, corresponding to the parameters
%   eta and delta
%  
% -------------------------------------------------------------------------
% Copyright (C) Xenia Miscouridou, University of Oxford
% xenia.miscouridou@spc.ox.ac.uk
% October 2018
%--------------------------------------------------------------------------


    times = train_data_forw;
    times_rec = train_data_forw;

    samples = [];

    mcmcparam = objmcmc.settings;
    modelparameters = objmodel.param;

    niter = mcmcparam.niter;
    nburn = mcmcparam.nburn;
    thin = mcmcparam.thin;
    nmh = mcmcparam.hyper.nmh;

    prop_eta.name = mcmcparam.hyper.prop_eta.name;
    prop_eta.param =mcmcparam.hyper.prop_eta.param;

    prop_delta.name = mcmcparam.hyper.prop_delta.name;
    prop_delta.param = mcmcparam.hyper.prop_delta.param;
    prop_delta.max = mcmcparam.hyper.prop_delta.max;
    prior_eta.param = mcmcparam.hyper.prior_eta.param;

    prior_delta.param = mcmcparam.hyper.prior_delta.param;

    max_T = modelparameters.max_T;
    N_e = size(times, 1);   % number of processes
    max_N = size(times, 2); % maximum number of events among processes
   
    mu = modelparameters.base;


    num_events = modelparameters.N_events_per_pair(:,1);
    num_events_rec = modelparameters.N_events_per_pair(:,2);

    all_times = unique([train_data_forw;train_data_backw]);   
    all_times = all_times(~isnan(all_times));

    % Initialize output
    nsamples = floor((niter - nburn)/thin);
    samples.eta = zeros( nsamples,1);
    samples.delta = zeros(nsamples,1);
   
    T = max_T;
   
    %% initialize the variables eta,delta
    delta = rand()*100;
    eta = rand()*delta; 
    %precomputations
    logR = -log(delta) + log(  sum(sum( 1 - exp(-delta.*(T-all_times)))));     
    S_multi = calculate_S_multi(delta, eta, times, times_rec, precomputed_diff, N_e, max_N, num_events,num_events_rec );   
    intensities = bsxfun(@plus, S_multi.*eta,mu);

    verbose=1;
    for iter=1:niter  
        %tic       
        if verbose
              if rem(iter, 20)==0
                fprintf('\n');
                fprintf('iter=%d/%d\n', iter, niter);
                fprintf('\n');
                fprintf('eta = %f, delta = %f', eta, delta)
               fprintf('\n');
               fprintf('\n');
              end
        end
   
        
         % Update eta
        [eta, intensities] =  update_eta_multi(eta, delta, mu, logR, S_multi, intensities, prop_eta, prior_eta, nmh);
        
        % Update delta
        [delta, logR, S_multi, intensities] = update_delta_multi(delta, max_T, logR, S_multi, precomputed_diff,intensities, eta, mu, times, times_rec, all_times, N_e, max_N,num_events,  num_events_rec, prop_delta, prior_delta, nmh);

       
        % Save estimates
        if (iter > nburn && rem((iter - nburn), thin) == 0)
                ind = ((iter-nburn)/thin);
                samples.eta(ind) = eta;
                samples.delta(ind) = delta;
        end
        

    end

    % time = toc;
    % fprintf('Time: %s', time);
end