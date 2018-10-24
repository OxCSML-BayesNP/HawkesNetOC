function sam = estimate_kernel_multi(train_data_forw, train_data_backw, BASE_2, precomputed_diff, num_events, modelparameters, mcmcparam)

% estimate_kernel_multi runs an MCMC sampler for posterior inference
% eta,delta.
%
% sam = estimate_kernel_multi(train_data_forw, train_data_backw, BASE_2, precomputed_diff, num_events, modelparameters, mcmcparam)
%
% -------------------------------------------------------------------------
% INPUTS
%   
%   - train_data_forw:  the forward event times of the processes
%   - train_data_backw: the backward event times of the processes
%   - BASE_2:
%   - precomputed_diff: a precomputed term to speed up inference
%   - num_events: 
%   - modelparameters: struct with the model parameters and hyperparameters
%   - mcmcparam: struct with the mcmc hyperparameters and parameters 
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
    times_rec = train_data_backw;
    mu = BASE_2;
  
    objmcmc=struct;
    objmodel=struct;
  
    objmcmc.settings=mcmcparam;
    objmodel.param = modelparameters;
    objmodel.param.base = mu;
    objmodel.param.N_events_per_pair = num_events;

    seed = randi(floor(intmax/12));
  
    fprintf(' Starting the MCMC chains')

    if getpoolsize
        parfor ( k=1:objmcmc.settings.nchains, getpoolsize() )
            rng(seed+k);
           [sam(:,k), sta(:,k)] = HPmcmc_multi_edge(times, times_rec,precomputed_diff, objmodel, objmcmc);
        end
    else
        
        for k=1:objmcmc.settings.nchains
            rng(seed+k);
           [sam(:,k)] = HPmcmc_multi_edge(times, times_rec,precomputed_diff, objmodel, objmcmc);
        end
        
    end
   
end 
