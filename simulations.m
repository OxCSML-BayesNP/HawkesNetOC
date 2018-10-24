% Simulate, visualise and perform inference on a bivariate process i.e. 
% a pair of mutually exciting Hawkes processes that share the base intensity 

% simulate_bivariate_process
addpath ./CGGP ./utils ./GGP ./Hawkes ./

clear all;
close all;

%% Simulate Data
addpath ./
K = 50; % num events
eta = .5; % jump size  [some references denote this by eta\delta]
delta = 1;% exponential decay
mu = 0.15; % base intensity (shared among the bivariate process 
[T_1, T_2, lambda_1_plus, lambda_2_plus, N_1, N_2] = simulate_bivariate_process(K, mu, eta, delta);
[T_3, T_4, lambda_3_plus, lambda_4_plus, N_3, N_4] = simulate_bivariate_process(K/2, mu/2, eta, delta);

T_3=[T_3;NaN(K/2,1)];
T_4=[T_4;NaN(K/2,1)];

BASE = [mu*ones(2,1);mu/2*ones(2,1)];
max_T = max([T_1;T_2;T_3;T_4]);
train_data_forw=   [T_1';T_2';T_3';T_4'];
train_data_backw = [T_2';T_1';T_4';T_3'];

% Visualize the process
    % counting process
plot_counting_process(T_1,T_2,max_T);
    % intensity process
plot_intensity(mu,eta,delta,T_1,T_2);


%% Inference
%% addpath
addpath ./CGGP ./utils ./GGP ./Hawkes ./
estimate_kernel = true;
if estimate_kernel
    istest = 0; keep_all=1;
    if istest
        niter = 500;
        nsamples = niter/2; % Nb of Monte Carlo samples to return
        nchains = 1;
        nburn = floor(1*niter/3); 
        thin = ceil((niter-nburn)/nsamples);

    else
        niter = 8000;
        if keep_all
            nsamples = niter;
        else
            nsamples = 500;
        end
        nchains = 1;
        nburn = floor(2*niter/3); 
        thin = ceil((niter-nburn)/nsamples);
    end   

    mcmcparam.niter = niter;
    mcmcparam.nburn = nburn;
    mcmcparam.thin = thin;
    mcmcparam.nchains =nchains;
    mcmcparam.hyper.nmh = 10;

    %proposals
    prop_eta.name = 'Normal';
    prop_eta.param = 1.5;
    mcmcparam.hyper.prop_eta = prop_eta;

    prop_delta.name = 'Normal';
    prop_delta.param = 2.5;
    mcmcparam.hyper.prop_delta = prop_delta;

    mcmcparam.hyper.prop_delta.max = 2000;

    %priors
    mcmcparam.hyper.prior_delta.param = 0.01;
    mcmcparam.hyper.prior_eta.param = 0.01;

    modelparameters.max_T = max_T;
    modelparameters.p = 1;
    
    num_events= [[K;K;K/2;K/2],[K;K;K/2;K/2]];
     N_e = size(train_data_forw,1);
        precomputed_diff = zeros(N_e,1) ;   
        for edge=1:N_e;
            l=1; 
            indices = (train_data_backw(edge,:)<train_data_forw(edge,l));
            first_term = -(train_data_forw(edge, l) - sum(train_data_backw(edge,indices))) ;


            for l = 2:num_events(edge,1)  
                l;
                indices = (train_data_backw(edge,:)<train_data_forw(edge,l))&(train_data_backw(edge,:)>=train_data_forw(edge,l-1));
                if ~sum(indices)
                    precomputed_diff(edge,1 )  =  exp(-(train_data_forw(edge, l) - sum(train_data_backw(edge,indices))+first_term)) ;
                end

            end
        end
        precomputed_diff(isnan(precomputed_diff))=1;

    [samples] =    estimate_kernel_multi(train_data_forw, train_data_backw, BASE, precomputed_diff,num_events,modelparameters,mcmcparam);

    params = [];
    params.nchains=nchains;
    params.niter=niter;
    params.nsamples =nsamples;
    figure;plot(samples(1).delta);
    hold all;
    plot(delta*ones(length(samples(1).delta),1))
    legend('Estimate','True')
    figure;plot(samples(1).eta);
    hold all;
    plot(eta*ones(length(samples(1).eta),1))
    legend('eta')
end



 


