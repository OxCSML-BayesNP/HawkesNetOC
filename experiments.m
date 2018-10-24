% Experiments
% college dataset
% load the dataset, preprocess data and perform the two stage inference 
% 1) to estimate the base intensity parameters muhat 
% 2) to estimate the kernel intensity parameters conditionally on muhat
%
%
%

clear all
close all
collegedata = true;
istest = 1;
addpath ./Hawkes ./CGGP ./GGP ./utils ./

if  collegedata
    % college
    addpath ./data 
    load data/college/college
    indices = randperm(length(meta.source),round(length(meta.source)));
    source = meta.source(indices);target = meta.target(indices);times = meta.times(indices);    
    max_T = max(times);T=max_T;max_T=1;
    times = times/max_T;
    p = 2;
end

%% Posterior Inference stage 1
%
%

% preprocess
estimate_base = 1;
[times_matrix_forw,times_matrix_backw,G,BASE, pairs] = preprocess(source,target,times,p,istest,estimate_base);

%split train - test
[train_ind,times_matrix_forw_2,times_matrix_backw_2,train_data_forw,train_data_backw,test_data_forw,test_data_backw,N_train_links_forw, N_train_links_backw,N_test_links_forw,N_test_links_backw, T_train,max_T] =   split_train_test(times_matrix_forw, times_matrix_backw, 0.6);

%proportion of training links to test links
sum(N_train_links_forw+ N_train_links_backw)./sum(N_test_links_forw+ N_train_links_forw + N_test_links_backw+ N_train_links_backw)

% prediction links estimate from CCRM
BASE_2 = BASE(train_ind,:);
pairs_2 = pairs(train_ind,:);
num_events = [pairs_2(:,3),pairs_2(:,4)];
t1 = T_train;
t2 = max_T;
mean_N = BASE_2*(t2-t1)*T;
sqrt(mean((N_test_links_backw+N_test_links_forw-2*mean_N).^2))


%% Posterior Inference stage 2
%
%
estimate_kernel = true ;
if estimate_kernel
    %% addpath
    keep_all=1;
    if istest
        niter = 100;
        nsamples = niter; % Nb of Monte Carlo samples to return
        nchains = 1;
        nburn = floor(1*niter/3); 
        thin = ceil((niter-nburn)/nsamples);

    else
        niter = 1000;
        if keep_all
            nsamples = niter;
        else
            nsamples = 500;
        end
        nchains = 1;
        nburn = floor(1*niter/3); 
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
    mcmcparam.hyper.prior_delta.param = .01;
    mcmcparam.hyper.prior_eta.param = 0.01;

    modelparameters.max_T = max_T;
    modelparameters.p = 1;
    
    N_e = size(train_data_forw,1);
    % precomputations
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
    precomputed_diff(isnan(precomputed_diff)) = 1;
       
    [samples] = estimate_kernel_multi(train_data_forw, train_data_backw,BASE_2,precomputed_diff,num_events, modelparameters, mcmcparam);
    params = [];
    params.nchains=nchains;
    params.niter=niter;
    params.nsamples =nsamples;

    %%%% plots and prediction
    plot_figures(samples, T, params)
    eta_estimate = mean(samples(1).eta);
    delta_estimate = mean(samples(1).delta);
    t1 = T_train;
    t2 = max_T;
    mean_N = predict(BASE_2*T,eta_estimate,delta_estimate,t1,t2);
    RMSE = sqrt(mean((N_test_links_forw+N_test_links_backw-2*mean_N).^2));


end

