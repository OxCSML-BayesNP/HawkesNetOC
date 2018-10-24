function plot_figures(samples, max_T, params)
% plot_figures plots the parameter estimates
% -------------------------------------------------------------------------
% INPUTS
%   
%   - samples: struct with the kernel parameter posterior estimates
%   - max_T: the maximum time for observed events for all processes
%   - params: model parameters
% 
% OUTPUTS
%
%   Returns the plot.
%  
% -------------------------------------------------------------------------
% Copyright (C) Xenia Miscouridou, University of Oxford
% xenia.miscouridou@spc.ox.ac.uk
% October 2018
%--------------------------------------------------------------------------


    niter = params.niter;
    nburn = floor(2*niter/3); 
    nsamples=params.nsamples;
    thin = ceil((niter-nburn)/nsamples);

    %get nchains
    figure;
    subplot(2,1,1);   
    for i =1:params.nchains       
        plot(thin:thin:(niter-nburn), squeeze(samples(i).eta(thin:thin:(niter-nburn)))/squeeze(samples(i).delta(thin:thin:(niter-nburn))));
        ylabel('eta/delta')
      %  xaxis
        hold all;
    end
    subplot(2,1,2);
    for i =1:params.nchains
        plot(thin:thin:(niter-nburn), squeeze(samples(i).eta(thin:thin:(niter-nburn))));
        ylabel('eta')
        hold all;
    end
    
    figure;
    subplot(2,1,1);
    for i =1:params.nchains
        plot(thin:thin:(niter-nburn),squeeze(samples(i).delta(thin:thin:(niter-nburn))));
        hold on;
       ylabel('delta')
        hold all;
    end
        subplot(2,1,2);
    for i =1:params.nchains
        plot(thin:thin:(niter-nburn), squeeze(samples(i).eta(thin:thin:(niter-nburn)).*samples(i).delta(thin:thin:(niter-nburn))/max_T));
        hold on;
        ylabel('eta/delta')
        hold all;
    end
    


end