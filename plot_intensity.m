function [lambda,lambda2] = plot_intensity(mu, eta, delta, T_1, T_2)
% plot_intensity plots the intensity for the processes in a pair of mutually exciting processes
% ie plots
% (lambda_{ij}, lambda_{ji}), where
% lambda_{ij}(t) = mu_{ij} + eta sum_{T_{ji}< t}  exp(-delta*(t-T_{ji}) ) 
% lambda_{ji}(t) = mu_{ji} + eta sum_{T_{ij}< t}  exp(-delta*(t-T_{ij}) )
%  mu_{ji} = mu_{ij}
% -------------------------------------------------------------------------
% INPUTS
%   
%   - mu: the base intensity of the process (it is shared within the pair)
%   - eta: Hawkes kernel parameter, the step size 
%   - delta: Hawkes kernel parameter, the exponential decay 
%   - T_1: the forward event times for the process
%   - T_2: the backward event times for the process
% 
% OUTPUTS
%
%   Returns the intensity plot.
%  
% -------------------------------------------------------------------------
% Copyright (C) Xenia Miscouridou, University of Oxford
% xenia.miscouridou@spc.ox.ac.uk
% October 2018
%--------------------------------------------------------------------------
    figure;
    max_T = max(max(T_1,T_2));
    t=0:.0001:max_T;
    timepoints = sort(unique(sort([T_1;t'])));
    t_back = T_2;
    lambda = zeros*length(timepoints);
    for j = 1:length(timepoints)
        s =(timepoints(j));
        lambda(j) = mu + eta*sum(exp(-delta*(s-t_back(t_back<s))));
    end
    
    timepoints = sort(unique(sort([T_2;t'])));
    lambda2=zeros*length(timepoints);
    t_back = T_1;
    for j = 1:length(timepoints)
        s =(timepoints(j));
        lambda2(j) = mu + eta*sum(exp(-delta*(s-t_back(t_back<s))));
    end
    
    subplot(2,1,1)
    plot(lambda,'r')
    ylim([0,0.5+max(lambda)]);
    legend('lambda_{12}')
    %xlim([0,max_T])
    hold all;
     subplot(2,1,2)
    plot(lambda2,'g')
    ylim([0,0.5+max(lambda2)]);
    legend('lambda_{21}')


end