function [T_1, T_2, lambda_1_plus, lambda_2_plus, N_1, N_2] = simulate_bivariate_process(K, mu, eta, delta)
% simulate_bivariate_process simulates a bivariate hawkes process with exponential kernel,
% i.e. a pair of mutually exciting Hawkes processes with intensities 
% (lambda_{ij}, lambda_{ji}), where
% lambda_{ij}(t) = mu_{ij} + eta sum_{T_{ji}< t}  exp(-delta*(t-T_{ji}) ) 
% lambda_{ji}(t) = mu_{ji} + eta sum_{T_{ij}< t}  exp(-delta*(t-T_{ij}) )
%  mu_{ji} = mu_{ij}
% the correspodning integer valued counting processes are denoted by 
%  N_{ij} and N_{ji}
% -------------------------------------------------------------------------
% INPUTS
%   
%   - K:  the number of event times to simulate
%   - mu: the base intensity of the process (it is shared within the pair mu_{ji}=mu_{ij} )
%   - eta: Hawkes kernel parameter, the step size 
%   - delta: Hawkes kernel parameter, the exponential decay 
% 
% OUTPUTS
%
%   - T_1: the forward event times for the process
%   - T_2: the backward event times for the process
%   - lambda_1_plus: the intensity function evaluated at the event timepoints from the right 
%   - lambda_2_plus: the intensity function evaluated at the event timepoints from the left 
%   - N_1: the counting process for T_1 
%   - N_2: the counting process for T_2
%  
% -------------------------------------------------------------------------
% Copyright (C) Xenia Miscouridou, University of Oxford
% xenia.miscouridou@spc.ox.ac.uk
% October 2018
%--------------------------------------------------------------------------
%
% Reference
% A. Dassios and H. Zhao. Exact simulation of hawkes process with exponentially decaying intensity.
% Electronic Communications in Probability, 18(62):1-13, 2013.

    T_1=zeros(K,1);
    T_2=zeros(K,1);
    lambda_1_plus=zeros(K,1);lambda_1_minus=lambda_1_plus;
    lambda_2_plus=zeros(K,1);lambda_2_minus = lambda_2_plus;
    N_1=zeros(K,1);
    N_2=zeros(K,1);
    S_1=zeros(K,1);
    S_2=zeros(K,1);
    D_1=zeros(K,1);
    D_2=zeros(K,1);


    % k=0
    k=0;
    T_1(k+1)=0;T_2(k+1)=0;
    lambda_1_plus(k+1)=  mu;
    lambda_2_plus(k+1) = mu;
    lambda_1_minus(k+1) = mu;
    lambda_2_minus(k+1) =mu;

    N_1(k+1) =0 ;
    N_2(k+1) =0;

    %simulation
    for k = 1:(K-1)

        %increment for process 2 that needs process 1 info
        D_1(k+1) = 1 + delta*log(rand)/(lambda_1_plus(k) - mu); 
        b = -log(rand)/mu;
        if  D_1(k+1)<0
            S_1(k+1) = b;
        else
            a=-log(D_1(k+1))/delta;
            S_1(k+1) = min(a,b);
        end


        %process 2
        T_2(k+1) = T_2(k)+S_1(k+1);
        N_2(k+1) = N_2(k)+1;
        lambda_2_minus(k+1) = (lambda_2_plus(k) - mu)*exp(-delta*(S_1(k+1))) + mu;
        lambda_2_plus(k+1)=lambda_2_minus(k+1)+eta;

        %increment for process 1 that needs process 2 info
        D_2(k+1) = 1 + delta*log(rand)/(lambda_2_plus(k) -mu);
        b = -log(rand)/mu;

        if  D_2(k+1) <0
            S_2(k+1) = b;
        else
            if D_2(k+1) >0
             a=-log(D_2(k+1))/delta;
             S_2(k+1) = min(a,b);
            end
        end

        %process 1
        T_1(k+1) = T_1(k)+S_2(k+1);
        N_1(k+1) = N_1(k)+1;
        lambda_1_minus(k+1) = (lambda_1_plus(k) - mu)*exp(-delta*( S_2(k+1)) ) + mu;
        lambda_1_plus(k+1)=lambda_1_minus(k+1)+eta;

    end
   
end
