function mean_N = predict(mu, eta, delta, t1, t2)
%  predict calculates the mean number of events per process
% -------------------------------------------------------------------------
% INPUTS
%   
%   - mu:  Hawkes base intensity
%   - eta: Hawkes kernel parameter for the step size
%   - delta: Hawkes kernel parameter for the exponential decay
%   - t1: forward event times of the process
%   - t2: backward event times of the process
% OUTPUTS
%   - mean_N: mean number of events
% 
% -------------------------------------------------------------------------
% Copyright (C) Xenia Miscouridou, University of Oxford
% xenia.miscouridou@spc.ox.ac.uk
% October 2018
%--------------------------------------------------------------------------
% Reference:
% A. Dassios and H. Zhao. Exact simulation of hawkes process with exponentially decaying intensity. 
% Electronic Communications in Probability, 18(62):1-13, 2013.


    k = delta-eta;
    mean_N = delta*(t2-t1)*(1/k) + (1-delta/k)/k*(exp(-k*t2) - exp(-k*t2));
    mean_N = mean_N.*mu;
end