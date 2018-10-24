function S_multi = calculate_S_multi(delta, eta, times, times_rec, precomputed_diff, N_e, max_N, num_events,num_events_rec )
%
% calculate_S_multi calculates recursively a term for the intensity of all
% the pairs of processes
%
% -------------------------------------------------------------------------
% INPUTS
%   
%   - delta: current estimate of delta
%   - eta: current estimate of delta
%   - times: forward event times for the processes
%   - times_rec: backward event times for the processes
%   - precomputed_diff: a precomputed term to speed up calculations
%   - N_e: number of pairs of processes
%   - max_N: maximum number of events among all processes
%   - num_events: num of forward events for every process
%   - num_events_rec: num of backward events for everyprocess
% 
% OUTPUTS
%   - S_multi: value of this term of the intensity
%  
% -------------------------------------------------------------------------
% Copyright (C) Xenia Miscouridou, University of Oxford
% xenia.miscouridou@spc.ox.ac.uk
% October 2018
%--------------------------------------------------------------------------


    S_multi = zeros(N_e, max_N);
    for edge=1:N_e
        S_multi(edge, 1:num_events(edge)) = calculate_S(delta, eta, times(edge,1:num_events(edge)), times_rec(edge,1:num_events_rec(edge)),precomputed_diff(edge), num_events(edge));

    end
        
end