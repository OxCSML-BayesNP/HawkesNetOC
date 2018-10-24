function [delta, logR, S_multi, intensities] = update_delta_multi(delta, max_T, logR, S_multi, precomputed_diff, intensities, eta, mu, times, times_rec, all_times, N_e,max_N, num_events, num_events_rec, prop_delta, prior_delta, nmh)
%
% update_delta_multi implements the MH step for delta 
%
%
% -------------------------------------------------------------------------
% INPUTS
%   
%   - delta: current estimate of delta
%   - max_T: the right limit for all event times
%   - logR: term involved in the intensity
%   - S_multi: term involved in the intensity
%   - precomputed_diff: a precomputed term to speed up calculations
%   - intensities: intensities evaluated at the current estimates of eta
%   and delta
%   - eta: the current estimate of eta
%   - mu: the base intensities
%   - times: forward event times for the processes
%   - times_rec: backward event times for the processes
%   - all_times: all event times
%   - N_e: number of pairs of processes
%   - max_N: maximum number of events among all processes
%   - num_events: num of forward events for every process
%   - num_events_rec: num of backward events for every process
%   - prop_delta: struct with the name and parameters for the proposal
%   distribution for delta for the MH step
%   - prior_delta: struct with the name and parameters for the prior
%   distribution for delta
%   - nmh: number of MH steps
% 
% OUTPUTS
%   - delta:  the new value of delta
%   - logR: new value of this term of the intensity
%   - S_multi: new value of this term of the intensity
%   - intensities: intensities evaluated at the new estimates of eta
%   and delta
% -------------------------------------------------------------------------
% Copyright (C) Xenia Miscouridou, University of Oxford
% xenia.miscouridou@spc.ox.ac.uk
% October 2018
%--------------------------------------------------------------------------

for nitermh=1:nmh
    %PROPOSAL  
    if strcmp(prop_delta.name,'Exponential')
        deltanew = exprnd(prop_delta.param);
        while (deltanew <eta) || (deltanew > prop_delta.max)
            deltanew = exprnd(prop_delta.param);
        end;
        proposal_term = - prop_delta.param*(deltanew - delta);

    end

    if strcmp(prop_delta.name,'Normal');
        deltanew = ( randn*prop_delta.param + delta );            
        while ((deltanew < eta)||(deltanew > prop_delta.max)) 
           deltanew = ( randn*prop_delta.param + delta); 
        end;
        proposal_term=0;
    end

    %prior
    prior_term = -prior_delta.param*(deltanew - delta);

    logR_new = -log(deltanew) + log(  sum(sum( 1 - exp(-deltanew.*(max_T-all_times)))));
    S_multi_new = calculate_S_multi(deltanew, eta, times, times_rec,precomputed_diff , N_e, max_N, num_events,num_events_rec );
    intensities_new = bsxfun(@plus, S_multi_new.*eta,mu);
    term1 = sum(sum(log(intensities_new) -log(intensities)));  
    term2 = - eta*(exp(logR_new) - exp(logR));

    % MH step
    if log(rand()) <  term1 + term2 + proposal_term + prior_term
        delta = deltanew;
        logR = logR_new;
        S_multi = S_multi_new;
        intensities = intensities_new;
    end
end

end
