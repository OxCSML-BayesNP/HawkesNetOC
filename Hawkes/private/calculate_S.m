function S = calculate_S(delta, eta, times,times_rec, precomp_diff, N_forw)  
%
% calculate_S calculates recursively a term for the intensity of a single
% pair of processes
%
% -------------------------------------------------------------------------
% INPUTS
%   
%   - delta: current estimate of delta
%   - eta: current estimate of delta
%   - times: event times for forward processes
%   - times_rec: event times for backward processes
%   - precomputed_diff: a precomputed term to speed up calculations
%   - N_forw: number of forward events for that process
% 
% OUTPUTS
%   - S : value of this term of the intensity
%  
% -------------------------------------------------------------------------
% Copyright (C) Xenia Miscouridou, University of Oxford
% xenia.miscouridou@spc.ox.ac.uk
% October 2018
%--------------------------------------------------------------------------


   if N_forw>0
    Z = zeros(1,N_forw);
    Z(1) =  exp( -delta*(times(1)))*sum( exp(- times_rec(times_rec<times(1)) )) ;
    if N_forw > 1
        for l = 2:N_forw          
          rec_term = Z(l-1) * exp(-delta*( times(l) - times(l-1) )) ;
          Z(l) = rec_term + precomp_diff.^delta ;
        end
    end
   else
       Z = 0;
   end
    S = Z;
    size(Z);
    S(S==0)=0;
    S(isnan(S))=0;
end
