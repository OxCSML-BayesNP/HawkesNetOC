function [times_matrix_forw, times_matrix_backw, pairs, G] = prepare_data(source,target,times)
% prepare_data prepares the data and brings it in the right format
%
% -------------------------------------------------------------------------
% INPUTS
%   
%   - source:  the sender nodes for each event timepoint
%   - target: the receiver nodes for each event timepoint
%   - times: a vector with all the event times/interactions
% 
% OUTPUTS
%
%   - times_matrix_forw: a matrix with all the forward event times of the processes
%   - times_matrix_backw: a matrix with all the backward event times of the processes
%   - pairs: source and target node pairs corresponding to each process (forward and backward)
%   - G: logical sparse matrix, the corresponding graph 
%
%  
% -------------------------------------------------------------------------
% Copyright (C) Xenia Miscouridou, University of Oxford
% xenia.miscouridou@spc.ox.ac.uk
% October 2018
%--------------------------------------------------------------------------


    non_self_index=(find(source~=target));
    source=source(non_self_index); target=target(non_self_index); times=times(non_self_index);   
    
    G = sparse(source,target, 1);

    
    [ST(:,1),ST(:,2), d_forw] = find(G);
    invmap = [1:max(size(G))]';
    %create times matrix
    n_HP = size(ST,1);
    max_numevents = max(d_forw);
    times_matrix_forw = zeros(n_HP,max_numevents)*NaN;
    times_matrix_backw = zeros(n_HP,max_numevents)*NaN;
    d_backw = zeros(n_HP,1);
    
    for l=1:n_HP
        if rem(l,1000)==0
            fprintf('...processing data for edge %d \n...', l)
        end
        %store the forward events
        index = (source == invmap(ST(l,1))&target == invmap(ST(l,2)));
        
             t = times(index);
             times_matrix_forw(l,1:d_forw(l)) = t;
         
        %store the backward events - may be empty
        index2 = (source == invmap(ST(l,2))&target == invmap(ST(l,1)));
        if ~sum(index2)==0
            d_backw(l) = nnz(index2);
            t2 = times(index2);
            times_matrix_backw(l,1:d_backw(l)) = t2;            
        end       
    end
       
    times_matrix_forw = bsxfun(@sort, times_matrix_forw,2);
    times_matrix_backw = bsxfun(@sort, times_matrix_backw,2);
    pairs = [ST,d_forw,d_backw];
  
end