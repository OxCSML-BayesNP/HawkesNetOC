function [times_matrix_forw, times_matrix_backw ,G, BASE, pairs] = preprocess(source, target, times, p,istest, estimate_base)
% preprocess preprocess data and estimates the base intensities if specified
%
% -------------------------------------------------------------------------
% INPUTS
%   
%   - source:  a vector with all the sender nodes for each event timepoint
%   - target: a vector with all the receiver nodes for each event timepoint
%   - times: a vector with all the event times/interactions
%   - p: number of latent communities for the base intensity
%   - istest: logical to indicate if this is a test
%   - estimate_base: logical to indicate if estimation of base intensities
%       is desired
% 
% OUTPUTS
%
%   - times_matrix_forw: a matrix with all the forward event times of the processes
%   - times_matrix_backw: a matrix with all the  backward event times of the processes
%   - G: logical sparse matrix, the corresponding graph 
%   - BASE: the base intensities for each process
%   - pairs: source and target node pairs corresponding to each process (forward and backward)
%
%  
% -------------------------------------------------------------------------
% Copyright (C) Xenia Miscouridou, University of Oxford
% xenia.miscouridou@spc.ox.ac.uk
% October 2018
%--------------------------------------------------------------------------



    %remove loops
    non_self_index = (find(source ~= target));
    source = source(non_self_index); target = target(non_self_index); 
    times = times(non_self_index);   

    %arrange data in the right form
    [times_matrix_forw, times_matrix_backw, pairs, G] = prepare_data(source,target,times);
    n_HP = nnz(G);
    figure;spy(G)


    sumexpoS = zeros(n_HP,max(pairs(:,3)));

    bins =zeros(size(times_matrix_backw));
    max_T = max(times);
    for l=1:n_HP

        tfor = [0,times_matrix_forw(l,1:pairs(l,3)),max_T];
        tback = times_matrix_backw(l,:);
        [~, y, bin] = histcounts(tback,tfor);
        bins(l,1:length(bin))=bin;

        y = bin(bin~=0);
        nRow = length(y);
        nCol = max(y);
        A = zeros(nRow,nCol);

        for i=1:nRow
        A(i,y(i)) = 1; 
        end
       sumexpoS(l,1:size(A,2)) =  exp(times_matrix_backw(l,1:pairs(l,4)))*A;

    end

    %run MCMC on the whols graph to estimate BASE parameters
    addpath ./Hawkes ./utils

    G_symm = sparse([pairs(:,1);pairs(:,2)],[pairs(:,2);pairs(:,1)],1);
    ind = any( G_symm );G_symm = G_symm(ind, ind);
    Gforbase = G_symm|G_symm';

    invmap = find( full(ind) )';
    if estimate_base
        if istest
            n_chains_base = 1;
            p=2;
            estimates = estimate_BASE(Gforbase, p, n_chains_base, istest);

        else
            n_chains_base = 2;
            estimates = estimate_BASE(Gforbase, p, n_chains_base, false);

        end


        BASE = zeros(n_HP,1);
        for l=1:n_HP     
            BASE(l)= sum(estimates.w((invmap==pairs(l,1)),:).*estimates.w((invmap==pairs(l,2)),:));
        end

        if size(BASE,1) ~= n_HP
            error('check dims')
        end
    else
           BASE = zeros(n_HP,1);
    end


end