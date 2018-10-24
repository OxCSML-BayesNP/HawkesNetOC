function [train_ind,times_matrix_forw_2,times_matrix_backw_2,train_data_forw,train_data_backw,test_data_forw,test_data_backw,N_train_links_forw, N_train_links_backw,N_test_links_forw,N_test_links_backw, T_train,max_T] = split_train_test(times_matrix_forw, times_matrix_backw, tr_perc)
% split_train_test splits the data set into train - test 
%
% [train_ind, times_matrix_forw_2, times_matrix_backw_2, train_data_forw, train_data_backw, test_data_forw, test_data_backw, N_train_links_forw, N_train_links_backw, N_test_links_forw, N_test_links_backw, T_train,max_T] =...
%       split_train_test(times_matrix_forw, times_matrix_backw, tr_perc)
%
% -------------------------------------------------------------------------
% INPUTS
%   
%   - train_data_forw:  the forward event times of the processes
%   - train_data_backw: the backward event times of the processes
%   - tr_perc: the training percentage for the training times (different to the proportion of 
%    training links to total links ) 
% 
% OUTPUTS
%   - train_ind: indices of training times
%   - times_matrix_forw_2: a matrix with all the forward event times of the processes
%   - times_matrix_backw_2: a matrix with all the backward event times of the processes
%   - train_data_forw: a matrix with the training forward  event times of the processes
%   - train_data_backw: a matrix with the training  backward event times of all the processes
%   - test_data_forw: a matrix with the test forward event times of all the processes
%   - test_data_backw: a matrix with the test backward event times of all the processes
%   - N_train_links_forw: number of forward  training links for each process
%   - N_train_links_backw: number of  backward training links for each process
%   - N_test_links_forw: number of  forward test links for each process
%   - N_test_links_backw: number of  backward test links for each process
%   - T_train: training (maximum) time
%   - max_T: overall maximum time
%
%  
% -------------------------------------------------------------------------
% Copyright (C) Xenia Miscouridou, University of Oxford
% xenia.miscouridou@spc.ox.ac.uk
% October 2018
%--------------------------------------------------------------------------




    max_T = min( max(max( times_matrix_forw)), max(max(times_matrix_backw)));
    T_train = tr_perc*max_T;


    [train_ind,~] = find(times_matrix_forw(:,1)<T_train);
    times_matrix_forw_2 = times_matrix_forw(train_ind,:);
    times_matrix_backw_2 = times_matrix_backw(train_ind,:);

    total_links_forw = sum((~isnan(times_matrix_forw_2)),2);
    train_data_forw = times_matrix_forw_2.*(times_matrix_forw_2 < T_train);
    train_data_forw(train_data_forw == 0) = NaN;
    N_train_links_forw = sum(~isnan(train_data_forw),2);
    
    
    total_links_backw = sum((~isnan(times_matrix_backw_2)),2);
    train_data_backw = times_matrix_forw_2.*(times_matrix_backw_2 < T_train);
    train_data_backw( train_data_backw == 0 ) = NaN;
    N_train_links_backw = sum(~isnan(train_data_backw),2);

    test_ind_forw = find( times_matrix_forw_2(:,1)>T_train );
    test_data_forw = times_matrix_forw_2.*(times_matrix_forw_2>=T_train);
    test_data_forw(test_data_forw == 0) = NaN;
    N_test_links_forw = sum(~isnan(test_data_forw),2);
    
    
    test_ind_backw = find( times_matrix_backw_2(:,1)>T_train );
    test_data_backw = times_matrix_backw_2.*(times_matrix_backw_2>=T_train);
    test_data_backw(test_data_backw == 0) = NaN;
    N_test_links_backw = sum(~isnan(test_data_backw),2);
      

    
end