function plot_counting_process(times_for, times_rec, max_T)
% plot_counting_process plots the counting process for a given pair of processes
% i e plots the integer valued processes 
% N_{ij} and N_{ji} that have intensities
% lambda_{ij} and lambda_{ji} respectively
%
% -------------------------------------------------------------------------
% INPUTS
%   
%   - times_for: the forward event times for the process
%   - times_rec: the backward event times for the process
%   - max_T: the maximum time for an observed event for both processes
% 
% OUTPUTS
%
%   Returns the counting process plot.
%  
% -------------------------------------------------------------------------
% Copyright (C) Xenia Miscouridou, University of Oxford
% xenia.miscouridou@spc.ox.ac.uk
% October 2018
%--------------------------------------------------------------------------

    figure;
    TIMES{1}=times_for;TIMES{2}=times_rec;
    for i=1:2
        subplot(2,1,i)
        times=TIMES{i};
        num_events = numel(times);
        t = (0:0.001:(max_T))'; 
        % Start with all zeros: 
        unitstep = zeros(size(t)); 
        unitstep(t<times(1)) = 0;unitstep(t>=times(num_events)) = num_events ; 

        for j=1:(num_events-1)
            % But make everything corresponding to t>=1 one:
            unitstep( (t>=times(j)) & (t<times(j+1))) = j; 
            hold on;
            plot(times(j),0,'r*')
        end
        hold on;
        plot(t,unitstep,'b','linewidth',3)
        plot(times(num_events),0,'r*')
        
        % Repeat, with everything shifted to the right by 1 unit:  
        hold on
        plot([times(num_events) times(num_events)],[0 num_events],'r:')        % Vertical Line
        box off
        ylabel('N')
        xlabel('t')
    end

end
