function h = plot_cost(C_st, rep, prefix, suffix)

% plot_cost plots the cost of each MCMC sample 
%  the cost function is a permutation-invariant absolute loss on the weights
% see munkres.m
%
% INPUT
%   - C_st: matrix of cost values for MCMC output for every chain
%   - rep: output directory to save the figure
%   - prefix: output filename prefix
%   - suffix: output filename suffix
% OUTPUT
%
% - h: figure handle vector


h = figure;
plot(C_st)
xlabel('MCMC iterations')
ylabel('Cost')
nchains = size(C_st, 2);
leg = cell(nchains, 1);
for k=1:nchains
    leg{k} = ['Chain' num2str(k)];
end
legend(leg, 'fontsize', 16, 'location', 'Best');
legend boxoff
box off
savefigs(gcf, [prefix 'cost' suffix], rep)