function h = plot_logpost(logp, rep, prefix, suffix)

% plot_logpost plots the given log posterior probability distributions 
%
% INPUT
%   - logp: log posterior probabilities
%   - rep:  output directory
%   - prefix: output filename prefix
%   - suffix: output filename suffix
%
% OUTPUT
% h = vector of handles to the plotted lines.

h = figure;
plot(logp)
xlabel('MCMC iterations')
ylabel('Log-posterior')
nchains = size(logp, 2);
leg = cell(nchains, 1);
for k=1:nchains
    leg{k} = ['Chain' num2str(k)];
end
legend(leg, 'fontsize', 16, 'location', 'Best');
legend boxoff
box off
savefigs(gcf, [prefix 'logpost' suffix], rep)