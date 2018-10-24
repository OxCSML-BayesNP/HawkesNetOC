function [h2, centerbins, freq] = plot_degree(G, linespec, step)

% plot_degree plots the degree of each node for the observed
% adgacency matrix G
%
% INPUT
%   - G: observed binary adjacency matrix
%   - linespec: line specification (determines line type, marker symbol, and color of the plotted lines.)
%   - step: step size for the logarithmic bin edges in the pdf of the degree distribution
%
%   OUTPUT
%   - h2: loglog degree distribution figure
%   - centerbins: bin centers
%   - freq: frequencies for the counts of the degrees in each bin


if nargin<2
    linespec = '*';
end
if nargin<3
    step = 1;
end

G(isnan(G)) = 0.5; % fill missing by 0.5
deg = full(sum(G));

% Uses logarithmic binning to get a less noisy estimate of the
% pdf of the degree distribution

edgebins = 2.^(0:step:16);
sizebins = edgebins;
sizebins = edgebins(2:end) - edgebins(1:end-1);

sizebins(end+1) = 1;
centerbins = edgebins;
counts = histc(deg, edgebins);
freq = counts./sizebins/size(G, 1);
h2 = loglog(centerbins, freq, linespec);
xlabel('Degree', 'fontsize', 16)
ylabel('Distribution', 'fontsize', 16)
