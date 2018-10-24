function [h, centerbins, out] = plot_loglog(x, y, linespec, step)

%plot_loglog gives a log-log scale plot of the vectors  x Vs y
%
% INPUT
%   - x: vector
%   - y: vector
%   - linespec: line specification (determines line type, marker symbol, and color of the plotted lines)
%   - step: step size for the logarithmic bin edges in the pdf of the degree distribution
%
% OUTPUT
%   - h:  column vector of handles to the plotted lines
%   - centerbins: bin centers
%   - out: output values for y

if nargin<3
    linespec = '*';
end
if nargin<4
    step = 1;
end

edgebins = 2.^(0:step:12);
sizebins = edgebins;
sizebins = edgebins(2:end) - edgebins(1:end-1);

sizebins(end+1) = 1;
centerbins = edgebins;
[counts, bins] = histc(x, edgebins);
out = zeros(size(counts));
for i=1:numel(counts)
    if counts(i)>0
        out(i) = mean(y(bins==i));
    else
        out(i) = NaN;
    end
end
h = loglog(centerbins, out, linespec);
