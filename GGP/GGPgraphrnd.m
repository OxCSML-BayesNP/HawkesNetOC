function [G, w, w_rem, alpha, sigma, tau, T] = GGPgraphrnd(alpha, sigma, tau, varargin)

% GGPgraphrnd samples a GGP graph.
% [G, w, w_rem, alpha, sigma, tau, T] = GGPgraphrnd(alpha, sigma, tau)
% [__] = GGPgraphrnd(alpha, sigma, tau, T)
% [__] = GGPgraphrnd(alpha, sigma, tau, T, maxiter)
%
% -------------------------------------------------------------------------
% INPUTS
%   - alpha: positive scalar
%   - sigma: real in (-Inf, 1)
%   - tau: positive scalar
% Optional input:
%   - T: truncation threshold; positive scalar
%   - maxiter: maximum number of iterations for the adaptive thinning
%     strategy (default=1e8); if nargin<5, it uses thinning
%
% OUTPUTS
%   - G: sparse binary matrix; graph
%   - w: sociability parameters of nodes with at least one connection
%   - wrem: sum of the sociability parameters of nodes with no connection
%   - alpha: Parameter alpha of the GGP
%   - sigma: Parameter sigma of the GGP
%   - tau: Parameter tau of the GGP
%   - T: threshold
% -------------------------------------------------------------------------
% EXAMPLE
% alpha = 100; sigma = 0.5; tau = 1;
% G = GGPgraphrnd(alpha, sigma, tau);

% Copyright (C) Francois Caron, University of Oxford
% caron@stats.ox.ac.uk
% April 2015
%--------------------------------------------------------------------------

% Sample the parameters if needed
[alpha, sigma, tau] = GGPgetparams(alpha, sigma, tau);

% Sample the weights w
[w, T] = GGPrnd(alpha, sigma, tau, varargin{:});

% Samples the graph conditional on the weights w using the conditional Poisson model
cumsum_w = [0; cumsum(w)];
W_star = cumsum_w(end);  % Total mass of the GGP
D_star = poissrnd(W_star^2); % Total number of directed edges

temp = W_star * rand(D_star, 2);
[~, bin] = histc(temp, cumsum_w);
[ind, ~, ib]  = unique(bin(:));
indlog = false(size(w));
indlog(ind) = true;
w_rem = sum(w(~indlog));
w = w(ind);
ib = reshape(ib, size(bin));
G = sparse(ib(:, 1), ib(:, 2), ones(size(ib, 1), 1), numel(ind), numel(ind));
G = logical(G + G');
