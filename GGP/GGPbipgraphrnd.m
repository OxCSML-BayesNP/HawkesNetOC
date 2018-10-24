function [G, w1, w1_rem, w2, w2_rem, alpha1, sigma1, tau1, alpha2, sigma2, tau2, T1, T2]...
    = GGPbipgraphrnd(alpha1, sigma1, tau1, alpha2, sigma2, tau2, varargin)

%GGPbipgraphrnd samples a GGP bipartite graph.
% [G, w1, w1_rem, w2, w2_rem, alpha1, sigma1, tau1, alpha2, sigma2, tau2, T1, T2]...
%    = GGPbipgraphrnd(alpha1, sigma1, tau1, alpha2, sigma2, tau2, varargin)
%
% -------------------------------------------------------------------------
% INPUTS
%   - alpha1: positive scalar
%   - sigma1: real in (-Inf, 1)
%   - tau1: positive scalar
%   - alpha2: positive scalar
%   - sigma2: real in (-Inf, 1)
%   - tau2: positive scalar
% Optional input:
%   - T: truncation threshold; positive scalar
%
% OUTPUTS
%   - w1: sociability parameters of type 1 nodes with at least one connection
%   - w1_rem: sum of the sociability parameters of type 1 nodes with no connection
%   - w2: sociability parameters of type 2 nodes with at least one connection
%   - w2_rem: sum of the sociability parameters of type 2 nodes with no connection
%   - alpha1: Parameter alpha of the GGP of type 1 nodes
%   - sigma1: Parameter sigma of the GGP of type 1 nodes
%   - tau1: Parameter tau of the GGP of type 1 nodes
%   - alpha2: Parameter alpha of the GGP of type 2 nodes
%   - sigma2: Parameter sigma of the GGP of type 2 nodes
%   - tau2: Parameter tau of the GGP of type 2 nodes
%   - T1: threshold of type 1 nodes
%   - T2: threshold of type 2 nodes
% -------------------------------------------------------------------------
% EXAMPLE
% alpha1 = 100; sigma1 = -1; tau1 = 10;
% alpha2 = 50; sigma2 = 0.5; tau2 = 1;
% G = GGPbipgraphrnd(alpha1, sigma1, tau1, alpha2, sigma2, tau2);

% Copyright (C) Francois Caron, University of Oxford
% caron@stats.ox.ac.uk
% April 2015
%--------------------------------------------------------------------------

% Sample the parameters if needed
[alpha1, sigma1, tau1] = GGPgetparams(alpha1, sigma1, tau1);
[alpha2, sigma2, tau2] = GGPgetparams(alpha2, sigma2, tau2);

% Sample the weights w1 and w2
[w1, T1] = GGPrnd(alpha1, sigma1, tau1, varargin{:});
[w2, T2] = GGPrnd(alpha2, sigma2, tau2, varargin{:});

% Samples the graph conditional on the weights w using the conditional Poisson model
cumsum_w1 = [0; cumsum(w1)];
cumsum_w2 = [0; cumsum(w2)];
W1_star = cumsum_w1(end);  % Total mass of the GGP of type 1 nodes
W2_star = cumsum_w2(end);  % Total mass of the GGP of type 2 nodes
D_star = poissrnd(W1_star*W2_star); % Total number of directed edges

temp1 = W1_star * rand(D_star, 1);
temp2 = W2_star * rand(D_star, 1);
[~, bin1] = histc(temp1, cumsum_w1);
[~, bin2] = histc(temp2, cumsum_w2);
G = sparse(bin1, bin2, ones(numel(bin1), 1), numel(w1), numel(w2));
G = logical(G);

ind = sum(G, 1)>0;
G = G(:, ind);
w2_rem = sum(w2(~ind));
w2 = w2(ind);

ind = sum(G, 2)>0;
G = G(ind, :);
w1_rem = sum(w1(~ind));
w1 = w1(ind);
