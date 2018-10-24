function [G, w, w_rem, alpha, sigma, tau, Fdist, gamma, T, w0, beta] = CGGPgraphrnd(p, alpha, sigma, tau, Fdist, gamma, observe_all, necessarily_infinite, varargin)

% CGGPgraphrnd samples a Compound GGP graph.
% [G, w, w_rem, alpha, sigma, tau, Fdist, gamma] = CGGPgraphrnd(p, alpha, sigma, tau, Fdist)
% [__] = CGGPgraphrnd(p, alpha, sigma, tau, Fdist, gamma)
% [__] = CGGPgraphrnd(p, alpha, sigma, tau, Fdist, gamma, observe_all, necessarily_infinite)
% [__] = CGGPgraphrnd(p, alpha, sigma, tau, Fdist, gamma, observe_all, necessarily_infinite, T)
% [__] = CGGPgraphrnd(p, alpha, sigma, tau, Fdist, gamma, observe_all, necessarily_infinite, T, maxiter)
%
% -------------------------------------------------------------------------
% INPUTS
%   - p: positive integer. number of features
%   - alpha: positive scalar
%   - sigma: real in (-Inf, 1)
%   - tau: positive scalar
%   - Fdist: distribution F
% Optional input:
%   - gamma: positive column vector of length p (default=0)
%   - observe_all: logical, true if all nodes are observed (even if no edge) and sigma<0
%   - necessarily_infinite: logical, true forces an necessarily_infinite number of nodes: sigma in [0,1)
%   - T: truncation threshold; positive scalar
%
% OUTPUTS
%   - G: sparse logical matrix. the generated graph
%   - w: feature interest parameters of nodes with at least one connection
%   - w_rem: sum of the sociability parameters of nodes with no connection
%   - alpha: Parameter alpha of the GGP
%   - sigma: Parameter sigma of the GGP
%   - tau: Parameter tau of the GGP
%   - Fdist: Structure with parameters of F
%   - gamma: Parameters gamma
%   - T: Threshold
% -------------------------------------------------------------------------
% EXAMPLE
% p = 3; alpha = 100; sigma = 0.5; tau = 1e-4;
% Fdist.name = 'gamma'; Fdist.param = 0.01*ones(p,1);
% G = CGGPgraphrnd(p, alpha, sigma, tau, Fdist);

% Copyright (c) F. Caron (University of Oxford), A. Todeschini (Inria), and 
% X. Miscouridou (University of Oxford)
% caron@stats.ox.ac.uk
% adrien.todeschini@gmail.com
% xenia.miscouridou@spc.ox.ac.uk
% September 2017
%--------------------------------------------------------------------------


if nargin < 6
    gamma = zeros(p,1);
end
if nargin < 7
    observe_all = false;
elseif ~islogical(observe_all)
    error('observe_all must be logical')
end
if nargin < 8
    necessarily_infinite = false;
elseif ~islogical(necessarily_infinite)
    error('necessarily_infinite must be logical')
end

% Sample the parameters if needed
[alpha, sigma, tau, Fdist, gamma] = CGGPgetparams(p, alpha, sigma, tau, Fdist, gamma, 'rnd', observe_all, necessarily_infinite);

% Sample the weights w
[w, T, w0, beta] = CGGPrnd(p, alpha, sigma, tau, Fdist, gamma, varargin{:});

% Sample the directed graph conditional on the weights w using the conditional Poisson model
cumsum_w = [zeros(p, 1), cumsum(w, 1)'];
W_star = cumsum_w(:, end);  % Total mass of the GGP
D_star = poissrnd(W_star.^2); % Total number of directed edges

nedges = sum(D_star); % total number of edges (all features)
bin = zeros(nedges, 2); % each row is an edge

for k=1:p
    % sample edges of kth feature
    temp = W_star(k) * rand(D_star(k), 2);
    [~, bink] = histc(temp, cumsum_w(k, :));
    
    % append them in bin
    bin(sum(D_star(1:(k-1)))+1:sum(D_star(1:k)), :) = bink;
end

clear temp bink

if observe_all
    K = size(w, 1); % number of nodes
    w_rem = zeros(1,p);
else
    % indices of nodes with at least one edge
    [ind, ~, ib]  = unique(bin(:));
    K = numel(ind); % number of nodes
    indlog = false(size(w, 1), 1);
    indlog(ind) = true;
    w_rem = sum(w(~indlog, :), 1);
    switch(Fdist.name)
        case 'gamma'
            if isnumeric(Fdist.param)
                a = Fdist.param';
                b = a;
            else
                a = Fdist.param.a';
                b = Fdist.param.b';
            end
        otherwise
            error('unknown distribution')
    end
    w_rem_small = CGGPsumsmallrnd(p, log(alpha), sigma, tau, a, b, gamma', T);
    w_rem = w_rem + w_rem_small;
    clear indlog
    w = w(ind, :);
    w0 = w0(ind);
    beta = beta(ind, :);
    bin = reshape(ib, size(bin));
end

G = sparse(bin(:, 1), bin(:, 2), 1, K, K);
G = logical(G + G');

end
