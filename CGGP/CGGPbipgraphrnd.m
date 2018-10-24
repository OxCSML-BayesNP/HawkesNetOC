function [G, w1, w1_rem, w2, w2_rem, alpha1, sigma1, tau1, Fdist1, gamma1, alpha2, sigma2, tau2, Fdist2, gamma2, T1, T2, w01, w02, beta1, beta2]...
    = CGGPbipgraphrnd(p, alpha1, sigma1, tau1, Fdist1, gamma1, alpha2, sigma2, tau2, Fdist2, gamma2, observe_all, necessarily_infinite, varargin)

% CGGPbipgraphrnd samples a GGP bipartite graph.
% [G, w1, w1_rem, w2, w2_rem, alpha1, sigma1, tau1, Fdist1, gamma1, alpha2, sigma2, tau2, Fdist2, gamma2, T1, T2]...
%    = CGGPbipgraphrnd(alpha1, sigma1, tau1, Fdist1, gamma1, alpha2, sigma2, tau2, Fdist2, gamma2, varargin)
%
% -------------------------------------------------------------------------
% INPUTS
%   - p: positive integer. number of features
%   - alpha1: positive scalar
%   - sigma1: real in (-Inf, 1)
%   - tau1: positive scalar
%   - Fdist2: distribution F
%   - gamma: positive column vector of length p
%   - alpha2: positive scalar
%   - sigma2: real in (-Inf, 1)
%   - tau2: positive scalar
%   - Fdist2: distribution F
%   - gamma2: positive column vector of length p
% Optional input:
%   - observe_all: logical, true if all nodes are observed (even if no edge) and sigma<0
%   - necessarily_infinite: logical, true forces an necessarily_infinite number of nodes: sigma in [0,1)
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
% p = 3; alpha1 = 100; sigma1 = -1; tau1 = 10;
% Fdist1.name = 'gamma'; Fdist1.param = [.03; .01; .01]; gamma1 = zeros(p,1);
% alpha2 = 50; sigma2 = 0.5; tau2 = 1;
% Fdist2.name = 'gamma'; Fdist2.param = [.01; .03; .01]; gamma2 = zeros(p,1);
% G = GGPbipgraphrnd(alpha1, sigma1, tau1, Fdist1, gamma1, alpha2, sigma2, tau2, Fdist2, gamma2);


% Copyright (c) F. Caron (University of Oxford), A. Todeschini (Inria), and 
% X. Miscouridou (University of Oxford)
% caron@stats.ox.ac.uk
% adrien.todeschini@gmail.com
% xenia.miscouridou@spc.ox.ac.uk
% September 2015
%--------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%
%%% TODO: return expected number of connections not generated
%%%%%%%%%%%%%%%%%

if nargin < 12
    observe_all = false(2,1);
elseif ~islogical(observe_all)
    error('observe_all must be logical')
elseif numel(observe_all)==1
    observe_all = observe_all * true(2,1);
end
if nargin < 13
    necessarily_infinite = false(2,1);
elseif ~islogical(necessarily_infinite)
    error('necessarily_infinite must be logical')
elseif numel(necessarily_infinite)==1
    necessarily_infinite = necessarily_infinite * true(2,1);
end

% Sample the parameters if needed
[alpha1, sigma1, tau1, Fdist1, gamma1] = CGGPgetparams(p, alpha1, sigma1, tau1, Fdist1, gamma1, 'rnd', observe_all(1), necessarily_infinite(2));
[alpha2, sigma2, tau2, Fdist2, gamma2] = CGGPgetparams(p, alpha2, sigma2, tau2, Fdist2, gamma2, 'rnd', observe_all(2), necessarily_infinite(2));

% Sample the weights w1 and w2
[w1, T1, w01, beta1] = CGGPrnd(p, alpha1, sigma1, tau1, Fdist1, gamma1, varargin{:});
[w2, T2, w02, beta2] = CGGPrnd(p, alpha2, sigma2, tau2, Fdist2, gamma2, varargin{:});

% Samples the graph conditional on the weights using the conditional Poisson model
cumsum_w1 = [zeros(p, 1), cumsum(w1, 1)'];
cumsum_w2 = [zeros(p, 1), cumsum(w2, 1)'];
W1_star = cumsum_w1(:,end);  % Total mass of the GGP of type 1 nodes
W2_star = cumsum_w2(:,end);  % Total mass of the GGP of type 2 nodes
D_star = poissrnd(W1_star.*W2_star); % Total number of directed edges

bin1 = zeros(sum(D_star), 1);
bin2 = zeros(sum(D_star), 1);

for k=1:p
    temp = W1_star(k) * rand(D_star(k), 1);
    [~, bink] = histc(temp, cumsum_w1(k, :));
    bin1(sum(D_star(1:(k-1)))+1:sum(D_star(1:k))) = bink;
    
    temp = W2_star(k) * rand(D_star(k), 1);
    [~, bink] = histc(temp, cumsum_w2(k, :));
    bin2(sum(D_star(1:(k-1)))+1:sum(D_star(1:k))) = bink;
end

clear temp bink

G = sparse(bin1, bin2, ones(numel(bin1), 1), size(w1,1), size(w2,1));
G = logical(G);

clear bin1 bin2

if observe_all(1)
    w1_rem = zeros(1,p);
else
ind = sum(G, 2)>0;
G = G(ind, :);
w1_rem = sum(w1(~ind,:));
w1 = w1(ind,:);
w01 = w01(ind);
beta1 = beta1(ind, :);
end

if observe_all(2)
    w2_rem = zeros(1,p);
else
    ind = sum(G, 1)>0;
    G = G(:, ind);
    w2_rem = sum(w2(~ind,:));
    w2 = w2(ind,:);
    w02 = w02(ind);
    beta2 = beta2(ind, :);
end
