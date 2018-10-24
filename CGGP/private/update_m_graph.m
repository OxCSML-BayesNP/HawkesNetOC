%% Update of the latent counts (Gibbs)
function [m] = update_m_graph(logw, ind1, ind2)

% update_m_graph updates the latent counts correspondning to the unobserved connection
% of each of the K nodes to each of the p communities 
% updates are from the full conditional
% rate = 2*w(ind1,:).*w(ind2,:) - w(ind1,:).^2.*repmat((ind1==ind2), 1, p);

%   - m: matrix of size [K,p] with latent counts
%   - logw: matrix of size [K,p] with weight parameters
%   - ind1: vector of indices of source nodes
%   - ind2: vector of indices of target nodes 

[K, p] = size(logw);

rate = log(2) + logw(ind1,:) + logw(ind2,:);
self = ind1==ind2;
rate(self,:) = rate(self,:) - log(2);
rate = exp(rate);
sum_rate = sum(rate, 2);

d = tpoissrnd(sum_rate);

% Sample n_ijk conditional on w1 and w2
% n = mnrnd(d, bsxfun(@rdivide, rate, sum_rate));
n = mnrnd2(d, bsxfun(@rdivide, rate, sum_rate)); %%% 100 times faster

% Compute feature counts m_ik
m = zeros(K, p);
for k=1:p
    nk = sparse(ind1, ind2, n(:,k), K, K);
    m(:,k) = sum(nk,2) + sum(nk,1)';
end
end