function [w, logw, w0, logw0, rate] = update_w0_graph(w, logw, ...
    w0, logw0, beta, logbeta, w_rem, m, L, epsilon, sigma, tau, issimple)


% update_w0_graph updates the base weight parameters w_0 for inference on a simple graph
% in a Markov chain Monte Carlo sampler conditional on the current values
% of the rest of the parameters and hyperparameters
%
%   - w: matrix of size [K, p] with each row correspondning to the nodes'
%     weights for wach of the p communities
%   - logw: log(w)
%   - w0: vector of size [K,1] with each node's base weights
%   - logwo : log(w0)
%   - beta : matrix of size [K, p] of beta scores
%   - logbeta : log(beta)
%   - w_rem : vector of size [1, p] with the sum of remaining mass/weights
%   for each community (these weights are not observed in the graph)
%   - m: matrix of size [K,p] with the number of latent connections of each
%   node for each of the latent communities
%   - L: positive integer, number of leapfrog steps
%   - epsilon: [positive scalar, the random walk standard deviation for the proposal in the HMC step  
%   - sigma: scalar < 1, parameter of the CGGP  
%   - tau: positive scalar, parameter of the CGGP
%   - Fdist: structurel specifies the distribution for the beta scores
%   - issimple: logical to indicate whther the graph is simple or not

%   - rate: the MH acceptance rate
% ----------------------------------------------------------------------


K = size(w, 1);
sum_w = sum(w,1) + w_rem;

logw0prop = logw0;
p_w0 = randn(K, 1);
pprop_w0 = p_w0 - epsilon/2 * grad_U_w0(m, w0, beta, sum_w, sigma, tau, issimple);
for lp=1:L
    logw0prop = logw0prop + epsilon*pprop_w0;
    w0prop = exp(logw0prop);
    logwprop = bsxfun(@plus, logw0prop, logbeta);
    wprop = exp(logwprop);
    sum_wprop = sum(wprop, 1) + w_rem;
    if lp<L
        pprop_w0 = pprop_w0  - epsilon * grad_U_w0(m, w0prop, beta, sum_wprop, sigma, tau, issimple);
    else
        pprop_w0 = - pprop_w0 + epsilon/2 * grad_U_w0(m, w0prop, beta, sum_wprop, sigma, tau, issimple);
    end
end

logaccept = -sum(sum_wprop.^2) + sum(sum_w.^2) ...
    + sum(sum(m.*(logwprop - logw)))...
    - sigma * sum(logw0prop - logw0)...
    - tau * sum(w0prop - w0);

logaccept = logaccept -.5*sum(pprop_w0.^2-p_w0.^2);
if issimple % If simple graph, do not take into account self-connections
    logaccept = logaccept + sum(wprop(:).^2 - w(:).^2);
end

if isnan(logaccept)
    logaccept = -Inf;
end

if log(rand)<logaccept
    w = wprop;
    logw = logwprop;
    w0 = w0prop;
    logw0 = logw0prop;
end
rate = min(1, exp(logaccept));

end


%% Gradient
function [out] = grad_U_w0(m, w0, beta, sum_w, sigma, tau, issimple)
if issimple
    out = - sum(m, 2) + sigma + w0 .* (tau + 2*(beta.*sum_w'-sum(beta.^2,2)));
else
    out = - sum(m, 2) + sigma + w0 .* (tau + 2*(beta*sum_w'));
end

end