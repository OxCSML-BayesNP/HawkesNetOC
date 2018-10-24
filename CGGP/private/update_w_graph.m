function [w, logw, w0, logw0, beta, logbeta, rate] = update_w_graph(w, logw, ...
    w0, logw0, beta, logbeta, w_rem, m, L, epsilon, sigma, tau, Fdist, issimple)

% update_w_graph updates the weight parameters w for inference on a simple graph
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
%     for each community (these weights are not observed in the graph)
%   - m: matrix of size [K,p] with the number of latent connections of each
%     node for each of the latent communities
%   - L: positive integer, number of leapfrog steps
%   - epsilon: [positive scalar, the random walk standard deviation for the proposal in the HMC step  
%   - sigma: scalar < 1, parameter of the CGGP  
%   - tau: positive scalar, parameter of the CGGP
%   - Fdist: structure; distribution for the beta scores
%   - issimple: logical to indicate whther the graph is simple or not

%   - rate: the HMC acceptance rate
% ----------------------------------------------------------------------


% Update w
[K, p] = size(w);
sum_w = sum(w,1) + w_rem;

p_w0 = randn(K, 1);
p_beta = randn(K, p);
logw0prop = logw0;
logbetaprop = logbeta;
pprop_w0 = p_w0 - epsilon/2 * grad_U_w0(m, w0, beta, sum_w, sigma, tau, issimple);
pprop_beta = p_beta - epsilon/2 * grad_U_beta(m, w0, beta, sum_w, Fdist, issimple);
for lp=1:L
    logw0prop = logw0prop + epsilon*pprop_w0;
    logbetaprop = logbetaprop + epsilon*pprop_beta;
    w0prop = exp(logw0prop);
    betaprop = exp(logbetaprop);
    logwprop = bsxfun(@plus, logw0prop, logbetaprop);
    wprop = exp(logwprop);
    sum_wprop = sum(wprop, 1) + w_rem;
    if lp<L
        pprop_w0 = pprop_w0  - epsilon * grad_U_w0(m, w0prop, betaprop, sum_wprop, sigma, tau, issimple);
        pprop_beta = pprop_beta  - epsilon * grad_U_beta(m, w0prop, betaprop, sum_wprop, Fdist, issimple);
    else
        pprop_w0 = - pprop_w0 + epsilon/2 * grad_U_w0(m, w0prop, betaprop, sum_wprop, sigma, tau, issimple);
        pprop_beta = - pprop_beta + epsilon/2 * grad_U_beta(m, w0prop, betaprop, sum_wprop, Fdist, issimple);
    end
end

switch Fdist.name
    case 'gamma'
        if isnumeric(Fdist.param)
            a = Fdist.param';
            b = a;
        else
            a = Fdist.param.a';
            b = Fdist.param.b';
        end
    otherwise
        error('unknown distribution %s', Fdist.name);
end

logaccept = -sum(sum_wprop.^2) + sum(sum_w.^2) ...
    + sum(sum(m.*(logwprop - logw)))...
    - sigma * sum(logw0prop - logw0)...
    - tau * sum(w0prop - w0)...
    + sum(a .* sum(logbetaprop - logbeta, 1))...
    - sum(b .* sum(betaprop - beta, 1));

logaccept = logaccept -.5*sum(pprop_w0.^2-p_w0.^2) -.5*sum(pprop_beta(:).^2-p_beta(:).^2);
if issimple % If simple graph, do not take into account self-connections
    logaccept = logaccept + sum(wprop(:).^2 - w(:).^2);
end

if isnan(logaccept)
    warning('logaccept is nan')
    logaccept = -Inf;
end

if log(rand)<logaccept
    w = wprop;
    logw = logwprop;
    w0 = w0prop;
    logw0 = logw0prop;
    beta = betaprop;
    logbeta = logbetaprop;
end
rate = min(1, exp(logaccept));

end

%% Gradient
function [out] = grad_U_w0(m, w0, beta, sum_w, sigma, tau, issimple)
if issimple
    out = - sum(m, 2) + sigma + w0 .* (tau + 2*(beta*sum_w'-sum(beta.^2, 2)));
else
    out = - sum(m, 2) + sigma + w0 .* (tau + 2*(beta*sum_w'));
end

end

function [out] = grad_U_beta(m, w0, beta, sum_w, Fdist, issimple)
switch Fdist.name
    case 'gamma'
        if isnumeric(Fdist.param)
            a = Fdist.param';
            b = a;
        else
            a = Fdist.param.a';
            b = Fdist.param.b';
        end
        if issimple
            out = - bsxfun(@plus, a, m) + beta .* bsxfun(@plus, b, 2*bsxfun(@minus, w0*(sum_w), w0.^2));
        else
            out = - bsxfun(@plus, a, m) + beta .* bsxfun(@plus, b, 2*w0*sum_w);
        end
    otherwise
        error('unknown distribution %s', Fdist.name);
end
end
