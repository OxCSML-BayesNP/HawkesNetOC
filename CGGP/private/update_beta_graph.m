function [w, logw, beta, logbeta, rate] = update_beta_graph(w, logw, ...
    w0, logw0, beta, logbeta, w_rem, m, L, epsilon, Fdist, issimple)

% update_beta_graph updates the vectors of betas for all nodes when performing inference 
% on a simple graph conditional on the current values of the rest of the parameters and hyperparameters 
%
%   - w: matrix of size [K, p]; node weight parameters
%   - logw: log(w)
%   - w0: vector of size [K,1]; node base weights
%   - logwo : log(w0)
%   - beta : matrix of size [K, p]; beta scores
%   - logbeta : log(beta)
%   - w_rem : vector of size [1, p]; sum of remaining mass/weights
%   for each community (those not observed in the graph)
%   - m: matrix of size [K,p]; number of latent connections of the nodes
%   - L: positive integer; number of leapfrog steps
%   - epsilon: positive scalar; random walk standard deviation for the proposal in the HMC step  
%   - Fdist: structure; specifies the distribution of the beta scores
%   - issimple: logical; indicate whether the graph is simple or not

%   - rate: scalar; the MH acceptance rate
% ----------------------------------------------------------------------


[K, p] = size(w);
sum_w = sum(w,1) + w_rem;

logbetaprop = logbeta;
p_beta = randn(K, p);
pprop_beta = p_beta - epsilon/2 * grad_U_beta(m, w0, beta, sum_w, Fdist, issimple);
for lp=1:L
    logbetaprop = logbetaprop + epsilon*pprop_beta;
    betaprop = exp(logbetaprop);
    logwprop = bsxfun(@plus, logw0, logbetaprop);
    wprop = exp(logwprop);
    sum_wprop = sum(wprop, 1) + w_rem;
    if lp<L
        pprop_beta = pprop_beta  - epsilon * grad_U_beta(m, w0, betaprop, sum_wprop, Fdist, issimple);
    else
        pprop_beta = - pprop_beta + epsilon/2 * grad_U_beta(m, w0, betaprop, sum_wprop, Fdist, issimple);
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
    + sum(a .* sum(logbetaprop - logbeta, 1))...
    - sum(b .* sum(betaprop - beta, 1));

logaccept = logaccept -.5*sum(pprop_beta(:).^2-p_beta(:).^2);
if issimple % If simple graph, do not take into account self-connections
    logaccept = logaccept + sum(wprop(:).^2 - w(:).^2);
end

if isnan(logaccept)
    logaccept = -Inf;
end

if log(rand)<logaccept
    w = wprop;
    logw = logwprop;
    beta = betaprop;
    logbeta = logbetaprop;
end
rate = min(1, exp(logaccept));

end

%% Gradient

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
