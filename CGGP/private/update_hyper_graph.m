function [w_rem, alpha, logalpha, sigma, tau, Fdist, gamma, rate] = ...
    update_hyper_graph(sum_w, beta, logbeta, w0, logw0, w_rem, alpha, logalpha, sigma,...
    tau, Fdist, gamma, T, estim, hyper, rw_std, rw_alpha, MH_nb, sigma_01, sigma_neg, sumsmall_arg, varargin)

% update_hyper_graph updates the hyperparameters using MH
% when performing inference on a simple graph conditional on the current values
% of the rest of the parameters and hyperparameters 
%
%   - sum_w2: vector of size [1, p] with the sum of weights of type 2 
%   - beta : matrix of size [K, p] of beta scores
%   - logbeta : log(beta)
%   - w0: vector of size [K, 1] with each node's base weights
%   - logwo : log(w0)
%   - w_rem : vector of size [1, p] with the sum of remaining mass/weights
%   - alpha: positive scalar define the graph's restriction
%   - logalpha : log(alpha)
%   - sigma: scalar < 1, parameter of the CGGP  
%   - tau: positive scalar, parameter of the CGGP
%   - Fdist: distribution for the beta scores
%   - gamma: vector of size [1,p] with tilting parameters gamma
%   - T: positive scalar, threshold above which w_0 are generated
%   - estim: struct with logical fields to indicate the estimation or not
%   of the parameters
%   - hyper: structure with the hyperparameters a and b for the Fdist
%   - rw_std: standard deviation for random walk for the MH update of sigma 
%   - rw_alpha: standard deviation for random walk for the MH update of alpha 
%   - MH_nb: number of MH steps
%   - sigma_01: logical indicating whether sigma = 1  
%   - sigma_neg: logical indicating whether sigma is negative  
%   - sumsmall_arg
%
%   - gamma:  laplace transform exponent
%   - rate: rate of acceptance for MH
% --------------------------------------------------------------


if nargin<21
    sumsmall_arg = {};
end

[N,p] = size(beta);

lograte = zeros(MH_nb,1);

alpha_prop = alpha;
logalpha_prop = logalpha;
sigma_prop = sigma;
tau_prop = tau;
Fdist_prop = Fdist;
gamma_prop = gamma;
w_rem_prop = w_rem;

if estim.gamma || any(gamma>0)
    error('*** FIXME not implemented yet ***')
end

for nn=1:MH_nb
    
    %% proposal
    if estim.sigma
        % random walk proposal
        if sigma_01
            sigma_prop = logistic(logit(sigma) + rw_std.sigma.*randn);
        elseif sigma_neg
            sigma_prop = -exp(log(-sigma) + rw_std.sigma*randn);
        else
            sigma_prop = 1-exp(log(1-sigma) + rw_std.sigma*randn);
        end
    end
    if estim.tau
        % random walk proposal
        tau_prop = exp(log(tau) + rw_std.tau*randn);
    end
    switch Fdist.name
        case 'gamma'
            if isnumeric(Fdist.param)
                a = Fdist.param;
                a_prop = Fdist_prop.param;
                estim.a = estim.Fparam;
                estim.b = false;
                hyper.a = hyper.Fdist.param;
                if estim.a
                    % random walk
                    a_prop = exp(log(a) + rw_std.a .* randn(size(a)));
                end
                b = a;
                b_prop = a_prop;
                Fdist_prop.param = a_prop;
            else
                a = Fdist.param.a;
                b = Fdist.param.b;
                a_prop = Fdist_prop.param.a;
                b_prop = Fdist_prop.param.b;
                estim.a = estim.Fparam.a;
                estim.b = estim.Fparam.b;
                hyper.a = hyper.Fdist.param.a;
                hyper.b = hyper.Fdist.param.b;
                if estim.a
                    % random walk proposal
                    a_prop = exp(log(a) + rw_std.a .* randn(size(a)));
                end
                if estim.b
                    % random walk proposal
                    b_prop = exp(log(b) + rw_std.b .* randn(size(b)));
                end
                Fdist_prop.param.a = a_prop;
                Fdist_prop.param.b = b_prop;
            end
        otherwise
            error('Unknown distribution F')
    end
    if estim.alpha
        if rw_alpha % random walk proposal
            logalpha_prop = logalpha + rw_std.alpha*randn;
            alpha_prop = exp(logalpha_prop);
        else % gamma proposal
            if sigma_neg && hyper.observe_all
                lambda_prop = gamrnd(hyper.alpha(1) + N, 1/(-hyper.alpha(2)*sigma_prop*tau_prop^(-sigma_prop) + 1));
                logalpha_prop = log(lambda_prop)-sigma_prop*log(tau_prop)+log(-sigma_prop);
                alpha_prop = exp(logalpha_prop);
            else
                psi_prop = CGGPpsi(2*sum_w + w_rem, sigma_prop, tau_prop, Fdist_prop, gamma_prop', varargin{:});
                alpha_prop = gamrnd(hyper.alpha(1) + N, 1/(hyper.alpha(2) + psi_prop ));
                logalpha_prop = log(alpha_prop);
            end
        end
    end
    
    
    if estim.w_rem
        w_rem_prop = CGGPsumrnd(p, logalpha_prop, sigma_prop, tau_prop, ...
            Fdist_prop, gamma_prop + 2*sum_w' + w_rem', T, sumsmall_arg{:}, varargin{:});
    end
    
    %% Compute the acceptance probability
    lar = 0; % log acceptance ratio
    
    if estim.sigma
        if sigma_01
            if hyper.sigma(1)>0
                lar = lar + hyper.sigma(1)*(log(sigma_prop)-log(sigma));
            end
            if hyper.sigma(2)>0
                lar = lar + hyper.sigma(2)*(log(1-sigma_prop)-log(1-sigma));
            end
        elseif sigma_neg
            if hyper.sigma(1)>0
                lar = lar + hyper.sigma(1)*(log(-sigma_prop)-log(-sigma));
            end
            lar = lar - hyper.sigma(2)*(-sigma_prop+sigma);
        else
            if hyper.sigma(1)>0
                lar = lar + hyper.sigma(1)*(log(1-sigma_prop)-log(1-sigma));
            end
            lar = lar - hyper.sigma(2)*(-sigma_prop+sigma);
        end
        lar = lar - N*(gammaln(1-sigma_prop)-gammaln(1-sigma))...
            - (sigma_prop-sigma)*sum(logw0);
    end
    if estim.tau
        if hyper.tau(1)>0
            lar = lar + hyper.tau(1)*(log(tau_prop)-log(tau));
        end
        lar = lar - (hyper.tau(2) + sum(w0))*(tau_prop-tau);
    end
    if estim.alpha
        if rw_alpha
            lar = lar + (hyper.alpha(1) + N)*(logalpha_prop-logalpha)...
                - hyper.alpha(2)*(alpha_prop-alpha);
        else
            if sigma_neg && hyper.observe_all
                lar = lar + (hyper.alpha(1) + N)*(log(hyper.alpha(2) - tau^sigma/sigma) - log(hyper.alpha(2) - tau_prop^sigma_prop/sigma_prop));
            else
                psi = CGGPpsi(2*sum_w + w_rem_prop, sigma, tau, Fdist, gamma', varargin{:});
                lar = lar + (hyper.alpha(1) + N)*(log(hyper.alpha(2) + psi) - log(hyper.alpha(2) + psi_prop));
            end
        end
    end
    if ~estim.alpha || rw_alpha
        if sigma_neg && hyper.observe_all
            lar = lar + alpha_prop/sigma_prop*(tau_prop^sigma_prop) - alpha/sigma*(tau^sigma);
        else
            psi_prop = CGGPpsi(2*sum_w + w_rem, sigma_prop, tau_prop, Fdist_prop, gamma_prop', varargin{:});
            psi = CGGPpsi(2*sum_w + w_rem_prop, sigma, tau, Fdist, gamma', varargin{:});
            lar = lar - alpha_prop*psi_prop + alpha*psi;
        end
    end
    if estim.a
        if hyper.a(1)>0
            lar = lar + sum(hyper.a(:,1).*(log(a_prop)-log(a)));
        end
        lar = lar - sum(hyper.a(:,2).*(a_prop-a))...
            + sum((a_prop-a).*sum(logbeta,1)')...
            + N*sum(fillmat(-gammaln(a_prop)+gammaln(a), [p, 1]));
    end
    if estim.b
        if hyper.b(1)>0
            lar = lar + sum(hyper.b(:,1).*(log(b_prop)-log(b)));
        end
        lar = lar - sum(hyper.b(:,2).*(b_prop-b))...
            - sum((b_prop-b).*sum(beta,1)');
    end
    if estim.a || estim.b
        lar = lar + N*sum(fillmat(a_prop.*log(b_prop)-a.*log(b), [p,1]));
    end
    if estim.w_rem
        lar = lar - sum(w_rem_prop.^2) + sum(w_rem.^2);
    end
    
    if isnan(lar)
        keyboard
    end
    
    %% Accept-reject step
    if log(rand) < lar
        w_rem = w_rem_prop;
        alpha = alpha_prop;
        logalpha = logalpha_prop;
        sigma = sigma_prop;
        tau = tau_prop;
        Fdist = Fdist_prop;
        gamma = gamma_prop;
    end
    
    lograte(nn) = lar;
end
rate = mean(min(1, exp(lograte)));
end

function x = logit(p)
    x = log(p) - log(1-p);
end

function p  = logistic(x)
    p = 1./(1+exp(-x));
end