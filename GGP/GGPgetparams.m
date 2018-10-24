function [alpha, sigma, tau] = GGPgetparams(alpha, sigma, tau, mode, sigma_neg, sigma_01)
%GGPgetparams samples parameters of a generalized gamma process.
% [alpha, sigma, tau] = GGPgetparams(alpha, sigma, tau)
% If the parameter is a scalar, its value is fixed
% else if it is a vector of size two, it is sampled from a gamma distribution
% -------------------------------------------------------------------------
% INPUTS
%   - alpha: scalar or vector of size 2
%   - sigma: scalar or vector of size 2
%   - tau: scalar or vector of size 2
% Optional inputs
%   - mode: string defining the action taken when hyperparameters are given.
%      - 'rnd' (default): samples parameters from gamma with hyperparameters given.
%      - 'default': return a default value. useful for checks.
%      - 'empty': return empty [].
%      - 'init': return initial value for mcmc.
%      - sigma_neg: logical indicator for 
%      - sigma_01
% OUTPUTS
%   - alpha: positive scalar
%   - sigma: real in (-Inf, 1)
%   - tau: positive scalar
% -------------------------------------------------------------------------
% EXAMPLE
% hyper_alpha = [100, 1]; hyper_sigma = .5; hyper_tau = [1, 0.1];
% N = GGPgetparams(hyper_alpha, hyper_sigma, hyper_tau);

% Copyright (c) F. Caron (University of Oxford), A. Todeschini (Inria), and 
% X. Miscouridou (University of Oxford)
% caron@stats.ox.ac.uk
% adrien.todeschini@gmail.com
% xenia.miscouridou@spc.ox.ac.uk
% September 2017
%--------------------------------------------------------------------------

if nargin<4
    mode = 'rnd';
end
if nargin < 5
    sigma_neg = false;
end
if nargin < 6
    sigma_01 = false;
end
if sigma_01 && sigma_neg
    error('sigma_neg and sigma_01 can''t be both true')
end

mode = validatestring(mode, {'rnd', 'default', 'empty', 'init'});

if numel(alpha)==2 % alpha is missing
    hyper_alpha = alpha;
    if hyper_alpha(1)<0 || hyper_alpha(2)<0
        error('Hyperparameters for alpha must be nonnegative');
    end
    switch mode
        case 'default'
            alpha = 1;
        case 'empty'
            alpha = [];
        case 'init'
            alpha = 100*rand;
%             alpha = 1000*rand; %%%% TEMP
%             alpha = 100; %%%% TEMP
        case 'rnd'
            alpha = gamrnd(hyper_alpha(1), 1./hyper_alpha(2));
    end
elseif numel(alpha)~=1
    error('alpha must be either scalar or vector of length 2')
end

if numel(sigma)==2 % sigma is missing
    hyper_sigma = sigma;
    if hyper_sigma(1)<0 || hyper_sigma(2)<0
        error('Hyperparameters for sigma must be nonnegative');
    end
    switch mode
        case 'default'
            sigma = 0;
        case 'empty'
            sigma = [];
        case 'init'
            if sigma_01
                sigma = .5;
            elseif sigma_neg
                sigma = -2*rand;
%                 sigma = -10*rand; %%%% TEMP
%                 sigma = -.5; %%%% TEMP
            else
%                 sigma = 2*rand - 1;
                sigma = 1.5*rand - 1; %%% not too high sigma
            end
        case 'rnd'
            if sigma_01
                sigma = betarnd(hyper_sigma(1), hyper_sigma(2));
            elseif sigma_neg
                sigma = -gamrnd(hyper_sigma(1), 1/hyper_sigma(2));
            else
                sigma = 1 - gamrnd(hyper_sigma(1), 1/hyper_sigma(2));
            end
    end
elseif numel(sigma)~=1
    error('sigma must be either scalar or vector of length 2')
elseif sigma_neg && sigma>=0
    error('sigma>=0 and sigma_neg=true not allowed')
elseif sigma_01 && sigma<0
    error('sigma<0 and sigma_01=true not allowed')
end

if numel(tau)==2 % tau is missing
    hyper_tau = tau;
    if hyper_tau(1)<0 || hyper_tau(2)<0
        error('Hyperparameters for tau must be nonnegative');
    end
    switch mode
        case 'default'
            tau = 1;
        case 'empty'
            tau = [];
        case 'init'
            tau = 10*rand;
%             tau = 1; %%% TEMP
        case 'rnd'
            tau = gamrnd(hyper_tau(1), 1./hyper_tau(2));
    end
elseif numel(tau)~=1
    error('tau must be either scalar or vector of length 2')
end


end
