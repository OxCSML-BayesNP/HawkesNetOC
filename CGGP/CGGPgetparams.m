function [alpha, sigma, tau, Fdist, gamma] = CGGPgetparams(p, alpha, sigma, tau, Fdist, gamma, mode, varargin)

% CGGPgetparams samples parameters of a generalized gamma process.
% [alpha, sigma, tau, Fdist, gamma] = CGGPgetparams(p, alpha, sigma, tau, Fdist, gamma, mode)
% If the parameter is a scalar, its value is fixed
% else if it is a vector of size two, it is sampled from a gamma distribution
% -------------------------------------------------------------------------
% INPUTS
%   - p: positive integer
%   - alpha: scalar or vector of size 2
%   - sigma: scalar or vector of size 2
%   - tau: scalar or vector of size 2
%   - Fdist:  structure; specifies the distribution of the beta vectors
%   - gamma: column vector of size p or matrix of size [p,2]
% Optional inputs
%   - mode: string defining the action taken when hyperparameters are given.
%      - 'rnd' (default): samples parameters from gamma with hyperparameters given.
%      - 'default': return a default value. useful for checks.
%      - 'empty': return empty [].
%      - 'init': return initial value for mcmc.
%     
% OUTPUTS
%   - alpha: positive scalar
%   - sigma: real in (-Inf, 1)
%   - tau: positive scalar
%   - Fdist:  structure; specifies the distribution of the beta vectors
%   - gamma: positive vector of size p
% -------------------------------------------------------------------------
% EXAMPLE
% p = 3; hyper_alpha = [100, 1]; hyper_sigma = .5; hyper_tau = [1, 0.1];
% hyper_Fdist.name = 'gamma'; hyper_Fdist.param = [.01, 1]; gamma = zeros(p, 1);
% N = CGGPgetparams(p, hyper_alpha, hyper_sigma, hyper_tau, hyper_Fdist, gamma);

% Copyright (c) F. Caron (University of Oxford), A. Todeschini (Inria), and 
% X. Miscouridou (University of Oxford)
% caron@stats.ox.ac.uk
% adrien.todeschini@gmail.com
% xenia.miscouridou@spc.ox.ac.uk
% September 2017
%--------------------------------------------------------------------------

if nargin<7
    mode = 'rnd';
end

[alpha, sigma, tau] = GGPgetparams(alpha, sigma, tau, mode, varargin{:});

if size(gamma, 2)==2 % gamma is missing
    hyper_gamma = gamma;
    if any(hyper_gamma(:,1)<0) || any(hyper_gamma(:,2)<0)
        error('Hyperparameters for gamma must be nonnegative');
    end
    switch mode
        case 'default'
            gamma = zeros(p,1);
        case 'empty'
            gamma = [];
        case 'init'
            gamma = zeros(p,1);  %%%% TEMP
        case 'rnd'
            gamma = gamrnd(hyper_gamma(:,1), 1./hyper_gamma(:,2));
    end
elseif size(gamma, 2)~=1
    error('gamma must have either one or two columns')
end

switch(Fdist.name)
    case 'gamma'
        b = [];
        if isnumeric(Fdist.param)
            a = Fdist.param;
        else
            a = Fdist.param.a;
            b = Fdist.param.b;
        end
        if size(b,2)==2 % b is missing
            hyper_b = b;
            if any(hyper_b(:,1)<0) || any(hyper_b(:,2)<0)
                error('Hyperparameters for Fdist.param must be nonnegative');
            end
            switch mode
                case 'default'
                    b = 1/p*ones(size(hyper_b,1),1);
                case 'empty'
                    b = [];
                case 'init'
                    b = sort(gamrnd(1,1,[size(hyper_b,1),1]));
                    b = b/sum(b);
                case 'rnd'
                    b = gamrnd(hyper_b(:,1), 1./hyper_b(:,2));
            end
            
        end
        if size(a,2)==2 % a is missing
            hyper_a = a;
            if any(hyper_a(:,1)<0) || any(hyper_a(:,2)<0)
                error('Hyperparameters for Fdist.param must be nonnegative');
            end
            switch mode
                case 'default'
                    a = fillmat(1/sqrt(p).*b, [size(hyper_a,1), 1]);
                case 'empty'
                    a = [];
                case 'init'
%                     a = fillmat(1/sqrt(p).*b, [size(hyper_a,1), 1]);
                    a = 2*rand(size(hyper_a,1), 1).*fillmat(1/sqrt(p).*b, [size(hyper_a,1), 1]);
                case 'rnd'
                    a = gamrnd(hyper_a(:,1), 1./hyper_a(:,2));
            end
        end
        if isnumeric(Fdist.param)
            Fdist.param = a;
        else
            Fdist.param.a = a;
            Fdist.param.b = b;
        end
    otherwise
        error('Unknown Distribution %s', Fdist.name)
end
