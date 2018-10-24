function [S_w_small, mu, Var] = CGGPsumsmallrnd(p, logalpha, sigma, tau, a, b, gamma, T, varargin)

% CGGPsumsmallrnd implements a Gaussian approximation  (N(mu,var)) for the sum of
% weights below a threshold T in a (tilted) Compound Generalized Gamma Process in which
% atom locations are restricted on [0,alpha] and 
% the intensity dictating the jumps depends on hyperparameters sigma, tau,
% a, b, gamma
% -------------------------------------------------------------------------
% INPUTS
%   - p: positive integer
%   - logalpha: scalar
%   - a: scalar or vector of size [1,p]
%   - b: scalar or vector of size [1,p]
%   - gamma: vector of size [1,p]
%   - T: scalar, weight threshold
%
% OPTIONAL INPUTS
%   - mufun: char that specifies the method to calculate the mean in the Gaussian approximation
%   - Varfun: char that specifies the method to calculate the variance in the Gaussian approximation
%   - cutoffs: vector with cutoffs relating to the numerical evaluations
%    involved in the calculations of the mean and variance
% ----------------------------------------------------------------------
% OUTPUTS
%   - S_w_small: sum of weights below T
%   - mu: mean of the Gaussian approximation
%   - Var: variance of the Gaussian approximation
% -------------------------------------------------------------------------

% Copyright (c) F. Caron (University of Oxford), A. Todeschini (Inria), and 
% X. Miscouridou (University of Oxford)
% caron@stats.ox.ac.uk
% adrien.todeschini@gmail.com
% xenia.miscouridou@spc.ox.ac.uk
% September 2015


%% parse input arguments
ip = inputParser;
ip.KeepUnmatched = true;

addOptional(ip, 'mufun', 'mu2', @ischar);
addOptional(ip, 'Varfun', 'Var1', @ischar);
addOptional(ip, 'cutoffs', 1e-1, @isnumeric);
parse(ip, varargin{:});

mufun = ip.Results.mufun;
Varfun = ip.Results.Varfun;
cutoffs = ip.Results.cutoffs;

varargin = varargin((end-numel(fieldnames(ip.Unmatched))+1):end);

cutoffs = sort(cutoffs, 'descend');
cutoffs = cutoffs(cutoffs<T);

if numel(a)~=p
    a = a*ones(1,p);
end
if numel(b)~=p
    b = b*ones(1,p);
end
if numel(gamma)~=p
    gamma = gamma*ones(1,p);
end

%% Compute the mean
switch mufun
    case 'mu0'
        mu = zeros(1,p);
    case 'mu1'
        mu = mu1_fun(logalpha, sigma, tau, a, b, gamma, T, cutoffs, varargin{:});
    case 'mu2'
        mu = mu2_fun(logalpha, sigma, tau, a, b, gamma, T, cutoffs, varargin{:});
    case 'mu3'
        mu = mu3_fun(logalpha, sigma, tau, a, b, gamma, T);
    otherwise
        error('Unknown function')
end

%% Compute the covariance
switch Varfun
    case 'Var0'
        Var = zeros(p);
    case 'Var1'
        % numerical integration
        Var = Var1_fun(logalpha, sigma, tau, a, b, gamma, T, cutoffs, varargin{:});
    otherwise
        error('Unknown function')
end

%% Sample
[~,err] = cholcov(Var);
if err~=0
    S_w_small = max(0, mu);
else
    S_w_small = max(0, mvnrnd(mu, Var, 1));
end
end

%% Mean

% integrand
function out = rho_mu(x, sigma, tau, a, logb, gamma_b)
log1p_gbx = log1p(gamma_b*x);
out = exp(bsxfun(@minus, bsxfun(@minus, -sigma*log(x)-tau*x, bsxfun(@plus, logb, log1p_gbx)),...
    sum(bsxfun(@times, a, log1p_gbx), 1)));
end

function mu = mu1_fun(logalpha, sigma, tau, a, b, gamma, T, cutoffs, varargin)
% numerical integration above cutoffs
rho_mu_fun = @(x) rho_mu(x, sigma, tau, a', log(b'), (gamma./b)');
mu = 0;
cutoff = T;
for i=1:numel(cutoffs)
    mu = mu + integral(rho_mu_fun, cutoffs(i), cutoff, 'ArrayValued', true, varargin{:})'; %%%% NOTE: very expensive
    if any(isinf(mu(:))) || any(isnan(mu(:))) || any(imag(mu(:))~=0)
        error('numerical integration invalid at %g', cutoffs(i));
    end
    cutoff = cutoffs(i);
end

% asymtotic limit below cutoffs
log_cst = log(a)+logalpha-gammaln(1-sigma);
% mu = exp(log_cst) .* mu + exp(log_cst-log(b)+(1-sigma)*log(cutoff)-log(1-sigma));
mu = exp(log_cst) .* mu + exp(log_cst-log(b)) .* (cutoff^(1-sigma)/(1-sigma)-tau*cutoff^(2-sigma)/(2-sigma)); % with Taylor expansion of exp
end

%% Mean, using integral by parts

% integrand
function out = rho_mu2(x, onemsigma, tau, a, logb, gamma_b)
gbx = gamma_b*x;
log1p_gbx = log1p(gbx);
logout1 = bsxfun(@minus, onemsigma*log(x)-tau*x-logb, log1p_gbx) - sum(bsxfun(@times, a, log1p_gbx), 1);

gb_1pgbx = bsxfun(@rdivide, gamma_b, 1+gbx);
logout2 = log(tau + gb_1pgbx + sum(bsxfun(@times, a, gb_1pgbx)));

out = exp(logout1 + logout2);
end

function mu = mu2_fun(logalpha, sigma, tau, a, b, gamma, T, cutoffs, varargin)
gamma_b = (gamma./b);
rho_mu_fun = @(x) rho_mu2(x, 1-sigma, tau, a', log(b'), gamma_b');
mu = 0;
cutoff = T;
for i=1:numel(cutoffs)
    mu = mu + integral(rho_mu_fun, cutoffs(i), cutoff, 'ArrayValued', true, varargin{:})';
    if any(isinf(mu(:))) || any(isnan(mu(:))) || any(imag(mu(:))~=0)
        error('numerical integration invalid at %g', cutoffs(i));
    end
    cutoff = cutoffs(i);
end
mu = mu + integral(rho_mu_fun, 0, cutoff, 'ArrayValued', true, varargin{:})';

A = exp((1-sigma)*log(T)-log(1-sigma)-tau*T-log(b)-log1p(gamma_b*T)-sum(a.*log1p(gamma_b*T)));
mu = exp(log(a)+logalpha-gammaln(1-sigma)) .* ( A + 1/(1-sigma) * mu);
end

%% Mean, using integral by parts and neglecting T^(3-sigma) terms

function mu = mu3_fun(logalpha, sigma, tau, a, b, gamma, T)
gamma_b = (gamma./b);
log1p_gbT = log1p(gamma_b*T);
logA = (1-sigma)*log(T)-log(1-sigma)-tau*T-log(b)-log1p_gbT-sum(a.*log1p_gbT);
gb_1pgbT = gamma_b./(1+gamma_b*T);
B = T/(2-sigma) * (tau+gb_1pgbT +sum(a.*gb_1pgbT));
mu = exp(log(a)+logalpha-gammaln(1-sigma) + logA + log1p(B));
end

%% Variance

% integrand
function out = rho_var1(x, onemsigma, tau, a, logb, gamma_b)
log1p_gbx = log1p(gamma_b*x);
log_l1pgbx = bsxfun(@plus, logb, log1p_gbx);
out = exp(bsxfun(@minus, bsxfun(@minus, onemsigma*log(x)-tau*x, bsxfun(@plus, repmat(log_l1pgbx, [1,1,numel(gamma_b)]), shiftdim(log_l1pgbx', -1))),...
    sum(bsxfun(@times, a, log1p_gbx), 1)));
end

function Var = Var1_fun(logalpha, sigma, tau, a, b, gamma, T, cutoffs, varargin)
rho_var_fun = @(x) rho_var1(x, 1-sigma, tau, a', log(b'), (gamma./b)');

Var = 0;
cutoff = T;
for i=1:numel(cutoffs)
    Var = Var + squeeze(integral(rho_var_fun, cutoffs(i), cutoff, 'ArrayValued', true, varargin{:}));
    if any(isinf(Var(:))) || any(isnan(Var(:))) || any(imag(Var(:))~=0)
        error('numerical integration invalid at %g', cutoffs(i));
    end
    cutoff = cutoffs(i);
end

Var = Var + squeeze(integral(rho_var_fun, 0, cutoff, 'ArrayValued', true, varargin{:}));

Var = exp(log(a'*a+diag(a))+logalpha-gammaln(1-sigma)) .* Var;
end
