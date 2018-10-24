function out = CGGPpsi(t, sigma, tau, Fdist, gamma, varargin)

% CGGPPSI returns the Laplace exponent of a CGGP.
% out = CGGPpsi(t, sigma, tau, Fdist, gamma)
%
% The Laplace exponent of a CGGP evaluated at t=(t_1,...,t_p) is
% psi(t_1,...,t_p) = int ( 1-exp(-t_1*w_1-...-t_p*w_p ) \rho(dw_1,...,dw_p))
% which under the CGGP model with Gamma distributed scores WE HAVE
% \rho(dw_1,...,dw_p)= exp(-gamma_1*w_1,...,gamma_p*w_p) int w_0^{-p} F(w_1/w_0,...,w_p/w_0) rho_0(dw_0)
% rho_0 is the measure of the which admites the form
% rho(dw) = alpha/Gamma(1-sigma) * w^{-1-sigma} * exp(-tau*w)dw
% -------------------------------------------------------------------------
% INPUTS
%   - t: vector of positive scalars of length p
%   - alpha: positive scalar
%   - sigma: real in (-inf, 1)
%   - tau: positive scalar
%   - Fdist: structure that defines the scores distribution F
%   - gamma: positive tilting parameters
%
% Further optional inputs are passed to the 'integral' matlab function.
%
% OUTPUT
%   - out: Laplace exponent evaluated at the values t
% -------------------------------------------------------------------------
% EXAMPLE
% t = .1:.1:10;
% p=length(t);
% alpha = 100; sigma = 0.5; tau = 1; gamma = zeros(p,1);
% Fdist.name = 'gamma'; Fdist.param = 0.01*ones(p,1);
% out = CGGPpsi(t, sigma, tau, Fdist, gamma);

% Copyright (c) F. Caron (University of Oxford), A. Todeschini (Inria), and 
% X. Miscouridou (University of Oxford)
% caron@stats.ox.ac.uk
% adrien.todeschini@gmail.com
% xenia.miscouridou@spc.ox.ac.uk
% September 2017
%--------------------------------------------------------------------------

if any(gamma>0)
    error('This function deal only for a CGGP with zero tilting parameters, i.e. gamma=(0,...,0)')
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
end

if abs(sigma)<1e-8 % gamma process
    sigma = (sigma>0)*1e-8; %%%% truncate
end



if numel(a)==1 && isinf(a)
    return
end

t_b = (t./b)';
rho_fun = @(x) rho2(x, t_b, sigma, tau, a');
out = exp(-gammaln(2-sigma))/sigma*integral(rho_fun, 0, Inf, varargin{:}) - tau^sigma/sigma;

if imag(out)~=0 || out<0 || isnan(out)
    keyboard
end
end


% %% stable for all
function out = rho(x, t_b, sigma, tau, a)
out =  exp( - x.*tau - sum(bsxfun(@times, a, log1p(t_b*x)), 1)  + (-sigma)*log(x)) .* (tau +sum(bsxfun(@times, a.*t_b, (1+t_b*x).^(-1) ), 1));
out_nan = isnan(out);
x_pos = x>0;
ind_naninf = x_pos & (out_nan | isinf(out));
if any(ind_naninf)
    keyboard
end
end

function out = rho2(x, t_b, sigma, tau, a)

out =  exp( - x.*tau - sum(bsxfun(@times, a, log1p(t_b*x)), 1)  + (1-sigma)*log(x))...
    .* ((tau +sum(bsxfun(@times, a.*t_b, (1+t_b*x).^(-1) ), 1)).^2  + sum(bsxfun(@times, a.*(t_b.^2), (1+t_b*x).^(-2) ), 1)  );
out_nan = isnan(out);
x_pos = x>0;
ind_naninf = x_pos & (out_nan | isinf(out));
if any(ind_naninf)
    keyboard
end
end
