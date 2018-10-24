function lp =  logpriorcggp(sample, param)

% logpriorcggp computes the log prior for the model parameters and
% hyperparameters under the model prior CGGP
%
% INPUT
%   - sample: structure; MCMC sample from the posterior 
%   - param: struct with the hyperparameter values
% OUTPUT
%   - lp: log prior

 
w = max(sample.w, 1e-30);
if isfield(sample, 'alpha')
    alpha = sample.alpha;
else
    alpha = param.alpha;
end
if isfield(sample, 'sigma')
    sigma = sample.sigma;
else
    sigma = param.sigma;
end
if isfield(sample, 'tau')
    tau = sample.tau;
else
    tau = param.tau;
end
if isfield(sample, 'gamma')
    gamma = sample.gamma;
else
    gamma = param.gamma;
end
switch(param.Fdist.name)
    case 'gamma'
        if isfield(sample, 'Fparam')
            if isnumeric(sample.Fparam)
                a = sample.Fparam;
                b = a;
            else
                a = sample.Fparam.a;
                b = sample.Fparam.b;
            end
        else
            if isnumeric(param.Fdist.param)
                a = param.Fdist.param;
                b = a;
            else
                a = param.Fdist.param.a;
                b = param.Fdist.param.b;
            end
        end
    otherwise
        error('Unknown Distribution %s', param.Fdist.name)
end

[K, p] = size(w);

% prior of weights
if length(b)==p
    sum_bw1 = w * b;
else
    sum_bw1 = sum(w * b, 2);
end

lp = sum(sum(bsxfun(@times, a'-1, log(w))))...
    - (sigma + sum(a))/2 * sum(log(sum_bw1));
z = 2 * sqrt(tau * sum_bw1);
temp = log(besselk(sigma+sum(a), z, 1)) - z;
% For small arguments
issmall = isinf(temp);
temp(issmall) = gammaln(abs(sigma+sum(a)))-log(2) + abs(sigma+sum(a))*(log(2)-log(z(issmall)));

lp = lp + sum(sum(temp));
lp = lp + K*sum(a .* log(b))- K*sum(gammaln(a)) + K*(sigma + sum(a))/2*log(tau);
lp = lp + K*log(-sigma) - K*gammaln(1-sigma) - K*sigma *log(tau);

%% p(N_alpha| rest)
lam = -alpha/sigma*tau^sigma;
lp = lp + K*log(lam)-lam-gammaln(K);

%% hyperparameters

if numel(param.alpha)==2
    hyper = param.alpha;
    if all(hyper>0)
        lp = lp + gamlogpdf(alpha, hyper(1), hyper(2));
    end
end

if numel(param.sigma)==2
    hyper = param.sigma;
    if param.observe_all
        if all(hyper>0)
            lp = lp + betalogpdf(sigma, hyper(1), hyper(2));
        end
    elseif param.infinite
        if all(hyper>0)
            lp = lp + gamlogpdf(-sigma, hyper(1), hyper(2));
        end
    else
        if all(hyper>0)
            lp = lp + gamlogpdf(1-sigma, hyper(1), hyper(2));
        end
    end
end

if numel(param.sigma)==2
    hyper = param.sigma;
    if param.observe_all
        if sigma>=0
            error('sigma>=0 and observe_all=true not allowed')
        end
        if all(hyper>0)
            lp = lp + betalogpdf(sigma, hyper(1), hyper(2));
        end
    elseif param.infinite
        if sigma<0
            error('sigma<0 and infinite=true not allowed')
        end
        if all(hyper>0)
            lp = lp + gamlogpdf(-sigma, hyper(1), hyper(2));
        end
    else
        if all(hyper>0)
            lp = lp + gamlogpdf(1-sigma, hyper(1), hyper(2));
        end
    end
end

if numel(param.tau)==2
    hyper = param.tau;
    if all(hyper>0)
        lp = lp + gamlogpdf(tau, hyper(1), hyper(2));
    end
end

if size(param.gamma, 2)==2
    hyper_gamma = param.gamma;
    if all(hyper>0)
        lp = lp + sum(gamlogpdf(gamma, hyper_gamma(:,1), hyper_gamma(:,2)));
    end
end

switch(param.Fdist.name)
    case 'gamma'
        if isnumeric(param.Fdist.param)
            if size(param.Fdist.param, 2)==2
                hyper = param.Fdist.param;
                if all(hyper>0)
                    lp = lp + sum(gamlogpdf(a, hyper(:,1), hyper(:,2)));
                end
            end
        else
            if size(param.Fdist.param.a, 2)==2
                hyper = param.Fdist.param.a;
                if all(hyper>0)
                    lp = lp + sum(gamlogpdf(a, hyper(:,1), hyper(:,2)));
                end
            end
            if size(param.Fdist.param.b, 2)==2
                hyper = param.Fdist.param.b;
                if all(hyper>0)
                    lp = lp + sum(gamlogpdf(b, hyper(:,1), hyper(:,2)));
                end
            end
        end
    otherwise
        error('Unknown Distribution %s', param.Fdist.name)
end

end


function out = gamlogpdf(x,a,b)
out = -b.*x;
if x>0
    out = out + (a-1).*log(x);
end
end

function out = betalogpdf(x,a,b)
out = 0;
if x>0
    out = out + (a-1).*log(x);
end
if x<1
    out = out + (b-1).*log(1-x);
end
end
