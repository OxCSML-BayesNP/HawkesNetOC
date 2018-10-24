function x = tpoissrnd(lambda)

% TPOISSRND Sample from a zero-truncated Poisson distribution
% lambda: scalar, vector or matrix 
%
x = ones(size(lambda));
ind = (lambda > 1e-5); % below this value, x=1 w. very high proba
lambda_ind = lambda(ind);
x(ind) = poissinv(exp(-lambda_ind) +rand(size(lambda_ind)).*(1 - exp(-lambda_ind)), lambda_ind);
end
