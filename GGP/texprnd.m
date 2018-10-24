function samples = texprnd(lambda, a, n, m)

% Samples from a left truncated exponential distribution
%   - m: positive integer
%   - n: positive integer
%   - lambda: matrix of size [n,m]
%   - a: positive real


samples = -log( - rand(n, m).*exp(-lambda*a) + exp(-lambda*a) )./lambda;

end