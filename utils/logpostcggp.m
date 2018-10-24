function [lp, ll] = logpostcggp(G, sample, prior)

% logpostcggp computes the log posterior p(w|G) (up to a constant)
% which is the sum of the log likelihood p(G|w) and the log prior p(w) (up to a constant)


% INPUTS:
%   - G: observed binary adjacency matrix
%   - sample: structure; sample from the MCMC posterior distribution
%   - prior: object of class graphmodel; 
% ---------------------------------------------
% OUTPUT:
%   - lp: logposterior
%   - ll: loglikelihood 


% Loglikelihood p(G|w)
ll = loglikcggp(G, sample, prior);
lp = ll;
% Prior
for i=1:length(sample)
    lp = lp + logpriorcggp(sample(i), prior.param(i)); 
end


