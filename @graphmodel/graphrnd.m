function [G, varargout] = graphrnd(obj, varargin)

% graphrnd samples a graph from a graph model
% [G, varargout] = GRAPHRND(obj, varargin)
%
% -------------------------------------------------------------------------
% INPUTS
%   - obj: an object of the graphmodel class
%
% Optional inputs:
%   - T: truncation threshold (for a GGP or CGGP graph model)
%   - maxiter: maximum number of iterations for the adaptive thinning
%     strategy (default=1e8)
% -------------------------------------------------------------------------
% OUTPUTS
%   - G: sparse logical adjacency matrix
%
% Specific outputs:
% For GGP non bipartite:
%   G = GRAPHRND(obj)
%   G = GRAPHRND(obj, T)
%   G = GRAPHRND(obj, T, maxiter)
%   [G, w, w_rem, alpha, sigma, tau, T] = GRAPHRND(__)
%   - w: sociability parameters of nodes with at least one connection
%   - w_rem: sum of the sociability parameters of nodes with no connection
%   - alpha: Parameter alpha of the GGP
%   - sigma: Parameter sigma of the GGP
%   - tau: Parameter tau of the GGP
%   - T: threshold
% For GGP bipartite:
%   G = GRAPHRND(obj)
%   G = GRAPHRND(obj, T)
%   G = GRAPHRND(obj, T, maxiter)
%   [G, w1, w1_rem, w2, w2_rem, alpha1, sigma1, tau1, alpha2, sigma2, tau2, T1, T2] = GRAPHRND(__)
%   - w1: sociability parameters of type 1 nodes with at least one connection
%   - w1_rem: sum of the sociability parameters of type 1 nodes with no connection
%   - w2: sociability parameters of type 2 nodes with at least one connection
%   - w2_rem: sum of the sociability parameters of type 2 nodes with no connection
%   - alpha1: Parameter alpha of the GGP for type 1 nodes
%   - sigma1: Parameter sigma of the GGP for type 1 nodes
%   - tau1: Parameter tau of the GGP for type 1 nodes
%   - alpha2: Parameter alpha of the GGP for type 2 nodes
%   - sigma2: Parameter sigma of the GGP for type 2 nodes
%   - tau2: Parameter tau of the GGP for type 2 nodes
%   - T1: threshold for type 1 nodes
%   - T2: threshold for type 2 nodes
% For CGGP non bipartite:
%   G = GRAPHRND(obj)
%   G = GRAPHRND(obj, T)
%   G = GRAPHRND(obj, T, maxiter)
%   [G, w, w_rem, alpha, sigma, tau, Fdist, gamma, T] = GRAPHRND(__)
%   - w: Matrix of size [n,p]. Feature interest parameters of nodes with at least one connection
%   - w_rem: Vector of length p. Sums of the feature interest parameters of nodes with no connection
%   - alpha: Parameter alpha of the GGP
%   - sigma: Parameter sigma of the GGP
%   - tau: Parameter tau of the GGP
%   - Fdist: Structure. Parameters of the p-dimensional distribution F of the scores
%     of the compound CRM for type 1 nodes
%   - gamma: Vector of length p. Tilting parameters gamma of the compound CRM
%   - T: threshold
% For CGGP bipartite:
%   [G, w1, w1_rem, w2, w2_rem, alpha1, sigma1, tau1, Fdist1, gamma1, alpha2, sigma2, tau2, Fdist2, gamma2] = GRAPHRND(obj, T)
%   - w1: Matrix of size [n,p]. Feature interest parameters of type 1 nodes with at least one connection
%   - w1_rem: Vector of length p. Sums of the feature interest parameters of type 1 nodes with no connection
%   - w2: Matrix of size nxp. Feature interest parameters of type 2 nodes with at least one connection
%   - w2_rem: Vector of length p. Sums of the feature interest parameters of type 2 nodes with no connection
%   - alpha1: Parameter alpha of the GGP for type 1 nodes
%   - sigma1: Parameter sigma of the GGP for type 1 nodes
%   - tau1: Parameter tau of the GGP for type 1 nodes
%   - Fdist1: Structure. Parameters of the p-dimensional distribution F of the scores beta
%     of the compound CRM for type 1 nodes
%   - gamma1: Vector of length p. Tilting parameters gamma of the compound CRM for type 1 nodes
%   - alpha2: Parameter alpha of the GGP for type 2 nodes
%   - sigma2: Parameter sigma of the GGP for type 2 nodes
%   - tau2: Parameter tau of the GGP for type 2 nodes
%   - Fdist2: Structure. Parameters of the p-dimensional distribution F of the scores beta
%     of the compound CRM for type 2 nodes
%   - gamma2: Vector of length p. Tilting parameters gamma of the compound CRM for type 2 nodes
%   - T1: threshold for type 1 nodes
%   - T2: threshold for type 2 nodes
% -------------------------------------------------------------------------
%
% See also GRAPHMODEL, GRAPHMODEL/GRAPHMODEL

% Copyright (c) F. Caron (University of Oxford), A. Todeschini (Inria), and 
% X. Miscouridou (University of Oxford)
% caron@stats.ox.ac.uk
% adrien.todeschini@gmail.com
% xenia.miscouridou@spc.ox.ac.uk
% September 2015
% --------------------------------------------------------------------------

varargout = [];
switch(obj.type)
    case 'ER'
        G = ERgraphrnd(obj.param.n, obj.param.p);
        varargout = {};
    case 'GGP'
        if nargin>2
            error('Only 2 possible inputs for objects of type GGP')
        end
        switch(obj.typegraph)
            case 'undirected'
                [G, w, w_rem, alpha, sigma, tau, T] = GGPgraphrnd(...
                    obj.param.alpha, obj.param.sigma, obj.param.tau, ...
                    varargin{:});
                varargout = {w, w_rem, alpha, sigma, tau, T};
            case 'simple'
                [G, w, w_rem, alpha, sigma, tau, T] = GGPgraphrnd(...
                    obj.param.alpha, obj.param.sigma, obj.param.tau, ...
                    varargin{:});
                varargout = {w, w_rem, alpha, sigma, tau, T};
                G = G - diag(diag(G));
            case 'bipartite'
                [G, w1, w1_rem, w2, w2_rem, alpha1, sigma1, tau1, alpha2, sigma2, tau2, T1, T2] ...
                    = GGPbipgraphrnd(obj.param(1).alpha, obj.param(1).sigma, obj.param(1).tau, ...
                    obj.param(2).alpha, obj.param(2).sigma, obj.param(2).tau, varargin{:});
                varargout = {w1, w1_rem, w2, w2_rem, alpha1, sigma1, tau1, ...
                    alpha2, sigma2, tau2, T1, T2};
        end
    case 'CGGP'
        if nargin>2
            error('Only 2 possible inputs for objects of type CGGP')
        end
        switch(obj.typegraph)
            case 'undirected'
                [G, w, w_rem, alpha, sigma, tau, Fdist, gamma, T, w0, beta] = ...
                    CGGPgraphrnd(obj.param.p, obj.param.alpha, obj.param.sigma, ...
                    obj.param.tau, obj.param.Fdist, obj.param.gamma, ...
                    obj.param.observe_all, obj.param.necessarily_sparse, varargin{:});
                
                varargout = {w, w_rem, alpha, sigma, tau, Fdist, gamma, T, w0, beta};
            case 'simple'
                [G, w, w_rem, alpha, sigma, tau, Fdist, gamma, T, w0, beta] = ...
                    CGGPgraphrnd(obj.param.p, obj.param.alpha, obj.param.sigma, ...
                    obj.param.tau, obj.param.Fdist, obj.param.gamma, ...
                    obj.param.observe_all, obj.param.necessarily_sparse, varargin{:});
                
                varargout = {w, w_rem, alpha, sigma, tau, Fdist, gamma, T, w0, beta};
                
                G = G - diag(diag(G));
            case 'bipartite'
                [G, w1, w1_rem, w2, w2_rem, alpha1, sigma1, tau1, Fdist1, ...
                    gamma1, alpha2, sigma2, tau2, Fdist2, gamma2, T1, T2, ...
                    w01, w02, beta1, beta2] = CGGPbipgraphrnd(obj.param(1).p, ...
                    obj.param(1).alpha, obj.param(1).sigma, obj.param(1).tau, ...
                    obj.param(1).Fdist, obj.param(1).gamma, obj.param(2).alpha, ...
                    obj.param(2).sigma, obj.param(2).tau, obj.param(2).Fdist, ...
                    obj.param(2).gamma, [obj.param.observe_all], [obj.param.necessarily_sparse], varargin{:});
                
                varargout = {w1, w1_rem, w2, w2_rem, alpha1, sigma1, tau1,...
                    Fdist1, gamma1, alpha2, sigma2, tau2, Fdist2, gamma2, ...
                    T1, T2, w01, w02, beta1, beta2};
        end
    case 'BA'
        G = BAgraphrnd(obj.param.n);
    case 'Lloyd'
        G = Lloydgraphrnd(obj.param.n, obj.param.sig, obj.param.c, obj.param.d);
    case 'MMSB'
        [G, s, pi] = mmsbrnd(obj.param.n, obj.param.alpha, obj.param.W, obj.param.rho);
        G = sparse(logical(G));
        varargout = {s, pi}; 
end

end
