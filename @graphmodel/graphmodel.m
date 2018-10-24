classdef graphmodel
    % graphmodel class containing the parameters of the graph model
    %
    % PROPERTIES
    %   - <a href="matlab: help graphmodel/name">name</a>
    %   - <a href="matlab: help graphmodel/type">type</a>
    %   - <a href="matlab: help graphmodel/param">param</a>
    %   - <a href="matlab: help graphmodel/typegraph">typegraph</a>
    %
    % CLASS CONSTRUCTOR
    %   - <a href="matlab: help graphmodel/graphmodel">graphmodel</a>
    %
    % METHODS
    %   - <a href="matlab: help graphmodel/graphrnd">graphrnd</a>: samples a graph

    properties
        name; % name of the graph model (string)
        type; % type of graph (string)
        param; % parameters (structure)
        typegraph = 'undirected'; % undirected, simple (undirected, no loops, no multiple edges), bipartite
    end
    %---------------------------------------------------------------------
    methods

        %% Class constructor
        function obj = graphmodel(type, varargin)
            %GRAPHMODEL creates a graph model.
            % obj = GRAPHMODEL(type, varargin)
            %
            % Possible type: 'ER', 'GGP', 'CGGP', 'BA', 'Lloyd', 'MMSB'
            % -------------------------------------------------------------
            % obj = GRAPHMODEL('ER', n, p) creates an Erdos-R?nyi graph
            % model with n nodes and probability of connection p in [0,1]
            % -------------------------------------------------------------
            % obj = GRAPHMODEL('GGP', alpha, sigma, tau, typegraph) creates
            % a GGP graph model where
            %   - optional input typegraph can be undirected (default),
            %     simple or bipartite
            %   - If typegraph is undirected or simple:
            %     - alpha: double or vector of length 2.
            %       In the first case, alpha is supposed to be fixed.
            %       Otherwise, it is random and drawn from a gamma
            %       distribution of parameters alpha(1) and alpha(2)
            %     - sigma: double or vector of length 2.
            %       In the first case, sigma is supposed to be fixed.
            %       Otherwise, it is random and (1-sigma) is drawn from
            %       a gamma distribution of parameters sigma(1) and sigma(2)
            %     - tau: double or vector of length 2.
            %       In the first case, tau is supposed to be fixed.
            %       Otherwise, it is random and drawn from a gamma
            %       distribution of parameters tau(1) and tau(2)
            %   - If typegraph is bipartite:
            %     - alpha: cell of length 2.
            %       - alpha{1} corresponds to the parameter alpha of
            %         the GGP associated to the first type of node.
            %         It may be a double or a vector of length 2.
            %         (See above for details)
            %       - alpha{2} corresponds to the parameter alpha of
            %         the GGP associated to the second type of node.
            %         It may be a double or a vector of length 2.
            %         (See above for details)
            %     - sigma: cell of length 2. Same as above.
            %     - tau: cell of length 2. Same as above.
            % -------------------------------------------------------------
            % obj = GRAPHMODEL('CGGP', p, alpha, sigma, tau, Fdist, gamma, 
            %       typegraph, observe_all, necessarily_sparse) 
            %       creates a Compound GGP graph model where
            %   - p is the number of features
            %   - optional string input typegraph can be 'undirected' (default),
            %     'simple' or 'bipartite'
            %   - If typegraph is undirected or simple:
            %     - alpha: double or vector of length 2.
            %       In the first case, alpha is supposed to be fixed.
            %       Otherwise, it is random and drawn from a gamma
            %       distribution of parameters alpha(1) and alpha(2)
            %     - sigma: double or vector of length 2.
            %       In the first case, sigma is supposed to be fixed.
            %       Otherwise, it is random and (1-sigma) is drawn from
            %       a gamma distribution of parameters sigma(1) and sigma(2)
            %     - tau: double or vector of length 2.
            %       In the first case, tau is supposed to be fixed.
            %       Otherwise, it is random and drawn from a gamma
            %       distribution of parameters tau(1) and tau(2)
            %     - Fdist: distribution F of the p-dimensional scores beta.
            %         A struct with fields
            %         - name: name of the distribution. Possible values:
            %           'gamma': product of independant gamma distributions
            %         - param: parameters of the distribution.
            %           Case name='gamma': if param is numeric then the shape and
            %           rate parameters are equal, else if param is a struct
            %           it must contain fields a and b.
            %           All parameters values are either one or two columns matrices.
            %           In the first case, the parameter is supposed to be fixed.
            %           Otherwise, its p components are random and drawn from a
            %           Gamma(first column, second columns).
            %           The number of rows is either one or p. If one row, the
            %           p components are supposed equal.
            %     - gamma: one or two columns matrix.
            %       In the first case, gamma is supposed to be fixed.
            %       Otherwise, its p components are random and drawn from a gamma
            %       distribution of parameters gamma(:,1) and gamma(:,2).
            %       If empty, gamma is fixed to zero.
            %     - observe_all: logical (default=false); if set to true,
            %     indicates that all nodes are observed, including nodes
            %     with connection; necesarily in the sparse regime (\sigma<0)
            %     - necessarily_sparse: logical (default=false)
            %       By default, sigma can take values in (-infty,1) and
            %       cover dense and sparse regime; if
            %       necessarily_sparse=true, it forces the model to be sparse (1>sigma>0)
            %   - If typegraph is bipartite:
            %     - alpha: cell of length 2.
            %       - alpha{1} corresponds to the parameter alpha of
            %         the GGP associated to the first type of node.
            %         (See above for details)
            %       - alpha{2} corresponds to the parameter alpha of
            %         the GGP associated to the second type of node.
            %         (See above for details)
            %     - sigma: cell of length 2. Same as above.
            %     - tau: cell of length 2. Same as above.
            %     - Fdist: struct array of length 2. Same as above.
            %     - gamma: cell of length 2. Same as above.
            % -------------------------------------------------------------
            % obj = GRAPHMODEL('BA', n) creates a Barabasi-Albert graph
            % model with n nodes
            % -------------------------------------------------------------
            % obj = GRAPHMODEL('Lloyd', n, sig, c, d) creates a Lloyd graph
            % model with parameters n, sig, c and d
            % -------------------------------------------------------------
            % obj = GRAPHMODEL('MMSB', n, p, alpha, W, rho, epsilon) creates a (MMSB)
            % mixed membership stochastic blockmodel with parameters n, p, alpha, W, rho
            %   - n: number of nodes
            %   - p: number of features
            %   - alpha: Dirichlet parameter for sampling the node membership vectors
            %   - W: matrix of probabilities for group interactions
            %   - rho: sparsity parameter
            %   - epsilon: optional parameter that is part of an extension of the MMSB,
            %     equals the small probability of connection when node
            %     indicators point to different communities
            % -------------------------------------------------------------
            % EXAMPLES
            % n = 1000; p = 0.01;
            % obj = graphmodel('ER', n, p);
            % alpha = 100; sigma = 0.1; tau = 1;
            % obj2 = graphmodel('GGP', alpha, sigma, tau);
            % obj3 = graphmodel('BA', n);
            % hyper_alpha = [100, 1]; sigma = 0.1; tau = 1;
            % obj4 = graphmodel('GGP', hyper_alpha, sigma, tau);
            % alpha = {100,50}; sigma = {-1,.5}; tau = {[1,.1],1};
            % obj5 = graphmodel('GGP', alpha, sigma, tau, 'bipartite');
            % p = 2; alpha = 20; sigma = .5; tau = 1;
            % Fdist = struct('name', 'gamma', 'param', struct('a', 0.2, 'b', 1/p));
            % gamma = zeros(p,1);
            % obj6 = graphmodel('CGGP', p, alpha, sigma, tau, Fdist, gamma);
            % p = 3; alpha = {100,50}; sigma = {-1,.5}; tau = {1, 1};
            % Fdist2(1) = struct('name', 'gamma', 'param', struct('a', 0.2, 'b', 1/p));
            % Fdist2(2) = struct('name', 'gamma', 'param', struct('a', 0.2, 'b', 1/p));
            % gamma = zeros(p,1);
            % obj7 = graphmodel('CGGP', p, alpha, sigma, tau, Fdist2, gamma, 'bipartite');
            % n = 100; p = 4; alpha = 2; W = 0.5*ones(p,p);; rho = 0.02;
            % obj8 = graphmodel('MMSB',n, p, alpha, W, rho);

            % Copyright (c) F. Caron (University of Oxford), A. Todeschini (Inria), and 
            % X. Miscouridou (University of Oxford)
            % caron@stats.ox.ac.uk
            % adrien.todeschini@gmail.com
            % xenia.miscouridou@spc.ox.ac.uk
            % August 2017
            %--------------------------------------------------------------------------

            % Creates a graph model
            obj.type = type;
            switch(type)   
                
                case 'ER' 
                    %% -------------   ERDOS-RENYI GRAPH --------------  %%
                    if numel(varargin)~=2
                        error('Erdos-Renyi graph must have two arguments n and p');
                    end
                    obj.param.n = varargin{1};
                    obj.param.p = varargin{2};
                    obj.name = 'Erdos-Renyi';
                    ERcheckparams( obj.param.n, obj.param.p );                     
                
                case 'GGP' 
                    %% -------------   GGP GRAPH --------------          %%
                    if numel(varargin)>4
                        error('GGP graph must have at most 4 arguments');
                    end
                   
                    obj.name = 'Generalized gamma process';
                    if numel(varargin)< 1
                         obj.param = struct('alpha', [.01,.01],...
                        'sigma',[.01,.01],...
                        'tau',[.01,.01]);
                    end
                    
                    if nargin>3
                    obj.param = struct('alpha', varargin{1},...
                        'sigma', varargin{2},...
                        'tau', varargin{3});
                    end
                    if numel(varargin)==4
                        obj.typegraph = validatestring(varargin{4}, {'undirected', 'simple', 'bipartite'});
                    end

                    if numel(varargin)<4 || strcmp(obj.typegraph, 'undirected') || strcmp(obj.typegraph, 'simple')
                        [alpha, sigma, tau] = GGPgetparams(obj.param.alpha, obj.param.sigma, obj.param.tau, 'default');
                        GGPcheckparams(alpha, sigma, tau);
                    elseif strcmp(obj.typegraph, 'bipartite')
                        if numel(obj.param)==1
                            obj.param(2) = obj.param(1);
                        end
                        [alpha1, sigma1, tau1] = GGPgetparams(obj.param(1).alpha, obj.param(1).sigma, obj.param(1).tau, 'default');
                        GGPcheckparams(alpha1, sigma1, tau1);
                        [alpha2, sigma2, tau2] = GGPgetparams(obj.param(2).alpha, obj.param(2).sigma, obj.param(2).tau, 'default');
                        GGPcheckparams(alpha2, sigma2, tau2);
                    end
                    
                
                case 'CGGP'
                    %% -------------   Compound GGP GRAPH -------------  %%
                    if numel(varargin)>9
                        error('CGGP graph must have at most 9 arguments');
                    end

                    obj.name = 'Compound generalized gamma process';

                    observe_all = false;
                    if numel(varargin)>=8
                        observe_all = varargin{8};
                    end
                    necessarily_sparse = false;
                    if numel(varargin)>=9
                        necessarily_sparse = varargin{9};
                    end
                    
                    if numel(varargin)>0 %default values which will be overwritten if given by user
                        alpha= [.01,.01]; 
                        sigma=[.01,.01];
                        tau = [1,1]; 
                        gamma= zeros(varargin{1},1);
                        Fdist = struct;
                        Fdist.name ='gamma';
                        Fdist.param.a = 0.01*ones(varargin{1}, 2);
                        Fdist.param.b = 0.01*ones(varargin{1}, 2);
                    else
                        error('Not enough input arguments')
                    end
                         
                    
                    if numel(varargin)>1
                        alpha = varargin{2};
                    end
                    if numel(varargin)>2
                        sigma = varargin{3};
                    end
                    if numel(varargin)>3
                        tau = varargin{4};
                    end
                    if numel(varargin)>4
                        Fdist = varargin{5};                        
                    end
                    if numel(varargin)>5
                        gamma = varargin{6};
                    end
                       
                    obj.param = struct('p', varargin{1},...
                        'alpha', alpha,...
                        'sigma', sigma,...
                        'tau', tau,...
                        'gamma', gamma,...
                        'observe_all', observe_all,...
                        'necessarily_sparse', necessarily_sparse); 
                    
                    if numel(varargin)>=7
                        obj.typegraph = validatestring(varargin{7}, {'undirected', 'simple', 'bipartite'});
                    end

                    if numel(varargin)<7 || strcmp(obj.typegraph, 'undirected') || strcmp(obj.typegraph, 'simple')
                        obj.param.Fdist = Fdist;
                        [alpha, sigma, tau, Fdist, gamma] = CGGPgetparams(obj.param.p, obj.param.alpha,...
                            obj.param.sigma, obj.param.tau, obj.param.Fdist, obj.param.gamma, 'default', ...
                            obj.param.observe_all, obj.param.necessarily_sparse);
                        CGGPcheckparams(obj.param.p, alpha, sigma, tau, Fdist, gamma);
                    elseif strcmp(obj.typegraph, 'bipartite')
                        if numel(obj.param)==1
                            obj.param(2) = obj.param(1);
                        end
                        numel(Fdist)
                        if numel(Fdist)==1
                            Fdist(2) = Fdist(1);
                        end
                        for i=1:2
                            obj.param(i).Fdist = Fdist(i);
                        end
                        if obj.param(1).p ~= obj.param(2).p
                            error('The parameter p must be the same for the two types of nodes');
                        end
                        [alpha1, sigma1, tau1, Fdist1, gamma1] = CGGPgetparams(obj.param(1).p, obj.param(1).alpha,...
                            obj.param(1).sigma, obj.param(1).tau, obj.param(1).Fdist, obj.param(1).gamma, 'default', ...
                            obj.param(1).observe_all, obj.param(1).necessarily_sparse);
                        CGGPcheckparams(obj.param(1).p, alpha1, sigma1, tau1, Fdist1, gamma1);
                        [alpha2, sigma2, tau2, Fdist2, gamma2] = CGGPgetparams(obj.param(2).p, obj.param(2).alpha,...
                            obj.param(2).sigma, obj.param(2).tau, obj.param(2).Fdist, obj.param(2).gamma, 'default', ...
                            obj.param(2).observe_all, obj.param(2).necessarily_sparse);
                        CGGPcheckparams(obj.param(2).p, alpha2, sigma2, tau2, Fdist2, gamma2);
                    end
                    
                case 'BA'
                    %% -----------   BARABASI-ALBERT GRAPH ------------  %%
                    if numel(varargin)~=1
                        error('BA graph must have one parameter n');
                    end
                    obj.param.n = varargin{1};
                    obj.name = 'Barabasi-Albert';
                    BAcheckparams( obj.param.n );
                    
                case 'Lloyd'
                    %% -----------   Lloyd et al. GRAPH ------------  %%
                    if numel(varargin)~=4
                        error('Lloyd graph must have 4 parameters n, sig, c and d');
                    end
                    obj.name = 'Lloyd';
                    obj.param.n = varargin{1};
                    obj.param.sig = varargin{2};
                    obj.param.c = varargin{3};
                    obj.param.d = varargin{4};
                    Lloydcheckparams(obj.param.n, obj.param.sig, obj.param.c, obj.param.d);
                    
                case 'MMSB'
                    %% ----------   Mixed Membership SBM GRAPH --------  %%
                    if numel(varargin)>6
                        error('MMSB graph must have at most 5 parameters n, p, alpha, W, and rho, epsilon');
                    end
                    obj.name = 'Mixed membership stochastic blockmodel';               
                    if numel(varargin) == 5
                        epsilon = [];
                    elseif numel(varargin) == 6
                        epsilon = varargin{6};
                    end    
                    obj.param = struct('n', varargin{1},...
                            'p', varargin{2},...
                            'alpha', varargin{3},...
                            'W', varargin{4},...
                            'rho', varargin{5},...
                            'epsilon', epsilon);                                    
                otherwise
                    %% -----------   UNRECOGNISED GRAPH ------------  %%
                    error('Unknown type %s', type);
            end
        end

        %% Sample a graph
        [G, varargout] = graphrnd(obj, varargin)
    end
end
