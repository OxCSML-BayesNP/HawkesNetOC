classdef graphmcmc
    % graphmcmc class of MCMC parameters and samples for inference on graph
    %
    % PROPERTIES
    %   - <a href="matlab: help graphmcmc/prior">prior</a>: prior on random graph
    %   - <a href="matlab: help graphmcmc/settings">settings</a>: settings of the MCMC algorithm
    %   - <a href="matlab: help graphmcmc/samples">samples</a>: MCMC samples
    %   - <a href="matlab: help graphmcmc/stats">stats</a>: Statistics of the MCMC algorithm
    %
    % CLASS CONSTRUCTOR
    %   - <a href="matlab: help graphmcmc/graphmcmc">graphmcmc</a>
    %
    % METHODS
    %   - <a href="matlab: help graphmcmc/graphmcmcsamples">graphmcmcsamples</a>: runs a MCMC algorithm for posterior inference on graphs
    %   - <a href="matlab: help graphmcmc/graphest">graphest</a>: returns median estimates of the graph parameters
    %   - <a href="matlab: help graphmcmc/graphinit">graphinit</a>: normalizes the hyperparameters for identifiability
    %   - <a href="matlab: help graphmcmc/graphnormalize">graphnormalize</a>: runs a MCMC algorithm for posterior inference on graphs

    properties
        % Object of class graphmodel
        prior = graphmodel('GGP', [0.01, 0.01], [0.01, 0.01], [0.01, 0.01]);

        % Settings of the MCMC algorithm
        settings = struct('niter', 1000, 'nburn', 500, 'thin', 10, 'nchains', 3, 'store_w', true,...
            'hyper', struct('rw_std', [.02, .02], 'MH_nb', 5));

        % MCMC Samples from the posterior
        samples;

        % Statistics on the MCMC algorithm
        stats;
    end

    %% METHODS
    methods

        % Class constructor
        function obj = graphmcmc(objmodel, niter, nburn, thin, nchains, nadapt, store_w)

        % GRAPHMCMC creates an object of class graphmcmc
        % obj = GRAPHMCMC(objmodel, niter, nburn, thin, nadapt, nchains, store_w)
        %
        % -------------------------------------------------------------------------
        % INPUTS
        %   - objmodel: an object of class graphmodel
        % Optional inputs
        %   - niter: number of MCMC iterations (default: 1000)
        %   - nburn: number of burn-in iterations (default:niter/2)
        %   - thin: thinning of the MCMC output (default: 1)
        %   - nadapt: number of iterations for adaptations (default:nburn/2)
        %   - nchains: number of MCMC chains (default: 1)
        %   - store_w: logical. true if we want to store and return samples for w
        %                   (default: true)
        %
        % OUTPUT
        %   - objmcmc: an object of the graphmcmc class
        %
        % See also GRAPHMCMC, GRAPHMODEL, GRAPHMODEL.graphmodel, GRAPHMCMCSAMPLES, GRAPHEST
        %
        % -------------------------------------------------------------------------
        % EXAMPLE
        % objmodel = graphmodel('GGP', 100, 0, 1);
        % niter = 10000; nburn = 1000; thin = 10; nadapt = 500; nchains = 3
        % objmcmc = graphmcmc(objmodel, niter, nburn, thin, nadapt, nchains)

        % Copyright (c) F. Caron (University of Oxford), A. Todeschini (Inria), and 
        % X. Miscouridou (University of Oxford)
        % caron@stats.ox.ac.uk
        % adrien.todeschini@gmail.com
        % xenia.miscouridou@spc.ox.ac.uk
        % August 2017

            if ~isa(objmodel, 'graphmodel')
                error('First argument must be a model of class graphmodel');
            end
            obj.prior = objmodel;
            switch(objmodel.type)
                case 'GGP'
                    switch(objmodel.typegraph)
                        case {'undirected','simple'}
                            % Settings of the MCMC algorithm
                            obj.settings.leapfrog = struct('L', 5, 'epsilon', .1, 'nadapt', 250);
%                             obj.settings.hyper = struct('rw_std', [.02, .02], 'MH_nb', 2);
                            obj.settings.latent = struct('MH_nb', 1);
                            % Samples
                            obj.samples = struct('w', [], 'w_rem', [], ...
                                'alpha', [], 'logalpha', [], 'sigma', [], 'tau', []);
                            % stats
                            obj.stats = struct('rate', [], 'rate2', []);
                        case 'bipartite'
                            % Samples
                            obj.samples = struct('w1', [], 'w2', [],...
                                'w1_rem', [], 'w2_rem', [],...
                                'alpha1', [], 'sigma1', [], 'tau1', [],...
                                'alpha2', [], 'sigma2', [], 'tau2', [],...
                                'logalpha1', [], 'logalpha2', []);
                            % stats
                            obj.stats = struct();
                        otherwise
                            error('Inference not supported for graph of type %s %s', objmodel.type, objmodel.typegraph);
                    end

                case 'CGGP'
                    switch(objmodel.typegraph)
                        case {'undirected', 'simple'}
                            % Settings of the MCMC algorithm
                            obj.settings.leapfrog = struct('L', 10, 'epsilon', .2,...
                                'nadapt', 250, 'adaptrate', 0.01, 'adaptwidth', 1000);
                            rw_std = struct('alpha', 0.02, 'sigma', .02, 'tau', .02, 'gamma', .02, 'a', .02, 'b', .02);
                            obj.settings.hyper = struct('rw_std', rw_std, 'MH_nb', 1, ...
                                'nadapt', 250, 'adaptrate', 0.005, 'adaptwidth', 200);
                            obj.settings.latent = struct('MH_nb', 1, 'nlatent', 20);
                            obj.settings.totalmass = struct('ntotalmass', 10);
                            % Samples
                            obj.samples = struct('w', [], 'w_rem', [], ...
                                'alpha', [], 'logalpha', [], 'sigma', [], 'tau', [], ...
                                'gamma', [], 'Fparam', []);
                            % stats
                            obj.stats = struct('rate', [], 'rate2', []);
                        case 'bipartite'
                            % Settings of the MCMC algorithm
                            rw_std = struct('alpha', 0.02, 'sigma', .02, 'tau', .02, 'gamma', .02, 'a', .02, 'b', .02);
                            rw_std(2) = rw_std;
                            obj.settings.hyper = struct('rw_std', rw_std, 'MH_nb', 2, ...
                                'nadapt', 250, 'adaptrate', 0.01, 'adaptwidth', 1000);
                            obj.settings.latent = struct('nlatent', 20);
                            obj.settings.ntotalmass = 20;
                            % Samples
                            obj.samples = struct('w', [], 'w_rem', [], ...
                                'alpha', [], 'logalpha', [], 'sigma', [], 'tau', [], 'gamma', [], 'Fparam', []);
                            obj.samples(2,1) = obj.samples;
                            % stats
                            obj.stats = struct('rate', [], 'rate2', []);
                        otherwise
                            error('Inference not supported for graph of type %s %s', objmodel.type, objmodel.typegraph);
                    end
                    case 'MMSB'
                        switch(objmodel.typegraph)
                            case {'undirected', 'simple'}
                                % Settings of the MCMC algorithm
                                obj.settings = struct;
                                obj.settings.hyper = struct( 'nmh', 10,'rw_std',0.02 );
                                obj.settings.store_s = 1;
                                % Samples
                                obj.samples = struct('alpha', [], 'W', [], ...
                                    'rho', [], 'pis',[], 's', [] );
                                % stats
                                obj.stats = struct('rateW', [], 'rate_alpha', [],'rate_rho', []);
                            otherwise
                                error('Inference not supported for graph of type %s %s', objmodel.type, objmodel.typegraph);
                        end
                otherwise
                    error('Unknown type of graph %s', objmodel.type);
            end

            if nargin>1
                obj.settings.niter = niter;
            end
            if nargin>2
                obj.settings.nburn = nburn;
            elseif nargin>1
                obj.settings.nburn = floor(niter/2);
            end
            if nargin>3
                obj.settings.thin = thin;
            else
                obj.settings.thin = ceil((obj.settings.niter-obj.settings.nburn)/500);
            end
            if nargin > 4
                obj.settings.nchains = nchains;
            end

            if ~strcmp(objmodel.type, 'MMSB')  
                if nargin>5
                    obj.settings.leapfrog.nadapt = nadapt;
                    obj.settings.hyper.nadapt = nadapt;
                else
                    nburn = obj.settings.nburn;
                    obj.settings.leapfrog.nadapt = floor(nburn/2);  
                    obj.settings.hyper.nadapt = floor(nburn/2);
                end
                if nargin>6
                    obj.settings.store_w = store_w;
                end
                obj.settings.leapfrog.adaptwidth = floor(obj.settings.leapfrog.nadapt/2);
            end
        end

        % Runs a MCMC sampler with 1 community to use as initial values
        init = graphinit(objmcmc, G, niter1, verbose)
        
        % Runs a MCMC sampler
        objmcmc = graphmcmcsamples(objmcmc, G, varargin);

        % Normalizes the hyperparameters for identifiability
        objmcmc = graphnormalize(objmcmc)

        % Returns estimates from the MCMC output
        [estimates, varargout] = graphest(objmcmc, nburn, verbose, cost);
    
    end
end
