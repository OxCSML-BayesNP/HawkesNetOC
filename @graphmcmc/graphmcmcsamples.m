function objmcmc = graphmcmcsamples(objmcmc, G, verbose, init, varargin)

% graphmcmcsamples runs a MCMC algorithm for posterior inference on graphs
%
% objmcmc = graphmcmcsamples(objmcmc, G, verbose, init)
% -------------------------------------------------------------------------
% INPUTS
%   - objmcmc: an object of the class graphmcmc, containing the graph model
%           specifications and the parameters of the MCMC algorithn
%   - G: sparse binary adjacency matrix
% Optional inputs:
%   - verbose: logical. If true (default), print information
%   - init: structure of initial values for the markov chain (CGGP only)
%
% OUTPUT
%   - objmcmc: Updated graphmcmc object with the set of samples
% -------------------------------------------------------------------------
%
% See also GRAPHMCMC, GRAPHMCMC.GRAPHMCMC, GRAPHEST

% Copyright (c) F. Caron (University of Oxford), A. Todeschini (Inria), and 
% X. Miscouridou (University of Oxford)
% caron@stats.ox.ac.uk
% adrien.todeschini@gmail.com
% xenia.miscouridou@spc.ox.ac.uk
% October 2017

isparallel = false; 
if ~isempty(ver('distcomp')) % Check if parallel toolbox installed
    if getpoolsize()>0 % If parallel pool initiated
        isparallel = true; % Will run MCMC chains in parallel using parallel computing toolbox
    end
end

if nargin<3
    verbose = true;
end
if nargin<4
    init = [];
end

objmodel = objmcmc.prior;
if ~strcmp(objmodel.type,'MMSB')   
    if nargin<4 || isempty(init) || (isstruct(init) && isempty(fieldnames(init)))
        init = struct;
        if strcmp(objmodel.typegraph, 'bipartite')
            init(2,1:objmcmc.settings.nchains) = struct;
        else
            init(1,1:objmcmc.settings.nchains) = struct;
        end
    elseif size(init, 2) ~= objmcmc.settings.nchains
        init = repmat(init, [1,objmcmc.settings.nchains]);
    end
end

% random seed for parallel chains
seed = randi(floor(intmax/10));

if verbose
    print_info(objmcmc, G, init);
end

% random seed for parallel chains
seed = randi(floor(intmax/10));

tic
switch (objmodel.type)
    case 'GGP'
        switch(objmodel.typegraph)
            case {'undirected', 'simple'}  
                % Run MCMC algorithms
                if isparallel                   
                    parfor (k=1:objmcmc.settings.nchains, getpoolsize()) 
                        if verbose
                            print_info_chain(k, objmcmc.settings.nchains);                            
                        end
                        rng(seed+k)
                        [samples(:,k), stats(:,k)] = GGPgraphmcmc(G, ...
                            objmodel.param, objmcmc.settings, objmodel.typegraph, verbose);
                    end
                else % not parallel
                    for k=1:objmcmc.settings.nchains 
                        if verbose
                            print_info_chain(k, objmcmc.settings.nchains);                            
                        end
                        rng(seed+k);
                        [samples(:,k), stats(:,k)] = GGPgraphmcmc(G, ...
                            objmodel.param, objmcmc.settings, objmodel.typegraph, verbose);
                    end
                end
                objmcmc.samples = samples;
                objmcmc.stats = stats;
            case 'bipartite'
                % Run MCMC algorithms
                if isparallel                    
                    parfor (k=1:objmcmc.settings.nchains, getpoolsize()) 
                        rng(seed+k);
                        if verbose
                            print_info_chain(k, objmcmc.settings.nchains);                            
                        end
                        [samples(:,k), stats(:,k)] = GGPbipgraphmcmc(G, ...
                            objmodel.param, objmcmc.settings, verbose);
                    end
                else % not parallel                    
                    for k=1:objmcmc.settings.nchains 
                        rng(seed+k);
                        if verbose
                            print_info_chain(k, objmcmc.settings.nchains);                            
                        end
                        [samples(:,k), stats(:,k)] = GGPbipgraphmcmc(G, ...
                            objmodel.param, objmcmc.settings, verbose);
                    end
                end                
                objmcmc.samples = samples;
                objmcmc.stats = stats;
            otherwise
                error('Unknown type of graph %s', objmodel.typegraph);
        end

    case 'CGGP'
        switch(objmodel.typegraph)
            case {'undirected', 'simple'}
                % Run MCMC algorithms
                if ~isparallel
                    [samples,stats] = for_CGGPgraphmcmc(seed, G, objmodel, objmcmc, objmodel.typegraph, verbose, init, varargin{:});
                else
                    [samples,stats] = parfor_CGGPgraphmcmc(seed, G, objmodel, objmcmc, objmodel.typegraph, verbose, init, varargin{:});
                end
                for ch=1:objmcmc.settings.nchains 
                    % Compute some additional summary
                    samples(ch).mean_w_rem = mean(samples(ch).w_rem, 2);
                    samples(ch).mean_w = mean(samples(ch).w, 2);
                end
                objmcmc.samples = samples;
                objmcmc.stats = stats;                
             case 'bipartite'
                 error('Sampler for CGGP for bipartite graphs not implemented')
            otherwise
                error('Unknown type of graph %s', objmodel.typegraph);
        end
        
        
    case 'MMSB'
        switch(objmodel.typegraph)
            case {'undirected', 'simple'}
                % Run MCMC algorithms
                if ~isparallel   
                    [samples, stats] = for_mmsbmcmc(seed, G, objmodel, objmcmc, verbose);
                else  
                    [samples,stats] = parfor_mmsbmcmc(seed, G, objmodel, objmcmc, verbose);
                end
                objmcmc.samples = samples;
                objmcmc.stats = stats;
            otherwise
                error('Unknown type of graph %s', objmodel.typegraph);
        end
    otherwise
        error('Inference not implemented for graph model of type %s', objmodel.type);
end
if verbose
    print_info_end;
end

end


function [samples,stats] = parfor_CGGPgraphmcmc(seed, G, objmodel, objmcmc, typegraph, verbose, init, varargin)
parfor (k=1:objmcmc.settings.nchains, getpoolsize())
    rng(seed+k)
    if verbose
        print_info_chain(k, objmcmc.settings.nchains);                            
    end
    [samples(:,k), stats(:,k)] = CGGPgraphmcmc(G, ...
        objmodel.param, objmcmc.settings, typegraph, verbose, 'init', init(k), varargin{:});
end
end

function [samples,stats] = for_CGGPgraphmcmc(seed, G, objmodel, objmcmc, typegraph, verbose, init, varargin)
for k=1:objmcmc.settings.nchains
    rng(seed+k);
    if verbose
        print_info_chain(k, objmcmc.settings.nchains);                            
    end
    [samples(:,k), stats(:,k)] = CGGPgraphmcmc(G, ...
        objmodel.param, objmcmc.settings, typegraph, verbose, 'init', init(k), varargin{:});
end
end


function [samples,stats] = for_mmsbmcmc(seed, G, objmodel, objmcmc,verbose)
for k=1:objmcmc.settings.nchains
    rng(seed+k);
    if verbose
        print_info_chain(k, objmcmc.settings.nchains);                            
    end
    [samples(:,k), stats(:,k)] = mmsbmcmc(G, objmodel, objmcmc, verbose);
end
end

function [samples,stats] = parfor_mmsbmcmc(seed, G, objmodel, objmcmc, verbose)
parfor (k=1:objmcmc.settings.nchains, getpoolsize())
    rng(seed+k);
    if verbose
        print_info_chain(k, objmcmc.settings.nchains);                            
    end
    [samples(:,k), stats(:,k)] = mmsbmcmc(G, objmodel, objmcmc, verbose);
end
end


function print_info_chain(k, nchains)
fprintf('-----------------------------------\n')
fprintf(' Markov chain %d/%d \n', k, nchains);
fprintf('-----------------------------------\n')
end

function print_info(objmcmc, G, init)

niter = objmcmc.settings.niter;
nchains = objmcmc.settings.nchains;
% Get an estimate of the computation time by running a short chain
niter0 = min(400, ceil(niter/20));
objmcmctemp = objmcmc;
objmcmctemp.settings.nchains = 1;
objmcmctemp.settings.niter = niter0;
tic
graphmcmcsamples(objmcmctemp, G, false, init); % Run short MCMC

time = toc * niter/niter0 * nchains; % extrapolate running time
hours = floor(time/3600);
minutes = (time - hours*3600)/60;

% Compute number of nodes, edges and missing data
nnodes = size(G, 1);
if strcmp(objmcmc.prior.typegraph, 'simple') % If no self-loops
    [ind1, ~] = find(triu(G, 1));
else
    [ind1, ~] = find(triu(G));
end
nedges = numel(ind1);
G_ismiss = sparse(isnan(G));
nmiss = nnz(G_ismiss);

% Print some info
fprintf('-----------------------------------\n')
fprintf('Start MCMC for %s graphs\n', objmcmc.prior.type)
fprintf('Nb of nodes: %d - Nb of edges: %d (%d missing)\n', nnodes, nedges-nmiss, nmiss);
fprintf('Nb of chains: %d - Nb of iterations: %d\n', objmcmc.settings.nchains, objmcmc.settings.niter)
fprintf('Nb of parallel workers: %d\n', max(getpoolsize(), 1))
fprintf('Estimated computation time: %.0f hour(s) %.0f minute(s)\n', hours, minutes);
fprintf('Estimated end of computation: %s \n', datestr(now + time/3600/24));
end

function print_info_end()
time = toc;
hours = floor(time/3600);
minutes = (time - hours*3600)/60;
fprintf('-----------------------------------\n')
fprintf('End MCMC\n')
fprintf('Computation time: %.0f hour(s) %.0f minute(s)\n', hours, minutes);
fprintf('-----------------------------------\n')
end