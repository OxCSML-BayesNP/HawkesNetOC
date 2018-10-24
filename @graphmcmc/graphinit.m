function init = graphinit(objmcmc, G, niter1, verbose)

% graphinit runs a MCMC algorithm with p=1 community for to be used as an
% initialization for a CGGP model with p>1 communities
%
% init = graphinit(objmcmc, G, p);
% -------------------------------------------------------------------------
% INPUTS
%   - objmcmc: an object of the class graphmcmc, containing the graph model
%     specifications (has to be a CGGP graph) and the parameters of the MCMC algorithm for the CGGP
%   - G: sparse binary adjacency matrix
%   - niter1: number of MCMC iterations
% Optional inputs:
%   - verbose: logical
%
% OUTPUT
%   - init: structure of initial values for the markov chain
%
% -------------------------------------------------------------------------
%
% See also graphmcmcsamples, graphmcmc, graphmcmc/graphmcmc

% Copyright (c) F. Caron (University of Oxford), A. Todeschini (Inria), and 
% X. Miscouridou (University of Oxford)
% caron@stats.ox.ac.uk
% adrien.todeschini@gmail.com
% xenia.miscouridou@spc.ox.ac.uk
% August 2017

if nargin<3
    niter1 = 10000;
end
if nargin < 4
    verbose = true;
end
if ~strcmp(objmcmc.prior.type, 'CGGP')
    error('function graphinit not implemented for model of type %s', objmcmc.prior.type)
end
p = objmcmc.prior.param.p;
nchains = objmcmc.settings.nchains;

objmodelcggp = objmcmc.prior; % CGGP model 
if strcmp(objmodelcggp.typegraph, 'bipartite')
    error('Initialisation for the bipartite graph not implemented yet')
end

% Create a GGP model and MCMC structure 
objmodelggp = graphmodel('GGP', objmodelcggp.param.alpha, objmodelcggp.param.sigma, [.01, .01], objmodelcggp.typegraph);
objmcmcggp = graphmcmc(objmodelggp, niter1, niter1-1, 1, nchains);% objmcmc;

if verbose
    fprintf('-----------------------------------\n');
    fprintf('Start initialisation of the MCMC algorithm for CGGP\n');
end
N = size(G, 1);
objmcmcggp = graphmcmcsamples(objmcmcggp, G, false);
init = getsample(objmcmcggp.samples, size(objmcmcggp.samples(1).w, 2));

for ch=1:objmcmc.settings.nchains
    init(ch).beta = 1/sqrt(p)/init(ch).tau*gamrnd(1, 1, N, p);
    init(ch).w0 = init(ch).w*init(ch).tau; 
    init(ch).logalpha = init(ch).logalpha + init(ch).sigma*log(init(ch).tau); 
    init(ch).alpha = exp( init(ch).logalpha); 
    init(ch).tau = 1; 
    b = ones(p, 1 );
    a = ones(p, 1 );
    init(ch).Fdist.name = 'gamma';
    init(ch).Fdist.param.a = a;
    init(ch).Fdist.param.b = b;
    if numel(objmcmc.prior.param.gamma)==p
        init(ch).gamma = objmcmc.prior.param.gamma;
    else
        error('Case gamma unknown not implemented');
    end
   
    if objmcmc.prior.param.observe_all
        init(ch).w_rem = zeros(1, p);
    else
        init(ch).w_rem = repmat(1/sqrt(p).*init(ch).w_rem, [1,p]);
    end
end

init = rmfield(init, 'w');
if verbose
    fprintf('-----------------------------------\n');
    fprintf('End initialisation\n');
    fprintf('-----------------------------------\n');    
end

end