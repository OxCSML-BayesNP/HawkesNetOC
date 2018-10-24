function [h_trace] = plot_trace(samples, settings, variables, namesvar, trueval, rep, prefix, suffix)
% Plots the trace of the variables for the different MCMC chains
%
% INPUT
%   - samples: MCMC output samples
%   - settings: struct with the mcmc settings
%   - variables: cell of strings with the variable names
%   - namesvar: cell of strings with variables tex names
%   - trueval: cell of strings with the true variable values (empty if not known)
%   - rep:  output directory
%   - prefix: output filename prefix
%   - suffix: output filename suffix
%
% OUTPUT
%   - h_trace: figure handle vector for the trace plot


if nargin<2 || isempty(settings)
    niter = size(samples(1).logalpha, ndims(samples(1).logalpha));
    nburn = 0;
    thin = 1;
else
    niter = settings.niter;
    nburn = settings.nburn;
    thin = settings.thin;
end
if nargin<3
    variables = {'logalpha', 'sigma', 'tau', 'Fparam.a', 'Fparam.b', 'w_rem'};
    namesvar = {'$\log \alpha$', '$\sigma$', '$\tau$', '$a$', '$b$', '$w_*$'};
    namesvar(2,:) = cellfun(@(x) [x '^\prime'], namesvar, 'UniformOutput', false);
elseif nargin<4
    namesvar = variables;
    namesvar(2,:) = cellfun(@(x) [x '^\prime'], variables, 'UniformOutput', false);
end
colour = get(gca,'ColorOrder');
nbcolour = size(colour, 1);

if nargin<5
    trueval = {};
end
if nargin<6
    rep = [];
end
if nargin<7
    prefix = '';
end
if nargin<8
    suffix = '';
end

[ntypes, nchains] = size(samples);
leg = cell(nchains, 1);
for k=1:nchains
    leg{k} = ['Chain ' num2str(k)];
end


%% Trace plots and posterior histograms for parameters
h_trace = gobjects(ntypes, ntypes) ;
for t=1:ntypes% loop over types of nodes
    for i=1:numel(variables)% loop over variables
        var = strsplit(variables{i},'.');
        h_trace(t, i) = figure;
        for k=1:nchains
            colourk = colour(mod(k-1, nbcolour)+1,:);
            if length(var)==1
                h(k,:) = plot(thin:thin:(niter-nburn), squeeze(samples(t, k).(var{1})), 'color', colourk);
            elseif length(var)==2
                h(k,:) = plot(thin:thin:(niter-nburn), squeeze(samples(t, k).(var{1}).(var{2})),'color', colourk);
            end
            hold on
        end
        if ~isempty(trueval)
            h_true = plot([thin; (niter-nburn)], [trueval{t,i}(:)'; trueval{t,i}(:)'], 'g--', 'linewidth', 3);
            legend([h(:,1); h_true(1)], [leg; 'True'], 'fontsize', 16, 'location', 'Best');
        else
            legend(h(:,1), leg, 'fontsize', 16, 'location', 'Best');
        end
        legend boxoff
        xlabel('MCMC iterations', 'fontsize', 16);
        ylabel(namesvar{t,i} , 'fontsize', 16, 'interpreter', 'latex');
        box off
        xlim([0, niter-nburn])
        clear h;
        if ~isempty(rep)
%            savefigs(gcf, [prefix 'trace_' var{end} num2str(t) suffix], rep);
        end
    end
end
