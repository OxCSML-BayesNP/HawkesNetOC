function h_hist = plot_hist(samples, variables, namesvar, trueval, ind_feat, col_feat, rep, prefix, suffix)

% plot_hist plots the trace and histograms of the variables for the MCMC
% output

% INPUT
%   - samples: MCMC output samples
%   - variables: cell of strings with the variable names
%   - namesvar: cell of strings with variables tex names
%   - trueval: cell of strings with the true variable values (empty if not known)
%   - ind_feat: vector of size p (p = number of features) with the
%      ordered features 1,..,p
%   - col_feat: matrix of size px3 with the color properties for the features 
%   - rep: output directory to save the figure
%   - prefix: output filename prefix
%   - suffix: output filename suffix
%
% OUTPUT
%   - h_hist: vector of handles to the plotted lines



if nargin<2
    variables = {'logalpha', 'sigma', 'tau', 'Fparam.a', 'Fparam.b', 'w_rem'};
    namesvar = {'$\log \alpha$', '$\sigma$', '$\tau$', '$a$', '$b$', '$w_*$'};
    namesvar(2,:) = cellfun(@(x) [x '^\prime'], namesvar, 'UniformOutput', false);
elseif nargin<3
    namesvar = variables;
    namesvar(2,:) = cellfun(@(x) [x '^\prime'], variables, 'UniformOutput', false);
end
if nargin<4
    trueval = {};
end
p = size(samples(1).w, 2);
if nargin<5 || isempty(ind_feat)
    ind_feat = 1:p;
end
if nargin<6 || isempty(col_feat)
    col_feat = hsv(p);
end
if nargin<7
    rep = [];
end
if nargin<8
    prefix = '';
end
if nargin<9
    suffix = '';
end


ntypes = size(samples, 1);

%% Trace posterior histograms for parameters
h_hist = gobjects(ntypes, ntypes);
for t=1:ntypes% loop over types of nodes
    for i=1:numel(variables)% loop over variables
        var = strsplit(variables{i}, '.');
        h_hist(t, i) = figure;
        hold on
        if length(var)==1
            data = combine_chains(samples(t,:), var{1});
        elseif length(var)==2
            data = combine_chains2(samples(t,:), var{1}, var{2});
        end
        sz = size(data);
        if prod(sz(1:2))==1
            hist(squeeze(data), 20)
        else
            data = squeeze(data)';
            if size(data,2)==p
                hist(data(:,ind_feat), 10)
            else
                hist(data, 10)
            end
        end
        
        if ~isempty(trueval)
            h_true = plot(trueval{t,i}(:), 0, 'g*', 'linewidth', 3, 'markersize', 12);
        end
        
        if size(data,2)==p
            h = findobj(gca, 'Type', 'patch');
            for k=1:numel(h)
                set(h(k), 'FaceColor', col_feat(k,:), 'EdgeColor', 'w')
            end
            leg = cellfun(@(x) ['Feat. ' num2str(x)], num2cell(1:p), 'uniformoutput', false);
            if ~isempty(trueval)
                leg = [leg, 'True'];
            end
            legend(leg, 'location', 'northwest')
            legend boxoff
        else
            h = findobj(gca, 'Type', 'patch');
            for k=1:numel(h)
                set(h(k), 'EdgeColor', 'w')
            end
            if ~isempty(trueval)
                legend(h_true(1), 'True', 'location', 'northwest')
                legend boxoff
            end
        end
        
        xlabel(namesvar{t,i}, 'fontsize', 16, 'interpreter', 'latex');
        ylabel('Nb MCMC samples', 'fontsize', 16);
            
        box off
        if ~isempty(rep)
            savefigs(gcf, [prefix 'hist_' var{end} num2str(t) suffix], rep);
        end
    end
end

end

function samples_all = combine_chains(samples, name)
nchains = size(samples, 2);
nsamples = size(samples(1).(name), ndims(samples(1).(name)));
for k = nchains:-1:1
    ind_all = (k-1)*nsamples+1:k*nsamples;
    if ismatrix(samples(1,k).(name))
        samples_all(:, ind_all) = samples(1,k).(name);
    else
        samples_all(:, :, ind_all) = samples(1,k).(name);
    end
end
end

function samples_all = combine_chains2(samples, name1, name2)
nchains = size(samples, 2);
nsamples = size(samples(1).(name1).(name2), ndims(samples(1).(name1).(name2)));
for k = nchains:-1:1
    ind_all = (k-1)*nsamples+1:k*nsamples;
    if ismatrix(samples(1,k).(name1).(name2))
        samples_all(:, ind_all) = samples(1,k).(name1).(name2);
    else
        samples_all(:, :, ind_all) = samples(1,k).(name1).(name2);
    end
end
end
