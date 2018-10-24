function [h_boxplot, h_names] = plot_groups(w, groups, row_names, ind_features, label_groups, featnames, color, rep, prefix, suffix)

% plot_groups plots the weights of the nodes in groups
% (according to the true underlying groups)
%
% INPUT
%   - w: matrix of weight parameters
%   - groups: vector with true affiliations for the nodes
%   - row_names: names of nodes
%   - ind_features: vector of size p (p = number of features) with the
%      ordered features 1,..,p
%   - label_groups: names of groups
%   - featnames: names of features 
%   - color: matrix of size p x 3 of the rgb values for the group colors
%   - rep:  output directory
%   - prefix: output filename prefix
%   - suffix: output filename suffix

%
% OUTPUT
% h_boxplot: figure handle vector for the box plot
% h_names: figure handle vector for the feature names


p = size(w, 2);
if nargin<4
    ind_features = 1:p;
end
if nargin<5 || isempty(label_groups)
    label_groups = arrayfun(@(x) ['Group ' num2str(x)], sort(unique(groups)), 'uniformoutput', false);
end
if nargin<6 || isempty(featnames)
    featnames = arrayfun(@(x) ['Feature ' num2str(x)], 1:p, 'uniformoutput', false);
end
if nargin<7 || isempty(color)
    color = autumn(numel(label_groups));
end
if nargin<8
    rep=[];
end
if nargin<9
    prefix = '';
end
if nargin<10
    suffix = '';
end

%% Boxplot
xloc = .25+.1*(1:numel(label_groups));%(.35:.1)[0.35, .45];
h_boxplot = figure;
for k=1:p
    ind_k = ind_features(k);
    subplot(2,ceil(p/2),k)
    title(featnames{k})
    hold on
    boxplot(w(:,ind_k), groups, 'labels', label_groups)
    h = findobj(gca,'Tag','Box');
    for j=1:length(h)
        patch(get(h(j),'XData'),get(h(j),'YData'),color(end-j+1, :));
    end
    box off
    ylabel(sprintf('$w_{%d}$', ind_features(k)), 'interpreter', 'latex', 'fontsize', 20)
end
savefigs(gcf, [prefix 'boxplots' suffix], rep)


%% Names
if ~isempty(row_names)
    h_names = figure;
    for k=1:p
        ind_k = ind_features(k);
        subplot(2,ceil(p/2),k)
        title(featnames{k})
        for j=1:length(label_groups)
            ind = (groups ==j);
            text(xloc(j)*ones(size(w(ind,k))), w(ind,ind_k), row_names(ind), 'color', color(j, :),'FontSize',10, 'interpreter', 'none');
            hold on
            text(xloc(j), -max(w(:,ind_k))/7, label_groups{j}, 'fontsize', 16)
        end
        xlim([min(xloc)-.01, max(xloc)+.08])
        ylim([0, max(w(:,ind_k))])
        set(gca,'xtick',[])
        ylabel(sprintf('$w_{%d}$', ind_features(k)), 'interpreter', 'latex', 'fontsize', 20)
    end
    savefigs(gcf, [prefix 'names' suffix], rep)
else
    h_names = [];
end
