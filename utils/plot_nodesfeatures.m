function h = plot_nodesfeatures(w, ind, ind_features, names, featnames, color, rep, prefix, suffix)

% plot_nodesfeatures gives a boxplot with the normalized vector of
% the weights for some given the nodes
%
%
% INPUT
%   - w: matrix of weight parameters
%   - ind: indices of nodes in a given order
%   - ind_features: vector of size p (p = number of features) with the
%      ordered features 1,..,p
%   - names: names of groups
%   - featnames: names of features 
%   - color: matrix of size p x 3 of the rgb values for the group colors
%   - rep:  output directory
%   - prefix: output filename prefix
%   - suffix: output filename suffix

%
% OUTPUT
% h: vector of handles to the plotted lines

p = size(w, 2);

if nargin<5 || isempty(featnames)
    featnames = arrayfun(@(x) ['Feature ' num2str(x)], 1:p, 'uniformoutput', false);
end
if nargin<6 || isempty(color)
    color = hsv(p);
end
if nargin<7
    rep=[];
end
if nargin<8
    prefix = '';
end
if nargin<9
    suffix = '';
end

% normalize
w = bsxfun(@rdivide, w(ind, ind_features), sum(w(ind,:), 2));

h = figure;
bbar = barh(w, 'stacked');
legend(featnames, 'location', 'BestOutside')
legend boxoff
for k=1:p
    set(bbar(k), 'FaceColor', color(k,:), 'EdgeColor', color(k,:));
end
if nargin>=4
    set(gca,'yticklabel', names, 'ytick', 1:length(ind))
end
set(gca,'xtick',[])
set(gca,'XColor','w')
xlim([0,1])
ylim([0.5,length(ind)+0.5])
box off

if ~isempty(rep)
    savefigs(gcf, [prefix 'namespct' suffix], rep)
end
