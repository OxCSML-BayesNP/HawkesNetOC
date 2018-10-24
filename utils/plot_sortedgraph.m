function [h, densmat] = plot_sortedgraph(G, nodefeat1, nodefeat2, ind_features, ...
    xylab, rep, prefix, suffix, formats, plot_graph, plot_blocks)

% PLOT_SORTEDGRAPH plots the adjacency matrix of the graph where the nodes are sorted by their max feature
% -------------------------------------------------------------------------
% INPUT
%   - G: adjacency matrix of the observed graph 
%   - nodefeat1: max feature of each source node
%   - nodefeat2: max feature of each target node
%   - ind_features: vector of size p (p = number of features) with the
%    ordered features 1,..,p
%   - xylab: cell of strings with the labels for the x and y axes respectively
%   - rep:  output directory
%   - prefix: output filename prefix
%   - suffix: output filename suffix
%   - formats: string or cell of strings for the file formats
%   - plot_blocks: logical to indicate whether to plot blocks in the
%     adjacency matrix
%
% OUTPUT
%   - h: figure handle vector
%   - densmat: matrix with the densities of the blocks
% -------------------------------------------------------------------------

p = max(ind_features);
if nargin<5
    xylab = {'Nodes type 1', 'Nodes type 2'};
end
if nargin<6
    rep=[];
end
if nargin<7
    prefix = '';
end
if nargin<8
    suffix = '';
end
if nargin<9
    formats = {'png', 'epsc2', 'fig'};
end
if nargin<10
    plot_graph = true;
end
if nargin<11
    plot_blocks = true;
end

%% relabel node features
[~,nodefeat1] = ismember(nodefeat1, ind_features);
[~,nodefeat2] = ismember(nodefeat2, ind_features);

%% reorder nodes by feature
[nodefeat1,ind1] = sort(nodefeat1);
[nodefeat2,ind2] = sort(nodefeat2);
G = G(ind1, ind2);

%% get cutlines coord
icut = NaN(p,1);
jcut = NaN(p,1);
for k=1:p
    i = find(nodefeat1==k, 1, 'last');
    if ~isempty(i)
        icut(k) = i;
    end
    j = find(nodefeat2==k, 1, 'last');
    if ~isempty(j)
        jcut(k) = j;
    end
end

%% density of the blocks
for k=p:-1:1
    ok1 = nodefeat1==k;
    for l=p:-1:1
        temp = G(ok1, nodefeat2==l);
        %%% TODO take missing entries NaN into account
        densmat(k,l) = nnz(temp)/numel(temp);
    end
end


%% Plot the sorted graph

if plot_graph
    % plot graph matrix
    h = figure; hold on
    axis ij
    axis equal
    ylim([0, size(G,1)])
    xlim([0, size(G,2)])
    dens = nnz(G)/numel(G);
    if dens>.1
        Gplot = 2*G - 1;
        Gplot(isnan(Gplot)) = 0;
        imagesc(Gplot);
        linecolor = 'w';
    else
        spy(G);
        linecolor = 'k';
    end
    
    % plot lines of blocks
    for k=1:p-1
        if ~isnan(icut(k))
            plot([0,size(G,2)], [icut(k)+.5, icut(k)+.5], linecolor, 'linewidth', 1.5)
        end
        if ~isnan(jcut(k))
            plot([jcut(k)+.5, jcut(k)+.5], [0,size(G,1)], linecolor, 'linewidth', 1.5);
        end
    end
    xlabel(xylab{1})
    ylabel(xylab{2})
    
    if ~isempty(rep)
        savefigs(gcf, [prefix 'sortedgraph' suffix], rep, formats);
    end
    
end

%% plot blocks
if plot_blocks
    densmat2 = densmat-min(densmat(:));
    densmat2 = densmat2./max(densmat2(:));
    
    h = figure; hold on
    axis ij
    axis equal
    ylim([0, size(G,1)])
    xlim([0, size(G,2)])
    
    % plot rectangles colored by density
    i = 0;
    for k=1:p
        j = 0;
        if ~isnan(icut(k))
            for l=1:p
                if ~isnan(jcut(l))
                    %                 fill([j j jcut(l) jcut(l)], [i icut(k) icut(k) i], 1-densmat2(k,l)*[.5 .5 0], 'edgecolor', 'none');
                    fill([j j jcut(l) jcut(l)], [i icut(k) icut(k) i], 1-densmat2(k,l)*[.8 .8 0], 'edgecolor', 'none');
                    j = jcut(l);
                end
            end
            i = icut(k);
        end
    end
    
    
    % plot lines of blocks
    for k=1:p-1
        if ~isnan(icut(k))
            plot([0,size(G,2)], [icut(k), icut(k)], 'k', 'linewidth', 1.5);
        end
        if ~isnan(jcut(k))
            plot([jcut(k), jcut(k)], [0,size(G,1)], 'k', 'linewidth', 1.5);
        end
    end
    plot([0, size(G,2)], [0, 0], 'k', 'linewidth', .5);
    plot([size(G,2), size(G,2)], [0, size(G,1)], 'k', 'linewidth', .1);
    
    xlabel(xylab{1})
    ylabel(xylab{2})
    
    if ~isempty(rep)
        savefigs(gcf, [prefix 'blocks' suffix], rep);
    end
    
end
