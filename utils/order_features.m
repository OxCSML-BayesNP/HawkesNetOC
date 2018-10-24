function [ind_features] = order_features(G, nodefeat1, nodefeat2)

% order_features returns a p-vector with the features 1,..,p 
% with order according to the maximum afilliations of the nodes
% ----------------------------------------------------
% INPUTS
%   - G: adjacency matrix of the observed graph 
%   - nodefeat1: max feature of each source node
%   - nodefeat2: max feature of each target node
%----------------------------------------------------
% OUTPUT
%   - ind_features: vector of size p (p = number of features) with the
%    ordered features 1,..,p


p = max([nodefeat1;nodefeat2]);

% block diagonalize
for k = p:-1:1
    for l = p:-1:1
        temp = G(nodefeat1==k, nodefeat2==l);
        densmat(k,l) = full(mean(temp(:)));
    end
end

[~,ind_features] = sort(sum(densmat,2), 'descend');

densmat = densmat(ind_features, ind_features);
densmat = densmat-min(densmat(:));
densmat = densmat./max(densmat(:));

ind = symrcm(densmat);
ind_features = ind_features(ind);
