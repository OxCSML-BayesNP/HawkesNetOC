function [confmat] = print_classif(filename, nodefeat, groups, ind_feat, label_groups)

% print_classif calculates the classification accuracy and prints them in a
% file
%
% INPUT
%   - filename: string with the name of the file to print in the
%   classification performance
%   - nodefeat: max feature of each node
%   - groups: vector with true affiliations for the nodes
%   - ind_features: vector of size p (p = number of features) with the
%   - label_groups: names of groups
%
% OUTPUT
%   - confmat: confusion matrix with the classification accuracy


N = numel(groups); % nb of nodes
ngroups = max(groups);
nfeat = max(nodefeat);

% confusion matrix
for i = ngroups:-1:1
    for j = nfeat:-1:1
        confmat(i, j) = sum(groups==i & nodefeat==ind_feat(j));
    end
end

fid = fopen(filename, 'w');
fids = [1, fid];

mfprintf(fids, 'Classification performance\n')
mfprintf(fids, '==========================\n')

%% print confusion matrix in counts
mfprintf(fids, 'Confusion matrix (counts)\n')
mfprintf(fids, '-------------------------\n')
mfprintf(fids, 'Group        : ')
mfprintf(fids, 'Feat%2d ', 1:nfeat)
mfprintf(fids, '| Total\n')
for i=1:ngroups
    mfprintf(fids, '%-13s: ', label_groups{i})
    for j=1:nfeat
        mfprintf(fids, '%6d ', confmat(i,j))
    end
    mfprintf(fids, '|%6d\n', sum(confmat(i,:), 2))
end
mfprintf(fids, 'Total        : ')
for j=1:nfeat
    mfprintf(fids, '%6d ', sum(confmat(:,j)))
end
mfprintf(fids, '|%6d\n', sum(confmat(:)))
mfprintf(fids, '-------------------------\n')

%% print confusion matrix in percent
confmat = confmat/N;
mfprintf(fids, 'Confusion matrix (%%)\n')
mfprintf(fids, '-------------------------\n')
mfprintf(fids, 'Group        : ')
mfprintf(fids, 'Feat%2d ', 1:nfeat)
mfprintf(fids, '| Total\n')
for i=1:ngroups
    mfprintf(fids, '%-13s: ',  label_groups{i})
    for j=1:nfeat
        mfprintf(fids, '%6.2f ', confmat(i,j)*100)
    end
    mfprintf(fids, '|%6.2f\n', sum(confmat(i,:), 2)*100)
end
mfprintf(fids, 'Total        : ')
for j=1:nfeat
    mfprintf(fids, '%6.2f ', sum(confmat(:,j))*100)
end
mfprintf(fids, '|%6.2f\n', sum(confmat(:))*100)
mfprintf(fids, '-------------------------\n')

%% print classif perf metrics
% assign a group to each feature
[~, groupfeat] = max(confmat);

mfprintf(fids, 'Group assignments of features\n')
mfprintf(fids, '-------------------------\n')
for j=1:nfeat
    mfprintf(fids, 'Feat%2d: %s\n', j, label_groups{groupfeat(j)})
end
mfprintf(fids, '-------------------------\n')

ind = sub2ind(size(confmat), groupfeat, 1:nfeat);
acc = sum(confmat(ind));
err = 1-acc;
mfprintf(fids, 'Accuracy = %.2f\n', acc*100)
mfprintf(fids, 'Error = %.2f\n', err*100)

mfprintf(fids, '==========================\n')

fclose(fid);
