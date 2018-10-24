function r = mnrnd2(n,p)

% mnrnd2 generates random vectors from the multinomial distribution.
%
% INPUTS
% - n: vector of numbers of trials
% - p: matrix of event probabilities with the same number of rows as n
% OUTPUT
% - r: matrix of sampled counts with the same size as p

ind = count2ind(n); % counts transformed into vector of duplicated indices
x = mn1rnd(p(ind,:)); % sample discrete variables

% reduce x by summing for each index
r = zeros(numel(n)+1, size(p,2));
binedges = 0.5:1:(numel(n)+0.5);
for k=1:size(p,2)
    if any(x(:,k))
        r(:,k) = histc(ind(x(:,k)), binedges)';
    end
end

r(end,:) = [];

end

function x = count2ind(n)
%COUNT2IND Inverse of HISTC
% - n: counts for categorical values
% OUTPUT
% - x: vector of values where each i is repeated n(i) times
ind = [1; cumsum(n)+1];
ind(end) = [];
c = false(sum(n),1);
c(ind) = true;
x = cumsum(c);
end

function r = mn1rnd(p)
%MN1RND Random vectors from the multinomial distribution with n=1
u = rand(size(p,1), 1); % uniform random var
r = indic(u, p); % logical indicator
end

function x = indic(u, p)
%INDIC Bin indicator vectors
cs = cumsum(p, 2);
under = bsxfun(@lt, u, cs);
cs(:,end) = [];
over = bsxfun(@ge, u, cs);
x = [under(:,1), under(:,2:end) & over];
end
