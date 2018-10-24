function G = Lloydgraphrnd(N,sig,c,d)

%Lloydgraphrnd samples a graph from Lloyd et al. network model
% G = Lloydgraphrnd(N,sig,c,d)
%
% -------------------------------------------------------------------------
% INPUTS
%   - N: number of nodes
%   - sig: nugget noise
%   - c: RBF length-scale
%   - d: RBF scale param
%

% Reference
% J.R. Lloyd, P. Orbanz, Z. Ghahramani, D.M. Roy. Random function priors
% for exchangeable arrays. NIPS, 2012.

% Copyright (c) F. Caron (University of Oxford), A. Todeschini (Inria), and 
% X. Miscouridou (University of Oxford)
% caron@stats.ox.ac.uk
% adrien.todeschini@gmail.com
% xenia.miscouridou@spc.ox.ac.uk

% Draw N nodes:
U = rand(N,1);
U = sort(U,'ascend');
U_array = zeros(2,N*(N-1)/2);
count = 1;
for ii=1:N
    for jj=1:ii-1
        U_array(:,count) = [U(ii);U(jj)];
        count = count + 1;
    end
end

%%
% Define GP kernel:
nn=size(U_array,2);
K_tmp = zeros(nn,nn);
for ii=1:nn
    for jj=1:ii
        K_tmp(ii,jj) = d*exp(-c*norm(U_array(:,ii) - U_array(:,jj))^2) + sig*min(U_array(:,ii)==U_array(:,jj));
    end
end


%%
% Draw a GP sample
Theta = chol(K_tmp,'lower')'*randn(nn,1);

W = 1./(1+exp(Theta));
links = rand(nn,1) < W;

G = zeros(N);
count = 1;
for ii=1:N
    for jj=1:ii-1
        G(ii,jj) = links(count);
        G(jj,ii) = links(count);
        count = count + 1;
    end
end