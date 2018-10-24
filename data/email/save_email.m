%% http://snap.stanford.edu/data/email-Eu-core.html

fid = fopen('email-Eu-core.txt');
A = textscan(fid, '%f%f%f');
fclose(fid);
times=A{3};
ind1=A{1};
ind2=A{2};
num_events=length(A{1});
num_nodes=numel(unique([A{2};unique(A{1})]));   %986
ind1 = A{1}+1;          % since nodes are labeled from zero
ind2 = A{2}+1;
G = sparse(ind1, ind2, ones(num_events,1));

times = times/60/60/24;


%%relabeling
distinct_nodes = unique([ind1;ind2]);
N_nodes= numel(distinct_nodes);
all_nodes=[ind1;ind2];
id_nodes = interp1(1:max(distinct_nodes),distinct_nodes);


IND1 = NaN*zeros(length(ind1),1);
IND2 = IND1;
for k = 1:140
IND1(ind1==k )=k;
IND2(ind2==k)=k;
end

for k=141:986
    IND1(ind1==id_nodes(k) )=k;
    IND2(ind2==id_nodes(k) )=k;
end


meta =[];
meta.timepoints = times;
meta.source =IND1;
meta.target = IND2;
N_edges=nnz(G);
N_temporal=length(times);
fprintf('There are %d nodes, %d undirected edges and %d temporal edges. \n', N_nodes, N_edges, N_temporal)
save email.mat meta G

 
