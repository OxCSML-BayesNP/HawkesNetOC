%% http://snap.stanford.edu/data/email-Eu-core.html

fid = fopen('sx-mathoverflow.txt');
A = textscan(fid, '%f%f%f');
fclose(fid);
times = A{3};
ind1 = A{1};
ind2 = A{2};
num_events = length(A{1});
num_nodes = numel(unique([A{2};unique(A{1})]));   %1899
G = sparse(ind1, ind2, ones(num_events,1));
G_symm = sparse([ind1;ind2],[ind2,ind1],ones(2*num_events,1));

times = times-times(1); %from unix to seconds
times = times/60/60/24; %from seconds to days

%%relabeling not needed
distinct_nodes = unique([ind1;ind2]);
numel(distinct_nodes)
max(distinct_nodes)
all_nodes=[ind1;ind2];

K=length(distinct_nodes);[~,~,conn]=find(triu(G_symm));
num_events = length(times);

meta =[];
meta.times =times;
meta.source =ind1;
meta.target = ind2;


%% save data
save forum.mat G meta


 
