%% http://snap.stanford.edu/data/

fid = fopen('sx-askubuntu.txt');
A = textscan(fid, '%f%f%f');
fclose(fid);
times=A{3};
ind1=A{1};
ind2=A{2};
num_events=length(A{1});
num_nodes=numel(unique([A{2};unique(A{1})]));   %1899
G = sparse(ind1, ind2, ones(num_events,1));
G_symm=sparse([ind1;ind2],[ind2,ind1],ones(2*num_events,1));
[a,b]=sort(times);

times=a;
ind1=ind1(b);
ind2=ind2(b);
times=times-times(1); %from unix to seconds
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


source = meta.source;

[source, order]=sort(source);
target = meta.target(order);
timepoints = meta.times(order);


%% FULL DATASET
% shift all times by 0.01
timepoints=timepoints+0.01;
max_T = floor(max(timepoints))+1; %IN DAYS

full_data=[source,target,timepoints];
full_max_T=max_T;

fprintf('Full data set has %d data points \n\n', size(full_data,1));
fprintf('The total maximum time is %d \n \n', (max_T))  ;  
%G=sparse(source,target,ones(length(source),1));

G = sparse([source;target],[target;source],ones(length([source;target]),1));

fn = fieldnames(meta);

% Remove nodes with no edge 
ind = any(G);
G = G(ind, ind);
for i=1:length(fn)
    meta.(fn{i}) = meta.(fn{i})(ind);
end
map=find(full(ind));
G_bin_symm = G_symm|G_symm';


N_und_edges=nnz(G_bin_symm);
[source_und, target_und, N_links_per_und_pair] = find(triu(G_symm)); % between source_und(1) AND target_und(1)there are N_links_per_und_pair(1) links
N_pairs = numel(source_und);N_nodes=size(G_symm,1);N_temporal_edges=numel(sum(N_links_per_und_pair));
max_N_links_among_pairs = max(N_links_per_und_pair);

% bring data in a different form
times_matrix = NaN*zeros(N_pairs, max_N_links_among_pairs);
max_event_per_pair=zeros(N_pairs,1);min_event_per_pair=zeros(N_pairs,1);
source_target_times_matrix=zeros(N_pairs,2);
for l=1:N_pairs
	if rem(l,1000)==0
		fprintf('Processing edge/pair %d \n', l)    
    end  
    event_times = timepoints(find((source==source_und(l)&target==target_und(l))));
    event_times2 = timepoints(find((source==target_und(l)&target==source_und(l))));
    events=[event_times;event_times2];
    N_links_per_und_pair(l)= numel(events);
    if ~isempty(events)
    times_matrix(l,1:N_links_per_und_pair(l)) = (events);
    max_event_per_pair(l)=max(events);
    min_event_per_pair(l)=min(events);
    source_target_times_matrix(l,1)=source_und(l);
    source_target_times_matrix(l,2)=target_und(l);
    end
end

%% save data
save askubuntu.mat G meta times_matrix source_target_times_matrix map


 
