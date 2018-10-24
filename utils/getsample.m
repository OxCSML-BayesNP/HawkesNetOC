function S = getsample(samples, ind)

% getsample accepts as input the array of structures samples 
% and extracts the data corresponding to the specific index for all 
% fields and across all dimensions of the size of the structure 
%
% output is S which is a structure of arrays with the extracted data from
% samples
%

% If samples are the MCMC samples output from an object of class graphmcmc,
% the function returns the samples of all parameters sat 
% the specific iteration ind .


names = fieldnames(samples);
sf = names';
sf(2,1:numel(names)) = {[]};
S(1:numel(samples)) = struct(sf{:});
S = reshape(S, size(samples));
for t=1:numel(samples) %% loop over types of nodes and chains
    for j=1:numel(names)
        name = names{j};
        if isnumeric(samples(t).(name))
            if ismatrix(samples(t).(name))
                S(t).(name) = samples(t).(name)(:,ind);
            else
                S(t).(name) = samples(t).(name)(:,:,ind);
            end
        else
            fn = fieldnames(samples(t).(name));
            for v=1:numel(fn)
                if ismatrix(samples(t).(name).(fn{v}))
                    S(t).(name).(fn{v}) = samples(t).(name).(fn{v})(:,ind);
                else
                    S(t).(name).(fn{v}) = samples(t).(name).(fn{v})(:,:,ind);
                end
            end
        end
    end
end
end
