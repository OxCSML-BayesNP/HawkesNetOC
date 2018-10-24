function [samples_comb] = combine(samples)

% combine accepts as input the variable samples, which is an
% array of structures
% and concatenates across the second dimension of the size of samples
% to return samples_comb which is a structure of arrays 

% If samples are the MCMC samples output from an object of class graphmcmc,
% the function combine, combines the output of several chains

names = fieldnames(samples);
nsamples = size(samples(1).(names{1}), ndims(samples(1).(names{1})));
[ntypes, nchains] = size(samples);

if nchains==1
    samples_comb = samples;
    return
end

% else (nchains>1)
for t=ntypes:-1:1 %% loop over types of nodes
    for j=1:numel(names)
        name = names{j};
        if isempty(samples(t,1).(name))
            samples_comb(t,1).(name) = [];
        else
            if isnumeric(samples(t,1).(name))
                for k=nchains:-1:1
                    ind = (k-1)*nsamples+1:k*nsamples;
                    if ismatrix(samples(t,k).(name))
                        samples_comb(t,1).(name)(:, ind) = samples(t,k).(name);
                    else
                        samples_comb(t,1).(name)(:, :, ind) = samples(t,k).(name);
                    end
                end
            else
                fn = fieldnames(samples(t,1).(name));
                for v=1:numel(fn)
                    for k=nchains:-1:1
                        ind = (k-1)*nsamples+1:k*nsamples;
                        if ismatrix(samples(t,k).(name).(fn{v}))
                            samples_comb(t,1).(name).(fn{v})(:, ind) = ...
                                samples(t,k).(name).(fn{v});
                        else
                            samples_comb(t,1).(name).(fn{v})(:, :, ind) = ...
                                samples(t,k).(name).(fn{v});
                        end
                    end
                end
            end
        end
    end
end
    