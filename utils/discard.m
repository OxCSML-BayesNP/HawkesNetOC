function [samples] = discard(samples, nburn)

% discard discards the first burnin data
% from each of the dimensions of the structure samples
% and returns samples
%
% If samples are the MCMC output from an object of class graphmcmc,
% the function discard, discards the burn in part of the MCMC output for
% each chain 


names = fieldnames(samples);
nsamples = size(samples(1).(names{1}), ndims(samples(1).(names{1})));
[ntypes, nchains] = size(samples);

ind = nburn+1:nsamples;

for t=ntypes:-1:1 %% loop over types of nodes
    for j=1:numel(names)
        name = names{j};
        if isnumeric(samples(t,1).(name))
            for k=nchains:-1:1
                if ismatrix(samples(t,k).(name))
                    samples(t,k).(name) = samples(t,k).(name)(:, ind);
                else
                    samples(t,k).(name) = samples(t,k).(name)(:, :, ind);
                end
            end
        else
            fn = fieldnames(samples(t,1).(name));
            for v=1:numel(fn)
                for k=nchains:-1:1
                    if ismatrix(samples(t,k).(name).(fn{v}))
                        samples(t,k).(name).(fn{v}) = ...
                            samples(t,k).(name).(fn{v})(:, ind);
                    else
                        samples(t,k).(name).(fn{v}) = ...
                            samples(t,k).(name).(fn{v})(:, :, ind);
                    end
                end
            end
        end
    end
end
    