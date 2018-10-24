function S = arraystruct2structarray(samples)

% arraystruct2structarray transforms samples, from a
% n x m x p array of structures 
% to a structure of  where each field contains n x m xp arrays

names = fieldnames(samples);
[ntypes, nchains] = size(samples);
nsamples = size(samples(1).(names{1}), ndims(samples(1).(names{1})));

for t=ntypes:-1:1 %% loop over types of nodes
    for ch=nchains:-1:1 %% loop over chains
        for j=1:numel(names)
            name = names{j};
            if isnumeric(samples(t,ch).(name))
                for i = nsamples:-1:1;
                    S(t,ch,i).(name) = samples(t,ch).(name)(:,:,i);
                end
            else
                fn = fieldnames(samples(t,ch).(name));
                for v = 1:numel(fn)
                    for i = nsamples:-1:1;
                        S(t,ch,i).(name).(fn{v}) = samples(t,ch).(name).(fn{v})(:,:,i);
                    end
                end
            end
        end
    end
end

end