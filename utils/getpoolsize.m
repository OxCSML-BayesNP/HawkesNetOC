function n = getpoolsize()

% getpoolsize finds the number of workers in the current parallel pool
% and returns them in the output variable n

p = gcp('nocreate');
if isempty(p)
    n = 0;
else
    n = p.NumWorkers;
end
end