function print_summary(filename, title, G, niter, nburn, nchains, thin, p, outpath, tinit)

% Summary file
fid = fopen(fullfile(outpath, filename), 'w');
fprintf(fid, '------------------------------------------------------\n');
fprintf(fid, '-- %s --\n', title);
fprintf(fid, '------------------------------------------------------\n');
fprintf(fid, 'Start: %s\n', datestr(tinit));
fprintf(fid, 'End: %s\n', datestr(clock));
fprintf(fid, '------------------------------------------------------\n');
fprintf(fid, 'nnodes = %d x %d\n', size(G));
nedges = nnz(G);
nmiss = full(sum(isnan(G(:))));
nones = nedges - nmiss;
nzeros = numel(G) - nedges;
fprintf(fid, 'nedges = %d\n', nedges);
fprintf(fid, 'nones = %d\n', nones);
fprintf(fid, 'nzeros = %d\n', nzeros);
fprintf(fid, 'nmissing = %d\n', nmiss);
fprintf(fid, '------------------------------------------------------\n');
fprintf(fid, 'niter = %d\nnburn = %d\nnchains = %d\nthin = %d\n',niter,nburn,nchains,thin);
fprintf(fid, '------------------------------------------------------\n');
fprintf(fid, 'nfeatures = %d\n',p);
fclose(fid);
