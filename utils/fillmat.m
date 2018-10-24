function B = fillmat(A, sz)

% fillmat Replicates the input matrix A to obtain a matrix B of size sz.
%
% sz is a vector 
% A is a matrix 

sz(1:ndims(A)) = sz./size(A);
B = repmat(A, sz);
end