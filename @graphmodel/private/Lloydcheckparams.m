function Lloydcheckparams( n, sig, c, d )
%LLOYDCHECKPARAMS checks the parameters of a Erdos-Renyi graph


% Copyright (c) F. Caron (University of Oxford), A. Todeschini (Inria), and 
% X. Miscouridou (University of Oxford)
% caron@stats.ox.ac.uk
% adrien.todeschini@gmail.com
% xenia.miscouridou@spc.ox.ac.uk
% September 2015
%--------------------------------------------------------------------------

if ~isnumeric(n) || ~isnumeric(sig) || ~isnumeric(c) || ~isnumeric(d)
    error('Parameters must be numeric')
end
if n<=0 || sig<=0 || c<=0 || d<=0
    error('Parameters must be strictly positive');
end

end

