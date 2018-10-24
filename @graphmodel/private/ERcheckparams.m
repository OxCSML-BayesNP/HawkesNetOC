function ERcheckparams( n, p )
%ERCHECKPARAMS checks the parameters of a Erdos-Renyi graph

% Copyright (c) F. Caron (University of Oxford), A. Todeschini (Inria), and 
% X. Miscouridou (University of Oxford)
% caron@stats.ox.ac.uk
% adrien.todeschini@gmail.com
% xenia.miscouridou@spc.ox.ac.uk
% September 2015
%--------------------------------------------------------------------------

if ~isnumeric(n) || (floor(n)-n)~=0
    error('First parameter n must be an integer');
end
if ~isnumeric(p) || p>1 || p<0
    error('Second parameter p must be a real in (0,1)')
end

end

