function BAcheckparams( n )
%BACHECKPARAMS checks the parameters of a Barabasi-Albert (BA) model

% Copyright (c) F. Caron (University of Oxford), A. Todeschini (Inria), and 
% X. Miscouridou (University of Oxford)
% caron@stats.ox.ac.uk
% adrien.todeschini@gmail.com
% xenia.miscouridou@spc.ox.ac.uk
% September 2015
%--------------------------------------------------------------------------

if ~isnumeric(n) || (floor(n)-n)~=0
    error('Parameter n must be an integer');
end

end

