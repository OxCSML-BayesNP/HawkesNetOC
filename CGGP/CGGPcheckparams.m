function CGGPcheckparams(p, alpha, sigma, tau, Fdist, gamma)

% CGGPcheckparams checks the parameters of a Compound GGP.
% CGGPcheckparams(p, alpha, sigma, tau, Fdist, gamma)
%
%   Valid range is
%       alpha>0, tau>0, sigma in (-infinity, 1)
%           or
%       alpha>0, tau=0, sigma in (0,1)
%   if Fdist.name is 'gamma':
%       Fdist.param>0
%   and
%       gamma_k>0 for k=1...p
%
%   The function returns an error if the parameters are outside of this
%   range.

% Copyright (c) F. Caron (University of Oxford), A. Todeschini (Inria), and 
% X. Miscouridou (University of Oxford)
% caron@stats.ox.ac.uk
% adrien.todeschini@gmail.com
% xenia.miscouridou@spc.ox.ac.uk
% September 2017
%--------------------------------------------------------------------------

if p<=0 || p ~= round(p)
    error('p must be a positive integer')
end

GGPcheckparams(alpha, sigma, tau);

if any(gamma)<0
    error('All gamma components must be nonnegative')
end
if any(size(gamma)~=[p, 1])
    error('gamma must be a column vector of length p');
end

if ~isstruct(Fdist)
    error('Fdist must me a struct')
end
if ~isfield(Fdist, 'name')
    error('Fdist must have a field name')
end
if ~ischar(Fdist.name)
    error('Fdist.name must be a character string')
end
if ~isfield(Fdist, 'param')
    error('Fdist must have a field param')
end
switch(Fdist.name)
    case 'gamma'
        b = [];
        if isnumeric(Fdist.param)
            a = Fdist.param;
        elseif isstruct(Fdist.param)
            a = Fdist.param.a;
            b = Fdist.param.b;
        else
            error('Fdist.param must be either numeric or struct');
        end

        if size(a,1)~=1 && size(a,1)~=p
            error('a must have either one or p rows')
        end
        if any(a)<0
            error('a must have striclty positive entries')
        end

        if ~isempty(b)
            if (size(b,1)~=1 && size(b,1)~=p)
                error('b must have either one or p rows')
            end
            if any(b)<0
                error('b must have striclty positive entries')
            end
        end

    otherwise
        error('Unknown distribution %s for F', Fdist.name);
end

end
