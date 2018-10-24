function beta = scoreCGGPrnd(p, w0, Fdist, gamma)

% scoreCGGPrnd generates the scores beta that tune the levels of
% affiliations of the nodes to the communities
% They are drawn from distribution F conditionally on the base weight w_0, 
% the tilting parameter vector gamma and the number of communities p

% INPUT
%   - p: positive integer; vectors beta have size [1,p]
%   - w0: positive scalar; base weights drawn from the base Levy measure rho_0
%   - F: structure; specifies the distribution function of betas
%   - gamma: vector of size [p,1]; the positive tilting parameters 
%
% OUTPUT
%   - beta: vectors of dimension [1,p]; 
%


N = numel(w0);

% Sample beta
if N==0
    beta = zeros(0,p);
    return;
end
    
switch(Fdist.name)
    case 'gamma'
        
        if isnumeric(Fdist.param)
            a = Fdist.param';
            b = a;
        else
            a = Fdist.param.a';
            b = Fdist.param.b';
        end
        
        a = fillmat(a, [N, p]);
        b = fillmat(b, [1, p]);
        if any(gamma>0)
            b = bsxfun(@plus, b, w0*gamma');
        else
            b = repmat(b, N, 1);
        end
        
        beta = gamrnd(a, 1./b);
    otherwise
        error('Unknown Distribution %s', Fdist.name)
end
    