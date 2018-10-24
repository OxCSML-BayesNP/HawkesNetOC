function samples = etstablernd(V0, alpha, tau, n)

%ETSTABLERND samples from the exponentially tilted stable distribution
% samples = etstablernd(V0, alpha, tau, n) returns a n*1 vector of numbers
% distributed fron an exponentially tilted stable distribution with Laplace
% transform (in z)
%           exp(-V0 * ((z + tau)^alpha - tau^alpha))
%
% Uses the algorithm proposed in (Devroye, 2009) with corrections pointed 
% out in (Hofert, 2011)
% -------------------------------------------------------------------------
% INPUTS
%   - V0: positive scalar
%   - alpha: real in (0,1)
%   - tau: positive scalar
%   - n: integer
%
% OUTPUT
%   - samples: vector of length n
% -------------------------------------------------------------------------
% EXAMPLE
% samples = etstablernd(100,0.5, 1, 100);

% References:
% - Luc Devroye. Random variate generation for exponentially and polynomially
% tilted stable distributions. ACM Transactions on Modeling and Computer
% Simulation, vol. 19(4), 2009.
% - Marius Hofert. Sampling exponentially tilted stable distributions. 
% ACM Transactions on Modeling and Computer Simulation, vol. 22(1), 2011.
%
% Copyright (C) Francois Caron, University of Oxford
% caron@stats.ox.ac.uk
% April 2015
%--------------------------------------------------------------------------

% Check parameters
if alpha<=0 || alpha>=1
    error('alpha must be in ]0,1[');
end
if tau <0
    error('tau must be >=0');
end
if V0<=0
    error('V0 must be >0');
end

if nargin<4
    n = 1;
end

% lambda = tau * V0^(1/alpha); % rescale
lambda_alpha = tau^alpha * V0; % lambda^a


% Now we sample from an exponentially tilted distribution of parameters
% sigma, lambda, as in (Devroye, 2009)
gamma = lambda_alpha * alpha * (1-alpha);

% xi = 1/pi *(2+sqrt(pi/2)) * sqrt(2*gamma) + 1;
xi = 1/pi *((2+sqrt(pi/2)) * sqrt(2*gamma) + 1); % Correction in Hofert
psi = 1/pi * exp(-gamma * pi^2/8) * (2 + sqrt(pi/2)) * sqrt(gamma * pi);
w1 = xi * sqrt(pi/2/gamma);
w2 = 2 * psi * sqrt(pi);
w3 = xi * pi;
b = (1-alpha)/alpha;

samples = zeros(n, 1);
for i=1:n
%     i
    while 1
        % generate U with density g*/G*
        while 1
            % Generate U with density proportional to g**
            U = gen_U(w1, w2, w3, gamma);

            W = rand;
            zeta = sqrt(ratio_B(U, alpha));
            z = 1/(1 - (1 + alpha*zeta/sqrt(gamma))^(-1/alpha));
            rho = pi * exp(-lambda_alpha * (1-zeta^(-2))) ...
                * (xi * exp(-gamma*U^2/2) * (U>=0)*(gamma>=1) + psi/sqrt(pi-U)* (U>0)*(U<pi) + xi *(U>=0)*(U<=pi)*(gamma<1))...
                /((1 + sqrt(pi/2)) *sqrt(gamma)/zeta + z);
            if (U<pi && W*rho<=1)
                break;
            end        
        end

        % Generate X with density proportional to g(x, U)
        a = zolotarev(U, alpha);
        m = (b/a)^alpha * lambda_alpha;
        delta = sqrt(m*alpha/a);
        a1 = delta * sqrt(pi/2);
%         a2 = delta;
        a2 = a1 + delta; % correction in Hofert
        a3 = z/a;
        s = a1 + delta + a3;% correction in Hofert
        V_p = rand;    
        N_p = randn;
        E_p = -log(rand);
        if V_p<a1/s
            X = m - delta*abs(N_p);
        elseif V_p<a2/s
            X = delta * rand + m;
        else
            X = m + delta + a3 * E_p;
        end
%         X
        E = -log(rand);
%         cond = (a*(X-m) + lambda*(X^(-b) - m^(-b)) - N_p^2/2 * (X<m) - E_p * (X>m+delta));
        cond = (a*(X-m) + exp(1/alpha*log(lambda_alpha)-b*log(m))*((m/X)^b - 1) - N_p^2/2 * (X<m) - E_p * (X>m+delta));
        if ((X>=0) && (cond <=E))
            break;
        end   
    end
    samples(i) = exp( 1/alpha* log(V0) -b*log(X));% more stable than V0^(1/alpha) * X^(-b);
end

end

function U = gen_U(w1, w2, w3, gamma)

V = rand;
W_p = rand;
if gamma>=1
    if (V < w1/(w1+w2))
        U = abs(randn) /sqrt(gamma);
    else
        U = pi * (1 - W_p^2);
    end
else
    if (V < w3/(w3 + w2))
        U = pi * W_p;
    else
        U = pi * (1 - W_p^2);
    end
end
end

function out = ratio_B(x, sigma)

out = sinc(x) / (sinc(sigma * x))^sigma / (sinc((1-sigma)*x))^(1-sigma);
end 

function out = sinc(x)
out = sin(x)/x;
end

function out = zolotarev(u, sigma)
% Zolotarev function, cf (Devroye, 2009)
out = ((sin(sigma*u))^sigma * (sin((1-sigma)*u))^(1-sigma) / sin(u))^(1/(1-sigma));

end
