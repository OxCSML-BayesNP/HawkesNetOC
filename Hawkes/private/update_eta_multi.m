function [eta, intensities] = update_eta_multi(eta, delta, mu, logR, S_multi, intensities, prop_eta, prior_eta, nmh)
% update_eta_multi implements the MH step for eta parameter
%
%
% -------------------------------------------------------------------------
% INPUTS
%   
%   - eta:  current estimate of eta
%   - delta: current estimate of delta
%   - mu: the base intensities
%   - logR: term involved in the intensity
%   - S_multi: term involved in the intensity
%   - intensities: intensities evaluated at the current estimates of eta
%   and delta
%   - prop_eta: struct with the name and parameters for the proposal
%   distribution for eta in the MH step
%   - prior_eta: struct with the name and parameters for the prior
%   distribution for eta
%   - nmh: number of MH steps
% 
% OUTPUTS
%   - eta:  the new value of eta
%   - intensities: intensities evaluated at the new estimates of eta
%   and delta
%
% -------------------------------------------------------------------------
% Copyright (C) Xenia Miscouridou, University of Oxford
% xenia.miscouridou@spc.ox.ac.uk
% October 2018
%--------------------------------------------------------------------------


    for nitermh=1:nmh

        %PROPOSAL  
        if strcmp(prop_eta.name, 'Exponential');
            etanew = exprnd(prop_eta.param);
            while etanew < 0 || etanew > delta    
                etanew = exprnd(prop_eta.param) ;
            end
            proposal_term = - prop_eta.param*(etanew-eta);
        
        else
            if strcmp(prop_eta.name, 'Normal');
                etanew =  (randn()*prop_eta.param + eta);
                while etanew < 0 || etanew>delta 
                % while etanew<0
                  etanew =  randn()*prop_eta.param + eta ;             
                end
                proposal_term =0;
            end
        end

    end
    
    prior_term = -prior_eta.param*(etanew-eta);
    %prior_term
    intensities_new =  bsxfun(@plus, S_multi.*etanew, mu);
    term1 = sum(sum(log(intensities_new) -log(intensities)));  
    term2 = -(etanew-eta)*exp(logR); 
    %MH step
    if  log(rand()) <  term1+term2 + proposal_term + prior_term
        eta = etanew;
        intensities = intensities_new;
    end

end