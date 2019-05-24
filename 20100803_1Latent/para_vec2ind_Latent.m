function [muQ, rhoQ, rho, delta1, delta0,SIGe] = para_vec2ind_Latent(para,mature,data)  
% para_vec2ind_Latent passes values from a vector 'para' to individual
% parameters, it supports single or multiple factor pure latent models
% INPUTS:   para:   a vector of parameters 
%           mature: a structure
%           data:   a structure
% OUTPUTS:  muQ, rhoQ, rho, delta1, delta0: ATSM parameters 
% last updated 8/3/2010

n1 = length(mature.exact);                                                  

nn = n1*(n1+1)/2;
rhoQ = veclt(para(1:nn));
nn_ = nn + n1; 
delta1   = para(nn + 1:nn_)/1e3;
nn = nn_ + n1^2;
rho      = reshape(para(nn_+1:nn),n1,n1);
nn_ = nn + n1;
muQ     = para(nn+1:nn_);
%muQ     = para(nn+1:nn_);

if nargout > 4
    if length(para) > nn_
        delta0 = para(nn_+1);
    elseif nargin < 3
        error('need more inputs')
    else
        mature_v  = [mature.exact, mature.error];
        obs_yield = [data.Y_1,data.Y_2];
        y         = obs_yield(:,(mature_v==1));                             % updated on Jul/07/09, because the first one is not necessary the short rate
        delta0  = mean(y);
    end
    if nargout > 5
        SIGe = diag(para(nn_+2:end))/1e4;
    end
end