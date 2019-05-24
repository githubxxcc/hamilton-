function [muQ, rhoQ, rho, delta1, delta0] = para_str2ind_Latent(para_str)
% para_str2ind_Latent converts a structure para_str to 
% individual parameters 
% last updated on April/5/10

delta0      = para_str.delta0;
delta1      = para_str.delta1;
muQ         = para_str.muQ;
rhoQ        = para_str.rhoQ;
rho         = para_str.rho;


