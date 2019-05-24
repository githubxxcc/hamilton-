function para = para_ind2vec_Latent(muQ, rhoQ, rho, delta1, delta0,SIGe)
% para_ind2vec_Latent converts individual parameters into a vector
% last updated on Jul/22/09

para = [ltvec(rhoQ);delta1*1e3;rho(:);muQ;delta0];
if nargin >5 
    para = [para; diag(SIGe)*1e4];
end