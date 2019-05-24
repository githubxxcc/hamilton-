function paraOLS_vec = para_str2OLS(paraSTR_vec,mature)

[cQ, rhoQ, rho, delta1, delta0,SIGe] = para_vec2ind_Latent(paraSTR_vec,mature);
para_str = para_ind2str_Latent(cQ, rhoQ, rho, delta1, delta0);
[A1,A2,B1,B2] = ABdiff_Latent(para_str, mature);

A1_star = A1 - B1*rho/B1*A1;
phi_star11 = B1*rho/B1;
Omega_star1 = B1*B1';
A2_star = A2 - B2/B1*A1;
phi_star21 = B2/B1;
Omega_star2 = SIGe*SIGe';

temp1 = [A1_star';phi_star11'];
temp2 = [A2_star';phi_star21'];
paraOLS_vec = [temp1(:);ltvec(Omega_star1);temp2(:);diag(Omega_star2)];

