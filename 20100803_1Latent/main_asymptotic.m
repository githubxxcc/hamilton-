% last updated aug/10/2009
clear
clc
cd('Z:\RESEARCH\HW1_20090706_identification and estimation\20100803_clean code\20100803_1Latent');

load mature;                           
load APdata
load parameters_global

%% --------------------------------------------------------------------------
% OLS regression to get regression parameters and variance covariance
% matrix
Y_1 = data.Y_1;
Y_2 = data.Y_2(:,end);
T = length(Y_1)-1;

y1 = Y_1(2:end,:);
x1 = [ones(T,1),Y_1(1:end-1,:)];
para1 = x1\y1;
A1_star  = para1(1,:)';
phi_star11 = para1(2:end,:)';
u1t_star = y1 - x1*para1;
var_1 = var_OLS(u1t_star);

y2 = Y_2(2:end,:);
x2 = [ones(T,1),Y_1(2:end,:)];
para2 = x2\y2;
A2_star = para2(1,:)';
phi_star21 = para2(2:end,:)';
SIGe = std(y2 - x2*para2);
SIGe =diag(SIGe);
%% information matrix and asymptotic std error
H1 = kron(inv(var_1),x1'*x1);
d = duplication(size(var_1,1));
H1v = T/2*d'*kron(inv(var_1),inv(var_1))*d;
H2 = kron(inv(diag(SIGe.^2)),x2'*x2);
d = duplication(size(SIGe,1));
H2v = T/2*d'*kron(inv(diag(SIGe.^2)),inv(diag(SIGe.^2)))*d;

R = blkdiag(H1,H1v,H2,H2v);
paraSTR_vec = para_ind2vec_Latent(cQ, rhoQ, rho, delta1, delta0,SIGe);
GAMMA = numgrad('para_str2OLS',paraSTR_vec,mature);

variance = inv(GAMMA'*R*GAMMA);
para_std = sqrt(diag(variance));
[cQ_std, rhoQ_std, rho_std, delta1_std, delta0_std,SIGe_std] = para_vec2ind_Latent(para_std,mature); 
save asympt_std rhoQ_std delta1_std rho_std cQ_std delta0_std SIGe_std 