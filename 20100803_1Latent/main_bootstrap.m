% last updated aug/10/2009
clear
clc
cd('Z:\RESEARCH\HW1_20090706_identification and estimation\20100803_clean code\20100803_1Latent');

options  =  optimset('fsolve');
options  =  optimset(options , 'Display'     , 'off');
%options  = optimset(options, 'Algorithm','levenberg-marquardt');           %updated 08/04/2009, faster for square system
options  = optimset(options, 'MaxFunEvals',100000);
options  = optimset(options, 'MaxIter',10000);
options  = optimset(options, 'TolFun',1e-20);

load mature;           
load APdata                
load parameters_global


%--------------------------------------------------------------------------
% OLS regression to get regression parameters and variance covariance
% matrix
Y_1 = data.Y_1;
Y_2 = data.Y_2(:,end);
T = length(Y_1)-1;

y1 = Y_1(2:end,:);
x1 = [ones(T,1),Y_1(1:end-1,:)];
para1 = x1\y1;
u1t_star = y1 - x1*para1;
var_1 = var_OLS(u1t_star);

y2 = Y_2(2:end,:);
x2 = [ones(T,1),Y_1(2:end,:)];
para2 = x2\y2;
bm = std(y2 - x2*para2);

%% information matrix
H1 = kron(inv(var_1),x1'*x1);
std_1 = sqrt(diag(inv(H1)));
d = duplication(size(var_1,1));
H1v = T/2*d'*kron(inv(var_1),inv(var_1))*d;
std_1v = sqrt(diag(inv(H1v)));
H2 = kron(inv(diag(bm.^2)),x2'*x2);
std_2 = sqrt(diag(inv(H2)));
d = duplication(1);
H2v = T/2*d'*kron(inv(diag(bm.^2)),inv(diag(bm.^2)))*d;
std_2v = sqrt(diag(inv(H2v)));

%[std_1';std(para1_vec)]
%[std_2';std(para2_vec)]
%[std_1v',std_2v';std(var_vec)]

%%
n = 1000;

bm_vec = zeros(n,1);

rhoQdelta1_vec = zeros(n,9);
rho_vec = zeros(n,9);
cQdelta0_vec = zeros(n,4);
SSE_vec = zeros(n,1);
disorder_vec = zeros(n,1);
lam0_vec = zeros(n,3);
lam1_vec = zeros(n,9);

%{
para1_vec = zeros(n,12);
para2_vec = zeros(n,4);
var_vec = zeros(n,7);
%}
tic
for j = 1:n
    j
    Y1_j = zeros(T+1,3);
    Y2_j = zeros(T+1,1);
    Y1_j(1,:) = Y_1(1,:);
    Y2_j(1,:) = Y_2(1,:);
    U = randn(T,4);
    U = U*blkdiag(chol(var_1),diag(bm));
    
    u1t_star_j = U(:,1:3);
    u2t_star_j = U(:,4:end);
    
    for t = 1:T
        Y1_j(t+1,:) = [1,Y1_j(t,:)]* para1 + u1t_star_j(t,:);
        Y2_j(t+1,:) = [1,Y1_j(t+1,:)]* para2 + u2t_star_j(t,:);
    end
    
    y1 = Y1_j(2:end,:);
    x1 = [ones(T,1),Y1_j(1:end-1,:)];
    para1_j = x1\y1;
    A1_star  = para1_j(1,:)';
    phi_star11 = para1_j(2:end,:)';
    u1t_star = y1 - x1*para1_j;
    var_1j = var_OLS(u1t_star);
    
    y2 = Y2_j(2:end,:);
    x2 = [ones(T,1),Y1_j(2:end,:)];
    para2_j = x2\y2;
    A2_star = para2_j(1,:)';
    phi_star21 = para2_j(2:end,:)';
    bm_j = std(y2 - x2*para2_j);
    %{
    para1_vec(j,:) = para1_j(:)';
    para2_vec(j,:) = para2_j(:)';
    var_vec(j,:) = [ltvec(var_1j);(bm_j'.^2)]';
    %}
    
    bm_vec(j,:) = bm_j; 
    B1B1 = var_1j;
    B2B1 = phi_star21*var_1j;
    
    initial = [ltvec(rhoQ);delta1*1000];
    para = fsolve(@(x) bn2rhoQdelta1(x,B1B1,B2B1,mature),initial,options);   
    rhoQ_j = veclt(para(1:6));
    delta1_j = para(7:end)/1e3;   
    [f,B1,B2] = bn2rhoQdelta1(para,B1B1,B2B1,mature);
    
    rhoQdelta1_vec(j,:) = para';
    
    if rhoQ_j(1,1)<rhoQ_j(2,2) ||rhoQ_j(2,2)<rhoQ_j(3,3)
        warning('rhoQ disorder!!!')
        disorder_vec(j) = 1;
    end
    
    rho_j = B1\phi_star11*B1;
    rho_vec(j,:) = rho_j(:)';
    SSE_vec(j) = sum(f.^2);
    
    A1 = (eye(3)-phi_star11)\A1_star;
    A2 = phi_star21*A1+A2_star;

    initial = [delta0;cQ];
    para = fsolve(@(x) an2cQdelta0(x,A1,A2,mature,rhoQ_j,delta1_j),initial,options);
    
    cQ_j   = para(2:end);
    delta0_j = para(1);
    cQdelta0_vec(j,:) = para';
    
    lam0_vec(j,:) = -cQ_j';
    lam1_vec(j,:) = rho_j(:)'-rhoQ_j(:)';
end
toc

max(disorder_vec)
max(SSE_vec)
rhoQdelta1_true = [ltvec(rhoQ);delta1*1000];
rhoQdelta1_std = sqrt(sum((rhoQdelta1_vec - repmat(rhoQdelta1_true',n,1)).^2)/n);
rhoQ_std = veclt(rhoQdelta1_std(1:6));
delta1_std = rhoQdelta1_std(7:end)/1e3;

rho_std = sqrt(sum((rho_vec - repmat(rho(:)',n,1)).^2)/n);
rho_std = reshape(rho_std,3,3);

cQdelta0_true = [delta0;cQ];
cQdelta0_std = sqrt(sum((cQdelta0_vec - repmat(cQdelta0_true',n,1)).^2)/n);
delta0_std = cQdelta0_std(1);
cQ_std = cQdelta0_std(2:end);
bm_std = sqrt(sum((bm_vec - repmat(bm',n,1)).^2)/n);

lam0_true = -cQ;
lam1_true = rho(:)'-rhoQ(:)';
lam0_std = sqrt(sum((lam0_vec - repmat(lam0_true',n,1)).^2)/n);
lam1_std = sqrt(sum((lam1_vec - repmat(lam1_true,n,1)).^2)/n);
lam1_std = reshape(lam1_std,3,3);

save bootstrap_std rhoQ_std delta1_std rho_std cQ_std delta0_std bm_std lam0_std lam1_std