% solves the just-identified pure-latent ATSM
% last updated aug/3/2010
clear
clc
load mature                       
load APdata
%% milti-stage MD
%--------------------------------------------------------------------------
% OLS regression to get reduced form parameters and variance covariance
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
bm = std(y2 - x2*para2);

%--------------------------------------------------------------------------
options  =  optimset('fsolve');
options  =  optimset(options , 'Display'     , 'iter');
%options  = optimset(options, 'Algorithm','levenberg-marquardt');           %updated 08/04/2009, faster for square system
options  = optimset(options, 'MaxFunEvals',100000);
options  = optimset(options, 'MaxIter',10000);
options  = optimset(options, 'TolFun',1e-20);

%--------------------------------------------------------------------------
% STEP 1: from OLS parameters to rhoQ and delta1
B1B1 = var_1;
B2B1 = phi_star21*var_1;

rhoQ = diag(sort(rand(3,1)/2+0.5,1,'descend'));
delta1 = ones(3,1)*1e-4;    
initial = [ltvec(rhoQ);delta1*1000];
para = fsolve(@(x) bn2rhoQdelta1(x,B1B1,B2B1,mature),initial,options);                                                                              
rhoQ = veclt(para(1:6));
delta1 = para(7:end)/1e3;
rhoQ(delta1 < 0,:) = -rhoQ(delta1 < 0,:);                                   % normalize delta1 to be positive
rhoQ(:,delta1 < 0) = -rhoQ(:,delta1 < 0);
delta1(delta1 < 0) = -delta1(delta1 < 0);
para = [ltvec(rhoQ);delta1*1000];
[f,B1,B2,B] = bn2rhoQdelta1(para,B1B1,B2B1,mature);
plot(B)

% STEP 2: solve rho analytically
rho = B1\phi_star11*B1;

% STEP 3: solve cQ and delta0 numerically
A1 = (eye(3)-phi_star11)\A1_star;
A2 = phi_star21*A1+A2_star;

delta0 = 0.0045;
cQ = [0;0;0];
initial = [delta0;cQ];
para = fsolve(@(x) an2cQdelta0(x,A1,A2,mature,rhoQ,delta1),initial,options);
cQ   = para(2:end);
delta0 = para(1);

save parameters_global cQ rhoQ rho delta0 delta1 bm
