function [F,A] = an2cQdelta0(x,A1_data,A2_data,mature,rhoQ,delta1)
% an2cQdelta0 compares reduced form A1 and A2 with those implied from ATSM
% INPUTS:   x:              a vector of ATSM parameters determining A's,
%                           i.e. delta0 and cQ
%           A1_data,A2_Data:reduced form A's from OLS
%           mature:         a structure containing information on maturity
%                           mature.exact and mature.error
%           rhoQ, delta1:   ATSM parameters determining B's
% OUTPUTS:  F:              a scaled measure of difference between reduced
%                           form and structural A1 and A2
%           A:              a_n's for n from 1 to largest in mature
% LAST UPDATED 8/3/2010

sig = eye(3);

cQ   = x(2:end);
delta0 = x(1);

mature_ = [mature.exact,mature.error];
nN = max(max(mature.exact),max(mature.error));

A_bar_temp = -delta0;  %initial value for A_bar and B_bar
B_bar_temp = -delta1;
A = zeros(nN,1);
A(1)   = -A_bar_temp;

for i = 2:nN
    A_bar_temp = A_bar_temp + B_bar_temp'*cQ + 1/2*(B_bar_temp'*sig*sig'*B_bar_temp)-delta0;
    B_bar_temp = rhoQ'* B_bar_temp - delta1;
    A(i)   = -A_bar_temp/i;
end             

F = (A(mature_) -[A1_data;A2_data])*1e8;