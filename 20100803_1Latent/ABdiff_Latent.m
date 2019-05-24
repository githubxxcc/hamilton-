function [A1,A2,B1,B2] = ABdiff_Latent(para_str, mature) 
% ABdiff_Latent calculates A's and B's with difference equations
% USAGE: [A1,A2,B1,B2] = ABdiff_Latent(para_str, nN) 
% INPUTS:   para_str:   a structure of parameters
%           mature:     a structure containing mature.exact and
%                       mature.error
% OUTPUTS:  A1:         3*1 vector
%           A2:         Ne*1 vector
%           B1:         3*3 matrix
%           B2:         Ne*3 matrix
% last updated on April/5/10

nN = max(max(mature.exact),max(mature.error));
[cQ, rhoQ, rho, delta1, delta0] = para_str2ind_Latent(para_str);
NN = length(rho);

A_bar_temp = -delta0;  %initial value for A_bar and B_bar
B_bar_temp = -delta1;

A = zeros(nN,1);
B = zeros(NN,nN);

A(1)   = -A_bar_temp;  %A1 and B1
B(:,1) = -B_bar_temp;

for i = 2:nN
    A_bar_temp = A_bar_temp + B_bar_temp'*cQ + 1/2*(B_bar_temp'*B_bar_temp)-delta0;
    B_bar_temp = rhoQ'* B_bar_temp - delta1;
    A(i)   = -A_bar_temp/i;
    B(:,i) = -B_bar_temp/i;      
end                                                    

B = B';

A1 = A(mature.exact);
A2 = A(mature.error);
B1 = B(mature.exact,:);
B2 = B(mature.error,:); 