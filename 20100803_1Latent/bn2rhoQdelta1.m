function [F,B1,B2,B] = bn2rhoQdelta1(x,B1B1,B2B1,mature)
% function bn2rhoQdelta1 compares B1B2, B2B1 from reduced form and those
% from ATSM
% INPUTS:   x:          a vector of parameters containing rhoQ and delta1
%           B1B1,B2B1:  OLS estimates
%           mature:     a structure containing information relating to
%                       maturity, includeing mature.exact and mature.error
% OUTPUTS:  F:          a scaled vector evaluating the distance between
%                       reduced form and ATSM implied B1B1 and B2B1
%           B:          a matrix containing b_n for n from 1 to largest
%                       maturity in mature
% last updated 08/03/2010

n1 = length(mature.exact);
nn = n1*(n1+1)/2;
rhoQ   = veclt(x(1:nn));
delta1 = x(nn+1:end)/1e3; 
    
mature_ = [mature.exact mature.error];
nN = max(mature_ );
B_bar_temp = -delta1;
kk = length(rhoQ);
B = zeros(kk,nN);
B(:,1) = -B_bar_temp;

for i = 2:nN
    B_bar_temp = rhoQ'* B_bar_temp - delta1;
    B(:,i) = -B_bar_temp/i;      
end                                                    
B = B';
B1 = B(mature.exact,:);
B2 = B(mature.error,:);

temp1 = B1*B1' - B1B1;
temp2 = B2*B1'  - B2B1;

F = [ltvec(temp1);temp2(:)]*1200*1e7;
