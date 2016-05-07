function [ Distance ] = ComputeCDFDistance( f_1, f_2 )
%COMPUTECDFDISTANCE Computes the L1Norm distance between the 
%CDFs of f_1, and f_2
%   Computes the L1 Norm distance between two distributions f_1 and f_2
% Inputs: 
%   f_1: 1st distribution 
%   f_2: 2nd distribution
%
% Author: Lewis Li (lewisli@stanford.edu)
% Based on file: L1cdfDiff.m by Addy Satija

[f1,x1]=ecdf(f_1);
[f2,x2]=ecdf(f_2);
[a1,b1]=unique(x1);



f22=interp1(x1(b1),f1(b1),x2,'linear','extrap');

f22(f22>1)=1;
f22(f22<0)=0;

Distance=sum(abs(f2-f22))/length(f22);

end

