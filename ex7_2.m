%--------------------------------------------------------------------------
%   
%          ADJUSTMENT THEORY I
%    Exercise 7: Adjustment Calculation - part II  
% 
%   Author         : Anastasia Pasioti
%   Version        : October 09, 2018
%   Last changes   : December 13, 2021
%
%--------------------------------------------------------------------------

clc;
clear all;
close all;

%--------------------------------------------------------------------------
%   Task 1: Adjustment of a straight line
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
%   Observations and redundancy
%--------------------------------------------------------------------------
%Vector of observations 
L = [1.112;.880;.768;.830;1.175];

%Number of observations
no_n =  length(L);

%Number of unknowns
no_u = 3;

%Redundancy
r =no_n-no_u;

%--------------------------------------------------------------------------
%  Stochastic model
%--------------------------------------------------------------------------
%VC Matrix of the observations
S_LL = diag([.02^2]*ones(no_n,1)); %eye(no_n) * 0.02^2;

%Theoretical reference standard deviation
sigma_0 =1 ;      %a priori

%Cofactor matrix of the observations
Q_LL = 1/(sigma_0^2) * S_LL;

%Weight matrix
P = inv(Q_LL);

%--------------------------------------------------------------------------
%  Adjustment
%--------------------------------------------------------------------------
%Design matrix
X=[1;2;3;4;5];
A = [X.^2 X ones(no_n ,1)];
    
%Normal matrix
N = A' * P * A;
        
%Vector of right hand side of normal equations
n = A' * P * L;
    
%Inversion of normal matrix / Cofactor matrix of the unknowns
Q_XX = inv(N);
    
%Solution of normal equation
X_hat = Q_XX * n;
    
%Estimated unknown parameters
a = X_hat(1);
b = X_hat(2);
c = X_hat(3);
 
%Vector of residuals
v = A * X_hat - L;

%Objective function
vTPv = v' * P * v;

%Vector of adjusted observations
L_hat = L + v;

%Final check

  
%Empirical reference standard deviation
s_0 = sqrt(vTPv/r);       %a posteriori

%VC matrix of adjusted unknowns
S_XX_hat =s_0^2 * Q_XX;

%Standard deviation of the adjusted unknowns
s_X = sqrt(diag(S_XX_hat));

%Cofactor matrix of adjusted observations
Q_LL_hat = A * Q_XX * A';

%VC matrix of adjusted observations
S_LL_hat = s_0^2 * Q_LL_hat;

%Standard deviation of the adjusted observations
s_L_hat = sqrt(diag(S_LL_hat));

%Cofactor matrix of the residuals
Q_vv = Q_LL - Q_LL_hat;

%VC matrix of residuals
S_vv = s_0^2 * Q_vv;

%Standard deviation of the residuals
s_v = sqrt(diag(S_vv));


a
b
c




