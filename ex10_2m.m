%--------------------------------------------------------------------------
%   
%          ADJUSTMENT THEORY I
%  Exercise 10: Adjustment Calculation - part V  
% 
%   Author         : Anastasia Pasioti
%   Version        : October 11, 2018
%   Last changes   : January 17, 2022
%
%--------------------------------------------------------------------------

clc;
clear all;
close all;

%--------------------------------------------------------------------------
%   Task 1
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%   Observations and redundancy
%--------------------------------------------------------------------------
%Height differences/observations
L =[10.509; 5.360; -8.523; -7.348; -3.167; 15.881];        %[m]

%Benchmarks - error free
ha=350.908;
  l1=[ha; 0; 0; -ha; 0; ha] ;               %[m]

%Vector of observations with benchmarks
L_dash = L+l1;

%Number of observations
no_n = length(L);

%Number of unknowns
no_u = 3;

%Redundancy
r =  no_n-no_u;

%--------------------------------------------------------------------------
%  Stochastic model
%--------------------------------------------------------------------------
%VC Matrix of the observations
sh=[6 4 5 3 4 12].*10^-3;
S_LL = diag(sh).^2;

%Theoretical standard deviation
sigma_0 = 1     %a priori

%Cofactor matrix of the observations
Q_LL = 1 / sigma_0^2 * S_LL;

%Weight matrix
P = inv(Q_LL);

%--------------------------------------------------------------------------
%  Adjustment
%--------------------------------------------------------------------------
%Design matrix
A = [1 0 0;-1 1 0; 0 -1 1;0 0 -1; -1 0 1; 0 1 0];
    
%Normal matrix
N = A' * P * A;
        
%Vector of right hand side of normal equations
n = A' * P * L_dash;
    
%Inversion of normal matrix / Cofactor matrix of the unknowns
Q_XX = inv(N);
    
%Solution of normal equation
X_hat = Q_XX * n;
    
%Estimated unknown parameters
hb=X_hat(1);
hc=X_hat(2);
hd=X_hat(3);
     
%Vector of residuals
v = A * X_hat - L_dash;

%Objective function
vTPv = v' * P * v;

%Vector of adjusted observations
L_hat = L + v;

L_hat2 = L_dash + v;
%Final check    (to check for computational errors)
Phi_X_hat = [hb;hc-hb;hd-hc;-hd;hd-hb;hc];

if L_hat-Phi_X_hat<10^-2
  disp('Everything is fine!')
else
  disp('Something is wrong.')
end

%Empirical reference standard deviation
s_0 = sqrt(vTPv / r);     %a posteriori

%VC matrix of adjusted unknowns
S_XX_hat = s_0^2 * Q_XX;

%Standard deviation of the adjusted unknowns
s_X =sqrt(diag(S_XX_hat)); 

%Cofactor matrix of adjusted observations
Q_LL_hat = A * Q_XX * A';

%VC matrix of adjusted observations
S_LL_hat = s_0^2 * Q_LL_hat;

%Standard deviation of the adjusted observations
s_L_hat =sqrt(diag(S_LL_hat)); 

%Cofactor matrix of the residuals
Q_vv = Q_LL - Q_LL_hat;

%VC matrix of residuals
S_vv = s_0^2 * Q_vv;

%Standard deviation of the residuals
s_v = sqrt(diag(S_vv));

hb
hc
hd