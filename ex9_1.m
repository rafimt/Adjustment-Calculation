%--------------------------------------------------------------------------
%   
%          ADJUSTMENT THEORY I
%    Exercise 8: Adjustment Calculation - part III  
% 
%   Author         : Anastasia Pasioti
%   Version        : October 09, 2018
%   Last changes   : January 03, 2022
%
%--------------------------------------------------------------------------

clc;
clear all;
close all;

%--------------------------------------------------------------------------
%   Task 1 - Non-linear equation system
%--------------------------------------------------------------------------
disp('Task 1 - Non-linear adjustment problem!')

%--------------------------------------------------------------------------
%   Observations and initial values for the unknowns
%--------------------------------------------------------------------------
data=load('Exercise9Task1.txt');
% error free value
t=data(:,1);

%Vector of observations
L = data(:,2);

%Number of observations
no_n =length(L); 

%Initial values for the unknowns
a=3;
f =1/pi;
phi =0.3;
%Vector of initial values for the unknowns
X_0 = [a;f;phi];

%Number of unknowns
no_u = length(X_0);

%Redundancy
r =  no_n-no_u;

%--------------------------------------------------------------------------
%  Stochastic model
%--------------------------------------------------------------------------
%VC Matrix of the observations
sy=0.15;
S_LL = eye(no_n)*sy^2;

%Theoretical standard deviation
sigma_0 =  1;    %a priori

%Cofactor matrix of the observations
Q_LL = 1/sigma_0^2 * S_LL;

%Weight matrix
P = inv(Q_LL);

%--------------------------------------------------------------------------
%  Adjustment
%--------------------------------------------------------------------------
%breakoff conditions
epsilon = 10^-8;
delta = 10^-12;
max_x_hat = 10^Inf;

%Number of iterations
iteration = 0;

 while max_x_hat>epsilon || Check2>delta          
    
     %Observations as functions of the approximations for the unknowns
     L_0 = a*sin(2*pi*f*t+phi);
   

     
     %Vector of reduced observations
     l = L - L_0;
     
    
     %Design matrix with the elements from the Jacobian matrix J
     A = [sin(2*pi*f*t+phi), a*cos(2*pi*f*t+phi)*2*pi.*t, a*cos(2*pi*f*t+phi)];
   
     %Normal matrix
     N = A' * P * A;
     
     %Vector of right hand side of normal equations
     n = A' * P * l;
    
     %Inversion of normal matrix / Cofactor matrix of the unknowns
     Q_xx = N^(-1);
    
     %Solution of the normal equations
     x_hat = Q_xx * n;
       
     %Update
     X_0 =  X_0 + x_hat; 
    
    a = X_0(1);
     f = X_0(2);
     phi = X_0(3);
    
     %Check 1
     max_x_hat =max(abs(x_hat)); 
     
     %Vector of residuals
     v = A*x_hat - l;
 
     %Vector of adjusted observations
     L_hat = L + v;
    
     %Objective function
     vTPv = v' * P * v;
     phi_X_hat = a*sin(2*pi*f*t+phi);
    


     %Check 2
     Check2 = max(abs(L_hat-phi_X_hat));
    
     %Update number of iterations
     iteration = iteration+1;
  
 end



%Empirical reference standard deviation
s_0 = sqrt(vTPv/r);

%VC matrix of adjusted unknowns
S_XX_hat = s_0^2 * Q_xx;

%Standard deviation of the adjusted unknowns
s_X = sqrt(diag(S_XX_hat));

%Cofactor matrix of adjusted observations
Q_LL_hat =  A * Q_xx * A';

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
f
phi

