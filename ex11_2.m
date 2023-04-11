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
%Coordinates - error free values

% error free value
X=[4316.175;4036.242;2136.262;1324.177];
Y=[935.411;2055.452;2331.535;1189.218];
%Vector of observations
L = [3491.901;2706.417;922.862;1819.298];


%Number of observations
no_n =length(L); 

%Initial values for the unknowns
xp = 3000;
yp = 1500;

%Vector of initial values for the unknowns
X_0 = [xp yp]';

%Number of unknowns
no_u = length(X_0);

%Redundancy
r =  no_n-no_u;

%--------------------------------------------------------------------------
%  Stochastic model
%--------------------------------------------------------------------------
s_dist = 2;  %[mm]
ppm =2;  %[parts per million

%VC Matrix of the observations
s_L = s_dist*ones(no_n,1) + ppm*L/10^3;  %[mm]
%VC Matrix of the observations
S_LL =  diag(s_L.^2);

%Theoretical standard deviation
sigma_0 =  1;    %a priori

%Cofactor matrix of the observations
Q_LL = 1/sigma_0^2 * S_LL;

%Weight matrix
P = Q_LL^(-1);

%--------------------------------------------------------------------------
%  Adjustment
%--------------------------------------------------------------------------
%break-off conditions
epsilon = 10^-8;
delta = 10^-12;
max_x_hat = 10^Inf;

%Number of iterations
iteration = 0;

while max_x_hat>epsilon || Check2>delta          
    
     %Observations as functions of the approximations for the unknowns
     L_0 = sqrt((X_0(1)-X).^2+(X_0(2)-Y).^2);

     
     %Vector of reduced observations
     l = L - L_0;
    
     %Design matrix with the elements from the Jacobian matrix J
     A = [(X_0(1)-X)./L_0 (X_0(2)-Y)./L_0];
    
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
    xp=X_0(1);
    yp=X_0(2);
     
    
     %Check 1
     max_x_hat =max(abs(x_hat)); 
     
     %Vector of residuals
     v = A*x_hat - l;
 
     %Vector of adjusted observations
     L_hat = L + v;
    
     %Objective function
     vTPv = v' * P * v;
     phi_X_hat = sqrt((X_0(1)-X).^2+(X_0(2)-Y).^2);
    


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

%Standard deviation of the residualsxp =

s_v = sqrt(diag(S_vv));

xp
yp
