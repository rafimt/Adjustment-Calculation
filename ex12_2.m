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
%error-free values
% for i = 1:size(coord,1)
%   eval(['y' num2str(coord(i,1)) '=' num2str(coord(i,2)) ';']);
%   eval(['x' num2str(coord(i,1)) '=' num2str(coord(i,3)) ';']);
% end
% X=[321.052;310.527;506.222;522.646;501.494];
% Y=[682.415;203.526;251.992;420.028;594.553];
% x1=321.052;x2=310.527;x4=506.222;x5=522.646;x6=501.494;
% y1=682.415;y2=203.526;y4=251.992;y5=420.028;y6=594.553;
x1=321.052;x2=310.527;x4=506.222;x5=522.646;x6=501.494;
y1=682.415;y2=203.526;y4=251.992;y5=420.028;y6=594.553;

%Vector of observations
L1 =206.9094;L2=46.5027;L3=84.6449;L4=115.5251;L5=155.5891; %[206.9094;46.5027;84.6449;115.5251;155.5891]*pi/200;
L_angle = [L3-L2;L4-L3;L5-L4;L1-L5]*pi/200;


%Number of observations
no_n =length(L_angle); 

%Initial values for the unknowns
x3 = 242.900;
y3 = 493.700;

%Vector of initial values for the unknowns
X_0 = [x3;y3];

%Number of unknowns
no_u = length(X_0);


%Redundancy
r = no_n - no_u ;

%--------------------------------------------------------------------------
%  Stochastic model
%--------------------------------------------------------------------------
s_dir = 10^(-3)*pi/200;


%VC Matrix of the observations
S_LL =eye(5)*s_dir^2;
F=[0 -1 1 0 0; 0 0 -1 1 0; 0 0 0 -1 1; 1 0 0 0 -1];
S_XX=F*S_LL*F'; %Theoretical standard deviation
sigma_0 =  1;    %a priori

%Cofactor matrix of the observations
Q_LL = 1/sigma_0^2 * S_XX;

%Weight matrix
p = Q_LL^(-1);

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
     L0 = [angle(x3,y3,x4,y4,x2,y2);
           angle(x3,y3,x5,y5,x4,y4);
           angle(x3,y3,x6,y6,x5,y5);
           angle(x3,y3,x1,y1,x6,y6)];
     
     %Vector of reduced observations
     l = L_angle - L0;
    
     %Design matrix with the elements from the Jacobian matrix J
     A(1,1) = alpha_dx(x3,y3,x4,y4,x2,y2);
     A(1,2) = alpha_dy(x3,y3,x4,y4,x2,y2);

     A(2,1) = alpha_dx(x3,y3,x5,y5,x4,y4);
     A(2,2) = alpha_dy(x3,y3,x5,y5,x4,y4);

     A(3,1) = alpha_dx(x3,y3,x6,y6,x5,y5);
     A(3,2) = alpha_dy(x3,y3,x6,y6,x5,y5);

     A(4,1) = alpha_dx(x3,y3,x1,y1,x6,y6);
     A(4,2) = alpha_dy(x3,y3,x1,y1,x6,y6);

    
     %Normal matrix
     N = A' * p * A;
     
     %Vector of right hand side of normal equations
     n = A' * p * l;
    
     %Inversion of normal matrix / Cofactor matrix of the unknowns
     Q_xx =N^(-1);
    
     %Solution of the normal equations
     x_hat = Q_xx * n;
       
     %Update
     X_0 =  X_0 + x_hat; 
     x3 = X_0(1);
     y3 = X_0(2);
     
     
    
     %Check 1
     max_x_hat =max(abs(x_hat)); 
     
     %Vector of residuals
     v = A*x_hat - l;
 
     %Vector of adjusted observations
     L_hat = L_angle + v; 

     %Objective function
     vTPv = v' * p * v;

     phi_X_hat = [angle(x3,y3,x4,y4,x2,y2);
           angle(x3,y3,x5,y5,x4,y4);
           angle(x3,y3,x6,y6,x5,y5);
           angle(x3,y3,x1,y1,x6,y6)];
    


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

x3
y3
