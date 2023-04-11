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
Y = [682.415; 203.526; 251.992; 420.028; 594.553];
X = [321.052; 310.527; 506.222; 522.646; 501.494];
x7=219.089;
y7=485.959;


%Vector of observations
L = [206.9094; 46.5027; 84.6449; 115.5251; 155.5891]*pi/200;


%Number of observations
no_n =length(L); 

%Initial values for the unknowns
x3 = 200;
y3 = 400;
w3 = 0;
%Vector of initial values for the unknowns
X_0 = [x3;y3;w3];

%Number of unknowns
no_u = length(X_0);

%Number of constraints
no_b = 1;
%Redundancy
r = no_n - no_u + no_b;

%--------------------------------------------------------------------------
%  Stochastic model
%--------------------------------------------------------------------------
s_L = 0.001*pi/200;
S_LL = diag(ones(no_n,1)*s_L^2);
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
epsilon = 10^-7;
delta = 10^-12;
max_x_hat = 10^Inf;

%Number of iterations
iteration = 0;

while max_x_hat>epsilon || Check2>delta          
    
     %Observations as functions of the approximations for the unknowns
     L_0 = atan2((Y-X_0(2)),(X-X_0(1))) - w3;

     
     %Vector of reduced observations
     l = L - L_0;
    
     %Design matrix with the elements from the Jacobian matrix J
    A = [(Y-X_0(2))./((X-X_0(1)).^2+(Y-X_0(2)).^2) (-X+X_0(1))./((X-X_0(1)).^2+(Y-X_0(2)).^2) ones(no_n,1)*(-1)];

     C = [(-x7+X_0(1))/sqrt((x7-X_0(1))^2+(y7-X_0(2))^2) (-y7+X_0(2))/sqrt((x7-X_0(1))^2+(y7-X_0(2))^2) 0];;
   
   
   %Normal matrix
    N = A'*P*A;
     
    %Extended normal matrix
    %N_ext = [N,C'; C,zeros(size(C)(1),size(C)(2)];
    N_ext = [N C'; C 0];
    
    %Vector of right hand side of normal equations
    n = A'*P*l;
     
    %Extended vector of right hand side of normal equations
     w1 = 25 - sqrt((x7-X_0(1))^2+(y7-X_0(2))^2);
    n_ext = [n;w1];
    %Inversion of normal matrix / Cofactor matrix of the unknowns
    Q_xx_ext = inv(N_ext);
    Q_xx = Q_xx_ext(1:no_u,1:no_u);  % for stdandard deviation
    
    %Solution of the normal equations
    x_hat = Q_xx_ext * n_ext;  % in fact including x_hat and k (Lagrange multiplier)
       
    %Update
    X_0 = X_0 + x_hat(1:no_u);
    x3 = X_0(1);
    y3 = X_0(2);
    w3 = X_0(3);

    %Vector of Lagrange multipliers
    k = x_hat(no_u+1:end); % end means the last element
    
    %Check 1
    max_x_hat = max(abs(x_hat(1:no_u)));
     
    %Vector of residuals
    v = A*x_hat(1:no_u)-l;
 
    %Vector of adjusted observations
    L_hat = L+v;
    
    %Objective function
    vTPv = v'*P*v;
    
    %Functional relationships without the observations
    phi_X_hat = atan2((Y-X_0(2)),(X-X_0(1))) - w3;
%     dir = atan2((Y-X_0(2)),(X-X_0(1)));
%     for i = 1:length(dir)
%         if dir(i)<0
%             dir(i)=dir(i)+2*pi;
%         end
%     end
%     phi_X_hat = dir - w3;

    %Check 2
    Check2 = max(abs(L_hat-phi_X_hat));
    
    %Update number of iterations
    iteration = iteration+1;
  
end

if Check2<=delta
    disp('Everything is fine!')
else
    disp('Something is wrong.')
end

x3
y3
w3 = (w3*200/pi) + 400

%Constraint
constraint = sqrt((x7-X_0(1))^2+(y7-X_0(2))^2)

%Empirical reference standard deviation
s_0 = sqrt(vTPv/r);

%VC matrix of adjusted unknowns
S_XX_hat = s_0^2*Q_xx;

%Standard deviation of the adjusted unknowns
s_X = sqrt(diag(S_XX_hat));        

%Cofactor matrix of adjusted observations
Q_LL_hat = A*Q_xx*A';

%VC matrix of adjusted observations
S_LL_hat = s_0^2*Q_LL_hat;

%Standard deviation of the adjusted observations
s_L_hat = sqrt(diag(S_LL_hat));     

%Cofactor matrix of the residuals
Q_vv = Q_LL-Q_LL_hat;

%VC matrix of residuals
S_vv = s_0^2*Q_vv;

%Standard deviation of the residuals
s_v = sqrt(diag(S_vv));

  