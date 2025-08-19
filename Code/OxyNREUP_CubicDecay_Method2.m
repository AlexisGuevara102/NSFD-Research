% Oxy NREUP Summer Research - Homework 5
% Alexis Guevara

%{
This script numerically solves the IVP u'(t) = u(t)^3, u(0) = u_0 > 0
Domain is a<=t<b
True Solution is u(t) = u(0)/sqrt(1+2*t*(u_0)^2)
%}

%BE SURE TO UNCOMMNENT WHICH (METHOD & ERROR) TO OBSERVE BEFORE RUNNING CODE%

format long

f1  = @(u) -u^3;    % function
a   = 0; b  = 1;    % endpoints
u_0 = 1;    % initial condition
N_t   = input('Enter the number of partitions along the t-axis: ');
delta_t = (b-a)/N_t;
%T = delta_t * N_t; % FIXME not sure about this part
T = zeros(1,N_t+1);
T = a:delta_t:b;

% TRUE SOLUTION
z = 1./sqrt(1+2*T);


% DISCRETE METHODS (WHEN X & Y ARE IN SECOND DEGREE)

% TWO POINT UNITY APPROXIMATIONS

% FIRST METHOD
A = zeros(1,N_t+1);
A(1) = 1;
for j = 1:N_t
    A(j+1) = (sqrt((4*delta_t*(A(j)^2)) + 1) - 1)/(2*delta_t*A(j));
end
A(1) = 1;

% SECOND METHOD
B = zeros(1,N_t+1);
B(1) = 1;
for j = 1:N_t
    B(j+1) = (1/3)*(((3*sqrt(3)*sqrt(27*(delta_t^2)*(B(j)^10) - 4*delta_t*(B(j)^8)) - 27*delta_t*(B(j)^5) + 2*(B(j)^3))^(1/3)/(2^(1/3))) + ...
        (((2^(1/3))*(B(j)^2))/((3*sqrt(3)*sqrt(27*(delta_t^2)*(B(j)^10) - 4*delta_t*(B(j)^8)) - 27*delta_t*(B(j)^5) + 2*(B(j)^3))^(1/3))) + ...
        (B(j)));
end
B(1) = 1;

% TWO POINT AVERAGING UNITY APPROXIMATIONS

% THIRD METHOD
C = zeros(1,N_t+1);
C(1) = 1;
for j = 1:N_t
    C(j+1) = (sqrt(-(delta_t^2)*(C(j)^4) + 2*delta_t*(C(j)^2) + 1) - 1)/(delta_t*C(j));
end
C(1) = 1;

% FOURTH METHOD
D = zeros(1,N_t+1);
D(1) = 1;
for j = 1:N_t
    D(j+1) = (delta_t*(D(j)^3) - 2*D(j))^2/ ... 
    (6*(-delta_t^3*D(j)^9 + 6*(delta_t^2)*(D(j)^7) + ...
    6*sqrt(3)*sqrt((delta_t^4)*(D(j)^14) - ...
    6*(delta_t^3)*(D(j)^12) + 39*(delta_t^2)*(D(j)^10) - ...
    8*delta_t*(D(j)^8)) - 66*delta_t*(D(j)^5) + 8*(D(j)^3))^(1/3)) + ...
    (1/6)*(-(delta_t^3)*(D(j)^9) + 6*(delta_t^2)*(D(j)^7) + ...
    6*sqrt(3)*sqrt((delta_t^4)*(D(j)^14) - 6*(delta_t^3)*(D(j)^12) + ...
    39*(delta_t^2)*(D(j)^10) - 8*delta_t*(D(j)^8)) - 66*delta_t*(D(j)^5) + ...
    8*(D(j)^3))^(1/3) + (1/6)*(2*D(j) - delta_t*(D(j)^3));
end
D(1) = 1;


% Plot of numerical and true solution (ALL INITIAL METHODS)
%{
plot(T,z,T,A,'o',T,B,'+',T,C,'*',T,D,'x')
title(['Graph of Discrete Square Methods y_{k+1} ' ...
    'versus exact solution y = u(0)/sqrt(1+2*t*(u_0)^2)'])
xlabel('x')
ylabel('y')
legend('Exact Solution','Discrete Method #1','Discrete Method #2','Discrete Method #3','Dicrete Method #4')
%}


% COMPUTE ERRORS FOR EACH INITIAL METHOD
err_inf_method1_2 = norm(A - z,Inf); % L infinity error for Method 1
err_inf_method2_2 = norm(B - z,Inf); % L infinity error for Method 2
err_inf_method3_2 = norm(C - z,Inf); % L infinity error for Method 3
err_inf_method4_2 = norm(D - z,Inf); % L infinity error for Method 4


% STORE MAXIMUM ERROR
fid = fopen('Homework_5_Method2_Error','a');
fprintf(fid, 'n=%f             \n', delta_t); 
fprintf(fid, 'err_inf_method1_2=%14.8e   \n', err_inf_method1_2);
fprintf(fid, 'err_inf_method2_2=%14.8e   \n', err_inf_method2_2);
fprintf(fid, 'err_inf_method3_2=%14.8e   \n', err_inf_method3_2);
fprintf(fid, 'err_inf_method4_2=%14.8e   \n', err_inf_method4_2);
fclose(fid);

