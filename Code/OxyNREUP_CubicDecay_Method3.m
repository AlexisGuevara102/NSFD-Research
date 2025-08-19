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


% DISCRETE METHODS (WHEN X & Y ARE IN THIRD DEGREE)

% FIRST METHOD (ORDER 1)
A = zeros(1,N_t+1);
A(1) = 1;
for j = 1:N_t
    A(j+1) = (2^(1/3)*(sqrt(3)*sqrt((delta_t^3)*(27*delta_t*(A(j)^2) + 4)) + ...
        9*(delta_t^2)*A(j))^(2/3) - ...
        2*(3^(1/3))*delta_t)/((6^(2/3))*delta_t*(sqrt(3)*sqrt((delta_t^3)*(27*delta_t*(A(j)^2) + 4)) + ...
        9*(delta_t^2)*A(j))^(1/3));
end
A(1) = 1;

% SECOND METHOD (ORDER 2)
B = zeros(1,N_t+1);
B(1) = 1;
for j = 1:N_t
    B(j+1) = (-27*(delta_t^3)*(B(j)^3) + 54*(delta_t^2)*(B(j)) + sqrt(864*(delta_t^3) + (54*(delta_t^2)*(B(j)) - 27*(delta_t^3)*(B(j)^3))^2))^(1/3)/(3*2^(1/3)*delta_t) - (2*2^(1/3))/(-27*(delta_t^3)*(B(j)^3) + 54*(delta_t^2)*(B(j)) + sqrt(864*(delta_t^3) + (54*(delta_t^2)*(B(j)) - 27*(delta_t^3)*(B(j)^3))^2))^(1/3);
end
B(1) = 1;


% COMPUTE ERRORS FOR EACH INITIAL METHOD
err_inf_method1_3 = norm(A - z,Inf); % L infinity error for Method 1
err_inf_method2_3 = norm(B - z,Inf); % L infinity error for Method 2


% STORE MAXIMUM ERROR
fid = fopen('Homework_5_Method3_Error','a');
fprintf(fid, 'n=%f             \n', delta_t); 
fprintf(fid, 'err_inf_method1_3 = %14.8e   \n', err_inf_method1_3);
fprintf(fid, 'err_inf_method1_3 = %14.8e   \n', err_inf_method2_3);
fclose(fid);

