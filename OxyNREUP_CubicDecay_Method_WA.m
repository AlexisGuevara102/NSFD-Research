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
u_0 = 1;            % initial condition
N_t   = input('Enter the number of partitions along the t-axis: ');
delta_t = (b-a)/N_t;
%T = delta_t * N_t; % FIXME not sure about this part
T = zeros(1,N_t+1);
T = a:delta_t:b;

% TRUE SOLUTION
z = 1./sqrt(1+2*T);


% DISCRETE METHODS (WEIGHTED AVERAGES)

% FIRST METHOD (-u^3 = -1/3*(2*x^3 + y^3))
% ORDER 1
A = zeros(1,N_t+1);
A(1) = 1;
for j = 1:N_t
    A(j+1) = ((-4*(delta_t^3)*(A(j)^3) + 6*(delta_t^2)*(A(j)) + 2*sqrt((delta_t^3)*(4*(delta_t^3)*(A(j)^6) - 12*(delta_t^2)*(A(j)^4) + 9*(delta_t)*(A(j)^2) + 4)))^(2/3) - 2*2^(1/3)*(delta_t))/(2*(delta_t)*(-2*(delta_t^3)*(A(j)^3) + 3*(delta_t^2)*A(j) + sqrt((delta_t^3)*(4*(delta_t^3)*(A(j)^6) - 12*(delta_t^2)*(A(j)^4) + 9*delta_t*(A(j)^2) + 4)))^(1/3));
end
A(1) = 1;

% FIRST METHOD (-u^3 = -(2*x^4 + y^4)/(3*x))
% 
B = zeros(1,N_t+1);
B(1) = 1;
for j = 1:N_t
    B(j+1) = (1/2)*sqrt((4*2^(1/3)*(2*(delta_t^2)*(B(j)^4) - 3*delta_t*(B(j)^2)))/(delta_t*(sqrt(59049*(delta_t^2)*(B(j)^4) - 4*(24*(delta_t^2)*(B(j)^4) - 36*delta_t*(B(j)^2))^3) + 243*delta_t*(B(j)^2))^(1/3)) + (sqrt(59049*(delta_t^2)*(B(j)^4) - 4*(24*(delta_t^2)*(B(j)^4) - 36*delta_t*(B(j)^2))^3) + 243*delta_t*(B(j)^2))^(1/3)/(3*2^(1/3)*delta_t)) - (1/2)*sqrt(-(6*B(j))/(delta_t*sqrt((4*2^(1/3)*(2*(delta_t^2)*(B(j)^4) - 3*delta_t*(B(j)^2)))/(delta_t*(sqrt(59049*(delta_t^2)*(B(j)^4) - 4*(24*(delta_t^2)*(B(j)^4) - 36*delta_t*(B(j)^2))^3) + 243*delta_t*(B(j)^2))^(1/3)) + (sqrt(59049*(delta_t^2)*(B(j)^4) - 4*(24*(delta_t^2)*(B(j)^4) - 36*delta_t*(B(j)^2))^3) + 243*delta_t*(B(j)^2))^(1/3)/(3*2^(1/3)*delta_t))) - (sqrt(59049*(delta_t^2)*(B(j)^4) - 4*(24*(delta_t^2)*(B(j)^4) - 36*delta_t*(B(j)^2))^3) + 243*delta_t*(B(j)^2))^(1/3)/(3*2^(1/3)*delta_t) - (4*2^(1/3)*(2*(delta_t^2)*(B(j)^4) - 3*delta_t*(B(j)^2)))/(delta_t*(sqrt(59049*(delta_t^2)*(B(j)^4) - 4*(24*(delta_t^2)*(B(j)^4) - 36*delta_t*(B(j)^2))^3) + 243*delta_t*(B(j)^2))^(1/3)));
end
B(1) = 1;

plot(T,z,T,B)

%{
% COMPUTE ERRORS FOR EACH INITIAL METHOD
err_inf_WA_1 = norm(A - z,Inf); % L infinity error for WA Method 1


% STORE MAXIMUM ERROR
fid = fopen('Cubic_Decay_WA_Error','a');
fprintf(fid, 'n=%f             \n', delta_t); 
fprintf(fid, 'err_inf_WA_1=%14.8e   \n', err_inf_WA_1);
fclose(fid);
%}
