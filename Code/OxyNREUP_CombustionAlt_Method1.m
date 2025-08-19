% Oxy NREUP Summer Research - Combustion Alt. Method
% Alexis Guevara

%{
This script numerically solves the IVP u'(t) = -u(t)^3 + u(t), u(0) = u_0 > 0
Domain is a<=t<b
True Solution is u(t) = (((u_0) / sqrt(1-(u_0^2))).*exp(T)) ./ sqrt(1 + (((u_0) / sqrt(1-(u_0^2)))^2).*(exp(2.*T)))
%}

format long

f1  = @(u) -u^3 + u;    % function
a   = 0; b  = 1;    % endpoints
u_0 = 0.5;    % initial condition
N_t   = input('Enter the number of partitions along the t-axis: ');
delta_t = (b-a)/N_t;
% where T = delta_t * N_t;
T = zeros(1,N_t+1);
T = a:delta_t:b;

% EXACT SOLUTION
D_z = (u_0) / sqrt(1-(u_0^2));
z = (D_z.*exp(T)) ./ sqrt(1 + (D_z^2).*(exp(2.*T)));

% UNITY APPROXIMATION SCHEMES

% METHOD 1 (u_{k+1} / u_{k}) 
% -> ORDER 1
A = zeros(1,N_t+1);
A(1) = 0.5;
for j = 1:N_t
    A(j+1) = A(j) / (1 + delta_t*(A(j)^2)-delta_t);
end
A(1) = 0.5;

% METHOD 2 (u_{k} / u_{k+1}) 
% -> ORDER 1
B = zeros(1,N_t+1);
B(1) = 0.5;
for j = 1:N_t
    B(j+1) = (1/2)*(sqrt((B(j)^2)*(1 - 4*delta_t*((B(j)^2) - 1))) + B(j));
end
B(1) = 0.5;

% METHOD 3 ((u_{k+1} - 1) / (u_{k} - 1)) 
% -> ORDER 1
C = zeros(1,N_t+1);
C(1) = 0.5;
for j = 1:N_t
    C(j+1) = (C(j)*(delta_t*C(j) + delta_t + 1))/(delta_t*C(j)*(C(j) + 1) + 1);
end
C(1) = 0.5;

% METHOD 4 ((u_{k} - 1) / (u_{k+1} - 1)) 
% -> ORDER 1 (IMAGINARY NUMBERS)
D = zeros(1,N_t+1);
D(1) = 0.5;
for j = 1:N_t
    D(j+1) = (1/2)*(-sqrt(-4*delta_t*(D(j)^4) + 4*delta_t*(D(j)^3) + 4*delta_t*(D(j)^2) - 4*delta_t*D(j) + (D(j)^2) - 2*D(j) + 1) + D(j) + 1);
end
D(1) = 0.5;

% METHOD 5 ((u_{k+1} + 1) / (u_{k} + 1))
% -> ORDER 1
E = zeros(1,N_t+1);
E(1) = 0.5;
for j = 1:N_t
    E(j+1) = (E(j)*(delta_t*(-E(j)) + delta_t + 1))/(delta_t*(E(j) - 1)*E(j) + 1);
end
E(1) = 0.5;

% METHOD 6 ((u_{k} + 1) / (u_{k+1} + 1))
% -> ORDER 1
F = zeros(1,N_t+1);
F(1) = 0.5;
for j = 1:N_t
    F(j+1) = (1/2)*(sqrt(-4*delta_t*(F(j)^4) - 4*delta_t*(F(j)^3) + 4*delta_t*(F(j)^2) + 4*delta_t*F(j) + (F(j)^2) + 2*F(j) + 1) + F(j) - 1);
end
F(1) = 0.5;

% METHOD 7 ( (u_{k+1}^{2} - 1) / (u_{k}^{2} - 1))
% -> ORDER 1
G = zeros(1,N_t+1);
G(1) = 0.5;
for j = 1:N_t
    G(j+1) = (sqrt(4*(delta_t^2)*(G(j)^2) + 4*delta_t*(G(j)^2) + 1) - 1)/(2*delta_t*G(j));
end
G(1) = 0.5;

% METHOD 8 ( (u_{k}^{2} - 1) / (u_{k+1}^{2} - 1))
% H

% METHOD 9 ( (u_{k+1}^{3} - u_{k+1}) / (u_{k}^{3} - u_{k}))
% -> ORDER 1   
I = zeros(1,N_t+1);
I(1) = 0.5;
for j = 1:N_t
    I(j+1) = ((27*(delta_t^2)*I(j) + sqrt(729*(delta_t^4)*(I(j)^2) + 108*((1 -delta_t)^3)*(delta_t^3)))^(1/3)/(3*2^(1/3)*delta_t)) - (2^(1/3)*(1 - delta_t))/(27*(delta_t^2)*I(j) + sqrt(729*(delta_t^4)*(I(j)^2) + 108*((1 - delta_t)^3)*(delta_t^3)))^(1/3);
end
I(1) = 0.5;


% PLOT OF EXACT SOLUTION & EACH INITIAL METHOD
plot(T,z,T,A,'o',T,B,'+',T,C,'x',T,D,'*',T,E,':',T,F,'--',T,G,'-.',T,I,'.')
title(['Graph of Discrete Methods y_{k+1}' ...
    'versus Exact Solution u = (D_z.*exp(T)) ./ sqrt(1 + (D_z^2).*(exp(2.*T)))'])
xlabel('t')
ylabel('u')
legend('Exact Solution','A','B','C','D','E','F','G','I')


%{
% PLOT ERRORS OF EACH INITIAL METHOD
error_stndrd  = abs(S - z); % STANDARD FINITE ERROR
plot(T,error_stndrd,'o')
title('Graph of Errors of Discrete Methods y_{k+1}')
xlabel('t')
ylabel('u')
legend('Standard')
%}


% COMPUTE ERRORS FOR EACH INITIAL METHOD
err_inf_method1 = norm(A - z,Inf); % L infinity error for Method 1
err_inf_method2 = norm(B - z,Inf); % L infinity error for Method 2
err_inf_method3 = norm(C - z,Inf); % L infinity error for Method 3
err_inf_method4 = norm(D - z,Inf); % L infinity error for Method 4
err_inf_method5 = norm(E - z,Inf); % L infinity error for Method 5
err_inf_method6 = norm(F - z,Inf); % L infinity error for Method 6
err_inf_method7 = norm(G - z,Inf); % L infinity error for Method 7

err_inf_method9 = norm(I - z,Inf); % L infinity error for Method 9


% STORE MAXIMUM ERROR
fid = fopen('CombustionAlt_Method1_Errors','a');
fprintf(fid, 'n=%f             \n', delta_t); 
fprintf(fid, 'err_inf_method1=%14.8e   \n', err_inf_method1);
fprintf(fid, 'err_inf_method2=%14.8e   \n', err_inf_method2);
fprintf(fid, 'err_inf_method3=%14.8e   \n', err_inf_method3);
fprintf(fid, 'err_inf_method4=%14.8e   \n', err_inf_method4);
fprintf(fid, 'err_inf_method5=%14.8e   \n', err_inf_method5);
fprintf(fid, 'err_inf_method6=%14.8e   \n', err_inf_method6);
fprintf(fid, 'err_inf_method7=%14.8e   \n', err_inf_method7);

fprintf(fid, 'err_inf_method9=%14.8e   \n', err_inf_method9);

fclose(fid);

