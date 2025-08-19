% Oxy NREUP Summer Research - Combustion Model Schemes
% Alexis Guevara

%{
This script numerically solves the IVP u'(t) = (u^2)(1-u), u(0) = u_0 > 0
Domain is a<=t<b
True Solution is u(t) = u(0)/sqrt(1+2*t*(u_0)^2)
%}

format long

f1  = @(u) (u^2)*(1-u);    % function
a   = 0; b  = 1;           % endpoints
u_0 = 2;                   % initial condition
N_t   = input('Enter the number of partitions along the t-axis: ');
delta_t = (b-a)/N_t;
%T = delta_t * N_t; % FIXME not sure about this part
T = zeros(1,N_t+1);
T = a:delta_t:b;

%{
PROPERTIES OF COMBUSTION MODEL TO CONSIDER
INTERVAL
(-INF,0)  -> positive slope, decreasing, approaching u=0
(0,1)     -> positive slope, increasing, approaching u=1
1         -> slope 0
(1,INF)   -> negative slope, decreasing
%}

% TRUE SOLUTION (TREAT AS TRUE SOLUTION)
tspan = [T];
y0 = u_0;
combustion = @(t,y) y.^2 - y.^3;
[t,y] = ode78(combustion, tspan, y0);

plot(t,y(:,1))

% STANDARD FINITE DIFFERENCE SCHEME
S = zeros(1,N_t+1);
S(1) = u_0;
for j = 1:N_t
    S(j+1) = -delta_t*(S(j)^3) + delta_t*(S(j)^2) + S(j);
end
S(1) = u_0;

% MICKENS METHOD
M = zeros(1,N_t+1);
M(1) = u_0;
for j = 1:N_t
    M(j+1) = ((1 + 2*(1 - exp(-delta_t))*M(j))*M(j)) / (1 + (1 - exp(-delta_t))*(M(j) + (M(j)^2)));
end
M(1) = u_0;

% LUBUMA & PATIDAR METHOD
alpha_LP = 1; beta_LP  = -1;
c = -(2*beta_LP + 1) / ((eps^(-1))*(alpha_LP^2));
phi_LP = c*(1 - exp(-delta_t/c));
LP = zeros(1,N_t+1);
LP(1) = u_0;
for j = 1:N_t
    LP(j+1) = LP(j) + (LP(j)^2 - LP(j)^3)*(((eps^-1)*phi_LP) / (1 + (eps^-1)*phi_LP*(alpha_LP-1)*LP(j) + (eps^-1)*phi_LP*(1-beta_LP)*(LP(j)^2)));
end
LP(1) = u_0;

% ROEGER METHOD
alpha_R = 0;
R = zeros(1,N_t+1);
R(1) = u_0;
for j = 1:N_t
    R(j+1) = (1 + alpha_R*delta_t*R(j) + delta_t*R(j))*R(j) / (1 + alpha_R*delta_t + delta_t*(R(j)^2));
end
R(1) = u_0; 

% IBIJOLA & OBAYOMI SCHEMES
alpha_IO   = 1; beta_IO   = -1;
psi_IO = delta_t; % it is possible to change psi

% SCHEME A1 
A1 = zeros(1,N_t+1);
A1(1) = u_0;
for j = 1:N_t
    A1(j+1) = (A1(j) + psi_IO*(A1(j)^2)*(1 - alpha_IO - beta_IO*A1(j))) / (1 - alpha_IO*psi_IO*A1(j) + (1-beta_IO)*psi_IO*(A1(j)^2));
end
A1(1) = u_0;

% SCHEME A2 
A2 = zeros(1,N_t+1);
A2(1) = u_0;
for j = 1:N_t
    A2(j+1) = (A2(j) + psi_IO*(A2(j)^2)*(alpha_IO - beta_IO*A2(j))) / (1 + (alpha_IO - 1)*psi_IO*A2(j) + (1 - beta_IO)*psi_IO*(A2(j)^2));
end
A2(1) = u_0;

% SCHEME A3
A3 = zeros(1,N_t+1);
A3(1) = u_0;
for j = 1:N_t
    A3(j+1) = A3(j) + psi_IO*(A3(j)^2) - psi_IO*(A3(j)^3);
end
A3(1) = u_0;

% SCHEME A5
A5 = zeros(1,N_t+1);
A5(1) = u_0;
for j = 1:N_t
    A5(j+1) = A5(j) / (1 - psi_IO*(A5(j) - (A5(j)^2)));
end
A5(1) = u_0;


% DISCRETE METHODS (UNITY APPROXIMATIONS)
% FIRST METHOD
% cubic average scheme + square average scheme
C = zeros(1,N_t+1);
C(1) = u_0;
for j = 1:N_t
    C(j+1) = -(2^(1/3)*(6*delta_t - delta_t^2))/(3*delta_t*(-27*(delta_t^3)*(C(j)^3) + 27*(delta_t^3)*(C(j)^2) + 2*(delta_t^3) + 54*(delta_t^2)*C(j) - 18*(delta_t^2) + sqrt(4*(6*delta_t - delta_t^2)^3 + (-27*(delta_t^3)*(C(j)^3) + 27*(delta_t^3)*(C(j)^2) + 2*(delta_t^3) + 54*(delta_t^2)*C(j) - 18*(delta_t^2))^2))^(1/3)) + (-27*(delta_t^3)*(C(j)^3) + 27*(delta_t^3)*(C(j)^2) + 2*(delta_t^3) + 54*(delta_t^2)*C(j) - 18*(delta_t^2) + sqrt(4*(6*delta_t - delta_t^2)^3 + (-27*(delta_t^3)*(C(j)^3) + 27*(delta_t^3)*(C(j)^2) + 2*(delta_t^3) + 54*(delta_t^2)*C(j) - 18*(delta_t^2))^2))^(1/3)/(3*2^(1/3)*delta_t) + 1/3;
end
C(1) = u_0;


% PLOT OF TRUE SOLUTION & NUMERICAL METHODS (ALL INITIAL METHODS)

plot(t,y(:,1),'-o',T,S,'-+',T,C,'-*',T,M,'-x',T,LP,'-s',T,R,'-d',T,A1,'-^',T,A2,'-v',T,A3,'-p',T,A5,'-h')

title(['Graph of Discrete Methods u_{k+1} ' ...
    'versus Runge Kutta Method'])
xlabel('t')
ylabel('u(t)')
legend('Runge Kutta','SFD','Unity Approx.','Mickens','Lubuma & Patidar','Roeger','I&O - A1', ...
    'I&O - A2','I&O - A3','I&O - A5')


y_tpose = y';

%{
% METHOD 1 DISCRETE ERRORS FOR VALUES H = 1/N & PLOT
err_roeger_combustion = abs(R-y_tpose);
err_unity_combustion  = abs(C-y_tpose);

plot(T,err_roeger_combustion,'-+', ...
    T,err_unity_combustion,'-*')
title('Graph of Errors of Discrete Methods u_{k+1}')
xlabel('t-axis')
ylabel(['Error ' char(949) '_{k}'])
legend('Roeger1','Unity Approx.')
%}

%{
order_unity_t = [log2(0.015625) log2(0.031250) log2(0.062500) log2(0.125000) log2(0.250000)];
order_unity_u = [log2(2.98363181e-04) log2(1.19981796e-03) log2(4.90232643e-03) log2(2.13858740e-02) log2(8.62138310e-02)];

order_roeger_t = [log2(0.015625) log2(0.031250) log2(0.062500) log2(0.125000) log2(0.250000)];
order_roeger_u = [log2(1.34892022e-03) log2(2.66710368e-03) log2(5.21304095e-03) log2(1.00015951e-02) log2(1.84375314e-02)];

% PROVE ITS ORDER 2
plot(order_unity_t,order_unity_u,'o-')
title(['Order of Convergence of Cubic + Square Unity Approximation for Combustion Equation'])
xlabel(['log2(\Delta t)'])
ylabel(['log2(|| u_{num} - u_{Kutta Method} ||_{L}^{\infty})'])

% (log2(1.19981796e-03) - log2(2.98363181e-04)) / (log2(0.031250) - log2(0.015625))
% (log2(2.66710368e-03) - log2(1.34892022e-03)) / (log2(0.031250) - log2(0.015625))
%}

%{
% COMPUTE ERRORS FOR EACH INITIAL METHOD
err_inf_standard = norm(S  - y_tpose,Inf); % L infinity error for Standard Mthd
err_inf_mickens  = norm(M  - y_tpose,Inf); % L infinity error for Mickens Mthd
err_inf_lp       = norm(LP - y_tpose,Inf); % L infinity error for LP Mthd
err_inf_roeger   = norm(R  - y_tpose,Inf); % L infinity error for Roeger Mthd
err_inf_a1       = norm(A1 - y_tpose,Inf); % L infinity error for A1 Mthd
err_inf_a2       = norm(A2 - y_tpose,Inf); % L infinity error for A2 Mthd
err_inf_a3       = norm(A3 - y_tpose,Inf); % L infinity error for A3 Mthd
err_inf_a5       = norm(A5 - y_tpose,Inf); % L infinity error for A5 Mthd
err_inf_unity    = norm(C  - y_tpose,Inf); % L infinity error for Unity Approx


% STORE MAXIMUM ERROR
fid = fopen('CombustionModel_Method1_Error','a');
fprintf(fid, 'n=%f             \n', delta_t); 
fprintf(fid, 'err_inf_standard=%14.8e   \n', err_inf_standard);
fprintf(fid, 'err_inf_mickens=%14.8e   \n', err_inf_mickens);
fprintf(fid, 'err_inf_lp=%14.8e   \n', err_inf_lp);
fprintf(fid, 'err_inf_roeger=%14.8e   \n', err_inf_roeger);
fprintf(fid, 'err_inf_a1=%14.8e   \n', err_inf_a1);
fprintf(fid, 'err_inf_a2=%14.8e   \n', err_inf_a2);
fprintf(fid, 'err_inf_a3=%14.8e   \n', err_inf_a3);
fprintf(fid, 'err_inf_a5=%14.8e   \n', err_inf_a5);
fprintf(fid, 'err_inf_unity=%14.8e   \n', err_inf_unity);

fclose(fid);
%}