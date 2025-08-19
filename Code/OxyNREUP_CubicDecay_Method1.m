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


% DISCRETE METHODS (WHEN X & Y ARE IN FIRST DEGREE)

% FIRST METHOD ___________________________________________________________

G = zeros(1,N_t+1);
G(1) = 1;
for j = 1:N_t
    G(j+1) = G(j)/(1+(G(j)^2)*delta_t);
end
G(1) = 1;

%{
% n = 20
n_20 = 20; 
h_20 = (b-a)/n_20;
T_20 = zeros(1,n_20+1);
T_20 = a:h_20:b;
G20 = zeros(1,n_20+1);
G20(1) = 1/sqrt(3);
for j = 1:n_20
    G20(j+1) = G20(j)/(1+(G20(j)^2)*h_20);
end
G20(1) = 1/sqrt(3);

% n = 40
n_40 = 40;
h_40 = (b-a)/n_40;
T_40 = zeros(1,n_40+1);
T_40 = a:h_40:b;
G40 = zeros(1,n_40+1);
G40(1) = 1/sqrt(3);
for j = 1:n_40
    G40(j+1) = G40(j)/(1+(G40(j)^2)*h_40);
end
G40(1) = 1/sqrt(3);

% n = 100
n_100 = 100; 
h_100 = (b-a)/100;
T_100 = zeros(1,100+1);
T_100 = a:h_100:b;
G100 = zeros(1,n_100+1);
G100(1) = 1/sqrt(3);
for j = 1:n_100
    G100(j+1) = G100(j)/(1+(G100(j)^2)*h_100);
end
G100(1) = 1/sqrt(3);

%plot(T,z,T,G,'o',T_20,G20,'*',T_40,G40,'+',T_100,G100,'x')
%}


% SECOND METHOD __________________________________________________________

I = zeros(1,N_t+1);
I(1) = 1;
for j = 1:N_t
    I(j+1) = (I(j) + sqrt((I(j)^2) - 4*delta_t*(I(j)^4)))/2;
end
I(1) = 1;

%{
% n = 20
n_20 = 20; 
h_20 = (b-a)/n_20;
T_20 = zeros(1,n_20+1);
T_20 = a:h_20:b;
I20 = zeros(1,n_20+1);
I20(1) = 1/sqrt(3);
for j = 1:n_20
    I20(j+1) = (I20(j) + sqrt((I20(j)^2) - 4*h_20*(I20(j)^4)))/2;
end
I20(1) = 1/sqrt(3);

% n = 40
n_40 = 40;
h_40 = (b-a)/n_40;
T_40 = zeros(1,n_40+1);
T_40 = a:h_40:b;
I40 = zeros(1,n_40+1);
I40(1) = 1/sqrt(3);
for j = 1:n_40
    I40(j+1) = (I40(j) + sqrt((I40(j)^2) - 4*h_40*(I40(j)^4)))/2;
end
I40(1) = 1/sqrt(3);

% n = 100
n_100 = 100; 
h_100 = (b-a)/100;
T_100 = zeros(1,100+1);
T_100 = a:h_100:b;
I100 = zeros(1,n_100+1);
I100(1) = 1/sqrt(3);
for j = 1:n_100
    I100(j+1) = (I100(j) + sqrt((I100(j)^2) - 4*h_100*(I100(j)^4)))/2;
end
I100(1) = 1/sqrt(3);

%plot(T,z,T,I,'o',T_20,I20,'*',T_40,I40,'+',T_100,I100,'x')
%}


% THIRD METHOD ___________________________________________________________

H = zeros(1,N_t+1);
H(1) = 1;
for j = 1:N_t
    H(j+1) = (2*H(j)-(delta_t*(H(j)^3)))/(2+(delta_t*(H(j)^2)));
end
H(1) = 1;

%{
% n = 20
n_20 = 20; 
h_20 = (b-a)/n_20;
T_20 = zeros(1,n_20+1);
T_20 = a:h_20:b;
H20 = zeros(1,n_20+1);
H20(1) = 1/sqrt(3);
for j = 1:n_20
    H20(j+1) = (2*H20(j)-(h_20*(H20(j)^3)))/(2+(h_20*(H20(j)^2)));
end
H20(1) = 1/sqrt(3);

% n = 40
n_40 = 40;
h_40 = (b-a)/n_40;
T_40 = zeros(1,n_40+1);
T_40 = a:h_40:b;
H40 = zeros(1,n_40+1);
H40(1) = 1/sqrt(3);
for j = 1:n_40
    H40(j+1) = (2*H40(j)-(h_40*(H40(j)^3)))/(2+(h_40*(H40(j)^2)));
end
H40(1) = 1/sqrt(3);

% n = 100
n_100 = 100; 
h_100 = (b-a)/100;
T_100 = zeros(1,100+1);
T_100 = a:h_100:b;
H100 = zeros(1,n_100+1);
H100(1) = 1/sqrt(3);
for j = 1:n_100
    H100(j+1) = (2*H100(j)-(h_100*(H100(j)^3)))/(2+(h_100*(H100(j)^2)));
end
H100(1) = 1/sqrt(3);

%plot(T,z,T,H,'o',T_20,H20,'*',T_40,H40,'+',T_100,H100,'x')
%}


% FOURTH METHOD __________________________________________________________

L = zeros(1,N_t+1);
L(1) = 1;
for j = 1:N_t
    L(j+1) = (1/4)*L(j)*(sqrt((delta_t^2)*(L(j)^4) - 12*delta_t*(L(j)^2) + 4) - delta_t*(L(j)^2) + 2);
end
L(1) = 1;
%plot(T,z,T,G,'o',T,I,'+',T,H,'x',T,L,'*')

%{
% n = 20
n_20 = 20; 
h_20 = (b-a)/n_20;
T_20 = zeros(1,n_20+1);
T_20 = a:h_20:b;
L20 = zeros(1,n_20+1);
L20(1) = 1/sqrt(3);
for j = 1:n_20
    L20(j+1) = (1/4)*L20(j)*(sqrt((h_20^2)*(L20(j)^4) - 12*h_20*(L20(j)^2) + 4) - h_20*(L20(j)^2) + 2);
end
L20(1) = 1/sqrt(3);

% n = 40
n_40 = 40;
h_40 = (b-a)/n_40;
T_40 = zeros(1,n_40+1);
T_40 = a:h_40:b;
L40 = zeros(1,n_40+1);
L40(1) = 1/sqrt(3);
for j = 1:n_40
    L40(j+1) = (1/4)*L40(j)*(sqrt((h_40^2)*(L40(j)^4) - 12*h_40*(L40(j)^2) + 4) - h_40*(L40(j)^2) + 2);
end
L40(1) = 1/sqrt(3);

% n = 100
n_100 = 100; 
h_100 = (b-a)/100;
T_100 = zeros(1,100+1);
T_100 = a:h_100:b;
L100 = zeros(1,n_100+1);
L100(1) = 1/sqrt(3);
for j = 1:n_100
    L100(j+1) = (1/4)*L100(j)*(sqrt((h_100^2)*(L100(j)^4) - 12*h_100*(L100(j)^2) + 4) - h_100*(L100(j)^2) + 2);
end
L100(1) = 1/sqrt(3);

%plot(T,z,T,H,'o',T_20,H20,'*',T_40,H40,'+',T_100,H100,'x')
%}


% Plot of numerical and true solution (ALL INITIAL METHODS)
plot(T,z,T,G,'o',T,I,'+',T,H,'*',T,L,'x')
title('Graph of Discrete Methods y_{k+1} versus exact solution y = u(0)/sqrt(1+2*t*(u_0)^2)')
xlabel('x')
ylabel('y')
legend('Exact Solution','Discrete Method #1','Discrete Method #2','Discrete Method #3','Discrete Method #4')


% COMPUTE ERRORS FOR EACH INITIAL METHOD
err_inf_method1 = norm(G - z,Inf); % L infinity error for Method 1
err_inf_method2 = norm(I - z,Inf); % L infinity error for Method 2
err_inf_method3 = norm(H - z,Inf); % L infinity error for Method 3
err_inf_method4 = norm(L - z,Inf); % L infinity error for Method 4


% METHOD 1 DISCRETE ERRORS FOR VALUES H = 1/N & PLOT
%{
err_n_method1    = abs(G-z);
err_n20_method1  = abs(G20-(1./sqrt(1+2*T_20)));
err_n40_method1  = abs(G40-(1./sqrt(1+2*T_40)));
err_n100_method1 = abs(G100-(1./sqrt(1+2*T_100)));

plot(T,err_n_method1,'o', ...
    T_20,err_n20_method1,'+', ...
    T_40,err_n40_method1,'x', ...
    T_100,err_n100_method1,'*')
title('Discrete Error for various values of H=1/N')
xlabel('x')
ylabel(['Error ' char(949) '_{k}'])
legend('N=10','N=20','N=40','N=100')
%}

% METHOD 2 DISCRETE ERRORS FOR VALUES H = 1/N & PLOT
%{
err_n_method2    = abs(I-z);
err_n20_method2  = abs(I20-(1./sqrt(1+2*T_20)));
err_n40_method2  = abs(I40-(1./sqrt(1+2*T_40)));
err_n100_method2 = abs(I100-(1./sqrt(1+2*T_100)));

plot(T,err_n_method2,'o', ...
    T_20,err_n20_method2,'+', ...
    T_40,err_n40_method2,'x', ...
    T_100,err_n100_method2,'*')
title('Discrete Error for various values of H=1/N')
xlabel('x')
ylabel(['Error ' char(949) '_{k}'])
legend('N=10','N=20','N=40','N=100')
%}


% METHOD 3 DISCRETE ERRORS FOR VALUES H = 1/N & PLOT
%{
err_n_method3    = abs(H-z);
err_n20_method3  = abs(H20-(1./sqrt(1+2*T_20)));
err_n40_method3  = abs(H40-(1./sqrt(1+2*T_40)));
err_n100_method3 = abs(H100-(1./sqrt(1+2*T_100)));

plot(T,err_n_method3,'o', ...
    T_20,err_n20_method3,'+', ...
    T_40,err_n40_method3,'x', ...
    T_100,err_n100_method3,'*')
title('Discrete Error for various values of H=1/N')
xlabel('x')
ylabel(['Error ' char(949) '_{k}'])
legend('N=10','N=20','N=40','N=100')
%}


% METHOD 4 DISCRETE ERRORS FOR VALUES H = 1/N & PLOT
%{
err_n_method4    = abs(L-z);
err_n20_method4  = abs(L20-(1./sqrt(1+2*T_20)));
err_n40_method4  = abs(L40-(1./sqrt(1+2*T_40)));
err_n100_method4 = abs(L100-(1./sqrt(1+2*T_100)));

plot(T,err_n_method4,'o', ...
    T_20,err_n20_method4,'+', ...
    T_40,err_n40_method4,'x', ...
    T_100,err_n100_method4,'*')
title('Discrete Error for various values of H=1/N')
xlabel('x')
ylabel(['Error ' char(949) '_{k}'])
legend('N=10','N=20','N=40','N=100')
%}


% STORE MAXIMUM ERROR
fid = fopen('Homework_5_Method1_Error','a');
fprintf(fid, 'n=%f             \n', delta_t); 
fprintf(fid, 'err_inf_method1=%14.8e   \n', err_inf_method1);
fprintf(fid, 'err_inf_method2=%14.8e   \n', err_inf_method2);
fprintf(fid, 'err_inf_method3=%14.8e   \n', err_inf_method3);
fprintf(fid, 'err_inf_method4=%14.8e   \n', err_inf_method4);
fclose(fid);

