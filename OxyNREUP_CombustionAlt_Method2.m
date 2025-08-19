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

% BERNOULLI SCHEME
% EXACT
Ber = zeros(1,N_t+1);
Ber(1) = 0.5;
for j = 1:N_t
    Ber(j+1) = (sqrt(exp(2*delta_t))*Ber(j)) / sqrt((exp(2*delta_t) - 1)*(Ber(j)^2) + 1);
end
Ber(1) = 0.5;

plot(T,z,T,Ber,'o')


% UNITY APPROXIMATION SCHEMES
%{ 
For each unity approximation, use the 'cubic average' scheme for 
-u(t)^3 and then create DIFFERENT unity approx schemes for u(t). 
Check txtbk for possible methods.
%}

% METHOD 1
% cubic average scheme + none
% ORDER 1
A = zeros(1,N_t+1);
A(1) = 0.5;
for j = 1:N_t
    A(j+1) = (-27*(delta_t^3)*(A(j)^3) + 54*(delta_t^3)*(A(j)) + 54*(delta_t^2)*(A(j)) + sqrt(864*(delta_t^3) + (-27*(delta_t^3)*(A(j)^3) + 54*(delta_t^3)*(A(j)) + 54*(delta_t^2)*(A(j)))^2))^(1/3)/(3*2^(1/3)*(delta_t)) - (2*2^(1/3))/(-27*(delta_t^3)*(A(j)^3) + 54*(delta_t^3)*(A(j)) + 54*(delta_t^2)*(A(j)) + sqrt(864*(delta_t^3) + (-27*(delta_t^3)*(A(j)^3) + 54*(delta_t^3)*(A(j)) + 54*(delta_t^2)*(A(j)))^2))^(1/3);
end
A(1) = 0.5;

% METHOD 2
% cubic average scheme + average scheme: (u_{k+1} + u_{k}) / (2u_{k})
% ORDER 2
B = zeros(1,N_t+1);
B(1) = 0.5;
for j = 1:N_t
    B(j+1) = (-27*(delta_t^3)*(B(j)^3) + 27*(delta_t^3)*(B(j)) + 54*(delta_t^2)*(B(j)) + sqrt(108*((2 - delta_t)^3)*(delta_t^3) + (-27*(delta_t^3)*(B(j)^3) + 27*(delta_t^3)*(B(j)) + 54*(delta_t^2)*(B(j)))^2))^(1/3)/(3*2^(1/3)*delta_t) - (2^(1/3)*(2 - delta_t))/(-27*(delta_t^3)*(B(j)^3) + 27*(delta_t^3)*(B(j)) + 54*(delta_t^2)*(B(j)) + sqrt(108*((2 - delta_t)^3)*(delta_t^3) + (-27*(delta_t^3)*(B(j)^3) + 27*(delta_t^3)*(B(j)) + 54*(delta_t^2)*(B(j)))^2))^(1/3);
end
B(1) = 0.5;


% PLOT OF EXACT SOLUTION & EACH INITIAL METHOD
plot(T,z,T,B,'o')
title(['Graph of Cubic Average + Average Scheme u_{k+1}' ...
    'versus Exact Solution u(t)'])
xlabel('t')
ylabel('u')
legend('Exact Solution','Cubic Average + Average Scheme')


%{
% PLOT ERRORS OF EACH INITIAL METHOD
error_stndrd  = abs(S - z); % STANDARD FINITE ERROR
plot(T,error_stndrd,'o')
title('Graph of Errors of Discrete Methods y_{k+1}')
xlabel('t')
ylabel('u')
legend('Standard')
%}


%{
% COMPUTE ERRORS FOR EACH INITIAL METHOD
err_inf_method1_2 = norm(A - z,Inf);   % L infinity error for Method 1
err_inf_method2_2 = norm(B - z,Inf);   % L infinity error for Method 2
err_inf_bernoulli = norm(Ber - z,Inf); % L infinity error for Bernoulli
err_1_bernoulli   = norm(Ber - z,1);   % L1 error for Bernoulli
err_2_bernoulli   = norm(Ber - z,2);   % L2 error for Bernoulli


% STORE MAXIMUM ERROR
fid = fopen('CombustionAlt_Method2_Errors','a');
fprintf(fid, 'n=%f             \n', delta_t); 
fprintf(fid, 'err_inf_method1=%14.8e   \n', err_inf_method1_2);
fprintf(fid, 'err_inf_method2=%14.8e   \n', err_inf_method2_2);
fprintf(fid, 'err_inf_bernoulli=%14.8e \n', err_inf_bernoulli);
fprintf(fid, 'err_1_bernoulli=%14.8e   \n', err_1_bernoulli);
fprintf(fid, 'err_2_bernoulli=%14.8e   \n', err_2_bernoulli);

fclose(fid);
%}
