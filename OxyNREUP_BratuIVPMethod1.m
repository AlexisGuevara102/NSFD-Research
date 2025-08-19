% Oxy NREUP Summer Research - Bratu IVP Method
% Alexis Guevara

%{
This script numerically solves the Bratu IVP 
True Solution is 
u(T) = -2*log(cos(T));
v(T) = 2*tan(T);
%}

format long

a   = 0; b  = 1;     % endpoints
u_0 = 0; u_1 = 0;    %initial value
N_t   = input('Enter the number of partitions along the t-axis: ');
delta_t = (b-a)/N_t;
T = zeros(1,N_t+1);
T = a:delta_t:b;

% VECTORS TO HOLD DIFFERENT SCHEMES
u_true = zeros(1,length(T));
v_true = zeros(1,length(T));
v_heun = zeros(1,length(T));
u_heun = zeros(1,length(T));

% SETTING INITIAL CONDITION
u_heun(1) = u_0;
v_heun(1) = u_1;

% SCHEMES FOR U & V
for i = 1:length(T)-1
    u_heun(i+1) = u_heun(i) + (delta_t/2)*(v_heun(i) + v_heun(i) + delta_t*2*exp(u_heun(i)));
    v_heun(i+1) = v_heun(i) + (delta_t/2)*(2*exp(u_heun(i)) + 2*exp(u_heun(i) + delta_t*v_heun(i)));
    u_w_
end

% TRUE SOLUTION
u_true = -2*log(cos(T));
v_true = 2*tan(T);

%{
% PLOT OF TRUE SOLUTION & NUMERICAL METHODS
plot(T,u_true,T,u_heun,'o-');
xlabel('t');
ylabel('u');
legend('True Solution','Heun')
%}

%{
% PLOT OF TRUE SOLUTION & NUMERICAL METHODS (U vs V)
plot(u_true,v_true,u_heun,v_heun,'o-');
xlabel('u');
ylabel('v');
legend('True Solution','Heun')
%}



% COMPUTE ERRORS FOR EACH INITIAL METHOD
%err_inf_u_heun = norm(u_true - u_heun,Inf); % L infinity error for u Heun Method
%err_inf_v_heun = norm(v_true - v_heun,Inf); % L infinity error for v Heun Method

%{
% STORE MAXIMUM ERROR
fid = fopen('BratuIVP_Method1_Error','a');
fprintf(fid, 'n=%f             \n', delta_t); 
fprintf(fid, 'err_inf_u_heun=%14.8e   \n', err_inf_u_heun);
fprintf(fid, 'err_inf_v_heun=%14.8e   \n', err_inf_v_heun);
fclose(fid);
%}
