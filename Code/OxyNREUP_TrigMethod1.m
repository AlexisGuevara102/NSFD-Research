% Oxy NREUP Summer Research - Predator-Prey Model Schemes
% Alexis Guevara

%{
This script numerically solves the Trig Model u'' = -u
True Solution is z= sin(T), Z = cos(T)
%}

format long

a   = 0; b  = 2*pi;  % endpoints
u_0 = 0; u_1 = 1;    %initial value
N_t   = input('Enter the number of partitions along the t-axis: ');
delta_t = (b-a)/N_t;
T = zeros(1,N_t+1);
T = a:delta_t:b;

% VECTORS TO HOLD DIFFERENT SCHEMES
z = zeros(1, length(T));
Z = zeros(1,length(T));
u_sfd = zeros(1,length(T));
v_sfd = zeros(1, length(T));
u_mickens = zeros(1, length(T));
v_mickens = zeros(1, length(T));
u_heun = zeros(1, length(T));
v_heun = zeros(1, length(T));

% SETTING INITIAL CONDITION
u_sfd(1) = u_0;
v_sfd(1) = u_1;
u_mickens(1) = u_0;
v_mickens(1) = u_1;
u_heun(1) = u_0;
v_heun(1) = u_1;

for i = 1:length(T)-1
    u_sfd(i+1)     = v_sfd(i) * delta_t + u_sfd(i);
    v_sfd(i+1)     = v_sfd(i) - u_sfd(i) * delta_t;
    u_mickens(i+1) = v_mickens(i)*sin(delta_t) + cos(delta_t) * u_mickens(i);
    v_mickens(i+1) = v_mickens(i) * cos(delta_t) - u_mickens(i)*sin(delta_t);
    u_heun(i+1)    = u_heun(i) + (delta_t/2)*(v_heun(i) + v_heun(i) - delta_t*u_heun(i));
    v_heun(i+1)    = v_heun(i) + (delta_t/2)*(-u_heun(i) - u_heun(i) - delta_t*v_heun(i)); 
end

% TRUE SOLUTION 
z = sin(T);
Z = cos(T);

% PLOT OF TRUE SOLUTION & NUMERICAL METHODS
plot(T,z,T,Z,T,u_heun,'o-',T,v_heun,'+-');
xlabel('t');
ylabel('u');
legend('True Sine', 'True Cosine', 'Heun u', 'Heun v');

%{
% PLOT OF TRUE SOLUTION & NUMERICAL METHODS
plot(z,Z,u_heun,v_heun,'o-');
xlabel('u');
ylabel('v');
legend('True Solution', 'Heun');
%}


%{
% COMPUTE ERRORS FOR EACH INITIAL METHOD
err_inf_heun_u  = norm(z  - u_heun,Inf); % L infinity error for Mickens Mthd
err_inf_heun_v  = norm(Z  - v_heun,Inf); % L infinity error for Mickens Mthd

% STORE MAXIMUM ERROR
fid = fopen('Trig_Method1_Error','a');
fprintf(fid, 'n=%f             \n', delta_t); 
fprintf(fid, 'err_inf_heun_u=%14.8e   \n', err_inf_heun_u);
fprintf(fid, 'err_inf_heun_v=%14.8e   \n', err_inf_heun_v);
fclose(fid);
%}
