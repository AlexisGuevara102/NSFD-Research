% Oxy NREUP Summer Research - CUBIC DECAY DIFFERENT METHODS
% Alexis Guevara

%{
This script numerically solves the IVP u'(t) = -u(t)^3, u(0) = u_0 > 0
Domain is a<=t<b
True Solution is u(t) = u(0)/sqrt(1+2*t*(u_0)^2)
%}

format long

f1  = @(u) -u^3;    % function
a   = 0; b  = 1;    % endpoints
u_0 = 1;    % initial condition
N_t   = input('Enter the number of partitions along the t-axis: ');
delta_t = (b-a)/N_t;
% where T = delta_t * N_t;
T = zeros(1,N_t+1);
T = a:delta_t:b;

% EXACT SOLUTION
z = 1./sqrt(2*T+1);

% STANDARD FINITE METHOD
S = zeros(1,N_t+1);
S(1) = 1;
for j = 1:N_t
    S(j+1) = -(S(j)^3)*delta_t + S(j);
end
S(1) = 1;

% MICKENS METHOD
M = zeros(1,N_t+1);
M(1) = 1;
for j = 1:N_t
    M(j+1) = sqrt((M(j)^2)/(1+2*(M(j)^2)*delta_t));
end
M(1) = 1;

% E&O METHOD
E = zeros(1,N_t+1);
E(1) = 1;
for j = 1:N_t
    E(j+1) = E(j)*(delta_t*(E(j)^2)+2) / (2+3*delta_t*(E(j)^2));
end
E(1) = 1;

% GURSKI METHOD
G = zeros(1,N_t+1);
G(1) = 1;
for j = 1:N_t
    G(j+1) = (G(j)/2)*(exp(-4*G(j)*delta_t)+1);
end
G(1) = 1;

% CUBIC AVERAGE METHOD
C = zeros(1,N_t+1);
C(1) = 1;
for j = 1:N_t
    C(j+1) = (-27*(delta_t^3)*(C(j)^3) + 54*(delta_t^2)*(C(j)) + sqrt(864*(delta_t^3) + (54*(delta_t^2)*(C(j)) - 27*(delta_t^3)*(C(j)^3))^2))^(1/3)/(3*2^(1/3)*delta_t) - (2*2^(1/3))/(-27*(delta_t^3)*(C(j)^3) + 54*(delta_t^2)*(C(j)) + sqrt(864*(delta_t^3) + (54*(delta_t^2)*(C(j)) - 27*(delta_t^3)*(C(j)^3))^2))^(1/3);
end
C(1) = 1;

%{
% PLOT OF EACH EXACT SOLUTION & EACH INITIAL METHOD
plot(T,z,T,S,'k+',T,C,'ro')
%plot(T,z,T,S,'o',T,M,'+',T,E,'*',T,G.^(1/2),'x',T,C,'p')
title(['Graph of Discrete Methods u_{k+1}' ...
    'versus Exact Solution of u(t)'])
xlabel('t')
ylabel('u')
%legend('Exact Solution','SDF','Mickens','E & O','Gurski','Cubic Average')
legend('Exact Solution','SFD Scheme','Cubic Average Scheme')
%}



% PLOT ERRORS OF EACH INITIAL METHOD
error_stndrd  = abs(S - z); % STANDARD FINITE ERROR
error_mickens = abs(M - z); % MICKENS ERROR
error_eo      = abs(E - z); % E&0 ERROR
error_gurski  = abs(G.^(1/2) - z); % GURSKI ERROR
error_cubic   = abs(C - z); % CUBIC AVERAGE ERROR

plot(T,error_stndrd,'ko',T,error_cubic,'ro')
title('Graph of Errors of Discrete Methods u_{k+1}')
xlabel('t-axis')
ylabel(['Error ' char(949) '_{k}'])
legend('SFD','Cubic Average Scheme')



%{
% COMPUTE ERRORS FOR EACH INITIAL METHOD
error_inf_stndrd  = norm(S - z,Inf);
error_inf_mickens = norm(M - z,Inf);
error_inf_eo      = norm(E - z,Inf); 
error_inf_gurski  = norm(G.^(1/2) - z,Inf);

% STORE MAXIMUM ERROR
fid = fopen('CubicDecay_Diff_Methods_Errors','a');
fprintf(fid, 'n=%f             \n', delta_t); 
fprintf(fid, 'error_inf_stndrd=%14.8e   \n', error_inf_stndrd);
fprintf(fid, 'error_inf_mickens=%14.8e  \n', error_inf_mickens);
fprintf(fid, 'error_inf_eo=%14.8e       \n', error_inf_eo);
fprintf(fid, 'error_inf_gurski=%14.8e   \n', error_inf_gurski);
fclose(fid);
%}

