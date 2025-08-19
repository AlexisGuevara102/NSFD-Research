% This program is to test out various schemes for SIR. 
% It uses the Symbolic Math Toolbox and the SIR_Euler function.
sympref('FloatingPointOutput',true); % Needed for the symbolic function

% Initial conditions
S_init = 130; % # of Suspectiable people
I_init= 50;   % # of Infected people 
R_init= 20;   % # of Recovered people
P = S_init+I_init+R_init; % Population number

% RATES & TIME 
a = 1/10;  % Recovery rate
b = 1/100; % Infection rate
h = 0.5;
N = 400; 
time = 0:h:N; 

% DEFINE NEEDED CONSTANTS 
a_bar = a/2;
b_bar = b/2;
S_star = (a/b)^2; % Threshold Value--if S goes over this then you have an epidemic
S_inf  = S_star*(1-sqrt((sqrt(S_init/S_star)-1)^2+(I_init/S_star)))^2;
t_star = (2*sqrt(I_init))/a;
t_c = (1/b_bar)*(pi-atan((sqrt(I_init))/(sqrt(S_init)-sqrt(S_star))));
r_0 = (b/a)^2*S_init;

% FUNCTION FOR S
u = @(t) (sqrt(S_init)-(a_bar/b_bar))*cos(b_bar*t)-sqrt(I_init)*sin(b_bar*t)+(a_bar/b_bar);
% FUNCTTON FOR I
v = @(t) sqrt(I_init)*cos(b_bar*t) + (sqrt(S_init) - (a_bar/b_bar))*sin(b_bar);

% TRUE SOLUTION FOR SIR (AS SYMBOLIC FUNCTIONS)
syms S(t) I(t) R(t)
S(t) = piecewise((0<t) & (t<t_c), (u(t))^2, t>=t_c, S_inf);
I    = @(t) I_init + S_init - 2 * sqrt(S_star * S_init) + 2 * sqrt(S_star * S(t)) - S(t);
R    = @(t) P - S(t) - I(t);

%{
Uses SirEuler to for approximate solution
SIR is a matrix, whose columns give values for
S, I, and R as col. 1,2 and 3, resp.
Row i represents week i-1, e.g. row 1 is week 0
and row 53 is week 52.
The inputs are initial S,I,R values and 
a (recov rate), b (infect rate), and N
%}

% PLOT TRUE SOLUTION & EULER SCHEME
SIR = SIR_Euler(S_init,I_init,R_init,a,b,N,h);
plot(time,SIR(:,1),'k+',time,S(time),'b') % S
%plot(time,SIR(:,2),'k+',time,I(time),'g') % I
%plot(time,SIR(:,3),'k+',time,R(time),'r') % R
plot(time,abs(SIR(:,1)-S(time)),'bs') % S

% this plots SIR together with true values
%plot(time,S(time),'b',time,I(time),'g',time,R(time))
%legend('S','I','R')