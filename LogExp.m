%% CLEAR VARIABLES AND SCREEN

clear;
clc;

%% INITIAL CODITIONS

n = 1;                      % Number of floors                  [-]
mass = 1;                   % Mass of each floor                [kg]
damping = 0;                % Damping constant                  [kg/s]
kons = 1;                   % Building elastic constant         [kg/s^2]
f = 0;                      % Magnitude of earthquake           [N]

for i = 1:n
    d0(i,1) = 0;
    v0(i,1) = 0;
    m(i,1) = mass;
    k(i,1) = kons;
    c(i,1) = damping;
end
d0(1) = 1;
v0(1) = 0;

%% Matrix declaration

M = diag(m);        % Mass matrix
C = zeros(n);       % Damping matrix
K = zeros(n);       % Elastic matrix size

% Elastic and Damping matrix definition
for i = 1:n
    for j = 1:n
        if i == n && j == n
            K(i,j) = k(i);
            C(i,j) = c(i);
        elseif i == j
            K(i,j) = (k(i)+k(i+1));
            C(i,j) = (c(i)+c(i+1));
        elseif j == i+1
            K(i,j) = -k(i+1);
            C(i,j) = -c(i+1);
        elseif j == i - 1
            K(i,j) = -k(i);
            C(i,j) = -c(i);
        end
    end
end
F = zeros(n,1);         % Forcing vector size
F(1) = f;               % Forcing vector (only ground floor)
%% TIME
t0 = 0;
tf = 10;
h0 = 0.001;
hstep = 0.001;
hf = 1;

%% ANALYT SOL

t = t0:0.001:tf;
y = cos(t); 

%% NM PARAMS

beta = 1/4;
gamma = 1/2;

%% NM SOL
h0001 = 0.001;
[tNewmark0001, dNewmark0001] = Newmark1( M, C, K, F, t0, tf, h0001, d0, v0, beta, gamma);
e0001 = y(5000) - dNewmark0001(5000);

h001 = 0.01; 
[tNewmark001, dNewmark001] = Newmark1( M, C, K, F, t0, tf, h001, d0, v0, beta, gamma);
e001 = y(5000) - dNewmark001(500);

h01 = 0.1;
[tNewmark01, dNewmark01] = Newmark1( M, C, K, F, t0, tf, h01, d0, v0, beta, gamma);
e01 = y(5000) - dNewmark01(50);

h1 = 1;
[tNewmark1, dNewmark1] = Newmark1( M, C, K, F, t0, tf, h1, d0, v0, beta, gamma);
e1 = y(5000) - dNewmark1(5);

%% PLOT
e = abs([e001, e01, e1])
h = [h001, h01, h1];




hold on;
plot(tNewmark0001, dNewmark0001, 'b');
plot(tNewmark001, dNewmark001, '--r');
plot(tNewmark01, dNewmark01, 'og');
plot(tNewmark1, dNewmark1, 'k');


%{
loglog(h, e, 'b');
%}
title('Newmark errors');
xlabel('Step Size');
ylabel('Error');
grid on;

