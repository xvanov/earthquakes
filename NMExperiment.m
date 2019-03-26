% Multistroy Building
%% Clear variables and screen

clear;
clc;

%% Initial conditions

n = 1;              % Number of floors                  [-]
mass = 1;       % Mass of each floor                    [kg]
damping = 0;        % Damping constant                  [kg/s]
kons = 1;       % Building elastic constant             [kg/s^2]
f = 0;              % Magnitude of earthquake           [N]
% Time
t0 = 0;
tf = 2*pi;
h = 0.1;

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

%% NEWMARK SOL

% UNCONDITIONALLY STABLE
beta = 1.5;
gamma = 1;
[tNewmark1, dNewmark1, vNewmark1, aNewmark1] = Newmark2( M, C, K, F, t0, tf, h, d0, v0, beta, gamma);

% CONDITIONALLY STABLE
beta = 0.5;
gamma = 2;
[tNewmark2, dNewmark2, vNewmark2, aNewmark2] = Newmark2( M, C, K, F, t0, tf, h, d0, v0, beta, gamma);

% UNSTABLE
beta = 1;
gamma = 1/4;
[tNewmark3, dNewmark3, vNewmark3, aNewmark3] = Newmark2( M, C, K, F, t0, tf, h, d0, v0, beta, gamma);

% OPTIMAL
beta = 1/4;
gamma = 1/2;
[tNewmark4, dNewmark4, vNewmark4, aNewmark4] = Newmark2( M, C, K, F, t0, tf, h, d0, v0, beta, gamma);

%% RK4 SOL

[tRK4, dRK4, vRK4, aRK4] = RK4( M, C, K, F, t0, tf, h, d0, v0 );

%% ANALYT SOL

t = 0:0.1:2*pi;
y = cos(t);

%% SOL PLOT

%{
beta = 1/4;
gamma = 1/2;
h = 1;
[tNewmark05, dNewmark05, vNewmark05, aNewmark05] = Newmark2( M, C, K, F, t0, tf, h, d0, v0, beta, gamma);
h = 0.5;
[tNewmark01, dNewmark01, vNewmark01, aNewmark01] = Newmark2( M, C, K, F, t0, tf, h, d0, v0, beta, gamma);
h = 0.1;
[tNewmark001, dNewmark001, vNewmark001, aNewmark001] = Newmark2( M, C, K, F, t0, tf, h, d0, v0, beta, gamma);
hold on;
grid on;

plot(t,y, 'g');
plot(tNewmark05, dNewmark05(1,:), '--b');
plot(tNewmark01, dNewmark01(1,:), 'k');
plot(tNewmark001, dNewmark001(1,:), 'or');

title('Newmark step size comparison')
xlabel('Time')
ylabel('Solution')
legend('Analytical solution', 'Newmark \Delta t = 1', 'Newmark \Delta t = 0.5', 'Newmark \Delta t = 0.1');
%}

%{
beta = 1/4;
gamma = 1/2;
ts = 0.01:0.01:1;
tNewmark = zeros(1,629);
for i = 1:length(ts)
    h = ts(i);
    [tNewmark1, dNewmark1] = Newmark2( M, C, K, F, t0, tf, h, d0, v0, beta, gamma);
    tNewmark = [tNewmark; tNewmark1];
end
%}

t0 = 0;
tf = 5;

beta = 1/4;
gamma = 1/2;

h = 0.01;
[tNewmark1, dNewmark1] = Newmark2( M, C, K, F, t0, tf, h, d0, v0, beta, gamma);
t = t0:h:tf;
y = cos(t);
e001 = abs(y - dNewmark1);

h = 0.1;
[tNewmark2, dNewmark2] = Newmark2( M, C, K, F, t0, tf, h, d0, v0, beta, gamma);
t = t0:h:tf;
y = cos(t);
e01 = abs(y - dNewmark2);

h = 1;
[tNewmark3, dNewmark3] = Newmark2( M, C, K, F, t0, tf, h, d0, v0, beta, gamma);
t = t0:h:tf;
y = cos(t);
e1 = abs(y - dNewmark3);
y(2)

h = 0.001;
[tNewmark4, dNewmark4] = Newmark2( M, C, K, F, t0, tf, h, d0, v0, beta, gamma);
t = t0:h:tf;
y = cos(t);
e0001 = abs(y - dNewmark4);
y(2/0.001)

Ha = [0.001, 0.01, 0.1, 1];
Er = [e0001(5/0.001) , e001(5/0.01), e01(5/0.1), e1(5)];
dNewmark4(5/0.001);
dNewmark3(5)
dNewmark2(5/0.1)
dNewmark1(5/0.01)


hold on;
grid on;
plot(Ha, Er);
title('Newmark errors')
xlabel('Step Size')
ylabel('Error')

