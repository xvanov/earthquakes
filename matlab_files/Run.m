% Multistory Building
%% Clear variables and screen

clear;
clc;

%% Initial conditions

n = 102;              % Number of floors              [-]
mass = 3600000;       % Mass of each floor            [kg]
damping = 0;        % Damping constant              [kg/s]
kons = 200*10^6;      % Building elastic constant     [kg/s^2]
f = 0;              % Magnitude of earthquake      [N]
% Time
t0 = 0;
tf = 15;
h = 0.1;

for i = 1:n
    d0(i,1) = 0;
    v0(i,1) = 0;
    m(i,1) = mass;
    k(i,1) = kons;
    c(i,1) = damping;
end
d0(1) = 0;
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


%% Numerical sol

[tEuler, dEuler, vEuler, aEuler] = Euler( M, C, K, F, t0, tf, h, d0, v0 );
[tEulerb, dEulerb, vEulerb, aEulerb] = Eulerb( M, C, K, F, t0, tf, h, d0, v0 );
[tRK4, dRK4, vRK4, aRK4] = RK4( M, C, K, F, t0, tf, h, d0, v0 );
beta = 1/4;
gamma = 1/2;
[tNewmark1, dNewmark1, vNewmark1, aNewmark1] = Newmark1( M, C, K, F, t0, tf, h, d0, v0, beta, gamma);

%% Anatical sol

n = length(M);                 
u(:,1) = [d0;v0];

O = zeros(n);
I = eye(n);
G = -inv(M)*C;
Q = -inv(M)*K;
V = [zeros(n,1); inv(M)*F];
H = [O, I; Q, G];

[EV,E]=eig(H);
2*pi./E;

%2-floor example
A = [real(EV(:,1)),imag(EV(:,1)),real(EV(:,3)),imag(EV(:,3))];

ZZ = rref([A,u]);

t=0:0.1:15;
analyt = -ZZ(1,5)*imag(EV(1,1))*sin(imag(E(1,1))*t)+ZZ(4,5)*real(EV(1,3))*sin(imag(E(3,3))*t);
analyt2 = -ZZ(1,5)*imag(EV(2,1))*sin(imag(E(1,1))*t)+ZZ(4,5)*real(EV(2,3))*sin(imag(E(3,3))*t);
%% Plot sols
%{
hold on;
grid on;

title('1st floor motion')
xlabel('Time [s]')
ylabel('Displacement [m]')

%plot(tRK4, dRK4(1,:), 'r');
%plot(tEuler, dEuler(1,:), 'b')
plot(t,analyt, 'r');
plot(tNewmark1, dNewmark1(1,:),'k');
plot(tEulerb, dEulerb(1,:),'--b');


legend('Analytical', 'Newmark', 'Backward Euler');
%}
%% El Centro

load('elcentro.mat');
%hold on;
%plot(eq(1,1:2501),0.1*eq(2,1:2501));
a = 9.8*0.1*eq(2,1:2501)';
t0 = 0;
tf = eq(1,2501);
h = eq(1,2)-eq(1,1);
m = length(a);
F = zeros(n,m);
F(1,:) = a*mass;
beta = 1/4;
gamma = 1/2;
[tDNewmark, dDNewmark, vDNewmark, aDNewmark] = DynamicNewmark( M, C, K, F, t0, tf, h, d0, v0, beta, gamma);

%% Plot El Centro

plot(tDNewmark,dDNewmark(1,:));
title('2-Storey Family House (1st Floor, no damping)')
xlabel('Time [s]')
ylabel('Displacement [m]')

E;
2*pi./E;



