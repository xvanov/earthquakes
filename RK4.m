function [t, d, v, a] = RK4( M, C, K, F, t0, tf, h, d0, v0 )
%RK4 Solves the equation of motion Mx" + Cx' + Kx = F using the RK4 method.
%   Inputs
%   M           Mass matrix                     [n,n]
%   C           Damping matrix                  [n,n]
%   K           Elastic matrix                  [n,n]
%   F           Forcing vector                  [n,1]
%   t0          Starting time                   [1,1]
%   tf          Final time                      [1,1]
%   h           Time step                       [1,1]
%   d0          Initial displacements vector    [n,1]
%   v0          Initial velocity vector         [n,1]
%
%   Outputs
%   t           Discrete time interval          [1,m]
%   d           Displacement function vector    [n,m]
%   v           Velocity function vector        [n,m]
%   a           Acceleration function vector    [n,m]
%

%% Reduction of order

n = length(M);                  % Number of stories
u(:,1) = [d0;v0];

O = zeros(n);
I = eye(n);
G = -inv(M)*C;
Q = -inv(M)*K;
V = [zeros(n,1); inv(M)*F];

H = [O, I; Q, G];               % Stifness matrix

du = @(t, u) H*u + V;           % First-order differential equation

%% Preliminary calculations

w1 = 1/6;
w2 = 1/3;
w3 = 1/3;
w4 = 1/6;

t = t0:h:tf;
m = length(t);
a(:,1) = H*u(:,1) + V;

%% Time-step calculations

for i = 1:m-1
    k1 = h*feval( du, t(i),     u(:,i)      );
    k2 = h*feval( du, t(i)+h/2, u(:,i)+k1/2 );
    k3 = h*feval( du, t(i)+h/2, u(:,i)+k2/2 );
    k4 = h*feval( du, t(i)+h,   u(:,i)+k3   );
    u(:,i+1) = u(:,i) + k1*w1 + k2*w2 + k3*w3 + k4*w4;
    a(:,i+1) = feval( du, t(i+1), u(:,i+1) );
end

%% Output vector declaration

d = u(1:n,:);
v = u(n+1:2*n,:);
a = a(n+1:2*n,:);

end