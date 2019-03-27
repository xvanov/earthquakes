function [t, d, v, a] = Newmark2( M, C, K, F, t0, tf, h, d0, v0, beta, gamma )
%Newmark2 Solves the equation of motion Mx" + Cx' + Kx = F using the Newmark method.
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

%% Preliminary calculations

t = t0:h:tf;
m = length(t);

d(:,1) = d0;
v(:,1) = v0;
a(:,1) = inv(M)*(-C*v0 - K*d0 + F);

c1 = 1/(beta*h^2);
c2 = 1/(beta*h);
c3 = 1/(2*beta) - 1;
c4 = gamma/(beta*h);
c5 = 1 - gamma/beta;
c6 = h*(1-gamma/(2*beta));


W = inv(M*c1 + C*c4 + K);

%% Time-step calculations

for i = 1:m-1
    A(:,i+1) = c1*d(:,i) + c2*v(:,i) + c3*a(:,i);
    V(:,i+1) = -c4*d(:,i) + c5*v(:,i) + c6*a(:,i);
    d(:,i+1) = W*( F + M*A(:,i+1) - C*V(:,i+1) );
    v(:,i+1) = c4*d(:,i+1) + V(:,i+1);
    a(:,i+1) = c1*d(:,i+1) - A(:,i+1);
end

end