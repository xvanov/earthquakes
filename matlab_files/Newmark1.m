function [t, d, v, a] = Newmark1( M, C, K, F, t0, tf, h, d0, v0, beta, gamma )
%Newmark1 Solves the equation of motion Mx" + Cx' + Kx = F using the Newmark method.
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
a(:,1) = inv(M)*(-C*v0-K*d0+F);

W = inv(M + C*gamma*h + K*beta*h^2);

%% Time-step calculations

for i = 1:m-1
    D(:,i+1) = d(:,i) + v(:,i)*h + a(:,i)*(1-2*beta)*h^2/2;
    V(:,i+1) = v(:,i) + a(:,i)*(1-gamma)*h;
    a(:,i+1) = W*( F - C*V(:,i+1) - K*D(:,i+1) );
    d(:,i+1) = D(:,i+1) + a(:,i+1)*beta*h^2;
    v(:,i+1) = V(:,i+1) + a(:,i+1)*gamma*h;
end

end

