clear;
clc;
%% Conditional convergence limits vs beta
beta = linspace(1/12,1/4-0.01,100);
gamma = 1/2;

Omega = 4./((gamma+1/2)^2-4*beta);

yyaxis left;
plot(beta, Omega);
title('Conditional convergence limits')
xlabel('\beta')
ylabel('Convergence limit \Omega')



%% Periodicity Errors

yyaxis right;
perror = 0.5*(beta - 1/12)*1
plot(beta, perror);
ylabel('Error with \Omega = 1')

