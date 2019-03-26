clear;
clc;

h = 0:0.001:1;

Omega = h.^2/12;
Omega2= h.^2/24;

loglog(h, Omega)
title('Newmark error \beta = 1/4, \gamma = 1/2');
xlabel('Step Size');
ylabel('Error');
grid on;
