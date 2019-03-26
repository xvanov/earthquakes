clear;
clc;

load('elcentro.mat');
hold on;
plot(eq(1,1:2501),9.8*0.1*eq(2,1:2501));
title('El Centro Earthquake Acceleration Graph');
xlabel('Time [s]');
ylabel('Acceleration [m/s^2]');
