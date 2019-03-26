clear all;

theta = 0:pi/128:2*pi;
hrho2 = 1-1./(exp(i*theta));
plot(real(hrho2),imag(hrho2));
grid on;
axis([-3 3 -3 3]);
text(-1.0,0,'Stable','HorizontalAlignment','center');
text(-1.0,2,'Unstable','HorizontalAlignment','center');
title('Forward Euler''s Method Stability Region')
xlabel('Re(h\rho)')
ylabel('Im(h\rho)')