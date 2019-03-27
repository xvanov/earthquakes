clear all;
clc;

theta = 0:pi/128:2*pi;
hrho2 = 1-1./(exp(i*theta));
plot(real(hrho2),imag(hrho2));
grid on;
axis([-3 3 -3 3]);
yL = ylim;
line([0 0], yL);  %y-axis
text(-1,0,'Stable','HorizontalAlignment','center');
text(1,2,'Overstable','HorizontalAlignment','center');
text(1,0,'Unstable','HorizontalAlignment','center');
title('Backward Euler''s Method Stability Region')
xlabel('Re(h\rho)')
ylabel('Im(h\rho)')