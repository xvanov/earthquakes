clear all;
clc;

gamma = 0.5:0.1:3;
beta = 1/4*(gamma+1/2).^2;

plot(gamma,beta);
grid on;
axis([0 3 0 3]);
line([0.5 0.5], [0 3]);  %gamma = 1/2
text(1,2,'Unconditionally Stable','HorizontalAlignment','center');
text(2,1,'Conditionally Stable','HorizontalAlignment','center');
text(0.25,1,'Unstable','HorizontalAlignment','center');
title('Newmark Stability Region - Undamped Case')
xlabel('\gamma')
ylabel('\beta')