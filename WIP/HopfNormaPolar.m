clear; close all

beta = -0.5;
omega = -1;


y = [0.1 0.5];
tend = 1000;
dt = 0.01;
t = 0;
for i = 1:tend
ydot(1) = -y(1)^5 + lambda*((1)^3) + beta*y(1);
ydot(2) = ;
y = y+(dt*y);
ystore(:,i) = y;
t = t+dt;
tvec(i) = t;
end

plot(ystore(1,:),ystore(2,:)); shg
% plot(tvec,ystore)