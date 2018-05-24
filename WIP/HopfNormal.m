clear; close all

beta = 1;
omega = -1;
z = 0;

y = [4; 0];
tend = 1000;
dt = 0.01;
t = 0;
for i = 1:tend
    dydt = hopfnormalform(t,y',omega,beta,z);
y = y+(dt*dydt);
ystore(:,i) = y;
t = t+dt;
tvec(i) = t;
end

plot(ystore(1,:),ystore(2,:)); shg
plot(tvec,ystore)