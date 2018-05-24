clear; close all

Mu = 1;

y = [2; 0];
tend = 1e6;
dt = 0.00001;
t = 0;
for i = 1:tend
dydt = vanderpoldemo(t,y,Mu);
y = y+(dt*y);
ystore(:,i) = y;
t = t+dt;
tvec(i) = t;
end

plot(ystore(1,:),ystore(2,:)); shg
plot(tvec,ystore)