clear; close all

beta = 1;
omega = -1;

y = [4 1;
     2 -1];
tend = 10000;
dt = 0.01;
t = 0;
A = [0 1;
     1 0];
eps = 0;
Nr = size(A,1);
for i = 1:tend
    for r = 1:Nr
        z = sum(eps.*A(:,r).*y(:,1));
        dydt(r,:) = hopfnormalform(t,y(r,:)',omega,beta,z);
    end
    y = y +(dt*dydt);
    ystore(:,i) = y(:);
    t = t+dt;
    tvec(i) = t;
end

subplot(3,1,1); plot(ystore(3,:),ystore(4,:)); shg
subplot(3,1,2); plot(tvec,ystore([1 3],:))
subplot(3,1,3); plot(ystore(1,:),ystore(3,:));