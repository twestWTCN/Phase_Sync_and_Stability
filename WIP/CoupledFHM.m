clear; close all

alpha = [0.7 0.6 0.4];
beta = [0.8 0.6 0.7];
tau = [12.5 14 10.2];
Iext = 1;

y = [2  4;
     0.5 -1
     0.2 -3];
tend = 50000;
dt = 0.1;
t = 0;
A = [0 1 1;
     1 0 1;
     1 1 0];
eps = 0.02;
Nr = size(A,2);
for i = 1:tend
    for r = 1:Nr
         I = Iext + sum(eps.*A(:,r).*y(:,1));
        dydt(r,:) = FHN(t,y(r,:)',alpha(r),beta(r),tau(r),I);
    end
    y = y +(dt*dydt);
    ystore(:,i) = y(:);
    t = t+dt;
    tvec(i) = t;
end

subplot(3,1,1); plot(ystore(1,:),ystore(3,:)); shg
subplot(3,1,2); plot(tvec,ystore([1 2],:))
subplot(3,1,3); plot(ystore(1,:),ystore(2,:)); hold on
plot(ystore(1,end-floor(500/dt):end),ystore(2,end-floor(500/dt):end),'r');