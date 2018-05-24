clear; close all

rng(2421)
y = [1*pi; 2*pi; -1*pi/2];
omega = randnbetween(3,0.1,3,1);
dt = 0.01;
tend = 50;
tt = tend./dt;
t = 0;
A = [0 1 1;
     1 0 1;
     1 1 0];
K = 0;
Nr = size(A,2);
for i = 1:tt
    for r = 1:Nr
         I = (K/Nr).*sum(A(r,:).*sin(y-y(r))');
        dydt(r) = omega(r) + I;
    end
    y = y +(dt*dydt');
    ystore(:,i) = cos(y);
    t = t+dt;
    tvec(i) = t;
end

subplot(2,1,1); plot(ystore(1,:),ystore(2,:)); shg
subplot(2,1,2); plot(tvec,ystore([1 2 3],:))
% subplot(3,1,3); plot(ystore(1,:),ystore(2,:)); hold on
% plot(ystore(1,end-floor(500/dt):end),ystore(2,end-floor(500/dt):end),'r');