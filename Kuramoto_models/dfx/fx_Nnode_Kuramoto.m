function [ystore tvec] = fx_Nnode_Kuramoto(dt,tt,Nr,K,A,omega,y,burn)
t = 0;
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

ystore = ystore(:,burn:end);
tvec = tvec(burn:end);