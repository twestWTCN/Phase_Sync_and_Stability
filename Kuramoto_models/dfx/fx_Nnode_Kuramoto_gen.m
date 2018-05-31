function [ystore tvec] = fx_Nnode_Kuramoto_gen(dt,tt,Nr,K,A,omega,sigma,y,D,burn)
t = 0; hn = max(D);
stvec  = normrnd(0,sigma,size(A,1),tt);
for i = size(y,2):tt
    for r = 1:Nr
        for d=1:Nr
            yp(d) = y(d,end-D(r,d));
        end
         I = (K/Nr).*sum( A(r,:).*sin(yp-y(r,i)));
        dydt(r) = omega(r) + I;
    end
    y(:,i+1) = y(:,i) +(dt*dydt') + (stvec(:,i).*dt);
    ystore(:,i) = cos(y(:,i));
    t = t+dt;
    tvec(i) = t;
end

ystore = ystore(:,burn:end);
tvec = tvec(burn:end);