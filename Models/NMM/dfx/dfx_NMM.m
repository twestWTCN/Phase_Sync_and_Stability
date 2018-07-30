function dydt = dfx_NMM(y,A,K,Nr,D,par,H,T,stvec,i)
for r = 1:Nr
    for d=1:Nr
        yp(d) = y(d,2,end-D(r,d));
    end
    U = sum(K.*A(r,:).*yp);
    dydt(r,1) = (H(r)/T(r))+(sigm(U-y(r,2,i),par(:,r)))   -   ((2*y(r,1,i))/T(r))    -  (y(r,2,i)/T(r)^2); 
    dydt(r,2) = y(r,1,i);
end