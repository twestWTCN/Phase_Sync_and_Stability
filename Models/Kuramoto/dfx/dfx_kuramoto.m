function dydt = dfx_kuramoto(y,A,K,Nr,D,omega,i)
for r = 1:Nr
    for d=1:Nr
        yp(d) = y(d,end-D(r,d));
    end
    I = (K/Nr).*sum( A(r,:).*(sin(yp-y(r,i)) + 0.25*sin(2*(yp-y(r,i)))) );
    dydt(r) = omega(r) + I;
end