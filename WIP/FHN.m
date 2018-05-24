function dydt = kuramoto(t,y,a,b,tau,I)

dydt(1) = y(1) -((y(1)^3)/3) - y(2) + I;
dydt(2) = (y(1)+a-(b*y(2)));
dydt(2) = (1/tau)*dydt(2);
dydt = dydt';
