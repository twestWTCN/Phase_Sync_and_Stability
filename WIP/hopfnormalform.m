function dydt = hopfnormalform(t,y,omega,beta,z)

dydt(1) = (beta*y(1)) - y(2) + (omega*y(1)*(y(1)^2 + y(2)^2 +z^2));
dydt(2) = (beta*y(2)) + y(1) + (omega*y(2)*(y(1)^2 + y(2)^2));
dydt = dydt';
