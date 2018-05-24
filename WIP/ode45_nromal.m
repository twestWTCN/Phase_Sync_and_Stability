tspan = [0, 20];
y0 = [2; 0];
omega = -1; beta = 1;
ode = @(t,y) hopfnormalform(t,y,omega,beta);
[t,y] = ode45(ode, tspan, y0);






% Plot of the solution
plot(t,y(:,1))
xlabel('t')
ylabel('solution y')
title('van der Pol Equation, \mu = 1')