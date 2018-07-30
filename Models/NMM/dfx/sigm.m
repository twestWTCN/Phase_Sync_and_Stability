function y = sigm(x,par)
% Sigmoid function
% par(1) = Ck1; connectivity strength
% par(2) = e0;
% par(3) = r;
% par(4) = V0;
% par(5) = Ck2;
y = ((par(1)*par(2))/(1+exp(par(3)*(par(4)-(par(5)*x)))));