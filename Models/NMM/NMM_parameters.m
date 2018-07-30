function [par,H,T] = NMM_parameters(Nr,EImap)
C = 135;
sigpar = 1;
% Excitatory
sigpar(1,1) = C; %Ck1 connectivity strength
sigpar(2,1) = 5; % e0
sigpar(3,1) = 0.56; % r
sigpar(4,1) = 6; %V0
sigpar(5,1) = 0.8*C; %Ck2

% Inhibitory
sigpar(1,2) = 0.25*C; %Ck1 connectivity strength
sigpar(2,2) = 5; % e0
sigpar(3,2) = 0.56; % r
sigpar(4,2) = 6; %V0
sigpar(5,2) = 0.25*C; %Ck2

rat = [3.25/10 22/10];
T = [13 25 4 22];
for i = 1:Nr
    if EImap(i) == 1
        par(:,i) = sigpar(:,1);
        T(i) = T(i);
        H(i) = T(i)*rat(1);
    elseif EImap(i) == -1
        par(:,i) = sigpar(:,2);
        T(i) = T(i);
        H(i) = T(i)*rat(2);
    end
end