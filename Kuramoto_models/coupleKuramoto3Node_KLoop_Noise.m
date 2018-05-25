function [PLV dRPvar MsKappa LHat LVar RPvar] = coupleKuramoto3Node_KLoop_Noise(Klist,sigvar)

dt = 0.005;
tend = 500;
tt = tend./dt;
fsamp = 1/dt;
burn = 25*fsamp;

y = [1*pi; 2*pi; -1*pi/2; pi/2];
cfreq = 3; period = (1/3)*(1/dt);
omega = randnbetween(cfreq,0.5,4,1);

A =[0 1 1 1;
    1 0 1 1
    1 1 0 1
    1 1 1 0];
Nr = size(A,2);


sppart = 5; L = 0;

for i=1:numel(Klist)
    [ystore{i} tvec{i}] = fx_Nnode_Kuramoto_Noise(dt,tt,Nr,Klist(i),A,omega,sigvar,y,burn);
    a(1,:) = ystore{i}(1,:)-ystore{i}(2,:);
    a(2,:) = ystore{i}(2,:)-ystore{i}(3,:);
    a(3,:) = ystore{i}(3,:)-ystore{i}(4,:);
    a(4,:) = ystore{i}(4,:)-ystore{i}(1,:);
    a(5,:) = ystore{i}(4,:)-ystore{i}(2,:);
    a(6,:) = ystore{i}(3,:)-ystore{i}(1,:);
    for p = 1:size(a,1); PLV(p,i) = abs(mean(exp(1i*a(p,:)),2)); end
    for p = 1:size(a,1); dRPvar(p,i) = std(exp(1i*a(p,:)),[],2); end
    for p = 1:size(a,1); MsKappa(p,i) = sum(SRP_Lengths(a(p,:),diff(a(p,:)),0.005,fsamp,1))./diff(tvec{i}([1 end])); end
    for p = 1:size(a,1); LHat(p,i) = mean(log(SRP_Lengths(a(p,:),diff(a(p,:)),0.005,fsamp,1))); end
    for p = 1:size(a,1); LVar(p,i) = std(log(SRP_Lengths(a(p,:),diff(a(p,:)),0.005,fsamp,1))); end
    for p = 1:size(a,1); [dum segRP] = SRP_Lengths(a(p,:),diff(a(p,:)),0.005,fsamp,1);segRP = segRP(segRP~=0);
        RPvar(p,i) = sqrt(circ_var(segRP'));    end
%         if rem(i,sppart) == 0
%             L = L+1;
%             figure(1)
%             subplot(size(Klist,2)./sppart,2,(2*L)-1); plot(ystore{i}(1,:),ystore{i}(2,:));
%             subplot(size(Klist,2)./sppart,2,(2*L)); plot(tvec{i},ystore{i}([1 2 3],:));
%         end
end
% a =1;