function [PLV dRPvar MsKappa LHat LVar RPvar] = coupledKuramoto_wrapper_gen(R,Klist,dvar,omvar,sigvar,dflag)

dt = R.dt;
tend = R.tend;
tt = tend./dt;
fsamp = 1/dt;
burn = R.burn;
cfreq = R.cfreq; 
omega = randnbetween(cfreq,omvar,4,1);

A =[0 1 1 1;
    1 0 1 1
    1 1 0 1
    1 1 1 0];
Nr = size(A,2);

if dflag == 1
    Da = floor(abs(randnbetween(0.25,dvar,sum(A(:)~=0))).*fsamp);
    D = [0      Da(1)   Da(2)   Da(3)
        Da(4)   0       Da(5)   Da(6)
        Da(7)   Da(8)   0       Da(9)
        Da(10)  Da(11)  Da(12)  0];
    D(D==0) = 1;
    y = [1*pi; 2*pi; -1*pi/2; pi/2];
    y = repmat(y,1,max(Da)+1);
else
    D = zeros(Nr);
    y = [1*pi; 2*pi; -1*pi/2; pi/2];
end

sppart = 5; L = 0;
for i=1:numel(Klist)
    [ystore{i} tvec{i}] = fx_Nnode_Kuramoto_gen(dt,tt,Nr,Klist(i),A,omega,sigvar,y,D,burn);
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
    %             subplot(size(Klist,2)./sppart,2,(2*L)); plot(tvec{i},ystore{i}([1 2 3 4],:));
    %         end
end
a = 1;