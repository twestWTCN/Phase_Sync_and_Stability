clear; close all

rng(2421)
dt = 0.01;
tend = 500;
tt = tend./dt;
fsamp = 1/dt;
burn = 10*fsamp;

y = [1*pi; 2*pi; -1*pi/2; pi/2];
cfreq = 3; period = (1/3)*(1/dt);
omega = randnbetween(cfreq,0.1,4,1);

A =[0 0 0 1;
    1 0 0 0;
    0 1 0 0
    0 0 1 0];
Klist = logspace(-1.5,1,30);
Nr = size(A,2);
sppart = 5; L = 0;
for i=1:numel(Klist)
    [ystore{i} tvec{i}] = fx_Nnode_Kuramoto(dt,tt,Nr,Klist(i),A,omega,y,burn);
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
    if rem(i,sppart) == 0
        L = L+1;
        figure(1)
        subplot(size(Klist,2)./sppart,2,(2*L)-1); plot(ystore{i}(1,:),ystore{i}(2,:));
        subplot(size(Klist,2)./sppart,2,(2*L)); plot(tvec{i},ystore{i}([1 2 3],:));
    end
end
set(gcf,'Position',[757 117 880 887])

figure(2)
cmap = linspecer(6);
subplot(3,2,1); grid on
boundedline(log10(Klist),mean(PLV),std(PLV),'cmap',cmap(1,:),'transparency', 0.4); dh.LineWidth =1;
xlabel('Coupling (log K)'); ylabel('PLV');
subplot(3,2,2); grid on
boundedline(log10(Klist),mean(dRPvar),std(dRPvar),'cmap',cmap(2,:),'transparency', 0.4); dh.LineWidth =1;
xlabel('Coupling (log K)'); ylabel('std d\Delta\phi/dt')
subplot(3,2,3); grid on
boundedline(log10(Klist),mean(MsKappa),std(MsKappa),'cmap',cmap(3,:),'transparency', 0.4); dh.LineWidth =1;
xlabel('Coupling (log K)'); ylabel('stable:unstable')
subplot(3,2,4); grid on
boundedline(log10(Klist),mean(LHat),std(LHat),'cmap',cmap(4,:),'transparency', 0.4); dh.LineWidth =1;
xlabel('Coupling (log K)'); ylabel('mean log Length')
subplot(3,2,5); grid on
boundedline(log10(Klist),mean(LVar),std(LVar),'cmap',cmap(5,:),'transparency', 0.4); dh.LineWidth =1;
xlabel('Coupling (log K)'); ylabel('std log Length')
subplot(3,2,6); grid on
[dh da] = boundedline(log10(Klist),mean(RPvar),std(RPvar),'cmap',cmap(6,:),'transparency', 0.4); dh.LineWidth =1;
xlabel('Coupling (log K)'); ylabel('std \Delta\phi')

set(gcf,'Position',[680         103        1017         875])

% subplot(3,1,3); plot(ystore(1,:),ystore(2,:)); hold on
% plot(ystore(1,end-floor(500/dt):end),ystore(2,end-floor(500/dt):end),'r');