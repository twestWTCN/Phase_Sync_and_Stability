function coupleKuramoto_Generic(R,varnamez,varlist)
rng(R.seed)
rN = R.rN;
Klist = R.Klist;

% set defaults
omvar = 0.5;
sigvar = pi/4;
dvar = 0.05;

for L = 1:3
    switch varnamez
        case 'omvar'
            omvar = varlist(L);
        case 'sigvar'
            sigvar = varlist(L);
        case 'dvar'
            dvar = varlist(L);
    end
    parfor rp = 1:rN
        [PLV(:,:,rp) dRPvar(:,:,rp) MsKappa(:,:,rp) LHat(:,:,rp) LVar(:,:,rp) RPvar(:,:,rp)] = coupledKuramoto_wrapper_gen(R,Klist,dvar,omvar,sigvar,0)
        disp(rp)
    end
    PLV_om{L} = PLV;
    dRPvar_om{L} = dRPvar;
    MsKappa_om{L} = MsKappa;
    LHat_om{L} = LHat;
    LVar_om{L} = LVar;
    RPvar_om{L} = RPvar;
end
omtrans = [0.4 0.4 0.4];
lstyle = {'-','--',':'};
lwid = [1 1 1];
figure(1)
cmap = linspecer(6);
xlimmer = [-2.5 1];
for L = 1:2
    set(gcf,'Position',[680   103   800   875])
    subplot(3,2,1); grid on
    [dh(1,L) da] = plotSimOut(Klist,PLV_om{L},cmap(1,:),omtrans{L},lstyle{L},L,xlimmer)
    xlabel('Coupling (log K)'); ylabel('PLV');
    
    subplot(3,2,2); grid on
    [dh(2,L) da] = plotSimOut(Klist,dRPvar_om{L},cmap(2,:),omtrans{L},lstyle{L},L,xlimmer)
    xlabel('Coupling (log K)'); ylabel('std d\Delta\phi/dt');
    
    subplot(3,2,3); grid on
    [dh(3,L) da] = plotSimOut(Klist,dRPvar_om{L},cmap(3,:),omtrans{L},lstyle{L},L,xlimmer)
    xlabel('Coupling (log K)'); ylabel('stable:unstable');
    
    subplot(3,2,5); grid on
    [dh(4,L) da] = plotSimOut(Klist,dRPvar_om{L},cmap(4,:),omtrans{L},lstyle{L},L,xlimmer)
    xlabel('Coupling (log K)'); ylabel('mean log Length');
    
    subplot(3,2,6); grid on
    [dh(5,L) da] = plotSimOut(Klist,dRPvar_om{L},cmap(5,:),omtrans{L},lstyle{L},L,xlimmer)
    xlabel('Coupling (log K)'); ylabel('std log Length');
    
    subplot(3,2,4); grid on
    [dh(6,L) da] = plotSimOut(Klist,dRPvar_om{L},cmap(6,:),omtrans{L},lstyle{L},L,xlimmer)
    xlabel('Coupling (log K)'); ylabel('std \Delta\phi');
end

function [dh da] = plotSimOut(Klist,x,cmap,omtrans,lstyle,L,xlimmer)
[dh da] = boundedline(log10(Klist),mean(mean(x,3),1),mean(std(x,[],3),1),...
    'cmap',cmap,'transparency', omtrans,'alpha');
dh.LineWidth =lwid(L); dh.LineStyle = lstyle{L};
da.LineStyle = lstyle{L}; da.EdgeColor = cmap;
if L==3; legend(dh(1,:),{sprintf('\\sigma %.2f',omvar(1)) ,sprintf('\\sigma %.2f',omvar(2))}); end %,sprintf('\\sigma %.2f',omvar(3))},'Location','best'); end
xlim(xlimmer)

% subplot(3,1,3); plot(ystore(1,:),ystore(2,:)); hold on
% plot(ystore(1,end-floor(500/dt):end),ystore(2,end-floor(500/dt):end),'r');