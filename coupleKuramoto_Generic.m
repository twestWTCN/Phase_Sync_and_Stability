function coupleKuramoto_Generic(R,varnamez,varlist,dflag,Tag,fighan)
rng(R.seed)
rN = R.rN;
Klist = R.Klist;
% set defaults
omvar = 0.5;
sigvar =0;
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
        [d{1} d{2} d{3} d{4} d{5} d{6} SRPeps] = coupledKuramoto_wrapper_gen(R,Klist,dvar,omvar,sigvar,dflag,1)
    R.SRPeps = SRPeps;
    parfor rp = 1:rN
        [PLV(:,:,rp) dRPvar(:,:,rp) MsKappa(:,:,rp) LHat(:,:,rp) LVar(:,:,rp) RPvar(:,:,rp)] = coupledKuramoto_wrapper_gen(R,Klist,dvar,omvar,sigvar,dflag)
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
figure(fighan)
cmap = linspecer(6);
xlimmer = [-2.5 2.5];
lwid = [1 1 1];

for L = 1:3
    set(gcf,'Position',[1002         216         862         847])
    subplot(3,2,1); grid on
    [dh(1,L) da] = plotSimOut(Klist,PLV_om{L},cmap(1,:),omtrans(L),lstyle{L},lwid(L),xlimmer,L,varlist)
    xlabel('Coupling (log K)'); ylabel('PLV');
    if L==3; legend(dh(1,:),{sprintf('\\sigma %.2f',varlist(1)) ,sprintf('\\sigma %.2f',varlist(2)),sprintf('\\sigma %.2f',varlist(3))},'Location','best'); end
    
    subplot(3,2,2); grid on
    [dh(2,L) da] = plotSimOut(Klist,dRPvar_om{L},cmap(2,:),omtrans(L),lstyle{L},lwid(L),xlimmer,L,varlist)
    xlabel('Coupling (log K)'); ylabel('std d\Delta\phi/dt');
    if L==3; legend(dh(2,:),{sprintf('\\sigma %.2f',varlist(1)) ,sprintf('\\sigma %.2f',varlist(2)),sprintf('\\sigma %.2f',varlist(3))},'Location','best'); end
    
    subplot(3,2,3); grid on
    [dh(3,L) da] = plotSimOut(Klist,dRPvar_om{L},cmap(3,:),omtrans(L),lstyle{L},lwid(L),xlimmer,L,varlist)
    xlabel('Coupling (log K)'); ylabel('stable:unstable');
    if L==3; legend(dh(3,:),{sprintf('\\sigma %.2f',varlist(1)) ,sprintf('\\sigma %.2f',varlist(2)),sprintf('\\sigma %.2f',varlist(3))},'Location','best'); end
    
    subplot(3,2,5); grid on
    [dh(4,L) da] = plotSimOut(Klist,dRPvar_om{L},cmap(4,:),omtrans(L),lstyle{L},lwid(L),xlimmer,L,varlist)
    xlabel('Coupling (log K)'); ylabel('mean log Length');
    if L==3; legend(dh(4,:),{sprintf('\\sigma %.2f',varlist(1)) ,sprintf('\\sigma %.2f',varlist(2)),sprintf('\\sigma %.2f',varlist(3))},'Location','best'); end
    
    subplot(3,2,6); grid on
    [dh(5,L) da] = plotSimOut(Klist,dRPvar_om{L},cmap(5,:),omtrans(L),lstyle{L},lwid(L),xlimmer,L,varlist)
    xlabel('Coupling (log K)'); ylabel('std log Length');
    if L==3; legend(dh(5,:),{sprintf('\\sigma %.2f',varlist(1)) ,sprintf('\\sigma %.2f',varlist(2)),sprintf('\\sigma %.2f',varlist(3))},'Location','best'); end
    
    subplot(3,2,4); grid on
    [dh(6,L) da] = plotSimOut(Klist,dRPvar_om{L},cmap(6,:),omtrans(L),lstyle{L},lwid(L),xlimmer,L,varlist)
    if L==3; legend(dh(6,:),{sprintf('\\sigma %.2f',varlist(1)) ,sprintf('\\sigma %.2f',varlist(2)),sprintf('\\sigma %.2f',varlist(3))},'Location','best'); end
    xlabel('Coupling (log K)'); ylabel('std \Delta\phi');
end
savefigure_v2([R.path '\Results\Kuramoto\'],[Tag],fighan.Number,[],'-r100'); 
function [dh da] = plotSimOut(Klist,x,cmap,omtrans,lstyle,lwid,xlimmer,L,varlist)
[dh da] = boundedline(log10(Klist),mean(mean(x,3),1),mean(std(x,[],3),1),...
    'cmap',cmap,'transparency', omtrans,'alpha');
dh.LineWidth =lwid; dh.LineStyle = lstyle;
da.LineStyle = lstyle; da.EdgeColor = cmap;
xlim(xlimmer)

% subplot(3,1,3); plot(ystore(1,:),ystore(2,:)); hold on
% plot(ystore(1,end-floor(500/dt):end),ystore(2,end-floor(500/dt):end),'r');