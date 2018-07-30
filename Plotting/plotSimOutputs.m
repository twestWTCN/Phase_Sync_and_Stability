function plotSimOutputs(R,Tag,varlist,fighan)
load([R.path '\Results\Kuramoto\' Tag '_simstats'],'PLV_om','dRPvar_om','MsKappa_om','LHat_om','LVar_om','RPvar_om','rlxtime_om')
omtrans = [0.4 0.4 0.4];
lstyle = {'-','--',':'};
lwid = [1 1 1];
figure(fighan)
cmap = linspecer(6);
xlimmer = [-2.5 2.5];
lwid = [1 1 1];
Klist = R.Klist;
for L = 1:3
    set(gcf,'Position',[1002         216         862         847])
    subplot(3,2,1); grid on
    [dh(1,L) da] = plotSimOut(Klist,PLV_om{L},cmap(1,:),omtrans(L),lstyle{L},lwid(L),xlimmer,L,varlist);
    xlabel('Coupling (log K)'); ylabel('PLV');
    if L==3; legend(dh(1,:),{sprintf('\\sigma %.2f',varlist(1)) ,sprintf('\\sigma %.2f',varlist(2)),sprintf('\\sigma %.2f',varlist(3))},'Location','best'); end
    
    subplot(3,2,2); grid on
    [dh(2,L) da] = plotSimOut(Klist,dRPvar_om{L},cmap(2,:),omtrans(L),lstyle{L},lwid(L),xlimmer,L,varlist);
    xlabel('Coupling (log K)'); ylabel('std d\Delta\phi/dt');
    if L==3; legend(dh(2,:),{sprintf('\\sigma %.2f',varlist(1)) ,sprintf('\\sigma %.2f',varlist(2)),sprintf('\\sigma %.2f',varlist(3))},'Location','best'); end
    
    subplot(3,2,3); grid on
    [dh(3,L) da] = plotSimOut(Klist,MsKappa_om{L},cmap(3,:),omtrans(L),lstyle{L},lwid(L),xlimmer,L,varlist);
    xlabel('Coupling (log K)'); ylabel('stable:unstable');
    if L==3; legend(dh(3,:),{sprintf('\\sigma %.2f',varlist(1)) ,sprintf('\\sigma %.2f',varlist(2)),sprintf('\\sigma %.2f',varlist(3))},'Location','best'); end
    
    subplot(3,2,5); grid on
    [dh(4,L) da] = plotSimOut(Klist,LHat_om{L},cmap(4,:),omtrans(L),lstyle{L},lwid(L),xlimmer,L,varlist);
    xlabel('Coupling (log K)'); ylabel('mean log Length');
    if L==3; legend(dh(4,:),{sprintf('\\sigma %.2f',varlist(1)) ,sprintf('\\sigma %.2f',varlist(2)),sprintf('\\sigma %.2f',varlist(3))},'Location','best'); end
    
    subplot(3,2,6); grid on
    [dh(5,L) da] = plotSimOut(Klist,LVar_om{L},cmap(5,:),omtrans(L),lstyle{L},lwid(L),xlimmer,L,varlist);
    xlabel('Coupling (log K)'); ylabel('std log Length');
    if L==3; legend(dh(5,:),{sprintf('\\sigma %.2f',varlist(1)) ,sprintf('\\sigma %.2f',varlist(2)),sprintf('\\sigma %.2f',varlist(3))},'Location','best'); end
    
%     subplot(3,2,4); grid on
%     [dh(6,L) da] = plotSimOut(Klist,RPvar_om{L},cmap(6,:),omtrans(L),lstyle{L},lwid(L),xlimmer,L,varlist)
%     if L==3; legend(dh(6,:),{sprintf('\\sigma %.2f',varlist(1)) ,sprintf('\\sigma %.2f',varlist(2)),sprintf('\\sigma %.2f',varlist(3))},'Location','best'); end
%     xlabel('Coupling (log K)'); ylabel('std \Delta\phi');
    
    subplot(3,2,4); grid on
    [dh(6,L) da] = plotSimOut(Klist,rlxtime_om{L},cmap(6,:),omtrans(L),lstyle{L},lwid(L),xlimmer,L,varlist);
    if L==3; legend(dh(6,:),{sprintf('\\sigma %.2f',varlist(1)) ,sprintf('\\sigma %.2f',varlist(2)),sprintf('\\sigma %.2f',varlist(3))},'Location','best'); end
    xlabel('Coupling (log K)'); ylabel('rlxtime');
end
savefigure_v2([R.path '\Results\Kuramoto\'],[Tag],fighan.Number,[],'-r100'); 

function [dh da] = plotSimOut(Klist,x,cmap,omtrans,lstyle,lwid,xlimmer,L,varlist)
[dh da] = boundedline(log10(Klist),nanmean(nanmean(x,3),1),nanmean(nanstd(x,[],3),1),...
    'cmap',cmap,'transparency', omtrans,'alpha','nan','gap');
dh.LineWidth =lwid; dh.LineStyle = lstyle;
da.LineStyle = lstyle; da.EdgeColor = cmap;
xlim(xlimmer)
