clear; close all

rng(2919)
rN = 32;
Klist = logspace(-2.5,1,75);
omvar = [0.5 0.1 0.01];
for om = 1:3
    parfor rp = 1:rN
        [PLV(:,:,rp) dRPvar(:,:,rp) MsKappa(:,:,rp) LHat(:,:,rp) LVar(:,:,rp) RPvar(:,:,rp)] = coupleKuramoto3Node_KLoop(Klist,omvar(om))
        disp(rp)
    end
    PLV_om{om} = PLV;
    dRPvar_om{om} = dRPvar;
    MsKappa_om{om} = MsKappa;
    LHat_om{om} = LHat;
    LVar_om{om} = LVar;
    RPvar_om{om} = RPvar;
end
omtrans = [0.4 0.4 0.4];
lstyle = {'-','--',':'};
lwid = [1 1 1];
    figure(1)
    cmap = linspecer(6);

for om = 1:2
    subplot(3,2,1); grid on
    [dh(1,om) da] = boundedline(log10(Klist),mean(mean(PLV_om{om},3),1),mean(std(PLV_om{om},[],3),1),...
        'cmap',cmap(1,:),'transparency', omtrans(om),'alpha');
    dh(1,om).LineWidth =lwid(om); dh(1,om).LineStyle = lstyle{om};
    da.LineStyle = lstyle{om}; da.EdgeColor = cmap(1,:);
    xlabel('Coupling (log K)'); ylabel('PLV'); xlim([-2.5 1])
    if om==3; legend(dh(1,:),{sprintf('\\sigma %.2f',omvar(1)) ,sprintf('\\sigma %.2f',omvar(2))}); end %,sprintf('\\sigma %.2f',omvar(3))},'Location','best'); end
    
    subplot(3,2,2); grid on
    [dh(2,om) da] = boundedline(log10(Klist),mean(mean(dRPvar_om{om},3),1),mean(std(dRPvar_om{om},[],3),1),...
        'cmap',cmap(2,:),'transparency', omtrans(om),'alpha');
    dh(2,om).LineWidth =lwid(om); dh(2,om).LineStyle = lstyle{om};
    da.LineStyle = lstyle{om}; da.EdgeColor = cmap(2,:);
    xlabel('Coupling (log K)'); ylabel('std d\Delta\phi/dt'); xlim([-2.5 1])
    if om==3; legend(dh(2,:),{sprintf('\\sigma %.2f',omvar(1)) ,sprintf('\\sigma %.2f',omvar(2))}); end %,sprintf('\\sigma %.2f',omvar(3))},'Location','best'); end

    subplot(3,2,3); grid on
    [dh(3,om) da] = boundedline(log10(Klist),mean(mean(MsKappa_om{om},3),1),mean(std(MsKappa_om{om},[],3),1),...
        'cmap',cmap(3,:),'transparency', omtrans(om),'alpha');
    dh(3,om).LineWidth =lwid(om); dh(3,om).LineStyle = lstyle{om};
    da.LineStyle = lstyle{om}; da.EdgeColor = cmap(3,:);
    xlabel('Coupling (log K)'); ylabel('stable:unstable'); xlim([-2.5 1])
    if om==3; legend(dh(3,:),{sprintf('\\sigma %.2f',omvar(1)) ,sprintf('\\sigma %.2f',omvar(2))}); end %,sprintf('\\sigma %.2f',omvar(3))},'Location','best'); end

    subplot(3,2,5); grid on
    [dh(4,om) da] =boundedline(log10(Klist),mean(mean(LHat_om{om},3),1),mean(std(LHat_om{om},[],3),1),...
        'cmap',cmap(4,:),'transparency', omtrans(om),'alpha');
    dh(4,om).LineWidth =lwid(om); dh(4,om).LineStyle = lstyle{om};
    da.LineStyle = lstyle{om}; da.EdgeColor = cmap(4,:);
    xlabel('Coupling (log K)'); ylabel('mean log Length'); xlim([-2.5 1])
    if om==3; legend(dh(4,:),{sprintf('\\sigma %.2f',omvar(1)) ,sprintf('\\sigma %.2f',omvar(2))}); end %,sprintf('\\sigma %.2f',omvar(3))},'Location','best'); end

    subplot(3,2,6); grid on
    [dh(5,om) da] = boundedline(log10(Klist),mean(mean(LVar_om{om},3),1),mean(std(LVar_om{om},[],3),1),...
        'cmap',cmap(5,:),'transparency', omtrans(om),'alpha');
    dh(5,om).LineWidth =lwid(om); dh(5,om).LineStyle = lstyle{om};
    da.LineStyle = lstyle{om}; da.EdgeColor = cmap(5,:);
    xlabel('Coupling (log K)'); ylabel('std log Length'); xlim([-2.5 1])
    if om==3; legend(dh(5,:),{sprintf('\\sigma %.2f',omvar(1)) ,sprintf('\\sigma %.2f',omvar(2))}); end %,sprintf('\\sigma %.2f',omvar(3))},'Location','best'); end

    subplot(3,2,4); grid on
    [dh(6,om) da] = boundedline(log10(Klist),mean(mean(RPvar_om{om},3),1),mean(std(RPvar_om{om},[],3),1),...
        'cmap',cmap(6,:),'transparency', omtrans(om),'alpha');
    dh(6,om).LineWidth =lwid(om); dh(6,om).LineStyle = lstyle{om};
    da.LineStyle = lstyle{om}; da.EdgeColor = cmap(6,:);
    xlabel('Coupling (log K)'); ylabel('std \Delta\phi'); xlim([-2.5 1])
    if om==3; legend(dh(6,:),{sprintf('\\sigma %.2f',omvar(1)) ,sprintf('\\sigma %.2f',omvar(2))}); end %,sprintf('\\sigma %.2f',omvar(3))},'Location','best'); end
end
    set(gcf,'Position',[680   103   800   875])

% subplot(3,1,3); plot(ystore(1,:),ystore(2,:)); hold on
% plot(ystore(1,end-floor(500/dt):end),ystore(2,end-floor(500/dt):end),'r');