%% Plot Model Output
clear; close all
load('C:\Users\twest\Documents\Work\GitHub\Phase_Sync_and_Stability\Models\NMM\ABC_Fitted\STN_GPe\bmod.mat')
R.Bcond = 2;
p = bmod;
R.obs.gainmeth = R.obs.gainmeth(1);
R = setSimTime(R,32);
% Nought Connectivity
Anought = repmat(-32,2);
Acon = Anought;
Cswitch = linspace(-2,6,64);

% plotting ops
Lss = {'-','--'}; %,'-.'};
cmap = linspecer(3);


for jr = 1:3
    for i = 1; %:length(Cswitch)
        Acon = Anought;
        pnew = p;

        if jr == 1
            pnew = pnew;
        elseif jr==2
            pnew.A{1}(1,2) = -32;
        elseif jr == 3
            pnew.A{2}(2,1) = -32;
        end
        namerz = 'none';
        rng(12312);
        % Simulate New Data
        u = innovate_timeseries(R,m);
        u{1} = u{1}.*sqrt(R.IntP.dt);
        xsims = R.IntP.intFx(R,m.x,u,pnew,m);
        % Run Observer function
        if isfield(R.obs,'obsFx')
            xsims = R.obs.obsFx(xsims,m,pnew,R);
        end
        % Run Data Transform
        if isfield(R.obs,'transFx')
            [~,feat_sim] = R.obs.transFx(xsims,R.chloc_name,R.chloc_name,1/R.IntP.dt,R.obs.SimOrd,R);
        else
            feat_sim = xsims; % else take raw time series
        end
        
        figure(100+jr)
        for ni = 1:2
            plot(linspace(0,length(xsims{1})*R.IntP.dt,length(xsims{1})),xsims{1}(ni,:)+((ni-1)*5)','color',cmap(ni,:)); hold on
        end
        xlim([5 8]); ylim([-5 12])
        set(gcf,'Position',[100   100   500   202])
        hAxes = gca;
        %         hAxes.XRuler.Axle.LineStyle = 'none';
        axis off
        feat_simjr{jr} = feat_sim;
%         figure
%         R.plot.outFeatFx({},{feat_sim},R.data.feat_xscale,R,1,[])
    end
end
figure
R.plot.outFeatFx({feat_simjr{1}},{feat_simjr{2}},R.data.feat_xscale,R,1,[])
R.plot.outFeatFx({feat_simjr{3}},{},R.data.feat_xscale,R,1,[])
for i = 1:4
    subplot(2,2,i)
    g = gca;
    delete(g.Children(2))
    g.Children(1).LineStyle = '--';
    g.Children(2).LineStyle = '-';
    g.Children(2).LineStyle = '-';
    g.Children(2).LineWidth = 3;
    g.Children(3).LineStyle = '-.';
    
    if i == 2 || i == 3
        g.Children(1).Color = cmap(i-1,:)
        g.Children(2).Color = cmap(i-1,:)
        g.Children(3).Color = cmap(i-1,:)
    end
end
legend({'STN -> GPe','STN <-> GPe','STN <- GPe'},'Location','SouthEast','FontSize',12)
set(gcf,'Position',[963   290   842   759])