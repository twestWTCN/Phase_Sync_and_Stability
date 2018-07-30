%% Plot Model Output
clear; close all
addpath('C:\Users\twest\Documents\Work\MATLAB ADDONS\DrosteEffect-BrewerMap-221b913')
load('C:\Users\twest\Documents\Work\GitHub\Phase_Sync_and_Stability\Models\NMM\ABC_Fitted\Full_network\FULL_OFF_CROSS.mat')

R.Bcond = 2;
p = bmod;
R.obs.gainmeth = {R.obs.gainmeth{1}};
R = setSimTime(R,32);

% Nought Connectivity
Anought = repmat(-32,3);
Acon = Anought;

% plotting ops
Lss = {'-','--'}; %,'-.'};
cmap = linspecer(6);

CswitchI = 0; %linspace(-4,4,64);
% CswitchJ = linspace(-3.5,2,64);
% {'MMC'  'STR'  'GPE'  'STN'  'GPI'  'THAL'}
for i = 1:length(CswitchI)
%     [i j]
    Acon = Anought;
    pnew = p;
    
            pnew.A{1}(2,1) = -32;
            pnew.A{1}(4,1) = -32;
    %         pnew.D(3,2) = 0.8;
    %         pnew.A{1}(3,2) = CswitchJ(j);
    
    namerz = 'none';
    rng(12312);
    % Simulate New Data
    u = innovate_timeseries(R,m);
    u{1} = u{1}.*sqrt(R.IntP.dt);
    [xsims dum wflag] = R.IntP.intFx(R,m.x,u,pnew,m);
    if wflag == 0
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
        
                    figure(100+i)
                    for ni = 1:size(xsims{1},1)
                        plot(linspace(0,length(xsims{1})*R.IntP.dt,length(xsims{1})),xsims{1}(ni,:)-((ni-1)*5)','color',cmap(ni,:)); hold on
                    end
                    xlim([5 8]); ylim([-35 5])
                    set(gcf,'Position',[100   100   500   202])
                    hAxes = gca;
                    %         hAxes.XRuler.Axle.LineStyle = 'none';
                    axis off
        
        feat_simjr{i} = feat_sim;
        beta_ind = find(R.data.feat_xscale > 14  & R.data.feat_xscale <24);
        %                     STN(i,j) = max(squeeze(feat_sim(1,1,1,1,beta_ind)));
        %                     GPe(i,j) = max(feat_sim(1,2,2,1,beta_ind));
    else
        %                     STN(i,j) = NaN;
        %                     GPe(i,j) = NaN;
    end
            figure
            R.plot.outFeatFx({R.data.feat_emp},{feat_sim},R.data.feat_xscale,R,1,[])
end

close all
[cmap] = brewermap(128,'RdYlBu');
colormap(cmap)

subplot(2,1,1)
colGridCon(CswitchI,CswitchJ,GPe,4)
xlabel('M2 \rightarrow STN'); ylabel('STN \rightarrow M2')
title('GPe Beta Power');
set(gcf,'Position',[850   437   950   425])

subplot(2,1,2)
colGridCon(CswitchI,CswitchJ,STN,4)
xlabel('M2 \rightarrow  STN'); ylabel('STN \rightarrow M2')
title('STN Beta Power');

set(gcf,'Position',[1286         129         514         858])
