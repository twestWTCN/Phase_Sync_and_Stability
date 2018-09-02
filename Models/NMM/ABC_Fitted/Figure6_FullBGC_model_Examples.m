%% Plot Model Output
clear; close all
R = simannealsetup_NPD_Can_060718()
addpath('C:\Users\twest\Documents\Work\MATLAB ADDONS\DrosteEffect-BrewerMap-221b913')
% load('C:\Users\twest\Documents\Work\GitHub\Phase_Sync_and_Stability\Models\NMM\ABC_Fitted\Full_network\FULL_OFF_CROSS.mat')
% p = bmod;

 R.out.dag = '180718_COND1_full_run2'; % cross only
R.plot.save = 'True';
load(['C:\Users\twest\Documents\Work\GitHub\Phase_Sync_and_Stability\Models\NMM\ABC_Fitted\Full_network\' R.out.dag '\modelspec_' R.out.tag '_' R.out.dag '.mat'])
m = varo;
load(['C:\Users\twest\Documents\Work\GitHub\Phase_Sync_and_Stability\Models\NMM\ABC_Fitted\Full_network\' R.out.dag '\modelfit_' R.out.tag '_' R.out.dag '.mat'])
xobs1 = varo;
load(['C:\Users\twest\Documents\Work\GitHub\Phase_Sync_and_Stability\Models\NMM\ABC_Fitted\Full_network\' R.out.dag '\parBank_' R.out.tag '_' R.out.dag '.mat'])
parBank = varo;
p = spm_unvec(parBank(1:end-1,1),xobs1.Mfit.Pfit);
load(['C:\Users\twest\Documents\Work\GitHub\Phase_Sync_and_Stability\Models\NMM\ABC_Fitted\Full_network\' R.out.dag '\R_' R.out.tag '_' R.out.dag '.mat'])

% R.obs.gainmeth = {'difference','unitvar'};
R = setSimTime(R,256);

% Nought Connectivity
Anought = repmat(-32,3);
Acon = Anought;

% plotting ops
Lss = {'-','--'}; %,'-.'};
cmap = linspecer(6);
cmapp = cmap;
cmap(1,:) = cmapp(3,:);
cmap(2,:) = cmapp(4,:);
cmap(3,:) = cmapp(2,:);
cmap(4,:) = cmapp(1,:);

CswitchI = 0; %linspace(-4,4,64);
% CswitchJ = linspace(-3.5,2,64);
% {'MMC'  'STR'  'GPE'  'STN'  'GPI'  'THAL'}
for i = 1:length(CswitchI)
%     [i j]
    Acon = Anought;
    pnew = p;
    
%             pnew.A{1}(2,1) = -32;
%             pnew.A{1}(4,1) = -32;
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
                    plot([7.5 8],[-30 -30],'k-','LineWidth',3)
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
