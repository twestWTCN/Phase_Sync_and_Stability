%% Plot Model Output
clear; close all
addpath('C:\Users\twest\Documents\Work\MATLAB ADDONS\DrosteEffect-BrewerMap-221b913')
load('C:\Users\twest\Documents\Work\GitHub\Phase_Sync_and_Stability\Models\NMM\ABC_Fitted\STN_GPe\bmod.mat')
load('C:\Users\twest\Documents\Work\GitHub\Phase_Sync_and_Stability\Models\NMM\ABC_Fitted\MMC_ON\MMC_ON.mat')

R.rootn = 'C:\Users\twest\Documents\Work\Github\SimAnneal_NeuroModel\Projects\Rat_NPD\';
R.out.tag = 'ModComp';
daglist = {'NPD_ModComp_M1','NPD_ModComp_M2','NPD_ModComp_M3',...
    'NPD_ModComp_M4','NPD_ModComp_M5','NPD_ModComp_M6','NPD_ModComp_M7'};
mnum = 7;
% Load Model
load([R.rootn 'outputs\' R.out.tag '\' daglist{mnum} '\modelspec_' R.out.tag '_' daglist{mnum} '.mat'])
m = varo;
% load modelfit
load([R.rootn 'outputs\' R.out.tag '\' daglist{mnum} '\modelfit_' R.out.tag '_' daglist{mnum} '.mat'])
R = varo;
BPfit =  R.Mfit.BPfit;
R.Bcond = 2;
p = bmod;
R.obs.gainmeth = R.obs.gainmeth(1);
R = setSimTime(R,256);
R.obs.brn = 5;
%% Conjoin Model
R.chsim_name = {'GPE' 'STN' 'MMC'};
R.chloc_name = {'GPE' 'STN' 'MMC'};

m.m = 3; % # of sources
m.x = {[0 0]  [0 0] [0 0 0 0 0 0 0 0]}; % Initial states
m.Gint = [1 1 14];
m.Tint = [1 1 4];
m.n =  size([m.x{:}],2); % Number of states
% These outline the models to be used in compile function
for i = 1:numel(R.chsim_name)
    m.dipfit.model(i).source = R.chsim_name{i};
end

m.outstates = {[1 0]  [1 0] [0 0 0 0 0 0 1 0]};
R.obs.outstates = find([m.outstates{:}]);
for i=1:numel(R.chloc_name)
    R.obs.obsstates(i) = find(strcmp(R.chloc_name{i},R.chsim_name));
end

% Precompute xinds to make things easier with indexing
% Compute X inds (removes need to spm_unvec which is slow)
xinds = zeros(size(m.x,2),2);
for i = 1:size(m.x,2)
    if i == 1
        xinds(i,1) = 1;
        xinds(i,2) = size(m.x{i},2);
    else
        xinds(i,1) = xinds(i-1,2)+1;
        xinds(i,2) = xinds(i,1) + (size(m.x{i},2)-1);
    end
end
m.xinds = xinds;

% setup exogenous noise
% m.uset.p = DCM.Ep;
m.uset.p.covar = eye(m.m);
m.uset.p.scale = 5e1; %.*R.InstP.dt;


p.int{3} = xcross_ON_MMC.Mfit.Pfit.int{1};
p.A{1} = [p.A{1} [-32; -32]; [-32 -32 -32]];
p.A{2} = [p.A{2} [-32; -32]; [-32 -32 -32]];
p.D = [p.D [0; 0]; [0 0 -32]];
p.C = [p.C; xcross_ON_MMC.Mfit.Pfit.C];
% Conjoin Models



% Nought Connectivity
Anought = repmat(-32,3);
Acon = Anought;
Cswitch = linspace(-6,6,32);

% plotting ops
Lss = {'-','--'}; %,'-.'};
cmap = linspecer(3);


for i = 1:4
    [i]
    Acon = Anought;
    pnew = p;
    if i == 2
        pnew.A{1}(2,3) = 1;
    elseif i == 3
        pnew.A{1}(2,3) = 2;
%         pnew.D(3,2) = 0.8;
        pnew.A{1}(3,2) = -2;
%         pnew.C(3) = p.C(3)*0.1;
    elseif i == 4
        pnew.A{1}(3,2) = 1;
    end
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
        R.obs.trans.norm = 0;
        if isfield(R.obs,'transFx')
            [~,feat_sim] = R.obs.transFx(xsims,R.chloc_name,R.chloc_name,1/R.IntP.dt,R.obs.SimOrd,R);
        else
            feat_sim = xsims; % else take raw time series
        end
        
        figure(100+i)
        for ni = 1:3
            plot(linspace(0,length(xsims{1})*R.IntP.dt,length(xsims{1})),xsims{1}(ni,:)+((ni-1)*5)','color',cmap(ni,:),'LineWidth',1.5); hold on
        end
        xlim([5 8]); ylim([-5 13.5])
        hold on
        plot([7.5 8],[-4 -4],'k','LineWidth',2)
        
        set(gcf,'Position',[100   100   500   202])
        hAxes = gca;
        %         hAxes.XRuler.Axle.LineStyle = 'none';
        axis off
        feat_simjr{i} = feat_sim;
        beta_ind = find(R.data.feat_xscale > 18  & R.data.feat_xscale <30);
%         STN(i,j) = max(squeeze(feat_sim(1,1,1,1,beta_ind)));
%         GPe(i,j) = max(feat_sim(1,2,2,1,beta_ind));
    else
%         STN(i,j) = NaN;
%         GPe(i,j) = NaN;
    end
    %         figure
    %         R.plot.outFeatFx({},{feat_sim},R.data.feat_xscale,R,1,[])
end

figure
R.plot.outFeatFx({feat_simjr{2}},{feat_simjr{1}},R.data.feat_xscale,R,1,[])
R.plot.outFeatFx({feat_simjr{4}},{feat_simjr{3}},R.data.feat_xscale,R,1,[])
for i = 1:9
    subplot(3,3,i)
    g = gca;
%     delete(g.Children(2))
g.Children(1).LineWidth = 3;
    g.Children(1).LineStyle = ':';
    
    g.Children(2).LineWidth = 3;
    g.Children(2).LineStyle = '-.';
    
    g.Children(3).LineWidth = 3;
    g.Children(3).LineStyle = '--';
    
    g.Children(4).LineWidth = 1.25;
    g.Children(4).LineStyle = '-';
    
    if any(intersect(i,[2 3 6]))
        ylim([0 1])
        g.Children(1).Color = cmap(1,:);
        g.Children(2).Color = cmap(1,:);
        g.Children(3).Color = cmap(1,:);
        g.Children(4).Color = cmap(1,:);
    elseif any(intersect(i,[4 7 8]))
        ylim([0 1])
        g.Children(1).Color = cmap(2,:);
        g.Children(2).Color = cmap(2,:);
        g.Children(3).Color = cmap(2,:);
        g.Children(4).Color = cmap(2,:);
    end
end
legend({'Model 1','Model 2','Model 3','Model 4'},'Location','SouthEast','FontSize',12)
set(gcf,'Position',[963   290   842   759])
