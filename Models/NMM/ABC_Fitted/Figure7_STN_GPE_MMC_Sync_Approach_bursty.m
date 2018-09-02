%% Plot Model Output
clear; close all
addpath('C:\Users\twest\Documents\Work\MATLAB ADDONS\DrosteEffect-BrewerMap-221b913')
load('C:\Users\twest\Documents\Work\GitHub\Phase_Sync_and_Stability\Models\NMM\ABC_Fitted\STN_GPe\bmod.mat')
load('C:\Users\twest\Documents\Work\GitHub\Phase_Sync_and_Stability\Models\NMM\ABC_Fitted\MMC_ON\MMC_ON.mat')

R.Bcond = 2;
p = bmod;
R.obs.gainmeth = {R.obs.gainmeth{1}};
R = setSimTime(R,68);
R.obs.brn = 5;
R.obs.trans.norm = 0;

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

% plotting ops
Lss = {'-','--'}; %,'-.'};
cmap = linspecer(3);

CswitchI = linspace(-2,4,64);
CswitchJ = linspace(-3,1,64);

burstFano = NaN(size(CswitchI,3),size(CswitchJ,3),2);
burstVar = NaN(size(CswitchI,3),size(CswitchJ,3),2);
burstMean = NaN(size(CswitchI,3),size(CswitchJ,3),2);

for j = 1:length(CswitchJ)
    parfor i = 1:length(CswitchI)
        [i j]
        Acon = Anought;
        pnew = p;
        
        pnew.A{1}(2,3) = CswitchI(i);
        pnew.D(3,2) = 0.8;
        pnew.A{1}(3,2) = CswitchJ(j);
        
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
            
            %             figure(100+i)
            %             for ni = 1:3
            %                 plot(linspace(0,length(xsims{1})*R.IntP.dt,length(xsims{1})),xsims{1}(ni,:)+((ni-1)*5)','color',cmap(ni,:)); hold on
            %             end
            %             xlim([5 8]); ylim([-5 12])
            %             set(gcf,'Position',[100   100   500   202])
            %             hAxes = gca;
            %             %         hAxes.XRuler.Axle.LineStyle = 'none';
            %             axis off
            
            feat_simjr{i} = feat_sim;
            beta_ind = find(R.data.feat_xscale > 14  & R.data.feat_xscale <24);
            STN(i,j) = max(squeeze(feat_sim(1,1,1,1,beta_ind)));
            GPe(i,j) = max(feat_sim(1,2,2,1,beta_ind));
            
            % Bursts
            K = 1;
            bF = NaN(1,2);
            bV = NaN(1,2);
            bM = NaN(1,2);
            for K = 1:2
                [fxsims] = ft_preproc_bandpassfilter(xsims{1}, 1/R.IntP.dt, [14 21],[], 'fir');
                fstn = abs(hilbert(fxsims(K,:)));
                burst = (fstn>2.5*std(fstn));
                burstI = find(burst);
                consecSegs = SplitVec(burstI,'consecutive');
                burstL = cellfun('length',consecSegs);
                bF(K) = var(burstL)/mean(burstL);
                bV(K) = var(burstL);
                bM(K) = mean(burstL);
            end
            
            burstFano(i,j,:) = bF;
            burstVar(i,j,:) = bV;
            burstMean(i,j,:) = bM;
        else
            STN(i,j) = NaN;
            GPe(i,j) = NaN;
            burstFano(i,j,:) = NaN(1,2);
            burstVar(i,j,:) = NaN(1,2);
            burstMean(i,j,:) = NaN(1,2);
        end
        %         figure
        %         R.plot.outFeatFx({},{feat_sim},R.data.feat_xscale,R,1,[])
    end
end

close all
[cmap] = brewermap(128,'RdYlBu');
colormap(cmap)

subplot(2,1,1)
colGridCon(CswitchI,CswitchJ,GPe,4)
xlabel('M2 \rightarrow STN','FontSize',14); ylabel('STN \rightarrow M2','FontSize',14)
title('GPe Beta Power','FontSize',14); xlim([-2 4]);ylim([-3 1]); caxis([-2 1.5])
a = gca;
a.FontSize = 16;


set(gcf,'Position',[850   437   950   425])

subplot(2,1,2)
colGridCon(CswitchI,CswitchJ,STN,4)
xlabel('M2 \rightarrow  STN'); ylabel('STN \rightarrow M2')
title('STN Beta Power'); xlim([-2 4]);ylim([-3 1]); caxis([-2 1.5])
a = gca;
a.FontSize = 16;
set(gcf,'Position',[1286         129         514         858])

load('STN_GPe_6464_burst.mat')
figure
colormap(cmap)

subplot(2,1,1)
colGridCon(CswitchI,CswitchJ,log10(squeeze(burstMean(:,:,1))),4)
xlabel('M2 \rightarrow STN','FontSize',14); ylabel('STN \rightarrow M2','FontSize',14)
title('GPe Beta Fano Factor','FontSize',14); xlim([-2 4]);ylim([-3 1]);% caxis([-2 1.5])
a = gca;
a.FontSize = 16;


set(gcf,'Position',[850   437   950   425])

subplot(2,1,2)
colGridCon(CswitchI,CswitchJ,log10(squeeze(burstVar(:,:,2))),4)
xlabel('M2 \rightarrow  STN'); ylabel('STN \rightarrow M2')
title('STN Beta Fano Factor'); xlim([-2 4]);ylim([-3 1]);% caxis([-2 1.5])
a = gca;
a.FontSize = 16;
set(gcf,'Position',[1286         129         514         858])
