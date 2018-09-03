%% Plot Model Output
clear; close all
addpath('C:\Users\twest\Documents\Work\MATLAB ADDONS\DrosteEffect-BrewerMap-221b913')
addpath('C:\Users\twest\Documents\Work\MATLAB ADDONS\PEB_DFA')
addpath('C:\Users\twest\Documents\Work\MATLAB ADDONS\ATvDFA-package')

load('C:\Users\twest\Documents\Work\GitHub\Phase_Sync_and_Stability\Models\NMM\ABC_Fitted\STN_GPe\bmod.mat')
R.Bcond = 2;
p = bmod;
R.obs.gainmeth = R.obs.gainmeth(1);
R = setSimTime(R,512);
R.obs.brn = 5;
R.obs.trans.norm = 0;

% Nought Connectivity
Anought = repmat(-32,2);
Acon = Anought;
Cswitch = linspace(-2,5,64);

% plotting ops
Lss = {'-','--'}; %,'-.'};
cmap = linspecer(3);

DFA_Evi = zeros(2,length(Cswitch));
DFA_Alpha = zeros(2,length(Cswitch));

parfor i = 1:length(Cswitch)
    [i]
    Acon = Anought;
    pnew = p;
    
    pnew.int{1}.T = Cswitch(i);
    pnew.int{2}.T = Cswitch(i);
    
    namerz = 'none';
    rng(12312);
    % Simulate New Data
    u = innovate_timeseries(R,m);
    u{1} = u{1}.*sqrt(R.IntP.dt);
    [xsims dum wflag J] = R.IntP.intFx(R,m.x,u,pnew,m);
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
        
        %         figure(100+jr)
        %         for ni = 1:2
        %             plot(linspace(0,length(xsims{1})*R.IntP.dt,length(xsims{1})),xsims{1}(ni,:)+((ni-1)*5)','color',cmap(ni,:)); hold on
        %         end
        %         xlim([5 8]); ylim([-5 12])
        %         set(gcf,'Position',[100   100   500   202])
        %         hAxes = gca;
        %         %         hAxes.XRuler.Axle.LineStyle = 'none';
        %         axis off
        %         feat_simjr{jr} = feat_sim;
        beta_ind = find(R.data.feat_xscale > 14  & R.data.feat_xscale <30);
        maxInd = []; Frq = [];
        [STN(i) maxInd(1)] = max(squeeze(feat_sim(1,1,1,1,beta_ind)));
        Frq(1) = R.data.feat_xscale(maxInd(1));
        
        [GPe(i) maxInd(2)]= max(feat_sim(1,2,2,1,beta_ind));
        Frq(2) = R.data.feat_xscale(maxInd(2));
        
        % Compute Estimators of Stability
        alpha_tmp = []; evi_tmp =[]; var_tmp = [];
        for h = 1:2
            blims = [Frq(h)-2.5 Frq(h)+2.5];
            fs = 1/R.IntP.dt;
                        X = abs(hilbert(ft_preproc_bandpassfilter(xsims{1}(h,:),fs,blims,[],'fir')));
%             X = abs(hilbert(xsims{1}(h,:)));
            DFAP = [];
            DFAP(1) = fs; DFAP(2) = 8/blims(2);  DFAP(3) = 10;  DFAP(4) = 50;  DFAP(5) = 0;
            [bmod win evi alpha] = peb_dfa_gen(X,DFAP,0)
            
            var_tmp(h) = std(X);
            evi_tmp(h) = evi(2);
            alpha_tmp(h) = alpha(1);
        end
        Jeig = real(eig(J{1}));
%         Jeig(abs(Jeig)<1e-12) = NaN; % Remove very small (numerical error)
        sysEig(:,i) = Jeig
        DFA_Evi(:,i) = evi_tmp;
        DFA_Alpha(:,i) =alpha_tmp;
        EnvVar(:,i) = var_tmp;
    else       
        STN(i) = NaN;
        GPe(i) = NaN;
        sysEig(:,i) = NaN(4,1);
        DFA_Evi(:,i) = nan(2,1);
        DFA_Alpha(:,i) = nan(2,1);
        EnvVar(:,i) = nan(2,1);
    end
    %         figure
    %         R.plot.outFeatFx({},{feat_sim},R.data.feat_xscale,R,1,[])
end

close all
cmap = linspecer(2);
figure(1)
subplot(4,1,1)
plot(Cswitch,STN,'color',cmap(1,:),'LineWidth',2)
hold on
plot(Cswitch,GPe,'color',cmap(2,:),'LineWidth',2)
    xlabel('STN/GPe Time Constants')
    ylim([-1 10])
    xlim([-2 4])

ylabel('\beta Power')
legend({'STN','GPe'})
for h = 1:2
    subplot(4,1,h+1)
    yyaxis left
    ylabel('DFA \alpha / Env. Var.')
    a = gca;
%     a.YLabel.String = 'DFA Alpha / Envelope Variance';
    a.YColor = cmap(h,:);
    ag(1) = plot(Cswitch,DFA_Alpha(h,:),'color',cmap(h,:));
    hold on
    ag(3) = scatter(Cswitch,EnvVar(h,:),'filled','MarkerEdgeColor',cmap(h,:),'MarkerFaceColor',cmap(h,:));
    yyaxis right
    ylabel('Linearity')
    ag(2) = plot(Cswitch,DFA_Evi(h,:),'LineStyle','--','color',cmap(h,:).*0.1);
    a = gca;
%     a.YLabel.String = 'Evidence for Linearity';
    a.YColor = cmap(h,:).*0.1;
%     if h ==2
        legend(ag,{'DFA \alpha','L Score','Variance'},'box','on','Location','SouthWest')
%     end
    xlim([-2 4])
    xlabel('STN/GPe Time Constants')
end

    subplot(4,1,4)
    X = repmat(Cswitch,4,1)
scatter(X(:),sysEig(:),'kx');
hold on
plot([-2 4],[0 0 ],'k--')
ylim([-0.25 0.25])
    xlim([-2 4])

ylabel('eigenspectra')
    xlabel('Connection Strength')
    set(gcf,'Position',[680   294   538   684])