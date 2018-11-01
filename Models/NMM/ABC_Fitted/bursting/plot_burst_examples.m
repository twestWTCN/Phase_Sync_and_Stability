clear; close all
load('burst_save.mat')
load('R_bursts.mat')

cmap = linspecer(2);

plot(R.IntP.tvec_obs,xsims{1}(1,1:end),'LineWidth',3,'color',cmap(1,:))
hold on

ximsbp = ft_preproc_bandpassfilter(xsims{1},2000,[14 21],[],'fir');
plot(R.IntP.tvec_obs,1.5*ximsbp(1,1:end)-5,'LineWidth',3,'color',cmap(2,:))


plot([3.5 5],[-9 -9],'k','LineWidth',5)
xlim([2 4]); ylim([-10 5])
box off

