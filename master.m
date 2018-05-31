% Master for Phase Sync and Stability Paper
close all
R = makeheader_phasesyncstab;
add_phasesyncstab_paths(R)

h1 = figure;
coupleKuramoto_Generic(R,'omvar',[5 1 0.1],0,'OmegaVar',h1)
h2 = figure;
coupleKuramoto_Generic(R,'omvar',[5 1 0.1],1,'OmegaVar_Delay',h2)
h3 = figure;
coupleKuramoto_Generic(R,'sigvar',[pi pi/4 0],0,'SigVar',h3)
h4 = figure;
coupleKuramoto_Generic(R,'sigvar',[pi pi/4 0],1,'SigVar_Delay',h4)



% To do plot relative phase vs time and show metastable sections
% figure
% set(gcf,'Position',[1277         625         573         473])
% plot(phidata(1,:)-phidata(2,:)); grid on; 