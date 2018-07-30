% Master for Phase Sync and Stability Paper
close all
R = makeheader_phasesyncstab;
% add_phasesyncstab_paths(R)
R.rN = 2;
coupleKuramoto_Generic(R,'omvar',[2.5 1 0.1],0,'OmegaVar')
h1 = figure; plotSimOutputs(R,'OmegaVar',[2.5 1 0.1],h1)

coupleKuramoto_Generic(R,'omvar',[2.5 1 0.1],1,'OmegaVar_Delay')
h2 = figure; plotSimOutputs(R,'OmegaVar_Delay',[2.5 1 0.1],h2)

coupleKuramoto_Generic(R,'sigvar',[pi/128 pi/256 0],0,'SigVar')
h3 = figure; plotSimOutputs(R,'SigVar',[pi/128 pi/256 0],h3)

coupleKuramoto_Generic(R,'sigvar',[pi/128 pi/256 0],1,'SigVar_Delay')
h4 = figure; plotSimOutputs(R,'SigVar_Delay',[pi/128 pi/256 0],h4)



% To do

% 1) plot relative phase vs time and show metastable sections
% figure
% set(gcf,'Position',[1277         625         573         473])
% plot(phidata(1,:)-phidata(2,:)); grid on; 

% 2) Adjust the sync leg estimator to account for neutrally stable -
% non-sync'd solutions. Maybe connect to pertubation? Or look at 2o
% derivative?
% 3) To phase shuffle - wrap the phase signal then divide into consecutive
% segments. Resample the periods and concatanate.

R = makeheader_phasesyncstab;
coupleNMM_Generic(R,'sigvar',[22 1 0.1],1,'OmegaVar')