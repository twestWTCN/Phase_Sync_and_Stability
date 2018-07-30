function R = makeheader_phasesyncstab
if strcmp(getenv('COMPUTERNAME'),'SFLAP-2')
    R.path = 'C:\Users\Tim\Documents\Work\GIT\Phase_Sync_and_Stability';
    R.toolpath = 'C:\Users\Tim\Documents\MATLAB_ADDONS';
    R.gitpath = 'C:\Users\Tim\Documents\Work\GIT';
elseif strcmp(getenv('COMPUTERNAME'),'FREE')
    R.path = 'C:\Users\twest\Documents\Work\GitHub\Phase_Sync_and_Stability';
    R.toolpath = 'C:\Users\twest\Documents\Work\MATLAB ADDONS';
    R.gitpath = 'C:\Users\twest\Documents\Work\GitHub';
end
R.seed = 2312;
R.rN = 8;
R.Klist = logspace(-2,4,25); %logspace(-2.5,2.5,25);
R.dt = 0.1; %0.001;
R.tend = 1500; % 1500
R.burn = 50*(1/R.dt); %50*(1/R.dt);
R.cfreq = 15;
R.period = floor((3/R.cfreq)*(1/R.dt)); % Minimum sync size