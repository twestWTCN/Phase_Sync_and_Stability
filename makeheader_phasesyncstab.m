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