function [] = add_phasesyncstab_paths(R)
addpath(R.path);
addpath(genpath([R.path '\Kuramoto_models']));
addpath(genpath([R.path '\Tools']));
addpath(genpath([R.path '\Plotting']));
addpath(genpath([R.gitpath '\Common_tools']));
addpath([R.toolpath '\SplitVec'])
addpath([R.toolpath '\Circular_Statistics_Toolbox'])
addpath([R.toolpath '\linspecer'])
addpath(genpath([R.toolpath '\boundedline-pkg']))