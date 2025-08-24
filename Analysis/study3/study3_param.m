function [v] = sleep_param(v)

v.homedir = [v.rootdir 'Analysis/study3']; v.bhvdir = [v.rootdir 'Behavioral'];
v.vardir = [v.homedir '/vars']; if ~exist(v.vardir); mkdir(v.vardir); end
cd(v.vardir);
load('colors.mat');load('bhv_vars.mat'); load('matrix_connectivity.mat');

v.ngroup = 1;  % number of experimental groups (between subject factor)
v.nroi = size(seed_info,2); % number of ROI used in analyses
v.fdr_thres = 0.05; % note that script below assumes two-tailed tests of significance

v.netwname = {'DMN','DAN','VAN','LANG','MOT','VIS'};
v.netw= {1:4, 5:8, 9:10, 11:12, 13:17, 18:21};
v.nnetw    = length(v.netwname);
v.seedname = {'lANG','rANG','PCC','MPFC','lIPS','rIPS','lFEF','rFEF','rTPJ','rIFG','lTPJ','lIFG','lSMA','lCS','rCS','lS2','rS2','lV1V2','rV1V2','lMT','rMT'};
v.bandname = {'Delta','Theta','Alpha','Sigma','Beta','Gamma'};
delta      = 1:4; theta = 4:8; alpha = 8:12; sigma = 12:16; beta = 16:30; gamma = 30:80; 
v.band     = {delta, theta, alpha, sigma, beta, gamma}; 
v.nband= numel(v.band);
v.bandfreq = {1:4, 4:8, 8:12, 12:16, 16:30, 30:80}; 
v.nf = 80;
v.stages = {'NREM2', 'NREM3', 'REM', 'RS'};
v.nseeds = length(v.seedname);
v.nnetws = length(v.netwname);
v.nbands = length(v.bandname);
v.stagediff = {'NREM2_NREM3', 'NREM2_REM' ,'NREM3_REM' ,'NREM2_RS' ,  'NREM3_RS', 'REM_RS'};
v.stagediff_title = {'NREM2 - NREM3', 'NREM2 - REM',  'NREM3 - REM' ,'NREM2 - RS' ,'NREM3 - RS' ,'REM - RS' };
v.stagecorr = {'NREM2xNREM3', 'NREM2xREM'  ,'NREM3xREM', 'NREM2xRS', 'NREM3xRS' ,'REMxRS'};
v.sess = {'CTRL', 'LRN'};
v.nstages = length(v.stages);
v.nsess = length(v.sess); % number of sessions / RS scans per subject
v.nstagecorr = length(v.stagecorr);
v.nstagediff = length(v.stagediff);
v.ncycles = 1; 
v.levels = v.nsess * v.nstages * v.ncycles; % total number of groups and sessions

v.bhv_vars = {'offline'};
v.nbhv_vars = 1;

v.colors = colors; 
v.colors_onesided = colors_onesided;
reds = 0:0.0005:1; v.j_colors = [1 1 1; ones(length(reds),1), flipud(reds'), flipud(reds'); 1 0 0];

cd(v.vardir); save('v.mat','v')
end