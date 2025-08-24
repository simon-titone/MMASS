function [v] = study2_param(v)

load([v.vardir '/colors.mat']); 
load([v.vardir '/matrix_connectivity.mat']);
v.DM = [v.vardir '/design.mat'];

v.nsubj = length(v.subj); % number of independent subjects
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
v.stages = {'NREM2', 'NREM3', 'REM'};
v.nseeds = length(v.seedname);
v.nnetws = length(v.netwname);
v.nbands = length(v.bandname);
v.stagediff = {'NREM2_NREM3', 'NREM2_REM'  ,'NREM3_REM' };
v.stagediff_title = {'NREM2 - NREM3', 'NREM2 - REM'  ,'NREM3 - REM' };
v.stagecorr = {'NREM3-NREM2', 'REM-NREM2', 'REM-NREM3'};
% n.stagecorr = {'NREM3-NREM2', 'REM-NREM2',  'WAKE-NREM2', 'REM-NREM3', 'WAKE-NREM3', 'WAKE-REM' };
% n.stagediff = {'NREM2_NREM3', 'NREM2_REM','NREM2_WAKE'  ,'NREM3_REM', 'NREM3_WAKE', 'REM_WAKE' };
% n.stagediff_title = {'NREM2 -- NREM3', 'NREM2 -- REM','NREM2 -- WAKE'  ,'NREM3 -- REM', 'NREM3 -- WAKE', 'REM -- WAKE' };
v.nstages = 3;
v.nstagecorr = length(v.stagecorr);
v.nstagediff = length(v.stagediff);
v.ncycles = 3; 
v.levels = v.ncycles * v.nstages; % total number of groups and sessions

v.colors = colors; 
v.colors_onesided = colors_onesided;
reds = 0:0.0005:1; v.j_colors = [1 1 1; ones(length(reds),1), flipud(reds'), flipud(reds'); 1 0 0];

if strcmp(v.comp, 'simontitone')
    v.FT_dir = '/Users/simontitone/Documents/study2/matrix_conn';
else
    v.FT_dir = '/Volumes/Backup2/Analysis/study2/study2_ft_roi';
end
%label ptp number onto scatter plot dots
for cc = 1:v.nsubj              
    c(cc) = v.subj(cc);
end
v.subj_label= c';
v.violin_subj_label = [c' ; c'];
save([v.vardir filesep 'v.mat'],'v')
end