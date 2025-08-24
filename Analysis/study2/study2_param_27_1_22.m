function [v] = study2_param(v)

v.rootdir = ['/Users/' v.comp '/Google_Drive/PhD/MMASS/'];
v.homedir = [v.rootdir 'Analysis/study2'];
v.vardir = [v.homedir '/vars'];
cd(v.vardir);
load('colors.mat'); load('matrix_connectivity.mat');

ex = [9 14 17];
v.NREM2_c2_subj = 1:29; v.NREM2_c2_subj(v.NREM2_c2_subj(ex)) = [];  % n = 26

v.NREM2_c3_excl = [ex];  % n = 26
v.NREM2_c3_subj = 1:29; v.NREM2_c3_subj(v.NREM2_c3_subj(v.NREM2_c3_excl)) = [];  % n = 26

v.NREM2_c4_excl = [ex 2 4 8 12 13];  % n = 20
v.NREM2_c4_subj = 1:29; v.NREM2_c4_subj(v.NREM2_c4_subj(v.NREM2_c4_excl)) = [];  % n = 26

v.NREM3_c2_excl = [ex 28 12 24];  % n = 23
v.NREM3_c2_subj = 1:29; v.NREM3_c2_excl = sort(v.NREM3_c2_excl, 'descend'); for i = 1:length(v.NREM3_c2_excl); v.NREM3_c2_subj(v.NREM3_c2_subj(v.NREM3_c2_excl(i))) = []; end

v.NREM3_c3_excl = unique([v.NREM3_c2_excl 2 3 5 7 11 13 15 19 28]);  % n = 15
v.NREM3_c3_subj = 1:29; v.NREM3_c3_excl = sort(v.NREM3_c3_excl, 'descend'); for i = 1:length(v.NREM3_c3_excl); v.NREM3_c3_subj(v.NREM3_c3_subj(v.NREM3_c3_excl(i))) = []; end

v.REM_c2_excl = [ex 3 4 6 10 11 12 18 28 ];  % n = 18
v.REM_c2_subj = 1:29; v.REM_c2_excl = sort(v.REM_c2_excl, 'descend'); for i = 1:length(v.REM_c2_excl); v.REM_c2_subj(v.REM_c2_subj(v.REM_c2_excl(i))) = []; end

v.REM_c3_excl = unique([v.REM_c2_excl 1 2 4 8 10 ]);  % n = 15
v.REM_c3_subj = 1:29; v.REM_c3_excl = sort(v.REM_c3_excl, 'descend'); for i = 1:length(v.REM_c3_excl); v.REM_c3_subj(v.REM_c3_subj(v.REM_c3_excl(i))) = []; end

v.REM_c4_excl = unique([v.REM_c3_excl 2 4 6 8 10 12 13 18 29]);  % n = 13
v.REM_c4_subj = 1:29; v.REM_c4_excl = sort(v.REM_c4_excl, 'descend'); for i = 1:length(v.REM_c4_excl); v.REM_c4_subj(v.REM_c4_subj(v.REM_c4_excl(i))) = []; end



v.nsubj = 29; % number of independent subjects
v.ngroup = 1;  % number of experimental groups (between subject factor)
v.nroi = size(seed_info,2); % number of ROI used in analyses
v.fdr_thres = 0.05; % note that script below assumes two-tailed tests of significance
v.netwname = {'dmn','dan','van','lang','mot','vis'};
v.netw= {1:4, 5:8, 9:10, 11:12, 13:17, 18:21};
v.nnetw    = length(v.netwname);
v.seedname = {'lANG','rANG','PCC','MPFC','lIPS','rIPS','lFEF','rFEF','rTPJ','rIFG','lTPJ','lIFG','lSMA','lCS','rCS','lS2','rS2','lV1V2','rV1V2','lMT','rMT'};
v.bandname = {'delta','theta','alpha','sigma','beta','gamma'};
delta      = 1:4; theta = 4:8; alpha = 8:12; sigma = 12:16; beta = 16:30; gamma = 30:80; 
v.band     = {delta, theta, alpha, sigma, beta, gamma}; 
v.nband= numel(v.band);
v.bandfreq = {1:4, 4:8, 8:12, 12:16, 16:30, 30:80}; 
v.nf = 80;
v.stage = {'NREM2', 'NREM3', 'REM'};
v.nseeds = length(v.seedname);
v.nnetws = length(v.netwname);
v.nbands = length(v.bandname);
v.stagediff = {'NREM2_NREM3', 'NREM2_REM'  ,'NREM3_REM' };
v.stagediff_title = {'NREM2 - NREM3', 'NREM2 - REM'  ,'NREM3 - REM' };
v.stagecorr = {'NREM3-NREM2', 'REM-NREM2', 'REM-NREM3'};
% n.stagecorr = {'NREM3-NREM2', 'REM-NREM2',  'WAKE-NREM2', 'REM-NREM3', 'WAKE-NREM3', 'WAKE-REM' };
% n.stagediff = {'NREM2_NREM3', 'NREM2_REM','NREM2_WAKE'  ,'NREM3_REM', 'NREM3_WAKE', 'REM_WAKE' };
% n.stagediff_title = {'NREM2 -- NREM3', 'NREM2 -- REM','NREM2 -- WAKE'  ,'NREM3 -- REM', 'NREM3 -- WAKE', 'REM -- WAKE' };
v.nstages = length(v.stage);
v.nstagecorr = length(v.stagecorr);
v.nstagediff = length(v.stagediff);
v.ncycles = 3; 
v.levels = v.ncycles * v.nstages; % total number of groups and sessions

v.colors = colors; 
v.colors_onesided = colors_onesided;
reds = 0:0.0005:1; v.j_colors = [1 1 1; ones(length(reds),1), flipud(reds'), flipud(reds'); 1 0 0];

%label ptp number onto scatter plot dots
for cc = 1:v.nsubj              
    c(cc) = v.subj(cc);
end
v.subj_label= c';
v.violin_subj_label = [c' ; c'];
cd(v.vardir); save('v.mat','v')
end