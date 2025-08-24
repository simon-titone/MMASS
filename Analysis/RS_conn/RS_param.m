function [n] = RS_param(n,subj)

n.homedir = ['/Users/simontitone/GoogleDrive/PhD/Analysis/study1/'];
n.vardir = [n.homedir 'vars']; if ~exist(n.vardir); mkdir(n.vardir); end
cd(n.vardir);
load('matrix_connectivity.mat');load('colors.mat');load('bhv_vars.mat');
dmn        = 1:4; dan = 5:8; van = 9:10; lang = 11:12; mot = 13:17; vis = 18:21;
n.netwrange= {1:4, 5:8, 9:10, 11:12, 13:17, 18:21};
n.netw     = {dmn, dan, van, lang, mot, vis};
n.nnetw    = numel(n.netw);
n.netwname = {'DMN','DAN','VAN','LANG','MOT','VIS'};
n.seedname = {'lANG','rANG','PCC','MPFC','lIPS','rIPS','lFEF','rFEF','rTPJ','rIFG','lTPJ','lIFG','lSMA','lCS','rCS','lS2','rS2','lV1V2','rV1V2','lMT','rMT'};
delta      = 1:4; theta = 4:8; alpha = 8:13; beta = 13:30; gamma = 30:80;
n.band     = {delta, theta, alpha, beta, gamma}; 
n.nband= numel(n.band);
n.bandname = {'Delta','Theta','Alpha','Beta','Gamma'}; 
n.bandfreq = {1:4, 4:8, 8:13, 13:30, 30:80};
n.nf = 80;
n.subj     = length(subj); % number of independent subjects
% if n.subj < 17
    n.bhv_vars = {'Offline','Online'};
% else
%     n.bhv_vars = {'online'};
% end
n.sess     = 2; % number of sessions / RS scans per subject
n.group    = 1;  % number of experimental groups (between subject factor)
n.levels   = n.sess * n.group; % total number of groups and sessions
n.roi      = size(seed_info,2); % number of ROI used in analyses
n.seeds    = 21; n.netws = length(n.netwname); n.num_bhv_vars = length(n.bhv_vars);
n.bands    = length(n.bandname);
n.fdr_thres= 0.05; % note that script below assumes two-tailed tests of significance
n.colors   = colors; 
n.colors_onesided = colors_onesided;
reds = 0:0.0005:1; n.j_colors = [1 1 1; ones(length(reds),1), flipud(reds'), flipud(reds'); 1 0 0];
for cc = 1:length(subj)
    c(cc) = subj(cc);
end
n.subj_label= transpose(c);

cd(n.vardir); save('n.mat','n')