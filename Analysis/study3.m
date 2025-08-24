clearvars;clc;close all;
%% ADJUST
v.scatter = []; v.boxplot = []; v.violin = [1]; v.out_cutoff = 3; v.excel = [];
v.subj_labels = [1]; v.pval_labels = [1];  v.fig_labels = [1]; v.fig_title = [1];
v.z_test = []; v.FDR = 'all'; 
%--------------------------------------------------------------------------
%% Parameters
v.rootdir = ['/Users/simontitone/GoogleDrive/PhD/'];
v.vardir = [v.rootdir 'Analysis/study3/vars']; var = [v.vardir 'v.mat'];v.drive = 'Portable_backup';
save(var,'v'); study3_param(v);   load('v.mat')       % Analysis parameters and variables defined
%% Outlier Analyses
% study3_conn_outliers(v, NETW);                      % Connectivity outliers

%% v.group Means
% study3_within_level_means(v)                  % Within-level means, needs updating
% study3_within_level_stats(v, NETW)                  % Within-level statistics (tcdf needs R2018+ to run properly), not finished editing yet

%% Primary Connectivity Analyses
%%%%%% ANOVA
study3_intersession_ANOVA(v)
% study3_intersession_stagediff_ANOVA(v)
% study3_interstage_indv_ANOVA(v)                             
% study3_stage_ANOVA_3way(v)
% study3_allstage_ratio_ANOVA(v)
% study3_stage_sess_av_ANOVA(v)
% study3_3stages_ANOVA_wake(v)
% study3_stage_ANOVA_all(v)
% study3_stageBYsess_ANOVA(v)

%%%% Correlations
% study3_within_lvl_corr(v)
% study3_intersession_corr_bhv(v)
% study3_intersession_corr(v)
study3_interstage_corr_bhv(v)                       % Now has z-test & partial correlation 
% study3_interstage_corr(v)
% study3_interstage_NREM2_corr(v)
% study3_interstage_corr_wake(v)
% study3_interstage_intersess_diff_corr_bhv(v)
% study3_ratio_bhv(v)

%% Behavioral
% study3_scatter_bhv_outliers(speed,n);                      % Creates scatter plots of behavioral outliers
% study3_bhv_group_plot(v);                                  % Plot behavior at the v.group level
% study3_bhv_group_stats(v);                                 % Behavioral data stats, saves .mat file
% SRT(n,subj);                                        % Plots SRT bhv data
% PVT(n,subj);                                        % Saves .csv file with PVT data

%% Power analysis
% v.analysis = 'MMASS_SLEEP'; v.FT_extract = 0;     
% v.scatter = 0; v.plot_outliers = 1; v.logplot = 0;
% study3_power_analysis(v);     

%%%%% %NETW: network, network, subject, stage, session, cycle, band
%%%%% BANDS = 1-delta 2-theta 3-alpha 4-sigma 5-beta 6-gamma