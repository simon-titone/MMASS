%%%%%% NETW(network, network, subject, run, band)
%%%%%% BANDS = 1-delta 2-theta 3-alpha 4-beta 5-gamma 6-bb
clearvars; close all;clc;
%% ADJUST:
n.comp = 'simontitone'; n.drive = 'Portable_Backup';
% n.comp = 'u0127719'; n.drive = 'Backup2';

% n.group = [1 2 10 12 13 17 18 22];                    % Excluded participants (ANOVA, n = 21)
n.group = [1 2 10 12 13 17 18 22 7 14 16];              % Excluded participants (online, n = 18)
% n.group = [1 2 10 12 13 17 18 22 7 14 16 3 6 20];     % Excluded participants (offline, n = 15)

% Optional parameters (1 = on, 0 = off)
n.scatter = 1; 
n.boxplot = 0; 
n.out_cutoff = 3; 
n.violin = 0; 
n.FDR = 'all'; % 'all' or 'mot'
n.subj_labels = 1; 
n.pval_labels = 1; 
n.title = 1;
% % --------------------------------------------------------------------------
%% Parameters
subj = 1:29; 
n.group = sort(n.group, 'descend'); for i = 1:length(n.group); subj(subj(n.group(i))) = []; end
RS_param(n,subj); 
load('bhv_vars.mat'); load('n.mat');
speed.online = on_ITI(subj);speed.offline = off_ITI(subj); % speed.mean = mean_ITI(subj);
RS_create_design(n); load('design.mat');

%% Load connectivity matrices
RS_load(n, subj) ; load('n.mat');

%% Outlier Analyses
% RS_scatter_conn_outliers(NETW,subj,n);              % Connectivity outliers
% RS_intersession_outliers(NETW,n,subj);              % Intersession connectivity outliers

%% CONNECTIVITY ANALYSES
% RS_ANOVA(NETW,design_matrix,n);                     % ANOVA - n.FDR option not set up yet, currently set to MOT only
% RS_within_corr_bhv(NETW,speed,n);                   % Within-level correlations w/ behaviour
% RS_bhv_intersession_corr(NETW,speed,n);             % Behavioral Correlations with inter-session changes in connectivity

%% Behavioral
% scatter_bhv_outliers(speed,n);                      % creates scatter plots of behavioral outliers
% bhv_group_plot(n, subj);                              % Plot behavior at the group level
% bhv_group_stats(n, subj);                           % stats
% SRT(n,subj);                                        % Plots SRT bhv data
% PVT(n,subj);                                        % Saves .csv file with PVT data

%% Power analysis
% n.analysis = 'RS_FINAL_SEG'; n.FT_extract = 0;  n.scatter = 0; n.plot_outliers = 1; n.logplot = 0;
% RS_power_analysis(n,subj);                          % performs power analysis 

%% Optional
% RS_indv_subj_matrix(NETW,n);                        %needs adjusting
% RS_t_test_sess(NETW,n);                             %optional alternative to ANOVA -shows directionality of differences btwn groups

% RS_within_vs_btw(NETW,subj,n,design_matrix);        % Within Vs. between network connectivity; option for boxplot- only usable with R2018+
% RS_within_level(NETW,design_matrix,n);              % Within-level statistics (tcdf needs R2018+ to run properly)
% RS_within_level_means(NETW, n);                     % Within-level means

%%%%%% ST May 2021 - added SRT & PVT, added individual groupings for subject selection
%%%%%% ST Feb 2021 - added Power analysis, behavioral group plots & ANOVA scatter plots
%%%%%% ST Jan 2021 - Changed data loading & implemented functions
%%%%%% ST May 2020 - created loops for behavioral and bands for 3F, 3E 
%%%%%% ADJUSTED BY GA FEBRUARY 2020 - FENS ABSTRACT                           