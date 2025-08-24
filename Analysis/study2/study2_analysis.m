%%%%% NETW: network, network, subject, stage, cycle, band
%%%%% BANDS: 1-delta 2-theta 3-alpha 4-sigma 5-beta 6-gamma
clearvars;clc;close all
%% ADJUST
% v.comp = 'simontitone'; v.drive = 'Portable_backup';
v.comp = 'u0127719'; v.drive = 'Backup2';
analysis = 's'; % c = cycle or s = stage

v.scatter = [1]; v.violin = []; v.fig_labels = [1]; v.fig_title = [1]; v.write = [];
v.out_cutoff = 3; v.plot_outliers = 1; v.logplot = 0;

%%%%% Excluded Subjects
study2_excl(v);
v.vardir = ['/Users/' v.comp '/Google_Drive/PhD/MMASS/Analysis/study2/vars']; v.var = [v.vardir '/v.mat'];load(v.var)
%% Inter-cycle analysis:
%%%% 1 = NREM2_c1-2;  2 = NREM2_c2-3; 3 = NREM2 c1-3; 4 = NREM3_c1_c2; 5 = NREM3_c2_c3; 6 = NREM3_c1_c3;  7 = REM_c1_c2; 8 = REM_c2_c3; 9 = REM_c1_c3
if strcmp(analysis ,'c')
    %     for g = 1:v.n_intercycle_groups
    for g = 6
        v.subj = v.subj_IC(g,:); v.subj= v.subj(~isnan(v.subj)); v.cycle = v.cycle_IC(g); v.stage = v.stage_IC(g); v.t_group = v.group_name_IC{g};
        %%%%% Parameters
        study2_param(v); load(v.var)                                  % Analysis parameters and variables defined
        study2_DM(v); load(v.DM);                                     % Creates design matrices
        study2_load(v); load(v.var)                                   % Loads connectivity matrices and groups into networks
        %%%%% Analyses
        study2_cycle_ANOVA(v, NETW, design_matrix)                  % compares across cycles
        % study2_conn_outliers(v, NETW);                              % Connectivity outliers
        % study2_within_level_means(v, NETW)                          % Within-level means
        % study2_within_level_stats(v, NETW)                          % Within-level statistics (tcdf needs R2018+ to run properly), not finished editing yet
        % study2_power_analysis(v);
        % study2_power_ttest(v);
        % study2_power_intercycle_ttest(v)
    end
end

%% Inter-stage Analysis:
%%%% 1- c1_N2_N3; 2- c1_N2_REM'; 3-c1_N3_REM; 4-c2_N2_N3; 5-c2_N2_REM; 6-c2_N3_REM; 7-c3_N2_N3; 8-c3_N2_REM; 9-c3_N3_REM
if strcmp(analysis ,'s')
%     for g = 1:v.n_interstage_groups
    for g = 7
        v.subj = v.subj_IS(g,:); v.subj= v.subj(~isnan(v.subj)); v.cycle = v.cycle_IS(g); v.stage = v.stage_IS(g); v.t_group = v.group_name_IS{g};
        %%%% Parameters
        study2_param(v); load(v.var)                                  % Analysis parameters and variables defined
        study2_DM(v); load(v.DM);                                     % Creates design matrices
        study2_load(v); load(v.var)                                   % Loads connectivity matrices and groups into networks
        %%%% Analyses
        study2_interstage_ANOVA(v, NETW, design_matrix)               % Compares across stages within each cycle
        % study2_conn_outliers(v, NETW);                              % Connectivity outliers
        % study2_within_level_means(v, NETW)                          % Within-level means
        % study2_within_level_stats(v, NETW)                          % Within-level statistics (tcdf needs R2018+ to run properly), not finished editing yet
        % study2_power_analysis(v);
        % study2_power_ttest(v);
        % study2_power_interstage_ttest(v)
    end
end