clearvars; warning off; clc;close all;
comp = 'u0127719'; vol = 'Backup2';
% comp = 'simontitone'; vol = 'Portable_Backup';
study = 'study2';

% source_dir = '/Volumes/Backup2/Analysis/study2/NET/study2_NET/'; a_output = 'study2';
%         source_dir = '/Volumes/Backup2/Analysis/study2/NET/study2_prelim_NET/'; a_output = 'study2_prelim';
%         source_dir = '/Volumes/Backup2/Analysis/study2/NET/study2_NET_badcap/'; a_output = 'study2_badcap';
% source_dir = '/Volumes/Backup2/Analysis/study2/NET/PSD_corrected/badcap/'; a_output = 'study2_badcap_corrected';
source_dir = '/Volumes/Backup2/Analysis/study2/NET/PSD_corrected/reg/'; a_output = 'study2_reg_corrected';

savedir = ['/Users/' comp '/Google_Drive/PhD/MMASS/Analysis/' study '/quality_check/PSD/' a_output];
if ~exist(savedir, 'dir'); mkdir(savedir);end % cd(savedir);

i = 35;
savefile = [savedir filesep 'PSD_dataset' num2str(i) '.png'];
raw = [source_dir  'dataset' num2str(i) filesep 'eeg_signal/raw_eeg.dat'];
processed = [source_dir  'dataset' num2str(i) filesep 'eeg_signal/processed_eeg.dat'];
net_plotPSD(raw, processed);
title = mtit(['dataset' num2str(i)]);

savefile = [savedir filesep 'PSD_dataset' num2str(i) '.png'];
saveas(gcf, savefile);
