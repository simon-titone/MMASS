clearvars; warning off; clc; close all;
comp = 'u0127719'; vol = 'Backup2';
% comp = 'simontitone'; vol = 'Portable_Backup';
study = 'study3';

source_dir = '/Volumes/Backup2/Analysis/study3/SM_NET/'; a_output = 'study3_v3'; n_ds = 67;
savedir = ['/Users/' comp '/Google_Drive/PhD/MMASS/Analysis/' study '/quality_checks/PSD/' a_output];
if ~exist(savedir, 'dir'); mkdir(savedir);end % cd(savedir);

% for i = 1:n_ds
for i = 22
    savefile = [savedir filesep 'PSD_dataset' num2str(i) '.png'];
    raw = [source_dir  'dataset' num2str(i) filesep 'eeg_signal/raw_eeg.dat'];
    processed = [source_dir  'dataset' num2str(i) filesep 'eeg_signal/processed_eeg.dat'];
%     if ~exist(savefile)
        net_plotPSD(raw, processed);
        title = mtit(['dataset' num2str(i)]);
        saveas(gcf, savefile);
%     end
end
close all;