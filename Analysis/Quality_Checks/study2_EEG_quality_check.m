clearvars; warning off;
comp = 'u0127719'; vol = 'Backup2';
% comp = 'simontitone'; vol = 'Portable_Backup';
study = 'study2';

for ii= 1:3
    if ii ==1
        source_dir = '/Volumes/Backup2/Analysis/study2/NET/study2_prelim_NET/'; a_output = 'study2_prelim'; n_ds = 52;
    elseif ii ==2
        source_dir = '/Volumes/Backup2/Analysis/study2/NET/study2_NET_badcap/'; a_output = 'study2_badcap'; n_ds = 43;
    elseif ii ==3
        source_dir = '/Volumes/Backup2/Analysis/study2/NET/study2_NET/'; a_output = 'study2'; n_ds = 159;
    end
    
    savedir = ['/Users/' comp '/Google_Drive/PhD/MMASS/Analysis/' study '/quality_check/PSD/' a_output];
    if ~exist(savedir, 'dir'); mkdir(savedir);end % cd(savedir);
    
    for i = 1:n_ds
        savefile = [savedir filesep 'PSD_dataset' num2str(i) '.png'];
        raw = [source_dir  'dataset' num2str(i) filesep 'eeg_signal/raw_eeg.dat'];
        processed = [source_dir  'dataset' num2str(i) filesep 'eeg_signal/processed_eeg.dat'];
        if ~exist(savefile)
            net_plotPSD(raw, processed);
            title = mtit(['dataset' num2str(i)]);
            saveas(gcf, savefile);
        end
    end
    close all;
end
