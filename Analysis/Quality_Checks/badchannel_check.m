clearvars;clc;
study = 'study2'; comp = 'u0127719';
% for i= 1:3
%     if i ==1
%         source_dir = '/Volumes/Backup2/Analysis/study2/NET/study2_prelim_NET/'; a_output = 'study2_prelim'; n_ds = 52;
%     elseif i ==2
        source_dir = '/Volumes/Backup2/Analysis/study2/NET/study2_NET_badcap/'; a_output = 'study2_badcap'; n_ds = 43;
%     elseif i ==3
%         source_dir = '/Volumes/Backup2/Analysis/study2/NET/study2_NET/'; a_output = 'study2'; n_ds = 159;
%     end

    savepath = ['/Users/' comp '/Google_Drive/PhD/MMASS/Analysis/study2/quality_check/badchannel/' a_output ];
    for i= 1:n_ds
        fpath = [source_dir 'dataset' num2str(i) '/eeg_signal/processed_eeg_badchannel.mat' ];
        load(fpath)
        if isempty(badchanind)
            badchanind = 0;
        end
        l = length(badchanind);
        bad_channels(i,1:l) = badchanind;
    end
    bad_channels(bad_channels == 0) = NaN;
    save(savepath, 'bad_channels')
% end