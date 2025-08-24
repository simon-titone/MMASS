analysis = 'study2'; 
comp = 'u0127719';
vol_path = ['/Volumes/Backup2/Analysis/' analysis '/MMASS_SLEEP/dataset' ];
homedir = ['/Users/' comp '/Google_Drive/PhD/MMASS/Experiment/Analysis/' analysis '/quality_checks/badchannel_check'];
outputdir = ['/Users/' comp '/Google_Drive/PhD/MMASS/Experiment/Analysis/' analysis '/quality_checks/badchannel_check'];
n_dat = 162;

% % move files from backup to google drive
% for i = 1:n_dat
% a = [vol_path num2str(i) '/eeg_signal/processed_eeg_badchannel.mat'];
% b = [homedir '/dataset' num2str(i)];
% copyfile(a, b)
% end

for i= 1:n_dat
    fpath = [vol_path num2str(i) '/eeg_signal/processed_eeg_badchannel.mat' ];
    load(fpath)
    if isempty(badchanind)
        badchanind = 0;
    end
    l = length(badchanind);
    bad_channels(i,1:l) = badchanind;
end

bad_channels(bad_channels == 0) = NaN;
if ~exist(outputdir); mkdir(outputdir); end; cd(outputdir); 
savename = ['bad_channels_' analysis '.mat'];
save(savename, 'bad_channels')
