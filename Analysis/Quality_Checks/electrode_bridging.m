clear all; close all; warning off
% comp = 'simontitone'; drive = 'Portable_backup';
comp = 'u0127719';drive = 'Backup2';
% ds = '/Sleep/MMASS_SLEEP/'; a_output = 'Sleep_conn';
% ds = '/RS/RS_FINAL_SEG/';a_output = 'RS_conn';

% for i= 1:3
%     if i ==1
%         source_dir = '/study2/NET/study2_prelim_NET/'; a_output = 'study2_prelim'; n_ds = 52;
%     elseif i ==2
%         source_dir = '/study2/NET/study2_NET_badcap/'; a_output = 'study2_badcap'; n_ds = 43;
%     elseif i ==3
        source_dir = '/study2/NET/study2_NET/'; a_output = 'study2'; n_ds = 159;
%     end

savepath = ['/Users/' comp '/Google_Drive/PhD/MMASS/Analysis/study2/quality_check/electrode_bridging/' a_output ];
if ~exist(savepath); mkdir(savepath); end

for d = 1:n_ds
    savetitle = [savepath filesep 'e_bridge_dataset' num2str(d) '.png'];
    if ~exist(savetitle)
        % for d = 1:162
        fpath = ['/Volumes/' drive '/Analysis' source_dir 'dataset' num2str(d) '/eeg_source/sources_eeg.mat'];
        load(fpath);
        channel = source.sensor_data'; nchan = size(channel,2); threshold = 0.98;
        for c = 1:nchan
            for c2 = c:nchan
                t(c, c2) = corr(channel(:,c),channel(:,c2));
            end
        end
        figure, imagesc(t), colorbar, caxis([-1, 1]), colormap(redblue)
        th = zeros(nchan);
        th(t<-threshold) = t(t<-threshold); th(t>threshold) = t(t>threshold);
        figure, imagesc(th), colorbar
        title(['Electrode Bridging d' num2str(d) ]);
        cd(savepath)
        saveas(gcf, savetitle);
    end
end
close all
% end