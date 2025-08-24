clear all; close all; warning off;clc
comp = 'u0127719';drive = 'Backup2';
source_dir = '/study3/SM_NET/'; a_output = 'study3_v3'; n_ds = 67;
savepath = ['/Users/' comp '/Google_Drive/PhD/MMASS/Analysis/study3/quality_checks/electrode_bridging/' a_output ];
if ~exist(savepath); mkdir(savepath); end

for d = 1:n_ds
    savetitle = [savepath filesep 'e_bridge_dataset' num2str(d) '.png'];
    if ~exist(savetitle)
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