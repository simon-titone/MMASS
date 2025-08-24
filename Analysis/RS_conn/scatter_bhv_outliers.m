function scatter_bhv_outliers(speed,n)
bhvdir = [n.homedir '/results/bhv'];
if ~exist(bhvdir); mkdir(bhvdir); end; cd(bhvdir);
for bh=1:n.num_bhv_vars
    if bh==1
        behav_ByLevel = speed.offline;
    elseif bh==2
        behav_ByLevel = speed.online;
    elseif bh==3
        behav_ByLevel = speed.mean;
    end
    cut = std(behav_ByLevel) * 3;
    lowcut = mean(behav_ByLevel) - cut;
    highcut = mean(behav_ByLevel) + cut;
    figure
    for i = 1:length(behav_ByLevel)
        scatter(repmat(1,1,length(behav_ByLevel)),behav_ByLevel,'o', 'b');
    end
    xlim = [0.5 1.5]; ylim = [(highcut + 0.1) (lowcut - 0.1)];
    line([0.5 1.5],[lowcut lowcut],'Color','k','LineStyle','--');
    line([0.5 1.5],[highcut highcut],'Color','k','LineStyle','--');
    title([n.bhv_vars{bh} ]);ylabel([n.bhv_vars{bh} ]);
    labelpoints(1,behav_ByLevel,n.subj_label,'E', .05);
    saveas(gcf, [n.bhv_vars{bh} '_scatter.png' ]);    
end
close all;