function study3_within_level_means(v, NETW)
%%%%% NETW: network, network, subject, stage, session, band
groupmeans = [v.homedir filesep 'results' filesep 'group_means_n' num2str(v.nsubj)];
if ~exist(groupmeans); mkdir(groupmeans); end; cd(groupmeans)

for se = 1:v.nsess
    for st = 1:v.nstages
        for b = 1:v.nbands
            for w = 1:v.nnetw
                for ww = 1:v.nnetw
                    groupmean_netw(w,ww,st,se,b) = nanmean(NETW(w,ww,:,st,se,b),3);
                end
            end
        end
    end
end
for se = 1:v.nsess
    for st = 1:v.nstages
        for b = 1:v.nbands
                    matrix = groupmean_netw(:,:,st,se,b);
        figure(); set(gcf,'Color','white'); box OFF;
        imagesc(matrix); axis square; %caxis([-mm mm]);
        colorbar; daspect([1 1 1]); colormap(v.j_colors);
        hold on;
        for b2=1:v.nnetws    line([b2+0.5 b2+0.5],[0.5 v.nnetw+0.5],'color','k');   line([0.5 v.nnetw+0.5],[b2+0.5 b2+0.5],'color','k');    end
        set(gca,'YTick',[1:v.nnetw],'YTickLabel',v.netwname,'Fontsize',16);
        set(gca,'XTick',[1:v.nnetw],'XTickLabel',v.netwname,'Fontsize',16);
        title(['Group Means ' v.sess{se} ' ' v.stages{st} ' ' v.bandname{b}  ' n= ' num2str(v.nsubj) ]);
        saveas(gcf, ['Group_means_' v.sess{se} '_' v.stages{st} '_' v.bandname{b} '.png']);

        end
    end
end
close all
end