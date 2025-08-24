function study2_within_level_means(v, NETW)
%%%% NETW: network, network, subject, stage, cycle, band
groupmeans = [v.homedir '/results/group_means/'  v.t_group];
if ~exist(groupmeans); mkdir(groupmeans); end; cd(groupmeans)

for c = v.cycle
    for st = v.stage
        for b = 1:v.nbands
            for w = 1:v.nnetw
                for ww = 1:v.nnetw
                    groupmean_netw(w,ww,st,c,b) = nanmean(NETW(w,ww,:,st,c,b),3);
                end
            end
        end
    end
end

for c = 1:v.cycle
    for st = 1:v.stage
        for b = 1:v.nbands
            matrix = groupmean_netw(:,:,st,c,b);
                        tri_matrix = triu(matrix);

            figure(); set(gcf,'Color','white'); box OFF;
            imagesc(tri_matrix); axis square; %caxis([-mm mm]);
            colorbar; daspect([1 1 1]); colormap(v.j_colors);
            hold on;
            for b2=1:v.nnetws    line([b2+0.5 b2+0.5],[0.5 v.nnetw+0.5],'color','k');   line([0.5 v.nnetw+0.5],[b2+0.5 b2+0.5],'color','k');    end
            set(gca,'YTick',[1:v.nnetw],'YTickLabel',v.netwname,'Fontsize',16);
            set(gca,'XTick',[1:v.nnetw],'XTickLabel',v.netwname,'Fontsize',16);
            title(['Group Means cycle' num2str(c) ' ' v.stages{st} ' ' v.bandname{b}  ' n= ' num2str(v.nsubj) ]);
            saveas(gcf, ['GM_c' num2str(c) '_' v.stages{st} '_' v.bandname{b} '.png']);
        end
    end
end
close all
end