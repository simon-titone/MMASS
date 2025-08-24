function RS_within_level_means(NETW, n)

groupmeans = [n.homedir filesep 'results' filesep 'group_means_n' num2str(n.subj)];
if ~exist(groupmeans); mkdir(groupmeans); end; cd(groupmeans)

for r = 1:n.levels
    for b = 1:n.bands
        groupmean_netw(:,:,r,b) = nanmean(NETW(:,:,:,r,b),3);
    end
end

for b = 1:n.bands
    for r = 1:n.levels
        matrix = groupmean_netw(:,:,r,b);
        figure(); set(gcf,'Color','white'); box OFF;
        imagesc(matrix); axis square; %caxis([-mm mm]);
        colorbar; daspect([1 1 1]); colormap(n.j_colors);
        
        hold on;
        for b2=1:n.netws    line([b2+0.5 b2+0.5],[0.5 n.netws+0.5],'color','k');   line([0.5 n.netws+0.5],[b2+0.5 b2+0.5],'color','k');    end
        set(gca,'YTick',[1:n.netws],'YTickLabel',n.netwname,'Fontsize',16);
        set(gca,'XTick',[1:n.netws],'XTickLabel',n.netwname,'Fontsize',16);
        title(['Group Means ' n.bandname{b} ' Netw run' num2str(r) ' n= ' num2str(n.subj) ]);
        
        %         [~,x] = max(matrix,[],2); % max of each row
        %         plot(x,1:n.netws,'ko','MarkerSize',8,'MarkerFaceColor','k');
        %         [~,x] = max(matrix); % max of each column
        %         plot(1:n.netws,x,'ko','MarkerSize',8,'MarkerFaceColor','k');
        %
        saveas(gcf, ['Group_means_' n.bandname{b} '_netw_run_' num2str(r) '.png']);

    end
end

close all;
end