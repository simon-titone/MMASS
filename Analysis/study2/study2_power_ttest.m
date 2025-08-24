function study2_power_ttest(v)
warning off; close all;
if v.comp == 'simontitone' FT_dir = '/Users/simontitone/Documents/study2/matrix_conn';else
    FT_dir = '/Volumes/Backup2/Analysis/study2/study2_ft_roi'; end
savepath = [v.homedir '/results/power/' v.t_group '/'];
if ~exist(savepath); mkdir(savepath);end
sourcepath = ['/Volumes/' v.drive '/Analysis/RS/' v.analysis];
PSD = nan(v.nnetw, v.nbands,v.nsubj,v.nstages,v.ncycles);
for i = 1:v.nsubj
    for st= 1:v.stage
        for c = 1:v.ncycles
            s = v.subj(i);
            file = [FT_dir filesep 'P' num2str(s) '_' v.stages{st} '_c' num2str(c) filesep 'FT_roi.mat'];
            if exist(file); load( file );end
            for ii = 1:21
                power_seed(ii,:,:) = FT_roi(ii,:,:).*conj(FT_roi(ii,:,:));
            end
            PSD_temp(:,:,i,st,c) = nanmean(power_seed(:,:,:),3); %calculates average across time series
            for nt = 1:v.nnetw
                PSD_netw(nt,:,i,st,c) = nanmean(PSD_temp(v.netw{nt},:,i,st,c),1); %averages across networks
            end
            for b = 1:v.nband
                PSD(:,b,i,st,c) = nanmean(PSD_netw(:,v.bandfreq{b},i,st,c),2); %averages across bands
            end
        end
    end
end
clear nt;clear b;clear i;clear st;clear c; clear s;clear ii

for b = 1:v.nbands
    for st= 1:v.nstages
        for c = 1:v.ncycles
            for nt = 1:v.nnetw
                internetw = setdiff(1:v.nnetw,nt);
                for n2 = 1:length(internetw)
                    nt2 = internetw(n2);
                    r1 = PSD(nt,b,:,st,c);
                    r2 = PSD(nt2,b,:,st,c);
                    [h,p,ci,stats] = ttest(r1,r2,'tail','both');
                    inter_h(nt,nt2,b,st,c) = h;
                    inter_p(nt,nt2,b,st,c) = p;
                    inter_t(nt,nt2,b,st,c) = stats.tstat;
                end
            end
        end
    end
end
for st= 1:v.nstages
    for c = 1:v.ncycles
        for b = 1:v.nbands
            figure(); set(gcf,'Color','white'); box OFF;
            imagesc(inter_t(:,:,b,st,c)); axis square;
            colorbar; daspect([1 1 1]); colormap(redblue); 
%             color_bounds = [-max(-max(squeeze(inter_t(:,:,b,st,c)))), max(max(inter_t(:,:,b,st,c)))];
%             caxis(color_bounds);
            
            hold on; for b2=1:v.nnetws    line([b2+0.5 b2+0.5],[0.5 v.nnetws+0.5],'color','k');   line([0.5 v.nnetws+0.5],[b2+0.5 b2+0.5],'color','k');    end
            set(gca,'YTick',[1:v.nnetws],'YTickLabel',v.netwname,'Fontsize',16); set(gca,'XTick',[1:v.nnetws],'XTickLabel',v.netwname,'Fontsize',16);
            dim = set(gcf, 'Position',  [1000, 900, 500, 450]); %left bottom width height
            h = colorbar; set(get(h,'title'),'string','t-val', 'fontsize', 12); %ylabel(h, 'f-val');

            title({[v.bandname{b} ' ' v.stages{st} ' cycle' num2str(c) ] [' Power t-test, n= ' num2str(v.nsubj) ]});
            
            [ivect_unc,jvect_unc]=find(inter_h(:,:,b,st,c)==1);
            scatter(jvect_unc,ivect_unc, 80, 'w', 'o' ,'ok','LineWidth',2.5); clear ivect_unc; clear jvect_unc;
            clear ivect_unc; clear jvect_unc;
            
            %     [ivect,jvect]=find(p_corr<=0.05);
            %     scatter(jvect,ivect,9,'ok','filled');
            %     clear ivect; clear jvect;
            
            save_name = [savepath v.bandname{b} '_' v.stages{st} '_c' num2str(c) '_power_n' num2str(v.nsubj) '.png'];
            saveas(gcf, save_name);
        end
    end
end
close all;
end