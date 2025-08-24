function study2_power_interstage_ttest(v)
warning off; close all;
if v.comp == 'simontitone' FT_dir = '/Users/simontitone/Documents/study2/matrix_conn';else
    FT_dir = '/Volumes/Backup2/Analysis/study2/study2_ft_roi'; end
savepath = [v.homedir '/results/power/interstage/'];
if ~exist(savepath); mkdir(savepath);end
PSD = nan(v.nnetw, v.nbands,v.nsubj,v.nstages,v.ncycles);
for i = 1:v.nsubj
    for st= 1:v.nstages
        for c = 1:v.ncycles
            s = v.subj(i);
            file = [FT_dir filesep 'P' num2str(s) '_' v.stages{st} '_c' num2str(c) filesep 'FT_roi.mat'];
            if exist(file); load( file );end
            for ii = 1:21
                power(ii,:,:) = FT_roi(ii,:,:).*conj(FT_roi(ii,:,:));
            end
            PSD_temp(:,:,i,st,c) = nanmean(power(:,:,:),3); %calculates average across time series
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
    for st= v.stage
        for c = v.cycle
            for nt = 1:v.nnetws
                if st ==1; st1 = 1; st2 = 2; 
                elseif st ==2; st1 = 1; st2 = 3; 
                elseif st ==3; st1 = 2; st2 = 3; 
                end
                r1 = squeeze(PSD(nt,b,:,st1,c));
                r2 = squeeze(PSD(nt,b,:,st2,c));
                [h,p,ci,stats] = ttest(r1,r2,'tail','both');
                inter_h(nt,b,st,c) = h;
                inter_p(nt,b,st,c) = p;
                inter_t(nt,b,st,c) = stats.tstat;
            end
        end
    end
end

for st= v.stage
    for c = v.cycle
        if st == 1; comp = 'N2_N3'; comp_t1 = 'NREM2' ; comp_t2 = 'NREM3' ;
        elseif st == 2; comp = 'N2_REM'; comp_t1 = 'NREM2' ; comp_t2 = 'REM' ;
        elseif st == 3; comp = 'N3_REM'; comp_t1 = 'NREM3' ; comp_t2 = 'REM' ;
        end
        figure(); set(gcf,'Color','white'); box OFF;
        imagesc(inter_t(:,:,st,c)); axis square;
        colorbar; daspect([1 1 1]); colormap(redblue);
        %             color_bounds = [-max(-max(squeeze(inter_t(:,:,b,st,c)))), max(max(inter_t(:,:,b,st,c)))];
        %             caxis(color_bounds);
        
        hold on; for b2=1:v.nnetws    line([b2+0.5 b2+0.5],[0.5 v.nnetws+0.5],'color','k');   line([0.5 v.nnetws+0.5],[b2+0.5 b2+0.5],'color','k');    end
        set(gca,'YTick',[1:v.nnetws],'YTickLabel',v.netwname,'Fontsize',16); set(gca,'XTick',[1:v.nbands],'XTickLabel',v.bandname,'Fontsize',16);
        dim = set(gcf, 'Position',  [1000, 900, 500, 450]); %left bottom width height
        h = colorbar; set(get(h,'title'),'string','t-val', 'fontsize', 12); %ylabel(h, 'f-val');
        
        title({['cycle ' num2str(c) ' ' comp_t1 ' x ' comp_t2 ] [' Power t-test, n= ' num2str(v.nsubj) ]});
        
        [ivect_unc,jvect_unc]=find(inter_h(:,:,st,c)==1);
        scatter(jvect_unc,ivect_unc, 80, 'w', 'o' ,'ok','LineWidth',2.5); clear ivect_unc; clear jvect_unc;
        clear ivect_unc; clear jvect_unc;
        
        %     [ivect,jvect]=find(p_corr<=0.05);
        %     scatter(jvect,ivect,9,'ok','filled');
        %     clear ivect; clear jvect;
        
        save_name = [savepath  'c' num2str(c) '_' comp '_power_n' num2str(v.nsubj) '.png'];
        saveas(gcf, save_name);
        
        if v.scatter == 1
            power_intersess_scatter = [savepath 'scatter'];
            if ~exist(power_intersess_scatter); mkdir(power_intersess_scatter); end; cd(power_intersess_scatter);
            for w=1:v.nnetws
                if inter_h(w,:,st,c)==1
                    x = squeeze(NETW(w,ww,:,st,cyc1,b)); y = squeeze(NETW(w,ww,:,st,cyc2,b));
                    f = figure; hold on; colormap summer;
                    plot([1 2], [x y])
                    h = bar([1 2], [nanmean(x) nanmean(y)], .5);
                    scatter(repmat(1,1,length(x)), x, 'o', 'r');
                    labelpoints(1,x,v.subj_label,'E', .05);
                    scatter(repmat(2,1,length(y)),y,'o', 'b');
                    labelpoints(2,y,v.subj_label,'E', .05);
                    title({[ v.bandname{b} ' ' v.stages{st} ] [ v.netwname{w} '-' v.netwname{ww}  ', n= ' num2str(nsubj)]}, 'FontSize', 16);
                    ylabel('connectivity (z-score)', 'FontSize', 16);
                    grouporder={comp_t1,comp_t2};
                    xticks([1 2]); xticklabels(grouporder);
                    set(gcf, 'Position',  [10, 10, 300, 500]);
                    text(0.4, 0.98,['f-value=' num2str(ANOVA_F(w,ww,st,c,b))],'Units','normalized' );
                    text(0.4, 0.96,[ 'p-value(uncorr)=  ' num2str(ANOVA_p(w,ww,st,c,b))],'Units','normalized' );
                    saveas(gcf, [comp v.bandname{b} '_' v.netwname{w} '_' v.netwname{ww}   '_bar_scatter.png' ]);
                    clear('x', 'y', 'f');
                end
            end
        end
        close all;
    end  
end
close all;
end