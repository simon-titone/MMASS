function RS_power_analysis(n,subj)
warning off; close all;
powerdir = ['/Users/' n.comp '/Google_Drive/PhD/MMASS/Experiment/Analysis/RS_conn/Power/' n.analysis ];
outlierdir = [powerdir '/outliers']; if ~exist(outlierdir); mkdir(outlierdir); end

FT_dir = [n.homedir 'matrix_conn'];
sourcepath = ['/Volumes/' n.drive '/Analysis/RS/' n.analysis];
n.spath=  ['/Users/' n.comp '/Google_Drive/PhD/MMASS/Experiment/Analysis/RS_conn/Power/options.mat'];
log_plot_path = [n.homedir 'Power/' n.analysis '/plot_log'];

%%%% extracts FT_roi.mat
if n.FT_extract == 1
    outdir = [n.homedir 'Power/' n.analysis '/data' ];
    dat = n.subj*n.levels;
    for i = 1:length(dat)
        fpath = [sourcepath '/dataset' num2str(i) '/eeg_source/sources_eeg.mat' ];
        cd([outdir '/dataset' num2str(dat(i))]);
        ST_net_seed_connectivity(fpath, n.spath,n);
    end
end

for i = 1:n.subj
    for r = 1:n.levels
        s = subj(i);
        load( [n.homedir filesep 'matrix_conn' filesep 'subj' num2str(s) '_run' num2str(r) filesep 'FT_roi.mat'] )
        
        for ii = 1:21
            power_seed(ii,:,:) = FT_roi(ii,:,:).*conj(FT_roi(ii,:,:));
        end
        PSD_temp(:,:,i,r) = nanmean(power_seed(:,:,:),3);
        for nt = 1:n.nnetw
            PSD_netw(nt,:,i,r) = nanmean(PSD_temp(n.netwrange{nt},:,i,r),1);
        end
        for b = 1:n.nband
            PSD(:,b,i,r) = nanmean(PSD_netw(:,n.bandfreq{b},i,r),2);
        end
        %     cd(indvdir)
        %     figure
        %     for nt = 1:n.nnetw
        %         x= squeeze(nanmean(PSD_netw(nt,:,:),3));
        %         plot(x(1:40)) %individual plot
        %     end
    end
end

% Plots PSD with Jessica's format
if n.logplot == 1
    if ~exist(log_plot_path); mkdir(log_plot_path); end; cd(log_plot_path);
    
    netwpowpre = nanmean(PSD_netw(:,:,:,1),3);
    stdnetwpowpre = std(PSD_netw(:,:,:,1),[],3,'omitnan');
    sup_pct = netwpowpre + stdnetwpowpre/sqrt(n.subj);
    slo_pct = netwpowpre - stdnetwpowpre/sqrt(n.subj);
    
    for nt = 1:n.nnetw
        figure
        x = 1:80; y1 = squeeze(sup_pct(nt,:)); y2 = squeeze(slo_pct(nt,:));
        hold on, plot(x,y1,x,y2)
        fill([x,fliplr(x)],[y1,fliplr(y2)], [0.8 0.8 0.8]), set(gca,'YScale','log')
        semilogy(netwpowpre(nt,:),'color',[0.2 0.1 0.2],'linewidth',2)
        
        xlabel('Hz'); xlim([1 48]);
        title(n.netwname{nt});
        ax = gca;
        ax.TickLength = [0 0];
        %     daspect([1 .015 1 ])
        saveas(gcf, [n.netwname{nt}  '_log_PSD.png' ]);
    end
end
clear xlim


cd(outlierdir);
for r = 1:n.levels
    for nt = 1:n.nnetw
        for b = 1:n.bands
            out = squeeze(PSD(nt,b,:,r));
            [PSD_max(:,1), PSD_max(:,2)] = max(out); [PSD_min(:,1), PSD_min(:,2)] = min(out);
            cut = std(out) * 3; lowcut = mean(out) - cut; highcut = mean(out) + cut;
            if PSD_max(:,1) > highcut
                outliers(nt,b,r) = subj(PSD_max(1,2)) ;
                outliers_raw(nt,b,r) = PSD_max(1,2) ;
            elseif  PSD_min < lowcut
                outliers(nt,b,r) = subj(PSD_min(1,2));
                outliers_raw(nt,b,r) = PSD_min(1,2);
            else
                outliers(nt,b,r) = 0;
                outliers(nt,b,r) = 0;
            end
            if n.scatter == 1
                if PSD_max > highcut | PSD_min < lowcut
                    figure
                    maxval(1,:) = out(:);
                    for i = 1:length(out)
                        scatter(repmat(1,1,length(out)),out,'o', 'b');
                    end
                    xlim([0.5 1.5]); ylim([(highcut + (highcut/2)) (lowcut - (lowcut/2))]);
                    line([0.5 1.5],[lowcut lowcut],'Color','k','LineStyle','--');
                    line([0.5 1.5],[highcut highcut],'Color','k','LineStyle','--');
                    title([n.bandname{b} ' '  n.netwname{nt}  ' PSD outliers' ]);
                    ylabel('power'); xlabel('participant');
                    labelpoints(1,out,subj,'E', .05);
                    saveas(gcf, [n.bandname{b} '_'  n.netwname{nt}  '_PSD_outliers.png' ]);
                end
            end
        end
    end
end
cd(outlierdir);
outliervar = ['RS_power_outliers_' n.analysis];
save(outliervar, 'outliers');

% Plot per network averages across participants
plotdir = [powerdir '/plot']; if ~exist(plotdir); mkdir(plotdir); end; cd(plotdir);
for nt = 1:n.nnetw
    for r = 1:n.levels
        y = outliers_raw(nt,:,r); y(y==0) = []; y = unique(y); % uses the order based on the count from the script
        o = outliers(nt,:,r); o(o==0) = []; o = unique(o); % uses the actual participant numbers
        stderr = (std(PSD_netw(nt,:,:,r),0,3)/sqrt(length(PSD_netw(nt,:,:,r))));
        x = nanmean(PSD_netw(nt,:,:,r),3);
        figure; hold on
        h= zeros(length(y),1);
        shadedErrorBar([],x(1:40),stderr(1:40));
        h(1) = plot(x(1:40));
        title([n.netwname{nt} ' run ' num2str(r)]);
        xlim([0 40]);
        if n.plot_outliers == 1
            cd(outlierdir);
            for iii = 1:length(y)
                b = PSD_netw(nt,:,y(iii),r);
                h(iii+1) = plot(b(1:40));
                name{iii} = {[ 'P' num2str(o(1,iii)) ]};
            end
            if length(y) == 1; l = legend(h, 'mean', char(name{1}));
            elseif length(y) == 2; l = legend(h, 'mean', char(name{1}),char(name{2}));
            elseif length(y) == 3; l = legend(h, 'mean', char(name{1}),char(name{2}),char(name{3}));
            elseif length(y) == 4; l = legend(h, 'mean', char(name{1}),char(name{2}),char(name{3}),char(name{4}));
            elseif length(y) == 5; l = legend(h, 'mean', char(name{1}),char(name{2}),char(name{3}),char(name{4}),char(name{5}));
            elseif length(y) == 6; l = legend(h, 'mean', char(name{1}),char(name{2}),char(name{3}),char(name{4}),char(name{5}),char(name{6}));
            end
            saveas(gcf,['PSD_' n.netwname{nt} '_run' num2str(r) '_outliers.png']);
        elseif n.plot_outliers == 0
            saveas(gcf,['PSD_' n.netwname{nt} '_run' num2str(r) '.png']);
        end
    end
    close all;
end
end
