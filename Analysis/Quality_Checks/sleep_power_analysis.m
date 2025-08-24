function sleep_power_analysis(v)
warning off; close all;
powerdir = ['/Users/' v.comp '/Google_Drive/PhD/MMASS/Experiment/Analysis/Sleep/quality_checks/Power/' v.analysis ];
outlierdir = [powerdir '/outliers'];if ~exist(outlierdir); mkdir(outlierdir); end;

FT_dir = [v.homedir 'matrix_conn'];
sourcepath = ['/Volumes/' v.drive '/Analysis/RS/' v.analysis];
% v.spath=  ['/Users/' v.comp '/Google_Drive/PhD/MMASS/Experiment/Analysis/Sleep/quality_checks/options.mat'];
log_plot_path = [v.homedir 'Power/' v.analysis '/plot_log'];

%%%% extracts FT_roi.mat
if v.FT_extract == 1
    outdir = [v.homedir 'Power/' v.analysis '/data' ];
    dat = v.subj*v.levels;
    for i = 1:length(dat)
        fpath = [sourcepath '/dataset' num2str(i) '/eeg_source/sources_eeg.mat' ];
        cd([outdir '/dataset' num2str(dat(i))]);
        ST_net_seed_connectivity(fpath, v.spath,n);
    end
end

for i = 1:v.nsubj
    for st= 1:v.nstages
        for se = 1:v.nsess
            s = v.subj(i);
            load( [v.homedir filesep 'matrix_conn' filesep 'subj' num2str(s) '_' v.sess{se} '_stage' num2str(st) filesep 'FT_roi.mat'] )
            for ii = 1:21
                power_seed(ii,:,:) = FT_roi(ii,:,:).*conj(FT_roi(ii,:,:));
            end
            PSD_temp(:,:,i,st,se) = nanmean(power_seed(:,:,:),3);
            for nt = 1:v.nnetw
                PSD_netw(nt,:,i,st,se) = nanmean(PSD_temp(v.netw{nt},:,i,st,se),1);
            end
            for b = 1:v.nband
                PSD(:,b,i,st,se) = nanmean(PSD_netw(:,v.bandfreq{b},i,st,se),2);
            end
            %     cd(indvdir)
            %     figure
            %     for nt = 1:v.nnetw
            %         x= squeeze(nanmean(PSD_netw(nt,:,:),3));
            %         plot(x(1:40)) %individual plot
            %     end
        end
    end
end

% Plots PSD with Jessica's format - Needs adjusting for sleep analysis
if v.logplot == 1
    if ~exist(log_plot_path); mkdir(log_plot_path); end; cd(log_plot_path);
    
    netwpowpre = nanmean(PSD_netw(:,:,:,1),3);
    stdnetwpowpre = std(PSD_netw(:,:,:,1),[],3,'omitnan');
    sup_pct = netwpowpre + stdnetwpowpre/sqrt(v.subj);
    slo_pct = netwpowpre - stdnetwpowpre/sqrt(v.subj);
    
    for nt = 1:v.nnetw
        figure
        x = 1:80; y1 = squeeze(sup_pct(nt,:)); y2 = squeeze(slo_pct(nt,:));
        hold on, plot(x,y1,x,y2)
        fill([x,fliplr(x)],[y1,fliplr(y2)], [0.8 0.8 0.8]), set(gca,'YScale','log')
        semilogy(netwpowpre(nt,:),'color',[0.2 0.1 0.2],'linewidth',2)
        
        xlabel('Hz'); xlim([1 48]);
        title(v.netwname{nt});
        ax = gca;
        ax.TickLength = [0 0];
        %     daspect([1 .015 1 ])
        saveas(gcf, [v.netwname{nt}  '_log_PSD.png' ]);
    end
end
clear xlim


cd(outlierdir);
for st= 1:v.nstages
    for se = 1:v.nsess
        for b = 1:v.nbands
            for nt = 1:v.nnetw
                out = squeeze(PSD(nt,b,:,st,se));
                [PSD_max(:,1), PSD_max(:,2)] = max(out); [PSD_min(:,1), PSD_min(:,2)] = min(out);
                cut = std(out) * 3; lowcut = mean(out) - cut; highcut = mean(out) + cut;
                if PSD_max(:,1) > highcut
                    outliers(nt,b,st,se) = v.subj(PSD_max(1,2)) ;
                    outliers_raw(nt,b,st,se) = PSD_max(1,2) ;
                elseif  PSD_min < lowcut
                    outliers(nt,b,st,se) = v.subj(PSD_min(1,2));
                    outliers_raw(nt,b,st,se) = PSD_min(1,2);
                else
                    outliers(nt,b,st,se) = 0;
                    outliers(nt,b,st,se) = 0;
                end
                if v.scatter == 1
                    if PSD_max(1,1) > highcut | PSD_min(1,1) < lowcut
                        figure
                        maxval(1,:) = out(:);
                        for i = 1:length(out)
                            scatter(repmat(1,1,length(out)),out,'o', 'b');
                        end
                        xlim([0.5 1.5]); %ylim([ (round(lowcut) - 1) (highcut + 1)]);
                        line([0.5 1.5],[lowcut lowcut],'Color','k','LineStyle','--');
                        line([0.5 1.5],[highcut highcut],'Color','k','LineStyle','--');
                        title([ v.sess{se} ' ' v.stages{st} ' ' v.bandname{b} ' '  v.netwname{nt}  ' PSD outliers' ]);
                        ylabel('power'); xlabel('participant');
                        labelpoints(1,out,v.subj,'E', .05);
                        saveas(gcf, [v.sess{se} '_' v.stages{st} '_' v.bandname{b} '_' v.netwname{nt}  '_PSD_outliers.png' ]);
                    end
                end
            end
        end
    end
    
    
    cd(outlierdir);
    outliervar = ['RS_power_outliers_' v.analysis];
    save(outliervar, 'outliers');
end

% Plot per network averages across participants
plotdir = [powerdir '/plot']; if ~exist(plotdir); mkdir(plotdir); end; cd(plotdir);
for se = 1:v.nsess
    for st = 1:v.nstages
        for nt = 1:v.nnetw
            y = outliers_raw(nt,:,st,se); y(y==0) = []; y = unique(y); % uses the order based on the count from the script
            o = outliers(nt,:,st,se); o(o==0) = []; o = unique(o); % uses the actual participant numbers
            stderr = (std(PSD_netw(nt,:,:,st,se),0,3)/sqrt(length(PSD_netw(nt,:,:,st,se))));
            x = nanmean(PSD_netw(nt,:,:,st,se),3);
            figure; hold on
            h= zeros(length(y),1);
            shadedErrorBar([],x(1:40),stderr(1:40));
            h(1) = plot(x(1:40));
            title([v.netwname{nt} ' ' v.sess{se} ' ' v.stages{st}]);
            xlim([0 40]);
            if v.plot_outliers == 1
                cd(outlierdir);
                for iii = 1:length(y)
                    b = PSD_netw(nt,:,y(iii),st,se);
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
                saveas(gcf,['PSD_'  v.sess{se} '_' v.stages{st} '_' v.netwname{nt}  '_outliers.png']);
            elseif v.plot_outliers == 0
                saveas(gcf,['PSD_'  v.sess{se} '_' v.stages{st} '_' v.netwname{nt} '.png']);
            end
        end
        close all;
    end
end

end
