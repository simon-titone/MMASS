function RS_within_corr_bhv(NETW,speed,n)

withinlevel_corr = [n.homedir  'results' filesep 'within_level_corr_bhv_n' num2str(n.subj)];
if ~exist(withinlevel_corr); mkdir(withinlevel_corr); end;
withinlevel_bhv_scatter = [withinlevel_corr filesep 'within_level_corr_scatter'];
if ~exist(withinlevel_bhv_scatter); mkdir(withinlevel_bhv_scatter); end;

cov = NaN(n.netws,n.netws,n.levels,n.bands);
for bh=1:n.num_bhv_vars
    if bh==1
        behav_ByLevel = speed.offline;
    elseif bh==2
        behav_ByLevel = speed.online;
    end
    for i=1:n.levels
        for b = 1:n.bands

            cov=NETW(:,:,:,i,b);
            for w=1:n.netws
                for ww=1:n.netws
                    val=zeros(1);
                    vect=not(isnan(behav_ByLevel));
                    if sum(vect)>1
                        [val,px]=corr(squeeze(cov(w,ww,:)),behav_ByLevel(vect)');
                    end
                    rval(w,ww,:)=val;
                    pval(w,ww,:)=px;
                end
            end
            clear val; clear px; clear w; clear ww; clear vect;

            if n.FDR == 'mot'
                pval_FDR  = pval(5,:);
            elseif n.FDR == 'all'
                pval_FDR = [ pval(1:n.netws)  pval(n.netws+2:n.netws*2)  pval(n.netws*2+3:n.netws*3)  pval(n.netws*3+4:n.netws*4)  pval(n.netws*4+5:n.netws*5)  pval(n.netws*5+6:n.netws*6)];
            end

            % FDR correction
            [~, dummy2, ~, dummy4]= fdr_bh(pval_FDR,n.fdr_thres,'pdep','no');
            crit_p.Behavcov(i) = dummy2;
            adj_p.Behavcov(1:length(dummy4)) = dummy4;
            clear dummy2; clear dummy4;

            if i == 1
                sess = 'Pre-Task';
            elseif i ==2
                sess = 'Post-Task';
            end
            cd(withinlevel_bhv_scatter);

            mm=0.8;
            [ivect,jvect]=find(pval(:,:)<=crit_p.Behavcov(i));
            [ivect_unc,jvect_unc]=find(pval(:,:)<=n.fdr_thres);
            figure(); set(gcf,'Color','white'); box OFF;

            matrix = rval(:,:);
            tri_matrix = triu(matrix);

            imagesc(tri_matrix); axis square; caxis([-mm mm]); hold on;
            colormap('redblue'); daspect([1 1 1]);
            h = colorbar; set(get(h,'title'),'string','r-val', 'fontsize', 12); %ylabel(h, 'f-val');
            dim = set(gcf, 'Position',  [1000, 900, 500, 450]); %left bottom width height

            if n.title == 1
                title([ n.bandname{b} ' ' sess ' x ' n.bhv_vars{bh} ', {\it n} = ' num2str(n.subj) ]); % generate corresponding rfx figure
            end

            for a2=1:n.netws-1
                line([a2+0.5 a2+0.5],[a2+1.5 0.5],'color','k','LineWidth',0.8);
                line([a2-0.5 n.netws+0.5],[a2+0.5 a2+0.5],'color','k','LineWidth',0.8);
            end
            clear a2;

            for b2=1:n.netws
                line([b2+0.5 b2+0.5],[0.5 n.netws+0.5],'color','k','LineWidth',0.001);
                line([0.5 n.netws+0.5],[b2+0.5 b2+0.5],'color','k','LineWidth',0.001);
            end

            set(gca,'YTick',[1:n.netws],'YTickLabel',n.netwname,'Fontsize',16);
            set(gca,'XTick',[1:n.netws],'XTickLabel',n.netwname,'Fontsize',16);

            scatter(jvect_unc,ivect_unc, 80,'w', 'o' ,'ok','LineWidth',2.5);
            clear ivect_unc; clear jvect_unc;
            scatter(jvect,ivect,80,'w','ok','filled' ,'LineWidth',2.5);
            clear ivect; clear jvect;

            cd(withinlevel_corr);
            saveas(gcf, [n.bandname{b} '_' sess '_corr_' n.bhv_vars{bh} '.png' ]);

            %%create scatter plots (added 24/8/20 ST)
            reset(gcf); reset(gca);
            for w=1:n.netws
                for ww=1:n.netws
                    if pval(w,ww) <0.051 == 1
                        withinlevel_bhv_scatter_band = [withinlevel_bhv_scatter filesep n.bandname{b}];
                        if ~exist(withinlevel_bhv_scatter_band); mkdir(withinlevel_bhv_scatter_band); end;

                        figure
                        x = squeeze(cov(w,ww,:)); y = behav_ByLevel';
                        scatter(x,y,[], 'b', 'filled');lsline;
                        [p,~] = polyfit(x,y,1);

                        mdl = fitlm(x, y, 'linear');
                        [y_pred,y_ci] = predict(mdl,x, 'Alpha', 0.05);
                        Y = [x, y_ci,y_pred];
                        Y = sortrows(Y,1);
                        eb = [Y(:,4)-Y(:,2)];
                        boundedline(Y(:,1),Y(:,4),eb, '-r','alpha');
                        clear mdl; clear y_pred; clear y_ci; clear Y; clear eb;
                        if n.title == 1
                            title({[n.bandname{b} ' ' sess ' ' n.netwname{w} '-' n.netwname{ww} ] [' Connectivity x '  n.bhv_vars{bh} ]});
                            xlabel('connectivity (z-score)');ylabel([n.bhv_vars{bh} ' gains speed (s)']);
                        end
                        set(gcf, 'Position',  [100, 100, 480, 320]); %left bottom width height
                        a = get(gca,'XTickLabel');
                        set(gca,'XTickLabel',a,'fontsize',14) %,'FontWeight','bold'
                        if n.pval_labels == 1
                            text(0.6, 0.98,['r-value= ' num2str(round(rval(w,ww), 3))],'Units','normalized' ,'FontSize',14);
                            text(0.6, 0.93,['p-value=  ' num2str(round(pval(w,ww),3))],'Units','normalized' ,'FontSize',14);
                            %                         text(0.81, 0.05,['*uncorrected'],'Units','normalized' ,'FontSize',12);
                        end
                        if n.subj_labels == 1
                            labelpoints(x,y,n.subj_label,'E', .05);
                        end
                        scatter_name = [ withinlevel_bhv_scatter_band filesep  n.bhv_vars{bh} '_'  sess '_' n.bandname{b} '_' n.netwname{w} '_' n.netwname{ww}   '_conn_scatter.png' ];
                        saveas(gcf, scatter_name);
                        clear('x','y');
                    end
                end
            end
        end
        close all;
    end
end
close all
clear CC_zscore_cov_ByLevel; clear behav_ByLevel; clear pval_Behav;
clear i;
