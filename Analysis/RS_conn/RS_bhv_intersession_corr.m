function RS_bhv_intersession_corr(NETW,speed,n)

btwsess_bhv =[n.homedir filesep 'results' filesep 'intersession_corr_n' num2str(n.subj)];
if ~exist(btwsess_bhv); mkdir(btwsess_bhv);end; 
btwsess_bhv_scatter = [btwsess_bhv filesep 'intersession_corr_scatter'];
if ~exist(btwsess_bhv_scatter); mkdir(btwsess_bhv_scatter); end;
%%%%%% NETW(network, network, subject, run, band)

InterSession = NaN(n.netws,n.netws,n.subj,n.bands);
for b = 1:n.bands
    for i = 1:1:n.subj
        for w=1:n.netws
            for ww = 1:n.netws
                InterSession(w,ww,i,b) = NETW(w,ww,i,2,b) - NETW(w,ww,i,1,b);
            end
        end
    end
end

for bh=1:n.num_bhv_vars
    if bh==1
        behav_ByLevel = speed.offline;
    elseif bh==2
        behav_ByLevel = speed.online;
        %     elseif bh==3
        %         behav_ByLevel = speed.mean;
    end
    
    clear i; clear w; clear ww;
    for w=1:n.netws
        for ww=1:n.netws
            for b = 1:n.bands
                val=zeros(1);
                vect=not(isnan(behav_ByLevel));
                if sum(vect)>1
                    [val,px]=corr(squeeze( InterSession(w,ww,:,b)),behav_ByLevel(vect)');
                end
                intersession_r(w,ww,b,bh)=val;
                intersession_pval(w,ww,b,bh)=px;
            end
        end
    end
    clear val; clear px; clear w; clear ww; clear vect;
    %     intersession_pval = squeeze(intersession_pval); intersession_r = squeeze(intersession_r);
%     for b=1:n.bands
    for b=3
        intersess_pval_bandtemp = intersession_pval(:,:,b,bh);
        % pval.behav.cov_InterSession_4cor.bandtemp{b} = [intersess_pval_bandtemp(1:n.netws)  intersess_pval_bandtemp(n.netws+2:n.netws*2) intersess_pval_bandtemp(n.netws*2+3:n.netws*3)  intersess_pval_bandtemp(n.netws*3+4:n.netws*4)  intersess_pval_bandtemp(n.netws*4+5:n.netws*5)  intersess_pval_bandtemp(n.netws*5+6:n.netws*6)];
        if n.FDR == 'mot'
            intersess_pval_vect = [intersess_pval_bandtemp(5, 1:n.netws)];
        elseif n.FDR == 'all'
            intersess_pval_vect = [intersess_pval_bandtemp(1:n.netws)  intersess_pval_bandtemp(n.netws+2:n.netws*2) intersess_pval_bandtemp(n.netws*2+3:n.netws*3)  intersess_pval_bandtemp(n.netws*3+4:n.netws*4)  intersess_pval_bandtemp(n.netws*4+5:n.netws*5)  intersess_pval_bandtemp(n.netws*5+6:n.netws*6)];
        end
        
        %%create scatter plots (added 24/8/20 ST)
        reset(gcf); reset(gca);
        cd(btwsess_bhv_scatter);
        for w=1:n.netws
            for ww=1:n.netws
%                 if intersess_pval_bandtemp(w,ww) <0.05 == 1
                    figure
                    x = squeeze(InterSession(w,ww,:,b)); y = behav_ByLevel';
                    scatter(x,y,[], 'b', 'filled');lsline;
                    
                    mdl = fitlm(x, y, 'linear');
                    [y_pred,y_ci] = predict(mdl,x, 'Alpha', 0.05);
                    Y = [x, y_ci,y_pred];
                    Y = sortrows(Y,1);
                    eb = [Y(:,4)-Y(:,2)];
                    boundedline(Y(:,1),Y(:,4),eb, '-r','alpha');
                    clear mdl; clear y_pred; clear y_ci; clear Y; clear eb;
                    
                    title({[n.bandname{b} ' Intersession ' n.netwname{w} '-' n.netwname{ww} ] [' Connectivity x '  n.bhv_vars{bh} ]});
                    xlabel('connectivity (z-score)');ylabel([n.bhv_vars{bh} ' gains speed (s)']);
                    set(gcf, 'Position',  [100, 100, 480, 320]); %left bottom width height
                    a = get(gca,'XTickLabel');
                    set(gca,'XTickLabel',a,'fontsize',14) %,'FontWeight','bold'
                    %                     legend(strcat('r-value=', num2str(intersession_r(w,ww,b,bh))), strcat('p-value(uncorr)=  ', num2str(intersess_pval_bandtemp(w,ww))), 'Location','northoutside');
                    if n.pval_labels == 1
                        text(0.6, 0.98,['r-value= ' num2str(round(intersession_r(w,ww,b,bh), 3))],'Units','normalized' ,'FontSize',14);
                        text(0.6, 0.93,['p-value*=  ' num2str(round(intersess_pval_bandtemp(w,ww),3))],'Units','normalized' ,'FontSize',14);
                        text(0.81, 0.05,['*uncorrected'],'Units','normalized' ,'FontSize',12);
                    end
                    if n.subj_labels == 1
                        labelpoints(x,y,n.subj_label,'E', .05);
                    end
                    saveas(gcf, [n.bandname{b} '_' n.bhv_vars{bh} '_' n.netwname{w} '_' n.netwname{ww} '_intersess_scatter.png' ]);
                    clear('x','y');
%                 end
            end
        end
        
        [~, dummy2, ~, dummy4]=fdr_bh(intersess_pval_vect,n.fdr_thres,'pdep','no');
        crit_p.Behavcov = dummy2;
        adj_p.Behavcov(1:length(dummy4)) = dummy4;
        clear dummy2; clear dummy4;
        
        mm=0.8;
        [ivect,jvect]=find(intersession_pval(:,:,b)<=crit_p.Behavcov);
        [ivect_unc,jvect_unc]=find(intersess_pval_bandtemp(:,:)<=n.fdr_thres);
        figure(); set(gcf,'Color','white'); box OFF;
        dim = set(gcf, 'Position',  [1000, 900, 500, 450]); %left bottom width height

        matrix = intersession_r(:,:,b,bh);
        tri_matrix = triu(matrix);
        
        imagesc(tri_matrix); axis square; caxis([-mm mm]); colormap('redblue');  daspect([1 1 1]);
        h = colorbar; set(get(h,'title'),'string','r-val', 'fontsize', 12); %ylabel(h, 'f-val');
        
        title({[ n.bandname{b} ' Inter-Session Changes'] [ 'Correlated w/ ' n.bhv_vars{bh} ', {\it n} = ' num2str(n.subj)]}); % generate corresponding rfx figure
        
        hold on; 

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
        
        scatter(jvect_unc,ivect_unc, 80, 'w', 'o' ,'ok','LineWidth',2.5);
        clear ivect_unc; clear jvect_unc;
        scatter(jvect,ivect,80,'w','ok','filled' ,'LineWidth',2.5);
        clear ivect; clear jvect; clear a2; clear i; clear mm;
        
        saveas(gcf, [btwsess_bhv filesep n.bandname{b} '_intersession_x_' n.bhv_vars{bh}  '.png' ]);
        
        clear pval_Behav;
        clear i; clear dim
    end
    close all;
end
end