function sleep_btw_corr(v,NETW)
%%%%% NETW: network, network, subject, stage, session, band

btw_session_corr_dir = [v.homedir filesep 'results' filesep 'btw_session_corr_bhv_' filesep 'btw_session_corr_bhv_'];

% calculates difference between sessions, per stage (e.g. EXP NREM2 - CTRL NREM2)
for st= 1:v.nstages
    for b = 1:v.nbands
        for i = 1:1:v.nsubj
            for w=1:v.nnetws
                for ww = 1:v.nnetws
                    intersess_conn(w,ww,i,st,b) = NETW(w,ww,i,st,2,b) - NETW(w,ww,i,st,1,b);
                end
            end
        end
    end
end

% calculates difference between stages, within session (e.g. EXP NREM2 - EXP NREM3)
for se = 1:v.nsess
    for b = 1:v.nbands
        for i = 1:1:v.nsubj
            for w=1:v.nnetws
                for ww = 1:v.nnetws
                    interstage_conn(w,ww,i,se,1,b) = NETW(w,ww,i,1,se,b) - NETW(w,ww,i,2,se,b);
                    interstage_conn(w,ww,i,se,2,b) = NETW(w,ww,i,1,se,b) - NETW(w,ww,i,3,se,b);
                    interstage_conn(w,ww,i,se,3,b) = NETW(w,ww,i,2,se,b) - NETW(w,ww,i,3,se,b);
                end
            end
        end
    end
end

if v.group == 'offline'
    z = 1;
else
    z = 2:3;
end
%% Interstage correlation with behavior
for bh= z
    behavior = v.speed(bh,:);
    for se = 1:v.nsess
        for s= 1:v.nstages
            btw_sess_corr_sess_dir =[btw_session_corr_dir v.stages{s}];
            if ~exist(btw_sess_corr_sess_dir); mkdir(btw_sess_corr_sess_dir); end; cd(btw_sess_corr_sess_dir)
            for b = 1:v.nbands
                for w=1:v.nnetws
                    for ww=1:v.nnetws
                        [val,px]=corr(squeeze(interstage_conn(w,ww,:,se,st,b)),behavior');
                        r_val(w,ww,se,st,b)=val;
                        p_val(w,ww,se,st,b)=px;
                    end
                end
            end
            
            clear val; clear px; clear w; clear ww; clear vect;
            for b=1:v.nbands
                            pval_fdr = p_val(:,:,se,st,b);

                if v.FDR == 'mot'
                    p_val_temp_FDR_correction  = pval_fdr(5,:);
                elseif v.FDR == 'all'
                    p_val_temp_FDR_correction = [ pval_fdr(1:v.nnetw)  pval_fdr(v.nnetw+2:v.nnetw*2)  pval_fdr(v.nnetw*2+3:v.nnetw*3)  pval_fdr(v.nnetw*3+4:v.nnetw*4)  pval_fdr(v.nnetw*4+5:v.nnetw*5)  pval_fdr(v.nnetw*5+6:v.nnetw*6)];
                end
                [~, dummy2, ~, dummy4]=fdr_bh(p_val_temp_FDR_correction, v.fdr_thres,'pdep','no');
                crit_p.Behavcov = dummy2;
                adj_p.Behavcov(1:length(dummy4)) = dummy4;
                clear dummy2; clear dummy4;
                
                
                mm=0.8;
                [ivect,jvect]=find(pmat_behav.cov_InterSession.bandtemp{b}(:,:)<=crit_p.Behavcov);
                [ivect_unc,jvect_unc]=find(pmat_behav.cov_InterSession.bandtemp{b}(:,:)<=fdr_thres);
                figure(); set(gcf,'Color','white'); box OFF;
                imagesc(corrmat_behav.cov_InterSession.bandtemp{b}(:,:)); axis square; caxis([-mm mm]); colormap(colors); colorbar; daspect([1 1 1]);
                
                title([ bandname{b} ' Netw: ' stages{s} ' Corr w/ Factor ' bhv_vars{bh} ]); % generate corresponding rfx figure
                
                hold on; for a2=1:v.nnetws-1    line([a2+0.5 a2+0.5],[0.5 v.nnetws+0.5],'color','k');   line([0.5 v.nnetws+0.5],[a2+0.5 a2+0.5],'color','k');    end
                clear a2;
                %         set(gca,'YTick',[1],'YTickLabel',Label.Seed,'Fontsize',16);
                for b2=1:v.nnetws    line([b2+0.5 b2+0.5],[0.5 v.nnetws+0.5],'color','k');   line([0.5 v.nnetws+0.5],[b2+0.5 b2+0.5],'color','k');    end
                set(gca,'YTick',[1:v.nnetws],'YTickLabel',netwname,'Fontsize',16);
                set(gca,'XTick',[1:v.nnetws],'XTickLabel',netwname,'Fontsize',16);
                scatter(jvect_unc,ivect_unc,36,'ok');
                clear ivect_unc; clear jvect_unc;
                scatter(jvect,ivect,9,'ok','filled');
                clear ivect; clear jvect;
                %         xticklabel_rotate([1:n.roi_MSL],90,Label.MSL,'Fontsize',16);
                %         print([homedir '/matrix/' 'Level' num2str(i) '_MSL_Corr_Factor ' num2str(ww) '.tif'],'-dtiff','-r150');
                
                saveas(gcf, [bandname{b} '_' stages{s} '_corr_' bhv_vars{bh} '.png' ]);
                
                
                % Scatter plots
                for w=1:v.nnetws
                    for ww=1:v.nnetws
                        if p_val(w,ww) <0.05 == 1
                            figure
                            x = squeeze(interstage_conn(w,ww,:,se,st,b)); y = behavior';
                            scatter(x,y,[], 'b', 'filled');lsline;
                            mdl = fitlm(x, y, 'linear');
                            [y_pred,y_ci] = predict(mdl,x, 'Alpha', 0.05);
                            Y = [x, y_ci,y_pred];
                            Y = sortrows(Y,1);
                            eb = [Y(:,4)-Y(:,2)];
                            boundedline(Y(:,1),Y(:,4),eb, '-r','alpha');
                            clear mdl; clear y_pred; clear y_ci; clear Y; clear eb;
                            
                            title([v.bandname{b} ' ' v.sess{se} ' ' v.stages{st} ' ' v.netwname{w} '-' v.netwname{ww} ' connectivity x '  v.bhv_vars{bh} ]);xlabel('connectivity (z-score)');ylabel([v.bhv_vars{bh} ' gains speed']);
                            set(gcf, 'Position',  [100, 100, 500, 300]); %left bottom width height
                            legend(strcat('r-value=', num2str(r_val(w,ww))), strcat('p-value(uncorr)=  ', num2str(p_val(w,ww))), 'Location','northoutside');
                            labelpoints(x,y,v.subj_label,'E', .05);
                            saveas(gcf, [v.bandname{b} '_' v.sess{se} '_' v.stages{st} '_' v.bhv_vars{bh} '_' v.netwname{w} '_' v.netwname{ww}   '_scatter.png' ]);
                            clear('x','y');
                        end
                    end
                end
                
                clear CC_zscore_cov_InterSession; clear pval_Behav;
                clear i;
            end
        end
    end
    
    close all
    
    % Between-stages correlations to behavior
    CC_zscore_cov_InterStage = NaN(v.nnetws,v.nnetws,v.subj,v.nbands);
    
    for w=1:v.nnetws
        for ww=1:v.nnetws
            for b = 1:v.nbands
                val=zeros(1);
                vect=not(isnan(behavior));
                if sum(vect)>1
                    [val,px]=corr(squeeze(CC_zscore_cov_InterStage(w,ww,:,b)),behavior(vect)');
                end
                corrmat_behav.cov_InterStage.bandtemp{b}(w,ww,:)=val;
                pmat_behav.cov_InterStage.bandtemp{b}(w,ww,:)=px;
            end
        end
    end
    clear val; clear px; clear w; clear ww; clear vect;
    
    for b=1:v.nbands
        pval.behav.cov_InterStage_4cor.bandtemp{b} = [ pmat_behav.cov_InterStage.bandtemp{b}(1:v.nnetws)  pmat_behav.cov_InterStage.bandtemp{b}(v.nnetws+2:v.nnetws*2)  pmat_behav.cov_InterStage.bandtemp{b}(v.nnetws*2+3:v.nnetws*3)  pmat_behav.cov_InterStage.bandtemp{b}(v.nnetws*3+4:v.nnetws*4)  pmat_behav.cov_InterStage.bandtemp{b}(v.nnetws*4+5:v.nnetws*5)  pmat_behav.cov_InterStage.bandtemp{b}(v.nnetws*5+6:v.nnetws*6)];
        % for ww=1:n.behav
        %     count = 1;
        %     for nMSL_col = 1:1:n.roi_MSL
        %         pval_Behav(w,ww) = pmat_behav.cov_WithinLevel(w,ww);
        %         count = count + 1;
        %     end
        % end
        [~, dummy2, ~, dummy4]=fdr_bh(pval.behav.cov_InterStage_4cor.bandtemp{b},fdr_thres,'pdep','no');
        crit_p.Behavcov = dummy2;
        adj_p.Behavcov(1:length(dummy4)) = dummy4;
        clear dummy2; clear dummy4;
        
        mm=0.8;
        [ivect,jvect]=find(pmat_behav.cov_InterStage.bandtemp{b}(:,:)<=crit_p.Behavcov);
        [ivect_unc,jvect_unc]=find(pmat_behav.cov_InterStage.bandtemp{b}(:,:)<=fdr_thres);
        figure(); set(gcf,'Color','white'); box OFF;
        imagesc(corrmat_behav.cov_InterStage.bandtemp{b}(:,:)); axis square; caxis([-mm mm]); colormap(colors); colorbar; daspect([1 1 1]);
        
        title([ bandname{b} ' Netw: ' stagecorr{s} ' Corr w/ Factor ' bhv_vars{bh} ]); % generate corresponding rfx figure
        
        hold on; for a2=1:v.nnetws-1    line([a2+0.5 a2+0.5],[0.5 v.nnetws+0.5],'color','k');   line([0.5 v.nnetws+0.5],[a2+0.5 a2+0.5],'color','k');    end
        clear a2;
        for b2=1:v.nnetws    line([b2+0.5 b2+0.5],[0.5 v.nnetws+0.5],'color','k');   line([0.5 v.nnetws+0.5],[b2+0.5 b2+0.5],'color','k');    end
        set(gca,'YTick',[1:v.nnetws],'YTickLabel',netwname,'Fontsize',16);
        set(gca,'XTick',[1:v.nnetws],'XTickLabel',netwname,'Fontsize',16);
        scatter(jvect_unc,ivect_unc,36,'ok');
        clear ivect_unc; clear jvect_unc;
        scatter(jvect,ivect,9,'ok','filled');
        clear ivect; clear jvect;
        
        saveas(gcf, [bandname{b} '_' stagecorr{s} '_corr_' bhv_vars{bh} '.png' ]);
        
        clear CC_zscore_cov_InterStage; clear pval_Behav;
        clear i;
    end
end

close all
end