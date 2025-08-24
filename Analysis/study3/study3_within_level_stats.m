function sleep_within_level_stats(v,NETW)
%%%%% NETW: network, network, subject, stage, session, band

w_lvl =[v.homedir filesep 'results' filesep 'within_level_n' num2str(v.nsubj)];
if ~exist(w_lvl); mkdir(w_lvl); end; cd(w_lvl)

for se = 1:v.nsess
    for st = 1:v.nstages
        for b = 1:v.nbands
            
            for w=1:v.nnetws
                for ww=1:v.nnetws
                    vect = NETW(w,ww,:,st,se,b);
                    [~,~,STATS]=glmfit([design_matrix.withinlevel],vect,'normal','constant','off');
                    betanetws.level_mat(w,ww,:)=STATS.beta;
                    senetws.level_mat(w,ww,:)=STATS.se;
                end
            end
            clear w; clear ww; clear STATS; clear vect;
            
            for i=1:v.levels
                datanetws_ave=squeeze(betanetws.level_mat(:,:,i));
                datanetws_se=squeeze(senetws.level_mat(:,:,i));
                datanetws_tscore=datanetws_ave./datanetws_se;
                count = 1;
                for nnetw_lin = 1:1:v.nnetws
                    for nnetw_col = 1:1:v.nnetws
                        pval.withinLevelnetws(i, count) = 2*tcdf(abs(datanetws_tscore(nnetw_lin,nnetw_col)),v.nsubj-1, 'upper'); % 2 -tailed p value
                        t.withinLevelnetws(i,count) = datanetws_tscore(nnetw_lin,nnetw_col);
                        count = count + 1;
                    end
                end
                clear count;clear nnetw_col;clear nnetw_lin
                
                pval.withinLevelnetws_4cor = [pval.withinLevelnetws(:,1:v.nnetws) pval.withinLevelnetws(:,v.nnetws+2:v.nnetws*2) pval.withinLevelnetws(:,v.nnetws*2+3:v.nnetws*3) pval.withinLevelnetws(:,v.nnetws*3+4:v.nnetws*4) pval.withinLevelnetws(:,v.nnetws*4+5:v.nnetws*5) pval.withinLevelnetws(:,v.nnetws*5+6:v.nnetws*6)];
                t.withinLevelnetws_4cor = [t.withinLevelnetws(:,1:v.nnetws) t.withinLevelnetws(:,v.nnetws+2:v.nnetws*2) t.withinLevelnetws(:,v.nnetws*2+3:v.nnetws*3) t.withinLevelnetws(:,v.nnetws*3+4:v.nnetws*4) t.withinLevelnetws(:,v.nnetws*4+5:v.nnetws*5) t.withinLevelnetws(:,v.nnetws*5+6:v.nnetws*6)];
                
                % FDR correction
                [~, dummy2, ~, dummy4]=fdr_bh(pval.withinLevelnetws_4cor(i,:),fdr_thres,'pdep','no');
                crit_p.withinLevelnetws(i) = dummy2;
                adj_p.withinLevelnetws(i,1:length(dummy4)) = dummy4;
                clear dummy2; clear dummy4;
                
                index = find(pval.withinLevelnetws_4cor(i,:) == crit_p.withinLevelnetws(i));
                if ~isempty(index)
                    crit_t = t.withinLevelnetws_4cor(i,index(1));
                end
                crit_t_unc = tinv(fdr_thres/2,v.nsubj-1); % 2 tailed
                
                mm =10;
                if ~isempty(index)
                    [ivect,jvect]=find(abs(datanetws_tscore(:,:))>=abs(crit_t));
                end
                [ivect_unc,jvect_unc]=find(abs(datanetws_tscore(:,:))>=abs(crit_t_unc));
                
                figure(); set(gcf,'Color','white'); box OFF;
                imagesc(datanetws_tscore); axis square; caxis([-mm mm]); colormap(colors); colorbar; daspect([1 1 1]); 
                
                title([ v.bandname{b} ' ' v.sess{se} ' ' v.stages{st} ' within-level']); % generate corresponding rfx figure
                
                hold on; for a2=1:v.nnetws-1    line([a2+0.5 a2+0.5],[0.5 v.nnetws+0.5],'color','k');   line([0.5 v.nnetws+0.5],[a2+0.5 a2+0.5],'color','k');    end
                clear a2;
                scatter(jvect_unc,ivect_unc,36,'ok');
                clear ivect_unc; clear jvect_unc;
                if ~isempty(index)
                    scatter(jvect,ivect,9,'ok','filled');
                    clear ivect; clear jvect;
                end                
                for b2=1:v.nnetws    line([b2+0.5 b2+0.5],[0.5 v.nnetws+0.5],'color','k');   line([0.5 v.nnetws+0.5],[b2+0.5 b2+0.5],'color','k');    end
                set(gca,'YTick',[1:v.nnetws],'YTickLabel',netwname,'Fontsize',16);
                set(gca,'XTick',[1:v.nnetws],'XTickLabel',netwname,'Fontsize',16);
                                
                saveas(gcf, [v.bandname{b} ' ' v.sess{se} ' ' v.stages{st} '_within_stats.png' ]);
                
                clear datanetws_ave; clear datanetws_se; clear datanetws_tscore;
                clear crit_t; clear crit_t_unc; clear index;
            end
            clear a2; clear i; clear mm;
        end
    end
end
end
