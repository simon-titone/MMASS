function RS_within_level(NETW,design_matrix,n)

% % Within each level (i.e., group/session combo), fit z-transformed
% % correlation values with GLM. Use corresponding BETA and SE estimates to
% % compute a t-score. Write the BETA  estimates and t-scores to text files
% % (ffx and rfx, respectively) and generate corresponding plots

betaseeds.level_mat=zeros(1,n.seeds,n.levels);
betanetwks.level_mat=zeros(1,n.netws,n.levels);
seseeds.level_mat=zeros(1,n.seeds,n.levels);
senetwks.level_mat=zeros(1,n.netws,n.levels);
withinlevel = [n.homedir filesep 'results' filesep 'within_level_n' num2str(n.subj) ];
if ~exist(withinlevel); mkdir(withinlevel); end; cd(withinlevel);

for i= 1:n.levels
    for b = 1:n.bands
        for w=1:n.netws
            for ww=1:n.netws
                ntw = NETW(w,ww,:,:,b);
                vect=squeeze(ntw(:));
                [~,~,STATS]=glmfit([design_matrix.withinlevel],vect,'normal','constant','off');
                betanetws.level_mat(w,ww,:)=STATS.beta;
                senetws.level_mat(w,ww,:)=STATS.se;
            end
        end
        clear w; clear ww; clear STATS; clear vect;
        
        for i=1:n.levels

            datanetws_ave=squeeze(betanetws.level_mat(:,:,i));
            datanetws_se=squeeze(senetws.level_mat(:,:,i));
            datanetws_tscore=datanetws_ave./datanetws_se;
            
            count = 1;
            for nnetw_lin = 1:1:n.netws
                for nnetw_col = 1:1:n.netws
                    pval.withinLevelnetws(i, count) = 2*tcdf(abs(datanetws_tscore(nnetw_lin,nnetw_col)),n.subj-1, 'upper'); % 2 -tailed p value
                    %                 pval.withinLevelnetws(i, count) = 2*tcdf(abs(datanetws_tscore(nnetw_lin,nnetw_col)),n.subj-1); % 2 -tailed p value
                    t.withinLevelnetws(i,count) = datanetws_tscore(nnetw_lin,nnetw_col);
                    count = count + 1;
                end
            end
            clear count;clear nnetw_col;clear nnetw_lin
            
            pval.withinLevelnetws_4cor = [pval.withinLevelnetws(:,1:n.netws) pval.withinLevelnetws(:,n.netws+2:n.netws*2) pval.withinLevelnetws(:,n.netws*2+3:n.netws*3) pval.withinLevelnetws(:,n.netws*3+4:n.netws*4) pval.withinLevelnetws(:,n.netws*4+5:n.netws*5) pval.withinLevelnetws(:,n.netws*5+6:n.netws*6)];
            t.withinLevelnetws_4cor = [t.withinLevelnetws(:,1:n.netws) t.withinLevelnetws(:,n.netws+2:n.netws*2) t.withinLevelnetws(:,n.netws*2+3:n.netws*3) t.withinLevelnetws(:,n.netws*3+4:n.netws*4) t.withinLevelnetws(:,n.netws*4+5:n.netws*5) t.withinLevelnetws(:,n.netws*5+6:n.netws*6)];
            
            [~, dummy2, ~, dummy4]=fdr_bh(pval.withinLevelnetws_4cor(i,:),n.fdr_thres,'pdep','no');
            crit_p.withinLevelnetws(i) = dummy2;
            adj_p.withinLevelnetws(i,1:length(dummy4)) = dummy4;
            clear dummy2; clear dummy4;
            
            index = find(pval.withinLevelnetws_4cor(i,:) == crit_p.withinLevelnetws(i));
            crit_t = t.withinLevelnetws(i,index(1));
            crit_t_unc = tinv(n.fdr_thres/2,n.subj-1); % 2 tailed
            clear index;
            
            [ivect,jvect]=find(abs(datanetws_tscore(:,:))>=abs(crit_t));
            [ivect_unc,jvect_unc]=find(abs(datanetws_tscore(:,:))>=abs(crit_t_unc));
            
            figure(); set(gcf,'Color','white'); box OFF;
            imagesc(datanetws_tscore); axis square;  %caxis([-mm mm]);
            colormap(n.j_colors);
            colorbar; daspect([1 1 1]);
            title([ n.bandname{b} ' Ntw: within level' num2str(i) ' n= ' num2str(n.subj) ]); % generate corresponding rfx figure
            hold on; for a2=1:n.netws-1    line([a2+0.5 a2+0.5],[0.5 n.netws+0.5],'color','k');   line([0.5 n.netws+0.5],[a2+0.5 a2+0.5],'color','k');    end
            clear a2;
            scatter(jvect_unc,ivect_unc,36,'ok');
            clear ivect_unc; clear jvect_unc;
            
            scatter(jvect,ivect,9,'ok','filled');
            clear ivect; clear jvect;
            %         xticklabel_rotate([1:n.netws],90,netwname,'Fontsize',16);
            for b2=1:n.netws    line([b2+0.5 b2+0.5],[0.5 n.netws+0.5],'color','k');   line([0.5 n.netws+0.5],[b2+0.5 b2+0.5],'color','k');    end
            set(gca,'YTick',[1:n.netws],'YTickLabel',n.netwname,'Fontsize',16);
            set(gca,'XTick',[1:n.netws],'XTickLabel',n.netwname,'Fontsize',16);
            
            saveas(gcf, [n.bandname{b} '_netw_level' num2str(i) '_RFX.png' ]);

            clear dataseeds_ave; clear dataseeds_se; clear dataseeds_tscore;
            clear datanetws_ave; clear datanetws_se; clear datanetws_tscore;
            clear crit_t; clear crit_t_unc;
            
        end
    end
end

close all
clear a2; clear i; clear mm;
end