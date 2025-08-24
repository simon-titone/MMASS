function RS_ANOVA(NETW,design_matrix,n)
%%%%%% NETW(network, network, subject, run, band)
ANOVA = [n.homedir filesep 'results' filesep 'ANOVA_n' num2str(n.subj)];
if ~exist(ANOVA); mkdir(ANOVA); end; 
intersess_scatter = [n.homedir filesep 'results' filesep 'intersession_scatter_n' num2str(n.subj)];

F_netw_anova_mat=zeros(n.netws,n.netws,size(design_matrix.SessANOVA,2)); % 3 terms; session effec
pval_netw_anova_mat=ones(n.netws,n.netws,size(design_matrix.SessANOVA,2));

for b = 1:n.bands
    for w=1:n.netws
        for ww=1:n.netws
            ntw = NETW(w,ww,:,:,b);
            vect=squeeze(ntw(:));
            
            t = table(num2str(design_matrix.SessANOVA(1:n.subj)), vect(1:n.subj), vect(n.subj+1:n.subj*n.sess), 'VariableNames', {'Group', 'pre', 'post'});
            Sess = table([1 2]', 'VariableNames', {'Sessions'});
            rm = fitrm(t,'pre-post~1','WithinDesign', Sess);
            clear t; clear Sess; clear vect;
            [ranovatbl] = ranova(rm);  % Get within subject (session; session by group) effects
            pval_netw_anova_mat(w,ww) = ranovatbl{'(Intercept):Sessions','pValue'};
            F_netw_anova_mat(w,ww) = ranovatbl{'(Intercept):Sessions','F'};
            clear ranovatbl;
            clear rm;
        end
    end
    clear w; clear ww;
    
    for i=1:size(design_matrix.SessANOVA,2)
        
        count = 1;
        for nnetw_col = 1:1:n.netws
            pval.anovanetw(i, count) = pval_netw_anova_mat(1,nnetw_col,i);
            F_anovanetw(i,count) = F_netw_anova_mat(1,nnetw_col,i);
            count = count + 1;
        end
        clear count; clear nnetw_col;
        
        [~, dummy2, ~, dummy4]=fdr_bh(pval.anovanetw(i,:),n.fdr_thres,'pdep','no');
        crit_p.anovanetw(i) = dummy2;
        adj_p.anovanetw(i,1:length(dummy4)) = dummy4;
        clear dummy2; clear dummy4;
        index = find(pval.anovanetw(i,:) == crit_p.anovanetw(i));
        if isempty(index)
            crit_F = 1000; % place holder
        else
            crit_F = F_anovanetw(i,index(1));
        end
        clear index;
        
        mm=6;
        [ivect,jvect]=find(abs(F_netw_anova_mat(:,:,i))>=abs(crit_F));
        [ivect_unc,jvect_unc]=find(pval_netw_anova_mat(:,:,i)<=0.05);
        clear crit_F;
        
        figure(); set(gcf,'Color','white'); box OFF;
        imagesc(F_netw_anova_mat(:,:,i)); axis square; caxis([0 mm]); colormap(n.colors_onesided); colorbar; daspect([1 1 1]);
        title([ n.bandname{b} ' Timepoint Main Effect n= ' num2str(n.subj) ]); % generate corresponding rfx figure
    end
    hold on; for a2=1:n.netws-1  line([a2+0.5 a2+0.5],[0.5 n.netws+0.5],'color','k');   line([0.5 n.netws+0.5],[a2+0.5 a2+0.5],'color','k');    end
    clear a2;
    for b2=1:n.netws    line([b2+0.5 b2+0.5],[0.5 n.netws+0.5],'color','k');   line([0.5 n.netws+0.5],[b2+0.5 b2+0.5],'color','k');    end
    set(gca,'YTick',[1:n.netws],'YTickLabel',n.netwname,'Fontsize',16); set(gca,'XTick',[1:n.netws],'XTickLabel',n.netwname,'Fontsize',16);

    scatter(jvect_unc,ivect_unc,36,'ok');
    clear ivect_unc; clear jvect_unc;
    scatter(jvect,ivect,9,'ok','filled');
    clear ivect; clear jvect; clear a2; clear i; clear mm;
    cd(ANOVA);
    saveas(gcf, [n.bandname{b} 'netw_session_main_effect.png' ]);
    
    if n.scatter == 1
        if ~exist(intersess_scatter); mkdir(intersess_scatter); end; cd(intersess_scatter);
        for w=1:n.netws
            for ww=1:n.netws
                if pval_netw_anova_mat(w,ww) <0.052 == 1
                    x = squeeze(NETW(w,ww,:,1,b)); y = squeeze(NETW(w,ww,:,2,b));
                    figure; hold on; colormap summer;
                    plot([1 2], [x y])
                    h = bar([1 2], [mean(x) mean(y)], .5);
                    scatter(repmat(1,1,length(x)),x,'o', 'r');
                    labelpoints(1,x,n.subj_label,'E', .05);
                    scatter(repmat(2,1,length(y)),y,'o', 'b');
                    labelpoints(2,y,n.subj_label,'E', .05);
                    title([n.bandname{b} ' run 1 x run 2 ' n.netwname{w} '-' n.netwname{ww} ', n= ' num2str(n.subj)]);
                    ylabel('connectivity (z-score)'); xlabel('Run');
                    xticks([1 2]);
                    set(gcf, 'Position',  [10, 10, 300, 500]); 
%                     legend(strcat('f-value=', num2str(F_netw_anova_mat(w,ww))), strcat('p-value(uncorr)=  ', num2str(pval_netw_anova_mat(w,ww))), 'Location','northoutside');
                    text(0.4, 0.98,['f-value=' num2str(F_netw_anova_mat(w,ww))],'Units','normalized' );
                    text(0.4, 0.96,[ 'p-value(uncorr)=  ' num2str(pval_netw_anova_mat(w,ww))],'Units','normalized' );
                    saveas(gcf, [n.bandname{b} '_' n.netwname{w} '_' n.netwname{ww}   '_bar_scatter.png' ]);
                    clear('x','y');
                end
            end
        end
        close all
    end
end

close all
end