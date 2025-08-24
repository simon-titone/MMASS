function study2_interstage_ANOVA(v,NETW,design_matrix)
%%%% NETW: network, network, subject, stage, cycle, band
study2_anova_dir = [v.homedir filesep 'results' filesep 'ANOVA' ];
study2_anova_dir_stage = [study2_anova_dir filesep 'interstage'  filesep v.t_group];

%% Inter-Stage ANOVA
for c = v.cycle
    for st= v.stage
        if ~exist(study2_anova_dir_stage); mkdir(study2_anova_dir_stage); end
        if st ==1
            st1 = 1; st2 = 2; % compares stage 1 to stage 2
            st1_title = 'NREM2'; st2_title = 'NREM3';
        elseif st ==2
            st1 = 1; st2 = 3; % compares stage 1 to stage 3
            st1_title = 'NREM2'; st2_title = 'REM';
        elseif st ==3
            st1 = 2; st2 = 3; % compares stage 2 to stage 3
            st1_title = 'NREM3'; st2_title = 'REM';
        end
        
        for b = 1:v.nbands
            cd(study2_anova_dir_stage)
            for w=1:v.nnetws
                for ww=1:v.nnetws
                    
                    vect1 = squeeze(NETW(w,ww,:,st1,c,b));
                    vect2 = squeeze(NETW(w,ww,:,st2,c,b));
                    t = table(num2str(design_matrix.SessANOVA(1:v.nsubj)), vect1, vect2, 'VariableNames', {'Group', 'stage1', 'stage2'});
                    Sess = table([1 2]', 'VariableNames', {'Sessions'});
                    rm = fitrm(t,'stage1-stage2~1','WithinDesign', Sess);
                    clear t; clear Sess; clear vect2;
                    
                    [ranovatbl] = ranova(rm);  % Get within subject (session; session by group) effects
                    ANOVA_p(w,ww,st,c,b) = ranovatbl{'(Intercept):Sessions','pValue'};
                    ANOVA_F(w,ww,st,c,b) = ranovatbl{'(Intercept):Sessions','F'};
                    
                    clear ranovatbl;
                    clear rm;
                end
            end
            clear w; clear ww;
            i=1;
            for ntw = 1:1:v.nnetws
                    pval_vect(:,:,st,c,b) = [ANOVA_p(1,1:v.nnetws,st,c,b) ANOVA_p(2,2:v.nnetws,st,c,b) ANOVA_p(3,3:v.nnetws,st,c,b) ANOVA_p(4,4:v.nnetws,st,c,b) ANOVA_p(5,5:v.nnetws,st,c,b)  ANOVA_p(6,6:v.nnetws,st,c,b)];
                    F_anovanetw(:,:,st,c,b) = [ANOVA_F(1,1:v.nnetws,st,c,b) ANOVA_F(2,2:v.nnetws,st,c,b) ANOVA_F(3,3:v.nnetws,st,c,b) ANOVA_F(4,4:v.nnetws,st,c,b) ANOVA_F(5,5:v.nnetws,st,c,b)  ANOVA_F(6,6:v.nnetws,st,c,b)];
                    ntw = ntw + 1;
            end
            clear count; clear nnetw_col;
            
            [~, dummy2, ~, dummy4]=fdr_bh(pval_vect(:,:,st,c,b),v.fdr_thres,'pdep','no');
            crit_p(:,:,st,c,b) = dummy2;
            adj_p(i,1:length(dummy4),st,c,b) = dummy4;
            clear dummy2; clear dummy4;
            index = find(pval_vect(i,:,st,c,b) == crit_p);

            adj_p_mat(1,1:6,st,c,b) = adj_p(1,1:6,st,c,b); adj_p_mat(2,2:6,st,c,b) = adj_p(1,7:11,st,c,b);
            adj_p_mat(3,3:6,st,c,b) = adj_p(1,12:15,st,c,b); adj_p_mat(4,4:6,st,c,b) = adj_p(1,16:18,st,c,b);
            adj_p_mat(5,5:6,st,c,b) = adj_p(1,19:20,st,c,b); adj_p_mat(6,6,st,c,b) = adj_p(1,21,st,c,b);

                  
            if isempty(index); crit_F = 1000; % place holder
            else; crit_F = F_anovanetw(index);
            end
            clear index;
            
            mm=6;
            [ivect,jvect]=find(abs(ANOVA_F(:,:,st,c,b))>=abs(crit_F));
            [ivect_unc,jvect_unc]=find(ANOVA_p(:,:,st,c,b)<=0.05);
            clear crit_F;
            
            figure(); set(gcf,'Color','white'); box OFF;
            
            matrix = ANOVA_F(:,:,st,c,b);
            tri_matrix = triu(matrix);

            colormap(flipud(pink(10)));
             
            imagesc(tri_matrix);axis square;  daspect([1 1 1]); 
            ANOVA_max=max(max(ANOVA_F(:,:,st,c,b)));
            if gt(ANOVA_max, 6) ==1 && lt(ANOVA_max, 15) ==1;mm = ANOVA_max; elseif gt(ANOVA_max, 15) ==1 ; mm = 15; else; mm= 6;end;caxis([0 mm]);
            
            set(gcf,'Color','white');set(gca,'YDir','reverse');hold on;
            dim = set(gcf, 'Position',  [1200, 900, 500, 450]); %left bottom width height
            h = colorbar; set(get(h,'title'),'string','f-val', 'fontsize', 12);
            if  v.fig_labels == 1
                title({[ v.bandname{b}  ' cycle ' num2str(v.cycle) ] [ v.stagediff_title{st}  ' ANOVA, n=' num2str(v.nsubj)]}); % generate corresponding rfx figure
            else
                title({[ v.bandname{b}  ' cycle ' num2str(v.cycle) ] [ v.stagediff_title{st}  ' ANOVA']}); % generate corresponding rfx figure
            end
            for a2=1:v.nnetws-1; line([a2+0.5 a2+0.5],[a2+1.5 0.5],'color','k','LineWidth',0.8);line([a2-0.5 v.nnetws+0.5],[a2+0.5 a2+0.5],'color','k','LineWidth',0.8); end            
            for b2=1:v.nnetws; line([b2-0.5 b2-0.5],[0.5 v.nnetws+0.5],'color','k','LineWidth',0.001); line([0.5 v.nnetws+0.5],[b2-0.5 b2-0.5],'color','k','LineWidth',0.001);end
            clear a2; clear b2;
            
            set(gca,'YTick',[1:v.nnetws],'YTickLabel', v.netwname,'Fontsize',16);
            set(gca,'XTick',[1:v.nnetws],'XTickLabel', v.netwname,'Fontsize',16);
            
            scatter(jvect_unc,ivect_unc, 80, 'w', 'o' ,'ok','LineWidth',2.5); clear ivect_unc; clear jvect_unc;
            scatter(jvect,ivect,80,'w','ok','filled' ,'LineWidth',2.5); clear ivect; clear jvect;

            saveas(gcf, [study2_anova_dir_stage filesep 'c' num2str(v.cycle) '_' v.stagediff{st} '_' v.bandname{b} '_ANOVA.png'])
            clear a2; clear i; clear mm;

            if v.scatter == 1
                intersess_scatter = [study2_anova_dir_stage filesep 'scatter' filesep v.bandname{b}];
                if ~exist(intersess_scatter); mkdir(intersess_scatter); end; 
                for w=1:v.nnetws
                    for ww=1:v.nnetws
                        if ANOVA_p(w,ww,st,c,b) <0.05 == 1
                            x = squeeze(NETW(w,ww,:,st1,c,b)); y = squeeze(NETW(w,ww,:,st2,c,b));
                            F = figure; hold on; colormap summer;
                            plot([1 2], [x y])
                            h = bar([1 2], [nanmean(x) nanmean(y)], .5);
                            scatter(repmat(1,1,length(x)),x,'o', 'r');
                            labelpoints(1,x,v.subj_label,'E', .05);
                            scatter(repmat(2,1,length(y)),y,'o', 'b');
                            labelpoints(2,y,v.subj_label,'E', .05);
                            if v.fig_title == 1
                            title(['Cycle ' num2str(v.cycle) ' ' v.bandname{b} ' ' v.netwname{w} '-' v.netwname{ww} ]);
                            end
                            ylabel('connectivity (z-score)');
                            grouporder={st1_title,st2_title};
                            xticks([1 2]); xticklabels(grouporder);
                            set(gcf, 'Position',  [10, 10, 300, 500]);
                            F = round(ANOVA_F(w,ww,st,c,b),3);
                            p = round(ANOVA_p(w,ww,st,c,b),3);
                            if v.fig_labels == 1
                                text(0.4, 0.98,['F = ' num2str(F)] ,'Units','normalized');
                                if p > 0.001 == 1
                                    text(0.4, 0.95,['p(UC) = ' num2str(p)],'Units','normalized');
                                else
                                    text(0.4, 0.96,[ 'p(UC) < 0.001'],'Units','normalized' );
                                end
                            end
                            scatter_savename = [ intersess_scatter  filesep 'c' num2str(v.cycle) '_' v.stagediff{st} '_' v.bandname{b} '_' v.netwname{w} '_' v.netwname{ww}   '_scatter.png' ];
                            saveas(gcf, scatter_savename);
                            clear('x', 'y', 'f');
                        end
                    end
                end
                close all
            end
            
            if v.violin == 1
            
                for w=1:v.nnetws
                    for ww=1:v.nnetws
                        if ANOVA_p(w,ww,st,c,b) <0.05 == 1
                            ANOVA_violin = [study2_anova_dir_stage filesep 'violin' filesep v.bandname{b}];
                            if ~exist(ANOVA_violin); mkdir(ANOVA_violin); end; cd(ANOVA_violin);
                            x = [squeeze(NETW(w,ww,:,st1,c,b)); squeeze(NETW(w,ww,:,st2,c,b))];
                            y = [repmat({[st1_title]}, v.nsubj,1) ; repmat({[st2_title]}, v.nsubj,1)];
                            z = [repmat(1,v.nsubj,1) ; repmat(2, v.nsubj,1)];
                            grouporder={st1_title,st2_title};
                            figure; v1 = violinplot(x, y,'GroupOrder',grouporder, 'ShowMean', true);
                            ylabel('connectivity (z-score)');
                            a = get(gca,'XTickLabel');
                            set(gca,'XTickLabel',a,'fontsize',16)
                            set(gcf, 'Position',  [10, 10, 350, 500]);
                            if v.fig_title == 1
                                title(['Cycle ' num2str(v.cycle) ' ' v.bandname{b} ' ' v.netwname{w} '-' v.netwname{ww} ], 'fontsize', 16);
                            end
                            if v.fig_labels == 1
                                F = num2str(round(ANOVA_F(w,ww,st,c,b),3)); p = round(ANOVA_p(w,ww,st,c,b),4);
                                text(0.4, 0.98,['F =' F],'Units','normalized' );
                                if p > 0.001 == 1
                                    text(0.4, 0.96,[ 'p(UC) =  ' num2str(p)],'Units','normalized' );
                                else
                                    text(0.4, 0.96,[ 'p(UC) < 0.001'],'Units','normalized' );
                                end
                            end
                            violin_savename = [ANOVA_violin filesep 'c' num2str(v.cycle) '_' v.bandname{b} '_' v.netwname{w} '_' v.netwname{ww}   '_violin.png' ];
                            saveas(gcf, violin_savename);
                        end
                    end
                end
                close all;
            end
        end
    end
    close all;
end 
for c = v.cycle
    for st= v.stage
        for b = 1:v.nbands
            for w=1:v.nnetws
                for ww=1:v.nnetws
                    if ANOVA_p(w,ww,st,c,b) < 0.05 && ANOVA_p(w,ww,st,c,b) > 0
                        ANOVA_sig_F(w,ww,st,c,b) = round(ANOVA_F(w,ww,st,c,b),3); 
                        ANOVA_sig_p(w,ww,st,c,b) = round(ANOVA_p(w,ww,st,c,b),3);
                    else
                        ANOVA_sig_F(w,ww,st,c,b) = NaN;
                        ANOVA_sig_p(w,ww,st,c,b) = NaN;
                    end
               
                end
            end
            ANOVA_sig_F(:,:,st,c,b) = triu(ANOVA_sig_F(:,:,st,c,b));
            ANOVA_sig_p(:,:,st,c,b) = triu(ANOVA_sig_p(:,:,st,c,b));
        end
    end
end

if v.write == 1
for c = v.cycle
    for st= v.stage
        for b = 1:v.nbands
            for w=1:v.nnetws
                for ww=1:v.nnetws
                    if  ANOVA_sig_F(w,ww,st,c,b) ~=0 && ~isnan(ANOVA_sig_F(w,ww,st,c,b)) && ANOVA_sig_p(w,ww,st,c,b)>0.001
                        result_uncorr(w,ww,st,c,b) = {['F(1, ' num2str(v.nsubj-1)  ') = ' num2str(round(ANOVA_sig_F(w,ww,st,c,b),2)) ', p = ' num2str(round(ANOVA_sig_p(w,ww,st,c,b),2))]};
                    elseif  ANOVA_sig_F(w,ww,st,c,b) ~=0 && ~isnan(ANOVA_sig_F(w,ww,st,c,b)) && ANOVA_sig_p(w,ww,st,c,b)<0.001
                        result_uncorr(w,ww,st,c,b) = {['F(1, ' num2str(v.nsubj-1)  ') = ' num2str(round(ANOVA_sig_F(w,ww,st,c,b),2)) ', p < 0.001']};
                    else
                        result_uncorr(w,ww,st,c,b) = {'nan'};
                    end
                end
            end
        end
    end
end
for c = v.cycle
    for st= v.stage
        for b = 1:v.nbands
            for w=1:v.nnetws
                for ww=1:v.nnetws
                    if  ANOVA_sig_F(w,ww,st,c,b) ~=0 && ~isnan(ANOVA_sig_F(w,ww,st,c,b)) && adj_p_mat(w,ww,st,c,b)<0.05
                        result_corr(w,ww,st,c,b) = {['F(1, ' num2str(v.nsubj-1)  ') = ' num2str(round(ANOVA_sig_F(w,ww,st,c,b),2)) ', p = ' num2str(round(adj_p_mat(w,ww,st,c,b),4))]};
                    elseif  ANOVA_sig_F(w,ww,st,c,b) ~=0 && ~isnan(ANOVA_sig_F(w,ww,st,c,b)) && adj_p_mat(w,ww,st,c,b)<0.001
                        result_corr(w,ww,st,c,b) = {['F(1, ' num2str(v.nsubj-1)  ') = ' num2str(round(ANOVA_sig_F(w,ww,st,c,b),2)) ', p < 0.001']};
                    else
                        result_corr(w,ww,st,c,b) = {'nan'};
                    end
                end
            end
        end
    end
end
count_title = 1;
count_label = 4;
for c = v.cycle
    for st= v.stage
        for b = 1:v.nbands
            r= result_uncorr(:,:,st,c,b);
            filename = [v.homedir filesep 'results/ANOVA/interstage/' 'study2_IS_results_uncorr.xlsx'];
            t = {[' cycle ' num2str(v.cycle) ' ' v.stagediff_title{st} ' ' v.bandname{b} ', n= ' num2str(v.nsubj) ]};
            t_cell = ['A' num2str(count_title) ];
            writecell(t, filename, 'Range', t_cell, 'Sheet', v.t_group )

            label_cell = ['A' num2str(count_label) ];
            writecell(v.netwname', filename, 'Range', label_cell , 'Sheet', v.t_group)
            label_b_cell = ['B' num2str(count_label-1) ];
            writecell(v.netwname, filename, 'Range', label_b_cell , 'Sheet', v.t_group)
            
            results_cell = ['B' num2str(count_label)];
            writecell(r, filename, 'Range', results_cell, 'Sheet', v.t_group)

            count_label = count_label + 10;
            count_title = count_title + 10;
        end
    end
end
count_title = 1;
count_label = 4;
for c = v.cycle
    for st= v.stage
        for b = 1:v.nbands
            r= result_corr(:,:,st,c,b);
            filename = [v.homedir filesep 'results/ANOVA/interstage/' 'study2_IS_results_corr.xlsx'];
            t = {[' cycle ' num2str(v.cycle) ' ' v.stagediff_title{st} ' ' v.bandname{b} ', n= ' num2str(v.nsubj) ]};
            t_cell = ['A' num2str(count_title) ];
            writecell(t, filename, 'Range', t_cell, 'Sheet', v.t_group )

            label_cell = ['A' num2str(count_label) ];
            writecell(v.netwname', filename, 'Range', label_cell , 'Sheet', v.t_group)
            label_b_cell = ['B' num2str(count_label-1) ];
            writecell(v.netwname, filename, 'Range', label_b_cell , 'Sheet', v.t_group)
            
            results_cell = ['B' num2str(count_label)];
            writecell(r, filename, 'Range', results_cell, 'Sheet', v.t_group)

            count_label = count_label + 10;
            count_title = count_title + 10;
        end
    end
end
end
% savefile = [study2_anova_dir_stage  '/ANOVA_results'];
% save(savefile,'ANOVA_F','ANOVA_p','ANOVA_sig_F', 'ANOVA_sig_p');
end