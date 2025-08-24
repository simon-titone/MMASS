function study2_stage_ANOVA_3way(v,NETW,design_matrix)
%%%%% NETW = (network, network, subject, stage, cycle, band)
%%%%% BANDS = 1-delta 2-theta 3-alpha 4-sigma 5-beta 6-gamma
anova_dir = [v.homedir filesep 'results/ANOVA/interstage/'  v.t_group ];
if ~exist(anova_dir); mkdir(anova_dir); end; cd(anova_dir)
%% Inter-cycle ANOVA
% ANOVA_p = zeros(v.nnetw,v.nnetw,v.nstages,v.ncycles,v.nbands);
% ANOVA_F = zeros(v.nnetw,v.nnetw,v.nstages,v.ncycles,v.nbands);
for c= v.cycle
    for b = 1:v.nbands
        for w=1:v.nnetws
            for ww=1:v.nnetws

                vect1 = squeeze(NETW(w,ww,:,1,c,b));
                vect2 = squeeze(NETW(w,ww,:,2,c,b));
                vect3 = squeeze(NETW(w,ww,:,3,c,b));
                nsubj = length(~isnan(vect1));
                mat(:,1) = vect1;mat(:,2) = vect2;mat(:,3) = vect3;

                %             t = table(num2str(design_matrix.SessANOVA(1:v.nsubj)), vect1, vect2, vect3, 'VariableNames', {'Group', 'vect1', 'vect2', 'vect3'});
                %             Sess = table([1 2]', 'VariableNames', {'Sessions'});
                %             rm = fitrm(t,'vect1-vect2~1','WithinDesign', Sess);
                %             clear t; clear Sess; clear vect2;
                %             [ranovatbl] = ranova(rm);  % Get within subject (session; session by group) effects
                [tbl,rm] = simple_mixed_anova(mat);
                f = table2array(tbl(3,4));
                p = table2array(tbl(3,5));
                %
                ANOVA_p(w,ww,c,b) =  p;  %ranovatbl{'(Intercept):Sessions','pValue'};
                ANOVA_F(w,ww,c,b) =  f;    %ranovatbl{'(Intercept):Sessions','F'};

                clear ranovatbl; clear rm;
            end
        end
        clear w; clear ww;
        i = 1;
        pval_vect(:,:,c,b) = [ANOVA_p(1,1:v.nnetws,c,b) ANOVA_p(2,2:v.nnetws,c,b) ANOVA_p(3,3:v.nnetws,c,b) ANOVA_p(4,4:v.nnetws,c,b) ANOVA_p(5,5:v.nnetws,c,b)  ANOVA_p(6,6:v.nnetws,c,b)];
        F_anovanetw(:,:,c,b) = [ANOVA_F(1,1:v.nnetws,c,b) ANOVA_F(2,2:v.nnetws,c,b) ANOVA_F(3,3:v.nnetws,c,b) ANOVA_F(4,4:v.nnetws,c,b) ANOVA_F(5,5:v.nnetws,c,b)  ANOVA_F(6,6:v.nnetws,c,b)];
        clear ntw;
        [~, dummy2, ~, dummy4]=fdr_bh(pval_vect(:,:,c,b),v.fdr_thres,'pdep','no');
        crit_p(:,:,c,b) = dummy2;
        adj_p(:,:,c,b) = dummy4;
        adj_p_mat(1,1:6,c,b) = adj_p(1,1:6,c,b); adj_p_mat(2,2:6,c,b) = adj_p(1,7:11,c,b);
        adj_p_mat(3,3:6,c,b) = adj_p(1,12:15,c,b); adj_p_mat(4,4:6,c,b) = adj_p(1,16:18,c,b);
        adj_p_mat(5,5:6,c,b) = adj_p(1,19:20,c,b); adj_p_mat(6,6,c,b) = adj_p(1,21,c,b);

        clear dummy2; clear dummy4;
        index = find(pval_vect(i,:,c,b) == crit_p);

        if isempty(index); crit_F = 1000; % place holder
        else; crit_F = F_anovanetw(index);
        end
        clear index;

        [ivect,jvect]=find(abs(ANOVA_F(:,:,c,b))>=abs(crit_F));
        [ivect_unc,jvect_unc]=find(ANOVA_p(:,:,c,b)<=0.05);
        clear crit_F;

        matrix = ANOVA_F(:,:,c,b); tri_matrix = triu(matrix);

        imagesc(tri_matrix); axis square;   daspect([1 1 1]);
        set(gcf,'Color','white');set(gca,'YDir','reverse');hold on;
        dim = set(gcf, 'Position',  [1000, 900, 500, 450]); %left bottom width height
        h = colorbar; set(get(h,'title'),'string','f-val', 'fontsize', 12); %ylabel(h, 'f-val');
        ANOVA_max=max(max(ANOVA_F(:,:,c,b)));
        if gt(ANOVA_max, 6) ==1 && lt(ANOVA_max, 15) ==1;mm = ANOVA_max; elseif gt(ANOVA_max, 15) ==1 ; mm = 15; else; mm= 6;end;caxis([0 mm]);
        colormap(flipud(pink(10)));
        
        if  v.fig_labels == 1
            title({  [ v.bandname{b} ' Cycle ' num2str(c)  ' Stage Effect' ' , n= ' num2str(nsubj) ]}); % generate corresponding rfx figure %
        else
            title({  [ v.bandname{b} ' Cycle ' num2str(c)  ' Stage Effect']}); % generate corresponding rfx figure %' , n= ' num2str(nsubj)
            
        end
        for a2=1:v.nnetws-1; line([a2+0.5 a2+0.5],[a2+1.5 0.5],'color','k','LineWidth',0.8);line([a2-0.5 v.nnetws+0.5],[a2+0.5 a2+0.5],'color','k','LineWidth',0.8); end
        for b2=1:v.nnetws; line([b2-0.5 b2-0.5],[0.5 v.nnetws+0.5],'color','k','LineWidth',0.001); line([0.5 v.nnetws+0.5],[b2-0.5 b2-0.5],'color','k','LineWidth',0.001);end
        clear a2; clear b2;

        set(gca,'YTick',[1:v.nnetws],'YTickLabel', v.netwname,'Fontsize',16); set(gca,'XTick',[1:v.nnetws],'XTickLabel', v.netwname,'Fontsize',16);

        scatter(jvect_unc,ivect_unc, 80, 'w', 'o' ,'ok','LineWidth',2.5); clear ivect_unc; clear jvect_unc;
        scatter(jvect,ivect,80,'w','ok','filled' ,'LineWidth',2.5); clear ivect; clear jvect;

        saveas(gcf, [anova_dir filesep  'cycle_' v.bandname{b} '_c' num2str(c) '_ANOVA.png'])
        clear a2; clear i; clear mm; clear f

        if v.scatter == 1
            intersess_scatter = [anova_dir filesep v.t_group '_scatter'];
            if ~exist(intersess_scatter); mkdir(intersess_scatter); end; cd(intersess_scatter);
            for w=1:v.nnetws
                for ww=1:v.nnetws
                    if ANOVA_p(w,ww,c,b) <0.05 == 1
                        scatter_band = [intersess_scatter filesep v.bandname{b}];
                        if ~exist(scatter_band); mkdir(scatter_band); end
                        x = squeeze(NETW(w,ww,:,1,c,b)); y = squeeze(NETW(w,ww,:,2,c,b));
                        f = figure; hold on; colormap summer;
                        plot([1 2], [x y])
                        h = bar([1 2], [nanmean(x) nanmean(y)], .5);
                        scatter(repmat(1,1,length(x)), x, 'o', 'r');
                        labelpoints(1,x,v.subj_label,'E', .05);
                        scatter(repmat(2,1,length(y)),y,'o', 'b');
                        labelpoints(2,y,v.subj_label,'E', .05);
                        title({[ v.bandname{b} ' cycle ' num2str(c) ] [ v.netwname{w} '-' v.netwname{ww}  ]}, 'FontSize', 16); % ', n= ' num2str(nsubj)
                        ylabel('connectivity (z-score)', 'FontSize', 16);
                        grouporder={comp_t1omp_t2};
                        xticks([1 2]); xticklabels(grouporder);
                        set(gcf, 'Position',  [10, 10, 300, 500]);
                        f = num2str(round(ANOVA_F(w,ww,c,b),3));
                        p = num2str(round(ANOVA_p(w,ww,c,b),3));
                        text(0.1, 0.98,['F = ' f] ,'Units','normalized');
                        text(0.1, 0.95,['p(UC) = ' p],'Units','normalized');
                        saveas(gcf, [scatter_band filesep comp '_' v.bandname{b} '_' v.netwname{w} '_' v.netwname{ww}   '_sc.png' ]);
                        clear('x', 'y', 'f');
                    end
                end
            end
            close all
        end

        %% Violin plots for each significant btw session comparison
        if v.violin == 1
            ANOVA_violin = [anova_dir filesep v.t_group '_violin'];
            if ~exist(ANOVA_violin); mkdir(ANOVA_violin); end; cd(ANOVA_violin);
            for w=1:v.nnetws
                for ww=1:v.nnetws
                    if ANOVA_p(w,ww,c,b) <0.05 == 1
                        violin_band = [ANOVA_violin filesep v.bandname{b}];
                        if ~exist(violin_band); mkdir(violin_band); end
                        x = [squeeze(NETW(w,ww,:,1,c,b)); squeeze(NETW(w,ww,:,2,c,b)); squeeze(NETW(w,ww,:,3,c,b))];
                        y = [repmat({'NREM2'}, v.nsubj,1) ; repmat({'NREM3'}, v.nsubj,1); repmat({'REM'}, v.nsubj,1)];
                        grouporder={'NREM2', 'NREM3', 'REM'};
                        figure; v1 = violinplot(x, y,'GroupOrder',grouporder, 'ShowMean', true);
                        a = get(gca,'XTickLabel');
                        set(gca,'XTickLabel',a,'fontsize',16)
                        set(gcf, 'Position',  [10, 10, 350, 500]);
                        if v.fig_title == 1
                            title({[ v.bandname{b} ' cycle ' num2str(c) ' ' v.netwname{w} '-' v.netwname{ww} ]}, 'FontSize', 16);
                            ylabel('connectivity (z-score)', 'FontSize', 16);
                        end
                        if v.fig_labels == 1
                            f = num2str(round(ANOVA_F(w,ww,c,b),3)); 
                            text(0.1, 0.98,['F = ' f] ,'Units','normalized');                                 
                            p = round(ANOVA_p(w,ww,c,b),4); 
                                if p < 0.001
                                text(0.1, 0.95,['p(UC) < 0.001 '],'Units','normalized');
                                else
                                text(0.1, 0.95,['p(UC) = ' num2str(p)],'Units','normalized');
                                end
                        end
                        saveas(gcf, [violin_band filesep  'cycle_' num2str(c) '_' v.bandname{b} '_' v.netwname{w} '_' v.netwname{ww}   '_violin.png' ]);
                    end
                end
            end
            close all;
        end
    end
    close all;
end



for c= v.cycle
    for b = 1:v.nbands
        for w=1:v.nnetws
            for ww=1:v.nnetws
                if ANOVA_p(w,ww,c,b) < 0.05 && ANOVA_p(w,ww,c,b) > 0
                    ANOVA_sig_F(w,ww,c,b) = ANOVA_F(w,ww,c,b);
                    ANOVA_sig_p(w,ww,c,b) = ANOVA_p(w,ww,c,b);
                else
                    ANOVA_sig_F(w,ww,c,b) = NaN;
                    ANOVA_sig_p(w,ww,c,b) = NaN;
                end
            end
        end
        ANOVA_sig_F(:,:,c,b) = triu(ANOVA_sig_F(:,:,c,b));
        ANOVA_sig_p(:,:,c,b) = triu(ANOVA_sig_p(:,:,c,b));
    end
end

if v.write == 1
    for c = v.cycle
        for b = 1:v.nbands
            for w=1:v.nnetws
                for ww=1:v.nnetws
                    if  ANOVA_sig_F(w,ww,c,b) ~=0 && ~isnan(ANOVA_sig_F(w,ww,c,b)) && ANOVA_sig_p(w,ww,c,b)>0.001
                        result(w,ww,c,b) = {['F(1, ' num2str(v.nsubj-1)  ') = ' num2str(round(ANOVA_sig_F(w,ww,c,b),3)) ', p = ' num2str(round(ANOVA_sig_p(w,ww,c,b),3))]};
                    elseif  ANOVA_sig_F(w,ww,c,b) ~=0 && ~isnan(ANOVA_sig_F(w,ww,c,b)) && ANOVA_sig_p(w,ww,c,b)<0.001
                        result(w,ww,c,b) = {['F(1, ' num2str(v.nsubj-1)  ') = ' num2str(round(ANOVA_sig_F(w,ww,c,b),3)) ', p < 0.001']};
                    else
                        result(w,ww,c,b) = {'nan'};
                    end
                    if  ANOVA_sig_F(w,ww,c,b) ~=0 && ~isnan(ANOVA_sig_F(w,ww,c,b)) && adj_p_mat(w,ww,c,b)<0.05
                        result_corr(w,ww,c,b) = {['F(1, ' num2str(v.nsubj-1)  ') = ' num2str(round(ANOVA_sig_F(w,ww,c,b),3)) ', p = ' num2str(round(adj_p_mat(w,ww,c,b),3))]};
                    else
                        result_corr(w,ww,c,b) = {'nan'};
                    end
                end
            end
        end
    end

    count_title = 1;
    count_label = 4;
    for c = v.cycle
            for b = 1:v.nbands
                r= result(:,:,c,b);
                filename = [v.homedir filesep 'results/ANOVA/intercycle/' 'study2_IC_results_uncorr.xlsx'];
                t = {[comp_t1 ' vs ' comp_t2 ' ' v.bandname{b} ' ' c ', n= ' num2str(v.nsubj) ]};
                t_cell = ['A' num2str(count_title) ];
                writecell(t, filename, 'Range', t_cell, 'Sheet', v.t_group )
                label_cell = ['A' num2str(count_label) ];
                writecell(v.netwname', filename, 'Range', label_cell , 'Sheet', v.t_group)
                label_b_cell = ['B' num2str(count_label-1) ];
                writecell(v.netwname, filename, 'Range', label_b_cell , 'Sheet', v.t_group)
                results_cell = ['B' num2str(count_label)];
                writecell(r, filename, 'Range', results_cell, 'Sheet', v.t_group)
                clear t; clear t_cell; clear filename; clear r; clear label_cell; clear label_b_cell; clear results_cell; clear t_cell

                r= result_corr(:,:,c,b);
                filename = [v.homedir filesep 'results/ANOVA/intercycle/' 'study2_IC_results_corr.xlsx'];
                t = {[comp_t1 ' vs ' comp_t2 ' ' v.bandname{b} ' ' c ', n= ' num2str(v.nsubj) ]};
                t_cell = ['A' num2str(count_title) ];
                writecell(t, filename, 'Range', t_cell, 'Sheet', v.t_group )
                label_cell = ['A' num2str(count_label) ];
                writecell(v.netwname', filename, 'Range', label_cell , 'Sheet', v.t_group)
                label_b_cell = ['B' num2str(count_label-1) ];
                writecell(v.netwname, filename, 'Range', label_b_cell , 'Sheet', v.t_group)
                results_cell = ['B' num2str(count_label)];
                writecell(r, filename, 'Range', results_cell, 'Sheet', v.t_group)
                clear t; clear t_cell; clear filename; clear r; clear label_cell; clear label_b_cell; clear results_cell; clear t_cell
                count_label = count_label + 10;
                count_title = count_title + 10;
            end
    end
end
end