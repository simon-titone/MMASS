function study3_intersession_ANOVA(v)
%%%%% %NETW: network, network, subject, stage, session, cycle, band
%%%%% BANDS = 1-delta 2-theta 3-alpha 4-sigma 5-beta 6-gamma
%% Load Subjects
% v.ANOVA = [4 9 12 14 16 17 ];
v.ANOVA_NREM2 = [4 9 12 14 15 16 17 21];    % n = 21, P5 & P21 connectivity outliers
v.ANOVA_NREM3 = [4 9 10 12 14 15 16 17 21 28];
v.ANOVA_REM = [1 3 4 6 7 9 10 11 12 13 14 16 17 18 21 24 27 28 29];     % n = 10

% v.group = v.ANOVA_NREM2; stage = 1;
v.group = v.ANOVA_NREM3; stage = 2;
% v.group = v.ANOVA_REM; stage = 3;

v.subj = 1:29; v.group = sort(v.group, 'descend');
for i = 1:length(v.group); v.subj(v.subj(v.group(i))) = []; end
v.nsubj = length(v.subj); % number of independent subjects
study3_load(v);

cd(v.vardir);load('v.mat');
for cc = 1:v.nsubj; c(cc) = v.subj(cc); end
v.subj_label= c'; v.violin_subj_label = [c' ; c'];
%% Inter-session ANOVA

sleep_anova_dir = [v.homedir filesep 'results' filesep 'ANOVA' ];
if ~exist(sleep_anova_dir); mkdir(sleep_anova_dir); end; cd(sleep_anova_dir)

% NETW(NETW == 0) = nan;
%% Inter-session ANOVA
for st= stage
    analysis = [ v.stages{st}];
    ANOVAdir_ISe = [sleep_anova_dir  '/Intersession/' analysis];
    if ~exist(ANOVAdir_ISe); mkdir(ANOVAdir_ISe); end
    for b = 1:v.nbands
        cd(ANOVAdir_ISe)
        for w=1:v.nnetws
            for ww=1:v.nnetws
                clear mat;
                vect1 = squeeze(NETW(w,ww,:,st,1,1,b)); vect1(isnan(vect1)) = [];
                vect2 = squeeze(NETW(w,ww,:,st,2,1,b)); vect2(isnan(vect2)) = [];
                if length(vect1) ~= length(vect2); break; end
                nsubj = length(vect1);
                mat(:,1) = vect1;mat(:,2) = vect2;
                [tbl,~] = simple_mixed_anova(mat);
                ANOVA_F(w,ww,b,st) = table2array(tbl(3,4));
                ANOVA_p(w,ww,b,st) = table2array(tbl(3,5));
            end
        end
        clear w; clear ww;

        pval_fdr_temp = ANOVA_p(:,:,b,st);
        F_val_fdr_temp = ANOVA_F(:,:,b,st);

        pval_vect_temp = pval_fdr_temp(1,:);
        F_val_vect_temp = F_val_fdr_temp(1,:);
        for i = 2:6
            pval_vect_temp = [pval_vect_temp pval_fdr_temp(i,i:6)];
            F_val_vect_temp = [F_val_vect_temp F_val_fdr_temp(i,i:6)];
        end; clear i;

        if  v.FDR == 'all' % runs FDR correction on all comparisons
            pval_vect(:,:,b,st) = pval_vect_temp;
            F_anovanetw(:,:,b,st) = F_val_vect_temp;
        elseif  v.FDR == 'mot'
            pval_vect(:,:,b,st) = pval_fdr_temp(5,:); % runs FDR correction ONLY on motor line
            F_anovanetw(:,:,b,st) = F_val_fdr_temp(5,:);
        end
        clear ntw;
        
        [~, dummy2, ~, dummy4]=fdr_bh(pval_vect(:,:,b,st),v.fdr_thres,'pdep','no');
        crit_p(:,:,b,st) = dummy2;

        adj_p_temp = zeros(6,6);
        index = 1;
        for i = 1:6
            for j = i:6
                adj_p_temp(i,j) = dummy4(index);
                adj_p_temp(j,i) = dummy4(index);
                index = index + 1;
            end
        end; clear i;
        adj_p(:,:,b,st) = triu(adj_p_temp);

        clear dummy2; clear dummy4;
        index = find(pval_vect(1,:,b,st) == crit_p);

        if isempty(index); crit_F = 1000; % place holder
        else; crit_F = F_anovanetw(1,index(1),b,st);
        end
        clear index;

        mm=6;
        [ivect,jvect]=find(abs(ANOVA_F(:,:,b,st))>=abs(crit_F));
        [ivect_unc,jvect_unc]=find(ANOVA_p(:,:,b,st)<=0.05);
        clear crit_F;

        f = figure; f.Position(3) = 490;
        set(gcf,'Color','white'); box OFF;

        matrix = ANOVA_F(:,:,b,st);
        tri_matrix = triu(matrix);

        imagesc(tri_matrix); axis square; caxis([0 mm]); daspect([1 1 1]); %colormap(v.colors_onesided);
        h = colorbar; set(get(h,'title'),'string','f-val', 'fontsize', 12); %ylabel(h, 'f-val');
        %             ANOVA_max=max(max(ANOVA_F(:,:,st,b)));
        %             if gt(ANOVA_max, 6) ==1 && lt(ANOVA_max, 15) ==1;mm = ANOVA_max; elseif gt(ANOVA_max, 15) ==1 ; mm = 15; else; mm= 6;end;caxis([0 mm]);
        colormap(flipud(pink(10)));

        if v.fig_labels == 1
            title({['LRN vs. CTRL ' v.bandname{b} ' ' v.stages{st}  ' {\it n} = ' num2str(nsubj)]}); % generate corresponding rfx figure
        else
            title({['LRN vs. CTRL ' v.bandname{b} ' ' v.stages{st}] }); % generate corresponding rfx figure
        end


        hold on; for a2=1:v.nnetws-1  line([a2+0.5 a2+0.5],[0.5 v.nnetws+0.5],'color','k');   line([0.5 v.nnetws+0.5],[a2+0.5 a2+0.5],'color','k');    end
        clear a2;
        for b2=1:v.nnetws    line([b2+0.5 b2+0.5],[0.5 v.nnetws+0.5],'color','k');   line([0.5 v.nnetws+0.5],[b2+0.5 b2+0.5],'color','k');    end
        set(gca,'YTick',[1:v.nnetws],'YTickLabel', v.netwname,'Fontsize',16);
        set(gca,'XTick',[1:v.nnetws],'XTickLabel', v.netwname,'Fontsize',16);

        scatter(jvect_unc,ivect_unc, 80, 'w', 'o' ,'ok','LineWidth',2.5); clear ivect_unc; clear jvect_unc;
        scatter(jvect,ivect,80,'w','ok','filled' ,'LineWidth',2.5); clear ivect; clear jvect;

        saveas(gcf, [ANOVAdir_ISe filesep v.bandname{b} '_' v.stages{st} '_ANOVA.png'])

        clear a2; clear i; clear mm; clear f


        if v.scatter == 1
            intersess_scatter = [ANOVAdir_ISe '/' analysis '_scatter'];
            if ~exist(intersess_scatter); mkdir(intersess_scatter); end; cd(intersess_scatter);
            for w=1:v.nnetws
                for ww=1:v.nnetws
                    if ANOVA_p(w,ww,b,st) <0.05 == 1
                        x = squeeze(NETW(w,ww,:,st,1,1,b));
                        y = squeeze(NETW(w,ww,:,st,2,1,b));
                        f = figure; hold on; colormap summer;
                        plot([1 2], [x y])
                        h = bar([1 2], [mean(x) mean(y)], .5);
                        scatter(repmat(1,1,length(x)),x,'o', 'r');
                        labelpoints(1,x,v.subj_label,'E', .05);
                        scatter(repmat(2,1,length(y)),y,'o', 'b');
                        labelpoints(2,y,v.subj_label,'E', .05);

                        if v.fig_labels == 1
                            title([v.bandname{b} ' ' v.stages{st} ' ' v.netwname{w} '-' v.netwname{ww}  ', n= ' num2str(nsubj)]);
                        else
                            title([v.bandname{b} ' ' v.stages{st} ' ' v.netwname{w} '-' v.netwname{ww}]);
                        end
                        ylabel('connectivity (z-score)');
                        grouporder={'Control','Learning'};
                        xticks([1 2]); xticklabels(grouporder);
                        set(gcf, 'Position',  [10, 10, 300, 500]);
                        text(0.4, 0.98,['f-value=' num2str(ANOVA_F(w,ww,b,st))],'Units','normalized' );
                        text(0.4, 0.96,[ 'p-value(uncorr)=  ' num2str(ANOVA_p(w,ww,b,st))],'Units','normalized' );
                        saveas(gcf, [intersess_scatter filesep v.bandname{b} '_' v.netwname{w} '_' v.netwname{ww}   '_bar_scatter.png' ]);
                        clear('x','y', 'f');
                    end
                end
            end
            close all
        end

        %% Violin plots for each significant btw session comparison
        if v.violin == 1
            ANOVA_violin = [ANOVAdir_ISe '/' analysis '_violin'];
            if ~exist(ANOVA_violin); mkdir(ANOVA_violin); end; cd(ANOVA_violin);
            for w=1:v.nnetws
                for ww=1:v.nnetws
%                     if ANOVA_p(w,ww,b,st) <0.05 == 1
                        x1 = squeeze(NETW(w,ww,:,st,1,1,b)); x1(isnan(x1)) = [];
                        x2 = squeeze(NETW(w,ww,:,st,2,1,b)); x2(isnan(x2)) = [];
                        x = [x1; x2];
                        y = [repmat({'Control'}, nsubj,1) ; repmat({'Learning'}, nsubj,1)];
                        grouporder={'Control','Learning'};
                        figure; v1 = violinplot(x, y,'GroupOrder',grouporder, 'ShowMean', true);
                        a = get(gca,'XTickLabel');
                        set(gca,'XTickLabel',a,'fontsize',14) %,'FontWeight','bold'
                        ylabel('connectivity (z-score)');
                        set(gcf, 'Position',  [10, 10, 350, 500]);
                        if v.fig_labels == 1
                            T = {[v.stages{st}  ' ' v.bandname{b} ' ' v.netwname{w} '-' v.netwname{ww} ', {\it n} = ' num2str(nsubj)]};
                        else
                            T = {[v.stages{st}  ' ' v.bandname{b} ' ' v.netwname{w} '-' v.netwname{ww}]};
                        end
                        title(T);
                        if v.pval_labels == 1
                            text(0.4, 0.98,['f = ' num2str(round(ANOVA_F(w,ww,b,st),3))],'Units','normalized' ,'FontSize',12);
                            if ANOVA_p(w,ww,b,st) <0.001 == 1
                                text(0.4, 0.96, 'p(UC) <  0.001', 'Units','normalized' ,'FontSize',10);
                            else
                                text(0.4, 0.96,[ 'p(UC) =  ' num2str(round(ANOVA_p(w,ww,b,st),3))],'Units','normalized' ,'FontSize',10);
                            end
                        end
                        saveas(gcf, [v.bandname{b} '_' v.netwname{w} '_' v.netwname{ww}   '_violin.png' ]);
%                     end
                end
            end
            close all;
        end
    end
    close all;
end
clear mat;




%% Writes excel file of results

if v.excel == 1
    for c = 1
        for st= stage
            for b = 1:v.nbands
                for w=1:v.nnetws
                    for ww=1:v.nnetws
                        if ANOVA_p(w,ww,b,st) < 0.05 && ANOVA_p(w,ww,b,st) > 0.000000001
                            ANOVA_sig_F(w,ww,b,st) = ANOVA_F(w,ww,b,st);
                            ANOVA_sig_p(w,ww,b,st) = ANOVA_p(w,ww,b,st);
                        else
                            ANOVA_sig_F(w,ww,b,st) = NaN;
                            ANOVA_sig_p(w,ww,b,st) = NaN;
                        end

                    end
                end
            end
        end
    end
    for st= stage
        for b = 1:v.nbands
            ANOVA_sig_F(:,:,b,st) = triu(ANOVA_sig_F(:,:,b,st));
            ANOVA_sig_p(:,:,b,st) = triu(ANOVA_sig_p(:,:,b,st));
        end
    end

    for c = 1
        for st= stage
            for b = 1:v.nbands
                for w=1:v.nnetws
                    for ww=1:v.nnetws
                        if  ANOVA_sig_F(w,ww,b,st) ~=0 && ~isnan(ANOVA_sig_F(w,ww,b,st)) && ANOVA_sig_p(w,ww,b,st)>0.001
                            result(w,ww,b,st) = {['F(1, ' num2str(v.nsubj-1)  ') = ' num2str(round(ANOVA_sig_F(w,ww,b,st),3)) ', p = ' num2str(round(ANOVA_sig_p(w,ww,b,st),3))]};
                        elseif  ANOVA_sig_F(w,ww,b,st) ~=0 && ~isnan(ANOVA_sig_F(w,ww,b,st)) && ANOVA_sig_p(w,ww,b,st)<0.001
                            result(w,ww,b,st) = {['F(1, ' num2str(v.nsubj-1)  ') = ' num2str(round(ANOVA_sig_F(w,ww,b,st),3)) ', p < 0.001']};
                        else
                            result(w,ww,b,st) = {'nan'};
                        end
                    end
                end
            end
        end
    end
    % for c = 1
    %     for st= stage
    %         for b = 1:v.nbands
    %             for w=1:v.nnetws
    %                 for ww=1:v.nnetws
    %                     if  ANOVA_sig_F(w,ww,b,st) ~=0 && ~isnan(ANOVA_sig_F(w,ww,b,st)) && adj_p_mat(w,ww,b,st)<0.05
    %                         result_corr(w,ww,b,st) = {['F(1, ' num2str(v.nsubj-1)  ') = ' num2str(round(ANOVA_sig_F(w,ww,b,st),3)) ', p = ' num2str(round(adj_p_mat(w,ww,b,st),4))]};
    %                     else
    %                         result_corr(w,ww,b,st) = {'nan'};
    %                     end
    %                 end
    %             end
    %         end
    %     end
    % end

    ANOVA_sig_F(:,:,b,st) = triu(ANOVA_sig_F(:,:,b,st));
    ANOVA_sig_p(:,:,b,st) = triu(ANOVA_sig_p(:,:,b,st));

    count_title = 1;
    count_label = 4;
    v.t_group = 1;
    for c = 1
        for st= stage
            for b = 1:v.nbands
                r= result(:,:,b,st);
                filename = [v.homedir filesep 'results/ANOVA/Intersession/' filesep v.stages{st} filesep 'study3_Intersession_uncorr.xlsx'];
                t = {[ v.stages{st} ' ' v.bandname{b} ', n= ' num2str(v.nsubj) ]};
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
    % count_title = 1;
    % count_label = 4;
    % for c = 1
    %     for st= v.stage
    %         for b = 1:v.nbands
    %             r= result_corr(:,:,b,st);
    %             filename = [v.homedir filesep 'results/ANOVA/interstage/' 'study3_IS_results_corr.xlsx'];
    %             t = {[' cycle ' num2str(v.cycle)  ' EXP vs. CTRL,' v.bandname{b} ', n= ' num2str(v.nsubj) ]};
    %             t_cell = ['A' num2str(count_title) ];
    %             writecell(t, filename, 'Range', t_cell, 'Sheet', v.t_group )
    %
    %             label_cell = ['A' num2str(count_label) ];
    %             writecell(v.netwname', filename, 'Range', label_cell , 'Sheet', v.t_group)
    %             label_b_cell = ['B' num2str(count_label-1) ];
    %             writecell(v.netwname, filename, 'Range', label_b_cell , 'Sheet', v.t_group)
    %
    %             results_cell = ['B' num2str(count_label)];
    %             writecell(r, filename, 'Range', results_cell, 'Sheet', v.t_group)
    %
    %             count_label = count_label + 10;
    %             count_title = count_title + 10;
    %         end
    %     end
    % end
    % savefile = [ANOVAdir_ISe filesep 'ANOVA_results'];
    % save(savefile,'ANOVA_F','ANOVA_p');
end
end