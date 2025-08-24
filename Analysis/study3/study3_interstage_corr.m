function study3_interstage_corr(v)
%%%%% NETW: network, network, subject, stage, session,cycle, band
NREM2_NREM3 = [4 12 16 17];
NREM2_REM = [1 3 4 7 9 11 12 13 16 17 18 21 24 27 28 29];
NREM3_REM = [1 3 4 7 9 11 12 13 16 17 18 21 24 27 28 29];
NREM2_RS = [1 2 4 10 12 13 16 17 18 22];
NREM3_RS = [1 2 4 10 12 13 16 17 18 22];
REM_RS= [1 2 3 4 7 9 10 11 12 13 16 17 18 21 22 24 27 28 29];

%% Interstage correlation
for se = 2
    for st= 4
        inter_stage_corr_dir = [v.homedir filesep 'results' filesep 'Correlation' filesep 'interstage_corr' filesep v.stagediff{st}];
        if ~exist(inter_stage_corr_dir); mkdir(inter_stage_corr_dir); end; cd(inter_stage_corr_dir)

        if st == 1; st1 = 1; st2 = 2; v.group = NREM2_NREM3;
        elseif st == 2; st1 = 1; st2 = 3; v.group = NREM2_REM;
        elseif st == 3; st1 = 2; st2 = 3; v.group = NREM3_REM;
        elseif st == 4; st1 = 1; st2 = 4; v.group = NREM2_RS;
        elseif st == 5; st1 = 2; st2 = 4; v.group = NREM3_RS;
        elseif st == 6; st1 = 3; st2 = 4; v.group = REM_RS;
        end

        v.subj = 1:29; v.group = sort(v.group, 'descend');
        for i = 1:length(v.group); v.subj(v.subj(v.group(i))) = []; end
        v.nsubj = length(v.subj);
        study3_load(v); cd(v.vardir);load('v.mat');
        for cc = 1:v.nsubj; c(cc) = v.subj(cc); end; v.subj_label= c'; v.violin_subj_label = [c' ; c'];

        for w=1:v.nnetws
            for ww=1:v.nnetws
                for b = 1:v.nbands
                    stage1 = squeeze(NETW(w,ww,:,st1,se,1,b));
                    stage2 = squeeze(NETW(w,ww,:,st2,se,1,b)) ;
                    [val,px]=corr(stage1,stage2);
                    r_val(w,ww,se,st,b)=val;
                    p_val(w,ww,se,st,b)=px;
                end
            end
        end
        clear val; clear px; clear w; clear ww; clear vect;

        for b=1:v.nbands

            pval_fdr = p_val(:,:,se,st,b);

            if v.FDR == 'mot'
                p_val_temp_FDR_correction  = p_val(5,:);
            elseif v.FDR == 'all'
                p_val_temp_FDR_correction = [ p_val(1:v.nnetw)  p_val(v.nnetw+2:v.nnetw*2)  p_val(v.nnetw*2+3:v.nnetw*3)  p_val(v.nnetw*3+4:v.nnetw*4)  p_val(v.nnetw*4+5:v.nnetw*5)  p_val(v.nnetw*5+6:v.nnetw*6)];
            end

            [~, dummy2, ~, dummy4]=fdr_bh(p_val_temp_FDR_correction, v.fdr_thres,'pdep','no');
            crit_p = dummy2;
            adj_p(1:length(dummy4)) = dummy4;
            clear dummy2; clear dummy4;

            [ivect,jvect]=find(p_val(:,:,se,st,b)<=crit_p);
            [ivect_unc,jvect_unc]=find(p_val(:,:,se,st,b)<=v.fdr_thres);

            matrix = r_val(:,:,se,st,b); tri_matrix = triu(matrix);

            imagesc(tri_matrix); axis square; % colormap(v.colors); daspect([1 1 1]); hold on
            h = colorbar; set(get(h,'title'),'string','r-val', 'fontsize', 12); %ylabel(h, 'f-val');
            dim = set(gcf, 'Position',  [1000, 900, 500, 450]); %left bottom width height
            colormap('redblue'); limit = max(max(abs(tri_matrix))); caxis([-limit limit]);

            T = {[v.sess{se} ' ' v.stagecorr{st} ] [ v.bandname{b}  ', {\it n} = ' num2str(v.nsubj)] };

            title(T); % generate corresponding rfx figure

            for a2=1:v.nnetws-1; line([a2+0.5 a2+0.5],[a2+1.5 0.5],'color','k','LineWidth',0.8);line([a2-0.5 v.nnetws+0.5],[a2+0.5 a2+0.5],'color','k','LineWidth',0.8); end
            for b2=1:v.nnetws; line([b2-0.5 b2-0.5],[0.5 v.nnetws+0.5],'color','k','LineWidth',0.001); line([0.5 v.nnetws+0.5],[b2-0.5 b2-0.5],'color','k','LineWidth',0.001);end
            clear a2; clear b2;

            set(gca,'YTick',[1:v.nnetws],'YTickLabel',v.netwname,'Fontsize',16);
            set(gca,'XTick',[1:v.nnetws],'XTickLabel',v.netwname,'Fontsize',16);
            hold on;
            scatter(jvect_unc,ivect_unc, 80, 'w', 'o' ,'ok','LineWidth',2.5); clear ivect_unc; clear jvect_unc;
            scatter(jvect,ivect,80,'w','ok','filled' ,'LineWidth',2.5); clear ivect; clear jvect;

            S = [ inter_stage_corr_dir filesep v.bandname{b} '_' v.sess{se} '_' v.stagediff{st} '.png' ];
            saveas(gcf, S);

            clear CC_zscore_cov_InterSession; clear pval_Behav;
            clear i;

            if v.scatter == 1
                for w=1:v.nnetws
                    for ww=1:v.nnetws
                        if p_val(w,ww,se,st,b) <0.05 == 1
                            scatter_band = [inter_stage_corr_dir  '/scatter/'  v.bandname{b}];
                            if ~exist(scatter_band); mkdir(scatter_band); end

                            figure
                            x = squeeze(NETW(w,ww,:,st1,se,1,b));
                            y = squeeze(NETW(w,ww,:,st2,se,1,b));
                            scatter(x,y,[], 'b', 'filled');lsline;
                            mdl = fitlm(x, y, 'linear');
                            [y_pred,y_ci] = predict(mdl,x, 'Alpha', 0.05);
                            Y = [x, y_ci,y_pred];
                            Y = sortrows(Y,1);
                            eb = [Y(:,4)-Y(:,2)];
                            boundedline(Y(:,1),Y(:,4),eb, '-r','alpha');
                            clear mdl; clear y_pred; clear y_ci; clear Y; clear eb;

                            T_scatter = {[ v.sess{se} ' ' v.stagecorr{st}  ] [ v.bandname{b} ' '   v.netwname{w} '-' v.netwname{ww}  ', {\it n} = ' num2str(v.nsubj)]};

                            title(T_scatter); xlabel(v.stages{st1}); ylabel(v.stages{st2});
                            set(gcf, 'Position',  [100, 100, 400, 320]); %left bottom width height

                            if v.subj_labels == 1
                                labelpoints(x,y,v.subj_label,'E', .05);
                            end
                            if v.pval_labels == 1
                                text(0.6, 0.98,['r = ' num2str(round(r_val(w,ww,se,st,b), 3))],'Units','normalized' ,'FontSize',14);
                                if p_val(w,ww,se,st,b) <0.001 == 1
                                    text(0.6, 0.93, 'p(UC) < 0.001 ' ,'Units','normalized' ,'FontSize',14);
                                else
                                    text(0.6, 0.93,['p(UC) =  ' num2str(round(p_val(w,ww,se,st,b),3))],'Units','normalized' ,'FontSize',14);
                                end
                            end
                            S_scatter = [ scatter_band filesep v.stagediff{st} '_corr'  v.netwname{w} '_' v.netwname{ww} '.png' ];
                            saveas(gcf, S_scatter);
                            clear('x','y');
                        end
                    end
                end
            end
        end
    end
end
close all
end