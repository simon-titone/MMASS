function sleep_interstage_ANOVA(v)
%%%%% NETW = (network, network, subject, stage, session,cycle, band)
%% Load Subjects
cd(v.vardir);load('v.mat');
v.stage_ANOVA = [4 9 12 14 16 17 ];
v.stage_ANOVA_REM = [1 3 4 7 9 11 12 13 14 16 17 18 21 24 27 28 29];   % n = 12
v.interstage_wake = [1 2 4 10 12 13 16 17 18 22];
v.interstage_wake_REM = [1 2 3 4 7 9 10 11 12 13 16 17 18 21 22 24 27 28 29 19];
sleep_anova_dir = [v.homedir  '/results/ANOVA'];
if ~exist(sleep_anova_dir); mkdir(sleep_anova_dir); end; cd(sleep_anova_dir)
%% Inter-Stage ANOVA
for se = 2
    for c= 1
        for stdiff= 1:3
            sleep_anova_dir_stage = [sleep_anova_dir  '/interstage/' v.stagediff{stdiff}];
            if ~exist(sleep_anova_dir_stage); mkdir(sleep_anova_dir_stage); end
            if stdiff ==1
                st1 = 1; st2 = 2; % compares NREM2 to NREM3
                st1_title = 'NREM2'; st2_title = 'NREM3';
                v.group = v.stage_ANOVA;
            elseif stdiff ==2
                st1 = 1; st2 = 3; % compares NREM2 to REM
                st1_title = 'NREM2'; st2_title = 'REM';
                v.group = v.stage_ANOVA_REM;
            elseif stdiff ==3
                st1 = 2; st2 = 3; % compares NREM3 to REM
                st1_title = 'NREM3'; st2_title = 'REM';
                v.group = v.stage_ANOVA_REM;
            elseif stdiff ==4
                st1 = 1; st2 = 4; % compares NREM2 to wake
                st1_title = 'NREM2'; st2_title = 'RS';
                v.group = v.interstage_wake;
            elseif stdiff ==5
                st1 = 2; st2 = 4; % compares NREM3 to wake
                st1_title = 'NREM3'; st2_title = 'RS';
                v.group = v.interstage_wake;
            end

            for b = 1:v.nbands
                cd(sleep_anova_dir_stage)
                for w=1:v.nnetws
                    for ww=1:v.nnetws
                        v.subj = 1:29; v.group = sort(v.group, 'descend');
                        for i = 1:length(v.group); v.subj(v.subj(v.group(i))) = []; end
                        v.nsubj = length(v.subj); % number of independent subjects
                        study3_load(v);

                        cd(v.vardir);load('v.mat');
                        for cc = 1:v.nsubj; label(cc) = v.subj(cc); end
                        v.subj_label= label'; v.violin_subj_label = [label' ; label'];

                        vect1 = squeeze(NETW(w,ww,:,st1,se,c,b));
                        vect2 = squeeze(NETW(w,ww,:,st2,se,c,b));
                        mat(:,1) = vect1;mat(:,2) = vect2;
                        [tbl,rm] = simple_mixed_anova(mat);
                        f = table2array(tbl(3,4));
                        p = table2array(tbl(3,5));

                        ANOVA_p(w,ww,b,c,1) =  p;    %ranovatbl{'(Intercept):Sessions','pValue'};
                        ANOVA_F(w,ww,b,c,1) =  f;    %ranovatbl{'(Intercept):Sessions','F'};

                        clear ranovatbl;
                        clear rm;
                    end
                end
                clear w; clear ww;

                i=1;
                for ntw = 1:1:v.nnetws
                    if  v.FDR == 'all'
                        pval.anovanetw(i, ntw) = ANOVA_p(1,ntw,i);
                        F_anovanetw(i,ntw) = ANOVA_F(1,ntw,i);
                    elseif  v.FDR == 'mot'
                        pval.anovanetw(i, ntw) = ANOVA_p(5,ntw,i); % runs FDR correction ONLY on motor line
                        F_anovanetw(i,ntw) = ANOVA_F(5,ntw,i);
                        ntw = ntw + 1;
                    end
                end
                clear count; clear nnetw_col;

                [~, dummy2, ~, dummy4]=fdr_bh(pval.anovanetw(i,:),v.fdr_thres,'pdep','no');
                crit_p.anovanetw(i) = dummy2;
                adj_p.anovanetw(i,1:length(dummy4)) = dummy4;
                clear dummy2; clear dummy4;
                index = find(pval.anovanetw(i,:) == crit_p.anovanetw(i));

                if isempty(index); crit_F = 1000; % place holder
                else; crit_F = F_anovanetw(i,index(1));
                end
                clear index;

                mm=6;
                [ivect,jvect]=find(abs(ANOVA_F(:,:,b,c,1))>=abs(crit_F));
                [ivect_unc,jvect_unc]=find(ANOVA_p(:,:,b,c,1)<=0.05);
                clear crit_F;

                figure(); set(gcf,'Color','white'); box OFF;

                matrix = ANOVA_F(:,:,b,c,1);
                tri_matrix = triu(matrix);

                imagesc(tri_matrix); axis square; caxis([-mm mm]);
                %                 colormap(v.colors);
                daspect([1 1 1]); hold on
                h = colorbar; set(get(h,'title'),'string','f-val', 'fontsize', 12); %ylabel(h, 'f-val');

                ANOVA_max=max(max(ANOVA_F(:,:,b,c,1)));
                if gt(ANOVA_max, 6) ==1 && lt(ANOVA_max, 15) ==1;mm = ANOVA_max; elseif gt(ANOVA_max, 15) ==1 ; mm = 15; else; mm= 6;end;caxis([0 mm]);
                colormap(flipud(pink(10)));


                if  v.fig_labels == 1
                    title({  ['EXP '  v.bandname{b} ' ' st1_title ' vs. ' st2_title ' , n= ' num2str(v.nsubj) ]}); % generate corresponding rfx figure %
                else
                    title({  ['EXP ' v.bandname{b} ' ' st1_title ' vs. ' st2_title]}); % generate corresponding rfx figure %' , n= ' num2str(nsubj)
                end
                for a2=1:v.nnetws-1; line([a2+0.5 a2+0.5],[a2+1.5 0.5],'color','k','LineWidth',0.8);line([a2-0.5 v.nnetws+0.5],[a2+0.5 a2+0.5],'color','k','LineWidth',0.8); end
                for b2=1:v.nnetws; line([b2-0.5 b2-0.5],[0.5 v.nnetws+0.5],'color','k','LineWidth',0.001); line([0.5 v.nnetws+0.5],[b2-0.5 b2-0.5],'color','k','LineWidth',0.001);end
                clear a2; clear b2;

                set(gca,'YTick',[1:v.nnetws],'YTickLabel', v.netwname,'Fontsize',16);
                set(gca,'XTick',[1:v.nnetws],'XTickLabel', v.netwname,'Fontsize',16);

                scatter(jvect_unc,ivect_unc, 80, 'w', 'o' ,'ok','LineWidth',2.5); clear ivect_unc; clear jvect_unc;
                scatter(jvect,ivect,80,'w','ok','filled' ,'LineWidth',2.5); clear ivect; clear jvect;

                savefile = [sleep_anova_dir_stage filesep v.sess{se} '_' v.stagediff{stdiff} '_' v.bandname{b} '_ANOVA.png'];
                saveas(gcf, savefile)

                clear i; clear mm;

                if v.scatter == 1
                    intersess_scatter = [sleep_anova_dir_stage '/scatter/' v.bandname{b}];
                    for w=1:v.nnetws
                        for ww=1:v.nnetws
                            if ANOVA_p(w,ww,b,c,1) <0.05 == 1
                                if ~exist(intersess_scatter); mkdir(intersess_scatter); end; cd(intersess_scatter);

                                x = squeeze(NETW(w,ww,:,st1,se,c,b));
                                y = squeeze(NETW(w,ww,:,st2,se,c,b));
                                f = figure; hold on; colormap summer;
                                plot([1 2], [x y])
                                h = bar([1 2], [mean(x) mean(y)], .5);
                                scatter(repmat(1,1,length(x)),x,'o', 'r');
                                labelpoints(1,x,v.subj_label,'E', .05);
                                scatter(repmat(2,1,length(y)),y,'o', 'b');
                                labelpoints(2,y,v.subj_label,'E', .05);

                                if  v.fig_labels == 1
                                    title([v.bandname{b} ' Inter-stage ' v.netwname{w} '-' v.netwname{ww}  ', n= ' num2str(v.nsubj)]);
                                else
                                    title([v.bandname{b} ' Inter-stage ' v.netwname{w} '-' v.netwname{ww} ]);
                                end
                                ylabel('connectivity (z-score)');
                                grouporder={st1_title,st2_title};
                                xticks([1 2]); xticklabels(grouporder);
                                set(gcf, 'Position',  [10, 10, 300, 500]);
                                text(0.4, 0.98,['f = ' num2str(round(ANOVA_F(w,ww,b,c,1),3))],'Units','normalized' ,'FontSize',12);
                                if ANOVA_p(w,ww,b,c,1) <0.001 == 1
                                    text(0.4, 0.96, 'p(UC) <  0.001', 'Units','normalized' ,'FontSize',10);
                                else
                                    text(0.4, 0.96,[ 'p(UC) =  ' num2str(round(ANOVA_p(w,ww,b,c,1),3))],'Units','normalized' ,'FontSize',10);
                                end
                                saveas(gcf, [intersess_scatter filesep st1_title '_' st2_title '_' v.bandname{b} '_' v.netwname{w} '_' v.netwname{ww}   '_scatter.png' ]);
                                clear('x','y', 'f');
                            end
                        end
                    end
                    close all
                end
                if v.violin == 1
                    ANOVA_violin = [sleep_anova_dir_stage '/violin/' v.bandname{b}];
                    for w=1:v.nnetws
                        for ww=1:v.nnetws
                            if ANOVA_p(w,ww,b,c,1) <0.05 == 1
                                if ~exist(ANOVA_violin); mkdir(ANOVA_violin); end; cd(ANOVA_violin);
                                x = [squeeze(NETW(w,ww,:,st1,se,c,b)); squeeze(NETW(w,ww,:,st2,se,c,b))];
                                y = [repmat({[st1_title]}, v.nsubj,1) ; repmat({[st2_title]}, v.nsubj,1)];
                                z = [repmat(1,v.nsubj,1) ; repmat(2, v.nsubj,1)];
                                grouporder={st1_title,st2_title};
                                figure; v1 = violinplot(x, y,'GroupOrder',grouporder, 'ShowMean', true);
                                ylabel('connectivity (z-score)');
                                set(gcf, 'Position',  [10, 10, 350, 500]);
                                if  v.fig_labels == 1
                                T = {[v.sess{se} ' Inter-stage ANOVA '] [v.bandname{b} ' ' v.netwname{w} '-' v.netwname{ww} ', n= ' num2str(v.nsubj)]};
                                else
                                T = {[v.sess{se} ' Inter-stage ANOVA '] [v.bandname{b} ' ' v.netwname{w} '-' v.netwname{ww} ]};
                                end
                                title(T);
                                if v.pval_labels == 1
                                    text(0.4, 0.98,['f = ' num2str(round(ANOVA_F(w,ww,b,c,1),3))],'Units','normalized' ,'FontSize',12);
                                    if ANOVA_p(w,ww,b,c,1) <0.001 == 1
                                        text(0.4, 0.96, 'p(UC) <  0.001', 'Units','normalized' ,'FontSize',10);
                                    else
                                        text(0.4, 0.96,[ 'p(UC) =  ' num2str(round(ANOVA_p(w,ww,b,c,1),3))],'Units','normalized' ,'FontSize',10);
                                    end
                                end
                                saveas(gcf, [st1_title '_' st2_title '_' v.bandname{b} '_' v.netwname{w} '_' v.netwname{ww}   '_violin.png' ]);
                            end
                        end
                    end
                    close all;
                end
            end
        end
        close all;
    end
end
end