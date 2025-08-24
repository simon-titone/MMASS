function sleep_ANOVA(v,NETW,design_matrix)
%%%%% NETW = (network, network, subject, stage, session, band)

sleep_anova_dir = [v.homedir filesep 'results' filesep 'ANOVA'];
if ~exist(sleep_anova_dir); mkdir(sleep_anova_dir); end; cd(sleep_anova_dir)

%% Inter-Stage ANOVA
for se = 1:v.nsess
    for st= 1:v.nstagecorr
        sleep_anova_dir_stage = [sleep_anova_dir filesep v.sess{se} filesep v.stagediff{st}];
        if ~exist(sleep_anova_dir_stage); mkdir(sleep_anova_dir_stage); end
        
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
            cd(sleep_anova_dir_stage)
            for w=1:v.nnetws
                for ww=1:v.nnetws
                    
                    vect1 = squeeze(NETW(w,ww,:,st1,se,b));
                    vect2 = squeeze(NETW(w,ww,:,st2,se,b));
                    
                    t = table(num2str(design_matrix.SessANOVA(1:v.nsubj)), vect1, vect2, 'VariableNames', {'Group', 'stage1', 'stage2'});
                    
                    Sess = table([1 2]', 'VariableNames', {'Sessions'});
                    rm = fitrm(t,'stage1-stage2~1','WithinDesign', Sess);
                    clear t; clear Sess; clear vect2;
                    
                    [ranovatbl] = ranova(rm);  % Get within subject (session; session by group) effects
                    ANOVA_p(w,ww,b,se,st) = ranovatbl{'(Intercept):Sessions','pValue'};
                    ANOVA_F(w,ww,b,se,st) = ranovatbl{'(Intercept):Sessions','F'};
                    
                    clear ranovatbl;
                    clear rm;
                end
            end
            clear w; clear ww;
            
            %             for i=1:size(design_matrix.SessANOVA,2)
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
            [ivect,jvect]=find(abs(ANOVA_F(:,:,b,se,st))>=abs(crit_F));
            [ivect_unc,jvect_unc]=find(ANOVA_p(:,:,b,se,st)<=0.05);
            clear crit_F;
            
            figure(); set(gcf,'Color','white'); box OFF;
            imagesc(ANOVA_F(:,:,b,se,st)); axis square; caxis([0 mm]); colormap(v.colors_onesided); colorbar; daspect([1 1 1]);
            title([ v.sess{se} ' ' v.stagediff_title{st} ' ' v.bandname{b} ' ANOVA']); % generate corresponding rfx figure
            
            %             end
            
            hold on; for a2=1:v.nnetws-1  line([a2+0.5 a2+0.5],[0.5 v.nnetws+0.5],'color','k');   line([0.5 v.nnetws+0.5],[a2+0.5 a2+0.5],'color','k');    end
            clear a2;
            for b2=1:v.nnetws    line([b2+0.5 b2+0.5],[0.5 v.nnetws+0.5],'color','k');   line([0.5 v.nnetws+0.5],[b2+0.5 b2+0.5],'color','k');    end
            set(gca,'YTick',[1:v.nnetws],'YTickLabel', v.netwname,'Fontsize',16);
            set(gca,'XTick',[1:v.nnetws],'XTickLabel', v.netwname,'Fontsize',16);
            
            scatter(jvect_unc,ivect_unc,36,'ok'); clear ivect_unc; clear jvect_unc;
            scatter(jvect,ivect,9,'ok','filled'); clear ivect; clear jvect;
            
            saveas(gcf, [ v.sess{se} '_' v.stagediff{st} '_' v.bandname{b} '_ANOVA.png'])
            
            clear a2; clear i; clear mm;
            
            
            if v.violin == 1
                ANOVA_violin = [sleep_anova_dir_stage 'ANOVA_violin'];
                if ~exist(ANOVA_violin); mkdir(ANOVA_violin); end; cd(ANOVA_violin);
                for w=1:v.nnetws
                    for ww=1:v.nnetws
                        if ANOVA_p(w,ww,b,se,st) <0.05 == 1
                            x = [squeeze(NETW(w,ww,:,1,b)); squeeze(NETW(w,ww,:,2,b))];
                            y = [repmat({[st1_title]}, v.nsubj,1) ; repmat({[st2_title]}, v.nsubj,1)];
                            z = [repmat(1,v.nsubj,1) ; repmat(2, v.nsubj,1)];
                            grouporder={st1_title,st2_title};
                            figure; v1 = violinplot(x, y,'GroupOrder',grouporder);
                            ylabel('connectivity (z-score)');
                            set(gcf, 'Position',  [10, 10, 350, 500]);
                            title([v.sess{se} ' ' v.bandname{b} ' ANOVA Violin ' v.netwname{w} '-' v.netwname{ww} ', n= ' num2str(v.nsubj)]);
                            text(0.4, 0.98,['f-value=' num2str(ANOVA_F(w,ww,b,se,st))],'Units','normalized' );
                            text(0.4, 0.96,[ 'p-value(uncorr)=  ' num2str(ANOVA_p(w,ww,b,se,st))],'Units','normalized' );
                            saveas(gcf, [v.bandname{b} '_' v.netwname{w} '_' v.netwname{ww}   '_violin.png' ]);
                        end
                    end
                end
                close all;
            end
        end
    end
    close all;
end

%% Inter-session ANOVA

for st= 1:v.nstagecorr
    sleep_anova_dir_sess = [sleep_anova_dir  filesep 'inter_sess_' v.stages{st}];
    if ~exist(sleep_anova_dir_sess); mkdir(sleep_anova_dir_sess); end;
    for b = 1:v.nbands
        cd(sleep_anova_dir_sess)
        for w=1:v.nnetws
            for ww=1:v.nnetws
                
                btw_sess_vect1 = squeeze(NETW(w,ww,:,st,1,b));
                btw_sess_vect2 = squeeze(NETW(w,ww,:,st,2,b));
                t = table(num2str(design_matrix.SessANOVA(1:v.nsubj)), btw_sess_vect1, btw_sess_vect2, 'VariableNames', {'Group', 'stage1', 'stage2'});
                
                Sess = table([1 2]', 'VariableNames', {'Sessions'});
                rm = fitrm(t,'stage1-stage2~1','WithinDesign', Sess);
                clear t; clear Sess; clear vect2;
                
                [ranovatbl] = ranova(rm);  % Get within subject (session; session by group) effects
                ANOVA_p(w,ww,b,st) = ranovatbl{'(Intercept):Sessions','pValue'};
                ANOVA_F(w,ww,b,st) = ranovatbl{'(Intercept):Sessions','F'};
                
                clear ranovatbl;
                clear rm;
            end
        end
        clear w; clear ww;
        
        %             for i=1:size(design_matrix.SessANOVA,2)
        i = 1;
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
        [ivect,jvect]=find(abs(ANOVA_F(:,:,b,st))>=abs(crit_F));
        [ivect_unc,jvect_unc]=find(ANOVA_p(:,:,b,st)<=0.05);
        clear crit_F;
        
        figure(); set(gcf,'Color','white'); box OFF;
        imagesc(ANOVA_F(:,:,b,st)); axis square; caxis([0 mm]); colormap(v.colors_onesided); colorbar; daspect([1 1 1]);
        title([v.bandname{b} ' ' v.stages{st} ' ANOVA']); % generate corresponding rfx figure
        
    end
    
    hold on; for a2=1:v.nnetws-1  line([a2+0.5 a2+0.5],[0.5 v.nnetws+0.5],'color','k');   line([0.5 v.nnetws+0.5],[a2+0.5 a2+0.5],'color','k');    end
    clear a2;
    for b2=1:v.nnetws    line([b2+0.5 b2+0.5],[0.5 v.nnetws+0.5],'color','k');   line([0.5 v.nnetws+0.5],[b2+0.5 b2+0.5],'color','k');    end
    set(gca,'YTick',[1:v.nnetws],'YTickLabel', v.netwname,'Fontsize',16);
    set(gca,'XTick',[1:v.nnetws],'XTickLabel', v.netwname,'Fontsize',16);
    
    scatter(jvect_unc,ivect_unc,36,'ok'); clear ivect_unc; clear jvect_unc;
    scatter(jvect,ivect,9,'ok','filled'); clear ivect; clear jvect;
    
    saveas(gcf, [v.bandname{b} '_' v.stages{st} '_ANOVA.png'])
    
    clear a2; clear i; clear mm;
    
    if v.violin == 1
        ANOVA_violin = [sleep_anova_dir_sess 'ANOVA_violin'];
        if ~exist(ANOVA_violin); mkdir(ANOVA_violin); end; cd(ANOVA_violin);
        for w=1:v.nnetws
            for ww=1:v.nnetws
                if ANOVA_p(w,ww,b,st) <0.05 == 1
                    
                    x = [squeeze(NETW(w,ww,:,1,b)); squeeze(NETW(w,ww,:,2,b))];
                    y = [repmat({'Control'}, v.nsubj,1) ; repmat({'Experimental'}, v.nsubj,1)];
                    grouporder={'Control','Experimental'};
                    figure; v1 = violinplot(x, y,'GroupOrder',grouporder);
                    ylabel('connectivity (z-score)');
                    set(gcf, 'Position',  [10, 10, 350, 500]);
                    title([v.bandname{b} ' ANOVA Violin ' v.netwname{w} '-' v.netwname{ww} ', n= ' num2str(v.nsubj)]);
                    text(0.4, 0.98,['f-value=' num2str(ANOVA_F(w,ww,b,st))],'Units','normalized' );
                    text(0.4, 0.96,[ 'p-value(uncorr)=  ' num2str(ANOVA_p(w,ww,b,st))],'Units','normalized' );
                    saveas(gcf, [v.bandname{b} '_' v.netwname{w} '_' v.netwname{ww}   '_violin.png' ]);
                end
            end
        end
        close all;
    end
end
end

