function study2_stage_ANOVA_ratio(v,NETW)
%%%%% NETW = (network, network, subject, stage, cycle, band)
%%%%% BANDS = 1-delta 2-theta 3-alpha 4-sigma 5-beta 6-gamma
anova_dir = [v.homedir filesep 'results/ANOVA/ratio/'  v.t_group ];
if ~exist(anova_dir); mkdir(anova_dir); end; cd(anova_dir)
%% Between-stage ANOVA of the ratio of between/within network conn, performed per cycle
% Calculates ratio of within/between network conn per ppt
for c = v.cycle
    for st = 1:3
        for b = 1:v.nbands
            for i = 1:v.nsubj
                for w = 1:v.nnetw
                    withinnetw(i,b,st,c) = NETW(w,w,i,st,c,b);
                    internetw = setdiff(1:v.nnetw,w);
                    for n2 = 1:length(internetw)
                        btwnetw(i,b,st,c) = mean(NETW(internetw,w,i,st,c,b));
                    end
                    ratio(i,b,st,c) = withinnetw(i,b,st,c)/btwnetw(i,b,st,c);
                end
            end
        end
    end
end

c= v.cycle;
for b = 1:v.nbands
    for w=1:v.nnetws
        vect1 = squeeze(ratio(:,b,1,c));
        vect2 = squeeze(ratio(:,b,2,c));
        vect3 = squeeze(ratio(:,b,3,c));
        nsubj = length(~isnan(mean([vect1 vect2 vect3], 2)));
        mat(:,1) = vect1;mat(:,2) = vect2;mat(:,3) = vect3;
        [tbl,rm] = simple_mixed_anova(mat,[],{'stage'});
        f = table2array(tbl(3,4));
        p = table2array(tbl(3,5));
        ANOVA_p(b,c) =  p;  %ranovatbl{'(Intercept):Sessions','pValue'};
        ANOVA_F(b,c) =  f;    %ranovatbl{'(Intercept):Sessions','F'};
        clear ranovatbl; clear rm;
    end
end

%% Violin plot for ratio of w/in and btw network connectivity
if v.violin == 1
    c= v.cycle;
    for b = 1:v.nbands

        ANOVA_violin = [anova_dir filesep 'violin_ratio'];
        if ~exist(ANOVA_violin); mkdir(ANOVA_violin); end; cd(ANOVA_violin);

        r1 = ratio(:,b,1,c);        r2 = ratio(:,b,2,c);         r3 = ratio(:,b,3,c);
        x = [r1 r2 r3];


        y = [repmat({'NREM2'}, v.nsubj,1) ; repmat({'NREM3'}, v.nsubj,1); repmat({'REM'}, v.nsubj,1)];
        z = [repmat(1,v.nsubj,1) ; repmat(2, v.nsubj,1)];
        grouporder={'NREM2', 'NREM3', 'REM'};
        figure;
        v1 = violinplot(x, y,'GroupOrder',grouporder, 'ShowMean', true);
        ylabel('connectivity (z-score)');
        a = get(gca,'XTickLabel');
        set(gca,'XTickLabel',a,'fontsize',16)
        set(gcf, 'Position',  [10, 10, 350, 500]);
        if v.fig_title == 1
            title({['Cycle ' num2str(v.cycle) ' ' v.bandname{b} ', n = ' num2str(v.nsubj)] [  'within vs. between network connectivity' ]}, 'fontsize', 16);
        end
        if v.fig_labels == 1
            F = num2str(round(ANOVA_F(b,c))); p = round(ANOVA_p(b,c),4);
            text(0.4, 0.98,['F =' F],'Units','normalized' );
            if p > 0.001 == 1
                text(0.4, 0.96,[ 'p(UC) =  ' num2str(p)],'Units','normalized' );
            else
                text(0.4, 0.96,[ 'p(UC) < 0.001'],'Units','normalized' );
            end
        end
        violin_savename = [ANOVA_violin filesep 'c' num2str(v.cycle) '_' v.bandname{b} '_ratio_violin.png' ];
        saveas(gcf, violin_savename);
        close all;
    end
end

%% Violin plot for individual w/in and btw network values
for b = 1:v.nbands

    ANOVA_violin = [anova_dir filesep 'violin_indv'];
    if ~exist(ANOVA_violin); mkdir(ANOVA_violin); end; cd(ANOVA_violin);

    within1 = withinnetw(:,b,1,c); within2 = withinnetw(:,b,2,c);  within3 = withinnetw(:,b,3,c);
    between1 = btwnetw(:,b,1,c);   between2 = btwnetw(:,b,2,c);    between3 = btwnetw(:,b,3,c);

    NREM2_within = {'NREM2'}; NREM2_between = {'NREM2'};
    NREM3_within = {'NREM3'}; NREM3_between = {'NREM3'};
    REM_within = {'REM'}; REM_between = {'REM'};

    y1= repmat(NREM2_within, v.nsubj,1);    y2= repmat(NREM2_between, v.nsubj,1);
    y3= repmat(NREM3_within, v.nsubj,1);    y4= repmat(NREM3_between, v.nsubj,1);
    y5= repmat(REM_within, v.nsubj,1);    y6= repmat(REM_between, v.nsubj,1);

    x = [within1 between1 within2 between2 within3 between3];
%     y = [y1; y2; y3; y4; y5; y6];
    z = [repmat(1,v.nsubj,1) ; repmat(2, v.nsubj,1)];
    grouporder={NREM2_within, NREM2_between, NREM3_within, NREM3_between, REM_within, REM_between};
    figure;
    v1 = violinplot(x, y,'GroupOrder',grouporder, 'ShowMean', true);
    ylabel('connectivity (z-score)');
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'fontsize',16)
%     set(gcf, 'Position',  [10, 10, 350, 500]);
    if v.fig_title == 1
        title({['Cycle ' num2str(v.cycle) ' ' v.bandname{b} ', n = ' num2str(v.nsubj)] }, 'fontsize', 16);
    end
    if v.fig_labels == 1
        F = num2str(round(ANOVA_F(b,c))); p = round(ANOVA_p(b,c),4);
        text(0.4, 0.98,['F =' F],'Units','normalized' );
        if p > 0.001 == 1
            text(0.4, 0.96,[ 'p(UC) =  ' num2str(p)],'Units','normalized' );
        else
            text(0.4, 0.96,[ 'p(UC) < 0.001'],'Units','normalized' );
        end
    end
    violin_savename = [ANOVA_violin filesep 'c' num2str(v.cycle) '_' v.bandname{b} '_ratio_violin.png' ];
    saveas(gcf, violin_savename);
    close all;
    close all;
end
end