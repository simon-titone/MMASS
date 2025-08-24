function study2_stage_power_ANOVA(v,NETW)
%%%%% NETW = (network, network, subject, stage, cycle, band)
%%%%% BANDS = 1-delta 2-theta 3-alpha 4-sigma 5-beta 6-gamma
anova_dir = [v.homedir filesep 'results/ANOVA/interstage/'  v.t_group ];
if ~exist(anova_dir); mkdir(anova_dir); end; cd(anova_dir)
%% Power ANOVA
savepath = [v.homedir '/results/power/stage/' ];
if ~exist(savepath); mkdir(savepath); end
% PSD = nan(v.nnetw, v.nbands,v.nsubj,v.nstages,v.cycle);
for st= 1:v.nstages
    for c = 1:v.ncycles
        for i = 1:v.nsubj
            s = v.subj(i);
            file = [v.FT_dir filesep 'P' num2str(s) '_' v.stages{st} '_c' num2str(c) filesep 'FT_roi.mat'];
            if exist(file); load( file );end
            for ii = 1:21
                power_seed(ii,:,:) = FT_roi(ii,:,:).*conj(FT_roi(ii,:,:));
            end
            PSD_temp(:,:,st,c,i) = nanmean(power_seed(:,:,:),3); %calculates average across time series
            for nt = 1:v.nnetw
                PSD_netw(nt,:,st,c,i) = nanmean(PSD_temp(v.netw{nt},:,st,c,i),1); %averages across networks
            end
            for b = 1:v.nband
                PSD(:,b,st,c,i) = nanmean(PSD_netw(:,v.bandfreq{b},st,c,i),2); %averages across bands
            end
        end
    end
end
%%%%% PSD = (network, network, subject, stage, cycle)
for c = v.cycle
    for w=1:v.nnetws
        for b=1:v.nbands
            vect1 = squeeze(PSD(w,b,1,c,:));
            vect2 = squeeze(PSD(w,b,2,c,:));
            vect3 = squeeze(PSD(w,b,3,c,:));
            nsubj = length(~isnan(vect1));
            mat(:,1) = vect1;mat(:,2) = vect2;mat(:,3) = vect3;
            
            [tbl,rm] = simple_mixed_anova(mat);
            f = table2array(tbl(3,4));
            p = table2array(tbl(3,5));
            
            ANOVA_p(w,b,c) =  p; ANOVA_F(w,b,c) =  f;
        end
    end
    clear w; clear ww;
    i = 1;
    pval_vect(:,:,c) = [ANOVA_p(1,1:v.nnetws,c) ANOVA_p(2,2:v.nnetws,c) ANOVA_p(3,3:v.nnetws,c) ANOVA_p(4,4:v.nnetws,c) ANOVA_p(5,5:v.nnetws,c)  ANOVA_p(6,6:v.nnetws,c)];
    F_anovanetw(:,:,c) = [ANOVA_F(1,1:v.nnetws,c) ANOVA_F(2,2:v.nnetws,c) ANOVA_F(3,3:v.nnetws,c) ANOVA_F(4,4:v.nnetws,c) ANOVA_F(5,5:v.nnetws,c)  ANOVA_F(6,6:v.nnetws,c)];
    clear ntw;
    [~, dummy2, ~, dummy4]=fdr_bh(pval_vect(:,:,c),v.fdr_thres,'pdep','no');
    crit_p(:,:,c) = dummy2;
    adj_p(:,:,c) = dummy4;
    adj_p_mat(1,1:6,c) = adj_p(1,1:6,c); adj_p_mat(2,2:6,c) = adj_p(1,7:11,c);
    adj_p_mat(3,3:6,c) = adj_p(1,12:15,c); adj_p_mat(4,4:6,c) = adj_p(1,16:18,c);
    adj_p_mat(5,5:6,c) = adj_p(1,19:20,c); adj_p_mat(6,6,c) = adj_p(1,21,c);
    
    clear dummy2; clear dummy4;
    index = find(pval_vect(i,:,c) == crit_p);
    
    if isempty(index); crit_F = 1000; % place holder
    else; crit_F = F_anovanetw(index);
    end
    clear index;
    
    [ivect,jvect]=find(abs(ANOVA_F(:,:,c))>=abs(crit_F));
    [ivect_unc,jvect_unc]=find(ANOVA_p(:,:,c)<=0.05);
    clear crit_F;
    
    matrix = ANOVA_F(:,:,c); 
%     tri_matrix = triu(matrix);
    
    imagesc(matrix); axis square;   daspect([1 1 1]);
    set(gcf,'Color','white');set(gca,'YDir','reverse');hold on;
    dim = set(gcf, 'Position',  [1000, 900, 500, 450]); %left bottom width height
    h = colorbar; set(get(h,'title'),'string','f-val', 'fontsize', 12); %ylabel(h, 'f-val');
    ANOVA_max=max(max(ANOVA_F(:,:,c)));
    if gt(ANOVA_max, 6) ==1 && lt(ANOVA_max, 15) ==1;mm = ANOVA_max; elseif gt(ANOVA_max, 15) ==1 ; mm = 15; else; mm= 6;end;caxis([0 mm]);
    colormap(flipud(pink(10)));
    
    title({  [ v.bandname{b} ' cycle ' num2str(c) ' Stage Effect , n= ' num2str(nsubj) ]}); % generate corresponding rfx figure
    
    for a2=1:v.nnetws-1; line([a2+0.5 a2+0.5],[a2+1.5 0.5],'color','k','LineWidth',0.8);line([a2-0.5 v.nnetws+0.5],[a2+0.5 a2+0.5],'color','k','LineWidth',0.8); end
    for b2=1:v.nnetws; line([b2-0.5 b2-0.5],[0.5 v.nnetws+0.5],'color','k','LineWidth',0.001); line([0.5 v.nnetws+0.5],[b2-0.5 b2-0.5],'color','k','LineWidth',0.001);end
    clear a2; clear b2;
    
    set(gca,'YTick',[1:v.nnetws],'YTickLabel', v.netwname,'Fontsize',16);
    set(gca,'XTick',[1:v.nbands],'XTickLabel', v.bandname,'Fontsize',16);
    
    scatter(jvect_unc,ivect_unc, 80, 'w', 'o' ,'ok','LineWidth',2.5); clear ivect_unc; clear jvect_unc;
    scatter(jvect,ivect,80,'w','ok','filled' ,'LineWidth',2.5); clear ivect; clear jvect;
    
    saveas(gcf, [savepath filesep v.t_group 'stage_power_ANOVA.png'])
    clear a2; clear i; clear mm; clear f
    
    
    %% Violin plots for each significant btw session comparison
    if v.violin == 1
        ANOVA_violin = [savepath filesep v.t_group '_violin'];
        if ~exist(ANOVA_violin); mkdir(ANOVA_violin); end; cd(ANOVA_violin);
        for w=1:v.nnetws
            for b=1:v.nbands
                if ANOVA_p(w,b,c) <0.05 == 1
                    violin_band = [ANOVA_violin filesep v.bandname{b}];
                    if ~exist(violin_band); mkdir(violin_band); end
                    x = [squeeze(PSD(w,b,1,c,:)); squeeze(PSD(w,b,2,c,:)); squeeze(PSD(w,b,3,c,:))];
                    y = [repmat({'NREM2'}, v.nsubj,1) ; repmat({'NREM3'}, v.nsubj,1); repmat({'REM'}, v.nsubj,1)];
                    grouporder={'NREM2', 'NREM3', 'REM'};
                    figure; v1 = violinplot(x, y,'GroupOrder',grouporder, 'ShowMean', true);
                    a = get(gca,'XTickLabel');
                    set(gca,'XTickLabel',a,'fontsize',16)
                    set(gcf, 'Position',  [10, 10, 350, 500]);
                    if v.fig_title == 1
                        title({[ v.bandname{b} ' cycle ' num2str(c) ' ' v.netwname{w} ]}, 'FontSize', 16);
                        ylabel('Power', 'FontSize', 16);
                    end
                    if v.fig_labels == 1
                        f = num2str(round(ANOVA_F(w,b,c),3));
                        text(0.1, 0.98,['F = ' f] ,'Units','normalized');
                        p = round(ANOVA_p(w,b,c),4);
                        if p <0.001
                            text(0.1, 0.95,['p(UC) < 0.001 '],'Units','normalized');
                        else
                            text(0.1, 0.95,['p(UC) = ' num2str(p)],'Units','normalized');
                        end
                    end
                    saveas(gcf, [violin_band filesep  'cycle_' v.bandname{b} '_' v.netwname{w}   '_violin.png' ]);
                end
            end
        end
        close all;
    end
end
end