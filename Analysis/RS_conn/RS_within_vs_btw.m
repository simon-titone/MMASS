function RS_within_vs_btw(NETW,subj,n,design_matrix)
%%%%%% NETW(network, network, subject, run, band)
intra_vs_inter = [n.homedir filesep 'results' filesep 'intra_vs_inter_n' num2str(n.subj) ];
if ~exist(intra_vs_inter); mkdir(intra_vs_inter); end;


for b = 1:n.bands
    for r = 1:n.levels
        netw=NETW(:,:,:,:,b);
        for s = 1:n.subj
            for nt = 1:n.netws
                btwnetw = setdiff(1:n.netws,nt);
                intra_v_inter(nt,1,s,r,b) = NETW(nt,nt,s,r,b); %within network conn
                intra_v_inter(nt,2,s,r,b) = nanmean(netw(nt,btwnetw,s,r)); % mean of between network conn
            end
        end
    end
end

for b = 1:n.bands
    for r = 1:n.levels
        for s = 1:n.subj
            for nt = 1:n.netws
                IvI_diff(nt,1,s,b) = intra_v_inter(nt,1,s,1,b) - intra_v_inter(nt,1,s,2,b);
                IvI_diff(nt,2,s,b) = intra_v_inter(nt,2,s,1,b) - intra_v_inter(nt,2,s,2,b);
            end
        end
    end
end
for b = 1:5
    for nt = 1:n.netws
        mean_diff(nt,b,1) = mean(IvI_diff(nt,1,:,b),3);
        mean_diff(nt,b,2) = mean(IvI_diff(nt,2,:,b),3);
    end
end
for r = 1:2
figure; set(gcf,'Color','white'); box OFF;
imagesc(mean_diff(:,:,r)); axis square; 
if r ==1; a = 'within-network'; else; a = 'between-network';end;
colorbar; daspect([1 1 1]); title([ a ' conn inter-session diff' ]);
hold on; for a2=1:n.netws-1    line([a2+0.5 a2+0.5],[0.5 n.netws+0.5],'color','k');   line([0.5 n.netws+0.5],[a2+0.5 a2+0.5],'color','k');    end
colormap(n.colors_onesided);
for b2=1:n.netws    line([b2+0.5 b2+0.5],[0.5 n.netws+0.5],'color','k');   line([0.5 n.netws+0.5],[b2+0.5 b2+0.5],'color','k');    end
set(gca,'YTick',[1:n.netws],'YTickLabel',n.netwname,'Fontsize',14);
set(gca,'XTick',[1:n.netws],'XTickLabel',n.bandname,'Fontsize',14);
cd(intra_vs_inter);
end
% ANOVA of pre- vs. post-task intra vs. intra network connectivity

for b = 1:5
    for nt = 1:n.netws
        vect1=squeeze(IvI_diff(nt,1,:,b));
        vect2=squeeze(IvI_diff(nt,2,:,b));
        
        t = table(num2str(design_matrix.SessANOVA(1:n.subj)), vect1(:), vect2(:), 'VariableNames', {'Group', 'within', 'between'});
        Sess = table([1 2]', 'VariableNames', {'Sessions'});
        rm = fitrm(t,'within-between~1','WithinDesign', Sess);
        clear t; clear Sess; clear vect;
        [ranovatbl] = ranova(rm);  % Get within subject (session; session by group) effects
        pval_netw_anova_mat(nt,b) = ranovatbl{'(Intercept):Sessions','pValue'};
        F_netw_anova_mat(nt,b) = ranovatbl{'(Intercept):Sessions','F'};
        clear ranovatbl;
        clear rm;
    end
end

figure(); set(gcf,'Color','white'); box OFF;
imagesc(F_netw_anova_mat); axis square;
colorbar; daspect([1 1 1]); title(['inter vs. intra conn ANOVA' ]); % generate corresponding rfx figure
hold on; for a2=1:n.netws-1    line([a2+0.5 a2+0.5],[0.5 n.netws+0.5],'color','k');   line([0.5 n.netws+0.5],[a2+0.5 a2+0.5],'color','k');    end
clear a2;

[ivect_unc,jvect_unc]=find(pval_netw_anova_mat<=0.05);
scatter(jvect_unc,ivect_unc,36,'ok');
clear ivect_unc; clear jvect_unc;
colormap(n.colors_onesided);
for b2=1:n.netws    line([b2+0.5 b2+0.5],[0.5 n.netws+0.5],'color','k');   line([0.5 n.netws+0.5],[b2+0.5 b2+0.5],'color','k');    end
set(gca,'YTick',[1:n.netws],'YTickLabel',n.netwname,'Fontsize',14);
set(gca,'XTick',[1:n.netws],'XTickLabel',n.bandname,'Fontsize',14);
cd(intra_vs_inter);
saveas(gcf, [ 'intra_vs_inter_ANOVA.png' ]);


%%%%%% intra_v_inter(netw, within/between, subject, run, band)

%Boxplot per Network
if n.boxplot == 1
    box_dir = [n.homedir '/results/intra_v_inter_n' num2str(n.subj) '/box_n' num2str(n.subj)];
    if ~exist(box_dir); mkdir(box_dir); end; cd(box_dir);
    
    for b = 1:n.nband
        for r = 1:n.levels
            for ntw = 1:n.netws
                for i = 1:2
                    temp(:,ntw,i,r,b) = squeeze(intra_v_inter(ntw,i,:,r,b));
                end
                bplot1(:,b) = temp(:, ntw,1,r,b); bplot1 = squeeze(bplot1);
                bplot2(:,b) = temp(:, ntw,2,r,b); bplot2 = squeeze(bplot2);
                bplot = {bplot1 bplot2}; %%%FINISH ADJUSTING FOR NEW VARIABLE STRUCTURE
            end
        end
        figure;
        boxplotGroup({bplot1 bplot2},'PrimaryLabels', {'intra' 'inter'},'SecondaryLabels', {'delta','theta','alpha','beta','gamma'}, 'GroupLabelType', 'Vertical'); % only usable with R2018+
        T = mtit([n.netwname{ntw}]);
        saveas(gcf, [n.netwname{ntw} '_btw_vs_wthin_box_run' num2str(r) '_n' num2str(n.subj) '.png' ]);
    end
    
    %%%%%% temp(subj, netw, within/between, run, band)
    % Boxplot of inter vs. intra
    inter_delta = squeeze(mean(intra_v_inter_run1(:,1,:,1),1)); intra_delta = squeeze(mean(intra_v_inter_run1(:,2,:,1),1));
    inter_theta = squeeze(mean(intra_v_inter_run1(:,1,:,2),1));intra_theta = squeeze(mean(intra_v_inter_run1(:,2,:,2),1));
    inter_alpha = squeeze(mean(intra_v_inter_run1(:,1,:,3),1));intra_alpha = squeeze(mean(intra_v_inter_run1(:,2,:,3),1));
    inter_beta = squeeze(mean(intra_v_inter_run1(:,1,:,4),1));intra_beta = squeeze(mean(intra_v_inter_run1(:,2,:,4),1));
    inter_gamma = squeeze(mean(intra_v_inter_run1(:,1,:,5),1));intra_gamma = squeeze(mean(intra_v_inter_run1(:,2,:,5),1));
    inter = [inter_delta inter_theta inter_alpha inter_beta inter_gamma];
    intra = [intra_delta intra_theta intra_alpha intra_beta intra_gamma];
    figure; hold on;
    data = {intra, inter};
    
    boxplotGroup(data,'PrimaryLabels', {'intra' 'inter'},'SecondaryLabels', {'delta','theta','alpha','beta','gamma'}, 'GroupLabelType', 'Vertical'); %
    T = mtit('Inter & Intra Network Conn');
    saveas(gcf, [ 'btw_vs_wthin_box_run1_n' num2str(n.subj) '.png' ]);
    close all
end
%%%%%% intra_v_inter(netw, within/between, subject, run, band)

for b = 1:5
    for nt = 1:n.netws
        for r =1:n.levels
            %         [h,p,stats] = chi2gof(reshape(intra_v_inter(nt,1,:,r,b),[1,length(intra_v_inter_run1)]));
            %         chi2_run1(b,1,nt)= h;
            %         [h,p,stats] = chi2gof(reshape(intra_v_inter_run1(nt,2,:,b),[1,length(intra_v_inter_run1)]));
            %         chi2_run1(b,2,nt)= h;
            %
            %         [h,p,stats] = chi2gof(reshape(intra_v_inter_run2(nt,1,:,b),[1,length(intra_v_inter_run2)]));
            %         chi2_run2(b,1,nt)= h;
            %         [h,p,stats] = chi2gof(reshape(intra_v_inter_run2(nt,2,:,b),[1,length(intra_v_inter_run2)]));
            %         chi2_run2(b,2,nt)= h;
            
            [~,p,~,stats] = ttest(squeeze(intra_v_inter(nt,1,:,r,b)),squeeze(intra_v_inter(nt,2,:,r,b)));
            t_val(nt,b,r) = stats.tstat;
            p_val(nt,b,r) = p;
            
        end
    end
end

for b = 1:5
    for nt = 1:n.netws
        for r =1:n.levels
            [~, dummy2, ~, dummy4]=fdr_bh(p_val(nt,b,r),n.fdr_thres,'pdep','no');
            crit_p(nt,b,r) = dummy2;
            adj_p(:,:, r) = dummy4;
            clear dummy2; clear dummy4;
        end
    end
end

%% Scatter plots for each network/band showing participant distribution
if n.scatter == 1
    wthin_vs_btw_scatter = [n.homedir filesep 'results' filesep 'intra_vs_inter_n' num2str(n.subj) filesep 'wthin_vs_btw_scatter_n' num2str(n.subj) ];
    if ~exist(wthin_vs_btw_scatter); mkdir(wthin_vs_btw_scatter); end;
    for nt = 1:n.netws
        for b= 1:n.bands
            for l=1:n.levels
                if l == 1
                    r=1:n.subj;
                elseif l==2
                    r=(n.subj+1):(n.subj*2);
                end
                y(:,1) = transpose(squeeze(intra_v_inter(nt,1,r,b)));
                y(:,2) = transpose(squeeze(intra_v_inter(nt,2,r,b)));
                [maxconn(:,1), maxconn(:,2)] = max(conn); [minconn(:,1), minconn(:,2)] = min(conn);
                
                cut1 = std(y1) * 3; lowcut1 = mean(y1) - cut1; highcut1 = mean(y1) + cut1;
                cut2 = std(y2) * 3; lowcut2 = mean(y2) - cut2; highcut2 = mean(y2) + cut2;
                figure;hold on;
                pos1 = [0.1 0.1 0.3 0.8];
                
                subplot('Position',pos1);
                x=repmat(1,1,length(subj));
                scatter(x, y1,'o', 'b');
                title('within conn')
                xlim([.5 1.5]);
                line([0.5 1.5],[lowcut1 lowcut1],'Color','k','LineStyle','--');
                line([0.5 1.5],[highcut1 highcut1],'Color','k','LineStyle','--');
                %             ylim([0 .13]);
                labelpoints(x,y1,n.subj_label,'E', .4);
                
                pos2 = [0.5 0.1 0.3 0.8];
                subplot('Position',pos2)
                scatter(x, y2,'o', 'b');
                title('btw conn')
                xlim([.5 1.5]);
                line([0.5 1.5],[lowcut2 lowcut2],'Color','k','LineStyle','--');
                line([0.5 1.5],[highcut2 highcut2],'Color','k','LineStyle','--');
                %             ylim([0 .13]);
                labelpoints(x,y2,n.subj_label,'E', .4);
                
                bigtitle = [n.netwname{nt} ' ' n.bandname{b} ' run ' num2str(l)];
                mtit(bigtitle,'fontsize',14,'color',[0 0 1]);
                
                saveas(gcf, [n.netwname{nt} '_' n.bandname{b} '_run_' num2str(l) '.png'])
            end
        end
    end
end
cd(intra_vs_inter)
for r=1:n.levels
    %     mm =1;
    %     [ivect,jvect]=find(abs(t(:,:,r))>=abs(crit_p(:,:,r)));
    %     [ivect_unc,jvect_unc]=find(abs(t(:,:,r))>=abs(p(:,:,r)));
    
    figure(); set(gcf,'Color','white'); box OFF;
    imagesc(t_val(:,:,r)); axis square;
    colorbar; daspect([1 1 1]); title(['inter vs. intra conn run' num2str(r)]); % generate corresponding rfx figure
    hold on; for a2=1:n.netws-1    line([a2+0.5 a2+0.5],[0.5 n.netws+0.5],'color','k');   line([0.5 n.netws+0.5],[a2+0.5 a2+0.5],'color','k');    end
    clear a2;
    
    [ivect_unc,jvect_unc]=find(p_val(:,:,r)<=0.05);
    scatter(jvect_unc,ivect_unc,36,'ok');
    clear ivect_unc; clear jvect_unc;
    
    [ivect,jvect]=find(adj_p(:,:,r)<=0.05);
    scatter(jvect,ivect,9,'ok','filled');
    clear ivect; clear jvect;
    colormap('redblue');
    %         xticklabel_rotate([1:n.netws],90,n.netwname,'Fontsize',16);
    
    for b2=1:n.netws    line([b2+0.5 b2+0.5],[0.5 n.netws+0.5],'color','k');   line([0.5 n.netws+0.5],[b2+0.5 b2+0.5],'color','k');    end
    set(gca,'YTick',[1:n.netws],'YTickLabel',n.netwname,'Fontsize',14);
    set(gca,'XTick',[1:n.netws],'XTickLabel',n.bandname,'Fontsize',14);
    
    saveas(gcf, [ 'intra_vs_inter_run' num2str(r) '.png' ]);
end
close all;
end
