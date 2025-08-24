%% STEP 1: DEFINE PARAMETERS OF INTEREST - ADJUSTED BY GA FEBRUARY 2020 - FENS ABSTRACT
% ST May 2020 - created loops for behavioral and bands for 3F, 3E

clear all;
homedir = '/Users/u0127719/Google_Drive/PhD/MMASS/Experiment/Analysis/RS_conn'; %desktop
% homedir = '/Users/simontitone/Google_Drive/PhD/MMASS/Experiment/Analysis/RS_conn'; %laptop
cd(homedir);
load('matrix_connectivity.mat');
load('colors.mat');
load('bhv_vars.mat');
netwname = {'dmn','dan','van','lang','mot','vis'};
seedname = {'lANG','rANG','PCC','MPFC','lIPS','rIPS','lFEF','rFEF','rTPJ','rIFG','lTPJ','lIFG','lSMA','lCS','rCS','lS2','rS2','lV1V2','rV1V2','lMT','rMT'};
bandname = {'delta','theta','alpha','beta','gamma','bb'}; 
bhv_vars = {'offline', 'online', 'mean'};

subj = [2:11 13:16 19:29]; % full sample, not excluding behavioral issue participants
% subj = [2:11 14:16 19:29]; % n24, only excluding those without 2 RS scans & outliers on connectivity (P13)
% subj = [2 4 5 8:11 13 15 19 21:29]; %Sample with only participants who have clean BHV data
% subj = [2 4 5 8:9 11 15 19 21:29]; %Excluding P13 & P10, outliers on conn

n.subj = length(subj); % number of independent subjects
n.sess = 2; % number of sessions / RS scans per subject
n.group = 1;  % number of experimental groups (between subject factor)
n.levels = n.sess * n.group; % total number of groups and sessions
n.roi = size(seed_info,2); % number of ROI used in analyses
n.seeds = 21;
n.netws = 6;
n.bhv_vars = 3;
n.bands = length(bandname);
fdr_thres = 0.05; % note that script below assumes two-tailed tests of significance

nnetw = length(netwname);
nsubj = num2str(length(subj));
nbands = length(bandname);
for cc = 1:length(subj)
    c(cc) = subj(cc);
end
c= transpose(c);

% Create Design matrices; adviseable to check this manually as loops below
% will not work for all experimental designs
design_matrix.withinlevel=zeros(n.subj,n.levels);
design_matrix.Session1 = zeros(n.subj,1);
design_matrix.SessANOVA=zeros(n.subj,1); % 3 vectors; first specifies session effect; 2nd group and 3rd interaction;

ind=0;
level_counter = 0;
for i = 1:n.sess
    %      for ii = 1:n.group
    level_counter = level_counter + 1;
    for j=1:n.subj
        ind = ind + 1;
        design_matrix.withinlevel(ind,level_counter)=1;
        
        if i == 1
            design_matrix.Session1(ind,1) = 1;
            design_matrix.SessANOVA(ind,1) = 1;
        elseif i == 2
            design_matrix.SessANOVA(ind,1) = -1;
        end
       
    end
end
clear ind; clear level_counter; clear i; clear ii; clear j;

% % Behavioral factors - RS Final analysis NOV 2020
speed.online = on_ITI(subj);
speed.mean = mean_ITI(subj);
speed.offline = off_ITI(subj);


reds = 0:0.0005:1; j_colors = [1 1 1; ones(length(reds),1), flipud(reds'), flipud(reds'); 1 0 0];


%% STEP 2: LOAD INDIVIDUAL DATA - ADJUSTED BY GA

% LOAD DATA ONE SUBJECT AT AT TIME; STORE THE CORRELATION COEFFICIENTS IN
% MATRIX TITLED CC_zscore. IN THE Z-TRANSFORMED MATRIX, RE-WRITES ALL INF
% VALUES (I.E., CORRELATION WITHIN A SEED) TO NANs




cd([homedir '/matrix_conn/']);
counter = 1;
for l = 1:1:n.levels
    for s = 1:1:n.subj
        cd (['subj' num2str(subj(s)) '_run' num2str(l)]);
        for i = 1:size(bandname,2)%loop F
            load([bandname{i} '_netw.mat']);
            load([bandname{i} '_seeds.mat']);
            if strcmp(bandname{i},'alpha') == 1
                CC_zscore_netw_alpha(:,:,counter) =  corr_netw_band;
                CC_zscore_seed_alpha(:,:,counter) = conn_seeds_band;
            elseif strcmp(bandname{i},'beta') == 1
                CC_zscore_netw_beta(:,:,counter) =  corr_netw_band;
                CC_zscore_seed_beta(:,:,counter) = conn_seeds_band;
            elseif strcmp(bandname{i},'delta') == 1
                CC_zscore_netw_delta(:,:,counter) =  corr_netw_band;
                CC_zscore_seed_delta(:,:,counter) = conn_seeds_band;
            elseif strcmp(bandname{i},'gamma') == 1
                CC_zscore_netw_gamma(:,:,counter) =  corr_netw_band;
                CC_zscore_seed_gamma(:,:,counter) = conn_seeds_band;
            elseif strcmp(bandname{i},'theta') == 1
                CC_zscore_netw_theta(:,:,counter) =  corr_netw_band;
                CC_zscore_seed_theta(:,:,counter) = conn_seeds_band;
            elseif strcmp(bandname{i},'bb') == 1
                CC_zscore_netw_bb(:,:,counter) =  corr_netw_band;
                CC_zscore_seed_bb(:,:,counter) = conn_seeds_band;
            end
        end
        counter = counter + 1;
        cd ..
    end
end

% CC_zscore_seed_alpha(isnan(CC_zscore_seed_alpha))= 0;
% CC_zscore_seed_beta(isnan(CC_zscore_seed_beta))= 0;
% CC_zscore_seed_delta(isnan(CC_zscore_seed_delta))= 0;
% CC_zscore_seed_gamma(isnan(CC_zscore_seed_gamma))= 0;
% CC_zscore_seed_theta(isnan(CC_zscore_seed_theta))= 0;
% CC_zscore_seed_bb(isnan(CC_zscore_seed_bb))= 0;

clear ilevel; clear isub; clear corr_netw_band; clear counter; clear indices;

% %% STEP 3: WITHIN MSL NETWORK CONNECTIVITY ANALYSES - ADJUTSED BY GA
% 

%% Scatterplot of behavioral outliers
% cd /Users/simontitone/Google_Drive/PhD/MMASS/Experiment/Behavioral
% for bh=1:n.bhv_vars
%     if bh==1
%         behav_ByLevel = speed.offline;
%     elseif bh==2
%         behav_ByLevel = speed.online;
%     elseif bh==3
%         behav_ByLevel = speed.mean;
%     end
%     
%     cut = std(behav_ByLevel) * 3;
%     lowcut = mean(behav_ByLevel) - cut;
%     highcut = mean(behav_ByLevel) + cut;
%     figure
%     for i = 1:length(behav_ByLevel)
%         scatter(repmat(1,1,length(behav_ByLevel)),behav_ByLevel,'o', 'b');
%     end
%     xlim = [0.5 1.5]; ylim = [(highcut + 0.1) (lowcut - 0.1)];
%     line([0.5 1.5],[lowcut lowcut],'Color','k','LineStyle','--');
%     line([0.5 1.5],[highcut highcut],'Color','k','LineStyle','--');
%     title([bhv_vars{bh} ]);ylabel([bhv_vars{bh} ]);
%     labelpoints(1,behav_ByLevel,c,'E', .05);
%     saveas(gcf, [bhv_vars{bh} '_scatter.png' ]);    
% end
% close all;

%% Scatterplot of connectivity outliers
% conn_scatter = [homedir filesep 'results' filesep 'conn_scatter_n' nsubj];
% if ~exist(conn_scatter); mkdir(conn_scatter); end; cd(conn_scatter);
% 
% % outliers = nan(n.netws,n.netws,n.subj, n.bands, n.sess);
% 
% for b = 1:n.bands
%     if strcmp(bandname{b},'alpha') == 1; netw=CC_zscore_netw_alpha;
%     elseif strcmp(bandname{b},'beta') == 1; netw=CC_zscore_netw_beta;
%     elseif strcmp(bandname{b},'delta') == 1;netw=CC_zscore_netw_delta;
%     elseif strcmp(bandname{b},'gamma') == 1; netw=CC_zscore_netw_gamma;
%     elseif strcmp(bandname{b},'theta') == 1;netw=CC_zscore_netw_theta;
%     elseif strcmp(bandname{b},'bb') == 1;netw=CC_zscore_netw_bb;
%     end
%     for r = 1:2
%         for nt = 1:n.netws
%             for nt2 = 1:n.netws
%                 if r == 1; run = 1:n.subj; end; if r ==2; run = n.subj +1:n.subj*2;end
%                 conn = netw(nt, nt2, run); [maxconn(:,1), maxconn(:,2)] = max(conn); [minconn(:,1), minconn(:,2)] = min(conn);
%                 cut = std(conn) * 3; lowcut = mean(conn) - cut; highcut = mean(conn) + cut;
% %                 if maxconn > highcut
% %                     outliers(nt,nt2,r,b) = subj(maxconn(1,2)) ;
% %                 elseif  minconn < lowcut
% %                     outliers(nt,nt2,r,b) = minconn(1,2);
% %                 end
%                                 if maxconn > highcut || minconn < lowcut
%                                     figure
%                                     maxval(1,:) = conn(:);
%                                 for i = 1:length(conn)
%                                     scatter(repmat(1,1,length(conn)),conn,'o', 'b');
%                                 end
%                                 xlim = [0.5 1.5]; ylim = [(highcut + 0.1) (lowcut - 0.1)];
%                                 line([0.5 1.5],[lowcut lowcut],'Color','k','LineStyle','--');
%                                 line([0.5 1.5],[highcut highcut],'Color','k','LineStyle','--');
%                                 title([bandname{b} ' '  netwname{nt} ' ' netwname{nt2} ' run ' num2str(r)  ' connectivity' ]);
%                                 ylabel('connectivity (z-score)');
%                                 labelpoints(1,conn,c,'E', .05);
%                                 saveas(gcf, [bandname{b} '_'  netwname{nt} '_' netwname{nt2} '_run' num2str(r)  '_conn_scatter.png' ]);
%                   end
%             end
%         end
%     end
% end
% clear r; clear nt; clear nt2; clear b;
% close all;
% cd([homedir '/results']);
% outliervar = ['RS_conn_outliers_n' num2str(n.subj)];
% save(outliervar, 'outliers');


%% T-TEST OF DIFFERENCE BETWEEN SESSIONS
% sesseffect = [homedir filesep 'results' filesep 'sess_n' nsubj];
% if ~exist(sesseffect); mkdir(sesseffect); end; cd(sesseffect);
% 
% diff = zeros(n.netws, n.netws, n.subj,n.bands);
% for b = 1:n.bands
%     if strcmp(bandname{b},'alpha') == 1; netw=CC_zscore_netw_alpha;
%     elseif strcmp(bandname{b},'beta') == 1; netw=CC_zscore_netw_beta;
%     elseif strcmp(bandname{b},'delta') == 1;netw=CC_zscore_netw_delta;
%     elseif strcmp(bandname{b},'gamma') == 1; netw=CC_zscore_netw_gamma;
%     elseif strcmp(bandname{b},'theta') == 1;netw=CC_zscore_netw_theta;
%     elseif strcmp(bandname{b},'bb') == 1;netw=CC_zscore_netw_bb;
%     end
%     diff(:,:,1:n.subj,b) = minus(netw(:,:,1:n.subj), netw(:,:,n.subj+1:n.subj*n.levels));
% end
% 
% for b = 1:n.bands
%     for nt = 1:n.netws
%         for nt2 = 1:n.netws
%             [h,p,ci,stats] = ttest(diff(nt, nt2, :, b));
%             diff_h(nt, nt2,b) = h;
%             diff_p(nt, nt2,b) = p;
%             diff_t(nt, nt2,b) = stats.tstat;
%         end
%     end
% end
% 
% % % % FDR Correction
% % for b = 1:n.bands
% %     for nt = 1:nnetw
% %         [~, dummy2, ~, dummy4]=fdr_bh(run1_p(nt,b),fdr_thres,'pdep','no');
% %         crit_p_run1(nt,b) = dummy2;
% %         adj_p_run1(:,:) = dummy4;
% %         clear dummy2; clear dummy4;
% %         
% %         [~, dummy2, ~, dummy4]=fdr_bh(run2_p(nt,b),fdr_thres,'pdep','no');
% %         crit_p_run2(nt,b) = dummy2;
% %         adj_p_run2(:,:) = dummy4;
% %         clear dummy2; clear dummy4;
% %     end
% % end
% 
% for b = 1:n.bands
%     figure(); set(gcf,'Color','white'); box OFF;
%     imagesc(diff_t(:,:,b)); axis square; %caxis([-mm mm]);
%     colorbar; daspect([1 1 1]); reds = 0:0.0005:1;
% %     j_colors = [1 1 1; ones(length(reds),1), flipud(reds'), flipud(reds'); 1 0 0];
% %     colormap((j_colors))
%     
%     hold on; for b2=1:n.netws    line([b2+0.5 b2+0.5],[0.5 n.netws+0.5],'color','k');   line([0.5 n.netws+0.5],[b2+0.5 b2+0.5],'color','k');    end
%     set(gca,'YTick',[1:n.netws],'YTickLabel',netwname,'Fontsize',16);
%     set(gca,'XTick',[1:n.netws],'XTickLabel',netwname,'Fontsize',16);
%     title(['Session Effect ' bandname{b} ' Netw n= ' nsubj ]);
%     
%      [ivect_unc,jvect_unc]=find(diff_p(:,:,b)<=0.05);
%     scatter(jvect_unc,ivect_unc,36,'ok');
%     clear ivect_unc; clear jvect_unc;
%     
% %     [ivect,jvect]=find(p_corr<=0.05);
% %     scatter(jvect,ivect,9,'ok','filled');
% %     clear ivect; clear jvect;
% %     
%     
%     saveas(gcf, ['session_effect_' bandname{b} '_netw.png']);
% end

%% INDV PARTICIPANT MATRICES
% if ~exist([homedir filesep 'results' filesep 'indv'])
%     mkdir([homedir filesep 'results' filesep 'indv']);
% end
% cd('../results/indv')
% 
% for s = 1:n.subj
%     for b = 1:n.bands
%         
%     if strcmp(bandname{b},'alpha') == 1; netw=CC_zscore_netw_alpha;
%     elseif strcmp(bandname{b},'beta') == 1; netw=CC_zscore_netw_beta;
%     elseif strcmp(bandname{b},'delta') == 1;netw=CC_zscore_netw_delta;
%     elseif strcmp(bandname{b},'gamma') == 1; netw=CC_zscore_netw_gamma;
%     elseif strcmp(bandname{b},'theta') == 1;netw=CC_zscore_netw_theta;
%     elseif strcmp(bandname{b},'bb') == 1;netw=CC_zscore_netw_bb;
%     end
%     
%     matrix = netw(:,:,s);
%     min1 = min(min(matrix,2));
%     indv_min(s,b) = min(min1);
%     
%     max1 = max(max(matrix));
%     indv_max(s,b,1) = max(max1);
%     end
% end
% 

%% STEP 3A: Within-level means
% COMPUTES MEANS OF THE  Z-TRANSFORMED CORRELATION
% COEFFICIENTS WITHIN THE MSL NETWORK AND WITHIN EACH LEVEL (E.G.,
% SESSION AND/OR GROUP) AND GENERATES A CORRESPONDING FIGURE.
% 
% % % NETWORKS
% groupmean_netw_alpha_run1 = nanmean(CC_zscore_netw_alpha(:,:,1:n.subj),3);
% groupmean_netw_alpha_run2 = nanmean(CC_zscore_netw_alpha(:,:,n.subj+1:n.subj*2),3);
% groupmean_netw_beta_run1 = nanmean(CC_zscore_netw_beta(:,:,1:n.subj),3);
% groupmean_netw_beta_run2 = nanmean(CC_zscore_netw_beta(:,:,n.subj+1:n.subj*2),3);
% groupmean_netw_delta_run1 = nanmean(CC_zscore_netw_delta(:,:,1:n.subj),3);
% groupmean_netw_delta_run2 = nanmean(CC_zscore_netw_delta(:,:,n.subj+1:n.subj*2),3);
% groupmean_netw_theta_run1 = nanmean(CC_zscore_netw_theta(:,:,1:n.subj),3);
% groupmean_netw_theta_run2 = nanmean(CC_zscore_netw_theta(:,:,n.subj+1:n.subj*2),3);
% groupmean_netw_gamma_run1 = nanmean(CC_zscore_netw_gamma(:,:,1:n.subj),3);
% groupmean_netw_gamma_run2 = nanmean(CC_zscore_netw_gamma(:,:,n.subj+1:n.subj*2),3);
% groupmean_netw_bb_run1 = nanmean(CC_zscore_netw_bb(:,:,1:n.subj),3);
% groupmean_netw_bb_run2 = nanmean(CC_zscore_netw_bb(:,:,n.subj+1:n.subj*2),3);
% 
% 
% % % SEEDS
% groupmean_seed_alpha_run1 = nanmean(CC_zscore_seed_alpha(:,:,1:n.subj),3);
% groupmean_seed_alpha_run2 = nanmean(CC_zscore_seed_alpha(:,:,n.subj+1:n.subj*2),3);
% groupmean_seed_beta_run1 = nanmean(CC_zscore_seed_beta(:,:,1:n.subj),3);
% groupmean_seed_beta_run2 = nanmean(CC_zscore_seed_beta(:,:,n.subj+1:n.subj*2),3);
% groupmean_seed_delta_run1 = nanmean(CC_zscore_seed_delta(:,:,1:n.subj),3);
% groupmean_seed_delta_run2 = nanmean(CC_zscore_seed_delta(:,:,n.subj+1:n.subj*2),3);
% groupmean_seed_theta_run1 = nanmean(CC_zscore_seed_theta(:,:,1:n.subj),3);
% groupmean_seed_theta_run2 = nanmean(CC_zscore_seed_theta(:,:,n.subj+1:n.subj*2),3);
% groupmean_seed_gamma_run1 = nanmean(CC_zscore_seed_gamma(:,:,1:n.subj),3);
% groupmean_seed_gamma_run2 = nanmean(CC_zscore_seed_gamma(:,:,n.subj+1:n.subj*2),3);
% groupmean_seed_bb_run1 = nanmean(CC_zscore_seed_bb(:,:,1:n.subj),3);
% groupmean_seed_bb_run2 = nanmean(CC_zscore_seed_bb(:,:,n.subj+1:n.subj*2),3);
% 
% % example figure
% groupmeans = [homedir filesep 'results' filesep 'group_means_n' nsubj];
% if ~exist(groupmeans); mkdir(groupmeans); end; cd(groupmeans)
% 
% for b = 1:n.bands
%     for r = 1:n.levels
%         if strcmp(bandname{b},'alpha') == 1 && r==1
%             matrix = groupmean_netw_alpha_run1;
%         elseif strcmp(bandname{b},'alpha') == 1 && r==2
%             matrix = groupmean_netw_alpha_run2;
%         elseif strcmp(bandname{b},'beta') == 1 && r==1
%             matrix = groupmean_netw_beta_run1;
%         elseif strcmp(bandname{b},'beta') == 1 && r==2
%             matrix = groupmean_netw_beta_run2;
%         elseif strcmp(bandname{b},'delta') == 1 && r==1
%             matrix = groupmean_netw_delta_run1;
%         elseif strcmp(bandname{b},'delta') == 1 && r==2
%             matrix = groupmean_netw_delta_run2;
%         elseif strcmp(bandname{b},'theta') == 1 && r==1
%             matrix = groupmean_netw_theta_run1;
%         elseif strcmp(bandname{b},'theta') == 1 && r==2
%             matrix = groupmean_netw_theta_run2;
%         elseif strcmp(bandname{b},'gamma') == 1 && r==1
%             matrix = groupmean_netw_gamma_run1;
%         elseif strcmp(bandname{b},'gamma') == 1 && r==2
%             matrix = groupmean_netw_gamma_run2;
%         elseif strcmp(bandname{b},'bb') == 1 && r==1
%             matrix = groupmean_netw_bb_run1;
%         elseif strcmp(bandname{b},'bb') == 1 && r==2
%             matrix = groupmean_netw_bb_run2;
%         end
%     
%         figure(); set(gcf,'Color','white'); box OFF;
%         imagesc(matrix); axis square; %caxis([-mm mm]);
%         colorbar; daspect([1 1 1]); %colormap autumn;
%         
%  
%         colormap((j_colors))
%         
%         hold on;
%         for b2=1:n.netws    line([b2+0.5 b2+0.5],[0.5 n.netws+0.5],'color','k');   line([0.5 n.netws+0.5],[b2+0.5 b2+0.5],'color','k');    end
%         set(gca,'YTick',[1:n.netws],'YTickLabel',netwname,'Fontsize',16);
%         set(gca,'XTick',[1:n.netws],'XTickLabel',netwname,'Fontsize',16);
%         title(['Group Means ' bandname{b} ' Netw run' num2str(r) ' n= ' nsubj ]);
% 
% %         [~,x] = max(matrix,[],2); % max of each row
% %         plot(x,1:n.netws,'ko','MarkerSize',8,'MarkerFaceColor','k');
% %         [~,x] = max(matrix); % max of each column
% %         plot(1:n.netws,x,'ko','MarkerSize',8,'MarkerFaceColor','k');
% %         
%         saveas(gcf, ['Group_means_' bandname{b} '_netw_run_' num2str(r) '.png']);
% %         
% %         figure(); set(gcf,'Color','white'); box OFF;
% %         imagesc(CC_zscore_seed_alpha_run1); axis square; %caxis([-mm mm]);
% %         colormap(colors); colorbar; daspect([1 1 1]);
% %         title(['Group means alpha seed run' num2str(r)]);
% %         hold on;
% %         for b2=1:n.seeds    line([b2+0.5 b2+0.5],[0.5 n.seeds+0.5],'color','k');   line([0.5 n.seeds+0.5],[b2+0.5 b2+0.5],'color','k');    end
% %         set(gca,'YTick',[1:n.seeds],'YTickLabel',seedname,'Fontsize',16);
% %         set(gca,'XTick',[1:n.seeds],'XTickLabel',seedname,'Fontsize',16);
% %         set(gca,'XTickLabelRotation',45)
% %                 saveas(gcf, ['Group_means_' bandname{i} '_seed_run_' num2str(r) '.png']);
% 
%         
%     end
% end
% 
% close all;


%__________________________________________________________________________
%% STEP 3B: Baseline Means
% % SAME AS ABOVE BUT AVERAGING ACROSS THE TWO EXPERIMENTAL GROUPS FOR THE
% % FIRST SESSION (I.E., BASELINE RS SCAN)
% % ST- Removed because it isn't relevent to this analysis
% 

%% Within Vs. Between network conn per frequency band
intra_v_inter = zeros(6,2);
for b = 1:n.bands
    if strcmp(bandname{b},'alpha') == 1; netw=CC_zscore_netw_alpha;
    elseif strcmp(bandname{b},'beta') == 1; netw=CC_zscore_netw_beta;
    elseif strcmp(bandname{b},'delta') == 1;netw=CC_zscore_netw_delta;
    elseif strcmp(bandname{b},'gamma') == 1; netw=CC_zscore_netw_gamma;
    elseif strcmp(bandname{b},'theta') == 1;netw=CC_zscore_netw_theta;
    elseif strcmp(bandname{b},'bb') == 1;netw=CC_zscore_netw_bb;
    end
    for s = 1:n.subj * n.levels
        for nt = 1:nnetw
            btwnetw = setdiff(1:nnetw,nt);
            intra_v_inter(nt,1,s,b) = netw(nt,nt,s);
            intra_v_inter(nt,2,s,b) = nanmean(netw(nt,btwnetw,s));
        end
    end
%     diff(:,:,b) = intra_v_inter(:,1,1:25,b) - intra_v_inter(:,2,1:25,b);
%     mean_diff(:,b) = nanmean(diff(:,:,b), 2);
%     std_diff(:,b) = std(diff(:,:,b),1,2);
end
intra_v_inter_run1 = intra_v_inter(:,:,1:n.subj,:); intra_v_inter_run2 = intra_v_inter(:,:,n.subj+1:n.subj*n.levels,:);
% 
% % Boxplot of inter vs. intra
% cd([homedir '/results/intra_v_inter'])
% % for i = 1:2
% %     for b = 1:5
% %         for ntw = 1:6
% %             temp(:,b,ntw,i) = squeeze(intra_v_inter_run1(ntw,i,:, b));
% %         end
% %     end
% % end
% % 
% % for n = 1:n.netws
% %     figure;
% %     bplot = {temp(:,:, n, 1) temp(:,:, n, 2) };
% %     boxplotGroup(bplot,'PrimaryLabels', {'intra' 'inter'},'SecondaryLabels', {'delta','theta','alpha','beta','gamma'}, 'GroupLabelType', 'Vertical'); % only usable with R2020
% %     t = mtit([netwname{n}]);
% %     saveas(gcf, [netwname{n} '_btw_vs_wthin_box_run1_n' nsubj '.png' ]);
% % end
% 
% 
% inter_delta = squeeze(mean(intra_v_inter_run1(:,1,:,1),1)); intra_delta = squeeze(mean(intra_v_inter_run1(:,2,:,1),1));
% inter_theta = squeeze(mean(intra_v_inter_run1(:,1,:,2),1));intra_theta = squeeze(mean(intra_v_inter_run1(:,2,:,2),1));
% inter_alpha = squeeze(mean(intra_v_inter_run1(:,1,:,3),1));intra_alpha = squeeze(mean(intra_v_inter_run1(:,2,:,3),1));
% inter_beta = squeeze(mean(intra_v_inter_run1(:,1,:,4),1));intra_beta = squeeze(mean(intra_v_inter_run1(:,2,:,4),1));
% inter_gamma = squeeze(mean(intra_v_inter_run1(:,1,:,5),1));intra_gamma = squeeze(mean(intra_v_inter_run1(:,2,:,5),1));
% inter = [inter_delta inter_theta inter_alpha inter_beta inter_gamma];
% intra = [intra_delta intra_theta intra_alpha intra_beta intra_gamma];
% figure; hold on;
% data = {intra, inter};
% 
% boxplotGroup(data,'PrimaryLabels', {'intra' 'inter'},'SecondaryLabels', {'delta','theta','alpha','beta','gamma'}, 'GroupLabelType', 'Vertical'); %
% t = mtit('Inter & Intra Network Conn');
% saveas(gcf, [ 'btw_vs_wthin_box_run1_n' nsubj '.png' ]);
% close all
% 
% boxplot([delta theta alpha beta gamma], bandname(1:5), 'colors','r','width',0.18 )
% boxplot(intra, bandname(1:5),'colors','b','positions',pos_intra,'width',0.18 )

for b = 1:n.bands
    for nt = 1:nnetw
        
        [h,p,stats] = chi2gof(reshape(intra_v_inter_run1(nt,1,:,b),[1,length(intra_v_inter_run1)]));
        chi2_run1(b,1,nt)= h;
        [h,p,stats] = chi2gof(reshape(intra_v_inter_run1(nt,2,:,b),[1,length(intra_v_inter_run1)]));
        chi2_run1(b,2,nt)= h;
        
        [h,p,stats] = chi2gof(reshape(intra_v_inter_run2(nt,1,:,b),[1,length(intra_v_inter_run2)]));
        chi2_run2(b,1,nt)= h;
        [h,p,stats] = chi2gof(reshape(intra_v_inter_run2(nt,2,:,b),[1,length(intra_v_inter_run2)]));
        chi2_run2(b,2,nt)= h;
        
        [h,p,ci,stats] = ttest(intra_v_inter_run1(nt,1,:,b),intra_v_inter_run1(nt,2,:,b));
        run1_t(nt,b) = stats.tstat;
        run1_p(nt,b) = p;
        
        [h,p,ci,stats] = ttest(intra_v_inter_run2(nt,1,:,b),intra_v_inter_run2(nt,2,:,b));
        run2_t(nt,b) = stats.tstat;
        run2_p(nt,b) = p;
    end
end

for b = 1:n.bands
    for nt = 1:nnetw
        [~, dummy2, ~, dummy4]=fdr_bh(run1_p(nt,b),fdr_thres,'pdep','no');
        crit_p_run1(nt,b) = dummy2;
        adj_p_run1(:,:) = dummy4;
        clear dummy2; clear dummy4;
        
        [~, dummy2, ~, dummy4]=fdr_bh(run2_p(nt,b),fdr_thres,'pdep','no');
        crit_p_run2(nt,b) = dummy2;
        adj_p_run2(:,:) = dummy4;
        clear dummy2; clear dummy4;
    end
end

wthin_vs_btw_scatter = [homedir filesep 'results' filesep 'wthin_vs_btw_scatter_n' nsubj ];
if ~exist(wthin_vs_btw_scatter); mkdir(wthin_vs_btw_scatter); end; cd(wthin_vs_btw_scatter);
% 
% %% Scatter plots for each network/band showing participant distribution
% for nt = 1:n.netws
%     for b= 1:n.bands
%         for l=1:n.levels
%             if l == 1
%                 r=1:25;
%             elseif l==2
%                 r=26:50;
%             end
%             for i = 1:length(subj)
%                 c(i) = subj(i);
%                 y1 = intra_v_inter(nt,1,r,b); y1 = transpose(y1(:));
%                 y2 = intra_v_inter(nt,2,r,b); y2 = transpose(y2(:));
%             end
%             figure;hold on;
%             pos1 = [0.1 0.1 0.3 0.8];
%             subplot('Position',pos1);
%             x=repmat(1,1,length(subj));
%             scatter(x, y1,'o', 'b');
%             title('within conn')
%             xlim([.5 1.5]);
% %             ylim([0 .13]);
%             labelpoints(x,y1,c,'E', .4);
%             
%             pos2 = [0.5 0.1 0.3 0.8];
%             subplot('Position',pos2)
%             scatter(x, y2,'o', 'b');
%             title('btw conn')
%             xlim([.5 1.5]);
% %             ylim([0 .13]);
%             labelpoints(x,y2,c,'E', .4);
%             
%             bigtitle = [netwname{nt} ' ' bandname{b} ' run ' num2str(l)];
%             mtit(bigtitle,'fontsize',14,'color',[0 0 1]);
%             
%             saveas(gcf, [netwname{nt} '_' bandname{b} '_run_' num2str(l) '.png'])
%         end
%     end
% end

cd([homedir filesep 'results'])
for r=1:n.levels
    if r == 1
        tstat= run1_t;
        p_uncorr = run1_p;
        p_corr = adj_p_run1;
    elseif r==2
        tstat= run2_t;
        p_uncorr = run2_p;
        p_corr = adj_p_run2;
    end
 
 
%     mm =1;
    [ivect,jvect]=find(abs(run1_t(:,:))>=abs(crit_p_run1));
    [ivect_unc,jvect_unc]=find(abs(run1_t(:,:))>=abs(run1_p));
    
    figure(); set(gcf,'Color','white'); box OFF;
    imagesc(tstat); axis square;   
    colorbar; daspect([1 1 1]); title(['inter vs. intra conn run' num2str(r)]); % generate corresponding rfx figure
    hold on; for a2=1:n.netws-1    line([a2+0.5 a2+0.5],[0.5 n.netws+0.5],'color','k');   line([0.5 n.netws+0.5],[a2+0.5 a2+0.5],'color','k');    end
    clear a2;
   
%          
    [ivect_unc,jvect_unc]=find(p_uncorr<=0.05);
    scatter(jvect_unc,ivect_unc,36,'ok');
    clear ivect_unc; clear jvect_unc;
    
    [ivect,jvect]=find(p_corr<=0.05);
    scatter(jvect,ivect,9,'ok','filled');
    clear ivect; clear jvect;
    colormap('redblue');
    %         xticklabel_rotate([1:n.netws],90,netwname,'Fontsize',16);
    
    for b2=1:n.netws    line([b2+0.5 b2+0.5],[0.5 n.netws+0.5],'color','k');   line([0.5 n.netws+0.5],[b2+0.5 b2+0.5],'color','k');    end
    set(gca,'YTick',[1:n.netws],'YTickLabel',netwname,'Fontsize',14);
    set(gca,'XTick',[1:n.netws],'XTickLabel',bandname,'Fontsize',14);

    saveas(gcf, [ 'intra_vs_inter_run' num2str(r) '.png' ]);
end
close all;

%% STEP 3C: Within-level Statistics - ADJUSTED BY GA
% % % Within each level (i.e., group/session combo), fit z-transformed
% % % correlation values with GLM. Use corresponding BETA and SE estimates to
% % % compute a t-score. Write the BETA  estimates and t-scores to text files
% % % (ffx and rfx, respectively) and generate corresponding plots
% 
% betaseeds.level_mat=zeros(1,n.seeds,n.levels);
% betanetwks.level_mat=zeros(1,n.netws,n.levels);
% seseeds.level_mat=zeros(1,n.seeds,n.levels);
% senetwks.level_mat=zeros(1,n.netws,n.levels);
% 
% withinlevel = [homedir filesep 'results' filesep 'within_level_n' nsubj ];
% if ~exist(withinlevel); mkdir(withinlevel); end; cd(withinlevel);
% 
% for i= 1:n.levels
%     for b = 1:n.bands
%         %     % %Seeds
%         %     for w=1:n.seeds
%         %         for ww=1:n.seeds
%         %             if strcmp(bandname{b},'alpha') == 1
%         %                 vect=squeeze(CC_zscore_seed_alpha(w,ww,:));
%         %             elseif strcmp(bandname{b},'beta') == 1
%         %                 vect=squeeze(CC_zscore_seed_beta(w,ww,:));
%         %             elseif strcmp(bandname{b},'delta') == 1
%         %                 vect=squeeze(CC_zscore_seed_delta(w,ww,:));
%         %             elseif strcmp(bandname{b},'gamma') == 1
%         %                 vect=squeeze(CC_zscore_seed_gamma(w,ww,:));
%         %             elseif strcmp(bandname{b},'theta') == 1
%         %                 vect=squeeze(CC_zscore_seed_theta(w,ww,:));
%         %             elseif strcmp(bandname{b},'bb') == 1
%         %                 vect=squeeze(CC_zscore_seed_bb(w,ww,:));
%         %             end
%         %             [~,~,STATS]=glmfit([design_matrix.withinlevel],vect,'normal','constant','off');
%         %             betaseeds.level_mat(w,ww,:)=STATS.beta;
%         %             seseeds.level_mat(w,ww,:)=STATS.se;
%         %         end
%         %     end
%         %     clear w; clear ww; clear STATS; clear vect;
%         
%         % %Netws
%         for w=1:n.netws
%             for ww=1:n.netws
%                 if strcmp(bandname{b},'alpha') == 1
%                     vect=squeeze(CC_zscore_netw_alpha(w,ww,:)); % vector of correlations between specific ROIs (as specified by  w; length of vect = num of subjects)
%                 elseif strcmp(bandname{b},'beta') == 1
%                     vect=squeeze(CC_zscore_netw_beta(w,ww,:));
%                 elseif strcmp(bandname{b},'delta') == 1
%                     vect=squeeze(CC_zscore_netw_delta(w,ww,:));
%                 elseif strcmp(bandname{b},'gamma') == 1
%                     vect=squeeze(CC_zscore_netw_gamma(w,ww,:));
%                 elseif strcmp(bandname{b},'theta') == 1
%                     vect=squeeze(CC_zscore_netw_theta(w,ww,:));
%                 elseif strcmp(bandname{b},'bb') == 1
%                     vect=squeeze(CC_zscore_netw_bb(w,ww,:));
%                 end
%                 [~,~,STATS]=glmfit([design_matrix.withinlevel],vect,'normal','constant','off');
%                 betanetws.level_mat(w,ww,:)=STATS.beta;
%                 senetws.level_mat(w,ww,:)=STATS.se;
%             end
%         end
%         clear w; clear ww; clear STATS; clear vect;
%         
%         for i=1:n.levels
%             
%             %             dataseeds_ave=squeeze(betaseeds.level_mat(:,:,i));
%             %             dataseeds_se=squeeze(seseeds.level_mat(:,:,i));
%             %             dataseeds_tscore=dataseeds_ave./dataseeds_se;
%             %
%             datanetws_ave=squeeze(betanetws.level_mat(:,:,i));
%             datanetws_se=squeeze(senetws.level_mat(:,:,i));
%             datanetws_tscore=datanetws_ave./datanetws_se;
%             
%             %             pathx=[homedir '/matrix/' 'Run' num2str(i) '_allsubj_ffx_seeds.txt'];
%             %             save(pathx,'dataseeds_ave','-ascii');                            % write text file (all subjects within a group) labeled ffx with beta estimates
%             
% %             pathx=[homedir '/matrix/' 'Run' num2str(i) '_allsubj_ffx_netws.txt'];
% %             save(pathx,'datanetws_ave','-ascii');                            % write text file (all subjects within a group) labeled ffx with beta estimates
% %             clear pathx;
%             
%             %             mm=7;
%             %             figure(); set(gcf,'Color','white'); box OFF;
%             %             imagesc(data_ave); axis square; caxis([-mm mm]); colormap(colors); colorbar; daspect([1 1 1]); title([ bandname{b} ' MSL: Level' num2str(i) ' - FFX']); % generate corresponding ffx figure
%             %             hold on; for a2=1:n.roi_MSL-1    line([a2+0.5 a2+0.5],[0.5 n.roi_MSL+0.5],'color','k');   line([0.5 n.roi_MSL+0.5],[a2+0.5 a2+0.5],'color','k');    end
%             %             clear a2;
%             %             set(gca,'YTick',[1],'YTickLabel',Label.Seed,'Fontsize',16);
%             %             xticklabel_rotate([1:n.roi_MSL],90,Label.MSL,'Fontsize',16);
%             % %             print([homedir '/matrix/' 'MSL: group' num2str(i) 'allsubj_ffx.tif'],'-dtiff','-r150'); close;
%             % %
%             %             pathx=[homedir '/matrix/' 'MSL_level' num2str(i) '_allsubj_rfx.txt'];          % write text file (all subjects within a group) labeled as rfx with t-scores
%             %             save(pathx,'data_tscore','-ascii'); clear pathx;
%             
%             %             %seeds
%             %             count = 1;
%             %             for nseed_lin = 1:1:n.seeds
%             %                 for nseed_col = 1:1:n.seeds
%             %                     pval.withinLevelseeds(i, count) = 2*tcdf(abs(dataseeds_tscore(nseed_lin,nseed_col)),n.subj-1, 'upper'); % 2 -tailed p value
%             %                     %             pval.withinLevelseeds(i, count) = 2*tcdf(abs(dataseeds_tscore(nseed_lin,nseed_col)),n.subj-1); % 2 -tailed p value
%             %                     t.withinLevelseeds(i,count) = dataseeds_tscore(nseed_lin,nseed_col);
%             %                     count = count + 1;
%             %                 end
%             %             end
%             %             clear count;clear nseed_col;clear nseed_lin;
%             %
%             %             [~, dummy2, ~, dummy4]=fdr_bh(pval.withinLevelseeds(i,:),fdr_thres,'pdep','no');
%             %             crit_p.withinLevelseeds(i) = dummy2;
%             %             adj_p.withinLevelseeds(i,1:length(dummy4)) = dummy4;
%             %             clear dummy2; clear dummy4;
%             %
%             %             index = find(pval.withinLevelseeds(i,:) == crit_p.withinLevelseeds(i));
%             %             crit_t = t.withinLevelseeds(i,index(1));
%             %             crit_t_unc = tinv(fdr_thres/2,n.subj-1); % 2 tailed
%             %             clear index;
%             
%             %%netws
%             count = 1;
%             for nnetw_lin = 1:1:n.netws
%                 for nnetw_col = 1:1:n.netws
%                     pval.withinLevelnetws(i, count) = 2*tcdf(abs(datanetws_tscore(nnetw_lin,nnetw_col)),n.subj-1, 'upper'); % 2 -tailed p value
%                     %                 pval.withinLevelnetws(i, count) = 2*tcdf(abs(datanetws_tscore(nnetw_lin,nnetw_col)),n.subj-1); % 2 -tailed p value
%                     t.withinLevelnetws(i,count) = datanetws_tscore(nnetw_lin,nnetw_col);
%                     count = count + 1;
%                 end
%             end
%             clear count;clear nnetw_col;clear nnetw_lin
%             
%             pval.withinLevelnetws_4cor = [pval.withinLevelnetws(:,1:n.netws) pval.withinLevelnetws(:,n.netws+2:n.netws*2) pval.withinLevelnetws(:,n.netws*2+3:n.netws*3) pval.withinLevelnetws(:,n.netws*3+4:n.netws*4) pval.withinLevelnetws(:,n.netws*4+5:n.netws*5) pval.withinLevelnetws(:,n.netws*5+6:n.netws*6)];
%             t.withinLevelnetws_4cor = [t.withinLevelnetws(:,1:n.netws) t.withinLevelnetws(:,n.netws+2:n.netws*2) t.withinLevelnetws(:,n.netws*2+3:n.netws*3) t.withinLevelnetws(:,n.netws*3+4:n.netws*4) t.withinLevelnetws(:,n.netws*4+5:n.netws*5) t.withinLevelnetws(:,n.netws*5+6:n.netws*6)];
%             
%             [~, dummy2, ~, dummy4]=fdr_bh(pval.withinLevelnetws_4cor(i,:),fdr_thres,'pdep','no');
%             crit_p.withinLevelnetws(i) = dummy2;
%             adj_p.withinLevelnetws(i,1:length(dummy4)) = dummy4;
%             clear dummy2; clear dummy4;
%             
%             index = find(pval.withinLevelnetws_4cor(i,:) == crit_p.withinLevelnetws(i));
%             crit_t = t.withinLevelnetws(i,index(1));
%             crit_t_unc = tinv(fdr_thres/2,n.subj-1); % 2 tailed
%             clear index;
%             
%             %seeds
%             
%             %             mm=7;
%             %             [ivect,jvect]=find(abs(dataseeds_tscore(:,:))>=abs(crit_t));
%             %             [ivect_unc,jvect_unc]=find(abs(dataseeds_tscore(:,:))>=abs(crit_t_unc));
%             %
%             %             figure(); set(gcf,'Color','white'); box OFF;
%             %             imagesc(dataseeds_tscore); axis square; caxis([-mm mm]); colormap(colors); colorbar; daspect([1 1 1]); title([ bandname{b} ' Seed: level' num2str(i) ' - RFX']); % generate corresponding rfx figure
%             %             hold on; for a2=1:n.seeds-1    line([a2+0.5 a2+0.5],[0.5 n.seeds+0.5],'color','k');   line([0.5 n.seeds+0.5],[a2+0.5 a2+0.5],'color','k');    end
%             %             clear a2;
%             %             set(gca,'YTick',[1],'YTickLabel',seedname,'Fontsize',16);
%             %             scatter(jvect_unc,ivect_unc,36,'ok');
%             %             clear ivect_unc; clear jvect_unc;
%             %
%             %             scatter(jvect,ivect,9,'ok','filled');
%             %             clear ivect; clear jvect;
%             %
%             %             xticklabel_rotate([1:n.seeds],90,seedname,'Fontsize',16);
%             %             print([homedir '/matrix/' 'level' num2str(i) '_seeds_allsubj_rfx.tif'],'-dtiff','-r150');
%             %
%             
%             %%netw
% %             mm =10;
%             [ivect,jvect]=find(abs(datanetws_tscore(:,:))>=abs(crit_t));
%             [ivect_unc,jvect_unc]=find(abs(datanetws_tscore(:,:))>=abs(crit_t_unc));
%             
%             figure(); set(gcf,'Color','white'); box OFF;
%             imagesc(datanetws_tscore); axis square;  %caxis([-mm mm]);
%             colormap(j_colors); 
%             colorbar; daspect([1 1 1]); 
%             title([ bandname{b} ' Ntw: within level' num2str(i) ' n= ' nsubj ]); % generate corresponding rfx figure
%             hold on; for a2=1:n.netws-1    line([a2+0.5 a2+0.5],[0.5 n.netws+0.5],'color','k');   line([0.5 n.netws+0.5],[a2+0.5 a2+0.5],'color','k');    end
%             clear a2;
%             %         set(gca,'YTick',[1],'YTickLabel',netwname,'Fontsize',16);
%             scatter(jvect_unc,ivect_unc,36,'ok');
%             clear ivect_unc; clear jvect_unc;
%             
%             scatter(jvect,ivect,9,'ok','filled');
%             clear ivect; clear jvect;
%             
%             %         xticklabel_rotate([1:n.netws],90,netwname,'Fontsize',16);
%             
%             for b2=1:n.netws    line([b2+0.5 b2+0.5],[0.5 n.netws+0.5],'color','k');   line([0.5 n.netws+0.5],[b2+0.5 b2+0.5],'color','k');    end
%             set(gca,'YTick',[1:n.netws],'YTickLabel',netwname,'Fontsize',16);
%             set(gca,'XTick',[1:n.netws],'XTickLabel',netwname,'Fontsize',16);
%             
% %             print([homedir '/matrix/' 'level' num2str(i) '_netws_allsubj_rfx.tif'],'-dtiff','-r150');
%             saveas(gcf, [bandname{b} '_netw_level' num2str(i) '_RFX.png' ]);
%             
%             
%             
%             clear dataseeds_ave; clear dataseeds_se; clear dataseeds_tscore;
%             clear datanetws_ave; clear datanetws_se; clear datanetws_tscore;
%             clear crit_t; clear crit_t_unc;
%             
%         end
%     end
% end
% 
% close all
% clear a2; clear i; clear mm;


% %-------------------------------------------------------------------------
% %% STEP 3D: Statistics for Baseline Run
% % SAME AS ABOVE BUT AVERAGING ACROSS THE TWO EXPERIMENTAL GROUPS FOR THE
% % FIRST SESSION (I.E., BASELINE RS SCAN)
% % ST- Removed because it isn't relevent to this analysis

% %-------------------------------------------------------------------------
%% STEP 3E: ANOVA %% ADJUSTED BY GA
% % Run an ANOVA on z-transformed correlation values withthree terms:
% % main effects of session, group and the group x session interaction.
% 
% ANOVA = [homedir filesep 'results' filesep 'ANOVA_n' nsubj];
% if ~exist(ANOVA); mkdir(ANOVA); end; cd(ANOVA);
% 
% F_netw_anova_mat=zeros(n.netws,n.netws,size(design_matrix.SessANOVA,2)); % 3 terms; session effec
% pval_netw_anova_mat=ones(n.netws,n.netws,size(design_matrix.SessANOVA,2));
% 
% for b = 1:n.bands
%     for w=1:n.netws
%         for ww=1:n.netws
%             if strcmp(bandname{b},'alpha') == 1
%                 vect=squeeze(CC_zscore_netw_alpha(w,ww,:));
%             elseif strcmp(bandname{b},'beta') == 1
%                 vect=squeeze(CC_zscore_netw_beta(w,ww,:));
%             elseif strcmp(bandname{b},'delta') == 1
%                 vect=squeeze(CC_zscore_netw_delta(w,ww,:));
%             elseif strcmp(bandname{b},'gamma') == 1
%                 vect=squeeze(CC_zscore_netw_gamma(w,ww,:));
%             elseif strcmp(bandname{b},'theta') == 1
%                 vect=squeeze(CC_zscore_netw_theta(w,ww,:));
%             elseif strcmp(bandname{b},'bb') == 1
%                 vect=squeeze(CC_zscore_netw_bb(w,ww,:));
%             end
% %             vector of correlations between specific ROIs (as specified by q and w; length of vect = num of subjects)
% 
%             t = table(num2str(design_matrix.SessANOVA(1:n.subj)), vect(1:n.subj), vect(n.subj+1:n.subj*n.sess), 'VariableNames', {'Group', 'pre', 'post'});
% %                     t = table(num2str(vect(1:n.subj), vect(n.subj+1:n.subj*n.sess), 'VariableNames', {'pre', 'post'}));
% 
% 
%             Sess = table([1 2]', 'VariableNames', {'Sessions'});
%             rm = fitrm(t,'pre-post~1','WithinDesign', Sess);
%             clear t; clear Sess; clear vect;
% 
% %                 [anovatbl] = ranova(rm, 'WithinModel', 'Sessions');  % Get between subject (Group) effect
% %                 pval.netw_anova_mat(1,w,2) = anovatbl{'Group','pValue'};
% % 
% %                 F.netw_anova_mat(1,w,2) = anovatbl{'Group','F'};
% %                 clear anovatbl;
% 
%             [ranovatbl] = ranova(rm);  % Get within subject (session; session by group) effects
%             pval_netw_anova_mat(w,ww) = ranovatbl{'(Intercept):Sessions','pValue'};
% 
%             F_netw_anova_mat(w,ww) = ranovatbl{'(Intercept):Sessions','F'};
% 
% %                 pval.netw_anova_mat(1,w,1) = ranovatbl{'(Intercept):Sessions','pValue'};
% %             
% %                 F.netw_anova_mat(1,w,1) = ranovatbl{'(Intercept):Sessions','F'};
%             clear ranovatbl;
%             clear rm;
%         end
%     end
%     clear w; clear ww;
% 
%     for i=1:size(design_matrix.SessANOVA,2)
% 
%         count = 1;
%         for nnetw_col = 1:1:n.netws
%             pval.anovanetw(i, count) = pval_netw_anova_mat(1,nnetw_col,i);
%             F_anovanetw(i,count) = F_netw_anova_mat(1,nnetw_col,i);
%             count = count + 1;
%         end
%         clear count; clear nnetw_col;
% 
%         [~, dummy2, ~, dummy4]=fdr_bh(pval.anovanetw(i,:),fdr_thres,'pdep','no');
%         crit_p.anovanetw(i) = dummy2;
%         adj_p.anovanetw(i,1:length(dummy4)) = dummy4;
%         clear dummy2; clear dummy4;
%         index = find(pval.anovanetw(i,:) == crit_p.anovanetw(i));
%         if isempty(index)
%             crit_F = 1000; % place holder
%         else
%             crit_F = F_anovanetw(i,index(1));
%         end
%         clear index;
% 
%         mm=6;
%         [ivect,jvect]=find(abs(F_netw_anova_mat(:,:,i))>=abs(crit_F));
%         [ivect_unc,jvect_unc]=find(pval_netw_anova_mat(:,:,i)<=0.05);
%         clear crit_F;
% 
%         figure(); set(gcf,'Color','white'); box OFF; 
%         imagesc(F_netw_anova_mat(:,:,i)); axis square; caxis([0 mm]); colormap(colors_onesided); colorbar; daspect([1 1 1]);
%         title([ bandname{b} ' Timepoint Main Effect n= ' nsubj ]); % generate corresponding rfx figure
% 
%     end
% 
%     hold on; for a2=1:n.netws-1  line([a2+0.5 a2+0.5],[0.5 n.netws+0.5],'color','k');   line([0.5 n.netws+0.5],[a2+0.5 a2+0.5],'color','k');    end
%     clear a2; 
%     for b2=1:n.netws    line([b2+0.5 b2+0.5],[0.5 n.netws+0.5],'color','k');   line([0.5 n.netws+0.5],[b2+0.5 b2+0.5],'color','k');    end
%     set(gca,'YTick',[1:n.netws],'YTickLabel',netwname,'Fontsize',16);
%     set(gca,'XTick',[1:n.netws],'XTickLabel',netwname,'Fontsize',16);
% 
% 
%     scatter(jvect_unc,ivect_unc,36,'ok');
%     clear ivect_unc; clear jvect_unc;
% 
%     scatter(jvect,ivect,9,'ok','filled');
%     clear ivect; clear jvect;
% 
%  
% %     print([homedir  '/matrix/' bandname{b} 'Netw_Session_Main_Effect.tif'],'-dtiff','-r150');
%     saveas(gcf, [bandname{b} 'netw_session_main_effect.png' ]);
% 
%     clear a2; clear i; clear mm;
% end
% close all
% % 
% -------------------------------------------------------------------------
%% STEP 3F: Within-level correlations w/ behaviour %% ADJUSTED BY GA
% CORRELATIONS W/ BEHAV FACTORS (AGE, BASELINE GABA, CHANGE IN
% GABA, LEARNING MAGNITUDE, AND OFFLINE GAINS). ONE FOR EACH LEVEL (I.E.,
% SESSION AND GROUP COMBO)

withinlevel_bhv_scatter = [homedir filesep 'results' filesep 'within_level_corr_bhv_scatter_n' nsubj];
if ~exist(withinlevel_bhv_scatter); mkdir(withinlevel_bhv_scatter); end; 
% withinlevel_bhv_FDR_scatter = [homedir filesep 'results' filesep 'within_level_corr_bhv_scatter_FDR_n' nsubj];
% if ~exist(withinlevel_bhv_FDR_scatter); mkdir(withinlevel_bhv_FDR_scatter); end; 
withinlevel_bhv = [homedir filesep 'results' filesep 'within_level_corr_bhv_n' nsubj];
if ~exist(withinlevel_bhv); mkdir(withinlevel_bhv); end; 


CC_zscore_cov_ByLevel = NaN(n.netws,n.netws,n.levels,n.bands);

for bh=1:n.bhv_vars
    if bh==1
        behav_ByLevel = speed.offline;
    elseif bh==2
        behav_ByLevel = speed.online;
    elseif bh==3
        behav_ByLevel = speed.mean;
    end
    for i=1:n.levels
        for b = 1:n.bands
            if strcmp(bandname{b},'alpha') == 1
                CC_zscore_cov_ByLevel = CC_zscore_netw_alpha(1:n.netws,1:n.netws,find(design_matrix.withinlevel(:,i)));
            elseif strcmp(bandname{b},'beta') == 1
                CC_zscore_cov_ByLevel = CC_zscore_netw_beta(1:n.netws,1:n.netws,find(design_matrix.withinlevel(:,i)));
            elseif strcmp(bandname{b},'delta') == 1
                CC_zscore_cov_ByLevel = CC_zscore_netw_delta(1:n.netws,1:n.netws,find(design_matrix.withinlevel(:,i)));
            elseif strcmp(bandname{b},'gamma') == 1
                CC_zscore_cov_ByLevel = CC_zscore_netw_gamma(1:n.netws,1:n.netws,find(design_matrix.withinlevel(:,i)));
            elseif strcmp(bandname{b},'theta') == 1
                CC_zscore_cov_ByLevel = CC_zscore_netw_theta(1:n.netws,1:n.netws,find(design_matrix.withinlevel(:,i)));
            elseif strcmp(bandname{b},'bb') == 1
                CC_zscore_cov_ByLevel = CC_zscore_netw_bb(1:n.netws,1:n.netws,find(design_matrix.withinlevel(:,i)));
            end
            for w=1:n.netws
                for ww=1:n.netws
                    val=zeros(1);
                    vect=not(isnan(behav_ByLevel));
                    if sum(vect)>1
                        [val,px]=corr(squeeze(CC_zscore_cov_ByLevel(w,ww,:)),behav_ByLevel(vect)');
                    end
                    corrmat_behav.cov_WithinLevel(w,ww,:)=val;
                    pmat_behav.cov_WithinLevel(w,ww,:)=px;
                end
            end
            clear val; clear px; clear w; clear ww; clear vect;

            pval.behav.cov_WithinLevel_4cor = [ pmat_behav.cov_WithinLevel(1:n.netws)  pmat_behav.cov_WithinLevel(n.netws+2:n.netws*2)  pmat_behav.cov_WithinLevel(n.netws*2+3:n.netws*3)  pmat_behav.cov_WithinLevel(n.netws*3+4:n.netws*4)  pmat_behav.cov_WithinLevel(n.netws*4+5:n.netws*5)  pmat_behav.cov_WithinLevel(n.netws*5+6:n.netws*6)];
   
            % FDR correction
            [~, dummy2, ~, dummy4]=fdr_bh(pval.behav.cov_WithinLevel_4cor,fdr_thres,'pdep','no');
            crit_p.Behavcov(i) = dummy2;
            adj_p.Behavcov(1:length(dummy4)) = dummy4;
            clear dummy2; clear dummy4;
            
            
            cd(withinlevel_bhv_scatter);
            %%create scatter plots (added 24/8/20 ST)
            reset(gcf); reset(gca);
            
            for w=1:n.netws
                for ww=1:n.netws
                    if pmat_behav.cov_WithinLevel(w,ww) <0.051 == 1
                    figure
                    x = squeeze(CC_zscore_cov_ByLevel(w,ww,:)); y = transpose(behav_ByLevel);
                    scatter(x,y,[], 'b', 'filled');lsline;
                    title([bandname{b} ' run ' num2str(i) ' ' netwname{w} '-' netwname{ww} ' connectivity x '  bhv_vars{bh} ]);xlabel('connectivity (z-score)');ylabel([bhv_vars{bh} ' gains speed']);
                    set(gcf, 'Position',  [100, 100, 500, 200]);legend(strcat('r-value=', num2str(corrmat_behav.cov_WithinLevel(w,ww))), strcat('p-value(uncorr)=  ', num2str(pmat_behav.cov_WithinLevel(w,ww))), 'Location','northoutside');
                    labelpoints(x,y,c,'E', .05);
                    saveas(gcf, [bandname{b} '_' bhv_vars{bh} '_' netwname{w} '_' netwname{ww} '_run' num2str(i)  '_conn_scatter.png' ]);
                    clear('x','y');
                    end
                end
            end
            
        

            mm=0.8;
            [ivect,jvect]=find(pmat_behav.cov_WithinLevel(:,:)<=crit_p.Behavcov(i));
            [ivect_unc,jvect_unc]=find(pmat_behav.cov_WithinLevel(:,:)<=fdr_thres);
            figure(); set(gcf,'Color','white'); box OFF; 
            imagesc(corrmat_behav.cov_WithinLevel(:,:)); axis square; caxis([-mm mm]); 
            colormap('redblue'); colorbar; daspect([1 1 1]); 
            
            title([ bandname{b} ' Netw: Level' num2str(i) ' Corr w/ Factor ' bhv_vars{bh} ' n= ' nsubj ]); % generate corresponding rfx figure
            hold on; for a2=1:n.netws-1    line([a2+0.5 a2+0.5],[0.5 n.netws+0.5],'color','k');   line([0.5 n.netws+0.5],[a2+0.5 a2+0.5],'color','k');    end
            clear a2;
            %         set(gca,'YTick',[1],'YTickLabel',Label.Seed,'Fontsize',16);
            for b2=1:n.netws    line([b2+0.5 b2+0.5],[0.5 n.netws+0.5],'color','k');   line([0.5 n.netws+0.5],[b2+0.5 b2+0.5],'color','k');    end
            set(gca,'YTick',[1:n.netws],'YTickLabel',netwname,'Fontsize',16);
            set(gca,'XTick',[1:n.netws],'XTickLabel',netwname,'Fontsize',16);

            scatter(jvect_unc,ivect_unc,36,'ok');
            clear ivect_unc; clear jvect_unc;
            scatter(jvect,ivect,9,'ok','filled');
            clear ivect; clear jvect;

            cd(withinlevel_bhv);
            saveas(gcf, [bandname{b} '_netw_level' num2str(i) '_corr_' bhv_vars{bh} '.png' ]);
        end
    end
end
close all
clear CC_zscore_cov_ByLevel; clear behav_ByLevel; clear pval_Behav;
clear i;

% %% STEP 3G: Baseline run correlations with behaviour
% % SAME AS ABOVE BUT AVERAGING ACROSS THE TWO EXPERIMENTAL GROUPS FOR THE
% % FIRST SESSION (I.E., BASELINE RS SCAN)
% % ST- Removed because it isn't relevent to this analysis
% 
% % %-------------------------------------------------------------------------
%% STEP 3H: Behavioral Correlations with inter-session changes in connectivity %% ADJUSTED BY GA
% 
% btwsess_bhv =[homedir filesep 'results' filesep 'btw_session_corr_bhv_n' nsubj];
% if ~exist(btwsess_bhv); mkdir(btwsess_bhv);end; cd(btwsess_bhv);
% 
% btwsess_bhv_scatter = [homedir filesep 'results' filesep 'btw_session_corr_bhv_scatter_n' nsubj];
% if ~exist(btwsess_bhv_scatter); mkdir(btwsess_bhv_scatter); end; 
% 
% CC_zscore_cov_InterSession = NaN(n.netws,n.netws,n.subj,n.bands);
% 
% for bh=1:n.bhv_vars
%     if bh==1
%         behav_ByLevel = speed.offline;
%     elseif bh==2
%         behav_ByLevel = speed.online;
%     elseif bh==3
%         behav_ByLevel = speed.mean;
%     end
%     
%     for i = 1:1:n.subj
%         for w=1:n.netws
%             for ww = 1:n.netws
%                 for b = 1:n.bands
%                     if strcmp(bandname{b},'alpha') == 1
%                         CC_zscore_cov_InterSession(w,ww,i,b) = CC_zscore_netw_alpha(w,ww,i+n.subj) - CC_zscore_netw_alpha(w,ww,i);
%                     elseif strcmp(bandname{b},'beta') == 1
%                         CC_zscore_cov_InterSession(w,ww,i,b) = CC_zscore_netw_beta(w,ww,i+n.subj) - CC_zscore_netw_beta(w,ww,i);
%                     elseif strcmp(bandname{b},'delta') == 1
%                         CC_zscore_cov_InterSession(w,ww,i,b) = CC_zscore_netw_delta(w,ww,i+n.subj) - CC_zscore_netw_delta(w,ww,i);
%                     elseif strcmp(bandname{b},'gamma') == 1
%                         CC_zscore_cov_InterSession(w,ww,i,b) = CC_zscore_netw_gamma(w,ww,i+n.subj) - CC_zscore_netw_gamma(w,ww,i);
%                     elseif strcmp(bandname{b},'theta') == 1
%                         CC_zscore_cov_InterSession(w,ww,i,b) = CC_zscore_netw_theta(w,ww,i+n.subj) - CC_zscore_netw_theta(w,ww,i);
%                     elseif strcmp(bandname{b},'bb') == 1
%                         CC_zscore_cov_InterSession(w,ww,i,b) = CC_zscore_netw_bb(w,ww,i+n.subj) - CC_zscore_netw_bb(w,ww,i);
%                     end
%                 end
%             end
%         end
%     end
%     clear i; clear w; clear ww;
%     
%     for w=1:n.netws
%         for ww=1:n.netws
%             for b = 1:n.bands
%                 val=zeros(1);
%                 vect=not(isnan(behav_ByLevel));
%                 if sum(vect)>1
%                     [val,px]=corr(squeeze( CC_zscore_cov_InterSession(w,ww,:,b)),behav_ByLevel(vect)');
%                 end
%                 corrmat_behav.cov_InterSession.bandtemp{b}(w,ww,:)=val;
%                 pmat_behav.cov_InterSession.bandtemp{b}(w,ww,:)=px;
%             end
%         end
%     end
%     clear val; clear px; clear w; clear ww; clear vect;
%     
%     for b=1:n.bands
%         pval.behav.cov_InterSession_4cor.bandtemp{b} = [ pmat_behav.cov_InterSession.bandtemp{b}(1:n.netws)  pmat_behav.cov_InterSession.bandtemp{b}(n.netws+2:n.netws*2)  pmat_behav.cov_InterSession.bandtemp{b}(n.netws*2+3:n.netws*3)  pmat_behav.cov_InterSession.bandtemp{b}(n.netws*3+4:n.netws*4)  pmat_behav.cov_InterSession.bandtemp{b}(n.netws*4+5:n.netws*5)  pmat_behav.cov_InterSession.bandtemp{b}(n.netws*5+6:n.netws*6)];
%         
%         %%create scatter plots (added 24/8/20 ST)
%         reset(gcf); reset(gca);
%         cd(btwsess_bhv_scatter);
%         for w=1:n.netws
%             for ww=1:n.netws
%                 if pmat_behav.cov_InterSession.bandtemp{b}(w,ww) <0.051 == 1
%                     figure
%                     x = squeeze(CC_zscore_cov_InterSession(w,ww,:,b)); y = transpose(behav_ByLevel);
%                     scatter(x,y,[], 'b', 'filled');lsline;
%                     title([bandname{b} ' Intersession ' netwname{w} '-' netwname{ww} ' connectivity x '  bhv_vars{bh} ]);xlabel('connectivity (z-score)');ylabel([bhv_vars{bh} ' speed']);
%                     set(gcf, 'Position',  [100, 100, 500, 200]);legend(strcat('r-value=', num2str(corrmat_behav.cov_InterSession.bandtemp{b}(w,ww))), strcat('p-value(uncorr)=  ', num2str(pmat_behav.cov_InterSession.bandtemp{b}(w,ww))), 'Location','northoutside');
%                     labelpoints(x,y,c,'E', .05);
%                     saveas(gcf, [bandname{b} '_' bhv_vars{bh} '_' netwname{w} '_' netwname{ww} '_run' num2str(i)  '_conn_scatter.png' ]);
%                     clear('x','y');
%                 end
%             end
%         end
%         
%         
%         %     for ww=1:n.behav
%         %         count = 1;
%         %         for nMSL_col = 1:1:n.roi_MSL
%         %             pval_Behav(w,ww) = pmat_behav.cov_WithinLevel(w,ww);
%         %             count = count + 1;
%         %         end
%         [~, dummy2, ~, dummy4]=fdr_bh(pval.behav.cov_InterSession_4cor.bandtemp{b},fdr_thres,'pdep','no');
%         crit_p.Behavcov = dummy2;
%         adj_p.Behavcov(1:length(dummy4)) = dummy4;
%         clear dummy2; clear dummy4;
%         
%         mm=0.8;
%         [ivect,jvect]=find(pmat_behav.cov_InterSession.bandtemp{b}(:,:)<=crit_p.Behavcov);
%         [ivect_unc,jvect_unc]=find(pmat_behav.cov_InterSession.bandtemp{b}(:,:)<=fdr_thres);
%         figure(); set(gcf,'Color','white'); box OFF;
%         imagesc(corrmat_behav.cov_InterSession.bandtemp{b}(:,:)); axis square; caxis([-mm mm]); colormap('redblue'); colorbar; daspect([1 1 1]);
%         title([ bandname{b} ' Netw: InterSession  Corr w/ ' bhv_vars{bh} ' n= ' nsubj ]); % generate corresponding rfx figure
%         
%         hold on; for a2=1:n.netws-1    line([a2+0.5 a2+0.5],[0.5 n.netws+0.5],'color','k');   line([0.5 n.netws+0.5],[a2+0.5 a2+0.5],'color','k');    end
%         clear a2;
%         %         set(gca,'YTick',[1],'YTickLabel',Label.Seed,'Fontsize',16);
%         for b2=1:n.netws    line([b2+0.5 b2+0.5],[0.5 n.netws+0.5],'color','k');   line([0.5 n.netws+0.5],[b2+0.5 b2+0.5],'color','k');    end
%         set(gca,'YTick',[1:n.netws],'YTickLabel',netwname,'Fontsize',16);
%         set(gca,'XTick',[1:n.netws],'XTickLabel',netwname,'Fontsize',16);
%         
%         scatter(jvect_unc,ivect_unc,36,'ok');
%         clear ivect_unc; clear jvect_unc;
%         
%         scatter(jvect,ivect,9,'ok','filled');
%         clear ivect; clear jvect;
%         %         xticklabel_rotate([1:n.roi_MSL],90,Label.MSL,'Fontsize',16);
%         %         print([homedir '/matrix/' 'Level' num2str(i) '_MSL_Corr_Factor ' num2str(ww) '.tif'],'-dtiff','-r150');
%         cd(btwsess_bhv);
%         saveas(gcf, [bandname{b} '_netw_btw_session_corr_' bhv_vars{bh} '.png' ]);
%         
%           clear pval_Behav; %clear CC_zscore_cov_InterSession;clear behav_ByLevel;
%         clear i;
%     end
% end
% close all


%% STEP 3I: Within-group correlations with behaviour
% CORRELATIONS W/ BEHAV FACTORS (AGE, BASELINE GABA, CHANGE IN
% GABA, LEARNING MAGNITUDE, AND OFFLINE GAINS).
% ST- Removed because it isn't relevent to this analysis

%-------------------------------------------------------------------------
%% STEP 3J: RELATING BASELINE CONNECTIVITY TO INTER-SESSION CHANGES IN CONNECTIVITY
% ST- Removed because it isn't relevent to this analysis

%-------------------------------------------------------------------------
%% SECTION 4: EXTRACTING DATA BASED ON VISUAL INSPECTION OF RESULTS ABOVE TO PRODUCE ADDITIONAL PLOTS %% GA: NOT ADJUSTED
% ST- Removed because it isn't relevent to this analysis
