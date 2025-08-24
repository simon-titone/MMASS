function sleep_within_lvl_corr(v)
%%%%% %NETW: network, network, subject, stage, session, cycle, band
%%%%% BANDS = 1-delta 2-theta 3-alpha 4-sigma 5-beta 6-gamma
%% Define groups
%individual sessions
v.EXP_offline = unique([4 12 16 17 3 6 7 20 14 ]);                            % n = 20
v.EXP_offline_REM = unique([1 3 4 7 9 10 11 12 13 16 17 18 21 24 27 28 29]);

v.CTRL_NREM2 = unique([9 14 17 3 6 7 20 14   21]);
v.CTRL_NREM3 = unique([9 14 17 28 3 6 7 20 14  15 22]);
v.CTRL_REM = unique([3 4 6 9 10 11 12 14 17 28 3 6 7 20 14    21]);
%comparing correlations between EXP x offline & CTRL x offline

v.ztest_NREM2 = unique([3 4 6 7 9 12 14 16 17 20]);
v.ztest_NREM3 = unique([3 4 6 7 9 12 14 16 17 20 28]);
v.ztest_REM = unique([1 3 4 6 7 9 10 11 12 13 14 16 17 18 20 21 24 27 28 29]);

%% Select group

% v.group = v.EXP_offline; stage = 1:2;  sess = 2;

% v.group = v.EXP_offline_REM; stage = 3; sess_name = 'EXP';sess = 2;

v.group = v.CTRL_NREM2; stage = 1; sess_name = 'CTRL';sess = 1;
% v.group = v.CTRL_NREM3; stage = 2; sess_name = 'CTRL';sess = 1;
% v.group = v.CTRL_REM; stage = 3; sess_name = 'CTRL';sess = 1;

% for stage_ = 1:v.nstages
%     if stage_ == 1
%         v.group = v.ztest_NREM2; stage = 1 ;sess = 1:2;
%     elseif stage_ == 2
%         v.group = v.ztest_NREM3; stage = 2;sess = 1:2;
%     else
%         v.group = v.ztest_REM; stage = 3 ;sess = 1:2;
%     end
clear speed; clear v.subj; clear subj_label;
%% Load subjects
v.subj = 1:29; v.group = sort(unique(v.group), 'descend');
for i = 1:length(v.group); v.subj(v.subj(v.group(i))) = []; end
v.nsubj = length(v.subj); % number of independent subjects
study3_load(v);

cd(v.vardir);load('v.mat');
load('bhv_vars.mat', 'off_ITI');
speed(1,:) = off_ITI(v.subj);
for c2 = 1:v.nsubj; c1(c2) = v.subj(c2); end
subj_label= c1'; v.violin_subj_label = [c1' ; c1'];
bh = 1;
%% Correlation between behavior and sleep during control or experimental nights

w_lvl_corr = [v.homedir filesep 'results/Correlation/Within_level_x_BHV' ]; % name dir
if ~exist(w_lvl_corr); mkdir(w_lvl_corr);end; cd(w_lvl_corr); % create dir

for se = sess
    for st = stage
        w_lvl_corr_sess = [w_lvl_corr filesep v.sess{se}  filesep v.stages{st}];
        if ~exist(w_lvl_corr_sess); mkdir(w_lvl_corr_sess); end
        withinlevel_bhv_scatter = [ w_lvl_corr_sess filesep 'scatter_n' num2str(v.nsubj)];
        if ~exist(withinlevel_bhv_scatter); mkdir(withinlevel_bhv_scatter); end

        for b=1:v.nbands
            conn_temp = NETW(:,:,:,st,se,b);
            behavior =  speed';
            for w=1:v.nnetw
                for ww=1:v.nnetw

                    conn_vect = squeeze(NETW(w,ww,:,st,2,1,b));
                    n_nan = find(isnan(conn_vect));
                    sort(n_nan, 'descend');
                    for nn = 1:length(n_nan)
                        behavior(n_nan(nn)) = [];
                        conn_vect(n_nan(nn)) = [];
                    end

                    conn = squeeze(conn_temp(w,ww,:));
                    [val,px]=corr(conn,behavior);
                    r_val(w,ww,b,bh,st,se)=val;
                    p_val(w,ww,b,bh,st,se)=px;
                end
            end
            clear val; clear px; clear w; clear ww; clear vect;

            pval_fdr_temp = p_val(:,:,b,bh,st,se);
            if v.FDR == 'mot'
                p_val_temp_FDR_correction  = pval_fdr_temp(5,:);
            elseif v.FDR == 'all'
                p_val_temp_FDR_correction = [ pval_fdr_temp(1:v.nnetw)  pval_fdr_temp(v.nnetw+2:v.nnetw*2)  pval_fdr_temp(v.nnetw*2+3:v.nnetw*3)  pval_fdr_temp(v.nnetw*3+4:v.nnetw*4)  pval_fdr_temp(v.nnetw*4+5:v.nnetw*5)  pval_fdr_temp(v.nnetw*5+6:v.nnetw*6)];
            end

            [~, dummy2, ~, dummy4]=fdr_bh(p_val_temp_FDR_correction,v.fdr_thres,'pdep','no');
            crit_p = dummy2; adj_p(1:length(dummy4)) = dummy4;
            a = dummy4;
            clear dummy2; clear dummy4;
            adj_p_temp = nan(6,6);
            if v.FDR == 'mot'
                adj_p_temp(1,5) = a(1); adj_p_temp(2,5) = a(2); adj_p_temp(3,5) = a(3); adj_p_temp(4,5) = a(4);adj_p_temp(5,5) = a(5); adj_p_temp(5,6) = a(6);adj_p_mat(:,:,st,b) = adj_p_temp;
            elseif  v.FDR == 'all'
                adj_p_temp(1,1:6) = a(1,1:6); adj_p_temp(2,2:6) = a(1,7:11);adj_p_temp(3,3:6) = a(1,12:15); adj_p_temp(4,4:6) = a(1,16:18);adj_p_temp(5,5:6) = a(1,19:20); adj_p_temp(6,6) = a(1,21);
                adj_p_mat(:,:,st,b) = adj_p_temp;
            end

            [ivect,jvect]=find(p_val(:,:,b,bh,st,se)<=crit_p);
            [ivect_unc,jvect_unc]=find(p_val(:,:,b,bh,st,se)<=v.fdr_thres);

            matrix = r_val(:,:,b,bh,st,se); tri_matrix = triu(matrix);

            imagesc(tri_matrix); axis square;   daspect([1 1 1]);
            set(gcf,'Color','white');set(gca,'YDir','reverse');hold on;
            dim = set(gcf, 'Position',  [1000, 900, 500, 450]); %left bottom width height
            h = colorbar; set(get(h,'title'),'string','r-val', 'fontsize', 12); %ylabel(h, 'f-val');
            colormap('redblue');    limit = max(max(abs(tri_matrix))); caxis([-limit limit]);

            T_n = {[v.sess{se} ' ' v.bandname{b} ' '  v.stages{st} ' x ' v.bhv_vars{bh} ' {\it n} = ' num2str(v.nsubj) ]};
            T = {[ v.sess{se} ' ' v.bandname{b} ' ' v.stages{st} ' x ' v.bhv_vars{bh} ]};
            if v.fig_title == 1
                if v.fig_labels == 1 ; title(T_n);  else; title(T); end
            end


            for a2=1:v.nnetws-1; line([a2+0.5 a2+0.5],[a2+1.5 0.5],'color','k','LineWidth',0.8);line([a2-0.5 v.nnetws+0.5],[a2+0.5 a2+0.5],'color','k','LineWidth',0.8); end
            for b2=1:v.nnetws; line([b2-0.5 b2-0.5],[0.5 v.nnetws+0.5],'color','k','LineWidth',0.001); line([0.5 v.nnetws+0.5],[b2-0.5 b2-0.5],'color','k','LineWidth',0.001);end
            clear a2; clear b2;

            set(gca,'YTick',[1:v.nnetws],'YTickLabel', v.netwname,'Fontsize',16); set(gca,'XTick',[1:v.nnetws],'XTickLabel', v.netwname,'Fontsize',16);

            scatter(jvect_unc,ivect_unc, 80, 'w', 'o' ,'ok','LineWidth',2.5); clear ivect_unc; clear jvect_unc;
            scatter(jvect,ivect,80,'w','ok','filled' ,'LineWidth',2.5); clear ivect; clear jvect;

            saveas(gcf, [w_lvl_corr_sess filesep  v.sess{se} '_' v.stages{st} '_' v.bandname{b} '_' v.bhv_vars{bh} '.png'])

            %%create scatter plots (added 24/8/20 ST)
            reset(gcf); reset(gca);
            if v.scatter == 1
                for w=1:v.nnetws
                    for ww=1:v.nnetws
                        if p_val(w,ww,b,bh,st,se) <0.05 == 1
                            %                             if b == 3
                            scatter_band = [withinlevel_bhv_scatter filesep v.bandname{b}];
                            if ~exist(scatter_band); mkdir(scatter_band); end
                            figure
                            x = squeeze(conn_temp(w,ww,:)); y = behavior;
                            scatter(x,y,[], 'b', 'filled');lsline;
                            mdl = fitlm(x, y, 'linear');
                            [y_pred,y_ci] = predict(mdl,x, 'Alpha', 0.05);
                            Y = [x, y_ci,y_pred];
                            Y = sortrows(Y,1);
                            eb = [Y(:,4)-Y(:,2)];
                            boundedline(Y(:,1),Y(:,4),eb, '-r','alpha');
                            clear mdl; clear y_pred; clear y_ci; clear Y; clear eb;
                            T_scatter_n = {[ v.sess{se} ' ' v.stages{st} ' x '  v.bhv_vars{bh} ', {\it n} =  ' num2str(v.nsubj)  ] [v.bandname{b} ' ' v.netwname{w} '-' v.netwname{ww} ]};
                            T_scatter = {[ v.sess{se} ' ' v.stages{st} ' x '  v.bhv_vars{bh} ] [v.bandname{b} ' ' v.netwname{w} '-' v.netwname{ww} ]};
                            xlabel('connectivity (z-score)');ylabel([v.bhv_vars{bh} ' gains speed']);
                            xx = (std(x)/2);
                            x_min = (min(x)-xx); x_max = (max(x)+xx);
                            xlim([x_min x_max]);
                            if v.fig_title == 1
                            if v.fig_labels == 1 ; title(T_scatter_n);  else; title(T_scatter); end
                            end

                            set(gcf, 'Position',  [100, 100, 300, 320]); %left bottom width height

                            if v.pval_labels == 1
                                text(0.6, 0.98,['r = ' num2str(round(r_val(w,ww,b,bh,st,se), 3))],'Units','normalized' ,'FontSize',14);
                                if p_val(w,ww,b,bh,st,se) <0.001 == 1
                                    text(0.6, 0.93, 'p(UC) < 0.001 ' ,'Units','normalized' ,'FontSize',14);
                                else
                                    text(0.6, 0.93,['p(UC) =  ' num2str(round(p_val(w,ww,b,bh,st,se),3))],'Units','normalized' ,'FontSize',14);
                                end
                            end
                            if v.subj_labels == 1
                                labelpoints(x,y,subj_label,'E', .05);
                            end
                            S_scatter = [scatter_band filesep '_' v.sess{se} '_' v.bandname{b} '_' v.stages{st} '_'  v.bhv_vars{bh} '_' v.netwname{w} '_' v.netwname{ww}   '_scatter.png' ];
                            saveas(gcf, S_scatter);
                            clear('x','y');
                        end
                    end
                end
                close all;
            end
        end

    end
    close all;
    clear i;
end
%% Write excel file with results
if v.excel == 1
    % separate variable for significant p-values
    for st= stage
        for b = 1:v.nbands
            for w=1:v.nnetws
                for ww=1:v.nnetws
                    if p_val(w,ww,b,bh,st,se) < 0.05 && p_val(w,ww,b,bh,st,se) > 0
                        sig_r(w,ww,st,b) = r_val(w,ww,b,bh,st,se);
                        sig_p(w,ww,st,b) = p_val(w,ww,b,bh,st,se);
                    else
                        sig_r(w,ww,st,b) = NaN;
                        sig_p(w,ww,st,b) = NaN;
                    end
                end
            end
            sig_r(:,:,st,b) = triu(sig_r(:,:,st,b));
            sig_p(:,:,st,b) = triu(sig_p(:,:,st,b));
        end
    end

    for st= stage
        for b = 1:v.nbands
            for w=1:v.nnetws
                for ww=1:v.nnetws
                    if  sig_r(w,ww,st,b) ~=0 && ~isnan(sig_r(w,ww,st,b)) && sig_p(w,ww,st,b)>0.001
                        result(w,ww,st,b) = {['r(' num2str(v.nsubj-2)  ') = ' num2str(round(sig_r(w,ww,st,b),3)) ', p = ' num2str(round(sig_p(w,ww,st,b),3))]};
                    elseif  sig_r(w,ww,st,b) ~=0 && ~isnan(sig_r(w,ww,st,b)) && sig_p(w,ww,st,b)<0.001
                        result(w,ww,st,b) = {['r(' num2str(v.nsubj-2)  ') = ' num2str(round(sig_r(w,ww,st,b),3)) ', p < 0.001']};
                    else
                        result(w,ww,st,b) = {'nan'};
                    end
                    if  sig_r(w,ww,st,b) ~=0 && ~isnan(sig_r(w,ww,st,b)) && adj_p_mat(w,ww,st,b)<0.05 && adj_p_mat(w,ww,st,b)> 0
                        result_corr(w,ww,st,b) = {['r(' num2str(v.nsubj-2)  ') = ' num2str(round(sig_r(w,ww,st,b),3)) ', p = ' num2str(adj_p_mat(w,ww,st,b),3)]};
                    else
                        result_corr(w,ww,st,b) = {'nan'};
                    end
                end
            end
        end
    end

    count_title = 1;
    count_label = 4;
    v.t_group = 1;
    for st= stage
        for b = 1:v.nbands
            r= result(:,:,st,b);
            filename = [w_lvl_corr_sess filesep 'study3_stage_x_bhv_uncorrected.xlsx'];
            t = {[v.stages{st} ' ' v.bandname{b} ', n= ' num2str(v.nsubj) ]};
            t_cell = ['A' num2str(count_title) ];
            writecell(t, filename, 'Range', t_cell, 'Sheet', v.t_group )
            label_cell = ['A' num2str(count_label) ];
            writecell(v.netwname', filename, 'Range', label_cell , 'Sheet', v.t_group)
            label_b_cell = ['B' num2str(count_label-1) ];
            writecell(v.netwname, filename, 'Range', label_b_cell , 'Sheet', v.t_group)
            results_cell = ['B' num2str(count_label)];
            writecell(r, filename, 'Range', results_cell, 'Sheet', v.t_group)
            clear t; clear t_cell; clear filename; clear r; clear label_cell; clear label_b_cell; clear results_cell; clear t_cell

            r= result_corr(:,:,st,b);
            filename = [w_lvl_corr_sess filesep 'study3_stage_x_bhv_FDR.xlsx'];
            t = {[v.stages{st} ' ' v.bandname{b} ', n= ' num2str(v.nsubj) ]};
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
% save r_val.mat r_val -mat
% end
end

