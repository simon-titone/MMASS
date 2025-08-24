function study3_interstage_intersess_diff_corr_bhv(v)
%%%%% NETW: network, network, subject, stage, session,cycle, band

Intersess_NREM2_NREM3 = [3 4 6 7 9 12 14 16 17 20 28];
Intersess_NREM2_REM = [1 3 4 6 7 9 10 11 12 13 14 16 17 18 20 21 24 27 28 29];
Intersess_NREM3_REM = [1 3 4 6 7 9 10 11 12 13 14 16 17 18 20 21 24 27 28 29];

z_test_NREM2_NREM3 = [3 4 6 7 9 12 14 16 17 20 28];

v.group = Intersess_NREM2_NREM3; stagecorr = 1;
% v.group = Intersess_NREM2_REM; stagecorr = 2;
% v.group = Intersess_NREM3_REM; stagecorr = 3;

v.subj = 1:29; v.group = sort(v.group, 'descend');
for i = 1:length(v.group); v.subj(v.subj(v.group(i))) = []; end
v.nsubj = length(v.subj); % number of independent subjects

study3_load(v);

cd(v.vardir);load('v.mat');
load('bhv_vars.mat', 'off_ITI');
v.speed(1,:) = off_ITI(v.subj);

v.nsubj = length(v.subj); % number of independent subjects
for cc = 1:v.nsubj; c(cc) = v.subj(cc); end
v.subj_label= c'; v.violin_subj_label = [c' ; c'];clear c;
bh= 1;
behavior = v.speed(bh,:)';
c = 1;
%% Difference between stages, minus difference between session (e.g. (EXP NREM2 - EXP NREM3) - (CTRL NREM2 - CTRL NREM3))
for b = 1:v.nbands
    for i = 1:1:v.nsubj
        for w=1:v.nnetws
            for ww = 1:v.nnetws
                interstage_intersess_conn(w,ww,i,1,b) = ((NETW(w,ww,i,1,2,c,b) - NETW(w,ww,i,2,2,c,b)) - (NETW(w,ww,i,1,1,c,b) - NETW(w,ww,i,2,1,c,b))); %NREM2 - NREM3
                interstage_intersess_conn(w,ww,i,2,b) = ((NETW(w,ww,i,1,2,c,b) - NETW(w,ww,i,3,2,c,b)) - (NETW(w,ww,i,1,1,c,b) - NETW(w,ww,i,3,1,c,b))); %NREM2 - REM
                interstage_intersess_conn(w,ww,i,3,b) = ((NETW(w,ww,i,2,2,c,b) - NETW(w,ww,i,3,2,c,b)) - (NETW(w,ww,i,2,1,c,b) - NETW(w,ww,i,3,1,c,b))); %NREM3 - REM
            end
        end
    end
end

%% Interstage correlation with behavior
btw_session_corr_dir = [v.homedir filesep 'results' filesep 'Correlation' filesep 'Interstage_intersess_BHV' ];
for st= stagecorr
    inter_stage_corr_stage_dir =[btw_session_corr_dir filesep  v.stagediff{st}];
    if ~exist(inter_stage_corr_stage_dir); mkdir(inter_stage_corr_stage_dir); end; cd(inter_stage_corr_stage_dir)
    for b = 1:v.nbands
        for w=1:v.nnetws
            for ww=1:v.nnetws
                x = squeeze(interstage_intersess_conn(w,ww,:,st,b));
                [val,px]=corr(x,behavior);
                r_val(w,ww,st,b)=val;
                p_val(w,ww,st,b)=px;

            end
        end
    end
    clear val; clear px; clear w; clear ww; clear vect;

    for b=1:v.nbands
        pval_fdr_temp = p_val(:,:,st,b);

        index = 1;
        for w = 1:v.nnetw
            for ww = 1:v.nnetw
                pval_fdr(index) = pval_fdr_temp(w,ww);
                index = index + 1;
            end
        end

        if v.FDR == 'mot'
            p_val_temp_FDR_correction  = pval_fdr_temp(5,:);
        elseif v.FDR == 'all'
            p_val_temp_FDR_correction = pval_fdr;
        end

        [~, dummy2, ~, dummy4]=fdr_bh(p_val_temp_FDR_correction, v.fdr_thres,'pdep','no');
        crit_p = dummy2;
        a(1:length(dummy4)) = dummy4;
        clear dummy2; clear dummy4;
        adj_p_temp = nan(6,6);
        if v.FDR == 'mot'
            adj_p_temp(1,5) = a(1); adj_p_temp(2,5) = a(2); adj_p_temp(3,5) = a(3); adj_p_temp(4,5) = a(4);adj_p_temp(5,5) = a(5); adj_p_temp(5,6) = a(6);adj_p_mat(:,:,st,b) = adj_p_temp;
        elseif  v.FDR == 'all'
            adj_p_temp(1,1:6) = a(1,1:6); adj_p_temp(2,2:6) = a(1,7:11);adj_p_temp(3,3:6) = a(1,12:15); adj_p_temp(4,4:6) = a(1,16:18);adj_p_temp(5,5:6) = a(1,19:20); adj_p_temp(6,6) = a(1,21);
            adj_p_mat(:,:,st,b) = adj_p_temp;
        end

        [ivect,jvect]=find(p_val(:,:,st,b)<=crit_p);
        [ivect_unc,jvect_unc]=find(p_val(:,:,st,b)<=v.fdr_thres);

        matrix = r_val(:,:,st,b); tri_matrix = triu(matrix);

        imagesc(tri_matrix); axis square; % colormap(v.colors); daspect([1 1 1]); hold on
        h = colorbar; set(get(h,'title'),'string','r-val', 'fontsize', 12); %ylabel(h, 'f-val');
        dim = set(gcf, 'Position',  [1000, 900, 500, 450]); %left bottom width height
        colormap('redblue'); limit = max(max(abs(tri_matrix))); caxis([-limit limit]);
        if v.fig_labels == 1
            T = { ['EXP - CTRL ' v.bandname{b}  ', {\it n} = ' num2str(v.nsubj)] [ v.stagediff_title{st}  ' x '  v.bhv_vars{bh} ]};
        else
            T = { ['EXP - CTRL ' v.bandname{b} ] [ v.stagediff_title{st}  ' x '  v.bhv_vars{bh} ]};
        end
        title(T); % generate corresponding rfx figure

        for a2=1:v.nnetws-1; line([a2+0.5 a2+0.5],[a2+1.5 0.5],'color','k','LineWidth',0.8);line([a2-0.5 v.nnetws+0.5],[a2+0.5 a2+0.5],'color','k','LineWidth',0.8); end
        for b2=1:v.nnetws; line([b2-0.5 b2-0.5],[0.5 v.nnetws+0.5],'color','k','LineWidth',0.001); line([0.5 v.nnetws+0.5],[b2-0.5 b2-0.5],'color','k','LineWidth',0.001);end
        clear a2; clear b2;

        set(gca,'YTick',[1:v.nnetws],'YTickLabel',v.netwname,'Fontsize',16);
        set(gca,'XTick',[1:v.nnetws],'XTickLabel',v.netwname,'Fontsize',16);
        hold on;
        scatter(jvect_unc,ivect_unc, 80, 'w', 'o' ,'ok','LineWidth',2.5); clear ivect_unc; clear jvect_unc;
        scatter(jvect,ivect,80,'w','ok','filled' ,'LineWidth',2.5); clear ivect; clear jvect;

        S = [v.bandname{b} '_intersess' v.stagediff{st} '_' v.bhv_vars{bh} '.png' ];
        saveas(gcf, S);

        clear CC_zscore_cov_InterSession; clear pval_Behav;
        clear i;

        if v.scatter == 1
            for w=1:v.nnetws
                for ww=1:v.nnetws
                    if p_val(w,ww,st,b) <0.05 == 1
                        scatter_band = [inter_stage_corr_stage_dir  '/scatter/'  v.bandname{b}];
                        if ~exist(scatter_band); mkdir(scatter_band); end

                        figure
                        x = squeeze(interstage_intersess_conn(w,ww,:,st,b)); y = behavior;
                        scatter(x,y,[], 'b', 'filled');lsline;
                        mdl = fitlm(x, y, 'linear');
                        [y_pred,y_ci] = predict(mdl,x, 'Alpha', 0.05);
                        Y = [x, y_ci,y_pred];
                        Y = sortrows(Y,1);
                        eb = [Y(:,4)-Y(:,2)];
                        boundedline(Y(:,1),Y(:,4),eb, '-r','alpha');
                        clear mdl; clear y_pred; clear y_ci; clear Y; clear eb;
                        if v.fig_labels == 1
                            T_scatter = {['EXP - CTRL ' v.stagediff_title{st} ' x ' v.bhv_vars{bh} ] [v.netwname{w} '-' v.netwname{ww} ' ' v.bandname{b}  ', {\it n} = ' num2str(v.nsubj)]};
                        else
                            T_scatter = {['EXP - CTRL ' v.stagediff_title{st} ' x ' v.bhv_vars{bh} ] [v.netwname{w} '-' v.netwname{ww} ' '  v.bandname{b}  ]};
                        end
                        title(T_scatter); xlabel('connectivity (z-score)'); ylabel([v.bhv_vars{bh} ' gains speed']);
                        set(gcf, 'Position',  [100, 100, 400, 320]); %left bottom width height

                        if v.subj_labels == 1
                            labelpoints(x,y,v.subj_label,'E', .05);
                        end
                        if v.pval_labels == 1
                            text(0.6, 0.98,['r = ' num2str(round(r_val(w,ww,st,b), 3))],'Units','normalized' ,'FontSize',14);
                            if p_val(w,ww,st,b) <0.001 == 1
                                text(0.6, 0.93, 'p(UC) < 0.001 ' ,'Units','normalized' ,'FontSize',14);
                            else
                                text(0.6, 0.93,['p(UC) =  ' num2str(round(p_val(w,ww,st,b),3))],'Units','normalized' ,'FontSize',14);
                            end
                        end

                        S_scatter = [ scatter_band filesep v.stagediff{st} '_x_BHV'  v.netwname{w} '_' v.netwname{ww} '.png' ];
                        saveas(gcf, S_scatter);
                        clear('x','y');
                    end
                end
            end
        end
    end
end

close all


%% Z-score comparison between two sessions to determine differences in correlations
if v.z_test == 1
    for s = 1 %stage
        for b = 1:v.nbands % freq
            for w = 1:6 %netw
                for ww = 1:6
                    n = v.nsubj;
                    ctrl = r_val(w,ww,1,st,b);
                    exp = r_val(w,ww,2,st,b);
                    [p,z] = corr_rtest(ctrl,exp,n,n);
                    p_1tailed(w,ww,b,s) = p(1,1);
                    p_2tailed(w,ww,b,s) = p(1,2);
                    z_all(w,ww,b,s) = z;
                end
            end
        end
    end

    for s = 1 %stage
        for b = 1:6 % freq
            z_pval(:,:,b,s) = triu(p_2tailed(:,:,b,s));
            z_score(:,:,b,s) = triu(z_all(:,:,b,s));
        end
    end
    z_test_savedir = [v.homedir filesep 'results' filesep 'Correlation' filesep 'Interstage_x_BHV'];
    cd(z_test_savedir);
    save("zscore.mat", "z_pval", "z_score");

    for s= 1
        for b = 1:v.nbands
            for w=1:v.nnetws
                for ww=1:v.nnetws
                    if  z_score(w,ww,b,s) ~=0
                        z_test_result(w,ww,b,s) = {['z = ' num2str(round(z_score(w,ww,b,s),3)) ', p = ' num2str(round(z_pval(w,ww,b,s),3))]};
                    end
                end
            end
        end
    end
    % Writes excel file with results matrix

    count_title = 1;
    count_label = 4;
    v.t_group = 1;
    for b = 1:v.nbands
        r =  triu(z_test_result(:,:,b,s));
        filename = [z_test_savedir filesep 'study3_interstage_x_bhv_z-test.xlsx'];
        t = {[ ' ' v.bandname{b} ', n= ' num2str(v.nsubj) ]};
        t_cell = ['A' num2str(count_title) ];
        writecell(t, filename, 'Range', t_cell, 'Sheet', v.t_group )
        label_cell = ['A' num2str(count_label) ];
        writecell(v.netwname', filename, 'Range', label_cell , 'Sheet', v.t_group)
        label_b_cell = ['B' num2str(count_label-1) ];
        writecell(v.netwname, filename, 'Range', label_b_cell , 'Sheet', v.t_group)
        results_cell = ['B' num2str(count_label)];
        writematrix(r, filename, 'Range', results_cell, 'Sheet', v.t_group)
        clear t; clear t_cell; clear filename; clear r; clear label_cell; clear label_b_cell; clear results_cell; clear t_cell

        count_label = count_label + 10;
        count_title = count_title + 10;
    end
end

%% Write excel file with results
if v.excel == 1
    % separate variable for significant p-values
    for st= 1
        for b = 1:v.nbands
            for w=1:v.nnetws
                for ww=1:v.nnetws
                    if p_val(w,ww,st,b) < 0.05 && p_val(w,ww,st,b) > 0
                        sig_r(w,ww,st,b) = r_val(w,ww,st,b);
                        sig_p(w,ww,st,b) = p_val(w,ww,st,b);
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

    for st= 1
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
    for st= 1
        for b = 1:v.nbands
            r= result(:,:,st,b);
            filename = [btw_session_corr_dir filesep 'study3_interstage_x_bhv_uncorrected.xlsx'];
            t = {[v.stagediff_title{st} ' ' v.bandname{b} ', n= ' num2str(v.nsubj) ]};
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
            filename = [btw_session_corr_dir filesep 'study3_interstage_x_bhv_FDR.xlsx'];
            t = {[v.stagediff_title{st} ' ' v.bandname{b} ', n= ' num2str(v.nsubj) ]};
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

