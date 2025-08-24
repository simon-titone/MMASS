function sleep_btw_sess_bhv_corr(v)
%%%%% %NETW: network, network, subject, stage, session, cycle, band
%%%%% BANDS = 1-delta 2-theta 3-alpha 4-sigma 5-beta 6-gamma
%% Load Subjects
v.EXP_offline = [4 12 16 17 3 6 7 20 14 ];                            % n = 20
v.EXP_offline_REM = [1 3 4 7 9 10 11 12 13 16 17 18 21 24 27 28 29];

v.CTRL_NREM2 = [9 14 17 3 6 7 20 14];
v.CTRL_NREM3 = [9 14 17 28 3 6 7 20 14];
v.CTRL_REM = [3 4 6 9 10 11 12 14 17 28 3 6 7 20 14];

sessBYstage_corr = unique([v.EXP_offline v.EXP_offline_REM v.CTRL_NREM2 v.CTRL_NREM3 v.CTRL_REM ]);

v.group = v.EXP_offline; stage = 1:2; sess_name = 'EXP'; sess = 2;
% v.group = v.EXP_offline_REM; stage = 3; sess = 'EXP';sess = 2;

% v.group = v.CTRL_NREM2; stage = 1; sess = 'CTRL';sess = 1;
% v.group = v.CTRL_NREM3; stage = 2; sess = 'CTRL';sess = 1;
% v.group = v.CTRL_REM; stage = 3; sess = 'CTRL';sess = 1;

v.subj = 1:29; v.group = sort(unique(v.group), 'descend');
for i = 1:length(v.group); v.subj(v.subj(v.group(i))) = []; end
v.nsubj = length(v.subj); % number of independent subjects
Sleep_load(v);

cd(v.vardir);load('v.mat');
load('bhv_vars.mat', 'off_ITI');
v.speed(1,:) = off_ITI(v.subj);
for c2 = 1:v.nsubj; c1(c2) = v.subj(c2); end
v.subj_label= c1'; v.violin_subj_label = [c1' ; c1'];
%% Correlation between behavior and sleep during control or experimental nights

w_lvl_corr = [v.homedir filesep 'results' filesep sess_name '_corr_n' num2str(v.nsubj)];
if ~exist(w_lvl_corr); mkdir(w_lvl_corr);end; cd(w_lvl_corr);

% for se = 1:v.nsess
%     for c= 1:v.ncycles
%         for st= 1:v.nstages
%             for b = 1:v.nbands
%                 for i = 1:v.nsubj
%                     for w=1:v.nnetws
%                         for ww=1:v.nnetws
%                             if NETW(w,ww,i,st,1,c,b) ~= 0 && NETW(w,ww,i,st,2,c,b) ~= 0
%                                 NETW_temp(w,ww,se,c,st,b,i) = NETW(w,ww,i,st,se,c,b);
%                             else
%                                 NETW_temp(w,ww,se,c,st,b,i) = NaN;
%                             end
%                         end
%                     end
%                 end
%             end
%         end
%     end
% end



for se = sess
    w_lvl_corr_sess = [w_lvl_corr filesep v.sess{se}];
    if ~exist(w_lvl_corr_sess); mkdir(w_lvl_corr_sess);end;
    withinlevel_bhv_scatter = [ w_lvl_corr_sess filesep 'scatter_n' num2str(v.nsubj)];
    if ~exist(withinlevel_bhv_scatter); mkdir(withinlevel_bhv_scatter); end;
    for st = stage
        for bh= 1
            for b = 1:v.nbands
                conn_temp = NETW(:,:,:,st,se,b);
                behavior = v.speed';
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

                if v.FDR == 'mot'
                    p_val_temp_FDR_correction  = p_val(5,:);
                elseif v.FDR == 'all'
                    p_val_temp_FDR_correction = [ p_val(1:v.nnetw)  p_val(v.nnetw+2:v.nnetw*2)  p_val(v.nnetw*2+3:v.nnetw*3)  p_val(v.nnetw*3+4:v.nnetw*4)  p_val(v.nnetw*4+5:v.nnetw*5)  p_val(v.nnetw*5+6:v.nnetw*6)];
                end

                [~, dummy2, ~, dummy4]=fdr_bh(p_val_temp_FDR_correction,v.fdr_thres,'pdep','no');
                crit_p = dummy2; adj_p(1:length(dummy4)) = dummy4; clear dummy2; clear dummy4;

                [ivect,jvect]=find(p_val(:,:,b,bh,st,se)<=crit_p);
                [ivect_unc,jvect_unc]=find(p_val(:,:,b,bh,st,se)<=v.fdr_thres);

                matrix = r_val(:,:,b,bh,st,se); tri_matrix = triu(matrix);

                imagesc(tri_matrix); axis square;   daspect([1 1 1]);
                set(gcf,'Color','white');set(gca,'YDir','reverse');hold on;
                dim = set(gcf, 'Position',  [1000, 900, 500, 450]); %left bottom width height
                h = colorbar; set(get(h,'title'),'string','r-val', 'fontsize', 12); %ylabel(h, 'f-val');
                colormap('redblue');    limit = max(max(abs(tri_matrix))); caxis([-limit limit]);
                title({[v.bandname{b} ' ' v.sess{se} ' ' v.stages{st} ' x ' v.bhv_vars{bh} ' {\it n} = ' num2str(v.nsubj) ]}); % generate corresponding rfx figure

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

                                title({[v.bandname{b} ' ' v.sess{se} ' ' v.stages{st} ' ' ] [v.netwname{w} '-' v.netwname{ww} ' x '  v.bhv_vars{bh} ', {\it n} =  ' num2str(v.nsubj) ]});xlabel('connectivity (z-score)');ylabel([v.bhv_vars{bh} ' gains speed']);
                                set(gcf, 'Position',  [100, 100, 420, 320]); %left bottom width height
                                if v.pval_labels == 1
                                    text(0.6, 0.98,['r = ' num2str(round(r_val(w,ww,b,bh,st,se), 3))],'Units','normalized' ,'FontSize',14);
                                    if p_val(w,ww,b,bh,st,se) <0.001 == 1
                                        text(0.6, 0.93, 'p(UC) < 0.001 ' ,'Units','normalized' ,'FontSize',14);
                                    else
                                        text(0.6, 0.93,['p(UC) =  ' num2str(round(p_val(w,ww,b,bh,st,se),3))],'Units','normalized' ,'FontSize',14);
                                    end
                                end
                                if v.subj_labels == 1
                                    labelpoints(x,y,v.subj_label,'E', .05);
                                end
                                saveas(gcf, [scatter_band filesep v.sess{se} '_' v.stages{st} '_' v.bandname{b} '_' v.bhv_vars{bh} '_' v.netwname{w} '_' v.netwname{ww}   '_scatter.png' ]);
                                clear('x','y');
                            end
                        end
                    end
                    close all;
                end
            end
        end
    end

    close all;
    clear i;
end
save(['r_' sess_name],r_val)
end