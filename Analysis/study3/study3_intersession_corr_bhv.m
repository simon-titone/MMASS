function study3_intersession_corr(v,NETW)
%%%%% NETW: network, network, subject, stage, session,cycle, band
%% Load Subjects
v.intersess_corr = [4 9 12 14 16 17 15 21   3 6 7 20];  % n = 20
v.intersess_corr_REM = unique([1 3 4 6 7 9 10 11 12 13 14 16 17 18 21 24 27 28 29]); 

% v.group = v.intersess_corr; stage = 1:2;
v.group = v.intersess_corr_REM; stage = 3;

v.subj = 1:29; v.group = sort(v.group, 'descend');
for i = 1:length(v.group); v.subj(v.subj(v.group(i))) = []; end
v.nsubj = length(v.subj); % number of independent subjects
study3_load(v);

cd(v.vardir);load('v.mat');
load('bhv_vars.mat', 'off_ITI');
v.speed(1,:) = off_ITI(v.subj);

for cc = 1:v.nsubj; c(cc) = v.subj(cc); end
v.subj_label= c'; v.violin_subj_label = [c' ; c'];
intersession_corr_dir = [v.homedir filesep 'results' filesep 'intersession_corr_n' num2str(v.nsubj)];

%% calculates difference between sessions, per stage (e.g. EXP NREM2 - CTRL NREM2)
for se = 1:v.nsess
    for c= 1:v.ncycles
        for st= 1:v.nstages
            for b = 1:v.nbands
                for i = 1:v.nsubj
                    for w=1:v.nnetws
                        for ww=1:v.nnetws
                            if NETW(w,ww,i,st,1,c,b) ~= 0 && NETW(w,ww,i,st,2,c,b) ~= 0
                                intersess_NETW(w,ww,se,c,st,b,i) = NETW(w,ww,i,st,se,c,b);
                            else
                                intersess_NETW(w,ww,se,c,st,b,i) = NaN;
                            end
                        end
                    end
                end
            end
        end
    end
end

for st= 1:v.nstages
    for c= 1:v.ncycles
        for b = 1:v.nbands
            for i = 1:1:v.nsubj
                for w=1:v.nnetws
                    for ww = 1:v.nnetws
                        intersess_conn(w,ww,i,st,c,b) = intersess_NETW(w,ww,2,c,st,b,i) - intersess_NETW(w,ww,1,c,st,b,i);
                    end
                end
            end
        end
    end
end

for st= stage
    for c= 1:v.ncycles
        for b = 1:v.nbands
            for w=1:v.nnetws
                for ww=1:v.nnetws
                    behavior = v.speed';
                    conn_vect = squeeze(intersess_conn(w,ww,:,st,c,b));
                    n_nan = find(isnan(conn_vect));
%                     if n_nan >1
%                         return
%                     end
                    n_nan = sort(n_nan, 'descend');
                    for nn = 1:length(n_nan)
                        behavior(n_nan(nn)) = [];
                        conn_vect(n_nan(nn)) = [];
                    end
                    
                    [val,px]=corr(conn_vect,behavior);
                    r_val(w,ww,st,c,b)=val;
                    p_val(w,ww,st,c,b)=px;
                    
                end
            end
            bhv(st,c,b,:) = behavior;
        end
    end
end

clear val; clear px; clear w; clear ww; clear vect;
for c= 1:v.ncycles
    for st= stage
        intersession_stage_dir =[intersession_corr_dir filesep v.stages{st}];
        if ~exist(intersession_stage_dir); mkdir(intersession_stage_dir); end; cd(intersession_stage_dir)
        
        for b=1:v.nbands
            pval_fdr = p_val(:,:,st,c,b);
            if v.FDR == 'mot'
                p_val_temp_FDR  = p_val(5,:);
            elseif v.FDR == 'all'
                p_val_temp_FDR = [ p_val(1:v.nnetw)  p_val(v.nnetw+2:v.nnetw*2)  p_val(v.nnetw*2+3:v.nnetw*3)  p_val(v.nnetw*3+4:v.nnetw*4)  p_val(v.nnetw*4+5:v.nnetw*5)  p_val(v.nnetw*5+6:v.nnetw*6)];
            end
            [~, dummy2, ~, dummy4]=fdr_bh(p_val_temp_FDR, v.fdr_thres,'pdep','no');
            crit_p.Behavcov = dummy2;
            adj_p.Behavcov(1:length(dummy4)) = dummy4;
            clear dummy2; clear dummy4;
            
            [ivect,jvect]=find(p_val(:,:,st,c,b)<=crit_p.Behavcov);
            [ivect_unc,jvect_unc]=find(p_val(:,:,st,c,b)<=v.fdr_thres);
            figure(); set(gcf,'Color','white'); box OFF;
            
            matrix = r_val(:,:,st,c,b); tri_matrix = triu(matrix);
            
            imagesc(tri_matrix); axis square; % colormap(v.colors); daspect([1 1 1]); hold on
            h = colorbar; set(get(h,'title'),'string','r-val', 'fontsize', 12); %ylabel(h, 'f-val');
            dim = set(gcf, 'Position',  [1000, 900, 500, 450]); %left bottom width height
            colormap('redblue'); limit = max(max(abs(tri_matrix))); caxis([-limit limit]);
            
            if v.fig_labels ==1
            T = {[  'EXP-CTRL connectivity x offline'   ] [v.stages{st} ' ' v.bandname{b}  ', {\it n} = ' num2str(v.nsubj)] };
            else
            T = {[  'EXP-CTRL connectivity x offline'   ] [v.stages{st} ' ' v.bandname{b} ] };
            end
            
            title(T);
            
            for a2=1:v.nnetws-1; line([a2+0.5 a2+0.5],[a2+1.5 0.5],'color','k','LineWidth',0.8);line([a2-0.5 v.nnetws+0.5],[a2+0.5 a2+0.5],'color','k','LineWidth',0.8); end
            for b2=1:v.nnetws; line([b2-0.5 b2-0.5],[0.5 v.nnetws+0.5],'color','k','LineWidth',0.001); line([0.5 v.nnetws+0.5],[b2-0.5 b2-0.5],'color','k','LineWidth',0.001);end
            clear a2; clear b2;
            
            set(gca,'YTick',[1:v.nnetws],'YTickLabel',v.netwname,'Fontsize',16);
            set(gca,'XTick',[1:v.nnetws],'XTickLabel',v.netwname,'Fontsize',16);
            hold on;
            scatter(jvect_unc,ivect_unc, 80, 'w', 'o' ,'ok','LineWidth',2.5); clear ivect_unc; clear jvect_unc;
            scatter(jvect,ivect,80,'w','ok','filled' ,'LineWidth',2.5); clear ivect; clear jvect;
            
            savename = [intersession_stage_dir filesep v.bandname{b} '_' v.stages{st} '_intersess_corr' '.png'];
            saveas(gcf, savename);
            
            if v.scatter == 1
                % Creates scatter plots
                reset(gcf); reset(gca);
                for w=1:v.nnetws
                    for ww=1:v.nnetws
                        if p_val(w,ww,st,c,b) <0.05 == 1
                            scatter_band = [intersession_stage_dir  '/scatter/'  v.bandname{b}];
                            if ~exist(scatter_band); mkdir(scatter_band); end
                            
                            figure
                            x = squeeze(intersess_conn(w,ww,:,st,c,b)); y = behavior;
                            scatter(x,y,[], 'b', 'filled');lsline;
                            mdl = fitlm(x, y, 'linear');
                            [y_pred,y_ci] = predict(mdl,x, 'Alpha', 0.05);
                            Y = [x, y_ci,y_pred];
                            Y = sortrows(Y,1);
                            eb = [Y(:,4)-Y(:,2)];
                            boundedline(Y(:,1),Y(:,4),eb, '-r','alpha');
                            clear mdl; clear y_pred; clear y_ci; clear Y; clear eb;
                            if v.fig_labels ==1
                                T_scatter = {[  'EXP-CTRL connectivity x offline' ', {\it n} = ' num2str(v.nsubj)] [v.stages{st} ' ' v.bandname{b} ' ' v.netwname{w} '-' v.netwname{ww} ] };
                            else
                                T_scatter = {[  'EXP-CTRL connectivity x offline' ] [v.stages{st} ' ' v.bandname{b} ' ' v.netwname{w} '-' v.netwname{ww} ] };
                            end
                            title(T_scatter); xlabel('connectivity (z-score)'); ylabel([ 'offline gains speed']);
                            set(gcf, 'Position',  [100, 100, 420, 320]); %left bottom width height
                            if v.subj_labels == 1
                                labelpoints(x,y,v.subj_label,'E', .05);
                            end
                            if v.pval_labels == 1
                                text(0.6, 0.98,['r = ' num2str(round(r_val(w,ww,st,c,b), 3))],'Units','normalized' ,'FontSize',14);
                                if p_val(w,ww,st,c,b) <0.001 == 1
                                    text(0.6, 0.93, 'p(UC) < 0.001 ' ,'Units','normalized' ,'FontSize',14);
                                else
                                    text(0.6, 0.93,['p(UC) =  ' num2str(round(p_val(w,ww,st,c,b),3))],'Units','normalized' ,'FontSize',14);
                                end
                            end
                            S_scatter = [ scatter_band filesep   v.stages{st} '_' v.bandname{b}  '_offline_' v.netwname{w} '_' v.netwname{ww}   '_scatter.png' ];
                            saveas(gcf, S_scatter);
                            clear('x','y');
                        end
                    end
                end
            end
        end
        close all
    end
end