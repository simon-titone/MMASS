function Sleep_conn_outliers(v,NETW)
%%%%% Scatterplot of connectivity outliers
%%%%% %NETW: network, network, subject, stage, session, cycle, band
%%%%% BANDS = 1-delta 2-theta 3-alpha 4-beta 5-gamma 6-bb
%%%%% outliers = (netw, netw, stage, night, band)
%% Load Subjects
v.subj = 1:29; v.group = sort(v.group, 'descend');
for i = 1:length(v.group); v.subj(v.subj(v.group(i))) = []; end
Sleep_load(v);

%% Outlier analysis
conn_scatter = [v.homedir filesep 'results/outliers/conn_scatter_n' num2str(v.nsubj)];
outlierdir = [v.homedir filesep 'results/outliers'];if ~exist(outlierdir); mkdir(outlierdir); end

for se = 1:v.nsess
    for c = 1:v.ncycles
        for st = 1:v.nstages
            for b = 1:v.nbands
                for nt = 1:v.nnetws
                    for nt2 = 1:v.nnetws
                        conn=NETW(nt,nt2,:,st,se,1,b);
                        [maxconn(:,1), maxconn(:,2)] = max(conn); [minconn(:,1), minconn(:,2)] = min(conn);
                        cut = std(conn) * v.out_cutoff; lowcut = mean(conn) - cut; highcut = mean(conn) + cut;
                        if maxconn(1,1) > highcut
                            outliers(nt,nt2,st,se,b) = v.subj(maxconn(1,2)) ;
                            outliers_high(nt,nt2,st,se,b) = v.subj(maxconn(1,2)) ;
                        elseif maxconn(1,1) < highcut
                            outliers(nt,nt2,st,se,b) = 0 ;
                            outliers_high(nt,nt2,st,se,b) = 0 ;
                        end
                        if minconn(1,1) < lowcut
                            outliers(nt,nt2,st,se,b) = v.subj(minconn(1,2));
                            outliers_low(nt,nt2,st,se,b) = v.subj(minconn(1,2));
                        elseif maxconn(1,1) < highcut
                            outliers(nt,nt2,st,se,b) = 0 ;
                            outliers_low(nt,nt2,st,se,b) = 0 ;
                        end
                        if v.scatter == 1
                            if ~exist(conn_scatter); mkdir(conn_scatter); end; cd(conn_scatter);
                            if maxconn(1,1) > highcut | minconn(1,1) < lowcut
                                figure
                                maxval(1,:) = conn(:);
                                for i = 1:length(conn)
                                    scatter(repmat(1,1,length(conn)),conn,'o', 'b');
                                end
                                xlim = [0.5 1.5]; ylim = [(highcut + (highcut/2)) (lowcut - (lowcut/2))];
                                line([0.5 1.5],[lowcut lowcut],'Color','k','LineStyle','--');
                                line([0.5 1.5],[highcut highcut],'Color','k','LineStyle','--');
                                title([ v.sess{se} ' cycle ' num2str(c) ' ' v.stages{st} ' ' v.bandname{b}  ' ' v.netwname{nt} ' ' v.netwname{nt2} ' connectivity' ]);
                                ylabel('connectivity (z-score)');
                                labelpoints(1,conn,v.subj_label,'E', .05);
                                saveas(gcf, [v.sess{se} '_c'  num2str(c) '_' v.stages{st} '_' v.bandname{b} '_'  v.netwname{nt} '_' v.netwname{nt2} '_conn_scatter.png' ]);
                            end
                        end
                    end
                end
            end
        end
    end
end
clear r; clear nt; clear nt2; clear b;
close all;
cd(outlierdir);
outliervar = ['sleep_conn_outliers_n' num2str(v.nsubj)];
% save(['RS_conn_corr_outliers_n' num2str(n.subj)], 'outlier' );
save(outliervar, 'outliers', 'outliers_low' , 'outliers_high' );
end
