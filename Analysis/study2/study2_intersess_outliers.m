function study2_intersess_outliers(NETW,v)
outlierdir = [v.homedir filesep 'results/outliers/']; 
conn_scatter = [v.homedir filesep 'results/outliers/IS_scatter_n' num2str(v.nsubj) '_' v.t_group];

InterSession = NaN(v.nnetws,v.nnetws,v.nsubj,v.nstages, v.ncycles,v.nbands);
for c = v.cycle
    for st= v.stage
        for b = 1:v.nbands
            for i = 1:1:v.nsubj
                for w=1:v.nnetws
                    for ww = 1:v.nnetws
                        if c == 1
                            InterSession(w,ww,i,st,c,b) = NETW(w,ww,i,st,1,b) - NETW(w,ww,i,st,2,b);% t_intersess = 'c1_c2';
                        elseif c == 2
                            InterSession(w,ww,i,st,c,b) = NETW(w,ww,i,st,2,b) - NETW(w,ww,i,st,3,b);% t_intersess = 'c2_c3';
                        elseif c == 3
                            InterSession(w,ww,i,st,c,b) = NETW(w,ww,i,st,1,b) - NETW(w,ww,i,st,3,b);% t_intersess = 'c1_c3';
                        end
                    end;
                end;
            end;
        end
    end
end

for c = v.cycle
    for st= v.stage
        for b = 1:v.nbands
            for w = 1:v.nnetws
                for ww = 1:v.nnetws
                    conn = squeeze(InterSession(w,ww,:,st,c,b)); [maxconn(:,1), maxconn(:,2)] = max(conn); [minconn(:,1), minconn(:,2)] = min(conn);
                    cut = std(conn) * 3; lowcut = mean(conn) - cut; highcut = mean(conn) + cut;
                    if maxconn(1,1) > highcut
                        outliers_high(w,ww,st,c,b) = v.subj(maxconn(1,2)) ;
                    elseif maxconn(1,1) < highcut
                        outliers_high(w,ww,st,c,b) = 0 ;
                    end
                    if  minconn(1,1) < lowcut
                        outliers_low(w,ww,st,c,b) = v.subj(minconn(1,2));
                    elseif  minconn(1,1) > lowcut
                        outliers_low(w,ww,st,c,b) = 0;
                    end
                    if v.scatter == 1
                        if ~exist(conn_scatter); mkdir(conn_scatter); end; cd(conn_scatter);
                        if maxconn > highcut | minconn < lowcut
                            figure
                            maxval(1,:) = conn(:);
                            for i = 1:length(conn)
                                scatter(repmat(1,1,length(conn)),conn,'o', 'b');
                            end
                            xlim = [0.5 1.5]; ylim = [(highcut + (highcut/2)) (lowcut - (lowcut/2))];
                            line([0.5 1.5],[lowcut lowcut],'Color','k','LineStyle','--');
                            line([0.5 1.5],[highcut highcut],'Color','k','LineStyle','--');
                            title([v.stages{st} ' c' num2str(c) ' ' v.bandname{b} ' '  v.netwname{w} ' ' v.netwname{ww}  ' connectivity' ]);
                            ylabel('connectivity (z-score)');
                            labelpoints(1,conn,v.subj_label,'E', .05);
                            saveas(gcf, [v.stages{st} '_c' num2str(c) '_' v.bandname{b} '_'  v.netwname{w} '_' v.netwname{ww}   '_conn_scatter.png' ]);
                        end
                    end
                end
            end
        end
    end
end
clear w; clear ww; clear b;
close all;
if ~exist(outlierdir); mkdir(outlierdir); end; cd(outlierdir);
outliervar = [  v.t_group '_IS_outliers_n' num2str(v.nsubj)];
save(outliervar, 'outliers_low' , 'outliers_high' ); 
end