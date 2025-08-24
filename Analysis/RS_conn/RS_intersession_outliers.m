function RS_intersession_outliers(NETW,n,subj)
outlierdir = [n.homedir filesep 'results/outliers/']; 
conn_scatter = [n.homedir filesep 'results/outliers/conn__intersess_scatter_n' num2str(n.subj)];

InterSession = NaN(n.netws,n.netws,n.subj,n.bands);
for b = 1:n.bands
    for i = 1:1:n.subj
        for w=1:n.netws
            for ww = 1:n.netws
                InterSession(w,ww,i,b) = NETW(w,ww,i,2,b) - NETW(w,ww,i,1,b);
            end;
        end;
    end;
end

for b = 1:n.bands
    for nt = 1:n.netws
        for nt2 = 1:n.netws
            conn = squeeze(InterSession(nt, nt2,:,b)); [maxconn(:,1), maxconn(:,2)] = max(conn); [minconn(:,1), minconn(:,2)] = min(conn);
            cut = std(conn) * 3; lowcut = mean(conn) - cut; highcut = mean(conn) + cut;
            if maxconn(1,1) > highcut
%                 outliers(nt,nt2,b) = subj(maxconn(1,2)) ;
                outliers_high(nt,nt2,b) = subj(maxconn(1,2)) ;
            elseif maxconn(1,1) < highcut
%                 outliers(nt,nt2,b) = 0 ;
                outliers_high(nt,nt2,b) = 0 ;
            end
            if  minconn(1,1) < lowcut
%                 outliers(nt,nt2,b) = subj(minconn(1,2));
                outliers_low(nt,nt2,b) = subj(minconn(1,2));
            elseif  minconn(1,1) > lowcut
%                 outliers(nt,nt2,b) = 0;
                outliers_low(nt,nt2,b) = 0;
            end
            if n.scatter == 1
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
                    title([n.bandname{b} ' '  n.netwname{nt} ' ' n.netwname{nt2}  ' connectivity' ]);
                    ylabel('connectivity (z-score)');
                    labelpoints(1,conn,n.subj_label,'E', .05);
                    saveas(gcf, [n.bandname{b} '_'  n.netwname{nt} '_' n.netwname{nt2}   '_conn_scatter.png' ]);
                end
            end
        end
    end
end

clear nt; clear nt2; clear b;
close all;
if ~exist(outlierdir); mkdir(outlierdir); end; cd(outlierdir);
outliervar = ['RS_intersession_conn_outliers_n' num2str(n.subj)];
save(outliervar, 'outliers_low' , 'outliers_high' ); 



end