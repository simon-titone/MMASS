function scatter_conn_outliers(NETW,subj,n)
%%%%% Scatterplot of connectivity outliers
%%%%% NETW: network, network, subject, stage, session, band
%%%%% BANDS = 1-delta 2-theta 3-alpha 4-beta 5-gamma 6-bb
%%%%% outliers = nan(n.netws,n.netws,n.subj, n.bands, n.sess);
conn_scatter = [n.homedir filesep 'results/outliers/conn_scatter_n' num2str(n.subj)];
outlierdir = [n.homedir 'results/outliers'];if ~exist(outlierdir); mkdir(outlierdir); end;

for b = 1:n.bands
    for r = 1:n.levels
        for s = 1:n.subj
            for nt = 1:n.netws
                for nt2 = 1:n.netws
                    adj_group = setdiff(1:n.subj, s);
                    groupmean(nt,nt2,s, r, b) = nanmean(squeeze(NETW(nt, nt2, adj_group, r, b)));
                end
            end
        end        
    end
end 

for b = 1:n.bands
    for r = 1:n.levels
        for s = 1:n.subj
            indv = NETW(:,:, s, r, b);
            vect = [indv(1,1:n.netws) indv(2,2:n.netws) indv(3,3:n.netws) indv(4,4:n.netws) indv(5,5:n.netws) indv(6,6:n.netws)];
            g = groupmean(:,:,s, r, b);
            groupvect = [g(1,1:n.netws) g(2,2:n.netws) g(3,3:n.netws) g(4,4:n.netws) g(5,5:n.netws) g(6,6:n.netws)];
            [rho, pval] = corr(vect', groupvect');
            corr_r(s, r, b) = rho;
            corr_p(s, r, b) = pval;
        end
    end
end

for b = 1:n.bands
    for r = 1:n.levels
        m(:,r,b) = mean(corr_r(:, r, b), 1); stdev(:,r,b) = std(corr_r(:, r, b),1, 1);
        highcut(r,b) = (m(:,r,b) + (stdev(:,r,b)*2.5));lowcut(r, b) = (m(:,r,b) - (stdev(:,r,b)*2.5));
    end
end

for b = 1:n.bands
    for s = 1:n.subj
        for r = 1:n.levels
            if corr_r(s, r, b) > highcut(r,b) || corr_r(s, r, b) < lowcut(r, b)
                outlier(s, 1, r, b) = subj(s);
                outlier(s, 2, r, b) = corr_p(s, r, b);
                outlier(s, 3, r, b) = corr_r(s, r, b);
                outlier(s, 4, r, b) = b;
                outlier(s, 5, r, b) = r;
            else 
                outlier(s, 1, r, b) = 0;
                outlier(s, 2, r, b) = 0;
                outlier(s, 3, r, b) = 0;
                outlier(s, 4, r, b) = 0;
                outlier(s, 5, r, b) = 0;
            end
        end
    end
end  

for b = 1:n.bands
    netw=NETW(:,:,:,:,b);
    for r = 1:2
        for nt = 1:n.netws
            for nt2 = 1:n.netws
                if r == 1; run = 1:n.subj; end; if r ==2; run = n.subj +1:n.subj*2;end
                conn = netw(nt, nt2, run); [maxconn(:,1), maxconn(:,2)] = max(conn); [minconn(:,1), minconn(:,2)] = min(conn);
                cut = std(conn) * n.out_cutoff; lowcut = mean(conn) - cut; highcut = mean(conn) + cut;
                if maxconn(1,1) > highcut
                    outliers(nt,nt2,r,b) = subj(maxconn(1,2)) ;
                    outliers_high(nt,nt2,r,b) = subj(maxconn(1,2)) ;
                elseif maxconn(1,1) < highcut
                    outliers(nt,nt2,r,b) = 0 ;
                    outliers_high(nt,nt2,r,b) = 0 ;
                end
                if minconn(1,1) < lowcut
                    outliers(nt,nt2,r,b) = subj(minconn(1,2));
                    outliers_low(nt,nt2,r,b) = subj(minconn(1,2));
                elseif maxconn(1,1) < highcut
                    outliers(nt,nt2,r,b) = 0 ;
                    outliers_low(nt,nt2,r,b) = 0 ;
                end
                if n.scatter == 1
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
                        title([n.bandname{b} ' '  n.netwname{nt} ' ' n.netwname{nt2} ' run ' num2str(r)  ' connectivity' ]);
                        ylabel('connectivity (z-score)');
                        labelpoints(1,conn,n.subj_label,'E', .05);
                        saveas(gcf, [n.bandname{b} '_'  n.netwname{nt} '_' n.netwname{nt2} '_run' num2str(r)  '_conn_scatter.png' ]);
                    end
                end
            end
        end
    end
end
clear r; clear nt; clear nt2; clear b;
close all;
cd(outlierdir);
outliervar = ['RS_conn_outliers_n' num2str(n.subj)];
% save(['RS_conn_corr_outliers_n' num2str(n.subj)], 'outlier' ); 
save(outliervar, 'outliers', 'outliers_low' , 'outliers_high' ); 
end