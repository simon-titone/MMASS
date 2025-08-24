function study2_conn_outliers(v,NETW)
%%%%% Scatterplot of connectivity outliers
%%%%% NETW: network, network, subject, stage, session, band
%%%%% BANDS = 1-delta 2-theta 3-alpha 4-beta 5-gamma 6-bb
%%%%% outliers = (netw, netw, stage, night, band)
conn_scatter = [v.homedir filesep 'results/outliers/conn_scatter_n' num2str(v.nsubj)];
outlierdir = [v.homedir filesep 'results/outliers'];if ~exist(outlierdir); mkdir(outlierdir); end

outliers_low = zeros(v.nnetws, v.nnetws, length(v.stage), length(v.cycle), v.nbands);
outliers_high = zeros(v.nnetws, v.nnetws, length(v.stage), length(v.cycle), v.nbands);
for c = 2:v.ncycles
    for st = v.stage
        for b = 1:v.nbands
            for nt = 1:v.nnetws
                for nt2 = 1:v.nnetws
                    conn=NETW(nt,nt2,:,st,c,b);
                    [maxconn(:,1), maxconn(:,2)] = max(conn); [minconn(:,1), minconn(:,2)] = min(conn);
                    cut = nanstd(conn) * v.out_cutoff; lowcut = nanmean(conn) - cut; highcut = nanmean(conn) + cut;
                    if maxconn(1,1) > highcut
                        outliers_high(nt,nt2,1,c,b) = v.subj(maxconn(1,2)) ;
                    elseif maxconn(1,1) < highcut
                        outliers_high(nt,nt2,1,c,b) = NaN ;
                    end
                    outliers_high(:,:,1,c,b) = triu(outliers_high(:,:,1,c,b));

                    if minconn(1,1) < lowcut
                        outliers_low(nt,nt2,1,c,b) = v.subj(minconn(1,2));
                    elseif minconn(1,1) > lowcut
                        outliers_low(nt,nt2,1,c,b) = NaN ;
                    end
                    outliers_low(:,:,1,c,b) = triu(outliers_low(:,:,1,c,b));

                    if v.scatter == 1
                        if ~exist(conn_scatter); mkdir(conn_scatter); end; cd(conn_scatter);
                        %                         if maxconn(1,1) > highcut | minconn(1,1) < lowcut
                        figure
                        maxval(1,:) = conn(:);
                        l = length(conn);
                        for i = 1:l
                            scatter(repmat(1,1,l),conn(:),'o', 'b');
                        end
                        xlim = [0.5 1.5]; ylim = [(highcut + (highcut/2)) (lowcut - (lowcut/2))];
                        line([0.5 1.5],[lowcut lowcut],'Color','k','LineStyle','--');
                        line([0.5 1.5],[highcut highcut],'Color','k','LineStyle','--');
%                         title([v.bandname{b} ' ' v.t_group ' ' v.netwname{nt} ' ' v.netwname{nt2} ' connectivity' ]);
                        title([v.bandname{b} ' REM cycle ' num2str(c) ' ' v.netwname{nt} ' - ' v.netwname{nt2} ' Conn' ]);
                        ylabel('connectivity (z-score)');
                        labelpoints(1,conn(:),v.subj_label,'E', .05);
                        saveas(gcf, [v.bandname{b} '_c' num2str(c) '_' v.stages{st} '_'  v.netwname{nt} '_' v.netwname{nt2} '_conn_scatter.png' ]);
                    end
                end
                close all;
            end
        end
    end
end
% end
clear r; clear nt; clear nt2; clear b;
close all;
cd(outlierdir);
outliervar = [v.t_group '_outliers_n' num2str(v.nsubj)];
% writematrix(outliers, [outlierdir outliervar '.xls']);
% save(['RS_conn_corr_outliers_n' num2str(n.subj)], 'outlier' );
save(outliervar, 'outliers_low' , 'outliers_high' );
end
