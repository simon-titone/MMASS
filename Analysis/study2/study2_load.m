function study2_load(v)
tmp = zeros(v.nseeds,v.nseeds,v.nf,v.nsubj,v.ncycles,v.nstages);
for i = 1:v.nsubj
    for st= 1:v.nstages
        for c = 1:v.ncycles
            s = v.subj(i);
            file = [v.homedir filesep 'matrix_conn' filesep 'P' num2str(s) '_' v.stages{st} '_c' num2str(c) filesep 'matrix_connectivity.mat'];
            if exist(file); load(file); else; continue;end
            corr_matrix(corr_matrix == 0) = nan;
            tmp(:,:,:,i,st,c) = corr_matrix;
        end
    end
end
clear s
% organize data for following sections
% ------------------------------------
% conn_band = zeros(v.nseeds,v.nseeds,v.nband,v.nsubj,v.nstages,v.ncycles); % average over bands
for b = 1:v.nband
    conn_band(:,:,:,:,:,b) = squeeze(nanmean(tmp(:,:,v.band{b},:,:,:),3));
end
%NETW: network, network, subject, stage, cycle, band
for i = 1:v.nsubj
    for st= 1:v.nstages
        for c = 1:v.ncycles
            for b = 1:v.nbands
                s = v.subj(i);
%                 conn_seeds_band = conn_band(:,:,i,st,c,b);
                for nt = 1:v.nnetw
                    NETW(nt,nt,i, st, c, b) = nanmean(nanmean(conn_band(v.netw{nt},v.netw{nt},i,st,c,b),1),2);
                    internetw = setdiff(1:v.nnetw,nt);
                    for n2 = 1:length(internetw)
                        NETW(nt,internetw(n2),i, st, c, b) = nanmean(nanmean(conn_band(v.netw{nt},v.netw{internetw(n2)},i,st, c, b),1),2);
                        NETW(internetw(n2),nt,i, st, c, b) = NETW(nt,internetw(n2),i,st,c,b);
                    end
                end
            end
        end
    end
end
NETW(NETW == 0) = NaN;

save(v.var, 'NETW','v')
end