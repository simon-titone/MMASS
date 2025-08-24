function RS_load(n, subj)
tmp = zeros(n.seeds,n.seeds,n.nf,n.subj,n.sess);
for i = 1:n.subj
    for r = 1:n.sess
        s = subj(i);
        file = [n.homedir  'matrix_conn' filesep 'subj' num2str(s) '_run' num2str(r) filesep 'matrix_connectivity.mat'];
        if exist(file) ; load(file);
        else
            corr_matrix = zeros(n.seeds,n.seeds,n.nf) ;
        end
        corr_matrix(corr_matrix == 0) = nan;
        tmp(:,:,:,i,r) = corr_matrix;
    end
end
clear s
% organize data for following sections
% ------------------------------------
conn_band = zeros(n.seeds,n.seeds,n.subj,n.sess,n.bands); % average over bands
for b = 1:n.bands
    conn_band(:,:,:,:,b) = squeeze(nanmean(tmp(:,:,n.band{b},:,:),3));
end
% NETW = zeros(n.nnetw,n.nnetw,i,r,b);
for i = 1:n.subj
    for r = 1:n.sess
        for b = 1:n.bands
            s = subj(i);
            conn_seeds_band = conn_band(:,:,i,r,b);
            for nt = 1:n.nnetw
                NETW(nt,nt,i, r, b) = nanmean(nanmean(conn_band(n.netw{nt},n.netw{nt},i,r,b),1),2);
                internetw = setdiff(1:n.nnetw,nt);
                for n2 = 1:length(internetw)
                    NETW(nt,internetw(n2),i, r, b) = nanmean(nanmean(conn_band(n.netw{nt},n.netw{internetw(n2)},i,r, b),1),2);
                    NETW(internetw(n2),nt,i, r, b) = NETW(nt,internetw(n2),i,r,b);
                end
            end
        end
    end
end
cd(n.vardir); save('n.mat','n', 'NETW')
end