function indv_subj_matrix(NETW,n)
indv_subj_mat = [n.homedir filesep 'results' filesep 'indv'];
if ~exist(indv_subj_mat); mkdir(indv_subj_mat); end; cd(indv_subj_mat);

for s = 1:n.subj
    for b = 1:n.bands
        netw=NETW(:,:,:,:,b);
        
        matrix = netw(:,:,s);
        min1 = min(min(matrix,2));
        indv_min(s,b) = min(min1);
        max1 = max(max(matrix));
        indv_max(s,b,1) = max(max1);
    end
end
