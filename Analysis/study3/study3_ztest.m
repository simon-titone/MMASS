%%% Load file "r_val.mat" for the relevent sample size/comparison 
for s = 1 %stage
    for i = 1:6 % freq
        for w = 1:6 %netw
            for ww = 1:6
                if s == 1; n = 19;elseif s == 2; n = 18; elseif s == 3; n = 9;end
                ctrl = r_val(w,ww,i,1,s,1);
                exp = r_val(w,ww,i,1,s,2);
                [p,z] = corr_rtest(ctrl,exp,n,n);
                p_1tailed(w,ww,i,s) = p(1,1);
                p_2tailed(w,ww,i,s) = p(1,2);
                z_all(w,ww,i,s) = z;
            end
        end
    end
end


for s = 1 %stage
    for i = 1:6 % freq
        pval(:,:,i,s) = triu(p_2tailed(:,:,i,s));
        z_score(:,:,i,s) = triu(z_all(:,:,i,s));
    end
end
save("zscore.mat", "pval", "z_score");