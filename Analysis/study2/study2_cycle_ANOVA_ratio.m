function study2_cycle_ANOVA_ratio(v,NETW)
%%%%% NETW = (network, network, subject, stage, cycle, band)
%%%%% BANDS = 1-delta 2-theta 3-alpha 4-sigma 5-beta 6-gamma
anova_dir = [v.homedir filesep 'results/ANOVA/intercycle/'  v.t_group ];
if ~exist(anova_dir); mkdir(anova_dir); end; cd(anova_dir)
%% Inter-cycle ANOVA
% Calculates ratio of within/between network conn per ppt
for c = 1:3
    for st = v.stage
        for b = 1:v.nbands
            for i = 1:v.nsubj
                for w = 1:v.nnetw
                    withinnetw = NETW(w,w,i,st,c,b);
                    internetw = setdiff(1:v.nnetw,w);
                    for n2 = 1:length(internetw)
                        btwnetw = mean(NETW(internetw,w,i,st,c,b));
                    end
                    ratio(i,b,st,c) = withinnetw/btwnetw;
                end
            end
        end
    end
end

for st= v.stage
    for b = 1:v.nbands
        for w=1:v.nnetws
            vect1 = squeeze(ratio(:,b,st,1));
            vect2 = squeeze(ratio(:,b,st,2));
            vect3 = squeeze(ratio(:,b,st,3));
            nsubj = length(~isnan(mean([vect1 vect2 vect3], 2)));
            mat(:,1) = vect1;mat(:,2) = vect2;mat(:,3) = vect3;
            [tbl,rm] = simple_mixed_anova(mat);
            f = table2array(tbl(3,4));
            p = table2array(tbl(3,5));
            ANOVA_p(b,st) =  p;  %ranovatbl{'(Intercept):Sessions','pValue'};
            ANOVA_F(b,st) =  f;    %ranovatbl{'(Intercept):Sessions','F'};
            clear ranovatbl; clear rm;
        end
    end
end

end