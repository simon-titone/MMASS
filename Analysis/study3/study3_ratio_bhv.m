function study3_ratio_bhv(v)
%%%%% %NETW: network, network, subject, stage, session, cycle, band
%%%%% BANDS = 1-delta 2-theta 3-alpha 4-sigma 5-beta 6-gamma
%% Define groups
%individual sessions
v.EXP_offline = unique([4 12 16 17 3 6 7 20 14 ]);                            % n = 20
v.EXP_offline_REM = unique([1 3 4 7 9 10 11 12 13 16 17 18 21 24 27 28 29]);

v.CTRL_NREM2 = unique([9 14 17 3 6 7 20 14   21]);
v.CTRL_NREM3 = unique([9 14 17 28 3 6 7 20 14  15 22]);
v.CTRL_REM = unique([3 4 6 9 10 11 12 14 17 28 3 6 7 20 14    21]);
%comparing correlations between EXP x offline & CTRL x offline

v.ztest_NREM2 = unique([3 4 6 7 9 12 14 16 17 20]);
v.ztest_NREM3 = unique([3 4 6 7 9 12 14 16 17 20 28]);
v.ztest_REM = unique([1 3 4 6 7 9 10 11 12 13 14 16 17 18 20 21 24 27 28 29]);

%% Select group

v.group = v.EXP_offline; stage = 1:2;  sess = 2;

% v.group = v.EXP_offline_REM; stage = 3; sess_name = 'EXP';sess = 2;

% v.group = v.CTRL_NREM2; stage = 1; sess_name = 'CTRL';sess = 1;
% v.group = v.CTRL_NREM3; stage = 2; sess_name = 'CTRL';sess = 1;
% v.group = v.CTRL_REM; stage = 3; sess_name = 'CTRL';sess = 1;

%% Load subjects
v.subj = 1:29; v.group = sort(unique(v.group), 'descend');
for i = 1:length(v.group); v.subj(v.subj(v.group(i))) = []; end
v.nsubj = length(v.subj); % number of independent subjects
study3_load(v);

cd(v.vardir);load('v.mat');
load('bhv_vars.mat', 'off_ITI');
speed(1,:) = off_ITI(v.subj);
for c2 = 1:v.nsubj; c1(c2) = v.subj(c2); end
subj_label= c1'; v.violin_subj_label = [c1' ; c1'];
bh = 1;
%% Correlation between behavior and sleep during control or experimental nights

w_lvl_corr = [v.homedir filesep 'results/Correlation/ratio_x_BHV' ]; % name dir
if ~exist(w_lvl_corr); mkdir(w_lvl_corr);end; cd(w_lvl_corr); % create dir

% Find ratio of within to between network connectivity per participant
for se = sess
    for st = 1:3
        for b = 1:v.nbands
            for i = 1:v.nsubj
                for w = 1:v.nnetw
                    withinnetw = NETW(w,w,i,st,se,1,b);
                    internetw = setdiff(1:v.nnetw,w);
                    for n2 = 1:length(internetw)
                        btwnetw = mean(NETW(internetw,w,i,st,se,1,b));
                    end
                    ratio(i,b,st,se) = withinnetw/btwnetw;
                end
            end
        end
    end
end

for se = sess
    for st = stage
        for b = 1:v.nbands
            behavior =  speed';
            conn = squeeze(ratio(:,b,st,se));
            [val,px]=corr(conn,behavior);
            r_val(st,b)=val;
            p_val(st,b)=px;
        end
    end
end

end

