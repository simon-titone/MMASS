function study2_excl(v)
base_excl = [9 14 17];
%%%% Inter-cycle groups
NREM2_c1_c2_excl = [base_excl unique([21 13 28])];                                           
NREM2_c2_c3_excl = [base_excl unique([13 28 5 23])];    
NREM2_c1_c3_excl = [base_excl unique([21 5 23])];  
NREM2_all_excl = unique([NREM2_c1_c2_excl,NREM2_c2_c3_excl,NREM2_c1_c3_excl]);

NREM3_c1_c2_excl = [base_excl unique([28 15 22 12 24 25])];     
NREM3_c2_c3_excl = [base_excl unique([12 24 25 2 3 5 7 11 13 15 19 28 27])]; 
NREM3_c1_c3_excl = [base_excl unique([2 3 5 7 11 13 15 19 28 28 15 22])];     
NREM3_all_excl = unique([NREM3_c1_c2_excl,NREM3_c2_c3_excl,NREM3_c1_c3_excl]);

REM_c1_c2_excl = [base_excl unique([3 4 6 10 11 12 28 21 4 18 28 21 6])];     
REM_c2_c3_excl = [base_excl unique([4 18 28 21 6 1 2 4 8 10 15])]; 
REM_c1_c3_excl = [base_excl unique([1 2 4 8 10 15 3 4 6 10 11 12 28 21])];      
REM_all_excl = unique([REM_c1_c2_excl,REM_c2_c3_excl,REM_c1_c3_excl]);

%%%% Inter-stage groups
c1_N2_N3_excl = [base_excl 15 21 22 28   ];
c1_N2_REM_excl = [base_excl 3 4 6 10 11 12 28 21];
c1_N3_REM_excl = [base_excl 3 4 6 10 11 12 28 15 21 22];
c1_all_excl = unique([c1_N2_N3_excl,c1_N2_REM_excl,c1_N3_REM_excl]);

c2_N2_N3_excl = [base_excl 12 24  13 1 20];
c2_N2_REM_excl = unique([base_excl 4 18 28 21 6 13]);
c2_N3_REM_excl = unique([base_excl 4 18 28 21 6 12 24 25]);
c2_all_excl = unique([c2_N2_N3_excl,c2_N2_REM_excl,c2_N3_REM_excl]);

c3_N2_N3_excl = unique([base_excl 2 3 5 7 11 13 15 19 28 27 23]);
c3_N2_REM_excl = unique([base_excl 1 2 4 8 10 15 5 23]);
c3_N3_REM_excl = unique([base_excl 2 3 5 7 11 13 15 19 28 27 1 2 4 8 10 15]);
c3_all_excl = unique([c3_N2_N3_excl,c3_N2_REM_excl,c3_N3_REM_excl]);



%%% Defining inter-cycle sample 
v.n_intercycle_groups = 12; v.n_interstage_groups = 12;
v.subj_IC = NaN(v.n_intercycle_groups,26); v.subj_IS = NaN(v.n_interstage_groups,26);

v.NREM2_c1_c2 = 1:29; NREM2_c1_c2_excl = sort(NREM2_c1_c2_excl, 'descend'); for s = 1:length(NREM2_c1_c2_excl); v.NREM2_c1_c2(v.NREM2_c1_c2(NREM2_c1_c2_excl(s))) = []; end
for i = 1:length(v.NREM2_c1_c2); v.subj_IC(1,i) = v.NREM2_c1_c2(1,i); end; clear i 
v.NREM2_c2_c3 = 1:29; NREM2_c2_c3_excl = sort(NREM2_c2_c3_excl, 'descend'); for s = 1:length(NREM2_c2_c3_excl); v.NREM2_c2_c3(v.NREM2_c2_c3(NREM2_c2_c3_excl(s))) = []; end
for i = 1:length(v.NREM2_c2_c3); v.subj_IC(2,i) = v.NREM2_c2_c3(1,i);end; clear i 
v.NREM2_c1_c3 = 1:29; NREM2_c1_c3_excl = sort(NREM2_c1_c3_excl, 'descend'); for s = 1:length(NREM2_c1_c3_excl); v.NREM2_c1_c3(v.NREM2_c1_c3(NREM2_c1_c3_excl(s))) = []; end
for i = 1:length(v.NREM2_c1_c3); v.subj_IC(3,i) = v.NREM2_c1_c3(1,i);end; clear i 


v.NREM3_c1_c2 = 1:29; NREM3_c1_c2_excl = sort(NREM3_c1_c2_excl, 'descend'); for s = 1:length(NREM3_c1_c2_excl); v.NREM3_c1_c2(v.NREM3_c1_c2(NREM3_c1_c2_excl(s))) = []; end
for i = 1:length(v.NREM3_c1_c2); v.subj_IC(4,i) = v.NREM3_c1_c2(1,i);end; clear i 
v.NREM3_c2_c3 = 1:29; NREM3_c2_c3_excl = sort(NREM3_c2_c3_excl, 'descend'); for s = 1:length(NREM3_c2_c3_excl); v.NREM3_c2_c3(v.NREM3_c2_c3(NREM3_c2_c3_excl(s))) = []; end
for i = 1:length(v.NREM3_c2_c3); v.subj_IC(5,i) = v.NREM3_c2_c3(1,i);end; clear i 
v.NREM3_c1_c3 = 1:29; NREM3_c1_c3_excl = sort(NREM3_c1_c3_excl, 'descend'); for s = 1:length(NREM3_c1_c3_excl); v.NREM3_c1_c3(v.NREM3_c1_c3(NREM3_c1_c3_excl(s))) = []; end
for i = 1:length(v.NREM3_c1_c3); v.subj_IC(6,i) = v.NREM3_c1_c3(1,i);end; clear i 


v.REM_c1_c2 = 1:29; REM_c1_c2_excl = sort(REM_c1_c2_excl, 'descend'); for s = 1:length(REM_c1_c2_excl); v.REM_c1_c2(v.REM_c1_c2(REM_c1_c2_excl(s))) = []; end
for i = 1:length(v.REM_c1_c2); v.subj_IC(7,i) = v.REM_c1_c2(1,i);end; clear i 
v.REM_c2_c3 = 1:29; REM_c2_c3_excl = sort(REM_c2_c3_excl, 'descend'); for s = 1:length(REM_c2_c3_excl); v.REM_c2_c3(v.REM_c2_c3(REM_c2_c3_excl(s))) = []; end
for i = 1:length(v.REM_c2_c3); v.subj_IC(8,i) = v.REM_c2_c3(1,i);end; clear i 
v.REM_c1_c3 = 1:29; REM_c1_c3_excl = sort(REM_c1_c3_excl, 'descend'); for s = 1:length(REM_c1_c3_excl); v.REM_c1_c3(v.REM_c1_c3(REM_c1_c3_excl(s))) = []; end
for i = 1:length(v.REM_c1_c3); v.subj_IC(9,i) = v.REM_c1_c3(1,i);end; clear i 

v.NREM2_all = 1:29; NREM2_all_excl = sort(NREM2_all_excl, 'descend'); 
for s = 1:length(NREM2_all_excl); v.NREM2_all(v.NREM2_all(NREM2_all_excl(s))) = []; end
for i = 1:length(v.NREM2_all); v.subj_IC(10,i) = v.NREM2_all(1,i);end; clear i 

v.NREM3_all = 1:29;  NREM3_all_excl = sort( NREM3_all_excl, 'descend'); 
for s = 1:length( NREM3_all_excl); v.NREM3_all(v.NREM3_all( NREM3_all_excl(s))) = []; end
for i = 1:length(v.NREM3_all); v.subj_IC(11,i) = v.NREM3_all(1,i);end; clear i 

v.REM_all = 1:29; REM_all_excl = sort(REM_all_excl, 'descend'); 
for s = 1:length(REM_all_excl); v.REM_all(v.REM_all(REM_all_excl(s))) = []; end
for i = 1:length(v.REM_all); v.subj_IC(12,i) = v.REM_all(1,i);end; clear i 


v.cycle_IC = [1 2 3 1 2 3 1 2 3 nan nan nan];
v.stage_IC = [1 1 1 2 2 2 3 3 3 1 2 3];
v.group_name_IC = {'NREM2_c1_c2', 'NREM2_c2_c3', 'NREM2_c1_c3', 'NREM3_c1_c2', 'NREM3_c2_c3', 'NREM3_c1_c3', 'REM_c1_c2', 'REM_c2_c3', 'REM_c1_c3', 'NREM2_all', 'NREM3_all', 'REM_all'};

%%%% Inter-stage groups
v.c1_N2_N3 = 1:29; c1_N2_N3_excl = sort(c1_N2_N3_excl, 'descend'); for s = 1:length(c1_N2_N3_excl); v.c1_N2_N3(v.c1_N2_N3(c1_N2_N3_excl(s))) = []; end
for i = 1:length(v.c1_N2_N3); v.subj_IS(1,i) = v.c1_N2_N3(1,i);end; clear i 

v.c1_N2_REM = 1:29; c1_N2_REM_excl = sort(c1_N2_REM_excl, 'descend'); for s = 1:length(c1_N2_REM_excl); v.c1_N2_REM(v.c1_N2_REM(c1_N2_REM_excl(s))) = []; end
for i = 1:length(v.c1_N2_REM); v.subj_IS(2,i) = v.c1_N2_REM(1,i);end; clear i 

v.c1_N3_REM = 1:29; c1_N3_REM_excl = sort(c1_N3_REM_excl, 'descend'); for s = 1:length(c1_N3_REM_excl); v.c1_N3_REM(v.c1_N3_REM(c1_N3_REM_excl(s))) = []; end
for i = 1:length(v.c1_N3_REM); v.subj_IS(3,i) = v.c1_N3_REM(1,i);end; clear i 

v.c2_N2_N3 = 1:29; c2_N2_N3_excl = sort(c2_N2_N3_excl, 'descend'); for s = 1:length(c2_N2_N3_excl); v.c2_N2_N3(v.c2_N2_N3(c2_N2_N3_excl(s))) = []; end
for i = 1:length(v.c2_N2_N3); v.subj_IS(4,i) = v.c2_N2_N3(1,i);end; clear i 

v.c2_N2_REM = 1:29; c2_N2_REM_excl = sort(c2_N2_REM_excl, 'descend'); for s = 1:length(c2_N2_REM_excl); v.c2_N2_REM(v.c2_N2_REM(c2_N2_REM_excl(s))) = []; end
for i = 1:length(v.c2_N2_REM); v.subj_IS(5,i) = v.c2_N2_REM(1,i);end; clear i 

v.c2_N3_REM = 1:29; c2_N3_REM_excl = sort(c2_N3_REM_excl, 'descend'); for s = 1:length(c2_N3_REM_excl); v.c2_N3_REM(v.c2_N3_REM(c2_N3_REM_excl(s))) = []; end
for i = 1:length(v.c2_N3_REM); v.subj_IS(6,i) = v.c2_N3_REM(1,i);end; clear i 

v.c3_N2_N3 = 1:29; c3_N2_N3_excl = sort(c3_N2_N3_excl, 'descend'); for s = 1:length(c3_N2_N3_excl); v.c3_N2_N3(v.c3_N2_N3(c3_N2_N3_excl(s))) = []; end
for i = 1:length(v.c3_N2_N3); v.subj_IS(7,i) = v.c3_N2_N3(1,i);end; clear i 

v.c3_N2_REM = 1:29; c3_N2_REM_excl = sort(c3_N2_REM_excl, 'descend'); for s = 1:length(c3_N2_REM_excl); v.c3_N2_REM(v.c3_N2_REM(c3_N2_REM_excl(s))) = []; end
for i = 1:length(v.c3_N2_REM); v.subj_IS(8,i) = v.c3_N2_REM(1,i);end; clear i 

v.c3_N3_REM = 1:29; c3_N3_REM_excl = sort(c3_N3_REM_excl, 'descend'); for s = 1:length(c3_N3_REM_excl); v.c3_N3_REM(v.c3_N3_REM(c3_N3_REM_excl(s))) = []; end
for i = 1:length(v.c3_N3_REM); v.subj_IS(9,i) = v.c3_N3_REM(1,i);end; clear i 

v.c1_all = 1:29; c1_all_excl = sort(c1_all_excl, 'descend'); for s = 1:length(c1_all_excl); v.c1_all(v.c1_all(c1_all_excl(s))) = []; end
for i = 1:length(v.c1_all); v.subj_IS(10,i) = v.c1_all(1,i);end; clear i 

v.c2_all = 1:29; c2_all_excl = sort(c2_all_excl, 'descend'); for s = 1:length(c2_all_excl); v.c2_all(v.c2_all(c2_all_excl(s))) = []; end
for i = 1:length(v.c2_all); v.subj_IS(11,i) = v.c2_all(1,i);end; clear i 

v.c3_all = 1:29; c3_all_excl = sort(c3_all_excl, 'descend'); for s = 1:length(c3_all_excl); v.c3_all(v.c3_all(c3_all_excl(s))) = []; end
for i = 1:length(v.c3_all); v.subj_IS(12,i) = v.c3_all(1,i);end; clear i 


v.cycle_IS = [1 1 1 2 2 2 3 3 3 1 2 3];
v.stage_IS = [1 2 3 1 2 3 1 2 3 nan nan nan];
v.group_name_IS = {'c1_N2_N3', 'c1_N2_REM', 'c1_N3_REM', 'c2_N2_N3', 'c2_N2_REM', 'c2_N3_REM', 'c3_N2_N3', 'c3_N2_REM', 'c3_N3_REM','c1_all','c2_all','c3_all'};

save([v.vardir '/v.mat'], 'v')

end