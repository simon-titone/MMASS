function Sleep_create_design(v)

design_matrix.withinlevel=zeros(v.nsubj,v.levels);
design_matrix.Session1 = zeros(v.nsubj,1);
design_matrix.SessANOVA=zeros(v.nsubj,1); % 3 vectors; first specifies session effect; 2nd group and 3rd interaction;
ind=0;level_counter = 0;
for i = 1:v.nsess
    level_counter = level_counter + 1;
    for j=1:v.nsubj
        ind = ind + 1;
        design_matrix.withinlevel(ind,level_counter)=1;
        if i == 1
            design_matrix.Session1(ind,1) = 1;
            design_matrix.SessANOVA(ind,1) = 1;
        elseif i == 2
            design_matrix.SessANOVA(ind,1) = -1;
        end
       
    end
end
v.design_matrix = design_matrix;
cd(v.vardir); save('v.mat');
end