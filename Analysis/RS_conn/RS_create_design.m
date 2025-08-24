function RS_create_design(n)

design_matrix.withinlevel=zeros(n.subj,n.levels);
design_matrix.Session1 = zeros(n.subj,1);
design_matrix.SessANOVA=zeros(n.subj,1); % 3 vectors; first specifies session effect; 2nd group and 3rd interaction;
ind=0;level_counter = 0;
for i = 1:n.sess
    %      for ii = 1:n.group
    level_counter = level_counter + 1;
    for j=1:n.subj
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
cd(n.vardir); save('design.mat','design_matrix');
end