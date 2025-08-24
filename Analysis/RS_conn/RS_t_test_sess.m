function t_test_sess(NETW,n)
%%%%%% NETW(network, network, subject, run, band);
% T-TEST OF DIFFERENCE BETWEEN SESSIONS

sesseffect = [n.homedir filesep 'results' filesep 'sess_n' num2str(n.subj)];
if ~exist(sesseffect); mkdir(sesseffect); end; cd(sesseffect);

diff = zeros(n.netws, n.netws, n.subj,n.bands);
for b = 1:n.bands
    netw=NETW(:,:,:,:,b);
    diff(:,:,1:n.subj,b) = minus(netw(:,:,1:n.subj), netw(:,:,n.subj+1:n.subj*n.levels));
end

% for b = 1:n.bands
%     for nt = 1:n.netws
%         for nt2 = 1:n.netws
%             [h,p,ci,stats] = ttest(diff(nt, nt2, :, b));
%             diff_h(nt, nt2,b) = h;
%             diff_p(nt, nt2,b) = p;
%             diff_t(nt, nt2,b) = stats.tstat;
%         end
%     end
% end

for b = 1:n.bands
    for nt = 1:n.netws
        for nt2 = 1:n.netws
            r1 = NETW(nt,nt2,:,1,b);r2 = NETW(nt,nt2,:,2,b);
            [h,p,ci,stats] = ttest(r1,r2,'tail','both');
            diff_h(nt, nt2,b) = h;
            diff_p(nt, nt2,b) = p;
            diff_t(nt, nt2,b) = stats.tstat;
        end
    end
end
% % FDR Correction
% for b = 1:n.bands
%     for nt = 1:n.netws
%         [~, dummy2, ~, dummy4]=fdr_bh(run1_p(nt,b),fdr_thres,'pdep','no');
%         crit_p_run1(nt,b) = dummy2;
%         adj_p_run1(:,:) = dummy4;
%         clear dummy2; clear dummy4;
%         
%         [~, dummy2, ~, dummy4]=fdr_bh(run2_p(nt,b),fdr_thres,'pdep','no');
%         crit_p_run2(nt,b) = dummy2;
%         adj_p_run2(:,:) = dummy4;
%         clear dummy2; clear dummy4;
%     end
% end

for b = 1:n.bands
    figure(); set(gcf,'Color','white'); box OFF;
    imagesc(diff_t(:,:,b)); axis square; %caxis([-mm mm]);
    colorbar; daspect([1 1 1]); reds = 0:0.0005:1;
%     j_colors = [1 1 1; ones(length(reds),1), flipud(reds'), flipud(reds'); 1 0 0];
%     colormap((j_colors))
    
    hold on; for b2=1:n.netws    line([b2+0.5 b2+0.5],[0.5 n.netws+0.5],'color','k');   line([0.5 n.netws+0.5],[b2+0.5 b2+0.5],'color','k');    end
    set(gca,'YTick',[1:n.netws],'YTickLabel',n.netwname,'Fontsize',16);
    set(gca,'XTick',[1:n.netws],'XTickLabel',n.netwname,'Fontsize',16);
    title(['Session Effect ' n.bandname{b} ' Netw n= ' num2str(n.subj) ]);
    
    [ivect_unc,jvect_unc]=find(diff_p(:,:,b)<=0.05);
    scatter(jvect_unc,ivect_unc,36,'ok');
    clear ivect_unc; clear jvect_unc;
    
%     [ivect,jvect]=find(p_corr<=0.05);
%     scatter(jvect,ivect,9,'ok','filled');
%     clear ivect; clear jvect;
    
    
    saveas(gcf, ['session_effect_' n.bandname{b} '_netw.png']);
end
close all;