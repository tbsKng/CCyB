oldfolder = cd('Figures\Chart_Policy');

style{1}='-';style{2}='-.'; style{3}='--';style{4}=':';style{5}='-+';style{6}='-o';style{7}='--o';
col{1}=[0 0 0];col{2}=[1 0 0];col{3}=[0 0 1];col{4}=[0.4 0.4 0.4];col{5}='g';col{6}=[1 0 1];col{7}='m';


%% Risk Shock to all Sectors

h1=figure;pause(2);frame_h = get(handle(gcf),'JavaFrame');set(frame_h,'Maximized',1);
%hh=suptitle({['Bank Risk Shock under calibrated and desired TCR policy ' country_names(country_nr,:)];''})
hh=suptitle({'10% Risk Shock to All Sectors'});
%hh=suptitle({num2str(risk_shock_mag) '% Risk Shock to NFCs'});
P = get(hh,'Position');
set(hh,'Position',[P(1) P(2)+0.02 P(3)]);
set(hh,'FontSize',16);
for zr = 1:size(data_Se_surprise{1},2)
    sh=subplot(4,5,zr);
    
    
    %plot(1:1:16,data_Se_CCB_surprise{2}(1:16,zr),style{1}, 'color',col{1},'LineWidth', 2,'MarkerSize',2);hold all;
    plot(1:1:16,data_Se_CCB_surprise{1}(1:16,zr),style{3}, 'color',col{3},'LineWidth', 2,'MarkerSize',2);hold all;
    %        plot(1:1:16,data_Se_CCB_surprise{3}(1:16,zr),style{4}, 'color',col{4},'LineWidth', 2,'MarkerSize',2);hold all;
    %        plot(1:1:16,data_Se_CCB_surprise{4}(1:16,zr),style{5}, 'color',col{5},'LineWidth', 2,'MarkerSize',2);hold all;
    ph=plot(1:1:16,data_Se_surprise{1}(1:16,zr),style{2}, 'color',col{2},'LineWidth', 1.5,'MarkerSize',2);hold all;
    % plot(1:1:16,data_Se_CCB_surprise{6}(1:16,zr),style{1}, 'color',col{1},'LineWidth', 2,'MarkerSize',2);hold all;
    plot(1:1:16,data_Se_CCB_surprise{6}(1:16,zr),style{4}, 'color',col{4},'LineWidth', 2,'MarkerSize',2);hold all;
    
    line([1 16],[0 0],'color',[0 0 0]); %axis tight;
    
    if zr==17
        l={'Sectoral CCyB w/ Corp. CR to Defaults','Benchmark','Sectoral both loans.'};
        lh=legend(l,'Orientation','horizontal');
        sp=get(sh,'position');
        set(lh,'position',[sp(1),sp(2)-.12,sp(3),.1]);
    end
    
    title([varnames9(zr,:)],'FontSize',12);
    
    
    set(gca,'XTick',[4 8 12 16],'FontSize',13)
    xlim([1 16]);
    
end


saveas(h1,'News_SeQ8_surprise_Policy_3DM','fig')
saveas(h1,'News_SeQ8_surprise_Policy_3DM','png')
%saveas(h1,'News_SeQ8_surprise_Policy_3DM','pdf')


%% HH Risk Shock


% h1=figure;pause(2);frame_h = get(handle(gcf),'JavaFrame');set(frame_h,'Maximized',1);
% %hh=suptitle({['Bank Risk Shock under calibrated and desired TCR policy ' country_names(country_nr,:)];''})
% hh=suptitle({'10% Mortgage Risk Shock'});
% P = get(hh,'Position');
% set(hh,'Position',[P(1) P(2)+0.02 P(3)]);
% set(hh,'FontSize',16);
% for zr = 1:size(data_Sm_surprise{1},2)
%     sh=subplot(4,4,zr);
%
%         plot(1:1:16,data_Sm_CCB_surprise{2}(1:16,zr),style{1}, 'color',col{1},'LineWidth', 2,'MarkerSize',2);hold all;
%         plot(1:1:16,data_Sm_CCB_surprise{1}(1:16,zr),style{3}, 'color',col{3},'LineWidth', 2,'MarkerSize',2);hold all;
%         plot(1:1:16,data_Sm_CCB_surprise{3}(1:16,zr),style{4}, 'color',col{4},'LineWidth', 2,'MarkerSize',2);hold all;
%         plot(1:1:16,data_Sm_CCB_surprise{4}(1:16,zr),style{5}, 'color',col{5},'LineWidth', 2,'MarkerSize',2);hold all;
%         plot(1:1:16,data_Sm_surprise{1}(1:16,zr),style{2}, 'color',col{2},'LineWidth', 1.5,'MarkerSize',2);hold all;
%
%
%     line([1 16],[0 0],'color',[0 0 0]); %axis tight;
%     if zr==15
%         l={'CCyB Total Loans','Sectoral CCyB','TR Strict Inflation Stabilization','TR Response to Loans','Benchmark'};
%         lh=legend(l,'Orientation','horizontal');
%         sp=get(sh,'position');
%         set(lh,'position',[sp(1),sp(2)-.12,sp(3),.1]);
%     end
%
%     title([varnames9(zr,:)],'FontSize',12);
%
%
%     set(gca,'XTick',[4 8 12 16],'FontSize',13)
%     xlim([1 16]);
%
% end
%
%
% saveas(h1,'News_SmQ8_surprise_Policy_3DM','fig')
% saveas(h1,'News_SmQ8_surprise_Policy_3DM','png')
% %saveas(h1,'News_SmQ8_surprise_Policy_3DM','pdf')


%% Bank Risk Shock

% h1=figure;pause(2);frame_h = get(handle(gcf),'JavaFrame');set(frame_h,'Maximized',1);
% %hh=suptitle({['Bank Risk Shock under calibrated and desired TCR policy ' country_names(country_nr,:)];''})
% hh=suptitle({'10% Bank Risk Shock'});
% P = get(hh,'Position');
% set(hh,'Position',[P(1) P(2)+0.02 P(3)]);
% set(hh,'FontSize',16);
% for zr = 1:size(data_SF_surprise{1},2)
%     sh=subplot(4,4,zr);
%
%
%         plot(1:1:16,data_SF_CCB_surprise{2}(1:16,zr),style{1}, 'color',col{1},'LineWidth', 2,'MarkerSize',2);hold all;
%         plot(1:1:16,data_SF_CCB_surprise{1}(1:16,zr),style{3}, 'color',col{3},'LineWidth', 2,'MarkerSize',2);hold all;
%  %       plot(1:1:16,data_SF_CCB_surprise{3}(1:16,zr),style{4}, 'color',col{4},'LineWidth', 2,'MarkerSize',2);hold all;
%  %       plot(1:1:16,data_SF_CCB_surprise{4}(1:16,zr),style{5}, 'color',col{5},'LineWidth', 2,'MarkerSize',2);hold all;
%         plot(1:1:16,data_SF_surprise{1}(1:16,zr),style{2}, 'color',col{2},'LineWidth', 1.5,'MarkerSize',2);hold all;
%
%     line([1 16],[0 0],'color',[0 0 0]); %axis tight;
%
%     if zr==15
%         l={'CCyB Total Loans','Sectoral CCyB','Benchmark'};
%         lh=legend(l,'Orientation','horizontal');
%         sp=get(sh,'position');
%         set(lh,'position',[sp(1),sp(2)-.12,sp(3),.1]);
%     end
%
%     title([varnames9(zr,:)],'FontSize',12);
%
%
%     set(gca,'XTick',[4 8 12 16],'FontSize',13)
%     xlim([1 16]);
%
% end
%
%
% saveas(h1,'News_SFQ8_surprise_Policy_3DM','fig')
% saveas(h1,'News_SFQ8_surprise_Policy_3DM','png')
% %saveas(h1,'News_SFQ8_surprise_Policy_3DM','pdf')


%% HH Risk Shock and Mortgage Bank Risk

h1=figure;pause(2);frame_h = get(handle(gcf),'JavaFrame');set(frame_h,'Maximized',1);
%hh=suptitle({['Bank Risk Shock under calibrated and desired TCR policy ' country_names(country_nr,:)];''})
hh=suptitle({'10% HH Risk and Mortgage Bank Risk Shock'});
P = get(hh,'Position');
set(hh,'Position',[P(1) P(2)+0.02 P(3)]);
set(hh,'FontSize',16);
for zr = 1:size(data_SmSH_surprise{1},2)
    sh=subplot(4,5,zr);
    
    
    % plot(1:1:16,data_SmSH_CCB_surprise{2}(1:16,zr),style{1}, 'color',col{1},'LineWidth', 2,'MarkerSize',2);hold all;
    plot(1:1:16,data_SmSH_CCB_surprise{1}(1:16,zr),style{3}, 'color',col{3},'LineWidth', 2,'MarkerSize',2);hold all;
    %      plot(1:1:16,data_SmSH_CCB_surprise{3}(1:16,zr),style{4}, 'color',col{4},'LineWidth', 2,'MarkerSize',2);hold all;
    %      plot(1:1:16,data_SmSH_CCB_surprise{4}(1:16,zr),style{5}, 'color',col{5},'LineWidth', 2,'MarkerSize',2);hold all;
    plot(1:1:16,data_SmSH_surprise{1}(1:16,zr),style{2}, 'color',col{2},'LineWidth', 1.5,'MarkerSize',2);hold all;
    plot(1:1:16,data_SmSH_CCB_surprise{5}(1:16,zr),style{1}, 'color',col{1},'LineWidth', 2,'MarkerSize',2);hold all;
    plot(1:1:16,data_SmSH_CCB_surprise{6}(1:16,zr),style{4}, 'color',col{4},'LineWidth', 2,'MarkerSize',2);hold all;
    line([1 16],[0 0],'color',[0 0 0]); %axis tight;
    
    if zr==17
        l={'Sectoral CCyB w/ NFC CR risk weight adj.','Benchmark','Corp. Sectoral CCyB only','Sectoral both loans.'};
        lh=legend(l,'Orientation','horizontal');
        sp=get(sh,'position');
        set(lh,'position',[sp(1),sp(2)-.12,sp(3),.1]);
    end
    
    title([varnames9(zr,:)],'FontSize',12);
    
    
    set(gca,'XTick',[4 8 12 16],'FontSize',13)
    xlim([1 16]);
    
end


saveas(h1,'News_SmSHQ8_surprise_Policy_3DM','fig')
saveas(h1,'News_SmSHQ8_surprise_Policy_3DM','png')
%saveas(h1,'News_SmSHQ8_surprise_Policy_3DM','pdf')

%% NFC Risk Shock and NFC Bank Risk Shock

h1=figure;pause(2);frame_h = get(handle(gcf),'JavaFrame');set(frame_h,'Maximized',1);
%hh=suptitle({['Bank Risk Shock under calibrated and desired TCR policy ' country_names(country_nr,:)];''})
hh=suptitle({'10% NFC Risk and NFC Bank Risk Shock'});
P = get(hh,'Position');
set(hh,'Position',[P(1) P(2)+0.02 P(3)]);
set(hh,'FontSize',16);
for zr = 1:size(data_SeSF_surprise{1},2)
    sh=subplot(4,5,zr);
    
    
    %  plot(1:1:16,data_SeSF_CCB_surprise{2}(1:16,zr),style{1}, 'color',col{1},'LineWidth', 2,'MarkerSize',2);hold all;
    plot(1:1:16,data_SeSF_CCB_surprise{1}(1:16,zr),style{3}, 'color',col{3},'LineWidth', 2,'MarkerSize',2);hold all;
    %     plot(1:1:16,data_SeSF_CCB_surprise{3}(1:16,zr),style{4}, 'color',col{4},'LineWidth', 2,'MarkerSize',2);hold all;
    %    plot(1:1:16,data_SeSF_CCB_surprise{4}(1:16,zr),style{5}, 'color',col{5},'LineWidth', 2,'MarkerSize',2);hold all;
    plot(1:1:16,data_SeSF_surprise{1}(1:16,zr),style{2}, 'color',col{2},'LineWidth', 1.5,'MarkerSize',2);hold all;
    %  plot(1:1:16,data_SeSF_surprise{5}(1:16,zr),style{2}, 'color',col{2},'LineWidth', 1.5,'MarkerSize',2);hold all;
    plot(1:1:16,data_SeSF_CCB_surprise{5}(1:16,zr),style{1}, 'color',col{1},'LineWidth', 1.5,'MarkerSize',2);hold all;
    plot(1:1:16,data_SeSF_CCB_surprise{6}(1:16,zr),style{4}, 'color',col{4},'LineWidth', 1.5,'MarkerSize',2);hold all;
    
    
    line([1 16],[0 0],'color',[0 0 0]); %axis tight;
    
    if zr==17
        l={'Sectoral CCyB w/ mort. CR risk weight adj.','Benchmark','Corp. Sectoral CCyB only','Sectoral both loans.'};
        lh=legend(l,'Orientation','horizontal');
        sp=get(sh,'position');
        set(lh,'position',[sp(1),sp(2)-.12,sp(3),.1]);
    end
    
    title([varnames9(zr,:)],'FontSize',12);
    
    
    set(gca,'XTick',[4 8 12 16],'FontSize',13)
    xlim([1 16]);
    
end


saveas(h1,'News_SeSFQ8_surprise_Policy_3DM_1','fig')
saveas(h1,'News_SeSFQ8_surprise_Policy_3DM_1','png')
%saveas(h1,'News_SeSFQ8_surprise_Policy_3DM_1','pdf')

cd(oldfolder)