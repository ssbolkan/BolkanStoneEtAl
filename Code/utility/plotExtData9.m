function plotExtData9(distDataOFF,distDataON,etDataOFF,etDataON,...
                        groups2plot,ylimits1,ylimits2,posthocStats)
                    
if nargin < 8 || isempty(posthocStats)
    posthocStats = 0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%
% update this section for plotting of measured vars above
% panel 9d left - distance
dataArray   = []; data = []; xIdx = []; 
% groups2plot = [1 2 3]; % indirect aoe,nd,pc tasks
dataArray   = cellfun(@minus,distDataON,distDataOFF,'UniformOutput',false);
data        = [dataArray{groups2plot(1)}'; ...
               dataArray{groups2plot(2)}'; ...
               dataArray{groups2plot(3)}'];
xIdx        = [ones(size(dataArray{groups2plot(1)}))';   ...
               ones(size(dataArray{groups2plot(2)}))'*2; ...
               ones(size(dataArray{groups2plot(3)}))'*3];
yaxName     = 'delta distance (cm, on-off)';
%%%%%%%%%%%%%%%%%%%%%%%%%

groupNames  = {'aoe';'ctrl#1';'ctrl#2'};
catMarker   = {'x','x', 'x'};
colorIdx    = {[0.8 0.8 0.8],[0.99 0.7 0.99], [0.7 0.99 0.99]};
xlimits     = [0.5 3.5];

figure
subplot(1,2,1)
hold on
plotSpread(data,'categoryIdx',xIdx,'distributionIdx',xIdx,...
    'categoryColors',colorIdx,'categoryMarkers',catMarker);

nGroup = 1;
h = errorbar(nGroup,nanmean(dataArray{groups2plot(nGroup)}),nanstd(dataArray{groups2plot(nGroup)})./ ...
    sqrt(numel(dataArray{groups2plot(nGroup)})-1), ...
    '.-','linewidth',4,'color','k','markersize',0.1); h.CapSize = 16;

nGroup = 2; 
h = errorbar(nGroup,nanmean(dataArray{groups2plot(nGroup)}),nanstd(dataArray{groups2plot(nGroup)})./ ...
    sqrt(numel(dataArray{groups2plot(nGroup)})-1), ...
    '.-','linewidth',4,'color','m','markersize',0.1); h.CapSize = 16;

nGroup = 3;
h = errorbar(nGroup,nanmean(dataArray{groups2plot(nGroup)}),nanstd(dataArray{groups2plot(nGroup)})./ ...
    sqrt(numel(dataArray{groups2plot(nGroup)})-1), ...
    '.-','linewidth',4,'color','c','markersize',0.1); h.CapSize = 16;

set(gca,'xtick',1:3,'xticklabel',groupNames)
rotateXLabels(gca,45)
ylabel(yaxName,'fontsize',12)
xlim(xlimits)
ylim(ylimits1)
plot(xlimits(1):.5:xlimits(2),zeros(1,numel(xlimits(1):.5:xlimits(2))),'k-');
box off
set(gca,'TickDir','out');

xIdx = categorical(xIdx);
[p,t,stats] = anova1(data,xIdx,'off');
text(0.6,ylimits1(1)*.85,['group p= ' num2str(p)])
text(0.6,ylimits1(1)*.9,['df= ' num2str(stats.df)])
text(0.6,ylimits1(1)*.95,['F= ' num2str((t{2,5}))])      

if posthocStats == 1
     
    [g1Vg2_p, ~ , g1Vg2_stats] = ranksum(dataArray{groups2plot(1)}, ...
        dataArray{groups2plot(2)},'method','approximate');    
    plot(1:2,ones(1,2)*ylimits1(2)*.8,'k-')
    text(0.5,ylimits1(2)*.85,['p= ' num2str(g1Vg2_p)],'FontSize',6)
    text(0.5,ylimits1(2)*.8,['z= ' num2str(g1Vg2_stats.zval)],'FontSize',6)
    text(0.5,ylimits1(2)*.75,['df= ' num2str(numel(dataArray{groups2plot(1)})+ ...
        numel(dataArray{groups2plot(2)})-1)],'FontSize',6)

    [g1Vg3_p, ~ , g1Vg3_stats] = ranksum(dataArray{groups2plot(1)}, ...
        dataArray{groups2plot(3)},'method','approximate');
    plot(1:3,ones(1,3)*ylimits1(2)*.9,'k-')
    text(1.5,ylimits1(2)*.95,['p= ' num2str(g1Vg3_p)],'FontSize',6)
    text(1.5,ylimits1(2)*.9,['z= ' num2str(g1Vg3_stats.zval)],'FontSize',6)
    text(1.5,ylimits1(2)*.85,['df= ' num2str(numel(dataArray{groups2plot(1)})+ ...
        numel(dataArray{groups2plot(3)})-1)],'FontSize',6)
    
   
    [g2Vg3_p, ~ , g2Vg3_stats] = ranksum(dataArray{groups2plot(2)}, ...
        dataArray{groups2plot(3)},'method','approximate');
    plot(2:3,ones(1,2)*ylimits1(2)*.7,'k-')
    text(2.5,ylimits1(2)*.75,['p= ' num2str(g2Vg3_p)],'FontSize',6)
    text(2.5,ylimits1(2)*.7,['z= ' num2str(g2Vg3_stats.zval)],'FontSize',6)
    text(2.5,ylimits1(2)*.65,['df= ' num2str(numel(dataArray{groups2plot(2)})+ ...
        numel(dataArray{groups2plot(3)})-2)],'FontSize',6)
else
end


%%%%%%%%%%%%%%%%%%%%%%%%%%
% update this section for plotting of measured vars above
% panel 9d right - excess travel
dataArray   = []; data = []; xIdx = []; 
% groups2plot = [1 2 3]; % indirect aoe,nd,pc tasks
dataArray   = cellfun(@minus,etDataON,etDataOFF,'UniformOutput',false);
data        = [dataArray{groups2plot(1)}'; ...
               dataArray{groups2plot(2)}'; ...
               dataArray{groups2plot(3)}'];
xIdx        = [ones(size(dataArray{groups2plot(1)}))';   ...
               ones(size(dataArray{groups2plot(2)}))'*2; ...
               ones(size(dataArray{groups2plot(3)}))'*3];
yaxName     = 'delta trials w/ excess travel (%, on-off)';
%%%%%%%%%%%%%%%%%%%%%%%%%

groupNames  = {'aoe';'ctrl#1';'ctrl#2'};
catMarker   = {'x','x', 'x'};
colorIdx    = {[0.8 0.8 0.8],[0.99 0.7 0.99], [0.7 0.99 0.99]};
xlimits     = [0.5 3.5];

subplot(1,2,2)
hold on
plotSpread(data,'categoryIdx',xIdx,'distributionIdx',xIdx,...
    'categoryColors',colorIdx,'categoryMarkers',catMarker);

nGroup = 1;
h = errorbar(nGroup,nanmean(dataArray{groups2plot(nGroup)}),nanstd(dataArray{groups2plot(nGroup)})./ ...
    sqrt(numel(dataArray{groups2plot(nGroup)})-1), ...
    '.-','linewidth',4,'color','k','markersize',0.1); h.CapSize = 16;

nGroup = 2; 
h = errorbar(nGroup,nanmean(dataArray{groups2plot(nGroup)}),nanstd(dataArray{groups2plot(nGroup)})./ ...
    sqrt(numel(dataArray{groups2plot(nGroup)})-1), ...
    '.-','linewidth',4,'color','m','markersize',0.1); h.CapSize = 16;

nGroup = 3;
h = errorbar(nGroup,nanmean(dataArray{groups2plot(nGroup)}),nanstd(dataArray{groups2plot(nGroup)})./ ...
    sqrt(numel(dataArray{groups2plot(nGroup)})-1), ...
    '.-','linewidth',4,'color','c','markersize',0.1); h.CapSize = 16;

set(gca,'xtick',1:3,'xticklabel',groupNames)
rotateXLabels(gca,45)
ylabel(yaxName,'fontsize',12)
xlim(xlimits)
ylim(ylimits2)
plot(xlimits(1):.5:xlimits(2),zeros(1,numel(xlimits(1):.5:xlimits(2))),'k-');
box off
set(gca,'TickDir','out');

xIdx = categorical(xIdx);
[p,t,stats] = anova1(data,xIdx,'off');
text(0.6,ylimits2(1)*.85,['group p= ' num2str(p)])
text(0.6,ylimits2(1)*.9,['df= ' num2str(stats.df)])
text(0.6,ylimits2(1)*.95,['F= ' num2str((t{2,5}))])  

if posthocStats == 1
     
    [g1Vg2_p, ~ , g1Vg2_stats] = ranksum(dataArray{groups2plot(1)}, ...
        dataArray{groups2plot(2)},'method','approximate');    
    plot(1:2,ones(1,2)*ylimits2(2)*.8,'k-')
    text(0.5,ylimits2(2)*.85,['p= ' num2str(g1Vg2_p)],'FontSize',6)
    text(0.5,ylimits2(2)*.8,['z= ' num2str(g1Vg2_stats.zval)],'FontSize',6)
    text(0.5,ylimits2(2)*.75,['df= ' num2str(numel(dataArray{groups2plot(1)})+ ...
        numel(dataArray{groups2plot(2)})-1)],'FontSize',6)

    [g1Vg3_p, ~ , g1Vg3_stats] = ranksum(dataArray{groups2plot(1)}, ...
        dataArray{groups2plot(3)},'method','approximate');
    plot(1:3,ones(1,3)*ylimits2(2)*.9,'k-')
    text(1.5,ylimits2(2)*.95,['p= ' num2str(g1Vg3_p)],'FontSize',6)
    text(1.5,ylimits2(2)*.9,['z= ' num2str(g1Vg3_stats.zval)],'FontSize',6)
    text(1.5,ylimits2(2)*.85,['df= ' num2str(numel(dataArray{groups2plot(1)})+ ...
        numel(dataArray{groups2plot(3)})-1)],'FontSize',6)
    
   
    [g2Vg3_p, ~ , g2Vg3_stats] = ranksum(dataArray{groups2plot(2)}, ...
        dataArray{groups2plot(3)},'method','approximate');
    plot(2:3,ones(1,2)*ylimits2(2)*.7,'k-')
    text(2.5,ylimits2(2)*.75,['p= ' num2str(g2Vg3_p)],'FontSize',6)
    text(2.5,ylimits2(2)*.7,['z= ' num2str(g2Vg3_stats.zval)],'FontSize',6)
    text(2.5,ylimits2(2)*.65,['df= ' num2str(numel(dataArray{groups2plot(2)})+ ...
        numel(dataArray{groups2plot(3)})-2)],'FontSize',6)
else
end



end
