function handle = plotTaskShaping(data, protocol)

if strcmp(protocol,'PoissonBlocksCondensed3m') 
    firstTestMaze = 10;
    lastShapeMaze = 9;
    aoeFlag       = 1;
elseif strcmp(protocol,'PoissonBlocksPermanentCues')
    firstTestMaze = 7;
    lastShapeMaze = 6;
    aoeFlag       = 0;
else
end

numMice = size(data.daysToMaze,1);

figure
subplot(1,4,1); hold on

% days to reach testing maze
for jj =1:size(data.currMaze,1)
    mouseLearn  = [];
    mouseLearn  = data.currMaze(jj,:);
    if all(isnan(mouseLearn)) || max(mouseLearn) < firstTestMaze
        continue
    else
        removeSess     = [];
        removeSess     = find(mouseLearn==firstTestMaze,1,'first');   
        maxSess(jj)    = removeSess+1;
        mouseLearn(removeSess+1:end) = [];
        plot(mouseLearn,'-','color',[.6 .6 .6],'linewidth',.5); hold on
    end
end

medianData = round(nanmedian(data.currMaze(:,1:max(maxSess))));
if aoeFlag == 1
    medianData(medianData>10) = 10;
else
    medianData(medianData>7) = 7;
end
plot(medianData,'k-','linewidth',2)

if aoeFlag == 1
    ylim([0.5 10.5])
    xlim([0 numel(medianData)])
else
    ylim([0.5 7.5])
    xlim([0 numel(medianData)])
end
box off; 
set(gca,'fontsize',10,'ytick',1:max(data.maxMazeID_xMouse))
set(gca,'TickDir','out');
xlabel('session #','fontsize',12); 
ylabel('maze ID','fontsize',12)
title('days to reach testing maze','fontsize',12)
set(gca,'fontsize',10,'xtick',0:10:numel(medianData))


% days to reach successive maze
subplot(1,4,2)
hold on
for jj=1:(lastShapeMaze+1)
    days          = [];
    days          = data.daysToMaze(:,jj); 
    days(days==0) = NaN;
    
    toPlot          = []; toPlot      = days;
    xIdx          = []; xIdx      = ones(size(days,1),1).*jj;
    catMarker     = []; catMarker = {'x'};
    colorIdx      = []; colorIdx  = {[0.8 0.8 0.8]};
    
    plotSpread(toPlot,'categoryIdx',xIdx,'distributionIdx',xIdx,'categoryColors',colorIdx,'categoryMarkers','x','binWidth',0.0001);   
end
errorbar(1:(lastShapeMaze+1),nanmean(data.daysToMaze(:,1:(lastShapeMaze+1))),...
    nanstd(data.daysToMaze(:,1:(lastShapeMaze+1)),0,1)./sqrt(numMice-1),...
    'ko-','linewidth',3,'markerfacecolor','k','MarkerSize',.1)
box off; 
set(gca,'fontsize',10)%,'xtick',[1 6 11])%,'xticklabel',[mazeLbl{1} mazeLbl{6} mazeLbl{11}])
xlabel('maze','fontsize',12); ylabel('# sessions','fontsize',12)
title('days to reach successive maze','fontsize',12)
if aoeFlag == 1
    xticks([1:9])
    xticklabels({'1','2','3','4','5','6','7','8','9'}) % 
    xlim([0.5 9.5])
else
    xticks([1:6])
    xticklabels({'1','2','3','4','5','6'}) % ,'7','8','9'
    xlim([0.5 6.5])
end
set(gca,'TickDir','out'); 


% days spent per maze
subplot(1,4,3)
hold on
if aoeFlag == 1
    lastShapeMaze = 9;
else
    lastShapeMaze = 6;
end

for jj=1:lastShapeMaze
    days          = [];
    days          = data.daysPerMaze(:,jj); 
    days(days==0) = NaN;
    
    toPlot        = []; toPlot    = days;
    xIdx          = []; xIdx      = ones(size(days,1),1).*jj;
    catMarker     = []; catMarker = {'x'};
    colorIdx      = []; colorIdx  = {[0.9 0.9 0.9]};
   
    plotSpread(toPlot,'categoryIdx',xIdx,'distributionIdx',xIdx,'categoryColors',colorIdx,'categoryMarkers','x','binWidth',0.0001);
end

errorbar(1:lastShapeMaze,nanmean(data.daysPerMaze(:,1:lastShapeMaze)),...
    nanstd(data.daysPerMaze(:,1:lastShapeMaze),0,1)./sqrt(numMice-1),...
    'k.-','linewidth',3,'markerfacecolor','k','MarkerSize',.1)
box off; 
set(gca,'fontsize',10)%,'xtick',[1 6 11])%,'xticklabel',[mazeLbl{1} mazeLbl{6} mazeLbl{11}])
xlabel('maze','fontsize',12); ylabel('# sessions','fontsize',12)
title('days spent per maze','fontsize',12)
if aoeFlag == 1
    xticks([1:9])
    xticklabels({'1','2','3','4','5','6','7','8','9'}) % 
    xlim([0.5 9.5])
else
    xticks([1:6])
    xticklabels({'1','2','3','4','5','6'}) % ,'7','8','9'
    xlim([0.5 6.5])
end
set(gca,'TickDir','out'); 




% mean perf/maze
subplot(1,4,4)
hold on
for jj= 1:lastShapeMaze
    
    days          = [];
    days          = data.perfByMaze(:,jj); 
    days(days==0) = NaN;
    
    toPlot        = []; toPlot    = days;
    xIdx          = []; xIdx      = ones(size(days,1),1).*jj;
    catMarker     = []; catMarker = {'x'};
    colorIdx      = []; colorIdx  = {[0.8 0.8 0.8]};
    
    plotSpread(toPlot,'categoryIdx',xIdx,'distributionIdx',xIdx,'categoryColors',colorIdx,'categoryMarkers','x','binWidth',0.0001);   

end
errorbar(1:lastShapeMaze,nanmean(data.perfByMaze(:,1:lastShapeMaze)),...
    nanstd(data.perfByMaze(:,1:lastShapeMaze),0,1)./sqrt(numMice-1),...
    'k.-','linewidth',3,'markerfacecolor','k','MarkerSize',.1)

box off; set(gca,'fontsize',10)%,'xtick',[1 6 11])%,'xticklabel',[mazeLbl{1} mazeLbl{6} mazeLbl{11}])?
xlabel('maze','fontsize',12); ylabel('% correct','fontsize',12)
ylim([0 100])
title('perf per maze','fontsize',12)
if aoeFlag == 1
    xticks([1:9])
    xticklabels({'1','2','3','4','5','6','7','8','9'}) % 
    xlim([0.5 9.5])
else
    xticks([1:6])
    xticklabels({'1','2','3','4','5','6'}) % ,'7','8','9'
    xlim([0.5 6.5])
end
set(gca,'TickDir','out'); 