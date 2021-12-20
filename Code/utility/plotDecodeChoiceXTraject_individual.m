function handle = plotDecodeChoiceXTraject_individual(predChoice)


counter = 1;
xIdx    = [];
data    = []; 
data    = predChoice.trainTOW.xMouseAccMG;
data    = data*100;
predChoiceBins = predChoice.binIDs+11.5;
for nBin = 1:numel(predChoiceBins)
        xIdx(counter:counter+size(data,2)-1,1) = ones(size(data,2),1)*(predChoiceBins(nBin)-2.5);
        xValues(counter:counter+size(data,2)-1,1) = ones(size(data,2),1)*(predChoiceBins(nBin)+7.5);
        counter = counter+ size(data,2);   
end
indivColor  = {[0.8 0.3 0.8]};    
groupColor  = 'm'; 
data        = data';
data        = data(:);
handle      = subplot(1,3,1); hold on
colorIdx    = repmat(indivColor,1, numel(predChoiceBins));
xNames      = {[25:25:300]};
plotSpread(data,'categoryIdx',xIdx,'distributionIdx',xIdx,'xNames',xNames,'categoryColors',colorIdx,'xyOri','flipped'); hold on

counter = 1;
xIdx    = [];
data    = []; 
data    = predChoice.trainTOW.xMouseAccTOW;
data    = data*100;
for nBin = 1:numel(predChoiceBins)
        xIdx(counter:counter+size(data,2)-1,1) = ones(size(data,2),1)*(predChoiceBins(nBin)-5);
        counter = counter+ size(data,2);   
end
indivColor  = {[0.3 0.3 0.3]};    
groupColor  = 'k'; 
data        = data';
data        = data(:);
handle         = subplot(1,3,1); hold on
colorIdx    = repmat(indivColor,1, numel(predChoiceBins));
plotSpread(data,'categoryIdx',xIdx,'distributionIdx',xIdx,'xNames',xNames,'categoryColors',colorIdx,'xyOri','flipped'); hold on

counter     = 1;
xIdx        = [];
data        = []; 
data        = predChoice.trainTOW.xMouseAccPC;
data        = data*100;
for nBin = 1:numel(predChoiceBins)
        xIdx(counter:counter+size(data,2)-1,1) = ones(size(data,2),1)*(predChoiceBins(nBin));
        counter = counter+ size(data,2);   
end
indivColor  = {[0.3 0.3 0.8]};    
groupColor  = 'c'; 
data        = data';
data        = data(:);
handle         = subplot(1,3,1); hold on
colorIdx    = repmat(indivColor,1, numel(predChoiceBins));
plotSpread(data,'categoryIdx',xIdx,'distributionIdx',xIdx,'xNames',xNames,'categoryColors',colorIdx,'xyOri','flipped'); hold on
xIdx        = xIdx+6.5; 
data        = NaN(numel(data),1);
plotSpread(data,'categoryIdx',xIdx,'distributionIdx',xIdx,'xNames',xNames,'categoryColors',colorIdx,'xyOri','flipped'); hold on
box off
title('train on AoE')
set(gca,'TickDir','out');
xlim([20 100])
xticks([20 40 60 80 100])
ylim([0 300])
ylabel('y position (cm)')

counter     = 1;
xIdx        = [];
test        = []; 
test        = predChoice.trainMG.xMouseAccMG;
test        = test*100;
for nBin = 1:numel(predChoiceBins)
        xIdx(counter:counter+size(test,2)-1,1) = ones(size(test,2),1)*(predChoiceBins(nBin)-2.5);
        xValues(counter:counter+size(test,2)-1,1) = ones(size(test,2),1)*(predChoiceBins(nBin)+7.5);
        counter = counter+ size(test,2);   
end
indivColor  = {[0.8 0.3 0.8]};    groupColor     = 'm'; 
test        = test';
test        = test(:);
handle         = subplot(1,3,2); hold on
colorIdx    = repmat(indivColor,1, numel(predChoiceBins));
xNames      = {[25:25:300]};
plotSpread(test,'categoryIdx',xIdx,'distributionIdx',xIdx,'xNames',xNames,'categoryColors',colorIdx,'xyOri','flipped'); hold on

counter     = 1;
xIdx        = [];
test        = []; 
test        = predChoice.trainMG.xMouseAccTOW;
test        = test*100;

for nBin = 1:numel(predChoiceBins)
        xIdx(counter:counter+size(test,2)-1,1) = ones(size(test,2),1)*(predChoiceBins(nBin)-5);
        counter = counter+ size(test,2);   
end
indivColor  = {[0.3 0.3 0.3]};
groupColor  = 'k'; 
test        = test';
test        = test(:);
handle         = subplot(1,3,2); hold on
colorIdx    = repmat(indivColor,1, numel(predChoiceBins));
plotSpread(test,'categoryIdx',xIdx,'distributionIdx',xIdx,'xNames',xNames,'categoryColors',colorIdx,'xyOri','flipped'); hold on

counter     = 1;
xIdx        = [];
test        = []; 
test        = predChoice.trainMG.xMouseAccPC;
test        = test*100;
for nBin = 1:numel(predChoiceBins)
        xIdx(counter:counter+size(test,2)-1,1) = ones(size(test,2),1)*(predChoiceBins(nBin));
        counter = counter+ size(test,2);   
end
indivColor  = {[0.3 0.3 0.8]};
groupColor  = 'c'; 
test        = test';
test        = test(:);
handle         = subplot(1,3,2); hold on
colorIdx    = repmat(indivColor,1, numel(predChoiceBins));
plotSpread(test,'categoryIdx',xIdx,'distributionIdx',xIdx,'xNames',xNames,'categoryColors',colorIdx,'xyOri','flipped'); hold on
xIdx        = xIdx+6.5; 
test        = NaN(numel(test),1);
plotSpread(test,'categoryIdx',xIdx,'distributionIdx',xIdx,'xNames',xNames,'categoryColors',colorIdx,'xyOri','flipped'); hold on
box off
title('train on ctrl#1')
set(gca,'TickDir','out');
xlim([20 100])
xticks([20 40 60 80 100])
ylim([0 300])
xlabel('frac. decode accuracy')


counter     = 1;
xIdx        = [];
test        = [];
test        = predChoice.trainPC.xMouseAccMG;
test        = test*100;
for nBin = 1:numel(predChoiceBins)
        xIdx(counter:counter+size(test,2)-1,1) = ones(size(test,2),1)*(predChoiceBins(nBin)-2.5);
        xValues(counter:counter+size(test,2)-1,1) = ones(size(test,2),1)*(predChoiceBins(nBin)+7.5);
        counter = counter+ size(test,2);   
end
indivColor  = {[0.8 0.3 0.8]};
groupColor  = 'm'; 
test        = test';
test        = test(:);
handle         = subplot(1,3,3); hold on
colorIdx    = repmat(indivColor,1, numel(predChoiceBins));
xNames      = {[25:25:300]};
plotSpread(test,'categoryIdx',xIdx,'distributionIdx',xIdx,'xNames',xNames,'categoryColors',colorIdx,'xyOri','flipped'); hold on

counter     = 1;
xIdx        = [];
test        = [];
test        = predChoice.trainPC.xMouseAccTOW;
test        = test*100;
for nBin = 1:numel(predChoiceBins)
        xIdx(counter:counter+size(test,2)-1,1) = ones(size(test,2),1)*(predChoiceBins(nBin)-5);
        counter = counter+ size(test,2);   
end
indivColor  = {[0.3 0.3 0.3]};
groupColor  = 'k'; 
test        = test';
test        = test(:);
handle         = subplot(1,3,3); hold on
colorIdx    = repmat(indivColor,1, numel(predChoiceBins));
plotSpread(test,'categoryIdx',xIdx,'distributionIdx',xIdx,'xNames',xNames,'categoryColors',colorIdx,'xyOri','flipped'); hold on

counter     = 1;
xIdx        = [];
test        = [];
test        = predChoice.trainPC.xMouseAccPC;
test        = test*100;
for nBin = 1:numel(predChoiceBins)
        xIdx(counter:counter+size(test,2)-1,1) = ones(size(test,2),1)*(predChoiceBins(nBin));
        counter = counter+ size(test,2);   
end
indivColor  = {[0.3 0.3 0.8]};
groupColor  = 'c'; 
test        = test';
test        = test(:);
handle         = subplot(1,3,3); hold on
colorIdx    = repmat(indivColor,1, numel(predChoiceBins));
plotSpread(test,'categoryIdx',xIdx,'distributionIdx',xIdx,'xNames',xNames,'categoryColors',colorIdx,'xyOri','flipped'); hold on
xIdx        = xIdx+6.5;
test        = NaN(numel(test),1);
plotSpread(test,'categoryIdx',xIdx,'distributionIdx',xIdx,'xNames',xNames,'categoryColors',colorIdx,'xyOri','flipped'); hold on
box off
title('train on ctrl#2')
set(gca,'TickDir','out');
xlim([20 100])
xticks([20 40 60 80 100])
ylim([0 300])
