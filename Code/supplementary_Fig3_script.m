%% Script to generate ExtData Fig 5B-C data plots
% real time place preference with DMS laser
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load conditioned place preference data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
cd(globalParams.dataPath)
load('rtPPdata.mat', 'allDataTestDMS', 'allDataBaselineDMS')

%%
% get relevant data from baseline and test sessions (according to
% compiled ethovision excel files saved as .mat cppData above)
mouseIdx           = allDataTestDMS(:,13);      % mouseID
mouseIdx_BL        = allDataBaselineDMS(:,13);
genotypeIdx        = allDataTestDMS(:,14);      % 1 is d1; 2 is d2/a2a
genotypeIdx_BL     = allDataBaselineDMS(:,14);
virusIdx           = allDataTestDMS(:,15);      % 0 is no opsin; 1 is NpHR
virusIdx_BL        = allDataBaselineDMS(:,15);  
velocityIdx        = allDataTestDMS(:,9);       % ethovision output velocity per sample (in mm/s)
velocityIdx_BL     = allDataBaselineDMS(:,9);
distanceIdx        = allDataTestDMS(:,8);       % ethovision output distance per sample (in mm)
distanceIdx_BL     = allDataBaselineDMS(:,8);
laserSideIdx       = allDataTestDMS(:,16);      % 0 is right zone laserON; 1 is left zone laserON
laserSideIdx_BL    = allDataBaselineDMS(:,16);  
currLocationIdx    = allDataTestDMS(:,10);      % 1 is in right zone; 0 in left zone
currLocationIdx_BL = allDataBaselineDMS(:,10);  
cumulTimeIdx       = allDataTestDMS(:,1);       % ~30 hz ethovision samples (adds up in seconds)
cumulTimeIdx_BL    = allDataBaselineDMS(:,1);  

% analyze delta time and delta velocity in future/current laser side for BL/test (stim-nonStim)
mouseIDs = unique(mouseIdx);
for nMouse = 1:numel(mouseIDs)
    
    % first get mouse genotype and virus and laser stim side
    indexIn = [];
    indexIn = find(mouseIdx== mouseIDs(nMouse),1,'first');   
    mouseGenotype(nMouse) = genotypeIdx(indexIn);
    mouseVirus(nMouse)    = virusIdx(indexIn);
    mouseID(nMouse)       = mouseIdx(indexIn);
    laserSide(nMouse)     = laserSideIdx(indexIn);
    
    % for time and vel on test session
    mouseTime = []; 
    mouseTime = cumulTimeIdx(mouseIdx == mouseIDs(nMouse)); mouseTime = [NaN; diff(mouseTime)];
    mouseLoc  = [];
    mouseLoc  = currLocationIdx(mouseIdx == mouseIDs(nMouse));
    mouseVel  = []; 
    mouseVel  = velocityIdx(mouseIdx == mouseIDs(nMouse));
     
    timeON(nMouse)  = nansum(mouseTime(mouseLoc == laserSide(nMouse)));
    timeOFF(nMouse) = nansum(mouseTime(mouseLoc ~= laserSide(nMouse)));
    velON(nMouse)   = nanmean(mouseVel(mouseLoc == laserSide(nMouse)));
    velOFF(nMouse)  = nanmean(mouseVel(mouseLoc ~= laserSide(nMouse)));
    
    % for time and vel on baseline session
    mouseTime = []; 
    mouseTime = cumulTimeIdx_BL(mouseIdx_BL == mouseIDs(nMouse)); mouseTime = [NaN; diff(mouseTime)];
    mouseLoc  = [];
    mouseLoc  = currLocationIdx_BL(mouseIdx_BL == mouseIDs(nMouse));
    mouseVel  = []; 
    mouseVel  = velocityIdx_BL(mouseIdx_BL == mouseIDs(nMouse));
    
    timeON_BL(nMouse)  = nansum(mouseTime(mouseLoc == laserSide(nMouse)));
    timeOFF_BL(nMouse) = nansum(mouseTime(mouseLoc ~= laserSide(nMouse)));
    velON_BL(nMouse)   = nanmean(mouseVel(mouseLoc == laserSide(nMouse)));
    velOFF_BL(nMouse)  = nanmean(mouseVel(mouseLoc ~= laserSide(nMouse)));
    
    
end

% get on-off delta measurements
diffTimeBL   = timeON_BL - timeOFF_BL;
diffTimeTest = timeON - timeOFF;
diffVelBL    = velON_BL - velOFF_BL;
diffVelTest  = velON - velOFF;

clearvars -except allDataTestDMS allDataBaselineDMS ...
                  timeON_BL timeOFF_BL velON_BL velOFF_BL ...
                  timeON timeOFF velON velOFF...
                  diffTimeBL diffTimeTest diffVelBL diffVelTest ...
                  mouseGenotype mouseVirus mouseID laserSide
              
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot panel 1B: delta time on laser side (on-off)       
clearvars -except allDataTestDMS allDataBaselineDMS ...
                  timeON_BL timeOFF_BL velON_BL velOFF_BL ...
                  timeON timeOFF velON velOFF...
                  diffTimeBL diffTimeTest diffVelBL diffVelTest ...
                  mouseGenotype mouseVirus mouseID laserSide
                                           
%%%%%%%%%%%%%%%%%%% update this section for plotting of measured vars above
%%%%%%%%%%%%%%%%%%% This plots Baseline delta Time on laser side (on-off)
% remove mice with baseline bias greater than 60s
removeBLbias = 1;
if removeBLbias == 1
    removeBLInd  = abs(diffTimeBL)<60;
else
    removeBLInd  = ones(1,numel(diffTimeBL));
end
% mouseGenotype: 2 indirect (a2a or d2r-cre); 1 direct (d1r-cre)
% mouseVirus: 1 (NpHR); 2 (no opsin)
dataArray        = []; data = []; xIdx = []; 
dataArray{1,1}   = diffTimeBL(mouseGenotype==2 & mouseVirus ==1 & removeBLInd==1);
dataArray{1,2}   = diffTimeBL(mouseGenotype==1 & mouseVirus ==1 & removeBLInd==1);
dataArray{1,3}   = diffTimeBL(mouseVirus ==0 & removeBLInd==1);
data        = [diffTimeBL(mouseGenotype==2 & mouseVirus ==1 & removeBLInd==1)'; ...
               diffTimeBL(mouseGenotype==1 & mouseVirus ==1 & removeBLInd==1)'; ...
               diffTimeBL(mouseVirus ==0 & removeBLInd==1)'];
xIdx        = [ones(size(dataArray{1,1}))';   ...
               ones(size(dataArray{1,2}))'*2; ...
               ones(size(dataArray{1,3}))'*3];
ylimits     = [-500 500];
yaxName     = 'delta x-position (cm, on-off)';
titlestr    = 'Baseline';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

groupNames  = {'indirect';'direct';'no opsin'};
catMarker   = {'x','x', 'x'};
colorIdx    = {[0.8 0.8 0.8],[0.8 0.8 0.8], [0.8 0.8 0.8]};
xlimits     = [0.5 3.5];

figure
subplot(1,2,1)
hold on
plotSpread(data,'categoryIdx',xIdx,'distributionIdx',xIdx,...
    'categoryColors',colorIdx,'categoryMarkers',catMarker);
nGroup = 1;
h = errorbar(nGroup,nanmean(dataArray{nGroup}),nanstd(dataArray{nGroup})./sqrt(numel(dataArray{nGroup})-1), ...
    '.-','linewidth',4,'color','k','markersize',0.1); h.CapSize = 16;
nGroup = 2; 
h = errorbar(nGroup,nanmean(dataArray{nGroup}),nanstd(dataArray{nGroup})./sqrt(numel(dataArray{nGroup})-1), ...
    '.-','linewidth',4,'color','k','markersize',0.1); h.CapSize = 16;
nGroup = 3;
h = errorbar(nGroup,nanmean(dataArray{nGroup}),nanstd(dataArray{nGroup})./sqrt(numel(dataArray{nGroup})-1), ...
    '.-','linewidth',4,'color','k','markersize',0.1); h.CapSize = 16;

set(gca,'xtick',1:3,'xticklabel',groupNames)
rotateXLabels(gca,45)
ylabel(yaxName,'fontsize',12)
xlim(xlimits)
ylim(ylimits)
plot(xlimits(1):.5:xlimits(2),zeros(1,numel(xlimits(1):.5:xlimits(2))),'k-');
box off
set(gca,'TickDir','out');
title(titlestr)

xIdx = categorical(xIdx);
[p,t,stats] = anova1(data,xIdx,'off');
text(0.6,ylimits(2)*.9,['group p= ' num2str(p)])
text(0.6,ylimits(2)*.8,['df= ' num2str(stats.df)])
text(0.6,ylimits(2)*.7,['F= ' num2str((t{2,5}))])

clearvars -except allDataTestDMS allDataBaselineDMS ...
                  timeON_BL timeOFF_BL velON_BL velOFF_BL ...
                  timeON timeOFF velON velOFF...
                  diffTimeBL diffTimeTest diffVelBL diffVelTest ...
                  mouseGenotype mouseVirus mouseID laserSide removeBLInd
                                           
%%%%%%%%%%%%%%%%%% update this section for plotting of measured vars above
%%%%%%%%%%%%%%%%%% This plots Test delta Time on laser side (on-off)
% remove mice with baseline bias greater than 60s if removeBLbias==1 above
% mouseGenotype: 2 indirect (a2a or d2r-cre); 1 direct (d1r-cre)
% mouseVirus: 1 (NpHR); 2 (no opsin)
dataArray        = []; data = []; xIdx = []; 
dataArray{1,1}   = diffTimeTest(mouseGenotype==2 & mouseVirus ==1 & removeBLInd ==1);
dataArray{1,2}   = diffTimeTest(mouseGenotype==1 & mouseVirus ==1 & removeBLInd ==1);
dataArray{1,3}   = diffTimeTest(mouseVirus ==0 & removeBLInd ==1);
data        = [diffTimeTest(mouseGenotype==2 & mouseVirus ==1 & removeBLInd ==1)'; ...
               diffTimeTest(mouseGenotype==1 & mouseVirus ==1 & removeBLInd ==1)'; ...
               diffTimeTest(mouseVirus ==0 & removeBLInd ==1)'];
xIdx        = [ones(size(dataArray{1,1}))';   ...
               ones(size(dataArray{1,2}))'*2; ...
               ones(size(dataArray{1,3}))'*3];
ylimits     = [-500 500];
yaxName     = 'delta x-position (cm, on-off)';
titlestr    = 'Test';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

groupNames  = {'indirect';'direct';'no opsin'};
catMarker   = {'x','x', 'x'};
colorIdx    = {[0.8 0.8 0.8],[0.8 0.8 0.8], [0.8 0.8 0.8]};
xlimits     = [0.5 3.5];

subplot(1,2,2)
hold on
plotSpread(data,'categoryIdx',xIdx,'distributionIdx',xIdx,...
    'categoryColors',colorIdx,'categoryMarkers',catMarker);
nGroup = 1;
h = errorbar(nGroup,nanmean(dataArray{nGroup}),nanstd(dataArray{nGroup})./sqrt(numel(dataArray{nGroup})-1), ...
    '.-','linewidth',4,'color','k','markersize',0.1); h.CapSize = 16;
nGroup = 2; 
h = errorbar(nGroup,nanmean(dataArray{nGroup}),nanstd(dataArray{nGroup})./sqrt(numel(dataArray{nGroup})-1), ...
    '.-','linewidth',4,'color','k','markersize',0.1); h.CapSize = 16;
nGroup = 3;
h = errorbar(nGroup,nanmean(dataArray{nGroup}),nanstd(dataArray{nGroup})./sqrt(numel(dataArray{nGroup})-1), ...
    '.-','linewidth',4,'color','k','markersize',0.1); h.CapSize = 16;
set(gca,'xtick',1:3,'xticklabel',groupNames)
rotateXLabels(gca,45)
ylabel(yaxName,'fontsize',12)
xlim(xlimits)
ylim(ylimits)
plot(xlimits(1):.5:xlimits(2),zeros(1,numel(xlimits(1):.5:xlimits(2))),'k-');
box off
set(gca,'TickDir','out');
title(titlestr)

xIdx = categorical(xIdx);
[p,t,stats] = anova1(data,xIdx,'off');
text(0.6,ylimits(2)*.9,['group p= ' num2str(p)])
text(0.6,ylimits(2)*.8,['df= ' num2str(stats.df)])
text(0.6,ylimits(2)*.7,['F= ' num2str((t{2,5}))])



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot panel 1C: delta velocity on laser side (on-off)       
clearvars -except allDataTestDMS allDataBaselineDMS ...
                  timeON_BL timeOFF_BL velON_BL velOFF_BL ...
                  timeON timeOFF velON velOFF...
                  diffTimeBL diffTimeTest diffVelBL diffVelTest ...
                  mouseGenotype mouseVirus mouseID laserSide
                                           
%%%%%%%%%%%%%%%%%%% update this section for plotting of measured vars above
%%%%%%%%%%%%%%%%%%% This plots Baseline delta VELOCITY on laser side (on-off)
% remove mice with baseline bias greater than 60s
removeBLbias = 1;
if removeBLbias == 1
    removeBLInd  = abs(diffTimeBL)<60;
else
    removeBLInd  = ones(1,numel(diffTimeBL));
end
% mouseGenotype: 2 indirect (a2a or d2r-cre); 1 direct (d1r-cre)
% mouseVirus: 1 (NpHR); 2 (no opsin)
dataArray        = []; data = []; xIdx = []; 
dataArray{1,1}   = diffVelBL(mouseGenotype==2 & mouseVirus ==1 & removeBLInd==1);
dataArray{1,2}   = diffVelBL(mouseGenotype==1 & mouseVirus ==1 & removeBLInd==1);
dataArray{1,3}   = diffVelBL(mouseVirus ==0 & removeBLInd==1);
data        = [diffVelBL(mouseGenotype==2 & mouseVirus ==1 & removeBLInd==1)'; ...
               diffVelBL(mouseGenotype==1 & mouseVirus ==1 & removeBLInd==1)'; ...
               diffVelBL(mouseVirus ==0 & removeBLInd==1)'];
xIdx        = [ones(size(dataArray{1,1}))';   ...
               ones(size(dataArray{1,2}))'*2; ...
               ones(size(dataArray{1,3}))'*3];
ylimits     = [-30 30];
yaxName     = 'delta speed (cm/s, on-off)';
titlestr    = 'Baseline';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

groupNames  = {'indirect';'direct';'no opsin'};
catMarker   = {'x','x', 'x'};
colorIdx    = {[0.8 0.8 0.8],[0.8 0.8 0.8], [0.8 0.8 0.8]};
xlimits     = [0.5 3.5];

figure
subplot(1,2,1)
hold on
plotSpread(data,'categoryIdx',xIdx,'distributionIdx',xIdx,...
    'categoryColors',colorIdx,'categoryMarkers',catMarker);
nGroup = 1;
h = errorbar(nGroup,nanmean(dataArray{nGroup}),nanstd(dataArray{nGroup})./sqrt(numel(dataArray{nGroup})-1), ...
    '.-','linewidth',4,'color','k','markersize',0.1); h.CapSize = 16;
nGroup = 2; 
h = errorbar(nGroup,nanmean(dataArray{nGroup}),nanstd(dataArray{nGroup})./sqrt(numel(dataArray{nGroup})-1), ...
    '.-','linewidth',4,'color','k','markersize',0.1); h.CapSize = 16;
nGroup = 3;
h = errorbar(nGroup,nanmean(dataArray{nGroup}),nanstd(dataArray{nGroup})./sqrt(numel(dataArray{nGroup})-1), ...
    '.-','linewidth',4,'color','k','markersize',0.1); h.CapSize = 16;

set(gca,'xtick',1:3,'xticklabel',groupNames)
rotateXLabels(gca,45)
ylabel(yaxName,'fontsize',12)
xlim(xlimits)
ylim(ylimits)
plot(xlimits(1):.5:xlimits(2),zeros(1,numel(xlimits(1):.5:xlimits(2))),'k-');
box off
set(gca,'TickDir','out');
title(titlestr)

xIdx = categorical(xIdx);
[p,t,stats] = anova1(data,xIdx,'off');
text(0.6,ylimits(2)*.9,['group p= ' num2str(p)])
text(0.6,ylimits(2)*.8,['df= ' num2str(stats.df)])
text(0.6,ylimits(2)*.7,['F= ' num2str((t{2,5}))])

clearvars -except allDataTestDMS allDataBaselineDMS ...
                  timeON_BL timeOFF_BL velON_BL velOFF_BL ...
                  timeON timeOFF velON velOFF...
                  diffTimeBL diffTimeTest diffVelBL diffVelTest ...
                  mouseGenotype mouseVirus mouseID laserSide removeBLInd
                                           
%%%%%%%%%%%%%%%%%% update this section for plotting of measured vars above
%%%%%%%%%%%%%%%%%% This plots Test delta Time on laser side (on-off)
% remove mice with baseline bias greater than 60s if removeBLbias==1 above
% mouseGenotype: 2 indirect (a2a or d2r-cre); 1 direct (d1r-cre)
% mouseVirus: 1 (NpHR); 2 (no opsin)
dataArray        = []; data = []; xIdx = []; 
dataArray{1,1}   = diffVelTest(mouseGenotype==2 & mouseVirus ==1 & removeBLInd ==1);
dataArray{1,2}   = diffVelTest(mouseGenotype==1 & mouseVirus ==1 & removeBLInd ==1);
dataArray{1,3}   = diffVelTest(mouseVirus ==0 & removeBLInd ==1);
data        = [diffVelTest(mouseGenotype==2 & mouseVirus ==1 & removeBLInd ==1)'; ...
               diffVelTest(mouseGenotype==1 & mouseVirus ==1 & removeBLInd ==1)'; ...
               diffVelTest(mouseVirus ==0 & removeBLInd ==1)'];
xIdx        = [ones(size(dataArray{1,1}))';   ...
               ones(size(dataArray{1,2}))'*2; ...
               ones(size(dataArray{1,3}))'*3];
ylimits     = [-30 30];
yaxName     = 'delta speed (cm/s, on-off)';
titlestr    = 'Test';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

groupNames  = {'indirect';'direct';'no opsin'};
catMarker   = {'x','x', 'x'};
colorIdx    = {[0.8 0.8 0.8],[0.8 0.8 0.8], [0.8 0.8 0.8]};
xlimits     = [0.5 3.5];

subplot(1,2,2)
hold on
plotSpread(data,'categoryIdx',xIdx,'distributionIdx',xIdx,...
    'categoryColors',colorIdx,'categoryMarkers',catMarker);
nGroup = 1;
h = errorbar(nGroup,nanmean(dataArray{nGroup}),nanstd(dataArray{nGroup})./sqrt(numel(dataArray{nGroup})-1), ...
    '.-','linewidth',4,'color','k','markersize',0.1); h.CapSize = 16;
nGroup = 2; 
h = errorbar(nGroup,nanmean(dataArray{nGroup}),nanstd(dataArray{nGroup})./sqrt(numel(dataArray{nGroup})-1), ...
    '.-','linewidth',4,'color','k','markersize',0.1); h.CapSize = 16;
nGroup = 3;
h = errorbar(nGroup,nanmean(dataArray{nGroup}),nanstd(dataArray{nGroup})./sqrt(numel(dataArray{nGroup})-1), ...
    '.-','linewidth',4,'color','k','markersize',0.1); h.CapSize = 16;
set(gca,'xtick',1:3,'xticklabel',groupNames)
rotateXLabels(gca,45)
ylabel(yaxName,'fontsize',12)
xlim(xlimits)
ylim(ylimits)
plot(xlimits(1):.5:xlimits(2),zeros(1,numel(xlimits(1):.5:xlimits(2))),'k-');
box off
set(gca,'TickDir','out');
title(titlestr)

xIdx = categorical(xIdx);
[p,t,stats] = anova1(data,xIdx,'off');
text(0.6,ylimits(2)*.9,['group p= ' num2str(p)])
text(0.6,ylimits(2)*.8,['df= ' num2str(stats.df)])
text(0.6,ylimits(2)*.7,['F= ' num2str((t{2,5}))])


    