%% Script to generate ExtData Fig 9 data plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% cross task comparison of delta (laser on-off) velocity, distance, excess travel,
%  per trial STD, x-position and view angle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
cd(globalParams.dataPath)
load('OffOn_TasksOrGroup_all.mat', 'lgDMSD2_aoe'  ,'lgDMSD2_nd'  ,'lgDMSD2_pc'   , ...
                                   'lgDMSD1_aoe'  ,'lgDMSD1_nd'  ,'lgDMSD1_pc'   , ...
                                   'lgDMSCtrl_aoe','lgDMSCtrl_nd','lgDMSCtrl_pc' );

%%
groupOrder = {'DMSD2_aoe'  ,'DMSD2_nd'  ,'DMSD2_pc'  , ...
              'DMSD1_aoe'  ,'DMSD1_nd'  ,'DMSD1_pc'  , ...
              'DMSCtrl_aoe','DMSCtrl_nd','DMSCtrl_pc'};
counter = 1;
for groupNum = 1:numel(groupOrder)
    clear lg lgClean lgOFF lgON lgCleanOFF lgCleanON
    lg          = eval(['lg' groupOrder{groupNum}]);
    
    % remove field not of same length
    if isfield(lg,'keyFrameLabels')
        lg = rmfield(lg,'keyFrameLabels');
    else
    end
%     lg          = selectLogSubset(lg,[],[],[],[],1,0.6); % drops blocks with perf <0.6
    [lgClean,~] = selectLogSubset(lg,[],[],[],[],0,0.6); % drops blocks with perf <0.6 and excessTravel>0.1
        
    % invert logs to ipsi contra space
    lg        = invertLogs(lg);
    lgClean   = invertLogs(lgClean);
    
    % subselect mazes in aoe and perm cues with matching cue stats
    if contains(groupOrder{groupNum},'DMS') % but only for DMS logs where there is x-task comparisons
        if contains(groupOrder{groupNum},'aoe')
            removeInd = [];
            removeInd = lgClean.currMaze == 11;
            lgClean  = structfun(@(x) x(~removeInd),lgClean,'UniformOutput',false);
            removeInd = [];
            removeInd = lg.currMaze == 11;
            lg       = structfun(@(x) x(~removeInd),lg,'UniformOutput',false);
        elseif contains(groupOrder{groupNum},'pc')
            removeInd = [];
            removeInd = lgClean.currMaze == 7;
            lgClean  = structfun(@(x) x(~removeInd),lgClean,'UniformOutput',false);
            removeInd = [];
            removeInd = lg.currMaze == 7;
            lg       = structfun(@(x) x(~removeInd),lg,'UniformOutput',false);
        else
        end
    else
    end
    
    % remove whole trial inhibition if present in log
    if any(contains(unique(lg.laserEpoch),"whole"))
        removeInd = [];
        removeInd = ismember(lg.laserEpoch,"whole");
        lg       = structfun(@(x) x(~removeInd),lg,'UniformOutput',false);
        removeInd = [];
        removeInd = ismember(lgClean.laserEpoch,"whole");
        lgClean  = structfun(@(x) x(~removeInd),lgClean,'UniformOutput',false);
    else
    end
    
    % create seperate logCleanOFF and logCleanON, which inclues only on "cue" or 0-200cm
    removeMouse = [];
    removeMouse = lgClean.laserON ==0;
    lgCleanOFF = structfun(@(x) x(removeMouse),lgClean,'UniformOutput',false);
    removeMouse = [];
    removeMouse = lgClean.laserON == 1;
    lgCleanON  = structfun(@(x) x(removeMouse),lgClean,'UniformOutput',false);
    removeMouse = [];
    removeMouse = lg.laserON ==0;
    lgOFF      = structfun(@(x) x(removeMouse),lg,'UniformOutput',false);
    removeMouse = [];
    removeMouse = lg.laserON == 1;
    lgON       = structfun(@(x) x(removeMouse),lg,'UniformOutput',false);
    
    % just get number of trials pre and post trial selection here
    nTrialsSubselect_OFF(counter) = numel(lgCleanOFF.choice); 
    nTrialsTotal_OFF(counter)     = numel(lgOFF.choice); 
    nTrialsSubselect_ON(counter)  = numel(lgCleanON.choice);
    nTrialsTotal_ON(counter)      = numel(lgON.choice);
          
    % panel 9A-C - calculation of Yvelocity in bins of 25cm for 300cm stem
    % for laser off and laser on trials seperately
    binIDs = 1:25:301;
    for nBin = 1:numel(binIDs)-1
        yVelOFF_xMouse{counter}(nBin,:)  = xMouseSpeedXY(lgCleanOFF,[binIDs(nBin) binIDs(nBin+1)], 2);
        yVelON_xMouse{counter}(nBin,:)   = xMouseSpeedXY(lgCleanON,[binIDs(nBin) binIDs(nBin+1)], 2);
    end
    
    % panel 9D-F,left - calculation of distance
    % for laser off and laser on trials seperately
    distance = [];
    distance = 330+(lg.excessTravel*330);
    mouseNums = [];
    mouseNums = unique(lg.mouseID);
    for nMouse = 1:numel(mouseNums)
        distanceOFF_xMouse{counter}(:,nMouse) = nanmean(distance(lg.mouseID == mouseNums(nMouse)...
            & lg.laserON==0));
        distanceON_xMouse{counter}(:,nMouse)  = nanmean(distance(lg.mouseID == mouseNums(nMouse)...
            & lg.laserON==1));
    end
    
    % panel 9D-F,right - calculation of % trials with excessTravel
    % for laser off and laser on trials seperately
    mouseNums = []; etOFF = []; etON = [];
    mouseNums = unique(lg.mouseID);
    for nMouse = 1:numel(mouseNums)
        etOFF(nMouse) = (sum(lg.excessTravel(lg.mouseID == mouseNums(nMouse) & lg.laserON == 0) > .1 )/...
            numel(lg.excessTravel(lg.mouseID == mouseNums(nMouse) & lg.laserON == 0)))*100;
        etON(nMouse)  = (sum(lg.excessTravel(lg.mouseID == mouseNums(nMouse) & lg.laserON == 1) > .1 )/...
            numel(lg.excessTravel(lg.mouseID == mouseNums(nMouse) & lg.laserON == 1)))*100;
    end
    excTravelOFF_xMouse{counter} = etOFF;
    excTravelON_xMouse{counter}  = etON;
    
    % panel 9G-I - calculation of % trials with excessTravel
    posBins = []; posBins = 1:5:301;
    vadOFF                      = [];
    vadOFF                      = xMouseViewAngSD(lgOFF, posBins);
    vaDevOFF_xMouse{counter}(:) = vadOFF.viewAngSD;
    vadON                       = [];
    vadON                       = xMouseViewAngSD(lgON, posBins);
    vaDevON_xMouse{counter}(:)  = vadON.viewAngSD;
    
    % panel 9J-K,left - calculation of x-position (cm) in maze stem 0-200cm    
    avgXpos = [];
    mazeSection = []; mazeSection = [1 201];
    pos2    = []; % next line inverts matrix columns in each cell array of pos
    % hack workaround as next function samples 3rd column based
    % on 2nd column
    % pos2 is view angle, y-position, x-position
    % pos is x-position, y-position, view angle
    pos2    = cellfun(@fliplr, lgClean.pos,'UniformOutput',false); %
    avgXpos = sampleViewAngleVsY_average(pos2,mazeSection,mazeSection(2)-mazeSection(1));
    mouseNums  = [];
    mouseNums  = unique(lgClean.mouseID);
    for nMouse = 1:numel(mouseNums)
        xPos_OFF_xMouse{counter}(:,nMouse) = ...
            nanmean(avgXpos(:,lgClean.mouseID == mouseNums(nMouse) & ...
            lgClean.laserON == 0),2);
        xPos_ON_xMouse{counter}(:,nMouse) = ...
            nanmean(avgXpos(:,lgClean.mouseID == mouseNums(nMouse) & ...
            lgClean.laserON == 1),2);
    end
    
    % panel 9J-K,right - calculation of view angle (deg) in maze stem 0-200cm
    avgViewAng  = [];
    avgViewAng  = sampleViewAngleVsY_average(lgClean.pos,mazeSection,mazeSection(2)-mazeSection(1));
    mouseNums       = [];
    mouseNums       = unique(lgClean.mouseID);
    for nMouse = 1:numel(mouseNums)
        viewAng_OFF_xMouse{counter}(:,nMouse) = ...
            nanmean(avgViewAng(:,lgClean.mouseID == mouseNums(nMouse) & ...
            lgClean.laserON == 0),2);
        viewAng_ON_xMouse{counter}(:,nMouse) = ...
            nanmean(avgViewAng(:,lgClean.mouseID == mouseNums(nMouse) & ...
            lgClean.laserON == 1),2);
    end
    
    % for the next loop of analysis given lgs in groupOrder
    counter = counter +1;
end

%% clear all vars except source data for plotting
clearvars -except groupOrder nTrialsSubselect_OFF nTrialsTotal_OFF ...
                  nTrialsSubselect_ON nTrialsTotal_ON ...
                  lgDMSCtrl_aoe lgDMSCtrl_nd lgDMSCtrl_pc ...
                  lgDMSD1_aoe lgDMSD1_nd lgDMSD1_pc ...
                  lgDMSD2_aoe lgDMSD2_nd lgDMSD2_pc ...
                  yVelOFF_xMouse yVelON_xMouse ...
                  distanceOFF_xMouse distanceON_xMouse ... 
                  excTravelOFF_xMouse excTravelON_xMouse...
                  vaDevOFF_xMouse vaDevON_xMouse ...
                  xPos_OFF_xMouse xPos_ON_xMouse ...
                  viewAng_OFF_xMouse viewAng_ON_xMouse
              
%% PLOTTING BASED ON ExtData Fig 9 source data above
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% cross task comparison of delta (laser on-off) velocity, distance, 
%  excess travel, per trial STD, x-position and view angle

%%%%%%%%%%%%%%%%%%%%%%%%%%
% update this section for plotting of measured vars above
% plot panel ext data 9a : y-velocity off/on indirect             
names2plot  = {'DMSD2_aoe','DMSD2_nd','DMSD2_pc'};                  
groups2plot = [1 2 3];  
ylimits     = [0 300];
xlimits     = [20 100];
yaxName     = 'y-velocity (cm/s, on-off)';
%%%%%%%%%%%%%%%%%%%%%%%%%
dataArrayOFF   = []; dataOFF = []; xIdx = []; 
dataArrayON    = []; dataON  = []; dataArray = [];
dataArray      = cellfun(@minus,yVelON_xMouse,yVelOFF_xMouse,'UniformOutput',false);
dataArrayOFF   = [yVelOFF_xMouse(groups2plot(1)) ...
                  yVelOFF_xMouse(groups2plot(2)) ...
                  yVelOFF_xMouse(groups2plot(3))];
dataArrayON    = [yVelON_xMouse(groups2plot(1)) ...
                  yVelON_xMouse(groups2plot(2)) ...
                  yVelON_xMouse(groups2plot(3))];              
dataOFF        = [yVelOFF_xMouse{groups2plot(1)}'; ...
                  yVelOFF_xMouse{groups2plot(2)}'; ...
                  yVelOFF_xMouse{groups2plot(3)}'];
dataON         = [yVelON_xMouse{groups2plot(1)}'; ...
                  yVelON_xMouse{groups2plot(2)}'; ...
                  yVelON_xMouse{groups2plot(3)}'];              
xIdx        = [ones(size(dataArray{groups2plot(1)},2),1);   ...
               ones(size(dataArray{groups2plot(2)},2),1)*2; ...
               ones(size(dataArray{groups2plot(3)},2),1)*3];
figure
binIDs    = 25:25:300;
offOffset = -7.5;
onOffset  = -12.5;
for nGroup = 1:numel(groups2plot)
    subplot(1,3,nGroup)
    mean2plot = []; err2plot = [];
    mean2plot = nanmean(dataArrayOFF{nGroup},2);  
    err2plot  = (nanstd(dataArrayOFF{nGroup},0,2)) / sqrt(size(dataArrayOFF{nGroup},2)-1);  
    errorbar(mean2plot, binIDs' + offOffset ,err2plot,'horizontal','k-'); hold on
    mean2plot = []; err2plot = [];
    mean2plot = nanmean(dataArrayON{nGroup},2);  
    err2plot  = (nanstd(dataArrayON{nGroup},0,2)) / sqrt(size(dataArrayON{nGroup},2)-1);  
    errorbar(mean2plot, binIDs' + onOffset ,err2plot,'horizontal','-','Color',[.4 1 .4]); hold on   
    xlim(xlimits) 
    set(gca,'XTick',[25 50 75 100])
    box off
    set(gca,'TickDir','out')
    if contains(groupOrder{groups2plot(nGroup)},'aoe')
        title('AoE')    
    elseif contains(groupOrder{groups2plot(nGroup)},'nd')
        title('ctrl#1')
    elseif contains(groupOrder{groups2plot(nGroup)},'pc')
        title('ctrl#2')
    else
    end
    if nGroup == 1
        ylabel('y position (cm)')
    else
        set(gca,'YTickLabel',[])
    end
    
    if nGroup == 2
        xlabel('y velocity (cm/s)')
    end
end

% plot panel ext data 9b : y-velocity off/on direct             
%%%%%%%%%%%%%%%%%%%%%%%%%%
% update this section for plotting of measured vars above
names2plot  = {'DMSD1_aoe','DMSD1_nd','DMSD1_pc'};                  
groups2plot = [4 5 6];  
ylimits     = [0 300];
xlimits     = [20 100];
yaxName     = 'y-velocity (cm/s, on-off)';
%%%%%%%%%%%%%%%%%%%%%%%%%
dataArrayOFF   = []; dataOFF = []; xIdx = []; 
dataArrayON    = []; dataON  = []; dataArray = [];
dataArray      = cellfun(@minus,yVelON_xMouse,yVelOFF_xMouse,'UniformOutput',false);
dataArrayOFF   = [yVelOFF_xMouse(groups2plot(1)) ...
                  yVelOFF_xMouse(groups2plot(2)) ...
                  yVelOFF_xMouse(groups2plot(3))];
dataArrayON    = [yVelON_xMouse(groups2plot(1)) ...
                  yVelON_xMouse(groups2plot(2)) ...
                  yVelON_xMouse(groups2plot(3))];              
dataOFF        = [yVelOFF_xMouse{groups2plot(1)}'; ...
                  yVelOFF_xMouse{groups2plot(2)}'; ...
                  yVelOFF_xMouse{groups2plot(3)}'];
dataON         = [yVelON_xMouse{groups2plot(1)}'; ...
                  yVelON_xMouse{groups2plot(2)}'; ...
                  yVelON_xMouse{groups2plot(3)}'];              
xIdx        = [ones(size(dataArray{groups2plot(1)},2),1);   ...
               ones(size(dataArray{groups2plot(2)},2),1)*2; ...
               ones(size(dataArray{groups2plot(3)},2),1)*3];

figure
binIDs    = 25:25:300;
offOffset = -7.5;
onOffset  = -12.5;
for nGroup = 1:numel(groups2plot)
    subplot(1,3,nGroup)
    mean2plot = []; err2plot = [];
    mean2plot = nanmean(dataArrayOFF{nGroup},2);  
    err2plot  = (nanstd(dataArrayOFF{nGroup},0,2)) / sqrt(size(dataArrayOFF{nGroup},2)-1);  
    errorbar(mean2plot, binIDs' + offOffset ,err2plot,'horizontal','k-'); hold on
    mean2plot = []; err2plot = [];
    mean2plot = nanmean(dataArrayON{nGroup},2);  
    err2plot  = (nanstd(dataArrayON{nGroup},0,2)) / sqrt(size(dataArrayON{nGroup},2)-1);  
    errorbar(mean2plot, binIDs' + onOffset ,err2plot,'horizontal','-','Color',[.4 1 .4]); hold on   
    xlim(xlimits) 
    set(gca,'XTick',[25 50 75 100])
    box off
    set(gca,'TickDir','out')
    if contains(groupOrder{groups2plot(nGroup)},'aoe')
        title('AoE')    
    elseif contains(groupOrder{groups2plot(nGroup)},'nd')
        title('ctrl#1')
    elseif contains(groupOrder{groups2plot(nGroup)},'pc')
        title('ctrl#2')
    else
    end
    if nGroup == 1
        ylabel('y position (cm)')
    else
        set(gca,'YTickLabel',[])
    end
    
    if nGroup == 2
        xlabel('y velocity (cm/s)')
    end
end

% plot panel ext data 9c : y-velocity off/on no opsin             
%%%%%%%%%%%%%%%%%%%%%%%%%%
% update this section for plotting of measured vars above
names2plot  = {'DMSCtrl_aoe','DMSCtrl_nd','DMSCtrl_pc'};                  
groups2plot = [7 8 9];  
ylimits     = [0 300];
xlimits     = [20 100];
yaxName     = 'y-velocity (cm/s, on-off)';
%%%%%%%%%%%%%%%%%%%%%%%%%
dataArrayOFF   = []; dataOFF = []; xIdx = []; 
dataArrayON    = []; dataON  = []; dataArray = [];
dataArray      = cellfun(@minus,yVelON_xMouse,yVelOFF_xMouse,'UniformOutput',false);
dataArrayOFF   = [yVelOFF_xMouse(groups2plot(1)) ...
                  yVelOFF_xMouse(groups2plot(2)) ...
                  yVelOFF_xMouse(groups2plot(3))];
dataArrayON    = [yVelON_xMouse(groups2plot(1)) ...
                  yVelON_xMouse(groups2plot(2)) ...
                  yVelON_xMouse(groups2plot(3))];              
dataOFF        = [yVelOFF_xMouse{groups2plot(1)}'; ...
                  yVelOFF_xMouse{groups2plot(2)}'; ...
                  yVelOFF_xMouse{groups2plot(3)}'];
dataON         = [yVelON_xMouse{groups2plot(1)}'; ...
                  yVelON_xMouse{groups2plot(2)}'; ...
                  yVelON_xMouse{groups2plot(3)}'];              
xIdx        = [ones(size(dataArray{groups2plot(1)},2),1);   ...
               ones(size(dataArray{groups2plot(2)},2),1)*2; ...
               ones(size(dataArray{groups2plot(3)},2),1)*3];

figure
binIDs    = 25:25:300;
offOffset = -7.5;
onOffset  = -12.5;
for nGroup = 1:numel(groups2plot)
    subplot(1,3,nGroup)
    mean2plot = []; err2plot = [];
    mean2plot = nanmean(dataArrayOFF{nGroup},2);  
    err2plot  = (nanstd(dataArrayOFF{nGroup},0,2)) / sqrt(size(dataArrayOFF{nGroup},2)-1);  
    errorbar(mean2plot, binIDs' + offOffset ,err2plot,'horizontal','k-'); hold on
    mean2plot = []; err2plot = [];
    mean2plot = nanmean(dataArrayON{nGroup},2);  
    err2plot  = (nanstd(dataArrayON{nGroup},0,2)) / sqrt(size(dataArrayON{nGroup},2)-1);  
    errorbar(mean2plot, binIDs' + onOffset ,err2plot,'horizontal','-','Color',[.4 1 .4]); hold on   
    xlim(xlimits) 
    set(gca,'XTick',[25 50 75 100])
    box off
    set(gca,'TickDir','out')
    if contains(groupOrder{groups2plot(nGroup)},'aoe')
        title('AoE')    
    elseif contains(groupOrder{groups2plot(nGroup)},'nd')
        title('ctrl#1')
    elseif contains(groupOrder{groups2plot(nGroup)},'pc')
        title('ctrl#2')
    else
    end
    if nGroup == 1
        ylabel('y position (cm)')
    else
        set(gca,'YTickLabel',[])
    end
    
    if nGroup == 2
        xlabel('y velocity (cm/s)')
    end
end

%% clear all vars except source data for plotting
clearvars -except groupOrder nTrialsSubselect_OFF nTrialsTotal_OFF ...
                  nTrialsSubselect_ON nTrialsTotal_ON ...
                  lgDMSCtrl_aoe lgDMSCtrl_nd lgDMSCtrl_pc ...
                  lgDMSD1_aoe lgDMSD1_nd lgDMSD1_pc ...
                  lgDMSD2_aoe lgDMSD2_nd lgDMSD2_pc ...
                  yVelOFF_xMouse yVelON_xMouse ...
                  distanceOFF_xMouse distanceON_xMouse ... 
                  excTravelOFF_xMouse excTravelON_xMouse...
                  vaDevOFF_xMouse vaDevON_xMouse ...
                  xPos_OFF_xMouse xPos_ON_xMouse ...
                  viewAng_OFF_xMouse viewAng_ON_xMouse
              
%% plot panels 9d-f delta distance and excess travel indirect direct no opsin
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

groups2plot = [1 2 3]; % indirect aoe,ctrl#1,ctrl#2
ylimits1 = [-200 200 ];
ylimits2 = [-15 15 ];
posthocStats = 0;
plotExtData9(distanceOFF_xMouse,distanceON_xMouse,...
                excTravelOFF_xMouse,excTravelON_xMouse, ...
                groups2plot,ylimits1,ylimits2,posthocStats);
groups2plot = [4 5 6]; % direct aoe,ctrl#1,ctrl#2
plotExtData9(distanceOFF_xMouse,distanceON_xMouse,...
                excTravelOFF_xMouse,excTravelON_xMouse, ...
                groups2plot,ylimits1,ylimits2,posthocStats);
groups2plot = [7 8 9]; % no opsin aoe,ctrl#1,ctrl#2
plotExtData9(distanceOFF_xMouse,distanceON_xMouse,...
                excTravelOFF_xMouse,excTravelON_xMouse, ...
                groups2plot,ylimits1,ylimits2,posthocStats);            
            
%% clear all vars except source data for plotting
clearvars -except groupOrder nTrialsSubselect_OFF nTrialsTotal_OFF ...
                  nTrialsSubselect_ON nTrialsTotal_ON ...
                  lgDMSCtrl_aoe lgDMSCtrl_nd lgDMSCtrl_pc ...
                  lgDMSD1_aoe lgDMSD1_nd lgDMSD1_pc ...
                  lgDMSD2_aoe lgDMSD2_nd lgDMSD2_pc ...
                  yVelOFF_xMouse yVelON_xMouse ...
                  distanceOFF_xMouse distanceON_xMouse ... 
                  excTravelOFF_xMouse excTravelON_xMouse...
                  vaDevOFF_xMouse vaDevON_xMouse ...
                  xPos_OFF_xMouse xPos_ON_xMouse ...
                  viewAng_OFF_xMouse viewAng_ON_xMouse
              
              
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot panels 9g-i delta per-trial deviation in view angle

%%%%%%%%%%%%%%%%%%%%%%%%%%
% update this section for plotting of measured vars above
% panel 9g-i
dataArray   = []; data = []; xIdx = []; 
groups2plot = [1 2 3]; % indirect aoe,nd,pc tasks
dataArray   = cellfun(@minus,vaDevON_xMouse,vaDevOFF_xMouse,'UniformOutput',false);
data        = [dataArray{groups2plot(1)}'; ...
               dataArray{groups2plot(2)}'; ...
               dataArray{groups2plot(3)}'];
xIdx        = [ones(size(dataArray{groups2plot(1)}))';   ...
               ones(size(dataArray{groups2plot(2)}))'*2; ...
               ones(size(dataArray{groups2plot(3)}))'*3];
yaxName     = 'delta trials w/ excess travel (%, on-off)';
ylimits     = [-150 150];
titlstr     = 'indirect';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

groupNames  = {'aoe';'ctrl#1';'ctrl#2'};
catMarker   = {'x','x', 'x'};
colorIdx    = {[0.8 0.8 0.8],[0.99 0.7 0.99], [0.7 0.99 0.99]};
xlimits     = [0.5 3.5];

figure
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
ylim(ylimits)
plot(xlimits(1):.5:xlimits(2),zeros(1,numel(xlimits(1):.5:xlimits(2))),'k-');
box off
set(gca,'TickDir','out');
title(titlstr)

xIdx = categorical(xIdx);
[p,t,stats] = anova1(data,xIdx,'off');
text(0.6,ylimits(1)*.85,['group p= ' num2str(p)])
text(0.6,ylimits(1)*.9,['df= ' num2str(stats.df)])
text(0.6,ylimits(1)*.95,['F= ' num2str((t{2,5}))])        

%%%%%%%%%%%%%%%%%%%%%%%%%%
% update this section for plotting of measured vars above
% panel 9g-i
dataArray   = []; data = []; xIdx = []; 
groups2plot = [4 5 6]; % indirect aoe,nd,pc tasks
dataArray   = cellfun(@minus,vaDevON_xMouse,vaDevOFF_xMouse,'UniformOutput',false);
data        = [dataArray{groups2plot(1)}'; ...
               dataArray{groups2plot(2)}'; ...
               dataArray{groups2plot(3)}'];
xIdx        = [ones(size(dataArray{groups2plot(1)}))';   ...
               ones(size(dataArray{groups2plot(2)}))'*2; ...
               ones(size(dataArray{groups2plot(3)}))'*3];
yaxName     = 'delta trials w/ excess travel (%, on-off)';
ylimits     = [-150 150];
titlstr     = 'direct';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

groupNames  = {'aoe';'ctrl#1';'ctrl#2'};
catMarker   = {'x','x', 'x'};
colorIdx    = {[0.8 0.8 0.8],[0.99 0.7 0.99], [0.7 0.99 0.99]};
xlimits     = [0.5 3.5];

figure
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
ylim(ylimits)
plot(xlimits(1):.5:xlimits(2),zeros(1,numel(xlimits(1):.5:xlimits(2))),'k-');
box off
set(gca,'TickDir','out');
title(titlstr)

xIdx = categorical(xIdx);
[p,t,stats] = anova1(data,xIdx,'off');
text(0.6,ylimits(1)*.85,['group p= ' num2str(p)])
text(0.6,ylimits(1)*.9,['df= ' num2str(stats.df)])
text(0.6,ylimits(1)*.95,['F= ' num2str((t{2,5}))])       

%%%%%%%%%%%%%%%%%%%%%%%%%%
% update this section for plotting of measured vars above
% panel 9g-i
dataArray   = []; data = []; xIdx = []; 
groups2plot = [7 8 9]; % no opsin aoe,nd,pc tasks
dataArray   = cellfun(@minus,vaDevON_xMouse,vaDevOFF_xMouse,'UniformOutput',false);
data        = [dataArray{groups2plot(1)}'; ...
               dataArray{groups2plot(2)}'; ...
               dataArray{groups2plot(3)}'];
xIdx        = [ones(size(dataArray{groups2plot(1)}))';   ...
               ones(size(dataArray{groups2plot(2)}))'*2; ...
               ones(size(dataArray{groups2plot(3)}))'*3];
yaxName     = 'delta trials w/ excess travel (%, on-off)';
ylimits     = [-150 150];
titlstr     = 'no opsin';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

groupNames  = {'aoe';'ctrl#1';'ctrl#2'};
catMarker   = {'x','x', 'x'};
colorIdx    = {[0.8 0.8 0.8],[0.99 0.7 0.99], [0.7 0.99 0.99]};
xlimits     = [0.5 3.5];

figure
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
ylim(ylimits)
plot(xlimits(1):.5:xlimits(2),zeros(1,numel(xlimits(1):.5:xlimits(2))),'k-');
box off
set(gca,'TickDir','out');
title(titlstr)

xIdx = categorical(xIdx);
[p,t,stats] = anova1(data,xIdx,'off');
text(0.6,ylimits(1)*.85,['group p= ' num2str(p)])
text(0.6,ylimits(1)*.9,['df= ' num2str(stats.df)])
text(0.6,ylimits(1)*.95,['F= ' num2str((t{2,5}))])       
              
%% clear all vars except source data for plotting
clearvars -except groupOrder nTrialsSubselect_OFF nTrialsTotal_OFF ...
                  nTrialsSubselect_ON nTrialsTotal_ON ...
                  lgDMSCtrl_aoe lgDMSCtrl_nd lgDMSCtrl_pc ...
                  lgDMSD1_aoe lgDMSD1_nd lgDMSD1_pc ...
                  lgDMSD2_aoe lgDMSD2_nd lgDMSD2_pc ...
                  yVelOFF_xMouse yVelON_xMouse ...
                  distanceOFF_xMouse distanceON_xMouse ... 
                  excTravelOFF_xMouse excTravelON_xMouse...
                  vaDevOFF_xMouse vaDevON_xMouse ...
                  xPos_OFF_xMouse xPos_ON_xMouse ...
                  viewAng_OFF_xMouse viewAng_ON_xMouse              
            
%% plot panels 9j-l delta xPos and view angle indirect direct no opsin
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

groups2plot = [1 2 3];
ylimits1 = [-0.5 0.5 ];
ylimits2 = [-15 15 ];
posthocStats = 1;
plotExtData9(xPos_OFF_xMouse,xPos_ON_xMouse, ...
                viewAng_OFF_xMouse,viewAng_ON_xMouse,...
                groups2plot,ylimits1,ylimits2,posthocStats);
groups2plot = [4 5 6];
plotExtData9(xPos_OFF_xMouse,xPos_ON_xMouse, ...
                viewAng_OFF_xMouse,viewAng_ON_xMouse,...
                groups2plot,ylimits1,ylimits2,posthocStats);
groups2plot = [7 8 9];
plotExtData9(xPos_OFF_xMouse,xPos_ON_xMouse, ...
                viewAng_OFF_xMouse,viewAng_ON_xMouse,...
                groups2plot,ylimits1,ylimits2,posthocStats);       

%% clear all vars except source data for plotting
clearvars -except groupOrder nTrialsSubselect_OFF nTrialsTotal_OFF ...
                  nTrialsSubselect_ON nTrialsTotal_ON ...
                  lgDMSCtrl_aoe lgDMSCtrl_nd lgDMSCtrl_pc ...
                  lgDMSD1_aoe lgDMSD1_nd lgDMSD1_pc ...
                  lgDMSD2_aoe lgDMSD2_nd lgDMSD2_pc ...
                  yVelOFF_xMouse yVelON_xMouse ...
                  distanceOFF_xMouse distanceON_xMouse ... 
                  excTravelOFF_xMouse excTravelON_xMouse...
                  vaDevOFF_xMouse vaDevON_xMouse ...
                  xPos_OFF_xMouse xPos_ON_xMouse ...
                  viewAng_OFF_xMouse viewAng_ON_xMouse














