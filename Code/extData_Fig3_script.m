%% Script to generate ExtData Fig 6B-J data plots
% related to Fig 2 motor comparisons - extra motor plots and on-off stats
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
cd(globalParams.dataPath)
load('xTask_allOff_matchCueStats.mat')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% analyzes log files concatenated across all DMS groups performing accumulation of
% evidence(aoe), no distractors (nd), or permanent cues(pc) tasks
% (1) cue statistics of aoe (log.currMaze==10) and pc (log.currMaze==7) pre-selected to match
% (2) only laser off trials included in logs. 
% (3) all analyses subselect trial blocks with performance >60%
% (4) all anlyses except distance subselect trials with excessTravel<0.1 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
warning('off','all')
groupOrder = {'_aoe', '_nd', '_pc'};
counter = 1;
for groupNum = 1:numel(groupOrder)
    clear lg lgClean 
    lg          = eval(['concatLog_OFF' groupOrder{groupNum}]); 
    lg          = selectLogSubset(lg,[],[],[],[],1,0.6); % drops blocks with perf <0.6
    [lgClean,~] = selectLogSubset(lg,[],[],[],[],0,0.6); % drops blocks with perf <0.6 and excessTravel>0.1
           
    % just get number of trials pre and post trial selection here
    ipsiCon = 0; cleanLog = 1; relaxTravel = 0;
    tempAcc = []; tempSumm = [];
    [~, tempSumm] = xMousePerfBias(lg, ipsiCon, cleanLog, relaxTravel);
    nTrialsSubselect_xTask(counter) = tempSumm.nTrialsPostSelection;
    nTrialsTotal_xTask(counter)     = tempSumm.nTrialsPreSelection;

    % panel ExtData 6B - calculation of Yvelocity (cm/s) in maze stem 0-200cm  
    mazeSection = []; mazeSection = [1 201];
    yVel_xMouse{counter}(:)   = xMouseSpeedXY(lgClean, mazeSection, 2);
               
    % panel ExtData 6C - calculation of x-position (cm) in maze stem 0-200cm
    % for left and right choice trials seperately
    avgXpos = [];
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
        xPosL_xMouse{counter}(:,nMouse) = ...
            nanmean(avgXpos(:,lgClean.mouseID == mouseNums(nMouse) & ...
                              lgClean.choice == 0),2);
        xPosR_xMouse{counter}(:,nMouse) = ...
            nanmean(avgXpos(:,lgClean.mouseID == mouseNums(nMouse) & ...
                              lgClean.choice == 1),2);
    end
    
    % panel ExtData 6D - calculation of view angle (deg) in maze stem 0-200cm 
    % for left and right choice trials seperately
    avgViewAng  = [];
    avgViewAng  = sampleViewAngleVsY_average(lgClean.pos,mazeSection,mazeSection(2)-mazeSection(1));
    
    mouseNums       = [];
    mouseNums       = unique(lgClean.mouseID);
    for nMouse = 1:numel(mouseNums)
        viewAngL_xMouse{counter}(:,nMouse) = ...
            nanmean(avgViewAng(:,lgClean.mouseID == mouseNums(nMouse) & ...
                              lgClean.choice == 0),2);
        viewAngR_xMouse{counter}(:,nMouse) = ...
            nanmean(avgViewAng(:,lgClean.mouseID == mouseNums(nMouse) & ...
                              lgClean.choice == 1),2);
    end
         
    % panel ExtData 6E - calculation of % trials with excessTravel         
    mouseNums = []; et = [];
    mouseNums = unique(lg.mouseID);
    for nMouse = 1:numel(mouseNums)
        et(nMouse) = (sum(lg.excessTravel(lg.mouseID == mouseNums(nMouse)) > 1.1)/...
                         numel(lg.excessTravel(lg.mouseID == mouseNums(nMouse))))*100;
        
    end   
    excTravel_xMouse{counter} = et;
    
    % panel ExtData 6F - calculation of % trials with excessTravel
    posBins = []; posBins = 1:5:301;  
    vad     = []; vad  = xMouseViewAngSD(lg, posBins);   
    vaDev_xMouse{counter}(:)  = vad.viewAngSD;       
    
    % panel ExtData 6H - choice decoding from viewAngle (train/test within mouse)
    % lg.pos is 3 columns of x position, y position, view angle - function
    % samples 3rd column based on 2nd column
    blockAvg      = 25;   
    minMaxPos     = [1 300];
    numFolds      = 5;
    numSamples    = 10;
    fitType       = 'fitglm';
    cleanTrials   = 1;   
    weighted      = 1;
    zPred         = 0;
    ipCon         = 0;
    optimalThresh = 1;
    predChoice    = [];
    % this takes a little while
    predChoice    = decodeChoiceFromTrajectory(lg,blockAvg,minMaxPos,numFolds,...
                    numSamples,fitType,cleanTrials,weighted,zPred,ipCon,optimalThresh);
    decodeChoiceFromVA_xMouse{counter} = predChoice;            
    
    % panel ExtData 6G - choice decoding from x-position (train/test within mouse)    
    pos2 = [];    pos2 = cellfun(@fliplr, lg.pos,'UniformOutput',false); 
    lg.pos       = pos2; % hack to flip pos to 3 columns of view ang, ypos, xpos (below samples 3rd column baed on 2nd)
    blockAvg      = 25;   
    minMaxPos     = [1 300];
    numFolds      = 5;
    numSamples    = 10;
    fitType       = 'fitglm';
    cleanTrials   = 1;   
    weighted      = 1;
    zPred         = 0;
    ipCon         = 0;
    optimalThresh = 1;
    predChoice    = [];
    % this takes a little while
    predChoice    = decodeChoiceFromTrajectory(lg,blockAvg,minMaxPos,numFolds,...
                    numSamples,fitType,cleanTrials,weighted,zPred,ipCon,optimalThresh);
   decodeChoiceFromXpos_xMouse{counter} = predChoice;      
        
    % for the next loop of analysis given logs in groupOrder       
    counter = counter +1;
end
       
%% analysis for panel ExtData 6J - decode choice from view angle x tasks
% concat all task logs into one masterLog first
fieldsToConcat  = fieldnames(concatLog_OFF_aoe);
masterLog = [];
for nField = 1:numel(fieldsToConcat)
    masterLog.(fieldsToConcat{nField}) = [concatLog_OFF_aoe.(fieldsToConcat{nField}), ...
                                          concatLog_OFF_nd.(fieldsToConcat{nField}) , ...
                                          concatLog_OFF_pc.(fieldsToConcat{nField})];
end       

blockAvg               = 25;
minMaxPos              = [1 300];
numFolds               = 2;
numSamples             = 50;
fitType                = 'fitglm';
cleanTrials            = 1;
weighted               = 1;
zPred                  = 0;
ipCon                  = 0;
optimalThresh          = 0;
[predChoiceXva_xTask]  = ...
    decodeChoiceFromTraject_xTasksAndMice(masterLog,...
                                          blockAvg,...
                                          minMaxPos,...
                                          numFolds,...
                                          numSamples,...
                                          fitType,...
                                          cleanTrials,...
                                          weighted,...
                                          zPred,...
                                          ipCon,...
                                          optimalThresh);
%% analysis for panel ExtData 6I - decode choice from x-position x tasks
pos2                   = [];    
pos2                   = cellfun(@fliplr, masterLog.pos,'UniformOutput',false); 
masterLog.pos          = pos2; % hack to flip pos to 3 columns of view ang, ypos, xpos (below samples 3rd column baed on 2nd)
blockAvg               = 25;
minMaxPos              = [1 300];
numFolds               = 2;
numSamples             = 50;
fitType                = 'fitglm';
cleanTrials            = 1;
weighted               = 1;
zPred                  = 0;
ipCon                  = 0;
optimalThresh          = 0;
[predChoiceXpos_xTask] = ...
    decodeChoiceFromTraject_xTasksAndMice(masterLog,...
                                          blockAvg,...
                                          minMaxPos,...
                                          numFolds,...
                                          numSamples,...
                                          fitType,...
                                          cleanTrials,...
                                          weighted,...
                                          zPred,...
                                          ipCon,...
                                          optimalThresh);

warning('on','all')

%% clear all vars except source data for plotting
clearvars -except groupOrder nTrialsSubselect_xTask nTrialsTotal_xTask ...
           yVel_xMouse xPosL_xMouse xPosR_xMouse viewAngL_xMouse viewAngR_xMouse excTravel_xMouse ...
           vaDev_xMouse decodeChoiceFromVA_xMouse decodeChoiceFromXpos_xMouse ...
           predChoiceXpos_xTask predChoiceXva_xTask ...
           concatLog_OFF_aoe concatLog_OFF_nd concatLog_OFF_pc masterLog
       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plotting  - here is panel 6B y-velocity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot panel 6B : y-velocity 0-200cm of maze stem          
clearvars -except groupOrder nTrialsSubselect_xTask nTrialsTotal_xTask ...
           yVel_xMouse xPosL_xMouse xPosR_xMouse viewAngL_xMouse viewAngR_xMouse excTravel_xMouse ...
           vaDev_xMouse decodeChoiceFromVA_xMouse decodeChoiceFromXpos_xMouse ...
           predChoiceXpos_xTask predChoiceXva_xTask ...
           concatLog_OFF_aoe concatLog_OFF_nd concatLog_OFF_pc masterLog
                                           
%%%%%%%%%%%%%%%%%%%%%%%%%%
% update this section for plotting of measured vars above
dataArray   = []; data = []; xIdx = []; 
dataArray   = yVel_xMouse;
data        = [yVel_xMouse{1}'; ...
               yVel_xMouse{2}'; ...
               yVel_xMouse{3}'];
xIdx        = [ones(size(yVel_xMouse{1}))';   ...
               ones(size(yVel_xMouse{2}))'*2; ...
               ones(size(yVel_xMouse{3}))'*3];
ylimits     = [20 120];
yaxName     = 'y-velocity (cm/s, 0-200cm)';
%%%%%%%%%%%%%%%%%%%%%%%%%

groupNames  = {'aoe';'ctrl#1';'ctrl#2'};
catMarker   = {'x','x', 'x'};
colorIdx    = {[0.8 0.8 0.8],[0.99 0.7 0.99], [0.7 0.99 0.99]};
xlimits     = [0.5 3.5];

figure
hold on
plotSpread(data,'categoryIdx',xIdx,'distributionIdx',xIdx,...
    'categoryColors',colorIdx,'categoryMarkers',catMarker);

nGroup = 1;
h = errorbar(nGroup,nanmean(dataArray{nGroup}),nanstd(dataArray{nGroup})./sqrt(numel(dataArray{nGroup})-1), ...
    '.-','linewidth',4,'color','k','markersize',0.1); h.CapSize = 16;

nGroup = 2; 
h = errorbar(nGroup,nanmean(dataArray{nGroup}),nanstd(dataArray{nGroup})./sqrt(numel(dataArray{nGroup})-1), ...
    '.-','linewidth',4,'color','m','markersize',0.1); h.CapSize = 16;

nGroup = 3;
h = errorbar(nGroup,nanmean(dataArray{nGroup}),nanstd(dataArray{nGroup})./sqrt(numel(dataArray{nGroup})-1), ...
    '.-','linewidth',4,'color','c','markersize',0.1); h.CapSize = 16;

set(gca,'xtick',1:3,'xticklabel',groupNames)
rotateXLabels(gca,45)
ylabel(yaxName,'fontsize',12)
xlim(xlimits)
ylim(ylimits)
% plot(xlimits(1):.5:xlimits(2),zeros(1,numel(xlimits(1):.5:xlimits(2))),'k-');
box off
set(gca,'TickDir','out');

xIdx = categorical(xIdx);
[p,t,stats] = anova1(data,xIdx,'off');
text(0.6,ylimits(2)*.95,['group p= ' num2str(p)])
text(0.6,ylimits(2)*.9,['df= ' num2str(stats.df)])
text(0.6,ylimits(2)*.85,['F= ' num2str((t{2,5}))])      
       
       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plotting - here is extData Fig 6C - xpos 0-200cm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot panel 6C : x-position 0-200cm of maze stem for left and right choice trials seperately         
clearvars -except groupOrder nTrialsSubselect_xTask nTrialsTotal_xTask ...
           yVel_xMouse xPosL_xMouse xPosR_xMouse viewAngL_xMouse viewAngR_xMouse excTravel_xMouse ...
           vaDev_xMouse decodeChoiceFromVA_xMouse decodeChoiceFromXpos_xMouse ...
           predChoiceXpos_xTask predChoiceXva_xTask ...
           concatLog_OFF_aoe concatLog_OFF_nd concatLog_OFF_pc masterLog
                                           
%%%%%%%%%%%%%%%%%%%%%%%%%%
% update this section for plotting of measured vars above
dataArray   = []; data = []; xIdx = []; 
dataArray   = xPosL_xMouse;
data        = [xPosL_xMouse{1}'; ...
               xPosL_xMouse{2}'; ...
               xPosL_xMouse{3}'];
xIdx        = [ones(size(xPosL_xMouse{1}))';   ...
               ones(size(xPosL_xMouse{2}))'*2; ...
               ones(size(xPosL_xMouse{3}))'*3];
ylimits     = [-0.5 0.5];
yaxName     = 'x-position (cm, 0-200cm)';
titlestr    = 'left choice';
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
h = errorbar(nGroup,nanmean(dataArray{nGroup}),nanstd(dataArray{nGroup})./sqrt(numel(dataArray{nGroup})-1), ...
    '.-','linewidth',4,'color','k','markersize',0.1); h.CapSize = 16;

nGroup = 2; 
h = errorbar(nGroup,nanmean(dataArray{nGroup}),nanstd(dataArray{nGroup})./sqrt(numel(dataArray{nGroup})-1), ...
    '.-','linewidth',4,'color','m','markersize',0.1); h.CapSize = 16;

nGroup = 3;
h = errorbar(nGroup,nanmean(dataArray{nGroup}),nanstd(dataArray{nGroup})./sqrt(numel(dataArray{nGroup})-1), ...
    '.-','linewidth',4,'color','c','markersize',0.1); h.CapSize = 16;

set(gca,'xtick',1:3,'xticklabel',groupNames)
rotateXLabels(gca,45)
ylabel(yaxName,'fontsize',12)
xlim(xlimits)
ylim(ylimits)
plot(xlimits(1):.5:xlimits(2),zeros(1,numel(xlimits(1):.5:xlimits(2))),'k-');
title(titlestr)
box off
set(gca,'TickDir','out');

xIdx = categorical(xIdx);
[p,t,stats] = anova1(data,xIdx,'off');
text(0.6,ylimits(2)*.95,['group p= ' num2str(p)])
text(0.6,ylimits(2)*.9,['df= ' num2str(stats.df)])
text(0.6,ylimits(2)*.85,['F= ' num2str((t{2,5}))])      

%%%%%%%%%%%%%%%%%%%%%%%%%%
% update this section for plotting of measured vars above
dataArray   = []; data = []; xIdx = []; 
dataArray   = xPosR_xMouse;
data        = [xPosR_xMouse{1}'; ...
               xPosR_xMouse{2}'; ...
               xPosR_xMouse{3}'];
xIdx        = [ones(size(xPosR_xMouse{1}))';   ...
               ones(size(xPosR_xMouse{2}))'*2; ...
               ones(size(xPosR_xMouse{3}))'*3];
ylimits     = [-0.5 0.5];
% yaxName     = 'x-position (cm, 0-200cm)';
titlestr    = 'right choice';
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
h = errorbar(nGroup,nanmean(dataArray{nGroup}),nanstd(dataArray{nGroup})./sqrt(numel(dataArray{nGroup})-1), ...
    '.-','linewidth',4,'color','k','markersize',0.1); h.CapSize = 16;

nGroup = 2; 
h = errorbar(nGroup,nanmean(dataArray{nGroup}),nanstd(dataArray{nGroup})./sqrt(numel(dataArray{nGroup})-1), ...
    '.-','linewidth',4,'color','m','markersize',0.1); h.CapSize = 16;

nGroup = 3;
h = errorbar(nGroup,nanmean(dataArray{nGroup}),nanstd(dataArray{nGroup})./sqrt(numel(dataArray{nGroup})-1), ...
    '.-','linewidth',4,'color','c','markersize',0.1); h.CapSize = 16;

set(gca,'xtick',1:3,'xticklabel',groupNames)
rotateXLabels(gca,45)
% ylabel(yaxName,'fontsize',12)
xlim(xlimits)
ylim(ylimits)
plot(xlimits(1):.5:xlimits(2),zeros(1,numel(xlimits(1):.5:xlimits(2))),'k-');
title(titlestr)
box off
set(gca,'TickDir','out');

xIdx = categorical(xIdx);
[p,t,stats] = anova1(data,xIdx,'off');
text(0.6,ylimits(2)*.95,['group p= ' num2str(p)])
text(0.6,ylimits(2)*.9,['df= ' num2str(stats.df)])
text(0.6,ylimits(2)*.85,['F= ' num2str((t{2,5}))])   


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plotting - here is ExtData Fig 6D - view angle 0-200cm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot panel 6D : delta view angle 0-200cm of maze stem for left and right choice trials seperately         
clearvars -except groupOrder nTrialsSubselect_xTask nTrialsTotal_xTask ...
           yVel_xMouse xPosL_xMouse xPosR_xMouse viewAngL_xMouse viewAngR_xMouse excTravel_xMouse ...
           vaDev_xMouse decodeChoiceFromVA_xMouse decodeChoiceFromXpos_xMouse ...
           predChoiceXpos_xTask predChoiceXva_xTask ...
           concatLog_OFF_aoe concatLog_OFF_nd concatLog_OFF_pc masterLog
                                           
%%%%%%%%%%%%%%%%%%%%%%%%%%
% update this section for plotting of measured vars above
dataArray   = []; data = []; xIdx = []; 
dataArray   = viewAngL_xMouse;
data        = [viewAngL_xMouse{1}'; ...
               viewAngL_xMouse{2}'; ...
               viewAngL_xMouse{3}'];
xIdx        = [ones(size(viewAngL_xMouse{1}))';   ...
               ones(size(viewAngL_xMouse{2}))'*2; ...
               ones(size(viewAngL_xMouse{3}))'*3];
ylimits     = [-30 30];
yaxName     = 'view angle (deg, 0-200cm)';
titlestr    = 'left choice';
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
h = errorbar(nGroup,nanmean(dataArray{nGroup}),nanstd(dataArray{nGroup})./sqrt(numel(dataArray{nGroup})-1), ...
    '.-','linewidth',4,'color','k','markersize',0.1); h.CapSize = 16;

nGroup = 2; 
h = errorbar(nGroup,nanmean(dataArray{nGroup}),nanstd(dataArray{nGroup})./sqrt(numel(dataArray{nGroup})-1), ...
    '.-','linewidth',4,'color','m','markersize',0.1); h.CapSize = 16;

nGroup = 3;
h = errorbar(nGroup,nanmean(dataArray{nGroup}),nanstd(dataArray{nGroup})./sqrt(numel(dataArray{nGroup})-1), ...
    '.-','linewidth',4,'color','c','markersize',0.1); h.CapSize = 16;

set(gca,'xtick',1:3,'xticklabel',groupNames)
rotateXLabels(gca,45)
ylabel(yaxName,'fontsize',12)
xlim(xlimits)
ylim(ylimits)
plot(xlimits(1):.5:xlimits(2),zeros(1,numel(xlimits(1):.5:xlimits(2))),'k-');
title(titlestr)
box off
set(gca,'TickDir','out');

xIdx = categorical(xIdx);
[p,t,stats] = anova1(data,xIdx,'off');
text(0.6,ylimits(2)*.95,['group p= ' num2str(p)])
text(0.6,ylimits(2)*.9,['df= ' num2str(stats.df)])
text(0.6,ylimits(2)*.85,['F= ' num2str((t{2,5}))])      

%%%%%%%%%%%%%%%%%%%%%%%%%%
% update this section for plotting of measured vars above
dataArray   = []; data = []; xIdx = []; 
dataArray   = viewAngR_xMouse;
data        = [viewAngR_xMouse{1}'; ...
               viewAngR_xMouse{2}'; ...
               viewAngR_xMouse{3}'];
xIdx        = [ones(size(viewAngR_xMouse{1}))';   ...
               ones(size(viewAngR_xMouse{2}))'*2; ...
               ones(size(viewAngR_xMouse{3}))'*3];
ylimits     = [-30 30];
% yaxName     = 'view angle (deg, 0-200cm)';
titlestr    = 'right choice';
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
h = errorbar(nGroup,nanmean(dataArray{nGroup}),nanstd(dataArray{nGroup})./sqrt(numel(dataArray{nGroup})-1), ...
    '.-','linewidth',4,'color','k','markersize',0.1); h.CapSize = 16;

nGroup = 2; 
h = errorbar(nGroup,nanmean(dataArray{nGroup}),nanstd(dataArray{nGroup})./sqrt(numel(dataArray{nGroup})-1), ...
    '.-','linewidth',4,'color','m','markersize',0.1); h.CapSize = 16;

nGroup = 3;
h = errorbar(nGroup,nanmean(dataArray{nGroup}),nanstd(dataArray{nGroup})./sqrt(numel(dataArray{nGroup})-1), ...
    '.-','linewidth',4,'color','c','markersize',0.1); h.CapSize = 16;

set(gca,'xtick',1:3,'xticklabel',groupNames)
rotateXLabels(gca,45)
% ylabel(yaxName,'fontsize',12)
xlim(xlimits)
ylim(ylimits)
plot(xlimits(1):.5:xlimits(2),zeros(1,numel(xlimits(1):.5:xlimits(2))),'k-');
title(titlestr)
box off
set(gca,'TickDir','out');

xIdx = categorical(xIdx);
[p,t,stats] = anova1(data,xIdx,'off');
text(0.6,ylimits(2)*.95,['group p= ' num2str(p)])
text(0.6,ylimits(2)*.9,['df= ' num2str(stats.df)])
text(0.6,ylimits(2)*.85,['F= ' num2str((t{2,5}))])   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plotting  - here is panel 6E trials with excess travel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot panel 6E :  % trials excess travel 0-200cm of maze stem          
clearvars -except groupOrder nTrialsSubselect_xTask nTrialsTotal_xTask ...
           yVel_xMouse xPosL_xMouse xPosR_xMouse viewAngL_xMouse viewAngR_xMouse excTravel_xMouse ...
           vaDev_xMouse decodeChoiceFromVA_xMouse decodeChoiceFromXpos_xMouse ...
           predChoiceXpos_xTask predChoiceXva_xTask ...
           concatLog_OFF_aoe concatLog_OFF_nd concatLog_OFF_pc masterLog
                                           
%%%%%%%%%%%%%%%%%%%%%%%%%%
% update this section for plotting of measured vars above
dataArray   = []; data = []; xIdx = []; 
dataArray   = excTravel_xMouse;
data        = [excTravel_xMouse{1}'; ...
               excTravel_xMouse{2}'; ...
               excTravel_xMouse{3}'];
xIdx        = [ones(size(excTravel_xMouse{1}))';   ...
               ones(size(excTravel_xMouse{2}))'*2; ...
               ones(size(excTravel_xMouse{3}))'*3];
ylimits     = [0 25];
yaxName     = 'trials w/ excess travel (%, >330cm)';
%%%%%%%%%%%%%%%%%%%%%%%%%

groupNames  = {'aoe';'ctrl#1';'ctrl#2'};
catMarker   = {'x','x', 'x'};
colorIdx    = {[0.8 0.8 0.8],[0.99 0.7 0.99], [0.7 0.99 0.99]};
xlimits     = [0.5 3.5];

figure
hold on
plotSpread(data,'categoryIdx',xIdx,'distributionIdx',xIdx,...
    'categoryColors',colorIdx,'categoryMarkers',catMarker);

nGroup = 1;
h = errorbar(nGroup,nanmean(dataArray{nGroup}),nanstd(dataArray{nGroup})./sqrt(numel(dataArray{nGroup})-1), ...
    '.-','linewidth',4,'color','k','markersize',0.1); h.CapSize = 16;

nGroup = 2; 
h = errorbar(nGroup,nanmean(dataArray{nGroup}),nanstd(dataArray{nGroup})./sqrt(numel(dataArray{nGroup})-1), ...
    '.-','linewidth',4,'color','m','markersize',0.1); h.CapSize = 16;

nGroup = 3;
h = errorbar(nGroup,nanmean(dataArray{nGroup}),nanstd(dataArray{nGroup})./sqrt(numel(dataArray{nGroup})-1), ...
    '.-','linewidth',4,'color','c','markersize',0.1); h.CapSize = 16;

set(gca,'xtick',1:3,'xticklabel',groupNames)
rotateXLabels(gca,45)
ylabel(yaxName,'fontsize',12)
xlim(xlimits)
ylim(ylimits)
% plot(xlimits(1):.5:xlimits(2),zeros(1,numel(xlimits(1):.5:xlimits(2))),'k-');
box off
set(gca,'TickDir','out');

xIdx = categorical(xIdx);
[p,t,stats] = anova1(data,xIdx,'off');
text(0.6,ylimits(2)*.95,['group p= ' num2str(p)])
text(0.6,ylimits(2)*.9,['df= ' num2str(stats.df)])
text(0.6,ylimits(2)*.85,['F= ' num2str((t{2,5}))])   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plotting  - here is panel 6F per-trial standard deviation in view angle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot panel 6F :            
clearvars -except groupOrder nTrialsSubselect_xTask nTrialsTotal_xTask ...
           yVel_xMouse xPosL_xMouse xPosR_xMouse viewAngL_xMouse viewAngR_xMouse excTravel_xMouse ...
           vaDev_xMouse decodeChoiceFromVA_xMouse decodeChoiceFromXpos_xMouse ...
           predChoiceXpos_xTask predChoiceXva_xTask ...
           concatLog_OFF_aoe concatLog_OFF_nd concatLog_OFF_pc masterLog
                                           
%%%%%%%%%%%%%%%%%%%%%%%%%%
% update this section for plotting of measured vars above
dataArray   = []; data = []; xIdx = []; 
dataArray   = vaDev_xMouse;
data        = [vaDev_xMouse{1}'; ...
               vaDev_xMouse{2}'; ...
               vaDev_xMouse{3}'];
xIdx        = [ones(size(vaDev_xMouse{1}))';   ...
               ones(size(vaDev_xMouse{2}))'*2; ...
               ones(size(vaDev_xMouse{3}))'*3];
ylimits     = [0 180];
yaxName     = 'per-trial view angle STD (degrees)';
%%%%%%%%%%%%%%%%%%%%%%%%%

groupNames  = {'aoe';'ctrl#1';'ctrl#2'};
catMarker   = {'x','x', 'x'};
colorIdx    = {[0.8 0.8 0.8],[0.99 0.7 0.99], [0.7 0.99 0.99]};
xlimits     = [0.5 3.5];

figure
hold on
plotSpread(data,'categoryIdx',xIdx,'distributionIdx',xIdx,...
    'categoryColors',colorIdx,'categoryMarkers',catMarker);

nGroup = 1;
h = errorbar(nGroup,nanmean(dataArray{nGroup}),nanstd(dataArray{nGroup})./sqrt(numel(dataArray{nGroup})-1), ...
    '.-','linewidth',4,'color','k','markersize',0.1); h.CapSize = 16;

nGroup = 2; 
h = errorbar(nGroup,nanmean(dataArray{nGroup}),nanstd(dataArray{nGroup})./sqrt(numel(dataArray{nGroup})-1), ...
    '.-','linewidth',4,'color','m','markersize',0.1); h.CapSize = 16;

nGroup = 3;
h = errorbar(nGroup,nanmean(dataArray{nGroup}),nanstd(dataArray{nGroup})./sqrt(numel(dataArray{nGroup})-1), ...
    '.-','linewidth',4,'color','c','markersize',0.1); h.CapSize = 16;

set(gca,'xtick',1:3,'xticklabel',groupNames)
rotateXLabels(gca,45)
ylabel(yaxName,'fontsize',12)
xlim(xlimits)
ylim(ylimits)
% plot(xlimits(1):.5:xlimits(2),zeros(1,numel(xlimits(1):.5:xlimits(2))),'k-');
box off
set(gca,'TickDir','out');

xIdx = categorical(xIdx);
[p,t,stats] = anova1(data,xIdx,'off');
text(0.6,ylimits(2)*.95,['group p= ' num2str(p)])
text(0.6,ylimits(2)*.9,['df= ' num2str(stats.df)])
text(0.6,ylimits(2)*.85,['F= ' num2str((t{2,5}))])   


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plotting  - here is panel 6G-H decode choice from x-pos/view angle trajectory (train/test WITHIN mouse)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars -except groupOrder nTrialsSubselect_xTask nTrialsTotal_xTask ...
           yVel_xMouse xPosL_xMouse xPosR_xMouse viewAngL_xMouse viewAngR_xMouse excTravel_xMouse ...
           vaDev_xMouse decodeChoiceFromVA_xMouse decodeChoiceFromXpos_xMouse ...
           predChoiceXpos_xTask predChoiceXva_xTask ...
           concatLog_OFF_aoe concatLog_OFF_nd concatLog_OFF_pc masterLog
       
% plot panel 6G :  decode choice from x-pos traject WITHIN mouse        
plotDecodeChoiceXTraject(decodeChoiceFromXpos_xMouse,groupOrder)
       
% plot panel 6H :  decode choice from view angle traject WITHIN mouse        
plotDecodeChoiceXTraject(decodeChoiceFromVA_xMouse,groupOrder)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plotting  - here is panel 6I decode choice from x-pos trajecotry (train/test between mice)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot panel 6I :  decode choice from x-pos traject train/test between mice mouse        
clearvars -except groupOrder nTrialsSubselect_xTask nTrialsTotal_xTask ...
           yVel_xMouse xPosL_xMouse xPosR_xMouse viewAngL_xMouse viewAngR_xMouse excTravel_xMouse ...
           vaDev_xMouse decodeChoiceFromVA_xMouse decodeChoiceFromXpos_xMouse ...
           predChoiceXpos_xTask predChoiceXva_xTask ...
           concatLog_OFF_aoe concatLog_OFF_nd concatLog_OFF_pc masterLog
       
% individual mouse accuracy plotting
figure
plotDecodeChoiceXTraject_individual(predChoiceXpos_xTask)

% group mouse accuracy plotting
temp = [];
temp = predChoiceXpos_xTask;
trainAOE    = []; trainND    = []; trainPC    = [];
trainAOEsem = []; trainNDsem = []; trainPCsem = [];

trainAOE    = [temp.trainTOW.xMouseAccAvgTOW; temp.trainTOW.xMouseAccAvgMG; temp.trainTOW.xMouseAccAvgPC]; trainAOE = trainAOE*100;
trainND     = [temp.trainMG.xMouseAccAvgTOW;  temp.trainMG.xMouseAccAvgMG;  temp.trainMG.xMouseAccAvgPC];  trainND  = trainND*100;
trainPC     = [temp.trainPC.xMouseAccAvgTOW;  temp.trainPC.xMouseAccAvgMG;  temp.trainPC.xMouseAccAvgPC];  trainPC  = trainPC*100;

trainAOEsem = [temp.trainTOW.xMouseAccSemTOW; temp.trainTOW.xMouseAccSemMG; temp.trainTOW.xMouseAccSemPC]; trainAOEsem = trainAOEsem*100;
trainNDsem  = [temp.trainMG.xMouseAccSemTOW;  temp.trainMG.xMouseAccSemMG;  temp.trainMG.xMouseAccSemPC];  trainNDsem  = trainNDsem*100;
trainPCsem  = [temp.trainPC.xMouseAccSemTOW;  temp.trainPC.xMouseAccSemMG;  temp.trainPC.xMouseAccSemPC];  trainPCsem  = trainPCsem*100;

binIDs = [];
binIDs = temp.binIDs;
gapFill = []; gapFill = diff(binIDs(1:2))/2;
plotOffsets = [-5 0 5];
tempBins = (binIDs+gapFill)- plotOffsets(1);
figure
subplot(1,3,1)
errorbar(trainAOE(1,:), tempBins+plotOffsets(1), -trainAOEsem(1,:), trainAOEsem(1,:), '.','horizontal', 'color', 'k', 'markersize',4, 'markerfacecolor', 'k'); hold on
errorbar(trainAOE(2,:), tempBins+plotOffsets(2), -trainAOEsem(2,:), trainAOEsem(2,:), '.','horizontal', 'color', 'm', 'markersize',4, 'markerfacecolor', 'm'); hold on
errorbar(trainAOE(3,:), tempBins+plotOffsets(3), -trainAOEsem(3,:), trainAOEsem(3,:), '.','horizontal', 'color', 'c', 'markersize',4, 'markerfacecolor', 'c'); hold on

title('train on AoE')
box off
set(gca,'TickDir','out');
xlim([20 100])
xticks([20 40 60 80 100])
ylabel('y position (cm)')

subplot(1,3,2)
errorbar(trainND(1,:), tempBins+plotOffsets(1), -trainNDsem(1,:), trainNDsem(1,:), '.','horizontal','color', 'k', 'markersize',4, 'markerfacecolor', 'k'); hold on
errorbar(trainND(2,:), tempBins+plotOffsets(2), -trainNDsem(2,:), trainNDsem(2,:), '.','horizontal','color', 'm', 'markersize',4, 'markerfacecolor', 'm'); hold on
errorbar(trainND(3,:), tempBins+plotOffsets(3), -trainNDsem(3,:), trainNDsem(3,:), '.','horizontal','color', 'c', 'markersize',4, 'markerfacecolor', 'c'); hold on
% plot(trainMG(1,:),tempBins,'k-'); hold on
% plot(trainMG(2,:),tempBins,'m-'); hold on
% plot(trainMG(3,:),tempBins,'c-'); hold on
title('train on Ctrl#1')
box off
set(gca,'TickDir','out');
xlim([20 100])
xticks([20 40 60 80 100])
xlabel('% decode accuracy')

subplot(1,3,3)
errorbar(trainPC(1,:), tempBins+plotOffsets(1), -trainPCsem(1,:), trainPCsem(1,:), '.','horizontal','color', 'k', 'markersize',4, 'markerfacecolor', 'k'); hold on
errorbar(trainPC(2,:), tempBins+plotOffsets(2), -trainPCsem(2,:), trainPCsem(2,:), '.','horizontal','color', 'm', 'markersize',4, 'markerfacecolor', 'm'); hold on
errorbar(trainPC(3,:), tempBins+plotOffsets(3), -trainPCsem(3,:), trainPCsem(3,:), '.','horizontal','color', 'c', 'markersize',4, 'markerfacecolor', 'c'); hold on
title('train on Ctrl#2')
box off
set(gca,'TickDir','out');
xlim([20 100])
xticks([20 40 60 80 100])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plotting  - here is panel 6J decode choice from view angle trajecotry (train/test between mice)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot panel 6J :  decode choice from view angle traject within mouse        
clearvars -except groupOrder nTrialsSubselect_xTask nTrialsTotal_xTask ...
           yVel_xMouse xPosL_xMouse xPosR_xMouse viewAngL_xMouse viewAngR_xMouse excTravel_xMouse ...
           vaDev_xMouse decodeChoiceFromVA_xMouse decodeChoiceFromXpos_xMouse ...
           predChoiceXpos_xTask predChoiceXva_xTask ...
           concatLog_OFF_aoe concatLog_OFF_nd concatLog_OFF_pc masterLog
       
% individual mouse accuracy plotting
figure
plotDecodeChoiceXTraject_individual(predChoiceXva_xTask) % update this

% group mouse accuracy plotting
temp = [];
temp = predChoiceXva_xTask; % update this
trainAOE    = []; trainND    = []; trainPC    = [];
trainAOEsem = []; trainNDsem = []; trainPCsem = [];

trainAOE    = [temp.trainTOW.xMouseAccAvgTOW; temp.trainTOW.xMouseAccAvgMG; temp.trainTOW.xMouseAccAvgPC]; trainAOE = trainAOE*100;
trainND     = [temp.trainMG.xMouseAccAvgTOW;  temp.trainMG.xMouseAccAvgMG;  temp.trainMG.xMouseAccAvgPC];  trainND  = trainND*100;
trainPC     = [temp.trainPC.xMouseAccAvgTOW;  temp.trainPC.xMouseAccAvgMG;  temp.trainPC.xMouseAccAvgPC];  trainPC  = trainPC*100;

trainAOEsem = [temp.trainTOW.xMouseAccSemTOW; temp.trainTOW.xMouseAccSemMG; temp.trainTOW.xMouseAccSemPC]; trainAOEsem = trainAOEsem*100;
trainNDsem  = [temp.trainMG.xMouseAccSemTOW;  temp.trainMG.xMouseAccSemMG;  temp.trainMG.xMouseAccSemPC];  trainNDsem  = trainNDsem*100;
trainPCsem  = [temp.trainPC.xMouseAccSemTOW;  temp.trainPC.xMouseAccSemMG;  temp.trainPC.xMouseAccSemPC];  trainPCsem  = trainPCsem*100;

binIDs = [];
binIDs = temp.binIDs;
gapFill = []; gapFill = diff(binIDs(1:2))/2;
plotOffsets = [-5 0 5];
tempBins = (binIDs+gapFill)- plotOffsets(1);
figure
subplot(1,3,1)
errorbar(trainAOE(1,:), tempBins+plotOffsets(1), -trainAOEsem(1,:), trainAOEsem(1,:), '.','horizontal', 'color', 'k', 'markersize',4, 'markerfacecolor', 'k'); hold on
errorbar(trainAOE(2,:), tempBins+plotOffsets(2), -trainAOEsem(2,:), trainAOEsem(2,:), '.','horizontal', 'color', 'm', 'markersize',4, 'markerfacecolor', 'm'); hold on
errorbar(trainAOE(3,:), tempBins+plotOffsets(3), -trainAOEsem(3,:), trainAOEsem(3,:), '.','horizontal', 'color', 'c', 'markersize',4, 'markerfacecolor', 'c'); hold on

title('train on AoE')
box off
set(gca,'TickDir','out');
xlim([20 100])
xticks([20 40 60 80 100])
ylabel('y position (cm)')

subplot(1,3,2)
errorbar(trainND(1,:), tempBins+plotOffsets(1), -trainNDsem(1,:), trainNDsem(1,:), '.','horizontal','color', 'k', 'markersize',4, 'markerfacecolor', 'k'); hold on
errorbar(trainND(2,:), tempBins+plotOffsets(2), -trainNDsem(2,:), trainNDsem(2,:), '.','horizontal','color', 'm', 'markersize',4, 'markerfacecolor', 'm'); hold on
errorbar(trainND(3,:), tempBins+plotOffsets(3), -trainNDsem(3,:), trainNDsem(3,:), '.','horizontal','color', 'c', 'markersize',4, 'markerfacecolor', 'c'); hold on
% plot(trainMG(1,:),tempBins,'k-'); hold on
% plot(trainMG(2,:),tempBins,'m-'); hold on
% plot(trainMG(3,:),tempBins,'c-'); hold on
title('train on Ctrl#1')
box off
set(gca,'TickDir','out');
xlim([20 100])
xticks([20 40 60 80 100])
xlabel('% decode accuracy')

subplot(1,3,3)
errorbar(trainPC(1,:), tempBins+plotOffsets(1), -trainPCsem(1,:), trainPCsem(1,:), '.','horizontal','color', 'k', 'markersize',4, 'markerfacecolor', 'k'); hold on
errorbar(trainPC(2,:), tempBins+plotOffsets(2), -trainPCsem(2,:), trainPCsem(2,:), '.','horizontal','color', 'm', 'markersize',4, 'markerfacecolor', 'm'); hold on
errorbar(trainPC(3,:), tempBins+plotOffsets(3), -trainPCsem(3,:), trainPCsem(3,:), '.','horizontal','color', 'c', 'markersize',4, 'markerfacecolor', 'c'); hold on
title('train on Ctrl#2')
box off
set(gca,'TickDir','out');
xlim([20 100])
xticks([20 40 60 80 100])

%% clear all vars except source data for plotting
clearvars -except groupOrder nTrialsSubselect_xTask nTrialsTotal_xTask ...
           yVel_xMouse xPosL_xMouse xPosR_xMouse viewAngL_xMouse viewAngR_xMouse excTravel_xMouse ...
           vaDev_xMouse decodeChoiceFromVA_xMouse decodeChoiceFromXpos_xMouse ...
           predChoiceXpos_xTask predChoiceXva_xTask ...
           concatLog_OFF_aoe concatLog_OFF_nd concatLog_OFF_pc masterLog

