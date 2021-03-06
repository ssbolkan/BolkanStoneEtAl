%% Script to generate ExtData Fig 4B-G data plots
% related to Fig 1 virtual corridor - extra motor plots and on-off stats
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load virtual corridor behavior log files 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
cd(globalParams.dataPath)
load('virtualCorridor.mat')

%% analyzes average velocity, x-position, view angle, distance, excess travel, view angle deviation
groupOrder = {'D2_vc', 'D1_vc', 'Ctrl_vc'};
counter = 1;
for groupNum = 1:numel(groupOrder)
    clear lg
    lg = eval(['lgDMS' groupOrder{groupNum}]);
    
    %  log pre-process and data selection for velocity, x-position and view angle    
    % invert x-position and view angle on leftHemisphere inhibition, which
    % was always on odd numbered dates. positive values reflect
    % contralateral movement relative to hemisphere
    leftHemiDay = [];
    leftHemiDay = mod(lg.date,2);  
    for nTrial = 1:numel(leftHemiDay)
        if leftHemiDay(nTrial)
            % pos N x 3column matrix of x-pos,y-pos,view angle --default
            lg.pos{1,nTrial}(:,1)  = -lg.pos{1,nTrial}(:,1);
            lg.pos{1,nTrial}(:,3)  = -lg.pos{1,nTrial}(:,3);
           
            % pos2 N x 3column matrix of view angle,y-pos,x-pos --swapped order for convenience
            lg.pos2{1,nTrial}(:,1) = -lg.pos2{1,nTrial}(:,1); 
            lg.pos2{1,nTrial}(:,3) = -lg.pos2{1,nTrial}(:,3);
        else
        end
    end
    
    % remove trials with excess travel >10% of a 300cm stem
    removeInd = []; lgClean  = [];
    removeInd = lg.excessTravel > 1.1;    
    lgClean  = structfun(@(x) x(~removeInd),lg,'UniformOutput',false);
    
    % remove mice with fewer than 150 trials after excess travel removal
    tempcount      = 1;
    trialsPerMouse = [];    
    mouseRemove    = []; 
    mouseNums      = [];
    mouseNums      = unique(lgClean.mouseID);
    for nMouse = 1:numel(mouseNums)
        trialsPerMouse(nMouse) = length(lgClean.mouseID(lgClean.mouseID == mouseNums(nMouse)));
        if trialsPerMouse(nMouse)<150
            mouseRemove(tempcount) = mouseNums(nMouse);
            tempcount = tempcount + 1;
        end
    end
    for nMouse = 1:numel(mouseNums)
        if ~ismember(mouseNums(nMouse),mouseRemove)
            continue
        else
            removeMouse = [];
            removeMouse = lgClean.mouseID == mouseNums(nMouse);
            lgClean  = structfun(@(x) x(~removeMouse),lgClean,'UniformOutput',false);
        end
    end
    
    % create seperate logCleanOFF and logCleanON, which inclues only on "cue" or 0-200cm
    removeMouse = [];
    removeMouse = lgClean.laserOFF ==1;
    lgCleanOFF = structfun(@(x) x(removeMouse),lgClean,'UniformOutput',false);
    removeMouse = [];
    removeMouse = lgClean.laserCue == 1;
    lgCleanON  = structfun(@(x) x(removeMouse),lgClean,'UniformOutput',false);
    
    nTrialsClean_OFF{counter} = numel(lgCleanOFF.distance);
    nTrialsClean_ON{counter}  = numel(lgCleanON.distance);
    
    % panel B - calculation of avg Yvelocity 11-211cm which is laser region for
    % cue light on log and light off log
    mazeSection = []; mazeSection = [11 211];
    yVelOFF{counter}(:)   = xMouseSpeedXY(lgCleanOFF, mazeSection, 2);
    yVelON{counter}(:)    = xMouseSpeedXY(lgCleanON, mazeSection, 2);
    
    
    % panel 1C - calculation of avg x-position 11-211cm which is laser region for
    % cue light on log and light off log
    avgXposON  = [];
    avgXposOFF = [];
    avgXposON  = sampleViewAngleVsY_average(lgCleanON.pos2,[11 211],200);
    avgXposOFF = sampleViewAngleVsY_average(lgCleanOFF.pos2,[11 211],200);
    
    mouseNums       = [];
    mouseNums       = unique(lgClean.mouseID);
    for nMouse = 1:numel(mouseNums)
        xPosOFF{counter}(:,nMouse) = nanmean(avgXposOFF(:,lgCleanOFF.mouseID == mouseNums(nMouse)),2);
        xPosON{counter}(:,nMouse)  = nanmean(avgXposON(:,lgCleanON.mouseID == mouseNums(nMouse)),2);
    end
    
    % panel D - calculation of view angle 11-211cm which is laser region for
    % cue light on log and light off log
    avgViewAngleON  = []; avgViewAngleOFF = [];
    avgViewAngleON  = sampleViewAngleVsY_average(lgCleanON.pos,[11 211],200);
    avgViewAngleOFF = sampleViewAngleVsY_average(lgCleanOFF.pos,[11 211],200);
    
    mouseNums       = [];
    mouseNums       = unique(lgClean.mouseID);
    for nMouse = 1:numel(mouseNums)
        viewAngOFF{counter}(:,nMouse) = nanmean(avgViewAngleOFF(:,lgCleanOFF.mouseID == mouseNums(nMouse)),2);
        viewAngON{counter}(:,nMouse)  = nanmean(avgViewAngleON(:,lgCleanON.mouseID == mouseNums(nMouse)),2);
    end
         
    
    % panel E - calculation of average total distance per trial - 
    % no trial or mouse subselection - on "cue" light on (11-211cm) trials and 
    % off trials 
    clear lg
    lg = eval(['lgDMS' groupOrder{groupNum}]);
    
    removeMouse = [];
    removeMouse = lg.laserOFF == 1;
    lgOFF      = structfun(@(x) x(removeMouse),lg,'UniformOutput',false);
    removeMouse = [];
    removeMouse = lg.laserCue == 1;
    lgON       = structfun(@(x) x(removeMouse),lg,'UniformOutput',false);
    removeMouse = [];
    removeMouse = lg.laserWhole == 1;
    lg          = structfun(@(x) x(~removeMouse),lg,'UniformOutput',false);
    
    nTrialsALL_OFF{counter} = numel(lgOFF.distance);
    nTrialsALL_ON{counter}  = numel(lgON.distance);
    
    posBins = []; posBins = 11:5:311;  
    vsdOFF  = []; vsdON   = [];
    
    vsdOFF  = xMouseViewAngSD(lgOFF, posBins);
    vsdON   = xMouseViewAngSD(lgON, posBins);
    
    vaDevOFF{counter}(:)  = vsdOFF.viewAngSD;    
    vaDevON{counter}(:)   = vsdON.viewAngSD;
    
    
    % panel F - calculation of average total distance per trial - 
    % no trial or mouse subselection - on "cue" light on (11-211cm) trials and 
    % off trials 
    clear distance mouseNums distOFF distON
    distance = cellfun(@(x) x(end), lg.distance);
    
    mouseNums = [];
    mouseNums = unique(lg.mouseID);
    for nMouse = 1:numel(mouseNums)
        distOFF(nMouse) = nanmean(distance(lg.mouseID == mouseNums(nMouse)...
            & lg.laserOFF==1));
        distON(nMouse)  = nanmean(distance(lg.mouseID == mouseNums(nMouse) ...
            & lg.laserCue==1));
    end
    
    distanceOFF{counter} = distOFF;
    distanceON{counter}  = distON;
    
    
    % panel G - calculation of % excess travel trials - 
    % no trial or mouse subselection - on "cue" light on (11-211cm) trials and 
    % off trials 
   
    mouseNums = [];
    mouseNums = unique(lg.mouseID);
    for nMouse = 1:numel(mouseNums)
        etOFF(nMouse) = (sum(lg.excessTravel(lg.mouseID == mouseNums(nMouse)...
                         & lg.laserOFF == 1) > 1.1)/...
                         numel(lg.excessTravel(lg.mouseID == mouseNums(nMouse)...
                         & lg.laserOFF ==1)))*100;
        
        etON(nMouse)  = (sum(lg.excessTravel(lg.mouseID == mouseNums(nMouse)...
                         & lg.laserCue == 1) > 1.1)/...
                         numel(lg.excessTravel(lg.mouseID == mouseNums(nMouse)...
                         & lg.laserCue ==1)))*100;
    end
    
    excessTravelOFF{counter} = etOFF;
    excessTravelON{counter}  = etON;  
    
    % count up for next log file
    counter = counter +1;
end
                      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plotting the data for extData Fig 4B-G
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot panel 1B : delta y-velocity 0-200cm of linear track           
clearvars -except yVelOFF yVelON viewAngOFF viewAngON ...
                  xPosOFF xPosON distanceOFF distanceON ...
                  excessTravelOFF excessTravelON ...
                  vaDevOFF vaDevON groupOrder ...
                  lgDMSCtrl_vc lgDMSD1_vc lgDMSD2_vc ...
                  nTrialsClean_OFF nTrialsClean_ON ...
                  nTrialsALL_OFF nTrialsALL_ON
                                           
%%%%%%%%%%%%%%%%%%%%%%%%%%
% update this section for plotting of measured vars above
dataArray   = []; data = []; xIdx = []; 
dataArray   = cellfun(@minus,yVelON,yVelOFF,'UniformOutput',false);
data        = [yVelON{1}'-yVelOFF{1}'; ...
               yVelON{2}'-yVelOFF{2}'; ...
               yVelON{3}'-yVelOFF{3}'];
xIdx        = [ones(size(yVelON{1}))';   ...
               ones(size(yVelON{2}))'*2; ...
               ones(size(yVelON{3}))'*3];
ylimits     = [-8 8];
yaxName     = 'delta y-velocity (cm/s, on-off)';
%%%%%%%%%%%%%%%%%%%%%%%%%

groupNames  = {'indirect';'direct';'no opsin'};
catMarker   = {'x','x', 'x'};
colorIdx    = {[0.8 0.8 0.8],[0.8 0.8 0.8], [0.8 0.8 0.8]};
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

xIdx = categorical(xIdx);
[p,t,stats] = anova1(data,xIdx,'off');
text(0.6,ylimits(2)*.9,['group p= ' num2str(p)])
text(0.6,ylimits(2)*.8,['df= ' num2str(stats.df)])
text(0.6,ylimits(2)*.7,['F= ' num2str((t{2,5}))])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot panel 1C : delta x-position 0-200cm of linear track           
clearvars -except yVelOFF yVelON viewAngOFF viewAngON ...
                  xPosOFF xPosON distanceOFF distanceON ...
                  excessTravelOFF excessTravelON ...
                  vaDevOFF vaDevON groupOrder ...
                  lgDMSCtrl_vc lgDMSD1_vc lgDMSD2_vc ...
                  nTrialsClean_OFF nTrialsClean_ON ...
                  nTrialsALL_OFF nTrialsALL_ON
                                           
%%%%%%%%%%%%%%%%%%%%%%%%%%
% update this section for plotting of measured vars above
dataArray   = []; data = []; xIdx = []; 
dataArray   = cellfun(@minus,xPosON,xPosOFF,'UniformOutput',false);
data        = [xPosON{1}'-xPosOFF{1}'; ...
               xPosON{2}'-xPosOFF{2}'; ...
               xPosON{3}'-xPosOFF{3}'];
xIdx        = [ones(size(xPosON{1}))';   ...
               ones(size(xPosON{2}))'*2; ...
               ones(size(xPosON{3}))'*3];
ylimits     = [-1 1];
yaxName     = 'delta x-position (cm, on-off)';
%%%%%%%%%%%%%%%%%%%%%%%%%

groupNames  = {'indirect';'direct';'no opsin'};
catMarker   = {'x','x', 'x'};
colorIdx    = {[0.8 0.8 0.8],[0.8 0.8 0.8], [0.8 0.8 0.8]};
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

xIdx = categorical(xIdx);
[p,t,stats] = anova1(data,xIdx,'off');
text(0.6,ylimits(2)*.9,['group p= ' num2str(p)])
text(0.6,ylimits(2)*.8,['df= ' num2str(stats.df)])
text(0.6,ylimits(2)*.7,['F= ' num2str((t{2,5}))])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot panel 1D : delta view angle 0-200cm of linear track           
clearvars -except yVelOFF yVelON viewAngOFF viewAngON ...
                  xPosOFF xPosON distanceOFF distanceON ...
                  excessTravelOFF excessTravelON ...
                  vaDevOFF vaDevON groupOrder ...
                  lgDMSCtrl_vc lgDMSD1_vc lgDMSD2_vc ...
                  nTrialsClean_OFF nTrialsClean_ON ...
                  nTrialsALL_OFF nTrialsALL_ON
                                           
%%%%%%%%%%%%%%%%%%%%%%%%%%
% update this section for plotting of measured vars above
dataArray   = []; data = []; xIdx = []; 
dataArray   = cellfun(@minus,viewAngON,viewAngOFF,'UniformOutput',false);
data        = [viewAngON{1}'-viewAngOFF{1}'; ...
               viewAngON{2}'-viewAngOFF{2}'; ...
               viewAngON{3}'-viewAngOFF{3}'];
xIdx        = [ones(size(viewAngON{1}))';   ...
               ones(size(viewAngON{2}))'*2; ...
               ones(size(viewAngON{3}))'*3];
ylimits     = [-8 8];
yaxName     = 'delta view angle (deg, on-off)';
%%%%%%%%%%%%%%%%%%%%%%%%%

groupNames  = {'indirect';'direct';'no opsin'};
catMarker   = {'x','x', 'x'};
colorIdx    = {[0.8 0.8 0.8],[0.8 0.8 0.8], [0.8 0.8 0.8]};
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

xIdx = categorical(xIdx);
[p,t,stats] = anova1(data,xIdx,'off');
text(0.6,ylimits(2)*.9,['group p= ' num2str(p)])
text(0.6,ylimits(2)*.8,['df= ' num2str(stats.df)])
text(0.6,ylimits(2)*.7,['F= ' num2str((t{2,5}))])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot panel 1E : delta per trial view angle STD 0-200cm of linear track           
clearvars -except yVelOFF yVelON viewAngOFF viewAngON ...
                  xPosOFF xPosON distanceOFF distanceON ...
                  excessTravelOFF excessTravelON ...
                  vaDevOFF vaDevON groupOrder ...
                  lgDMSCtrl_vc lgDMSD1_vc lgDMSD2_vc ...
                  nTrialsClean_OFF nTrialsClean_ON ...
                  nTrialsALL_OFF nTrialsALL_ON
                                           
%%%%%%%%%%%%%%%%%%%%%%%%%%
% update this section for plotting of measured vars above
dataArray   = []; data = []; xIdx = []; 
dataArray   = cellfun(@minus,vaDevON,vaDevOFF,'UniformOutput',false);
data        = [vaDevON{1}'-vaDevOFF{1}'; ...
               vaDevON{2}'-vaDevOFF{2}'; ...
               vaDevON{3}'-vaDevOFF{3}'];
xIdx        = [ones(size(vaDevON{1}))';   ...
               ones(size(vaDevON{2}))'*2; ...
               ones(size(vaDevON{3}))'*3];
ylimits     = [-48 48];
yaxName     = 'delta per-trial view angle STD (deg, on-off)';
%%%%%%%%%%%%%%%%%%%%%%%%%

groupNames  = {'indirect';'direct';'no opsin'};
catMarker   = {'x','x', 'x'};
colorIdx    = {[0.8 0.8 0.8],[0.8 0.8 0.8], [0.8 0.8 0.8]};
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

xIdx = categorical(xIdx);
[p,t,stats] = anova1(data,xIdx,'off');
text(0.6,ylimits(2)*.9,['group p= ' num2str(p)])
text(0.6,ylimits(2)*.8,['df= ' num2str(stats.df)])
text(0.6,ylimits(2)*.7,['F= ' num2str((t{2,5}))])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot panel 1F : delta per-trial distance of linear track           
clearvars -except yVelOFF yVelON viewAngOFF viewAngON ...
                  xPosOFF xPosON distanceOFF distanceON ...
                  excessTravelOFF excessTravelON ...
                  vaDevOFF vaDevON groupOrder ...
                  lgDMSCtrl_vc lgDMSD1_vc lgDMSD2_vc ...
                  nTrialsClean_OFF nTrialsClean_ON ...
                  nTrialsALL_OFF nTrialsALL_ON
                                           
%%%%%%%%%%%%%%%%%%%%%%%%%%
% update this section for plotting of measured vars above
dataArray   = []; data = []; xIdx = []; 
dataArray   = cellfun(@minus,distanceON,distanceOFF,'UniformOutput',false);
data        = [distanceON{1}'-distanceOFF{1}'; ...
               distanceON{2}'-distanceOFF{2}'; ...
               distanceON{3}'-distanceOFF{3}'];
xIdx        = [ones(size(distanceON{1}))';   ...
               ones(size(distanceON{2}))'*2; ...
               ones(size(distanceON{3}))'*3];
ylimits     = [-150 150];
yaxName     = 'delta distance (cm, on-off)';
%%%%%%%%%%%%%%%%%%%%%%%%%

groupNames  = {'indirect';'direct';'no opsin'};
catMarker   = {'x','x', 'x'};
colorIdx    = {[0.8 0.8 0.8],[0.8 0.8 0.8], [0.8 0.8 0.8]};
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

xIdx = categorical(xIdx);
[p,t,stats] = anova1(data,xIdx,'off');
text(0.6,ylimits(2)*.9,['group p= ' num2str(p)])
text(0.6,ylimits(2)*.8,['df= ' num2str(stats.df)])
text(0.6,ylimits(2)*.7,['F= ' num2str((t{2,5}))])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot panel 1G : delta trials with excess travel in linear track         
clearvars -except yVelOFF yVelON viewAngOFF viewAngON ...
                  xPosOFF xPosON distanceOFF distanceON ...
                  excessTravelOFF excessTravelON ...
                  vaDevOFF vaDevON groupOrder ...
                  lgDMSCtrl_vc lgDMSD1_vc lgDMSD2_vc ...
                  nTrialsClean_OFF nTrialsClean_ON ...
                  nTrialsALL_OFF nTrialsALL_ON
                                           
%%%%%%%%%%%%%%%%%%%%%%%%%%
% update this section for plotting of measured vars above
dataArray   = []; data = []; xIdx = []; 
dataArray   = cellfun(@minus,excessTravelON,excessTravelOFF,'UniformOutput',false);
data        = [excessTravelON{1}'-excessTravelOFF{1}'; ...
               excessTravelON{2}'-excessTravelOFF{2}'; ...
               excessTravelON{3}'-excessTravelOFF{3}'];
xIdx        = [ones(size(excessTravelON{1}))';   ...
               ones(size(excessTravelON{2}))'*2; ...
               ones(size(excessTravelON{3}))'*3];
ylimits     = [-10 10];
yaxName     = 'delta trials w/ excess travel (%, on-off)';
%%%%%%%%%%%%%%%%%%%%%%%%%

groupNames  = {'indirect';'direct';'no opsin'};
catMarker   = {'x','x', 'x'};
colorIdx    = {[0.8 0.8 0.8],[0.8 0.8 0.8], [0.8 0.8 0.8]};
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

xIdx = categorical(xIdx);
[p,t,stats] = anova1(data,xIdx,'off');
text(0.6,ylimits(2)*.9,['group p= ' num2str(p)])
text(0.6,ylimits(2)*.8,['df= ' num2str(stats.df)])
text(0.6,ylimits(2)*.7,['F= ' num2str((t{2,5}))])

%% clear all except source data for plots
clearvars -except yVelOFF yVelON viewAngOFF viewAngON ...
                  xPosOFF xPosON distanceOFF distanceON ...
                  excessTravelOFF excessTravelON ...
                  vaDevOFF vaDevON groupOrder ...
                  lgDMSCtrl_vc lgDMSD1_vc lgDMSD2_vc ...
                  nTrialsClean_OFF nTrialsClean_ON ...
                  nTrialsALL_OFF nTrialsALL_ON
