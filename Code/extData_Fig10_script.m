% Script to generate ExtData Fig 14 data plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% cross GLM-HMM state comparison of motor variables
% and effects of delta (laser on-off) on motor variables:
% velocity, x-position and view angle,per trial STD, distance, excess
% travel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
cd(globalParams.dataPath)
load('DMS_AoE_GLMHMM_states.mat')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
groupOrder = {'DMSD2'  ,'DMSD1'  };
counter = 1;

for groupNum = 1:numel(groupOrder)
    clear lg lgClean lg_state1 lg_state2 lg_state3 lgTemp lgTempOFF lgTempON
    lg      = eval(['lg' groupOrder{groupNum}]);
    stateIdx = {'_state1', '_state2', '_state3'};
    
    % remove field not of same length
    if isfield(lg,'keyFrameLabels')
        lg = rmfield(lg,'keyFrameLabels');
    else
    end
    
    % just log total number of off/on trials per state
    for ss = 1:numel(stateIdx)        
        nTrialsALL_OFF(counter,ss) = sum(lg.laserON == 0 & lg.stateIndex == ss); 
        nTrialsALL_ON(counter,ss)  = sum(lg.laserON == 1 & lg.stateIndex == ss); 
    end
    
    % remove excess travel >0.1 for analysis of velocity, x-pos, view angle
    [lgClean,~] = selectLogSubset(lg,[],[],[],[],0,0.0); % excessTravel>0.1
        
    % invert logs to ipsi contra space
    lg        = invertLogs(lg);
    lgClean   = invertLogs(lgClean);
    
    % make seperate logs for each state
    removeMouse  = [];
    removeMouse  = lgClean.stateIndex == 1;
    lg_state1    = structfun(@(x) x(removeMouse),lgClean,'UniformOutput',false);
    removeMouse  = [];
    removeMouse  = lgClean.stateIndex == 2;
    lg_state2    = structfun(@(x) x(removeMouse),lgClean,'UniformOutput',false);
    removeMouse  = [];
    removeMouse  = lgClean.stateIndex == 3;
    lg_state3    = structfun(@(x) x(removeMouse),lgClean,'UniformOutput',false);
    
    % see if after trial selection a mouse does not occupy all 3 states   
    s1s3diff = []; s1s3diff = setdiff(lg_state1.mouseID,  lg_state3.mouseID);
    s2s3diff = []; s2s3diff = setdiff(lg_state2.mouseID, lg_state3.mouseID);
    s1s2diff = []; s1s2diff = setdiff(lg_state1.mouseID, lg_state2.mouseID);   
    nonUniqueMice = sort([s1s3diff s2s3diff s1s2diff]);
    idx = [false diff(nonUniqueMice)<1];  
    nonUniqueMice(idx) = [];
    
    % remove mice that don't occupy all 3 states after trial selection
    if ~isempty(nonUniqueMice)
        for nMouse = 1:numel(nonUniqueMice)
            removeMouse  = [];
            removeMouse  = lg_state1.mouseID == nonUniqueMice(nMouse);
            lg_state1    = structfun(@(x) x(~removeMouse),lg_state1,'UniformOutput',false);
            removeMouse  = [];
            removeMouse  = lg_state2.mouseID == nonUniqueMice(nMouse);
            lg_state2    = structfun(@(x) x(~removeMouse),lg_state2,'UniformOutput',false);
            removeMouse  = [];
            removeMouse  = lg_state3.mouseID == nonUniqueMice(nMouse);
            lg_state3    = structfun(@(x) x(~removeMouse),lg_state3,'UniformOutput',false);
        end
    else
    end
    
       
    for ss = 1:numel(stateIdx)
        clearvars lgTemp
        lgTemp = eval(['lg' stateIdx{ss}]);  
        
        % this adds field pos2 - hack for convenience of function handling
        % pos is x-position, y-position, view angle
        % pos2 is view angle, y-position, x-position        
        pos2    = cellfun(@fliplr, lgTemp.pos,'UniformOutput',false); %
        lgTemp.pos2 = pos2;
        
        % create seperate logs for light off and on trials
        removeMouse  = [];
        removeMouse  = lgTemp.laserON==1;
        lgTempOFF    = structfun(@(x) x(~removeMouse),lgTemp,'UniformOutput',false);
        lgTempON     = structfun(@(x) x(removeMouse),lgTemp,'UniformOutput',false);
        % just log total number of off/on trials per state post selection here for convenience
        nTrialsSUBSELECT_OFF(counter,ss) = numel(lgTempOFF.choice); 
        nTrialsSUBSELECT_ON(counter,ss)  = numel(lgTempON.choice);       
        
        % calculation of Yvelocity in bins of 25cm for 300cm stem
        % for each state laser off for cross state comp panel 14b,p
        % also for laser off and on for within state comp of laser effects this is panel i,w 
        binIDs = 1:25:301; 
        for nBin = 1:numel(binIDs)-1
            yVelOFF_xState{counter,ss}(nBin,:)  = xMouseSpeedXY(lgTempOFF,[binIDs(nBin) binIDs(nBin+1)], 2);
            yVelON_xState{counter,ss}(nBin,:)   = xMouseSpeedXY(lgTempON,[binIDs(nBin) binIDs(nBin+1)], 2);
        end       
        
        
        % calculation of per-trial average x-pos 0-300cm in 25cm bins only for laser off trials
        avgXPos = [];
        avgXpos = sampleViewAngleVsY_average(lgTempOFF.pos2,[1 301],25);   
        
        % gets per mouse average x-position in each spatial bin on left or right choice trials 
        % panels c and q - only for laser off trials
        mouseNums  = [];
        mouseNums  = unique(lgTempOFF.mouseID);
        for nMouse = 1:numel(mouseNums)
            xPosL_xState{counter,ss}(:,nMouse) = ...
                nanmean(avgXpos(:,lgTempOFF.mouseID == mouseNums(nMouse) & ...
                lgTempOFF.choice == 0),2);
            xPosR_xState{counter,ss}(:,nMouse) = ...
                nanmean(avgXpos(:,lgTempOFF.mouseID == mouseNums(nMouse) & ...
                lgTempOFF.choice == 1),2);
        end
        
        % same as above but for view angle
        % panels d and r - only for laser off trials
        avgViewAng  = [];
        avgViewAng  = sampleViewAngleVsY_average(lgTempOFF.pos,[1 301],25);
        
        mouseNums       = [];
        mouseNums       = unique(lgTempOFF.mouseID);
        for nMouse = 1:numel(mouseNums)
            viewAngL_xState{counter,ss}(:,nMouse) = ...
                nanmean(avgViewAng(:,lgTempOFF.mouseID == mouseNums(nMouse) & ...
                lgTempOFF.choice == 0),2);
            viewAngR_xState{counter,ss}(:,nMouse) = ...
                nanmean(avgViewAng(:,lgTempOFF.mouseID == mouseNums(nMouse) & ...
                lgTempOFF.choice == 1),2);
        end
        
        % calculation of delta (on-off) x-position (cm) in maze stem 0-200cm
        % panel j and x
        avgXpos = [];
        mazeSection = []; mazeSection = [1 201];
       
        avgXpos = sampleViewAngleVsY_average(lgTemp.pos2,mazeSection,mazeSection(2)-mazeSection(1));
        mouseNums  = [];
        mouseNums  = unique(lgTemp.mouseID);
        for nMouse = 1:numel(mouseNums)
            xPos_OFF_xState{counter,ss}(:,nMouse) = ...
                nanmean(avgXpos(:,lgTemp.mouseID == mouseNums(nMouse) & ...
                lgTemp.laserON == 0),2);
            xPos_ON_xState{counter,ss}(:,nMouse) = ...
                nanmean(avgXpos(:,lgTemp.mouseID == mouseNums(nMouse) & ...
                lgTemp.laserON == 1),2);
        end
        
        % calculation of delta (on-off) view angle (cm) in maze stem 0-200cm
        % panel k and y
        avgViewAng  = [];
        avgViewAng  = sampleViewAngleVsY_average(lgTemp.pos,mazeSection,mazeSection(2)-mazeSection(1));
        mouseNums       = [];
        mouseNums       = unique(lgTemp.mouseID);
        for nMouse = 1:numel(mouseNums)
            viewAng_OFF_xState{counter,ss}(:,nMouse) = ...
                nanmean(avgViewAng(:,lgTemp.mouseID == mouseNums(nMouse) & ...
                lgTemp.laserON == 0),2);
            viewAng_ON_xState{counter,ss}(:,nMouse) = ...
                nanmean(avgViewAng(:,lgTemp.mouseID == mouseNums(nMouse) & ...
                lgTemp.laserON == 1),2);
        end
        
        
    end
        
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % load log file again and don't remove excess travel for analysis
    % of distance, per-trial deviation in view angle and excess travel
    clear lg lgClean lg_state1 lg_state2 lg_state3 lgTemp lgTempOFF lgTempON
    lg      = eval(['lg' groupOrder{groupNum}]);
    stateIdx = {'_state1', '_state2', '_state3'};
    
    % remove field not of same length
    if isfield(lg,'keyFrameLabels')
        lg = rmfield(lg,'keyFrameLabels');
    else
    end
    
    % invert logs to ipsi contra space
    lg        = invertLogs(lg);
    
    % make seperate logs for each state
    removeMouse  = [];
    removeMouse  = lg.stateIndex == 1;
    lg_state1    = structfun(@(x) x(removeMouse),lg,'UniformOutput',false);
    removeMouse  = [];
    removeMouse  = lg.stateIndex == 2;
    lg_state2    = structfun(@(x) x(removeMouse),lg,'UniformOutput',false);
    removeMouse  = [];
    removeMouse  = lg.stateIndex == 3;
    lg_state3    = structfun(@(x) x(removeMouse),lg,'UniformOutput',false);
    
    % see if a mouse does not occupy all 3 states
    s1s3diff = []; s1s3diff = setdiff(lg_state1.mouseID,  lg_state3.mouseID);
    s2s3diff = []; s2s3diff = setdiff(lg_state2.mouseID, lg_state3.mouseID);
    s1s2diff = []; s1s2diff = setdiff(lg_state1.mouseID, lg_state2.mouseID);
    nonUniqueMice = sort([s1s3diff s2s3diff s1s2diff]);
    idx = [false diff(nonUniqueMice)<1];
    nonUniqueMice(idx) = [];
    
    % remove mice that don't occupy all 3 states
    if ~isempty(nonUniqueMice)
        for nMouse = 1:numel(nonUniqueMice)
            removeMouse  = [];
            removeMouse  = lg_state1.mouseID == nonUniqueMice(nMouse);
            lg_state1    = structfun(@(x) x(~removeMouse),lg_state1,'UniformOutput',false);
            removeMouse  = [];
            removeMouse  = lg_state2.mouseID == nonUniqueMice(nMouse);
            lg_state2    = structfun(@(x) x(~removeMouse),lg_state2,'UniformOutput',false);
            removeMouse  = [];
            removeMouse  = lg_state3.mouseID == nonUniqueMice(nMouse);
            lg_state3    = structfun(@(x) x(~removeMouse),lg_state3,'UniformOutput',false);
        end
    else
    end
    
    for ss = 1:numel(stateIdx)
        clearvars lgTemp
        lgTemp = eval(['lg' stateIdx{ss}]);
        
        % create seperate logs for light off and on trials
        removeMouse  = [];
        removeMouse  = lgTemp.laserON==1;
        lgTempOFF    = structfun(@(x) x(~removeMouse),lgTemp,'UniformOutput',false);
        lgTempON     = structfun(@(x) x(removeMouse),lgTemp,'UniformOutput',false);
        
        % this adds field pos2 - hack for convenience of function handling
        % pos is x-position, y-position, view angle
        % pos2 is view angle, y-position, x-position
        pos2    = cellfun(@fliplr, lgTemp.pos,'UniformOutput',false); %
        lgTemp.pos2 = pos2;
        
        % panel e and l and s and z- calculation of distance
        % for laser off and laser on trials seperately
        distance = [];
        distance = 330+(lgTemp.excessTravel*330);
        mouseNums = [];
        mouseNums = unique(lgTemp.mouseID);
        for nMouse = 1:numel(mouseNums)
            distanceOFF_xState{counter,ss}(:,nMouse) = nanmean(distance(lgTemp.mouseID == mouseNums(nMouse)...
                & lgTemp.laserON==0));
            distanceON_xState{counter,ss}(:,nMouse)  = nanmean(distance(lgTemp.mouseID == mouseNums(nMouse)...
                & lgTemp.laserON==1));           
        end
        
        % panel f and m and t and aa - calculation of % trials with excessTravel
        % for laser off and laser on trials seperately
        mouseNums = []; etOFF = []; etON = [];
        mouseNums = unique(lgTemp.mouseID);
        for nMouse = 1:numel(mouseNums)
            etOFF(nMouse) = (sum(lgTemp.excessTravel(lgTemp.mouseID == mouseNums(nMouse) & lgTemp.laserON == 0) > .1 )/...
                numel(lgTemp.excessTravel(lgTemp.mouseID == mouseNums(nMouse) & lgTemp.laserON == 0)))*100;
            etON(nMouse)  = (sum(lgTemp.excessTravel(lgTemp.mouseID == mouseNums(nMouse) & lgTemp.laserON == 1) > .1 )/...
                numel(lgTemp.excessTravel(lgTemp.mouseID == mouseNums(nMouse) & lgTemp.laserON == 1)))*100;
        end
        excTravelOFF_xState{counter,ss} = etOFF;
        excTravelON_xState{counter,ss}  = etON;
        
        % panel g n uu bb - calculation of view angle deviation
        posBins = []; posBins = 1:5:301;
        vadOFF                      = [];
        vadOFF                      = xMouseViewAngSD(lgTempOFF, posBins);
        vaDevOFF_xState{counter,ss}(:) = vadOFF.viewAngSD;
        vadON                       = [];
        vadON                       = xMouseViewAngSD(lgTempON, posBins);
        vaDevON_xState{counter,ss}(:)  = vadON.viewAngSD;
              
    end
       
    % for the next loop of analysis given lgs in groupOrder
    counter = counter +1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% clear all vars except source data for plotting
clearvars -except groupOrder stateIdx lgDMSD1 lgDMSD2 ...
                  nTrialsSUBSELECT_OFF nTrialsSUBSELECT_ON ...
                  nTrialsALL_OFF nTrialsALL_ON ...
                  nTrialsMATCHMICE_OFF nTrialsMATCHMICE_ON ...
                  viewAng_OFF_xState viewAng_ON_xState ...
                  viewAngL_xState viewAngR_xState ...
                  xPos_OFF_xState xPos_ON_xState ...
                  xPosL_xState xPosR_xState ... 
                  yVelOFF_xState yVelON_xState...                 
                  vaDevOFF_xState vaDevON_xState ...
                  excTravelOFF_xState excTravelON_xState ...
                  distanceOFF_xState distanceON_xState
              
%% PLOTTING BASED ON ExtData Fig 14 source data above
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% update this section for plotting of measured vars above
% plot panel ext data b and p : y-velocity laser off across states           
dataArray   = yVelOFF_xState;
ylimits     = [0 300];
xlimits     = [20 100];
xAxTicks    = [25 50 75 100];
xAxName     = 'y-velocity (cm/s, on-off)';
yAxName     = 'y-position (cm)';
%%%%%%%%%%%%%%%%%%%%%%%%%

binIDs     = 25:25:300;
plotOffset = [-17.5 -12.5 -7.5];
for nGroup = 1:size(dataArray,1)
    figure
    for nState = 1:size(dataArray,2)  
        mean2plot = []; err2plot = [];
        mean2plot = nanmean(dataArray{nGroup,nState},2);
        err2plot  = (nanstd(dataArray{nGroup,nState},0,2)) /...
                     sqrt(size(dataArray{nGroup,nState},2)-1);                 
        errorbar(mean2plot, binIDs' + plotOffset(nState) ,err2plot,'horizontal','Color',globalParams.stateColors{nState}); hold on
    end
        
    xlim(xlimits)
    set(gca,'XTick',xAxTicks)
    box off
    set(gca,'TickDir','out')
        
    if contains(groupOrder{nGroup},'D2')
        title('indirect')    
    elseif contains(groupOrder{nGroup},'D1')
        title('direct')
    else
    end
    
    ylabel(yAxName)
    xlabel(xAxName)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% clear all vars except source data for plotting
clearvars -except groupOrder lgDMSD1 lgDMSD2 ...
                  nTrialsSUBSELECT_OFF nTrialsSUBSELECT_ON ...
                  nTrialsALL_OFF nTrialsALL_ON ...
                  nTrialsMATCHMICE_OFF nTrialsMATCHMICE_ON ...
                  viewAng_OFF_xState viewAng_ON_xState ...
                  viewAngL_xState viewAngR_xState ...
                  xPos_OFF_xState xPos_ON_xState ...
                  xPosL_xState xPosR_xState ... 
                  yVelOFF_xState yVelON_xState...                 
                  vaDevOFF_xState vaDevON_xState ...
                  excTravelOFF_xState excTravelON_xState ...
                  distanceOFF_xState distanceON_xState

%% PLOTTING BASED ON ExtData Fig 14 source data above
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% update this section for plotting of measured vars above
% plot panel ext data c and q : laser off binned x-position across states           
dataArray1   = xPosL_xState;
dataArray2   = xPosR_xState;
ylimits     = [0 300];
xlimits1     = [-1.5 0.25];
xAxTicks1    = [-1 0];
xlimits2     = [-0.25 1.5];
xAxTicks2    = [0 1];
xAxName     = 'x-position (cm)';
yAxName     = 'y-position (cm)';
%%%%%%%%%%%%%%%%%%%%%%%%%

binIDs     = 25:25:300;
plotOffset = [-17.5 -12.5 -7.5];
for nGroup = 1:size(dataArray1,1)
    figure
    for nState = 1:size(dataArray1,2)  
        subplot(1,2,1)
        mean2plot = []; err2plot = [];
        mean2plot = nanmean(dataArray1{nGroup,nState},2);
        err2plot  = (nanstd(dataArray1{nGroup,nState},0,2)) /...
                     sqrt(size(dataArray1{nGroup,nState},2)-1);                 
        errorbar(mean2plot, binIDs' + plotOffset(nState) ,err2plot,'horizontal','Color',globalParams.stateColors{nState}); hold on
        plot(zeros(300,1),[1:300],'--','Color',[0.8 0.8 0.8])
        xlim(xlimits1)
        set(gca,'XTick',xAxTicks1)
        box off
        set(gca,'TickDir','out')
        ylabel(yAxName)
        xlabel(xAxName)
        
        subplot(1,2,2)
        mean2plot = []; err2plot = [];
        mean2plot = nanmean(dataArray2{nGroup,nState},2);
        err2plot  = (nanstd(dataArray2{nGroup,nState},0,2)) /...
                     sqrt(size(dataArray2{nGroup,nState},2)-1);                 
        errorbar(mean2plot, binIDs' + plotOffset(nState) ,err2plot,'horizontal','Color',globalParams.stateColors{nState}); hold on
        xlim(xlimits2)
        set(gca,'XTick',xAxTicks2)
        box off
        set(gca,'TickDir','out')
        xlabel(xAxName)
        plot(zeros(300,1),[1:300],'--','Color',[0.8 0.8 0.8])
        set(gca,'YAxisLocation','right');
    end
        
     subplot(1,2,1)
    if contains(groupOrder{nGroup},'D2')
        title('indirect')    
    elseif contains(groupOrder{nGroup},'D1')
        title('direct')
    else
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% clear all vars except source data for plotting
clearvars -except groupOrder lgDMSD1 lgDMSD2 ...
                  nTrialsSUBSELECT_OFF nTrialsSUBSELECT_ON ...
                  nTrialsALL_OFF nTrialsALL_ON ...
                  nTrialsMATCHMICE_OFF nTrialsMATCHMICE_ON ...
                  viewAng_OFF_xState viewAng_ON_xState ...
                  viewAngL_xState viewAngR_xState ...
                  xPos_OFF_xState xPos_ON_xState ...
                  xPosL_xState xPosR_xState ... 
                  yVelOFF_xState yVelON_xState...                 
                  vaDevOFF_xState vaDevON_xState ...
                  excTravelOFF_xState excTravelON_xState ...
                  distanceOFF_xState distanceON_xState              

              
%% PLOTTING BASED ON ExtData Fig 14 source data above
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% update this section for plotting of measured vars above
% plot panel ext data d and r : laser off binned view angle across states           
dataArray1   = viewAngL_xState;
dataArray2   = viewAngR_xState;
ylimits     = [0 300];
xlimits1     = [-45 10];
xAxTicks1    = [-30 -20 -10 0 10];
xlimits2     = [-10 45];
xAxTicks2    = [-10 0 10 20 30];
xAxName     = 'view angle (deg)';
yAxName     = 'y-position (cm)';
%%%%%%%%%%%%%%%%%%%%%%%%%

binIDs     = 25:25:300;
plotOffset = [-17.5 -12.5 -7.5];
for nGroup = 1:size(dataArray1,1)
    figure
    for nState = 1:size(dataArray1,2)  
        subplot(1,2,1)
        mean2plot = []; err2plot = [];
        mean2plot = nanmean(dataArray1{nGroup,nState},2);
        err2plot  = (nanstd(dataArray1{nGroup,nState},0,2)) /...
                     sqrt(size(dataArray1{nGroup,nState},2)-1);                 
        errorbar(mean2plot, binIDs' + plotOffset(nState) ,err2plot,'horizontal','Color',globalParams.stateColors{nState}); hold on
        plot(zeros(300,1),[1:300],'--','Color',[0.8 0.8 0.8])
        xlim(xlimits1)
        set(gca,'XTick',xAxTicks1)
        box off
        set(gca,'TickDir','out')
        ylabel(yAxName)
        xlabel(xAxName)
        
        subplot(1,2,2)
        mean2plot = []; err2plot = [];
        mean2plot = nanmean(dataArray2{nGroup,nState},2);
        err2plot  = (nanstd(dataArray2{nGroup,nState},0,2)) /...
                     sqrt(size(dataArray2{nGroup,nState},2)-1);                 
        errorbar(mean2plot, binIDs' + plotOffset(nState) ,err2plot,'horizontal','Color',globalParams.stateColors{nState}); hold on
        xlim(xlimits2)
        set(gca,'XTick',xAxTicks2)
        box off
        set(gca,'TickDir','out')
        xlabel(xAxName)
        plot(zeros(300,1),[1:300],'--','Color',[0.8 0.8 0.8])
        set(gca,'YAxisLocation','right');
    end
        
     subplot(1,2,1)
    if contains(groupOrder{nGroup},'D2')
        title('indirect')    
    elseif contains(groupOrder{nGroup},'D1')
        title('direct')
    else
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% clear all vars except source data for plotting
clearvars -except groupOrder lgDMSD1 lgDMSD2 ...
                  nTrialsSUBSELECT_OFF nTrialsSUBSELECT_ON ...
                  nTrialsALL_OFF nTrialsALL_ON ...
                  nTrialsMATCHMICE_OFF nTrialsMATCHMICE_ON ...
                  viewAng_OFF_xState viewAng_ON_xState ...
                  viewAngL_xState viewAngR_xState ...
                  xPos_OFF_xState xPos_ON_xState ...
                  xPosL_xState xPosR_xState ... 
                  yVelOFF_xState yVelON_xState...                 
                  vaDevOFF_xState vaDevON_xState ...
                  excTravelOFF_xState excTravelON_xState ...
                  distanceOFF_xState distanceON_xState                  

              
%% PLOTTING BASED ON ExtData Fig 14 source data above
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% update this section for plotting of measured vars above
% plot panel ext data d and r : laser off binned view angle across states           
dataArray1   = vaDevOFF_xState;
ylimits1     = [0 800];
yAxName1     = 'per-trial view angle STD (deg)';

dataArray2   = distanceOFF_xState;
ylimits2     = [300 950];
yAxName2     = 'distance(cm)';

dataArray3   = excTravelOFF_xState;
ylimits3     = [0 40];
yAxName3     = '% trials with excess travel';

xlimits     = [0.5 3.5];
xTickLabels = {'s1','s2','s3'};

%%%%%%%%%%%%%%%%%%%%%%%%%

for nGroup = 1:size(dataArray1,1)
    figure
    for nState = 1:size(dataArray1,2)  
        subplot(1,3,1)
        mean2plot = []; err2plot = [];
        mean2plot = nanmean(dataArray1{nGroup,nState},2);
        err2plot  = (nanstd(dataArray1{nGroup,nState},0,2)) /...
                     sqrt(size(dataArray1{nGroup,nState},2)-1);                                     
        errorbar(nState, mean2plot,err2plot,'Color',globalParams.stateColors{nState}); hold on
        ylim(ylimits1)
        box off
        set(gca,'TickDir','out')
        ylabel(yAxName1)      
        xlim(xlimits)
        
        subplot(1,3,2)
        mean2plot = []; err2plot = [];
        mean2plot = nanmean(dataArray2{nGroup,nState},2);
        err2plot  = (nanstd(dataArray2{nGroup,nState},0,2)) /...
                     sqrt(size(dataArray2{nGroup,nState},2)-1);                                     
        errorbar(nState, mean2plot,err2plot,'Color',globalParams.stateColors{nState}); hold on
        ylim(ylimits2)
        box off
        set(gca,'TickDir','out')
        ylabel(yAxName2)      
        xlim(xlimits)
        
        subplot(1,3,3)
        mean2plot = []; err2plot = [];
        mean2plot = nanmean(dataArray3{nGroup,nState},2);
        err2plot  = (nanstd(dataArray3{nGroup,nState},0,2)) /...
                     sqrt(size(dataArray3{nGroup,nState},2)-1);                                     
        errorbar(nState, mean2plot,err2plot,'Color',globalParams.stateColors{nState}); hold on
        ylim(ylimits3)
        box off
        set(gca,'TickDir','out')
        ylabel(yAxName3)      
        xlim(xlimits)    
    end
        
    subplot(1,3,1)
    plot(1:3,[dataArray1{nGroup,1}' ...
              dataArray1{nGroup,2}' ...
              dataArray1{nGroup,3}'], ...
              '-', 'Color', [0.8 0.8 0.8]);
    subplot(1,3,2)
    plot(1:3,[dataArray2{nGroup,1}' ...
              dataArray2{nGroup,2}' ...
              dataArray2{nGroup,3}'], ...
              '-', 'Color', [0.8 0.8 0.8]);
    subplot(1,3,3)
    plot(1:3,[dataArray3{nGroup,1}' ...
              dataArray3{nGroup,2}' ...
              dataArray3{nGroup,3}'], ...
              '-', 'Color', [0.8 0.8 0.8]);
    
    data1 = []; data2 = []; data3 = [];    
    data1   = [dataArray1{nGroup,1}'; ...
               dataArray1{nGroup,2}'; ...
               dataArray1{nGroup,3}'];              
    xIdx1   = [ones(size(dataArray1{nGroup,1},2),1);   ...
               ones(size(dataArray1{nGroup,2},2),1)*2; ...
               ones(size(dataArray1{nGroup,3},2),1)*3];
    subj1   = [[1:numel(dataArray1{nGroup,1})]';   ...
               [1:numel(dataArray1{nGroup,2})]'; ...
               [1:numel(dataArray1{nGroup,3})]'];
     
    data2   = [dataArray2{nGroup,1}'; ...
               dataArray2{nGroup,2}'; ...
               dataArray2{nGroup,3}'];              
    xIdx2   = [ones(size(dataArray2{nGroup,1},2),1);   ...
               ones(size(dataArray2{nGroup,2},2),1)*2; ...
               ones(size(dataArray2{nGroup,3},2),1)*3];
    subj2   = [[1:numel(dataArray2{nGroup,1})]';   ...
               [1:numel(dataArray2{nGroup,2})]'; ...
               [1:numel(dataArray2{nGroup,3})]'];
          
    data3   = [dataArray3{nGroup,1}'; ...
               dataArray3{nGroup,2}'; ...
               dataArray3{nGroup,3}'];              
    xIdx3   = [ones(size(dataArray3{nGroup,1},2),1);   ...
               ones(size(dataArray3{nGroup,2},2),1)*2; ...
               ones(size(dataArray3{nGroup,3},2),1)*3];
    subj3   = [[1:numel(dataArray3{nGroup,1})]';   ...
               [1:numel(dataArray3{nGroup,2})]'; ...
               [1:numel(dataArray3{nGroup,3})]'];       
              
    X = []; pState = []; fState = [];
    X(:,1) = data1; 
    X(:,2) = xIdx1; 
    X(:,3) = subj1;
    [pState,~,fState] = RMAOV1(X,.05);   
    subplot(1,3,1)
    text(0.6,ylimits1(2)*.95,['state p= ' num2str(pState)])
    text(0.6,ylimits1(2)*.90,['F= ' num2str(fState)])
    set(gca,'XTickLabel',{'s1','s2','s3'});     
    if contains(groupOrder{nGroup},'D2')
        title('indirect')    
    elseif contains(groupOrder{nGroup},'D1')
        title('direct')
    else
    end
    
    X = []; pState = []; fState = [];
    X(:,1) = data2; 
    X(:,2) = xIdx2; 
    X(:,3) = subj2;
    [pState,~,fState] = RMAOV1(X,.05);   
    subplot(1,3,2)
    text(0.6,ylimits2(2)*.95,['state p= ' num2str(pState)])
    text(0.6,ylimits2(2)*.90,['F= ' num2str(fState)])
    set(gca,'XTickLabel',{'s1','s2','s3'});     
    if contains(groupOrder{nGroup},'D2')
        title('indirect')    
    elseif contains(groupOrder{nGroup},'D1')
        title('direct')
    else
    end
    
    X = []; pState = []; fState = [];
    X(:,1) = data3; 
    X(:,2) = xIdx3; 
    X(:,3) = subj3;
    [pState,~,fState] = RMAOV1(X,.05);   
    subplot(1,3,3)
    text(0.6,ylimits3(2)*.95,['state p= ' num2str(pState)])
    text(0.6,ylimits3(2)*.90,['F= ' num2str(fState)])
    set(gca,'XTickLabel',{'s1','s2','s3'});     
    if contains(groupOrder{nGroup},'D2')
        title('indirect')    
    elseif contains(groupOrder{nGroup},'D1')
        title('direct')
    else
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% clear all vars except source data for plotting
clearvars -except groupOrder lgDMSD1 lgDMSD2 ...
                  nTrialsSUBSELECT_OFF nTrialsSUBSELECT_ON ...
                  nTrialsALL_OFF nTrialsALL_ON ...
                  nTrialsMATCHMICE_OFF nTrialsMATCHMICE_ON ...
                  viewAng_OFF_xState viewAng_ON_xState ...
                  viewAngL_xState viewAngR_xState ...
                  xPos_OFF_xState xPos_ON_xState ...
                  xPosL_xState xPosR_xState ... 
                  yVelOFF_xState yVelON_xState...                 
                  vaDevOFF_xState vaDevON_xState ...
                  excTravelOFF_xState excTravelON_xState ...
                  distanceOFF_xState distanceON_xState

%% PLOTTING BASED ON ExtData Fig 14 source data above
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% update this section for plotting of measured vars above
% plot panel ext data i and w : laser off/on binned y-vel across states           
dataArrayOFF   = yVelOFF_xState;
dataArrayON    = yVelON_xState;
ylimits     = [0 300];
xlimits     = [20 100];
xAxTicks    = [25 50 75 100];
xAxName     = 'y-velocity (cm/s, on-off)';
yAxName     = 'y-position (cm)';
%%%%%%%%%%%%%%%%%%%%%%%%%

binIDs     = 25:25:300;
plotOffset = [-17.5 -12.5 -7.5];
for nGroup = 1:size(dataArrayOFF,1)
    figure
    for nState = 1:size(dataArrayOFF,2)  
        subplot(1,3,nState)
        mean2plot = []; err2plot = [];
        mean2plot = nanmean(dataArrayOFF{nGroup,nState},2);
        err2plot  = (nanstd(dataArrayOFF{nGroup,nState},0,2)) /...
                     sqrt(size(dataArrayOFF{nGroup,nState},2)-1);                 
        errorbar(mean2plot, binIDs' + plotOffset(1) ,err2plot,'horizontal','Color','k'); hold on
        xlim(xlimits)
        set(gca,'XTick',xAxTicks)
        box off
        set(gca,'TickDir','out')
        title(['state_' num2str(nState)])
        if nState == 1
            ylabel(yAxName)            
        elseif nState == 2
            xlabel(xAxName)
        else
        end
        
        subplot(1,3,nState)
        mean2plot = []; err2plot = [];
        mean2plot = nanmean(dataArrayON{nGroup,nState},2);
        err2plot  = (nanstd(dataArrayON{nGroup,nState},0,2)) /...
                     sqrt(size(dataArrayON{nGroup,nState},2)-1);                 
        errorbar(mean2plot, binIDs' + plotOffset(2) ,err2plot,'horizontal','Color',globalParams.laserColor); hold on
        xlim(xlimits)
        set(gca,'XTick',xAxTicks)
        box off
        set(gca,'TickDir','out')
        title(['state_' num2str(nState)])
        if nState == 1
            ylabel(yAxName)            
        elseif nState == 2
            xlabel(xAxName)
        else
        end
        
    end
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% clear all vars except source data for plotting
clearvars -except groupOrder lgDMSD1 lgDMSD2 ...
                  nTrialsSUBSELECT_OFF nTrialsSUBSELECT_ON ...
                  nTrialsALL_OFF nTrialsALL_ON ...
                  nTrialsMATCHMICE_OFF nTrialsMATCHMICE_ON ...
                  viewAng_OFF_xState viewAng_ON_xState ...
                  viewAngL_xState viewAngR_xState ...
                  xPos_OFF_xState xPos_ON_xState ...
                  xPosL_xState xPosR_xState ... 
                  yVelOFF_xState yVelON_xState...                 
                  vaDevOFF_xState vaDevON_xState ...
                  excTravelOFF_xState excTravelON_xState ...
                  distanceOFF_xState distanceON_xState     
              
%% PLOTTING BASED ON ExtData Fig 14 source data above
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% update this section for plotting of measured vars above
% plot panel ext data j and k and x and y
% laser off binned view angle across states           
dataArray1   = cellfun(@minus,xPos_ON_xState,xPos_OFF_xState,'UniformOutput',false);
ylimits1     = [-0.5 0.5];
yAxName1     = 'delta x-position (cm, on-off)';

dataArray2   = cellfun(@minus,viewAng_ON_xState,viewAng_OFF_xState,'UniformOutput',false);;
ylimits2     = [-20 20];
yAxName2     = 'delta view angle (deg, on-off)';

xlimits      = [0.5 3.5];
xTickLabels  = {'s1','s2','s3'}; 

%%%%%%%%%%%%%%%%%%%%%%%%%

for nGroup = 1:size(dataArray1,1)
    figure
    for nState = 1:size(dataArray1,2)  
        subplot(1,2,1)
        mean2plot = []; err2plot = [];
        mean2plot = nanmean(dataArray1{nGroup,nState},2);
        err2plot  = (nanstd(dataArray1{nGroup,nState},0,2)) /...
                     sqrt(size(dataArray1{nGroup,nState},2)-1);                                     
        errorbar(nState, mean2plot,err2plot,'Color',globalParams.stateColors{nState}); hold on
        ylim(ylimits1)
        box off
        set(gca,'TickDir','out')
        ylabel(yAxName1)      
        xlim(xlimits)
        
        subplot(1,2,2)
        mean2plot = []; err2plot = [];
        mean2plot = nanmean(dataArray2{nGroup,nState},2);
        err2plot  = (nanstd(dataArray2{nGroup,nState},0,2)) /...
                     sqrt(size(dataArray2{nGroup,nState},2)-1);                                     
        errorbar(nState, mean2plot,err2plot,'Color',globalParams.stateColors{nState}); hold on
        ylim(ylimits2)
        box off
        set(gca,'TickDir','out')
        ylabel(yAxName2)      
        xlim(xlimits)
        
    end
        
    subplot(1,2,1)
    plot(1:3,[dataArray1{nGroup,1}' ...
              dataArray1{nGroup,2}' ...
              dataArray1{nGroup,3}'], ...
              '-', 'Color', [0.8 0.8 0.8]);
    subplot(1,2,2)
    plot(1:3,[dataArray2{nGroup,1}' ...
              dataArray2{nGroup,2}' ...
              dataArray2{nGroup,3}'], ...
              '-', 'Color', [0.8 0.8 0.8]);
    
    data1 = []; data2 = []; data3 = [];    
    data1   = [dataArray1{nGroup,1}'; ...
               dataArray1{nGroup,2}'; ...
               dataArray1{nGroup,3}'];              
    xIdx1   = [ones(size(dataArray1{nGroup,1},2),1);   ...
               ones(size(dataArray1{nGroup,2},2),1)*2; ...
               ones(size(dataArray1{nGroup,3},2),1)*3];
    subj1   = [[1:numel(dataArray1{nGroup,1})]';   ...
               [1:numel(dataArray1{nGroup,2})]'; ...
               [1:numel(dataArray1{nGroup,3})]'];
     
    data2   = [dataArray2{nGroup,1}'; ...
               dataArray2{nGroup,2}'; ...
               dataArray2{nGroup,3}'];              
    xIdx2   = [ones(size(dataArray2{nGroup,1},2),1);   ...
               ones(size(dataArray2{nGroup,2},2),1)*2; ...
               ones(size(dataArray2{nGroup,3},2),1)*3];
    subj2   = [[1:numel(dataArray2{nGroup,1})]';   ...
               [1:numel(dataArray2{nGroup,2})]'; ...
               [1:numel(dataArray2{nGroup,3})]'];
                           
    X = []; pState = []; fState = [];
    X(:,1) = data1; 
    X(:,2) = xIdx1; 
    X(:,3) = subj1;
    [pState,~,fState] = RMAOV1(X,.05);   
    subplot(1,2,1)
    text(0.6,ylimits1(1)*.90,['state p= ' num2str(pState)])
    text(0.6,ylimits1(1)*.95,['F= ' num2str(fState)])
    plot(0.5:3.5,zeros(1,4),'-','Color', [0.8 0.8 0.8])
    set(gca,'XTickLabel',{'s1','s2','s3'});     
    if contains(groupOrder{nGroup},'D2')
        title('indirect')    
    elseif contains(groupOrder{nGroup},'D1')
        title('direct')
    else
    end
    
    X = []; pState = []; fState = [];
    X(:,1) = data2; 
    X(:,2) = xIdx2; 
    X(:,3) = subj2;
    [pState,~,fState] = RMAOV1(X,.05);   
    subplot(1,2,2)
    text(0.6,ylimits2(1)*.90,['state p= ' num2str(pState)])
    text(0.6,ylimits2(1)*.95,['F= ' num2str(fState)])
    plot(0.5:3.5,zeros(1,4),'-','Color', [0.8 0.8 0.8])    
    set(gca,'XTickLabel',{'s1','s2','s3'});     
    if contains(groupOrder{nGroup},'D2')
        title('indirect')    
    elseif contains(groupOrder{nGroup},'D1')
        title('direct')
    else
    end
    
    subplot(1,2,1)
   [s1Vs2_p, ~ , s1Vs2_stats] = signrank(dataArray1{nGroup,1}, ...
        dataArray1{nGroup,2},'method','approximate');    
    plot(1:2,ones(1,2)*ylimits1(2)*.8,'k-')
    text(0.5,ylimits1(2)*.85,['p= ' num2str(s1Vs2_p)],'FontSize',6)
    text(0.5,ylimits1(2)*.8,['z= ' num2str(s1Vs2_stats.zval)],'FontSize',6)
    text(0.5,ylimits1(2)*.75,['df= ' num2str(numel(dataArray1{nGroup,1})+ ...
        numel(dataArray1{nGroup,2})-1)],'FontSize',6)

    [s1Vs3_p, ~ , s1Vs3_stats] = signrank(dataArray1{nGroup,1}, ...
        dataArray1{nGroup,3},'method','approximate');
    plot(1:3,ones(1,3)*ylimits1(2)*.9,'k-')
    text(1.5,ylimits1(2)*.95,['p= ' num2str(s1Vs3_p)],'FontSize',6)
    text(1.5,ylimits1(2)*.9,['z= ' num2str(s1Vs3_stats.zval)],'FontSize',6)
    text(1.5,ylimits1(2)*.85,['df= ' num2str(numel(dataArray1{nGroup,1})+ ...
        numel(dataArray1{nGroup,3})-1)],'FontSize',6)
       
    [s2Vs3_p, ~ , s2Vs3_stats] = signrank(dataArray1{nGroup,2}, ...
        dataArray1{nGroup,3},'method','approximate');
    plot(2:3,ones(1,2)*ylimits1(2)*.7,'k-')
    text(2.5,ylimits1(2)*.75,['p= ' num2str(s2Vs3_p)],'FontSize',6)
    text(2.5,ylimits1(2)*.7,['z= ' num2str(s2Vs3_stats.zval)],'FontSize',6)
    text(2.5,ylimits1(2)*.65,['df= ' num2str(numel(dataArray1{nGroup,2})+ ...
        numel(dataArray1{nGroup,3})-2)],'FontSize',6)
    
    subplot(1,2,2)
   [s1Vs2_p, ~ , s1Vs2_stats] = signrank(dataArray2{nGroup,1}, ...
        dataArray2{nGroup,2},'method','approximate');    
    plot(1:2,ones(1,2)*ylimits2(2)*.8,'k-')
    text(0.5,ylimits2(2)*.85,['p= ' num2str(s1Vs2_p)],'FontSize',6)
    text(0.5,ylimits2(2)*.8,['z= ' num2str(s1Vs2_stats.zval)],'FontSize',6)
    text(0.5,ylimits2(2)*.75,['df= ' num2str(numel(dataArray2{nGroup,1})+ ...
        numel(dataArray2{nGroup,2})-1)],'FontSize',6)

    [s1Vs3_p, ~ , s1Vs3_stats] = signrank(dataArray2{nGroup,1}, ...
        dataArray2{nGroup,3},'method','approximate');
    plot(1:3,ones(1,3)*ylimits2(2)*.9,'k-')
    text(1.5,ylimits2(2)*.95,['p= ' num2str(s1Vs3_p)],'FontSize',6)
    text(1.5,ylimits2(2)*.9,['z= ' num2str(s1Vs3_stats.zval)],'FontSize',6)
    text(1.5,ylimits2(2)*.85,['df= ' num2str(numel(dataArray2{nGroup,1})+ ...
        numel(dataArray2{nGroup,3})-1)],'FontSize',6)
       
    [s2Vs3_p, ~ , s2Vs3_stats] = signrank(dataArray2{nGroup,2}, ...
        dataArray2{nGroup,3},'method','approximate');
    plot(2:3,ones(1,2)*ylimits2(2)*.7,'k-')
    text(2.5,ylimits2(2)*.75,['p= ' num2str(s2Vs3_p)],'FontSize',6)
    text(2.5,ylimits2(2)*.7,['z= ' num2str(s2Vs3_stats.zval)],'FontSize',6)
    text(2.5,ylimits2(2)*.65,['df= ' num2str(numel(dataArray2{nGroup,2})+ ...
        numel(dataArray2{nGroup,3})-2)],'FontSize',6)
   
end              
  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% clear all vars except source data for plotting
clearvars -except groupOrder lgDMSD1 lgDMSD2 ...
                  nTrialsSUBSELECT_OFF nTrialsSUBSELECT_ON ...
                  nTrialsALL_OFF nTrialsALL_ON ...
                  nTrialsMATCHMICE_OFF nTrialsMATCHMICE_ON ...
                  viewAng_OFF_xState viewAng_ON_xState ...
                  viewAngL_xState viewAngR_xState ...
                  xPos_OFF_xState xPos_ON_xState ...
                  xPosL_xState xPosR_xState ... 
                  yVelOFF_xState yVelON_xState...                 
                  vaDevOFF_xState vaDevON_xState ...
                  excTravelOFF_xState excTravelON_xState ...
                  distanceOFF_xState distanceON_xState     
              
%% PLOTTING BASED ON ExtData Fig 14 source data above
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% update this section for plotting of measured vars above
% plot panel ext data l,m,n and z,aa,bb 
% delta (on-off) view angle deviation, distance, % trials with excess travel           
dataArray1   = cellfun(@minus,vaDevON_xState,vaDevOFF_xState,'UniformOutput',false);
ylimits1     = [-400 800];
yAxName1     = 'delta per-trial view angle STD (deg, on-off)';

dataArray2   = cellfun(@minus,distanceON_xState,distanceOFF_xState,'UniformOutput',false);
ylimits2     = [-400 1200];
yAxName2     = 'delta distance (cm, on-off)';

dataArray3   = cellfun(@minus,excTravelON_xState,excTravelOFF_xState,'UniformOutput',false);
ylimits3     = [-25 25];
yAxName3     = 'delta trials with excess travel (%, on-off)';

xlimits     = [0.5 3.5];
xTickLabels = {'s1','s2','s3'};

%%%%%%%%%%%%%%%%%%%%%%%%%

for nGroup = 1:size(dataArray1,1)
    figure
    for nState = 1:size(dataArray1,2)  
        subplot(1,3,1)
        mean2plot = []; err2plot = [];
        mean2plot = nanmean(dataArray1{nGroup,nState},2);
        err2plot  = (nanstd(dataArray1{nGroup,nState},0,2)) /...
                     sqrt(size(dataArray1{nGroup,nState},2)-1);                                     
        errorbar(nState, mean2plot,err2plot,'Color',globalParams.stateColors{nState}); hold on
        ylim(ylimits1)
        box off
        set(gca,'TickDir','out')
        ylabel(yAxName1)      
        xlim(xlimits)
        
        subplot(1,3,2)
        mean2plot = []; err2plot = [];
        mean2plot = nanmean(dataArray2{nGroup,nState},2);
        err2plot  = (nanstd(dataArray2{nGroup,nState},0,2)) /...
                     sqrt(size(dataArray2{nGroup,nState},2)-1);                                     
        errorbar(nState, mean2plot,err2plot,'Color',globalParams.stateColors{nState}); hold on
        ylim(ylimits2)
        box off
        set(gca,'TickDir','out')
        ylabel(yAxName2)      
        xlim(xlimits)
        
        subplot(1,3,3)
        mean2plot = []; err2plot = [];
        mean2plot = nanmean(dataArray3{nGroup,nState},2);
        err2plot  = (nanstd(dataArray3{nGroup,nState},0,2)) /...
                     sqrt(size(dataArray3{nGroup,nState},2)-1);                                     
        errorbar(nState, mean2plot,err2plot,'Color',globalParams.stateColors{nState}); hold on
        ylim(ylimits3)
        box off
        set(gca,'TickDir','out')
        ylabel(yAxName3)      
        xlim(xlimits)    
    end
        
    subplot(1,3,1)
    plot(1:3,[dataArray1{nGroup,1}' ...
              dataArray1{nGroup,2}' ...
              dataArray1{nGroup,3}'], ...
              '-', 'Color', [0.8 0.8 0.8]);
    subplot(1,3,2)
    plot(1:3,[dataArray2{nGroup,1}' ...
              dataArray2{nGroup,2}' ...
              dataArray2{nGroup,3}'], ...
              '-', 'Color', [0.8 0.8 0.8]);
    subplot(1,3,3)
    plot(1:3,[dataArray3{nGroup,1}' ...
              dataArray3{nGroup,2}' ...
              dataArray3{nGroup,3}'], ...
              '-', 'Color', [0.8 0.8 0.8]);
    
    data1 = []; data2 = []; data3 = [];    
    data1   = [dataArray1{nGroup,1}'; ...
               dataArray1{nGroup,2}'; ...
               dataArray1{nGroup,3}'];              
    xIdx1   = [ones(size(dataArray1{nGroup,1},2),1);   ...
               ones(size(dataArray1{nGroup,2},2),1)*2; ...
               ones(size(dataArray1{nGroup,3},2),1)*3];
    subj1   = [[1:numel(dataArray1{nGroup,1})]';   ...
               [1:numel(dataArray1{nGroup,2})]'; ...
               [1:numel(dataArray1{nGroup,3})]'];
     
    data2   = [dataArray2{nGroup,1}'; ...
               dataArray2{nGroup,2}'; ...
               dataArray2{nGroup,3}'];              
    xIdx2   = [ones(size(dataArray2{nGroup,1},2),1);   ...
               ones(size(dataArray2{nGroup,2},2),1)*2; ...
               ones(size(dataArray2{nGroup,3},2),1)*3];
    subj2   = [[1:numel(dataArray2{nGroup,1})]';   ...
               [1:numel(dataArray2{nGroup,2})]'; ...
               [1:numel(dataArray2{nGroup,3})]'];
          
    data3   = [dataArray3{nGroup,1}'; ...
               dataArray3{nGroup,2}'; ...
               dataArray3{nGroup,3}'];              
    xIdx3   = [ones(size(dataArray3{nGroup,1},2),1);   ...
               ones(size(dataArray3{nGroup,2},2),1)*2; ...
               ones(size(dataArray3{nGroup,3},2),1)*3];
    subj3   = [[1:numel(dataArray3{nGroup,1})]';   ...
               [1:numel(dataArray3{nGroup,2})]'; ...
               [1:numel(dataArray3{nGroup,3})]'];       
              
    X = []; pState = []; fState = [];
    X(:,1) = data1; 
    X(:,2) = xIdx1; 
    X(:,3) = subj1;
    [pState,~,fState] = RMAOV1(X,.05);   
    subplot(1,3,1)
    text(0.6,ylimits1(2)*.95,['state p= ' num2str(pState)])
    text(0.6,ylimits1(2)*.90,['F= ' num2str(fState)])
    set(gca,'XTickLabel',{'s1','s2','s3'});     
    if contains(groupOrder{nGroup},'D2')
        title('indirect')    
    elseif contains(groupOrder{nGroup},'D1')
        title('direct')
    else
    end
    
    X = []; pState = []; fState = [];
    X(:,1) = data2; 
    X(:,2) = xIdx2; 
    X(:,3) = subj2;
    [pState,~,fState] = RMAOV1(X,.05);   
    subplot(1,3,2)
    text(0.6,ylimits2(2)*.95,['state p= ' num2str(pState)])
    text(0.6,ylimits2(2)*.90,['F= ' num2str(fState)])
    set(gca,'XTickLabel',{'s1','s2','s3'});     
    if contains(groupOrder{nGroup},'D2')
        title('indirect')    
    elseif contains(groupOrder{nGroup},'D1')
        title('direct')
    else
    end
    
    X = []; pState = []; fState = [];
    X(:,1) = data3; 
    X(:,2) = xIdx3; 
    X(:,3) = subj3;
    [pState,~,fState] = RMAOV1(X,.05);   
    subplot(1,3,3)
    text(0.6,ylimits3(2)*.95,['state p= ' num2str(pState)])
    text(0.6,ylimits3(2)*.90,['F= ' num2str(fState)])
    set(gca,'XTickLabel',{'s1','s2','s3'});     
    if contains(groupOrder{nGroup},'D2')
        title('indirect')    
    elseif contains(groupOrder{nGroup},'D1')
        title('direct')
    else
    end
    
end          

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% clear all vars except source data for plotting
clearvars -except groupOrder lgDMSD1 lgDMSD2 ...
                  nTrialsSUBSELECT_OFF nTrialsSUBSELECT_ON ...
                  nTrialsALL_OFF nTrialsALL_ON ...
                  nTrialsMATCHMICE_OFF nTrialsMATCHMICE_ON ...
                  viewAng_OFF_xState viewAng_ON_xState ...
                  viewAngL_xState viewAngR_xState ...
                  xPos_OFF_xState xPos_ON_xState ...
                  xPosL_xState xPosR_xState ... 
                  yVelOFF_xState yVelON_xState...                 
                  vaDevOFF_xState vaDevON_xState ...
                  excTravelOFF_xState excTravelON_xState ...
                  distanceOFF_xState distanceON_xState    