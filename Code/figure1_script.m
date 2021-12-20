%% Script to generate Fig 1 data plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTE THAT Code/Utility subfolder contains compiled versions of
% SUPPORTING MEX FUNCTIONS from 'Matlab Offline Files SDK' AND
% 'TankMouseVR' REPOS. If a local recompile is necessary... 
% (1) CLONE REPOS AND ADD TO YOUR MATLAB PATH
% (2) MAKE SURE YOU HAVE A COMPATABLE VISUAL STUDIO (OR OTHER) C++ COMPILER FOR
% YOUR VERSION OF MATLAB
% (2) RUN compile_utilities_pipeline() TO COMPILE MEX FILES IN
% 'TankMouseVR-master'
% (3) SET YOUR DIRECTORY TO '[your local path]\Matlab Offline Files
% SDK\mexPlex\' AND RUN build_and_verify_mexPlex() (ALSO SEE 'HOW TO BUILD
% mexPlex.txt' IN THE SAME FOLDER FOR MORE DETAIL)
% NOTE: mexPlex only necessary if loading .plx data, which is not the case
% if using spikeTS, laserTS, and unit waveforms saved in spikeDataAll.mat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plots PSTH and raster for example units in panel 1B-C
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
cd(globalParams.dataPath)
load('runningWheelSpikes.mat')

xlims   = []; xlims   = [-5 10];
binsize = []; binsize = 0.25;
plotWaveform = 1; % set to 1 for seperate fig of unit waveform placed as inset

% plot panel 1B
exampleFile = spikeDataAll(17);
subsamp = []; subsamp = 1;
psthRasterExample_metaFile(exampleFile,xlims,binsize,subsamp,[],plotWaveform)
% plot panel 1C
exampleFile = spikeDataAll(70);
subsamp = []; subsamp = 2;
psthRasterExample_metaFile(exampleFile,xlims,binsize,subsamp,[],plotWaveform)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load virtual corridor behavior lg files for fig 1G-J
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
cd(globalParams.dataPath)
load('virtualCorridor.mat')

%% analyzes the lgs for y-velocity, x-position, view angle and distance
groupOrder = {'D2_vc', 'D1_vc', 'Ctrl_vc'};
counter = 1;
for groupNum = 1:numel(groupOrder)
    
    lg = eval(['lgDMS' groupOrder{groupNum}]);
    
    %  lg pre-process and data selection for velocity, x-position and view angle    
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
    
    % create seperate lgCleanOFF and lgCleanON, which inclues only on "cue" or 0-200cm
    removeMouse = [];
    removeMouse = lgClean.laserOFF ==1;
    lgCleanOFF = structfun(@(x) x(removeMouse),lgClean,'UniformOutput',false);
    removeMouse = [];
    removeMouse = lgClean.laserCue == 1;
    lgCleanON  = structfun(@(x) x(removeMouse),lgClean,'UniformOutput',false);
    
    % panel 1G - calculation of Yvelocity in bins of 25cm for 300cm stem
    binIDs = []; binIDs = 11:25:311;
    for nBin = 1:numel(binIDs)-1
        yVelOFF{counter}(nBin,:)   = xMouseSpeedXY(lgCleanOFF,[binIDs(nBin) binIDs(nBin+1)], 2);
        yVelON{counter}(nBin,:)    = xMouseSpeedXY(lgCleanON,[binIDs(nBin) binIDs(nBin+1)], 2);
    end
    
    % panel 1H - calculation of x-position in bins of 25cm for 300cm stem
    avgXposON  = sampleViewAngleVsY_average(lgCleanON.pos2,[11 311],25);
    avgXposOFF = sampleViewAngleVsY_average(lgCleanOFF.pos2,[11 311],25);
    
    mouseNums       = [];
    mouseNums       = unique(lgClean.mouseID);
    for nMouse = 1:numel(mouseNums)
        xPosOFF{counter}(:,nMouse) = nanmean(avgXposOFF(:,lgCleanOFF.mouseID == mouseNums(nMouse)),2);
        xPosON{counter}(:,nMouse)  = nanmean(avgXposON(:,lgCleanON.mouseID == mouseNums(nMouse)),2);
    end
    
    % panel 1I - calculation of view angle in bins of 25cm for 300cm stem
    avgViewAngleON  = sampleViewAngleVsY_average(lgCleanON.pos,[11 311],25);
    avgViewAngleOFF = sampleViewAngleVsY_average(lgCleanOFF.pos,[11 311],25);
    
    mouseNums       = [];
    mouseNums       = unique(lgClean.mouseID);
    for nMouse = 1:numel(mouseNums)
        viewAngOFF{counter}(:,nMouse) = nanmean(avgViewAngleOFF(:,lgCleanOFF.mouseID == mouseNums(nMouse)),2);
        viewAngON{counter}(:,nMouse)  = nanmean(avgViewAngleON(:,lgCleanON.mouseID == mouseNums(nMouse)),2);
    end
         
    % panel 1J - calculation of average total distance per trial - no trial or mouse subselection
    clear lg distance mouseNums distOFF distON
    lg = eval(['lgDMS' groupOrder{groupNum}]);
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
    
    counter = counter +1;
end

%% clear all except the key processed variables to plot (source data)
clearvars -except yVelOFF yVelON viewAngOFF viewAngON ...
                  xPosOFF xPosON distanceOFF distanceON ...
                  lgDMSCtrl_vc lgDMSD1_vc lgDMSD2_vc groupOrder
              
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plotting the data for fig 1G-J
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plot panel 1G : y-velocity             
figure
binIDs    = 25:25:300;
offOffset = -7.5;
onOffset  = -12.5;
for nGroup = 1:numel(groupOrder)
    subplot(1,3,nGroup)
    mean2plot = []; err2plot = [];
    mean2plot = nanmean(yVelOFF{nGroup},2);  
    err2plot  = (nanstd(yVelOFF{nGroup},0,2)) / sqrt(size(yVelOFF{nGroup},2)-1);  
    errorbar(mean2plot, binIDs' + offOffset ,err2plot,'horizontal','k-'); hold on
    mean2plot = []; err2plot = [];
    mean2plot = nanmean(yVelON{nGroup},2);  
    err2plot  = (nanstd(yVelON{nGroup},0,2)) / sqrt(size(yVelON{nGroup},2)-1);  
    errorbar(mean2plot, binIDs' + onOffset ,err2plot,'horizontal','-','Color',[.4 1 .4]); hold on   
    xlim([20 90]) 
    box off
    set(gca,'TickDir','out')
    if contains(groupOrder(nGroup),'D2')
        title('indirect')    
    elseif contains(groupOrder(nGroup),'D1')
        title('direct')
    elseif contains(groupOrder(nGroup),'Ctrl')
        title('no opsin')
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

%% plot panel 1H - x position          
figure
binIDs    = 25:25:300;
offOffset = -7.5;
onOffset  = -12.5;
for nGroup = 1:numel(groupOrder)
    subplot(1,3,nGroup)
    mean2plot = []; err2plot = [];
    mean2plot = nanmean(xPosOFF{nGroup},2);  
    err2plot  = (nanstd(xPosOFF{nGroup},0,2)) / sqrt(size(xPosOFF{nGroup},2)-1);  
    errorbar(mean2plot, binIDs' + offOffset ,err2plot,'horizontal','k-'); hold on
    mean2plot = []; err2plot = [];
    mean2plot = nanmean(xPosON{nGroup},2);  
    err2plot  = (nanstd(xPosON{nGroup},0,2)) / sqrt(size(xPosON{nGroup},2)-1);  
    errorbar(mean2plot, binIDs' + onOffset ,err2plot,'horizontal','-','Color',[.4 1 .4]); hold on   
    xlim([-2.9 2.9]) 
    box off
    set(gca,'TickDir','out')
    if contains(groupOrder(nGroup),'D2')
        title('indirect')    
    elseif contains(groupOrder(nGroup),'D1')
        title('direct')
    elseif contains(groupOrder(nGroup),'Ctrl')
        title('no opsin')
    else
    end
    if nGroup == 1
        ylabel('y position (cm)')
    else
        set(gca,'YTickLabel',[])
    end
    
    if nGroup == 2
        xlabel('x position (cm)')
    end
end


%% plot panel 1I - view angle              
figure
binIDs    = 25:25:300;
offOffset = -7.5;
onOffset  = -12.5;
for nGroup = 1:numel(groupOrder)
    subplot(1,3,nGroup)
    mean2plot = []; err2plot = [];
    mean2plot = nanmean(viewAngOFF{nGroup},2);  
    err2plot  = (nanstd(viewAngOFF{nGroup},0,2)) / sqrt(size(viewAngOFF{nGroup},2)-1);  
    errorbar(mean2plot, binIDs' + offOffset ,err2plot,'horizontal','k-'); hold on
    mean2plot = []; err2plot = [];
    mean2plot = nanmean(viewAngON{nGroup},2);  
    err2plot  = (nanstd(viewAngON{nGroup},0,2)) / sqrt(size(viewAngON{nGroup},2)-1);  
    errorbar(mean2plot, binIDs' + onOffset ,err2plot,'horizontal','-','Color',[.4 1 .4]); hold on   
    xlim([-20 20]) 
    box off
    set(gca,'TickDir','out')
    if contains(groupOrder(nGroup),'D2')
        title('indirect')    
    elseif contains(groupOrder(nGroup),'D1')
        title('direct')
    elseif contains(groupOrder(nGroup),'Ctrl')
        title('no opsin')
    else
    end
    if nGroup == 1
        ylabel('y position (cm)')
    else
        set(gca,'YTickLabel',[])
    end
    
    if nGroup == 2
        xlabel('view angle (deg.)')
    end
end

%% plot panel 1J - distance

figure
for nGroup = 1:numel(groupOrder)
    subplot(1,3,nGroup)
    mean2plot = []; err2plot = [];
    mean2plot = nanmean(distanceOFF{nGroup},2);  
    err2plot  = (nanstd(distanceOFF{nGroup},0,2)) / sqrt(size(distanceOFF{nGroup},2)-1);
    errorbar(1, mean2plot ,err2plot,'Color', 'k','LineStyle','none'); hold on
    mean2plot = []; err2plot = [];
    mean2plot = nanmean(distanceON{nGroup},2);  
    err2plot  = (nanstd(distanceON{nGroup},0,2)) / sqrt(size(distanceON{nGroup},2)-1);
    errorbar(2, mean2plot ,err2plot,'Color', [0.5 1 0.5],'LineStyle','none'); hold on
    plot(1:2, [distanceOFF{nGroup}' distanceON{nGroup}'],'Color', [0.8 0.8 0.8]); hold on
    xlim([0.5 2.5])
    ylim([340 700])
    box off; set(gca,'TickDir','out');
    set(gca,'xtick',1:2,'xticklabel',{'off';'on'})
    
     if contains(groupOrder(nGroup),'D2')
        title('indirect')    
    elseif contains(groupOrder(nGroup),'D1')
        title('direct')
    elseif contains(groupOrder(nGroup),'Ctrl')
        title('no opsin')
    else
    end
    if nGroup == 1
        ylabel('y position (cm)')
    else
        set(gca,'YTickLabel',[])
    end
    
    if nGroup == 2
        xlabel('view angle (deg.)')
    end
end

    
 %% scrap   
    




