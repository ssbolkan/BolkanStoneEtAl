%% Script to generate ExtData Fig 8 data plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% cross task comparison of delta bias with stats text
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
cd(globalParams.dataPath)
load('OffOn_TasksOrGroup_all.mat', 'lgDMSD2_aoe'  ,'lgDMSD2_nd'  ,'lgDMSD2_pc'  , ...
                                   'lgDMSD1_aoe'  ,'lgDMSD1_nd'  ,'lgDMSD1_pc'  , ...
                                   'lgDMSCtrl_aoe','lgDMSCtrl_nd','lgDMSCtrl_pc', ...
                                   'lgNACD2_aoe'  ,'lgNACD1_aoe' );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% analyzes log files for opto induced choice bias
% (1) cue statistics of aoe (log.currMaze==10) and pc (log.currMaze==7) pre-selected to match
% (2) all analyses subselect trial blocks with performance >60% and trials with excessTravel<0.1 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

groupOrder = {'DMSD2_aoe'  ,'DMSD2_nd'  ,'DMSD2_pc'  , ...
              'DMSD1_aoe'  ,'DMSD1_nd'  ,'DMSD1_pc'  , ...
              'DMSCtrl_aoe','DMSCtrl_nd','DMSCtrl_pc', ...
              'NACD2_aoe'  ,'NACD1_aoe' };
counter = 1;
for groupNum = 1:numel(groupOrder)
    clear lg lgClean 
    lg          = eval(['lg' groupOrder{groupNum}]); 
       
    % remove field not of same length
    if isfield(lg,'keyFrameLabels')
        lg = rmfield(lg,'keyFrameLabels');
    else
    end
    
    % drops blocks with perf <0.6 and excessTravel>0.1
    [lgClean,~] = selectLogSubset(lg,[],[],[],[],0,0.6); 
    
    % subselect mazes in aoe and perm cues with matching cue stats
    % did this for DMS/NAC comparison cause could seem odd if DMS data
    % differed in repeated plots of the data
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
        
    % overall bias, delta bias, log of total off/on trials here
    ipsiCon = 1; 
    tempPerf = []; tempSumm = [];
    [tempPerf, tempSumm] = xMousePerfBias(lgClean, ipsiCon);
    biasOFF_xMouse{counter}(1,:)    = [tempPerf(:).contraBiasOFF]*100;   
    biasON_xMouse{counter}(1,:)     = [tempPerf(:).contraBiasON]*100;     
    biasDelta_xMouse{counter}(1,:)  = [tempPerf(:).contraBiasDiff]*100;         
    nTrialsOFF_xGroup(counter)       = sum([tempPerf(:).nTrialsOFF]);
    nTrialsON_xGroup(counter)        = sum([tempPerf(:).nTrialsON]);   
    counter = counter + 1; 
end

%% 
clearvars -except lgDMSD2_aoe lgDMSD1_aoe lgDMSCtrl_aoe ...
                  lgDMSD2_nd lgDMSD1_nd lgDMSCtrl_nd ...
                  lgDMSD2_pc lgDMSD1_pc lgDMSCtrl_pc ...
                  lgNACD2_aoe lgNACD1_aoe ...
                  groupOrder biasOFF_xMouse biasON_xMouse biasDelta_xMouse ...
                  nTrialsOFF_xGroup nTrialsON_xGroup
                        
%% panel 8b-f %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
% plot off vs on bias and delta bias for indirect/direct/no opsin DMS
% during evidence accumulation 
displayStats = 1; % 0 to suppress display of 
                  % (1) 1-way ANOVA of group on delta bias
                  % (2) post-hoc ranksum unpaired comparison of delta bias
                  %     across tasks                  
      
% update these 2(3 for title) variables to determine group to plot comparisons, each must much
% the order defined in groupOrder to be accurate.%%%%%%%%%%%%%%%%%%%%%%%%%%
titlestr    = 'DMSindirect'; % DMSindirect DMSdirect DMSnoOpsin indirect direct
names2plot  = {'DMSD2_aoe'  ,'DMSD2_nd'   ,'DMSD2_pc'};   %panel 8B - indirect DMS: aoe, ctrl#1, ctrl#2
%             {'DMSD2_aoe'  ,'DMSD2_nd'   ,'DMSD2_pc'};   %panel 8B - indirect DMS: aoe, ctrl#1, ctrl#2
%             {'DMSD1_aoe'  ,'DMSD1_nd'   ,'DMSD1_pc'};   %panel 8C - direct DMS: aoe, ctrl#1, ctrl#2
%             {'DMSCtrl_aoe','DMSCtrl_nd' ,'DMSCtrl_pc'}; %panel 8D - no opsin DMS: aoe, ctrl#1, ctrl#2
%             {'DMSD2_aoe'  ,'NACD2_aoe'};                %panel 8E - AoE indirect pathway: DMS vs NAc
%             {'DMSD1_aoe'  ,'NACD1_aoe'};                %panel 8F - AoE direct pathway: DMS vs NAc
groups2plot = [1 2 3]; %panel 8B - indirect DMS: aoe, ctrl#1, ctrl#2
%             [1 2 3]; %panel 8B - indirect DMS: aoe, ctrl#1, ctrl#2
%             [4 5 6]; %panel 8C - direct DMS: aoe, ctrl#1, ctrl#2
%             [7 8 9]; %panel 8D - no opsin DMS: aoe, ctrl#1, ctrl#2
%             [1 10];  %panel 8E - AoE indirect pathway: DMS vs NAc
%             [4 11];  %panel 8F - AoE direct pathway: DMS vs NAc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if numel(groups2plot) == 3
    xNames     = {'AoE';'ctrl#1';'ctrl#2'};
    barColor   = {'k','m','c'};     
    colorIdx   = {[0.8 0.8 0.8],[0.99 0.7 0.99], [0.7 0.99 0.99]}; 
    xLimits    = [0.5 3.5];  
    yLimits    = [-90 90];
    catMarker  = {'x','x', 'x'};
else
    xNames     = {'DMS';'NAc'};
    barColor   = {'k','k'};     
    colorIdx   = {[0.8 0.8 0.8],[0.8 0.8 0.8]}; 
    xLimits    = [0.5 2.5];
    yLimits    = [-90 90];
    catMarker  = {'x','x'};
end

dataArrayDelta = [];
dataDelta      = []; 
dataArrayDelta = biasDelta_xMouse; 
if numel(groups2plot) == 3
    dataDelta      = [biasDelta_xMouse{groups2plot(1)}'; ...
                      biasDelta_xMouse{groups2plot(2)}'; ...
                      biasDelta_xMouse{groups2plot(3)}'];
    xIdx           = [ones(size(biasDelta_xMouse{groups2plot(1)}))'*1;   ...
                      ones(size(biasDelta_xMouse{groups2plot(2)}))'*2;  ...
                      ones(size(biasDelta_xMouse{groups2plot(3)}))'*3];
else
    dataDelta      = [biasDelta_xMouse{groups2plot(1)}'; ...
                      biasDelta_xMouse{groups2plot(2)}'];
    xIdx           = [ones(size(biasDelta_xMouse{groups2plot(1)}))'*1;   ...
                      ones(size(biasDelta_xMouse{groups2plot(2)}))'*2];
end


figure
plot(xLimits(1):xLimits(2),zeros(1,numel(xLimits(1):xLimits(2))),'-','Color',[0.6 0.6 0.6]); hold on

plotSpread(dataDelta,'categoryIdx',xIdx,'distributionIdx',xIdx,...
    'categoryColors',colorIdx,'categoryMarkers',catMarker); hold on

nTask = find(contains(groupOrder,names2plot{1}));
h = errorbar(1,nanmean(dataArrayDelta{nTask}),nanstd(dataArrayDelta{nTask})./sqrt(numel(dataArrayDelta{nTask})-1), ...
    '.-','linewidth',4,'color',barColor{1},'markersize',0.1); h.CapSize = 14; hold on
text(0.5,yLimits(1)*.85, [num2str(numel(dataArrayDelta{nTask})) ' mice' ],'FontSize',6)
text(0.5,yLimits(1)*.9, [num2str(nTrialsOFF_xGroup(nTask)) ' off' ],'FontSize',6)
text(0.5,yLimits(1)*.95, [num2str(nTrialsON_xGroup(nTask)) ' on' ],'FontSize',6)

nTask = find(contains(groupOrder,names2plot{2}));
h = errorbar(2,nanmean(dataArrayDelta{nTask}),nanstd(dataArrayDelta{nTask})./sqrt(numel(dataArrayDelta{nTask})-1), ...
    '.-','linewidth',4,'color',barColor{2},'markersize',0.1); h.CapSize = 14; hold on
text(1.5,yLimits(1)*.85, [num2str(numel(dataArrayDelta{nTask})) ' mice' ],'FontSize',6)
text(1.5,yLimits(1)*.9, [num2str(nTrialsOFF_xGroup(nTask)) ' off' ],'FontSize',6)
text(1.5,yLimits(1)*.95, [num2str(nTrialsON_xGroup(nTask)) ' on' ],'FontSize',6)

if numel(groups2plot) == 3
    nTask = find(contains(groupOrder,names2plot{3}));
    h = errorbar(3,nanmean(dataArrayDelta{nTask}),nanstd(dataArrayDelta{nTask})./sqrt(numel(dataArrayDelta{nTask})-1), ...
        '.-','linewidth',4,'color',barColor{3},'markersize',0.1); h.CapSize = 14; hold on
    text(2.5,yLimits(1)*.85, [num2str(numel(dataArrayDelta{nTask})) ' mice' ],'FontSize',6)
    text(2.5,yLimits(1)*.9, [num2str(nTrialsOFF_xGroup(nTask)) ' off' ],'FontSize',6)
    text(2.5,yLimits(1)*.95, [num2str(nTrialsON_xGroup(nTask)) ' on' ],'FontSize',6)
else
end

ylabel('delta bias (%, on-off)','fontsize',12)
box off
ylim(yLimits)
xlim(xLimits)
set(gca,'TickDir','out')
set(gca,'xtick',1:numel(names2plot),'xticklabel',xNames)
rotateXLabels(gca,45)
title(titlestr)
            
if displayStats == 1
    % unpaired stats for subplot 4 - delta bias across groups
    groupVec   = categorical(xIdx);
    [p,t,stats] = anova1(dataDelta,groupVec,'off');
    text(xLimits(2)-.5,yLimits(1)*0.85,['Group p= ' num2str(p)],'FontSize',6)
    text(xLimits(2)-.5,yLimits(1)*0.9,['df= ' num2str(stats.df)],'FontSize',6)
    text(xLimits(2)-.5,yLimits(1)*0.95,['F= ' num2str((t{2,5}))],'FontSize',6)
    
    if numel(groups2plot) == 3
        [g1Vg2_p, ~ , g1Vg2_stats] = ranksum(dataArrayDelta{groups2plot(1)}, ...
                                     dataArrayDelta{groups2plot(2)},'method','approximate');
        plot(1:2,[76,76],'k-')
        text(0.75,78,['p= ' num2str(g1Vg2_p)],'FontSize',6)
        text(1.25,78,['z= ' num2str(g1Vg2_stats.zval)],'FontSize',6)
        text(1.75,78,['df= ' num2str(numel(dataArrayDelta{groups2plot(1)})+ ...
             numel(dataArrayDelta{groups2plot(2)})-2)],'FontSize',6)  
         
         
        [g1Vg3_p, ~ , g1Vg3_stats] = ranksum(dataArrayDelta{groups2plot(1)}, ...
                                     dataArrayDelta{groups2plot(3)},'method','approximate');                                 
        plot(1:3,[80,80,80],'k-')
        text(1.5,82,['p= ' num2str(g1Vg3_p)],'FontSize',6)
        text(2,82,['z= ' num2str(g1Vg3_stats.zval)],'FontSize',6)
        text(2.5,82,['df= ' num2str(numel(dataArrayDelta{groups2plot(1)})+ ...
             numel(dataArrayDelta{groups2plot(3)})-2)],'FontSize',6)                        
         
        [g2Vg3_p, ~ , g2Vg3_stats] = ranksum(dataArrayDelta{groups2plot(2)}, ...
                                     dataArrayDelta{groups2plot(3)},'method','approximate');
        plot(2:3,[72,72],'k-')
        text(2,74,['p= ' num2str(g2Vg3_p)],'FontSize',6)
        text(2.5,74,['z= ' num2str(g2Vg3_stats.zval)],'FontSize',6)
        text(3,74,['df= ' num2str(numel(dataArrayDelta{groups2plot(2)})+ ...
             numel(dataArrayDelta{groups2plot(3)})-2)],'FontSize',6)  
         
    else
        
        [g1Vg2_p, ~ , g1Vg2_stats] = ranksum(dataArrayDelta{groups2plot(1)}, ...
                                     dataArrayDelta{groups2plot(2)},'method','approximate');         
        plot(1:2,[80,80],'k-')
        text(0.75,82,['p= ' num2str(g1Vg2_p)],'FontSize',6)
        text(1.25,82,['z= ' num2str(g1Vg2_stats.zval)],'FontSize',6)
        text(1.75,82,['df= ' num2str(numel(dataArrayDelta{groups2plot(1)})+ ...
             numel(dataArrayDelta{groups2plot(2)})-2)],'FontSize',6)                          
    end                                                           
else
end              
     

%% 
clearvars -except lgDMSD2_aoe lgDMSD1_aoe lgDMSCtrl_aoe ...
                  lgDMSD2_nd lgDMSD1_nd lgDMSCtrl_nd ...
                  lgDMSD2_pc lgDMSD1_pc lgDMSCtrl_pc ...
                  lgNACD2_aoe lgNACD1_aoe ...
                  groupOrder biasOFF_xMouse biasON_xMouse biasDelta_xMouse ...
                  nTrialsOFF_xGroup nTrialsON_xGroup
                        
%% panel 8b-f %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
% plot off vs on bias and delta bias for indirect/direct/no opsin DMS
% during evidence accumulation 
displayStats = 1; % 0 to suppress display of 
                  % (1) 1-way ANOVA of group on delta bias
                  % (2) post-hoc ranksum unpaired comparison of delta bias
                  %     across tasks                  
      
% update these 2(3 for title) variables to determine group to plot comparisons, each must much
% the order defined in groupOrder to be accurate.%%%%%%%%%%%%%%%%%%%%%%%%%%
titlestr    = 'DMSdirect'; % DMSindirect DMSdirect DMSnoOpsin indirect direct
names2plot  = {'DMSD1_aoe'  ,'DMSD1_nd'   ,'DMSD1_pc'};   %panel 8C - direct DMS: aoe, ctrl#1, ctrl#2
%             {'DMSD2_aoe'  ,'DMSD2_nd'   ,'DMSD2_pc'};   %panel 8B - indirect DMS: aoe, ctrl#1, ctrl#2
%             {'DMSD1_aoe'  ,'DMSD1_nd'   ,'DMSD1_pc'};   %panel 8C - direct DMS: aoe, ctrl#1, ctrl#2
%             {'DMSCtrl_aoe','DMSCtrl_nd' ,'DMSCtrl_pc'}; %panel 8D - no opsin DMS: aoe, ctrl#1, ctrl#2
%             {'DMSD2_aoe'  ,'NACD2_aoe'};                %panel 8E - AoE indirect pathway: DMS vs NAc
%             {'DMSD1_aoe'  ,'NACD1_aoe'};                %panel 8F - AoE direct pathway: DMS vs NAc
groups2plot = [4 5 6]; %panel 8C - direct DMS: aoe, ctrl#1, ctrl#2
%             [1 2 3]; %panel 8B - indirect DMS: aoe, ctrl#1, ctrl#2
%             [4 5 6]; %panel 8C - direct DMS: aoe, ctrl#1, ctrl#2
%             [7 8 9]; %panel 8D - no opsin DMS: aoe, ctrl#1, ctrl#2
%             [1 10];  %panel 8E - AoE indirect pathway: DMS vs NAc
%             [4 11];  %panel 8F - AoE direct pathway: DMS vs NAc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if numel(groups2plot) == 3
    xNames     = {'AoE';'ctrl#1';'ctrl#2'};
    barColor   = {'k','m','c'};     
    colorIdx   = {[0.8 0.8 0.8],[0.99 0.7 0.99], [0.7 0.99 0.99]}; 
    xLimits    = [0.5 3.5];  
    yLimits    = [-90 90];
    catMarker  = {'x','x', 'x'};
else
    xNames     = {'DMS';'NAc'};
    barColor   = {'k','k'};     
    colorIdx   = {[0.8 0.8 0.8],[0.8 0.8 0.8]}; 
    xLimits    = [0.5 2.5];
    yLimits    = [-90 90];
    catMarker  = {'x','x'};
end

dataArrayDelta = [];
dataDelta      = []; 
dataArrayDelta = biasDelta_xMouse; 
if numel(groups2plot) == 3
    dataDelta      = [biasDelta_xMouse{groups2plot(1)}'; ...
                      biasDelta_xMouse{groups2plot(2)}'; ...
                      biasDelta_xMouse{groups2plot(3)}'];
    xIdx           = [ones(size(biasDelta_xMouse{groups2plot(1)}))'*1;   ...
                      ones(size(biasDelta_xMouse{groups2plot(2)}))'*2;  ...
                      ones(size(biasDelta_xMouse{groups2plot(3)}))'*3];
else
    dataDelta      = [biasDelta_xMouse{groups2plot(1)}'; ...
                      biasDelta_xMouse{groups2plot(2)}'];
    xIdx           = [ones(size(biasDelta_xMouse{groups2plot(1)}))'*1;   ...
                      ones(size(biasDelta_xMouse{groups2plot(2)}))'*2];
end


figure
plot(xLimits(1):xLimits(2),zeros(1,numel(xLimits(1):xLimits(2))),'-','Color',[0.6 0.6 0.6]); hold on

plotSpread(dataDelta,'categoryIdx',xIdx,'distributionIdx',xIdx,...
    'categoryColors',colorIdx,'categoryMarkers',catMarker); hold on

nTask = find(contains(groupOrder,names2plot{1}));
h = errorbar(1,nanmean(dataArrayDelta{nTask}),nanstd(dataArrayDelta{nTask})./sqrt(numel(dataArrayDelta{nTask})-1), ...
    '.-','linewidth',4,'color',barColor{1},'markersize',0.1); h.CapSize = 14; hold on
text(0.5,yLimits(1)*.85, [num2str(numel(dataArrayDelta{nTask})) ' mice' ],'FontSize',6)
text(0.5,yLimits(1)*.9, [num2str(nTrialsOFF_xGroup(nTask)) ' off' ],'FontSize',6)
text(0.5,yLimits(1)*.95, [num2str(nTrialsON_xGroup(nTask)) ' on' ],'FontSize',6)

nTask = find(contains(groupOrder,names2plot{2}));
h = errorbar(2,nanmean(dataArrayDelta{nTask}),nanstd(dataArrayDelta{nTask})./sqrt(numel(dataArrayDelta{nTask})-1), ...
    '.-','linewidth',4,'color',barColor{2},'markersize',0.1); h.CapSize = 14; hold on
text(1.5,yLimits(1)*.85, [num2str(numel(dataArrayDelta{nTask})) ' mice' ],'FontSize',6)
text(1.5,yLimits(1)*.9, [num2str(nTrialsOFF_xGroup(nTask)) ' off' ],'FontSize',6)
text(1.5,yLimits(1)*.95, [num2str(nTrialsON_xGroup(nTask)) ' on' ],'FontSize',6)

if numel(groups2plot) == 3
    nTask = find(contains(groupOrder,names2plot{3}));
    h = errorbar(3,nanmean(dataArrayDelta{nTask}),nanstd(dataArrayDelta{nTask})./sqrt(numel(dataArrayDelta{nTask})-1), ...
        '.-','linewidth',4,'color',barColor{3},'markersize',0.1); h.CapSize = 14; hold on
    text(2.5,yLimits(1)*.85, [num2str(numel(dataArrayDelta{nTask})) ' mice' ],'FontSize',6)
    text(2.5,yLimits(1)*.9, [num2str(nTrialsOFF_xGroup(nTask)) ' off' ],'FontSize',6)
    text(2.5,yLimits(1)*.95, [num2str(nTrialsON_xGroup(nTask)) ' on' ],'FontSize',6)
else
end

ylabel('delta bias (%, on-off)','fontsize',12)
box off
ylim(yLimits)
xlim(xLimits)
set(gca,'TickDir','out')
set(gca,'xtick',1:numel(names2plot),'xticklabel',xNames)
rotateXLabels(gca,45)
title(titlestr)
            
if displayStats == 1
    % unpaired stats for subplot 4 - delta bias across groups
    groupVec   = categorical(xIdx);
    [p,t,stats] = anova1(dataDelta,groupVec,'off');
    text(xLimits(2)-.5,yLimits(1)*0.85,['Group p= ' num2str(p)],'FontSize',6)
    text(xLimits(2)-.5,yLimits(1)*0.9,['df= ' num2str(stats.df)],'FontSize',6)
    text(xLimits(2)-.5,yLimits(1)*0.95,['F= ' num2str((t{2,5}))],'FontSize',6)
    
    if numel(groups2plot) == 3
        [g1Vg2_p, ~ , g1Vg2_stats] = ranksum(dataArrayDelta{groups2plot(1)}, ...
                                     dataArrayDelta{groups2plot(2)},'method','approximate');
        plot(1:2,[76,76],'k-')
        text(0.75,78,['p= ' num2str(g1Vg2_p)],'FontSize',6)
        text(1.25,78,['z= ' num2str(g1Vg2_stats.zval)],'FontSize',6)
        text(1.75,78,['df= ' num2str(numel(dataArrayDelta{groups2plot(1)})+ ...
             numel(dataArrayDelta{groups2plot(2)})-2)],'FontSize',6)  
         
         
        [g1Vg3_p, ~ , g1Vg3_stats] = ranksum(dataArrayDelta{groups2plot(1)}, ...
                                     dataArrayDelta{groups2plot(3)},'method','approximate');                                 
        plot(1:3,[80,80,80],'k-')
        text(1.5,82,['p= ' num2str(g1Vg3_p)],'FontSize',6)
        text(2,82,['z= ' num2str(g1Vg3_stats.zval)],'FontSize',6)
        text(2.5,82,['df= ' num2str(numel(dataArrayDelta{groups2plot(1)})+ ...
             numel(dataArrayDelta{groups2plot(3)})-2)],'FontSize',6)                        
         
        [g2Vg3_p, ~ , g2Vg3_stats] = ranksum(dataArrayDelta{groups2plot(2)}, ...
                                     dataArrayDelta{groups2plot(3)},'method','approximate');
        plot(2:3,[72,72],'k-')
        text(2,74,['p= ' num2str(g2Vg3_p)],'FontSize',6)
        text(2.5,74,['z= ' num2str(g2Vg3_stats.zval)],'FontSize',6)
        text(3,74,['df= ' num2str(numel(dataArrayDelta{groups2plot(2)})+ ...
             numel(dataArrayDelta{groups2plot(3)})-2)],'FontSize',6)  
         
    else
        
        [g1Vg2_p, ~ , g1Vg2_stats] = ranksum(dataArrayDelta{groups2plot(1)}, ...
                                     dataArrayDelta{groups2plot(2)},'method','approximate');         
        plot(1:2,[80,80],'k-')
        text(0.75,82,['p= ' num2str(g1Vg2_p)],'FontSize',6)
        text(1.25,82,['z= ' num2str(g1Vg2_stats.zval)],'FontSize',6)
        text(1.75,82,['df= ' num2str(numel(dataArrayDelta{groups2plot(1)})+ ...
             numel(dataArrayDelta{groups2plot(2)})-2)],'FontSize',6)                          
    end                                                           
else
end     

%% 
clearvars -except lgDMSD2_aoe lgDMSD1_aoe lgDMSCtrl_aoe ...
                  lgDMSD2_nd lgDMSD1_nd lgDMSCtrl_nd ...
                  lgDMSD2_pc lgDMSD1_pc lgDMSCtrl_pc ...
                  lgNACD2_aoe lgNACD1_aoe ...
                  groupOrder biasOFF_xMouse biasON_xMouse biasDelta_xMouse ...
                  nTrialsOFF_xGroup nTrialsON_xGroup
                        
%% panel 8b-f %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
% plot off vs on bias and delta bias for indirect/direct/no opsin DMS
% during evidence accumulation 
displayStats = 1; % 0 to suppress display of 
                  % (1) 1-way ANOVA of group on delta bias
                  % (2) post-hoc ranksum unpaired comparison of delta bias
                  %     across tasks                  
      
% update these 2(3 for title) variables to determine group to plot comparisons, each must much
% the order defined in groupOrder to be accurate.%%%%%%%%%%%%%%%%%%%%%%%%%%
titlestr    = 'DMSnoOpsin'; % DMSindirect DMSdirect DMSnoOpsin indirect direct
names2plot  = {'DMSCtrl_aoe','DMSCtrl_nd' ,'DMSCtrl_pc'}; %panel 8D - no opsin DMS: aoe, ctrl#1, ctrl#2
%             {'DMSD2_aoe'  ,'DMSD2_nd'   ,'DMSD2_pc'};   %panel 8B - indirect DMS: aoe, ctrl#1, ctrl#2
%             {'DMSD1_aoe'  ,'DMSD1_nd'   ,'DMSD1_pc'};   %panel 8C - direct DMS: aoe, ctrl#1, ctrl#2
%             {'DMSCtrl_aoe','DMSCtrl_nd' ,'DMSCtrl_pc'}; %panel 8D - no opsin DMS: aoe, ctrl#1, ctrl#2
%             {'DMSD2_aoe'  ,'NACD2_aoe'};                %panel 8E - AoE indirect pathway: DMS vs NAc
%             {'DMSD1_aoe'  ,'NACD1_aoe'};                %panel 8F - AoE direct pathway: DMS vs NAc
groups2plot = [7 8 9]; %panel 8D - no opsin DMS: aoe, ctrl#1, ctrl#2 
%             [1 2 3]; %panel 8B - indirect DMS: aoe, ctrl#1, ctrl#2
%             [4 5 6]; %panel 8C - direct DMS: aoe, ctrl#1, ctrl#2
%             [7 8 9]; %panel 8D - no opsin DMS: aoe, ctrl#1, ctrl#2
%             [1 10];  %panel 8E - AoE indirect pathway: DMS vs NAc
%             [4 11];  %panel 8F - AoE direct pathway: DMS vs NAc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if numel(groups2plot) == 3
    xNames     = {'AoE';'ctrl#1';'ctrl#2'};
    barColor   = {'k','m','c'};     
    colorIdx   = {[0.8 0.8 0.8],[0.99 0.7 0.99], [0.7 0.99 0.99]}; 
    xLimits    = [0.5 3.5];  
    yLimits    = [-90 90];
    catMarker  = {'x','x', 'x'};
else
    xNames     = {'DMS';'NAc'};
    barColor   = {'k','k'};     
    colorIdx   = {[0.8 0.8 0.8],[0.8 0.8 0.8]}; 
    xLimits    = [0.5 2.5];
    yLimits    = [-90 90];
    catMarker  = {'x','x'};
end

dataArrayDelta = [];
dataDelta      = []; 
dataArrayDelta = biasDelta_xMouse; 
if numel(groups2plot) == 3
    dataDelta      = [biasDelta_xMouse{groups2plot(1)}'; ...
                      biasDelta_xMouse{groups2plot(2)}'; ...
                      biasDelta_xMouse{groups2plot(3)}'];
    xIdx           = [ones(size(biasDelta_xMouse{groups2plot(1)}))'*1;   ...
                      ones(size(biasDelta_xMouse{groups2plot(2)}))'*2;  ...
                      ones(size(biasDelta_xMouse{groups2plot(3)}))'*3];
else
    dataDelta      = [biasDelta_xMouse{groups2plot(1)}'; ...
                      biasDelta_xMouse{groups2plot(2)}'];
    xIdx           = [ones(size(biasDelta_xMouse{groups2plot(1)}))'*1;   ...
                      ones(size(biasDelta_xMouse{groups2plot(2)}))'*2];
end


figure
plot(xLimits(1):xLimits(2),zeros(1,numel(xLimits(1):xLimits(2))),'-','Color',[0.6 0.6 0.6]); hold on

plotSpread(dataDelta,'categoryIdx',xIdx,'distributionIdx',xIdx,...
    'categoryColors',colorIdx,'categoryMarkers',catMarker); hold on

nTask = find(contains(groupOrder,names2plot{1}));
h = errorbar(1,nanmean(dataArrayDelta{nTask}),nanstd(dataArrayDelta{nTask})./sqrt(numel(dataArrayDelta{nTask})-1), ...
    '.-','linewidth',4,'color',barColor{1},'markersize',0.1); h.CapSize = 14; hold on
text(0.5,yLimits(1)*.85, [num2str(numel(dataArrayDelta{nTask})) ' mice' ],'FontSize',6)
text(0.5,yLimits(1)*.9, [num2str(nTrialsOFF_xGroup(nTask)) ' off' ],'FontSize',6)
text(0.5,yLimits(1)*.95, [num2str(nTrialsON_xGroup(nTask)) ' on' ],'FontSize',6)

nTask = find(contains(groupOrder,names2plot{2}));
h = errorbar(2,nanmean(dataArrayDelta{nTask}),nanstd(dataArrayDelta{nTask})./sqrt(numel(dataArrayDelta{nTask})-1), ...
    '.-','linewidth',4,'color',barColor{2},'markersize',0.1); h.CapSize = 14; hold on
text(1.5,yLimits(1)*.85, [num2str(numel(dataArrayDelta{nTask})) ' mice' ],'FontSize',6)
text(1.5,yLimits(1)*.9, [num2str(nTrialsOFF_xGroup(nTask)) ' off' ],'FontSize',6)
text(1.5,yLimits(1)*.95, [num2str(nTrialsON_xGroup(nTask)) ' on' ],'FontSize',6)

if numel(groups2plot) == 3
    nTask = find(contains(groupOrder,names2plot{3}));
    h = errorbar(3,nanmean(dataArrayDelta{nTask}),nanstd(dataArrayDelta{nTask})./sqrt(numel(dataArrayDelta{nTask})-1), ...
        '.-','linewidth',4,'color',barColor{3},'markersize',0.1); h.CapSize = 14; hold on
    text(2.5,yLimits(1)*.85, [num2str(numel(dataArrayDelta{nTask})) ' mice' ],'FontSize',6)
    text(2.5,yLimits(1)*.9, [num2str(nTrialsOFF_xGroup(nTask)) ' off' ],'FontSize',6)
    text(2.5,yLimits(1)*.95, [num2str(nTrialsON_xGroup(nTask)) ' on' ],'FontSize',6)
else
end

ylabel('delta bias (%, on-off)','fontsize',12)
box off
ylim(yLimits)
xlim(xLimits)
set(gca,'TickDir','out')
set(gca,'xtick',1:numel(names2plot),'xticklabel',xNames)
rotateXLabels(gca,45)
title(titlestr)
            
if displayStats == 1
    % unpaired stats for subplot 4 - delta bias across groups
    groupVec   = categorical(xIdx);
    [p,t,stats] = anova1(dataDelta,groupVec,'off');
    text(xLimits(2)-.5,yLimits(1)*0.85,['Group p= ' num2str(p)],'FontSize',6)
    text(xLimits(2)-.5,yLimits(1)*0.9,['df= ' num2str(stats.df)],'FontSize',6)
    text(xLimits(2)-.5,yLimits(1)*0.95,['F= ' num2str((t{2,5}))],'FontSize',6)
    
    if numel(groups2plot) == 3
        [g1Vg2_p, ~ , g1Vg2_stats] = ranksum(dataArrayDelta{groups2plot(1)}, ...
                                     dataArrayDelta{groups2plot(2)},'method','approximate');
        plot(1:2,[76,76],'k-')
        text(0.75,78,['p= ' num2str(g1Vg2_p)],'FontSize',6)
        text(1.25,78,['z= ' num2str(g1Vg2_stats.zval)],'FontSize',6)
        text(1.75,78,['df= ' num2str(numel(dataArrayDelta{groups2plot(1)})+ ...
             numel(dataArrayDelta{groups2plot(2)})-2)],'FontSize',6)  
         
         
        [g1Vg3_p, ~ , g1Vg3_stats] = ranksum(dataArrayDelta{groups2plot(1)}, ...
                                     dataArrayDelta{groups2plot(3)},'method','approximate');                                 
        plot(1:3,[80,80,80],'k-')
        text(1.5,82,['p= ' num2str(g1Vg3_p)],'FontSize',6)
        text(2,82,['z= ' num2str(g1Vg3_stats.zval)],'FontSize',6)
        text(2.5,82,['df= ' num2str(numel(dataArrayDelta{groups2plot(1)})+ ...
             numel(dataArrayDelta{groups2plot(3)})-2)],'FontSize',6)                        
         
        [g2Vg3_p, ~ , g2Vg3_stats] = ranksum(dataArrayDelta{groups2plot(2)}, ...
                                     dataArrayDelta{groups2plot(3)},'method','approximate');
        plot(2:3,[72,72],'k-')
        text(2,74,['p= ' num2str(g2Vg3_p)],'FontSize',6)
        text(2.5,74,['z= ' num2str(g2Vg3_stats.zval)],'FontSize',6)
        text(3,74,['df= ' num2str(numel(dataArrayDelta{groups2plot(2)})+ ...
             numel(dataArrayDelta{groups2plot(3)})-2)],'FontSize',6)  
         
    else
        
        [g1Vg2_p, ~ , g1Vg2_stats] = ranksum(dataArrayDelta{groups2plot(1)}, ...
                                     dataArrayDelta{groups2plot(2)},'method','approximate');         
        plot(1:2,[80,80],'k-')
        text(0.75,82,['p= ' num2str(g1Vg2_p)],'FontSize',6)
        text(1.25,82,['z= ' num2str(g1Vg2_stats.zval)],'FontSize',6)
        text(1.75,82,['df= ' num2str(numel(dataArrayDelta{groups2plot(1)})+ ...
             numel(dataArrayDelta{groups2plot(2)})-2)],'FontSize',6)                          
    end                                                           
else
end     

%% 
clearvars -except lgDMSD2_aoe lgDMSD1_aoe lgDMSCtrl_aoe ...
                  lgDMSD2_nd lgDMSD1_nd lgDMSCtrl_nd ...
                  lgDMSD2_pc lgDMSD1_pc lgDMSCtrl_pc ...
                  lgNACD2_aoe lgNACD1_aoe ...
                  groupOrder biasOFF_xMouse biasON_xMouse biasDelta_xMouse ...
                  nTrialsOFF_xGroup nTrialsON_xGroup
                        
%% panel 8b-f %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
% plot off vs on bias and delta bias for indirect/direct/no opsin DMS
% during evidence accumulation 
displayStats = 1; % 0 to suppress display of 
                  % (1) 1-way ANOVA of group on delta bias
                  % (2) post-hoc ranksum unpaired comparison of delta bias
                  %     across tasks                  
      
% update these 2(3 for title) variables to determine group to plot comparisons, each must much
% the order defined in groupOrder to be accurate.%%%%%%%%%%%%%%%%%%%%%%%%%%
titlestr    = 'indirect'; % DMSindirect DMSdirect DMSnoOpsin indirect direct
names2plot  = {'DMSD2_aoe'  ,'NACD2_aoe'};                %panel 8E - AoE indirect pathway: DMS vs NAc
%             {'DMSD2_aoe'  ,'DMSD2_nd'   ,'DMSD2_pc'};   %panel 8B - indirect DMS: aoe, ctrl#1, ctrl#2
%             {'DMSD1_aoe'  ,'DMSD1_nd'   ,'DMSD1_pc'};   %panel 8C - direct DMS: aoe, ctrl#1, ctrl#2
%             {'DMSCtrl_aoe','DMSCtrl_nd' ,'DMSCtrl_pc'}; %panel 8D - no opsin DMS: aoe, ctrl#1, ctrl#2
%             {'DMSD2_aoe'  ,'NACD2_aoe'};                %panel 8E - AoE indirect pathway: DMS vs NAc
%             {'DMSD1_aoe'  ,'NACD1_aoe'};                %panel 8F - AoE direct pathway: DMS vs NAc
groups2plot = [1 10];  %panel 8E - AoE indirect pathway: DMS vs NAc
%             [1 2 3]; %panel 8B - indirect DMS: aoe, ctrl#1, ctrl#2
%             [4 5 6]; %panel 8C - direct DMS: aoe, ctrl#1, ctrl#2
%             [7 8 9]; %panel 8D - no opsin DMS: aoe, ctrl#1, ctrl#2
%             [1 10];  %panel 8E - AoE indirect pathway: DMS vs NAc
%             [4 11];  %panel 8F - AoE direct pathway: DMS vs NAc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if numel(groups2plot) == 3
    xNames     = {'AoE';'ctrl#1';'ctrl#2'};
    barColor   = {'k','m','c'};     
    colorIdx   = {[0.8 0.8 0.8],[0.99 0.7 0.99], [0.7 0.99 0.99]}; 
    xLimits    = [0.5 3.5];  
    yLimits    = [-90 90];
    catMarker  = {'x','x', 'x'};
else
    xNames     = {'DMS';'NAc'};
    barColor   = {'k','k'};     
    colorIdx   = {[0.8 0.8 0.8],[0.8 0.8 0.8]}; 
    xLimits    = [0.5 2.5];
    yLimits    = [-90 90];
    catMarker  = {'x','x'};
end

dataArrayDelta = [];
dataDelta      = []; 
dataArrayDelta = biasDelta_xMouse; 
if numel(groups2plot) == 3
    dataDelta      = [biasDelta_xMouse{groups2plot(1)}'; ...
                      biasDelta_xMouse{groups2plot(2)}'; ...
                      biasDelta_xMouse{groups2plot(3)}'];
    xIdx           = [ones(size(biasDelta_xMouse{groups2plot(1)}))'*1;   ...
                      ones(size(biasDelta_xMouse{groups2plot(2)}))'*2;  ...
                      ones(size(biasDelta_xMouse{groups2plot(3)}))'*3];
else
    dataDelta      = [biasDelta_xMouse{groups2plot(1)}'; ...
                      biasDelta_xMouse{groups2plot(2)}'];
    xIdx           = [ones(size(biasDelta_xMouse{groups2plot(1)}))'*1;   ...
                      ones(size(biasDelta_xMouse{groups2plot(2)}))'*2];
end


figure
plot(xLimits(1):xLimits(2),zeros(1,numel(xLimits(1):xLimits(2))),'-','Color',[0.6 0.6 0.6]); hold on

plotSpread(dataDelta,'categoryIdx',xIdx,'distributionIdx',xIdx,...
    'categoryColors',colorIdx,'categoryMarkers',catMarker); hold on

nTask = find(contains(groupOrder,names2plot{1}));
h = errorbar(1,nanmean(dataArrayDelta{nTask}),nanstd(dataArrayDelta{nTask})./sqrt(numel(dataArrayDelta{nTask})-1), ...
    '.-','linewidth',4,'color',barColor{1},'markersize',0.1); h.CapSize = 14; hold on
text(0.5,yLimits(1)*.85, [num2str(numel(dataArrayDelta{nTask})) ' mice' ],'FontSize',6)
text(0.5,yLimits(1)*.9, [num2str(nTrialsOFF_xGroup(nTask)) ' off' ],'FontSize',6)
text(0.5,yLimits(1)*.95, [num2str(nTrialsON_xGroup(nTask)) ' on' ],'FontSize',6)

nTask = find(contains(groupOrder,names2plot{2}));
h = errorbar(2,nanmean(dataArrayDelta{nTask}),nanstd(dataArrayDelta{nTask})./sqrt(numel(dataArrayDelta{nTask})-1), ...
    '.-','linewidth',4,'color',barColor{2},'markersize',0.1); h.CapSize = 14; hold on
text(1.5,yLimits(1)*.85, [num2str(numel(dataArrayDelta{nTask})) ' mice' ],'FontSize',6)
text(1.5,yLimits(1)*.9, [num2str(nTrialsOFF_xGroup(nTask)) ' off' ],'FontSize',6)
text(1.5,yLimits(1)*.95, [num2str(nTrialsON_xGroup(nTask)) ' on' ],'FontSize',6)

if numel(groups2plot) == 3
    nTask = find(contains(groupOrder,names2plot{3}));
    h = errorbar(3,nanmean(dataArrayDelta{nTask}),nanstd(dataArrayDelta{nTask})./sqrt(numel(dataArrayDelta{nTask})-1), ...
        '.-','linewidth',4,'color',barColor{3},'markersize',0.1); h.CapSize = 14; hold on
    text(2.5,yLimits(1)*.85, [num2str(numel(dataArrayDelta{nTask})) ' mice' ],'FontSize',6)
    text(2.5,yLimits(1)*.9, [num2str(nTrialsOFF_xGroup(nTask)) ' off' ],'FontSize',6)
    text(2.5,yLimits(1)*.95, [num2str(nTrialsON_xGroup(nTask)) ' on' ],'FontSize',6)
else
end

ylabel('delta bias (%, on-off)','fontsize',12)
box off
ylim(yLimits)
xlim(xLimits)
set(gca,'TickDir','out')
set(gca,'xtick',1:numel(names2plot),'xticklabel',xNames)
rotateXLabels(gca,45)
title(titlestr)
            
if displayStats == 1
    % unpaired stats for subplot 4 - delta bias across groups
    groupVec   = categorical(xIdx);
    [p,t,stats] = anova1(dataDelta,groupVec,'off');
    text(xLimits(2)-.5,yLimits(1)*0.85,['Group p= ' num2str(p)],'FontSize',6)
    text(xLimits(2)-.5,yLimits(1)*0.9,['df= ' num2str(stats.df)],'FontSize',6)
    text(xLimits(2)-.5,yLimits(1)*0.95,['F= ' num2str((t{2,5}))],'FontSize',6)
    
    if numel(groups2plot) == 3
        [g1Vg2_p, ~ , g1Vg2_stats] = ranksum(dataArrayDelta{groups2plot(1)}, ...
                                     dataArrayDelta{groups2plot(2)},'method','approximate');
        plot(1:2,[76,76],'k-')
        text(0.75,78,['p= ' num2str(g1Vg2_p)],'FontSize',6)
        text(1.25,78,['z= ' num2str(g1Vg2_stats.zval)],'FontSize',6)
        text(1.75,78,['df= ' num2str(numel(dataArrayDelta{groups2plot(1)})+ ...
             numel(dataArrayDelta{groups2plot(2)})-2)],'FontSize',6)  
         
         
        [g1Vg3_p, ~ , g1Vg3_stats] = ranksum(dataArrayDelta{groups2plot(1)}, ...
                                     dataArrayDelta{groups2plot(3)},'method','approximate');                                 
        plot(1:3,[80,80,80],'k-')
        text(1.5,82,['p= ' num2str(g1Vg3_p)],'FontSize',6)
        text(2,82,['z= ' num2str(g1Vg3_stats.zval)],'FontSize',6)
        text(2.5,82,['df= ' num2str(numel(dataArrayDelta{groups2plot(1)})+ ...
             numel(dataArrayDelta{groups2plot(3)})-2)],'FontSize',6)                        
         
        [g2Vg3_p, ~ , g2Vg3_stats] = ranksum(dataArrayDelta{groups2plot(2)}, ...
                                     dataArrayDelta{groups2plot(3)},'method','approximate');
        plot(2:3,[72,72],'k-')
        text(2,74,['p= ' num2str(g2Vg3_p)],'FontSize',6)
        text(2.5,74,['z= ' num2str(g2Vg3_stats.zval)],'FontSize',6)
        text(3,74,['df= ' num2str(numel(dataArrayDelta{groups2plot(2)})+ ...
             numel(dataArrayDelta{groups2plot(3)})-2)],'FontSize',6)  
         
    else
        
        [g1Vg2_p, ~ , g1Vg2_stats] = ranksum(dataArrayDelta{groups2plot(1)}, ...
                                     dataArrayDelta{groups2plot(2)},'method','approximate');         
        plot(1:2,[80,80],'k-')
        text(0.75,82,['p= ' num2str(g1Vg2_p)],'FontSize',6)
        text(1.25,82,['z= ' num2str(g1Vg2_stats.zval)],'FontSize',6)
        text(1.75,82,['df= ' num2str(numel(dataArrayDelta{groups2plot(1)})+ ...
             numel(dataArrayDelta{groups2plot(2)})-2)],'FontSize',6)                          
    end                                                           
else
end     

%% clear all vars except source data for plots
clearvars -except lgDMSD2_aoe lgDMSD1_aoe lgDMSCtrl_aoe ...
                  lgDMSD2_nd lgDMSD1_nd lgDMSCtrl_nd ...
                  lgDMSD2_pc lgDMSD1_pc lgDMSCtrl_pc ...
                  lgNACD2_aoe lgNACD1_aoe ...
                  groupOrder biasOFF_xMouse biasON_xMouse biasDelta_xMouse ...
                  nTrialsOFF_xGroup nTrialsON_xGroup
                        
%% panel 8b-f %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
% plot off vs on bias and delta bias for indirect/direct/no opsin DMS
% during evidence accumulation 
displayStats = 1; % 0 to suppress display of 
                  % (1) 1-way ANOVA of group on delta bias
                  % (2) post-hoc ranksum unpaired comparison of delta bias
                  %     across tasks                  
      
% update these 2(3 for title) variables to determine group to plot comparisons, each must much
% the order defined in groupOrder to be accurate.%%%%%%%%%%%%%%%%%%%%%%%%%%
titlestr    = 'direct'; % DMSindirect DMSdirect DMSnoOpsin indirect direct
names2plot  = {'DMSD1_aoe'  ,'NACD1_aoe'};                %panel 8F - AoE direct pathway: DMS vs NAc
%             {'DMSD2_aoe'  ,'DMSD2_nd'   ,'DMSD2_pc'};   %panel 8B - indirect DMS: aoe, ctrl#1, ctrl#2
%             {'DMSD1_aoe'  ,'DMSD1_nd'   ,'DMSD1_pc'};   %panel 8C - direct DMS: aoe, ctrl#1, ctrl#2
%             {'DMSCtrl_aoe','DMSCtrl_nd' ,'DMSCtrl_pc'}; %panel 8D - no opsin DMS: aoe, ctrl#1, ctrl#2
%             {'DMSD2_aoe'  ,'NACD2_aoe'};                %panel 8E - AoE indirect pathway: DMS vs NAc
%             {'DMSD1_aoe'  ,'NACD1_aoe'};                %panel 8F - AoE direct pathway: DMS vs NAc
groups2plot = [4 11];  %panel 8F - AoE direct pathway: DMS vs NAc
%             [1 2 3]; %panel 8B - indirect DMS: aoe, ctrl#1, ctrl#2
%             [4 5 6]; %panel 8C - direct DMS: aoe, ctrl#1, ctrl#2
%             [7 8 9]; %panel 8D - no opsin DMS: aoe, ctrl#1, ctrl#2
%             [1 10];  %panel 8E - AoE indirect pathway: DMS vs NAc
%             [4 11];  %panel 8F - AoE direct pathway: DMS vs NAc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if numel(groups2plot) == 3
    xNames     = {'AoE';'ctrl#1';'ctrl#2'};
    barColor   = {'k','m','c'};     
    colorIdx   = {[0.8 0.8 0.8],[0.99 0.7 0.99], [0.7 0.99 0.99]}; 
    xLimits    = [0.5 3.5];  
    yLimits    = [-90 90];
    catMarker  = {'x','x', 'x'};
else
    xNames     = {'DMS';'NAc'};
    barColor   = {'k','k'};     
    colorIdx   = {[0.8 0.8 0.8],[0.8 0.8 0.8]}; 
    xLimits    = [0.5 2.5];
    yLimits    = [-90 90];
    catMarker  = {'x','x'};
end

dataArrayDelta = [];
dataDelta      = []; 
dataArrayDelta = biasDelta_xMouse; 
if numel(groups2plot) == 3
    dataDelta      = [biasDelta_xMouse{groups2plot(1)}'; ...
                      biasDelta_xMouse{groups2plot(2)}'; ...
                      biasDelta_xMouse{groups2plot(3)}'];
    xIdx           = [ones(size(biasDelta_xMouse{groups2plot(1)}))'*1;   ...
                      ones(size(biasDelta_xMouse{groups2plot(2)}))'*2;  ...
                      ones(size(biasDelta_xMouse{groups2plot(3)}))'*3];
else
    dataDelta      = [biasDelta_xMouse{groups2plot(1)}'; ...
                      biasDelta_xMouse{groups2plot(2)}'];
    xIdx           = [ones(size(biasDelta_xMouse{groups2plot(1)}))'*1;   ...
                      ones(size(biasDelta_xMouse{groups2plot(2)}))'*2];
end


figure
plot(xLimits(1):xLimits(2),zeros(1,numel(xLimits(1):xLimits(2))),'-','Color',[0.6 0.6 0.6]); hold on

plotSpread(dataDelta,'categoryIdx',xIdx,'distributionIdx',xIdx,...
    'categoryColors',colorIdx,'categoryMarkers',catMarker); hold on

nTask = find(contains(groupOrder,names2plot{1}));
h = errorbar(1,nanmean(dataArrayDelta{nTask}),nanstd(dataArrayDelta{nTask})./sqrt(numel(dataArrayDelta{nTask})-1), ...
    '.-','linewidth',4,'color',barColor{1},'markersize',0.1); h.CapSize = 14; hold on
text(0.5,yLimits(1)*.85, [num2str(numel(dataArrayDelta{nTask})) ' mice' ],'FontSize',6)
text(0.5,yLimits(1)*.9, [num2str(nTrialsOFF_xGroup(nTask)) ' off' ],'FontSize',6)
text(0.5,yLimits(1)*.95, [num2str(nTrialsON_xGroup(nTask)) ' on' ],'FontSize',6)

nTask = find(contains(groupOrder,names2plot{2}));
h = errorbar(2,nanmean(dataArrayDelta{nTask}),nanstd(dataArrayDelta{nTask})./sqrt(numel(dataArrayDelta{nTask})-1), ...
    '.-','linewidth',4,'color',barColor{2},'markersize',0.1); h.CapSize = 14; hold on
text(1.5,yLimits(1)*.85, [num2str(numel(dataArrayDelta{nTask})) ' mice' ],'FontSize',6)
text(1.5,yLimits(1)*.9, [num2str(nTrialsOFF_xGroup(nTask)) ' off' ],'FontSize',6)
text(1.5,yLimits(1)*.95, [num2str(nTrialsON_xGroup(nTask)) ' on' ],'FontSize',6)

if numel(groups2plot) == 3
    nTask = find(contains(groupOrder,names2plot{3}));
    h = errorbar(3,nanmean(dataArrayDelta{nTask}),nanstd(dataArrayDelta{nTask})./sqrt(numel(dataArrayDelta{nTask})-1), ...
        '.-','linewidth',4,'color',barColor{3},'markersize',0.1); h.CapSize = 14; hold on
    text(2.5,yLimits(1)*.85, [num2str(numel(dataArrayDelta{nTask})) ' mice' ],'FontSize',6)
    text(2.5,yLimits(1)*.9, [num2str(nTrialsOFF_xGroup(nTask)) ' off' ],'FontSize',6)
    text(2.5,yLimits(1)*.95, [num2str(nTrialsON_xGroup(nTask)) ' on' ],'FontSize',6)
else
end

ylabel('delta bias (%, on-off)','fontsize',12)
box off
ylim(yLimits)
xlim(xLimits)
set(gca,'TickDir','out')
set(gca,'xtick',1:numel(names2plot),'xticklabel',xNames)
rotateXLabels(gca,45)
title(titlestr)
            
if displayStats == 1
    % unpaired stats for subplot 4 - delta bias across groups
    groupVec   = categorical(xIdx);
    [p,t,stats] = anova1(dataDelta,groupVec,'off');
    text(xLimits(2)-.5,yLimits(1)*0.85,['Group p= ' num2str(p)],'FontSize',6)
    text(xLimits(2)-.5,yLimits(1)*0.9,['df= ' num2str(stats.df)],'FontSize',6)
    text(xLimits(2)-.5,yLimits(1)*0.95,['F= ' num2str((t{2,5}))],'FontSize',6)
    
    if numel(groups2plot) == 3
        [g1Vg2_p, ~ , g1Vg2_stats] = ranksum(dataArrayDelta{groups2plot(1)}, ...
                                     dataArrayDelta{groups2plot(2)},'method','approximate');
        plot(1:2,[76,76],'k-')
        text(0.75,78,['p= ' num2str(g1Vg2_p)],'FontSize',6)
        text(1.25,78,['z= ' num2str(g1Vg2_stats.zval)],'FontSize',6)
        text(1.75,78,['df= ' num2str(numel(dataArrayDelta{groups2plot(1)})+ ...
             numel(dataArrayDelta{groups2plot(2)})-2)],'FontSize',6)  
         
         
        [g1Vg3_p, ~ , g1Vg3_stats] = ranksum(dataArrayDelta{groups2plot(1)}, ...
                                     dataArrayDelta{groups2plot(3)},'method','approximate');                                 
        plot(1:3,[80,80,80],'k-')
        text(1.5,82,['p= ' num2str(g1Vg3_p)],'FontSize',6)
        text(2,82,['z= ' num2str(g1Vg3_stats.zval)],'FontSize',6)
        text(2.5,82,['df= ' num2str(numel(dataArrayDelta{groups2plot(1)})+ ...
             numel(dataArrayDelta{groups2plot(3)})-2)],'FontSize',6)                        
         
        [g2Vg3_p, ~ , g2Vg3_stats] = ranksum(dataArrayDelta{groups2plot(2)}, ...
                                     dataArrayDelta{groups2plot(3)},'method','approximate');
        plot(2:3,[72,72],'k-')
        text(2,74,['p= ' num2str(g2Vg3_p)],'FontSize',6)
        text(2.5,74,['z= ' num2str(g2Vg3_stats.zval)],'FontSize',6)
        text(3,74,['df= ' num2str(numel(dataArrayDelta{groups2plot(2)})+ ...
             numel(dataArrayDelta{groups2plot(3)})-2)],'FontSize',6)  
         
    else
        
        [g1Vg2_p, ~ , g1Vg2_stats] = ranksum(dataArrayDelta{groups2plot(1)}, ...
                                     dataArrayDelta{groups2plot(2)},'method','approximate');         
        plot(1:2,[80,80],'k-')
        text(0.75,82,['p= ' num2str(g1Vg2_p)],'FontSize',6)
        text(1.25,82,['z= ' num2str(g1Vg2_stats.zval)],'FontSize',6)
        text(1.75,82,['df= ' num2str(numel(dataArrayDelta{groups2plot(1)})+ ...
             numel(dataArrayDelta{groups2plot(2)})-2)],'FontSize',6)                          
    end                                                           
else
end     


%% clear all vars except source data for plots
clearvars -except lgDMSD2_aoe lgDMSD1_aoe lgDMSCtrl_aoe ...
                  lgDMSD2_nd lgDMSD1_nd lgDMSCtrl_nd ...
                  lgDMSD2_pc lgDMSD1_pc lgDMSCtrl_pc ...
                  lgNACD2_aoe lgNACD1_aoe ...
                  groupOrder biasOFF_xMouse biasON_xMouse biasDelta_xMouse ...
                  nTrialsOFF_xGroup nTrialsON_xGroup