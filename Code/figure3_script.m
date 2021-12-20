%% Script to generate Fig 3 data plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% cross task comparison of indirect and direct pathway inhibition DMS and NAc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
cd(globalParams.dataPath)
load('OffOn_TasksOrGroup_all.mat')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% analyzes log files for opto induced choice bias
% (1) cue statistics of aoe (log.currMaze==10) and pc (log.currMaze==7) pre-selected to match
% (2) all analyses subselect trial blocks with performance >60% and trials with excessTravel<0.1 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

groupOrder = {'DMSD2_aoe','DMSD1_aoe','DMSCtrl_aoe', ...
              'DMSD2_nd' ,'DMSD1_nd' ,'DMSCtrl_nd', ...
              'DMSD2_pc' ,'DMSD1_pc' ,'DMSCtrl_pc', ...
              'NACD2_aoe','NACD1_aoe', 'NACCtrl_aoe'};
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

%% clear all except the key processed variables to plot (source data)
clearvars -except lgDMSD2_aoe lgDMSD1_aoe lgDMSCtrl_aoe ...
                  lgDMSD2_nd lgDMSD1_nd lgDMSCtrl_nd ...
                  lgDMSD2_pc lgDMSD1_pc lgDMSCtrl_pc ...
                  lgNACD2_aoe lgNACD1_aoe lgNACCtrl_aoe ...
                  groupOrder biasOFF_xMouse biasON_xMouse biasDelta_xMouse ...
                  nTrialsOFF_xGroup nTrialsON_xGroup

%% panel 2D-E %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
% plot off vs on bias and delta bias for indirect/direct/no opsin DMS
% during evidence accumulation 
displayStats = 1; % 0 to suppress display of 
                  % (1) 1-way ANOVA of group on delta bias
                  % (2) post-hoc ranksum unpaired comparison of delta bias
                  %     with no opsin group                  
      
% update these 4 variables to determine group to plot and color                  
names2plot  = {'DMSD2_aoe','DMSD1_aoe','DMSCtrl_aoe'};                  
groups2plot = [1 2 3];  
barColor    = 'k';     
colorIdx    = {[0.8 0.8 0.8],[0.8 0.8 0.8],[0.8 0.8 0.8]}; % [0.8 0.8 0.8],[0.99 0.7 0.99], [0.7 0.99 0.99]

figure
subplot(1,4,1)
dataArrayOFF = []; dataArrayON = []; dataArrayDelta = [];
data = []; dataDelta = []; xIdx = [];
dataArrayOFF   = biasOFF_xMouse;
dataArrayON    = biasON_xMouse; 
dataArrayDelta = biasDelta_xMouse;

data(:,1)    = [biasOFF_xMouse{1}'; ...
                biasOFF_xMouse{2}'; ...
                biasOFF_xMouse{3}'];
data(:,2)    = [biasON_xMouse{1}'; ...
                biasON_xMouse{2}'; ...
                biasON_xMouse{3}'];     
dataDelta    = [biasDelta_xMouse{1}'; ...
                biasDelta_xMouse{2}'; ...
                biasDelta_xMouse{3}'];   
xIdx         = [ones(size(biasOFF_xMouse{1}))'*groups2plot(1);   ...
                ones(size(biasOFF_xMouse{2}))'*groups2plot(2); ...
                ones(size(biasOFF_xMouse{3}))'*groups2plot(3)];
catMarker   = {'x','x', 'x'};

% subplot 4 - delta bias indirect, direct, no opsin
subplot(1,4,4);
plot([.5:3.5],zeros(1,4),'-','Color',[0.6 0.6 0.6]); hold on
plotSpread(dataDelta,'categoryIdx',xIdx,'distributionIdx',xIdx,...
    'categoryColors',colorIdx,'categoryMarkers',catMarker); 

nTask = find(contains(groupOrder,names2plot{1}));
h = errorbar(1,nanmean(dataArrayDelta{nTask}),nanstd(dataArrayDelta{nTask})./sqrt(numel(dataArrayDelta{nTask})-1), ...
    '.-','linewidth',4,'color',barColor,'markersize',0.1); h.CapSize = 14;

nTask = find(contains(groupOrder,names2plot{2}));
h = errorbar(2,nanmean(dataArrayDelta{nTask}),nanstd(dataArrayDelta{nTask})./sqrt(numel(dataArrayDelta{nTask})-1), ...
    '.-','linewidth',4,'color',barColor,'markersize',0.1); h.CapSize = 14;

nTask = find(contains(groupOrder,names2plot{3}));
h = errorbar(3,nanmean(dataArrayDelta{nTask}),nanstd(dataArrayDelta{nTask})./sqrt(numel(dataArrayDelta{nTask})-1), ...
    '.-','linewidth',4,'color',barColor,'markersize',0.1); h.CapSize = 14;

ylabel('delta contra. bias (%, on-off)','fontsize',12)
xlim([.5 3.5])
ylim([-100 100])
box off
set(gca,'TickDir','out')
set(gca,'xtick',1:3,'xticklabel',{'indirect';'direct';'no opsin'})
rotateXLabels(gca,45)
            
% subplot 1 - off vs on bias indirect pathway                      
nTask = find(contains(groupOrder,names2plot{1}));
subplot(1,4,1)
plot([.5:2.5],zeros(1,3),'-','Color',[0.6 0.6 0.6]); hold on
plot(1:2, data(xIdx==find(nTask==groups2plot),:),'-','Color', [0.9 0.9 0.9]); hold on
h = errorbar(1,nanmean(dataArrayOFF{nTask}),nanstd(dataArrayOFF{nTask})./sqrt(numel(dataArrayOFF{nTask})-1), ...
    '.-','linewidth',4,'color','k','markersize',0.1); h.CapSize = 16;
h = errorbar(2,nanmean(dataArrayON{nTask}),nanstd(dataArrayON{nTask})./sqrt(numel(dataArrayON{nTask})-1), ...
    '.-','linewidth',4,'color','g','markersize',0.1); h.CapSize = 16;
set(gca,'xtick',1:2,'xticklabel',{'off';'on'})
ylabel('% bias (contra-ipsi)','fontsize',12)
xlim([.5 2.5])
ylim([-80 80])
box off
set(gca,'TickDir','out')
title('indirect')
text(0.5,-55, [num2str(numel(dataArrayOFF{nTask})) ' mice' ])
text(0.5,-65, [num2str(nTrialsOFF_xGroup(nTask)) ' off' ])
text(0.5,-75, [num2str(nTrialsON_xGroup(nTask)) ' on' ])

% subplot 2 - off vs on bias direct pathway                      
nTask = find(contains(groupOrder,names2plot{2}));
subplot(1,4,2)
plot([.5:2.5],zeros(1,3),'-','Color',[0.6 0.6 0.6]); hold on
plot(1:2, data(xIdx==find(nTask==groups2plot),:),'-','Color', [0.9 0.9 0.9]); hold on
h = errorbar(1,nanmean(dataArrayOFF{nTask}),nanstd(dataArrayOFF{nTask})./sqrt(numel(dataArrayOFF{nTask})-1), ...
    '.-','linewidth',4,'color','k','markersize',0.1); h.CapSize = 16;
h = errorbar(2,nanmean(dataArrayON{nTask}),nanstd(dataArrayON{nTask})./sqrt(numel(dataArrayON{nTask})-1), ...
    '.-','linewidth',4,'color','g','markersize',0.1); h.CapSize = 16;
set(gca,'xtick',1:2,'xticklabel',{'off';'on'})
ylabel('% bias (contra-ipsi)','fontsize',12)
xlim([.5 2.5])
ylim([-80 80])
box off
set(gca,'TickDir','out')
title('direct')
text(0.5,-55, [num2str(numel(dataArrayOFF{nTask})) ' mice' ])
text(0.5,-65, [num2str(nTrialsOFF_xGroup(nTask)) ' off' ])
text(0.5,-75, [num2str(nTrialsON_xGroup(nTask)) ' on' ])

% subplot 3 - off vs on bias no opsin group                     
nTask = find(contains(groupOrder,names2plot{3}));
subplot(1,4,3)
plot([.5:2.5],zeros(1,3),'-','Color',[0.6 0.6 0.6]); hold on
plot(1:2, data(xIdx==find(nTask==groups2plot),:),'-','Color', [0.9 0.9 0.9]); hold on
h = errorbar(1,nanmean(dataArrayOFF{nTask}),nanstd(dataArrayOFF{nTask})./sqrt(numel(dataArrayOFF{nTask})-1), ...
    '.-','linewidth',4,'color','k','markersize',0.1); h.CapSize = 16;
h = errorbar(2,nanmean(dataArrayON{nTask}),nanstd(dataArrayON{nTask})./sqrt(numel(dataArrayON{nTask})-1), ...
    '.-','linewidth',4,'color','g','markersize',0.1); h.CapSize = 16;
set(gca,'xtick',1:2,'xticklabel',{'off';'on'})
ylabel('% bias (contra-ipsi)','fontsize',12)
xlim([.5 2.5])
ylim([-80 80])
box off
set(gca,'TickDir','out')
title('no opsin')
text(0.5,-55, [num2str(numel(dataArrayOFF{nTask})) ' mice' ])
text(0.5,-65, [num2str(nTrialsOFF_xGroup(nTask)) ' off' ])
text(0.5,-75, [num2str(nTrialsON_xGroup(nTask)) ' on' ])

if displayStats == 1
    % unpaired stats for subplot 4 - delta bias across groups
    subplot(1,4,4)     
    groupVec   = categorical(xIdx);
    [p,t,stats] = anova1(dataDelta,groupVec,'off');
    text(0.6,-80,['Group p= ' num2str(p)],'FontSize',6)
    text(0.6,-86,['df= ' num2str(stats.df)],'FontSize',6)
    text(0.6,-92,['F= ' num2str((t{2,5}))],'FontSize',6)
    
    [d2Vno_p, ~ , d2Vno_stats] = ranksum(dataArrayDelta{groups2plot(1)}, ...
                                 dataArrayDelta{groups2plot(3)},'method','approximate');
    [d1Vno_p, ~ , d1Vno_stats] = ranksum(dataArrayDelta{groups2plot(2)}, ...
                                 dataArrayDelta{groups2plot(3)},'method','approximate');
    
    plot(1:3,[90,90,90],'k-')
    text(1,98,['p= ' num2str(d2Vno_p)],'FontSize',6)
    text(1,95,['z= ' num2str(d2Vno_stats.zval)],'FontSize',6)
    text(1,92,['df= ' num2str(numel(dataArrayDelta{groups2plot(1)})+numel(dataArrayDelta{groups2plot(3)})-2)],'FontSize',6)
    
    plot(2:3,[75,75],'k-')
    text(2,83,['p= ' num2str(d1Vno_p)],'FontSize',6)
    text(2,80,['z= ' num2str(d1Vno_stats.zval)],'FontSize',6)
    text(2,77,['df= ' num2str(numel(dataArrayDelta{groups2plot(2)})+numel(dataArrayDelta{groups2plot(3)})-2)],'FontSize',6)

%     % NOT USED -- See ExtendedDataFig8 for related stats
%     % paired stats for subplot 1-3 - off vs on signrank for each group
%     nTask = find(contains(groupOrder,'DMSD2_aoe'));
%     subplot(1,4,nTask)
%     [d2_offVon_p, ~ , d2_offVon_stats] = signrank(dataArrayOFF{groups2plot(nTask)}, ...
%                                  dataArrayON{groups2plot(nTask)},'method','approximate');
%     plot(1:2,[70,70],'k-')
%     text(1,78,['p= ' num2str(d2_offVon_p)],'FontSize',6)
%     text(1,75,['z= ' num2str(d2_offVon_stats.zval)],'FontSize',6)
%     text(1,72,['df= ' num2str(numel(dataArrayOFF{groups2plot(nTask)})-1)],'FontSize',6)
                                                           
else
end


%% panel 2G-H %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
% plot off vs on bias and delta bias for indirect/direct/no opsin DMS
% during no distractors task 
clearvars -except lgDMSD2_aoe lgDMSD1_aoe lgDMSCtrl_aoe ...
                  lgDMSD2_nd lgDMSD1_nd lgDMSCtrl_nd ...
                  lgDMSD2_pc lgDMSD1_pc lgDMSCtrl_pc ...
                  lgNACD2_aoe lgNACD1_aoe lgNACCtrl_aoe ...
                  groupOrder biasOFF_xMouse biasON_xMouse biasDelta_xMouse ...
                  nTrialsOFF_xGroup nTrialsON_xGroup

displayStats = 1; % 0 to suppress display of 
                  % (1) 1-way ANOVA of group on delta bias
                  % (2) post-hoc ranksum unpaired comparison of delta bias
                  %     with no opsin group

% update these 4 variables to detrmine group to plot and color                  
names2plot  = {'DMSD2_nd','DMSD1_nd','DMSCtrl_nd'};                  
groups2plot = [4 5 6];  
barColor    = 'm';
colorIdx    = {[0.99 0.7 0.99],[0.99 0.7 0.99],[0.99 0.7 0.99]}; % [0.8 0.8 0.8],[0.99 0.7 0.99], [0.7 0.99 0.99]

figure
dataArrayOFF = []; dataArrayON = []; dataArrayDelta = [];
data = []; dataDelta = []; xIdx = []; 
dataArrayOFF   = biasOFF_xMouse;
dataArrayON    = biasON_xMouse; 
dataArrayDelta = biasDelta_xMouse;

data(:,1)    = [biasOFF_xMouse{groups2plot(1)}'; ...
                biasOFF_xMouse{groups2plot(2)}'; ...
                biasOFF_xMouse{groups2plot(3)}'];
data(:,2)    = [biasON_xMouse{groups2plot(1)}'; ...
                biasON_xMouse{groups2plot(2)}'; ...
                biasON_xMouse{groups2plot(3)}'];     
dataDelta    = [biasDelta_xMouse{groups2plot(1)}'; ...
                biasDelta_xMouse{groups2plot(2)}'; ...
                biasDelta_xMouse{groups2plot(3)}'];   
xIdx         = [ones(size(biasOFF_xMouse{groups2plot(1)}))';   ...
                ones(size(biasOFF_xMouse{groups2plot(2)}))'*2; ...
                ones(size(biasOFF_xMouse{groups2plot(3)}))'*3];
catMarker   = {'x','x', 'x'};

% subplot 4 - delta bias indirect, direct, no opsin
subplot(1,4,4);
plot([.5:3.5],zeros(1,4),'-','Color',[0.6 0.6 0.6]); hold on
plotSpread(dataDelta,'categoryIdx',xIdx,'distributionIdx',xIdx,...
    'categoryColors',colorIdx,'categoryMarkers',catMarker); 

nTask = find(contains(groupOrder,names2plot{1}));
h = errorbar(1,nanmean(dataArrayDelta{nTask}),nanstd(dataArrayDelta{nTask})./sqrt(numel(dataArrayDelta{nTask})-1), ...
    '.-','linewidth',4,'color',barColor,'markersize',0.1); h.CapSize = 14;

nTask = find(contains(groupOrder,names2plot{2}));
h = errorbar(2,nanmean(dataArrayDelta{nTask}),nanstd(dataArrayDelta{nTask})./sqrt(numel(dataArrayDelta{nTask})-1), ...
    '.-','linewidth',4,'color',barColor,'markersize',0.1); h.CapSize = 14;

nTask = find(contains(groupOrder,names2plot{3}));
h = errorbar(3,nanmean(dataArrayDelta{nTask}),nanstd(dataArrayDelta{nTask})./sqrt(numel(dataArrayDelta{nTask})-1), ...
    '.-','linewidth',4,'color',barColor,'markersize',0.1); h.CapSize = 14;

ylabel('delta contra. bias (%, on-off)','fontsize',12)
xlim([.5 3.5])
ylim([-100 100])
box off
set(gca,'TickDir','out')
set(gca,'xtick',1:3,'xticklabel',{'indirect';'direct';'no opsin'})
rotateXLabels(gca,45)
            
% subplot 1 - off vs on bias indirect pathway                      
nTask = find(contains(groupOrder,names2plot{1}));
subplot(1,4,1)
plot([.5:2.5],zeros(1,3),'-','Color',[0.6 0.6 0.6]); hold on
plot(1:2, data(xIdx==find(nTask==groups2plot),:),'-','Color', [0.9 0.9 0.9]); hold on
h = errorbar(1,nanmean(dataArrayOFF{nTask}),nanstd(dataArrayOFF{nTask})./sqrt(numel(dataArrayOFF{nTask})-1), ...
    '.-','linewidth',4,'color','k','markersize',0.1); h.CapSize = 16;
h = errorbar(2,nanmean(dataArrayON{nTask}),nanstd(dataArrayON{nTask})./sqrt(numel(dataArrayON{nTask})-1), ...
    '.-','linewidth',4,'color','g','markersize',0.1); h.CapSize = 16;
set(gca,'xtick',1:2,'xticklabel',{'off';'on'})
ylabel('% bias (contra-ipsi)','fontsize',12)
xlim([.5 2.5])
ylim([-80 80])
box off
set(gca,'TickDir','out')
title('indirect')
text(0.5,-55, [num2str(numel(dataArrayOFF{nTask})) ' mice' ])
text(0.5,-65, [num2str(nTrialsOFF_xGroup(nTask)) ' off' ])
text(0.5,-75, [num2str(nTrialsON_xGroup(nTask)) ' on' ])

% subplot 2 - off vs on bias direct pathway                      
nTask = find(contains(groupOrder,names2plot{2}));
subplot(1,4,2)
plot([.5:2.5],zeros(1,3),'-','Color',[0.6 0.6 0.6]); hold on
plot(1:2, data(xIdx==find(nTask==groups2plot),:),'-','Color', [0.9 0.9 0.9]); hold on
h = errorbar(1,nanmean(dataArrayOFF{nTask}),nanstd(dataArrayOFF{nTask})./sqrt(numel(dataArrayOFF{nTask})-1), ...
    '.-','linewidth',4,'color','k','markersize',0.1); h.CapSize = 16;
h = errorbar(2,nanmean(dataArrayON{nTask}),nanstd(dataArrayON{nTask})./sqrt(numel(dataArrayON{nTask})-1), ...
    '.-','linewidth',4,'color','g','markersize',0.1); h.CapSize = 16;
set(gca,'xtick',1:2,'xticklabel',{'off';'on'})
ylabel('% bias (contra-ipsi)','fontsize',12)
xlim([.5 2.5])
ylim([-80 80])
box off
set(gca,'TickDir','out')
title('direct')
text(0.5,-55, [num2str(numel(dataArrayOFF{nTask})) ' mice' ])
text(0.5,-65, [num2str(nTrialsOFF_xGroup(nTask)) ' off' ])
text(0.5,-75, [num2str(nTrialsON_xGroup(nTask)) ' on' ])

% subplot 3 - off vs on bias no opsin group                     
nTask = find(contains(groupOrder,names2plot{3}));
subplot(1,4,3)
plot([.5:2.5],zeros(1,3),'-','Color',[0.6 0.6 0.6]); hold on
plot(1:2, data(xIdx==find(nTask==groups2plot),:),'-','Color', [0.9 0.9 0.9]); hold on
h = errorbar(1,nanmean(dataArrayOFF{nTask}),nanstd(dataArrayOFF{nTask})./sqrt(numel(dataArrayOFF{nTask})-1), ...
    '.-','linewidth',4,'color','k','markersize',0.1); h.CapSize = 16;
h = errorbar(2,nanmean(dataArrayON{nTask}),nanstd(dataArrayON{nTask})./sqrt(numel(dataArrayON{nTask})-1), ...
    '.-','linewidth',4,'color','g','markersize',0.1); h.CapSize = 16;
set(gca,'xtick',1:2,'xticklabel',{'off';'on'})
ylabel('% bias (contra-ipsi)','fontsize',12)
xlim([.5 2.5])
ylim([-80 80])
box off
set(gca,'TickDir','out')
title('no opsin')
text(0.5,-55, [num2str(numel(dataArrayOFF{nTask})) ' mice' ])
text(0.5,-65, [num2str(nTrialsOFF_xGroup(nTask)) ' off' ])
text(0.5,-75, [num2str(nTrialsON_xGroup(nTask)) ' on' ])

if displayStats == 1
    % unpaired stats for subplot 4 - delta bias across groups
    subplot(1,4,4)     
    groupVec   = categorical(xIdx);
    [p,t,stats] = anova1(dataDelta,groupVec,'off');
    text(0.6,-80,['Group p= ' num2str(p)],'FontSize',6)
    text(0.6,-86,['df= ' num2str(stats.df)],'FontSize',6)
    text(0.6,-92,['F= ' num2str((t{2,5}))],'FontSize',6)
    
    [d2Vno_p, ~ , d2Vno_stats] = ranksum(dataArrayDelta{groups2plot(1)}, ...
                                 dataArrayDelta{groups2plot(3)},'method','approximate');
    [d1Vno_p, ~ , d1Vno_stats] = ranksum(dataArrayDelta{groups2plot(2)}, ...
                                 dataArrayDelta{groups2plot(3)},'method','approximate');
    
    plot(1:3,[90,90,90],'k-')
    text(1,98,['p= ' num2str(d2Vno_p)],'FontSize',6)
    text(1,95,['z= ' num2str(d2Vno_stats.zval)],'FontSize',6)
    text(1,92,['df= ' num2str(numel(dataArrayDelta{groups2plot(1)})+numel(dataArrayDelta{groups2plot(3)})-2)],'FontSize',6)
    
    plot(2:3,[75,75],'k-')
    text(2,83,['p= ' num2str(d1Vno_p)],'FontSize',6)
    text(2,80,['z= ' num2str(d1Vno_stats.zval)],'FontSize',6)
    text(2,77,['df= ' num2str(numel(dataArrayDelta{groups2plot(2)})+numel(dataArrayDelta{groups2plot(3)})-2)],'FontSize',6)

%     % NOT USED -- See ExtendedDataFig8 for related stats
%     % paired stats for subplot 1-3 - off vs on signrank for each group
%     nTask = find(contains(groupOrder,'DMSD2_aoe'));
%     subplot(1,4,nTask)
%     [d2_offVon_p, ~ , d2_offVon_stats] = signrank(dataArrayOFF{groups2plot(nTask)}, ...
%                                  dataArrayON{groups2plot(nTask)},'method','approximate');
%     plot(1:2,[70,70],'k-')
%     text(1,78,['p= ' num2str(d2_offVon_p)],'FontSize',6)
%     text(1,75,['z= ' num2str(d2_offVon_stats.zval)],'FontSize',6)
%     text(1,72,['df= ' num2str(numel(dataArrayOFF{groups2plot(nTask)})-1)],'FontSize',6)
                                                           
else
end


%% panel 2J-K %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
% plot off vs on bias and delta bias for indirect/direct/no opsin DMS
% during permanent cues task
clearvars -except lgDMSD2_aoe lgDMSD1_aoe lgDMSCtrl_aoe ...
                  lgDMSD2_nd lgDMSD1_nd lgDMSCtrl_nd ...
                  lgDMSD2_pc lgDMSD1_pc lgDMSCtrl_pc ...
                  lgNACD2_aoe lgNACD1_aoe lgNACCtrl_aoe ...
                  groupOrder biasOFF_xMouse biasON_xMouse biasDelta_xMouse ...
                  nTrialsOFF_xGroup nTrialsON_xGroup

displayStats = 1; % 0 to suppress display of 
                  % (1) 1-way ANOVA of group on delta bias
                  % (2) post-hoc ranksum unpaired comparison of delta bias
                  %     with no opsin group
% changing these three fields determines group data to plot and color
names2plot  = {'DMSD2_pc','DMSD1_pc','DMSCtrl_pc'};                  
groups2plot = [7 8 9];  
barColor    = 'c';
colorIdx    = {[0.7 0.99 0.99],[0.7 0.99 0.99],[0.7 0.99 0.99]}; % [0.8 0.8 0.8],[0.99 0.7 0.99], [0.7 0.99 0.99]

figure
dataArrayOFF = []; dataArrayON = []; dataArrayDelta = [];
data = []; dataDelta = []; xIdx = []; 

dataArrayOFF   = biasOFF_xMouse;
dataArrayON    = biasON_xMouse; 
dataArrayDelta = biasDelta_xMouse;

data(:,1)    = [biasOFF_xMouse{groups2plot(1)}'; ...
                biasOFF_xMouse{groups2plot(2)}'; ...
                biasOFF_xMouse{groups2plot(3)}'];
data(:,2)    = [biasON_xMouse{groups2plot(1)}'; ...
                biasON_xMouse{groups2plot(2)}'; ...
                biasON_xMouse{groups2plot(3)}'];     
dataDelta    = [biasDelta_xMouse{groups2plot(1)}'; ...
                biasDelta_xMouse{groups2plot(2)}'; ...
                biasDelta_xMouse{groups2plot(3)}'];   
xIdx         = [ones(size(biasOFF_xMouse{groups2plot(1)}))';   ...
                ones(size(biasOFF_xMouse{groups2plot(2)}))'*2; ...
                ones(size(biasOFF_xMouse{groups2plot(3)}))'*3];
catMarker   = {'x','x', 'x'};

% subplot 4 - delta bias indirect, direct, no opsin
subplot(1,4,4);
plot([.5:3.5],zeros(1,4),'-','Color',[0.6 0.6 0.6]); hold on
plotSpread(dataDelta,'categoryIdx',xIdx,'distributionIdx',xIdx,...
    'categoryColors',colorIdx,'categoryMarkers',catMarker); 

nTask = find(contains(groupOrder,names2plot{1}));
h = errorbar(1,nanmean(dataArrayDelta{nTask}),nanstd(dataArrayDelta{nTask})./sqrt(numel(dataArrayDelta{nTask})-1), ...
    '.-','linewidth',4,'color',barColor,'markersize',0.1); h.CapSize = 14;

nTask = find(contains(groupOrder,names2plot{2}));
h = errorbar(2,nanmean(dataArrayDelta{nTask}),nanstd(dataArrayDelta{nTask})./sqrt(numel(dataArrayDelta{nTask})-1), ...
    '.-','linewidth',4,'color',barColor,'markersize',0.1); h.CapSize = 14;

nTask = find(contains(groupOrder,names2plot{3}));
h = errorbar(3,nanmean(dataArrayDelta{nTask}),nanstd(dataArrayDelta{nTask})./sqrt(numel(dataArrayDelta{nTask})-1), ...
    '.-','linewidth',4,'color',barColor,'markersize',0.1); h.CapSize = 14;

ylabel('delta contra. bias (%, on-off)','fontsize',12)
xlim([.5 3.5])
ylim([-100 100])
box off
set(gca,'TickDir','out')
set(gca,'xtick',1:3,'xticklabel',{'indirect';'direct';'no opsin'})
rotateXLabels(gca,45)
            
% subplot 1 - off vs on bias indirect pathway                      
nTask = find(contains(groupOrder,names2plot{1}));
subplot(1,4,1)
plot([.5:2.5],zeros(1,3),'-','Color',[0.6 0.6 0.6]); hold on
plot(1:2, data(xIdx==find(nTask==groups2plot),:),'-','Color', [0.9 0.9 0.9]); hold on
h = errorbar(1,nanmean(dataArrayOFF{nTask}),nanstd(dataArrayOFF{nTask})./sqrt(numel(dataArrayOFF{nTask})-1), ...
    '.-','linewidth',4,'color','k','markersize',0.1); h.CapSize = 16;
h = errorbar(2,nanmean(dataArrayON{nTask}),nanstd(dataArrayON{nTask})./sqrt(numel(dataArrayON{nTask})-1), ...
    '.-','linewidth',4,'color','g','markersize',0.1); h.CapSize = 16;
set(gca,'xtick',1:2,'xticklabel',{'off';'on'})
ylabel('% bias (contra-ipsi)','fontsize',12)
xlim([.5 2.5])
ylim([-80 80])
box off
set(gca,'TickDir','out')
title('indirect')
text(0.5,-55, [num2str(numel(dataArrayOFF{nTask})) ' mice' ])
text(0.5,-65, [num2str(nTrialsOFF_xGroup(nTask)) ' off' ])
text(0.5,-75, [num2str(nTrialsON_xGroup(nTask)) ' on' ])

% subplot 2 - off vs on bias direct pathway                      
nTask = find(contains(groupOrder,names2plot{2}));
subplot(1,4,2)
plot([.5:2.5],zeros(1,3),'-','Color',[0.6 0.6 0.6]); hold on
plot(1:2, data(xIdx==find(nTask==groups2plot),:),'-','Color', [0.9 0.9 0.9]); hold on
h = errorbar(1,nanmean(dataArrayOFF{nTask}),nanstd(dataArrayOFF{nTask})./sqrt(numel(dataArrayOFF{nTask})-1), ...
    '.-','linewidth',4,'color','k','markersize',0.1); h.CapSize = 16;
h = errorbar(2,nanmean(dataArrayON{nTask}),nanstd(dataArrayON{nTask})./sqrt(numel(dataArrayON{nTask})-1), ...
    '.-','linewidth',4,'color','g','markersize',0.1); h.CapSize = 16;
set(gca,'xtick',1:2,'xticklabel',{'off';'on'})
ylabel('% bias (contra-ipsi)','fontsize',12)
xlim([.5 2.5])
ylim([-80 80])
box off
set(gca,'TickDir','out')
title('direct')
text(0.5,-55, [num2str(numel(dataArrayOFF{nTask})) ' mice' ])
text(0.5,-65, [num2str(nTrialsOFF_xGroup(nTask)) ' off' ])
text(0.5,-75, [num2str(nTrialsON_xGroup(nTask)) ' on' ])

% subplot 3 - off vs on bias no opsin group                     
nTask = find(contains(groupOrder,names2plot{3}));
subplot(1,4,3)
plot([.5:2.5],zeros(1,3),'-','Color',[0.6 0.6 0.6]); hold on
plot(1:2, data(xIdx==find(nTask==groups2plot),:),'-','Color', [0.9 0.9 0.9]); hold on
h = errorbar(1,nanmean(dataArrayOFF{nTask}),nanstd(dataArrayOFF{nTask})./sqrt(numel(dataArrayOFF{nTask})-1), ...
    '.-','linewidth',4,'color','k','markersize',0.1); h.CapSize = 16;
h = errorbar(2,nanmean(dataArrayON{nTask}),nanstd(dataArrayON{nTask})./sqrt(numel(dataArrayON{nTask})-1), ...
    '.-','linewidth',4,'color','g','markersize',0.1); h.CapSize = 16;
set(gca,'xtick',1:2,'xticklabel',{'off';'on'})
ylabel('% bias (contra-ipsi)','fontsize',12)
xlim([.5 2.5])
ylim([-80 80])
box off
set(gca,'TickDir','out')
title('no opsin')
text(0.5,-55, [num2str(numel(dataArrayOFF{nTask})) ' mice' ])
text(0.5,-65, [num2str(nTrialsOFF_xGroup(nTask)) ' off' ])
text(0.5,-75, [num2str(nTrialsON_xGroup(nTask)) ' on' ])

if displayStats == 1
    % unpaired stats for subplot 4 - delta bias across groups
    subplot(1,4,4)     
    groupVec   = categorical(xIdx);
    [p,t,stats] = anova1(dataDelta,groupVec,'off');
    text(0.6,-80,['Group p= ' num2str(p)],'FontSize',6)
    text(0.6,-86,['df= ' num2str(stats.df)],'FontSize',6)
    text(0.6,-92,['F= ' num2str((t{2,5}))],'FontSize',6)
    
    [d2Vno_p, ~ , d2Vno_stats] = ranksum(dataArrayDelta{groups2plot(1)}, ...
                                 dataArrayDelta{groups2plot(3)},'method','approximate');
    [d1Vno_p, ~ , d1Vno_stats] = ranksum(dataArrayDelta{groups2plot(2)}, ...
                                 dataArrayDelta{groups2plot(3)},'method','approximate');
    
    plot(1:3,[90,90,90],'k-')
    text(1,98,['p= ' num2str(d2Vno_p)],'FontSize',6)
    text(1,95,['z= ' num2str(d2Vno_stats.zval)],'FontSize',6)
    text(1,92,['df= ' num2str(numel(dataArrayDelta{groups2plot(1)})+numel(dataArrayDelta{groups2plot(3)})-2)],'FontSize',6)
    
    plot(2:3,[75,75],'k-')
    text(2,83,['p= ' num2str(d1Vno_p)],'FontSize',6)
    text(2,80,['z= ' num2str(d1Vno_stats.zval)],'FontSize',6)
    text(2,77,['df= ' num2str(numel(dataArrayDelta{groups2plot(2)})+numel(dataArrayDelta{groups2plot(3)})-2)],'FontSize',6)

%     % NOT USED -- See ExtendedDataFig8 for related stats
%     % paired stats for subplot 1-3 - off vs on signrank for each group
%     nTask = find(contains(groupOrder,'DMSD2_aoe'));
%     subplot(1,4,nTask)
%     [d2_offVon_p, ~ , d2_offVon_stats] = signrank(dataArrayOFF{groups2plot(nTask)}, ...
%                                  dataArrayON{groups2plot(nTask)},'method','approximate');
%     plot(1:2,[70,70],'k-')
%     text(1,78,['p= ' num2str(d2_offVon_p)],'FontSize',6)
%     text(1,75,['z= ' num2str(d2_offVon_stats.zval)],'FontSize',6)
%     text(1,72,['df= ' num2str(numel(dataArrayOFF{groups2plot(nTask)})-1)],'FontSize',6)
                                                           
else
end

%% panel 2O-P %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
% plot off vs on bias and delta bias for indirect/direct/no opsin NAC
% during evidence accumulation 
displayStats = 1; % 0 to suppress display of 
                  % (1) 1-way ANOVA of group on delta bias
                  % (2) post-hoc ranksum unpaired comparison of delta bias
                  %     with no opsin group                  
      
% update these 4 variables to determine group to plot and color                  
names2plot  = {'NACD2_aoe','NACD1_aoe','NACCtrl_aoe'};                  
groups2plot = [10 11 12];  
barColor    = 'k';     
colorIdx    = {[0.8 0.8 0.8],[0.8 0.8 0.8],[0.8 0.8 0.8]}; % [0.8 0.8 0.8],[0.99 0.7 0.99], [0.7 0.99 0.99]

figure
subplot(1,4,1)
dataArrayOFF = []; dataArrayON = []; dataArrayDelta = [];
data = []; dataDelta = []; xIdx = [];

dataArrayOFF   = biasOFF_xMouse;
dataArrayON    = biasON_xMouse; 
dataArrayDelta = biasDelta_xMouse;

data(:,1)    = [biasOFF_xMouse{groups2plot(1)}'; ...
                biasOFF_xMouse{groups2plot(2)}'; ...
                biasOFF_xMouse{groups2plot(3)}'];
data(:,2)    = [biasON_xMouse{groups2plot(1)}'; ...
                biasON_xMouse{groups2plot(2)}'; ...
                biasON_xMouse{groups2plot(3)}'];     
dataDelta    = [biasDelta_xMouse{groups2plot(1)}'; ...
                biasDelta_xMouse{groups2plot(2)}'; ...
                biasDelta_xMouse{groups2plot(3)}'];   
xIdx         = [ones(size(biasOFF_xMouse{groups2plot(1)}))';   ...
                ones(size(biasOFF_xMouse{groups2plot(2)}))'*2; ...
                ones(size(biasOFF_xMouse{groups2plot(3)}))'*3];
catMarker   = {'x','x', 'x'};

% subplot 4 - delta bias indirect, direct, no opsin
subplot(1,4,4);
plot([.5:3.5],zeros(1,4),'-','Color',[0.6 0.6 0.6]); hold on
plotSpread(dataDelta,'categoryIdx',xIdx,'distributionIdx',xIdx,...
    'categoryColors',colorIdx,'categoryMarkers',catMarker); 

nTask = find(contains(groupOrder,names2plot{1}));
h = errorbar(1,nanmean(dataArrayDelta{nTask}),nanstd(dataArrayDelta{nTask})./sqrt(numel(dataArrayDelta{nTask})-1), ...
    '.-','linewidth',4,'color',barColor,'markersize',0.1); h.CapSize = 14;

nTask = find(contains(groupOrder,names2plot{2}));
h = errorbar(2,nanmean(dataArrayDelta{nTask}),nanstd(dataArrayDelta{nTask})./sqrt(numel(dataArrayDelta{nTask})-1), ...
    '.-','linewidth',4,'color',barColor,'markersize',0.1); h.CapSize = 14;

nTask = find(contains(groupOrder,names2plot{3}));
h = errorbar(3,nanmean(dataArrayDelta{nTask}),nanstd(dataArrayDelta{nTask})./sqrt(numel(dataArrayDelta{nTask})-1), ...
    '.-','linewidth',4,'color',barColor,'markersize',0.1); h.CapSize = 14;

ylabel('delta contra. bias (%, on-off)','fontsize',12)
xlim([.5 3.5])
ylim([-100 100])
box off
set(gca,'TickDir','out')
set(gca,'xtick',1:3,'xticklabel',{'indirect';'direct';'no opsin'})
rotateXLabels(gca,45)
            
% subplot 1 - off vs on bias indirect pathway                      
nTask = find(contains(groupOrder,names2plot{1}));
subplot(1,4,1)
plot([.5:2.5],zeros(1,3),'-','Color',[0.6 0.6 0.6]); hold on
plot(1:2, data(xIdx==find(nTask==groups2plot),:),'-','Color', [0.9 0.9 0.9]); hold on
h = errorbar(1,nanmean(dataArrayOFF{nTask}),nanstd(dataArrayOFF{nTask})./sqrt(numel(dataArrayOFF{nTask})-1), ...
    '.-','linewidth',4,'color','k','markersize',0.1); h.CapSize = 16;
h = errorbar(2,nanmean(dataArrayON{nTask}),nanstd(dataArrayON{nTask})./sqrt(numel(dataArrayON{nTask})-1), ...
    '.-','linewidth',4,'color','g','markersize',0.1); h.CapSize = 16;
set(gca,'xtick',1:2,'xticklabel',{'off';'on'})
ylabel('% bias (contra-ipsi)','fontsize',12)
xlim([.5 2.5])
ylim([-80 80])
box off
set(gca,'TickDir','out')
title('indirect')
text(0.5,-55, [num2str(numel(dataArrayOFF{nTask})) ' mice' ])
text(0.5,-65, [num2str(nTrialsOFF_xGroup(nTask)) ' off' ])
text(0.5,-75, [num2str(nTrialsON_xGroup(nTask)) ' on' ])

% subplot 2 - off vs on bias direct pathway                      
nTask = find(contains(groupOrder,names2plot{2}));
subplot(1,4,2)
plot([.5:2.5],zeros(1,3),'-','Color',[0.6 0.6 0.6]); hold on
plot(1:2, data(xIdx==find(nTask==groups2plot),:),'-','Color', [0.9 0.9 0.9]); hold on
h = errorbar(1,nanmean(dataArrayOFF{nTask}),nanstd(dataArrayOFF{nTask})./sqrt(numel(dataArrayOFF{nTask})-1), ...
    '.-','linewidth',4,'color','k','markersize',0.1); h.CapSize = 16;
h = errorbar(2,nanmean(dataArrayON{nTask}),nanstd(dataArrayON{nTask})./sqrt(numel(dataArrayON{nTask})-1), ...
    '.-','linewidth',4,'color','g','markersize',0.1); h.CapSize = 16;
set(gca,'xtick',1:2,'xticklabel',{'off';'on'})
ylabel('% bias (contra-ipsi)','fontsize',12)
xlim([.5 2.5])
ylim([-80 80])
box off
set(gca,'TickDir','out')
title('direct')
text(0.5,-55, [num2str(numel(dataArrayOFF{nTask})) ' mice' ])
text(0.5,-65, [num2str(nTrialsOFF_xGroup(nTask)) ' off' ])
text(0.5,-75, [num2str(nTrialsON_xGroup(nTask)) ' on' ])

% subplot 3 - off vs on bias no opsin group                     
nTask = find(contains(groupOrder,names2plot{3}));
subplot(1,4,3)
plot([.5:2.5],zeros(1,3),'-','Color',[0.6 0.6 0.6]); hold on
plot(1:2, data(xIdx==find(nTask==groups2plot),:),'-','Color', [0.9 0.9 0.9]); hold on
h = errorbar(1,nanmean(dataArrayOFF{nTask}),nanstd(dataArrayOFF{nTask})./sqrt(numel(dataArrayOFF{nTask})-1), ...
    '.-','linewidth',4,'color','k','markersize',0.1); h.CapSize = 16;
h = errorbar(2,nanmean(dataArrayON{nTask}),nanstd(dataArrayON{nTask})./sqrt(numel(dataArrayON{nTask})-1), ...
    '.-','linewidth',4,'color','g','markersize',0.1); h.CapSize = 16;
set(gca,'xtick',1:2,'xticklabel',{'off';'on'})
ylabel('% bias (contra-ipsi)','fontsize',12)
xlim([.5 2.5])
ylim([-80 80])
box off
set(gca,'TickDir','out')
title('no opsin')
text(0.5,-55, [num2str(numel(dataArrayOFF{nTask})) ' mice' ])
text(0.5,-65, [num2str(nTrialsOFF_xGroup(nTask)) ' off' ])
text(0.5,-75, [num2str(nTrialsON_xGroup(nTask)) ' on' ])

if displayStats == 1
    % unpaired stats for subplot 4 - delta bias across groups
    subplot(1,4,4)     
    groupVec   = categorical(xIdx);
    [p,t,stats] = anova1(dataDelta,groupVec,'off');
    text(0.6,-80,['Group p= ' num2str(p)],'FontSize',6)
    text(0.6,-86,['df= ' num2str(stats.df)],'FontSize',6)
    text(0.6,-92,['F= ' num2str((t{2,5}))],'FontSize',6)
    
    [d2Vno_p, ~ , d2Vno_stats] = ranksum(dataArrayDelta{groups2plot(1)}, ...
                                 dataArrayDelta{groups2plot(3)},'method','approximate');
    [d1Vno_p, ~ , d1Vno_stats] = ranksum(dataArrayDelta{groups2plot(2)}, ...
                                 dataArrayDelta{groups2plot(3)},'method','approximate');
    
    plot(1:3,[90,90,90],'k-')
    text(1,98,['p= ' num2str(d2Vno_p)],'FontSize',6)
    text(1,95,['z= ' num2str(d2Vno_stats.zval)],'FontSize',6)
    text(1,92,['df= ' num2str(numel(dataArrayDelta{groups2plot(1)})+numel(dataArrayDelta{groups2plot(3)})-2)],'FontSize',6)
    
    plot(2:3,[75,75],'k-')
    text(2,83,['p= ' num2str(d1Vno_p)],'FontSize',6)
    text(2,80,['z= ' num2str(d1Vno_stats.zval)],'FontSize',6)
    text(2,77,['df= ' num2str(numel(dataArrayDelta{groups2plot(2)})+numel(dataArrayDelta{groups2plot(3)})-2)],'FontSize',6)

%     % NOT USED -- See ExtendedDataFig8 for related stats
%     % paired stats for subplot 1-3 - off vs on signrank for each group
%     nTask = find(contains(groupOrder,'DMSD2_aoe'));
%     subplot(1,4,nTask)
%     [d2_offVon_p, ~ , d2_offVon_stats] = signrank(dataArrayOFF{groups2plot(nTask)}, ...
%                                  dataArrayON{groups2plot(nTask)},'method','approximate');
%     plot(1:2,[70,70],'k-')
%     text(1,78,['p= ' num2str(d2_offVon_p)],'FontSize',6)
%     text(1,75,['z= ' num2str(d2_offVon_stats.zval)],'FontSize',6)
%     text(1,72,['df= ' num2str(numel(dataArrayOFF{groups2plot(nTask)})-1)],'FontSize',6)
                                                           
else
end

