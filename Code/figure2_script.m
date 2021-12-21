%% Script to generate Fig 2 data plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% cross task comparison of performance and motor indicators 2C-G
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

groupOrder = {'_aoe', '_nd', '_pc'};
counter = 1;
for groupNum = 1:numel(groupOrder)
    clear lg lgClean 
    lg          = eval(['concatLog_OFF' groupOrder{groupNum}]); 
    lg          = selectLogSubset(lg,[],[],[],[],1,0.6); % drops blocks with perf <0.6
    [lgClean,~] = selectLogSubset(lg,[],[],[],[],0,0.6); % drops blocks with perf <0.6 and excessTravel>0.1
    
    % overall accuracy here panel 2C
    ipsiCon = 0; cleanLog = 1; relaxTravel = 0;
    tempAcc = []; tempSumm = [];
    [tempAcc, tempSumm] = xMousePerfBias(lg, ipsiCon, cleanLog, relaxTravel);
    accuracy_xMouse{counter}(1,:)   = [tempAcc(:).perfOFF]*100;     
    nTrialsSubselect_xTask(counter) = tempSumm.nTrialsPostSelection;
    nTrialsTotal_xTask(counter)     = tempSumm.nTrialsPreSelection;

    % panel 2D - calculation of Yvelocity in bins of 25cm for 300cm stem
    binIDs = 1:25:301;    
    for nBin = 1:numel(binIDs)-1
        yVel_xMouse{counter}(nBin,:)   = xMouseSpeedXY(lgClean,[binIDs(nBin) binIDs(nBin+1)], 2);
    end
    
    % panel 2E - calculation of x-position in bins of 25cm for 300cm stem 
    % for left and right choice trials seperately
    avgXpos = [];
    pos2    = []; % next line inverts matrix columns in each cell array of pos
                  % hack workaround as next function samples 3rd column based
                  % on 2nd column
                  % pos2 is view angle, y-position, x-position
                  % pos is x-position, y-position, view angle
    pos2    = cellfun(@fliplr, lgClean.pos,'UniformOutput',false); %       
    avgXpos = sampleViewAngleVsY_average(pos2,[1 301],25);
    
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
    
    % panel 2F - calculation of view angle in bins of 25cm for 300cm stem 
    % for left and right choice trials seperately
    avgViewAng  = [];
    avgViewAng  = sampleViewAngleVsY_average(lgClean.pos,[1 301],25);
    
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
         
    % panel 2G - calculation of average total distance per trial - 
    % no trial subselection
    distance = [];
    distance = 330+(lg.excessTravel*330);
    
    mouseNums = [];
    mouseNums = unique(lg.mouseID);
    for nMouse = 1:numel(mouseNums)
        distance_xMouse{counter}(:,nMouse) = nanmean(distance(lg.mouseID == mouseNums(nMouse)));
    end
        
    % for the next loop of analysis given logs in groupOrder       
    counter = counter +1;
end

%% clear all except the key processed variables to plot (source data)
clearvars -except groupOrder accuracy_xMouse nTrialsSubselect_xTask nTrialsTotal_xTask ...
           yVel_xMouse xPosL_xMouse xPosR_xMouse viewAngL_xMouse viewAngR_xMouse distance_xMouse ...
           concatLog_OFF_aoe concatLog_OFF_nd concatLog_OFF_pc binIDs
       
%% panel 2C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
% plot accuracy by task
displayStats = 1; % 0 to suppress display of 1-way ANOVA of task on accuracy
                  % and ranksum unpaired comparisons between each task
                  
figure
dataArray   = []; data = []; xIdx = []; colorIdx = [];
dataArray   = accuracy_xMouse;
data        = [accuracy_xMouse{1}'; ...
               accuracy_xMouse{2}'; ...
               accuracy_xMouse{3}'];
xIdx        = [ones(size(accuracy_xMouse{1}))';   ...
               ones(size(accuracy_xMouse{2}))'*2; ...
               ones(size(accuracy_xMouse{3}))'*3];
catMarker   = {'x','x', 'x'};
colorIdx    = {[0.8 0.8 0.8],[0.99 0.7 0.99], [0.7 0.99 0.99]};
hold on
plotSpread(data,'categoryIdx',xIdx,'distributionIdx',xIdx,...
    'categoryColors',colorIdx,'categoryMarkers',catMarker);

nTask = 1;
h = errorbar(nTask,nanmean(dataArray{nTask}),nanstd(dataArray{nTask})./sqrt(numel(dataArray{nTask})-1), ...
    '.-','linewidth',4,'color','k','markersize',0.1); h.CapSize = 16;
text(0.8, 64, [num2str(numel(dataArray{nTask})) ' mice' ])
text(0.8, 62, [num2str(nTrialsSubselect_xTask(nTask)) ' trials' ])

nTask = 2; 
h = errorbar(nTask,nanmean(dataArray{nTask}),nanstd(dataArray{nTask})./sqrt(numel(dataArray{nTask})-1), ...
    '.-','linewidth',4,'color','m','markersize',0.1); h.CapSize = 16;
text(1.8, 64, [num2str(numel(dataArray{nTask})) ' mice' ])
text(1.8, 62, [num2str(nTrialsSubselect_xTask(nTask)) ' trials' ])

nTask = 3;
h = errorbar(nTask,nanmean(dataArray{nTask}),nanstd(dataArray{nTask})./sqrt(numel(dataArray{nTask})-1), ...
    '.-','linewidth',4,'color','c','markersize',0.1); h.CapSize = 16;
text(2.8, 64, [num2str(numel(dataArray{nTask})) ' mice' ])
text(2.8, 62, [num2str(nTrialsSubselect_xTask(nTask)) ' trials' ])

set(gca,'xtick',1:3,'xticklabel',{'AoE';'ctrl#1';'ctrl#2'})
rotateXLabels(gca,45)
ylabel('accuracy (%)','fontsize',12)
xlim([.65 3.35])
ylim([60 100])
box off
set(gca,'TickDir','out');

if displayStats == 1
    xIdx = categorical(xIdx);
    [p,t,stats] = anova1(data,xIdx,'off');
    text(1.5,74,['ANOVA p= ' num2str(p)])
    text(1.5,72,['df= ' num2str(stats.df)])
    text(1.5,70,['F= ' num2str((t{2,5}))])
    
    [aoeVpc_p, ~ , aoeVpc_stats] = ranksum(dataArray{1}, dataArray{3},'method','approximate');
    [aoeVnd_p, ~ , aoeVnd_stats] = ranksum(dataArray{1}, dataArray{2},'method','approximate');
    [ndVpc_p , ~ , ndVpc_stats]  = ranksum(dataArray{2} , dataArray{3},'method','approximate');
    
    plot(1:3,[98,98,98],'k-')
    text(1,99,['p= ' num2str(aoeVpc_p)])
    text(1.8,99,['z= ' num2str(aoeVpc_stats.zval)])
    text(2.8,99,['df= ' num2str(numel(dataArray{1})+numel(dataArray{3})-1)])
    
    plot(1:2,[96,96],'k-')
    text(1,97,['p= ' num2str(aoeVnd_p)])
    text(1,95,['z= ' num2str(aoeVnd_stats.zval)])
    text(1,93,['df= ' num2str(numel(dataArray{1})+numel(dataArray{2})-1)])
    
    plot(2:3,[94,94],'k-')
    text(2.8,97,['p= ' num2str(ndVpc_p)])
    text(2.8,95,['z= ' num2str(ndVpc_stats.zval)])
    text(2.8,93,['df= ' num2str(numel(dataArray{2})+numel(dataArray{3})-1)])
else
end

%% panel 2E %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot y-velocity across maze stem
clearvars -except groupOrder accuracy_xMouse nTrialsSubselect_xTask nTrialsTotal_xTask ...
           yVel_xMouse xPosL_xMouse xPosR_xMouse viewAngL_xMouse viewAngR_xMouse distance_xMouse ...
           concatLog_OFF_aoe concatLog_OFF_nd concatLog_OFF_pc binIDs
displayStats = 1; % 0 to suppress display of 1-way ANOVA of task on accuracy
                  % and ranksum unpaired comparisons between each task
ANOVA_sig    = 0.05; 
postHoc_sig  = 0.05;                   
       
bins2plot   = binIDs(2:end)-13;
plotOffsets = [0 5 -5];
dataArray   = []; data = []; xIdx = []; colorIdx = [];
dataArray   = yVel_xMouse;
data        = [yVel_xMouse{1}'; ...
               yVel_xMouse{2}'; ...
               yVel_xMouse{3}'];
xIdx        = [ones(size(yVel_xMouse{1},2),1);   ...
               ones(size(yVel_xMouse{2},2),1)*2; ...
               ones(size(yVel_xMouse{3},2),1)*3];           
                  
figure
nTask = 1;
meanPlot = []; meanPlot = nanmean(dataArray{nTask},2);
errPlot  = []; errPlot  = (nanstd(dataArray{nTask},0,2)) / sqrt(size(dataArray{nTask},2)-1);
errorbar(meanPlot, bins2plot'+ plotOffsets(nTask) ,errPlot,'horizontal','k-'); hold on
text(25,80,['AoE n= ' num2str(size(dataArray{nTask},2)) ],'Color','k','FontSize',10) 

nTask = 2;
meanPlot = []; meanPlot = nanmean(dataArray{nTask},2);
errPlot  = []; errPlot  = (nanstd(dataArray{nTask},0,2)) / sqrt(size(dataArray{nTask},2)-1);
errorbar(meanPlot, bins2plot'+ plotOffsets(nTask) ,errPlot,'horizontal','m-');
text(25,70,['ctrl#1 n= ' num2str(size(dataArray{nTask},2)) ],'Color','m','FontSize',10) 

nTask = 3;
meanPlot = []; meanPlot = nanmean(dataArray{nTask},2);
errPlot  = []; errPlot  = (nanstd(dataArray{nTask},0,2)) / sqrt(size(dataArray{nTask},2)-1);
errorbar(meanPlot, bins2plot'+ plotOffsets(nTask) ,errPlot,'horizontal','c-');
text(25,60,['ctrl#2 n= ' num2str(size(dataArray{nTask},2)) ],'Color','c','FontSize',10) 

ylim([-0.5 300.5])
xlim([20 80])
xlabel('y-velocity (cm/s)')
set(gca,'TickDir','out');
ylabel('y-position (cm)')
box off

if displayStats == 1    
    xIdx = categorical(xIdx);      
    for nBin = 1:numel(bins2plot)
        p = []; t = []; stats = [];
        [p,t,stats]      = anova1(data(:,nBin),xIdx,'off');
        pTask(nBin)      = p; 
        Fstat(nBin)      = t{2,5}; 
        degFreedom(nBin) = stats.df;
        
        if pTask(nBin) < ANOVA_sig
            plot(22,bins2plot(nBin),'k*')
        else
        end
              
        aoeVpc_stats   = [];
        [aoeVpc_p(nBin), ~ , aoeVpc_stats] = ranksum(dataArray{1}(nBin,:), ...
                                       dataArray{3}(nBin,:),'method','approximate');
        aoeVpc_z(nBin) = aoeVpc_stats.zval;  
        
        if aoeVpc_p(nBin) < postHoc_sig
            plot(78,bins2plot(nBin),'c*')
        else
        end     
                                   
        aoeVnd_stats   = [];                           
        [aoeVnd_p(nBin), ~ , aoeVnd_stats] = ranksum(dataArray{1}(nBin,:), ...
                                       dataArray{2}(nBin,:),'method','approximate');
        aoeVnd_z(nBin) = aoeVnd_stats.zval; 
        
        if aoeVnd_p(nBin) < postHoc_sig
            plot(79,bins2plot(nBin),'m*')
        else
        end
                    
        ndVpc_stats    = [];
        [ndVpc_p(nBin), ~ , ndVpc_stats]  = ranksum(dataArray{2}(nBin,:) , ...
                                       dataArray{3}(nBin,:),'method','approximate');
        ndVpc_z(nBin)  = ndVpc_stats.zval;
        
        if ndVpc_p(nBin) < postHoc_sig
            plot(80,bins2plot(nBin),'k*')
        else
        end
        
    end
else
end
       
%% panel 2E %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot x-position trajectory on left and right trials by task    
clearvars -except groupOrder accuracy_xMouse nTrialsSubselect_xTask nTrialsTotal_xTask ...
           yVel_xMouse xPosL_xMouse xPosR_xMouse viewAngL_xMouse viewAngR_xMouse distance_xMouse ...
           concatLog_OFF_aoe concatLog_OFF_nd concatLog_OFF_pc binIDs  
displayStats = 1; % 0 to suppress display of 1-way ANOVA of task on accuracy
                  % and ranksum unpaired comparisons between each task
ANOVA_sig    = 0.05; 
postHoc_sig  = 0.05;
    
       
bins2plot   = binIDs(2:end)-13;
plotOffsets = [0 5 -5];
dataArray   = []; data = []; xIdx = []; colorIdx = [];
dataArray   = xPosL_xMouse;
data        = [xPosL_xMouse{1}'; ...
               xPosL_xMouse{2}'; ...
               xPosL_xMouse{3}'];
xIdx        = [ones(size(xPosL_xMouse{1},2),1);   ...
               ones(size(xPosL_xMouse{2},2),1)*2; ...
               ones(size(xPosL_xMouse{3},2),1)*3];           
                  
figure
subplot(1,2,1)
nTask = 1;
meanPlot = []; meanPlot = nanmean(dataArray{nTask},2);
errPlot  = []; errPlot  = (nanstd(dataArray{nTask},0,2)) / sqrt(size(dataArray{nTask},2)-1);
errorbar(meanPlot, bins2plot'+ plotOffsets(nTask) ,errPlot,'horizontal','k-'); hold on
text(-1.4,80,['AoE n= ' num2str(size(dataArray{nTask},2)) ],'Color','k','FontSize',10) 

nTask = 2;
meanPlot = []; meanPlot = nanmean(dataArray{nTask},2);
errPlot  = []; errPlot  = (nanstd(dataArray{nTask},0,2)) / sqrt(size(dataArray{nTask},2)-1);
errorbar(meanPlot, bins2plot'+ plotOffsets(nTask) ,errPlot,'horizontal','m-');
text(-1.4,70,['ctrl#1 n= ' num2str(size(dataArray{nTask},2)) ],'Color','m','FontSize',10) 

nTask = 3;
meanPlot = []; meanPlot = nanmean(dataArray{nTask},2);
errPlot  = []; errPlot  = (nanstd(dataArray{nTask},0,2)) / sqrt(size(dataArray{nTask},2)-1);
errorbar(meanPlot, bins2plot'+ plotOffsets(nTask) ,errPlot,'horizontal','c-');
text(-1.4,60,['ctrl#2 n= ' num2str(size(dataArray{nTask},2)) ],'Color','c','FontSize',10) 

plot(zeros(300,1),[1:300],'--','Color',[0.8 0.8 0.8])
title('left choice')
ylim([-0.5 300.5])
xlim([-1.5 0.25])
xlabel('x-position (cm)')
set(gca,'TickDir','out');
ylabel('y-position (cm)')
box off

if displayStats == 1    
    xIdx = categorical(xIdx);  
    pTask    = [];  Fstat    = [];  degFreedom = [];
    aoeVnd_p = [];  aoeVnd_z = []; 
    aoeVpc_p = [];  aoeVpc_z = []; 
    ndVpc_p  = [];  ndVpc_z  = [];  
    for nBin = 1:numel(bins2plot)
        p = []; t = []; stats = [];
        [p,t,stats]      = anova1(data(:,nBin),xIdx,'off');
        pTask(nBin)      = p; 
        Fstat(nBin)      = t{2,5}; 
        degFreedom(nBin) = stats.df;
        
        if pTask(nBin) < ANOVA_sig 
            plot(0.2,bins2plot(nBin),'k*')
        else
        end
              
        aoeVpc_stats   = [];
        [aoeVpc_p(nBin), ~ , aoeVpc_stats] = ranksum(dataArray{1}(nBin,:), ...
                                       dataArray{3}(nBin,:),'method','approximate');
        aoeVpc_z(nBin) = aoeVpc_stats.zval;  
        
        if aoeVpc_p(nBin) < postHoc_sig
            plot(-1.5,bins2plot(nBin),'c*')
        else
        end     
                                   
        aoeVnd_stats   = [];                           
        [aoeVnd_p(nBin), ~ , aoeVnd_stats] = ranksum(dataArray{1}(nBin,:), ...
                                       dataArray{2}(nBin,:),'method','approximate');
        aoeVnd_z(nBin) = aoeVnd_stats.zval; 
        
        if aoeVnd_p(nBin) < postHoc_sig
            plot(-1.45,bins2plot(nBin),'m*')
        else
        end
                    
        ndVpc_stats    = [];
        [ndVpc_p(nBin), ~ , ndVpc_stats]  = ranksum(dataArray{2}(nBin,:) , ...
                                       dataArray{3}(nBin,:),'method','approximate');
        ndVpc_z(nBin)  = ndVpc_stats.zval;
        
        if ndVpc_p(nBin) < postHoc_sig
            plot(-1.4,bins2plot(nBin),'k*')
        else
        end
        
    end   
    pTask_L      = pTask;
    Fstat_L      = Fstat;
    degFreedom_L = degFreedom;
    aoeVnd_p_L   = aoeVnd_p;
    aoeVnd_z_L   = aoeVnd_z;
    aoeVpc_p_L   = aoeVpc_p;
    aoeVpc_z_L   = aoeVpc_z;
    ndVpc_p_L    = ndVpc_p;
    ndVpc_z_L    = ndVpc_z;   
else
end


% right choice trials
bins2plot   = binIDs(2:end)-13;
plotOffsets = [0 5 -5];
dataArray   = []; data = []; xIdx = []; colorIdx = [];
dataArray   = xPosR_xMouse;
data        = [xPosR_xMouse{1}'; ...
               xPosR_xMouse{2}'; ...
               xPosR_xMouse{3}'];
xIdx        = [ones(size(xPosR_xMouse{1},2),1);   ...
               ones(size(xPosR_xMouse{2},2),1)*2; ...
               ones(size(xPosR_xMouse{3},2),1)*3];           
                  
subplot(1,2,2)
nTask = 1;
meanPlot = []; meanPlot = nanmean(dataArray{nTask},2);
errPlot  = []; errPlot  = (nanstd(dataArray{nTask},0,2)) / sqrt(size(dataArray{nTask},2)-1);
errorbar(meanPlot, bins2plot'+ plotOffsets(nTask) ,errPlot,'horizontal','k-'); hold on

nTask = 2;
meanPlot = []; meanPlot = nanmean(dataArray{nTask},2);
errPlot  = []; errPlot  = (nanstd(dataArray{nTask},0,2)) / sqrt(size(dataArray{nTask},2)-1);
errorbar(meanPlot, bins2plot'+ plotOffsets(nTask) ,errPlot,'horizontal','m-');

nTask = 3;
meanPlot = []; meanPlot = nanmean(dataArray{nTask},2);
errPlot  = []; errPlot  = (nanstd(dataArray{nTask},0,2)) / sqrt(size(dataArray{nTask},2)-1);
errorbar(meanPlot, bins2plot'+ plotOffsets(nTask) ,errPlot,'horizontal','c-');

plot(zeros(300,1),[1:300],'--','Color',[0.8 0.8 0.8])
title('right choice')
ylim([-0.5 300.5])
xlim([-0.25 1.5])
xlabel('x-position (cm)')
box off
set(gca,'YAxisLocation','right');

if displayStats == 1    
    xIdx = categorical(xIdx);  
    pTask    = [];  Fstat    = [];  degFreedom = [];
    aoeVnd_p = [];  aoeVnd_z = []; 
    aoeVpc_p = [];  aoeVpc_z = []; 
    ndVpc_p  = [];  ndVpc_z  = [];    
    
    for nBin = 1:numel(bins2plot)
        p = []; t = []; stats = [];
        [p,t,stats]      = anova1(data(:,nBin),xIdx,'off');
        pTask(nBin)      = p; 
        Fstat(nBin)      = t{2,5}; 
        degFreedom(nBin) = stats.df;
        
        if pTask(nBin) < ANOVA_sig 
            plot(-0.2,bins2plot(nBin),'k*')
        else
        end
              
        aoeVpc_stats   = [];
        [aoeVpc_p(nBin), ~ , aoeVpc_stats] = ranksum(dataArray{1}(nBin,:), ...
                                       dataArray{3}(nBin,:),'method','approximate');
        aoeVpc_z(nBin) = aoeVpc_stats.zval;  
        
        if aoeVpc_p(nBin) < postHoc_sig
            plot(1.5,bins2plot(nBin),'c*')
        else
        end     
                                   
        aoeVnd_stats   = [];                           
        [aoeVnd_p(nBin), ~ , aoeVnd_stats] = ranksum(dataArray{1}(nBin,:), ...
                                       dataArray{2}(nBin,:),'method','approximate');
        aoeVnd_z(nBin) = aoeVnd_stats.zval; 
        
        if aoeVnd_p(nBin) < postHoc_sig
            plot(1.45,bins2plot(nBin),'m*')
        else
        end
                    
        ndVpc_stats    = [];
        [ndVpc_p(nBin), ~ , ndVpc_stats]  = ranksum(dataArray{2}(nBin,:) , ...
                                       dataArray{3}(nBin,:),'method','approximate');
        ndVpc_z(nBin)  = ndVpc_stats.zval;
        
        if ndVpc_p(nBin) < postHoc_sig
            plot(1.4,bins2plot(nBin),'k*')
        else
        end
        
    end   
    pTask_R      = pTask;
    Fstat_R      = Fstat;
    degFreedom_R = degFreedom;
    aoeVnd_p_R   = aoeVnd_p;
    aoeVnd_z_R   = aoeVnd_z;
    aoeVpc_p_R   = aoeVpc_p;
    aoeVpc_z_R   = aoeVpc_z;
    ndVpc_p_R    = ndVpc_p;
    ndVpc_z_R    = ndVpc_z;   
else
end




%% panel 2F %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot view angle trajectory on left and right trials by task 
clearvars -except groupOrder accuracy_xMouse nTrialsSubselect_xTask nTrialsTotal_xTask ...
           yVel_xMouse xPosL_xMouse xPosR_xMouse viewAngL_xMouse viewAngR_xMouse distance_xMouse ...
           concatLog_OFF_aoe concatLog_OFF_nd concatLog_OFF_pc binIDs  
displayStats = 1; % 0 to suppress display of 1-way ANOVA of task on accuracy
                  % and ranksum unpaired comparisons between each task
ANOVA_sig    = 0.05; 
postHoc_sig  = 0.05;
         
bins2plot   = binIDs(2:end)-13;
plotOffsets = [0 5 -5];
dataArray   = []; data = []; xIdx = []; colorIdx = [];
dataArray   = viewAngL_xMouse;
data        = [viewAngL_xMouse{1}'; ...
               viewAngL_xMouse{2}'; ...
               viewAngL_xMouse{3}'];
xIdx        = [ones(size(viewAngL_xMouse{1},2),1);   ...
               ones(size(viewAngL_xMouse{2},2),1)*2; ...
               ones(size(viewAngL_xMouse{3},2),1)*3];           
                  
figure
subplot(1,2,1)
nTask = 1;
meanPlot = []; meanPlot = nanmean(dataArray{nTask},2);
errPlot  = []; errPlot  = (nanstd(dataArray{nTask},0,2)) / sqrt(size(dataArray{nTask},2)-1);
errorbar(meanPlot, bins2plot'+ plotOffsets(nTask) ,errPlot,'horizontal','k-'); hold on
text(-40,80,['AoE n= ' num2str(size(dataArray{nTask},2)) ],'Color','k','FontSize',10) 

nTask = 2;
meanPlot = []; meanPlot = nanmean(dataArray{nTask},2);
errPlot  = []; errPlot  = (nanstd(dataArray{nTask},0,2)) / sqrt(size(dataArray{nTask},2)-1);
errorbar(meanPlot, bins2plot'+ plotOffsets(nTask) ,errPlot,'horizontal','m-');
text(-40,70,['ctrl#1 n= ' num2str(size(dataArray{nTask},2)) ],'Color','m','FontSize',10) 

nTask = 3;
meanPlot = []; meanPlot = nanmean(dataArray{nTask},2);
errPlot  = []; errPlot  = (nanstd(dataArray{nTask},0,2)) / sqrt(size(dataArray{nTask},2)-1);
errorbar(meanPlot, bins2plot'+ plotOffsets(nTask) ,errPlot,'horizontal','c-');
text(-40,60,['ctrl#2 n= ' num2str(size(dataArray{nTask},2)) ],'Color','c','FontSize',10) 

plot(zeros(300,1),[1:300],'--','Color',[0.8 0.8 0.8])
title('left choice')
ylim([-0.5 300.5])
xlim([-45 5])
xlabel('view angle (deg)')
set(gca,'TickDir','out');
ylabel('y-position (cm)')
xticks([-45 -30 -15 0])
box off

if displayStats == 1    
    xIdx = categorical(xIdx);  
    pTask    = [];  Fstat    = [];  degFreedom = [];
    aoeVnd_p = [];  aoeVnd_z = []; 
    aoeVpc_p = [];  aoeVpc_z = []; 
    ndVpc_p  = [];  ndVpc_z  = [];  
    for nBin = 1:numel(bins2plot)
        p = []; t = []; stats = [];
        [p,t,stats]      = anova1(data(:,nBin),xIdx,'off');
        pTask(nBin)      = p; 
        Fstat(nBin)      = t{2,5}; 
        degFreedom(nBin) = stats.df;
        
        if pTask(nBin) < ANOVA_sig 
            plot(4,bins2plot(nBin),'k*')
        else
        end
              
        aoeVpc_stats   = [];
        [aoeVpc_p(nBin), ~ , aoeVpc_stats] = ranksum(dataArray{1}(nBin,:), ...
                                       dataArray{3}(nBin,:),'method','approximate');
        aoeVpc_z(nBin) = aoeVpc_stats.zval;  
        
        if aoeVpc_p(nBin) < postHoc_sig
            plot(-45,bins2plot(nBin),'c*')
        else
        end     
                                   
        aoeVnd_stats   = [];                           
        [aoeVnd_p(nBin), ~ , aoeVnd_stats] = ranksum(dataArray{1}(nBin,:), ...
                                       dataArray{2}(nBin,:),'method','approximate');
        aoeVnd_z(nBin) = aoeVnd_stats.zval; 
        
        if aoeVnd_p(nBin) < postHoc_sig
            plot(-44,bins2plot(nBin),'m*')
        else
        end
                    
        ndVpc_stats    = [];
        [ndVpc_p(nBin), ~ , ndVpc_stats]  = ranksum(dataArray{2}(nBin,:) , ...
                                       dataArray{3}(nBin,:),'method','approximate');
        ndVpc_z(nBin)  = ndVpc_stats.zval;
        
        if ndVpc_p(nBin) < postHoc_sig
            plot(-43,bins2plot(nBin),'k*')
        else
        end
        
    end   
    pTask_L      = pTask;
    Fstat_L      = Fstat;
    degFreedom_L = degFreedom;
    aoeVnd_p_L   = aoeVnd_p;
    aoeVnd_z_L   = aoeVnd_z;
    aoeVpc_p_L   = aoeVpc_p;
    aoeVpc_z_L   = aoeVpc_z;
    ndVpc_p_L    = ndVpc_p;
    ndVpc_z_L    = ndVpc_z;   
else
end


% right choice trials
bins2plot   = binIDs(2:end)-13;
plotOffsets = [0 5 -5];
dataArray   = []; data = []; xIdx = []; colorIdx = [];
dataArray   = viewAngR_xMouse;
data        = [viewAngR_xMouse{1}'; ...
               viewAngR_xMouse{2}'; ...
               viewAngR_xMouse{3}'];
xIdx        = [ones(size(viewAngR_xMouse{1},2),1);   ...
               ones(size(viewAngR_xMouse{2},2),1)*2; ...
               ones(size(viewAngR_xMouse{3},2),1)*3];           
                  
subplot(1,2,2)
nTask = 1;
meanPlot = []; meanPlot = nanmean(dataArray{nTask},2);
errPlot  = []; errPlot  = (nanstd(dataArray{nTask},0,2)) / sqrt(size(dataArray{nTask},2)-1);
errorbar(meanPlot, bins2plot'+ plotOffsets(nTask) ,errPlot,'horizontal','k-'); hold on

nTask = 2;
meanPlot = []; meanPlot = nanmean(dataArray{nTask},2);
errPlot  = []; errPlot  = (nanstd(dataArray{nTask},0,2)) / sqrt(size(dataArray{nTask},2)-1);
errorbar(meanPlot, bins2plot'+ plotOffsets(nTask) ,errPlot,'horizontal','m-');

nTask = 3;
meanPlot = []; meanPlot = nanmean(dataArray{nTask},2);
errPlot  = []; errPlot  = (nanstd(dataArray{nTask},0,2)) / sqrt(size(dataArray{nTask},2)-1);
errorbar(meanPlot, bins2plot'+ plotOffsets(nTask) ,errPlot,'horizontal','c-');

plot(zeros(300,1),[1:300],'--','Color',[0.8 0.8 0.8])
title('right choice')
ylim([-0.5 300.5])
xlim([-5 45])
xlabel('view angle (deg)')
box off
set(gca,'YAxisLocation','right');
xticks([0 15 30 45])

if displayStats == 1    
    xIdx = categorical(xIdx);  
    pTask    = [];  Fstat    = [];  degFreedom = [];
    aoeVnd_p = [];  aoeVnd_z = []; 
    aoeVpc_p = [];  aoeVpc_z = []; 
    ndVpc_p  = [];  ndVpc_z  = [];    
    
    for nBin = 1:numel(bins2plot)
        p = []; t = []; stats = [];
        [p,t,stats]      = anova1(data(:,nBin),xIdx,'off');
        pTask(nBin)      = p; 
        Fstat(nBin)      = t{2,5}; 
        degFreedom(nBin) = stats.df;
        
        if pTask(nBin) < ANOVA_sig 
            plot(-4,bins2plot(nBin),'k*')
        else
        end
              
        aoeVpc_stats   = [];
        [aoeVpc_p(nBin), ~ , aoeVpc_stats] = ranksum(dataArray{1}(nBin,:), ...
                                       dataArray{3}(nBin,:),'method','approximate');
        aoeVpc_z(nBin) = aoeVpc_stats.zval;  
        
        if aoeVpc_p(nBin) < postHoc_sig
            plot(1.5,bins2plot(nBin),'c*')
        else
        end     
                                   
        aoeVnd_stats   = [];                           
        [aoeVnd_p(nBin), ~ , aoeVnd_stats] = ranksum(dataArray{1}(nBin,:), ...
                                       dataArray{2}(nBin,:),'method','approximate');
        aoeVnd_z(nBin) = aoeVnd_stats.zval; 
        
        if aoeVnd_p(nBin) < postHoc_sig
            plot(44,bins2plot(nBin),'m*')
        else
        end
                    
        ndVpc_stats    = [];
        [ndVpc_p(nBin), ~ , ndVpc_stats]  = ranksum(dataArray{2}(nBin,:) , ...
                                       dataArray{3}(nBin,:),'method','approximate');
        ndVpc_z(nBin)  = ndVpc_stats.zval;
        
        if ndVpc_p(nBin) < postHoc_sig
            plot(43,bins2plot(nBin),'k*')
        else
        end
        
    end   
    pTask_R      = pTask;
    Fstat_R      = Fstat;
    degFreedom_R = degFreedom;
    aoeVnd_p_R   = aoeVnd_p;
    aoeVnd_z_R   = aoeVnd_z;
    aoeVpc_p_R   = aoeVpc_p;
    aoeVpc_z_R   = aoeVpc_z;
    ndVpc_p_R    = ndVpc_p;
    ndVpc_z_R    = ndVpc_z;   
else
end                


%% panel 2G %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot distance by task 
clearvars -except groupOrder accuracy_xMouse nTrialsSubselect_xTask nTrialsTotal_xTask ...
           yVel_xMouse xPosL_xMouse xPosR_xMouse viewAngL_xMouse viewAngR_xMouse distance_xMouse ...
           concatLog_OFF_aoe concatLog_OFF_nd concatLog_OFF_pc binIDs 
displayStats = 1; % 0 to suppress display of 1-way ANOVA of task on accuracy
                  % and ranksum unpaired comparisons between each task
                 
                  
figure
dataArray   = []; data = []; xIdx = []; colorIdx = [];
dataArray   = distance_xMouse;
data        = [distance_xMouse{1}'; ...
               distance_xMouse{2}'; ...
               distance_xMouse{3}'];
xIdx        = [ones(size(distance_xMouse{1}))';   ...
               ones(size(distance_xMouse{2}))'*2; ...
               ones(size(distance_xMouse{3}))'*3];
catMarker   = {'x','x', 'x'};
colorIdx    = {[0.8 0.8 0.8],[0.99 0.7 0.99], [0.7 0.99 0.99]};
hold on
plotSpread(data,'categoryIdx',xIdx,'distributionIdx',xIdx,...
    'categoryColors',colorIdx,'categoryMarkers',catMarker);

nTask = 1;
h = errorbar(nTask,nanmean(dataArray{nTask}),nanstd(dataArray{nTask})./sqrt(numel(dataArray{nTask})-1), ...
    '.-','linewidth',4,'color','k','markersize',0.1); h.CapSize = 16;
text(0.8, 545, [num2str(numel(dataArray{nTask})) ' mice' ])
text(0.8, 535, [num2str(nTrialsTotal_xTask(nTask)) ' trials' ])

nTask = 2;
h = errorbar(nTask,nanmean(dataArray{nTask}),nanstd(dataArray{nTask})./sqrt(numel(dataArray{nTask})-1), ...
    '.-','linewidth',4,'color','m','markersize',0.1); h.CapSize = 16;
text(1.8, 545, [num2str(numel(dataArray{nTask})) ' mice' ])
text(1.8, 535, [num2str(nTrialsTotal_xTask(nTask)) ' trials' ])

nTask = 3;
h = errorbar(nTask,nanmean(dataArray{nTask}),nanstd(dataArray{nTask})./sqrt(numel(dataArray{nTask})-1), ...
    '.-','linewidth',4,'color','c','markersize',0.1); h.CapSize = 16;
text(2.8, 545, [num2str(numel(dataArray{nTask})) ' mice' ])
text(2.8, 535, [num2str(nTrialsTotal_xTask(nTask)) ' trials' ])
set(gca,'xtick',1:3,'xticklabel',{'AoE';'ctrl#1';'ctrl#2'})
rotateXLabels(gca,45)
ylabel('distance (cm)','fontsize',12)
xlim([.65 3.35])
ylim([300 550])
box off
set(gca,'TickDir','out');

if displayStats == 1
    xIdx = categorical(xIdx);
    [p,t,stats] = anova1(data,xIdx,'off');
    text(1.2,400,['ANOVA p= ' num2str(p)])
    text(1.2,390,['df= ' num2str(stats.df)])
    text(1.2,380,['F= ' num2str((t{2,5}))])
    
    [aoeVpc_p, ~ , aoeVpc_stats] = ranksum(dataArray{1}, dataArray{3},'method','approximate');
    [aoeVnd_p, ~ , aoeVnd_stats] = ranksum(dataArray{1}, dataArray{2},'method','approximate');
    [ndVpc_p , ~ , ndVpc_stats]  = ranksum(dataArray{2} , dataArray{3},'method','approximate');
    
    plot(1:3,[500,500,500],'k-')
    text(1,505,['p= ' num2str(aoeVpc_p)])
    text(1.8,505,['z= ' num2str(aoeVpc_stats.zval)])
    text(2.8,505,['df= ' num2str(numel(dataArray{1})+numel(dataArray{3})-1)])
    
    plot(1:2,[480,480],'k-')
    text(1,485,['p= ' num2str(aoeVnd_p)])
    text(1,475,['z= ' num2str(aoeVnd_stats.zval)])
    text(1,465,['df= ' num2str(numel(dataArray{1})+numel(dataArray{2})-1)])
    
    plot(2:3,[460,460],'k-')
    text(2.5,465,['p= ' num2str(ndVpc_p)])
    text(2.5,455,['z= ' num2str(ndVpc_stats.zval)])
    text(2.5,445,['df= ' num2str(numel(dataArray{2})+numel(dataArray{3})-1)])
else
end

%%
clearvars -except groupOrder accuracy_xMouse nTrialsSubselect_xTask nTrialsTotal_xTask ...
           yVel_xMouse xPosL_xMouse xPosR_xMouse viewAngL_xMouse viewAngR_xMouse distance_xMouse ...
           concatLog_OFF_aoe concatLog_OFF_nd concatLog_OFF_pc binIDs