%% Script to generate ExtData Fig 7C,F,I,L data plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% cross task psychometrics with indirect and direct pathway inhibition DMS and NAc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
cd(globalParams.dataPath)
load('OffOn_TasksOrGroup_all.mat')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% analyzes log files for opto induced choice bias
% (1) cue statistics of aoe (log.currMaze==10) and pc (log.currMaze==7) pre-selected to match
% (2) all analyses subselect trial blocks with performance >60% and trials with excessTravel<0.1 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars psychPerfByMouse_OFF psychPerfBinsByMouse_OFF ...
          psychPerfByMouse_ON psychPerfBinsByMouse_ON

groupOrder = {'DMSD2_aoe'  ,'DMSD2_nd'  ,'DMSD2_pc', ...
              'DMSD1_aoe'  ,'DMSD1_nd'  ,'DMSD1_pc', ...
              'DMSCtrl_aoe','DMSCtrl_nd','DMSCtrl_pc', ...
              'NACD2_aoe'  ,'NACD1_aoe' ,'NACCtrl_aoe'};
          
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

    % remove whole trial inhibition if present in lg
    if any(contains(unique(lg.laserEpoch),"whole"))
        removeInd = []; 
        removeInd = ismember(lg.laserEpoch,"whole");
        lg       = structfun(@(x) x(~removeInd),lg,'UniformOutput',false);
        removeInd = []; 
        removeInd = ismember(lgClean.laserEpoch,"whole");
        lgClean  = structfun(@(x) x(~removeInd),lgClean,'UniformOutput',false);
    else
    end
    
    % invert logs to ipsi contra space
    lg        = invertLogs(lg);
    lgClean   = invertLogs(lgClean);
    
    % get binned performance per mouse of data using psychometricFit 
    binnedPerfOFF = [];
    binnedPerfON  = [];
    xBinsOFF   = [];
    xBinsON    = [];
    mouseName  = []; 
    mouseName  = unique(lgClean.mouseID);
    psychBins  = -15:5:15;
    for nMouse = 1:numel(mouseName)
        
        psychOFF   = [];
        psychON    = [];
        psychOFF                = psychometricFit(lgClean.choice...
                                 (lgClean.mouseID==mouseName(nMouse) & lgClean.laserON==0),...
                                  lgClean.nCues_RminusL(lgClean.mouseID==mouseName(nMouse)& lgClean.laserON==0),...
                                  1,psychBins,1,1);
        % if less than 10 trials in a bin fill perf with NAN
        psychOFF.perfPsych_binned(psychOFF.nTrials_binned<10) = NaN;                    
        psychOFF.perfPsych_xaxisBins(psychOFF.nTrials_binned<10) = NaN; 
        binnedPerfOFF(nMouse,:) = psychOFF.perfPsych_binned;
        xBinsOFF(nMouse,:)      = psychOFF.perfPsych_xaxisBins;
        
        
        psychON                 = psychometricFit(lgClean.choice...
                                  (lgClean.mouseID==mouseName(nMouse) & lgClean.laserON==1),...
                                   lgClean.nCues_RminusL(lgClean.mouseID==mouseName(nMouse)& lgClean.laserON==1),...
                                   1,psychBins,1,1);
        % if less than 10 trials in a bin fill perf with NAN                       
        psychON.perfPsych_binned(psychON.nTrials_binned<10) = NaN;                    
        psychON.perfPsych_xaxisBins(psychON.nTrials_binned<10) = NaN; 
        binnedPerfON(nMouse,:) = psychON.perfPsych_binned;
        xBinsON(nMouse,:)      = psychON.perfPsych_xaxisBins;                       
                                
    end
    
    % save per mouse binned psychometric performance (actual data average, not fits)
    psychPerfByMouse_OFF{counter}     = binnedPerfOFF *100;
    psychPerfBinsByMouse_OFF{counter} = xBinsOFF;
    
    psychPerfByMouse_ON{counter}     = binnedPerfON*100;
    psychPerfBinsByMouse_ON{counter} = xBinsON;
    
    % count up for next log
    counter = counter+1;
    
end

%% clearvars except source data for plotting
clearvars -except lgDMSD2_aoe lgDMSD1_aoe lgDMSCtrl_aoe ...
                  lgDMSD2_nd lgDMSD1_nd lgDMSCtrl_nd ...
                  lgDMSD2_pc lgDMSD1_pc lgDMSCtrl_pc ...
                  lgNACD2_aoe lgNACD1_aoe lgNACCtrl_aoe ...
                  groupOrder psychPerfByMouse_OFF psychPerfBinsByMouse_OFF ...
                  psychPerfByMouse_ON psychPerfBinsByMouse_ON
              
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%              
%% PLOTTING of Ext Data Fig 7C,F,I,L
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
% update these 2 variables to determine groups to plot - the name and group 
% numbers should match the ordering in groupOrder variable to be accurate
% commented pairings are correct match for plots in panels ExtData Fig 7
names2plot  = {'DMSD2_aoe'  ,'DMSD2_nd'   ,'DMSD2_pc'};    % panel C- indirect DMS: aoe, ctrl#1, ctrl#2
%             {'DMSD2_aoe'  ,'DMSD2_nd'   ,'DMSD2_pc'};    % panel C- indirect DMS: aoe, ctrl#1, ctrl#2   
%             {'DMSD1_aoe'  ,'DMSD1_nd'   ,'DMSD1_pc'};    % panel F- direct DMS: aoe, ctrl#1, ctrl#2
%             {'DMSCtrl_aoe','DMSCtrl_nd' ,'DMSCtrl_pc'};  % panel I- no opsin DMS: aoe, ctrl#1, ctrl#2  
%             {'NACD2_aoe'  ,'NACD1_aoe'  ,'NACCtrl_aoe'}; % panel L- NAc: indirect, direct, no opsin    

groups2plot = [1 2 3];    % panel C- indirect DMS: aoe, ctrl#1, ctrl#2
%             [1 2 3];    % panel C- indirect DMS: aoe, ctrl#1, ctrl#2
%             [4 5 6];    % panel F- direct DMS: aoe, ctrl#1, ctrl#2
%             [7 8 9];    % panel I- no opsin DMS: aoe, ctrl#1, ctrl#2  
%             [10 11 12]; % panel L- NAc: indirect, direct, no opsin              

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
barOFF      = 'k';   
barON       = [.1 .9 .1];
indivOFF    = [0.8 0.8 0.8]; 
indivON     = [.6 0.9 0.6]; 
xLimits     = [-14.5 14.5];
yLimits     = [0 100];

for nGroup =1:numel(groups2plot)
    subplot(1,3,nGroup)
    
    plot(psychPerfBinsByMouse_OFF{1,groups2plot(nGroup)}',psychPerfByMouse_OFF{1,groups2plot(nGroup)}','Color',indivOFF); hold on
    plot(psychPerfBinsByMouse_ON{1,groups2plot(nGroup)}',psychPerfByMouse_ON{1,groups2plot(nGroup)}','Color',indivON)
    xlim(xLimits)
    ylim(yLimits)
    
    % plot laser on cross mouse average
    x = []; y = []; ySEM = []; l = []; u = [];
    x    = nanmean(psychPerfBinsByMouse_ON{1,groups2plot(nGroup)},1);
    y    = nanmean(psychPerfByMouse_ON{1,groups2plot(nGroup)},1);
    ySEM = nanstd(psychPerfByMouse_ON{1,groups2plot(nGroup)},1)/(sqrt(size(psychPerfByMouse_ON{1,groups2plot(nGroup)},1)-1));
    l    = -ySEM;
    u    = ySEM;
    h    = errorbar(x, y, l, u,'.-','linewidth',2,'color',barON,'markersize',0.1); h.CapSize = 8; hold on
    
    % plot laser off cross mouse average
    x = []; y = []; ySEM = []; l = []; u = [];
    x    = nanmean(psychPerfBinsByMouse_OFF{1,groups2plot(nGroup)},1);
    y    = nanmean(psychPerfByMouse_OFF{1,groups2plot(nGroup)},1);
    ySEM = nanstd(psychPerfByMouse_OFF{1,groups2plot(nGroup)},1)/(sqrt(size(psychPerfByMouse_OFF{1,groups2plot(nGroup)},1)-1));
    l    = -ySEM;
    u    = ySEM;
    h    = errorbar(x, y, l, u,'.-','linewidth',2,'color',barOFF,'markersize',0.1); h.CapSize = 8; hold on
    
    plot(xLimits(1):xLimits(2),ones(1,numel(xLimits(1):xLimits(2)))*50,'--','Color',[0.5 0.5 0.5])
    plot(zeros(1,numel(yLimits(1):yLimits(2))),yLimits(1):yLimits(2),'--','Color',[0.5 0.5 0.5])
    set(gca,'TickDir','out');
    box off
    xticks([-10 0 10])
    yticks([0 20 40 60 80 100])
    xlabel('#Contra-Ipsi')
    ylabel('went contra (%)')
    
    if contains(names2plot(nGroup),'DMS')
        if contains(names2plot(nGroup),'aoe')
            title('AoE')
        elseif contains(names2plot(nGroup),'nd')
            title('ctrl#1')
        elseif contains(names2plot(nGroup),'pc')
            title('ctrl#2')
        else
        end
    elseif contains(names2plot(nGroup),'NAC')
        if contains(names2plot(nGroup),'D2')
            title('indirect')
        elseif contains(names2plot(nGroup),'D1')
            title('direct')
        elseif contains(names2plot(nGroup),'Ctrl')
            title('no opsin')
        else
        end
    else
    end
end

%% clearvars except source data for plotting
clearvars -except lgDMSD2_aoe lgDMSD1_aoe lgDMSCtrl_aoe ...
                  lgDMSD2_nd lgDMSD1_nd lgDMSCtrl_nd ...
                  lgDMSD2_pc lgDMSD1_pc lgDMSCtrl_pc ...
                  lgNACD2_aoe lgNACD1_aoe lgNACCtrl_aoe ...
                  groupOrder psychPerfByMouse_OFF psychPerfBinsByMouse_OFF ...
                  psychPerfByMouse_ON psychPerfBinsByMouse_ON
              
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%              
%% PLOTTING of Ext Data Fig 7C,F,I,L
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
% update these 2 variables to determine groups to plot - the name and group 
% numbers should match the ordering in groupOrder variable to be accurate
% commented pairings are correct match for plots in panels ExtData Fig 7
names2plot  = {'DMSD1_aoe'  ,'DMSD1_nd'   ,'DMSD1_pc'};    % panel F- direct DMS: aoe, ctrl#1, ctrl#2
%             {'DMSD2_aoe'  ,'DMSD2_nd'   ,'DMSD2_pc'};    % panel C- indirect DMS: aoe, ctrl#1, ctrl#2   
%             {'DMSD1_aoe'  ,'DMSD1_nd'   ,'DMSD1_pc'};    % panel F- direct DMS: aoe, ctrl#1, ctrl#2
%             {'DMSCtrl_aoe','DMSCtrl_nd' ,'DMSCtrl_pc'};  % panel I- no opsin DMS: aoe, ctrl#1, ctrl#2  
%             {'NACD2_aoe'  ,'NACD1_aoe'  ,'NACCtrl_aoe'}; % panel L- NAc: indirect, direct, no opsin    

groups2plot = [4 5 6];    % panel F- direct DMS: aoe, ctrl#1, ctrl#2
%             [1 2 3];    % panel C- indirect DMS: aoe, ctrl#1, ctrl#2
%             [4 5 6];    % panel F- direct DMS: aoe, ctrl#1, ctrl#2
%             [7 8 9];    % panel I- no opsin DMS: aoe, ctrl#1, ctrl#2  
%             [10 11 12]; % panel L- NAc: indirect, direct, no opsin              

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
barOFF      = 'k';   
barON       = [.1 .9 .1];
indivOFF    = [0.8 0.8 0.8]; 
indivON     = [.6 0.9 0.6]; 
xLimits     = [-14.5 14.5];
yLimits     = [0 100];

for nGroup =1:numel(groups2plot)
    subplot(1,3,nGroup)
    
    plot(psychPerfBinsByMouse_OFF{1,groups2plot(nGroup)}',psychPerfByMouse_OFF{1,groups2plot(nGroup)}','Color',indivOFF); hold on
    plot(psychPerfBinsByMouse_ON{1,groups2plot(nGroup)}',psychPerfByMouse_ON{1,groups2plot(nGroup)}','Color',indivON)
    xlim(xLimits)
    ylim(yLimits)
    
    % plot laser on cross mouse average
    x = []; y = []; ySEM = []; l = []; u = [];
    x    = nanmean(psychPerfBinsByMouse_ON{1,groups2plot(nGroup)},1);
    y    = nanmean(psychPerfByMouse_ON{1,groups2plot(nGroup)},1);
    ySEM = nanstd(psychPerfByMouse_ON{1,groups2plot(nGroup)},1)/(sqrt(size(psychPerfByMouse_ON{1,groups2plot(nGroup)},1)-1));
    l    = -ySEM;
    u    = ySEM;
    h    = errorbar(x, y, l, u,'.-','linewidth',2,'color',barON,'markersize',0.1); h.CapSize = 8; hold on
    
    % plot laser off cross mouse average
    x = []; y = []; ySEM = []; l = []; u = [];
    x    = nanmean(psychPerfBinsByMouse_OFF{1,groups2plot(nGroup)},1);
    y    = nanmean(psychPerfByMouse_OFF{1,groups2plot(nGroup)},1);
    ySEM = nanstd(psychPerfByMouse_OFF{1,groups2plot(nGroup)},1)/(sqrt(size(psychPerfByMouse_OFF{1,groups2plot(nGroup)},1)-1));
    l    = -ySEM;
    u    = ySEM;
    h    = errorbar(x, y, l, u,'.-','linewidth',2,'color',barOFF,'markersize',0.1); h.CapSize = 8; hold on
    
    plot(xLimits(1):xLimits(2),ones(1,numel(xLimits(1):xLimits(2)))*50,'--','Color',[0.5 0.5 0.5])
    plot(zeros(1,numel(yLimits(1):yLimits(2))),yLimits(1):yLimits(2),'--','Color',[0.5 0.5 0.5])
    set(gca,'TickDir','out');
    box off
    xticks([-10 0 10])
    yticks([0 20 40 60 80 100])
    xlabel('#Contra-Ipsi')
    ylabel('went contra (%)')
    
    if contains(names2plot(nGroup),'DMS')
        if contains(names2plot(nGroup),'aoe')
            title('AoE')
        elseif contains(names2plot(nGroup),'nd')
            title('ctrl#1')
        elseif contains(names2plot(nGroup),'pc')
            title('ctrl#2')
        else
        end
    elseif contains(names2plot(nGroup),'NAC')
        if contains(names2plot(nGroup),'D2')
            title('indirect')
        elseif contains(names2plot(nGroup),'D1')
            title('direct')
        elseif contains(names2plot(nGroup),'Ctrl')
            title('no opsin')
        else
        end
    else
    end
end

%% clearvars except source data for plotting
clearvars -except lgDMSD2_aoe lgDMSD1_aoe lgDMSCtrl_aoe ...
                  lgDMSD2_nd lgDMSD1_nd lgDMSCtrl_nd ...
                  lgDMSD2_pc lgDMSD1_pc lgDMSCtrl_pc ...
                  lgNACD2_aoe lgNACD1_aoe lgNACCtrl_aoe ...
                  groupOrder psychPerfByMouse_OFF psychPerfBinsByMouse_OFF ...
                  psychPerfByMouse_ON psychPerfBinsByMouse_ON
              
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%              
%% PLOTTING of Ext Data Fig 7C,F,I,L
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
% update these 2 variables to determine groups to plot - the name and group 
% numbers should match the ordering in groupOrder variable to be accurate
% commented pairings are correct match for plots in panels ExtData Fig 7
names2plot  = {'DMSCtrl_aoe','DMSCtrl_nd' ,'DMSCtrl_pc'};  % panel I- no opsin DMS: aoe, ctrl#1, ctrl#
%             {'DMSD2_aoe'  ,'DMSD2_nd'   ,'DMSD2_pc'};    % panel C- indirect DMS: aoe, ctrl#1, ctrl#2   
%             {'DMSD1_aoe'  ,'DMSD1_nd'   ,'DMSD1_pc'};    % panel F- direct DMS: aoe, ctrl#1, ctrl#2
%             {'DMSCtrl_aoe','DMSCtrl_nd' ,'DMSCtrl_pc'};  % panel I- no opsin DMS: aoe, ctrl#1, ctrl#2  
%             {'NACD2_aoe'  ,'NACD1_aoe'  ,'NACCtrl_aoe'}; % panel L- NAc: indirect, direct, no opsin    

groups2plot = [7 8 9];    % panel I- no opsin DMS: aoe, ctrl#1, ctrl#2  
%             [1 2 3];    % panel C- indirect DMS: aoe, ctrl#1, ctrl#2
%             [4 5 6];    % panel F- direct DMS: aoe, ctrl#1, ctrl#2
%             [7 8 9];    % panel I- no opsin DMS: aoe, ctrl#1, ctrl#2  
%             [10 11 12]; % panel L- NAc: indirect, direct, no opsin              

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
barOFF      = 'k';   
barON       = [.1 .9 .1];
indivOFF    = [0.8 0.8 0.8]; 
indivON     = [.6 0.9 0.6]; 
xLimits     = [-14.5 14.5];
yLimits     = [0 100];

for nGroup =1:numel(groups2plot)
    subplot(1,3,nGroup)
    
    plot(psychPerfBinsByMouse_OFF{1,groups2plot(nGroup)}',psychPerfByMouse_OFF{1,groups2plot(nGroup)}','Color',indivOFF); hold on
    plot(psychPerfBinsByMouse_ON{1,groups2plot(nGroup)}',psychPerfByMouse_ON{1,groups2plot(nGroup)}','Color',indivON)
    xlim(xLimits)
    ylim(yLimits)
    
    % plot laser on cross mouse average
    x = []; y = []; ySEM = []; l = []; u = [];
    x    = nanmean(psychPerfBinsByMouse_ON{1,groups2plot(nGroup)},1);
    y    = nanmean(psychPerfByMouse_ON{1,groups2plot(nGroup)},1);
    ySEM = nanstd(psychPerfByMouse_ON{1,groups2plot(nGroup)},1)/(sqrt(size(psychPerfByMouse_ON{1,groups2plot(nGroup)},1)-1));
    l    = -ySEM;
    u    = ySEM;
    h    = errorbar(x, y, l, u,'.-','linewidth',2,'color',barON,'markersize',0.1); h.CapSize = 8; hold on
    
    % plot laser off cross mouse average
    x = []; y = []; ySEM = []; l = []; u = [];
    x    = nanmean(psychPerfBinsByMouse_OFF{1,groups2plot(nGroup)},1);
    y    = nanmean(psychPerfByMouse_OFF{1,groups2plot(nGroup)},1);
    ySEM = nanstd(psychPerfByMouse_OFF{1,groups2plot(nGroup)},1)/(sqrt(size(psychPerfByMouse_OFF{1,groups2plot(nGroup)},1)-1));
    l    = -ySEM;
    u    = ySEM;
    h    = errorbar(x, y, l, u,'.-','linewidth',2,'color',barOFF,'markersize',0.1); h.CapSize = 8; hold on
    
    plot(xLimits(1):xLimits(2),ones(1,numel(xLimits(1):xLimits(2)))*50,'--','Color',[0.5 0.5 0.5])
    plot(zeros(1,numel(yLimits(1):yLimits(2))),yLimits(1):yLimits(2),'--','Color',[0.5 0.5 0.5])
    set(gca,'TickDir','out');
    box off
    xticks([-10 0 10])
    yticks([0 20 40 60 80 100])
    xlabel('#Contra-Ipsi')
    ylabel('went contra (%)')
    
    if contains(names2plot(nGroup),'DMS')
        if contains(names2plot(nGroup),'aoe')
            title('AoE')
        elseif contains(names2plot(nGroup),'nd')
            title('ctrl#1')
        elseif contains(names2plot(nGroup),'pc')
            title('ctrl#2')
        else
        end
    elseif contains(names2plot(nGroup),'NAC')
        if contains(names2plot(nGroup),'D2')
            title('indirect')
        elseif contains(names2plot(nGroup),'D1')
            title('direct')
        elseif contains(names2plot(nGroup),'Ctrl')
            title('no opsin')
        else
        end
    else
    end
end

%% clearvars except source data for plotting
clearvars -except lgDMSD2_aoe lgDMSD1_aoe lgDMSCtrl_aoe ...
                  lgDMSD2_nd lgDMSD1_nd lgDMSCtrl_nd ...
                  lgDMSD2_pc lgDMSD1_pc lgDMSCtrl_pc ...
                  lgNACD2_aoe lgNACD1_aoe lgNACCtrl_aoe ...
                  groupOrder psychPerfByMouse_OFF psychPerfBinsByMouse_OFF ...
                  psychPerfByMouse_ON psychPerfBinsByMouse_ON
              
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%              
%% PLOTTING of Ext Data Fig 7C,F,I,L
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
% update these 2 variables to determine groups to plot - the name and group 
% numbers should match the ordering in groupOrder variable to be accurate
% commented pairings are correct match for plots in panels ExtData Fig 7
names2plot  = {'NACD2_aoe'  ,'NACD1_aoe'  ,'NACCtrl_aoe'}; % panel L- NAc: indirect, direct, no opsin    
%             {'DMSD2_aoe'  ,'DMSD2_nd'   ,'DMSD2_pc'};    % panel C- indirect DMS: aoe, ctrl#1, ctrl#2   
%             {'DMSD1_aoe'  ,'DMSD1_nd'   ,'DMSD1_pc'};    % panel F- direct DMS: aoe, ctrl#1, ctrl#2
%             {'DMSCtrl_aoe','DMSCtrl_nd' ,'DMSCtrl_pc'};  % panel I- no opsin DMS: aoe, ctrl#1, ctrl#2  
%             {'NACD2_aoe'  ,'NACD1_aoe'  ,'NACCtrl_aoe'}; % panel L- NAc: indirect, direct, no opsin    

groups2plot = [10 11 12]; % panel L- NAc: indirect, direct, no opsin  
%             [1 2 3];    % panel C- indirect DMS: aoe, ctrl#1, ctrl#2
%             [4 5 6];    % panel F- direct DMS: aoe, ctrl#1, ctrl#2
%             [7 8 9];    % panel I- no opsin DMS: aoe, ctrl#1, ctrl#2  
%             [10 11 12]; % panel L- NAc: indirect, direct, no opsin              

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
barOFF      = 'k';   
barON       = [.1 .9 .1];
indivOFF    = [0.8 0.8 0.8]; 
indivON     = [.6 0.9 0.6]; 
xLimits     = [-14.5 14.5];
yLimits     = [0 100];

for nGroup =1:numel(groups2plot)
    subplot(1,3,nGroup)
    
    plot(psychPerfBinsByMouse_OFF{1,groups2plot(nGroup)}',psychPerfByMouse_OFF{1,groups2plot(nGroup)}','Color',indivOFF); hold on
    plot(psychPerfBinsByMouse_ON{1,groups2plot(nGroup)}',psychPerfByMouse_ON{1,groups2plot(nGroup)}','Color',indivON)
    xlim(xLimits)
    ylim(yLimits)
    
    % plot laser on cross mouse average
    x = []; y = []; ySEM = []; l = []; u = [];
    x    = nanmean(psychPerfBinsByMouse_ON{1,groups2plot(nGroup)},1);
    y    = nanmean(psychPerfByMouse_ON{1,groups2plot(nGroup)},1);
    ySEM = nanstd(psychPerfByMouse_ON{1,groups2plot(nGroup)},1)/(sqrt(size(psychPerfByMouse_ON{1,groups2plot(nGroup)},1)-1));
    l    = -ySEM;
    u    = ySEM;
    h    = errorbar(x, y, l, u,'.-','linewidth',2,'color',barON,'markersize',0.1); h.CapSize = 8; hold on
    
    % plot laser off cross mouse average
    x = []; y = []; ySEM = []; l = []; u = [];
    x    = nanmean(psychPerfBinsByMouse_OFF{1,groups2plot(nGroup)},1);
    y    = nanmean(psychPerfByMouse_OFF{1,groups2plot(nGroup)},1);
    ySEM = nanstd(psychPerfByMouse_OFF{1,groups2plot(nGroup)},1)/(sqrt(size(psychPerfByMouse_OFF{1,groups2plot(nGroup)},1)-1));
    l    = -ySEM;
    u    = ySEM;
    h    = errorbar(x, y, l, u,'.-','linewidth',2,'color',barOFF,'markersize',0.1); h.CapSize = 8; hold on
    
    plot(xLimits(1):xLimits(2),ones(1,numel(xLimits(1):xLimits(2)))*50,'--','Color',[0.5 0.5 0.5])
    plot(zeros(1,numel(yLimits(1):yLimits(2))),yLimits(1):yLimits(2),'--','Color',[0.5 0.5 0.5])
    set(gca,'TickDir','out');
    box off
    xticks([-10 0 10])
    yticks([0 20 40 60 80 100])
    xlabel('#Contra-Ipsi')
    ylabel('went contra (%)')
    
    if contains(names2plot(nGroup),'DMS')
        if contains(names2plot(nGroup),'aoe')
            title('AoE')
        elseif contains(names2plot(nGroup),'nd')
            title('ctrl#1')
        elseif contains(names2plot(nGroup),'pc')
            title('ctrl#2')
        else
        end
    elseif contains(names2plot(nGroup),'NAC')
        if contains(names2plot(nGroup),'D2')
            title('indirect')
        elseif contains(names2plot(nGroup),'D1')
            title('direct')
        elseif contains(names2plot(nGroup),'Ctrl')
            title('no opsin')
        else
        end
    else
    end
end

%% clearvars except source data for plotting
clearvars -except lgDMSD2_aoe lgDMSD1_aoe lgDMSCtrl_aoe ...
                  lgDMSD2_nd lgDMSD1_nd lgDMSCtrl_nd ...
                  lgDMSD2_pc lgDMSD1_pc lgDMSCtrl_pc ...
                  lgNACD2_aoe lgNACD1_aoe lgNACCtrl_aoe ...
                  groupOrder psychPerfByMouse_OFF psychPerfBinsByMouse_OFF ...
                  psychPerfByMouse_ON psychPerfBinsByMouse_ON
