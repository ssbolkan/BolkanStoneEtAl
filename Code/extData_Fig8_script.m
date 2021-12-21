%% Script to generate ExtData Fig 12A-I data plots
% related to Fig 6 and 7
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load DMS logs with state probability and most likely state indices 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
cd(globalParams.dataPath)
load('DMS_AoE_GLMHMM_states.mat')

%% analysis of reward and rewardRate at all within-session transitions
reward = []; rewardRate = [];
groupOrder = {'DMSD1', 'DMSD2' };
counter = 1;
for groupNum = 1:numel(groupOrder)
    clear lg  
    lg       = eval(['lg' groupOrder{groupNum}]); 
    
    % initialize and get reward delivered on trial, if choice was correct
    rewardGiven = []; 
    rewardGiven = .004*lg.rewardScale; % converted to mL. 4ul*rewardScale reflects actual size of reward per trial
    accuracy    = []; 
    accuracy    = lg.choice == lg.trialType; % if choice matches reward side
    
    % find all trial indices where most likely state changes
    %  and then remove all indicies that are between session transitions
    removeInd = []; transInd = [];
    transInd = find(diff(lg.stateIndex)~=0);
    for transID = 1:numel(transInd)
        if lg.sessionID(transInd(transID)) ~= lg.sessionID(transInd(transID)+1) || ...
                lg.sessionID(transInd(transID)) ~= lg.sessionID(transInd(transID)-1)
            removeInd = [removeInd transID];
        else
        end
    end
    transInd(removeInd) = [];
    
    % initialize counters and vars to save for 6 possible transition types
    s1TOs3_COUNTER = 0; s1TOs3_reward = []; s1TOs3_rewardRate = [];
    s1TOs2_COUNTER = 0; s1TOs2_reward = []; s1TOs2_rewardRate = [];
    s2TOs3_COUNTER = 0; s2TOs3_reward = []; s2TOs3_rewardRate = [];
    s2TOs1_COUNTER = 0; s2TOs1_reward = []; s2TOs1_rewardRate = [];
    s3TOs1_COUNTER = 0; s3TOs1_reward = []; s3TOs1_rewardRate = [];
    s3TOs2_COUNTER = 0; s3TOs2_reward = []; s3TOs2_rewardRate = [];
    
    % loop through all within session transitions
    for transID = 1:numel(transInd)
        % if transition from state 1-->2
        if lg.stateIndex(transInd(transID)) == 1 && lg.stateIndex(transInd(transID)+1) == 2
            startInd = []; startInd = find(diff(lg.sessionID)~=0) ; % get index of all between sess transitions
            startInd = startInd(startInd<transInd(transID)); % select all between sess transitions before current state transition
            startInd = startInd(end) + 1; % the closest of these (+1 trial) is the first trial of session with current state transition
            s1TOs2_COUNTER = s1TOs2_COUNTER+1;
            toSumRewardTemp = []; toSumAccuracyTemp = []; toSumTrialDurTemp = [];
            toSumRewardTemp    = rewardGiven(startInd:transInd(transID)); % scaled rewardGiven per trial up to transition (mL)
            toSumAccuracyTemp  = accuracy(startInd:transInd(transID)); % whether trials up to transition was rewarded
            toSumTrialDurTemp  = lg.trialDurFull(startInd:transInd(transID)); % per trial durations including ITI up to transition (in seconds)
            s1TOs2_reward(s1TOs2_COUNTER) = sum(toSumRewardTemp(toSumAccuracyTemp==1)); % total reward (mL) from sess start to transition
            s1TOs2_rewardRate(s1TOs2_COUNTER) = sum(toSumRewardTemp(toSumAccuracyTemp==1))/(sum(toSumTrialDurTemp)/60)*1000; % reward (uL) / sess duration (mL) up to transition
            
            % if transition from state 1-->3
        elseif lg.stateIndex(transInd(transID)) == 1 && lg.stateIndex(transInd(transID)+1) == 3
            startInd = []; startInd = find(diff(lg.sessionID)~=0) ; % get index of all between sess transitions
            startInd = startInd(startInd<transInd(transID)); % select all between sess transitions before current state transition
            startInd = startInd(end) + 1; % the closest of these (+1 trial) is the first trial of session with current state transition
            s1TOs3_COUNTER = s1TOs3_COUNTER+1;
            toSumRewardTemp = []; toSumAccuracyTemp = []; toSumTrialDurTemp = [];
            toSumRewardTemp    = rewardGiven(startInd:transInd(transID)); % scaled rewardGiven per trial up to transition (mL)
            toSumAccuracyTemp  = accuracy(startInd:transInd(transID)); % whether trials up to transition was rewarded
            toSumTrialDurTemp  = lg.trialDurFull(startInd:transInd(transID)); % per trial durations including ITI up to transition (in seconds)
            s1TOs3_reward(s1TOs3_COUNTER) = sum(toSumRewardTemp(toSumAccuracyTemp==1)); % total reward (mL) from sess start to transition
            s1TOs3_rewardRate(s1TOs3_COUNTER) = sum(toSumRewardTemp(toSumAccuracyTemp==1))/(sum(toSumTrialDurTemp)/60)*1000; % reward (uL) / sess duration (mL) up to transition
            
            % if transition from state 2-->3
        elseif lg.stateIndex(transInd(transID)) == 2 && lg.stateIndex(transInd(transID)+1) == 3
            startInd = []; startInd = find(diff(lg.sessionID)~=0) ; % get index of all between sess transitions
            startInd = startInd(startInd<transInd(transID)); % select all between sess transitions before current state transition
            startInd = startInd(end) + 1; % the closest of these (+1 trial) is the first trial of session with current state transition
            s2TOs3_COUNTER = s2TOs3_COUNTER+1;
            toSumRewardTemp = []; toSumAccuracyTemp = []; toSumTrialDurTemp = [];
            toSumRewardTemp    = rewardGiven(startInd:transInd(transID)); % scaled rewardGiven per trial up to transition (mL)
            toSumAccuracyTemp  = accuracy(startInd:transInd(transID)); % whether trials up to transition was rewarded
            toSumTrialDurTemp  = lg.trialDurFull(startInd:transInd(transID)); % per trial durations including ITI up to transition (in seconds)
            s2TOs3_reward(s2TOs3_COUNTER) = sum(toSumRewardTemp(toSumAccuracyTemp==1)); % total reward (mL) from sess start to transition
            s2TOs3_rewardRate(s2TOs3_COUNTER) = sum(toSumRewardTemp(toSumAccuracyTemp==1))/(sum(toSumTrialDurTemp)/60)*1000; % reward (uL) / sess duration (mL) up to transition
            
            % if transition from state 2-->1
        elseif lg.stateIndex(transInd(transID)) == 2 && lg.stateIndex(transInd(transID)+1) == 1
            startInd = []; startInd = find(diff(lg.sessionID)~=0) ; % get index of all between sess transitions
            startInd = startInd(startInd<transInd(transID)); % select all between sess transitions before current state transition
            startInd = startInd(end) + 1; % the closest of these (+1 trial) is the first trial of session with current state transition
            s2TOs1_COUNTER = s2TOs1_COUNTER+1;
            toSumRewardTemp    = rewardGiven(startInd:transInd(transID)); % scaled rewardGiven per trial up to transition (mL)
            toSumAccuracyTemp  = accuracy(startInd:transInd(transID)); % whether trials up to transition was rewarded
            toSumTrialDurTemp  = lg.trialDurFull(startInd:transInd(transID)); % per trial durations including ITI up to transition (in seconds)
            s2TOs1_reward(s2TOs1_COUNTER) = sum(toSumRewardTemp(toSumAccuracyTemp==1)); % total reward (mL) from sess start to transition
            s2TOs1_rewardRate(s2TOs1_COUNTER) = sum(toSumRewardTemp(toSumAccuracyTemp==1))/(sum(toSumTrialDurTemp)/60)*1000; % reward (uL) / sess duration (mL) up to transition
            
            % if transition from state 3-->1
        elseif lg.stateIndex(transInd(transID)) == 3 && lg.stateIndex(transInd(transID)+1) == 1
            startInd = []; startInd = find(diff(lg.sessionID)~=0) ; % get index of all between sess transitions
            startInd = startInd(startInd<transInd(transID)); % select all between sess transitions before current state transition
            startInd = startInd(end) + 1; % the closest of these (+1 trial) is the first trial of session with current state transition
            s3TOs1_COUNTER = s3TOs1_COUNTER+1;
            toSumRewardTemp = []; toSumAccuracyTemp = []; toSumTrialDurTemp = [];
            toSumRewardTemp    = rewardGiven(startInd:transInd(transID)); % scaled rewardGiven per trial up to transition (mL)
            toSumAccuracyTemp  = accuracy(startInd:transInd(transID)); % whether trials up to transition was rewarded
            toSumTrialDurTemp  = lg.trialDurFull(startInd:transInd(transID)); % per trial durations including ITI up to transition (in seconds)
            s3TOs1_reward(s3TOs1_COUNTER) = sum(toSumRewardTemp(toSumAccuracyTemp==1)); % total reward (mL) from sess start to transition
            s3TOs1_rewardRate(s3TOs1_COUNTER) = sum(toSumRewardTemp(toSumAccuracyTemp==1))/(sum(toSumTrialDurTemp)/60)*1000; % reward (uL) / sess duration (mL) up to transition
            
            % if transition from state 3-->2
        elseif lg.stateIndex(transInd(transID)) == 3 && lg.stateIndex(transInd(transID)+1) == 2
            startInd = []; startInd = find(diff(lg.sessionID)~=0) ; % get index of all between sess transitions
            startInd = startInd(startInd<transInd(transID)); % select all between sess transitions before current state transition
            startInd = startInd(end) + 1; % the closest of these (+1 trial) is the first trial of session with current state transition
            s3TOs2_COUNTER = s3TOs2_COUNTER+1;
            toSumRewardTemp = []; toSumAccuracyTemp = []; toSumTrialDurTemp = [];
            toSumRewardTemp    = rewardGiven(startInd:transInd(transID)); % scaled rewardGiven per trial up to transition (mL)
            toSumAccuracyTemp  = accuracy(startInd:transInd(transID)); % whether trials up to transition was rewarded
            toSumTrialDurTemp  = lg.trialDurFull(startInd:transInd(transID)); % per trial durations including ITI up to transition (in seconds)
            s3TOs2_reward(s3TOs2_COUNTER) = sum(toSumRewardTemp(toSumAccuracyTemp==1)); % total reward (mL) from sess start to transition
            s3TOs2_rewardRate(s3TOs2_COUNTER) = sum(toSumRewardTemp(toSumAccuracyTemp==1))/(sum(toSumTrialDurTemp)/60)*1000; % reward (uL) / sess duration (mL) up to transition
        else
        end
    end
    
    % save vars
    reward{counter}.s1TOs2     = s1TOs2_reward;
    reward{counter}.s1TOs3     = s1TOs3_reward;
    reward{counter}.s2TOs1     = s2TOs1_reward;
    reward{counter}.s2TOs3     = s2TOs3_reward;
    reward{counter}.s3TOs1     = s3TOs1_reward;
    reward{counter}.s3TOs2     = s3TOs2_reward;
    
    rewardRate{counter}.s1TOs2 = s1TOs2_rewardRate;
    rewardRate{counter}.s1TOs3 = s1TOs3_rewardRate;
    rewardRate{counter}.s2TOs1 = s2TOs1_rewardRate;
    rewardRate{counter}.s2TOs3 = s2TOs3_rewardRate;
    rewardRate{counter}.s3TOs1 = s3TOs1_rewardRate;
    rewardRate{counter}.s3TOs2 = s3TOs2_rewardRate;
    
    % count up for next lg
    counter = counter +1;
    
end

%% clearvars except source data used for plotting
clearvars -except lgDMSD2 lgDMSD1 ...
                  groupOrder reward rewardRate

%% PLOTTING panel 12C,D,E,F
% update these variables to match for
d2orD1     = 2; % 2 for D2 and 1 or D1 (following groupOrder)
dataToplot = reward{d2orD1}; % reward or rewardRate   
yNames     = 'cumulative reward at transition (mL)'; 
%            'cumulative reward at transition (mL)'
%            'reward rate at transition (uL/min)'
yLimits    = [0 1.5];
%            [0 1.5] for reward
%            [0 40]  for rewardRate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if d2orD1 == 1
    titlestr = 'direct';
else
    titlestr = 'indirect';
end

sxTOs1 = []; sxTOs2 = []; sxTOs3 = [];
sxTOs1 = [dataToplot.s2TOs1';dataToplot.s3TOs1'];
sxTOs2 = [dataToplot.s1TOs2';dataToplot.s3TOs2'];
sxTOs3 = [dataToplot.s1TOs3';dataToplot.s2TOs3'];

figure
data = []; xIdx = []; 
dataArray = []; colorIdx = [];
dataArray = {sxTOs1,sxTOs2,sxTOs3};
data      = [sxTOs1; sxTOs2; sxTOs3];
xIdx      = [ones(numel(sxTOs1),1); ones(numel(sxTOs2),1)*2; ones(numel(sxTOs3),1)*3 ];
catMarker = {'x','x', 'x'};
colorIdx  = {[0 0.5 0.8]  [0.9 0.7 0.2] [0.9 0 0]};
hold on
plotSpread(data,'categoryIdx',xIdx,'distributionIdx',xIdx,'categoryColors',colorIdx,'categoryMarkers',catMarker);
h = errorbar(1,nanmean(sxTOs1),nanstd(sxTOs1)./sqrt(numel(sxTOs1)-1), ...
    '.-','linewidth',4,'color','b','markersize',0.1);  h.CapSize = 40;
h = errorbar(2,nanmean(sxTOs2),nanstd(sxTOs2)./sqrt(numel(sxTOs2)-1), ...
    '.-','linewidth',4,'color','y','markersize',0.1);  h.CapSize = 40;
h = errorbar(3,nanmean(sxTOs3),nanstd(sxTOs3)./sqrt(numel(sxTOs3)-1), ...
    '.-','linewidth',4,'color','r','markersize',0.1);  h.CapSize = 40;

set(gca,'xtick',1:3,'xticklabel',{'s1';'s2';'s3'})
ylabel(yNames,'fontsize',12)
xlim([.5 3.5])
ylim(yLimits)
title(titlestr)
box off
xIdx2 = categorical(xIdx);
[p,t,stats] = anova1(data,xIdx2,'off');

[s1_s2_p,~,s1_s2_stats] = ranksum(data(xIdx==1),data(xIdx==2),'method','approximate');
[s1_s3_p,~,s1_s3_stats] = ranksum(data(xIdx==1),data(xIdx==3),'method','approximate');
[s2_s3_p,~,s2_s3_stats] = ranksum(data(xIdx==2),data(xIdx==3),'method','approximate');

text(1.5,yLimits(2)*.75,['state p=' num2str(p)],'Color', 'k')
text(1.5,yLimits(2)*.70,['f=' num2str(cell2mat(t(2,5)))],'Color', 'k')
text(1.5,yLimits(2)*.65,['df=' num2str(stats.df)],'Color', 'k')

plot(1:0.5:2,ones(1,3)*yLimits(2)*.9,'k-')
text(1,yLimits(2)*.92,['p=' num2str(round(s1_s2_p,3))],'FontSize',8)
text(1.5,yLimits(2)*.92,['z=' num2str(round(s1_s2_stats.zval,3))],'FontSize',8)
plot(2:0.5:3,ones(1,3)*yLimits(2)*.85,'k-')
text(2,yLimits(2)*.87,['p=' num2str(round(s2_s3_p,3))],'FontSize',8)
text(2.5,yLimits(2)*.87,['z=' num2str(round(s2_s3_stats.zval,3))],'FontSize',8)
plot(1:0.5:3,ones(1,5)*yLimits(2)*.95,'k-')
text(1.2,yLimits(2)*.97,['p=' num2str(round(s1_s3_p,3))],'FontSize',8)
text(1.8,yLimits(2)*.97,['z=' num2str(round(s1_s3_stats.zval,3))],'FontSize',8)

%% clearvars except source data used for plotting
clearvars -except lgDMSD2 lgDMSD1 ...
                  groupOrder reward rewardRate

%% PLOTTING panel 12C,D,E,F
% update these variables to match for
d2orD1     = 2; % 2 for D2 and 1 or D1 (following groupOrder)
dataToplot = rewardRate{d2orD1}; % reward or rewardRate   
yNames     = 'reward rate at transition (uL/min)'; 
%            'cumulative reward at transition (mL)'
%            'reward rate at transition (uL/min)'
yLimits    = [0 40];
%            [0 1.5] for reward
%            [0 40]  for rewardRate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if d2orD1 == 1
    titlestr = 'direct';
else
    titlestr = 'indirect';
end

sxTOs1 = []; sxTOs2 = []; sxTOs3 = [];
sxTOs1 = [dataToplot.s2TOs1';dataToplot.s3TOs1'];
sxTOs2 = [dataToplot.s1TOs2';dataToplot.s3TOs2'];
sxTOs3 = [dataToplot.s1TOs3';dataToplot.s2TOs3'];

figure
data = []; xIdx = []; 
dataArray = []; colorIdx = [];
dataArray = {sxTOs1,sxTOs2,sxTOs3};
data      = [sxTOs1; sxTOs2; sxTOs3];
xIdx      = [ones(numel(sxTOs1),1); ones(numel(sxTOs2),1)*2; ones(numel(sxTOs3),1)*3 ];
catMarker = {'x','x', 'x'};
colorIdx  = {[0 0.5 0.8]  [0.9 0.7 0.2] [0.9 0 0]};
hold on
plotSpread(data,'categoryIdx',xIdx,'distributionIdx',xIdx,'categoryColors',colorIdx,'categoryMarkers',catMarker);
h = errorbar(1,nanmean(sxTOs1),nanstd(sxTOs1)./sqrt(numel(sxTOs1)-1), ...
    '.-','linewidth',4,'color','b','markersize',0.1);  h.CapSize = 40;
h = errorbar(2,nanmean(sxTOs2),nanstd(sxTOs2)./sqrt(numel(sxTOs2)-1), ...
    '.-','linewidth',4,'color','y','markersize',0.1);  h.CapSize = 40;
h = errorbar(3,nanmean(sxTOs3),nanstd(sxTOs3)./sqrt(numel(sxTOs3)-1), ...
    '.-','linewidth',4,'color','r','markersize',0.1);  h.CapSize = 40;

set(gca,'xtick',1:3,'xticklabel',{'s1';'s2';'s3'})
ylabel(yNames,'fontsize',12)
xlim([.5 3.5])
ylim(yLimits)
title(titlestr)
box off
xIdx2 = categorical(xIdx);
[p,t,stats] = anova1(data,xIdx2,'off');

[s1_s2_p,~,s1_s2_stats] = ranksum(data(xIdx==1),data(xIdx==2),'method','approximate');
[s1_s3_p,~,s1_s3_stats] = ranksum(data(xIdx==1),data(xIdx==3),'method','approximate');
[s2_s3_p,~,s2_s3_stats] = ranksum(data(xIdx==2),data(xIdx==3),'method','approximate');

text(1.5,yLimits(2)*.75,['state p=' num2str(p)],'Color', 'k')
text(1.5,yLimits(2)*.70,['f=' num2str(cell2mat(t(2,5)))],'Color', 'k')
text(1.5,yLimits(2)*.65,['df=' num2str(stats.df)],'Color', 'k')

plot(1:0.5:2,ones(1,3)*yLimits(2)*.9,'k-')
text(1,yLimits(2)*.92,['p=' num2str(round(s1_s2_p,3))],'FontSize',8)
text(1.5,yLimits(2)*.92,['z=' num2str(round(s1_s2_stats.zval,3))],'FontSize',8)
plot(2:0.5:3,ones(1,3)*yLimits(2)*.85,'k-')
text(2,yLimits(2)*.87,['p=' num2str(round(s2_s3_p,3))],'FontSize',8)
text(2.5,yLimits(2)*.87,['z=' num2str(round(s2_s3_stats.zval,3))],'FontSize',8)
plot(1:0.5:3,ones(1,5)*yLimits(2)*.95,'k-')
text(1.2,yLimits(2)*.97,['p=' num2str(round(s1_s3_p,3))],'FontSize',8)
text(1.8,yLimits(2)*.97,['z=' num2str(round(s1_s3_stats.zval,3))],'FontSize',8)

%% clearvars except source data used for plotting
clearvars -except lgDMSD2 lgDMSD1 ...
                  groupOrder reward rewardRate

%% PLOTTING panel 12C,D,E,F
% update these variables to match for
d2orD1     = 1; % 2 for D2 and 1 or D1 (following groupOrder)
dataToplot = reward{d2orD1}; % reward or rewardRate   
yNames     = 'cumulative reward at transition (mL)'; 
%            'cumulative reward at transition (mL)'
%            'reward rate at transition (uL/min)'
yLimits    = [0 1.5];
%            [0 1.5] for reward
%            [0 40]  for rewardRate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if d2orD1 == 1
    titlestr = 'direct';
else
    titlestr = 'indirect';
end

sxTOs1 = []; sxTOs2 = []; sxTOs3 = [];
sxTOs1 = [dataToplot.s2TOs1';dataToplot.s3TOs1'];
sxTOs2 = [dataToplot.s1TOs2';dataToplot.s3TOs2'];
sxTOs3 = [dataToplot.s1TOs3';dataToplot.s2TOs3'];

figure
data = []; xIdx = []; 
dataArray = []; colorIdx = [];
dataArray = {sxTOs1,sxTOs2,sxTOs3};
data      = [sxTOs1; sxTOs2; sxTOs3];
xIdx      = [ones(numel(sxTOs1),1); ones(numel(sxTOs2),1)*2; ones(numel(sxTOs3),1)*3 ];
catMarker = {'x','x', 'x'};
colorIdx  = {[0 0.5 0.8]  [0.9 0.7 0.2] [0.9 0 0]};
hold on
plotSpread(data,'categoryIdx',xIdx,'distributionIdx',xIdx,'categoryColors',colorIdx,'categoryMarkers',catMarker);
h = errorbar(1,nanmean(sxTOs1),nanstd(sxTOs1)./sqrt(numel(sxTOs1)-1), ...
    '.-','linewidth',4,'color','b','markersize',0.1);  h.CapSize = 40;
h = errorbar(2,nanmean(sxTOs2),nanstd(sxTOs2)./sqrt(numel(sxTOs2)-1), ...
    '.-','linewidth',4,'color','y','markersize',0.1);  h.CapSize = 40;
h = errorbar(3,nanmean(sxTOs3),nanstd(sxTOs3)./sqrt(numel(sxTOs3)-1), ...
    '.-','linewidth',4,'color','r','markersize',0.1);  h.CapSize = 40;

set(gca,'xtick',1:3,'xticklabel',{'s1';'s2';'s3'})
ylabel(yNames,'fontsize',12)
xlim([.5 3.5])
ylim(yLimits)
title(titlestr)
box off
xIdx2 = categorical(xIdx);
[p,t,stats] = anova1(data,xIdx2,'off');

[s1_s2_p,~,s1_s2_stats] = ranksum(data(xIdx==1),data(xIdx==2),'method','approximate');
[s1_s3_p,~,s1_s3_stats] = ranksum(data(xIdx==1),data(xIdx==3),'method','approximate');
[s2_s3_p,~,s2_s3_stats] = ranksum(data(xIdx==2),data(xIdx==3),'method','approximate');

text(1.5,yLimits(2)*.75,['state p=' num2str(p)],'Color', 'k')
text(1.5,yLimits(2)*.70,['f=' num2str(cell2mat(t(2,5)))],'Color', 'k')
text(1.5,yLimits(2)*.65,['df=' num2str(stats.df)],'Color', 'k')

plot(1:0.5:2,ones(1,3)*yLimits(2)*.9,'k-')
text(1,yLimits(2)*.92,['p=' num2str(round(s1_s2_p,3))],'FontSize',8)
text(1.5,yLimits(2)*.92,['z=' num2str(round(s1_s2_stats.zval,3))],'FontSize',8)
plot(2:0.5:3,ones(1,3)*yLimits(2)*.85,'k-')
text(2,yLimits(2)*.87,['p=' num2str(round(s2_s3_p,3))],'FontSize',8)
text(2.5,yLimits(2)*.87,['z=' num2str(round(s2_s3_stats.zval,3))],'FontSize',8)
plot(1:0.5:3,ones(1,5)*yLimits(2)*.95,'k-')
text(1.2,yLimits(2)*.97,['p=' num2str(round(s1_s3_p,3))],'FontSize',8)
text(1.8,yLimits(2)*.97,['z=' num2str(round(s1_s3_stats.zval,3))],'FontSize',8)

%% clearvars except source data used for plotting
clearvars -except lgDMSD2 lgDMSD1 ...
                  groupOrder reward rewardRate

%% PLOTTING panel 12C,D,E,F
% update these variables to match for
d2orD1     = 1; % 2 for D2 and 1 or D1 (following groupOrder)
dataToplot = rewardRate{d2orD1}; % reward or rewardRate   
yNames     = 'reward rate at transition (uL/min)'; 
%            'cumulative reward at transition (mL)'
%            'reward rate at transition (uL/min)'
yLimits    = [0 40];
%            [0 1.5] for reward
%            [0 40]  for rewardRate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if d2orD1 == 1
    titlestr = 'direct';
else
    titlestr = 'indirect';
end

sxTOs1 = []; sxTOs2 = []; sxTOs3 = [];
sxTOs1 = [dataToplot.s2TOs1';dataToplot.s3TOs1'];
sxTOs2 = [dataToplot.s1TOs2';dataToplot.s3TOs2'];
sxTOs3 = [dataToplot.s1TOs3';dataToplot.s2TOs3'];

figure
data = []; xIdx = []; 
dataArray = []; colorIdx = [];
dataArray = {sxTOs1,sxTOs2,sxTOs3};
data      = [sxTOs1; sxTOs2; sxTOs3];
xIdx      = [ones(numel(sxTOs1),1); ones(numel(sxTOs2),1)*2; ones(numel(sxTOs3),1)*3 ];
catMarker = {'x','x', 'x'};
colorIdx  = {[0 0.5 0.8]  [0.9 0.7 0.2] [0.9 0 0]};
hold on
plotSpread(data,'categoryIdx',xIdx,'distributionIdx',xIdx,'categoryColors',colorIdx,'categoryMarkers',catMarker);
h = errorbar(1,nanmean(sxTOs1),nanstd(sxTOs1)./sqrt(numel(sxTOs1)-1), ...
    '.-','linewidth',4,'color','b','markersize',0.1);  h.CapSize = 40;
h = errorbar(2,nanmean(sxTOs2),nanstd(sxTOs2)./sqrt(numel(sxTOs2)-1), ...
    '.-','linewidth',4,'color','y','markersize',0.1);  h.CapSize = 40;
h = errorbar(3,nanmean(sxTOs3),nanstd(sxTOs3)./sqrt(numel(sxTOs3)-1), ...
    '.-','linewidth',4,'color','r','markersize',0.1);  h.CapSize = 40;

set(gca,'xtick',1:3,'xticklabel',{'s1';'s2';'s3'})
ylabel(yNames,'fontsize',12)
xlim([.5 3.5])
ylim(yLimits)
title(titlestr)
box off
xIdx2 = categorical(xIdx);
[p,t,stats] = anova1(data,xIdx2,'off');

[s1_s2_p,~,s1_s2_stats] = ranksum(data(xIdx==1),data(xIdx==2),'method','approximate');
[s1_s3_p,~,s1_s3_stats] = ranksum(data(xIdx==1),data(xIdx==3),'method','approximate');
[s2_s3_p,~,s2_s3_stats] = ranksum(data(xIdx==2),data(xIdx==3),'method','approximate');

text(1.5,yLimits(2)*.75,['state p=' num2str(p)],'Color', 'k')
text(1.5,yLimits(2)*.70,['f=' num2str(cell2mat(t(2,5)))],'Color', 'k')
text(1.5,yLimits(2)*.65,['df=' num2str(stats.df)],'Color', 'k')

plot(1:0.5:2,ones(1,3)*yLimits(2)*.9,'k-')
text(1,yLimits(2)*.92,['p=' num2str(round(s1_s2_p,3))],'FontSize',8)
text(1.5,yLimits(2)*.92,['z=' num2str(round(s1_s2_stats.zval,3))],'FontSize',8)
plot(2:0.5:3,ones(1,3)*yLimits(2)*.85,'k-')
text(2,yLimits(2)*.87,['p=' num2str(round(s2_s3_p,3))],'FontSize',8)
text(2.5,yLimits(2)*.87,['z=' num2str(round(s2_s3_stats.zval,3))],'FontSize',8)
plot(1:0.5:3,ones(1,5)*yLimits(2)*.95,'k-')
text(1.2,yLimits(2)*.97,['p=' num2str(round(s1_s3_p,3))],'FontSize',8)
text(1.8,yLimits(2)*.97,['z=' num2str(round(s1_s3_stats.zval,3))],'FontSize',8)

%% clearvars except source data used for plotting
clearvars -except lgDMSD2 lgDMSD1 ...
                  groupOrder reward rewardRate