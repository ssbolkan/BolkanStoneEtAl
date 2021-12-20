function [bias, dataSumm] = xMousePerfBias(lg,ipsiCon, cleanLog, relaxTravel)

if nargin < 2 || isempty(ipsiCon)
    ipsiCon = 0;
end

if nargin < 3 || isempty(cleanLog)
    cleanLog = false; % of true removes trial blocks >0.6 performance
end

if nargin < 4 || isempty(relaxTravel)
    relaxTravel = false;  % if true then allows excess travel trials >0.1
end

if ipsiCon
    lg = invertLogs(lg);
else
end

dataSumm.nTrialsPreSelection = numel(lg.choice);
dataSumm.nMicePreSelection   = numel(unique(lg.mouseID));

if cleanLog == 1   
    [lg, ~]   = selectLogSubset(lg,[],[],[],[],relaxTravel,0.6); 
    % lg, minNT, metaMouse, minHistory, selectTrialsWithHistory, relaxCriteria, perfTh)
    dataSumm.nTrialsPostSelection = numel(lg.choice);
    dataSumm.nMicePostSelection = numel(unique(lg.mouseID));
else
    dataSumm.nTrialsPostSelection = numel(lg.choice);
    dataSumm.nMicePostSelection = numel(unique(lg.mouseID));
end

stdInterval = normcdf(1, 0, 1) - normcdf(-1, 0, 1);

mouseName = unique(lg.mouseID);
for nMouse = 1:numel(mouseName)
        
  [bias(nMouse).perfOFF,bias(nMouse).perfOFF_binoErr(:)] = ...
    binofit(sum(lg.trialType(~lg.laserON & lg.mouseID == mouseName(nMouse)) == ...
    lg.choice(~lg.laserON & lg.mouseID == mouseName(nMouse))),...
    sum(~lg.laserON & lg.mouseID == mouseName(nMouse)),1-stdInterval);

  [bias(nMouse).perfON,bias(nMouse).perfON_binoErr(:)] = ...
    binofit(sum(lg.trialType(lg.laserON & lg.mouseID == mouseName(nMouse) ) == ...
    lg.choice(lg.laserON & lg.mouseID == mouseName(nMouse) )),...
    sum(lg.laserON & lg.mouseID == mouseName(nMouse)),1-stdInterval);   
  
    bias(nMouse).nTrialsOFF = numel(lg.choice(~lg.laserON & lg.mouseID == mouseName(nMouse)));
    bias(nMouse).nTrialsON  = numel(lg.choice(lg.laserON & lg.mouseID == mouseName(nMouse)));
    
  [bias(nMouse).perfOFF_contra,bias(nMouse).perfOFF_binoErr_contra(:)] = ...
    binofit(sum(lg.trialType(~lg.laserON & lg.mouseID == mouseName(nMouse) & lg.trialType == globalParams.rightCode) == ...
    lg.choice(~lg.laserON & lg.mouseID == mouseName(nMouse) & lg.trialType == globalParams.rightCode)),...
    sum(~lg.laserON & lg.mouseID == mouseName(nMouse) & lg.trialType == globalParams.rightCode),1-stdInterval);

  [bias(nMouse).perfON_contra,bias(nMouse).perfON_binoErr_contra(:)] = ...
    binofit(sum(lg.trialType(lg.laserON & lg.mouseID == mouseName(nMouse) & lg.trialType == globalParams.rightCode) == ...
    lg.choice(lg.laserON & lg.mouseID == mouseName(nMouse) & lg.trialType == globalParams.rightCode)),...
    sum(lg.laserON & lg.mouseID == mouseName(nMouse) & lg.trialType == globalParams.rightCode),1-stdInterval);
  
  [bias(nMouse).perfOFF_ipsi,bias(nMouse).perfOFF_binoErr_ipsi(:)] = ...
    binofit(sum(lg.trialType(~lg.laserON & lg.mouseID == mouseName(nMouse) & lg.trialType == globalParams.leftCode) == ...
    lg.choice(~lg.laserON & lg.mouseID == mouseName(nMouse) & lg.trialType == globalParams.leftCode)),...
    sum(~lg.laserON & lg.mouseID == mouseName(nMouse) & lg.trialType == globalParams.leftCode),1-stdInterval);

  [bias(nMouse).perfON_ipsi,bias(nMouse).perfON_binoErr_ipsi(:)] = ...
    binofit(sum(lg.trialType(lg.laserON & lg.mouseID == mouseName(nMouse) & lg.trialType == globalParams.leftCode) == ...
    lg.choice(lg.laserON & lg.mouseID == mouseName(nMouse) & lg.trialType == globalParams.leftCode)),...
    sum(lg.laserON & lg.mouseID == mouseName(nMouse) & lg.trialType == globalParams.leftCode),1-stdInterval);
  
	bias(nMouse).contraBiasOFF    = bias(nMouse).perfOFF_contra - bias(nMouse).perfOFF_ipsi;
	bias(nMouse).contraBiasON     = bias(nMouse).perfON_contra - bias(nMouse).perfON_ipsi; 
    
    bias(nMouse).contraBiasDiff   = bias(nMouse).contraBiasON - bias(nMouse).contraBiasOFF;   
    bias(nMouse).perfDiff         = bias(nMouse).perfON - bias(nMouse).perfOFF;
end