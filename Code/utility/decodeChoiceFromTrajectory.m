function [predChoiceVA] = decodeChoiceFromTrajectory(lg,blockAvg,minMaxPos,numFolds,numSamples,fitType,cleanLog,weighted,zPred,ipCon,optimalThresh)

if nargin < 2 || isempty(blockAvg)
    blockAvg = 25; %  [25 1]
end

if nargin < 3 || isempty(minMaxPos)
    minMaxPos = [1 300]; % [1 300]
end

if nargin < 4  || isempty(numFolds)
    numFolds = 5;
end

if nargin < 5  || isempty(numSamples)
    numSamples = 10;
end

if nargin < 6 || isempty(fitType)
    fitType = 'fitglm'; % previous opts ridgeSeperate ridgeTogether removed from this
end

if nargin < 7 || isempty(cleanLog)
    cleanLog = true; % removes below 0.6 perf and 0.1 excess travel
end

if nargin < 8 || isempty(weighted)
    weighted = true; % option within fitglmRegression function
end

if nargin < 9 || isempty(zPred)
    zPred = false; % zScore predMat or not
end

if nargin < 10 || isempty(ipCon)
    ipCon = false; % convert to ipsicon or not
end

if nargin < 11 || isempty(optimalThresh)
    optimalThresh = false; % convert to ipsicon or not
end

predChoiceVA.blockAvg      = blockAvg;
predChoiceVA.maxPos        = minMaxPos;
predChoiceVA.numFolds      = numFolds;
predChoiceVA.numSamples    = numSamples;
predChoiceVA.fitType       = fitType;
predChoiceVA.weighted      = weighted;
predChoiceVA.zPred         = zPred;
predChoiceVA.ipCon         = ipCon;
predChoiceVA.optimalThresh = optimalThresh;


%% must remove NAN and -1 choice trials in log fields for later fit
removeMouse = [];
removeMouse = isnan(lg.choice) | lg.choice==-1;
if isfield(lg,'keyFrameLabels')
    keyFrameLabelsTemp = []; keyFrameLabelsTemp = lg.keyFrameLabels;
    lg = rmfield(lg,'keyFrameLabels');
end
lg  = structfun(@(x) x(~removeMouse),lg,'UniformOutput',false);

%% remove perf below 0.6 and excess travel > 0.1
if cleanLog == 1
    relaxedCriteria = false;
    [lg, ~]   = selectLogSubset(lg,[],[],[],[],relaxedCriteria,0);  % lg, minNT, metaMouse, minHistory, selectTrialsWithHistory, relaxCriteria, perfTh)
else
end

%% invert to ipsi contra laser hemisphere
if ipCon == 1
    lg = invertLogs(lg);
else
end
%% get VA from each position 0:1cm:300cm and then average in increments of input blockAvg

% viewAngle     = sampleViewAngleVsY(lg.pos, 0:1:330); % sampleViewAngleVsY(lg.pos, posBins,1);
viewAngle     = sampleViewAngleVsY(lg.pos, minMaxPos(1):1:minMaxPos(end)); % sampleViewAngleVsY(lg.pos, posBins,1);
sz            = size(viewAngle);
binAvg        = blockAvg(1);
n             = size(viewAngle, 1);             % Length of first dimension
nc            = n - mod(n, binAvg);             % Multiple of binAvg
np            = nc / binAvg;                    % Length of result / number of trials
xx            = reshape(viewAngle(1:nc, :), binAvg, np, []); % [binAvg x np(trials) x size(viewAngles,2)]
avgViewAngle  = sum(xx, 1) / binAvg;            % Mean over 1st dim
avgViewAngle  = reshape(avgViewAngle, np, []);  % Remove leading dim of length 1


% used for posLogic to select maze region
% binIDs               = 0:blockAvg:329;
binIDs               = minMaxPos(1):blockAvg:minMaxPos(end);
% binIDs               = (binIDs - 1) + blockAvg/2;
predChoiceVA.binIDs  = binIDs;


%% loop through mice
mice             = unique(lg.mouseID);
for iMouse = 1:numel(mice)
    %% clear and initialize variables used in mouse loops
    laserON               = [];    choice                = [];    sel                   = [];
    viewAng               = [];    viewAngOFF            = [];    viewAngON             = [];
    choiceOFF             = [];    choiceON              = [];    newLog                = [];
    newLogOFF             = [];    newLogON              = [];    viewAngOFF_R          = [];
    viewAngOFF_L          = [];    viewAngON_R           = [];    viewAngON_L           = [];
    nCues_RminusL         = [];    pos                   = [];    nCues_RminusL_OFF     = [];
    nCues_RminusL_ON      = [];    posOFF                = [];    posON                 = [];
    prevChoice            = [];    nextChoice            = [];    date                  = [];
    dateOFF               = [];    dateON                = [];
    
    %% extract relevant log variables for fitting
    newLog                = selectMouseTrials(lg, mice(iMouse));
    [choice,sel]          = selectMouseTrials(lg, mice(iMouse), 'choice');
    prevChoice            = [0 choice(1:end-1)];
    nextChoice            = [choice(2:end) 0];
    nCues_RminusL         = selectMouseTrials(lg, mice(iMouse), 'nCues_RminusL');
    pos                   = selectMouseTrials(lg, mice(iMouse), 'pos');
    laserON               = selectMouseTrials(lg, mice(iMouse), 'laserON');
    date                  = selectMouseTrials(lg, mice(iMouse), 'date');
    dateOFF               = date(laserON==0);
    dateON                = date(laserON==1);
    viewAng               = avgViewAngle(:,sel);
    viewAng               = viewAng + (rand(size(viewAng)) - 0.5)*1e-4;     % prevent decoding degeneracies
    viewAngOFF            = viewAng(:,laserON==0);
    viewAngON             = viewAng(:,laserON==1);
    choiceOFF             = choice(laserON==0);
    choiceON              = choice(laserON==1);
    viewAngOFF_R          = viewAng(:,laserON==0 & choice == 1);
    viewAngOFF_L          = viewAng(:,laserON==0 & choice == 0);
    viewAngON_R           = viewAng(:,laserON==1 & choice == 1);
    viewAngON_L           = viewAng(:,laserON==1 & choice == 0);
    posOFF                = pos(laserON==1);
    posON                 = pos(laserON==0);
    nCues_RminusL_OFF     = nCues_RminusL(laserON==0);
    nCues_RminusL_ON      = nCues_RminusL(laserON==1);
    
    
    %% extract and save light off and light on log files in case useful for plotting
    try keyFrameLabels = []; keyFrameLabels = newLog.keyFrameLabels; catch; end
    try newLog = rmfield(newLog,'keyFrameLabels'); catch; end
    removeMouse = []; removeMouse = newLog.laserON==0; % remove light on
    newLogOFF   = structfun(@(x) x(removeMouse),newLog,'UniformOutput',false);
    
    removeMouse = []; removeMouse = newLog.laserON==1; % remove light on
    newLogON    = structfun(@(x) x(removeMouse),newLog,'UniformOutput',false);
    
    predChoiceVA.mouse(iMouse).newLog         = newLog;
    predChoiceVA.mouse(iMouse).newLogOFF      = newLogOFF;
    predChoiceVA.mouse(iMouse).newLogON       = newLogON;
    predChoiceVA.mouse(iMouse).viewAngOFF     = viewAngOFF;
    predChoiceVA.mouse(iMouse).viewAngON      = viewAngON;
    predChoiceVA.mouse(iMouse).viewAngOFF_R   = viewAngOFF_R;
    predChoiceVA.mouse(iMouse).viewAngOFF_L   = viewAngOFF_L;
    predChoiceVA.mouse(iMouse).viewAngON_R    = viewAngON_R;
    predChoiceVA.mouse(iMouse).viewAngON_L    = viewAngON_L;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% early versions include ridge regression but dropped from here
    
    switch fitType                
        case 'fitglm'
            %% define predictors and outputs and save for quick reference
            rawPredMatOFF            = double(viewAngOFF)';
            output                   = double(choiceOFF)';
            predChoiceVA.mouse(iMouse).rawPredMatOFF ...
                =  rawPredMatOFF;
            predChoiceVA.mouse(iMouse).choiceOFF ...
                =  output;
            
            rawPredMatON             = double(viewAngON)';
            predChoiceVA.mouse(iMouse).rawPredMatON ...
                =  rawPredMatON;
            predChoiceVA.mouse(iMouse).choiceON ...
                =  double(choiceON)';
            
            %% zscore prediction matrices and save if opt selected
            if zPred == 1
                predMatOFF           = zscore(rawPredMatOFF,0,2);
                predMatON            = zscore(rawPredMatON,0,2);
            elseif zPred == 0
                predMatOFF           = rawPredMatOFF;
                predMatON            = rawPredMatON;
            end
            
            predChoiceVA.mouse(iMouse).predMatOFF ...
                =  predMatOFF;
            predChoiceVA.mouse(iMouse).predMatON ...
                =  predMatON;
            
            %% initialize fit variables
            bestCoeff       = [];     accuracyOFF     = [];    accIntervalOFF  = [];
            accuracyON      = [];     accIntervalON   = [];    coeff           = [];
            fitInfo         = [];     allCoeffs       = [];    thresh          = [];
            pseudoExp       = [];     yON             = [];
            yON             = choiceON';
            
            %% set up cross-validated index for off data
            warning('off', 'stats:cvpartition:KFoldMissingGrp');
            randGenerator           = RandStream('mt19937ar', 'Seed', 723176);
            pseudoExp               = cvpartition(output(:,end), 'KFold', numFolds, randGenerator);
            for iMC = 2:numSamples
                pseudoExp(iMC)      = repartition(pseudoExp(1));
            end
            warning('on', 'stats:cvpartition:KFoldMissingGrp');
                        
            for nbin = 1:numel(binIDs)
                %% do fit to all bins train/test light off data according to pseudoExp
                % note can add lasso here instead
                coeff = []; fitInfo = [];
                [coeff, fitInfo] = fitglmRegression(predMatOFF(:,nbin), output, pseudoExp,numSamples,weighted);
                                
                %% appply train set coeff according to pseudoExp to predict light on data
                for iMC = 1:numel(pseudoExp)
                    for iFold = 1:numFolds
                        % Select train/test sets with user-specified pruning if so desired
                        if isfield(pseudoExp(iMC).training(iFold),'idx')
                            iTrain              = pseudoExp(iMC).training(iFold).idx;
                        else
                            iTrain              = pseudoExp(iMC).training(iFold);
                        end
                        
                        % coeff fits
                        mcCoeff               = fitInfo.mcCoeff{iMC,iFold};                        
                        % test light off prediction
                        fitInfo.CVTestONPrediction{iMC,iFold} = glmval(mcCoeff,predMatON(:,nbin), 'logit');                        
                    end
                end
                
                %% Compute accuracy as the fraction of correct predictions, with spreads across CV samples
                % Cross-validated prediction accuracy
                for iMC = 1:size(fitInfo.CVTestPrediction,1)
                    for iFold = 1:size(fitInfo.CVTestPrediction,2)
                        if optimalThresh == 1
                            % Compute optimal threshold for separating training set populations according to truth
                            trainSel          = fitInfo.CVExperiments(iMC).training(iFold);
                            predictionTrain   = fitInfo.CVTrainPrediction{iMC,iFold};
                            try
                                threshold     = optimalClassThreshold(predictionTrain, output(trainSel)>0.5);
                            catch
                                threshold     = getThreshold(predictionTrain);
                            end
                            thresh(iMC,iFold) = threshold;
                        else
                            threshold = 0.5;
                        end
                        thresh(iMC,iFold)                        = threshold;
                        % Compute accuracy by applying threshold to prediction in test set and comparing to truth
                        testSel                                  = fitInfo.CVExperiments(iMC).test(iFold);
                        predictionOFF                            = fitInfo.CVTestPrediction{iMC,iFold};
                        predTruthOFF                             = output(testSel);
                        predictionOFF(predictionOFF < threshold) = 0;
                        predictionOFF(predictionOFF > threshold) = 1;
                        predAccuracyOFF(iMC,iFold)               = sum(predictionOFF == predTruthOFF)./numel(predTruthOFF);
                        
                        predictionON                             = fitInfo.CVTestONPrediction{iMC,iFold};
                        predTruthON                              = yON;
                        predictionON(predictionON < threshold)   = 0;
                        predictionON(predictionON > threshold)   = 1;
                        predAccuracyON(iMC,iFold)                = sum(predictionON == predTruthON)./numel(predTruthON);
                        
                    end
                end
                
                cl             = normcdf([0 1], 0, 1); % confidence level
                accuracyOFF    = nanmean(predAccuracyOFF(:));
                accIntervalOFF = quantile(predAccuracyOFF(:), cl);
                
                accuracyON     = nanmean(predAccuracyON(:));
                accIntervalON  = quantile(predAccuracyON(:), cl);
                
                %% save relevant variables
                predChoiceVA.mouse(iMouse).fitInfo{nbin}            = fitInfo;
                predChoiceVA.mouse(iMouse).bestCoeff(nbin,:)        = coeff';
                %                 predChoiceVA.mouse(iMouse).allCoeffs{nbin}          = allCoeffs;
                predChoiceVA.mouse(iMouse).accuracyOFF(nbin)        = accuracyOFF;
                predChoiceVA.mouse(iMouse).accIntervalOFF(nbin,:)   = accIntervalOFF';
                predChoiceVA.mouse(iMouse).accuracyON(nbin)         = accuracyON;
                predChoiceVA.mouse(iMouse).accIntervalON(nbin,:)    = accIntervalON';
                predChoiceVA.mouse(iMouse).threshold{nbin}          = thresh;
                predChoiceVA.mouse(iMouse).fitInfo{nbin}.pseudoExp  = pseudoExp;
                          
            end

    end
    
end

end