function [predChoiceVA] = decodeChoiceFromTraject_xTasksAndMice(lg,blockAvg,minMaxPos,numFolds,numSamples,fitType,cleanLog,weighted,zPred,ipCon,optimalThresh)

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
    fitType = 'fitglm'; % opts ridgeSeperate ridgeTogether fitglm
end

if nargin < 7 || isempty(cleanLog)
    cleanLog = true; % removes below 0.6 perf and 0.1 excess travel
end

if nargin < 8 || isempty(weighted)
    weighted = true; % currentlly not supported for ridge
end

if nargin < 9 || isempty(zPred)
    zPred = false; % zScore predMat or not
end

if nargin < 10 || isempty(ipCon)
    ipCon = false; % convert to ipsicon or not
end

if nargin < 11 || isempty(optimalThresh)
    optimalThresh = false; %
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

%%
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

binIDs               = minMaxPos(1):blockAvg:minMaxPos(end);
predChoiceVA.binIDs  = binIDs;

%% clear and initialize variables used in mouse loops
choice                = [];    sel                   = [];    taskID                = [];
viewAng               = [];    viewAng_R             = [];    viewAng_L             = [];
viewAngTOW            = [];    viewAngTOW_R          = [];    viewAngTOW_L          = [];
viewAngMG             = [];    viewAngMG_R           = [];    viewAngMG_L           = [];
viewAngPC             = [];    viewAngPC_R           = [];    viewAngPC_L           = [];
choiceTOW             = [];    choiceMG              = [];    choicePC              = [];
mouseID               = [];    mouseIDTOW            = [];    mouseIDMG             = [];
mouseIDPC             = [];

%% extract relevant log variables for fitting
choice                = lg.choice;
taskID                = lg.taskID;
mouseID               = lg.mouseID;
viewAng               = avgViewAngle;
viewAng               = viewAng + (rand(size(viewAng)) - 0.5)*1e-4;     % prevent decoding degeneracies
viewAng_R             = viewAng(:,choice == 1);
viewAng_L             = viewAng(:,choice == 0);

viewAngTOW            = viewAng(:,ismember(taskID, 'aoe'));
viewAngTOW_R          = viewAng(:,ismember(taskID, 'aoe') & choice == 0);
viewAngTOW_L          = viewAng(:,ismember(taskID, 'aoe') & choice == 1);
choiceTOW             = choice(ismember(taskID,'aoe'));
mouseIDTOW            = mouseID(ismember(taskID,'aoe'));

viewAngMG             = viewAng(:,ismember(taskID,'ct1'));
viewAngMG_R           = viewAng(:,ismember(taskID,'ct1') & choice == 0);
viewAngMG_L           = viewAng(:,ismember(taskID,'ct1') & choice == 1);
choiceMG              = choice(ismember(taskID,'ct1'));
mouseIDMG             = mouseID(ismember(taskID,'ct1'));

viewAngPC             = viewAng(:,ismember(taskID,'ct2'));
viewAngPC_R           = viewAng(:,ismember(taskID,'ct2') & choice == 0);
viewAngPC_L           = viewAng(:,ismember(taskID,'ct2') & choice == 1);
choicePC              = choice(ismember(taskID,'ct2'));
mouseIDPC             = mouseID(ismember(taskID,'ct2'));


%% extract and save light off and light on log files in case useful for plotting

predChoiceVA.viewAng        = viewAng;
predChoiceVA.viewAng_R      = viewAng_R;
predChoiceVA.viewAng_L      = viewAng_L;
predChoiceVA.choice         = choice;

predChoiceVA.viewAngTOW     = viewAngTOW;
predChoiceVA.mouseIDTOW     = mouseIDTOW;
predChoiceVA.viewAngTOW_R   = viewAngTOW_R;
predChoiceVA.viewAngTOW_L   = viewAngTOW_L;
predChoiceVA.choiceTOW      = choiceTOW;

predChoiceVA.viewAngMG      = viewAngMG;
predChoiceVA.mouseIDMG      = mouseIDMG;
predChoiceVA.viewAngMG_R    = viewAngMG_R;
predChoiceVA.viewAngMG_L    = viewAngMG_L;
predChoiceVA.choiceMG       = choiceMG;

predChoiceVA.viewAngPC      = viewAngPC;
predChoiceVA.mouseIDPC      = mouseIDPC;
predChoiceVA.viewAngPC_R    = viewAngPC_R;
predChoiceVA.viewAngPC_L    = viewAngPC_L;
predChoiceVA.choicePC       = choicePC;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

switch fitType
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    case 'fitglm'
        
        predMatTOW               = double(viewAngTOW)';
        outputTOW                = double(choiceTOW)';
        predMatMG                = double(viewAngMG)';
        outputMG                 = double(choiceMG)';
        predMatPC                = double(viewAngPC)';
        outputPC                 = double(choicePC)';
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %  train on AoE first, test on MG and PC
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% initialize fit variables for
        bestCoeff        = [];     coeff           = [];   fitInfo         = [];
        allCoeffs        = [];     thresh          = [];   pseudoExp       = [];
        accuracyTOW      = [];     accIntervalTOW  = [];   stdTOW          = [];
        accuracyMG       = [];     accIntervalMG   = [];   stdMG           = [];
        accuracyPC       = [];     accIntervalPC   = [];   stdPC           = [];
        predictionTOW    = [];     predAccuracyTOW = [];   predTruthTOW    = [];
        predictionMG     = [];     predAccuracyMG  = [];   predTruthMG     = [];
        predictionPC     = [];     predAccuracyPC  = [];   predTruthPC     = [];
        xMousePredAccTOW = [];     xMousePredAccMG = [];   xMousePredAccPC = [];
        xMouseAccTOW     = [];     xMouseAccMG     = [];   xMouseAccPC     = [];
        xMouseAccSemTOW  = [];     xMouseAccSemMG  = [];   xMouseAccSemPC  = [];
        xMouseAccStdTOW  = [];     xMouseAccStdMG  = [];   xMouseAccStdPC  = [];
        
        %% set up cross-validated index for training on TOW data
        warning('off', 'stats:cvpartition:KFoldMissingGrp');
        randGenerator         = RandStream('mt19937ar', 'Seed', 723176);
        pseudoExp             = cvpartition(outputTOW(:,end), 'KFold', numFolds, randGenerator);
        for iMC = 2:numSamples
            pseudoExp(iMC)      = repartition(pseudoExp(1));
        end
        warning('on', 'stats:cvpartition:KFoldMissingGrp');
        
        for nbin = 1:numel(binIDs)
            %% do fit to all bins train/test light off data according to pseudoExp
            % note can add lasso here instead
            coeff = []; fitInfo = [];
            [coeff, fitInfo] = fitglmRegression(predMatTOW(:,nbin), outputTOW, pseudoExp,numSamples,weighted);
            
            %% appply train set coeff according to pseudoExp to predict alternate task data
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
                    fitInfo.CVTestMGPrediction{iMC,iFold} = glmval(mcCoeff,predMatMG(:,nbin), 'logit');
                    fitInfo.CVTestPCPrediction{iMC,iFold} = glmval(mcCoeff,predMatPC(:,nbin), 'logit');
                    
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
                            threshold     = optimalClassThreshold(predictionTrain, outputTOW(trainSel)>0.5);
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
                    predictionTOW                            = fitInfo.CVTestPrediction{iMC,iFold};
                    predTruthTOW                             = outputTOW(testSel);
                    predictionTOW(predictionTOW < threshold) = 0;
                    predictionTOW(predictionTOW > threshold) = 1;
                    predAccuracyTOW(iMC,iFold)               = sum(predictionTOW == predTruthTOW)./numel(predTruthTOW);
                    
                    mouseTOW = unique(mouseIDTOW);
                    for nMice = 1:numel(mouseTOW)
                        xMousePred      = []; xMousePred      = predictionTOW(mouseIDTOW(testSel) == mouseTOW(nMice));
                        xMousePredTruth = []; xMousePredTruth = outputTOW(testSel==1 & (mouseIDTOW == mouseTOW(nMice))');
                        xMousePred(xMousePred < threshold)    = 0;
                        xMousePred(xMousePred > threshold)    = 1;
                        xMousePredAccTOW(iMC,iFold,nMice)           = sum(xMousePred == xMousePredTruth)./numel(xMousePredTruth);
                    end
                    
                    predictionMG                             = fitInfo.CVTestMGPrediction{iMC,iFold};
                    predTruthMG                              = outputMG;
                    predictionMG(predictionMG < threshold)   = 0;
                    predictionMG(predictionMG > threshold)   = 1;
                    predAccuracyMG(iMC,iFold)                = sum(predictionMG == predTruthMG)./numel(predTruthMG);
                    
                    mouseMG = unique(mouseIDMG);
                    for nMice = 1:numel(mouseMG)
                        xMousePred      = []; xMousePred      = predictionMG(mouseIDMG' == mouseMG(nMice));
                        xMousePredTruth = []; xMousePredTruth = outputMG(mouseIDMG' == mouseMG(nMice));
                        xMousePred(xMousePred < threshold)    = 0;
                        xMousePred(xMousePred > threshold)    = 1;
                        xMousePredAccMG(iMC,iFold,nMice)           = sum(xMousePred == xMousePredTruth)./numel(xMousePredTruth);
                    end
                    
                    predictionPC                             = fitInfo.CVTestPCPrediction{iMC,iFold};
                    predTruthPC                              = outputPC;
                    predictionPC(predictionPC < threshold)   = 0;
                    predictionPC(predictionPC > threshold)   = 1;
                    predAccuracyPC(iMC,iFold)                = sum(predictionPC == predTruthPC)./numel(predTruthPC);
                    
                    mousePC = unique(mouseIDPC);
                    for nMice = 1:numel(mousePC)
                        xMousePred      = []; xMousePred      = predictionPC(mouseIDPC' == mousePC(nMice));
                        xMousePredTruth = []; xMousePredTruth = outputPC(mouseIDPC' == mousePC(nMice));
                        xMousePred(xMousePred < threshold)    = 0;
                        xMousePred(xMousePred > threshold)    = 1;
                        xMousePredAccPC(iMC,iFold,nMice)      = sum(xMousePred == xMousePredTruth)./numel(xMousePredTruth);
                    end
                    
                end
            end
            
            cl             = normcdf([0 1], 0, 1); % confidence level
            accuracyTOW    = nanmean(predAccuracyTOW(:));
            stdTOW         = nanstd(predAccuracyTOW(:));
            accIntervalTOW = quantile(predAccuracyTOW(:), cl);
            xMouseAccTOW   = squeeze(nanmean(nanmean(xMousePredAccTOW,1)));
            xMouseAccSemTOW = nanstd(xMouseAccTOW,0,1)/(sqrt(size(xMouseAccTOW,1)-1));
            xMouseAccStdTOW = nanstd(xMouseAccTOW,0,1);
            
            accuracyMG     = nanmean(predAccuracyMG(:));
            stdMG          = nanstd(predAccuracyMG(:));
            accIntervalMG  = quantile(predAccuracyMG(:), cl);
            xMouseAccMG    = squeeze(nanmean(nanmean(xMousePredAccMG,1)));
            xMouseAccSemMG = nanstd(xMouseAccMG,0,1)/(sqrt(size(xMouseAccMG,1)-1));
            xMouseAccStdMG = nanstd(xMouseAccMG,0,1);
            
            accuracyPC     = nanmean(predAccuracyPC(:));
            stdPC          = nanstd(predAccuracyPC(:));
            accIntervalPC  = quantile(predAccuracyPC(:), cl);
            xMouseAccPC    = squeeze(nanmean(nanmean(xMousePredAccPC,1)));
            xMouseAccSemPC = nanstd(xMouseAccPC,0,1)/(sqrt(size(xMouseAccPC,1)-1));
            xMouseAccStdPC = nanstd(xMouseAccPC,0,1);
            
            %% save relevant variables
            predChoiceVA.trainTOW.fitInfo{nbin}            = fitInfo;
            predChoiceVA.trainTOW.bestCoeff(nbin,:)        = coeff';
            %                 predChoiceVA.mouse(iMouse).allCoeffs{nbin}          = allCoeffs;
            predChoiceVA.trainTOW.accuracyTOW(nbin)        = accuracyTOW;
            predChoiceVA.trainTOW.stdTOW(nbin)             = stdTOW;
            predChoiceVA.trainTOW.accIntervalTOW(nbin,:)   = accIntervalTOW';
            predChoiceVA.trainTOW.xMouseAccTOW(nbin,:)     = xMouseAccTOW;
            predChoiceVA.trainTOW.xMouseAccAvgTOW(nbin)    = nanmean(xMouseAccTOW);
            predChoiceVA.trainTOW.xMouseAccSemTOW(nbin)    = xMouseAccSemTOW;
            predChoiceVA.trainTOW.xMouseAccStdTOW(nbin)    = xMouseAccStdTOW;
            
            predChoiceVA.trainTOW.accuracyMG(nbin)         = accuracyMG;
            predChoiceVA.trainTOW.stdMG(nbin)              = stdMG;
            predChoiceVA.trainTOW.accIntervalMG(nbin,:)    = accIntervalMG';
            predChoiceVA.trainTOW.xMouseAccMG(nbin,:)      = xMouseAccMG;
            predChoiceVA.trainTOW.xMouseAccAvgMG(nbin)     = nanmean(xMouseAccMG);
            predChoiceVA.trainTOW.xMouseAccSemMG(nbin)     = xMouseAccSemMG;
            predChoiceVA.trainTOW.xMouseAccStdMG(nbin)     = xMouseAccStdMG;
            
            predChoiceVA.trainTOW.accuracyPC(nbin)         = accuracyPC;
            predChoiceVA.trainTOW.stdPC(nbin)              = stdPC;
            predChoiceVA.trainTOW.accIntervalPC(nbin,:)    = accIntervalPC';
            predChoiceVA.trainTOW.xMouseAccPC(nbin,:)      = xMouseAccPC;
            predChoiceVA.trainTOW.xMouseAccAvgPC(nbin)     = nanmean(xMouseAccPC);
            predChoiceVA.trainTOW.xMouseAccSemPC(nbin)     = xMouseAccSemPC;
            predChoiceVA.trainTOW.xMouseAccStdPC(nbin)     = xMouseAccStdPC;
            
            predChoiceVA.trainTOW.threshold{nbin}          = thresh;
            predChoiceVA.trainTOW.fitInfo{nbin}.pseudoExp  = pseudoExp;
            
        end
               
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% now train on MG, and test on AoE and PC
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% initialize fit variables for train on MG
        bestCoeff        = [];     coeff           = [];   fitInfo         = [];
        allCoeffs        = [];     thresh          = [];   pseudoExp       = [];
        accuracyTOW      = [];     accIntervalTOW  = [];   stdTOW          = [];
        accuracyMG       = [];     accIntervalMG   = [];   stdMG           = [];
        accuracyPC       = [];     accIntervalPC   = [];   stdPC           = [];
        predictionTOW    = [];     predAccuracyTOW = [];   predTruthTOW    = [];
        predictionMG     = [];     predAccuracyMG  = [];   predTruthMG     = [];
        predictionPC     = [];     predAccuracyPC  = [];   predTruthPC     = [];
        xMousePredAccTOW = [];     xMousePredAccMG = [];   xMousePredAccPC = [];
        xMouseAccTOW     = [];     xMouseAccMG     = [];   xMouseAccPC     = [];
        xMouseAccSemTOW  = [];     xMouseAccSemMG  = [];   xMouseAccSemPC  = [];
        xMouseAccStdTOW  = [];     xMouseAccStdMG  = [];   xMouseAccStdPC  = [];
        
        %% set up cross-validated index for training on MG data
        warning('off', 'stats:cvpartition:KFoldMissingGrp');
        randGenerator         = RandStream('mt19937ar', 'Seed', 723176);
        pseudoExp             = cvpartition(outputMG(:,end), 'KFold', numFolds, randGenerator);
        for iMC = 2:numSamples
            pseudoExp(iMC)      = repartition(pseudoExp(1));
        end
        warning('on', 'stats:cvpartition:KFoldMissingGrp');
        
        
        for nbin = 1:numel(binIDs)
            %% do fit to all bins train/test light off data according to pseudoExp
            % note can add lasso here instead
            coeff = []; fitInfo = [];
            [coeff, fitInfo] = fitglmRegression(predMatMG(:,nbin), outputMG, pseudoExp,numSamples,weighted);
            
            
            %% appply train set coeff according to pseudoExp to predict alternate task data
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
                    fitInfo.CVTestTOWPrediction{iMC,iFold} = glmval(mcCoeff,predMatTOW(:,nbin), 'logit');
                    fitInfo.CVTestPCPrediction{iMC,iFold} = glmval(mcCoeff,predMatPC(:,nbin), 'logit');
                    
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
                            threshold     = optimalClassThreshold(predictionTrain, outputMG(trainSel)>0.5);
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
                    predictionMG                            = fitInfo.CVTestPrediction{iMC,iFold};
                    predTruthMG                             = outputMG(testSel);
                    predictionMG(predictionMG < threshold) = 0;
                    predictionMG(predictionMG > threshold) = 1;
                    predAccuracyMG(iMC,iFold)               = sum(predictionMG == predTruthMG)./numel(predTruthMG);
                    
                    mouseMG = unique(mouseIDMG);
                    for nMice = 1:numel(mouseMG)
                        xMousePred      = []; xMousePred      = predictionMG(mouseIDMG(testSel) == mouseMG(nMice));
                        xMousePredTruth = []; xMousePredTruth = outputMG(testSel==1 & (mouseIDMG == mouseMG(nMice))');
                        xMousePred(xMousePred < threshold)    = 0;
                        xMousePred(xMousePred > threshold)    = 1;
                        xMousePredAccMG(iMC,iFold,nMice)      = sum(xMousePred == xMousePredTruth)./numel(xMousePredTruth);
                    end
                    
                    
                    predictionTOW                             = fitInfo.CVTestTOWPrediction{iMC,iFold};
                    predTruthTOW                              = outputTOW;
                    predictionTOW(predictionTOW < threshold)   = 0;
                    predictionTOW(predictionTOW > threshold)   = 1;
                    predAccuracyTOW(iMC,iFold)                = sum(predictionTOW == predTruthTOW)./numel(predTruthTOW);
                    
                    mouseTOW = unique(mouseIDTOW);
                    for nMice = 1:numel(mouseTOW)
                        xMousePred      = []; xMousePred      = predictionTOW(mouseIDTOW' == mouseTOW(nMice));
                        xMousePredTruth = []; xMousePredTruth = outputTOW(mouseIDTOW' == mouseTOW(nMice));
                        xMousePred(xMousePred < threshold)    = 0;
                        xMousePred(xMousePred > threshold)    = 1;
                        xMousePredAccTOW(iMC,iFold,nMice)     = sum(xMousePred == xMousePredTruth)./numel(xMousePredTruth);
                    end                   
                    
                    predictionPC                             = fitInfo.CVTestPCPrediction{iMC,iFold};
                    predTruthPC                              = outputPC;
                    predictionPC(predictionPC < threshold)   = 0;
                    predictionPC(predictionPC > threshold)   = 1;
                    predAccuracyPC(iMC,iFold)                = sum(predictionPC == predTruthPC)./numel(predTruthPC);
                    
                    mousePC = unique(mouseIDPC);
                    for nMice = 1:numel(mousePC)
                        xMousePred      = []; xMousePred      = predictionPC(mouseIDPC' == mousePC(nMice));
                        xMousePredTruth = []; xMousePredTruth = outputPC(mouseIDPC' == mousePC(nMice));
                        xMousePred(xMousePred < threshold)    = 0;
                        xMousePred(xMousePred > threshold)    = 1;
                        xMousePredAccPC(iMC,iFold,nMice)      = sum(xMousePred == xMousePredTruth)./numel(xMousePredTruth);
                    end
                    
                end
            end
            
            cl             = normcdf([0 1], 0, 1); % confidence level
            accuracyMG    = nanmean(predAccuracyMG(:));
            stdMG         = nanstd(predAccuracyMG(:));
            accIntervalMG = quantile(predAccuracyMG(:), cl);
            xMouseAccMG    = squeeze(nanmean(nanmean(xMousePredAccMG,1)));
            xMouseAccSemMG = nanstd(xMouseAccMG,0,1)/(sqrt(size(xMouseAccMG,1)-1));
            xMouseAccStdMG = nanstd(xMouseAccMG,0,1);
            
            accuracyTOW     = nanmean(predAccuracyTOW(:));
            stdTOW          = nanstd(predAccuracyTOW(:));
            accIntervalTOW  = quantile(predAccuracyTOW(:), cl);
            xMouseAccTOW    = squeeze(nanmean(nanmean(xMousePredAccTOW,1)));
            xMouseAccSemTOW = nanstd(xMouseAccTOW,0,1)/(sqrt(size(xMouseAccTOW,1)-1));
            xMouseAccStdTOW = nanstd(xMouseAccTOW,0,1);
            
            accuracyPC     = nanmean(predAccuracyPC(:));
            stdPC          = nanstd(predAccuracyPC(:));
            accIntervalPC  = quantile(predAccuracyPC(:), cl);
            xMouseAccPC    = squeeze(nanmean(nanmean(xMousePredAccPC,1)));
            xMouseAccSemPC = nanstd(xMouseAccPC,0,1)/(sqrt(size(xMouseAccPC,1)-1));
            xMouseAccStdPC = nanstd(xMouseAccPC,0,1);
            
            %% save relevant variables
            predChoiceVA.trainMG.fitInfo{nbin}            = fitInfo;
            predChoiceVA.trainMG.bestCoeff(nbin,:)        = coeff';
            %                 predChoiceVA.mouse(iMouse).allCoeffs{nbin}          = allCoeffs;
            predChoiceVA.trainMG.accuracyTOW(nbin)        = accuracyTOW;
            predChoiceVA.trainMG.stdTOW(nbin)             = stdTOW;
            predChoiceVA.trainMG.accIntervalTOW(nbin,:)   = accIntervalTOW';
            predChoiceVA.trainMG.xMouseAccTOW(nbin,:)     = xMouseAccTOW;
            predChoiceVA.trainMG.xMouseAccAvgTOW(nbin)    = nanmean(xMouseAccTOW);
            predChoiceVA.trainMG.xMouseAccSemTOW(nbin)    = xMouseAccSemTOW;
            predChoiceVA.trainMG.xMouseAccStdTOW(nbin)    = xMouseAccStdTOW;
            
            predChoiceVA.trainMG.accuracyMG(nbin)         = accuracyMG;
            predChoiceVA.trainMG.stdMG(nbin)              = stdMG;
            predChoiceVA.trainMG.accIntervalMG(nbin,:)    = accIntervalMG';
            predChoiceVA.trainMG.xMouseAccMG(nbin,:)      = xMouseAccMG;
            predChoiceVA.trainMG.xMouseAccAvgMG(nbin)     = nanmean(xMouseAccMG);
            predChoiceVA.trainMG.xMouseAccSemMG(nbin)     = xMouseAccSemMG;
            predChoiceVA.trainMG.xMouseAccStdMG(nbin)     = xMouseAccStdMG;
            
            predChoiceVA.trainMG.accuracyPC(nbin)         = accuracyPC;
            predChoiceVA.trainMG.stdPC(nbin)              = stdPC;
            predChoiceVA.trainMG.accIntervalPC(nbin,:)    = accIntervalPC';
            predChoiceVA.trainMG.xMouseAccPC(nbin,:)      = xMouseAccPC;
            predChoiceVA.trainMG.xMouseAccAvgPC(nbin)     = nanmean(xMouseAccPC);
            predChoiceVA.trainMG.xMouseAccSemPC(nbin)     = xMouseAccSemPC;
            predChoiceVA.trainMG.xMouseAccStdPC(nbin)     = xMouseAccStdPC;
            
            predChoiceVA.trainMG.threshold{nbin}          = thresh;
            predChoiceVA.trainMG.fitInfo{nbin}.pseudoExp  = pseudoExp;
            
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% train on PC and test on MG and AoE
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% initialize fit variables for
        bestCoeff        = [];     coeff           = [];   fitInfo         = [];
        allCoeffs        = [];     thresh          = [];   pseudoExp       = [];
        accuracyTOW      = [];     accIntervalTOW  = [];   stdTOW          = [];
        accuracyMG       = [];     accIntervalMG   = [];   stdMG           = [];
        accuracyPC       = [];     accIntervalPC   = [];   stdPC           = [];
        predictionTOW    = [];     predAccuracyTOW = [];   predTruthTOW    = [];
        predictionMG     = [];     predAccuracyMG  = [];   predTruthMG     = [];
        predictionPC     = [];     predAccuracyPC  = [];   predTruthPC     = [];
        xMousePredAccTOW = [];     xMousePredAccMG = [];   xMousePredAccPC = [];
        xMouseAccTOW     = [];     xMouseAccMG     = [];   xMouseAccPC     = [];
        xMouseAccSemTOW  = [];     xMouseAccSemMG  = [];   xMouseAccSemPC  = [];
        xMouseAccStdTOW  = [];     xMouseAccStdMG  = [];   xMouseAccStdPC  = [];
        
        %% set up cross-validated index for training on PC data FIRST
        warning('off', 'stats:cvpartition:KFoldMissingGrp');
        randGenerator         = RandStream('mt19937ar', 'Seed', 723176);
        pseudoExp             = cvpartition(outputPC(:,end), 'KFold', numFolds, randGenerator);
        for iMC = 2:numSamples
            pseudoExp(iMC)      = repartition(pseudoExp(1));
        end
        warning('on', 'stats:cvpartition:KFoldMissingGrp');
        
        
        for nbin = 1:numel(binIDs)
            %% do fit to all bins train/test light off data according to pseudoExp
            % note can add lasso here instead
            coeff = []; fitInfo = [];
            [coeff, fitInfo] = fitglmRegression(predMatPC(:,nbin), outputPC, pseudoExp,numSamples,weighted);
            
            
            %% appply train set coeff according to pseudoExp to predict alternate task data
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
                    fitInfo.CVTestTOWPrediction{iMC,iFold} = glmval(mcCoeff,predMatTOW(:,nbin), 'logit');
                    fitInfo.CVTestMGPrediction{iMC,iFold}  = glmval(mcCoeff,predMatMG(:,nbin), 'logit');
                    
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
                            threshold     = optimalClassThreshold(predictionTrain, outputPC(trainSel)>0.5);
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
                    predictionPC                            = fitInfo.CVTestPrediction{iMC,iFold};
                    predTruthPC                             = outputPC(testSel);
                    predictionPC(predictionPC < threshold) = 0;
                    predictionPC(predictionPC > threshold) = 1;
                    predAccuracyPC(iMC,iFold)               = sum(predictionPC == predTruthPC)./numel(predTruthPC);
                    
                    mousePC = unique(mouseIDPC);
                    for nMice = 1:numel(mousePC)
                        xMousePred      = []; xMousePred      = predictionPC(mouseIDPC(testSel) == mousePC(nMice));
                        xMousePredTruth = []; xMousePredTruth = outputPC(testSel==1 & (mouseIDPC == mousePC(nMice))');
                        xMousePred(xMousePred < threshold)    = 0;
                        xMousePred(xMousePred > threshold)    = 1;
                        xMousePredAccPC(iMC,iFold,nMice)      = sum(xMousePred == xMousePredTruth)./numel(xMousePredTruth);
                    end
                    
                    
                    predictionTOW                             = fitInfo.CVTestTOWPrediction{iMC,iFold};
                    predTruthTOW                              = outputTOW;
                    predictionTOW(predictionTOW < threshold)   = 0;
                    predictionTOW(predictionTOW > threshold)   = 1;
                    predAccuracyTOW(iMC,iFold)                = sum(predictionTOW == predTruthTOW)./numel(predTruthTOW);
                    
                    mouseTOW = unique(mouseIDTOW);
                    for nMice = 1:numel(mouseTOW)
                        xMousePred      = []; xMousePred      = predictionTOW(mouseIDTOW' == mouseTOW(nMice));
                        xMousePredTruth = []; xMousePredTruth = outputTOW(mouseIDTOW' == mouseTOW(nMice));
                        xMousePred(xMousePred < threshold)    = 0;
                        xMousePred(xMousePred > threshold)    = 1;
                        xMousePredAccTOW(iMC,iFold,nMice)     = sum(xMousePred == xMousePredTruth)./numel(xMousePredTruth);
                    end
                    
                    
                    predictionMG                             = fitInfo.CVTestMGPrediction{iMC,iFold};
                    predTruthMG                              = outputMG;
                    predictionMG(predictionMG < threshold)   = 0;
                    predictionMG(predictionMG > threshold)   = 1;
                    predAccuracyMG(iMC,iFold)                = sum(predictionMG == predTruthMG)./numel(predTruthMG);
                    
                    mouseMG = unique(mouseIDMG);
                    for nMice = 1:numel(mouseMG)
                        xMousePred      = []; xMousePred      = predictionMG(mouseIDMG' == mouseMG(nMice));
                        xMousePredTruth = []; xMousePredTruth = outputMG(mouseIDMG' == mouseMG(nMice));
                        xMousePred(xMousePred < threshold)    = 0;
                        xMousePred(xMousePred > threshold)    = 1;
                        xMousePredAccMG(iMC,iFold,nMice)     = sum(xMousePred == xMousePredTruth)./numel(xMousePredTruth);
                    end
                    
                end
            end
            
            cl             = normcdf([0 1], 0, 1); % confidence level
            accuracyMG     = nanmean(predAccuracyMG(:));
            stdMG          = nanstd(predAccuracyMG(:));
            accIntervalMG  = quantile(predAccuracyMG(:), cl);
            xMouseAccMG    = squeeze(nanmean(nanmean(xMousePredAccMG,1)));
            xMouseAccSemMG = nanstd(xMouseAccMG,0,1)/(sqrt(size(xMouseAccMG,1)-1));
            xMouseAccStdMG = nanstd(xMouseAccMG,0,1);
            
            accuracyPC     = nanmean(predAccuracyPC(:));
            stdPC          = nanstd(predAccuracyPC(:));
            accIntervalPC  = quantile(predAccuracyPC(:), cl);
            xMouseAccPC    = squeeze(nanmean(nanmean(xMousePredAccPC,1)));
            xMouseAccSemPC = nanstd(xMouseAccPC,0,1)/(sqrt(size(xMouseAccPC,1)-1));
            xMouseAccStdPC = nanstd(xMouseAccPC,0,1);
            
            accuracyTOW     = nanmean(predAccuracyTOW(:));
            stdTOW          = nanstd(predAccuracyTOW(:));
            accIntervalTOW  = quantile(predAccuracyTOW(:), cl);
            xMouseAccTOW    = squeeze(nanmean(nanmean(xMousePredAccTOW,1)));
            xMouseAccSemTOW = nanstd(xMouseAccTOW,0,1)/(sqrt(size(xMouseAccTOW,1)-1));
            xMouseAccStdTOW = nanstd(xMouseAccTOW,0,1);
            
            %% save relevant variables
            predChoiceVA.trainPC.fitInfo{nbin}            = fitInfo;
            predChoiceVA.trainPC.bestCoeff(nbin,:)        = coeff';
            %                 predChoiceVA.mouse(iMouse).allCoeffs{nbin}          = allCoeffs;
            predChoiceVA.trainPC.accuracyTOW(nbin)        = accuracyTOW;
            predChoiceVA.trainPC.stdTOW(nbin)             = stdTOW;
            predChoiceVA.trainPC.accIntervalTOW(nbin,:)   = accIntervalTOW';
            predChoiceVA.trainPC.xMouseAccTOW(nbin,:)     = xMouseAccTOW;
            predChoiceVA.trainPC.xMouseAccAvgTOW(nbin)    = nanmean(xMouseAccTOW);
            predChoiceVA.trainPC.xMouseAccSemTOW(nbin)    = xMouseAccSemTOW;
            predChoiceVA.trainPC.xMouseAccStdTOW(nbin)    = xMouseAccStdTOW;
            
            predChoiceVA.trainPC.accuracyMG(nbin)         = accuracyMG;
            predChoiceVA.trainPC.stdMG(nbin)              = stdMG;
            predChoiceVA.trainPC.accIntervalMG(nbin,:)    = accIntervalMG';
            predChoiceVA.trainPC.xMouseAccMG(nbin,:)      = xMouseAccMG;
            predChoiceVA.trainPC.xMouseAccAvgMG(nbin)     = nanmean(xMouseAccMG);
            predChoiceVA.trainPC.xMouseAccSemMG(nbin)     = xMouseAccSemMG;
            predChoiceVA.trainPC.xMouseAccStdMG(nbin)     = xMouseAccStdMG;
            
            predChoiceVA.trainPC.accuracyPC(nbin)         = accuracyPC;
            predChoiceVA.trainPC.stdPC(nbin)              = stdPC;
            predChoiceVA.trainPC.accIntervalPC(nbin,:)    = accIntervalPC';
            predChoiceVA.trainPC.xMouseAccPC(nbin,:)      = xMouseAccPC;
            predChoiceVA.trainPC.xMouseAccAvgPC(nbin)     = nanmean(xMouseAccPC);
            predChoiceVA.trainPC.xMouseAccSemPC(nbin)     = xMouseAccSemPC;
            predChoiceVA.trainPC.xMouseAccStdPC(nbin)     = xMouseAccStdPC;
            
            predChoiceVA.trainPC.threshold{nbin}          = thresh;
            predChoiceVA.trainPC.fitInfo{nbin}.pseudoExp  = pseudoExp;
            
        end
        
end

end