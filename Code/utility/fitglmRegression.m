function  [coeff, FitInfo] = fitglmRegression(X, y, numCVFolds, numCVSamples,biasWeighted) 

% [coeff, fitInfo] = ridgeRegression(X, y, lambdaRatoo, numCVFolds, numCVSamples) 
% wrapper to perform ridge regression using a CV object and having similar
% output structure to lasso. Calls the matlab function ridge.m
% X is predictor matrix, y is response vector, lambdaRatio is vector
% containing the lambda hyperparameter. numCVfolds and numCVsamples are to
% set up k-fold CV object. alternatively numCVFolds can be the CV object
% itself
% coeff is n params x n lambda matrix for full data solution, fitInfo
% contains CV info similar to lasso and elasticNetRegression

 %% Default arguments

if nargin < 3 || isempty(numCVFolds)
  numCVFolds            = 5;
end
if nargin < 4 || isempty(numCVSamples)
  numCVSamples          = 10;
end

if nargin < 5 || isempty(biasWeighted)
  biasWeighted          = false;
end

%% Cross-validation setup
if isnumeric(numCVFolds)
  warning('off', 'stats:cvpartition:KFoldMissingGrp');
  pseudoExp           = cvpartition(y, 'KFold', numCVFolds);
  for iMC = 2:numCVSamples
    pseudoExp(iMC)    = repartition(pseudoExp(1));
  end
  warning('on', 'stats:cvpartition:KFoldMissingGrp');
else
  pseudoExp           = numCVFolds;
  numCVFolds          = pseudoExp(1).NumTestSets;
end

%% central value fits
FitInfo.CVExperiments     = pseudoExp;
FitInfo.Coefficients      = glmfit(X,y,'binomial','logit');
coeff                     = FitInfo.Coefficients;

%% Cross validated goodness-of-fit
FitInfo.CVTrainPrediction = cell(numel(pseudoExp), numCVFolds);
FitInfo.CVTestPrediction  = cell(numel(pseudoExp), numCVFolds);

for iMC = 1:numel(pseudoExp)
  for iFold = 1:numCVFolds
    % Select train/test sets with user-specified pruning if so desired
    if isfield(pseudoExp(iMC).training(iFold),'idx')
      iTrain              = pseudoExp(iMC).training(iFold).idx;
      iTest               = pseudoExp(iMC).test(iFold).idx;
    else
      iTrain              = pseudoExp(iMC).training(iFold);
      iTest               = pseudoExp(iMC).test(iFold);
    end
    % Perform fit for this CV experiment and generate/test prediction
    
    if biasWeighted == 1
        wts=zeros(numel(y(iTrain)),1);
        class1=find(y(iTrain)==1);  %right
        class2=find(y(iTrain)==0);  %left
        
        wts(class1)=1/numel(class1);
        wts(class2)=1/numel(class2);
        wts=wts./max(wts);
        offset = zeros(numel(y(iTrain)),1);
        
        mcCoeff               = glmfit(X(iTrain,:), y(iTrain),'binomial','logit','off',offset,wts);
    else
        mcCoeff               = glmfit(X(iTrain,:), y(iTrain),'binomial','logit');
    end
    
%         nR     = []; nR        = sum(y(iTrain)==1);               % number of "R" choices
%         nL     = []; nL        = sum(y(iTrain)==0);              % number of "L" choices
%         choWts = []; choWts    = zeros(numel(y(iTrain)),1);       % compute weights
%         choWts(y(iTrain)==0)   = nR/nR;                           % weights for L trials
%         choWts(y(iTrain)==1)   = nL/nR;                           % weights for R trials
%         %     x2 = x - wts'*x./(sum(wts));                        % re-center x so weighted mean is zero
%         temp1  = []; temp1     = X(iTrain) - ((choWts.*X(iTrain))./(sum(choWts)));
%         X(iTrain,:)              = temp1;
%         
%         XX            = [ones(numel(y(iTrain)),1) X(iTrain,:)]; %add those ones
%         nfpwnllfun    = @(theta)(-(choWts.*y(iTrain))'*(XX*theta) + choWts'*log(1+exp(XX*theta)));
%         theta0        = [0;0]; % initial parameters
%         choice_betas  = fminunc(nfpwnllfun, theta0);     
    
    FitInfo.mcCoeff{iMC,iFold}           = mcCoeff;
    % train prediction
    FitInfo.CVTrainPrediction{iMC,iFold} = glmval(mcCoeff,X(iTrain,:),'logit');
    FitInfo.CVTestPrediction{iMC,iFold}  = glmval(mcCoeff,X(iTest,:), 'logit');

    % fitting without bias term
%     temp = []; temp = glmval(mcCoeff(end),X(iTrain,:),'logit');
%     FitInfo.CVTrainPrediction{iMC,iFold} = temp(:,end);
%     temp = []; temp = glmval(mcCoeff(end),X(iTest,:),'logit');
%     FitInfo.CVTestPrediction{iMC,iFold}  = temp(:,end);
    

  end
end

end