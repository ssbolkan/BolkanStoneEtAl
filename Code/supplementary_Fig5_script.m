%% Script to generate ExtData Fig 15 data plots
% learning in permanent cues/no distractors shaping or aoe/no distractors
% shaping
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load learning metalogs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
cd(globalParams.dataPath)
load('taskShapingData.mat')

%%
shapingData     = aoeShaping;
shapingProtocol = 'PoissonBlocksCondensed3m'; % PoissonBlocksCondensed3m PoissonBlocksPermanentCues
aoeLearning     = getTaskShaping(shapingData, shapingProtocol, 0);

shapingData     = pcShaping;
shapingProtocol = 'PoissonBlocksPermanentCues'; % PoissonBlocksCondensed3m PoissonBlocksPermanentCues
pcLearning      = getTaskShaping(shapingData, shapingProtocol, 0); 

%% PLOTTING
shapingProtocol = 'PoissonBlocksCondensed3m'; % PoissonBlocksCondensed3m PoissonBlocksPermanentCues
plotTaskShaping(aoeLearning, shapingProtocol)

shapingProtocol = 'PoissonBlocksPermanentCues'; % PoissonBlocksCondensed3m PoissonBlocksPermanentCues
plotTaskShaping(pcLearning, shapingProtocol)