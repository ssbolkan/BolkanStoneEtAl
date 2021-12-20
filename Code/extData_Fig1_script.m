%% Script to generate Extended Data Fig 1C-H data plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% effects of indirect or direct pathway inhibition in mice on a running wheel 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTE THAT Code/Utility subfolder contains compiled versions of
% SUPPORTING MEX FUNCTIONS from 'Matlab Offline Files SDK' AND
% 'TankMouseVR' REPOS. If a local recompile is necessary... 
% (1) CLONE REPOS AND ADD TO YOUR MATLAB PATH
% (2) MAKE SURE YOU HAVE A COMPATABLE VISUAL STUDIO (OR OTHER) C++ COMPILER FOR
% YOUR VERSION OF MATLAB
% (2) RUN compile_utilities_pipeline() TO COMPILE MEX FILES IN
% 'TankMouseVR-master'
% (3) SET YOUR DIRECTORY TO '[your local path]\Matlab Offline Files
% SDK\mexPlex\' AND RUN build_and_verify_mexPlex() (ALSO SEE 'HOW TO BUILD
% mexPlex.txt' IN THE SAME FOLDER FOR MORE DETAIL)
% NOTE: mexPlex only necessary if loading .plx data, which is not the case
% if using spikeTS, laserTS, and unit waveforms saved in spikeDataAll.mat

clear all
cd(globalParams.dataPath)
load('runningWheelSpikes.mat')

%% plots PSTH and raster for example units in panel ExtData 1C and 1F
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% relies on plexon matlab SDK and compiling mex functions
% Under Home tab --> add ons --> search MATLAB Support for MinGW-w64 C/C++ Compiler 

xlims   = []; xlims   = [-5 10];
binsize = []; binsize = 0.25;
plotWaveform = 1; % set to 1 for seperate fig of unit waveform
% for rasters with too many spikes, spikes were subsampled in paper figs to make editable in 
% illustrator for formatting outside matlab. make 'doSubSamp == 1' to get exact fig raster
% for all significantly light modulated units exact subsample is also
% indicated in raster titles (parenthetical) of ExtData Fig3
doSubSamp = 1;

% plot panel 1C left
if ~doSubSamp; subsamp = 1; else; subsamp = 1; end
exampleFile = spikeDataAll(15);
psthRasterExample_metaFile(exampleFile,xlims,binsize,subsamp,[],plotWaveform)
% plot panel 1C right
if ~doSubSamp; subsamp = 1; else; subsamp = 8; end
exampleFile = spikeDataAll(19);
psthRasterExample_metaFile(exampleFile,xlims,binsize,subsamp,[],plotWaveform)

% plot panel 1F left
if ~doSubSamp; subsamp = 1; else; subsamp = 1; end
exampleFile = spikeDataAll(109);
psthRasterExample_metaFile(exampleFile,xlims,binsize,subsamp,[],plotWaveform)
% plot panel 1F right
if ~doSubSamp; subsamp = 1; else; subsamp = 2; end
exampleFile = spikeDataAll(70);
psthRasterExample_metaFile(exampleFile,xlims,binsize,subsamp,[],plotWaveform)

cd(globalParams.dataPath)
clearvars -except spikeDataAll

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plots delta firing rate histograms and pie charts for ExtData 1D and 1G
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first 1D : sets plot and title as indirect pathway on vs pre
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mouse = []; neuronIdx = []; alpha = []; data2comp = []; sig2comp = [];
mouse      = 'a2a'; % this with next line sets index to indirect pathway
neuronIdx  = contains([spikeDataAll.animalName],mouse);
alpha      = 0.05 / sum(neuronIdx);    % pvalue corrected for # neurons
data2Comp  = [spikeDataAll.diffOnOff]; % data to plot
sig2Comp   = [spikeDataAll.pOffOn];    % p value to use
titlestr   = 'indirect: on vs. pre';
xlabelstr  = 'firing rate (Hz, on-pre)';

% gets data index
decIndex  = neuronIdx == 1    & ...
            sig2Comp  < alpha & ...
            data2Comp < 0 ;
incIndex  = neuronIdx == 1    & ...
            sig2Comp  < alpha & ...
            data2Comp > 0 ;    
nsIndex   = neuronIdx == 1    & ...
            sig2Comp  > alpha;   
        
% initialize histogram variables
n  = [];  edges = []; bin  = [];
n2 = []; edges2 = []; bin2 = [];
n3 = []; edges3 = []; bin3 = [];
nsData  = []; incData = []; decData = [];

% gets actual data values
nsData  = data2Comp(nsIndex == 1);
incData = data2Comp(incIndex == 1);
decData = data2Comp(decIndex == 1);

% makes histogram counts
[n,edges,bin]    = histcounts(incData, -20:1:20);
[n2,edges2,bin2] = histcounts(decData, -20:1:20);
[n3,edges3,bin3] = histcounts(nsData, -20:1:20);

% organizes data for plot
nPlot = []; edgePlot = [];
nPlot(:,1)    = n;     nPlot(:,2)    = n2;     nPlot(:,3)    = n3;    
edgePlot(:,1) = edges; edgePlot(:,2) = edges2; edgePlot(:,3) = edges3;
edgePlot = edgePlot';

% plots histogram stacked
figure
subplot(2,2,1)
H = bar(edgePlot(1,1:end-1),nPlot,'stacked');
hold on
H(1).FaceColor = 'g';
H(2).FaceColor = 'r';
H(3).FaceColor = 'none';
hold on
box off
ylabel('# neurons')
xlabel(xlabelstr)
set(gca,'TickDir','out')
title(titlestr)
xlim([-21 21])
ylim([0 30])
plot(zeros(1,30),1:30,'k:')

% plots pie chart
nonsig = []; sigInc = []; sigDec = []; labels = [];
nonsig = length(nsData);
sigInc = length(incData);
sigDec = length(decData);
totalN = sum([nonsig sigInc sigDec]);
explode = [0 0 1];
subplot(2,2,3)
labels = {num2str(sigInc),num2str(sigDec),num2str(nonsig)};
pie([sigInc/totalN sigDec/totalN nonsig/totalN],explode,labels)
colormap([0 1 0; 1 0 0; 0.9 0.9 0.9])
title(titlestr)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sets plot and title as indirect pathway POST vs pre
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mouse = []; neuronIdx = []; alpha = []; data2comp = []; sig2comp = [];
mouse      = 'a2a'; % still select indirect pathway
neuronIdx  = contains([spikeDataAll.animalName],mouse);
alpha      = 0.05 / sum(neuronIdx); % correct pval for number of neurons
data2Comp  = [spikeDataAll.diffPostPre]; % now selecting post vs pre
sig2Comp   = [spikeDataAll.pPrePost];    % update pval for post vs pre
titlestr   = 'indirect: post vs. pre';
xlabelstr  = 'firing rate (Hz, post-pre)';

% gets data index
decIndex  = neuronIdx == 1    & ...
            sig2Comp  < alpha & ...
            data2Comp < 0 ;
incIndex  = neuronIdx == 1    & ...
            sig2Comp  < alpha & ...
            data2Comp > 0 ;    
nsIndex   = neuronIdx == 1    & ...
            sig2Comp  > alpha;   
        
% initialize histogram variables        
n  = [];  edges = []; bin  = [];
n2 = []; edges2 = []; bin2 = [];
n3 = []; edges3 = []; bin3 = [];
nsData  = []; incData = []; decData = [];

% gets actual data values
nsData  = data2Comp(nsIndex == 1);
incData = data2Comp(incIndex == 1);
decData = data2Comp(decIndex == 1);

% makes histogram counts
[n,edges,bin]    = histcounts(incData, -20:1:20);
[n2,edges2,bin2] = histcounts(decData, -20:1:20);
[n3,edges3,bin3] = histcounts(nsData, -20:1:20);

% organizes data for plot
nPlot = []; edgePlot = [];
nPlot(:,1)    = n;     nPlot(:,2)    = n2;     nPlot(:,3)    = n3;    
edgePlot(:,1) = edges; edgePlot(:,2) = edges2; edgePlot(:,3) = edges3;
edgePlot = edgePlot';

% plots histogram stacked
subplot(2,2,2)
H = bar(edgePlot(1,1:end-1),nPlot,'stacked');
hold on
H(1).FaceColor = 'g';
H(2).FaceColor = 'r';
H(3).FaceColor = 'none';
hold on
box off
ylabel('# neurons')
xlabel(xlabelstr)
set(gca,'TickDir','out')
title(titlestr)
xlim([-21 21])
ylim([0 30])
plot(zeros(1,30),1:30,'k:')

% plots pie chart
nonsig = []; sigInc = []; sigDec = []; labels = [];
nonsig = length(nsData);
sigInc = length(incData);
sigDec = length(decData);
totalN = sum([nonsig sigInc sigDec]);
explode = [0 0 1];
subplot(2,2,4)
labels = {num2str(sigInc),num2str(sigDec),num2str(nonsig)};
pie([sigInc/totalN sigDec/totalN nonsig/totalN],explode,labels)
colormap([0 1 0; 1 0 0; 0.9 0.9 0.9])
title(titlestr)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1G next : sets plot and title as DIRECT pathway on vs pre
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mouse = []; neuronIdx = []; alpha = []; data2comp = []; sig2comp = [];
mouse      = 'd1'; % this with next line sets index to direct pathway
neuronIdx  = contains([spikeDataAll.animalName],mouse);
alpha      = 0.05 / sum(neuronIdx);    % pvalue corrected for # neurons
data2Comp  = [spikeDataAll.diffOnOff]; % data to plot
sig2Comp   = [spikeDataAll.pOffOn];    % p value to use
titlestr   = 'direct: on vs. pre';
xlabelstr  = 'firing rate (Hz, on-pre)';

% gets data index
decIndex  = neuronIdx == 1    & ...
            sig2Comp  < alpha & ...
            data2Comp < 0 ;
incIndex  = neuronIdx == 1    & ...
            sig2Comp  < alpha & ...
            data2Comp > 0 ;    
nsIndex   = neuronIdx == 1    & ...
            sig2Comp  > alpha;   
        
% initialize histogram variables
n  = [];  edges = []; bin  = [];
n2 = []; edges2 = []; bin2 = [];
n3 = []; edges3 = []; bin3 = [];
nsData  = []; incData = []; decData = [];

% gets actual data values
nsData  = data2Comp(nsIndex == 1);
incData = data2Comp(incIndex == 1);
decData = data2Comp(decIndex == 1);

% makes histogram counts
[n,edges,bin]    = histcounts(incData, -20:1:20);
[n2,edges2,bin2] = histcounts(decData, -20:1:20);
[n3,edges3,bin3] = histcounts(nsData, -20:1:20);

% organizes data for plot
nPlot = []; edgePlot = [];
nPlot(:,1)    = n;     nPlot(:,2)    = n2;     nPlot(:,3)    = n3;    
edgePlot(:,1) = edges; edgePlot(:,2) = edges2; edgePlot(:,3) = edges3;
edgePlot = edgePlot';

% plots histogram stacked
figure
subplot(2,2,1)
H = bar(edgePlot(1,1:end-1),nPlot,'stacked');
hold on
H(1).FaceColor = 'g';
H(2).FaceColor = 'r';
H(3).FaceColor = 'none';
hold on
box off
ylabel('# neurons')
xlabel(xlabelstr)
set(gca,'TickDir','out')
title(titlestr)
xlim([-21 21])
ylim([0 30])
plot(zeros(1,30),1:30,'k:')

% plots pie chart
nonsig = []; sigInc = []; sigDec = []; labels = [];
nonsig = length(nsData);
sigInc = length(incData);
sigDec = length(decData);
totalN = sum([nonsig sigInc sigDec]);
explode = [0 0 1];
subplot(2,2,3)
labels = {num2str(sigInc),num2str(sigDec),num2str(nonsig)};
pie([sigInc/totalN sigDec/totalN nonsig/totalN],explode,labels)
colormap([0 1 0; 1 0 0; 0.9 0.9 0.9])
title(titlestr)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sets plot and title as DIRECT pathway POST vs pre
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mouse = []; neuronIdx = []; alpha = []; data2comp = []; sig2comp = [];
mouse      = 'd1'; % still select indirect pathway
neuronIdx  = contains([spikeDataAll.animalName],mouse);
alpha      = 0.05 / sum(neuronIdx); % correct pval for number of neurons
data2Comp  = [spikeDataAll.diffPostPre]; % now selecting post vs pre
sig2Comp   = [spikeDataAll.pPrePost];    % update pval for post vs pre
titlestr   = 'direct: post vs. pre';
xlabelstr  = 'firing rate (Hz, post-pre)';

% gets data index
decIndex  = neuronIdx == 1    & ...
            sig2Comp  < alpha & ...
            data2Comp < 0 ;
incIndex  = neuronIdx == 1    & ...
            sig2Comp  < alpha & ...
            data2Comp > 0 ;    
nsIndex   = neuronIdx == 1    & ...
            sig2Comp  > alpha;   
        
% initialize histogram variables        
n  = [];  edges = []; bin  = [];
n2 = []; edges2 = []; bin2 = [];
n3 = []; edges3 = []; bin3 = [];
nsData  = []; incData = []; decData = [];

% gets actual data values
nsData  = data2Comp(nsIndex == 1);
incData = data2Comp(incIndex == 1);
decData = data2Comp(decIndex == 1);

% makes histogram counts
[n,edges,bin]    = histcounts(incData, -20:1:20);
[n2,edges2,bin2] = histcounts(decData, -20:1:20);
[n3,edges3,bin3] = histcounts(nsData, -20:1:20);

% organizes data for plot
nPlot = []; edgePlot = [];
nPlot(:,1)    = n;     nPlot(:,2)    = n2;     nPlot(:,3)    = n3;    
edgePlot(:,1) = edges; edgePlot(:,2) = edges2; edgePlot(:,3) = edges3;
edgePlot = edgePlot';

% plots histogram stacked
subplot(2,2,2)
H = bar(edgePlot(1,1:end-1),nPlot,'stacked');
hold on
H(1).FaceColor = 'g';
H(2).FaceColor = 'r';
H(3).FaceColor = 'none';
hold on
box off
ylabel('# neurons')
xlabel(xlabelstr)
set(gca,'TickDir','out')
title(titlestr)
xlim([-21 21])
ylim([0 30])
plot(zeros(1,30),1:30,'k:')

% plots pie chart
nonsig = []; sigInc = []; sigDec = []; labels = [];
nonsig = length(nsData);
sigInc = length(incData);
sigDec = length(decData);
totalN = sum([nonsig sigInc sigDec]);
explode = [0 0 1];
subplot(2,2,4)
labels = {num2str(sigInc),num2str(sigDec),num2str(nonsig)};
pie([sigInc/totalN sigDec/totalN nonsig/totalN],explode,labels)
colormap([0 1 0; 1 0 0; 0.9 0.9 0.9])
title(titlestr)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plots z scored firing rate for ExtData 1E and 1H
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
zFRmat = []; totNeurons = []; numBins = []; binIDs = [];
totNeurons = size(spikeDataAll,1);
numBins    = size(spikeDataAll(1).zFR,1);
zFRmat     = reshape([spikeDataAll.zFR],[numBins,totNeurons]); zFRmat = zFRmat';
binIDs     = spikeDataAll.binIDs; binIDs = binIDs+0.25;

% first 1E left: sets plot and title as indirect pathway on vs pre
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mouse = []; neuronIdx = []; alpha = []; data2comp = []; sig2comp = [];
decIndex = []; incIndex = []; nsIndex = [];
mouse      = 'a2a'; % still select indirect pathway
neuronIdx  = contains([spikeDataAll.animalName],mouse);
alpha      = 0.05 / sum(neuronIdx);    % correct pval for number of neurons
data2Comp  = [spikeDataAll.diffOnOff]; % now selecting on vs pre
sig2Comp   = [spikeDataAll.pOffOn];    % update pval for on vs pre
titlestr   = 'indirect: on vs. pre';
ylabelstr  = 'z-scored firing rate';
xlabelstr  = 'time (s)';

% gets data index
decIndex   = neuronIdx == 1    & ...
             sig2Comp  < alpha & ...
             data2Comp < 0 ;
incIndex   = neuronIdx == 1    & ...
             sig2Comp  < alpha & ...
             data2Comp > 0 ;    
nsIndex    = neuronIdx == 1    & ...
             sig2Comp  > alpha;   

% plots on vs pre non-sig and sig decreasers         
figure
ax1 = subplot(1,2,1);
errorbarplot_sb(binIDs, zFRmat(nsIndex == 1,:))
hold on
errorbarplot_sb(binIDs, zFRmat(decIndex,:),'r','r')
hold on
lineLength = -10:1:2;
plot(zeros(1,length(lineLength)), lineLength,'k--')
plot(ones(1,length(lineLength))*5, lineLength,'k--')
xlabel(xlabelstr)
ylabel(ylabelstr)
title(titlestr)
set(gca,'TickDir','out')
ylim([-10 2])
xlim([-5 10])
text(-4,-6,{'n=' num2str(sum(nsIndex))})
text(-4,-8,{'n=' num2str(sum(decIndex))},'Color','r')
hold on

% next 1E right : sets plot and title as indirect pathway POST vs pre
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mouse = []; neuronIdx = []; alpha = []; data2comp = []; sig2comp = [];
decIndex = []; incIndex = []; nsIndex = [];
mouse      = 'a2a'; % still select indirect pathway
neuronIdx  = contains([spikeDataAll.animalName],mouse);
alpha      = 0.05 / sum(neuronIdx);    % correct pval for number of neurons
data2Comp  = [spikeDataAll.diffPostPre]; % now selecting post vs pre
sig2Comp   = [spikeDataAll.pPrePost];    % update pval for post vs pre
titlestr   = 'indirect: post vs. pre';
ylabelstr  = 'z-scored firing rate';
xlabelstr  = 'time (s)';

% gets data index
decIndex   = neuronIdx == 1    & ...
             sig2Comp  < alpha & ...
             data2Comp < 0 ;
incIndex   = neuronIdx == 1    & ...
             sig2Comp  < alpha & ...
             data2Comp > 0 ;    
nsIndex    = neuronIdx == 1    & ...
             sig2Comp  > alpha;   

ax2 = subplot(1,2,2);
errorbarplot_sb(binIDs, zFRmat(nsIndex == 1,:))
hold on
errorbarplot_sb(binIDs, zFRmat(decIndex==1,:),'r','r')
hold on
plot(zeros(1,length(lineLength)), lineLength,'k--')
plot(ones(1,length(lineLength))*5, lineLength,'k--')
xlabel(xlabelstr)
ylabel(ylabelstr)
title(titlestr)
set(gca,'TickDir','out')
ylim([-10 2])
xlim([-5 10])
text(-4,-6,{'n=' num2str(sum(nsIndex))})
text(-4,-8,{'n=' num2str(sum(decIndex))},'Color','r')
linkaxes([ax1, ax2],'xy')

%% next 1H left: sets plot and title as direct pathway on vs pre
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mouse = []; neuronIdx = []; alpha = []; data2comp = []; sig2comp = [];
decIndex = []; incIndex = []; nsIndex = [];
mouse      = 'd1'; % still select direct pathway
neuronIdx  = contains([spikeDataAll.animalName],mouse);
alpha      = 0.05 / sum(neuronIdx);    % correct pval for number of neurons
data2Comp  = [spikeDataAll.diffOnOff]; % now selecting on vs pre
sig2Comp   = [spikeDataAll.pOffOn];    % update pval for on vs pre
titlestr   = 'direct: on vs. pre';
ylabelstr  = 'z-scored firing rate';
xlabelstr  = 'time (s)';

% gets data index
decIndex   = neuronIdx == 1    & ...
             sig2Comp  < alpha & ...
             data2Comp < 0 ;
incIndex   = neuronIdx == 1    & ...
             sig2Comp  < alpha & ...
             data2Comp > 0 ;    
nsIndex    = neuronIdx == 1    & ...
             sig2Comp  > alpha;   

% plots on vs pre non-sig and sig decreasers         
figure
ax1 = subplot(1,2,1);
errorbarplot_sb(binIDs, zFRmat(nsIndex == 1,:))
hold on
errorbarplot_sb(binIDs, zFRmat(decIndex,:),'r','r')
hold on
lineLength = -10:1:2;
plot(zeros(1,length(lineLength)), lineLength,'k--')
plot(ones(1,length(lineLength))*5, lineLength,'k--')
xlabel(xlabelstr)
ylabel(ylabelstr)
title(titlestr)
set(gca,'TickDir','out')
ylim([-10 2])
xlim([-5 10])
text(-4,-6,{'n=' num2str(sum(nsIndex))})
text(-4,-8,{'n=' num2str(sum(decIndex))},'Color','r')
hold on

% next 1H right : sets plot and title as direct pathway POST vs pre
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mouse = []; neuronIdx = []; alpha = []; data2comp = []; sig2comp = [];
decIndex = []; incIndex = []; nsIndex = [];
mouse      = 'd1'; % still select direct pathway
neuronIdx  = contains([spikeDataAll.animalName],mouse);
alpha      = 0.05 / sum(neuronIdx);    % correct pval for number of neurons
data2Comp  = [spikeDataAll.diffPostPre]; % now selecting post vs pre
sig2Comp   = [spikeDataAll.pPrePost];    % update pval for post vs pre
titlestr   = 'direct: post vs. pre';
ylabelstr  = 'z-scored firing rate';
xlabelstr  = 'time (s)';

% gets data index
decIndex   = neuronIdx == 1    & ...
             sig2Comp  < alpha & ...
             data2Comp < 0 ;
incIndex   = neuronIdx == 1    & ...
             sig2Comp  < alpha & ...
             data2Comp > 0 ;    
nsIndex    = neuronIdx == 1    & ...
             sig2Comp  > alpha;   

ax2 = subplot(1,2,2);
errorbarplot_sb(binIDs, zFRmat(nsIndex == 1,:))
hold on
errorbarplot_sb(binIDs, zFRmat(decIndex==1,:),'r','r')
hold on
plot(zeros(1,length(lineLength)), lineLength,'k--')
plot(ones(1,length(lineLength))*5, lineLength,'k--')
xlabel(xlabelstr)
ylabel(ylabelstr)
title(titlestr)
set(gca,'TickDir','out')
ylim([-10 2])
xlim([-5 10])
text(-4,-6,{'n=' num2str(sum(nsIndex))})
text(-4,-8,{'n=' num2str(sum(decIndex))},'Color','r')
linkaxes([ax1, ax2],'xy')

