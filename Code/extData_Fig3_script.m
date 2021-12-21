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

cd(globalParams.dataPath)
load('runningWheelSpikes.mat')
clearvars -except spikeDataAll
%% plots PSTH and raster for example units in panel ExtData 1C and 1F
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% relies on plexon matlab SDK and compiling mex functions
% Under Home tab --> add ons --> search MATLAB Support for MinGW-w64 C/C++ Compiler 
fileName = [];
for nn = 1:size(spikeDataAll,1)
    fileName{nn,1} = spikeDataAll(nn).fileName;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first ExtData Fig 1A indirect pathway neurons with sig dec during on vs pre
mouse = []; neuronIdx = []; alpha = []; data2comp = []; sig2comp = [];
files2plot = [];  decIndex = [];
mouse      = 'a2a'; % this with next line sets index to indirect pathway
neuronIdx  = contains([spikeDataAll.animalName],mouse);
alpha      = 0.05 / sum(neuronIdx);    % pvalue corrected for # neurons
data2Comp  = [spikeDataAll.diffOnOff]; % data to plot
sig2Comp   = [spikeDataAll.pOffOn];    % p value to use
doSubSamp  = 1; % logical to subsample to reduce size of figures (see note in ExtData Fig1)
suppressLabels = 1; % to see neuron details set this to 0 (although more difficult to see rasters this way)

% gets data index
decIndex  = neuronIdx == 1    & ...
            sig2Comp  < alpha & ...
            data2Comp < 0 ;

files2plot  = find(decIndex==1);
% files2plot  = fileName(decIndex==1);
figure; hold on
xlims       = []; xlims   = [-5 10];
binsize     = []; binsize = 0.25;

% for some reason, neurons ordered differently in fig than in index obtained here 
% 'reOrderedD2' matches to fig. subSampD2 matches subsampling in fig. set
% doSubSamp==0 above to not do subsampling
reOrderedD2 = [18 17 16 14 15 13 12 11 10 9 8 7 6 5 4 3  2 1];
subSampD2   = [3  22  1  1  1  2  1  2  8 5 1 4 1 5 3 15 1 1];
for nn = 1:length(files2plot)
    subplot(11,2,nn);
    if doSubSamp == 1
        rasterExample_metaFile(spikeDataAll(files2plot(reOrderedD2(nn))),xlims,binsize,subSampD2(nn),[],[],suppressLabels);hold on
    else
        rasterExample_metaFile(spikeDataAll(files2plot(reOrderedD2(nn))),xlims,binsize,[],[],[],suppressLabels);
    end
end        


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% next extData Fig 3B direct pathway neurons with sig dec during on vs pre
mouse = []; neuronIdx = []; alpha = []; data2comp = []; sig2comp = [];
files2plot = [];  decIndex = [];
mouse      = 'd1'; % this with next line sets index to indirect pathway
neuronIdx  = contains([spikeDataAll.animalName],mouse);
alpha      = 0.05 / sum(neuronIdx);    % pvalue corrected for # neurons
data2Comp  = [spikeDataAll.diffOnOff]; % data to plot
sig2Comp   = [spikeDataAll.pOffOn];    % p value to use
doSubSamp  = 1; % logical to subsample to reduce size of figures (see note in ExtData Fig1)
suppressLabels = 1; % to see neuron details set this to 0 (although more difficult to see rasters this way)

% gets data index
decIndex  = neuronIdx == 1    & ...
            sig2Comp  < alpha & ...
            data2Comp < 0 ;

files2plot  = find(decIndex==1);        
% files2plot  = fileName(decIndex==1);
figure; hold on
xlims       = []; xlims   = [-5 10];
binsize     = []; binsize = 0.25;


% for some reason, neurons ordered differently in fig than in index obtained here 
% 'reOrderedD1' matches to fig. subSampD1 matches subsampling in fig. set
% doSubSamp==0 above to not do subsampling
reOrderedD1 = [21 20 18 19 17 16 15 14 13 12 11 10 9 8 7 6 5 4 3 2 1];
subSampD1   = [3  1  1  3  1  2  2  3  2  1  4  2  1 3 1 1 1 1 2 2 1];
for nn = 1:length(files2plot)
    subplot(11,2,nn);
    if doSubSamp == 1
        rasterExample_metaFile(spikeDataAll(files2plot(reOrderedD1(nn))),xlims,binsize,subSampD1(nn),[],[],suppressLabels);hold on
    else
        rasterExample_metaFile(spikeDataAll(files2plot(reOrderedD1(nn))),xlims,binsize,[],[],[],suppressLabels);
    end
end

