spikeDataAll = % spikeDataAll is the processed, spike-sorted, single-unit 
                % output from 32-ch neuronexus silicon optrode recordings
                % in NpHR-expressing A2a- or D1R-Cre mice receiving 5s laser 
                % (5mW,532nm) sweeps while ambulating on a running wheel
                % each row is a single-unit with the following fields (columns):
    filename       : character array      % local raw data folder/plexon files
    animalName     : cell array of        % animal ID - used for determining genotype, e.g.: 
                     characters           %  <contains(spikeDataAll(nNeuron).animalName,'a2a')> --> indirect
    chanNum        : 1X1 numeric          % silicon optrode channel
    numSpikes      : 1X1 numeric          % total number of spike in recording (~15min)
    waveformPoints : 1X1 numeric          % number of samples for each waveform (acquisition sample rate was 20kHz)
                                          % output of plx_waves_v from OmniPlex
                                          % and MAP Offline SDK Bundle (plexon)
                                          % https://plexon.com/software-downloads/#software-downloads-SDKs
    spikeTS        : numSpikesX1 double   % spike timestamps in seconds
    waveforms1     : numSikesXwaveform-   % waveform of spike on channel 1 
                     Points double        % spike-sorting used neighboring silicon sites
    waveforms2     : numSikesXwaveform-   % waveform of spike on channel 2 
                     Points double        % 
    avgWaveform1   : 1X20 double          % average of waveforms1 20 samples at 20kHz
    avgWaveform2   : 1X20 double          % average of waveforms2 20 samples at 20kHz
    laserTS        : 1XnumLaserSweeps dbl % laser onset in seconds. sweeps were 5s square pulse, followed by 15s ITI
    binIDs         : 1XspikeBins double   % .25s bins used for summing spikes in PSTHs around laserTS (i.e. -5:0.25:10 around laser onset)
    zFR            : 1XspikeBins double   % z-scored firing rate in 0.25s bins -5 to 10s around laser onset
                                          % zFR(nBin) = (firingRate(nBin) - bsMeanFR)./ bsStdFR
                                          % bsMeanFR = mean firing rate -5:0s pre laser
                                          % bsStdFR  = standard deviation in firing rate across trials -5:0s pre laser
    FRhz           : 1XspikeBins double   % binned firing rate in Hz -5:0.25:10s around laser onset 
    bsAvgHzTrials  : 1XnumLaserSweeps dbl % baseline (-5 to 0s of laser onset) average firing rate in Hz across all laser sweeps
    onAvgHzTrials  : 1XnumLaserSweeps dbl % laser on (0 to 5s of laser onset) average firing rate in Hz across all laser sweeps
    rebAvgHzTrials : 1XnumLaserSweeps dbl % post laser (5 to 7s of laser onset) average firing rate in Hz across all laser sweeps
    pPrePost       : 1x1 numeric          % cross laser sweep trials paired two-tailed t-Test pvalue of rebAvgHzTrials and bsAvgHzTrials
    prePostStats   : 1x1 structure        % paired two-tailed t-Test stats for above (tStat, degrees of freedom, sd) 
    pOffOn         : 1x1 numeric          % cross laser sweep trials paired two-tailed t-Test pvalue of bsAvgHzTrials and onAvgHzTrials
    offOnStats     : 1x1 structure        % paired two-tailed t-Test stats for above (tStat, degrees of freedom, sd) 
    bsAverageFR    : 1x1 numeric          % cross laser sweep trial average of bsAvgHzTrials
    onAverageFR    : 1x1 numeric          % cross laser sweep trial average of onAvgHzTrials
    rebAverageFR   : 1x1 numeric          % cross laser sweep trial average of rebAvgHzTrials
    diffOnOff      : 1x1 numeric          % simply onAverageFR-bsAverageFR
                                          % used to ID decreases, increasers, e.g. spikeDataAll(1).diffOnOff<0 --> decreaser
    diffPostPre    : 1x1 numeric          % simply onAverageFR-bsAverageFR
                                          % used to ID decreases, increasers, e.g. spikeDataAll(1).diffOnOff<0 --> decreaser
    modIndex       : 1x1 numeric          % (onMean - preMean) / (onMean + preMean); 
    rebIndex       : 1x1 numeric          % (postMean - preMean) / (postMean + preMean);
                                          % last two fields used for comparison to Owen, Liu, Kreitzer (2019)
    
    
    