function psthRasterExample_metaFile(unitData, xlims, binsize, subsamp, color, plotWaveform)

if nargin < 2
    xlims  = [-5 10];
    binsize = 0.25;
    subsamp = 1;
end

if nargin < 3
    binsize = 0.25;
    subsamp = 1;
end

if nargin < 4 || isempty(subsamp)
    subsamp = 1;
end

if nargin < 5 || isempty(color)
    color = [0 0 0];
end

if nargin < 6
    plotWaveform = 1;
end

[fileName, ~, ~] = fileparts(unitData.fileName);

ch1wave = unitData.avgWaveform1;
ch2wave = unitData.avgWaveform2;

if max(abs(ch1wave)) > max(abs(ch2wave))
    wave = unitData.waveforms1;
else
    wave = unitData.waveforms2;
end

laserTimeStamp = [];
laserTimeStamp = unitData.laserTS;
spikes = [];
spikes = unitData.spikeTS;

figure
h=zeros(1,2);
h(1)=subplot(2,1,1);
binnedspikedata = [];
Avs = []; StdErr = [];
[binnedspikedata]=histc(spikes, min(spikes):binsize:max(spikes));
[Avs, StdErr] = TimeTriggeredAv(binnedspikedata, min(spikes):binsize:max(spikes), 1./binsize, abs(xlims(1))*1e3, xlims(2)*1e3,laserTimeStamp);
b=bar([xlims(1):binsize:xlims(2)],Avs/binsize, 'histc'); hold on
b.EdgeColor   = 'none';
b.FaceColor   = color;
ylim_curr     = get(gca,'ylim');
xtick_curr    = get(gca,'xtick');
xlim_curr     = get(gca,'xlim');
xlim_curr(1)  = xlims(1)-2;
xlim_curr(2)  = xlims(2)+2;
set(gca,'xlim',xlim_curr)
plot(zeros(1,numel(ylim_curr(1):ylim_curr(2))),ylim_curr(1):ylim_curr(2),'k:')
plot(ones(1,numel(ylim_curr(1):ylim_curr(2)))*5,ylim_curr(1):ylim_curr(2),'k:')
plot(ones(1,numel(ylim_curr(1):ylim_curr(2)))*7,ylim_curr(1):ylim_curr(2),'k:')
box off
set(gca,'TickDir','out')
ylabel('Firing Rate (Hz)')
xlabel('time (s)')
deleteStr = 'Z:\Scott\intan\NpHR DMS\';
titleStr = erase(fileName,deleteStr);
title(titleStr)


h(2)=subplot(2,1,2);
for trialnum=1:length(laserTimeStamp)
    zeroedraster=spikes(1:subsamp:end)-laserTimeStamp(trialnum);
    rasterplot(trialnum,zeroedraster(zeroedraster<xlims(2) & zeroedraster>xlims(1))');
end
box off
set(gca,'TickDir','out')
set(gca,'YTick',[])
ylim_curr = get(gca,'ylim');
xtick_curr = get(gca,'xtick');
hold on
plot(zeros(1,numel(ylim_curr(1):ylim_curr(2))),ylim_curr(1):ylim_curr(2),'k:')
plot(ones(1,numel(ylim_curr(1):ylim_curr(2)))*5,ylim_curr(1):ylim_curr(2),'k:')
plot(ones(1,numel(ylim_curr(1):ylim_curr(2)))*7,ylim_curr(1):ylim_curr(2),'k:')
linkaxes(h,'x')

% figure
% for trialnum= 1:length(laserTimeStamp)
%     zeroedraster=spikes(1:subsamp:end)-laserTimeStamp(trialnum);
%     trialSpikes = []; trialSpikes = zeroedraster(zeroedraster<xlims(2) & zeroedraster>xlims(1))';
%     
%     [trialspikedata]=histc(trialSpikes, xlims(1):binsize:xlims(2));
%     trialFR = trialspikedata./binsize;
%     trialFR = smoothdata(trialFR,2,'movmean',4);
%     plot(xlims(1):binsize:xlims(2),trialFR,'-','Color',[0.4 0.4 0.4]); hold on   
% end
%     plot(xlims(1):binsize:xlims(2),Avs/binsize,'k-','LineWidth',3);
    

if plotWaveform == 1
figure

% laserTimeStampOff = laserTimeStamp +5;
% laserTimeStampPre = laserTimeStamp -5;
% offWaves = wave
% onWaves  = wave
try
randInd   = randsample(1:size(wave,1),100);
randWaves = smoothdata(wave(randInd,:),2,'movmean',3);
meanWaveform = smoothdata(nanmean(wave,1),2,'movmean',3);
plot(1:unitData.waveformPoints, randWaves,'Color', [0.4 0.4 0.4]); hold on
plot(1:unitData.waveformPoints, meanWaveform,'k-','LineWidth',3)

set(gca,'TickDir','out')
xlabel('samples (20kHz)')
ylabel('milliVolt')
title(titleStr)
box off
catch
randWaves = smoothdata(wave,2,'movmean',3);
plot(1:unitData.waveformPoints, randWaves,'Color', [0.4 0.4 0.4])
set(gca,'TickDir','out')
xlabel('samples (20kHz)')
ylabel('milliVolt')
title(titleStr)
box off
end
else
end

end
