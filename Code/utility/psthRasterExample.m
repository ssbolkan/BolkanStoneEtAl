function psthRasterExample(fileName, xlims, binsize, subsamp, color, plotWaveform)

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
    plotWaveform = 0;
end

[pathstr, ~, ~] = fileparts(fileName);
cd(pathstr)
load('dataSync.mat', 'dig', 'dig_ts')
% load('data.mat', 'filters')

ch       = 1;
u        = 1;
n = []; npw = []; ts = []; wave = []; wave1 = []; wave2 = [];
[n, npw, ts, wave1] = plx_waves_v(fileName, ch, u);

ch2      = 2;
[~, ~, ~, wave2] = plx_waves_v(fileName, ch2, u);

ch1wave = nanmean(wave1,1);
ch2wave = nanmean(wave2,1);

if max(abs(ch1wave)) > max(abs(ch2wave))
    wave = wave1;
else
    wave = wave2;
end

laserTime      = [];
laserTime      = dig_ts(dig==1);

laserTimeStamp = [];
test = []; test = diff(laserTime)>1;
laserTimeStamp = laserTime(find(diff(laserTime)>1) ) -1; % - (fs-1)
spikes = [];
spikes = ts;

figure
h=zeros(1,2);
h(1)=subplot(2,1,1);
binnedspikedata = [];
Avs = []; StdErr = [];
[binnedspikedata]=histc(spikes, min(spikes):binsize:max(spikes));
[Avs, StdErr] = TimeTriggeredAv(binnedspikedata, min(spikes):binsize:max(spikes), 1./binsize, abs(xlims(1))*1e3, xlims(2)*1e3,laserTimeStamp);
b=bar([xlims(1):binsize:xlims(2)],Avs/binsize, 'histc')
b.EdgeColor = 'none';
b.FaceColor = color;
ylim_curr = get(gca,'ylim');
xtick_curr = get(gca,'xtick');
hold on
plot(zeros(1,numel(ylim_curr(1):ylim_curr(2))),ylim_curr(1):ylim_curr(2),'k:')
plot(ones(1,numel(ylim_curr(1):ylim_curr(2)))*5,ylim_curr(1):ylim_curr(2),'k:')
plot(ones(1,numel(ylim_curr(1):ylim_curr(2)))*7,ylim_curr(1):ylim_curr(2),'k:')

box off
set(gca,'TickDir','out')
ylabel('Firing Rate (Hz)')
xlabel('time (s)')
% deleteStr = 'Z:\Scott\intan\NpHR DMS\';
% titleStr = erase(fileName,deleteStr);
title(fileName)
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


if plotWaveform == 1
    figure
    try
        randInd   = randsample(1:size(wave,1),100);
        randWaves = smoothdata(wave(randInd,:),2,'movmean',3);
        meanWaveform = smoothdata(nanmean(wave,1),2,'movmean',3);
        plot(1:npw, randWaves,'Color', [0.4 0.4 0.4]); hold on
        plot(1:npw, meanWaveform,'k-','LineWidth',3)
        
        set(gca,'TickDir','out')
        xlabel('samples (20kHz)')
        ylabel('milliVolt')
        box off
    catch
        randWaves = smoothdata(wave,2,'movmean',3);
        plot(1:npw, randWaves,'Color', [0.4 0.4 0.4])
        set(gca,'TickDir','out')
        xlabel('samples (20kHz)')
        ylabel('milliVolt')
        box off
    end
else
end

end
