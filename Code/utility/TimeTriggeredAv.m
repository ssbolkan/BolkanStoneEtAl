%USAGE:[Avs, StdErr] = timeTriggeredAv(Trace, trace_ts, sf, timeBefore, timeAfter, triggers_ts, groups)
% Written by Kenneth D. Harris 
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html
% any comments, or if you make any extensions
% let me know at harris@axon.rutgers.edu
% 
%
% ALTERED 7/5/2005 by Joshua Gordon
%   To compute triggered averages based on time, rather than samples.
%   timeBefore & timeAfter given in ms; triggers in timestamps
%
% computes triggered averages from Trace at the times given by T.
% Trace may be 2D, in which case columns are averaged separately.
% The output will then be of the form Avs(Time, Column)
% nBefore and nAfter give the number of samples before and after
% to use.
%
% groups is a group label for the trigger points T.  In this case the
% output will be Avs(Time, Column, Group)
%
% StdErr gives standard error in the same arrangement

function [Avs, StdErr, Save] = timeTriggeredAv(Trace, trace_ts, sf, timeBefore, timeAfter, triggers_ts, groups)

if min(triggers_ts)<min(trace_ts) | max(triggers_ts)>max(trace_ts)
    warning('triggers before or after trace period are being ignored')
    ind = find(triggers_ts<min(trace_ts)| triggers_ts>max(trace_ts));
    triggers_ts(ind) = [];
end

if (nargin<7 | length(groups) == 1)
	groups = ones(length(triggers_ts), 1);
end

% convert before and after to samples
nBefore = ceil(timeBefore*sf/1000);
nAfter = ceil(timeAfter*sf/1000);

% convert triggers_ts to trace indices
triggers = ones(length(triggers_ts),1);
for i = 1:length(triggers_ts)
    subtractions = abs(trace_ts - triggers_ts(i));
    [tmp,triggers(i)] = nanmin(subtractions);
end 

if size(Trace,2)>16  %if trace has more than 16 columns, then
    Trace = Trace';  %change orientation (make each vector column vector)
end
nColumns = size(Trace,2);
nSamples = nBefore + nAfter + 1;
nGroups = nanmax(groups);
maxTime = size(Trace, 1);

Avs = zeros(nSamples, nColumns, nGroups);

BlockSize = floor(2000000/nSamples); % memory saving parameter

for grp = 1:nGroups

	Sum = zeros(nSamples, nColumns);
	SumSq = zeros(nSamples, nColumns);
	MyTriggers = find(groups==grp & triggers > nBefore & triggers <= maxTime-nAfter);
	nTriggers = length(MyTriggers);
	
	% go through triggers in groups of BlockSize to save memory
	for Block = 1:ceil(nTriggers/BlockSize)
		BlockTriggers = MyTriggers(1+(Block-1)*BlockSize:min(Block*BlockSize,nTriggers));
		nBlockTriggers = length(BlockTriggers);
		
		TimeOffsets = repmat(-nBefore:nAfter, nBlockTriggers, 1);
		TimeCenters = repmat(triggers(BlockTriggers), 1, nSamples);
		TimeIndex = TimeOffsets + TimeCenters;
		
		Waves = Trace(TimeIndex,:);
		Waves = reshape(Waves, [nBlockTriggers, nSamples, nColumns]);
		Sum = Sum + reshape(nansum(Waves, 1), [nSamples, nColumns]);
        if nargout>2
            Save{Block} = Waves;
        end
		SumSq = SumSq + reshape(nansum(Waves.^2,1), [nSamples, nColumns]);
	end
	
	Avs(:,:,grp) = Sum/nTriggers;
	StdErr(:,:,grp) = sqrt(SumSq/nTriggers - Avs(:,:,grp).^2) / sqrt(nTriggers);
end

                                                                                                                                                                                                                                                                                                                                                                            