function [speedXY, speedXY_trials] = xMouseSpeedXY(lg,mazeSection, XorY)

% average per-mouse speed in x or y dimensions
% lg is flattened behavior log
% mazeSection is range for calculation in cm, default [250 300]

if nargin < 2 || isempty(mazeSection)
    mazeSection = [250 300];
end
if nargin < 3 || isempty(XorY)
    XorY = 1:2; 
end

mice        = unique(lg.mouseID);
speedXY     = nan(numel(mice), 1);

for iMouse = 1:numel(mice)
  [pos,time]    = selectMouseTrials(lg, mice(iMouse), 'pos', 'time');
  ntrials       = numel(pos);
  veloTrials    = nan(ntrials,1);
  
  for itrial = 1:ntrials
    startidx    = find(pos{itrial}(:,2) >= mazeSection(1),   1, 'first');
    endidx      = find(pos{itrial}(:,2) <  mazeSection(2), 1, 'last');
    displ       = diff(pos{itrial}(startidx:endidx,XorY));
    
% Old method -- summing total displacement across all iterations in bin    
%     if numel(XorY) > 1
%         displ   = sum(sqrt(sum(displ.^2,2)));
%     else
%         displ   = sqrt(sum(displ.^2));
%     end
    
% calculate speed at each virmen iteration in a spatial bin and then average across
% iterations
    elt         = time{itrial}(startidx:endidx);
    speed       = displ ./ diff(elt);    
    speed(isinf(speed)) = []; 
    veloTrials(itrial) = nanmean(speed);

  end
  
% [n,edges,bin] = histcounts(veloTrials, min(veloTrials):1:max(veloTrials));
% figure
% histogram(n,edges)

  speedXY_trials{iMouse} = veloTrials;
  speedXY(iMouse) = nanmean(veloTrials);

end