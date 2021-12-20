function threshold = getThreshold(prediction,cut)

threshold = min(prediction) + (max(prediction)-min(prediction))/2;

% if cut == 0
%   if abs(threshold - cut) > .1;     threshold = cut; end
% else
%   if abs(threshold - cut) > .1*cut; threshold = cut; end
% end

end