function avgViewAngle = sampleViewAngleVsY_average(pos,minMaxPos,blockAvg)

if nargin < 2 || isempty(minMaxPos)
    minMaxPos = [1 300];
end

if nargin < 3 || isempty(blockAvg)
    blockAvg  = [25];
end

viewAngle     = sampleViewAngleVsY(pos, minMaxPos(1):1:minMaxPos(end)); % sampleViewAngleVsY(lg.pos, posBins,1);
viewAngle     = sampleViewAngleVsY(pos, minMaxPos(1):1:minMaxPos(end)); % sampleViewAngleVsY(lg.pos, posBins,1);
sz            = size(viewAngle);
binAvg        = blockAvg(1);
n             = size(viewAngle, 1);             % Length of first dimension
nc            = n - mod(n, binAvg);             % Multiple of binAvg
np            = nc / binAvg;                    % Length of result / number of trials
xx            = reshape(viewAngle(1:nc, :), binAvg, np, []); % [binAvg x np(trials) x size(viewAngles,2)]
avgViewAngle  = sum(xx, 1) / binAvg;            % Mean over 1st dim
avgViewAngle  = reshape(avgViewAngle, np, []);

end