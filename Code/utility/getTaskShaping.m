function dataOut = getTaskShaping(shapingData, shapingProtocol, removeIncomplete)

if nargin < 3 || isempty(removeIncomplete)
    removeIncomplete = 1;
else
end

% loop just to get max mazeID and max number of shaping/testing sessions
% per mouse
nmice = numel(shapingData);
for jj = 1:nmice
    totSess_xMouse(jj)    = numel(shapingData{jj}.mainMazeID);
    maxMazeID_xMouse(jj)  = max(shapingData{jj}.mainMazeID);
end

% initialze vars to save as NaN matrix 
% (total number of mice x max number of shaping sessions)
perfByMaze  = nan(nmice,max(maxMazeID_xMouse));
daysPerMaze = nan(nmice,max(totSess_xMouse));
daysToMaze  = nan(nmice,max(totSess_xMouse));
currMaze    = nan(nmice,max(totSess_xMouse));
nmazes      = max(maxMazeID_xMouse);

for jj = 1:nmice
    
    % hack to drop mice that did not have an initial session on maze 1-3
    % inidcating that they did not go through full shaping progression
    % (1-3 because some mice can progress to 2 or 3 on first session)
    % otherwise fills currMaze with mainMazeID on all shaping sessions
    if  shapingData{jj}.mainMazeID(1) < 4
        currMaze(jj,1:totSess_xMouse(jj)) = shapingData{jj}.mainMazeID; 
    else
        currMaze(jj,1:totSess_xMouse(jj)) = nan(1,max(totSess_xMouse(jj)));
    end
    
    % sometimes this fails (ignore these mice for now- data remains NaN)
    try
        for mz = 1:nmazes
             
            idx       = [];
            idx       = find(shapingData{jj}.mainMazeID == mz);
            
            % rarely a mouse steps up a maze and later is bumped back
            % this selects only the first block of sessions at mazeID
            diffIndex = find(diff(idx)>1);
            if ~isempty(diffIndex)
                idx = idx(1:diffIndex);
            end
            
            % somtimes a mouse 
            if ~isempty(idx)
                perfByMaze(jj,mz)  = nanmean(shapingData{jj}.perfOverall(idx));
                daysPerMaze(jj,mz) = length(idx);
            else 
                % rarely, mostly for maze 1-3, a mouse can advance across
                % shaping mazes within a session, therefore daysPerMaze
                % is defined as maze, perfByMaze ignored for now
                daysPerMaze(jj,mz) = 1;
            end
            
            % ignore maze 1 as all mice start here
            % and up to final shaping maze sum the number of sessions to
            % reach
            if mz>1 
                daysToMaze(jj,mz)  = sum(daysPerMaze(jj,1:mz-1));
            % for given shaping procedure, ignore multiple sessions after reaching 
            % the first testing maze (7/8 for permCues, 10/11 for AoE 
            elseif mz>6 && strcmp(shapingProtocol, 'PoissonBlocksPermanentCues')
                daysToMaze(jj,mz)  = idx(1);
            elseif mz>9 && strcmp(shapingProtocol, 'PoissonBlocksCondensed3m')
                daysToMaze(jj,mz)  = idx(1);    
            else
            end
            
            
        end
    catch
%         keyboard
    end
end

if removeIncomplete == 1
    % remove test sessions from vars (and mice that didnt learn/have complete shaping data)
    if strcmp(shapingProtocol, 'PoissonBlocksPermanentCues')
        removeInd = maxMazeID_xMouse==7 | maxMazeID_xMouse ==8;
    elseif strcmp(shapingProtocol, 'PoissonBlocksCondensed3m')
        removeInd = maxMazeID_xMouse==10 | maxMazeID_xMouse ==11;
    else
    end
    
    maxMazeID_xMouse(~removeInd) = [];
    totSess_xMouse(~removeInd)   = [];
    daysToMaze(~removeInd',:)    = [];
    perfByMaze(~removeInd',:)    = [];
    daysPerMaze(~removeInd',:)   = [];
    currMaze(~removeInd',:)      = [];
    
    removeInd                    = [];
    removeInd                    = isnan(currMaze(:,1));
    maxMazeID_xMouse(removeInd') = [];
    totSess_xMouse(removeInd')   = [];
    daysToMaze(removeInd,:)      = [];
    perfByMaze(removeInd,:)      = [];
    daysPerMaze(removeInd,:)     = [];
    currMaze(removeInd,:)        = [];
else
end

dataOut.maxMazeID_xMouse     = maxMazeID_xMouse;
dataOut.totSess_xMouse       = totSess_xMouse;
dataOut.daysToMaze           = daysToMaze;
dataOut.perfByMaze           = perfByMaze;
dataOut.daysPerMaze          = daysPerMaze;
dataOut.currMaze             = currMaze;
