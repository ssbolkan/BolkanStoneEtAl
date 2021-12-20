function lg_ipCon = invertLogs(lg)
lg_ipCon                         = lg;
choice                           = lg.choice;
date                             = lg.date;
trialType                        = lg.trialType;
rightIndex                       = mod(lg.date,2)==0 ;
lg_ipCon.choice(rightIndex)      = ~choice(rightIndex);
lg_ipCon.trialType(rightIndex)   = ~trialType(rightIndex);
 
pos_invert = [];    % make contra positive and ipsi negative
for numTrial = 1:length(choice)
    if mod(lg.date(numTrial),2)==0 && choice(numTrial) == 1 % right hemi, right choice -- ipsi
       pos_invert{numTrial}(:,1) = -lg.pos{numTrial}(:,1); % ipsi should be neg (right is normally positive, so invert)
       pos_invert{numTrial}(:,2) = lg.pos{numTrial}(:,2);
       pos_invert{numTrial}(:,3) = -lg.pos{numTrial}(:,3); 
    elseif mod(lg.date(numTrial),2)==0 && choice(numTrial) == 0 % right hemi, left choice -- contra
       pos_invert{numTrial}(:,1) = -lg.pos{numTrial}(:,1); % contra should be pos (left is normally negative, so invert)
       pos_invert{numTrial}(:,2) = lg.pos{numTrial}(:,2);
       pos_invert{numTrial}(:,3) = -lg.pos{numTrial}(:,3);  
    elseif mod(lg.date(numTrial),2)==1 && choice(numTrial) == 1 % left hemi, right choice -- contra
       pos_invert{numTrial}(:,1) = lg.pos{numTrial}(:,1);  % contra should be pos (right is normally positive, so keep)
       pos_invert{numTrial}(:,2) = lg.pos{numTrial}(:,2);
       pos_invert{numTrial}(:,3) = lg.pos{numTrial}(:,3); 
    elseif mod(lg.date(numTrial),2)==1 && choice(numTrial) == 0 % left hemi, left choice -- ipsi
       pos_invert{numTrial}(:,1) = lg.pos{numTrial}(:,1); % ipsi should be neg (left is normally negative, so keep)
       pos_invert{numTrial}(:,2) = lg.pos{numTrial}(:,2);
       pos_invert{numTrial}(:,3) = lg.pos{numTrial}(:,3);
    else
       pos_invert{numTrial}(:,1) = lg.pos{numTrial}(:,1); % nan and -1 choice
       pos_invert{numTrial}(:,2) = lg.pos{numTrial}(:,2);
       pos_invert{numTrial}(:,3) = lg.pos{numTrial}(:,3); 
    end
end
lg_ipCon.pos                     = pos_invert;

% pos numbers mean more right cues -- on even dates/right hemi, this is
% ipsi, so invert -- on left hemi days this is contra, so keep
nCues_RminusL_invert             = lg.nCues_RminusL; 
nCues_RminusL_invert(rightIndex) = -lg.nCues_RminusL(rightIndex);
lg_ipCon.nCues_RminusL           = nCues_RminusL_invert;

% on even dates/right hemi -- right cues are ipsi (become left ones) and left
% cues are contra (become right ones)
cueDur_L_invert                  = lg.cueDur_L;
cueDur_R_invert                  = lg.cueDur_R;
cueDur_L_invert(rightIndex)      = lg.cueDur_R(rightIndex);
cueDur_R_invert(rightIndex)      = lg.cueDur_L(rightIndex);
lg_ipCon.currDur_L               = cueDur_L_invert;
lg_ipCon.currDur_R               = cueDur_R_invert;

cueOffset_L_invert               = lg.cueOffset_L;
cueOffset_R_invert               = lg.cueOffset_R;
cueOffset_L_invert(rightIndex)   = lg.cueOffset_R(rightIndex);
cueOffset_R_invert(rightIndex)   = lg.cueOffset_L(rightIndex);
lg_ipCon.cueOffset_L             = cueOffset_L_invert;
lg_ipCon.cueOffset_R             = cueOffset_R_invert;

cueOnset_L_invert                = lg.cueOnset_L;
cueOnset_R_invert                = lg.cueOnset_R;
cueOnset_L_invert(rightIndex)    = lg.cueOnset_R(rightIndex);
cueOnset_R_invert(rightIndex)    = lg.cueOnset_L(rightIndex);
lg_ipCon.cueOnset_L              = cueOnset_L_invert;
lg_ipCon.cueOnset_R              = cueOnset_R_invert;

cueOrder_invert                  = lg.cueOrder;
for ncell= 1:numel(cueOrder_invert)
  if rightIndex(ncell)==1
       cueOrder_invert{ncell} = ~lg.cueOrder{ncell};
  else
  end
end
lg_ipCon.cueOrder              = cueOrder_invert;

cuePos_L_invert                = lg.cuePos_L;
cuePos_R_invert                = lg.cuePos_R;
cuePos_L_invert(rightIndex)    = lg.cuePos_R(rightIndex);
cuePos_R_invert(rightIndex)    = lg.cuePos_L(rightIndex);
lg_ipCon.cuePos_L              = cuePos_L_invert;
lg_ipCon.cuePos_R              = cuePos_R_invert;

nCues_L_invert                = lg.nCues_L;
nCues_R_invert                = lg.nCues_R;
nCues_L_invert(rightIndex)    = lg.nCues_R(rightIndex);
nCues_R_invert(rightIndex)    = lg.nCues_L(rightIndex);
lg_ipCon.nCues_L              = nCues_L_invert;
lg_ipCon.nCues_R              = nCues_R_invert;


end