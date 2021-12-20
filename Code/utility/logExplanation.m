lg = % lg is a x-mouse concatenated data structure with many fields that reflect trial by trial task information 
     % each log may have modest differences in the included fields, but
     % below is a summary of nearly all possible (and relevant) fields
               choice: [1 x nTrials single]  % choice that the animal made. 
                                               % 1 right, 0 left, -1 choice timed out, NaN manually aborted trials
             cuePos_L: [1 x nTrials cell]    % position (cm) in maze of left cues (note in AoE cues visible 10cm prior)
             cuePos_R: [1 x nTrials cell]    % position (cm) in maze of right cues 
             currMaze: [1 x nTrials uint8]   % mazeID (see extendedData Fig 15 for maze descriptions across tasks)
         excessTravel: [1 x nTrials single]  % maze is ~330 cm long. excessTravel = (330-ydistanceTravelled)/330. 
                                               % negative values indicate mouse cut corner of maze arm. 
                                               % positive values mean mouse did some travel up and back down maze. 
                                               % bad trials thresholded out at excessTravel>0.1
    firstTrialofBlock: [1 x nTrials logical] % task follows a block structure of different mazes. when this ==1 it is first block of a new maze
                                               % in AoE mice warm up on currMaze==4 with a continuous visual guide, no distractors and no delay period 
                                               % upon performance criterion mice progress to main test maze (e.g. currMaze==10 or 11)
                                               % if performance drops at main test maze, mice get a block of 10 easy trials on currMaze ==7, no distractors, cues don't dissapear after appearance  
          isLaserSess: [1 x nTrials logical] % 1 for all trials if laser session (note - on rare occassions, isLaserSess == 1 (true) yet no laserON trials b/c mice did not complete enough main test maze trials
                 lCue: [1 x nTrials single]  % length of cue region in cm
              lMemory: [1 x nTrials single]  % length of memory region in cm
               lStart: [1 x nTrials single]  % length of startbox in cm                                  
              laserON: [1 x nTrials logical] % 1 if laser on, 0 if laser off
    meanPerfBlockCtrl: [1 x nTrials double]  % mean performance in a trial block (i.e. same maze) during light off trials. 
                                               % in code typically bad trial blocks thresholded out at meanPerfBlockCtrl<0.6               
     meanPerfMainMaze: [1 x nTrials double]  % mean performance in a session only considering main test maze (e.g. for AoE currMaze == 10 or 11)          
              nCues_L: [1 x nTrials single]  % total number of cues on left
              nCues_R: [1 x nTrials single]  % total number of cues on right
        nCues_RminusL: [1 x nTrials single]  % delta (right-left) cues
          nCues_total: [1 x nTrials single]  % total number of all cues
                  pos: [1 x nTrials cell]    % a cell is a 3 column x VirMEn iterations in trial to reward (~60Hz) matrix. 
                                               % column1->x-position (cm); column2->y-position(cm); column3->view angle(degrees)
                                               % negative values are leftward, positive rightward
          rewardScale: [1 x nTrials single]  % size of the reward which is 4ul x rewardScale
           sensorDots: [1 x nTrials cell]    % a cell is a 5 column x VirMEn iterations in trial through reward/ITI when VR world is blank (~60Hz). note this is why more iterations than e.g. pos 
                                               % column3->x-axis dots; column4->y-axis dots; 
                                               % column5-> duration in ms of virmen iteration
                                               % this is the raw ball motion measurement for translating ball->VR movement 
                                               % dots is calibrated for each ball on each rig - reflects detected dots for 1 revolution of 8cm circumference 
             speedCue: [1 x nTrials single]  % average x-y speed during cue region
            speedStem: [1 x nTrials single]  % average x-y speed during maze stem (cue+memory regions)
            stemDispl: [1 x nTrials single]  % x-y displacement during maze stem (cue+memory region)
        stemDisplNorm: [1 x nTrials single]  % above but normalized to maze stem length
                 time: [1 x nTrials cell]    % time elapsed per trial at each virmen iteration
             trialDur: [1 x nTrials single]  % total trial duration including reward/ITI
          trialDurCue: [1 x nTrials single]  % total trial duration in cue region
         trialDurFull: [1 x nTrials single]  % total trial duration to reward/no reward
          trialDurMem: [1 x nTrials single]  % total trial duration in memory region
            trialProb: [1 x nTrials cell]    % related to de-biasing algorithm - probability of drawing right (1,1) or left (2,1) trial 
            trialType: [1 x nTrials single]  % 1 if right is rewarded arm, 0 if left is rewarded arm. 
                                               % on a small percent of trials, nCues_L = nCues_R, and trialType (ie rewarded arm) is randomly selected. 
                                               % mouse receives reward when choice==trialType 
              mouseID: [1 x nTrials double]  % mouse entered ID
           genotypeID: [1 x nTrials double]  % D1R-Cre is 1, D2R-Cre is 2, A2A-cre is 3 
            sessionID: [1 x nTrials double]  % sessionIDs are specific to each mouse's training history, so IDs can repeat across mice
                 date: [1 x nTrials double]  % if even date (right hemisphere laser), odd date (left hemisphere laser); 
                                             % can use mod(lg.date,2)==0 for rightIndex, and mod(lg.date,2)==1 for leftIndex  
           laserEpoch: [1 x nTrials cell]    % 'cue' or 'whole'
               taskID: [1 x nTrials cell]    % 'aoe' or 'nd' or 'pc'
           stateIndex: [1 x nTrials cell]    % most likely glm-hmm state
           state1prob: [1 x nTrials cell]    % postieror probability of glm-hmm state 1
           state2prob: [1 x nTrials cell]    % postieror probability of glm-hmm state 2
           state3prob: [1 x nTrials cell]    % postieror probability of glm-hmm state 3
           
           
               
               