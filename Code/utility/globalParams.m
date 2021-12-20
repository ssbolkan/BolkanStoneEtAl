classdef globalParams
    
    properties (Constant)        
        % paths
        dataPath        = 'Z:\Paper_Code\Bolkan_Stone_NN_2021\Data';
        rightCode       = 1
        leftCode        = 0
        nilCode         = -1
        abortCode       = NaN
        
        % analysis parameters
        perfTh          = .6; % default performance cutoff
        psychbins       = -15:5:15; % -15:2:15;%psychometric curve x axis
        
        % plotting 
        stateColors = {[[39 110 167]./255], [[237 177 32]./255], [[237 34 36]./255]};
        laserColor  = [[25 225 25]./255];
      
        towersCl            = 'k';
        visGuideCl          = 'm';
        memGuideCl          = 'c';
                      
                        
    end
    

    
    
end