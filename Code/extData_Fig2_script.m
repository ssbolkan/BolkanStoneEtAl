%% Script to generate Extended Data Fig 2C, 2F, 2K, 2I, 2L data plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% stereology manually quantified within Leica LASX software
% actual images can be found in cup Scott\confocal\RNAscope
% this script calls count values saved in excel docs for each line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all

prompt   = 'Enter one of following: d1_cre or a2a_cre or d2_cre or d1_d2overlap or d2_d1overlap   ';
plotWhat = input(prompt,'s');

% plotWhat    = 'd1_cre'; 
% change this line to following opts to obtain:
% ExtData 2C - d1_cre
% ExtData 2F - a2a_cre
% ExtData 2I - d2_cre
% ExtData 2K - d1_d2overlap
% ExtData 2L - d2_d1overlap

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pathname = globalParams.dataPath;
cd(pathname)
excelFiles  = dir([pathname,'/*.xlsx']);
colorIdx    = {[.6 .1 .1],[0 .4 .8], [0.5 0.1 0.6]};
catMarker   = {'x','x', 'x'};
xNames      = {'DMS','NAcC','NAcSh'};
xLimits     = [0.5 3.5];
if contains(plotWhat,'overlap')
    yLimits = [0 10];
else
    yLimits = [60 100];
end


switch plotWhat
%%%%%%%%%%%%%%%%%%%% plots ExtData 2C or 2F %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case {'d1_cre','a2a_cre'} % only have specificity for these mice B4:B10
        if contains(plotWhat,'d1')
            nFile = 1;
        else
            nFile = 2;
        end
        filename = excelFiles(nFile).name;
        xlRange = 'B4:B10'; 
        for nSheet = 1:3
            dataArray{nSheet} = xlsread(filename,nSheet,xlRange)    ;
        end
        
        data        = [dataArray{1}; ...
            dataArray{2}; ...
            dataArray{3}];
        xIdx        = [ones(size(dataArray{1}));   ...
            ones(size(dataArray{2}))*2; ...
            ones(size(dataArray{3}))*3];
        figure
        hold on
        plotSpread(data,'categoryIdx',xIdx,'distributionIdx',xIdx,...
            'categoryColors',colorIdx,'categoryMarkers',catMarker,'showMM',4,...
            'distributionColors',colorIdx);
        set(gca,'xtick',1:3,'xticklabel',xNames)
        xtickangle(60)
        set(gca,'TickDir','out')
        box off
        xlim(xLimits)
        ylim(yLimits)
        deleteStr = '_cre'; lineStr = erase(plotWhat,deleteStr);
        titlstr = {'specificity' lineStr ' / GFP' };
        title(titlstr)
        ylabel('% overlap')
                      
%%%%%%%%%%%%%%%%%%%% plots ExtData 2I %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
    case 'd2_cre' % have specificity B4:B10 and penetrance D4:D10
        nFile = 3;
        filename = excelFiles(nFile).name;
        xlRange1 = 'B4:B10'; xlRange2 = 'D4:D10';
        for nSheet = 1:3
            dataArray1{nSheet} = xlsread(filename,nSheet,xlRange1);
            dataArray2{nSheet} = xlsread(filename,nSheet,xlRange2);
        end
        
        data1        = [dataArray1{1}; ...
            dataArray1{2}; ...
            dataArray1{3}];
        xIdx1        = [ones(size(dataArray1{1}));   ...
            ones(size(dataArray1{2}))*2; ...
            ones(size(dataArray1{3}))*3];
        
        data2        = [dataArray2{1}; ...
            dataArray2{2}; ...
            dataArray2{3}];
        xIdx2        = [ones(size(dataArray2{1}));   ...
            ones(size(dataArray2{2}))*2; ...
            ones(size(dataArray2{3}))*3];
        
        figure
        subplot(1,2,1)
        hold on
        plotSpread(data1,'categoryIdx',xIdx1,'distributionIdx',xIdx1,...
            'categoryColors',colorIdx,'categoryMarkers',catMarker,'showMM',4,...
            'distributionColors',colorIdx);
        set(gca,'xtick',1:3,'xticklabel',xNames)
        xtickangle(60)
        set(gca,'TickDir','out')
        box off
        xlim(xLimits)
        ylim(yLimits)
        deleteStr = '_cre'; lineStr = erase(plotWhat,deleteStr);
        titlstr = {'specificity' lineStr ' / Cre' };
        title(titlstr)
        ylabel('% overlap')
        
        subplot(1,2,2)
        hold on
        plotSpread(data2,'categoryIdx',xIdx2,'distributionIdx',xIdx2,...
            'categoryColors',colorIdx,'categoryMarkers',catMarker,'showMM',4,...
            'distributionColors',colorIdx);
        set(gca,'xtick',1:3,'xticklabel',xNames)
        xtickangle(60)
        set(gca,'TickDir','out')
        box off
        xlim(xLimits)
        ylim(yLimits)
        deleteStr = '_cre'; lineStr = erase(plotWhat,deleteStr);
        titlstr = {'penetrance' lineStr ' / Cre' };
        title(titlstr)

%%%%%%%%%%%%%%%%%%%% plots ExtData 2K %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
    case 'd1_d2overlap' % did this in d2-cre and d1-tom mice F4:F10
        nCount = 1;
        for nFile = 3:4
            filename = excelFiles(nFile).name;
            xlRange = 'F4:F10';
            for nSheet = 1:3
                dataToAdd = [];
                dataToAdd =  xlsread(filename,nSheet,xlRange);
                dataArray{nSheet}(nCount:(nCount+numel(dataToAdd)-1),1) = xlsread(filename,nSheet,xlRange);
            end
            nCount = nCount + numel(dataToAdd);
        end
        
        data        = [dataArray{1}; ...
                       dataArray{2}; ...
                       dataArray{3}];
        xIdx        = [ones(size(dataArray{1}));   ...
                       ones(size(dataArray{2}))*2; ...
                       ones(size(dataArray{3}))*3];
        figure
        hold on
        plotSpread(data,'categoryIdx',xIdx,'distributionIdx',xIdx,...
            'categoryColors',colorIdx,'categoryMarkers',catMarker,'showMM',4,...
            'distributionColors',colorIdx);
        set(gca,'xtick',1:3,'xticklabel',xNames)
        xtickangle(60)
        set(gca,'TickDir','out')
        box off
        xlim(xLimits)
        ylim(yLimits)      
        ylabel('% D1R / D2R')
        
        
%%%%%%%%%%%%%%%%%%%% plots ExtData 2L %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                        
    case 'd2_d1overlap' % only did this in d1-tom mice H4:H10
        nFile = 4; 
        filename = excelFiles(nFile).name;
        xlRange = 'H4:H10';
        for nSheet = 1:3
            dataArray{nSheet} = xlsread(filename,nSheet,xlRange);
        end
        
        data        = [dataArray{1}; ...
                       dataArray{2}; ...
                       dataArray{3}];
        xIdx        = [ones(size(dataArray{1}));   ...
                       ones(size(dataArray{2}))*2; ...
                       ones(size(dataArray{3}))*3];
        figure
        hold on
        plotSpread(data,'categoryIdx',xIdx,'distributionIdx',xIdx,...
            'categoryColors',colorIdx,'categoryMarkers',catMarker,'showMM',4,...
            'distributionColors',colorIdx);
        set(gca,'xtick',1:3,'xticklabel',xNames)
        xtickangle(60)
        set(gca,'TickDir','out')
        box off
        xlim(xLimits)
        ylim(yLimits)      
        ylabel('% D2R / D1R')
        
    otherwise
        disp('not an option')
end
