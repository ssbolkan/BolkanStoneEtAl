function handle = plotDecodeChoiceXTraject(predChoice,groupOrder)

plotOffsets = [-5 0 5];
binIDs      = predChoice{1,1}.binIDs + 11.5;
counter2    = 1;
xNames      = {[25:25:300]};
figure

for ff = 1:numel(groupOrder)
    accOff = [];         
    for nMouse = 1:size(predChoice{1,ff}.mouse,2)
        accOff(nMouse,:) = predChoice{1,ff}.mouse(nMouse).accuracyOFF * 100;
    end
        
    if contains(groupOrder{ff},'_aoe')
        indivColor     = {[0.3 0.3 0.3]};    
        groupColor     = 'k';  
    elseif contains(groupOrder{ff},'_nd')
        indivColor     = {[0.3 0.6 0.6]};   
        groupColor     = 'm';     
    elseif contains(groupOrder{ff},'_pc')
        indivColor     = {[0.6 0.3 0.6]};    
        groupColor     = 'c';            
    else
    end
    
    handle = subplot(1,2,2);
    x       = []; x       = binIDs + plotOffsets(ff);
    y       = []; y       = nanmean(accOff,1);
    ySEM    = []; ySEM    = nanstd(accOff,0,1)/(sqrt(size(accOff,1)-1));
    low     = []; low     = -ySEM;
    hi      = []; hi      = ySEM;
    errorbar(y, x, low, hi, '.','horizontal', 'color', groupColor, 'markersize',4, 'markerfacecolor', groupColor); 
    hold on
    box off
    set(gca,'TickDir','out');
    ylabel('y position (cm)');
    xlabel('% decode accuracy')    
    title('group')
    xlim([40 100])
       
    data      = [];
    data      = reshape(accOff,size(accOff,1)*size(accOff,2),1);
    colorIdx  = [];
    colorIdx  = repmat(indivColor,1, numel(binIDs));
    counter = 1;
    xIdx    = [];
    for nBin = 1:numel(binIDs)
        xIdx(counter:counter+size(accOff,1)-1,1) = ones(size(accOff,1),1)*(binIDs(nBin)+ plotOffsets(ff));
        counter = counter+ size(accOff,1);
    end
    
    handle       = subplot(1,2,1);
    plotSpread(data,'categoryIdx',xIdx,'distributionIdx',xIdx,'xNames',xNames,'categoryColors',colorIdx,'xyOri','flipped'); hold on     
    box off
    ylabel('y position (cm)');
    xlabel('% decode accuracy')    
    title('individual')
    xlim([40 100])
    
    % hack to make axes match
    if ff == 3
        xIdx        = xIdx+6.5; 
        data        = NaN(numel(data),1);
        plotSpread(data,'categoryIdx',xIdx,'distributionIdx',xIdx,'xNames',xNames,'categoryColors',colorIdx,'xyOri','flipped'); hold on
        box off
        set(gca,'TickDir','out');
        xlim([40 100])
        xticks([20 40 60 80 100])
        ylim([0 300])
        ylabel('y position (cm)')
        xlabel('% decode accuracy')  
    end
   
    % used for stats below - not implemented in useable function form ATM
%     offDataTask(counter2:(counter2+size(accOff,1))-1,1:numel(binIDs)) = accOff;
%     offGroup(counter2:(counter2+size(accOff,1))-1,1) = ones(size(accOff,1),1)*ff;
%     counter2 = counter2 + size(accOff,1);
%     groupedData{ff} = accOff;
end

% for getting stats and reporting on fig - not in function useable form
% for nBins = 1:size(offDataTask,2)
%     datavec  = []; datavec = offDataTask(:,nBins);
%     groupvec = offGroup;
%     [p_task(nBins),tableOut_task{nBins},ANOVA_task{nBins}] ...
%         = anovan(datavec,{groupvec},'varnames',{'group'},'display','off');
%     if p_task(nBins) < 0.05
%         if lillietest(datavec(groupvec==1)) || lillietest(datavec(groupvec==3))
%             [tow_pc_p,~,tow_pc_stats]   = ranksum(datavec(groupvec==1),datavec(groupvec==3));
%         else
%             [~,tow_pc_p,~,tow_pc_stats]    = ttest2(datavec(groupvec==1),datavec(groupvec==3));
%         end
%         if tow_pc_p < 0.05/3; ax1 =subplot(1,2,1); plot(0.4,(binIDs(nBins)+gapFill)-1,'c*'); else; end
%         if lillietest(datavec(groupvec==1)) || lillietest(datavec(groupvec==2))
%             [tow_mg_p,~,tow_mg_stats]   = ranksum(datavec(groupvec==1),datavec(groupvec==2));
%         else
%             [~,tow_mg_p,~,tow_mg_stats]    = ttest2(datavec(groupvec==1),datavec(groupvec==2));
%         end
%         if tow_mg_p < 0.05/3; ax1 =subplot(1,2,1); plot(0.425,(binIDs(nBins)+gapFill)-1,'m*'); else; end
%         if lillietest(datavec(groupvec==2)) || lillietest(datavec(groupvec==3))
%             [pc_mg_p,~,pc_mg_stats]   = ranksum(datavec(groupvec==2),datavec(groupvec==3));
%         else
%             [~,pc_mg_p,~,pc_mg_stats]    = ttest2(datavec(groupvec==2),datavec(groupvec==3));
%         end
%         if pc_mg_p < 0.05/3; ax1 =subplot(1,2,1); plot(0.45,(binIDs(nBins)+gapFill)-1,'k*'); else; end
%     else
%     end
% end
