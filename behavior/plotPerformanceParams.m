% deltaBias %difference between left and right blocks at 0% contrast
% meanBias %mean of left and right blocks at 0% contrast
% meanLapse %mean of sum(leftBlockLapses, rightBlockLapses)

colors = lines(5);
for m = 1:length(mouseList)
    
    mouseName = char(mouseList{m});
    expDate = char(expList{m}{1});
    expNum = expList{m}{2};
    hemisphere = hemList(m);
    [expRef, ~] = data.constructExpRef(mouseName,expDate,expNum);
    fprintf('Loading %s...',expRef)
    [behavioralData, expInfo, neuralData] = data.loadDataset(mouseName, expDate, expNum);
    
    if hemisphere < 0
        [~, iBlockTrials0] = selectCondition(expInfo, [-.05 0 .05], behavioralData, ...
            initTrialConditions('highRewardSide','left'));
        [~, cBlockTrials0] = selectCondition(expInfo, [-.05 0 .05], behavioralData, ...
            initTrialConditions('highRewardSide','right'));
        meanContraChoice = mean(expInfo.block.events.responseValues(cBlockTrials0) == 1);
        meanIpsiChoice = mean(expInfo.block.events.responseValues(iBlockTrials0) == 1);
    else
        [~, cBlockTrials0] = selectCondition(expInfo, [-.05 0 .05], behavioralData, ...
            initTrialConditions('highRewardSide','left'));
        [~, iBlockTrials0] = selectCondition(expInfo, [-.05 0 .05], behavioralData, ...
            initTrialConditions('highRewardSide','right'));
        meanContraChoice = mean(expInfo.block.events.responseValues(cBlockTrials0) == -1);
        meanIpsiChoice = mean(expInfo.block.events.responseValues(iBlockTrials0) == -1);
    end
    
    [~, lapseTrials] = selectCondition(expInfo, [-1 1], behavioralData, ...
            initTrialConditions());
    meanLapse(m) = mean(expInfo.block.events.feedbackValues(lapseTrials) == 1);   
    deltaBias(m) = meanContraChoice - meanIpsiChoice;
    meanBias(m) = (meanContraChoice + meanIpsiChoice)/2;
    
    
end

for m = 1:length(mouseList)
    mouseName = char(mouseList{m});
    if contains(mouseName,'031')
        plotColors(m,:) = colors(1,:);
    elseif contains(mouseName,'032')
        plotColors(m,:) = colors(2,:);
    elseif contains(mouseName,'038')
        plotColors(m,:) = colors(3,:);
    elseif contains(mouseName,'046')
        plotColors(m,:) = colors(4,:);
    elseif contains(mouseName,'047')
        plotColors(m,:) = colors(5,:);
    end
end
%%
figure;
hold on;
% line([0 0],[0 1],'LineStyle','--','Color',[.5 .5 .5]);
scatter(deltaBias,meanLapse,60,plotColors,'filled');
xlim([-.03 .6]);
ylim([.5 1]);
set(gca,'tickdir','out')
xlabel('block bias')
ylabel('1 – total lapse')
