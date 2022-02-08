switches = find(diff(expInfo.block.events.highRewardSideValues  ));
for is = 1:length(switches)
    range(is,:) = [switches(is) - 99:switches(is),switches(is) + 1:switches(is)+100];
    rhs(is,:) = expInfo.block.events.highRewardSideValues(range(is,end));
end

for is = 1:length(switches)
    if rhs(is,:) == 1
        earlyHighTrials = selectCondition(expInfo, getUniqueContrasts(expInfo), behavioralData, ...
            initTrialConditions('whichTrials',range(is,:)','movementDir','ccw','movementTime','early'));
        rangeEarly(is,:) = earlyHighTrials(range(is,:));
        allHighTrials = selectCondition(expInfo, getUniqueContrasts(expInfo), behavioralData, ...
            initTrialConditions('whichTrials',range(is,:)','movementDir','ccw','movementTime','all'));
        rangeAll(is,:) = allHighTrials(range(is,:));
    elseif rhs(is,:) == -1
        earlyHighTrials = selectCondition(expInfo, getUniqueContrasts(expInfo), behavioralData, ...
            initTrialConditions('whichTrials',range(is,:)','movementDir','cw','movementTime','early'));
        rangeEarly(is,:) = earlyHighTrials(range(is,:));
        allHighTrials = selectCondition(expInfo, getUniqueContrasts(expInfo), behavioralData, ...
            initTrialConditions('whichTrials',range(is,:)','movementDir','cw','movementTime','all'));
        rangeAll(is,:) = allHighTrials(range(is,:));
    end
end

drange = 1:10:200 ;
for d = 1:length(drange)
    a(d) = sum(sum(rangeEarly(:,drange(d):drange(d)+9)))/sum(sum(rangeAll(:,drange(d):drange(d)+9)));
end
%%


for is = 1:length(switches)

    earlyHighTrials = selectCondition(expInfo, getUniqueContrasts(expInfo), behavioralData, ...
            initTrialConditions('whichTrials',range(is,:)','movementTime','early'));
    rangeEarly(is,:) = earlyHighTrials(range(is,:));
    allHighTrials = selectCondition(expInfo, getUniqueContrasts(expInfo), behavioralData, ...
            initTrialConditions('whichTrials',range(is,:)','movementTime','all'));
        rangeAll(is,:) = allHighTrials(range(is,:));
end
    
drange = 1:10:200 ;
for d = 1:length(drange)
    a(d) = sum(sum(rangeEarly(:,drange(d):drange(d)+9)))/sum(sum(rangeAll(:,drange(d):drange(d)+9)));
end

%%
iX = 52:61;
earlyHighTrials_L = selectCondition(expInfo(iX), getUniqueContrasts(expInfo(iX)), behavioralData(iX), ...
    initTrialConditions('highRewardSide','left','movementDir','cw','movementTime','early','responseType','all'));
earlyHighTrials_R = selectCondition(expInfo(iX), getUniqueContrasts(expInfo(iX)), behavioralData(iX), ...
    initTrialConditions('highRewardSide','right','movementDir','ccw','movementTime','early','responseType','all'));

allHighTrials_L = selectCondition(expInfo(iX), getUniqueContrasts(expInfo(iX)), behavioralData(iX), ...
    initTrialConditions('highRewardSide','left','movementDir','cw','movementTime','all','responseType','all'));
allHighTrials_R = selectCondition(expInfo(iX), getUniqueContrasts(expInfo(iX)), behavioralData(iX), ...
    initTrialConditions('highRewardSide','left','movementDir','cw','movementTime','all','responseType','all'));

earlyHighTrials = earlyHighTrials_R + earlyHighTrials_L;
allHighTrials = allHighTrials_R + allHighTrials_L;

earlyLowTrials_L = selectCondition(expInfo(iX), getUniqueContrasts(expInfo(iX)), behavioralData(iX), ...
    initTrialConditions('highRewardSide','left','movementDir','ccw','movementTime','early','responseType','all'));
earlyLowTrials_R = selectCondition(expInfo(iX), getUniqueContrasts(expInfo(iX)), behavioralData(iX), ...
    initTrialConditions('highRewardSide','right','movementDir','cw','movementTime','early','responseType','all'));
allLowTrials_L = selectCondition(expInfo(iX), getUniqueContrasts(expInfo(iX)), behavioralData(iX), ...
    initTrialConditions('highRewardSide','left','movementDir','ccw','movementTime','all','responseType','all'));
allLowTrials_R = selectCondition(expInfo(iX), getUniqueContrasts(expInfo(iX)), behavioralData(iX), ...
    initTrialConditions('highRewardSide','right','movementDir','cw','movementTime','all','responseType','all'));

earlyLowTrials = earlyLowTrials_R + earlyLowTrials_L;
allLowTrials = allLowTrials_R + allLowTrials_L;

sum(earlyHighTrials)/sum(allHighTrials)
sum(earlyLowTrials)/sum(allLowTrials)

%%
for iX = 1:61
    early_early = selectCondition(expInfo(iX), getUniqueContrasts(expInfo(iX)), behavioralData(iX), ...
    initTrialConditions('preStimMovement','quiescent','movementTime','early','responseType','all','pastMovementTime','early'));

late_early = selectCondition(expInfo(iX), getUniqueContrasts(expInfo(iX)), behavioralData(iX), ...
    initTrialConditions('preStimMovement','quiescent','movementTime','early','responseType','all','pastMovementTime','late'));

early_all = selectCondition(expInfo(iX), getUniqueContrasts(expInfo(iX)), behavioralData(iX), ...
    initTrialConditions('preStimMovement','quiescent','movementTime','all','responseType','all','pastMovementTime','early'));

late_all = selectCondition(expInfo(iX), getUniqueContrasts(expInfo(iX)), behavioralData(iX), ...
    initTrialConditions('preStimMovement','quiescent','movementTime','all','responseType','all','pastMovementTime','late'));


propEarly_followingEarly(iX,1) = nansum(early_early)/nansum(early_all);
propEarly_followingLate(iX,1) = nansum(late_early)/nansum(late_all);
end
%%
animalID = [ones(15,1); 2*ones(15,1); 3*ones(21,1); 4*ones(10,1)];
for iX = 1:61
correct_early = selectCondition(expInfo(iX), getUniqueContrasts(expInfo(iX)), behavioralData(iX), ...
    initTrialConditions('preStimMovement','quiescent','movementTime','early','pastResponseType','correct','pastMovementTime','all'));

incorrect_early = selectCondition(expInfo(iX), getUniqueContrasts(expInfo(iX)), behavioralData(iX), ...
    initTrialConditions('preStimMovement','quiescent','movementTime','early','pastResponseType','incorrect','pastMovementTime','all'));

correct_all = selectCondition(expInfo(iX), getUniqueContrasts(expInfo(iX)), behavioralData(iX), ...
    initTrialConditions('preStimMovement','quiescent','movementTime','all','pastResponseType','correct','pastMovementTime','all'));

incorrect_all = selectCondition(expInfo(iX), getUniqueContrasts(expInfo(iX)), behavioralData(iX), ...
    initTrialConditions('preStimMovement','quiescent','movementTime','all','pastResponseType','incorrect','pastMovementTime','all'));


propEarly_followingCorrect(iX,1) = nansum(correct_early)/nansum(correct_all);
propEarly_followingIncorrect(iX,1) = nansum(incorrect_early)/nansum(incorrect_all);

end

%%

for iX = 1:61
active_early = selectCondition(expInfo(iX), getUniqueContrasts(expInfo(iX)), behavioralData(iX), ...
    initTrialConditions('preStimMovement','active','movementTime','early','specificRTs',[.1 Inf]));

quiesc_early = selectCondition(expInfo(iX), getUniqueContrasts(expInfo(iX)), behavioralData(iX), ...
    initTrialConditions('preStimMovement','quiescent','movementTime','early','specificRTs',[.1 Inf]));

active_all = selectCondition(expInfo(iX), getUniqueContrasts(expInfo(iX)), behavioralData(iX), ...
    initTrialConditions('preStimMovement','active','movementTime','all','specificRTs',[.1 Inf]));

[~, quiesc_all] = selectCondition(expInfo(iX), getUniqueContrasts(expInfo(iX)), behavioralData(iX), ...
    initTrialConditions('preStimMovement','quiescent','movementTime','all','specificRTs',[.1 3]));


propEarly_followingActive(iX,1) = nansum(active_early)/nansum(active_all);
propEarly_followingQuiescent(iX,1) = nansum(quiesc_early)/nansum(quiesc_all);

end

%%
figure;
hold on;
subplot(1,3,1);
for iX = 1:61
    %plot the prop(early) for each session, paired between categories
    p = plot([1, 2],[propEarly_followingCorrect(iX),propEarly_followingIncorrect(iX)], '-', 'LineWidth', 2, 'Color',[.5 .5 .5],'MarkerSize', 10);
    %transparency
    p.Color(4) = .5;
    %don't erase
    hold on
end
%set x axis limits
xlim([.75 2.25])
%turn the box off
box off
%set tick direction
set(gca,'tickdir','out')
%label only 1 and 2 on the x axis...
xticks([1 2])
%...and rename them
set(gca, 'XTickLabels', {'was correct', 'was incorrect'})
%label axes
xlabel('Previous trial')
ylabel('P(impulsive moves)')

subplot(1,3,2);
for iX = 1:61
    p= plot([1, 2],[propEarly_followingEarly(iX),propEarly_followingLate(iX)], '-', 'LineWidth', 2, 'Color',[.5 .5 .5],'MarkerSize', 10);
    p.Color(4) = .5;
    hold on
end
xlim([.75 2.25])
box off
set(gca,'tickdir','out')
xticks([1 2])
set(gca, 'XTickLabels', {'was early', 'was late'})
xlabel('Previous trial')

subplot(1,3,3);
for iX = 1:61
    p= plot([1, 2],[propEarly_followingActive(iX),propEarly_followingQuiescent(iX)], '-', 'LineWidth', 2, 'Color',[.5 .5 .5],'MarkerSize', 10);
    p.Color(4) = .5;
    hold on
end
xlim([.75 2.25])
box off
set(gca,'tickdir','out')
xticks([1 2])
set(gca, 'XTickLabels', {'was active', 'was quiescent'})
xlabel('Prestimulus activity')




