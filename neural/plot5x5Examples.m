

trialTypes = getTrialTypes(expInfo,behavioralData,'late');
nt = length(behavioralData.eventTimes(1).daqTime);

% extract stimulus, choice, feedback, value, and block values
trueStimuli = expInfo.block.events.contrastValues(1:nt);
trialCorrectChoice = expInfo.block.events.correctResponseValues(1:nt);
trueBlocks = expInfo.block.events.highRewardSideValues(1:nt);
trueChoices = expInfo.block.events.responseValues(1:nt);
trueFeedback = double(expInfo.block.events.feedbackValues(1:nt));

% assign the 0% stimuli as either 'left' or 'right' depending on the
% preassigned correct choice (not the mouse's choice)
trueStimuli(trueStimuli == 0) = eps;
trueStimuli(abs(trueStimuli) < .05) = ...
    trueStimuli(abs(trueStimuli) < .05).* trialCorrectChoice(abs(trueStimuli) < .05);
trueSide = trueStimuli > 0;

% low rewards are possible on sign-mismatched block and stimulus
% high rewards are possible on sign-matched block and stimulus
% 1 = high, 0 = low
trueValue(trueBlocks.*sign(trueStimuli) == -1) = 0;
trueValue(trueBlocks.*sign(trueStimuli) == 1) = 1;

[~, whichTrials] = selectCondition(expInfo, getUniqueContrasts(expInfo), behavioralData, ...
        initTrialConditions('movementTime','late'));

trialTypes.singleVar.value{1} = intersect(whichTrials,find(trueValue == 1));
trialTypes.singleVar.value{2} = intersect(whichTrials,find(trueValue == 0));

trueFeedback = double(expInfo.block.events.feedbackValues(1:nt));
trialTypes.singleVar.value{1} = intersect(whichTrials,find(trueFeedback == 1));
trialTypes.singleVar.value{2} = intersect(whichTrials,find(trueFeedback == 0));

%%

figure;
cellList = [768 661 520 623 50];
cellList = [390 1827 520 623 1617];

yl = [.81 .6 .6 1.3 .6]

variables = {'contrast' 'direction' 'outcome' 'block' 'value'};
colors{1} = [0 0 .5; 0 0 1; 0 .4 1; .6 .8 1; .75 .75 .75; .8 .45 .45; 1 0 0; .5 0 0; .25 0 0];
colors{2} = [.5 0 1; 1 0 .5];
colors{3} = [.1 .7 .1; .75 0 0];
colors{4} = [0.1 .7 .1; 1 .6 0];
colors{5} = [0.8 0 1; 0 0.6 0.1];
for c = 1:length(cellList)
    for v = 1:length(variables)
        subplot(length(cellList),length(variables),(c-1)*length(cellList)+v)
        conditions = trialTypes.singleVar.(matlab.lang.makeValidName(variables{v}));
        if v == 5
            nc = 2;
        else
            nc = length(conditions);
        end
        hold on
        line([0 0],[-.1 1.5],'Color',[.5 .5 .5],'LineStyle',':')
        line([0.8 0.8],[-.1 1.5],'Color',[.5 .5 .5],'LineStyle',':')
        for n = 1:nc
            plot(neuralData.eta.eventWindow,mean(neuralData.eta.alignedResps{1}(conditions{n},:,cellList(c))),...
            'color',colors{v}(n,:),...
            'linewidth',2)
        end
        xlim([-.55 2])
        ylim([-.05 yl(c)])
        xticks([0 .8])
        set(gca, 'XTickLabels', {'stim' 'go'})
        prettyPlot(gca)
    end
end


%%

imagesc(nan(length(isort1), 2))

cellList = [623 1392 363 153 781];

%%
for c = 1:length(cellList)
    cellno = cellList(c);
    text(1,1958-find(isort1==cellno),strcat({'cell '},num2str(cellno)),'Color','y','FontSize',8,'horizontalalignment','right')
    line([1 2],[1958-find(isort1==cellno) 1958-find(isort1==cellno)],'Color','y','LineStyle','-')
end