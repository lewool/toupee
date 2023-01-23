for m = 1:length(mouseList)
    
    mouseName = char(mouseList{m});
    expDate = char(expList{m}{1});
    expNum = expList{m}{2};
    hemisphere = hemList(m);
    [expRef, ~] = data.constructExpRef(mouseName,expDate,expNum);
    fprintf('Loading %s...\n',expRef)
    [behavioralData, expInfo] = data.loadBehavioralDataset(mouseName, expDate, expNum);
    

%%

    et = behavioralData;
    contrasts = getUniqueContrasts(expInfo);
    nt = length(et.eventTimes(1).daqTime);

    trialCorrectChoice = expInfo.block.events.correctResponseValues(1:nt);
%     trueStimuli = expInfo.block.events.contrastValues(1:nt);
    sidedStimuli = expInfo.block.events.contrastValues(1:nt);
    sidedStimuli(sidedStimuli == 0) = eps;
    sidedStimuli(abs(sidedStimuli) < .05) = ...
        sidedStimuli(abs(sidedStimuli) < .05).* ...
        trialCorrectChoice(abs(sidedStimuli) < .05);
    trueChoices = et.wheelMoves.epochs(5).moveDir;
    trueBlock = expInfo.block.events.highRewardSideValues(1:nt);
    trueFeedback = zeros(1,nt);
    trueFeedback(trueChoices .* trialCorrectChoice > 0) = 1;
    trueFeedback(trueChoices .* trialCorrectChoice < 0) = 0;
    trueValue(sidedStimuli .* trueBlock > 0) = 1;
    trueValue(sidedStimuli .* trueBlock < 0) = 0;
    
    [impTrials, ~] = selectCondition(expInfo, contrasts, et, ...
        initTrialConditions('movementTime','early'));

    maxVels = behavioralData.wheelMoves.epochs(5).peakVel;
    RTs = behavioralData.wheelMoves.epochs(5).onsetTimes - behavioralData.eventTimes(1).daqTime;
    trials = intersect(...
        find((RTs < prctile(RTs,97.5)) & (RTs > prctile(RTs,2.5))),...
        find(~isnan(trueChoices)));
    
    fitData{m,1} = impTrials(trials)';
    fitData{m,2} = maxVels(trials)';
    fitData{m,3} = RTs(trials)';
    fitData{m,4} = sidedStimuli(trials)';
    fitData{m,5} = trueChoices(trials)';
    fitData{m,6} = trueBlock(trials)';
    fitData{m,7} = trueFeedback(trials)';
    fitData{m,8} = trueValue(trials)';


    clearvars -except mouseList expList hemList fitData
    
end

%% see how much each predictor contributes to the model

ml = cellfun(@char,mouseList,'UniformOutput',false);
um = unique(ml);
conditions = {'wheel all', 'rt all', 'wheel imp', 'rt imp', 'wheel pat', 'rt pat'};

addedValue = [];
pValues = [];
for c = 1:length(conditions)
    for u = 1:length(um)
        whichSessions = strcmp(ml,um{u});

        fitArray = cell2mat(fitData(whichSessions,:));
        fitArray(fitArray(:,5) == 0,:) = [];
        fitArray(isnan(fitArray(:,5)),:) = [];
        
        if strcmp(conditions{c},'wheel all')
            x = fitArray(:,4:end);
            y = fitArray(:,2);
        elseif strcmp(conditions{c},'rt all')
            x = fitArray(:,4:end);
            y = fitArray(:,3);
        elseif strcmp(conditions{c},'wheel imp')
            x = fitArray(logical(fitArray(:,1)),4:end);
            y = fitArray(logical(fitArray(:,1)),2);
        elseif strcmp(conditions{c},'rt imp')
            x = fitArray(logical(fitArray(:,1)),4:end);
            y = fitArray(logical(fitArray(:,1)),3);
        elseif strcmp(conditions{c},'wheel pat')
            x = fitArray(logical(~fitArray(:,1)),4:end);
            y = fitArray(logical(~fitArray(:,1)),2);
        elseif strcmp(conditions{c},'rt pat')
            x = fitArray(logical(~fitArray(:,1)),4:end);
            y = fitArray(logical(~fitArray(:,1)),3);
        end

        x = [ones(length(x),1) x];
        x(:,2) = abs(x(:,2));
        y = abs(y);

        nfold = 40;
        cvp = cvpartition(size(x,1),'Kfold',nfold);
        for n = 1:nfold
            c_full = x(cvp.training(n),:)\y(cvp.training(n));
            SS_res = sum((y(cvp.test(n)) - x(cvp.test(n),:)*c_full).^2);
            SS_tot = sum((y(cvp.test(n)) - mean(y(cvp.test(n)))).^2);
            ev_full(n,:) = 1-SS_res/SS_tot;

            for v = 2:size(x,2)
                redc = 1:size(x,2);
                redc(v) = [];
                c_red = x(cvp.training(n),redc)\y(cvp.training(n));
                SS_red = sum((y(cvp.test(n)) - x(cvp.test(n),redc)*c_red).^2);
                ev_red(n,v) = 1-SS_red/SS_tot;
            end
        end
        addedValue(u,:,c) = (mean(ev_full) - mean(ev_red))*100;
        pValues(u,:,c) = anovan(y,{x(:,2) x(:,3) x(:,4) x(:,5) x(:,6)},'display','off');

    end
end
numSig = squeeze(sum(pValues < .05))';

%% plot
figure;
for s = 1:length(conditions)
    mean_av = nanmean(addedValue(:,2:end,s));
    sem_av = nanstd(addedValue(:,2:end,s))/sqrt(size(addedValue,1));
    subplot(3,2,s)
    h1 = bar(mean_av,'FaceColor',[.5 .5 .5]);
    hold on
    h2 = errorbar(mean_av,sem_av,'Color','k','LineStyle','none');
    prettyPlot(gca)
    set(gca, 'XTickLabels', {'contrast' 'choice' 'block' 'outcome' 'value'})
    xlabel('Predictor')
    ylabel('Explained variance (%)')
    ylim([-1 12])
    txtHeight = 1+(h1.YData + h2.YPositiveDelta);
    for i = 1:5
        text(i,txtHeight(i),strcat(num2str(numSig(s,i)),'/',num2str(size(addedValue,1))),'HorizontalAlignment','center')
    end
    
    if s == 1
        title('Wheel velocity (all trials)')
    elseif s == 2
        title('Response time (all trials)')
    elseif s == 3
        title('Wheel velocity (imp. trials)')
    elseif s == 4
        title('Response time (imp. trials)')
    elseif s == 5
        title('Wheel velocity (pat. trials)')
    elseif s == 6
        title('Response time (pat. trials)')
    end
    
end