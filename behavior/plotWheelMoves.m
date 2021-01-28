% clear all
% expInfo = initExpInfo({{'LEW031'}},{{'2020-02-03',1,[1]}});
% matched = 0;
% 
% % expInfo = initExpInfo('LEW031');
% % matched = 1;
% 
% %% load data & extract some variables to make this code still work
% 
% if matched == 1
%     [expInfo, neuralData, behavioralData] = processExperiment(expInfo,'matched');
%     [neuralData] = alignResps(expInfo, neuralData, behavioralData);
%     [neuralData] = getSignificantActivity(expInfo, behavioralData, neuralData);
%     combinedNeuralData = combineNeuralData(expInfo, behavioralData, neuralData,'matched');
% 
%     alignedResps = combinedNeuralData.matched.eta.alignedResps;
%     eventWindow = combinedNeuralData.matched.eta.eventWindow;
%     bfcH = combinedNeuralData.matched.stats.bfcH;
%     pLabels = combinedNeuralData.matched.stats.labels;
%     events = combinedNeuralData.matched.eta.events;
% elseif matched == 0 
%     [expInfo, neuralData, behavioralData] = processExperiment(expInfo);
%     
% end

%% set up trial conditions to compare

clear contrastConditions trialConditions labels condIdx
contrasts = getUniqueContrasts(expInfo);
allContrasts = getAllContrasts(expInfo);

%set up trial conditions for hi-L and hi-R blocks
trialConditions{1} = initTrialConditions('movementDir','cw','movementTime','early');
trialConditions{2} = initTrialConditions('movementDir','ccw','movementTime','early');
trialConditions{3} = initTrialConditions('movementDir','cw','movementTime','late');
trialConditions{4} = initTrialConditions('movementDir','ccw','movementTime','late');

trialLabels{1} = 'leftEarly_';
trialLabels{2} = 'rightEarly_';
trialLabels{3} = 'leftLate_';
trialLabels{4} = 'rightLate_';

contrastConditions{1} = contrasts;
contrastLabels{1} = 'sAll';

testTrials = 1:2:size(behavioralData.eventTimes(1).daqTime,2);
trainTrials = 2:2:size(behavioralData.eventTimes(1).daqTime,2);

d = 1;
for c = 1:length(contrastConditions)
    for t = 1:length(trialConditions)
        [~, condIdx{d,:}.all] = selectCondition(expInfo, contrastConditions{c}, behavioralData, trialConditions{t});
        condIdx{d,:}.test = intersect(testTrials,condIdx{d}.all);
        condIdx{d,:}.train = intersect(trainTrials,condIdx{d}.all);
        labels{d,1} = strcat(trialLabels{t},contrastLabels{c});
        d = d+1;
    end
end

%% all traces

figure;
hold on;

colors = [0 .4 1; 1 0 0];
titles = {'choose left (early)' 'choose right (early)' 'choose left( late)' 'choose right(late)'};

for c = 1:length(condIdx)
    subplot(2,2,c)
    hold on
    if mod(c,2) == 1
        color = colors(1,:);
    else
        color = colors(2,:);
    end
    testIdx = condIdx{c, 1}.all;
    for t = testIdx
        rawTraceTimes = []; rawTraceValues = [];
        rawTraceTimes = behavioralData.wheelMoves.traces.time{t};
        rawTraceValues = behavioralData.wheelMoves.traces.pos{t};

        eventTime = interp1(rawTraceTimes,rawTraceTimes,behavioralData.eventTimes(1).daqTime(t),'nearest','extrap');
        feedbackTime = interp1(rawTraceTimes,rawTraceTimes,behavioralData.eventTimes(5).daqTime(t),'nearest','extrap');
        eventIdx = find(rawTraceTimes == eventTime - 0.5);
        feedbackIdx = find(rawTraceTimes == feedbackTime);
        relativeTraceTimes{t,:} = rawTraceTimes - eventTime;

        relativeTraceValues{t,:} = rawTraceValues;

        p = plot(relativeTraceTimes{t,:}(1:end),relativeTraceValues{t,:}(1:end),'color',color,'linewidth',1);
        p.Color(4) = 0.1;

    end
    
    line([0 0],[-100 100],'LineStyle','--','color',[.5 .5 .5])
    xlim([-0.5 2]);
    ylim([-100 100]);
    xlabel('time from stimulus onset (s)')
    ylabel('wheel position (mm)')
    title(titles{c});
end

%% add mean traces

for iTrial = 1:length(behavioralData.wheelMoves.traces.pos)
    rawTraceTimes = behavioralData.wheelMoves.traces.time{iTrial};
        rawTraceValues = behavioralData.wheelMoves.traces.pos{iTrial};

        eventTime = interp1(rawTraceTimes,rawTraceTimes,behavioralData.eventTimes(1).daqTime(iTrial),'nearest','extrap');
        feedbackTime = interp1(rawTraceTimes,rawTraceTimes,behavioralData.eventTimes(5).daqTime(iTrial),'nearest','extrap');
        eventIdx = find(rawTraceTimes == eventTime - 0.5);
        feedbackIdx = find(rawTraceTimes == feedbackTime);
        traceTime = rawTraceTimes - eventTime;
        traceValue = rawTraceValues;
        

        startIdx = find(traceTime == -.5);
        endIdx = find(traceTime == 1.5);
        timeSeries = eventIdx:eventIdx+2000;
        
        truncValues(iTrial,:) = traceValue(timeSeries);
end
%%

% figure;
% hold on;

colors = [0 .4 1; 1 0 0];
titles = {'choose left (early)' 'choose right (early)' 'choose left (late)' 'choose right (late)'};

for c = 1:length(condIdx)
    subplot(2,2,c)
    hold on
    if mod(c,2) == 1
        color = colors(1,:);
    else
        color = colors(2,:);
    end
    testIdx = condIdx{c, 1}.all;
    
    meanTrace = nanmean(truncValues(testIdx,:))';
    semTrace = nanstd(truncValues(testIdx,:),[],1)';
    window = -.5:0.001:1.5;
    
    plot(window,meanTrace,'color',color,'linewidth',2);
    hold on
%     plot(window,meanTrace + semTrace,'color',color,'linewidth',2);
    plotCI = fill([window';flipud(window')],[(meanTrace-semTrace);flipud(meanTrace+semTrace)],color, 'LineStyle', 'none');
    alpha(0.2);
    
    line([0 0],[-50 50],'LineStyle','--','color',[.5 .5 .5])
    xlim([-0.5 1.5]);
    ylim([-50 50]);
    xlabel('time from stimulus onset (s)')
    ylabel('wheel position (mm)')
    title(titles{c});
end

