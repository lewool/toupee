%% initialize experiment details

expInfo = initExpInfo({{'LEW032'}},{{'2020-02-14',3,[3]}});
matched = 0;

% expInfo = initExpInfo('LEW031');
% matched = 1;

%% load data & extract some variables to make this code still work

if matched == 1
    [expInfo, neuralData, behavioralData] = processExperiment(expInfo,'matched');
    [neuralData] = alignResps(expInfo, neuralData, behavioralData);
    [neuralData] = getSignificantActivity(expInfo, behavioralData, neuralData);
    combinedNeuralData = combineNeuralData(expInfo, behavioralData, neuralData,'matched');

    alignedResps = combinedNeuralData.matched.eta.alignedResps;
    eventWindow = combinedNeuralData.matched.eta.eventWindow;
    bfcH = combinedNeuralData.matched.stats.bfcH;
    pLabels = combinedNeuralData.matched.stats.labels;
    events = combinedNeuralData.matched.eta.events;
elseif matched == 0 
    [expInfo, neuralData, behavioralData] = processExperiment(expInfo);
    [neuralData] = alignResps(expInfo, neuralData, behavioralData);
    [neuralData] = getSignificantActivity(expInfo, behavioralData, neuralData);
    combinedNeuralData = combineNeuralData(expInfo, behavioralData, neuralData);

    alignedResps = neuralData.eta.alignedResps;
    eventWindow = neuralData.eta.eventWindow;
    bfcH = neuralData.stats.bfcH;
    pLabels = neuralData.stats.labels;
    events = neuralData.eta.events;
end

%% choose cells
whichCells = 'rightStim'; %choose from 'labels' array
if strcmp(whichCells, 'all')
    plotCells = 1:size(alignedResps{1},3);
else
    plotCells = find(bfcH(:,strcmp(pLabels,whichCells)) > 0);
end

%% set up trial conditions to compare

clear contrastConditions trialConditions labels condIdx
contrasts = getUniqueContrasts(expInfo);

%set up trial conditions for hi-L and hi-R blocks
trialConditions{1} = initTrialConditions('highRewardSide','left','movementDir','cw','movementTime','late');
trialConditions{2} = initTrialConditions('highRewardSide','left','movementDir','ccw','movementTime','late');
trialConditions{3} = initTrialConditions('highRewardSide','right','movementDir','cw','movementTime','late');
trialConditions{4} = initTrialConditions('highRewardSide','right','movementDir','ccw','movementTime','late');
trialConditions{5} = initTrialConditions('movementDir','cw','movementTime','late');
trialConditions{6} = initTrialConditions('movementDir','ccw','movementTime','late');

trialLabels{1} = 'bL_mL_';
trialLabels{2} = 'bL_mR_';
trialLabels{3} = 'bR_mL_';
trialLabels{4} = 'bR_mR_';
trialLabels{5} = 'bAll_mL_';
trialLabels{6} = 'bAll_mR_';

contrastConditions{1} = contrasts(contrasts<0);
contrastConditions{2} = contrasts(contrasts>0);
contrastConditions{3} = contrasts(contrasts==0);
contrastLabels{1} = 'sL';
contrastLabels{2} = 'sR';
contrastLabels{3} = 's0';

% for o = 1:length(contrasts)
%     contrastConditions{o} = contrasts(o);
%     if contrasts(o) < 0
%         contrastLabels{o} = num2str(contrasts(o));
%     elseif contrasts(o) > 0
%         contrastLabels{o} = strcat('+',num2str(contrasts(o)));
%     end
% end

testTrials = 1:2:size(alignedResps{1},1);
trainTrials = 2:2:size(alignedResps{1},1);

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

%% set up some plot values

condList = {'bAll_mL_sL' 'bAll_mR_sL' 'bAll_mL_s0' 'bAll_mR_s0' 'bAll_mL_sR' 'bAll_mR_sR'};

color_stimL = [0 .4 1];
color_stimR = [1 0 0];
color_stim0 = [0 0 0];
color_moveL = [0 0 1];
color_moveR = [0.5 0 0];

stimTimeColors = [...
    0 .4 1;
    0 .4 1;
    .4 .4 .4;
    .4 .4 .4;
    1 0 0;
    1 0 0 ...
];

moveTimeColors = [...
    0 0 1;
    .5 0 0;
    0 0 1;
    .5 0 0;
    0 0 1;
    .5 0 0 ...
];  

for iCond = 1:length(condList)
    tl(iCond) = size(condIdx{strcmp(labels,condList{iCond})}.all,2);
end
buffer = 5;
total_length = sum(tl)+ buffer*(length(condList));

%choose which events to plot
plotEvents = events(1:2);
k = 1;

%% plot
fig = figure;
set(fig, 'Position', [680 76 876 912]);
hold on;

if ~exist('k') == 1
    k = 1;
end

max_k = length(plotCells);

while k <= max_k

for e = 1:length(plotEvents)
    for iCond = 1:length(condList)
    spidx1 = sub2ind([length(plotEvents) total_length], e, sum(tl(1:iCond))-tl(iCond)+buffer*iCond);
        spidx2 = sub2ind([length(plotEvents) total_length], e, buffer*iCond+sum(tl(1:iCond)));
        subplot(total_length,length(plotEvents),[spidx1 spidx2])
        cla;
    end
end

    for iCond = 1:length(condList)
    whichTrials = condIdx{strcmp(labels,condList{iCond})}.all;
    numTrials = size(whichTrials,2);
    %find time diff between stimOn and movOn
    relativeTimes = et(7).daqTime(whichTrials) - et(1).daqTime(whichTrials);

    %always sort trials from min to max difference between stimOn and movOn,
    %irrespective of ETA alignment
    [~,sortIdx] = sort(relativeTimes,'ascend');

    trialTimes = [...
        et(1).daqTime(whichTrials(sortIdx))',...
        et(7).daqTime(whichTrials(sortIdx))',...
        et(5).daqTime(whichTrials(sortIdx))',...
    ];

    at(:,:,1) = trialTimes - trialTimes(:,1);
    at(:,:,2) = trialTimes - trialTimes(:,2);
    at(:,:,3) = trialTimes - trialTimes(:,3);
    
    for e = 1:length(plotEvents)
        spidx1 = sub2ind([length(plotEvents) total_length], e, sum(tl(1:iCond))-tl(iCond)+buffer*iCond);
        spidx2 = sub2ind([length(plotEvents) total_length], e, buffer*iCond+sum(tl(1:iCond)));
        subplot(total_length,length(plotEvents),[spidx1 spidx2])
        f = imagesc(eventWindow,1:numTrials,alignedResps{e}(whichTrials(sortIdx),:,plotCells(k)));
        colormap(flipud(gray));
        hold on;
        st = plot(at(:,1,e),1:numTrials,'k.');
        set(st, 'MarkerEdgeColor',stimTimeColors(iCond,:));
        mt = plot(at(:,2,e),1:numTrials,'b.');
        set(mt, 'MarkerEdgeColor',moveTimeColors(iCond,:));
%         rt = plot(at(:,3,e),1:numTrials,'k.');
        box off
        set(gca,'ytick',[]);
        
        iMax(iCond) = max(max(f.CData));
        iMin(iCond) = min(min(f.CData));
        
        if e == 1 && iCond == 1
            title('stimulus onset')
        elseif e == 2 && iCond == 1
            title('movement onset')
        elseif e == 3 && iCond == 1
            title('reward onset')
        end
        if e == 1
            text(-2,5,num2str(numTrials))
            text(-3,5,condList{iCond}(6:end),'Interpreter','none')
        end
        if iCond ~=length(condList)
            set(gca,'xtick',[]);
        end
        if iCond == length(condList)
            xlabel('time (s)')
        end
    end
    clear trialTimes at
end

    for e = 1:length(plotEvents)
        for iCond = 1:length(condList)
        spidx1 = sub2ind([length(plotEvents) total_length], e, sum(tl(1:iCond))-tl(iCond)+buffer*iCond);
            spidx2 = sub2ind([length(plotEvents) total_length], e, buffer*iCond+sum(tl(1:iCond)));
            subplot(total_length,length(plotEvents),[spidx1 spidx2])
            caxis([min(iMin) max(iMax)]);
        end
    end

    was_a_key = waitforbuttonpress;
    
    if was_a_key && strcmp(get(fig, 'CurrentKey'), 'leftarrow')
      k = max(1, k - 1);
    elseif was_a_key && strcmp(get(fig, 'CurrentKey'), 'rightarrow')
      k = min(max_k, k + 1);
    end
    
end
