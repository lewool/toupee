%% initialize experiment details

clear all
expInfo = initExpInfo({{'LEW031'}},{{'2020-02-03',1,[1]}});
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

whichCells = 'leftStim'; %choose from 'pLabels' array
if strcmp(whichCells, 'all')
    plotCells = 1:size(alignedResps{1},3);
else
    plotCells = find(bfcH(:,strcmp(pLabels,whichCells)) > 0);
end

%% set up trial conditions to compare

clear contrastConditions trialConditions labels condIdx
contrasts = getUniqueContrasts(expInfo);
allContrasts = getAllContrasts(expInfo);

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

%% equalize contrasts between correct and incorrect trials

%decide whether to use this
balanced = 1;
% balanced = 0;

matchList = {...
    'bAll_mL_sL' , 'bAll_mR_sL';...
    'bAll_mL_sR' , 'bAll_mR_sR';...
    'bL_mL_sL' , 'bL_mR_sL';...
    'bL_mL_sR' , 'bL_mR_sR';...
    'bR_mL_sL' , 'bR_mR_sL';...
    'bR_mL_sR' , 'bR_mR_sR'...
    };

for iM = 1:size(matchList,1)
    mLTrials = condIdx{strcmp(labels,matchList{iM,1})}.all; 
    mRTrials = condIdx{strcmp(labels,matchList{iM,2})}.all; 
    mLContrasts = allContrasts(mLTrials);
    mRContrasts = allContrasts(mRTrials);
    
    subsetmLTrials = [];
    subsetmRTrials = [];
    
    for c = 1:length(contrasts)
        nmLT = sum(mLContrasts == contrasts(c));
        nmRT = sum(mRContrasts == contrasts(c));
        minShared = min([nmLT nmRT]);
        subsetmLTrials = [subsetmLTrials randsample(mLTrials(mLContrasts == contrasts(c)),minShared)];
        subsetmRTrials = [subsetmRTrials randsample(mRTrials(mRContrasts == contrasts(c)),minShared)];
    end
    
    whichMatchTrials{iM,1} = subsetmLTrials;
    whichMatchTrials{iM,2} = subsetmRTrials;
end

%% set up some plot values

et = behavioralData.eventTimes;

%which conditions from 'condIdx'/'labels' you want to plot (each gets a subplot)
condList = {...
    'bL_mL_sL' ... 
    'bL_mR_sL' ...
    'bL_mL_s0' ...
    'bL_mR_s0' ...
    'bL_mL_sR' ...
    'bL_mR_sR' ...
    'bR_mL_sL' ... 
    'bR_mR_sL' ...
    'bR_mL_s0' ...
    'bR_mR_s0' ...
    'bR_mL_sR' ...
    'bR_mR_sR' ...
};

% which colors you want to use (length should = condList)
stimTimeColors = [...
    0 .4 1;
    0 .4 1;
    .4 .4 .4;
    .4 .4 .4;
    1 0 0;
    1 0 0;
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
    .5 0 0;
    0 0 1;
    .5 0 0;
    0 0 1;
    .5 0 0;
    0 0 1;
    .5 0 0 ...
];  

%compute the total size of the figure based on trials
for iCond = 1:length(condList)
    if balanced && sum(sum(strcmp(condList(iCond),matchList)))
         tl(iCond) = size(whichMatchTrials{strcmp(condList(iCond),matchList)},2);
    else
        tl(iCond) = size(condIdx{strcmp(labels,condList{iCond})}.all,2);
    end
end
psth = 40;
buffer = 5;
total_length = sum(tl)+ buffer*(length(condList)+1) + psth;

%choose which events to plot
plotEvents = events(1:2);

%reset plotter to the beginning
k = 1;

%% plot (all trials)
fig = figure;
set(fig, 'Position', [680 76 876 912]);
hold on;

if ~exist('k') == 1
    k = 1;
end

max_k = length(plotCells);

while k <= max_k

    %clear subplots
    for e = 1:length(plotEvents)
        for iCond = 1:length(condList)
            spidx1 = sub2ind([length(plotEvents) total_length], e, sum(tl(1:iCond))-tl(iCond)+buffer*iCond);
            spidx2 = sub2ind([length(plotEvents) total_length], e, buffer*iCond+sum(tl(1:iCond)));
            subplot(total_length,length(plotEvents),[spidx1 spidx2])
            cla;
        end
        spidxA = sub2ind([length(plotEvents) total_length], e, total_length-psth);
        spidxB = sub2ind([length(plotEvents) total_length], e, total_length);    
        subplot(total_length,length(plotEvents),[spidxA spidxB])
        cla;
    end
    
    for e = 1:length(plotEvents)
        
        % for each condition in 'condList' (above)
        for iCond = 1:length(condList)
            
            % extract the relevant trials for that condition
            if balanced && sum(sum(strcmp(condList(iCond),matchList)))
                whichTrials = whichMatchTrials{strcmp(condList(iCond),matchList)};
            else
                whichTrials = condIdx{strcmp(labels,condList{iCond})}.all;
            end
            numTrials = size(whichTrials,2);

            %find time diff between stimOn and movOn, sort trials by this
            %difference
            timeDiffs = et(7).daqTime(whichTrials) - et(1).daqTime(whichTrials);
            [~,sortIdx] = sort(timeDiffs,'ascend');

            %record stimOn, movOn, rewardOn times per trial
            trialTimes = [...
                et(1).daqTime(whichTrials(sortIdx))',...
                et(7).daqTime(whichTrials(sortIdx))',...
                et(5).daqTime(whichTrials(sortIdx))',...
            ];

            %0-align all times to either stimOn, moveOn, or rewardOn
            relativeTimes(:,:,1) = trialTimes - trialTimes(:,1);
            relativeTimes(:,:,2) = trialTimes - trialTimes(:,2);
            relativeTimes(:,:,3) = trialTimes - trialTimes(:,3);

            spidx1 = sub2ind([length(plotEvents) total_length], e, sum(tl(1:iCond))-tl(iCond)+buffer*iCond);
            spidx2 = sub2ind([length(plotEvents) total_length], e, buffer*iCond+sum(tl(1:iCond)));
            subplot(total_length,length(plotEvents),[spidx1 spidx2])
            f = imagesc(eventWindow,1:numTrials,alignedResps{e}(whichTrials(sortIdx),:,plotCells(k)));
            colormap(flipud(gray));
            hold on;
            st = plot(relativeTimes(:,1,e),1:numTrials,'k.');
            set(st, 'MarkerEdgeColor',stimTimeColors(iCond,:));
            mt = plot(relativeTimes(:,2,e),1:numTrials,'b.');
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
                text(-.95,5,num2str(numTrials))
                text(-1.95,5,condList{iCond}(6:end),'Interpreter','none')
                xlim([-1 2]);
            elseif e == 2
                xlim([-2 1]);
            end
            
            set(gca,'xtick',[]);
            clear trialTimes relativeTimes
            %psth plots
        
            psth_mean(iCond,:,e) = mean(alignedResps{e}(whichTrials(sortIdx),:,plotCells(k)));

            spidxA = sub2ind([length(plotEvents) total_length], e, total_length-psth);
            spidxB = sub2ind([length(plotEvents) total_length], e, total_length);    
            subplot(total_length,length(plotEvents),[spidxA spidxB])
            plotPSTH = plot(eventWindow, psth_mean(iCond,:,e));
            set(plotPSTH, 'LineWidth',2,'Color',stimTimeColors(iCond,:));
            if mod(iCond,2) == 0 
                set(plotPSTH, 'LineStyle','--')
            else
                set(plotPSTH, 'LineStyle','-')
            end
            hold on;
        
        end
     end   
    yMin = min(min(min(psth_mean)));
    yMax = max([.1 max(max(max(psth_mean)))]);
    
    for e = 1:length(plotEvents)
        for iCond = 1:length(condList)
            spidx1 = sub2ind([length(plotEvents) total_length], e, sum(tl(1:iCond))-tl(iCond)+buffer*iCond);
            spidx2 = sub2ind([length(plotEvents) total_length], e, buffer*iCond+sum(tl(1:iCond)));
            subplot(total_length,length(plotEvents),[spidx1 spidx2])
            caxis([min(iMin) max(iMax)]);
        end
        spidxA = sub2ind([length(plotEvents) total_length], e, total_length-psth);
        spidxB = sub2ind([length(plotEvents) total_length], e, total_length);    
        subplot(total_length,length(plotEvents),[spidxA spidxB])
        box off
        xlabel('time (s)')
        if e == 1
            xlim([-1 2])
            ylim([yMin yMax]);
            ylabel('z-scored activity')
        elseif e == 2
            xlim([-2 1])
            ylim([yMin yMax]);
        end
    end

    was_a_key = waitforbuttonpress;
    
    if was_a_key && strcmp(get(fig, 'CurrentKey'), 'leftarrow')
      k = max(1, k - 1);
    elseif was_a_key && strcmp(get(fig, 'CurrentKey'), 'rightarrow')
      k = min(max_k, k + 1);
    end
    
end
    