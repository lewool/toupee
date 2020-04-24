%% initialize experiment details
clear all;
expInfo = initExpInfo({{'LEW032'}},{{'2020-02-17',2,[2]}});
% expInfo = initExpInfo('LEW031');
matched = 0;
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

%% get cell responses at a particular timepoint

Fs = 0.1;
%%%%%%%%%%%%%%%% compute baseline activity

% align traces to stim onset
event = 'stimulusOnTimes';
stim_alignedTraces = alignedResps{strcmp(events,event)};
stim_eventWindow = eventWindow;

%designate a baseline window
stim_eventIdx = find(stim_eventWindow == 0);
stim_preTime = [-0.5 0] / Fs;
baselineIdx = stim_eventIdx + stim_preTime(1) : stim_eventIdx;

%compute the mean baseline activity per cell, per trial (trials x neurons)
baselineResps = squeeze(mean(stim_alignedTraces(:,baselineIdx,:),2));

%%%% compute peristimulus activity

%designate a peristimulus window
stimTime = [0 0.3] / Fs;
stimIdx = stim_eventIdx + stimTime(1) :stim_eventIdx + stimTime(2);

%compute the mean peristimulus activity per cell, per trial (trials x neurons)
stimResps = squeeze(mean(stim_alignedTraces(:,stimIdx,:),2));

%%%%%%%%%%%%%%%% compute perimovement activity

% align traces to movement onset
event = 'prestimulusQuiescenceEndTimes';
mov_alignedTraces = alignedResps{strcmp(events,event)};
mov_eventWindow = eventWindow;

%designate a movement window
mov_eventIdx = find(mov_eventWindow == 0);
movTime = [-0.2 0.1] / Fs;
movIdx = mov_eventIdx + movTime(1) : mov_eventIdx + movTime(2);

%compute the mean perimovement activity per cell, per trial (trials x neurons)
movResps = squeeze(mean(mov_alignedTraces(:,movIdx,:),2));

%%%%%%%%%%%%%%%% compute premovement activity

%designate a movement window
pmov_eventIdx = find(mov_eventWindow == 0);
pmovTime = [-0.7 -0.1] / Fs;
pmovIdx = pmov_eventIdx + pmovTime(1) : pmov_eventIdx + pmovTime(2);

%compute the mean perimovement activity per cell, per trial (trials x neurons)
pmovResps = squeeze(mean(mov_alignedTraces(:,pmovIdx,:),2));

%%%%%%%%%%%%%%%% compute perireward activity

% align traces to movement onset
event = 'feedbackTimes';
rew_alignedTraces = alignedResps{strcmp(events,event)};
rew_eventWindow = eventWindow;

%designate a movement window
rew_eventIdx = find(rew_eventWindow == 0);
rewTime = [0 0.2] / Fs;
rewIdx = rew_eventIdx + rewTime(1) : rew_eventIdx + rewTime(2);

%compute the mean perireward activity per cell, per trial (trials x neurons)
rewResps = squeeze(mean(rew_alignedTraces(:,rewIdx,:),2));

%% set up trial conditions to compare

clear contrastConditions trialConditions labels condIdx
contrasts = getUniqueContrasts(expInfo);
allContrasts = getAllContrasts(expInfo);

%set up trial conditions for hi-L and hi-R blocks
trialConditions{1} = initTrialConditions('highRewardSide','left','movementDir','cw','movementTime','late');
trialConditions{2} = initTrialConditions('highRewardSide','left','movementDir','ccw','movementTime','late');
trialConditions{3} = initTrialConditions('highRewardSide','right','movementDir','cw','movementTime','late');
trialConditions{4} = initTrialConditions('highRewardSide','right','movementDir','ccw','movementTime','late');
trialConditions{5} = initTrialConditions('highRewardSide','left','movementTime','late');
trialConditions{6} = initTrialConditions('highRewardSide','right','movementTime','late');
trialConditions{7} = initTrialConditions('movementDir','cw','movementTime','late');
trialConditions{8} = initTrialConditions('movementDir','ccw','movementTime','late');

trialLabels{1} = 'bL_mL_';
trialLabels{2} = 'bL_mR_';
trialLabels{3} = 'bR_mL_';
trialLabels{4} = 'bR_mR_';
trialLabels{5} = 'bL_mAll_';
trialLabels{6} = 'bR_mAll_';
trialLabels{7} = 'bAll_mL_';
trialLabels{8} = 'bAll_mR_';

contrastConditions{1} = contrasts(contrasts<0);
contrastConditions{2} = contrasts(contrasts>0);
contrastConditions{3} = contrasts(contrasts~=0);
contrastLabels{1} = 'sL';
contrastLabels{2} = 'sR';
contrastLabels{3} = 'sAll';

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

%% equalize contrasts between trial types

%decide whether to use this
balanced = 1;
% balanced = 0;
if balanced == 1
matchList = {...
    'bAll_mL_sL' , 'bAll_mR_sL';...
    'bAll_mL_sR' , 'bAll_mR_sR';...
    'bL_mL_sL' , 'bR_mL_sL';...
    'bL_mR_sL' , 'bR_mR_sL';...
    'bL_mL_sR' , 'bR_mL_sR';...
    'bL_mR_sR' , 'bR_mR_sR';...
    };

for iM = 1:size(matchList,1)
    mLTrials = condIdx{strcmp(labels,matchList{iM,1})}.all; 
    mRTrials = condIdx{strcmp(labels,matchList{iM,2})}.all; 
    mLContrasts = allContrasts(mLTrials);
    mRContrasts = allContrasts(mRTrials);
    
    subsetLTrials = [];
    subsetRTrials = [];
    
    for c = 1:length(contrasts)
        nmLT = sum(mLContrasts == contrasts(c));
        nmRT = sum(mRContrasts == contrasts(c));
        minShared = min([nmLT nmRT]);
        subsetLTrials = [subsetLTrials randsample(mLTrials(mLContrasts == contrasts(c)),minShared)];
        subsetRTrials = [subsetRTrials randsample(mRTrials(mRContrasts == contrasts(c)),minShared)];
    end
    
    whichMatchTrials{iM,1} = subsetLTrials;
    whichMatchTrials{iM,2} = subsetRTrials;
end
end
%%
close all
figMatrix = figure;
colors = [0.1 0.7 0.1; 1 .6 0; 0 .4 1; 1 0 0];

set(figMatrix,'position',[40 90 1230 900]);
respPeriod = stimResps;

whichCells = 'stim'; %choose from 'pLabels' array
if strcmp(whichCells, 'all')
    plotCells = 1:size(alignedResps{1},3);
else
    plotCells = find(bfcH(:,strcmp(pLabels,whichCells)) > 0);
end

% condVector = {'bR_mL_sL' 'bR_mR_sL' 'bR_mL_sR' 'bR_mR_sR'};
% color = colors(2,:);
% condVector = {'bL_mL_sL' 'bL_mR_sL' 'bL_mL_sR' 'bL_mR_sR'};
% color = colors(1,:);
condVector = {'bAll_mL_sL' 'bAll_mR_sL' 'bAll_mL_sR' 'bAll_mR_sR'};
color = [0 0 0];
mSize = 12;

clear X;
for iCond = 1:length(condVector)
    if balanced && sum(sum(strcmp(condVector(iCond),matchList)))
        whichTrials = whichMatchTrials{strcmp(condVector(iCond),matchList)};
    else
        whichTrials = condIdx{strcmp(labels,condVector{iCond})}.all;
    end
    numTrials = size(whichTrials,2);
    X(:,iCond) = mean(respPeriod(whichTrials,plotCells),1);
end

axLim = 1.05 * [min(min(X)) max(max(X)) min(min(X)) max(max(X))];
for r = 1:4
    x1 = X(:,r);
    for c = 1:4
        x2 = X(:,c);
        subplot(4,4,sub2ind([4 4],r,c));
        hold on;
        if r == c
%             cla;
            h = histogram(x1,linspace(axLim(1),axLim(2),11),'Normalization','probability');
            set(h,'FaceColor',color,...
                'FaceAlpha',.2...
                );
            set(gca,'Color','none');
            box off;
            axis square;
            xlim([axLim(1) axLim(2)]);
        else
            p = scatter(x1,x2,mSize);
            l = line([axLim(1) axLim(2)],[axLim(1) axLim(2)]);
            set(p,...
                'MarkerEdgeColor','none',...
                'Marker','o',...
                'MarkerFaceColor',color,...
                'MarkerFaceAlpha',.2...
                );
            set(l,...
                'LineStyle', '--',...
                'LineWidth',1,...
                'Color',[.5 .5 .5]...
                );
            axis(axLim);
        end
        axis square
        
        if c == 4
            xlabel(condVector{r},'Interpreter','none')
        end
        
        if r == 1
            ylabel(condVector{c},'Interpreter','none')
        end
        
    end
end


condVectorY = {'bR_mL_sL' 'bR_mR_sL' 'bR_mL_sR' 'bR_mR_sR'};
condVectorX = {'bL_mL_sL' 'bL_mR_sL' 'bL_mL_sR' 'bL_mR_sR'};

mSize = 12;
color = [0 0 0];

clear X Y;
for iCond = 1:length(condVectorX)
    if balanced && sum(sum(strcmp(condVectorX(iCond),matchList)))
        whichTrials = whichMatchTrials{strcmp(condVectorX(iCond),matchList)};
    else
        whichTrials = condIdx{strcmp(labels,condVectorX{iCond})}.all;
    end
    numTrials = size(whichTrials,2);
    X(:,iCond) = mean(respPeriod(whichTrials,plotCells),1);
end

for iCond = 1:length(condVectorY)
    if balanced && sum(sum(strcmp(condVectorY(iCond),matchList)))
        whichTrials = whichMatchTrials{strcmp(condVectorY(iCond),matchList)};
    else
        whichTrials = condIdx{strcmp(labels,condVectorY{iCond})}.all;
    end
    numTrials = size(whichTrials,2);
    Y(:,iCond) = mean(respPeriod(whichTrials,plotCells),1);
end

axLim = 1.05 * [min(min([X Y])) max(max([X Y])) min(min([X Y])) max(max([X Y]))];
figVector = figure;
set(figVector,'position',[1280 90 260 900]);

for s = 1:4
    subplot(4,1,s)
    x = X(:,s);
    y = Y(:,s);
    p = scatter(x,y,mSize);
        l = line([axLim(1) axLim(2)],[axLim(1) axLim(2)]);
    title(condVectorX{s}(end-4:end),'Interpreter','none');  
    xlabel('bL')
    ylabel('bR')
    set(p,...
        'MarkerEdgeColor','none',...
        'Marker','o',...
        'MarkerFaceColor',color,...
        'MarkerFaceAlpha',.2...
        );
    set(l,...
        'LineStyle', '--',...
        'LineWidth',1,...
        'Color',[.5 .5 .5]...
        );
    axis(axLim);
    axis square
end

