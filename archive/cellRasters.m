%% load experiment details

% close all;
clear all;

mouseName = 'LEW008';
expDate = '2019-02-07';
expNum = 1;
expSeries = [1];

%% load data

[block, Timeline] = data.loadData(mouseName, expDate, expNum);

%% get event timings and wheel trajectories

signalsNames = {'stimulusOnTimes' 'interactiveOnTimes' 'stimulusOffTimes'};
[eventTimes, wheelTrajectories] = getEventTimes(block, Timeline, signalsNames);

%% load traces
[allFcell, ops] = loadExpTraces(mouseName, expDate, expSeries);

%%
stimEvent = 'stimulusOnTimes';
[stim_alignedTraces, eventWindow] = getExpTraces_noUpsample(mouseName, expDate, expNum, expSeries, allFcell, eventTimes, ops, stimEvent);

movEvent = 'prestimulusQuiescenceEndTimes';
[mov_alignedTraces, eventWindow] = getExpTraces_noUpsample(mouseName, expDate, expNum, expSeries, allFcell, eventTimes, ops, movEvent);

%% select cells with the properties you want
plotAll = chooseCellType(cellLabel, mouseName, expDate, expNum, expSeries, block, allFcell, eventTimes, ops);

%% plot raster browser - equal size subplots, move-aligned
fig = figure;
hold on
set(fig,'Position',[290   100   740   210])
startIdx = find(eventWindow == -1);
endIdx = find(eventWindow == 2);

if ~exist('k') == 1
    k = 1;
end
max_k = length(plotAll);

% select the conditions
contrasts = unique(block.events.contrastValues);

[~, condIdx_highStimL_highL] = selectCondition(block, contrasts(contrasts<0), eventTimes, 'all', 'all', 'all', 'left', 'all', 'all', 'all', 'all', 'all');
[~, condIdx_highStimL_highR] = selectCondition(block, contrasts(contrasts<0), eventTimes, 'all', 'all', 'all', 'right', 'all', 'all', 'all', 'all', 'all');

[~, condIdx_highStimR_highL] = selectCondition(block, contrasts(contrasts>0), eventTimes, 'all', 'all', 'all', 'left', 'all', 'all', 'all', 'all', 'all');
[~, condIdx_highStimR_highR] = selectCondition(block, contrasts(contrasts>0), eventTimes, 'all', 'all', 'all', 'right', 'all', 'all', 'all', 'all', 'all');

highStimLHighL_axes = subplot(1,4,1);
highStimLHighR_axes = subplot(1,4,2);
highStimRHighL_axes = subplot(1,4,3);
highStimRHighR_axes = subplot(1,4,4);

while k <= max_k
    
    subplot(highStimLHighL_axes);cla;
    subplot(highStimLHighR_axes);cla;
    subplot(highStimRHighL_axes);cla;
    subplot(highStimRHighR_axes);cla; 

    subplot(highStimLHighL_axes);
    title(strcat('Left high (',num2str(length(condIdx_highStimL_highL)),')'));
    relativeMovTimes = zeros(length(condIdx_highStimL_highL),1)';
    relativeStimTimes = eventTimes(1).daqTime(condIdx_highStimL_highL)-eventTimes(7).daqTime(condIdx_highStimL_highL);
    [~,sortIdx] = sort(relativeStimTimes,'descend');
    f = imagesc(eventWindow(startIdx:endIdx), 1:length(condIdx_highStimL_highL),mov_alignedTraces{plotAll(k,1)}.eventSpikes(condIdx_highStimL_highL(sortIdx),startIdx:endIdx,plotAll(k,2)));
    hold on;
    plot(relativeStimTimes(sortIdx),1:length(condIdx_highStimL_highL),'g.');
    etaLine = line([0 0],[1 length(condIdx_highStimL_highL)]);
    set(etaLine,'Color',[1 0 1]);
    colormap(flipud(gray));
    xlim([eventWindow(startIdx) eventWindow(endIdx)]);
    box off
    iMax(1) = max(max(f.CData));
    iMin(1) = min(min(f.CData));
    ax = gca;
    ax.TickDir = 'out';
    ax.YTick = [];

    subplot(highStimRHighL_axes);
    title(strcat('Right low (',num2str(length(condIdx_highStimR_highL)),')'));
    relativeMovTimes = zeros(length(condIdx_highStimR_highL),1)';
    relativeStimTimes = eventTimes(1).daqTime(condIdx_highStimR_highL)-eventTimes(7).daqTime(condIdx_highStimR_highL);
    [~,sortIdx] = sort(relativeStimTimes,'descend');
    f = imagesc(eventWindow(startIdx:endIdx), 1:length(condIdx_highStimR_highL),mov_alignedTraces{plotAll(k,1)}.eventSpikes(condIdx_highStimR_highL(sortIdx),startIdx:endIdx,plotAll(k,2)));
    hold on;
    plot(relativeStimTimes(sortIdx),1:length(condIdx_highStimR_highL),'g.');
    etaLine = line([0 0],[1 length(condIdx_highStimR_highL)]);
    set(etaLine,'Color',[1 0 1]);
    colormap(flipud(gray));
    xlim([eventWindow(startIdx) eventWindow(endIdx)]);
    box off
    iMax(2) = max(max(f.CData));
    iMin(2) = min(min(f.CData));
    ax = gca;
    ax.TickDir = 'out';
    ax.YTick = [];

    subplot(highStimLHighR_axes);
    title(strcat('Left low (',num2str(length(condIdx_highStimL_highR)),')'));
    relativeMovTimes = zeros(length(condIdx_highStimL_highR),1)';
    relativeStimTimes = eventTimes(1).daqTime(condIdx_highStimL_highR)-eventTimes(7).daqTime(condIdx_highStimL_highR);
    [~,sortIdx] = sort(relativeStimTimes,'descend');
    f = imagesc(eventWindow(startIdx:endIdx), 1:length(condIdx_highStimL_highR),mov_alignedTraces{plotAll(k,1)}.eventSpikes(condIdx_highStimL_highR(sortIdx),startIdx:endIdx,plotAll(k,2)));
    hold on;
    plot(relativeStimTimes(sortIdx),1:length(condIdx_highStimL_highR),'g.');
    etaLine = line([0 0],[1 length(condIdx_highStimL_highR)]);
    set(etaLine,'Color',[1 0 1]);
    colormap(flipud(gray));
    xlim([eventWindow(startIdx) eventWindow(endIdx)]);
    box off
    iMax(3) = max(max(f.CData));
    iMin(3) = min(min(f.CData));
    ax = gca;
    ax.TickDir = 'out';
    ax.YTick = [];
    
    subplot(highStimRHighR_axes);
    title(strcat('Right high (',num2str(length(condIdx_highStimR_highR)),')'));
    relativeMovTimes = zeros(length(condIdx_highStimR_highR),1)';
    relativeStimTimes = eventTimes(1).daqTime(condIdx_highStimR_highR)-eventTimes(7).daqTime(condIdx_highStimR_highR);
    [~,sortIdx] = sort(relativeStimTimes,'descend');
    f = imagesc(eventWindow(startIdx:endIdx), 1:length(condIdx_highStimR_highR),mov_alignedTraces{plotAll(k,1)}.eventSpikes(condIdx_highStimR_highR(sortIdx),startIdx:endIdx,plotAll(k,2)));
    hold on;
    plot(relativeStimTimes(sortIdx),1:length(condIdx_highStimR_highR),'g.');
    etaLine = line([0 0],[1 length(condIdx_highStimR_highR)]);
    set(etaLine,'Color',[1 0 1]);
    colormap(flipud(gray));
    xlim([eventWindow(startIdx) eventWindow(endIdx)]);
    box off
    iMax(4) = max(max(f.CData));
    iMin(4) = min(min(f.CData));
    ax = gca;
    ax.TickDir = 'out';
    ax.YTick = [];
    
    subplot(highStimLHighL_axes);
    caxis([min(iMin) max(iMax)]);
    subplot(highStimLHighR_axes);
    caxis([min(iMin) max(iMax)]);
    subplot(highStimRHighL_axes);
    caxis([min(iMin) max(iMax)]);
    subplot(highStimRHighR_axes);
    caxis([min(iMin) max(iMax)]);
    
    was_a_key = waitforbuttonpress;
    if was_a_key && strcmp(get(fig, 'CurrentKey'), 'leftarrow')
      k = max(1, k - 1);
    elseif was_a_key && strcmp(get(fig, 'CurrentKey'), 'rightarrow')
      k = min(max_k, k + 1);
    end
end