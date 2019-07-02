%% load experiment details

% close all;
clear all;

mouseName = 'LEW005';
expDate = '2018-06-10';
expNum = 2;
expSeries = [2 3];

%% load data

[block, Timeline] = loadData(mouseName, expDate, expNum);

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
plotAll = chooseCellType('mov', mouseName, expDate, expNum, expSeries, block, allFcell, eventTimes, ops);

%% plot raster browser
fig = figure(203);
hold on
set(fig,'Position',[290   100   380   1220])
startIdx = find(eventWindow == -1);
endIdx = find(eventWindow == 2);

if ~exist('k') == 1
    k = 1;
end
max_k = length(plotAll);

% select the conditions
contrasts = unique(block.events.contrastValues);

[~, condIdx_hitL_highL] = selectCondition(block, contrasts(contrasts<0), eventTimes, 'all', 'all', 'all', 'left', 'correct', 'all', 'all', 'all', 'all');
[~, condIdx_hitL_highR] = selectCondition(block, contrasts(contrasts<0), eventTimes, 'all', 'all', 'all', 'right', 'correct', 'all', 'all', 'all', 'all');
[~, condIdx_errorL_highL] = selectCondition(block, contrasts(contrasts<0), eventTimes, 'all', 'all', 'all', 'left', 'incorrect', 'all', 'all', 'all', 'all');
[~, condIdx_errorL_highR] = selectCondition(block, contrasts(contrasts<0), eventTimes, 'all', 'all', 'all', 'right', 'incorrect', 'all', 'all', 'all', 'all');

[~, condIdx_hitR_highL] = selectCondition(block, contrasts(contrasts>0), eventTimes, 'all', 'all', 'all', 'left', 'correct', 'all', 'all', 'all', 'all');
[~, condIdx_hitR_highR] = selectCondition(block, contrasts(contrasts>0), eventTimes, 'all', 'all', 'all', 'right', 'correct', 'all', 'all', 'all', 'all');
[~, condIdx_errorR_highL] = selectCondition(block, contrasts(contrasts>0), eventTimes, 'all', 'all', 'all', 'left', 'incorrect', 'all', 'all', 'all', 'all');
[~, condIdx_errorR_highR] = selectCondition(block, contrasts(contrasts>0), eventTimes, 'all', 'all', 'all', 'right', 'incorrect', 'all', 'all', 'all', 'all');

% set up the subplot sizes based on trial numbers
maxLength_hitL = max([length(condIdx_hitL_highL) length(condIdx_hitL_highR)]);
maxLength_errorL = max([length(condIdx_errorL_highL) length(condIdx_errorL_highR)]);
maxLength_hitR = max([length(condIdx_hitR_highL) length(condIdx_hitR_highR)]);
maxLength_errorR = max([length(condIdx_errorR_highL) length(condIdx_errorR_highR)]);
subplotLength = maxLength_hitL + maxLength_errorL + maxLength_hitR + maxLength_errorR;
marginSize = round(subplotLength*.05);

figureLength = subplotLength + 3*marginSize;
figureSize = 2*figureLength;
startsubplot1 = 1;
startsubplot2 = (maxLength_hitL + marginSize)*2-1;
startsubplot3 = (maxLength_hitL + maxLength_errorL + 2*marginSize)*2 - 1;
startsubplot4 = (maxLength_hitL + maxLength_errorL + maxLength_hitR + 3*marginSize)*2 -1;

hitLHighL_axes = subplot(figureLength,2,startsubplot1:2:2*length(condIdx_hitL_highL));
hitLHighR_axes = subplot(figureLength,2,startsubplot1+1:2:2*length(condIdx_hitL_highR));

errorLHighL_axes = subplot(figureLength,2,startsubplot2:2:startsubplot2+2*length(condIdx_errorL_highL));
errorLHighR_axes = subplot(figureLength,2,startsubplot2+1:2:startsubplot2+2*length(condIdx_errorL_highR));

hitRHighL_axes = subplot(figureLength,2,startsubplot3:2:startsubplot3+2*length(condIdx_hitR_highL));
hitRHighR_axes = subplot(figureLength,2,startsubplot3+1:2:startsubplot3+2*length(condIdx_hitR_highR));

errorRHighL_axes = subplot(figureLength,2,startsubplot4:2:startsubplot4+2*length(condIdx_errorR_highL));
errorRHighR_axes = subplot(figureLength,2,startsubplot4+1:2:startsubplot4+2*length(condIdx_errorR_highR));

while k <= max_k
    
    subplot(hitLHighL_axes);cla;
    subplot(errorLHighL_axes);cla;
    subplot(hitRHighL_axes);cla;
    subplot(errorRHighL_axes);cla;
    subplot(hitLHighR_axes);cla;
    subplot(errorLHighR_axes);cla;
    subplot(hitRHighR_axes);cla;
    subplot(errorRHighR_axes);cla

    subplot(hitLHighL_axes);
    relativeStimTimes = zeros(length(condIdx_hitL_highL),1)';
    relativeMovTimes = eventTimes(7).daqTime(condIdx_hitL_highL)-eventTimes(1).daqTime(condIdx_hitL_highL);
    [~,sortIdx] = sort(relativeMovTimes);
    imagesc(eventWindow(startIdx:endIdx), 1:length(condIdx_hitL_highL),stim_alignedTraces{plotAll(k,1)}.eventSpikes(condIdx_hitL_highL(sortIdx),startIdx:endIdx,plotAll(k,2)));
    hold on;
    plot(relativeMovTimes(sortIdx),1:length(condIdx_hitL_highL),'k.');
    stimLine = line([0 0],[1 length(condIdx_hitL_highL)]);
    set(stimLine,'Color','k');
    colormap(flipud(gray));
    xlim([-1 3]);
    box off
    set(gca,'ytick',[])
    set(gca,'xtick',[])
    title('highL')
    ylabel('hitL')

    subplot(hitLHighR_axes);
    relativeStimTimes = zeros(length(condIdx_hitL_highR),1)';
    relativeMovTimes = eventTimes(7).daqTime(condIdx_hitL_highR)-eventTimes(1).daqTime(condIdx_hitL_highR);
    [~,sortIdx] = sort(relativeMovTimes);
    imagesc(eventWindow(startIdx:endIdx), 1:length(condIdx_hitL_highR),stim_alignedTraces{plotAll(k,1)}.eventSpikes(condIdx_hitL_highR(sortIdx),startIdx:endIdx,plotAll(k,2)));
    hold on;
    plot(relativeMovTimes(sortIdx),1:length(condIdx_hitL_highR),'k.');
    stimLine = line([0 0],[1 length(condIdx_hitL_highR)]);
    set(stimLine,'Color','k');
    colormap(flipud(gray));
    xlim([-1 3]);
    box off
    set(gca,'ytick',[])
    set(gca,'xtick',[])
    title('highR')

    subplot(errorLHighL_axes);
    relativeStimTimes = zeros(length(condIdx_errorL_highL),1)';
    relativeMovTimes = eventTimes(7).daqTime(condIdx_errorL_highL)-eventTimes(1).daqTime(condIdx_errorL_highL);
    [~,sortIdx] = sort(relativeMovTimes);
    imagesc(eventWindow(startIdx:endIdx), 1:length(condIdx_errorL_highL),stim_alignedTraces{plotAll(k,1)}.eventSpikes(condIdx_errorL_highL(sortIdx),startIdx:endIdx,plotAll(k,2)));
    hold on;
    plot(relativeMovTimes(sortIdx),1:length(condIdx_errorL_highL),'k.');
    stimLine = line([0 0],[1 length(condIdx_errorL_highL)]);
    set(stimLine,'Color','k');
    colormap(flipud(gray));
    xlim([-1 3]);
    box off
    set(gca,'ytick',[])
    set(gca,'xtick',[])
    ylabel('errorL')

    subplot(errorLHighR_axes);
    relativeStimTimes = zeros(length(condIdx_errorL_highR),1)';
    relativeMovTimes = eventTimes(7).daqTime(condIdx_errorL_highR)-eventTimes(1).daqTime(condIdx_errorL_highR);
    [~,sortIdx] = sort(relativeMovTimes);
    imagesc(eventWindow(startIdx:endIdx), 1:length(condIdx_errorL_highR),stim_alignedTraces{plotAll(k,1)}.eventSpikes(condIdx_errorL_highR(sortIdx),startIdx:endIdx,plotAll(k,2)));
    hold on;
    plot(relativeMovTimes(sortIdx),1:length(condIdx_errorL_highR),'k.');
    stimLine = line([0 0],[1 length(condIdx_errorL_highR)]);
    set(stimLine,'Color','k');
    colormap(flipud(gray));
    xlim([-1 3]);
    box off
    set(gca,'ytick',[])
    set(gca,'xtick',[])

    subplot(hitRHighL_axes);
    relativeStimTimes = zeros(length(condIdx_hitR_highL),1)';
    relativeMovTimes = eventTimes(7).daqTime(condIdx_hitR_highL)-eventTimes(1).daqTime(condIdx_hitR_highL);
    [~,sortIdx] = sort(relativeMovTimes);
    imagesc(eventWindow(startIdx:endIdx), 1:length(condIdx_hitR_highL),stim_alignedTraces{plotAll(k,1)}.eventSpikes(condIdx_hitR_highL(sortIdx),startIdx:endIdx,plotAll(k,2)));
    hold on;
    plot(relativeMovTimes(sortIdx),1:length(condIdx_hitR_highL),'k.');
    stimLine = line([0 0],[1 length(condIdx_hitR_highL)]);
    set(stimLine,'Color','k');
    colormap(flipud(gray));
    xlim([-1 3]);
    box off
    set(gca,'ytick',[])
    set(gca,'xtick',[])
    ylabel('hitR')

    subplot(hitRHighR_axes);
    relativeStimTimes = zeros(length(condIdx_hitR_highR),1)';
    relativeMovTimes = eventTimes(7).daqTime(condIdx_hitR_highR)-eventTimes(1).daqTime(condIdx_hitR_highR);
    [~,sortIdx] = sort(relativeMovTimes);
    imagesc(eventWindow(startIdx:endIdx), 1:length(condIdx_hitR_highR),stim_alignedTraces{plotAll(k,1)}.eventSpikes(condIdx_hitR_highR(sortIdx),startIdx:endIdx,plotAll(k,2)));
    hold on;
    plot(relativeMovTimes(sortIdx),1:length(condIdx_hitR_highR),'k.');
    stimLine = line([0 0],[1 length(condIdx_hitR_highR)]);
    set(stimLine,'Color','k');
    colormap(flipud(gray));
    xlim([-1 3]);
    box off
    set(gca,'ytick',[])
    set(gca,'xtick',[])

    subplot(errorRHighL_axes);
    relativeStimTimes = zeros(length(condIdx_errorR_highL),1)';
    relativeMovTimes = eventTimes(7).daqTime(condIdx_errorR_highL)-eventTimes(1).daqTime(condIdx_errorR_highL);
    [~,sortIdx] = sort(relativeMovTimes);
    imagesc(eventWindow(startIdx:endIdx), 1:length(condIdx_errorR_highL),stim_alignedTraces{plotAll(k,1)}.eventSpikes(condIdx_errorR_highL(sortIdx),startIdx:endIdx,plotAll(k,2)));
    hold on;
    plot(relativeMovTimes(sortIdx),1:length(condIdx_errorR_highL),'k.');
    stimLine = line([0 0],[1 length(condIdx_errorR_highL)]);
    set(stimLine,'Color','k');
    colormap(flipud(gray));
    xlim([-1 3]);
    set(gca,'ytick',[])
    box off
    ylabel('errorR')

    subplot(errorRHighR_axes);
    relativeStimTimes = zeros(length(condIdx_errorR_highR),1)';
    relativeMovTimes = eventTimes(7).daqTime(condIdx_errorR_highR)-eventTimes(1).daqTime(condIdx_errorR_highR);
    [~,sortIdx] = sort(relativeMovTimes);
    imagesc(eventWindow(startIdx:endIdx), 1:length(condIdx_errorR_highR),stim_alignedTraces{plotAll(k,1)}.eventSpikes(condIdx_errorR_highR(sortIdx),startIdx:endIdx,plotAll(k,2)));
    hold on;
    plot(relativeMovTimes(sortIdx),1:length(condIdx_errorR_highR),'k.');
    stimLine = line([0 0],[1 length(condIdx_errorR_highR)]);
    set(stimLine,'Color','k');
    colormap(flipud(gray));
    xlim([-1 3]);
    set(gca,'ytick',[])
    box off

    was_a_key = waitforbuttonpress;
    if was_a_key && strcmp(get(fig, 'CurrentKey'), 'leftarrow')
      k = max(1, k - 1);
    elseif was_a_key && strcmp(get(fig, 'CurrentKey'), 'rightarrow')
      k = min(max_k, k + 1);
    end
end
