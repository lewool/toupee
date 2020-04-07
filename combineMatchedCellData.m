
[allFcell, expInfo] = loadMatchedCellData(mouseName);
alignedRespsInd = cell(length(expInfo),length(events));
%% loop through experiments and align cells/events per session

for ex = 1:length(expInfo)

    % get event timings and wheel trajectories
    [eventTimes{ex}, wheelTrajectories{ex}] = getEventTimes(expInfo(ex), {'stimulusOnTimes' 'interactiveOnTimes' 'stimulusOffTimes'});

    % load traces
    [cellResps{ex}, respTimes{ex}] = getCellResps(expInfo(ex), allFcell(ex,:));
    cellResps{ex} = zscore(cellResps{ex});

    % align calcium traces to the event you want
    events = {'stimulusOnTimes' 'prestimulusQuiescenceEndTimes' 'rewardOnTimes'};
    
    for ev = 1:length(events)
        [alignedRespsInd{ex,ev}, eventWindow] = alignResps(expInfo(ex), cellResps{ex}, respTimes{ex}, eventTimes{ex}, events{ev});
    end
end

%%
alignedResps = cell(1,length(events));
for ev = 1:length(events)
    alignedResps{ev} = cat(1,alignedRespsInd{:,ev});
end

