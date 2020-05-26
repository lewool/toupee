function [expInfo, neuralData, behavioralData] = processExperiment(expInfo)


%% load cell data

matchTag = expInfo(1).cellMatched;
expInfo = data.loadExpData(expInfo);

if matchTag
    [allFcell, expInfo] = loadMatchedCellData(expInfo);
else
    [allFcell, expInfo] = loadCellData(expInfo);
end
 
%% get event times

allEventTimes = getEventTimes(expInfo, {'stimulusOnTimes' 'interactiveOnTimes' 'stimulusOffTimes'});
allWheelMoves = getWheelMoves(expInfo, allEventTimes);

if length(allEventTimes) == length(allWheelMoves)
    behavioralData = struct('eventTimes',allEventTimes,'wheelMoves', allWheelMoves);
else
    error('Event/wheel experiment lengths do not match')
end

%% collate cell responses across planes

[cellResps, respTimes] = getCellResps(expInfo, allFcell);

%% assemble into relevant structs, for tidiness

neuralData = struct('allFcell',allFcell,'cellResps',cellResps,'respTimes',respTimes);

%% align resps

[neuralData] = alignResps(expInfo, neuralData, behavioralData);
[neuralData] = getSignificantActivity(expInfo, behavioralData, neuralData,0);

%%
