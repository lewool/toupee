function [expInfo, neuralData, behavioralData] = processExperiment(expInfo)


%% load cell data

matchTag = expInfo(1).cellMatched;
expInfo = data.loadExpData(expInfo);

if matchTag
    [allFcell, expInfo] = loadMatchedCellData(expInfo);
else
    try
        [allFcell, expInfo] = loadCellData(expInfo);
    catch
        [allFcell, expInfo] = loadCellData_multiExp(expInfo);
    end
end
 
%% get event times

try %task
    try %old tl
        allEventTimes = getEventTimes(expInfo, {'stimulusOnTimes' 'interactiveOnTimes' 'stimulusOffTimes'});
        allWheelMoves = getWheelMoves(expInfo, allEventTimes);

        if length(allEventTimes) == length(allWheelMoves)
            behavioralData = struct('eventTimes',allEventTimes,'wheelMoves', allWheelMoves);
        else
            error('Event/wheel experiment lengths do not match')
        end
    catch %new tl
        allEventTimes = newtl.getEventTimes(expInfo, {'stimulusOnTimes' 'interactiveOnTimes' 'stimulusOffTimes'});
        allWheelMoves = newtl.getWheelMoves(expInfo, allEventTimes);

        if length(allEventTimes) == length(allWheelMoves)
            behavioralData = struct('eventTimes',allEventTimes,'wheelMoves', allWheelMoves);
        else
            error('Event/wheel experiment lengths do not match')
        end
    end
catch %passive
    allEventTimes = passive.getEventTimes(expInfo, {'stimulusOnTimes' 'stimulusOffTimes'});
    allWheelMoves = passive.getWheelMoves(expInfo, allEventTimes);
    behavioralData = struct('eventTimes',allEventTimes,'wheelMoves', allWheelMoves);
end
%%


%% collate cell responses across planes

[cellResps, respTimes] = getCellResps(expInfo, allFcell);

%% assemble into relevant structs, for tidiness

neuralData = struct('allFcell',allFcell,'cellResps',cellResps,'respTimes',respTimes);

%% align resps

try
    [neuralData] = alignResps(expInfo, neuralData, behavioralData);
%     [neuralData] = alignNeuropil(expInfo, neuralData, behavioralData);
catch
    [neuralData] = alignResps(expInfo, neuralData, behavioralData,2,{'stimulusOnTimes'});
end

try
    [neuralData] = getSignificantActivity(expInfo, behavioralData, neuralData,0);
catch
    [neuralData] = passive.getSignificantActivity(expInfo, behavioralData, neuralData,0);
end
%%