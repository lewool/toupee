function [expInfo, behavioralData] = processBehavior(expInfo)

expInfo = data.loadExpData(expInfo);
 
%% get event times

allEventTimes = getEventTimes(expInfo, {'stimulusOnTimes' 'interactiveOnTimes' 'stimulusOffTimes'});
allWheelMoves = getWheelMoves(expInfo, allEventTimes);

if length(allEventTimes) == length(allWheelMoves)
    behavioralData = struct('eventTimes',allEventTimes,'wheelMoves', allWheelMoves);
else
    error('Event/wheel experiment lengths do not match')
end
