function [expInfo, behavioralData] = processBehavior(expInfo)

expInfo = data.loadExpData(expInfo);
 
%% get event times

% allEventTimes = getEventTimes(expInfo, {'stimulusOnTimes' 'interactiveOnTimes' 'stimulusOffTimes'});
% allWheelMoves = getWheelMoves(expInfo, allEventTimes);
% 
% if length(allEventTimes) == length(allWheelMoves)
%     behavioralData = struct('eventTimes',allEventTimes,'wheelMoves', allWheelMoves);
% else
%     error('Event/wheel experiment lengths do not match')
% end

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