function [expInfo, neuralData, behavioralData] = processExperiment(expInfo)

fprintf('(1/6) loading session data...')

matchTag = expInfo(1).cellMatched;
expInfo = data.loadExpData(expInfo);

fprintf('done\n')
%% load cell data

fprintf('(2/6) loading cell data...')
if matchTag
    [allFcell, expInfo] = loadMatchedCellData(expInfo);
else
    try
        [allFcell, expInfo] = loadCellData(expInfo);
    catch
        [allFcell, expInfo] = loadCellData_multiExp(expInfo);
    end
end
fprintf('done\n')

%% get event times

fprintf('(3/6) extracting events/wheel...')
for iX = 1:length(expInfo)
    try %task
        try %old tl
            allEventTimes{iX} = getEventTimes(expInfo(iX), {'stimulusOnTimes' 'interactiveOnTimes' 'stimulusOffTimes'});
            allWheelMoves{iX} = getWheelMoves(expInfo(iX), allEventTimes{iX});

            if length(allEventTimes{iX}) == length(allWheelMoves{iX})
                behavioralData(iX) = struct('eventTimes',allEventTimes{iX},'wheelMoves', allWheelMoves{iX});
            else
                error('Event/wheel experiment lengths do not match')
            end
        catch %new tl
            allEventTimes{iX} = newtl.getEventTimes(expInfo(iX), {'stimulusOnTimes' 'interactiveOnTimes' 'stimulusOffTimes'});
            allWheelMoves{iX} = newtl.getWheelMoves(expInfo(iX), allEventTimes{iX});

            if length(allEventTimes{iX}) == length(allWheelMoves{iX})
                behavioralData(iX) = struct('eventTimes',allEventTimes{iX},'wheelMoves', allWheelMoves{iX});
            else
                error('Event/wheel experiment lengths do not match')
            end
        end
    catch %passive
        allEventTimes{iX} = passive.getEventTimes(expInfo(iX), {'stimulusOnTimes' 'stimulusOffTimes'});
        allWheelMoves{iX} = passive.getWheelMoves(expInfo(iX), allEventTimes{iX});
        behavioralData(iX) = struct('eventTimes',allEventTimes{iX},'wheelMoves', allWheelMoves{iX});
    end
end
fprintf('done\n')
    
%% collate cell responses across planes

fprintf('(4/6) extracting cell responses...')

[cellResps, respTimes] = getCellResps(expInfo, allFcell);
    

%% assemble into relevant structs, for tidiness

neuralData = struct('allFcell',allFcell,'cellResps',cellResps,'respTimes',respTimes);
fprintf('done\n')

%% align resps
fprintf('(5/6) computing ETAs...')
try
    [neuralData] = alignResps(expInfo, neuralData, behavioralData);
%     [neuralData] = alignNeuropil(expInfo, neuralData, behavioralData);
catch
    [neuralData] = alignResps(expInfo, neuralData, behavioralData,2,{'stimulusOnTimes' 'stimulusOffTimes'});
end
fprintf('done\n')

%% get significance
fprintf('(6/6) computing stats...')
try
    [neuralData] = getSignificantActivity(expInfo, behavioralData, neuralData,0);
catch
    [neuralData] = passive.getSignificantActivity(expInfo, behavioralData, neuralData,0);
end
fprintf('done!\n')

    %%