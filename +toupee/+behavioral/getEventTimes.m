function [expInfo, eventTimes] =...
    getEventTimes(expInfo, eventNames, varargin)
% Gets event times according to multiple clocks
%
% Event times for experiments are recorded in signals times (the times that
% code updated the value of an event on the stimulus computer) and in
% daq times (the times of values the ni-daq i/o device read from a sensor
% (e.g. the wheel or a photodiode) or triggered an actuator (e.g. the
% reward valve). The daq times can be thought of as rig times (i.e. the 
% time an event was experienced by a subject on the rig). For many signals
% events, there is no direct corresponding rig time (e.g. when the 
% interactive period of a trial begins), so this code estimates what the 
% rig time would be for those events.
%
%
% Inputs:
% -------
% expInfo : table
%   A table containing relevant information and data variables (columns) 
%   for particular experiment sessions (rows).
%
% eventNames : char array OR cell array
%   Names of events whose times to return. Each cell in a cell array
%   contains the name of one event. Must match an event name in a block
%   file's 'events', 'inputs', or 'outputs' struct, and must be an event
%   that occurs exactly once per trial.
%
% 'sessions' : int scalar OR int array OR char array OR cell array 
% (optional name-value pair)
%   Specific sessions for which to return the events' times, instead of
%   from all sessions. (default = all)
%
% 'phdFlipThresh' : double array (optional name-value pair)
%   The percentage change of the photodiode value (wrt to the max 
%   photodiode value) over a specified amount time to qualify as a visual
%   stimulus flip event. (default = 5% change over 10 ms)
%
%
% Outputs:
% --------
% expInfo : table
%   The updated `expInfo`.
%
% eventTimes : table
%   A 1-row table containing as columns (one for each session in `expInfo`)
%   tables with three columns: event names, 'signalsTimes', and 'rigTimes'.
%
%
% Examples:
% ---------
% 1) For a single session, get the stimulusOn, interactiveOn, stimulusOff,
% and reward event rig times.
%   deats = {'LEW031', '2020-02-03', 1};
%   files = {'block', 'timeline'};
%   expInfo = toupee.meta.processExperiment(deats, files);
%   eventNames = {'stimulusOn', 'interactiveOn', 'stimulusOff', 'reward'};
%   [expInfo, eventTimes] =...
%       toupee.behavioral.getEventTimes(expInfo, eventNames);
%
% 2) For multiple sessions, get the stimulusOn, interactiveOn, stimulusOff,
% and reward event rig times, and specify estimating a phd flip when the
% signal changes at least 7.5% of its max value of 10 ms.
%   deats = {{'LEW031', '2020-02-03', 1},...
%            {'LEW032', '2020-02-28', 1, [1, 2]}};
%   files = {'block', 'timeline'};
%   expInfo = toupee.meta.processExperiment(deats, files);
%   eventNames = {'stimulusOn', 'interactiveOn', 'stimulusOff', 'reward'};
%   [expInfo, eventTimes] =...
%       toupee.behavioral.getEventTimes(expInfo, eventNames,...
%                                       'phdFlipThresh', [0.075, 0.01]);
%
%
% See Also:
% ---------
% toupee.behavioral.getWheelMoves
% toupee.behavioral.getTrials
%
% @todo getting some event times requires special cases (e.g.
% 'interactiveOn')

% Check if eventTimes is already a column/table within behavioralData.
% If it is, merge tables.

%% Prerun checks.
% Import `iif`.
import toupee.misc.iif
% Ensure input args are of proper type.
p = inputParser;
% must be char or cell array of chars
isValidEvents = @(x) ischar(x) || (iscell(x) && ischar(x{1}));
% must be char or cell array of chars or scalar int or array of ints
isValidSessions = @(x) isValidEvents(x) || all(mod(x, 1) == 0);
% must be double 2 element array
isValidPhdFlipThresh = @(x) isnumeric(x) && numel(x) == 2;
addRequired(p, 'expInfo');
addRequired(p, 'events', isValidEvents);
% default value is all sessions
addParameter(p, 'sessions', expInfo.('Row'), isValidSessions);
% default value is 0.05
addParameter(p, 'phdFlipThresh', [0.05, 0.01], isValidPhdFlipThresh);
parse(p, expInfo, eventNames, varargin{:});
expInfo = p.Results.expInfo;
eventNames = p.Results.events;
if ischar(eventNames), eventNames = {eventNames}; end  % cellify
sessions = p.Results.sessions;
phdFlipThresh = p.Results.phdFlipThresh;
% Initialize table for all sessions, and table for individual sesisons:
eventTimes = table();
% If 'reward' was specified as an event, add an additional 'estReward' row
% to contain "woulda been" estimated reward times for incorrect trials.
if any(contains(eventNames, 'reward'))
    rowNames = [eventNames, {'estReward'}];
else
    rowNames = eventNames;
end
colNames = {'signalsTimes', 'rigTimes'};
colTypes = {'cell', 'cell'};
nR = numel(rowNames);  % number of rows
nC = numel(colNames);  % number of columns
%% Get event times.
nE = numel(sessions);  % number of experiment sessions
% Get event times for each session
for iE = 1:nE
    % Initialize table for current session.
    curEventTimes = table('Size', [nR, nC], 'VariableNames', colNames,...
                          'VariableTypes', colTypes, 'RowNames', rowNames);
    % Extract relevant data from this session.
    expRef = sessions{iE};  % session expRef
    block = expInfo.BlockFile{expRef};  % block data
    evts = block.events;  % signals events data
    ins = block.inputs;  % signals inputs data
    outs = block.outputs;  % signals outputs data
    allEvts = [evts, ins, outs];
    timeline = expInfo.TimelineFile{expRef};  % timeline data
    timelineTimes = timeline.rawDAQTimestamps{1}';
    nS = timeline.rawDAQSampleCount{1};  % total number of daq samples
    daqFs = diff(timelineTimes(1:2));  % daq sampling rate
    % Use photodiode flip times (i.e. a change in photodiode current value
    % that can be inferred as due to a visual stimulus flip-to-screen 
    % event) to estimate daq times of events from those events' signals 
    % times.
    % raw photodiode data from daq
    phdRaw = timeline.rawDAQData{1}...
             (:, strcmp('photoDiode', timeline.hw.inputs{1}.name));
    phdRawNorm = phdRaw / max(phdRaw);  % unit normalized photodiode vals
    sThresh = phdFlipThresh(2) / daqFs;  % n samples threshold for flip
    vThresh = phdFlipThresh(1);  % V % change threshold for flip
    flipMask = false(nS, 1);  % mask for flip at sample
    % Get photodiode flip times: changes in `phdRawNorm` of more than 
    % `vThresh` in less than `sThresh` indicate a flip
    % @todo gotta be a faster way to compute running peak2peak
    for iS = 1:(nS - sThresh)
        if peak2peak(phdRawNorm(iS:(iS + sThresh))) > vThresh
           flipMask(iS) = true;
        end
    end
    flipIdxsLast = find(flipMask);
    % get just a single time for each flip event (at the end of the flip)
    flipIdxsLast = flipIdxsLast(diff(flipIdxsLast) > 1);
    flipTimes = timelineTimes(flipIdxsLast);
    % For each time of an event, take the first photodiode flip time after
    % the signals time as an estimate of the rig time.
    nT = numel(evts.endTrialValues{1});  % number of completed trials
    for iN = 1:numel(eventNames)
        switch eventNames
            % Special case for reward: no need to use photodiode, instead
            % use raw reward trace from daq and estimate what the reward
            % time would have been for incorrect trials.
            case contains('reward')
                % Get all reward times.
                % reward trace from daq
                rewardTrace = timeline.rawDAQData{1}...
                    (:, strcmp('rewardCommand',...
                               timeline.hw.inputs{1}.name));
                % threshold for reward delivery
                rewardThresh = max(rewardTrace) / 2;
                % indices of reward trace where reward valve was open
                % (includes hotkey delivered rewards)
                allRewardIdxs = find(rewardTrace > rewardThresh);
                % get just a single time for each reward valve opening
                allRewardIdxsFirst =...
                    [allRewardIdxs(1);...
                     allRewardIdxs(find(diff(allRewardIdxs) > 1) + 1)];
                allRewardTimes = timelineTimes(allRewardIdxsFirst);
                % @todo:
                % Remove hotkey rewards to filter down to task rewards.
                hotkeyRewardMask = strcmp('r', ins.keyboardValues{1});
                signalsHotkeyRewardTimes =...
                    ins.keyboardValues{1}(hotkeyRewardMask);
                % interpolate signals hotkey reward times to daq reward
                % times
                rigHotkeyRewardTimes =...
                    interp1(allRewardTimes, allRewardTimes,...
                            signalsHotkeyRewardTimes, 'nearest');
                [~, hkIdx] = intersect(allRewardTimes,...
                                       rigHotkeyRewardTimes);
                taskRewardTimes = allRewardTimes;
                taskRewardTimes(hkIdx) = [];
                % Estimate "woulda been" reward times for incorrect trials:
                % this is the signals incorrect feedback times + (the mean
                % diff between signals output reward times and signals
                % correct feedback times) + (the mean diff between rig
                % output reward times and signals output reward times).
                signalsIncorrectFeedbackTimes =...
                    evts.feedbackTimes{1}(~evnts.feedbackValues{1});
                meanSignalsOutputFeedbackDiff =...
                    mean(outs.rewardTimes{1}(:)...
                         - signalsCorrectFeedbackTimes(:));
                signalsCorrectFeedbackTimes =...
                    evts.feedbackTimes{1}(evnts.feedbackValues{1});
                meanSignalsRigRewDiff =...
                    mean(taskRewardTimes(:)...
                         - signalsCorrectFeedbackTimes(:));
                estIncorrectRewardTimes = signalsIncorrectFeedbackTimes(:)...
                                          + meanSignalsOutputFeedbackDiff...
                                          + meanSignalsRigRewDiff;
                estTaskRewardTimes = sort([estIncorrectRewardTimes;...
                                           taskRewardTimes]);
                % Assign times to table.
                curEventTimes.rigTimes{'reward'} = taskRewardTimes;
                curEventTimes.rigTimes{'estReward'} = estTaskRewardTimes;
                curEventTimes.signalsTimes{'reward'} =...
                    outs.rewardTimes{1}(:);
                curEventTimes.signalsTimes{'estReward'} =...
                    evts.feedbackTimes{1}(:) + meanSignalsOutputFeedbackDiff;
            % For all other events, use photodiode flip event to estimate
            % rig time from signals time. @todo this may not be best prac..
            otherwise
                name = iif(contains(eventNames{iN}, 'Times'),...
                           eventNames{iN}, strcat(eventNames{iN}));
                signalsTimes = allEvts.(name){1}(1:nT);
                % First find closest 'stimWindowUpdate' time to event time.
                swut = block.stimWindowUpdateTimes{1};
                updateTimes = interp1(swut, swut, signalsTimes(:),...
                                      'nearest', 'extrap');
                % Then find closest phd flip to 'stimWindowUpdate' time.
                rigTimes = zeros(nT, 1);
                for iT = 1:nT
                    % must be at least a 10 ms delay between 'stimWindowUpdate'
                    % and nearest phd flip
                    rigTimes(iT) = find((updateTimes(iT) - flipTimes)...
                                        < -0.01, 1, 'first');
                end
                % Assign times to table.
                curEventTimes.signalsTimes{(name)} = signalsTimes;
                curEventTimes.rigTimes{(name)} = rigTimes;  
        end
    end

    % Assign table for current session to `eventTimes` & `expInfo` tables.
    eventTimes.(expRef) = curEventTimes;
    expInfo.behavioralData{expRef, 'eventTimes'} = {curEventTimes};
    
end

end



