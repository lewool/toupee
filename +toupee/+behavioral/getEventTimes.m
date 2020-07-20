function [expInfo, eventTimes] =...
    getEventTimes(expInfo, eventNames, varargin)
% Gets event times for multiple clocks
%
% Event times for experiments are recorded in signals times (the times that
% code updated the value of an event on the stimulus computer) and in
% daq times (the times of values the ni-daq i/o device read from a sensor
% (e.g. the wheel or a photodiode) or triggered an actuator (e.g. the
% reward valve). The daq times can be thought of as rig times (i.e. the 
% time an event was experienced by a subject on the rig). For many signals
% events, there is no direct corresponding rig time (e.g. when the 
% interactive period of a trial begins), so this code estimates what the 
% rig time would be for those events, assuming those events coincide with
% a flip to screen of the visual stimulus.
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
%   from all sessions. (Default: all)
%
% 'phdFlipThresh' : double array (optional name-value pair)
%   The percentage change of the photodiode value (wrt to the max 
%   photodiode value) over a specified amount time to qualify as a visual
%   stimulus flip event. (Default: 5% change over 10 ms)
% 
% 'scrnRefreshRate' : int scalar (optional name-value pair)
%   The refresh rate (in hz) of the physical screen that displays visual
%   stimuli.
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
%   eventNames = {'stimulusOn', 'interactiveOn', 'stimulusOff', 'response',...
%                 'reward'};
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
% @todo add support for event times not associated with stim flip
% @todo check to merge final table assignments

%% Prerun checks.
% Import `iif`.
import toupee.misc.iif
% Ensure input args are of proper type.
p = inputParser;
% must be char or cell array of chars
isValidEventNames = @(x) ischar(x) || (iscell(x) && ischar(x{1}));
% must be char or cell array of chars or scalar int or array of ints
isValidSessions = @(x) isValidEventNames(x) ...
                       || (isnumeric(x) && all(mod(x, 1) == 0));
% must be double 2 element array
isValidPhdFlipThresh = @(x) isnumeric(x) && numel(x) == 2;
isValidScrnRefreshRate = @(x) isscalar(x) && isnumeric(x) ...
                              && mod(x, 1) == 0;

addRequired(p, 'expInfo');
addRequired(p, 'eventNames', isValidEventNames);
% default value is all sessions
addParameter(p, 'sessions', expInfo.('Row'), isValidSessions);
% default value is 5% change over 10 ms
addParameter(p, 'phdFlipThresh', [0.05, 0.01], isValidPhdFlipThresh);
% default value is 60 hz
addParameter(p, 'scrnRefreshRate', 60, isValidScrnRefreshRate);

parse(p, expInfo, eventNames, varargin{:});
expInfo = p.Results.expInfo;
eventNames = p.Results.eventNames;
if ischar(eventNames), eventNames = {eventNames}; end  % cellify
sessions = p.Results.sessions;
phdFlipThresh = p.Results.phdFlipThresh;
scrnRefreshRate = p.Results.scrnRefreshRate;
% Initialize table for all sessions, and table for individual sesisons:
eventTimes = table();
% If 'reward' was specified as an event, add an additional 'estReward' row
% to contain "woulda been" estimated reward times for incorrect trials.
if any(contains(eventNames, 'reward'))
    rowNames = [eventNames, {'estReward'}];
else
    rowNames = eventNames;
end
colNames = {'signalsTimes', 'rigTimes', 'swu-FlipTimes'};
colTypes = {'cell', 'cell', 'cell'};
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
    timelineTimes = timeline.rawDAQTimestamps{1}(:);
    nS = timeline.rawDAQSampleCount{1};  % total number of daq samples
    % Ensure daq timestamps make sense.
    daqFsErrId = 'toupee:meta:processExperiment:badInput';
    daqFsErrMsg = 'Daq times were found to be non-linearly sampled.';
    %#ok<*COMPNOT>
    assert(all(diff(timelineTimes, 2)) == 0, daqFsErrId, daqFsErrMsg);
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
    vThresh = phdFlipThresh(1);  % V change threshold for flip
    flipMask = false(nS, 1);  % mask for flip at sample
    % Get photodiode flip times: change in `phdRawNorm` of more than 
    % `vThresh` in less than `sThresh` indicates a flip.
    % For each sample, see if it belongs to a flip by looking forward in
    % time within `sThresh` to see if there is a significant change in V.
    for iS = 1:(nS - sThresh)
        deltaV = max(abs(phdRawNorm(iS) - phdRawNorm(iS:(iS + sThresh))));
        if deltaV > vThresh
            flipMask(iS) = true;
        end
    end
    flipIdxsLast = find(flipMask);
    % Get just a single time for each flip event (at the end of the flip).
    % get minimum number of samples required between two flips, as a
    % function of physical screen refresh rate
    minFlipGap = ceil(1 / scrnRefreshRate * daqFs);
    flipIdxsLast = flipIdxsLast(diff(flipIdxsLast) > minFlipGap);
    flipTimes = timelineTimes(flipIdxsLast);
    % For each time of an event, take the first photodiode flip time after
    % the signals time as an estimate of the rig time.
    nT = numel(evts.endTrialValues{1});  % number of completed trials
    for iN = 1:numel(eventNames)
        if contains(eventNames{iN}, 'reward')
            % Special case for reward: no need to use photodiode, instead
            % use raw reward trace from daq and estimate what the reward
            % time would have been for incorrect trials.
            % Get all reward times.
            % reward trace from daq
            rewardTrace =...
                timeline.rawDAQData{1}...
                (:, strcmp('rewardCommand', timeline.hw.inputs{1}.name));
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
            % Remove hotkey rewards to filter down to task rewards.
            hotkeyRewardMask = strcmp('r', ins.keyboardValues{1});
            signalsHotkeyRewardTimes =...
                ins.keyboardValues{1}(hotkeyRewardMask);
            % interpolate signals hotkey reward times to daq reward
            % times
            rigHotkeyRewardTimes =...
                interp1(allRewardTimes, allRewardTimes,...
                        signalsHotkeyRewardTimes, 'nearest');
            [~, hkIdx] = intersect(allRewardTimes, rigHotkeyRewardTimes);
            taskRewardTimes = allRewardTimes;
            taskRewardTimes(hkIdx) = [];
            % Estimate "woulda been" reward times for incorrect trials:
            % this is the signals incorrect feedback times + (the mean
            % diff between signals output reward times and signals
            % correct feedback times) + (the mean diff between rig
            % output reward times and signals output reward times).
            signalsIncorrectFeedbackTimes =...
                evts.feedbackTimes{1}(~evts.feedbackValues{1});
            signalsCorrectFeedbackTimes =...
                evts.feedbackTimes{1}(evts.feedbackValues{1});
            meanSignalsOutputFeedbackDiff =...
                mean(outs.rewardTimes{1}(:)...
                - signalsCorrectFeedbackTimes(:));
            meanSignalsRigRewDiff =...
                mean(taskRewardTimes(:)...
                - outs.rewardTimes{1}(:));
            estIncorrectRewardTimes =...
                signalsIncorrectFeedbackTimes(:)...
                + meanSignalsOutputFeedbackDiff...
                + meanSignalsRigRewDiff;
            % includes reward times for correct trials, and estimated
            % reward times for incorrect trials
            estTaskRewardTimes = sort([taskRewardTimes;...
                                       estIncorrectRewardTimes]);
            % Assign times to table.
            curEventTimes.rigTimes{'reward'} = taskRewardTimes;
            curEventTimes.rigTimes{'estReward'} = estTaskRewardTimes;
            curEventTimes.signalsTimes{'reward'} =...
                outs.rewardTimes{1}(:);
            % includes signals reward times for correct trials, and
            % estimated signals reward times for incorrect trials
            curEventTimes.signalsTimes{'estReward'} = sort(...
                [outs.rewardTimes{1},...
                 (evts.feedbackTimes{1}((~evts.feedbackValues{1}))...
                  + meanSignalsOutputFeedbackDiff)])';
            % Compute and assign mean signals -> daq out delay to table.
            expInfo.behavioralData.daqDelay{iE} = meanSignalsRigRewDiff;
        else
            % For all other events, use photodiode flip event to estimate
            % rig time from signals time. @todo this may not be best prac..
            name = iif(contains(eventNames{iN}, 'Times'),...
                       eventNames{iN}, strcat(eventNames{iN}, 'Times'));
            signalsTimes = allEvts.(name){1}(1:nT);
            % First find closest 'stimWindowUpdate' time to event time.
            swut = block.stimWindowUpdateTimes{1};
            updateTimes = interp1(swut, swut, signalsTimes(:),...
                                  'nearest', 'extrap');
            % Then find closest phd flip to 'stimWindowUpdate' time.
            rigTimes = zeros(nT, 1);
            for iT = 1:nT
                rigTimes(iT) =...
                    flipTimes(find((updateTimes(iT) - flipTimes) < 0,...
                                   1, 'first'));
            end
            % estimated time diff between swu and flip
            swuFlipTimes = rigTimes - updateTimes;
            % Assign times to table.
            curEventTimes.('signalsTimes'){(name(1:(end - 5)))} =...
                signalsTimes(:);
            curEventTimes.('rigTimes'){(name(1:(end - 5)))} =...
                rigTimes(:);
            curEventTimes.('swu-FlipTimes'){(name(1:(end - 5)))} =...
                swuFlipTimes(:);
        end
    end

    % Assign table for current session to `eventTimes` & `expInfo` tables.
    eventTimes.(expRef) = {curEventTimes};
    expInfo.behavioralData{expRef, 'eventTimes'} = {curEventTimes};
    
    % Add trial durations.
    trialDurs = evts.endTrialTimes{1}(1:nT) - evts.newTrialTimes{1}(1:nT);
    expInfo.behavioralData{expRef, 'eventTimes'}{1}...
        .('signalsTimes'){'trialDur'} = trialDurs;
    
    % Get durations of 'interactiveDelay' periods if got 'stimulusOn' and
    % 'interactiveOn' times: 'interactiveOn' - 'stimulusOn' for each trial
    if sum(contains(eventTimes{1, expRef}{1}.Properties.RowNames,...
                    ["stimulusOn", "interactiveOn"])) == 2
        eventTimes{1, expRef}{1}.('signalsTimes'){'interactiveDelayDur'} =...
            eventTimes{1, expRef}{1}.('signalsTimes'){'interactiveOn'}...
            -  eventTimes{1, expRef}{1}.('signalsTimes'){'stimulusOn'};
        eventTimes{1, expRef}{1}.('rigTimes'){'interactiveDelayDur'} =...
            eventTimes{1, expRef}{1}.('rigTimes'){'interactiveOn'}...
            -  eventTimes{1, expRef}{1}.('rigTimes'){'stimulusOn'};
    end
    if sum(contains(expInfo.behavioralData{expRef, 'eventTimes'}{1}...
                    .Properties.RowNames, ["stimulusOn", "interactiveOn"]))...
                    == 2
        expInfo.behavioralData{expRef, 'eventTimes'}{1}...
            .('signalsTimes'){'interactiveDelayDur'} =...
            expInfo.behavioralData{expRef, 'eventTimes'}{1}...
                .('signalsTimes'){'interactiveOn'}...
            -  expInfo.behavioralData{expRef, 'eventTimes'}{1}...
                .('signalsTimes'){'stimulusOn'};
        expInfo.behavioralData{expRef, 'eventTimes'}{1}...
            .('rigTimes'){'interactiveDelayDur'} =...
            expInfo.behavioralData{expRef, 'eventTimes'}{1}...
                .('rigTimes'){'interactiveOn'}...
            -  expInfo.behavioralData{expRef, 'eventTimes'}{1}...
                .('rigTimes'){'stimulusOn'};
    end
    
    % Get durations of 'interactiveOn' periods if got 'interactiveOn' and
    % 'response' times: 'response' - 'interactiveOn' for each trial
    if sum(contains(eventTimes{1, expRef}{1}.Properties.RowNames,...
                    ["interactiveOn", "response"])) == 2
        eventTimes{1, expRef}{1}.('signalsTimes'){'interactiveDur'} =...
            eventTimes{1, expRef}{1}.('signalsTimes'){'response'}...
            -  eventTimes{1, expRef}{1}.('signalsTimes'){'interactiveOn'};
        eventTimes{1, expRef}{1}.('rigTimes'){'interactiveDur'} =...
            eventTimes{1, expRef}{1}.('rigTimes'){'response'}...
            -  eventTimes{1, expRef}{1}.('rigTimes'){'interactiveOn'};
    end
    if sum(contains(expInfo.behavioralData{expRef, 'eventTimes'}{1}...
                    .Properties.RowNames, ["interactiveOn", "response"]))...
                    == 2
        expInfo.behavioralData{expRef, 'eventTimes'}{1}...
            .('signalsTimes'){'interactiveDur'} =...
            expInfo.behavioralData{expRef, 'eventTimes'}{1}...
                .('signalsTimes'){'response'}...
            -  expInfo.behavioralData{expRef, 'eventTimes'}{1}...
                .('signalsTimes'){'interactiveOn'};
        expInfo.behavioralData{expRef, 'eventTimes'}{1}...
            .('rigTimes'){'interactiveDur'} =...
            expInfo.behavioralData{expRef, 'eventTimes'}{1}...
                .('rigTimes'){'response'}...
            -  expInfo.behavioralData{expRef, 'eventTimes'}{1}...
                .('rigTimes'){'interactiveOn'};
    end

end

end
