function [expInfo, wheelMoves] = getWheelMoves(expInfo, varargin)
% Gets per-trial wheel information
%
%
% Inputs:
% -------
% expInfo : table
%   A table containing relevant information and data variables (columns) 
%   for particular experiment sessions (rows).
%
% 'sessions' : int scalar OR int array OR char array OR cell array 
% (optional name-value pair)
%   Specific sessions for which to filter trials from, instead of from all
%   sessions. (Default: all)
%
% 'eventNames' : char array OR cell array (optional name-value pair)
%   Specified events (used as time-markers) around which wheel info is
%   returned. Each cell contains the name of an event. (Default: none)
%
% 'eventWindows' : numeric array OR cell array (optional name-value pair)
%   The time windows (in s) around each event specified in `eventNames`.
%   Each cell should contain a two-element numeric array containing the
%   time window to return: e.g. [-0.5 0.5] represents a half second before
%   and after the event. (Default: none)
%
% 'fs' : double scalar (optional name-value pair)
%   The sampling frequency (in hz) of the data acquisition device that
%   sampled the wheel. (Default: 1000 hz).
%
% 'gradFn' : function handle (optional name-value pair)
%   A function handle containing the function to use on the wheel position
%   data to compute the wheel velocity. (Default: The numerical gradient of
%   a single-pass, moving-average with a 10 ms window)
%
% 'wheelSpecs' : struct (optional name-value pair)
%   Contains specifications on wheel parameters. The possible fields and
%   their values are:
%       radius : double scalar
%           The radius of the wheel in m. (Default: 0.031 m)
%       res : int scalar
%           The total tick resolution of a full wheel revolution, including
%           the encoding. (Default: 400 (100 segments per revolution, X4
%           encoding for a kulber 100 ppr encoder))
%       gain : double scalar
%           The gain of the wheel displacement when used as a signal in
%           signals. (Default: 5)
%
%
% Outputs:
% --------
% expInfo : table
%   The updated `expInfo`.
%
% wheelMoves : table
%   A 1-row table of tables, where each column contains a wheel moves table
%   for a session. Each subtable contains rows for windows around specified
%   events, and the columns are continuous position, continuous velocity,
%   peak velocity, continuous acceleration, peak acceleration, continuous
%   direction, initial movement direction, final movement direction, and
%   continuous movement classification.
%
%
% Examples:
% ---------
% 1) For multiple sessions, for all trials, get all wheel moves before and
% after 0.5 seconds of the following events: 'stimulusOn', 'interactiveOn',
% 'stimulusOff', 'response', 'estReward'
%   deats = {{'LEW031', '2020-02-03', 1},...
%            {'LEW032', '2020-02-28', 1, [1, 2]}};
%   files = {'block', 'timeline'};
%   expInfo = toupee.meta.processExperiment(deats, files);
%   eventNames = {'stimulusOn', 'interactiveOn', 'stimulusOff', 'response',...
%                 'reward'};
%   expInfo =...
%       toupee.behavioral.getEventTimes(expInfo, eventNames,...
%                                       'phdFlipThresh', [0.075, 0.01]);
%   eventNames = {'stimulusOn', 'interactiveOn', 'stimulusOff', 'response',...
%                 'estReward'};
%   eventWindows = {[-0.5 0.5], [-0.5 0.5], [-0.5 0.5], [-0.5 0.5],...
%                   [-0.5 0.5]};
%   fs = 1000;
%   gradFn = @(x) gradient(movmean(x, 10));
%   wheelSpecs.radius = 0.031; wheelSpecs.res = 400; wheelSpecs.gain = 5;
%   expInfo = ...
%       toupee.behavioral.getWheelMoves(expInfo, 'eventNames', eventNames, ...
%                                       'eventWindows', eventWindows, ...
%                                       'fs', fs, 'gradFn', gradFn, ...
%                                       'wheelSpecs', wheelSpecs); 
%
% get event times, get wheel moves, get trials
% 
% 
% See Also:
% ---------
% toupee.behavioral.getTrials
% toupee.behavioral.getTrialBlocks
% toupee.behavioral.getEventTimes
% toupee.behavioral.wheel.getMoves
% toupee.behavioral.wheel.computeVelAcc

%% Prerun checks.
% Imports
import toupee.behavioral.*
% Turn off warning for assigning to a subset of rows of a table at a time.
warning('off', 'MATLAB:table:RowsAddedNewVars')
warning('off', 'MATLAB:table:RowsAddedExistingVars')
% Ensure input args are of proper type.
p = inputParser;
% must be char or cell array of chars or scalar int or array of ints
isValidSessions = @(x) isValidEvents(x) || all(mod(x, 1) == 0);
% must be char or cell array of chars
isValidEventNames = @(x) ischar(x) || (iscell(x) && ischar(x{1}));
% must be a two-elem vector with second elem greater than first, or a cell
% array of such vectors
isValidEventWindows = @(x) isnumeric(x) && numel(x) == 2 && x(2) > x(1) ...
                      || (iscell(x) && all(cellfun(@(y) ...
                                           isnumeric(y) ...
                                           && numel(y) == 2 ...
                                           && y(2) > y(1), x)));
% must be an int scalar
isValidFs = @(x) isnumeric(x) && isscalar(x) && mod(x, 1) == 0;
% must be a function handle
isValidGradFn = @(x) isa(x, 'function_handle');
% must be a struct with proper fields
isValidWheelSpecs = ...
    @(x) isstruct(x) ...
         && all(contains(fieldnames(x), ["radius", "res", "gain"])) ...
         && (isscalar(x.radius) && isnumeric(x.radius) && x.radius > 0) ...
         && (isscalar(x.res) && isnumeric(x.res) && x.res > 0) ...
         && (isscalar(x.gain) && isnumeric(x.gain) && x.gain > 0);

addRequired(p, 'expInfo');
addParameter(p, 'sessions', expInfo.('Row'), isValidSessions);
addParameter(p, 'eventNames', {}, isValidEventNames);
addParameter(p, 'eventWindows', {}, isValidEventWindows);
addParameter(p, 'fs', 1000, isValidFs);
addParameter(p, 'gradFn', @(x) gradient(movmean(x, 10)), isValidGradFn);
addParameter(p, 'wheelSpecs', ...
             struct('radius', 0.031, 'res', 400, 'gain', 5), ...
             isValidWheelSpecs);

parse(p, expInfo, varargin{:});
sessions = p.Results.sessions;
sessions = expInfo.Row(sessions);
eventNames = p.Results.eventNames;
eventWindows = p.Results.eventWindows;
fs = p.Results.fs;
gradFn = p.Results.gradFn;
wheelSpecs = p.Results.wheelSpecs;

%% Get wheel moves.
% Initialize `wheelMoves` table.
rowNames = ['fullTrials', ...
            cellfun(@(x, y) [x, ': [', num2str(y(1)), ', ', num2str(y(2)), ...
                             ']'], eventNames, eventWindows, 'uni', 0)]';
nR = numel(rowNames);  % number of rows
colNames = {'rigTimes', 'position', 'velocity', 'acceleration', 'nMoves', ...
            'moveOn', 'moveOff', 'moveDisplacement', 'moveDirection', ...
            'moveClass', 'movePeakVelocity', 'movePeakAcceleration'};
nC = numel(colNames);  % number of columns
colTypes = cellstr(repmat('cell', nC, 1));
wheelMoves = table('Size', [nR, nC], 'VariableNames', colNames, ...
                   'VariableTypes', colTypes, 'RowNames', rowNames);
% Create placeholder elements in `eventNames` and `eventWindows` for
% returning wheel moves for each full trial
eventNames = [{'fullTrials'}, eventNames];
eventWindows = [{[0, 0]}, eventWindows];
% Go experiment-by-experiment, event-by-event, trial-by-trial.
% the resolution of wheel movement in m
wheelTick2MeterX = (2 * wheelSpecs.radius * pi) / wheelSpecs.res;
nE = numel(sessions);  % number of experiment sessions 
for iE = 1:nE
    % Extract relevant data from this session.
    expRef = sessions{iE};  % session expRef
    block = expInfo.BlockFile{expRef};  % block data
    timeline = expInfo.TimelineFile{expRef};  % timeline data
    evts = block.events;  % events data
    nT = numel(evts.endTrialValues{1});  % number of completed trials
    eventTimes = expInfo.behavioralData.eventTimes;  % event times
    % block wheel position converted to m from counter values
    xRaw = block.inputs.wheelValues{1}(2:end)' * wheelTick2MeterX;
    % block time in s
    tRaw = block.inputs.wheelTimes{1}(2:end)';
    % Ensure timestamps of wheel data are monotonically increasing.
    tErrId = 'toupee:behavioral:getWheelMoves:nonmonotonicT';
    tErrMsg = ...
        ['The timestamps of the wheel data (block.inputs.wheelTimes) ' ...
         'are not monotonically increasing.'];
    assert(all(diff(tRaw) > 0), tErrId, tErrMsg);
    % Estimate wheel time in rig time from block time: estimate the delay 
    % time for timeline input -> block input (wheel data) as double the 
    % delay time as block output -> timeline output (reward data), and
    % subtract it.  @todo acknowledge `daqDelay` comes from `getEventTimes`
    daqDelay = expInfo.behavioralData.daqDelay{iE};
    tRaw = tRaw - (daqDelay * 2);
    t = timeline.rawDAQTimestamps{1}(:);  % timeline time
    % Interpolate the wheel trace to linearly sample in rig time.
    x = interp1(tRaw, xRaw, t, 'nearest', 'extrap');
    % Compute instantaneous velocity and acceleration for entire session.
    [v, a] = wheel.computeVelAcc(x, t, 'fs', fs, 'gradFn', gradFn);
    % Per specified event, per trial, get:
    % continuous position, continuous velocity, peak velocity, continuous
    % acceleration, peak acceleration, continuous direction, initial 
    % direction, final direction, and continuous movement classification.
    tic
    for iN = 1:numel(eventNames)
        rowName = rowNames{iN};
        evt = eventNames{iN};
        % Initialize cells to get all wheel move info per event.
        tEvt = cell(nT, 1);
        xEvt = cell(nT, 1);
        vEvt = cell(nT, 1);
        aEvt = cell(nT, 1);
        nMovesEvt = cell(nT, 1);
        moveOnEvt = cell(nT, 1);
        moveOffEvt = cell(nT, 1);
        moveDisplacementEvt = cell(nT, 1);
        moveDirectionEvt = cell(nT, 1);
        moveClassEvt = cell(nT, 1);
        movePeakVelocityEvt = cell(nT, 1);
        movePeakAccelerationEvt = cell(nT, 1);
        % For each trial, get wheel move info for specified time windows
        % around event.
        for iT = 1:nT
            switch evt
                % Special case to get indices for full trials: don't have
                % rigTimes for newTrial + endTrial, so estimate them from
                % block data.
                case 'fullTrials'  
                startTime = evts.newTrialTimes{1}(iT) + daqDelay;
                endTime = evts.endTrialTimes{1}(iT) + daqDelay;
                % For all other events, use `rigTimes` column from
                % `eventTimes` table
                otherwise
                startTime = eventTimes{iE}{evt, 'rigTimes'}{1}(iT) + eventWindows{iN}(1);
                endTime = eventTimes{iE}{evt, 'rigTimes'}{1}(iT) + eventWindows{iN}(2);
            end
            iStartEvt = round(startTime * fs);
            iEndEvt = round(endTime * fs);
            evtIdxs = iStartEvt:1:iEndEvt;
            tEvt{iT} = t(evtIdxs);
            xEvt{iT} = x(evtIdxs);
            vEvt{iT} = v(evtIdxs);
            aEvt{iT} = a(evtIdxs);
            [moveOn, moveOff, moveDisplacement, moveDirection, moveClass, ...
             movePeakVelocity, movePeakAcceleration] = ...
                wheel.getMoves(xEvt{iT}, tEvt{iT}, 'fs', fs, ...
                               'gradFn', gradFn);
            nMovesEvt{iT} = numel(moveOn);
            moveOnEvt{iT} = moveOn;
            moveOffEvt{iT} = moveOff;
            moveDisplacementEvt{iT} = moveDisplacement;
            moveDirectionEvt{iT} = moveDirection;
            moveClassEvt{iT} = moveClass;
            movePeakVelocityEvt{iT} = movePeakVelocity;
            movePeakAccelerationEvt{iT} = movePeakAcceleration;
        end
        % For each event, assign cells with wheel move info for all trials
        % to `wheelMoves` table
        wheelMoves{rowName, 'rigTimes'} = {tEvt};
        wheelMoves{rowName, 'position'} = {xEvt};
        wheelMoves{rowName, 'velocity'} = {vEvt};
        wheelMoves{rowName, 'acceleration'} = {aEvt};
        wheelMoves{rowName, 'nMoves'} = {nMovesEvt};
        wheelMoves{rowName, 'moveOn'} = {moveOnEvt};
        wheelMoves{rowName, 'moveOff'} = {moveOffEvt};
        wheelMoves{rowName, 'moveDisplacement'} = {moveDisplacementEvt};
        wheelMoves{rowName, 'moveDirection'} = {moveDirectionEvt};
        wheelMoves{rowName, 'moveClass'} = {moveClassEvt};
        wheelMoves{rowName, 'movePeakVelocity'} = {movePeakVelocityEvt};
        wheelMoves{rowName, 'movePeakAcceleration'} = ...
            {movePeakAccelerationEvt};
    end
    toc
    % Assign `wheelMoves` to proper expRef in `expInfo.behavioralData`
    expInfo.behavioralData{expRef, 'wheelMoves'} = {wheelMoves};
end

end









