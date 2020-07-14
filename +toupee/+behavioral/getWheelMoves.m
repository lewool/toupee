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
% sessions : int scalar OR int array OR char array OR cell array 
% (optional name-value pair)
%   Specific sessions for which to filter trials from, instead of from all
%   sessions. (Default: all)
%
% 'eventNames' : char array OR cell array (optional name-value pair)
%   Specified events (used as time-markers) around which wheel info is
%   returned. Each cell contains the name of an event. (Default: none)
%
% 'eventWindows' : numeric array OR cell array (optional name-value pair)
%   The time windows around each event specified in `eventNames`. Each cell
%   should contain a two-element numeric array containing the time window
%   to return: e.g. [-0.5 0.5] represents a half second before and after
%   the event. (Default: none)
%
% 'fs' : double scalar (optional name-value pair)
%   The sampling frequency of the data acquisition device that sampled the
%   wheel. (Default: 1000).
%
% 'gradFn' : function handle (optional name-value pair)
%   A function handle containing the function to use on the wheel position
%   data to compute the wheel velocity. (Default: The numerical gradient of
%   a single-pass, moving-average with a 5 ms window)
%
% 'wheelSpecs' : struct (optional name-value pair)
%   Contains specifications on wheel parameters. The possible fields and
%   their values are:
%       radius : double scalar
%           The radius of the wheel in m. (Default: .031)
%       res : int scalar
%           The total tick resolution of a full wheel revolution. 
%           (Default: 400)
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
%
%
% See Also:
% ---------
% toupee.behavioral.getTrials
% toupee.behavioral.getTrialBlocks
% toupee.behavioral.getEventTimes
%

%% Prerun checks.
% Imports
import toupee.behavioral.wheel.computeVelAcc
% Ensure input args are of proper type.
p = inputParser;
% must be char or cell array of chars or scalar int or array of ints
isValidSessions = @(x) isValidEvents(x) || all(mod(x, 1) == 0);
% must be char or cell array of chars
isValidEventNames = @(x) ischar(x) || (iscell(x) && ischar(x{1}));
% must be a two-elem vector with second elem greater than first, or a cell
% array of such vectors
isValidEventWindows = @(x) isnumeric(x) && numel(x) == 2 && x(2) > x(1)...
                      || (iscell(x) && all(cellfun(@(y)...
                                           isnumeric(y)...
                                           && numel(y) == 2 ...
                                           && y(2) > y(1), x)));
% must be an int scalar
isValidFs = @(x) isscalar(x) && mod(x, 1) == 0;
% must be a function handle
isValidGradFn = @(x) isa(x, 'function_handle');
% must be a struct with proper fields
isValidWheelSpecs =...
    @(x) isstruct(x)...
         && all(contains(fieldnames(x), ["radius", "res", "gain"]))...
         && (isscalar(x.radius) && isnumeric(x.radius) && x.radius > 0)...
         && (isscalar(x.res) && isnumeric(x.res) && x.res > 0)...
         && (isscalar(x.gain) && isnumeric(x.gain) && x.gain > 0);

addParameter(p, 'sessions', expInfo.('Row'), isValidSessions);
addParameter(p, 'eventNames', isValidEventNames, {});
addParameter(p, 'eventWindows', isValidEventWindows, {});
addParameter(p, 'fs', 1000, isValidFs);
addParameter(p, 'gradFn', @(x) gradient(movmean(x, 10)), isValidGradFn);
addParameter(p, 'wheelSpecs',...
             struct('radius', .031, 'res', 400, 'gain', 5),...
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
nE = numel(sessions);  % number of experiment sessions
% Go experiment-by-experiment, event-by-event, trial-by-trial.
for iE = 1:nE
    % Extract relevant data from this session.
    expRef = sessions{iE};  % session expRef
    block = expInfo.BlockFile{expRef};  % block data
    timeline = expInfo.TimelineFile{expRef};  % timeline data
    evts = block.events;  % events data
    nT = numel(evts.endTrialValues{1});  % number of completed trials
    % block wheel position in m
    x = block.inputs.wheelMMValues{1}' ./ 1000;  
    % block time in s
    tRaw = block.inputs.wheelMMTimes{1}';
    % Estimate wheel time in rig time from block time: estimate the delay 
    % time for timeline input -> block input (wheel data) as double the 
    % delay time as block output -> timeline output (reward data), and
    % subtract it.  @todo acknowledge `daqDelay` comes from `getEventTimes`
    tRaw = tRaw - (expInfo.behavioralData.daqDelay{iE} * 2);
    t = timeline.rawDAQTimestamps{1}(:);  % timeline time
    % Interpolate the wheel trace to evenly sample in rig time.
    x = interp1(tRaw, x, t, 'linear', 'extrap');
    % compute instantaneous velocity and acceleration for entire session
    [v, a] = computeVelAcc(x, t);
    % get trial durations (add daqDelay to block times for new trial and end trial)
    
    % Per specified event, per trial, get:
    % continuous position, continuous velocity, peak velocity, continuous
    % acceleration, peak acceleration, continuous direction, initial 
    % direction, final direction, and continuous movement classification.

end

end









