function [expInfo] = ...
    processExperiment(details, behavioralSpecs, neuralSpecs)
% Gets experiment information and data for given session(s)
% 
% Inputs:
% -------
% details : cell array
%   The session(s) details. Use a nested cell for each session, with four
%   elements in each nested cell: subject name (char array), date (char 
%   array in datestr format), session (int scalar), and series (int array).
% behavioralSpecs : cell array of name-value pairs
%   An array of function handles of functions in `+toupee\+behavioral\`.
%   These functions will be called in order to set up the 'behavioralData'
%   field of the returned `expInfo` struct.
% neuralFns : cell array of name-value pairs
%   An array of function handles of functions in `+toupee\+neural\`. These
%   functions will be called in order to set up the 'neuralData' field of
%   the returned `expInfo` struct.
% 
% Outputs:
% --------
% expInfo : struct
%   A struct array with each element containing fields with information for
%   a particular session. The fields for each struct element are:
%   'subject', 'expDate', 'expNum', 'expSeries', 'block', 'timeline',
%   'numPlanes', 'numChannels', 'behavioralData', 'neuralData'.
%
% Examples:
% ---------
% 1) Return experiment info for a single session from a single subject.
%   details = {{'LEW031', '2020-02-03', 1, 1}};
%   expInfo = toupee.meta.processExperiment(details);
% 2) Return experiment info for multiple sessions from multiple subjects.
%   details = {{'LEW031', '2020-02-03', 1, 1},... 
%              {'LEW015', '2019-03-21', 1, 1},...
%              {'LEW005', '2018-06-10', 2, [2 3]}};
%   expInfo = toupee.meta.processExperiment(details);

% import all other functions in this subpackage, and `iif`
import toupee.meta.*
import toupee.misc.iif

% Do some checks on input args

% For each experiment session...
for e = 1:length(details)
    % Initialize expInfo struct.
    expInfo(e) = struct(...
        'subject', details{e}{1},...
        'expDate', details{e}{2},...
        'expNum', details{e}{3},...
        'expSeries', details{e}{4},...
        'block', [],...
        'timeline', [],...
        'numPlanes', [],...
        'numChannels', [],...
        'behavioralData', [],...
        'neuralData', []);  %#ok<AGROW>
    % Load exp data (block + timeline files)
    subject = expInfo(e).subject;
    expDate = expInfo(e).expDate;
    expNum = expInfo(e).expNum;
    expRef = constructExpRef(subject, expDate, expNum);
    serverPaths = getPaths().server;
    for s = 1:length(serverPaths)  % check each server for exp data
        p = serverPaths{s};
        blockFilePath = fullfile(p, subject, expDate, num2str(expNum),...
                                 strcat(expRef, '_Block.mat'));
        timelineFilePath = fullfile(p, subject, expDate, num2str(expNum),...
                                    strcat(expRef, '_Timeline.mat'));
        % Load block file
        if isempty(expInfo.block) && isfile(blockFilePath)
            expInfo.block = load(blockFilePath);
        end
        % Load timeline file
        if isempty(expInfo.timeline) && isfile(timelineFilePath)
            expInfo.timeline = load(timelineFilePath);
        end
    end
end

%% load behavioral data

%% load neural data
 
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

expInfo.behavioralData = behavioralData;
expInfo.neuralData = neuralData;

end

