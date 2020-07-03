function [expInfo, trialBlocks] = getTrialBlocks(expInfo, sessions)
% Separates trials into the block they came from
%
%
% Inputs:
% -------
% expInfo : table
%   A table containing relevant information and data variables (columns) 
%   for particular experiment sessions (rows).
%
% sessions : int scalar OR int array OR char array OR cell array (optional)
%   Specific sessions for which to filter trials from, instead of from all
%   sessions.
%
% Outputs:
% --------
% expInfo : table
%   The updated `expInfo`.
%
% trialBlocks : cell array
%   Contains `numel(expInfo)` number of cells, where each of these cells is
%   a [nBlocks x 2] cell array, where each row representes a block: the
%   first element in each row contains an array of the trial indices of
%   that block, and the second element contains info on the block type 
%   (left high reward side == -1, right high reward side == 1)
%
%
% Examples:
% ---------
% 1) For a single session, get the trial blocks.
%   deats = {'LEW031', '2020-02-03', 1};
%   files = 'block';
%   expInfo = toupee.meta.processExperiment(deats, files);
%   [expInfo, trialBlocks] = toupee.behavioral.getTrialBlocks(expInfo);
%
% 2) For multiple sessions, get each session's trial blocks.
%   deats = {{'LEW031', '2020-02-03', 1},...
%            {'LEW032', '2020-02-28', 1, [1, 2]}};
%   files = 'block';
%   expInfo = toupee.meta.processExperiment(deats, files);
%   expInfo = toupee.behavioral.getTrialBlocks(expInfo);
%
%
% See Also:
% ---------
% toupee.behavioral.getTrials
%
% @todo adopt for high probability side experiments
%

%% Prerun checks.
% Get sessions info
if nargin < 2  % then use all sessions
    sessions = expInfo.('Row');
else
    sessions = expInfo.Row(sessions);
end

%% Get trial blocks.
nE = size(expInfo, 1);  % number of experiment sessions
trialBlocks = cell(1, nE);  % initialize `blocks`
% Get trial blocks for each session
for iE = 1:nE
    % Extract relevant data from this session.
    expRef = sessions{iE};  % session expRef
    b = expInfo.BlockFile{expRef};  % block data
    e = b.events;  % events data
    nT = numel(e.endTrialValues{1});  % number of trials
    try
        req = e.highRewardSideValues{1}(1:nT);  % contains block type info
    catch ex
        warning(ex.identifier,...
                strcat(ex.message(1:(end-1)),[...
                ' in the saved events for %s; cannot compute ''%s''. ',...
                'Continuing to the next session.']), expRef, mfilename);
    end
    switchIdx = find(diff(req) ~= 0) + 1;  % trial index at block switch
    startIdx = [1, switchIdx];  % trial idx at block start
    endIdx = [switchIdx - 1, nT];  % trial index at block end
    curTrialBlocks = arrayfun(@(x, y) [x:1:y], startIdx, endIdx,...
                              'uni', 0)';  %#ok<*NBRAK>
    curTrialBlocks(:,2) = num2cell(req(startIdx));
    % Assign to `expInfo` and `trialBlocks`.
    expInfo.behavioralData{expRef, 'trialBlocks'} = {curTrialBlocks};
    trialBlocks{iE} = curTrialBlocks;
end
