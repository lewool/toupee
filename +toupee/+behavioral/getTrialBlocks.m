function [expInfo, blocks] = getTrialBlocks(expInfo)
% Separates trials into the block they came from
%
%
% Inputs:
% -------
% expInfo : struct array
%   A struct containing relevant information and data for particular
%   experiment sessions.
%
%
% Outputs:
% --------
% expInfo : struct array
%   A struct containing relevant information and data for particular
%   experiment sessions.
%
% blocks : cell array
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
%   details = {'LEW031', '2020-02-03', 1}
%   files = 'block';
%   expInfo = toupee.meta.processExperiment(details, files);
%   expInfo = toupee.behavioral.getTrialBlocks(expInfo);
%
% 2) For multiple sessions, get each session's trial blocks.
%   details = {{'LEW031', '2020-02-03', 1},... 
%              {'LEW037', '2020-03-13', 1},...
%              {'LEW005', '2018-06-10', 2, [2 3]}};
%   files = 'block';
%   expInfo = toupee.meta.processExperiment(details, files);
%   expInfo = toupee.behavioral.getTrialBlocks(expInfo);
%
%
% See Also:
% ---------
% toupee.behavioral.filterTrials
%

ne = numel(expInfo);  % number of experiment sessions
blocks = cell(1, ne);  % initialize `blocks`
% Get trial blocks for each session
for e = 1:ne
    % Extract relevant data from this session.
    b = expInfo(e).block;  % Rigbox block file
    nt = numel(b.events.endTrialValues);  % number of trials
    hrsv = b.events.highRewardSideValues(1:nt);  % block type
    switchIdx = find(diff(hrsv) ~= 0) + 1;  % trial index at block switch
    startIdx = [1, switchIdx];  % trial idx at block start
    endIdx = [switchIdx - 1, nt];  % trial index at block end
    curTrialBlocks = arrayfun(@(x, y) [x:1:y], startIdx, endIdx,...
                              'UniformOutput', false)';  %#ok<*NBRAK>
    curTrialBlocks(:,2) = num2cell(hrsv(startIdx));
    % Assign.
    expInfo(e).behavioralData.trials.blocks = curTrialBlocks;
    blocks{e} = curTrialBlocks;
end
