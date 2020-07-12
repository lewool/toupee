function [expInfo, wheelMoves] = getWheelMoves(expInfo, sessions)
% Gets per-trial wheel information
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
%
% Outputs:
% --------
% expInfo : table
%   The updated `expInfo`.
%
% wheelMoves : table
%   Contains continuous position, continuous velocity, peak velocity, 
%   final movement direction, movement classifications per binned time, and 
%   movement relative to specified signals events for each trial in a 
%   session. The number of cells is equal to the number of specified
%   sessions.
%
%
% Examples:
% ---------
% 1) For a single session: find trials where stim was presented on the
%   'high-reward' side:
%   deats = {'LEW031', '2020-02-03', 1};
%   files = {'timeline', 'block'};
%   expInfo = toupee.meta.processExperiment(deats, files);
%   conditions = struct('highRewardSideConcordance', 'concordant');
%   [expInfo, mask] =...
%       toupee.behavioral.getTrials(expInfo, conditions,...
%                                   'concordantHigh');
%
% 2) For multiple sessions: find 1) trials where stim was presented on the
% 'high-reward' side, 2) correct trials where stim was presented on the
% 'high-reward' side, and 3) correct trials where stim was presented on the
% 'low-reward' side.
%   deats = {{'LEW031', '2020-02-03', 1},...
%            {'LEW032', '2020-02-28', 1, [1, 2]}};
%   files = {'timeline', 'block'};
%   expInfo = toupee.meta.processExperiment(deats, files);
%   conditions1 = struct('highRewardSideConcordance', 'concordant');
%   conditions2 = struct('highRewardSideConcordance', 'concordant',...
%                        'response', 'correct');
%   conditions3 = struct('highRewardSideConcordance', 'discordant',...
%                        'response', 'correct');
%   expInfo =...
%       toupee.behavioral.getTrials(expInfo, conditions1,...
%                                   'concordantHigh');
%   expInfo =...
%       toupee.behavioral.getTrials(expInfo, conditions2, ...
%                                   'concordantHighCorrect');
%   expInfo =...
%       toupee.behavioral.getTrials(expInfo, conditions3, ...
%                                   'discordantHighCorrect');
% 
% 3) For multiple sessions: for the first session find trials where stim
% was presented on the 'high-reward' side, and for the second session find
% correct response trials.
%   deats = {{'LEW031', '2020-02-03', 1},...
%            {'LEW032', '2020-02-28', 1, [1, 2]}};
%   files = {'timeline', 'block'};
%   expInfo = toupee.meta.processExperiment(deats, files);
%   conditions1 = struct('highRewardSideConcordance', 'concordant');
%   session1 = 1;
%   conditions2 = struct('response', 'correct');
%   session2 = 2;
%   expInfo =...
%       toupee.behavioral.getTrials(expInfo, conditions1,...
%                                   'concordantHigh', session1);
%   expInfo =...
%       toupee.behavioral.getTrials(expInfo, conditions2,...
%                                   'correct', session2);
%
%
% See Also:
% ---------
% toupee.behavioral.getWheelMoves
% toupee.behavioral.getTrialBlocks
% toupee.behavioral.getEventTimes
%
% @todo add explanations for each field in `conditions`
% @todo add more documentation
% @todo add code for 'past' conditions
%
end