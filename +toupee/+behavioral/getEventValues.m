function [expInfo, eventValues] = getEventValues(expInfo, names)
% Gets signals' event values from a block file in an `expInfo` struct
%
%
% Inputs:
% -------
% expInfo : struct array
%   A struct containing relevant information and data for particular
%   experiment sessions.
%
% names : char array OR cell array
%   The names of the signals saved in `block.events` of a block file. Use a
%   cell array if returning the values for 2 or more events. Use a nested
%   nested cell array to get event values for multiple sessions.
%
%
% Outputs:
% --------
% expInfo : struct array
%   A struct containing relevant information and data for particular
%   experiment sessions.
%
% eventValues : struct array
%   An array where each row represents a different event, and the 2nd
%   column in each row contains the signals values for that event in a cell
%   array, where each cell in that array corresponds to an experiment
%   session.
%
%
% Examples:
% ---------
% 1) For a single session: get the contrast, action, and response values
% for each trial.
%   details = {'LEW031', '2020-02-03', 1};
%   files = {'block'};
%   expInfo = toupee.meta.processExperiment(details, files);
%   names = 
%   [expInfo, eventValues] =...
%       toupee.behavioral.getEventValues(expInfo, names);
%
%
% 2) For multiple sessions: find 1) trials where stim was presented on the
% 'high-reward' side, 2) correct trials where stim was presented on the
% 'high-reward' side, and 3) correct trials where stim was presented on the
% 'low-reward' side.
%   details = {{'LEW031', '2020-02-03', 1},...
%              {'LEW032', '2020-03-12', 1, [1, 2]}};
%   files = {'timeline', 'block'};
%   expInfo = toupee.meta.processExperiment(details, files);
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
%       toupee.behavioral.getTrials(expInfo, conditions1, ...
%                                   'discordantHighCorrect');
% 
% 3) For multiple sessions: for the first session find trials where stim
% was presented on the 'high-reward' side, and for the second session find
% correct response trials.
%   conditions(1) = struct('highRewardSideConcordance', 'concordant');
%   conditions(2) = struct('response', 'correct');
%
% @todo add explanations for each field in `conditions`
% @todo add more documentation
% @todo add code for 'past' conditions
end