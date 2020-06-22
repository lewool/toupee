function [expInfo, mask, idx] = selectTrials(expInfo, name, conditions)
% Selects trials from an experiment session that match some condition(s).
%
%
% Inputs:
% -------
% expInfo : struct array
%   A struct containing relevant information and data for particular
%   experiment sessions.
%
% name : cell array
%   User-defined name for the field in which the trials info will be saved
%   into `expInfo.behavioralData.trials.`
%
% conditions : struct array
%   A struct whose fields specify which trials to select. The possible
%   fields and their values are:
%       response : 'correct', 'incorrect'
%       reaction : 'early', 'normal'
%       outcome : 'rewarded', 'unrewarded'
%       repeatType : 'random', 'baited'
%       contrasts : [array of values between 0-1]
%       movement : 'ccw', 'cw'
%       highRewardSide : 'left', 'right'
%       concordant : 'high', 'low', 'discordant'
%       pastStimulus : 'right', 'left', 'zero'
%       pastMovement : 'ccw', 'cw'
%       pastResponse : 'correct', 'incorrect'
%       trialsBack : [any integer]
%       switchBlocks : 'beforeLeft', 'afterLeft', 'beforeRight',
%                      'afterRight'
%       whichTrials : [boolean mask array OR numeric int idx array]
%
%
% Outputs:
% --------
% expInfo : struct
%   A struct containing relevant information and data for particular
%   experiment sessions.
%
% mask : cell array
%   Contains logical arrays in each cell, with the length of the logical 
%   arrays equal to the number of trials for the given experiment session.
% 
% idx : cell array
%   Contains arrays in each cell, with each array containing the trial
%   indices that matched the specified conditions for the specified
%   experiment sessions.
%
% Examples:
% ---------
% 1) For a single session: find trials where stim was presented on the
%   'high-reward' side:
%
%
% 2) For multiple sessions: find 1) trials where stim was presented on the
% 'high-reward' side, 2) correct trials where stim was presented on the
% 'high-reward' side, and 3) correct trials where stim was presented on the
% 'low-reward' side.
%   details = {{'LEW031', '2020-02-03', 1},...
%              {'LEW032', '2020-03-12', 1, [1, 2]}};
%   specs = {'timeline', 'block'};
%   expInfo = processExperiment(details, specs);
%   conditions1 = struct('concordant', 'high');
%   conditions2 = struct('concordant', 'high', 'response', 'correct');
%   conditions3 = struct('concordant', 'low', 'response', 'correct');
% 
% 3) For multiple sessions, for the first session, find trials that
% occurred on the 'high-reward' side, and for the second session, find
% correct response trials.
%
%
% @todo add explanations for each field in `conditions`
%

% Do some checks on input args
if ~iscell(name) || ~ischar(name{1}) || ~isstruct(conditions)
    error('toupee:meta:selectTrials:badInput',...
          ['The "name" input arg should be a char array, and the '...
           '"conditions" input arg should be a struct']);
end

% Get all the fieldnames in `conditions`
fields = fieldnames(conditions);

for e = 1:numel(expInfo)
    
end

end