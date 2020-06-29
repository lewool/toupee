function [expInfo, mask, idx] =...
    getTrials(expInfo, conditions, colName, sessions)
% Filters trials from experiment session(s) that match some condition(s)
%
%
% Inputs:
% -------
% expInfo : table
%   A table containing relevant information and data variables (columns) 
%   for particular experiment sessions (rows).
%
% conditions : struct array
%   A struct whose fields specify which trials to select. The possible
%   fields and their values are:
%       reaction : 'early', 'preStimOn', 'preGoCue', 'normal', 'late'
%       action : 'right', 'left', 'timeout'
%       response : 'correct', 'incorrect'
%       outcome : 'rewarded', 'unrewarded'
%       repeatType : 'random', 'baited'
%       stimulusSide : 'right', 'left', 'zero'
%       contrasts : <array of values between [-1, 1]>
%       movementDir : 'right', 'left'
%       highRewardSide : 'right', 'left' 
%       highRewardSideConcordance : 'concordant', 'discordant'
%       pastResponse : 'correct', 'incorrect'
%       pastStimulusSide : 'right', 'left', 'zero'
%       pastMovementDir : 'right', 'left'
%       nTrialsPast : <integer greater than 0, less than nTrials>
%       circaBlockSwitch : 'beforeRight', 'afterRight', 'beforeLeft',
%                          'afterLeft'
%       nTrialsCirca : <integer greater than 0, less than nTrials>
%       whichTrials : <boolean mask array OR numeric int idx array>
%
% colName : char array OR cell array
%   User-defined name for the column in which the trials' info will be
%   saved into in `expInfo.behavioralData`.
%
%
% Outputs:
% --------
% expInfo : table
%   The updated `expInfo`.
%
% mask : cell array
%   Contains logical arrays in each cell, with the number of cells equal to
%   the number of sessions in `expInfo`, and the length of the logical 
%   arrays equal to the number of completed trials for the given session.
% 
% idx : cell array
%   Contains arrays in each cell, with the number of cells equal to
%   the number of sessions in `expInfo`, and each logical array within a 
%   cell containing the indices that matched the specified conditions for 
%   the specified experiment sessions.
%
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
% @todo add optionality to filter for specific, not all, sessions
%

% Do some checks on input args.
if ~ischar(colName) && (~iscell(colName) || ~ischar(colName{1}))...
    || ~isstruct(conditions)  % make sure input args are correct types
    error('toupee:meta:getTrials:badInput',...
          ['The "name" input arg must be a char array, and the '...
           '"conditions" input arg must be a struct']);
end
% See if all provided fieldnames are valid.
validNames = {'reaction', 'action', 'resoonse', 'outcome', 'repeatType', ....
              'stimulusSide', 'contrasts', 'movementDir', ...
              'highRewardSide', 'highRewardSideConcordance', ...
              'pastResponse', 'pastStimulusSide', 'pastMovementDir', ...
              'nTrialsBack', 'circaBlockSwitch', 'nTrialsCirca', ...
              'whichTrials'};
givenNames = fieldnames(conditions);          
if numel(find(strcmp(givenNames, validNames))) ~= numel(givenNames) 
    error('toupee:meta:getTrials:badInput',...
          ['At least one of the fields of the "conditions" input arg ',...
           'does not have a valid name.']);
end

% Get specified trials for each experiment session.
for e = 1:numel(expInfo)
    % Extract relevant data from this session.
    b = expInfo(e).block;  % block
    nt = numel(b.events.endTrialValues);  % number of trials
    % @todo wm = getWheelMoves;
    % Preassign the mask as all trials
    mask = true(1, nt);
    
    % Start selecting trials from conditions:
    
    % reaction @todo need wheel moves
    if isfield(conditions, 'reaction')  
        switch conditions.reaction
            case 'preStimOn'
            case 'preGoCue'
            case 'early'
            case 'normal'
            case 'late'
        end
    end
    
    % action
    if isfield(conditions, 'action')
        switch conditions.action
            case 'right'
                mask2 = b.events.responseValues(1:nt) == 1;
            case 'left'
                mask2 = b.events.responseValues(1:nt) == -1;
            case 'timeout'
                mask2 = b.events.responseValues(1:nt) == 0;
            otherwise
                error('toupee:behavioral:getTrials:badInput',...
                      ['''%s'' is not a valid value for "action". The ',...
                       'value must be ''right'', ''left'', or ',...
                       '''timeout''.'], conditions.action);
        end
        mask = mask & mask2;
    end
    
    % response
    if isfield(conditions, 'response')
        switch conditions.response
            case 'correct'
                mask2 = b.events.feedbackValues(1:nt);
            case 'incorrect'
                mask2 = ~b.events.feedbackValues(1:nt);
            otherwise
                error('toupee:behavioral:getTrials:badInput',...
                      ['''%s'' is not a valid value for "response". The ',...
                       'value must be ''correct'', or ''incorrect''.'],...
                       conditions.response);
        end
        mask = mask & mask2;
    end
    
    % outcome
    if isfield(conditions, 'outcome')
        rewardedTrials = b.events.rewardSizeValues(1:nt) > 0;
        switch conditions.outcome
            case 'rewarded'
                mask2 = rewardedTrials;
            case 'unrewarded'
                mask2 = ~rewardedTrials;
            otherwise
                error('toupee:behavioral:getTrials:badInput',...
                      ['''%s'' is not a valid value for "outcome". The ',...
                       'value must be ''rewarded'', or ''unrewarded''.'],...
                       conditions.outcome);
        end
        mask = mask & mask2;
    end
    
    % repeatType
    if isfield(conditions, 'repeatType')
        notRepeat = b.events.repeatNumValues(1:nt) == 1;
        switch conditions.repeatType
            case 'random'
                mask2 = notRepeat;
            case 'baited'
                mask2 = ~notRepeat;
            otherwise
                error('toupee:behavioral:getTrials:badInput',...
                      ['''%s'' is not a valid value for "repeatType". ',...
                       'The value must be ''random'', or ''baited''.'],...
                       conditions.repeatType);
        end
        mask = mask & mask2;
    end
    
    % stimulusSide
    if isfield(conditions, 'stimulusSide')
        switch conditions.stimulusSide
            case 'right'
                mask2 = b.events.stimulusSide(1:nt) == 1;
            case 'left'
                mask2 = b.events.stimulusSide(1:nt) == -1;
            case 'zero'
                mask2 = b.events.stimulusSide(1:nt) == 0;
            otherwise
                error('toupee:behavioral:getTrials:badInput',...
                      ['''%s'' is not a valid value for "stimulusSide". ',...
                       'The value must be ''right'', ''left'', or ',...
                       '''zero''.'], conditions.stimulusSide);
        end
        mask = mask & mask2;
    end
    
    % contrasts
    if isfield(conditions, 'contrasts')
        if any(conditions.contrasts < -1 | conditions.contrasts > 1)
            error('toupee:behavioral:getTrials:badInput',...
                  ['''%s'' is not a valid value for "contrasts". ',...
                   'The value(s) must be between [-1, 1].'],...
                   num2str(conditions.contrasts));
        end
        mask2 = any(b.events.contrastValues(1:nt)...
                    == conditions.contrasts(:), 1);
        mask = mask & mask2;
    end
    
    % movementDir @todo need wheel moves
    if isfield(conditions, 'movementDir')
        switch conditions.movement
            case 'right'
            case 'left'
        end
    end

    % highRewardSide
    if isfield(conditions, 'highRewardSide')
        switch conditions.highRewardSide
            case 'right'
                mask2 = b.events.highRewardSideValues(1:nt) == 1;
            case 'left'
                mask2 = b.events.highRewardSideValues(1:nt) == -1;
            otherwise
                error('toupee:behavioral:getTrials:badInput',...
                      ['''%s'' is not a valid value for ',...
                       '"highRewardSide". The value must be ''right'' or ',...
                       '''left''.'], conditions.highRewardSide);
        end
        mask = mask & mask2;
    end
    
    % highRewardSideConcordance
    if isfield(conditions, 'highRewardSideConcordance')
        % when stimulus appears on high reward side, or stimulus contrast
        % is 0 (we include this latter case for ease of plotting
        % psychometrics)
        concordance = sign(b.events.contrastValues(1:nt))...
                      == sign(b.events.highRewardSideValues(1:nt))...
                      |...
                      sign(b.events.contrastValues(1:nt)) == 0;
        switch conditions.highRewardSideConcordance
            case 'concordant'
                mask2 = concordance;
            case 'discordant'
                mask2 = ~concordance;
            otherwise
                error('toupee:behavioral:getTrials:badInput',...
                      ['''%s'' is not a valid value for ',...
                       '"highRewardSideConcordance". The value must be ',...
                       '''concordant'', or ''discordant''.'],...
                       conditions.highRewardSideConcordance);
        end
        mask = mask & mask2;
    end
    
    % pastResponse
    if isfield(conditions, 'pastResponse')
        switch conditions.pastResponse
            case 'left'
            case 'right'
        end
    end
    
    % pastStimulusSide
    if isfield(conditions, 'pastStimulusSide')
        switch conditions.pastStimulusSide
            case 'left'
            case 'right'
            case 'zero'
        end
    end
    
    % pastMovementDir
    if isfield(conditions, 'pastMovementDir')
        switch conditions.pastMovementDir
            case 'left'
            case 'right'
        end
    end
    
    % circaBlockSwitch
    if isfield(conditions, 'circaBlockSwitch')
        if ~isfield(conditions, 'nTrialsCirca')
            warning('toupee:behavioral:getTrials:missingInput',... 
                    ['Using "circaBlockSwitch" without specifying ',...
                     '"nTrialsCirca". Will set "nTrialsCirca" to 50 as ',...
                     'default.'])
            ntc = 50;
        else
            ntc = conditions.nTrialsCirca;
        end
        hrsv = b.events.highRewardSideValues(1:nt);
        switch conditions.switchBlocks
            case 'beforeRight'
                switchIdx = find(diff(hrsv) > 0) + 1; 
                maskIdx = arrayfun(@(x) [(x - ntc):1:x], switchIdx,...
                                   'UniformOutput', false); %#ok<*NBRAK>
                maskIdx = horzcat(maskIdx{:});
            case 'afterRight'
                switchIdx = find(diff(hrsv) > 0) + 1; 
                maskIdx = arrayfun(@(x) [(x + ntc):1:x], switchIdx,...
                                   'UniformOutput', false);
                maskIdx = horzcat(maskIdx{:});
            case 'beforeLeft'
                switchIdx = find(diff(hrsv) < 0) + 1; 
                maskIdx = arrayfun(@(x) [(x - ntc):1:x], switchIdx,...
                                   'UniformOutput', false);
                maskIdx = horzcat(maskIdx{:});
            case 'afterLeft'
                switchIdx = find(diff(hrsv) < 0) + 1; 
                maskIdx = arrayfun(@(x) [(x + ntc):1:x], switchIdx,...
                                   'UniformOutput', false);
                maskIdx = horzcat(maskIdx{:});
            otherwise
                error('toupee:behavioral:getTrials:badInput',...
                      ['''%s'' is not a valid value for ',...
                       '"circaBlockSwitch". The value must be ',...
                       '''beforeRight'', ''afterRight'', ''beforeLeft'', ',...
                       'or ''afterLeft''.'], conditions.circaBlockSwitch);
        end
        % Remove too high and too low values in mask.
        if any(maskIdx < 0) || any(maskIdx > nt)
            warning(['In "circaBlockSwitch", some indices circa the block ',...
                    'switch either exceeded the number of trials, or ',...
                    'deceeded 0. These trials will be removed.']);
            maskIdx(maskIdx < 0) = [];
            maskIdx(maskIdx > nt) = [];
        end
        mask2 = false(1, nt);
        mask2(maskIdx) = 1;
        mask = mask & mask2;
    end
    
    % whichTrials
    if isfield(conditions, 'whichTrials')
        if ~isnumeric(conditions.whichTrials)...
           || ~islogical(conditions.whichTrials)
            error('toupee:behavioral:getTrials:badInput',...
                      ['"whichTrials" must be a logical array of length ',...
                       'nTrials, or it must be a numeric array ',...
                       'containing the indices of the trials to pass the ',...
                       'filter.']);
        elseif isnumeric(conditions.whichTrials)
            mask2 = false(1, nt);
            mask2(whichTrials) = 1;
        elseif islogical(conditions.whichTrials)
            mask2 = conditions.whichTrials;
        end
        mask = mask & mask2;
    end
    
    % Finalize `idx`, and add to `expInfo`
    idx = find(mask);
    maskName = strcat(colName, 'Mask');
    idxName = strcat(colName, 'Idx');
    expInfo(e).behavioralData.trials.(maskName) = mask;
    expInfo(e).behavioralData.trials.(idxName) = idx;
    
end

end