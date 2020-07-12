function [expInfo, mask] =...
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
%       choice : 'right', 'left', 'timeout'
%       response : 'correct', 'incorrect'
%       outcome : 'rewarded', 'unrewarded'
%       repeatType : 'random', 'baited'
%       stimulusSide : 'right', 'left', 'zero'
%       contrasts : <array of values between [-1, 1]>
%       movementDir : 'right', 'left'
%       highRewardSide : 'right', 'left'
%       highRewardSideConcordance : 'concordant', 'discordant'
%       quiescent : 'true', 'false'
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
% mask : cell array
%   Contains logical arrays in each cell, with the number of cells equal to
%   the number of specified sessions, and the length of the logical 
%   arrays equal to the number of completed trials for the given session.
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

%% Prerun checks.
% Import all functions in `+misc`.
import toupee.misc.*
% Ensure input args are of proper type.
if ~ischar(colName) && (~iscell(colName) || ~ischar(colName{1}))...
    || ~isstruct(conditions)
    error('toupee:behavioral:getTrials:badInput',...
          ['The "colName" input arg must be a char array, and the '...
           '"conditions" input arg must be a struct']);
end
% See if all provided fieldnames are valid.
validNames = {...
    'reaction', 'choice', 'response', 'outcome', 'repeatType',...
    'stimulusSide', 'contrasts', 'movementDir', 'quiescent',...
    'highRewardSide', 'highRewardSideConcordance', 'pastResponse',...
    'pastStimulusSide', 'pastMovementDir', 'nTrialsBack',...
    'circaBlockSwitch', 'nTrialsCirca', 'whichTrials'};
givenNames = fieldnames(conditions);
matchedNames = cellfun(@(x) find(strcmpi(x, validNames)), givenNames,...
                       'uni', 0);
% If can't find a match for one of `conditions`' fieldnames, throw error.
badConds = find(cell2mat(cellfun(@(x) isempty(x), matchedNames,...
                         'uni', 0)));
if ~isempty(badConds)
    error('toupee:behavioral:getTrials:badInput', strcat(['The following ',...
          'provided fields of the "conditions" input arg do not have ',...
          'valid names: '], sprintf(' %s', givenNames{badConds})));
end
% Get sessions info
if nargin < 4  % then use all sessions
    sessions = expInfo.('Row');
else  % use specified sessions
    sessions = expInfo.Row(sessions);
end
%% Filter trials.
% initialize `maskCell`; each cell will contain mask for that session
nE = numel(sessions);  % number of experiment sessions
maskCell = cell(1, nE);  
% Get specified trials for each specified experiment session.
for iE = 1:nE
    % Extract relevant data from this session.
    expRef = sessions{iE};  % session expRef
    block = expInfo.BlockFile{expRef};  % block data
    evts = block.events;  % events data
    nT = numel(evts.endTrialValues{1});  % number of completed trials
    % @todo wm = getWheelMoves;
    % Preassign the mask to return all trials
    mask = true(1, nT);
    
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
    
    % choice
    if isfield(conditions, 'choice')
        try
            req = evts.responseValues{1}(1:nT);  % required data for condition
        catch ex
            warning(ex.identifier,...
                    strcat(ex.message(1:(end-1)),[...
                    ' in the saved events for %s; cannot compute the ',...
                    '''choice'' condition. Continuing to the next ',...
                    'session.']), expRef);
            continue
        end
        switch conditions.choice
            case 'right'
                mask2 = req == 1;
            case 'left'
                mask2 = req == -1;
            case 'timeout'
                mask2 = req == 0;
            otherwise
                error('toupee:behavioral:getTrials:badInput',...
                      ['''%s'' is not a valid value for "choice". The ',...
                       'value must be ''right'', ''left'', or ',...
                       '''timeout''.'], conditions.choice);
        end
        mask = mask & mask2;
    end
    
    % response
    if isfield(conditions, 'response')
        try
            req = evts.feedbackValues{1}(1:nT);
        catch ex
            warning(ex.identifier,...
                    strcat(ex.message(1:(end-1)),[...
                    ' in the saved events for %s; cannot compute the ',...
                    '''response'' condition. Continuing to the next ',...
                    'session.']), expRef);
            continue
        end
        switch conditions.response
            case 'correct'
                mask2 = req;
            case 'incorrect'
                mask2 = ~req;
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
        try
            req = evts.rewardSizeValues{1}(1:nT);
        catch ex
            warning(ex.identifier,...
                    strcat(ex.message(1:(end-1)),[...
                    ' in the saved events for %s; cannot compute the ',...
                    '''outcome'' condition. Continuing to the next ',...
                    'session.']), expRef);
            continue
        end
        rewardedTrials = req > 0;
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
        try
            req = evts.repeatNumValues{1}(1:nT);
        catch ex
            warning(ex.identifier,...
                    strcat(ex.message(1:(end-1)),[...
                    ' in the saved events for %s; cannot compute the ',...
                    '''repeatType'' condition. Continuing to the next ',...
                    'session.']), expRef);
            continue
        end
        notRepeat = req == 1;
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
        try
            req = sign(evts.contrastValues{1}(1:nT));
        catch ex
            warning(ex.identifier,...
                    strcat(ex.message(1:(end-1)),[...
                    ' in the saved events for %s; cannot compute the ',...
                    '''stimulusSide'' condition. Continuing to the next ',...
                    'session.']), expRef);
            continue
        end
        switch conditions.stimulusSide
            case 'right'
                mask2 = req == 1;
            case 'left'
                mask2 = req == -1;
            case 'zero'
                mask2 = req == 0;
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
        try
            req = evts.contrastValues{1}(1:nT);
        catch ex
            warning(ex.identifier,...
                    strcat(ex.message(1:(end-1)),[...
                    ' in the saved events for %s; cannot compute the ',...
                    '''contrasts'' condition. Continuing to the next ',...
                    'session.']), expRef);
            continue
        end
        if any(conditions.contrasts < -1 | conditions.contrasts > 1)
            error('toupee:behavioral:getTrials:badInput',...
                  ['''%s'' is not a valid value for "contrasts". ',...
                   'The value(s) must be between [-1, 1].'],...
                   num2str(conditions.contrasts));
        end
        mask2 = any(req == conditions.contrasts(:));
        mask = mask & mask2;
    end
    
    % movementDir @todo need wheel moves
    if isfield(conditions, 'movementDir')
        switch conditions.movementDir
            case 'right'
            case 'left'
        end
    end
    
    % quiescent @todo need wheel moves
    if isfield(conditions, 'quiescent')
        try
        catch ex
            warning(ex.identifier,...
                    strcat(ex.message(1:(end-1)),[...
                    ' in the saved events for %s; cannot compute the ',...
                    '''quiescent'' condition. Continuing ',...
                    'to the next session.']), expRef);
            continue
        end
        % when stimulus appears on high reward side, or stimulus contrast
        % is 0 (we include this latter case for ease of plotting
        % psychometrics)
        switch conditions.quiescent
            case 'true'
            case 'false'
            otherwise
                error('toupee:behavioral:getTrials:badInput',...
                      ['''%s'' is not a valid value for ',...
                       '"quiescent". The value must be '...
                       '''true'', or ''false''.'],...
                       conditions.quiescent);
        end
        mask = mask & mask2;
    end

    % highRewardSide
    if isfield(conditions, 'highRewardSide')
        try
            req = evts.highRewardSideValues{1}(1:nT);
        catch ex
            warning(ex.identifier,...
                    strcat(ex.message(1:(end-1)),[...
                    ' in the saved events for %s; cannot compute the ',...
                    '''highRewardSide'' condition. Continuing to the ',...
                    'next session.']), expRef);
            continue
        end
        switch conditions.highRewardSide
            case 'right'
                mask2 = req == 1;
            case 'left'
                mask2 = req == -1;
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
        try
            req = evts.contrastValues{1}(1:nT);
            req2 = evts.highRewardSideValues{1}(1:nT);
        catch ex
            warning(ex.identifier,...
                    strcat(ex.message(1:(end-1)),[...
                    ' in the saved events for %s; cannot compute the ',...
                    '''highRewardSideConcordance'' condition. Continuing ',...
                    'to the next session.']), expRef);
            continue
        end
        % when stimulus appears on high reward side, or stimulus contrast
        % is 0 (we include this latter case for ease of plotting
        % psychometrics)
        concordance = sign(req) == sign(req2) | sign(req) == 0;
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
        try
            req = evts.highRewardSideValues{1}(1:nT);
        catch ex
            warning(ex.identifier,...
                    strcat(ex.message(1:(end-1)),[...
                    ' in the saved events for %s; cannot compute the ',...
                    '''circaBlockSwitch'' condition. Continuing ',...
                    'to the next session.']), expRef);
              continue
        end
        switch conditions.circaBlockSwitch
            case 'beforeRight'
                switchIdx = find(diff(req) > 0) + 1; 
                maskIdx = arrayfun(@(x) [(x - ntc):1:x], switchIdx,...
                                   'uni', 0); %#ok<*NBRAK>
            case 'afterRight'
                switchIdx = find(diff(req) > 0) + 1; 
                maskIdx = arrayfun(@(x) [x:1:(x + ntc)], switchIdx,...
                                   'uni', 0);
            case 'beforeLeft'
                switchIdx = find(diff(req) < 0) + 1; 
                maskIdx = arrayfun(@(x) [(x - ntc):1:x], switchIdx,...
                                   'uni', 0);
            case 'afterLeft'
                switchIdx = find(diff(req) < 0) + 1; 
                maskIdx = arrayfun(@(x) [x:1:(x + ntc)], switchIdx,...
                                   'uni', 0);
            otherwise
                error('toupee:behavioral:getTrials:badInput',...
                      ['''%s'' is not a valid value for ',...
                       '"circaBlockSwitch". The value must be ',...
                       '''beforeRight'', ''afterRight'', ''beforeLeft'', ',...
                       'or ''afterLeft''.'], conditions.circaBlockSwitch);
        end
        % Remove too high and too low values in mask.
        maskIdx = horzcat(maskIdx{:});
        if any(maskIdx < 0) || any(maskIdx > nT)
            warning(['In "circaBlockSwitch", some indices circa the block ',...
                    'switch either exceeded the number of trials, or ',...
                    'deceeded 0. These trials will be removed.']);
            maskIdx(maskIdx < 0) = [];
            maskIdx(maskIdx > nT) = [];
        end
        mask2 = false(1, nT);
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
            mask2 = false(1, nT);
            mask2(whichTrials) = 1;
        elseif islogical(conditions.whichTrials)
            mask2 = conditions.whichTrials;
        end
        mask = mask & mask2;
    end
    
    % Finalize `mask` and add to `expInfo`
    if ~contains(colName, 'trials', 'IgnoreCase', true)
        maskName = strcat(colName, 'Trials');
    end 
    expInfo.behavioralData{expRef, (maskName)} = {mask};
    maskCell{iE} = mask;
end

mask = maskCell;  % assign `maskCell` to `mask` to return output

end