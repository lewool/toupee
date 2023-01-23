function [condLogical, condIdx] = selectCondition(expInfo, contrast, behavioralData, trialConditions)
% A wee function for picking the trial types you want to look at when analyzing your cells
% INPUTS: block struct, a contrast/side (this can be a vector), and trialConditions, a struct
% that specifies the types of trials you want to include. Both 'expInfo' and eventTimes' can 
% be 1xn structs containing exp data from multiple days for concatenation
% OUTPUT: an index of all the relevant trials from the experiment(s)

% 8 March 2018 Written by LEW
% 23 July 2018 Updated to accept multiple blocks/eventTimes
% 19 November 2018 Added optional trial-history selections
% 13 May 2019 Added rewarded/unrewarded outcome selection (likelihood experiments)
% 20 Nov 2019 Changed so that trial conditions are input from a single
% struct, initialized in 'initTrialConditions.m'
% 27 Mar 2020 Updated to take expInfo directly, fixed concatenation

%% unpack trialConditions struct
repeatType = trialConditions.repeatType;
movementDir = trialConditions.movementDir;
movementTime = trialConditions.movementTime;
highRewardSide = trialConditions.highRewardSide;
preStimMovement = trialConditions.preStimMovement;
responseType = trialConditions.responseType;
rewardOutcome = trialConditions.rewardOutcome;
pastStimulus = trialConditions.pastStimulus;
pastMovementDir = trialConditions.pastMovementDir;
pastMovementTime = trialConditions.pastMovementTime;
pastResponseType = trialConditions.pastResponseType;
trialsBack = trialConditions.trialsBack;
switchBlocks = trialConditions.switchBlocks;
trimBlocks = trialConditions.trimBlocks;
whichTrials = trialConditions.whichTrials;
specificRTs = trialConditions.specificRTs;

%% 

% get length of expInfo blocks
numExps = length(expInfo);
condLogical = [];
for iExp = 1:numExps
    
    % extract data from single experiment
    b = expInfo(iExp).block;
    et = behavioralData(iExp).eventTimes;
    wm = behavioralData(iExp).wheelMoves;

    nt = numel(b.events.endTrialTimes);
    contrastVal = b.events.contrastValues(1:nt); 
    
    firstMoveDirs = zeros(1,nt);
    firstMoveDirs((~isnan(wm.epochs(3).onsetTimes))) = wm.epochs(3).moveDir(~isnan(wm.epochs(3).onsetTimes));
    firstMoveDirs((~isnan(wm.epochs(2).onsetTimes))) = wm.epochs(2).moveDir(~isnan(wm.epochs(2).onsetTimes));
    
%%%%%%%%%%%% ASSIGNMENT SEGMENT (FROM TRIALCONDITIONS INPUT) %%%%%%%%%%%%%%
    
    switch repeatType
        case 'random' %randomly presented trials
            idxRepeat = b.events.repeatNumValues == 1;
            idxRepeat = idxRepeat(1:nt);
        case 'baited' %trials repeated after an incorrect response
            idxRepeat = b.events.repeatNumValues > 1;
            idxRepeat = idxRepeat(1:nt);
        case 'all'
            idxRepeat = true(1,nt);
        otherwise
            error('Choose a valid repeatType: "random", "baited", or "all"');
    end

    switch movementDir
        case 'ccw' %moved the stimulus from right to left
            if strcmp(movementTime,'early')
                idxDirection = wm.epochs(2).moveDir == 1;
%                 idxDirection = b.events.responseValues == 1;
            elseif strcmp(movementTime,'late')
%                 idxDirection = wm.epochs(3).moveDir == 1;
                idxDirection = b.events.responseValues == 1;
            elseif strcmp(movementTime,'all')
%                 idxDirection = firstMoveDirs == 1;
                idxDirection = b.events.responseValues == 1;
            end
            idxDirection = idxDirection(1:nt);
        case 'cw' %moved the stimulus from left to right
            if strcmp(movementTime,'early')
                idxDirection = wm.epochs(2).moveDir == -1;
%                 idxDirection = b.events.responseValues == -1;
            elseif strcmp(movementTime,'late')
%                 idxDirection = wm.epochs(3).moveDir == -1;
                idxDirection = b.events.responseValues == -1;
            elseif strcmp(movementTime,'all')
%                 idxDirection = firstMoveDirs == -1;
                idxDirection = b.events.responseValues == -1;
            end
            idxDirection = idxDirection(1:nt);
        case 'all'
            idxDirection = true(1,nt);
        otherwise
            error('Choose a valid movementType: "ccw", "cw", or "all"');
    end

    switch movementTime
        case 'early'
            idxMovement = wm.epochs(2).isMoving;
            idxMovement = idxMovement(1:nt);
        case 'late'
            idxMovement = ~wm.epochs(2).isMoving .* wm.epochs(3).isMoving;
        case 'all'
            idxMovement = true(1,nt);
        otherwise
            error('Choose a valid movementTime: "early", "late", or "all"');
    end

    switch highRewardSide
        case 'left' %left-stimulus trials had larger value/prob
            try
                idxRewardSide = b.events.highRewardSideValues   == -1;
                idxRewardSide = idxRewardSide(1:nt);
            catch
                try
                    idxRewardSide = b.events.likelyRewardSideValues   == -1;
                    idxRewardSide = idxRewardSide(1:nt);
                catch
                    warning('No biased-side parameter found; including all trials by default')
                    idxRewardSide = true(1,nt);
                end
            end

        case 'right' %right-stimulus trials had larger value/prob
            try
                idxRewardSide = b.events.highRewardSideValues   == 1;
                idxRewardSide = idxRewardSide(1:nt);
            catch
                try
                    idxRewardSide = b.events.likelyRewardSideValues   == 1;
                    idxRewardSide = idxRewardSide(1:nt);
                catch
                    warning('No biased-side parameter found; including all trials by default')
                    idxRewardSide = true(1,nt);
                end
            end

        case 'all'
            idxRewardSide = true(1,nt);
        otherwise
            error('Choose a valid highRewardSide: "left", "right", or "all"');
    end
    
    switch preStimMovement
        case 'active'
            idxPreStim = wm.epochs(1).isMoving;
            idxPreStim = idxPreStim(1:nt);
        case 'quiescent'
            idxPreStim = ~wm.epochs(1).isMoving;
            idxPreStim = idxPreStim(1:nt);
        case 'all'
            idxPreStim = true(1,nt);
        otherwise
            error('Choose a valid preStimMovement: "active", "quiescent", or "all"');
    end
    
    switch responseType
        case 'correct'
%             idxCorrect = firstMoveDirs == b.events.correctResponseValues(1:nt);
            idxCorrect = b.events.feedbackValues == 1;
            idxCorrect = idxCorrect(1:nt);
        case 'incorrect'
%             idxCorrect = firstMoveDirs ~= b.events.correctResponseValues(1:nt);
            idxCorrect = b.events.feedbackValues == 0;
            idxCorrect = idxCorrect(1:nt);
        case 'all'
            idxCorrect = true(1,nt);
        otherwise
            error('Choose a valid responseType: "correct", "incorrect", or "all"');
    end

    switch rewardOutcome
        case 'rewarded'
            try 
                idxRewarded = ~isnan(et(5).daqTime(:))';
                idxRewarded = idxRewarded(1:nt);
            catch 
                warning('No biased-side parameter found; including all trials by default')
                idxRewarded = true(1,nt);
            end

        case 'unrewarded'
            try 
                idxRewarded = (isnan(et(4).daqTime(:)))'.*(b.events.feedbackValues == 1);
                idxRewarded = idxRewarded(1:nt);
            catch 
                warning('No biased-side parameter found; including all trials by default')
                idxRewarded = true(1,nt);
            end

        case 'all'
            idxRewarded = true(1,nt);
        otherwise
            error('Choose a valid rewardOutcome: "rewarded", "unrewarded", or "all"');
    end        

    %track the stimulus 1 trial in the past
    switch pastStimulus 
        case 'right' %stimulus was on the right
            idxPastStimulus = double(b.events.contrastValues > 0);
            idxPastStimulus = [nan(1,trialsBack), idxPastStimulus(1:nt-trialsBack)];
        case 'left' %stimulus was on the left
            idxPastStimulus = double(b.events.contrastValues < 0);
            idxPastStimulus = [nan(1,trialsBack), idxPastStimulus(1:nt-trialsBack)];
        case 'zero'
            idxPastStimulus = double(b.events.contrastValues == 0);
            idxPastStimulus = [nan(1,trialsBack), idxPastStimulus(1:nt-trialsBack)];
        case 'all'
            idxPastStimulus = true(1,nt);
        otherwise
            error('Choose a valid pastStimulus: "right", "left", "zero", or "all"');
    end

    %track the choice 1 trial in the past
    switch pastMovementDir 
        case 'ccw' %moved the stimulus left
            idxPastDirection = double(b.events.responseValues == 1);
            idxPastDirection = [nan(1,trialsBack), idxPastDirection(1:nt-trialsBack)];
        case 'cw' %moved the stimulus right
            idxPastDirection = double(b.events.responseValues == -1);
            idxPastDirection = [nan(1,trialsBack), idxPastDirection(1:nt-trialsBack)];
        case 'all'
            idxPastDirection = true(1,nt);
        otherwise
            error('Choose a valid pastMovementDir: "ccw", "cw", or "all"');
    end
    
    %track the timing 1 trial in the past
    switch pastMovementTime 
        case 'early'
            idxPastMovement = wm.epochs(2).isMoving;
            idxPastMovement = [nan(1,trialsBack), idxPastMovement(1:nt-trialsBack)];
        case 'late'
            idxPastMovement = ~wm.epochs(2).isMoving .* wm.epochs(3).isMoving;
            idxPastMovement = [nan(1,trialsBack), idxPastMovement(1:nt-trialsBack)];
        case 'all'
            idxPastMovement = true(1,nt);
        otherwise
            error('Choose a valid pastMovementTime: "early", "late", or "all"');
    end

    %track the outcome 1 trial in the past
    switch pastResponseType 
        case 'correct'
            idxPastCorrect = double(b.events.feedbackValues == 1);
            idxPastCorrect = [nan(1,trialsBack), idxPastCorrect(1:nt-trialsBack)];
        case 'incorrect'
            idxPastCorrect = double(b.events.feedbackValues == 0);
            idxPastCorrect = [nan(1,trialsBack), idxPastCorrect(1:nt-trialsBack)];
        case 'all'
            idxPastCorrect = true(1,nt);
        otherwise
            error('Choose a valid pastResponseType: "correct", "incorrect", or "all"');
    end 

    switch switchBlocks
        case 'beforeLeft'
            bl = [];
            switchLength = 50;
            switchPoints = find(diff(b.events.highRewardSideValues) < 0) + 1;
            for s = 1:length(switchPoints)
                bl = [bl (switchPoints(s) - switchLength):(switchPoints(s) - 1)];
            end
            idxSwitch = ismember(1:nt,bl);
        case 'afterLeft'
            al = [];
            switchLength = 50;
            switchPoints = find(diff(b.events.highRewardSideValues) < 0) + 1;
            for s = 1:length(switchPoints)
                al = [al switchPoints(s):(switchPoints(s) + switchLength)];
            end
            idxSwitch = ismember(1:nt,al);
        case 'beforeRight'
            br = [];
            switchLength = 50;
            switchPoints = find(diff(b.events.highRewardSideValues) > 0) + 1;
            for s = 1:length(switchPoints)
                br = [br (switchPoints(s) - switchLength):(switchPoints(s) - 1)];
            end
            idxSwitch = ismember(1:nt,br);
        case 'afterRight'
            ar = [];
            switchLength = 50;
            switchPoints = find(diff(b.events.highRewardSideValues) > 0) + 1;
            for s = 1:length(switchPoints)
                ar = [ar switchPoints(s):(switchPoints(s) + switchLength)];
            end
            idxSwitch = ismember(1:nt,ar);    
        case 'all'
            idxSwitch = true(1,nt);
        otherwise
            error('Choose a valid switchType: "before", "after", or "all"');
    end
    
    if isnumeric(trimBlocks)
        idxTrim = false(1,nt);
        switchPoints = find(diff(b.events.highRewardSideValues) ~= 0);
        for s = 1:length(switchPoints)
            if s == 1
                idxTrim(trimBlocks+1:switchPoints(s)) = true;
            elseif s < length(switchPoints)
            	idxTrim(switchPoints(s-1)+trimBlocks+1:switchPoints(s)) = true;
            elseif s == length(switchPoints)
                idxTrim(switchPoints(s-1)+trimBlocks+1:switchPoints(s)) = true;
                idxTrim(switchPoints(s)+trimBlocks+1:nt) = true;
            end
        end
    elseif ~isnumeric(trimBlocks) && strcmp(trimBlocks,'all')
        idxTrim = true(1,nt);
    else
        error('Choose a scalar value, or "all"');
    end

    if isnumeric(whichTrials)
        idxWhich = false(1,nt);
        idxWhich(whichTrials) = true;
    elseif ~isnumeric(whichTrials)
        idxWhich = true(1,nt);
    else
        error('Choose a valid indexing range, or "all"');
    end

    if isnumeric(specificRTs)
        RTs = wm.epochs(5).onsetTimes - et(1).daqTime;
        idxRT = RTs >= specificRTs(1) & RTs <= specificRTs(2);
        idxRT = idxRT(1:nt);
    elseif ~isnumeric(specificRTs)
        idxRT = true(1,nt);
    else
        error('Choose a valid indexing range, or "all"');
    end

%%%%%%%%%%%%%%%%%%%%%% END ASSIGNMENT SEGMENT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    
    % find the trials that pass all conditions (logical vector)
    cl = ...
        ismember(contrastVal, contrast).* ...
        idxRepeat .* ...
        idxDirection .* ...
        idxMovement .* ...
        idxRewardSide .* ...
        idxPreStim .* ...
        idxCorrect .* ...
        idxRewarded .* ...
        idxPastStimulus .* ...
        idxPastDirection .* ...
        idxPastMovement .* ...
        idxPastCorrect .* ...
        idxSwitch .* ...
        idxTrim .* ...
        idxWhich .* ...
        idxRT;
    
    % add logical vector to previous ones
    condLogical = [condLogical,cl];
end
    % index the whole logical vector (length = GRAND TOTAL of trials)
    [~, condIdx] = find(condLogical > 0);
end