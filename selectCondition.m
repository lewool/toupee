function [condLogical, condIdx] = selectCondition(block, contrast, eventTimes, repeatType, movementDir, movementTime, highRewardSide, responseType, rewardOutcome, pastStimulus, pastMovementDir, pastResponseType)
% A wee function for picking the trial types you want to look at when analyzing your cells
% INPUTS: block struct, a contrast/side (this can be a vector), and some flags for
% the types of trials you want to include. If you don't care about a flag,
% choose 'all'. Both 'block' and eventTimes' can be nx1 structs containing
% exp data from multiple days for concatenation
% OUTPUT: an index of all the relevant trials from the experiment(s)

% 8 March 2018 Written by LEW
% 23 July 2018 Updated to accept multiple blocks/eventTimes
% 19 November 2018 Added optional trial-history selections
% 13 May 2019 Added rewarded/unrewarded outcome selection (likelihood experiments)

trialsBack = 1;

if nargin < 10
    warning('1-trial-back history was not analyzed');
    pastStimulus = 'all';
    pastMovementDir = 'all';
    pastResponseType = 'all';
end

condLogical = [];
for iBlock = 1:length(block)
    if length(block) == 1
        b = block;
        et = eventTimes;
    else
        b = block{iBlock};
        et = eventTimes{iBlock};
    end
    
    nt = numel(b.events.endTrialTimes);
    contrastVal = b.events.contrastValues(1:nt);

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
        case 'ccw' %moved the stimulus left
            idxDirection = b.events.responseValues == 1;
            idxDirection = idxDirection(1:nt);
        case 'cw' %moved the stimulus right
            idxDirection = b.events.responseValues == -1;
            idxDirection = idxDirection(1:nt);
        case 'all'
            idxDirection = true(1,nt);
        otherwise
            error('Choose a valid movementType: "ccw", "cw", or "all"');
    end

    switch movementTime
        case 'early'
            idxMovement = et(7).daqTime(1:nt) - et(2).daqTime(1:nt) <= 0;
            idxMovement = idxMovement(1:nt);
        case 'late'
            idxMovement = et(7).daqTime(1:nt) - et(2).daqTime(1:nt) > 0;
            idxMovement = idxMovement(1:nt);
        case 'all'
            idxMovement = true(1,nt);
        otherwise
            error('Choose a valid movementTime: "early", "late", or "all"');
    end

    switch highRewardSide
        case 'left' %left-stimulus trials had larger reward
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

        case 'right' %right-stimulus trials had larger reward
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

    switch responseType
        case 'correct'
            idxCorrect = b.events.feedbackValues == 1;
            idxCorrect = idxCorrect(1:nt);
        case 'incorrect'
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
                likelihoodTest = block.events.likelyRewardSideValues;
                idxRewarded = ~isnan(eventTimes(5).daqTime(:))';
                idxRewarded = idxRewarded(1:nt);
            catch 
                warning('No biased-side parameter found; including all trials by default')
                idxRewarded = true(1,nt);
            end
            
        case 'unrewarded'
            try 
                likelihoodTest = block.events.likelyRewardSideValues;
                idxRewarded = (isnan(eventTimes(4).daqTime(:)))'.*(block.events.feedbackValues == 1);
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
            error('Choose a valid movementType: "ccw", "cw", or "all"');
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
            error('Choose a valid movementType: "ccw", "cw", or "all"');
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
            error('Choose a valid responseType: "correct", "incorrect", or "all"');
    end     

    cl = ismember(contrastVal, contrast).* idxRepeat .* idxDirection .* idxMovement .* idxRewardSide .* idxCorrect .*idxRewarded .* idxPastStimulus .* idxPastDirection .* idxPastCorrect;
    condLogical = [condLogical,cl];
    
end
[~, condIdx] = find(condLogical > 0);
end