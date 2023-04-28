function [predictors, windows] = getPredictors(expInfo, behavioralData, featureList, Fs)

% set up default struct
predictors = struct;

structNames = fieldnames(predictors);
contrasts = unique(expInfo.block.events.contrastValues);
nt = length(behavioralData.eventTimes(1).daqTime);
hemisphere = expInfo.hemisphere;

% where to start kernel window
stimStart = -0.5 * 1/Fs;
stimWindow = 2.5 * 1/Fs;
moveStart = -0.5 * 1/Fs;
moveWindow = 1.5 * 1/Fs;
blockStart = -.5 * 1/Fs;
blockWindow = 2.5 * 1/Fs;
rewardStart = -0.5 * 1/Fs;
rewardWindow = 1.5 * 1/Fs;
interactionStart = -0.5 * 1/Fs;
interactionWindow = 1.5 * 1/Fs;
valueStart = -.5 * 1/Fs;
valueWindow = 2.5 * 1/Fs;
RTStart = -1 * 1/Fs;
RTWindow = 2 *1/Fs;


% set up different trial conditions
[whichTrials, ~] = selectCondition(expInfo, contrasts, behavioralData, ...
        initTrialConditions('movementTime','all'));
whichTrials = logical(whichTrials);

[~, ppl, ppr, ~, ~] = getPsychometric(expInfo, behavioralData, find(whichTrials), contrasts);
pp = mean([ppl;ppr]);

if hemisphere < 0
    trueStimuli = expInfo.block.events.contrastValues(1:nt);
    trialCorrectChoice = expInfo.block.events.correctResponseValues(1:nt);
    trueBlocks = expInfo.block.events.highRewardSideValues(1:nt);
    trueChoices = expInfo.block.events.responseValues(1:nt);
    pp = pp;
else
    trueStimuli = -expInfo.block.events.contrastValues(1:nt);
    trialCorrectChoice = -expInfo.block.events.correctResponseValues(1:nt);
    trueBlocks = -expInfo.block.events.highRewardSideValues(1:nt);
    trueChoices = -expInfo.block.events.responseValues(1:nt);
    pp = fliplr(1 - pp);
end
  
RTs = behavioralData.wheelMoves.epochs(5).onsetTimes - behavioralData.eventTimes(1).daqTime;
velocities = abs(behavioralData.wheelMoves.epochs(5).peakVel);
velTrace = cat(2,behavioralData.wheelMoves.traces.vel{:});
timeTrace = cat(2,behavioralData.wheelMoves.traces.time{:});

trueFeedback = double(expInfo.block.events.feedbackValues(1:nt));
signedFeedback = (trueFeedback*2)-1;

% assign the 0% stimuli as either 'left' or 'right' depending on the
% preassigned correct choice (not the mouse's choice)
trueStimuli(trueStimuli == 0) = eps;
trueStimuli(abs(trueStimuli) < .05) = ...
    trueStimuli(abs(trueStimuli) < .05).* trialCorrectChoice(abs(trueStimuli) < .05);

% low rewards are possible on sign-mismatched block and stimulus
% high rewards are possible on sign-matched block and stimulus
% 2 = high, 1 = low
actualValue(trueBlocks.*sign(trueStimuli) == -1) = 1;
actualValue(trueBlocks.*sign(trueStimuli) == 1) = 2;
actualValue(trueFeedback == 0) = 0;
actualSignedValue = actualValue .* sign(trueStimuli);

trueValue = trueBlocks.*sign(trueStimuli);

perf = abs([1 1 1 1 nan 0 0 0 0] - pp);
trueSignedExpReward = nan(1,nt);
for c = 1:length(contrasts)
    trueSignedExpReward(trueStimuli == contrasts(c)) = perf(c) * trueValue(trueStimuli == contrasts(c)) .* sign(trueStimuli(trueStimuli == contrasts(c)));
    trueExpReward(trueStimuli == contrasts(c)) = perf(c) * trueValue(trueStimuli == contrasts(c));
end
expRewardLevels = unique(trueExpReward);
expRewardLevels(isnan(expRewardLevels)) = [];
signedExpRewardLevels = unique(trueSignedExpReward);
signedExpRewardLevels(isnan(signedExpRewardLevels)) = [];


% interactions take into account the direction of the choice and the
% availability of high or low reward (dictated by block)
% 2 = correct contra turn & high reward
% 1 = correct contra turn & low reward
% 0 = incorrect
% -1 = correct ipsi low
% -2 = correct ipsi high

trueInteraction((trueBlocks == -1)  & (sign(trueFeedback) == 1) & (trueChoices == -1)) = -2;
trueInteraction((trueBlocks == 1)  & (sign(trueFeedback) == 1) & (trueChoices == -1)) = -1;
trueInteraction((trueBlocks == -1)  & (sign(trueFeedback) == 0) & (trueChoices == 1)) = eps;
trueInteraction((trueBlocks == 1)  & (sign(trueFeedback) == 0) & (trueChoices == 1)) = eps;

trueInteraction((trueBlocks == 1)  & (sign(trueFeedback) == 1) & (trueChoices == 1)) = 2;
trueInteraction((trueBlocks == -1)  & (sign(trueFeedback) == 1) & (trueChoices == 1)) = 1;
trueInteraction((trueBlocks == 1)  & (sign(trueFeedback) == 0) & (trueChoices == -1)) = -eps;
trueInteraction((trueBlocks == -1)  & (sign(trueFeedback) == 0) & (trueChoices == -1)) = -eps;

%% fetch trials

RTs1 = find((RTs <= prctile(RTs, 25)));
RTs2 = find((RTs > prctile(RTs, 25)) & (RTs <= prctile(RTs, 50)));
RTs3 = find((RTs > prctile(RTs, 50)) & (RTs <= prctile(RTs, 75)));
RTs4 = find((RTs > prctile(RTs, 75)));

lowStimIpsiTrials = find((trueStimuli > -0.25) .* (trueStimuli < 0) .* whichTrials);
highStimIpsiTrials = find((trueStimuli < -0.25) .* whichTrials);
lowStimContraTrials = find((trueStimuli < 0.25) .* (trueStimuli > 0) .* whichTrials);
highStimContraTrials = find((trueStimuli > 0.25) .* whichTrials);
zeroStimTrials = find((trueStimuli > -0.005) .* (trueStimuli < 0.005) .* whichTrials);

actionTrials = find(whichTrials);
leftChoiceTrials = find((trueChoices < 0) .* whichTrials);
rightChoiceTrials = find((trueChoices > 0) .* whichTrials);
leftBlockTrials = find((trueBlocks < 0) .* whichTrials);
rightBlockTrials = find((trueBlocks > 0) .* whichTrials);
correctTrials = find((trueFeedback > 0) .* whichTrials);
incorrectTrials = find((trueFeedback == 0) .* whichTrials);
highRewardTrials = find((actualValue == 2) .* whichTrials);
lowRewardTrials = find((actualValue == 1) .* whichTrials);

highValueTrials = find((trueValue == 1) .* whichTrials);
lowValueTrials = find((trueValue == -1) .* whichTrials);

% ipsiHighValTrials = find(trueSignedExpReward <= signedExpRewardLevels(4));
% ipsiLowValTrials = find((trueSignedExpReward <= signedExpRewardLevels(8)) .* (trueSignedExpReward > signedExpRewardLevels(4)));
% contraLowValTrials = find((trueSignedExpReward <= signedExpRewardLevels(12)) .* (trueSignedExpReward > signedExpRewardLevels(8)));
% contraHighValTrials = find((trueSignedExpReward <= signedExpRewardLevels(16)) .* (trueSignedExpReward > signedExpRewardLevels(12)));

% lowestExpValTrials = find(trueExpReward <= expRewardLevels(4));
% lowExpValTrials =find((trueExpReward <= expRewardLevels(8)) .* (trueExpReward > expRewardLevels(4)));
% highExpValTrials =find((trueExpReward <= expRewardLevels(12)) .* (trueExpReward > expRewardLevels(8)));
% highestExpValTrials =find((trueExpReward <= expRewardLevels(16)) .* (trueExpReward > expRewardLevels(12)));

% highContraTrials = find(trueInteraction == 2);
% lowContraTrials = find(trueInteraction == 1);
% noneContraTrials = find(trueInteraction == eps);
% noneIpsiTrials = find(trueInteraction == -eps);
% lowIpsiTrials = find(trueInteraction == -1);
% highIpsiTrials =  find(trueInteraction == -2);


%%


prestimulusTimes = behavioralData.eventTimes(1).daqTime - 1;
% prestimulusTimes(isnan(prestimulusTimes)) = [];
stimulusTimes = behavioralData.eventTimes(1).daqTime;
% stimulusTimes(isnan(stimulusTimes)) = [];
gocueTimes = behavioralData.eventTimes(2).daqTime;
movementTimes = behavioralData.wheelMoves.epochs(5).onsetTimes;
% movementTimes(isnan(movementTimes)) = [];
outcomeTimes = behavioralData.eventTimes(5).daqTime;
% outcomeTimes(isnan(outcomeTimes)) = [];
blockTimes = behavioralData.wheelMoves.epochs(5).onsetTimes;

% replace defaults with values
for p = featureList   
   switch char(p)
       
       case 'stimulus_all'
%            nonzeroContrasts = contrasts(abs(contrasts) > 0);
           for c = 1:length(contrasts)
               [~, si] = selectCondition(expInfo, contrasts(c), behavioralData, ...
                   initTrialConditions('movementTime','late'));
               predictors.(matlab.lang.makeValidName(strcat('stimulus',num2str(c)))).times = stimulusTimes(si)';
               predictors.(matlab.lang.makeValidName(strcat('stimulus',num2str(c)))).values = ones(length(si),1);
           end
           
           windows.stimulus = linspace(stimStart, stimStart+stimWindow-1, stimWindow);

        case 'stimulus'
           %find high-contrast contra trials
           predictors.(matlab.lang.makeValidName('stimulusContraHigh')).times = stimulusTimes(highStimContraTrials)';
           predictors.(matlab.lang.makeValidName('stimulusContraHigh')).values = ones(length(highStimContraTrials),1);
           predictors.(matlab.lang.makeValidName('stimulusContraHigh')).color = [0 .0 1];
           
           %find low-contrast contra trials
           predictors.(matlab.lang.makeValidName('stimulusContraLow')).times = stimulusTimes(lowStimContraTrials)';
           predictors.(matlab.lang.makeValidName('stimulusContraLow')).values = ones(length(lowStimContraTrials),1);
           predictors.(matlab.lang.makeValidName('stimulusContraLow')).color = [0 .4 1];
           
           %find zero-contrast trials
           predictors.(matlab.lang.makeValidName('stimulusZero')).times = stimulusTimes(zeroStimTrials)';
           predictors.(matlab.lang.makeValidName('stimulusZero')).values = ones(length(zeroStimTrials),1);
           predictors.(matlab.lang.makeValidName('stimulusZero')).color = [.75 .75 .75];
           
           %find low-contrast ipsi trials
           predictors.(matlab.lang.makeValidName('stimulusIpsiLow')).times = stimulusTimes(lowStimIpsiTrials)';
           predictors.(matlab.lang.makeValidName('stimulusIpsiLow')).values = ones(length(lowStimIpsiTrials),1);
           predictors.(matlab.lang.makeValidName('stimulusIpsiLow')).color = [1 0 0];
           
           %find high-contrast ipsi trials
           predictors.(matlab.lang.makeValidName('stimulusIpsiHigh')).times = stimulusTimes(highStimIpsiTrials)';
           predictors.(matlab.lang.makeValidName('stimulusIpsiHigh')).values = ones(length(highStimIpsiTrials),1);
           predictors.(matlab.lang.makeValidName('stimulusIpsiHigh')).color = [.5 0 0];
            
           windows.stimulus = linspace(stimStart, stimStart+stimWindow-1, stimWindow);
       
       case 'RT'
           %find 1st fastest RT trials
           predictors.(matlab.lang.makeValidName('RT1')).times = stimulusTimes(RTs1)';
           predictors.(matlab.lang.makeValidName('RT1')).values = ones(length(RTs1),1);
           predictors.(matlab.lang.makeValidName('RT1')).color = [1 0 1];
           
           %find 2nd fastest RT trials
           predictors.(matlab.lang.makeValidName('RT2')).times = stimulusTimes(RTs2)';
           predictors.(matlab.lang.makeValidName('RT2')).values = ones(length(RTs2),1);
           predictors.(matlab.lang.makeValidName('RT2')).color = [0.86 .5 .86];
           
           %find 3rd fastest RT trials
           predictors.(matlab.lang.makeValidName('RT3')).times = stimulusTimes(RTs3)';
           predictors.(matlab.lang.makeValidName('RT3')).values = ones(length(RTs3),1);
           predictors.(matlab.lang.makeValidName('RT3')).color = [.5 .86 .5];
           
           %find 4th fastest RT trials
           predictors.(matlab.lang.makeValidName('RT4')).times = stimulusTimes(RTs4)';
           predictors.(matlab.lang.makeValidName('RT4')).values = ones(length(RTs4),1);
           predictors.(matlab.lang.makeValidName('RT4')).color = [0 1 0];
            
           windows.RT = linspace(RTStart, RTStart+RTWindow-1, RTWindow);
            
        case 'action'
           predictors.(matlab.lang.makeValidName('action')).times = movementTimes(actionTrials)';
           predictors.(matlab.lang.makeValidName('action')).values = ones(length(actionTrials),1);        
           predictors.(matlab.lang.makeValidName('action')).color = [0 0 0];
           
           windows.action = linspace(moveStart, moveStart+moveWindow-1, moveWindow);
           
       case 'velocity'
           
           predictors.(matlab.lang.makeValidName('velocity')).times = timeTrace';
           predictors.(matlab.lang.makeValidName('velocity')).values = abs(zscore(velTrace))';        
           predictors.(matlab.lang.makeValidName('velocity')).color = [0 1 0];
           
%            windows.velocity = linspace(moveStart, moveStart+moveWindow-1, moveWindow);
           windows.velocity = 0;
        
        case 'signedvelocity'
           
           predictors.(matlab.lang.makeValidName('signedvelocity')).times = timeTrace';
           predictors.(matlab.lang.makeValidName('signedvelocity')).values = zscore(velTrace)';        
           predictors.(matlab.lang.makeValidName('signedvelocity')).color = [0 1 0];
           
           windows.velocity = linspace(moveStart, moveStart+moveWindow-1, moveWindow);
       
       case 'choice'
           predictors.(matlab.lang.makeValidName('choice')).times = movementTimes(actionTrials)';
           predictors.(matlab.lang.makeValidName('choice')).values = trueChoices(actionTrials)';
           predictors.(matlab.lang.makeValidName('choice')).color = [.5 .5 .5];
           
           windows.choice = linspace(moveStart, moveStart+moveWindow-1, moveWindow);
            
       case 'outcome'
           
           %find correct choices
           predictors.(matlab.lang.makeValidName('outcomeCorrect')).times = outcomeTimes(correctTrials)';
           predictors.(matlab.lang.makeValidName('outcomeCorrect')).values = ones(length(correctTrials),1);
           predictors.(matlab.lang.makeValidName('outcomeCorrect')).color = [0 .5 0];
           
           %find incorrect choices
           predictors.(matlab.lang.makeValidName('outcomeIncorrect')).times = outcomeTimes(incorrectTrials)';
           predictors.(matlab.lang.makeValidName('outcomeIncorrect')).values = ones(length(incorrectTrials),1);
           predictors.(matlab.lang.makeValidName('outcomeIncorrect')).color = [.7 .1 0];
           
           windows.outcome = linspace(rewardStart, rewardStart+rewardWindow-1, rewardWindow);
       
%        case 'outcome'
%            predictors.(matlab.lang.makeValidName('outcome')).times = outcomeTimes(actionTrials)';
%            predictors.(matlab.lang.makeValidName('outcome')).values = signedFeedback(actionTrials)';
%            predictors.(matlab.lang.makeValidName('outcome')).color = [0 .5 0];
%            
%            windows.outcome = linspace(rewardStart, rewardStart+rewardWindow-1, rewardWindow);
       
       case 'reward'
           
           %find correct high choices
           predictors.(matlab.lang.makeValidName('rewardHigh')).times = outcomeTimes(highRewardTrials)';
           predictors.(matlab.lang.makeValidName('rewardHigh')).values = ones(length(highRewardTrials),1);
           predictors.(matlab.lang.makeValidName('rewardHigh')).color = [0 .5 0];
           
           %find correct low choices
           predictors.(matlab.lang.makeValidName('rewardLow')).times = outcomeTimes(lowRewardTrials)';
           predictors.(matlab.lang.makeValidName('rewardLow')).values = ones(length(lowRewardTrials),1);
           predictors.(matlab.lang.makeValidName('rewardLow')).color = [.5 1 .5];
           
           %find incorrect choices
           predictors.(matlab.lang.makeValidName('rewardNone')).times = outcomeTimes(incorrectTrials)';
           predictors.(matlab.lang.makeValidName('rewardNone')).values = ones(length(incorrectTrials),1);
           predictors.(matlab.lang.makeValidName('rewardNone')).color = [.7 .1 0];
            
           windows.outcome = linspace(rewardStart, rewardStart+rewardWindow-1, rewardWindow);
           
       case 'interaction'
           %find high  choices
           predictors.(matlab.lang.makeValidName('intHigh')).times = outcomeTimes(highRewardTrials)';
           predictors.(matlab.lang.makeValidName('intHigh')).values = trueChoices(highRewardTrials);
           predictors.(matlab.lang.makeValidName('intHigh')).color = [0 .5 0];
           
           %find low  choices
           predictors.(matlab.lang.makeValidName('intLow')).times = outcomeTimes(lowRewardTrials)';
           predictors.(matlab.lang.makeValidName('intLow')).values = trueChoices(lowRewardTrials);
           predictors.(matlab.lang.makeValidName('intLow')).color = [.5 1 .5];
           
           %find incorrect  choices
           predictors.(matlab.lang.makeValidName('intNone')).times = outcomeTimes(incorrectTrials)';
           predictors.(matlab.lang.makeValidName('intNone')).values = trueChoices(incorrectTrials);
           predictors.(matlab.lang.makeValidName('intNone')).color = [.7 .1 0];
            
           windows.outcome = linspace(rewardStart, rewardStart+rewardWindow-1, rewardWindow);
           
        case 'block'
           predictors.(matlab.lang.makeValidName('block')).times = stimulusTimes(whichTrials)';
           predictors.(matlab.lang.makeValidName('block')).values = trueBlocks(whichTrials)';
           predictors.(matlab.lang.makeValidName('block')).color = [1 .75 0];
           
           windows.block = linspace(blockStart, blockStart+blockWindow-1, blockWindow);
           
       case 'value'
           predictors.(matlab.lang.makeValidName('value')).times = stimulusTimes(whichTrials)';
           predictors.(matlab.lang.makeValidName('value')).values = trueValue(whichTrials)';
           predictors.(matlab.lang.makeValidName('value')).color = [1 .5 0];
           
           windows.value = linspace(valueStart, valueStart+valueWindow-1, valueWindow);
       
       case 'signedvalue'
           predictors.(matlab.lang.makeValidName('signedvalue')).times = gocueTimes(highValueTrials)';
           predictors.(matlab.lang.makeValidName('signedvalue')).values = sign(trueSignedValue(highValueTrials));
           predictors.(matlab.lang.makeValidName('signedvalue')).color = [.75 0 .5];
           
           windows.value = linspace(valueStart, valueStart+valueWindow-1, valueWindow);
           
       case 'expvalue'
           %find high-contrast contra trials
           predictors.(matlab.lang.makeValidName('highestValue')).times = gocueTimes(highestExpValTrials)';
           predictors.(matlab.lang.makeValidName('highestValue')).values = ones(length(highestExpValTrials),1);
           predictors.(matlab.lang.makeValidName('highestValue')).color = [0 .0 1];
           
           %find low-contrast left trials
           predictors.(matlab.lang.makeValidName('highValue')).times = gocueTimes(highExpValTrials)';
           predictors.(matlab.lang.makeValidName('highValue')).values = ones(length(highExpValTrials),1);
           predictors.(matlab.lang.makeValidName('highValue')).color = [0 .4 1];
           %find low-contrast right trials
           predictors.(matlab.lang.makeValidName('lowValue')).times = gocueTimes(lowExpValTrials)';
           predictors.(matlab.lang.makeValidName('lowValue')).values = ones(length(lowExpValTrials),1);
           predictors.(matlab.lang.makeValidName('lowValue')).color = [1 0 0];
           
           %find high-contrast right trials
           predictors.(matlab.lang.makeValidName('lowestValue')).times = gocueTimes(lowestExpValTrials)';
           predictors.(matlab.lang.makeValidName('lowestValue')).values = ones(length(lowestExpValTrials),1);
           predictors.(matlab.lang.makeValidName('lowestValue')).color = [.5 0 0];
           
           windows.value = linspace(valueStart, valueStart+valueWindow-1, valueWindow);
       
       case 'signedexpvalue'
           %find high-contrast contra trials
           predictors.(matlab.lang.makeValidName('valueContraHigh')).times = gocueTimes(contraHighValTrials)';
           predictors.(matlab.lang.makeValidName('valueContraHigh')).values = ones(length(contraHighValTrials),1);
           predictors.(matlab.lang.makeValidName('valueContraHigh')).color = [0 .0 1];
           
           %find low-contrast left trials
           predictors.(matlab.lang.makeValidName('valueContraLow')).times = gocueTimes(contraLowValTrials)';
           predictors.(matlab.lang.makeValidName('valueContraLow')).values = ones(length(contraLowValTrials),1);
           predictors.(matlab.lang.makeValidName('valueContraLow')).color = [0 .4 1];
           %find low-contrast right trials
           predictors.(matlab.lang.makeValidName('valueIpsiLow')).times = gocueTimes(ipsiLowValTrials)';
           predictors.(matlab.lang.makeValidName('valueIpsiLow')).values = ones(length(ipsiLowValTrials),1);
           predictors.(matlab.lang.makeValidName('valueIpsiLow')).color = [1 0 0];
           
           %find high-contrast right trials
           predictors.(matlab.lang.makeValidName('valueIpsiHigh')).times = gocueTimes(ipsiHighValTrials)';
           predictors.(matlab.lang.makeValidName('valueIpsiHigh')).values = ones(length(ipsiHighValTrials),1);
           predictors.(matlab.lang.makeValidName('valueIpsiHigh')).color = [.5 0 0];
           
           windows.value = linspace(valueStart, valueStart+valueWindow-1, valueWindow);
       
       otherwise
           error('"%s" is not a recognized feature name',char(p));
   end
end



%%%% workbench

% [~,lmi] = selectCondition(block, contrasts, eventTimes, initTrialConditions('movementTime',patience,'movementDir','cw'));
%             predictors.(matlab.lang.makeValidName('movementLeft')).times = ...
%                 eventTimes(strcmp({eventTimes.event},'prestimulusQuiescenceEndTimes')).daqTime(lmi)';
%             predictors.(matlab.lang.makeValidName('movementLeft')).times = 
%             [~,rmi] = selectCondition(block, contrasts, eventTimes, initTrialConditions('movementTime',patience,'movementDir','ccw'));
%             predictors.(matlab.lang.makeValidName('movementRight')).times = ...
%                 eventTimes(strcmp({eventTimes.event},'prestimulusQuiescenceEndTimes')).daqTime(rmi)';

%  %find high left trials
%             predictors.(matlab.lang.makeValidName('highLeftReward')).times = ...
%                 eventTimes(strcmp({eventTimes.event},'stimulusOnTimes')).daqTime(lhr)';
%             predictors.(matlab.lang.makeValidName('highLeftReward')).values = ones(length(lhr),1);
%             
%             %find high right trials
%             predictors.(matlab.lang.makeValidName('highRightReward')).times = ...
%                 eventTimes(strcmp({eventTimes.event},'stimulusOnTimes')).daqTime(rhr)';
%             predictors.(matlab.lang.makeValidName('highRightReward')).values = ones(length(rhr),1);