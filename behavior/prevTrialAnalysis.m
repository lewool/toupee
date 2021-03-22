%set up structure called prevTrials that will be called when plotting
prevTrial = struct('prevTrialCond',[],'fraction',[]);
condList = {{'Early-hit'} {'Late-hit'} {'Early-miss'} {'Late-miss'} {'All'}};
for a =1:5
    prevTrial(a).prevTrialCond = condList{a};
end

%get trials of all different conditions, and convert them from logicals to list of all trials that meet the condictions  
for iX = 1:length(expInfo) 
    %find lists of trials with diff conditions happening in previous trial 
    [prevEarlyHit{iX}, ~] = selectCondition(expInfo(iX), [-1 -.5 -.12 -.05 0 .05 .12 .5 1],...
        behavioralData(iX), initTrialConditions('preStimMovement','quiescent',...
        'specificRTs',[.1 3],'pastMovementTime','early','pastResponseType','correct','movementTime','early'));
    [prevEarlyHit{iX}] = find(prevEarlyHit{iX}==1);
    
    [prevLateHit{iX}, ~] = selectCondition(expInfo(iX), [-1 -.5 -.12 -.05 0 .05 .12 .5 1],...
        behavioralData(iX), initTrialConditions('movementTime','late','preStimMovement','quiescent',...
        'specificRTs',([.1 3]),'pastMovementTime','late','pastResponseType','correct','movementTime','early'));
    [prevLateHit{iX}] = find(prevLateHit{iX}==1);

    [prevEarlyMiss{iX}, ~] = selectCondition(expInfo(iX), [-1 -.5 -.12 -.05 0 .05 .12 .5 1],...
        behavioralData(iX), initTrialConditions('preStimMovement','quiescent',...
        'specificRTs',([.1 3]),'pastMovementTime','early','pastResponseType','incorrect','movementTime','early'));
    [prevEarlyMiss{iX}] = find(prevEarlyMiss{iX}==1);

    [prevLateMiss{iX}, ~] = selectCondition(expInfo(iX), [-1 -.5 -.12 -.05 0 .05 .12 .5 1],...
        behavioralData(iX), initTrialConditions('movementTime','late','preStimMovement','quiescent',...
        'specificRTs',([.1 3]),'pastMovementTime','late','pastResponseType','incorrect','movementTime','early'));
    [prevLateMiss{iX}] = find(prevLateMiss{iX}==1);

    [earlyTrials{iX}, ~] = selectCondition(expInfo(iX), [-1 -.5 -.12 -.05 0 .05 .12 .5 1],...
        behavioralData(iX), initTrialConditions('movementTime','early','preStimMovement','quiescent',...
        'specificRTs',([.1 3])));
    [earlyTrials{iX}] = find(earlyTrials{iX}==1);

    [allTrials{iX}, ~] = selectCondition(expInfo(iX), [-1 -.5 -.12 -.05 0 .05 .12 .5 1],...
        behavioralData(iX), initTrialConditions('movementTime','all','preStimMovement','quiescent',...
        'specificRTs',([.1 3])));
    [allTrials{iX}] = find(allTrials{iX}==1);

    
    %calculate the fraction of early to late trials for all conditions 
    prevTrial(1).fraction{iX} = length(prevEarlyHit{iX}) ./ length(allTrials{iX});
    prevTrial(2).fraction{iX} = length(prevLateHit{iX}) ./ length(allTrials{iX});
    prevTrial(3).fraction{iX} = length(prevEarlyMiss{iX}) ./ length(allTrials{iX});
    prevTrial(4).fraction{iX} = length(prevLateMiss{iX}) ./ length(allTrials{iX});
    prevTrial(5).fraction{iX} = length(earlyTrials{iX}) ./ length(allTrials{iX});
  

end

%%
boxplot(prevTrial.prevTrialCond, prevTrial.fraction)





