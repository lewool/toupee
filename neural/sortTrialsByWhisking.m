function [relativeTimes, sortIdxWhisk] = sortTrialsByWhisking(whichTrials, eyeData, et, wm)

%find time diff between stimOn and moveOn, sort trials by this
%difference
timeDiffs = wm.epochs(5).onsetTimes(whichTrials) - et(1).daqTime(whichTrials);
[~,sortIdx] = sort(timeDiffs,'ascend');

meanPreStimWhisk = mean(eyeData.eta.alignedFace{1}(whichTrials,95:101,2),2);
[~,sortIdxWhisk] = sort(meanPreStimWhisk,'descend'); 

%record stimOn, movOn, rewardOn times per trial
trialTimes = [...
    et(1).daqTime(whichTrials(sortIdxWhisk))',...
    wm.epochs(5).onsetTimes(whichTrials(sortIdxWhisk))',...
    et(5).daqTime(whichTrials(sortIdxWhisk))',...
];

%0-align all times to either stimOn, moveOn, or rewardOn
relativeTimes(:,:,1) = trialTimes - trialTimes(:,1);
relativeTimes(:,:,2) = trialTimes - trialTimes(:,2);
relativeTimes(:,:,3) = trialTimes - trialTimes(:,3);
