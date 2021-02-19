function [relativeTimes, sortIdWhisk] = sortTrialByWhisk(whichTrials,eyeData, et, wm)

%find time diff between stimOn and moveOn, sort trials by this
%difference
timeDiffs = wm.epochs(5).onsetTimes(whichTrials) - et(1).daqTime(whichTrials);
alignedFace = eyeData.eta.alignedFace;
timeRange = alignedFace{1}(:,95:101,2);
[~,sortIdx] = sort(timeDiffs,'ascend');
[~,sortIdWhisk] = sort(timeRange,'ascend');



%record stimOn, movOn, rewardOn times per trial
trialTimes = [...
    et(1).daqTime(whichTrials(sortIdx))',...
    wm.epochs(5).onsetTimes(whichTrials(sortIdx))',...
    et(5).daqTime(whichTrials(sortIdx))',...
];

%0-align all times to either stimOn, moveOn, or rewardOn
relativeTimes(:,:,1) = trialTimes - trialTimes(:,1);
relativeTimes(:,:,2) = trialTimes - trialTimes(:,2);
relativeTimes(:,:,3) = trialTimes - trialTimes(:,3);
