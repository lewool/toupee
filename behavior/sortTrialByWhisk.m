function [relativeTimes, sortIdWhisk] = sortTrialByWhisk(whichTrials,eyeData, et, wm)

%find time diff between stimOn and moveOn, sort trials by this difference
timeDiffs = wm.epochs(5).onsetTimes(whichTrials{1})- et(1).daqTime(whichTrials{1});

%sort trials by time to plot 1st move % reward time 
[~,sortIdx] = sort(timeDiffs,'ascend');

%sort trials by amount of whisking pre-stim
alignedFace = eyeData.eta.alignedFace;
timeRange = alignedFace{1}(:,95:101,2);
[~,sortIdWhisk] = sort(timeRange,'ascend');



%record stimOn, moveOn, rewardOn times per trial
trialTimes = [...
    et(1).daqTime(whichTrials(sortIdx))',...
    wm.epochs(5).onsetTimes(whichTrials(sortIdx))',...
    et(5).daqTime(whichTrials(sortIdx))',...
];

%0-align all times to either stimOn, moveOn, or rewardOn
relativeTimes(:,:,1) = trialTimes - trialTimes(:,1);
relativeTimes(:,:,2) = trialTimes - trialTimes(:,2);
relativeTimes(:,:,3) = trialTimes - trialTimes(:,3);
