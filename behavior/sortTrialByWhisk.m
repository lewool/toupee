function [sortIdx] = sortTrialByWhisk(whichTrials,eyeData)


%sort trials by time to plot 1st move % reward time 

%sort trials by amount of whisking pre-stim
alignedFace = eyeData.eta.alignedFace;
meanPreStimWhisk = mean(alignedFace{1}(:,95:101,2),2);
[~,sortIdx] = sort(meanPreStimWhisk,'descend');


%{
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
%}