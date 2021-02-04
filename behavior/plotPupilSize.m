%% initiation
%specify which session to process 
%load data and cretae the data structure
expInfo = initExpInfo({{'LEW031'}},{{'2020-02-14',1,[1]}});
[expInfo, neuralData, behavioralData] = processExperiment(expInfo);
eyeData = getEyeData(expInfo);
[eyeData] = alignFace(expInfo, eyeData, behavioralData);


%% plot pupil size aligned to stim 
%plot(eyeData.eta.eventWindow,eyeData.eta.alignedFace{1}(:,:,1))
plot(eyeData.eta.eventWindow,mean(eyeData.eta.alignedFace{1}(:,:,1)))
xlabel('time, stim=0')
ylabel('pupil size')

%%
%aligned to 1st move
[~, trialsEarly] = ...
	selectCondition(expInfo, ...
	getUniqueContrasts(expInfo), ...
	behavioralData, ...
	initTrialConditions('specificRTs',[0 5]));


%plot(eyeData.eta.eventWindow,eyeData.eta.alignedFace{2}(:,:,1))
plot(eyeData.eta.eventWindow,nanmean(eyeData.eta.alignedFace{2}(trialsEarly,:,1)))
xlabel('time, move=0')
ylabel('pupil size')

%%
%aligned to reward
%plot(eyeData.eta.eventWindow,eyeData.eta.alignedFace{3}(:,:,1))
plot(eyeData.eta.eventWindow,mean(eyeData.eta.alignedFace{3}(:,:,1)))
xlabel('time, reward=0')
ylabel('pupil size')

%%
%plot early trials at 1st move
[~, trialsEarly] = ...
	selectCondition(expInfo, ...
	getUniqueContrasts(expInfo), ...
	behavioralData, ...
	initTrialConditions('movementTime','early'));


%plot(eyeData.eta.eventWindow,eyeData.eta.alignedFace{2}(:,:,1))
plot(eyeData.eta.eventWindow, ...
    nanmean(eyeData.eta.alignedFace{2}(trialsEarly,:,1)))
xlabel('time, move=0')
ylabel('pupil size')
title('Early trials')
 
