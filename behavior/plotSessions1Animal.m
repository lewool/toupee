%% initiation
%specify which session to process 
%load data and cretae the data structure
mouseName={{'LEW031'}};
expDate = {{'2020-02-17',1,[1]}};
expInfo = initExpInfo(mouseName,expDate);
[expInfo, neuralData, behavioralData] = processExperiment(expInfo);
eyeData = getEyeData(expInfo);
[eyeData] = alignFace(expInfo, eyeData, behavioralData);

%%
%plot early trials at 1st move
[~, trialsEarly] = ...
	selectCondition(expInfo, ...
	getUniqueContrasts(expInfo), ...
	behavioralData, ...
	initTrialConditions('movementTime','early'));

%log the pupil data vector into your data array
   allSessionsArray(3,:)...
   = nanmean(eyeData.eta.alignedFace{2}(trialsEarly,:,1));

%%
%plot 
plot(eyeData.eta.eventWindow,allSessionsArray)
xlabel('time, move=0')
ylabel('pupil size')
%title(str(mouseName)+str(expDate)+'Early trials')
