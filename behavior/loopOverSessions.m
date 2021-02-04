%% make a list of mice/experiments you want to analyze

%{
mouseList = {...
    {'LEW031'}...
    {'LEW032'}};
%}

mouseName = {{'LEW031'}};
expList = { ...
    {'2020-02-03',1,[1]},...
    {'2020-02-14',1,[1]}...
    {'2020-02-17',1,[1]}...
    {'2020-02-18',1,[1]}...
    {'2020-02-25',1,[1]}...    
    };
    
%% load all the experiments into expInfo
expInfo = initExpInfo(mouseName,expList);

%% process the usual data
% the script knows to loop over all the experiments you listed above
% this will take a while but the command line will print progress
[expInfo, neuralData, behavioralData] = processExperiment(expInfo);
eyeData = getEyeData(expInfo);
[eyeData] = alignFace(expInfo, eyeData, behavioralData);

%% loop through the experiments for chosen analysis
%let's say you want to collect a mean Time-aligned trace from each
%experiment and put it into a persistent variable that logs every loop
%initialize the data array in advance
%nb: lengthOfTimeWindow will be different whether you are grabbing eye data
%(length = 201) or neural data (length = 41)

%make arrays of correct size
earlySessionsPupil = zeros(length(expInfo),201);
lateSessionsPupil = zeros(length(expInfo),201);
earlySessionsWhisk = zeros(length(expInfo),201);
lateSessionsWhisk = zeros(length(expInfo),201);
allSessionsPaw = zeros(length(expInfo),201);

%index EARLY pupil size and add to vector 
for iX = 1:length(expInfo)
    [~,trialsEarly_iX] = selectCondition(expInfo(iX), ...
	getUniqueContrasts(expInfo(iX)),behavioralData(iX), ...
	initTrialConditions('movementTime','early'));
    earlySessionsPupil(iX,:) = nanmean(eyeData(iX).eta.alignedFace{2}(trialsEarly_iX,:,1));

    %index LATE pupil size
    [~,trialsLate_iX] = selectCondition(expInfo(iX), ...
        getUniqueContrasts(expInfo(iX)),behavioralData(iX), ...
        initTrialConditions('movementTime','late'));
    lateSessionsPupil(iX,:) = nanmean(eyeData(iX).eta.alignedFace{2}(trialsLate_iX,:,1));
   
    %index EARLY whisking (ROI 2)
    earlySessionsWhisk(iX,:) = ...
        nanmean(eyeData(iX).eta.alignedFace{2}(trialsEarly_iX,:,2));
    
    %late whisking
    lateSessionsWhisk(iX,:) = ...
        nanmean(eyeData(iX).eta.alignedFace{2}(trialsLate_iX,:,2));
    
     %Paw movement (all sessions)
    allSessionsPaw(iX,:) = ...
        nanmean(eyeData(iX).eta.alignedFace{2}(:,:,4));
    
    %index neural data EARLY trials
    cellMeansE{iX} = squeeze(...
        nanmean(neuralData(iX).eta.alignedResps{2}(trialsEarly_iX,:,:)));
    
    %index neural data LATE trials
    cellMeansL{iX} = squeeze(...
        nanmean(neuralData(iX).eta.alignedResps{2}(trialsLate_iX,:,:)));
     
end

%% plot the experiment-by-experiment vectors in your persistent variable
%you can also do computations on this array (e.g., compute the SEM)
%plot early trial pupil
figure;
plot(eyeData(iX).eta.eventWindow,earlySessionsPupil)
xlabel('Time(s) Movement aligned')
ylabel('pupil size')
title(mouseName+": Early trials")

%plot late trial pupil
figure;
plot(eyeData(iX).eta.eventWindow,lateSessionsPupil)
xlabel('Time(s) Movement aligned')
ylabel('pupil size')
title(mouseName+": Late trials")

%plot early trial whisking 
figure;
plot(eyeData(iX).eta.eventWindow,earlySessionsWhisk)
xlabel('Time(s) Movement aligned')
ylabel('Whisker Stimulus')
title(mouseName+": Early trials")
%xlim([-0.6, 0.1]) 

%plot late trial whisking
figure;
plot(eyeData(iX).eta.eventWindow,lateSessionsWhisk)
xlabel('Time(s) Movement aligned')
ylabel('Whisker Stimulus')
title(mouseName+": Late trials")
%xlim([-0.6, 0.1]) 

%%
%overview of mean of all cells in all sessions for early and late trials 
%imagesc(neuralData(iX).eta.eventWindow,1:size(neuralData.eta.alignedResps{2},3),cellMeans{iX}')
figure;
animalMeanNeuralE = cat(2,cellMeansE{1},cellMeansE{2},cellMeansE{3},cellMeansE{4},cellMeansE{5});
[~,meanPeaksE] = max(animalMeanNeuralE);
[~,sIdxE] = sort(meanPeaksE);
imagesc(neuralData(1).eta.eventWindow,1:size(animalMeanNeuralE,3),animalMeanNeuralE(:,sIdxE))
xlabel('Time(s) Movement aligned')
title(mouseName+": Early trials")

figure;
animalMeanNeuralL = cat(2,cellMeansL{1},cellMeansL{2},cellMeansL{3},cellMeansL{4},cellMeansL{5});
[~,meanPeaksL] = max(animalMeanNeuralL);
[~,sIdxL] = sort(meanPeaksL);
imagesc(neuralData(1).eta.eventWindow,1:size(animalMeanNeuralL,3),animalMeanNeuralL(:,sIdxL))
xlabel('Time(s) Movement aligned')
title(mouseName+": Late trials")

%%
%simple line plot early trial neural
figure;
plot(neuralData(iX).eta.eventWindow,cellMeansE{iX})
xlabel('Time(s) Movement aligned')
ylabel('Neural activity')
title(mouseName+": Early trials")
xlim([-0.3 0.2])

%simple line plot late trial neural
figure;
plot(neuralData(iX).eta.eventWindow,cellMeansL{iX})
xlabel('Time(s) Movement aligned')
ylabel('Neural activity')
title(mouseName+": Late trials")
xlim([-0.3 0.2])

%%
%plot paw (& wheel?) movement
figure;
plot(eyeData(iX).eta.eventWindow,allSessionsPaw)
hold on
xlabel('Time(s) Movement aligned')
ylabel('Paw Movement ')
title(mouseName+": All trials")
%overlay wheel movements? below code needs fixing 
%plot(eyeData(iX).eta.eventWindow,behavioralData(iX).wheelMoves.traces.time)










