function  [earlySessionsWhisk,lateSessionsWhisk] = earlyVsLate(expInfo, behavioralData, neuralData,eyeData)

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
    earlySessionsPupil(iX,:) = nanmean(eyeData(iX).eta.alignedFace{1}(trialsEarly_iX,:,1));

    %index LATE pupil size
    [~,trialsLate_iX] = selectCondition(expInfo(iX), ...
        getUniqueContrasts(expInfo(iX)),behavioralData(iX), ...
        initTrialConditions('movementTime','late'));
    lateSessionsPupil(iX,:) = nanmean(eyeData(iX).eta.alignedFace{1}(trialsLate_iX,:,1));
 
   
    %index EARLY whisking (ROI 2)
    earlySessionsWhisk(iX,:) = ...
        nanmean(eyeData(iX).eta.alignedFace{1}(trialsEarly_iX,:,2));
    
    %late whisking
    lateSessionsWhisk(iX,:) = ...
        nanmean(eyeData(iX).eta.alignedFace{1}(trialsLate_iX,:,2));
     
    
     %Paw movement (all trials)
    allSessionsPaw(iX,:) = ...
        nanmean(eyeData(iX).eta.alignedFace{1}(:,:,4));
    
    %index neural data EARLY trials
    cellMeansE{iX} = squeeze(...
        nanmean(neuralData(iX).eta.alignedResps{1}(trialsEarly_iX,:,:)));
    
    %index neural data LATE trials
    cellMeansL{iX} = squeeze(...
        nanmean(neuralData(iX).eta.alignedResps{1}(trialsLate_iX,:,:)));
     
end

%% plot the experiment-by-experiment vectors from the sessions-arrays
%plot early trial pupil
figure;
plot(eyeData(iX).eta.eventWindow,earlySessionsPupil)
xlabel('Time(s) Stimulus aligned')
ylabel('pupil size')
title(mouseName+": Early trials")

%plot late trial pupil
figure;
plot(eyeData(iX).eta.eventWindow,lateSessionsPupil)
xlabel('Time(s) Stimulus aligned')
ylabel('pupil size')
title(mouseName+": Late trials")

%plot early trial whisking 
figure;
plot(eyeData(iX).eta.eventWindow,earlySessionsWhisk)
xlabel('Time(s) Stimulus aligned')
ylabel('Whisker Stimulus')
title(mouseName+": Early trials")
%xlim([-0.6, 0.1]) 

%plot late trial whisking
figure;
plot(eyeData(iX).eta.eventWindow,lateSessionsWhisk)
xlabel('Time(s) Stimulus aligned')
ylabel('Whisker Stimulus')
title(mouseName+": Late trials")
%xlim([-0.6, 0.1]) 

%%
%snake plot of neural activity
%mean of all trials per cell in all sessions for early and late trials
%1 horizontal line is 1 cell, colour coded for activity  
%SHOULD BE MADE NORMALIZED!!
%cells are sorted by the time of their max peak
%early trials 
animalMeanNeuralE = cat(2,cellMeansE{1},cellMeansE{1},cellMeansE{3},cellMeansE{4},cellMeansE{5});
animalMeanNeuralEt = transpose(animalMeanNeuralE);
[~,meanPeaksE] = max(animalMeanNeuralE);
[~,sIdxE] = sort(meanPeaksE);
sIdxE = transpose(sIdxE);
figure;
imagesc(neuralData(iX).eta.eventWindow,1:size(animalMeanNeuralEt,1),animalMeanNeuralEt(sIdxE,:))
xlabel('Time(s) Stimulus aligned')
ylabel('#cells')
title(mouseName+": Early trials")
line([0,0],[0,70000])
xlim([-0.5 0])

%late trials 
figure;
animalMeanNeuralL = cat(2,cellMeansL{1},cellMeansL{1},cellMeansL{3},cellMeansL{4},cellMeansL{5});
animalMeanNeuralLt = transpose(animalMeanNeuralL);
[~,meanPeaksL] = max(animalMeanNeuralL);
[~,sIdxL] = sort(meanPeaksL);
sIdxL = transpose(sIdxL);
imagesc(neuralData(1).eta.eventWindow,1:size(animalMeanNeuralLt,3),animalMeanNeuralLt(sIdxL,:))
xlabel('Time(s) Stimulus aligned')
ylabel('#cells')
title(mouseName+": Late trials")
line([0,0],[0,70000])
xlim([-0.5 0])

%%
%compare mean neural activity at baseline -0.5-0 sec pre-stim 
%create array for the overall population mean activity in each session
populationActivityE = zeros(41,length(cellMeansE));
populationActivityL = zeros(41,length(cellMeansE));
%errE = std(cellMeansE{:})/sqrt(length(cellMeansE{:}));
%(irow,:)
figure;
for iexp = 1:length(cellMeansE)
    for irow = 1:41
        populationActivityE(irow,iexp) = nanmean(cellMeansE{iexp}(irow,:));
        populationActivityL(irow,iexp) = nanmean(cellMeansL{iexp}(irow,:));
    end
    plot(neuralData(iX).eta.eventWindow,populationActivityE(:,iexp),...
         'Color', [0, 0.4470, 0.7410],'Linewidth',2)
    plot(neuralData(iX).eta.eventWindow,populationActivityL(:,iexp),...
        'Color', [105/255, 165/255, 131/255],'Linewidth',2)
    xlim([-0.5 0])
    xlabel('Time(s) Stimulus aligned')
    ylabel('Neural activity')
    hold on;
end
legend('Late', 'Early','Location', 'Northwest')
%%
%plot paw (& wheel?) movement
figure;
plot(eyeData(iX).eta.eventWindow,allSessionsPaw)
hold on
xlabel('Time(s) Stimulus aligned')
ylabel('Paw Movement ')
title(mouseName+": All trials")
%overlay wheel movements? below code needs fixing 
%plot(eyeData(iX).eta.eventWindow,behavioralData(iX).wheelMoves.traces.time)


