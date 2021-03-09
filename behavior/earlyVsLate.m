function  [earlyTrialsWhisk,lateTrialsWhisk] = earlyVsLate(expInfo, behavioralData, neuralData,eyeData)

%% loop through the experiments for chosen analysis
%let's say you want to collect a mean Time-aligned trace from each
%experiment and put it into a persistent variable that logs every loop
%initialize the data array in advance
%nb: lengthOfTimeWindow will be different whether you are grabbing eye data
%(length = 201) or neural data (length = 41)

%make arrays of correct size
earlyTrialsPupil = zeros(length(expInfo),201);
lateTrialsPupil = zeros(length(expInfo),201);
earlyTrialsWhisk = zeros(length(expInfo),201);
lateTrialsWhisk = zeros(length(expInfo),201);
allTrialsPaw = zeros(length(expInfo),201);

%index EARLY pupil size and add to vector 
for iX = 1:length(expInfo)
    [~,trialsEarly_iX] = selectCondition(expInfo(iX), ...
	getUniqueContrasts(expInfo(iX)),behavioralData(iX), ...
	initTrialConditions('movementTime','early'));
    earlyTrialsPupil(iX,:) = nanmean(eyeData(iX).eta.alignedFace{1}(trialsEarly_iX,:,1));

    %index LATE pupil size
    [~,trialsLate_iX] = selectCondition(expInfo(iX), ...
        getUniqueContrasts(expInfo(iX)),behavioralData(iX), ...
        initTrialConditions('movementTime','late'));
    lateTrialsPupil(iX,:) = nanmean(eyeData(iX).eta.alignedFace{1}(trialsLate_iX,:,1));
 
   
    %index EARLY whisking (ROI 2)
    earlyTrialsWhisk(iX,:) = ...
        nanmean(eyeData(iX).eta.alignedFace{1}(trialsEarly_iX,:,2));
    
    %late whisking
    lateTrialsWhisk(iX,:) = ...
        nanmean(eyeData(iX).eta.alignedFace{1}(trialsLate_iX,:,2));
     
    
     %Paw movement (all trials)
    allTrialsPaw(iX,:) = ...
        nanmean(eyeData(iX).eta.alignedFace{1}(:,:,4));
    
    %index neural data EARLY trials
    cellMeansE{iX} = squeeze(...
        nanmean(neuralData(iX).eta.alignedResps{1}(trialsEarly_iX,:,:)));
    
    %index neural data LATE trials
    cellMeansL{iX} = squeeze(...
        nanmean(neuralData(iX).eta.alignedResps{1}(trialsLate_iX,:,:)));
     
end

%% plot the experiment-by-experiment vectors from the Trials-arrays
%plot early vs late trial pupil
figure;
for iX = 1: length(expInfo)
    plot(eyeData(iX).eta.eventWindow,earlyTrialsPupil,'Color',[235/255, 108/255, 0/255])
    plot(eyeData(iX).eta.eventWindow,mean(earlyTrialsPupil),'Color',[255/255, 128/255, 0/255],'Linewidth', 3)
    xlabel('Time(s) Stimulus aligned')
    ylabel('Pupil size')
    title("Early vs late trial pupil")
    hold on
    %plot late trial whisking
    plot(eyeData(iX).eta.eventWindow,lateTrialsPupil,'Color',[0/255, 133/255, 133/255])
    plot(eyeData(iX).eta.eventWindow,mean(lateTrialsPupil),'Color',[0/255, 153/255, 153/255], 'Linewidth', 3)
end

%%
%plot early vs late whisking 
figure;
for iX = 1: length(expInfo)
    plot(eyeData(iX).eta.eventWindow,earlyTrialsWhisk,'Color',[235/255, 108/255, 0/255])
    plot(eyeData(iX).eta.eventWindow,mean(earlyTrialsWhisk),'Color',[255/255, 128/255, 0/255],'Linewidth', 3)
    xlabel('Time(s) Stimulus aligned')
    ylabel('Whisking')
    title("Early vs late trial whisking")
    hold on
    %plot late trial whisking
    plot(eyeData(iX).eta.eventWindow,lateTrialsWhisk,'Color',[0/255, 133/255, 133/255])
    plot(eyeData(iX).eta.eventWindow,mean(lateTrialsWhisk),'Color',[0/255, 153/255, 153/255], 'Linewidth', 3)
end
%%
%snake plot of neural activity
%mean of all trials per cell in all Trials for early and late trials
%1 horizontal line is 1 cell, colour coded for activity  
%SHOULD BE MADE NORMALIZED!!
%cells are sorted by the time of their max peak
%early trials 
animalMeanNeuralE = cat(2,cellMeansE{1},cellMeansE{1},cellMeansE{3},cellMeansE{4},cellMeansE{5});
animalMeanNeuralEt = transpose(animalMeanNeuralE);
animalMeanNeuralEt = zscore(animalMeanNeuralEt);
[~,meanPeaksE] = max(animalMeanNeuralE);
[~,sIdxE] = sort(meanPeaksE);
sIdxE = transpose(sIdxE);
figure;
imagesc(neuralData(iX).eta.eventWindow,1:size(animalMeanNeuralEt,1),animalMeanNeuralEt(sIdxE,:))
xlabel('Time(s) Stimulus aligned')
ylabel('#cells')
title("Early trials")
line([0,0],[0,70000])
%xlim([-0.5 0])

%late trials 
figure;
animalMeanNeuralL = cat(2,cellMeansL{1},cellMeansL{1},cellMeansL{3},cellMeansL{4},cellMeansL{5});
animalMeanNeuralLt = transpose(animalMeanNeuralL);
animalMeanNeuralLt = zscore(animalMeanNeuralLt);
[~,meanPeaksL] = max(animalMeanNeuralL);
[~,sIdxL] = sort(meanPeaksL);
sIdxL = transpose(sIdxL);
imagesc(neuralData(1).eta.eventWindow,1:size(animalMeanNeuralLt,3),animalMeanNeuralLt(sIdxL,:))
xlabel('Time(s) Stimulus aligned')
ylabel('#cells')
title("Late trials")
line([0,0],[0,70000])
%xlim([-0.5 0])

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
plot(eyeData(iX).eta.eventWindow,allTrialsPaw)
hold on
xlabel('Time(s) Stimulus aligned')
ylabel('Paw Movement ')
title("All trials")
%overlay wheel movements? below code needs fixing 
%plot(eyeData(iX).eta.eventWindow,behavioralData(iX).wheelMoves.traces.time)


