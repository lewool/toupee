function  [earlyTrialsWhisk,lateTrialsWhisk,earlyTrialsPupil,lateTrialsPupil] = earlyVsLate(expInfo, behavioralData, neuralData,eyeData)

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

%choose our conditions
%remove trials when mouse was moving wheel before stim onset with if
%statement

for iX = 1:length(expInfo) 
    
    [earlyTrials{iX}, ~] = selectCondition(expInfo(iX), [-1 -.5 -.12 -.05 0 .05 .12 .5 1],...
        behavioralData(iX), initTrialConditions('movementTime','early','pastMovementTime','early','preStimMovement','quiescent','specificRTs',([.1 3])));
    [lateTrials{iX}, ~] = selectCondition(expInfo(iX), [-1 -.5 -.12 -.05 0 .05 .12 .5 1],...
        behavioralData(iX), initTrialConditions('movementTime','late','pastMovementTime','late','preStimMovement','quiescent','specificRTs',([.1 3])));
    %remove trials when wheel moves pre-stim
    %earlyTrialsWithoutPrestimWheel{iX} = earlyLogical{iX} .* ~behavioralData(iX).wheelMoves.epochs(1).isMoving;
    %lateTrialsWithoutPrestimWheel{iX} = lateLogical{iX} .* ~behavioralData(iX).wheelMoves.epochs(1).isMoving;
      
    earlyTrials{iX} = find(earlyTrials{iX} == 1);
    lateTrials{iX} = find(lateTrials{iX} == 1);


    %place our different indexed ROIs in vectors 
    %index EARLY pupil
    earlyTrialsPupil(iX,:) = nanmean(eyeData(iX).eta.alignedFace{3}(earlyTrials{iX},:,1));

    %index LATE pupil size
    lateTrialsPupil(iX,:) = nanmean(eyeData(iX).eta.alignedFace{3}(lateTrials{iX},:,1));
 
   
    %index EARLY whisking (ROI 2)
    earlyTrialsWhisk(iX,:) = ...
        nanmean(eyeData(iX).eta.alignedFace{3}(earlyTrials{iX},:,2));
    
    %late whisking
    lateTrialsWhisk(iX,:) = ...
        nanmean(eyeData(iX).eta.alignedFace{3}(lateTrials{iX},:,2));
     
    
     %Paw movement 
    earlyTrialsPaw(iX,:) = ...
        nanmean(eyeData(iX).eta.alignedFace{3}(earlyTrials{iX},:,4));
    lateTrialsPaw(iX,:) = ...
        nanmean(eyeData(iX).eta.alignedFace{3}(lateTrials{iX},:,4));
    
    %index neural data EARLY trials
    cellMeansE{iX} = squeeze(...
        nanmean(neuralData(iX).eta.alignedResps{1}(earlyTrials{iX},:,:)));
    
    %index neural data LATE trials
    cellMeansL{iX} = squeeze(...
        nanmean(neuralData(iX).eta.alignedResps{1}(lateTrials{iX},:,:)));
     
end

%% plot the experiment-by-experiment vectors from the Trials-arrays
%plot early vs late trial pupil
figure;
for iX = 1: length(expInfo)
    plot(eyeData(iX).eta.eventWindow,earlyTrialsPupil,'Color',[235/255, 108/255, 0/255])
    plot(eyeData(iX).eta.eventWindow,mean(earlyTrialsPupil),'Color',[255/255, 128/255, 0/255],'Linewidth', 3)
    xlabel('Time(s) reward aligned','Fontsize',14)
    ylabel('Pupil size','Fontsize',14)
    %title("Early vs late trial pupil")
    hold on
    plot(eyeData(iX).eta.eventWindow,lateTrialsPupil,'Color',[0/255, 133/255, 133/255])
    plot(eyeData(iX).eta.eventWindow,mean(lateTrialsPupil),'Color',[0/255, 153/255, 153/255], 'Linewidth', 3)
    line([0,0],[-1,2],'Color','k','Linestyle','--');
    xlim([-0.1 1])
    box off
end

%%
%plot early vs late whisking 
figure; 
for iX = 1: length(expInfo)
    plot(eyeData(iX).eta.eventWindow,earlyTrialsWhisk,'Color',[235/255, 108/255, 0/255])
    plot(eyeData(iX).eta.eventWindow,mean(earlyTrialsWhisk),'Color',[255/255, 128/255, 0/255],'Linewidth', 3)
    xlabel('Time(s) reward aligned','Fontsize',14)
    ylabel('Whisking','Fontsize',14)
    %title("Early vs late trial whisking")
    hold on
    plot(eyeData(iX).eta.eventWindow,lateTrialsWhisk,'Color',[0/255, 133/255, 133/255])
    plot(eyeData(iX).eta.eventWindow,mean(lateTrialsWhisk),'Color',[0/255, 153/255, 153/255], 'Linewidth', 3)
    line([0,0],[-1,1.5],'Color','k','Linestyle','--')
    box off
end
%%
%plot early vs late whisking 
figure; 
for iX = 1: length(expInfo)
    plot(eyeData(iX).eta.eventWindow,earlyTrialsWhisk,'Color',[235/255, 108/255, 0/255])
    plot(eyeData(iX).eta.eventWindow,mean(earlyTrialsWhisk),'Color',[255/255, 128/255, 0/255],'Linewidth', 3)
    xlabel('Time(s) reward aligned','Fontsize',14)
    ylabel('Whisking','Fontsize',14)
    %title("Early vs late trial whisking")
    hold on
    plot(eyeData(iX).eta.eventWindow,lateTrialsWhisk,'Color',[0/255, 133/255, 133/255])
    plot(eyeData(iX).eta.eventWindow,mean(lateTrialsWhisk),'Color',[0/255, 153/255, 153/255], 'Linewidth', 3)
    xlim([-0.1 1])
    line([0,0],[-1,2],'Color','k','Linestyle','--')
    box off
end
%%
%plot paw movement
figure;
for iX = 1: length(expInfo)
    plot(eyeData(iX).eta.eventWindow,earlyTrialsPaw,'Color',[235/255, 108/255, 0/255])
    plot(eyeData(iX).eta.eventWindow,mean(earlyTrialsPaw),'Color',[255/255, 128/255, 0/255],'Linewidth', 3)
    xlabel('Time(s) reward aligned','Fontsize',14)
    ylabel('Paw movement','Fontsize',14)
    %title("Early vs late trial paw movement")
    hold on
    %plot late trial paw
    plot(eyeData(iX).eta.eventWindow,lateTrialsPaw,'Color',[0/255, 133/255, 133/255])
    plot(eyeData(iX).eta.eventWindow,mean(lateTrialsPaw),'Color',[0/255, 153/255, 153/255], 'Linewidth', 3)
    %xlim([-0.5 0])
    line([0,0],[-1,2.5],'Color','k','Linestyle','--')
    box off
end
%overlay wheel movements? below code needs fixing 
%plot(eyeData(iX).eta.eventWindow,behavioralData(iX).wheelMoves.traces.time)
%%
%plot paw movement
figure;
for iX = 1: length(expInfo)
    plot(eyeData(iX).eta.eventWindow,earlyTrialsPaw,'Color',[235/255, 108/255, 0/255])
    plot(eyeData(iX).eta.eventWindow,mean(earlyTrialsPaw),'Color',[255/255, 128/255, 0/255],'Linewidth', 3)
    xlabel('Time(s) reward aligned','Fontsize',14)
    ylabel('Paw movement','Fontsize',14)
    %title("Early vs late trial paw movement")
    hold on
    %plot late trial paw
    plot(eyeData(iX).eta.eventWindow,lateTrialsPaw,'Color',[0/255, 133/255, 133/255])
    plot(eyeData(iX).eta.eventWindow,mean(lateTrialsPaw),'Color',[0/255, 153/255, 153/255], 'Linewidth', 3)
    xlim([-0.5 0.5])
    line([0,0],[-1,2.5],'Color','k','Linestyle','--')
    box off
end
%%
%{
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
xlabel('Time(s) stimulus aligned')
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
xlabel('Time(s) stimulus aligned')
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
    xlabel('Time(s) stimulus aligned')
    ylabel('Neural activity')
    hold on;
end
legend('Late', 'Early','Location', 'Northwest')

%}
