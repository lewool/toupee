
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
%%
%plot whisking off all trials aligned to stim for early and late trials 
 [~, trialsE] = ...
	selectCondition(expInfo, ...
	getUniqueContrasts(expInfo), ...
	behavioralData, ...
	initTrialConditions('movementTime','early'));

[~, trialsL] = ...
	selectCondition(expInfo, ...
	getUniqueContrasts(expInfo), ...
	behavioralData, ...
	initTrialConditions('movementTime','late'));

%plot(eyeData.eta.eventWindow,eyeData.eta.alignedFace{2}(:,:,1))
figure;
for t =1:length(trialsL)
    plot(eyeData.eta.eventWindow, ...
        smooth(eyeData.eta.alignedFace{1}(trialsE(t),:,2)),'Color',[255/255, 128/255, 0/255])
    hold off
    alpha(0.5)
    hold on
    plot(eyeData.eta.eventWindow, ...
        mean(eyeData.eta.alignedFace{1}(trialsE,:,2)),'Color',[230/255, 108/255, 0/255],'Linewidth',3)
    xlabel('time, stim=0')
    ylabel('whisking')
    ylim([-2 2])
    hold on 
    plot(eyeData.eta.eventWindow, ...
        smooth(eyeData.eta.alignedFace{1}(t,:,2)),'Color',[0/255, 153/255, 153/255])
    alpha(0.5)
    plot(eyeData.eta.eventWindow, ...
        mean(eyeData.eta.alignedFace{1}(trialsL,:,2)),'Color',[0/255, 133/255, 133/255],'Linewidth',3)
    xlabel('time, stim=0')
    ylabel('whisking')
    hold on
end   

   



