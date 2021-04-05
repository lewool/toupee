whichSessions = 1:9; %LEW031
clear pupilEarly pupilLate whiskEarly whiskLate pawEarly pawLate
for iX = whichSessions
    [~, earlyTrials] = selectCondition(expInfo(iX), getUniqueContrasts(expInfo(iX)), behavioralData(iX), ...
    initTrialConditions('preStimMovement','quiescent','movementTime','early','specificRTs',[.1 Inf]));
    [~, lateTrials] = selectCondition(expInfo(iX), getUniqueContrasts(expInfo(iX)), behavioralData(iX), ...
    initTrialConditions('preStimMovement','quiescent','movementTime','late','specificRTs',[.1 Inf]));
    pupilEarly(iX,:) = smooth(mean(eyeData(iX).eta.alignedFace{2}(earlyTrials,:,1),1));
    pupilLate(iX,:) = smooth(mean(eyeData(iX).eta.alignedFace{2}(lateTrials,:,1),1));
    whiskEarly(iX,:) = smooth(mean(eyeData(iX).eta.alignedFace{3}(earlyTrials,:,2),1));
    whiskLate(iX,:) = smooth(mean(eyeData(iX).eta.alignedFace{3}(lateTrials,:,2),1));
    pawEarly(iX,:) = smooth(mean(eyeData(iX).eta.alignedFace{2}(earlyTrials,:,4),1));
    pawLate(iX,:) = smooth(mean(eyeData(iX).eta.alignedFace{2}(lateTrials,:,4),1));
end
pupilEarly = pupilEarly(whichSessions,:);
pupilLate = pupilLate(whichSessions,:);
whiskEarly = whiskEarly(whichSessions,:);
whiskLate = whiskLate(whichSessions,:);
pawEarly = pawEarly(whichSessions,:);
pawLate = pawLate(whichSessions,:);

%% plot pupil
allpupil = [pupilEarly; pupilLate];
pupilArray{1,1} = 1:size(pupilEarly,1);
pupilArray{2,1} = [1:size(pupilEarly,1)] + size(pupilEarly,1);
figure;
[meanPSTH, semPSTH, ~] = computePSTHs(allpupil,pupilArray);
plotPSTHs(eyeData(1).eta.eventWindow, cell2mat(meanPSTH), cell2mat(semPSTH), [1 .5 0; 0 .6 .6],'-')
set(gca,'tickdir','out', 'Fontsize',14);
xlabel('Time from movement onset (s)')
ylabel('Pupil area (zscored)')
box off
y=ylim;
line([0 0],[y(1) y(2)],'LineStyle','--','Linewidth',1,'Color',[.7 .7 .7]);
%xlim([-1.5 .1])
%legend()


%% plot whisking
allwhisk = [whiskEarly; whiskLate];
whiskArray{1,1} = 1:size(whiskEarly,1);
whiskArray{2,1} = [1:size(whiskEarly,1)] + size(whiskEarly,1);
figure;
[meanPSTH, semPSTH, ~] = computePSTHs(allwhisk,whiskArray);
plotPSTHs(eyeData(1).eta.eventWindow, cell2mat(meanPSTH), cell2mat(semPSTH), [1 .5 0; 0 .6 .6],'-')
set(gca,'tickdir','out','Fontsize',14);
xlabel('Time from reward (s)')
ylabel('Whisking (zscored)')
box off
y=ylim;
line([0 0],[y(1) y(2)],'LineStyle','--','Linewidth',1,'Color',[.7 .7 .7]);
xlim([-.2 2])

%% plot paw
allpaw = [pawEarly; pawLate];
pawArray{1,1} = 1:size(pawEarly,1);
pawArray{2,1} = [1:size(pawEarly,1)] + size(pawEarly,1);
figure;
[meanPSTH, semPSTH, ~] = computePSTHs(allpaw,pawArray);
plotPSTHs(eyeData(1).eta.eventWindow, cell2mat(meanPSTH), cell2mat(semPSTH), [1 .5 0; 0 .6 .6],'-')
set(gca,'tickdir','out','Fontsize',14);
xlabel('Time from movement onset (s)')
ylabel('Paw movement (zscored)')
box off
y=ylim;
line([0 0],[y(1) y(2)],'LineStyle','--','Linewidth',1,'Color',[.7 .7 .7]);
%xlim([-.2 2])



