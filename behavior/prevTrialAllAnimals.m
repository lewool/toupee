%% Import the data
[~, ~, raw] = xlsread('C:\Users\Ella Svahn\Documents\eyedata\PrevTrialEffect_allAnimals.xlsx','Sheet1','A2:G123');
raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};
% Replace non-numeric cells with NaN
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % Find non-numeric cells
raw(R) = {NaN}; % Replace non-numeric cells

% Create output variable
data = reshape([raw{:}],size(raw));
% Create table
PT = table;

% Allocate imported array to column variable names
PT.animalID = data(:,1);
PT.propEarly_followingLate = data(:,2);
PT.propEarly_followingEarly = data(:,3);
PT.propEarly_followingCorrect = data(:,4);
PT.propEarly_followingIncorrect = data(:,5);
PT.propEarly_followingActive = data(:,6);
PT.propEarly_followingQuiescent = data(:,7);
% Clear temporary variables
clearvars data raw R;
%% plot previous correct vs incorrect 
colors = [ 0 0 1; 0 1 0; 0 .7 .7; 0 0 0];
figure;
for iX = 1:61
    pickColor = PT.animalID(iX);
    p = plot([1, 2],[PT.propEarly_followingCorrect(iX),PT.propEarly_followingIncorrect(iX)],...
        '-', 'LineWidth', 2, 'Color',colors(pickColor,:),'MarkerSize', 10);
    p.Color(4) = .5;
    hold on
end
xlim([.75 2.25])
ylim([0 1])
box off
set(gca,'tickdir','out', 'Fontsize', 12)
xticks([1 2])
set(gca, 'XTickLabels', {'Rewarded', 'Non-rewarded'})
xlabel('Previous trial condition')
ylabel('Fraction of early moves')
line([1 2],[.92 .92],'Color','k')
text(1.5,.94,'*','Fontsize',16)

%% plot previous earlry vs late 
colors = [ 0 0 1; 0 1 0; 0 .7 .7; 0 0 0];
figure;
for iX = 1:61
    pickColor = PT.animalID(iX);
    p = plot([1, 2],[PT.propEarly_followingEarly(iX),PT.propEarly_followingLate(iX)],...
        '-', 'LineWidth', 2, 'Color',colors(pickColor,:),'MarkerSize', 10);
    p.Color(4) = .5;
    hold on
end
hold on
xlim([.75 2.25])
ylim([0 1])
box off
set(gca,'tickdir','out','Fontsize', 12)
xticks([1 2])
set(gca, 'XTickLabels', {'Early', 'Late'})
xlabel('Previous trial condition')
ylabel('Fraction of early moves')
line([1 2],[.92 .92],'Color','k')
text(1.5,.94,'*','Fontsize',16)


%% plot pre stim period in current trial 
colors = [ 0 0 1; 0 1 0; 0 .7 .7; 0 0 0];
figure;
for iX = 1:61
    pickColor = PT.animalID(iX);
    p = plot([1, 2],[PT.propEarly_followingActive(iX),PT.propEarly_followingQuiescent(iX)],...
        '-', 'LineWidth', 2, 'Color',colors(pickColor,:),'MarkerSize', 10);
    p.Color(4) = .5;
    hold on
end
hold on
xlim([.75 2.25])
ylim([0 1])
box off
set(gca,'tickdir','out','Fontsize', 12)
xticks([1 2])
set(gca, 'XTickLabels', {'Active pre-stim', 'Quiescent pre-stim'})
ylabel('Fraction of impulsive moves')
xlabel('Pre-stimulus state')
line([1 2],[.92 .92],'Color','k')
text(1.5,.94,'*','Fontsize',16)

%% do a Two-sample Kolmogorov-Smirnov test
disp('movetime')
[h,p] = kstest2(PT.propEarly_followingEarly,PT.propEarly_followingLate)
disp('reward')
[h,p] = kstest2(PT.propEarly_followingCorrect,PT.propEarly_followingIncorrect)
disp('pre-stimmove')
[h,p] = kstest2(PT.propEarly_followingActive,PT.propEarly_followingQuiescent)


%% beeswarm plot (not yet working!)
%{
figure;
beeswarm(PT-animalID,[PT.propEarly_followingEarly,PT.propEarly_followingLate]',...
   'sort_style','hex','dot_size',1,'overlay_style','sd','colormap','winter');
hold on
%xlim([0 5])
%ylim([0 30])
ylabel('Fraction of early trials of all trials','Fontsize',14)
%}
%% plot line mean+SEM whisking and pupil for prev Trial conditions 
whichSessions = 10:18; %LEW032
clear pupilEarly pupilLate pupilHit pupilMiss whiskEarly whiskLate whiskHit whiskMiss
for iX = whichSessions
   [~,prevEarly] = selectCondition(expInfo(iX), [-1 -.5 -.12 -.05 0 .05 .12 .5 1],...
        behavioralData(iX), initTrialConditions('preStimMovement','quiescent',...
        'specificRTs',[.1 3],'pastMovementTime','early','movementTime','early'));
    [~,prevLate] = selectCondition(expInfo(iX), [-1 -.5 -.12 -.05 0 .05 .12 .5 1],...
        behavioralData(iX), initTrialConditions('preStimMovement','quiescent',...
        'specificRTs',([.1 3]),'pastMovementTime','late','movementTime','early'));
    [~,prevHit] = selectCondition(expInfo(iX), [-1 -.5 -.12 -.05 0 .05 .12 .5 1],...
        behavioralData(iX), initTrialConditions('preStimMovement','quiescent',...
        'specificRTs',([.1 3]),'pastResponseType','correct','movementTime','early'));
    [~,prevMiss] = selectCondition(expInfo(iX), [-1 -.5 -.12 -.05 0 .05 .12 .5 1],...
        behavioralData(iX), initTrialConditions('preStimMovement','quiescent',...
        'specificRTs',([.1 3]),'pastResponseType','incorrect','movementTime','early'));    
        
    pupilEarly(iX,:) = smooth(mean(eyeData(iX).eta.alignedFace{1}(prevEarly,:,1),1));
    pupilLate(iX,:) = smooth(mean(eyeData(iX).eta.alignedFace{1}(prevLate,:,1),1));
    whiskEarly(iX,:) = smooth(mean(eyeData(iX).eta.alignedFace{1}(prevEarly,:,2),1));
    whiskLate(iX,:) = smooth(mean(eyeData(iX).eta.alignedFace{1}(prevLate,:,2),1));
    pupilHit(iX,:) = smooth(mean(eyeData(iX).eta.alignedFace{1}(prevHit,:,1),1));
    pupilMiss(iX,:) = smooth(mean(eyeData(iX).eta.alignedFace{1}(prevMiss,:,1),1));
    whiskHit(iX,:) = smooth(mean(eyeData(iX).eta.alignedFace{1}(prevHit,:,2),1));
    whiskMiss(iX,:) = smooth(mean(eyeData(iX).eta.alignedFace{1}(prevMiss,:,2),1));   
end

pupilEarly = pupilEarly(whichSessions,:);
pupilLate = pupilLate(whichSessions,:);
whiskEarly = whiskEarly(whichSessions,:);
whiskLate = whiskLate(whichSessions,:);
pupilHit = pupilHit(whichSessions,:);
pupilMiss = pupilMiss(whichSessions,:);
whiskHit = whiskHit(whichSessions,:);
whiskMiss = whiskMiss(whichSessions,:);


%% plot pupil
clear allpupil pupilArray
allpupil = [pupilEarly; pupilLate];
pupilArray{1,1} = 1:size(pupilEarly,1);
pupilArray{2,1} = [1:size(pupilEarly,1)] + size(pupilEarly,1);
figure;
[meanPSTH, semPSTH, ~] = computePSTHs(allpupil,pupilArray);
plotPSTHs(eyeData(1).eta.eventWindow, cell2mat(meanPSTH), cell2mat(semPSTH), [.8 .4 0; 0 .5 .5],'-')
set(gca,'tickdir','out', 'Fontsize',14);
xlabel('Time from stimulus onset (s)')
ylabel('Pupil area (zscored)')
box off
y=ylim;
line([0 0],[y(1) y(2)],'LineStyle','--','Linewidth',1,'Color',[.7 .7 .7]);
%xlim([-1 1])
%legend({'','Prev early','','Prev late'},'Location','northwest',...
%    'Orientation','horizontal','Fontsize',12)
%legend boxoff

clear allpupil pupilArray
allpupil = [pupilHit; pupilMiss];
pupilArray{1,1} = 1:size(pupilHit,1);
pupilArray{2,1} = [1:size(pupilHit,1)] + size(pupilHit,1);
figure;
[meanPSTH, semPSTH, ~] = computePSTHs(allpupil,pupilArray);
plotPSTHs(eyeData(1).eta.eventWindow, cell2mat(meanPSTH), cell2mat(semPSTH), [.9 .2 0; 0 .6 .2],'-')
set(gca,'tickdir','out', 'Fontsize',14);
xlabel('Time from stimulus onset (s)')
ylabel('Pupil area (zscored)')
box off
y=ylim;
line([0 0],[y(1) y(2)],'LineStyle','--','Linewidth',1,'Color',[.7 .7 .7]);
%legend({'','Prev Reward','','Prev No Reward'},'Location','northwest',...
%    'Orientation','horizontal','Fontsize',12)
%legend boxoff
%xlim([-1 1])

%% plot whisking
clear allwhisk whiskArray
allwhisk = [whiskEarly; whiskLate];
whiskArray{1,1} = 1:size(whiskEarly,1);
whiskArray{2,1} = [1:size(whiskEarly,1)] + size(whiskEarly,1);
figure;
[meanPSTH, semPSTH, ~] = computePSTHs(allwhisk,whiskArray);
plotPSTHs(eyeData(1).eta.eventWindow, cell2mat(meanPSTH), cell2mat(semPSTH), [.8 .4 0; 0 .5 .5],'-')
set(gca,'tickdir','out','Fontsize',14);
xlabel('Time from stimulus onset (s)')
ylabel('Whisking (zscored)')
box off
y=ylim;
line([0 0],[y(1) y(2)],'LineStyle','--','Linewidth',1,'Color',[.7 .7 .7]);
%xlim([-1 1])
%legend({'','Prev early','','Prev late'},'Location','northwest',...
%    'Orientation','horizontal','Fontsize',12)
%legend boxoff

clear allwhisk whiskArray
allwhisk = [whiskHit; whiskMiss];
whiskArray{1,1} = 1:size(whiskHit,1);
whiskArray{2,1} = [1:size(whiskHit,1)] + size(whiskHit,1);
figure;
[meanPSTH, semPSTH, ~] = computePSTHs(allwhisk,whiskArray);
plotPSTHs(eyeData(1).eta.eventWindow, cell2mat(meanPSTH), cell2mat(semPSTH), [.9 .2 0; 0 .6 .2],'-')
set(gca,'tickdir','out','Fontsize',14);
xlabel('Time from stimulus onset (s)')
ylabel('Whisking (zscored)')
box off
y=ylim;
line([0 0],[y(1) y(2)],'LineStyle','--','Linewidth',1,'Color',[.7 .7 .7]);
%legend({'','Prev Reward','','Prev No Reward'},'Location','northwest',...
%    'Orientation','horizontal','Fontsize',12)
%legend boxoff
%xlim([-1 1])