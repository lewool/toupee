%% overview

% 

%% load experiment details
cd('C:\Users\Wool\Documents\GitHub\toupee');
% close all;
clearvars -except hiR_histoValues hiL_histoValues stR_histoValues stL_histoValues 

mouseName = 'LEW008';
expDate = '2019-01-29';
expNum = 1;
expSeries = [1];

%% load data

[block, Timeline] = loadData(mouseName, expDate, expNum);

%% get event timings and wheel trajectories
cd('C:\Users\Wool\Documents\GitHub\toupee\behavior');

signalsNames = {'stimulusOnTimes' 'interactiveOnTimes' 'stimulusOffTimes'};
[eventTimes, wheelTrajectories] = getEventTimes(block, Timeline, signalsNames);

%% load traces
cd('C:\Users\Wool\Documents\GitHub\toupee');

[allFcell, ops] = loadExpTraces(mouseName, expDate, expSeries);

%% align calcium traces to the event you want

% cut the trace into trial-by-trial traces, aligned to a particular event
% event = 'prestimulusQuiescenceEndTimes';
event = 'stimulusOnTimes';
% event = 'rewardOnTimes';
[alignedTraces, eventWindow] = getExpTraces(mouseName, expDate, expNum, expSeries, allFcell, eventTimes, ops, event);
% [alignedTraces, alignedSpikes, eventWindow] = alignSpikes(mouseName, expDate, expNum, expSeries, allFcell, eventTimes, ops, event);

%% select cells with the properties you want
plotAll = chooseCellType('all', mouseName, expDate, expNum, expSeries, block, allFcell, eventTimes, ops);

%% initialize some plot values

%imaging rate
numPlanes = length(alignedTraces);
Fs = 15;%/ numPlanes;

%event window
eventIdx = find(eventWindow == 0);

%stimulus response index
crfTime = 0.1;
vPreIdx = eventIdx - 1;
vPeriIdx = eventIdx + ceil(crfTime*Fs);

responseWindow = eventIdx+6;
contrasts = unique(block.events.contrastValues);

%% report the mean response of each cell to your event, per trial
%output is a matrix of size trials x cells

numCompleteTrials = size(alignedTraces{1}.eventSpikes,1);
onsetResps = zeros(numCompleteTrials,size(plotAll,1));
for k = 1:size(plotAll,1)
    onsetResps(:,k) = squeeze(nanmean(alignedTraces{plotAll(k,1)}.eventSpikes(:,responseWindow,plotAll(k,2)),2));
end

someNaNs = find(isnan(onsetResps));
badCells = unique(ceil(someNaNs/numCompleteTrials));
for r = fliplr(badCells)
    onsetResps(:,r) = [];
end
%% 
% cd('C:\Users\Wool\Documents\GitHub\toupee\behavior');
% [condLogical, condIdxLeft] = selectCondition(block,contrasts, eventTimes,'all','all','all','left','all','all','all','all','all');
% [condLogical, condIdxRight] = selectCondition(block,contrasts, eventTimes,'all','all','all','right','all','all','all','all','all');
% 
% onsetResps_hiLeft = onsetResps(condIdxLeft,:);
% onsetResps_Ltrain = onsetResps_hiLeft(2:2:end,:);
% onsetResps_Ltest = onsetResps_hiLeft(1:2:end,:);
% 
% onsetResps_hiRight = onsetResps(condIdxRight,:);
% onsetResps_Rtrain = onsetResps_hiRight(2:2:end,:);
% onsetResps_Rtest = onsetResps_hiRight(1:2:end,:);
% 
% onsetResps_meanLeft = nanmean(onsetResps_Ltrain,1);
% onsetResps_meanRight = nanmean(onsetResps_Rtrain,1);
% 
% for iTrial = 1:size(onsetResps_Ltest,1)
%     LdotsL(iTrial,:) = dot(onsetResps_Ltest(iTrial,:),onsetResps_meanLeft);
%     LdotsR(iTrial,:) = dot(onsetResps_Ltest(iTrial,:),onsetResps_meanRight);
% end
% for iTrial = 1:size(onsetResps_Rtest,1)
%     RdotsR(iTrial,:) = dot(onsetResps_Rtest(iTrial,:),onsetResps_meanRight);
%     RdotsL(iTrial,:) = dot(onsetResps_Rtest(iTrial,:),onsetResps_meanLeft);
% end
% 
% maxLim = 10 * round((max([max(RdotsR) max(RdotsL) max(LdotsR) max(LdotsL)]) * 1.1)/10);
% minLim = 10 * round((min([min(RdotsR) min(RdotsL) min(LdotsR) min(LdotsL)]) * 1.1)/10);
% minLim = -150;
% maxLim = 800;
% stepLim = round((maxLim-minLim)/40)*10;
% figure(100);
% ax = gca;
% hold on;
% axis square;
% plotLine = line([minLim maxLim],[minLim maxLim]);
% set(plotLine,'Color', [.5 .5 .5],'LineStyle','--')
% plotL = scatter(LdotsL,LdotsR,8,'b');
% plotR = scatter(RdotsL,RdotsR,8,'ro');
% set(plotL,'Marker','o','MarkerEdgeColor',[0 102 255]./255);
% axis([minLim maxLim minLim maxLim]);
% plotL.MarkerEdgeAlpha = 0.75;
% hold on;
% plotR.MarkerEdgeAlpha = 0.75;
% ax.TickDir = 'out';
% ylabel('x_i · xbar_R')
% xlabel('x_i · xbar_R')
% % xticks([minLim:stepLim:maxLim])
% % yticks([minLim:stepLim:maxLim])
% title(strcat(mouseName,{' '},expDate,{' '},num2str(eventWindow(responseWindow))));

%%

h =  findobj('type','figure');
n = length(h);

cd('C:\Users\Wool\Documents\GitHub\toupee\behavior');
block = expInfo.block;
numCompleteTrials = size(onsetResps,1);
[condLogical, condIdxLeftHigh] = selectCondition(block,contrasts, eventTimes,trialConditions{1});
[condLogical, condIdxRightHigh] = selectCondition(block,contrasts, eventTimes,trialConditions{2});
[condLogical, condIdxLeftStim] = selectCondition(block,contrasts(contrasts<0), eventTimes,trialConditions{3});
[condLogical, condIdxRightStim] = selectCondition(block,contrasts(contrasts>0), eventTimes,trialConditions{3});

[condLogical, condIdxStimL_highL] = selectCondition(block,contrasts(contrasts<0), eventTimes,trialConditions{1});
[condLogical, condIdxStimR_highL] = selectCondition(block,contrasts(contrasts>0), eventTimes,trialConditions{1});
[condLogical, condIdxStimL_highR] = selectCondition(block,contrasts(contrasts<0), eventTimes,trialConditions{2});
[condLogical, condIdxStimR_highR] = selectCondition(block,contrasts(contrasts>0), eventTimes,trialConditions{2});

% [condLogical, condIdxStim0_highL] = selectCondition(block,contrasts(contrasts==0), eventTimes,'all','all','all','left','correct','all','all','all','all');
% [condLogical, condIdxStim0_highR] = selectCondition(block,contrasts(contrasts==0), eventTimes,'all','all','all','right','correct','all','all','all','all');


trainTrials = 2:2:numCompleteTrials;
testTrials = 1:2:numCompleteTrials;

highLeft_trainTrials = intersect(trainTrials,condIdxLeftHigh);
highRight_trainTrials = intersect(trainTrials,condIdxRightHigh);
stimLeft_trainTrials = intersect(trainTrials,condIdxLeftStim);
stimRight_trainTrials = intersect(trainTrials,condIdxRightStim);

stimL_highL = intersect(testTrials,condIdxStimL_highL);
stimR_highL = intersect(testTrials,condIdxStimR_highL);
stimL_highR = intersect(testTrials,condIdxStimL_highR);
stimR_highR = intersect(testTrials,condIdxStimR_highR);

onsetResps_meanLeft_high = nanmean(onsetResps(highLeft_trainTrials,:),1);
onsetResps_meanRight_high = nanmean(onsetResps(highRight_trainTrials,:),1);

clear LLdotsL LLdotsR LRdotsL LRdotsR RLdotsL RLdotsR RRdotsL RRdotsR

for iTrial = stimL_highL
    LLdotsL(iTrial,:) = dot(onsetResps(iTrial,:),onsetResps_meanLeft_high)/(norm(onsetResps(iTrial,:))*norm(onsetResps_meanLeft_high));
    LLdotsR(iTrial,:) = dot(onsetResps(iTrial,:),onsetResps_meanRight_high)/(norm(onsetResps(iTrial,:))*norm(onsetResps_meanRight_high));
end

for iTrial = stimR_highL
    RLdotsR(iTrial,:) = dot(onsetResps(iTrial,:),onsetResps_meanRight_high)/(norm(onsetResps(iTrial,:))*norm(onsetResps_meanRight_high));
    RLdotsL(iTrial,:) = dot(onsetResps(iTrial,:),onsetResps_meanLeft_high)/(norm(onsetResps(iTrial,:))*norm(onsetResps_meanLeft_high));
end

for iTrial = stimL_highR
    LRdotsL(iTrial,:) = dot(onsetResps(iTrial,:),onsetResps_meanLeft_high)/(norm(onsetResps(iTrial,:))*norm(onsetResps_meanLeft_high));
    LRdotsR(iTrial,:) = dot(onsetResps(iTrial,:),onsetResps_meanRight_high)/(norm(onsetResps(iTrial,:))*norm(onsetResps_meanRight_high));
end

for iTrial = stimR_highR
    RRdotsR(iTrial,:) = dot(onsetResps(iTrial,:),onsetResps_meanRight_high)/(norm(onsetResps(iTrial,:))*norm(onsetResps_meanRight_high));
    RRdotsL(iTrial,:) = dot(onsetResps(iTrial,:),onsetResps_meanLeft_high)/(norm(onsetResps(iTrial,:))*norm(onsetResps_meanLeft_high));
end

mSize = 15;
% hiR_L = [0 .5 1];
% hiL_R = [1 .35 .35];
% hiL_L = [0 0 1];
% hiR_R = [.65 0 0 ];

reward = [0 0 0];
axMax = 1.5*max([max(LLdotsL) max(LLdotsR) max(RLdotsL) max(RLdotsR) max(LRdotsL) max(LRdotsR) max(RRdotsL) max(RRdotsR)]);
axMin = 1.5*min([max(LLdotsL) min(LLdotsR) min(RLdotsL) min(RLdotsR) min(LRdotsL) min(LRdotsR) min(RRdotsL) min(RRdotsR)]);

LLdotsL(LLdotsL == 0) = [];
LLdotsR(LLdotsR == 0) = [];
RLdotsL(RLdotsL == 0) = [];
RLdotsR(RLdotsR == 0) = [];
LRdotsL(LRdotsL == 0) = [];
LRdotsR(LRdotsR == 0) = [];
RRdotsL(RRdotsL == 0) = [];
RRdotsR(RRdotsR == 0) = [];

figure(n+2);
subplot(2,2,1);
hiR_L = [1 0 0];
hiL_R = [0 .4 1];
hiL_L = [0 .4 1];
hiR_R = [1 0 0];
ax1 = gca;
hold on;
axis square;
plotLine = line([axMin axMax],[axMin axMax]);
set(plotLine,'Color', [.5 .5 .5],'LineStyle','--')

plotLL = scatter(LLdotsL,LLdotsR,mSize,'bo');
set(plotLL,'MarkerEdgeColor',hiL_L, 'Marker','o','MarkerFaceColor',hiL_L,'MarkerFaceAlpha',.4)
plotRL = scatter(RLdotsL,RLdotsR,mSize,'bo');
set(plotRL,'MarkerEdgeColor',hiL_R, 'Marker','o','MarkerFaceColor',hiL_R,'MarkerFaceAlpha',.4)
plotLR = scatter(LRdotsL,LRdotsR,mSize,'ro');
set(plotLR,'MarkerEdgeColor',hiR_L, 'Marker','o','MarkerFaceColor',hiR_L,'MarkerFaceAlpha',.4)
plotRR = scatter(RRdotsL,RRdotsR,mSize,'ro');
set(plotRR,'MarkerEdgeColor',hiR_R, 'Marker','o','MarkerFaceColor',hiR_R,'MarkerFaceAlpha',.4)

axis([axMin axMax axMin axMax]);
hold on;
ax1.TickDir = 'out';
ylabel('similarity to high-right mean activity')
xlabel('similarity to high-left mean activity')
% xticks([minLim:stepLim:maxLim])
% yticks([minLim:stepLim:maxLim])
title('trials split by reward block');

subplot(2,2,2);
hiL_R = [1 0 0];
hiR_L = [0 .4 1];
hiL_L = [0 .4 1];
hiR_R = [1 0 0];
ax1 = gca;
hold on;
axis square;
plotLine = line([axMin axMax],[axMin axMax]);
set(plotLine,'Color', [.5 .5 .5],'LineStyle','--')

plotLL = scatter(LLdotsL,LLdotsR,mSize,'bo');
set(plotLL,'MarkerEdgeColor',hiL_L, 'Marker','o','MarkerFaceColor',hiL_L,'MarkerFaceAlpha',.4)
plotRL = scatter(RLdotsL,RLdotsR,mSize,'bo');
set(plotRL,'MarkerEdgeColor',hiL_R, 'Marker','o','MarkerFaceColor',hiL_R,'MarkerFaceAlpha',.4)
plotLR = scatter(LRdotsL,LRdotsR,mSize,'ro');
set(plotLR,'MarkerEdgeColor',hiR_L, 'Marker','o','MarkerFaceColor',hiR_L,'MarkerFaceAlpha',.4)
plotRR = scatter(RRdotsL,RRdotsR,mSize,'ro');
set(plotRR,'MarkerEdgeColor',hiR_R, 'Marker','o','MarkerFaceColor',hiR_R,'MarkerFaceAlpha',.4)

axis([axMin axMax axMin axMax]);
hold on;
ax1.TickDir = 'out';
ylabel('similarity to high-right mean activity')
xlabel('similarity to high-left mean activity')
% xticks([minLim:stepLim:maxLim])
% yticks([minLim:stepLim:maxLim])
title('trials split by stimulus side');

figure(n+3);
subplot(2,2,1);
ax2 = gca;
% set(ax2,'xtick',[0])
set(ax2,'ytick',[])
hold on
plotHiR = histfit([RRdotsL;LRdotsL]-[RRdotsR;LRdotsR],[],'kernel');
plotHiR(1).FaceAlpha = 0;
plotHiR(1).EdgeAlpha = 0;
plotHiL = histfit([LLdotsL;RLdotsL]-[LLdotsR;RLdotsR],[],'kernel');
plotHiL(1).FaceAlpha = 0;
plotHiL(2).Color = hiL_L;
plotHiL(1).EdgeAlpha = 0;

hiR_histoValues{length(hiR_histoValues)+1} = [RRdotsL;LRdotsL]-[RRdotsR;LRdotsR];
hiL_histoValues{length(hiL_histoValues)+1} = [LLdotsL;RLdotsL]-[LLdotsR;RLdotsR];

subplot(2,2,2);
ax2 = gca;
% set(ax2,'xtick',[0])
set(ax2,'ytick',[])
hold on
plotHiR = histfit([RRdotsL;RLdotsL]-[RRdotsR;RLdotsR],[],'kernel');
plotHiR(1).FaceAlpha = 0;
plotHiR(1).EdgeAlpha = 0;
plotHiL = histfit([LLdotsL;LRdotsL]-[LLdotsR;LRdotsR],[],'kernel');
plotHiL(1).FaceAlpha = 0;
plotHiL(2).Color = hiL_L;
plotHiL(1).EdgeAlpha = 0;


onsetResps_meanLeft_stim = nanmean(onsetResps(stimLeft_trainTrials,:),1);
onsetResps_meanRight_stim = nanmean(onsetResps(stimRight_trainTrials,:),1);

clear LLdotsL LLdotsR LRdotsL LRdotsR RLdotsL RLdotsR RRdotsL RRdotsR

for iTrial = stimL_highL
    LLdotsL(iTrial,:) = dot(onsetResps(iTrial,:),onsetResps_meanLeft_stim)/(norm(onsetResps(iTrial,:))*norm(onsetResps_meanLeft_stim));
    LLdotsR(iTrial,:) = dot(onsetResps(iTrial,:),onsetResps_meanRight_stim)/(norm(onsetResps(iTrial,:))*norm(onsetResps_meanRight_stim));
end

for iTrial = stimR_highL
    RLdotsR(iTrial,:) = dot(onsetResps(iTrial,:),onsetResps_meanRight_stim)/(norm(onsetResps(iTrial,:))*norm(onsetResps_meanRight_stim));
    RLdotsL(iTrial,:) = dot(onsetResps(iTrial,:),onsetResps_meanLeft_stim)/(norm(onsetResps(iTrial,:))*norm(onsetResps_meanLeft_stim));
end

for iTrial = stimL_highR
    LRdotsL(iTrial,:) = dot(onsetResps(iTrial,:),onsetResps_meanLeft_stim)/(norm(onsetResps(iTrial,:))*norm(onsetResps_meanLeft_stim));
    LRdotsR(iTrial,:) = dot(onsetResps(iTrial,:),onsetResps_meanRight_stim)/(norm(onsetResps(iTrial,:))*norm(onsetResps_meanRight_stim));
end

for iTrial = stimR_highR
    RRdotsR(iTrial,:) = dot(onsetResps(iTrial,:),onsetResps_meanRight_stim)/(norm(onsetResps(iTrial,:))*norm(onsetResps_meanRight_stim));
    RRdotsL(iTrial,:) = dot(onsetResps(iTrial,:),onsetResps_meanLeft_stim)/(norm(onsetResps(iTrial,:))*norm(onsetResps_meanLeft_stim));
end

LLdotsL(LLdotsL == 0) = [];
LLdotsR(LLdotsR == 0) = [];
RLdotsL(RLdotsL == 0) = [];
RLdotsR(RLdotsR == 0) = [];
LRdotsL(LRdotsL == 0) = [];
LRdotsR(LRdotsR == 0) = [];
RRdotsL(RRdotsL == 0) = [];
RRdotsR(RRdotsR == 0) = [];

figure(n+2);
subplot(2,2,3);
hiR_L = [1 0 0];
hiL_R = [0 .4 1];
hiL_L = [0 .4 1];
hiR_R = [1 0 0];
ax1 = gca;
hold on;
axis square;
plotLine = line([axMin axMax],[axMin axMax]);
set(plotLine,'Color', [.5 .5 .5],'LineStyle','--')

plotLL = scatter(LLdotsL,LLdotsR,mSize,'bo');
set(plotLL,'MarkerEdgeColor',hiL_L, 'Marker','o','MarkerFaceColor',hiL_L,'MarkerFaceAlpha',.4)
plotRL = scatter(RLdotsL,RLdotsR,mSize,'bo');
set(plotRL,'MarkerEdgeColor',hiL_R, 'Marker','o','MarkerFaceColor',hiL_R,'MarkerFaceAlpha',.4)
plotLR = scatter(LRdotsL,LRdotsR,mSize,'ro');
set(plotLR,'MarkerEdgeColor',hiR_L, 'Marker','o','MarkerFaceColor',hiR_L,'MarkerFaceAlpha',.4)
plotRR = scatter(RRdotsL,RRdotsR,mSize,'ro');
set(plotRR,'MarkerEdgeColor',hiR_R, 'Marker','o','MarkerFaceColor',hiR_R,'MarkerFaceAlpha',.4)

axis([axMin axMax axMin axMax]);
hold on;
ax1.TickDir = 'out';
ylabel('similarity to stim-right mean activity')
xlabel('similarity to stim-left mean activity')
% xticks([minLim:stepLim:maxLim])
% yticks([minLim:stepLim:maxLim])
title('trials split by reward block');

subplot(2,2,4);
hiL_R = [1 0 0];
hiR_L = [0 .4 1];
hiL_L = [0 .4 1];
hiR_R = [1 0 0];
ax1 = gca;
hold on;
axis square;
plotLine = line([axMin axMax],[axMin axMax]);
set(plotLine,'Color', [.5 .5 .5],'LineStyle','--')

plotLL = scatter(LLdotsL,LLdotsR,mSize,'bo');
set(plotLL,'MarkerEdgeColor',hiL_L, 'Marker','o','MarkerFaceColor',hiL_L,'MarkerFaceAlpha',.4)
plotRL = scatter(RLdotsL,RLdotsR,mSize,'bo');
set(plotRL,'MarkerEdgeColor',hiL_R, 'Marker','o','MarkerFaceColor',hiL_R,'MarkerFaceAlpha',.4)
plotLR = scatter(LRdotsL,LRdotsR,mSize,'ro');
set(plotLR,'MarkerEdgeColor',hiR_L, 'Marker','o','MarkerFaceColor',hiR_L,'MarkerFaceAlpha',.4)
plotRR = scatter(RRdotsL,RRdotsR,mSize,'ro');
set(plotRR,'MarkerEdgeColor',hiR_R, 'Marker','o','MarkerFaceColor',hiR_R,'MarkerFaceAlpha',.4)

axis([axMin axMax axMin axMax]);
hold on;
ax1.TickDir = 'out';
ylabel('similarity to stim-right mean activity')
xlabel('similarity to stim-left mean activity')
% xticks([minLim:stepLim:maxLim])
% yticks([minLim:stepLim:maxLim])
title('trials split by stimulus side');

figure(n+3);
subplot(2,2,3);
ax2 = gca;
% set(ax2,'xtick',[0])
set(ax2,'ytick',[])
hold on
plotHiR = histfit([RRdotsL;LRdotsL]-[RRdotsR;LRdotsR],[],'kernel');
plotHiR(1).FaceAlpha = 0;
plotHiR(1).EdgeAlpha = 0;
plotHiL = histfit([LLdotsL;RLdotsL]-[LLdotsR;RLdotsR],[],'kernel');
plotHiL(1).FaceAlpha = 0;
plotHiL(2).Color = hiL_L;
plotHiL(1).EdgeAlpha = 0;

subplot(2,2,4);
ax2 = gca;
% set(ax2,'xtick',[0])
set(ax2,'ytick',[])
hold on
plotHiR = histfit([RRdotsL;RLdotsL]-[RRdotsR;RLdotsR],[],'kernel');
plotHiR(1).FaceAlpha = 0;
plotHiR(1).EdgeAlpha = 0;
plotHiL = histfit([LLdotsL;LRdotsL]-[LLdotsR;LRdotsR],[],'kernel');
plotHiL(1).FaceAlpha = 0;
plotHiL(2).Color = hiL_L;
plotHiL(1).EdgeAlpha = 0;

stR_histoValues{length(stR_histoValues)+1} = [RRdotsL;RLdotsL]-[RRdotsR;RLdotsR];
stL_histoValues{length(stL_histoValues)+1} = [LLdotsL;LRdotsL]-[LLdotsR;LRdotsR];


%% trajectories
cd('C:\Users\Wool\Documents\GitHub\toupee\behavior');
block = expInfo.block;
numCompleteTrials = size(onsetResps,1);

colorLightRed = [1 0 0];
colorLightBlue = [0 .4 1];
colorDarkRed = [.5 0 0];
colorDarkBlue = [0 0 1];

% 1. get trial-type indices 
trialConditions{1} = initTrialConditions('highRewardSide','left','responseType','correct','movementTime','all');
trialConditions{2} = initTrialConditions('highRewardSide','right','responseType','correct','movementTime','all');
trialConditions{3} = initTrialConditions('responseType','correct','movementTime','all');


[condLogical, condIdxLeftHigh] = selectCondition(block,contrasts(contrasts~=0), eventTimes,trialConditions{1});
[condLogical, condIdxRightHigh] = selectCondition(block,contrasts(contrasts~=0), eventTimes,trialConditions{2});
[condLogical, condIdxLeftStim] = selectCondition(block,contrasts(contrasts<0), eventTimes,trialConditions{3});
[condLogical, condIdxRightStim] = selectCondition(block,contrasts(contrasts>0), eventTimes,trialConditions{3});

[condLogical, condIdxStimL_highL] = selectCondition(block,contrasts(contrasts<0), eventTimes,trialConditions{1});
[condLogical, condIdxStimR_highL] = selectCondition(block,contrasts(contrasts>0), eventTimes,trialConditions{1});
[condLogical, condIdxStimL_highR] = selectCondition(block,contrasts(contrasts<0), eventTimes,trialConditions{2});
[condLogical, condIdxStimR_highR] = selectCondition(block,contrasts(contrasts>0), eventTimes,trialConditions{2});

trainTrials = 2:2:numCompleteTrials;
testTrials = 1:2:numCompleteTrials;


% 2. first pick a 'mother time' amd report the mean response of each cell @it, per trial
% output is a matrix of size trials x cells
% eventIdx = find(eventWindow == 0);
% motherTime = 1.3;
% motherIdx = eventIdx + ceil(motherTime*Fs);

numCompleteTrials = size(alignedResps{1},1);
motherResps = zeros(numCompleteTrials,size(plotCells,2));
for k = 1:size(plotCells,2)
    motherResps(:,k) = squeeze(nanmean(alignedResps{2}(:,motherIdx,plotCells(k)),2));
end

someNaNs = find(isnan(motherResps));
badCells = unique(ceil(someNaNs/numCompleteTrials));
for r = fliplr(badCells)
    motherResps(:,r) = [];
end

% 3. compute the mean vectors upon which to project trial-by-trial vectors
motherResps_reward = nanmean(motherResps(intersect(trainTrials,condIdxLeftHigh),:),1) - nanmean(motherResps(intersect(trainTrials,condIdxRightHigh),:),1);
motherResps_stim = nanmean(motherResps(intersect(trainTrials,condIdxLeftStim),:),1) - nanmean(motherResps(intersect(trainTrials,condIdxRightStim),:),1);

% 4. pick another time window to evaluate the responses of cells on
% individual trials (this will be projected on the mother vector)

% 4a. plot a few time points as scatterplots

% indices = [1 5 10 15 20];
% figure;
% set(gcf,'Position',[100 100 1795 645]);
% set(gcf,'renderer','Painters')
% hold on
% for p = 1:length(indices)
%     projIdx = eventIdx + indices(p);
% 
% projResps = zeros(numCompleteTrials,size(plotCells,2));
% for k = 1:size(plotCells,2)
%     projResps(:,k) = squeeze(nanmean(alignedResps{2}(:,projIdx,plotCells(k)),2));
% end
% someNaNs = find(isnan(projResps));
% badCells = unique(ceil(someNaNs/numCompleteTrials));
% for r = fliplr(badCells)
%     projResps(:,r) = [];
% end
% 
% for iTrial = testTrials
%     projDotReward(iTrial,:) = dot(projResps(iTrial,:),motherResps_reward)/(norm(projResps(iTrial,:))*norm(motherResps_reward));
%     projDotStim(iTrial,:) = dot(projResps(iTrial,:),motherResps_stim)/(norm(projResps(iTrial,:))*norm(motherResps_stim));
% end
% 
% subplot(2,5,p)
% hold on
% plotStimLTrials = scatter(projDotStim(intersect(testTrials,condIdxLeftStim)),projDotReward(intersect(testTrials,condIdxLeftStim)));
% set(plotStimLTrials,'MarkerEdgeColor',colorLightBlue, 'Marker','o','MarkerFaceColor',colorLightBlue,'MarkerFaceAlpha',.4)
% plotStimRTrials = scatter(projDotStim(intersect(testTrials,condIdxRightStim)),projDotReward(intersect(testTrials,condIdxRightStim)));
% set(plotStimRTrials,'MarkerEdgeColor',colorLightRed, 'Marker','o','MarkerFaceColor',colorLightRed,'MarkerFaceAlpha',.4)
% xlim([-.4 .4]);
% ylim([-.4 .4]);
% yticks([-.4 -.2 0 .2 .4])
% xticks([-.4 -.2 0 .2 .4])
% axis square
% if p == 1
%     ylabel('highL - highR')
%     legend({'left stim','right stim'})
% end
% xlabel('stimL - stimR')
% 
% subplot(2,5,p+5)
% hold on
% plotStimLTrials = scatter(projDotStim(intersect(testTrials,condIdxLeftHigh)),projDotReward(intersect(testTrials,condIdxLeftHigh)));
% set(plotStimLTrials,'MarkerEdgeColor',colorLightBlue, 'Marker','o','MarkerFaceColor',colorLightBlue,'MarkerFaceAlpha',.4)
% plotStimRTrials = scatter(projDotStim(intersect(testTrials,condIdxRightHigh)),projDotReward(intersect(testTrials,condIdxRightHigh)));
% set(plotStimRTrials,'MarkerEdgeColor',colorLightRed, 'Marker','o','MarkerFaceColor',colorLightRed,'MarkerFaceAlpha',.4)
% xlim([-.4 .4]);
% ylim([-.4 .4]);
% yticks([-.4 -.2 0 .2 .4])
% xticks([-.4 -.2 0 .2 .4])
% axis square
% if p == 1
%     ylabel('highL - highR')
%     legend({'left high','right high'})
% end
% xlabel('stimL - stimR')
% 
% clear projDotReward projDotStim
% end


% 4b. plot every timepoint and collect a trajectory vector

for p = 1:length(eventWindow)

projResps =  squeeze(alignedResps{2}(:,p,plotCells));

someNaNs = find(isnan(projResps));
badCells = unique(ceil(someNaNs/numCompleteTrials));
projResps(:,badCells) = [];

for iTrial = testTrials
    projDotReward(iTrial,:) = dot(projResps(iTrial,:),motherResps_reward)/(norm(projResps(iTrial,:))*norm(motherResps_reward));
    projDotStim(iTrial,:) = dot(projResps(iTrial,:),motherResps_stim)/(norm(projResps(iTrial,:))*norm(motherResps_stim));
end

leftRewardsTrajectory(p,1) = mean(projDotStim(intersect(testTrials,condIdxLeftHigh)));
leftRewardsTrajectory(p,2) = mean(projDotReward(intersect(testTrials,condIdxLeftHigh)));

rightRewardsTrajectory(p,1) = mean(projDotStim(intersect(testTrials,condIdxRightHigh)));
rightRewardsTrajectory(p,2) = mean(projDotReward(intersect(testTrials,condIdxRightHigh)));

leftStimsTrajectory(p,1) = mean(projDotStim(intersect(testTrials,condIdxLeftStim)));
leftStimsTrajectory(p,2) = mean(projDotReward(intersect(testTrials,condIdxLeftStim)));

rightStimsTrajectory(p,1) = mean(projDotStim(intersect(testTrials,condIdxRightStim)));
rightStimsTrajectory(p,2) = mean(projDotReward(intersect(testTrials,condIdxRightStim)));

end
%
figure;
hold on
set(gcf,'Position',[100 100 1025 420]);
set(gcf,'renderer','Painters');
subplot(1,2,1)
hold on;
plot(leftStimsTrajectory(:,1),leftStimsTrajectory(:,2),'Color',colors(1,:),'LineWidth',1);
plot(rightStimsTrajectory(:,1),rightStimsTrajectory(:,2),'Color','r','LineWidth',1);
% plot(leftStimsTrajectory(1,1),leftStimsTrajectory(1,2),'ko','MarkerFaceColor','k')
% plot(leftStimsTrajectory(eventIdx,1),leftStimsTrajectory(eventIdx,2),'ko','MarkerFaceColor',[0 .6 .12])
% plot(leftStimsTrajectory(end,1),leftStimsTrajectory(end,2),'ko','MarkerFaceColor',colorDarkRed)
% plot(rightStimsTrajectory(1,1),rightStimsTrajectory(1,2),'ko','MarkerFaceColor','k')
% plot(rightStimsTrajectory(eventIdx,1),rightStimsTrajectory(eventIdx,2),'ko','MarkerFaceColor',[0 .6 .12])
% plot(rightStimsTrajectory(end,1),rightStimsTrajectory(end,2),'ko','MarkerFaceColor',colorDarkRed)
xlim([-.25 .25]);
xlim([-.25 .25]);
ylim([-.25 .25]);
axis square
ylabel('highL - highR')
% legend({'left stim','right stim'})
xlabel('stimL - stimR')

subplot(1,2,2);
hold on;
plot(leftRewardsTrajectory(:,1),leftRewardsTrajectory(:,2),'Color',colors(3,:),'LineWidth',1);
plot(rightRewardsTrajectory(:,1),rightRewardsTrajectory(:,2),'Color',colors(4,:),'LineWidth',1);
% plot(leftRewardsTrajectory(1,1),leftRewardsTrajectory(1,2),'ko','MarkerFaceColor','k')
% plot(leftRewardsTrajectory(eventIdx,1),leftRewardsTrajectory(eventIdx,2),'ko','MarkerFaceColor',[0 .6 .12])
% plot(leftRewardsTrajectory(end,1),leftRewardsTrajectory(end,2),'ko','MarkerFaceColor',colorDarkRed)
% plot(rightRewardsTrajectory(1,1),rightRewardsTrajectory(1,2),'ko','MarkerFaceColor','k')
% plot(rightRewardsTrajectory(eventIdx,1),rightRewardsTrajectory(eventIdx,2),'ko','MarkerFaceColor',[0 .6 .12])
% plot(rightRewardsTrajectory(end,1),rightRewardsTrajectory(end,2),'ko','MarkerFaceColor',colorDarkRed)
xlim([-.25 .25]);
xlim([-.25 .25]);
ylim([-.25 .25]);
axis square
ylabel('highL - highR')
% legend({'left high','right high'})
xlabel('stimL - stimR')

%%

colors = [0 .4 1; 1 0 0; 0.1 0.7 0.1; 1 .6 0];

leftRewardsTrajectory_int(:,1) = interp1(eventWindow, leftRewardsTrajectory(:,1)', linspace(eventWindow(1),eventWindow(end), 50));
leftRewardsTrajectory_int(:,2) = interp1(eventWindow, leftRewardsTrajectory(:,2)', linspace(eventWindow(1),eventWindow(end), 50));
rightRewardsTrajectory_int(:,1) = interp1(eventWindow, rightRewardsTrajectory(:,1)', linspace(eventWindow(1),eventWindow(end), 50));
rightRewardsTrajectory_int(:,2) = interp1(eventWindow, rightRewardsTrajectory(:,2)', linspace(eventWindow(1),eventWindow(end), 50));
leftStimsTrajectory_int(:,1) = interp1(eventWindow, leftStimsTrajectory(:,1)', linspace(eventWindow(1),eventWindow(end), 50));
leftStimsTrajectory_int(:,2) = interp1(eventWindow, leftStimsTrajectory(:,2)', linspace(eventWindow(1),eventWindow(end), 50));
rightStimsTrajectory_int(:,1) = interp1(eventWindow, rightStimsTrajectory(:,1)', linspace(eventWindow(1),eventWindow(end), 50));
rightStimsTrajectory_int(:,2) = interp1(eventWindow, rightStimsTrajectory(:,2)', linspace(eventWindow(1),eventWindow(end), 50));
%%
figure(4);
set(gcf,'color','w');

subplot(1,2,2)
ax1 = gca;
% ax1.TickDir = 'out';
plot(0,0)
hold on
xlim([-.25 .25]);
xlim([-.25 .25]);
ylim([-.25 .25]);
axis square
box off

cd('C:\Users\Wool\Desktop\tempFigs')
vidObj1 = VideoWriter('movie.avi');
open(vidObj1);
loops = 50;
for f = 1:loops-1
    
    plot(leftRewardsTrajectory_int(f:f+1,1),leftRewardsTrajectory_int(f:f+1,2),'Color',colors(3,:),'LineWidth',1);
    plot(rightRewardsTrajectory_int(f:f+1,1),rightRewardsTrajectory_int(f:f+1,2),'Color',colors(4,:),'LineWidth',1);
    if f == 17
        plot(leftRewardsTrajectory_int(f,1),leftRewardsTrajectory_int(f,2),'ko','MarkerFaceColor',[0 .6 .12])
        plot(rightRewardsTrajectory_int(f,1),rightRewardsTrajectory_int(f,2),'ko','MarkerFaceColor',[0 .6 .12])
    end

    hold on
    F1(f) = getframe(gcf);
    writeVideo(vidObj1,F1(f));
end

subplot(1,2,1)
plot(0,0)
hold on
xlim([-.25 .25]);
xlim([-.25 .25]);
ylim([-.25 .25]);
axis square
box off
ax2 = gca;
% ax2.TickDir = 'out';
% cd('C:\Users\Wool\Desktop\tempFigs')
% vidObj2 = VideoWriter('stim.avi');
% open(vidObj2);
loops = 50;
for f = 1:loops-1
    plot(leftStimsTrajectory_int(f:f+1,1),leftStimsTrajectory_int(f:f+1,2),'Color',colors(1,:),'LineWidth',1);
    plot(rightStimsTrajectory_int(f:f+1,1),rightStimsTrajectory_int(f:f+1,2),'Color',colors(2,:),'LineWidth',1);
    if f == 17
        plot(leftStimsTrajectory_int(f,1),leftStimsTrajectory_int(f,2),'ko','MarkerFaceColor',[0 .6 .12])
        plot(rightStimsTrajectory_int(f,1),rightStimsTrajectory_int(f,2),'ko','MarkerFaceColor',[0 .6 .12])
    end

    hold on
    F1(f) = getframe(gcf);
    writeVideo(vidObj1,F1(f));
end


%%
fig = figure;
movie(fig,F1,1)
cd('C:\Users\Wool\Desktop\tempFigs')
VideoWriter(F1,'rewardMovie.avi')

fig = figure;
movie(fig,F2,1)
cd('C:\Users\Wool\Desktop\tempFigs')
VideoWriter(F1,'stimMovie.avi')

% Prepare the new file.
    %   vidObj = VideoWriter('peaks.avi');
    %   open(vidObj);
    %
    %   % Create an animation.
    %   Z = peaks; surf(Z);
    %   axis tight
    %   set(gca,'nextplot','replacechildren');
    %
    %   for k = 1:20
    %      surf(sin(2*pi*k/20)*Z,Z)
    %
    %      % Write each frame to the file.
    %      currFrame = getframe;
    %      writeVideo(vidObj,currFrame);
    %   end
    % 
    %   % Close the file.
    %   close(vidObj);




