%% load experiment details

% close all;
clc;
clear all;
cd('C:\Users\Wool\Documents\MATLAB\expPipeline\expPipeline');

mouseName = 'LEW008';
expDate = '2019-04-11';
expNum = 1;
expSeries = [1];

%% load data
cd('C:\Users\Wool\Documents\GitHub\toupee\behavior')
[block, Timeline] = loadData(mouseName, expDate, expNum);

%% get event timings and wheel trajectories
cd('C:\Users\Wool\Documents\MATLAB\expPipeline\expPipeline');

signalsNames = {'stimulusOnTimes' 'stimulusOffTimes'};
[eventTimes, wheelTrajectories] = getEventTimes_passive(block, Timeline, signalsNames);

%% load traces
[allFcell, ops] = loadExpTraces(mouseName, expDate, expSeries);

%% align calcium traces to the event you want
cd('C:\Users\Wool\Documents\MATLAB\expPipeline\expPipeline');

% cut the trace into trial-by-trial traces, aligned to a particular event
event = 'stimulusOnTimes';
[alignedTraces, alignedSpikes, eventWindow] = alignSpikes(mouseName, expDate, expNum, expSeries, allFcell, eventTimes, ops, event);

%% initialize some plot values

%imaging rate
numPlanes = length(alignedTraces);
Fs = 15;%/ numPlanes;

%event window
eventIdx = find(eventWindow == 0);

%stimulus response index
crfTime = 0.2;
vPreIdx = eventIdx - 1;
vPeriIdx = eventIdx + ceil(crfTime*Fs);

responseWindow = eventIdx+3;
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

%% do PCA on the 3d matrix (cells x time x trials)

%reshape the matrix to cells x time x trials
M = permute(alignedSpikes,[3 2 1]);
[i,~] = find(isnan(M));
nanCells = unique(i);

for c = fliplr(nanCells')
    M(c,:,:) = [];
end

% center the data
Mcent = zeros(size(M,1), size(M,2),size(M,3));
for iCell = 1:size(M,1)
    Mslice = squeeze(M(iCell,:,:));
    MsliceMean = mean(Mslice,1);
    Mcent(iCell,:,:) = Mslice - MsliceMean;
end
% Mcent = M;

Mresh = reshape(Mcent,size(Mcent,1), size(Mcent,2)*size(Mcent,3));
[U,S,V] = svd(Mresh,'econ');
newM = S*abs(V');
newMresh = reshape(newM,size(M,1),size(M,2), size(M,3));

%% report the mean response of each cell to your event, per trial
%output is a matrix of size trials x cells

numPCs = 3;
responseWindow = 40;

onsetPCs = squeeze(newMresh(1:numPCs,responseWindow,:))';


%% 
cd('C:\Users\Wool\Documents\GitHub\toupee\behavior');
[condLogical, condIdxLeft] = selectCondition(block,contrasts, eventTimes,'all','all','all','left','all','all','all','all','all');
[condLogical, condIdxRight] = selectCondition(block,contrasts, eventTimes,'all','all','all','right','all','all','all','all','all');

onsetPCs_hiLeft = onsetPCs(condIdxLeft,:);
onsetPCs_Ltrain = onsetPCs_hiLeft(2:2:end,:);
onsetPCs_Ltest = onsetPCs_hiLeft(1:2:end,:);

onsetPCs_hiRight = onsetPCs(condIdxRight,:);
onsetPCs_Rtrain = onsetPCs_hiRight(2:2:end,:);
onsetPCs_Rtest = onsetPCs_hiRight(1:2:end,:);

onsetPCs_meanLeft = nanmean(onsetPCs_Ltrain,1);
onsetPCs_meanRight = nanmean(onsetPCs_Rtrain,1);

for iTrial = 1:size(onsetPCs_Ltest,1)
    LdotsL(iTrial,:) = dot(onsetPCs_Ltest(iTrial,:),onsetPCs_meanLeft);
    LdotsR(iTrial,:) = dot(onsetPCs_Ltest(iTrial,:),onsetPCs_meanRight);
end
for iTrial = 1:size(onsetPCs_Rtest,1)
    RdotsR(iTrial,:) = dot(onsetPCs_Rtest(iTrial,:),onsetPCs_meanRight);
    RdotsL(iTrial,:) = dot(onsetPCs_Rtest(iTrial,:),onsetPCs_meanLeft);
end

maxLim = 10 * round((max([max(RdotsR) max(RdotsL) max(LdotsR) max(LdotsL)]) * 1.1)/10);
minLim = 10 * round((min([min(RdotsR) min(RdotsL) min(LdotsR) min(LdotsL)]) * 1.1)/10);
stepLim = round((maxLim-minLim)/40)*10;
figure;
ax = gca;
hold on;
axis square;
plotLine = line([minLim maxLim],[minLim maxLim]);
set(plotLine,'Color', [.5 .5 .5],'LineStyle','--')
plotL = scatter(LdotsL,LdotsR,8,'b');
plotR = scatter(RdotsL,RdotsR,8,'ro');
set(plotL,'Marker','o','MarkerEdgeColor',[0 102 255]./255);
axis([minLim maxLim minLim maxLim]);
plotL.MarkerEdgeAlpha = 0.75;
hold on;
plotR.MarkerEdgeAlpha = 0.75;
ax.TickDir = 'out';
ylabel('x_i · xbar_R')
xlabel('x_i · xbar_R')
xticks([minLim:stepLim:maxLim])
yticks([minLim:stepLim:maxLim])
title(strcat(mouseName,{' '},expDate,{' '},num2str(expNum)));
