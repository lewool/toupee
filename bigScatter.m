%% initialize figs
balanced = 1;
axLimA = [0 .25 0 .25];
axLimB = [0 .25 0 .25];
mSize = 20;
color = [0 0 0];

scatterFig = figure(1);
set(gcf,'position', [98         254        1366         657]);
for s = 1:2
    subplot(2,4,s)
    l = line([axLimA(1) axLimA(2)],[axLimA(1) axLimA(2)]);
    set(l,...
        'LineStyle', '--',...
        'LineWidth',1,...
        'Color',[.5 .5 .5]...
        );
    hold on
end
for s = 1:4
    subplot(2,4,s+4)
    l = line([axLimA(1) axLimA(2)],[axLimA(1) axLimA(2)]);
    set(l,...
        'LineStyle', '--',...
        'LineWidth',1,...
        'Color',[.5 .5 .5]...
        );
    hold on
end

hold on;

timelineFig = figure(2);
hold on;
%%
for ex = 1:9
%% import data

nd = neuralData(ex);
[baselineResps, stimResps, pmovResps, movResps, rewResps] = getEpochResps(nd.eta);
alignedResps = nd.eta.alignedResps;
bfcH = nd.stats.bfcH;
pLabels = nd.stats.labels;
eventWindow = nd.eta.eventWindow;

%% set up trial conditions to compare

clear contrastConditions trialConditions labels condIdx
contrasts = getUniqueContrasts(expInfo(ex));
allContrasts = getAllContrasts(expInfo(ex));

%set up trial conditions for hi-L and hi-R blocks
trialConditions{1} = initTrialConditions('highRewardSide','left','movementDir','cw','movementTime','late');
trialConditions{2} = initTrialConditions('highRewardSide','left','movementDir','ccw','movementTime','late');
trialConditions{3} = initTrialConditions('highRewardSide','right','movementDir','cw','movementTime','late');
trialConditions{4} = initTrialConditions('highRewardSide','right','movementDir','ccw','movementTime','late');
trialConditions{5} = initTrialConditions('highRewardSide','left','movementTime','late');
trialConditions{6} = initTrialConditions('highRewardSide','right','movementTime','late');
trialConditions{7} = initTrialConditions('movementDir','cw','movementTime','late');
trialConditions{8} = initTrialConditions('movementDir','ccw','movementTime','late');

trialLabels{1} = 'bL_mL_';
trialLabels{2} = 'bL_mR_';
trialLabels{3} = 'bR_mL_';
trialLabels{4} = 'bR_mR_';
trialLabels{5} = 'bL_mAll_';
trialLabels{6} = 'bR_mAll_';
trialLabels{7} = 'bAll_mL_';
trialLabels{8} = 'bAll_mR_';

contrastConditions{1} = contrasts(contrasts<0);
contrastConditions{2} = contrasts(contrasts>0);
contrastConditions{3} = contrasts(contrasts~=0);
contrastLabels{1} = 'sL';
contrastLabels{2} = 'sR';
contrastLabels{3} = 'sAll';

testTrials = 1:2:size(alignedResps{1},1);
trainTrials = 2:2:size(alignedResps{1},1);

d = 1;
for c = 1:length(contrastConditions)
    for t = 1:length(trialConditions)
        [~, condIdx{d,:}.all] = selectCondition(expInfo(ex), contrastConditions{c}, behavioralData(ex), trialConditions{t});
        condIdx{d,:}.test = intersect(testTrials,condIdx{d}.all);
        condIdx{d,:}.train = intersect(trainTrials,condIdx{d}.all);
        labels{d,1} = strcat(trialLabels{t},contrastLabels{c});
        d = d+1;
    end
end

%% equalize contrasts between trial types
clear whichMatchTrials;

matchList = {...
    'bAll_mL_sL' , 'bAll_mR_sL';...
    'bAll_mL_sR' , 'bAll_mR_sR';...
    'bL_mL_sL' , 'bR_mR_sL';...
    'bL_mL_sR' , 'bR_mR_sR';...
    'bL_mL_sL' , 'bR_mR_sL';...
    'bL_mL_sR' , 'bR_mR_sR';...
    };

for iM = 1:size(matchList,1)
    mLTrials = condIdx{strcmp(labels,matchList{iM,1})}.all; 
    mRTrials = condIdx{strcmp(labels,matchList{iM,2})}.all; 
    mLContrasts = allContrasts(mLTrials);
    mRContrasts = allContrasts(mRTrials);
    
    subsetLTrials = [];
    subsetRTrials = [];
    
    for c = 1:length(contrasts)
        nmLT = sum(mLContrasts == contrasts(c));
        nmRT = sum(mRContrasts == contrasts(c));
        minShared = min([nmLT nmRT]);
        subsetLTrials = [subsetLTrials randsample(mLTrials(mLContrasts == contrasts(c)),minShared)];
        subsetRTrials = [subsetRTrials randsample(mRTrials(mRContrasts == contrasts(c)),minShared)];
    end
    
    whichMatchTrials{iM,1} = subsetLTrials;
    whichMatchTrials{iM,2} = subsetRTrials;
end

%%
% close all
colors = [0.1 0.7 0.1; 1 .6 0; 0 .4 1; 1 0 0];

respPeriod = stimResps;

if ex < 9
    whichCells = 'leftStim'; %choose from 'pLabels' array
    if strcmp(whichCells, 'all')
        plotCells = 1:size(alignedResps{1},3);
    else
        plotCells = find(bfcH(:,strcmp(pLabels,whichCells)) > 0);
    end
else
        whichCells = 'rightStim'; %choose from 'pLabels' array
    if strcmp(whichCells, 'all')
        plotCells = 1:size(alignedResps{1},3);
    else
        plotCells = find(bfcH(:,strcmp(pLabels,whichCells)) > 0);
    end
end

%% stim versus move
condVector = {'bAll_mL_sL' 'bAll_mR_sL' 'bAll_mL_sR' 'bAll_mR_sR'};

clear X;
for iCond = 1:length(condVector)
    if balanced && sum(sum(strcmp(condVector(iCond),matchList)))
        whichTrials = whichMatchTrials{strcmp(condVector(iCond),matchList)};
    else
        whichTrials = condIdx{strcmp(labels,condVector{iCond})}.all;
    end
    numTrials = size(whichTrials,2);
    X(:,iCond) = mean(respPeriod(whichTrials,plotCells),1);
end

%% high versus low
condVectorY = {'bR_mL_sL' 'bR_mR_sL' 'bR_mL_sR' 'bR_mR_sR'};
condVectorZ = {'bL_mL_sL' 'bL_mR_sL' 'bL_mL_sR' 'bL_mR_sR'};

clear Z Y;
for iCond = 1:length(condVectorZ)
    if balanced && sum(sum(strcmp(condVectorZ(iCond),matchList)))
        whichTrials = whichMatchTrials{strcmp(condVectorZ(iCond),matchList)};
    else
        whichTrials = condIdx{strcmp(labels,condVectorZ{iCond})}.all;
    end
    numTrials = size(whichTrials,2);
    Z(:,iCond) = mean(respPeriod(whichTrials,plotCells),1);
end

for iCond = 1:length(condVectorY)
    if balanced && sum(sum(strcmp(condVectorY(iCond),matchList)))
        whichTrials = whichMatchTrials{strcmp(condVectorY(iCond),matchList)};
    else
        whichTrials = condIdx{strcmp(labels,condVectorY{iCond})}.all;
    end
    numTrials = size(whichTrials,2);
    Y(:,iCond) = mean(respPeriod(whichTrials,plotCells),1);
end

%% plot all

if ex < 9
    xCC = X(:,1);
    xCI = X(:,2);
    xII = X(:,4);
    xIC = X(:,3);
else
    xCC = X(:,4);
    xCI = X(:,3);
    xII = X(:,1);
    xIC = X(:,2);
end
%%


figure(1);
ax1 = subplot(2,4,1);
hold on;

p = scatter(xCC,xCI,mSize);
xlabel('chose contra')
ylabel('chose ipsi')
title('stim contra')
set(p,...
    'MarkerEdgeColor','none',...
    'Marker','o',...
    'MarkerFaceColor',[0 .4 1],...
    'MarkerFaceAlpha',.3...
    );
axis(axLimA);
axis square

ax2 = subplot(2,4,2);
hold on;
p = scatter(xIC,xII,mSize);
xlabel('chose contra')
ylabel('chose ipsi')
title('stim ipsi')
set(p,...
    'MarkerEdgeColor','none',...
    'Marker','o',...
    'MarkerFaceColor',[1 0 0],...
    'MarkerFaceAlpha',.3...
    );
axis(axLimA);
axis square

for s = 1:4
    ax = subplot(2,4,s+4);
    hold on;
    x = Z(:,s);
    y = Y(:,s);
    p = scatter(x,y,mSize);
    title(condVectorZ{s}(end-4:end),'Interpreter','none');  
    xlabel('high reward contra')
    ylabel('high reward ipsi')
    if s < 3
    set(p,...
        'MarkerEdgeColor','none',...
        'Marker','o',...
        'MarkerFaceColor',[0 .4 1],...
        'MarkerFaceAlpha',.3...
        );
    elseif s >= 3
        set(p,...
        'MarkerEdgeColor','none',...
        'Marker','o',...
        'MarkerFaceColor',[1 0 0],...
        'MarkerFaceAlpha',.3...
        );
    end
    axis(axLimB);
    axis square
end

%% 
condVector = {'bAll_mL_sL' 'bAll_mR_sL' 'bAll_mL_sR' 'bAll_mR_sR'};
color = [0 0 0];

clear XX;
for iCond = 1:length(condVector)
    if balanced && sum(sum(strcmp(condVector(iCond),matchList)))
        whichTrials = whichMatchTrials{strcmp(condVector(iCond),matchList)};
    else
        whichTrials = condIdx{strcmp(labels,condVector{iCond})}.all;
    end
    numTrials = size(whichTrials,2);
    XX(:,:,iCond) = squeeze(mean(alignedResps{1}(whichTrials,:,plotCells),1));
end

if ex > 9
    meanPrefIdx_C = mean(XX(:,:,1) - XX(:,:,2),2);
    meanPrefIdx_I = mean(XX(:,:,3) - XX(:,:,4),2);
    stdPrefIdx = std(XX(:,:,1) - XX(:,:,2),[],2)/sqrt(size(XX,2));
else
   meanPrefIdx_C = mean(XX(:,:,4) - XX(:,:,3),2);
   meanPrefIdx_I = mean(XX(:,:,2) - XX(:,:,1),2);
stdPrefIdx = std(XX(:,:,4) - XX(:,:,3),[],2)/sqrt(size(XX,2)); 
end

windowSize = 3;
b = (1/windowSize)*ones(1,windowSize);
a = 1;
meanC_smoothed(:,ex) = filter(b, a, meanPrefIdx_C);
meanI_smoothed(:,ex) = filter(b, a, meanPrefIdx_I);
std_smoothed(:,ex) = filter(b, a, stdPrefIdx);
%%
figure(2);
hold on;
subplot(1,2,1)
hold on;
plotLine = plot(eventWindow,meanC_smoothed(:,ex));
set(plotLine,'LineWidth',2,'Color',[0 .4 1]);
% plotCI = fill([eventWindow';flipud(eventWindow')],[(LPI_smoothed-std_smoothed);flipud((LPI_smoothed+std_smoothed))],[0 0 0], 'LineStyle', 'none');
plotLine.Color(4) = (0.15);
xlim([-1 2])
ylabel('Contra choice – Ipsi choice')
xlabel('Time from stimulus onset (s)')
set(gca,'tickdir','out')

subplot(1,2,2)
hold on;
plotLine = plot(eventWindow,meanI_smoothed(:,ex));
set(plotLine,'LineWidth',2,'Color',[1 0 0]);
% plotCI = fill([eventWindow';flipud(eventWindow')],[(LPI_smoothed-std_smoothed);flipud((LPI_smoothed+std_smoothed))],[0 0 0], 'LineStyle', 'none');
plotLine.Color(4) = (0.15);
xlim([-1 2])
ylabel('Contra choice – Ipsi choice')
xlabel('Time from stimulus onset (s)')
set(gca,'tickdir','out')

end

figure(2);
set(gcf,'position',[816   606   781   269])

hold on;
subplot(1,2,1)
plotLine = plot(eventWindow,mean(meanC_smoothed,2));
hold on;
set(plotLine,'LineWidth',3,'Color',[0 .4 1]);
subplot(1,2,2)
hold on;
plotLine = plot(eventWindow,mean(meanI_smoothed,2));
set(plotLine,'LineWidth',3,'Color',[1 0 0]);


%% WORKSHOP

% allX = [X(:,:,1) X(:,:,2)];
% [M,N] = size(allX);
% % Preserve the row indices
% rowIndex = repmat((1:M)',[1 N]);
% for i = 1:38
% % Get randomized column indices by sorting a second random array
% [~,randomizedColIndex] = sort(rand(M,N),2);
% % Need to use linear indexing to create B
% newLinearIndex = sub2ind([M,N],rowIndex,randomizedColIndex);
% shuffleX(:,:,i) = allX(newLinearIndex);
% end
% diffShuffleX = squeeze(mean(shuffleX(:,1:19,:) - shuffleX(:,20:38,:),2));
% sortedDiffShuffle = sort(diffShuffleX,2);
% meanShuffle = mean(diffShuffleX,2);
% stdShuffle = std(diffShuffleX,[],2)/sqrt(size(diffShuffleX,2));
% 
% plotMS = plot(eventWindow, meanShuffle);
% set(plotMS,'LineWidth',2,'Color','k');
% plotSCI = fill([eventWindow';flipud(eventWindow')],[(meanShuffle-stdShuffle);flipud((meanShuffle+stdShuffle))],'k', 'LineStyle', 'none');
% alpha(0.2);

% figMatrix = figure;
% set(figMatrix,'position',[40 90 1230 900]);
% 
% for r = 1:4
%     x1 = X(:,r);
%     for c = 1:4
%         x2 = X(:,c);
%         subplot(4,4,sub2ind([4 4],r,c));
%         hold on;
%         if r == c
% %             cla;
%             h = histogram(x1,linspace(axLimA(1),axLimA(2),11),'Normalization','probability');
%             set(h,'FaceColor',color,...
%                 'FaceAlpha',.2...
%                 );
%             set(gca,'Color','none');
%             box off;
%             axis square;
%             xlim([axLimA(1) axLimA(2)]);
%         else
%             p = scatter(x1,x2,mSize);
%             l = line([axLimA(1) axLimA(2)],[axLimA(1) axLimA(2)]);
%             set(p,...
%                 'MarkerEdgeColor','none',...
%                 'Marker','o',...
%                 'MarkerFaceColor',color,...
%                 'MarkerFaceAlpha',.2...
%                 );
%             set(l,...
%                 'LineStyle', '--',...
%                 'LineWidth',1,...
%                 'Color',[.5 .5 .5]...
%                 );
%             axis(axLimA);
%         end
%         axis square
%         
%         if c == 4
%             xlabel(condVector{r},'Interpreter','none')
%         end
%         
%         if r == 1
%             ylabel(condVector{c},'Interpreter','none')
%         end
%         
%     end
% end

% figVector = figure;
% set(figVector,'position',[1280 90 260 900]);
% for s = 1:4
%     subplot(4,1,s)
%     x = Z(:,s);
%     y = Y(:,s);
%     p = scatter(x,y,mSize);
%         l = line([axLimB(1) axLimB(2)],[axLimB(1) axLimB(2)]);
%     title(condVectorZ{s}(end-4:end),'Interpreter','none');  
%     xlabel('bL')
%     ylabel('bR')
%     set(p,...
%         'MarkerEdgeColor','none',...
%         'Marker','o',...
%         'MarkerFaceColor',color,...
%         'MarkerFaceAlpha',.2...
%         );
%     set(l,...
%         'LineStyle', '--',...
%         'LineWidth',1,...
%         'Color',[.5 .5 .5]...
%         );
%     axis(axLimB);
%     axis square
% end
