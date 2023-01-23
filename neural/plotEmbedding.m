%%
cd('C:\Users\Wool\Documents\GitHub\rastermap\matlab')
[isort1, isort2, ~] = mapTmap(neuralData.cellResps');
% [isort2, isort1, ~] = mapTmap(neuralData.cellResps);

%%

[neuralData] = getSignificantActivity(expInfo, behavioralData, neuralData,0);
hemisphere = 1;
eventWindow = neuralData.eta.eventWindow;
trialTypes = getTrialTypes(expInfo,behavioralData,'late');

%% some trial types 
nt = length(expInfo.block.events.endTrialValues);
trueStimuli = expInfo.block.events.contrastValues(1:nt);
trialCorrectChoice = expInfo.block.events.correctResponseValues(1:nt);
trueBlocks = expInfo.block.events.highRewardSideValues(1:nt);
trueChoices = expInfo.block.events.responseValues(1:nt);
trueFeedback = expInfo.block.events.feedbackValues(1:nt);

% assign the 0% stimuli as either 'left' or 'right' depending on the
% preassigned correct choice (not the mouse's choice)
trueStimuli(trueStimuli == 0) = eps;
trueStimuli(abs(trueStimuli) < .05) = ...
    trueStimuli(abs(trueStimuli) < .05).* trialCorrectChoice(abs(trueStimuli) < .05);

% low rewards are possible on sign-mismatched block and stimulus
% high rewards are possible on sign-matched block and stimulus
% 1 = high, 0 = low
trueValue(trueBlocks.*sign(trueStimuli) == -1) = -1;
trueValue(trueBlocks.*sign(trueStimuli) == 1) = 1;


%% plot sample of activity
snapshot = [110 140];
figure;
% subplot(50,1,[1 35])
imagesc(neuralData.respTimes,1:length(isort1),neuralData.cellResps(:,isort1)');
hold on

contrasts = getUniqueContrasts(expInfo);
plotColors = [.25 0 0; .5 0 0;  1 0 0;  .8 .45 .45; .75 .75 .75; .6 .8 1; 0 .4 1; 0 0 1; 0 0 .5];
for t = 1:length(expInfo.block.events.endTrialTimes)
    scol = plotColors(find(contrasts == expInfo.block.events.contrastValues(t)),:);
    line([expInfo.block.events.newTrialTimes(t) expInfo.block.events.newTrialTimes(t)],[0 length(isort1)],'Color','k','LineStyle','--')
    line([behavioralData.eventTimes(1).daqTime(t) behavioralData.eventTimes(1).daqTime(t)],[-70 0],'Color','k','LineStyle','-','LineWidth',2)
    if expInfo.block.events.feedbackValues(t) == 1
        line([behavioralData.eventTimes(5).daqTime(t) behavioralData.eventTimes(5).daqTime(t)],[-50 0],'Color',[0 .6 .1],'LineStyle','-','LineWidth',2)
    else
        line([behavioralData.eventTimes(5).daqTime(t) behavioralData.eventTimes(5).daqTime(t)],[-50 0],'Color',[.7 .1 0],'LineStyle','-','LineWidth',2)
    end
end
ylim([-50 length(isort1)])
axis xy
colormap(flipud(gray));
caxis([0 .4])
prettyPlot(gca)
set(gca,'TickDir','out','TickLength', [.00 .00])
set(gca, 'YTickLabels', {})
xlim(snapshot)
set(gca, 'XTickLabels', {'' '' ''})
axis off
figure;
subplot(5,1,1)
plot(expInfo.block.inputs.wheelTimes(1:end-1),smooth(diff(expInfo.block.inputs.wheelValues),8),'Color',[.5 0 .65],'LineWidth',1.5)
axis off
ylim([min(smoothdata(diff(expInfo.block.inputs.wheelValues))) max(smoothdata(diff(expInfo.block.inputs.wheelValues)))])
xlim(snapshot)

subplot(5,1,5)
plot(eyeData.timeInterp,eyeData.proc.pupil.area_smooth,'Color',[.9 0 0],'LineWidth',1.5)
xlim(snapshot)
axis off

subplot(5,1,4)
plot(eyeData.timeInterp,eyeData.proc.face{1, 1}.motion,'Color',[1 0.5 0],'LineWidth',1.5)
xlim(snapshot)
axis off

subplot(5,1,3)
plot(eyeData.timeInterp,eyeData.proc.face{1, 2}.motion,'Color',[.25 .75 0],'LineWidth',1.5)
xlim(snapshot)
axis off

subplot(5,1,2)
plot(eyeData.timeInterp,eyeData.proc.face{1, 3}.motion,'Color',[0 .25 .75],'LineWidth',1.5)
xlim(snapshot)
axis off
%%
whichComparison = 1;
switch whichComparison
    case 1
        plotFreqs = {'leftStim' 'rightStim'};
        plotCols = [0 .4 1; 1 0 0];
    case 2
        plotFreqs = {'leftMov' 'rightMov'};
        plotCols = [0.5 0 1; 1 0 0.5];
    case 3
        plotFreqs = {'stim' 'leftStim' 'rightStim'};
        plotCols = [0.5 .5 0.5;0 .4 1; 1 0 0];
    case 4
        plotFreqs = {'correct' 'incorrect'};
        plotCols = [0 0.5 0; .7 .1 0];
    case 5
        plotFreqs = {'mov' 'leftMov' 'rightMov'};
        plotCols = [.5 .5 .5; 0.5 0 1; 1 0 0.5];
    case 6
        plotFreqs = {'leftVal' 'rightVal'};
        plotCols = [0 .4 1; 1 0 0];
end
clear tx;
for pf = 1:length(plotFreqs)
    whichCells = plotFreqs{pf};
    plotCells = find(neuralData.stats.pValues(:,strcmp(neuralData.stats.labels,whichCells)) < 0.01);
    clear idx
    for p = 1:length(plotCells)
        idx(p) = find(isort1 == plotCells(p));
    end
    logidx = false(1,length(isort1));
    logidx(idx) = true;
    typeFrequency(:,pf) = smooth(movmean(logidx,30),30);
    [h,hh] = max(typeFrequency(:,pf));
    tx(pf,:) = [h hh];
end

contrasts = getUniqueContrasts(expInfo);    
if hemisphere > 0 
    con = contrasts;
else
    con = fliplr(contrasts);
end

fig2 = figure;
set(fig2,'position',[106 476 2408 1150])
subplot(1,length(con)+1,1)

if ~strcmp(whichCells,'all')
    hold on;
%     plot(zeros(1,length(isort1)),1:length(isort1),'k');
    for pf = 1:length(plotFreqs)
        plotSignal(1:length(isort1),zeros(1,length(isort1)),typeFrequency(:,pf)',zeros(1,length(isort1)),plotCols(pf,:),'none')
        text(length(isort1)-pf*30,.85,plotFreqs{pf},'Color',plotCols(pf,:),'HorizontalAlignment','left');
    end
%     fill([length(isort1) length(isort1) length(isort1)-200 length(isort1)-200],[.075 .025 .025 .075],'k')
%     text(length(isort1)-100,-.045,'200 cells','Color','k','HorizontalAlignment','center','rotation',90);
    set(gca,'view',[-90 90])
    xlim([0 length(isort1)])
    ylim([0 .85])
%     set(gca, 'XDir','reverse')
    axis off
end

for c = 1: length(con)
    [~, whichTrials] = selectCondition(expInfo, con(c), behavioralData, initTrialConditions('movementTime','late','responseType','all'));
    meanTrialActivity = (squeeze(nanmean(neuralData.eta.alignedResps{1}(whichTrials,:,isort1),1))');
    subplot(1,length(con)+1,c+1)
    imagesc(eventWindow,1:size(neuralData.eta.alignedResps{1},3),smoothdata(meanTrialActivity));
    colormap(flipud(gray));
    hold on;
    line([0 0],[0 size(neuralData.eta.alignedResps{1},3)],'LineStyle','--','Color','k');
    xlim([-.5 2]);
    caxis([-0 .5])
    if c == 1
%         ylabel('Cells')
        set(gca, 'YTickLabels', {})
        xlabel('Time from stimulus onset')
    else axis off
    end
    
    box off 
    set(gca,'tickdir','out')
    title(strcat(num2str(contrasts(c)*100),'%'))
    
%     if ~strcmp(whichCells,'all')
%         for p = 1:length(plotCells)
%             lh = line([-.5 2],[idx(p) idx(p)],'Color',overlay);
%             lh.Color(4) = .3;
%         end
%     end
    
    axis xy
end

%%
    figure;
    [~, allTrials] = selectCondition(expInfo, contrasts, behavioralData, initTrialConditions('movementTime','late','responseType','all'));
    meanTrialActivity = (squeeze(nanmean(neuralData.eta.alignedResps{1}(allTrials,:,isort1),1))');
    imagesc(neuralData.eta.eventWindow,1:size(neuralData.eta.alignedResps{1},3),smoothdata(meanTrialActivity));
    colormap(flipud(gray));
    hold on;
    line([0 0],[0 size(neuralData.eta.alignedResps{1},3)],'LineStyle','--','Color','k');
    line([0.8 0.8],[0 size(neuralData.eta.alignedResps{1},3)],'LineStyle','--','Color','k');
    xlim([-.5 2]);
    caxis([0.05 .25])
    xticks([-1 0 1])
    ax = gca;
    ax.YColor = 'w';
    
    box off 
    set(gca,'tickdir','out')
    
%     if ~strcmp(whichCells,'all')
%         for p = 1:length(plotCells)
%             lh = line([-.5 2],[idx(p) idx(p)],'Color',overlay);
%             lh.Color(4) = .3;
%         end
%     end
    
    axis xy
%%
ca = [-.15 .15];
gamma = .7;
et = behavioralData;
contrasts = getUniqueContrasts(expInfo);    

figure;
ax = subplot(1,6,1);
[~, allTrials] = selectCondition(expInfo, contrasts, behavioralData, initTrialConditions('movementTime','late','responseType','all'));
meanTrialActivity = (squeeze(nanmean(neuralData.eta.alignedResps{1}(allTrials,:,isort1),1))');
imagesc(neuralData.eta.eventWindow,1:size(neuralData.eta.alignedResps{1},3),smoothdata(meanTrialActivity));
hold on;
line([0 0],[0 size(neuralData.eta.alignedResps{1},3)],'LineStyle','--','Color','k');
line([0.8 0.8],[0 size(neuralData.eta.alignedResps{1},3)],'LineStyle','--','Color','k');
xlim([-.5 2]);
xticks([0 .8])
% set(gca, 'XTickLabels', {'' '' ''})
caxis([0.05 .3])
axis xy
prettyPlot(gca)
ax.YColor = 'w';

ax2 = subplot(1,6,2);
sLmL = trialTypes.intVar.all.side_direction{1,1};
sLmR = trialTypes.intVar.all.side_direction{1,2};
sRmR = trialTypes.intVar.all.side_direction{3,2};
sRmL = trialTypes.intVar.all.side_direction{3,1};
% meanTrialActivity = ...
%     ((squeeze(nanmean(neuralData.eta.alignedResps{1}(sLmL,:,isort1),1))') + (squeeze(nanmean(neuralData.eta.alignedResps{1}(sLmR,:,isort1),1))')) - ...
%     ((squeeze(nanmean(neuralData.eta.alignedResps{1}(sRmR,:,isort1),1))') + (squeeze(nanmean(neuralData.eta.alignedResps{1}(sRmL,:,isort1),1))'));
meanTrialActivity = ...
    ((squeeze(nanmean(neuralData.eta.alignedResps{1}([sLmL sLmR],:,isort1),1))')) - ...
    ((squeeze(nanmean(neuralData.eta.alignedResps{1}([sRmR sRmL],:,isort1),1))'));

imagesc(neuralData.eta.eventWindow,1:size(neuralData.eta.alignedResps{1},3),smoothdata(meanTrialActivity));
hold on;
line([0 0],[0 size(neuralData.eta.alignedResps{1},3)],'LineStyle','--','Color','k');
line([0.8 0.8],[0 size(neuralData.eta.alignedResps{1},3)],'LineStyle','--','Color','k');
axis xy
caxis(ca)
prettyPlot(gca)
cm2 = colormapThruWhite([1 0 0],[0 .0 1],100,gamma);
ax2.YColor = 'w';
xlim([-.5 2]);
xticks([0 .8])
% set(gca, 'XTickLabels', {'' '' ''})

ax5 = subplot(1,6,3);
leftTrials = cat(2,trialTypes.intVar.cb3D.direction_block{1,1},trialTypes.intVar.cb3D.direction_block{2,1});
rightTrials = cat(2,trialTypes.intVar.cb3D.direction_block{1,2},trialTypes.intVar.cb3D.direction_block{2,2});
meanTrialActivity = (squeeze(nanmean(neuralData.eta.alignedResps{1}(leftTrials,:,isort1),1))') - (squeeze(nanmean(neuralData.eta.alignedResps{1}(rightTrials,:,isort1),1))');
imagesc(neuralData.eta.eventWindow,1:size(neuralData.eta.alignedResps{1},3),smoothdata(meanTrialActivity));
hold on;
line([0 0],[0 size(neuralData.eta.alignedResps{1},3)],'LineStyle','--','Color','k');
line([0.8 0.8],[0 size(neuralData.eta.alignedResps{1},3)],'LineStyle','--','Color','k');
axis xy
caxis(ca)
prettyPlot(gca) %0.1 .7 .1; 1 .6 0
cm5 = colormapThruWhite([1 .6 0],[.1 .7 .1],100,gamma);
ax5.YColor = 'w';
xlim([-.5 2]);
xticks([0 .8])
% set(gca, 'XTickLabels', {'' '' ''})

ax6 = subplot(1,6,4);
[~, hiVal] = find(trueValue == 1);
[~, lowVal] = find(trueValue == 0);
meanTrialActivity = (squeeze(nanmean(neuralData.eta.alignedResps{1}(hiVal,:,isort1),1))') - (squeeze(nanmean(neuralData.eta.alignedResps{1}(lowVal,:,isort1),1))');
imagesc(neuralData.eta.eventWindow,1:size(neuralData.eta.alignedResps{1},3),smoothdata(meanTrialActivity));
hold on;
line([0 0],[0 size(neuralData.eta.alignedResps{1},3)],'LineStyle','--','Color','k');
line([0.8 0.8],[0 size(neuralData.eta.alignedResps{1},3)],'LineStyle','--','Color','k');
axis xy
caxis(ca)
prettyPlot(gca) %0.1 .7 .1; 1 .6 0
cm6 = colormapThruWhite([1 0 1],[0 .75 0],100,gamma);
ax6.YColor = 'w';
xlim([-.5 2]);
xticks([0 .8])
% set(gca, 'XTickLabels', {'' '' ''})

ax3 = subplot(1,6,5);
[~, leftStimTrials] = selectCondition(expInfo, contrasts, et, initTrialConditions('movementTime','late','responseType','all','movementDir','cw'));
[~, rightStimTrials] = selectCondition(expInfo, contrasts, et, initTrialConditions('movementTime','late','responseType','all','movementDir','ccw'));
meanTrialActivity = (squeeze(nanmean(neuralData.eta.alignedResps{1}(leftStimTrials,:,isort1),1))') - (squeeze(nanmean(neuralData.eta.alignedResps{1}(rightStimTrials,:,isort1),1))');
imagesc(neuralData.eta.eventWindow,1:size(neuralData.eta.alignedResps{1},3),smoothdata(meanTrialActivity));
hold on;
line([0 0],[0 size(neuralData.eta.alignedResps{1},3)],'LineStyle','--','Color','k');
line([0.8 0.8],[0 size(neuralData.eta.alignedResps{1},3)],'LineStyle','--','Color','k');
axis xy
caxis(ca)
prettyPlot(gca)
cm3 = colormapThruWhite([1 0 .5],[.5 0 1],100,gamma);
ax3.YColor = 'w';
xlim([-.5 2]);
xticks([0 .8])
% set(gca, 'XTickLabels', {'' '' ''})

ax4 = subplot(1,6,6);
[~, leftStimTrials] = selectCondition(expInfo, contrasts, et, initTrialConditions('movementTime','late','responseType','correct','movementDir','all'));
[~, rightStimTrials] = selectCondition(expInfo, contrasts, et, initTrialConditions('movementTime','late','responseType','incorrect','movementDir','all'));
meanTrialActivity = (squeeze(nanmean(neuralData.eta.alignedResps{1}(leftStimTrials,:,isort1),1))') - (squeeze(nanmean(neuralData.eta.alignedResps{1}(rightStimTrials,:,isort1),1))');
imagesc(neuralData.eta.eventWindow,1:size(neuralData.eta.alignedResps{1},3),smoothdata(meanTrialActivity));
hold on;
line([0 0],[0 size(neuralData.eta.alignedResps{1},3)],'LineStyle','--','Color','k');
line([0.8 0.8],[0 size(neuralData.eta.alignedResps{1},3)],'LineStyle','--','Color','k');
axis xy
caxis(ca)
prettyPlot(gca)
cm4 = colormapThruWhite([.7 0 .1],[0 0.5 0],100,gamma);
ax4.YColor = 'w';
xlim([-.5 2]);
xticks([0 .8])
% set(gca, 'XTickLabels', {'' '' ''})

colormap(ax,flipud(gray));
colormap(ax2,cm2);
colormap(ax3,cm3);
colormap(ax4,cm4);
colormap(ax5,cm5);
colormap(ax6,cm6);

%%

maxVels = abs(behavioralData.wheelMoves.epochs(5).peakVel);
maxVels(isnan(maxVels)) = 0;
RTs = behavioralData.wheelMoves.epochs(5).onsetTimes(1:nt) - behavioralData.eventTimes(1).daqTime(1:nt);
RTs(isnan(RTs)) = 0;

contrasts = getUniqueContrasts(expInfo);
[~, whichTrials] = selectCondition(expInfo, contrasts, behavioralData, ...
    initTrialConditions('movementTime','late','specificRTs',[.1 3]));

for t = 1:size(neuralData.eta.alignedResps{1},2)
    array = squeeze(neuralData.eta.alignedResps{1}(whichTrials,t,:));
    velCorr(:,t) = corr(array,maxVels(whichTrials)');
    RTCorr(:,t) = corr(array,RTs(whichTrials)');
    stimCorr(:,t) = corr(array,trueStimuli(whichTrials)');
    choiceCorr(:,t) = corr(array,trueChoices(whichTrials)');
    feedbackCorr(:,t) = corr(array,trueFeedback(whichTrials)');
    valCorr(:,t) = corr(array,trueValue(whichTrials)');
    blockCorr(:,t) = corr(array,trueBlocks(whichTrials)');
end
    
ca = [-.2 .2];
xl = [-.5 2];
cm2 = colormapThruWhite([1 0 0],[0 .0 1],100,gamma);
cm5 = colormapThruWhite([1 0 .5],[.5 0 1],100,gamma);
cm6 = colormapThruWhite([.7 0 .1],[0 0.5 0],100,gamma);
cm4 = colormapThruWhite([1 .6 0],[.1 .7 .1],100,gamma);
cm3 = colormapThruWhite([1 0 1],[0 .75 0],100,gamma);
cm7 = colormapThruWhite([1 .6 0],[.1 .7 .1],100,gamma);
cm8 = colormapThruWhite([1 0 1],[0 .75 0],100,gamma);

figure;
ax = subplot(1,8,1);
[~, allTrials] = selectCondition(expInfo, contrasts, behavioralData, initTrialConditions('movementTime','late','responseType','all'));
meanTrialActivity = (squeeze(nanmean(neuralData.eta.alignedResps{1}(allTrials,:,isort1),1))');
imagesc(neuralData.eta.eventWindow,1:size(neuralData.eta.alignedResps{1},3),smoothdata(meanTrialActivity));
hold on;
line([0 0],[0 size(neuralData.eta.alignedResps{1},3)],'LineStyle','--','Color','k');
line([0.8 0.8],[0 size(neuralData.eta.alignedResps{1},3)],'LineStyle','--','Color','k');
xlim(xl);
xticks([0 .8])
% set(gca, 'XTickLabels', {'' '' ''})
caxis([0.05 .3])
axis xy
prettyPlot(gca)
ax.YColor = 'w';

ax2 = subplot(1,8,2);
imagesc(neuralData.eta.eventWindow,1:length(isort1),smoothdata(stimCorr(isort1,:),2,'gaussian',.2))
colormap(BlueWhiteRed)
hold on;
line([0 0],[0 size(neuralData.eta.alignedResps{1},3)],'LineStyle','--','Color','k');
line([0.8 0.8],[0 size(neuralData.eta.alignedResps{1},3)],'LineStyle','--','Color','k');
xlim([-.5 2]);
xticks([0 .8])
axis xy
prettyPlot(gca)
caxis(ca)
xlim(xl)
ax2.YColor = 'w';

ax3 = subplot(1,8,3);
imagesc(neuralData.eta.eventWindow,1:length(isort1),smoothdata(valCorr(isort1,:),2,'gaussian',.2))
colormap(BlueWhiteRed)
hold on;
line([0 0],[0 size(neuralData.eta.alignedResps{1},3)],'LineStyle','--','Color','k');
line([0.8 0.8],[0 size(neuralData.eta.alignedResps{1},3)],'LineStyle','--','Color','k');
xlim([-.5 2]);
xticks([0 .8])
axis xy
prettyPlot(gca)
caxis(ca)
xlim(xl)
ax3.YColor = 'w';

ax4 = subplot(1,8,4);
imagesc(neuralData.eta.eventWindow,1:length(isort1),smoothdata(blockCorr(isort1,:),2,'gaussian',.2))
colormap(BlueWhiteRed)
hold on;
line([0 0],[0 size(neuralData.eta.alignedResps{1},3)],'LineStyle','--','Color','k');
line([0.8 0.8],[0 size(neuralData.eta.alignedResps{1},3)],'LineStyle','--','Color','k');
xlim([-.5 2]);
xticks([0 .8])
axis xy
prettyPlot(gca)
caxis(ca)
xlim(xl)
ax4.YColor = 'w';

ax5 = subplot(1,8,5);
imagesc(neuralData.eta.eventWindow,1:length(isort1),smoothdata(choiceCorr(isort1,:),2,'gaussian',.2))
colormap(BlueWhiteRed)
hold on;
line([0 0],[0 size(neuralData.eta.alignedResps{1},3)],'LineStyle','--','Color','k');
line([0.8 0.8],[0 size(neuralData.eta.alignedResps{1},3)],'LineStyle','--','Color','k');
xlim([-.5 2]);
xticks([0 .8])
axis xy
prettyPlot(gca)
caxis(ca)
xlim(xl)
ax5.YColor = 'w';

ax6 = subplot(1,8,6);
imagesc(neuralData.eta.eventWindow,1:length(isort1),smoothdata(feedbackCorr(isort1,:),2,'gaussian',.2))
colormap(BlueWhiteRed)
hold on;
line([0 0],[0 size(neuralData.eta.alignedResps{1},3)],'LineStyle','--','Color','k');
line([0.8 0.8],[0 size(neuralData.eta.alignedResps{1},3)],'LineStyle','--','Color','k');
xlim([-.5 2]);
xticks([0 .8])
axis xy
prettyPlot(gca)
caxis(ca)
xlim(xl)
ax6.YColor = 'w';

ax7 = subplot(1,8,7);
imagesc(neuralData.eta.eventWindow,1:length(isort1),smoothdata(velCorr(isort1,:),2,'gaussian',.2))
colormap(BlueWhiteRed)
hold on;
line([0 0],[0 size(neuralData.eta.alignedResps{1},3)],'LineStyle','--','Color','k');
line([0.8 0.8],[0 size(neuralData.eta.alignedResps{1},3)],'LineStyle','--','Color','k');
xlim([-.5 2]);
xticks([0 .8])
axis xy
prettyPlot(gca)
caxis(ca)
xlim(xl)
ax7.YColor = 'w';

ax8 = subplot(1,8,8);
imagesc(neuralData.eta.eventWindow,1:length(isort1),smoothdata(RTCorr(isort1,:),2,'gaussian',.2))
colormap(BlueWhiteRed)
hold on;
line([0 0],[0 size(neuralData.eta.alignedResps{1},3)],'LineStyle','--','Color','k');
line([0.8 0.8],[0 size(neuralData.eta.alignedResps{1},3)],'LineStyle','--','Color','k');
xlim([-.5 2]);
xticks([0 .8])
axis xy
prettyPlot(gca)
caxis(ca)
xlim(xl)
ax8.YColor = 'w';

colormap(ax,flipud(gray));
colormap(ax2,cm2);
colormap(ax3,cm3);
colormap(ax4,cm4);
colormap(ax5,cm5);
colormap(ax6,cm6);
colormap(ax7,cm7);
colormap(ax8,cm8);

    
%% WORKBENCH

% whichCells = 'rightStim'; %choose from 'pLabels' array
% if strcmp(whichCells, 'all')
%     plotCells = (1:size(neuralData.eta.alignedResps{1},3))';
% else
%     plotCells = find(neuralData.stats.pValues(:,strcmp(neuralData.stats.labels,whichCells)) < 0.01);
% end
% 
% clear idx
% for p = 1:length(plotCells)
%     idx(p) = find(isort1 == plotCells(p));
% end

% if strcmp(whichCells,'leftStim') || strcmp(whichCells,'rightStim')
%     overlay = [0 1 .2];
% elseif strcmp(whichCells,'leftMov')
%     if hemisphere > 0
%         overlay = [0 .4 1];
%     else
%         overlay = [1 0 0];
%     end
% elseif strcmp(whichCells,'rightMov')
%     if hemisphere < 0
%         overlay = [0 .4 1];
%     else
%         overlay = [1 0 0];
%     end
% elseif strcmp(whichCells,'value')
%     overlay = [1 0 1];
% elseif strcmp(whichCells,'hit')
%     overlay = [0 .7 0];
% elseif strcmp(whichCells,'miss')
%     overlay = [.5 0 0];
% elseif strcmp(whichCells,'impulsive')
%     overlay = [1 .5 0];
% elseif strcmp(whichCells,'patient')
%     overlay = [0 .5 .5];
% end
