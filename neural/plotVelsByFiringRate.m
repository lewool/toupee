for m = 1:length(mouseList)
    
    mouseName = char(mouseList{m});
    expDate = char(expList{m}{1});
    expNum = expList{m}{2};
    hemisphere = hemList(m);
    [expRef, ~] = data.constructExpRef(mouseName,expDate,expNum);
    fprintf('Loading %s...\n',expRef)
    [behavioralData, expInfo, neuralData] = data.loadDataset(mouseName, expDate, expNum);
    [baselineResps, stimResps, pmovResps, movResps, rewResps, preCueResps] = getEpochResps(neuralData.eta);
%%

    et = behavioralData;
    contrasts = getUniqueContrasts(expInfo);
    nt = length(et.eventTimes(1).daqTime);

    trialCorrectChoice = expInfo.block.events.correctResponseValues(1:nt);
    trueStimuli = expInfo.block.events.contrastValues(1:nt);
    sidedStimuli = expInfo.block.events.contrastValues(1:nt);
    sidedStimuli(sidedStimuli == 0) = eps;
    sidedStimuli(abs(sidedStimuli) < .05) = ...
        sidedStimuli(abs(sidedStimuli) < .05).* ...
        trialCorrectChoice(abs(sidedStimuli) < .05);
    trueChoices = et.wheelMoves.epochs(5).moveDir;
    trueBlock = expInfo.block.events.highRewardSideValues(1:nt);
    trueFeedback = zeros(1,nt);
    trueFeedback(trueChoices .* trialCorrectChoice > 0) = 1;
    trueFeedback(trueChoices .* trialCorrectChoice < 0) = 0;
    trueValue(sidedStimuli .* trueBlock > 0) = 1;
    trueValue(sidedStimuli .* trueBlock < 0) = 0;
    
    [impTrials, ~] = selectCondition(expInfo, contrasts, et, ...
        initTrialConditions('movementTime','early'));

    maxVels = behavioralData.wheelMoves.epochs(5).peakVel;
    RTs = behavioralData.wheelMoves.epochs(5).onsetTimes - behavioralData.eventTimes(1).daqTime;
    trials = intersect(...
        find((RTs < prctile(RTs,97.5)) & (RTs > prctile(RTs,2.5))),...
        find(~isnan(trueChoices)));
    
    fitData{m,1} = impTrials(trials)';
    fitData{m,2} = maxVels(trials)';
    fitData{m,3} = RTs(trials)';
    fitData{m,4} = sidedStimuli(trials)';
    fitData{m,5} = trueChoices(trials)';
    fitData{m,6} = trueBlock(trials)';
    fitData{m,7} = trueFeedback(trials)';
    fitData{m,8} = trueValue(trials)';

    clearvars -except mouseList expList hemList fitData
    
end

%%


behavioralData = getTrialVelocities(expInfo, behavioralData);
for t = 1:length(neuralData.eta.eventWindow)
    [~,dsew_idx(t)] = min(abs(behavioralData.eta.eventWindow - neuralData.eta.eventWindow(t)));
end

%%
colors = [0 0 .5; 0 0 1; 0 .4 1; .6 .8 1; .75 .75 .75; .75 .75 .75; .8 .45 .45; 1 0 0; .5 0 0; .25 0 0];

epsContrasts = unique(fitData{m,4});

whichCells = find(neuralData.stats.hTest(:,strcmp(neuralData.stats.labels,'choice')) > 0);

ETA = 2;
% figure;
selectTs = [15 17 19 21 23 25];
for t = 1:length(selectTs)
    for c = 1:length(epsContrasts)
        subplot(length(selectTs),length(epsContrasts),(t-1)*length(epsContrasts) + c)
        cla
        hiTrials = trials((fitData{m,4} == epsContrasts(c)) & (fitData{m,8} == 1) & (fitData{m,5} == -1));
        loTrials = trials((fitData{m,4} == epsContrasts(c)) & (fitData{m,8} == 0) & (fitData{m,5} == -1));
        scatter(behavioralData.eta.alignedVels{ETA}(hiTrials,dsew_idx(selectTs(t))), ...
            nanmean(neuralData.eta.alignedResps{ETA}(hiTrials,selectTs(t),whichCells),3),...
            25,'MarkerFaceColor',colors(c,:),'MarkerEdgeColor','none','MarkerFaceAlpha',.7);
        correlation(t,c,1) = corr(...
            behavioralData.eta.alignedVels{ETA}(hiTrials,dsew_idx(selectTs(t))),...
            nanmean(neuralData.eta.alignedResps{ETA}(hiTrials,selectTs(t),whichCells),3));
        hold on;
        scatter(behavioralData.eta.alignedVels{ETA}(loTrials,dsew_idx(selectTs(t))), ...
            nanmean(neuralData.eta.alignedResps{ETA}(loTrials,selectTs(t),whichCells),3),...
            25,'MarkerFaceColor','none','MarkerEdgeColor',colors(c,:),'MarkerEdgeAlpha',.7);
        correlation(t,c,2) = corr(...
            behavioralData.eta.alignedVels{ETA}(loTrials,dsew_idx(selectTs(t))),...
            nanmean(neuralData.eta.alignedResps{ETA}(loTrials,selectTs(t),whichCells),3));
        xlim([-200 200])
        ylim([.0 .5])
        line([0 0],[0 .5],'Color',[.5 .5 .5],'linestyle',':')
    end
end

%%
preMoveIdx = 18:19;
postMoveIdx = 21:22;
figure;


for c = 1:length(epsContrasts)
    subplot(1,length(epsContrasts),c)
    cla
    hiTrials = trials((fitData{m,4} == epsContrasts(c)) & (fitData{m,8} == 1) & (fitData{m,5} == 1));
    loTrials = trials((fitData{m,4} == epsContrasts(c)) & (fitData{m,8} == 0) & (fitData{m,5} == 1));
    scatter(nanmean(behavioralData.eta.alignedVels{ETA}(hiTrials,dsew_idx(postMoveIdx)),2), ...
        nanmean(nanmean(neuralData.eta.alignedResps{ETA}(hiTrials,postMoveIdx,whichCells),3),2),...
        25,'MarkerFaceColor',colors(c,:),'MarkerEdgeColor','none','MarkerFaceAlpha',.7);
    hold on;
    scatter(nanmean(behavioralData.eta.alignedVels{ETA}(loTrials,dsew_idx(postMoveIdx)),2), ...
        nanmean(nanmean(neuralData.eta.alignedResps{ETA}(loTrials,postMoveIdx,whichCells),3),2),...
        25,'MarkerFaceColor','none','MarkerEdgeColor',colors(c,:),'MarkerEdgeAlpha',.7);
    xlim([-200 200])
    ylim([.0 .5])
    line([0 0],[0 .5],'Color',[.5 .5 .5],'linestyle',':')
end

%%
whichCells = find(neuralData.stats.hTest(:,strcmp(neuralData.stats.labels,'choice')) == 1);
hiTrials = trials((fitData{m,8} == 1) & (fitData{m,5} == 1));
loTrials = trials((fitData{m,8} == 0) & (fitData{m,5} == 1));

resps1 = nanmean(nanmean(neuralData.eta.alignedResps{ETA}(hiTrials,postMoveIdx,whichCells),3),2);
vels1 = (nanmean(behavioralData.eta.alignedVels{ETA}(hiTrials,dsew_idx(postMoveIdx)),2));

resps2 = nanmean(nanmean(neuralData.eta.alignedResps{ETA}(loTrials,postMoveIdx,whichCells),3),2);
vels2 = (nanmean(behavioralData.eta.alignedVels{ETA}(loTrials,dsew_idx(postMoveIdx)),2));

figure;
scatter(resps1,vels1)
hold on
scatter(resps2,vels2)

corr(resps1,vels1)
corr(resps2,vels2)
