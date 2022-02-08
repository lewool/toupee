subDim = ceil(sqrt(length(mouseList)));
figure(1);
for m = 1:length(mouseList)
    
    mouseName = char(mouseList{m});
    expDate = char(expList{m}{1});
    expNum = expList{m}{2};
    
    [expRef, ~] = data.constructExpRef(mouseName,expDate,expNum);
    fprintf('fetching %s...\n',expRef)
    [behavioralData, expInfo, neuralData] = data.loadDataset(mouseName, expDate, expNum);
    expInfo.hemisphere = hemList(m);
    [baselineResps, stimResps, pmovResps, movResps, rewResps, preCueResps] = getEpochResps(neuralData.eta);


    %% discard trials with early movements
    contrasts = getUniqueContrasts(expInfo);
    [~, whichTrials] = selectCondition(expInfo, contrasts, behavioralData, initTrialConditions('movementTime','late'));
    trialTypes = getTrialTypes(expInfo, behavioralData, 'late');

    [~, earlyTrials] = selectCondition(expInfo, getUniqueContrasts(expInfo), behavioralData, ...
        initTrialConditions('preStimMovement','quiescent','movementTime','early','specificRTs',[.1 inf]));

    [~, lateTrials] = selectCondition(expInfo, getUniqueContrasts(expInfo), behavioralData, ...
        initTrialConditions('preStimMovement','quiescent','movementTime','late','specificRTs',[.1 inf]));

    [correct_highSide, correct_lowSide, incorrect] = getHighLoRewardTrials(expInfo, behavioralData);
    whichCells = 1:length(neuralData.stats.bfcH(:,5));


%% compute action vs stimulus coding for each cell

    action = nanmean((movResps(whichTrials,:)));%-baselineResps(whichTrials,:)),1);
    stim = nanmean((stimResps(whichTrials,:)));%-baselineResps(whichTrials,:)),1);
    rew = nanmean((rewResps(whichTrials,:)));%-baselineResps(whichTrials,:)),1);
    relMovResps = movResps;
    relStimResps = stimResps;
    relRewResps = rewResps;

    if expInfo.hemisphere < 0

        dir = ... 
            (nanmean(relMovResps(trialTypes.intVar.all.contrast_direction{5,1},whichCells)) - ...
            nanmean(relMovResps(trialTypes.intVar.all.contrast_direction{5,2},whichCells)));

        side = ...
            (nanmean(relStimResps(trialTypes.singleVar.side{1},whichCells)) - ...
            nanmean(relStimResps(trialTypes.singleVar.side{3},whichCells)));

        feed = ...
            (nanmean(relRewResps(trialTypes.intVar.all.contrast_outcome{5,1},whichCells)) - ...
            nanmean(relRewResps(trialTypes.intVar.all.contrast_outcome{5,2},whichCells)));

        timing = ...
            (nanmean(relMovResps(earlyTrials,whichCells)) - ...
            nanmean(relMovResps(lateTrials,whichCells)));

        hilo = ...
            (nanmean(relRewResps(correct_highSide,whichCells)) - ...
            nanmean(relRewResps(correct_lowSide,whichCells)));
    else

        dir = ... 
            (nanmean(relMovResps(trialTypes.intVar.all.contrast_direction{5,2},whichCells)) - ...
            nanmean(relMovResps(trialTypes.intVar.all.contrast_direction{5,1},whichCells)));

        side = ...
            (nanmean(relStimResps(trialTypes.singleVar.side{3},whichCells)) - ...
            nanmean(relStimResps(trialTypes.singleVar.side{1},whichCells)));

        feed = ...
            (nanmean(relRewResps(trialTypes.intVar.all.contrast_outcome{5,1},whichCells)) - ...
            nanmean(relRewResps(trialTypes.intVar.all.contrast_outcome{5,2},whichCells)));

        timing = ...
            (nanmean(relMovResps(earlyTrials,whichCells)) - ...
            nanmean(relMovResps(lateTrials,whichCells)));

        hilo = ...
            (nanmean(relRewResps(correct_highSide,whichCells)) - ...
            nanmean(relRewResps(correct_lowSide,whichCells)));
    end

%%
    meanColor = [1 0 0];
    ms = 4;
    mm = 50;
    ticksmax = .3;
    % close all;
    dirIndex = dir;
    sideIndex = side;
    rewIdx = feed;
    timeIdx = timing;
    hiloIdx = hilo;
    [~,sortDirIdx] = sort(dirIndex);
    [~,sortRewIdx] = sort(rewIdx);
    
    figure(1);
    subplot(subDim,subDim,m)
    lim = max([abs(sideIndex), abs(dirIndex)])*1.1;
    line([0 0],[-lim lim],'LineStyle',':','Color',[.5 .5 .5])
    hold on;
    line([-lim lim],[0 0],'LineStyle',':','Color',[.5 .5 .5])
    scatter(dirIndex,sideIndex,ms,...
        'MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha', .7);
    xlim([-lim lim]);
    ylim([-lim lim]);
    xticks([-ticksmax 0 ticksmax])
    yticks([-ticksmax 0 ticksmax])
    axis square
    box off
    set(gca,'tickdir','out')
    ylabel('stim selectivity')
    xlabel('choice selectivity')

%% WORKBENCH
figure;
set(gcf,'position', [1000 1050 948 290])
subplot(1,4,1)
lim = max([abs(sideIndex), abs(dirIndex)])*1.1;
line([0 0],[-lim lim],'LineStyle',':','Color',[.5 .5 .5])
hold on;
line([-lim lim],[0 0],'LineStyle',':','Color',[.5 .5 .5])
scatter(dirIndex,sideIndex,ms,...
    'MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha', .7);
xlim([-lim lim]);
ylim([-lim lim]);
xticks([-.2 0 .2])
yticks([-.2 0 .2])
axis square
box off
set(gca,'tickdir','out')
ylabel('stim selectivity')
xlabel('choice selectivity')

subplot(1,4,2)
lim = max([abs(rewIdx), abs(dirIndex)])*1.1;
line([0 0],[-lim lim],'LineStyle',':','Color',[.5 .5 .5])
hold on;
line([-lim lim],[0 0],'LineStyle',':','Color',[.5 .5 .5])
scatter(dirIndex,rewIdx,ms,...
    'MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha', .7);
xlim([-lim lim]);
ylim([-lim lim]);
xticks([-.2 0 .2])
yticks([-.2 0 .2])
axis square
box off
set(gca,'tickdir','out')
ylabel('reward selectivity')
xlabel('choice selectivity')

subplot(1,4,3)
lim = max([abs(sideIndex), abs(rewIdx)])*1.1;
line([0 0],[-lim lim],'LineStyle',':','Color',[.5 .5 .5])
hold on;
line([-lim lim],[0 0],'LineStyle',':','Color',[.5 .5 .5])
scatter(sideIndex,rewIdx,ms,...
    'MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha', .7);
xlim([-lim lim]);
ylim([-lim lim]);
xticks([-.2 0 .2])
yticks([-.2 0 .2])
axis square
box off
set(gca,'tickdir','out')
ylabel('reward selectivity')
xlabel('stim selectivity')

subplot(1,4,4)
lim = max([abs(rewIdx), abs(timeIdx)])*1.1;
line([0 0],[-lim lim],'LineStyle',':','Color',[.5 .5 .5])
hold on;
line([-lim lim],[0 0],'LineStyle',':','Color',[.5 .5 .5])
scatter(rewIdx,timeIdx,ms,...
    'MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha', .7);
xlim([-lim lim]);
ylim([-lim lim]);
xticks([-.2 0 .2])
yticks([-.2 0 .2])
axis square
box off
set(gca,'tickdir','out')
ylabel('impulsivity selectivity')
xlabel('reward selectivity')

fit{m}.dir_side = polyfit(dirIndex,sideIndex,1);
fit{m}.dir_reward = polyfit(dirIndex,rewIdx,1);
fit{m}.side_reward = polyfit(sideIndex,rewIdx,1);
fit{m}.imp_reward = polyfit(timeIdx,rewIdx,1);

correlation{m}.dir_side = corr(dirIndex',sideIndex');
correlation{m}.dir_reward = corr(dirIndex',rewIdx');
correlation{m}.side_reward = corr(sideIndex',rewIdx');
correlation{m}.imp_reward = corr(timeIdx',rewIdx');

clearvars -except mouseList expList hemList subDim fit correlation
end


%%
for m = 1:length(mouseList)
    corr_ds(m) = correlation{m}.dir_side;
    corr_dr(m) = correlation{m}.dir_reward;
    corr_sr(m) = correlation{m}.side_reward;
    corr_ir(m) = correlation{m}.imp_reward;

    
    
    
end

corr_ds(isnan(corr_ds))=[];
corr_dr(isnan(corr_dr))=[];
corr_sr(isnan(corr_sr))=[];
corr_ir(isnan(corr_ir))=[];
%%
figure;
beeswarm([...
    ones(1,length(corr_ds))'; ...
    2*ones(1,length(corr_dr))'; ...
    3*ones(1,length(corr_sr))';
    4*ones(1,length(corr_ir))'],[...
    corr_ds'; 
    corr_dr'; 
    corr_sr';
    corr_ir'],...
    'sort_style','hex','dot_size',1,'overlay_style','sd');

ylim([-.8 .8]);
set(gca,'tickdir','out')
xticks([1 2 3 4])
set(gca, 'XTickLabels', {'Stim vs. choice' 'Choice vs. feedback' 'Stim vs. feedback' 'Impulsivity vs. feedback' })
ylabel('Correlation (r)')
ylim([-.8 .800001]);
ll = line([.5 4.5],[0 0],'LineStyle','--','Color',[.5 .5 .5]);
uistack(ll,'bottom');










% clearvars -except mouseList expList hemList
