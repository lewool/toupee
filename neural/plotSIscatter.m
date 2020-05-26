allmID = [];
allbID = [];
allsID = [];
allrID = [];
for ex = 1:15


% nd = combinedNeuralData.matched;
nd = neuralData(ex);
[baselineResps, stimResps, pmovResps, movResps, rewResps] = getEpochResps(nd.eta);

%%
contrasts = getUniqueContrasts(expInfo);
[~, stimConds{1}] = selectCondition(expInfo(ex), contrasts(contrasts<0), behavioralData(ex), initTrialConditions('movementTime','late'));
[~, stimConds{2}] = selectCondition(expInfo(ex), contrasts(contrasts==0), behavioralData(ex), initTrialConditions('movementTime','late'));
[~, stimConds{3}] = selectCondition(expInfo(ex), contrasts(contrasts>0), behavioralData(ex), initTrialConditions('movementTime','late'));

[~, moveConds{1}] = selectCondition(expInfo(ex), contrasts, behavioralData(ex), initTrialConditions('movementDir','cw','movementTime','late'));
[~, moveConds{2}] = selectCondition(expInfo(ex), contrasts, behavioralData(ex), initTrialConditions('movementDir','ccw','movementTime','late'));

[~, rewConds{1}] = selectCondition(expInfo(ex), contrasts, behavioralData(ex), initTrialConditions('responseType','incorrect','movementTime','late'));
[~, rewConds{2}] = selectCondition(expInfo(ex), contrasts, behavioralData(ex), initTrialConditions('responseType','correct','movementTime','late'));

[~, blockConds{1}] = selectCondition(expInfo(ex), contrasts, behavioralData(ex), initTrialConditions('highRewardSide','left','movementTime','late'));
[~, blockConds{2}] = selectCondition(expInfo(ex), contrasts, behavioralData(ex), initTrialConditions('highRewardSide','right','movementTime','late'));

mLbL = intersect(moveConds{1},blockConds{1});
mRbL = intersect(moveConds{2},blockConds{1});
mLbR = intersect(moveConds{1},blockConds{2});
mRbR = intersect(moveConds{2},blockConds{2});

ns = min([numel(mLbL) numel(mLbR) numel(mRbL) numel(mRbR)]);

moveConds{1} = [randsample(mLbL,ns) randsample(mLbR,ns)];
moveConds{2} = [randsample(mRbL,ns) randsample(mRbR,ns)];
blockConds{1} = [randsample(mLbL,ns) randsample(mRbL,ns)];
blockConds{2} = [randsample(mLbR,ns) randsample(mRbR,ns)];

%%
if ex < 9
    whichCells = 'rightMov'; %choose from 'pLabels' array
    if strcmp(whichCells, 'all')
        plotCells = 1:size(nd.eta.alignedResps{1},3);
    else
        plotCells = find(nd.stats.bfcH(:,strcmp(nd.stats.labels,whichCells)) > 0);
    end
else
     whichCells = 'leftMov'; %choose from 'pLabels' array
    if strcmp(whichCells, 'all')
        plotCells = 1:size(nd.eta.alignedResps{1},3);
    else
        plotCells = find(nd.stats.bfcH(:,strcmp(nd.stats.labels,whichCells)) > 0);
    end
end
%% 
resps = movResps(:,plotCells);

if ex > 9
movDirIdx = (nanmean(resps(moveConds{2},:)) - nanmean(resps(moveConds{1},:)))./...
    (nanmean(resps(moveConds{2},:)) + nanmean(resps(moveConds{1},:)));
blockDirIdx = (nanmean(resps(blockConds{2},:)) - nanmean(resps(blockConds{1},:)))./...
    (nanmean(resps(blockConds{2},:)) + nanmean(resps(blockConds{1},:)));
stimDirIdx = (nanmean(resps(stimConds{3},:)) - nanmean(resps(stimConds{1},:)))./...
    (nanmean(resps(stimConds{3},:)) + nanmean(resps(stimConds{1},:)));
rewDirIdx = (nanmean(resps(rewConds{1},:)) - nanmean(resps(rewConds{2},:)))./...
    (nanmean(resps(rewConds{1},:)) + nanmean(resps(rewConds{2},:)));
else
    movDirIdx = (nanmean(resps(moveConds{1},:)) - nanmean(resps(moveConds{2},:)))./...
    (nanmean(resps(moveConds{1},:)) + nanmean(resps(moveConds{2},:)));
blockDirIdx = (nanmean(resps(blockConds{1},:)) - nanmean(resps(blockConds{2},:)))./...
    (nanmean(resps(blockConds{1},:)) + nanmean(resps(blockConds{2},:)));
stimDirIdx = (nanmean(resps(stimConds{1},:)) - nanmean(resps(stimConds{3},:)))./...
    (nanmean(resps(stimConds{1},:)) + nanmean(resps(stimConds{3},:)));
rewDirIdx = (nanmean(resps(rewConds{1},:)) - nanmean(resps(rewConds{2},:)))./...
    (nanmean(resps(rewConds{1},:)) + nanmean(resps(rewConds{2},:)));
end

allmID = [allmID movDirIdx];
allbID = [allbID blockDirIdx];
allsID = [allsID stimDirIdx];
allrID = [allrID rewDirIdx];

end

%%
axlim = 1;
hmax = 600;
yp = allbID;
xp = allmID;
color = [0 .4 1];
% color = [1 0 0];
% color = [0 0 0];
scatterFig = figure;
set(scatterFig,'position',[579 633 464 350]);
hold on;
line([-1 1],[0 0],'LineStyle','--','Color',[0.5 0.5 0.5]);
line([0 0],[-1 1],'LineStyle','--','Color',[0.5 0.5 0.5]);
scatter(xp,yp,20,'MarkerEdgeColor','none','MarkerFaceColor', color,'MarkerFaceAlpha', 0.2)
axis square
axis([-axlim axlim -axlim axlim])
set(gca,'ytick',[-1 -.5 0 .5 1]);
set(gca,'xtick',[-1 -.5 0 .5 1]);
set(gca,'tickdir','out')
ylabel('Block side index')
xlabel('Side index')
xpl = axes('Position',[.223 .11 .591 .15]);
xh = histogram(xp,linspace(-1,1,35));
ylim([0 hmax])
xlim([-axlim axlim]);
set(xh,'FaceColor',color,'FaceAlpha',.2);
axis off

ypl = axes('Position',[.21 .11 .11 .813]);
yh = histogram(yp,linspace(-1,1,35));
set(ypl,'view',[90 -90])
set(yh,'FaceColor',color,'FaceAlpha',.2);
ylim([0 hmax])
xlim([-axlim axlim]);

axis off

figure
hold on
yh = histogram(yp,linspace(-1,1,35),'normalization','probability');
set(yh,'FaceColor',color,'FaceAlpha',.2);
ylim([0 1])
xlim([-axlim axlim]);
