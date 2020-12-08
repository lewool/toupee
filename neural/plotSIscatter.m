allmID = [];
allbID = [];
allsID = [];
allrID = [];
for ex = 1


% nd = combinedNeuralData.matched;
nd = neuralData(ex);
[baselineResps, stimResps, pmovResps, movResps, rewResps] = getEpochResps(nd.eta);

trialTypes = getTrialTypes(expInfo(ex), behavioralData(ex), 'late');

%%
contrasts = getUniqueContrasts(expInfo);
stimConds{1} = trialTypes.intVar.cb2D.side{1};
stimConds{2} = trialTypes.intVar.cb2D.side{2};
stimConds{3} = trialTypes.intVar.cb2D.side{3};

moveConds{1} = trialTypes.intVar.cb2D.direction{1};
moveConds{2} = trialTypes.intVar.cb2D.direction{2};

rewConds{1} = trialTypes.intVar.cb2D.outcome{2};
rewConds{2} = trialTypes.intVar.cb2D.outcome{1};

blockConds{1} = trialTypes.intVar.cb3D.block{1};
blockConds{2} = trialTypes.intVar.cb3D.block{2};


%%
if ex < 9
    whichCells = 'stim'; %choose from 'pLabels' array
    if strcmp(whichCells, 'all')
        plotCells = 1:size(nd.eta.alignedResps{1},3);
    else
        plotCells = find(nd.stats.bfcH(:,strcmp(nd.stats.labels,whichCells)) > 0);
    end
else
     whichCells = 'advanceMov'; %choose from 'pLabels' array
    if strcmp(whichCells, 'all')
        plotCells = 1:size(nd.eta.alignedResps{1},3);
    else
        plotCells = find(nd.stats.bfcH(:,strcmp(nd.stats.labels,whichCells)) > 0);
    end
end
%% 
resps = stimResps(:,plotCells) - baselineResps(:,plotCells);

if ex > 0
movDirIdx = (nanmean(resps(moveConds{2},:)) - nanmean(resps(moveConds{1},:)))./...
    (nanmean(resps(moveConds{2},:)) + nanmean(resps(moveConds{1},:)));
blockDirIdx = (nanmean(resps(blockConds{2},:)) - nanmean(resps(blockConds{1},:)))./...
    (nanmean(resps(blockConds{2},:)) + nanmean(resps(blockConds{1},:)));
stimDirIdx = (nanmean(resps(stimConds{3},:)) - nanmean(resps(stimConds{1},:)))./...
    (nanmean(resps(stimConds{3},:)) + nanmean(resps(stimConds{1},:)));
rewDirIdx = (nanmean(resps(rewConds{2},:)) - nanmean(resps(rewConds{1},:)))./...
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
hmax = 10;
yp = allsID;
xp = allmID;
% color = [0 .4 1];
% color = [1 0 0];
color = [0 0 0];
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
ylabel('Left stim   <-------->   Right stim')
xlabel('Left choice   <-------->   Right choice')
% xlabel('Left stimulus   <-------->   Right stimulus')
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
set(gcf,'renderer','Painters');
% figure
% hold on
% yh = histogram(yp,linspace(-1,1,35),'normalization','probability');
% set(yh,'FaceColor',color,'FaceAlpha',.2);
% ylim([0 1])
% xlim([-axlim axlim]);
