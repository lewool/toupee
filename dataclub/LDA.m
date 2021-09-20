numExps = 1;
stimDim = nan(9, 41, 2, length(numExps));
movDim = nan(9, 41, 2, length(numExps));
clear whichTrials;

for iX = numExps

%% initialize cells & response arrays

[baselineResps, stimResps, pmovResps, movResps, rewResps] = getEpochResps(neuralData(iX).eta);
whichETA = 1; %use stimulus-aligned
plotCells = getWhichCells('all',neuralData(iX));

%% set up trial conditions to compare
trialTypes = getTrialTypes(expInfo(iX), behavioralData(iX), 'late');
contrasts = getUniqueContrasts(expInfo(iX));

sL = [trialTypes.singleVar.contrast{1} trialTypes.singleVar.contrast{2}];
sR = [trialTypes.singleVar.contrast{8} trialTypes.singleVar.contrast{9}];
mL = cat(2,trialTypes.intVar.cb2D.contrast_direction{:,1});
mR = cat(2,trialTypes.intVar.cb2D.contrast_direction{:,2});

testTrials = 1:2:size(baselineResps,1);
trainTrials = 2:2:size(baselineResps,1);

%% some color stuff
zeroIdx = find(contrasts == 0);
walkup = length(contrasts) - zeroIdx;
walkback = zeroIdx - 1;
allColors = [0 0 .5; 0 0 1; 0 .4 1; .6 .8 1; .75 .75 .75; .8 .45 .45; 1 0 0; .5 0 0; .25 0 0];
zeroGray = find(allColors(:,1) == .75);
colors = allColors(zeroGray-walkback:zeroGray + walkup,:);
ms = 15;


%% compute dot products
clear whichResps
[baselineResps, stimResps, pmovResps, movResps, rewResps] = getEpochResps(neuralData(iX).eta);

eventWindow = neuralData(iX).eta.eventWindow;

dotProd_mL = nan(size(baselineResps,1),length(eventWindow));
dotProd_mR = nan(size(baselineResps,1),length(eventWindow));
dotProd_sL = nan(size(baselineResps,1),length(eventWindow));
dotProd_sR = nan(size(baselineResps,1),length(eventWindow));

meanMoveLeft = nan(1,length(baselineResps));
meanMoveRight = nan(1,length(baselineResps));
meanStimLeft = nan(1,length(baselineResps));
meanStimRight = nan(1,length(baselineResps));


meanMoveLeft = nanmean(movResps(intersect(mL,trainTrials),:),1);
meanMoveRight = nanmean(movResps(intersect(mR,trainTrials),:),1);
meanStimLeft = nanmean(stimResps(intersect(sL,trainTrials),:),1);
meanStimRight = nanmean(stimResps(intersect(sR,trainTrials),:),1);

for iT = 1:size(neuralData(iX).eta.alignedResps{1},2)
    whichResps = neuralData(iX).eta.alignedResps{whichETA}(:,iT,plotCells);

    fullNorm = 1; % 1 if just angle, 0 if angle and relative magnitude
    for iTrial = testTrials

        if fullNorm == 1
            tNorm = norm(whichResps(iTrial,:));
        else
            tNorm = 1;
        end

        dotProd_mL(iTrial,iT) = dot(whichResps(iTrial,:),meanMoveLeft)/(tNorm*norm(meanMoveLeft));
        dotProd_mR(iTrial,iT) = dot(whichResps(iTrial,:),meanMoveRight)/(tNorm*norm(meanMoveRight));
        dotProd_sL(iTrial,iT) = dot(whichResps(iTrial,:),meanStimLeft)/(tNorm*norm(meanStimLeft));
        dotProd_sR(iTrial,iT) = dot(whichResps(iTrial,:),meanStimRight)/(tNorm*norm(meanStimRight));

    end

end

dotProd_mL(dotProd_mL == 0) = NaN;
dotProd_mR(dotProd_mR == 0) = NaN;
dotProd_sL(dotProd_sL == 0) = NaN;
dotProd_sR(dotProd_sR == 0) = NaN;

%%

for c = 1:length(trialTypes.intVar.all.contrast_direction)
    whichTrials{c,1} = trialTypes.intVar.all.contrast_direction{c,1};
    whichTrials{c,2} = trialTypes.intVar.all.contrast_direction{c,2};
end
eventWindow = neuralData(iX).eta.eventWindow;

for iT = 1:length(eventWindow)
    for c = 1:length(whichTrials)
        if expInfo(iX).hemisphere < 0
            stimDim(c,iT,2,iX) = nanmean(dotProd_sL(whichTrials{c,1},iT) - dotProd_sR(whichTrials{c,1},iT));
            stimDim(c,iT,1,iX) = nanmean(dotProd_sL(whichTrials{c,2},iT) - dotProd_sR(whichTrials{c,2},iT));
            movDim(c,iT,2,iX) = nanmean(dotProd_mL(whichTrials{c,1},iT) - dotProd_mR(whichTrials{c,1},iT));
            movDim(c,iT,1,iX) = nanmean(dotProd_mL(whichTrials{c,2},iT) - dotProd_mR(whichTrials{c,2},iT));
        else
            stimDim(c,iT,1,iX) = nanmean(dotProd_sR(whichTrials{c,1},iT) - dotProd_sL(whichTrials{c,1},iT));
            stimDim(c,iT,2,iX) = nanmean(dotProd_sR(whichTrials{c,2},iT) - dotProd_sL(whichTrials{c,2},iT));
            movDim(c,iT,1,iX) = nanmean(dotProd_mR(whichTrials{c,1},iT) - dotProd_mL(whichTrials{c,1},iT));
            movDim(c,iT,2,iX) = nanmean(dotProd_mR(whichTrials{c,2},iT) - dotProd_mL(whichTrials{c,2},iT));
        end
    end
end
end

%%

%correct for hemisphere

for iX = numExps
    if expInfo(iX).hemisphere < 0
        stimDim(:,:,1,iX) = flipud(stimDim(:,:,1,iX));
        stimDim(:,:,2,iX) = flipud(stimDim(:,:,2,iX));
        movDim(:,:,2,iX) = flipud(movDim(:,:,2,iX));
        movDim(:,:,1,iX) = flipud(movDim(:,:,1,iX));
    end
end
%%
lim = .15;
tbi = 11:41;
xi = 1
figure;
hold on
subplot(1,2,1)
hold on
line([-lim lim],[0 0],'LineStyle','--','Color',[.5 .5 .5])
line([0 0],[-lim lim],'LineStyle','--','Color',[.5 .5 .5])
xlim([-lim lim]);
ylim([-lim lim]);
xlabel('stimulus dimension')
ylabel('movement dimension')
hold on
    box off
    axis square
subplot(1,2,2)
hold on
line([-lim lim],[0 0],'LineStyle','--','Color',[.5 .5 .5])
line([0 0],[-lim lim],'LineStyle','--','Color',[.5 .5 .5])
xlim([-lim lim]);
ylim([-lim lim]);
xlabel('stimulus dimension')
ylabel('movement dimension')
hold on
    box off
    axis square
%
for it = length(tbi)
for c = 1:4
    subplot(1,2,1)
    plot(smooth(nanmean(stimDim(c,tbi(1:it),1,xi),4)),smooth(nanmean(movDim(c,tbi(1:it),1,xi),4)),'k-','Color',colors(c,:),'LineWidth',2)
    hold on
    box off
    axis square
    set(gca,'tickdir','out')
    subplot(1,2,2)
    plot(smooth(nanmean(stimDim(c,tbi(1:it),2,xi),4)),smooth(nanmean(movDim(c,tbi(1:it),2,xi),4)),'k-','Color',colors(c,:),'LineWidth',2)
    hold on
    box off
    axis square
    set(gca,'tickdir','out')
end
    
for c = 6:9
    subplot(1,2,2)
    plot(smooth(nanmean(stimDim(c,tbi(1:it),1,xi),4)),smooth(nanmean(movDim(c,tbi(1:it),1,xi),4)),'k-','Color',colors(c,:),'LineWidth',2)
    hold on
    box off
    axis square
    set(gca,'tickdir','out')
    subplot(1,2,1)
    plot(smooth(nanmean(stimDim(c,tbi(1:it),2,xi),4)),smooth(nanmean(movDim(c,tbi(1:it),2,xi),4)),'k-','Color',colors(c,:),'LineWidth',2)
    hold on
    box off
    axis square
    set(gca,'tickdir','out')
end
% printfig(gcf,strcat('LEW movie',num2str(xi)))
end
%% plot stim resps
% epIdx = 15:20;
epIdx = 20:25;
% epIdx = 28:36;
% epIdx = 37:41;
for iX = 1:length(expInfo)
    for c = 1:length(whichTrials)
        meanContraMove_stim(c,iX) = nanmean(stimDim(c,epIdx,1,iX));
        meanIpsiMove_stim(c,iX) = nanmean(stimDim(c,epIdx,2,iX));
        meanContraMove_mov(c,iX) = nanmean(movDim(c,epIdx,1,iX));
        meanIpsiMove_mov(c,iX) = nanmean(movDim(c,epIdx,2,iX));
    end
end
%
lim=.15;
figure;
hold on
subplot(1,2,1)
hold on
line([-lim lim],[0 0],'LineStyle','--','Color',[.5 .5 .5])
line([0 0],[-lim lim],'LineStyle','--','Color',[.5 .5 .5])
xlim([-lim lim]);
ylim([-lim lim]);
axis square
box off
set(gca,'tickdir','out')
xlabel('stimulus dimension')
ylabel('movement dimension')
subplot(1,2,2)
hold on
line([-lim lim],[0 0],'LineStyle','--','Color',[.5 .5 .5])
line([0 0],[-lim lim],'LineStyle','--','Color',[.5 .5 .5])
xlim([-lim lim]);
ylim([-lim lim]);
axis square
box off
set(gca,'tickdir','out')
xlabel('stimulus dimension')
ylabel('movement dimension')

sessions = 1;
for c=1:4
    
    subplot(1,2,1)
    hold on
    plot(meanContraMove_stim(c,sessions),meanContraMove_mov(c,sessions),'ko','MarkerSize',8,'MarkerEdgeColor','none','MarkerFaceColor',colors(c,:),'LineStyle','none')
    subplot(1,2,2)
    hold on
    plot(meanIpsiMove_stim(c,sessions),meanIpsiMove_mov(c,sessions),'ko','MarkerSize',8,'MarkerEdgeColor','none','MarkerFaceColor',colors(c,:),'LineStyle','none')
end

for c=6:9
    
    subplot(1,2,2)
    hold on
    plot(meanContraMove_stim(c,sessions),meanContraMove_mov(c,sessions),'ko','MarkerSize',8,'MarkerEdgeColor','none','MarkerFaceColor',colors(c,:),'LineStyle','none')
    subplot(1,2,1)
    hold on
    plot(meanIpsiMove_stim(c,sessions),meanIpsiMove_mov(c,sessions),'ko','MarkerSize',8,'MarkerEdgeColor','none','MarkerFaceColor',colors(c,:),'LineStyle','none')
end
%%
lim=.06;
figure;
hold on
subplot(1,2,1)
hold on
line([-lim lim],[0 0],'LineStyle','--','Color',[.5 .5 .5])
line([0 0],[-lim lim],'LineStyle','--','Color',[.5 .5 .5])
xlim([-lim lim]);
ylim([-lim lim]);
axis square
box off
set(gca,'tickdir','out')
xlabel('stimulus dimension')
ylabel('movement dimension')

for c = 1:4
    if c < 5
        P2 = nanmean([meanContraMove_stim(c,:)',meanContraMove_mov(c,:)'],1);
        P1 = nanmean([meanIpsiMove_stim(c,:)',meanIpsiMove_mov(c,:)'],1);
        U =  P2(:,1)-P1(:,1);
        V =  P2(:,2)-P1(:,2);
        quiver(P1(:,1),P1(:,2),U,V,0,'LineWidth',3,'Color',colors(c,:),'MaxHeadSize',1)
    else
        P1 = nanmean([meanContraMove_stim(c,:)',meanContraMove_mov(c,:)'],1);
        P2 = nanmean([meanIpsiMove_stim(c,:)',meanIpsiMove_mov(c,:)'],1);
        U =  P2(:,1)-P1(:,1);
        V =  P2(:,2)-P1(:,2);
        quiver(P1(:,1),P1(:,2),U,V,0,'LineWidth',3,'Color',colors(c,:),'MaxHeadSize',1)
    end
end
subplot(1,2,2)
hold on
line([-lim lim],[0 0],'LineStyle','--','Color',[.5 .5 .5])
line([0 0],[-lim lim],'LineStyle','--','Color',[.5 .5 .5])
xlim([-lim lim]);
ylim([-lim lim]);
axis square
box off
set(gca,'tickdir','out')
xlabel('stimulus dimension')
ylabel('movement dimension')

for c = 6:9
    if c < 10
        P2 = nanmean([meanContraMove_stim(c,:)',meanContraMove_mov(c,:)'],1);
        P1 = nanmean([meanIpsiMove_stim(c,:)',meanIpsiMove_mov(c,:)'],1);
        U =  P2(:,1)-P1(:,1);
        V =  P2(:,2)-P1(:,2);
        quiver(P1(:,1),P1(:,2),U,V,0,'LineWidth',3,'Color',colors(c,:),'MaxHeadSize',1)
    else
     P1 = nanmean([meanContraMove_stim(c,:)',meanContraMove_mov(c,:)']);
        P2 = nanmean([meanIpsiMove_stim(c,:)',meanIpsiMove_mov(c,:)']);
        U =  P2(:,1)-P1(:,1);
        V =  P2(:,2)-P1(:,2);
        quiver(P1(:,1),P1(:,2),U,V,0,'LineWidth',3,'Color',colors(c,:),'MaxHeadSize',1)
    end
end

%% choose a subset of trials to visualize
iX = 1;
[baselineResps, stimResps, pmovResps, movResps, rewResps] = getEpochResps(neuralData(iX).eta);
whichETA = 1; %use stimulus-aligned
plotCells = getWhichCells('all',neuralData(iX));
trialTypes = getTrialTypes(expInfo(iX), behavioralData(iX), 'late');
contrasts = getUniqueContrasts(expInfo(iX));
axMin = -.2;
axMax = .6;
clear D1 D2

whichTrials{1} = trialTypes.intVar.all.side_direction{1,1};
whichTrials{2} = trialTypes.intVar.all.side_direction{3,2};
color{1} = colors(3,:);
color{2} = colors(7,:);

trialTypes = getTrialTypes(expInfo(iX), behavioralData(iX), 'late');
contrasts = getUniqueContrasts(expInfo(iX));

sL = trialTypes.singleVar.side{1};
sR = trialTypes.singleVar.side{3};
mL = trialTypes.singleVar.direction{1};
mR = trialTypes.singleVar.direction{2};
sLmL = trialTypes.intVar.all.side_direction{1,1};
sLmR = trialTypes.intVar.all.side_direction{1,2};
sRmL = trialTypes.intVar.all.side_direction{3,1};
sRmR = trialTypes.intVar.all.side_direction{3,2};

testTrials = 1:2:size(baselineResps,1);
trainTrials = 2:2:size(baselineResps,1);

dotProd_mL = nan(size(baselineResps,1),length(eventWindow));
dotProd_mR = nan(size(baselineResps,1),length(eventWindow));
dotProd_sL = nan(size(baselineResps,1),length(eventWindow));
dotProd_sR = nan(size(baselineResps,1),length(eventWindow));

meanMoveLeft = nan(1,length(baselineResps));
meanMoveRight = nan(1,length(baselineResps));
meanStimLeft = nan(1,length(baselineResps));
meanStimRight = nan(1,length(baselineResps));

meanMoveLeft = nanmean(movResps(intersect(mL,trainTrials),:),1);
meanMoveRight = nanmean(movResps(intersect(mR,trainTrials),:),1);
meanStimLeft = nanmean(stimResps(intersect(sL,trainTrials),:),1);
meanStimRight = nanmean(stimResps(intersect(sR,trainTrials),:),1);


for iT = 1:size(neuralData(iX).eta.alignedResps{1},2)
    whichResps = neuralData(iX).eta.alignedResps{whichETA}(:,iT,plotCells);

%     % compute means from 'train' trails
%     meanMoveLeft = nanmean(whichResps(intersect(mL,trainTrials),:),1);
%     meanMoveRight = nanmean(whichResps(intersect(mR,trainTrials),:),1);
%     meanStimLeft = nanmean(whichResps(intersect(sL,trainTrials),:),1);
%     meanStimRight = nanmean(whichResps(intersect(sR,trainTrials),:),1);

    fullNorm = 1; % 1 if just angle, 0 if angle and relative magnitude
    for iTrial = testTrials

        if fullNorm == 1
            tNorm = norm(whichResps(iTrial,:));
        else
            tNorm = 1;
        end

        dotProd_mL(iTrial,iT) = dot(whichResps(iTrial,:),meanMoveLeft)/(tNorm*norm(meanMoveLeft));
        dotProd_mR(iTrial,iT) = dot(whichResps(iTrial,:),meanMoveRight)/(tNorm*norm(meanMoveRight));
        dotProd_sL(iTrial,iT) = dot(whichResps(iTrial,:),meanStimLeft)/(tNorm*norm(meanStimLeft));
        dotProd_sR(iTrial,iT) = dot(whichResps(iTrial,:),meanStimRight)/(tNorm*norm(meanStimRight));

    end

end

dotProd_mL(dotProd_mL == 0) = NaN;
dotProd_mR(dotProd_mR == 0) = NaN;
dotProd_sL(dotProd_sL == 0) = NaN;
dotProd_sR(dotProd_sR == 0) = NaN;

for iT = 1:length(eventWindow)
    D1(1,iT) = nanmean(dotProd_sR(whichTrials{1},iT) - dotProd_sL(whichTrials{1},iT));
    D1(2,iT) = nanmean(dotProd_sR(whichTrials{2},iT) - dotProd_sL(whichTrials{2},iT));
    D2(1,iT) = nanmean(dotProd_mR(whichTrials{1},iT) - dotProd_mL(whichTrials{1},iT));
    D2(2,iT) = nanmean(dotProd_mR(whichTrials{2},iT) - dotProd_mL(whichTrials{2},iT));
end

for iT = 1:length(eventWindow)
    D1_all{1,iT} = (dotProd_sR(whichTrials{1},iT) - dotProd_sL(whichTrials{1},iT));
    D1_all{2,iT} = (dotProd_sR(whichTrials{2},iT) - dotProd_sL(whichTrials{2},iT));
    D2_all{1,iT} = (dotProd_mR(whichTrials{1},iT) - dotProd_mL(whichTrials{1},iT));
    D2_all{2,iT} = (dotProd_mR(whichTrials{2},iT) - dotProd_mL(whichTrials{2},iT));
end

%% figure



eventWindow = neuralData(iX).eta.eventWindow;

for iT = 1:length(eventWindow)
subplot(1,3,1)
line([axMin axMax],[axMin axMax],'LineStyle','--','Color',[.5 .5 .5])
hold on
scatter(dotProd_sL(whichTrials{1},iT),dotProd_sR(whichTrials{1},iT),ms*2,'MarkerEdgeColor','none','MarkerFaceColor',color{1},'MarkerFaceAlpha', .7)
scatter(dotProd_sL(whichTrials{2},iT),dotProd_sR(whichTrials{2},iT),ms*2,'MarkerEdgeColor','none','MarkerFaceColor',color{2},'MarkerFaceAlpha', .7)
xlim([axMin axMax]);
ylim([axMin axMax]);
axis square
box off
set(gca,'tickdir','out')
xlabel('mean R_{leftStim}')
ylabel('mean R_{rightStim}')

subplot(1,3,2)
line([axMin axMax],[axMin axMax],'LineStyle','--','Color',[.5 .5 .5])
hold on
scatter(dotProd_mL(whichTrials{1},iT),dotProd_mR(whichTrials{1},iT),ms*2,'MarkerEdgeColor','none','MarkerFaceColor',color{1},'MarkerFaceAlpha', .7)
scatter(dotProd_mL(whichTrials{2},iT),dotProd_mR(whichTrials{2},iT),ms*2,'MarkerEdgeColor','none','MarkerFaceColor',color{2},'MarkerFaceAlpha', .7)
xlim([axMin axMax]);
ylim([axMin axMax]);
axis square
box off
set(gca,'tickdir','out')
xlabel('mean R_{leftMov}')
ylabel('mean R_{rightMov}')

subplot(1,3,3)
hold on
line([-.4 .4],[0 0],'LineStyle','--','Color',[.5 .5 .5])
line([0 0],[-.4 .4],'LineStyle','--','Color',[.5 .5 .5])
xlim([-.4 .4]);
ylim([-.4 .4]);
plot(D1(1,iT),D2(1,iT),'o','MarkerSize',ms/2,'MarkerEdgeColor','none','MarkerFaceColor',color{1})
plot(D1(2,iT),D2(2,iT),'o','MarkerSize',ms/2,'MarkerEdgeColor','none','MarkerFaceColor',color{2})
axis square
box off
set(gca,'tickdir','out')
ylabel('R_{rightMov} – R_{leftMov}')
xlabel('R_{rightStim} – R_{leftStim}')
printfig(gcf,strcat('LDA scatter movie',num2str(iT)))
subplot(1,3,1)
cla
subplot(1,3,2)
cla
subplot(1,3,3)
cla
end

%%
figure;
subplot(1,3,1)
iT = 15;
scatter(D1_all{1,iT},D2_all{1,iT},ms*2,'MarkerEdgeColor','none','MarkerFaceColor',color{1},'MarkerFaceAlpha', .2)
hold on
scatter(D1_all{2,iT},D2_all{2,iT},ms*2,'MarkerEdgeColor','none','MarkerFaceColor',color{2},'MarkerFaceAlpha', .2)
plot(D1(1,iT),D2(1,iT),'o','MarkerSize',ms/2,'MarkerEdgeColor','k','MarkerFaceColor',color{1})
plot(D1(2,iT),D2(2,iT),'o','MarkerSize',ms/2,'MarkerEdgeColor','k','MarkerFaceColor',color{2})

line([-.4 .4],[0 0],'LineStyle','--','Color',[.5 .5 .5])
line([0 0],[-.4 .4],'LineStyle','--','Color',[.5 .5 .5])
xlim([-.2 .2]);
ylim([-.2 .2]);
axis square
box off
set(gca,'tickdir','out')
ylabel('choice dimension')
xlabel('stimulus dimension')
yticks([-.2 -.1 0 .1 .2])
set(gca, 'YTickLabels', {'-0.2', '-0.1', '0','0.1', '0.2'})

subplot(1,3,2)
iT = 25;
scatter(D1_all{1,iT},D2_all{1,iT},ms*2,'MarkerEdgeColor','none','MarkerFaceColor',color{1},'MarkerFaceAlpha', .2)
hold on
scatter(D1_all{2,iT},D2_all{2,iT},ms*2,'MarkerEdgeColor','none','MarkerFaceColor',color{2},'MarkerFaceAlpha', .2)
plot(D1(1,iT),D2(1,iT),'o','MarkerSize',ms/2,'MarkerEdgeColor','k','MarkerFaceColor',color{1})
plot(D1(2,iT),D2(2,iT),'o','MarkerSize',ms/2,'MarkerEdgeColor','k','MarkerFaceColor',color{2})
line([-.4 .4],[0 0],'LineStyle','--','Color',[.5 .5 .5])
line([0 0],[-.4 .4],'LineStyle','--','Color',[.5 .5 .5])
xlim([-.2 .2]);
ylim([-.2 .2]);
yticks([-.2 -.1 0 .1 .2])
set(gca, 'YTickLabels', {'-0.2', '-0.1', '0','0.1', '0.2'})

axis square
box off
set(gca,'tickdir','out')
ylabel('choice dimension')
xlabel('stimulus dimension')

subplot(1,3,3)
iT = 35;
scatter(D1_all{1,iT},D2_all{1,iT},ms*2,'MarkerEdgeColor','none','MarkerFaceColor',color{1},'MarkerFaceAlpha', .2)
hold on
scatter(D1_all{2,iT},D2_all{2,iT},ms*2,'MarkerEdgeColor','none','MarkerFaceColor',color{2},'MarkerFaceAlpha', .2)
plot(D1(1,iT),D2(1,iT),'o','MarkerSize',ms/2,'MarkerEdgeColor','k','MarkerFaceColor',color{1})
plot(D1(2,iT),D2(2,iT),'o','MarkerSize',ms/2,'MarkerEdgeColor','k','MarkerFaceColor',color{2})
line([-.4 .4],[0 0],'LineStyle','--','Color',[.5 .5 .5])
line([0 0],[-.4 .4],'LineStyle','--','Color',[.5 .5 .5])
xlim([-.2 .2]);
ylim([-.2 .2]);
axis square
box off
set(gca,'tickdir','out')
ylabel('choice dimension')
xlabel('stimulus dimension')
yticks([-.2 -.1 0 .1 .2])
set(gca, 'YTickLabels', {'-0.2', '-0.1', '0','0.1', '0.2'})

