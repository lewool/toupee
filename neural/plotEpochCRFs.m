function plotEpochCRFs(expInfo, neuralData, CRF, subpop, expIdx)

if ~exist('expIdx')
    expIdx = 1:length(expInfo);
end
contrasts = getUniqueContrasts(expInfo);
eventWindow = neuralData(1).eta.eventWindow;

%% initialize the plot

zeroIdx = find(contrasts == 0);
walkup = length(contrasts) - zeroIdx;
walkback = zeroIdx - 1;
% allColors = [.25 0 0;.5 0 0 ;1 0 0;.8 .45 .45;.75 .75 .75;.55 .55 .55;.35 .35 .35;.15 .15 .15;0 0 0];
allColors = [0 0 .5; 0 0 1; 0 .4 1; .6 .8 1; .75 .75 .75; .8 .45 .45; 1 0 0; .5 0 0; .25 0 0];

zeroGray = find(allColors(:,1) == .75);
colors = allColors(zeroGray-walkback:zeroGray + walkup,:);
dirColors = [0 .4 1; 1 0 0];

%% get the data

bCRF_mC = CRF.(matlab.lang.makeValidName(subpop)).mC.bCRF(:,expIdx);
bCRF_mI = CRF.(matlab.lang.makeValidName(subpop)).mI.bCRF(:,expIdx);
sCRF_mC = CRF.(matlab.lang.makeValidName(subpop)).mC.sCRF(:,expIdx);
sCRF_mI = CRF.(matlab.lang.makeValidName(subpop)).mI.sCRF(:,expIdx);
pCRF_mC = CRF.(matlab.lang.makeValidName(subpop)).mC.pCRF(:,expIdx);
pCRF_mI = CRF.(matlab.lang.makeValidName(subpop)).mI.pCRF(:,expIdx);
mCRF_mC = CRF.(matlab.lang.makeValidName(subpop)).mC.mCRF(:,expIdx);
mCRF_mI = CRF.(matlab.lang.makeValidName(subpop)).mI.mCRF(:,expIdx);

%% plotting

maxY = 1.1*max(max([nanmean(bCRF_mC,2) nanmean(bCRF_mI,2) nanmean(sCRF_mC,2) nanmean(sCRF_mI,2) nanmean(mCRF_mC,2) mean(mCRF_mI,2)]));
minY = 1.1*min(min([nanmean(bCRF_mC,2) nanmean(bCRF_mI,2) nanmean(sCRF_mC,2) nanmean(sCRF_mI,2) nanmean(mCRF_mC,2) mean(mCRF_mI,2)]));
minX = contrasts(1)*1.1;
maxX = contrasts(end)*1.1;

figure;hold on
set(gcf,'position',[79 576 1708 321])

subplot(1,4,1)
errorbar(contrasts, nanmean(bCRF_mC,2),nanstd(bCRF_mC,[],2)/sqrt(size(bCRF_mC,2)),...
    'o','MarkerFaceColor','k','LineStyle','-','color','k','MarkerEdgeColor','none','capsize',0);
hold on;
errorbar(contrasts, nanmean(bCRF_mI,2),nanstd(bCRF_mI,[],2)/sqrt(size(bCRF_mI,2)),...
    'o','MarkerFaceColor','k','LineStyle',':','color','k','MarkerEdgeColor','none','capsize',0);
for pC = 1:length(contrasts)
    plot(contrasts(pC),nanmean(bCRF_mC(pC,:),2),'o','MarkerFaceColor',colors(pC,:),'MarkerEdgeColor','none')
    plot(contrasts(pC),nanmean(bCRF_mI(pC,:),2),'o','MarkerFaceColor',colors(pC,:),'MarkerEdgeColor','none')
end
ylim([minY maxY])
xlim([minX maxX])
box off
set(gca,'tickdir','out')
title('Baseline')
ylabel('Neural activity')
xlabel('Contrast')

subplot(1,4,2)
errorbar(contrasts, nanmean(sCRF_mC,2),nanstd(sCRF_mC,[],2)/sqrt(size(sCRF_mC,2)),...
    'LineStyle','-','color','k','MarkerEdgeColor','none','capsize',0);
hold on;
errorbar(contrasts, nanmean(sCRF_mI,2),nanstd(sCRF_mI,[],2)/sqrt(size(sCRF_mI,2)),...
    'LineStyle',':','color','k','MarkerEdgeColor','none','capsize',0);
for pC = 1:length(contrasts)
    plot(contrasts(pC),nanmean(sCRF_mC(pC,:),2),'o','MarkerFaceColor',colors(pC,:),'MarkerEdgeColor','none')
    plot(contrasts(pC),nanmean(sCRF_mI(pC,:),2),'o','MarkerFaceColor',colors(pC,:),'MarkerEdgeColor','none')
end
ylim([minY maxY])
xlim([minX maxX])
box off
set(gca,'tickdir','out')
title('Stimulus')

subplot(1,4,3)
errorbar(contrasts, nanmean(pCRF_mC,2),nanstd(pCRF_mC,[],2)/sqrt(size(pCRF_mC,2)),...
    'o','MarkerFaceColor','k','LineStyle','-','color','k','MarkerEdgeColor','none','capsize',0);
hold on;
errorbar(contrasts, nanmean(pCRF_mI,2),nanstd(pCRF_mI,[],2)/sqrt(size(pCRF_mI,2)),...
    'o','MarkerFaceColor','k','LineStyle',':','color','k','MarkerEdgeColor','none','capsize',0);
for pC = 1:length(contrasts)
    plot(contrasts(pC),nanmean(pCRF_mC(pC,:),2),'o','MarkerFaceColor',colors(pC,:),'MarkerEdgeColor','none')
    plot(contrasts(pC),nanmean(pCRF_mI(pC,:),2),'o','MarkerFaceColor',colors(pC,:),'MarkerEdgeColor','none')
end
ylim([minY maxY])
box off
xlim([minX maxX])
set(gca,'tickdir','out')
title('Pre-movement')

subplot(1,4,4)
errorbar(contrasts, nanmean(mCRF_mC,2),nanstd(mCRF_mC,[],2)/sqrt(size(mCRF_mC,2)),...
    'o','MarkerFaceColor','k','LineStyle','-','color','k','MarkerEdgeColor','none','capsize',0);
hold on;
errorbar(contrasts, nanmean(mCRF_mI,2),nanstd(mCRF_mI,[],2)/sqrt(size(mCRF_mI,2)),...
    'o','MarkerFaceColor','k','LineStyle',':','color','k','MarkerEdgeColor','none','capsize',0);
for pC = 1:length(contrasts)
    plot(contrasts(pC),nanmean(mCRF_mC(pC,:),2),'o','MarkerFaceColor',colors(pC,:),'MarkerEdgeColor','none')
    plot(contrasts(pC),nanmean(mCRF_mI(pC,:),2),'o','MarkerFaceColor',colors(pC,:),'MarkerEdgeColor','none')
end
ylim([minY maxY])
xlim([minX maxX])
box off
set(gca,'tickdir','out')
title('Movement')
