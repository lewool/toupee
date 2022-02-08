zeroIdx = find(contrasts == 0);
walkup = length(contrasts) - zeroIdx;
walkback = zeroIdx - 1;
allColors = [0 0 .5; 0 0 1; 0 .4 1; .6 .8 1; .75 .75 .75; .8 .45 .45; 1 0 0; .5 0 0; .25 0 0];

zeroGray = find(allColors(:,1) == .75);
colors = allColors(zeroGray-walkback:zeroGray + walkup,:);

%%

clear R_mL R_mR vel_mR vel_mL
iX_range = 1:14;
for iX = 1:14
    % plot responses to contrast x movement
    [baselineResps, stimResps, pmovResps, movResps, rewResps] = getEpochResps(neuralData(iX).eta);

%     whichNeurons = 1;
    whichNeurons = neuralData(iX).stats.bfcH(:,6) == 1 ;
%     whichNeurons = 1:length(neuralData(iX).stats.bfcH);

    trialTypes = getTrialTypes(expInfo(iX), behavioralData(iX), 'late');
    whichTrials = trialTypes.intVar.all.contrast_direction;
    vels = ((behavioralData(iX).wheelMoves.epochs(5).peakVel));
    for c = 1:size(whichTrials,1)
        mL = whichTrials{c,2};
        mR = whichTrials{c,2};

        R_mL(c,iX) = nanmean(nanmean(movResps(mL,whichNeurons),2));
        R_mR(c,iX) = nanmean(nanmean(movResps(mR,whichNeurons),2));
        vel_mL(c,iX) = nanmean(vels(mL))';
        vel_mR(c,iX) = nanmean(vels(mR))';

    end   


end

figure(1);
hold on
set(gca,'tickdir','out')
box off
xlabel('Mean wheel velocity')
ylabel('Mean neural activity')

for c = 1:size(whichTrials,1)
    errorbar(...
        nanmean(vel_mL(c,iX_range)),...
        nanmean(R_mL(c,iX_range)),...
        -nanstd(R_mL(c,iX_range))/sqrt(size(R_mL(c,iX_range),2)),...
        nanstd(R_mL(c,iX_range))/sqrt(size(R_mL(c,iX_range),2)),...
        -nanstd(vel_mL(c,iX_range))/sqrt(size(vel_mL(c,iX_range),2)),...
        nanstd(vel_mL(c,iX_range))/sqrt(size(vel_mL(c,iX_range),2)),...
        'o','MarkerFaceColor','k','LineStyle','none','color','k','MarkerEdgeColor','none','MarkerFaceColor',colors(c,:),'capsize',0);
    p2 = plot(...
    nanmean(vel_mL(c,iX_range)),...
    nanmean(R_mL(c,iX_range)),...
    'o','MarkerSize',8,'MarkerFaceColor','k','LineStyle','none','color','k','MarkerEdgeColor','none','MarkerFaceColor',colors(c,:));
end

figure(2);
hold on
set(gca,'tickdir','out')
box off
xlabel('Mean wheel velocity')
ylabel('Mean neural activity')

for c = 1:size(whichTrials,1)
plot(...
    (vel_mL(c,iX_range)),...
    (R_mL(c,iX_range)),...
    'o','MarkerSize',8,'MarkerFaceColor','k','LineStyle','none','color','k','MarkerEdgeColor','none','MarkerFaceColor',colors(c,:));
% %     errorbar(...
%         nanmean(vel_mR{c}),...
%         nanmean(R_mR{c}),...
%         -nanstd(R_mR{c})/sqrt(size(R_mR{c},1)),...
%         nanstd(R_mR{c})/sqrt(size(R_mR{c},1)),...
%         -nanstd(vel_mR{c})/sqrt(size(vel_mR{c},1)),...
%         nanstd(vel_mR{c})/sqrt(size(vel_mR{c},1)),...
%         'o','MarkerFaceColor','k','LineStyle','-','color','k','MarkerEdgeColor','none','MarkerFaceColor',colors(c,:),'capsize',0);

end

%%
figure;
hold on;

for iX = 1:5
    % plot responses to contrast x movement
    [baselineResps, stimResps, pmovResps, movResps, rewResps] = getEpochResps(neuralData(iX).eta);

    whichNeurons = 1;
    whichNeurons = neuralData(iX).stats.bfcH(:,6) == 1 ;
%     whichNeurons = 1:length(neuralData(iX).stats.bfcH);

    trialTypes = getTrialTypes(expInfo(iX), behavioralData(iX), 'late');
    whichTrials = trialTypes.intVar.all.contrast_direction;
    clear R_mL R_mR vel_mR vel_mL
    for c = 1:size(whichTrials,1)
        mL = whichTrials{c,1};
        mR = whichTrials{c,2};

        R_mL{c} = nanmean(baselineResps(mL,whichNeurons),2);
        R_mR{c} = nanmean(baselineResps(mR,whichNeurons),2);
        vel_mL{c} = -1;
        vel_mR{c} = 1;

    end   

    for c = 1:size(whichTrials,1)
        line([nanmean(vel_mL{c}) nanmean(vel_mR{c})],[nanmean(R_mL{c}) nanmean(R_mR{c})],'LineStyle','-','Color',[.5 .5 .5])
        plot(nanmean(vel_mL{c}),nanmean(R_mL{c}),'ko','MarkerFaceColor',colors(c,:),'MarkerEdgeColor','none')
        plot(nanmean(vel_mR{c}),nanmean(R_mR{c}),'ko','MarkerFaceColor',colors(c,:),'MarkerEdgeColor','none')
    end

end










%%
figure; 
hold on

for c = 1:size(whichTrials,1)
    scatter(R_mL{c},vel_mL{c},40,'MarkerEdgeColor','none','MarkerFaceColor',colors(c,:),'MarkerFaceAlpha', 1)
    scatter(R_mR{c},vel_mR{c},40,'MarkerEdgeColor','none','MarkerFaceColor',colors(c,:),'MarkerFaceAlpha', 1)
end
%%
figure;
hold on;
for c = 1:size(whichTrials,1)
    errorbar(...
        nanmean(vel_mL{c}),...
        nanmean(R_mL{c}),...
        -nanstd(R_mL{c})/sqrt(size(R_mL{c},1)),...
        nanstd(R_mL{c})/sqrt(size(R_mL{c},1)),...
        -nanstd(vel_mL{c})/sqrt(size(vel_mL{c},1)),...
        nanstd(vel_mL{c})/sqrt(size(vel_mL{c},1)),...
        'o','MarkerFaceColor','k','LineStyle','-','color','k','MarkerEdgeColor','none','MarkerFaceColor',colors(c,:),'capsize',0);

%     errorbar(...
%         nanmean(vel_mR{c}),...
%         nanmean(R_mR{c}),...
%         -nanstd(R_mR{c})/sqrt(size(R_mR{c},1)),...
%         nanstd(R_mR{c})/sqrt(size(R_mR{c},1)),...
%         -nanstd(vel_mR{c})/sqrt(size(vel_mR{c},1)),...
%         nanstd(vel_mR{c})/sqrt(size(vel_mR{c},1)),...
%         'o','MarkerFaceColor','k','LineStyle','-','color','k','MarkerEdgeColor','none','MarkerFaceColor',colors(c,:),'capsize',0);

end
