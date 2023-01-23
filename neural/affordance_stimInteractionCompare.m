
for m = 1:length(mouseList)
    mouseName = char(mouseList{m});
    expDate = char(expList{m}{1});
    expNum = expList{m}{2};
    hemisphere = hemList(m);
    [expRef, ~] = data.constructExpRef(mouseName,expDate,expNum);
    fprintf('fetching %s...\n',expRef)
    [behavioralData, expInfo, neuralData] = data.loadDataset(mouseName, expDate, expNum);
    [baselineResps, stimResps, pmovResps, movResps, rewResps, preCueResps] = getEpochResps(neuralData.eta);

    %%
    ETA = 1;
    contrasts = getUniqueContrasts(expInfo);
    trialTypes = getTrialTypes(expInfo, behavioralData, 'late');
    
    trialsChoice = trialTypes.intVar.all.contrast_direction;
    trialsBlock = trialTypes.intVar.all.contrast_block;
    
%     trialsChoice = trialTypes.intVar.cb3D.contrast_direction;
%     trialsBlock = trialTypes.intVar.cb3D.contrast_block;

    if hemisphere > 0
        trialsChoice = trialsChoice;
        trialsBlock = trialsBlock;
    else
        trialsChoice = rot90(rot90(trialsChoice));
        trialsBlock = rot90(rot90(trialsBlock));
    end
    
    cellResps = neuralData.eta.alignedResps{ETA};
    
    whichResps = baselineResps;
    clear respsChoice respsBlock  
    for c = 1:length(contrasts)
        
        clear trialsChoiceA trialsChoiceB trialsBlockA trialsBlockB
        
        trialsChoiceA = trialsChoice{c,1};
        trialsChoiceB = trialsChoice{c,2};
        respsChoice(c,:) = squeeze(nanmean(whichResps(trialsChoiceA,:,:),1)) - squeeze(nanmean(whichResps(trialsChoiceB,:,:),1));
%         exampleRespChoice(c,1,:) = squeeze(nanmean(whichResps(trialsChoiceA,:,:),1));
%         exampleRespChoice(c,2,:) = squeeze(nanmean(whichResps(trialsChoiceB,:,:),1));
        trialsBlockA = trialsBlock{c,1};
        trialsBlockB = trialsBlock{c,2};
        respsBlock(c,:) = squeeze(nanmean(whichResps(trialsBlockA,:,:),1)) - squeeze(nanmean(whichResps(trialsBlockB,:,:),1));
%         exampleRespBlock(c,1,:) = squeeze(nanmean(whichResps(trialsBlockA,:,:),1));
%         exampleRespBlock(c,2,:) = squeeze(nanmean(whichResps(trialsBlockB,:,:),1));
    end
    
    
%     saveResps{m,:} = {exampleRespChoice exampleRespBlock};
    medianContrastChoiceDiff(m,:) = median(respsChoice,2);
    medianContrastBlockDiff(m,:) = median(respsBlock,2);
end


%%

plotColors = flipud([.25 0 0; .5 0 0;  1 0 0;  .8 .45 .45; .75 .75 .75; .6 .8 1; 0 .4 1; 0 0 1; 0 0 .5]);
figure;
subplot(1,2,1)
hold on
errorbar(contrasts, nanmedian(exampleRespChoice(:,1,:),3),nanstd(exampleRespChoice(:,1,:),[],3)/sqrt(length(exampleRespChoice)),...
    'Marker','none','LineStyle','-','MarkerSize',8,'MarkerFaceColor','k',...
    'MarkerEdgeColor','none','LineWidth',1,'Color','k','capsize',0)
for c = 1:length(contrasts)
    pp = plot(contrasts(c),nanmedian(exampleRespChoice(c,1,:)),...
        'MarkerSize',8,'Marker','o','MarkerFaceColor',plotColors(c,:),'MarkerEdgeColor','k');
    uistack(pp,'top');
end
hold on
errorbar(contrasts, nanmedian(exampleRespChoice(:,2,:),3),nanstd(exampleRespChoice(:,2,:),[],3)/sqrt(length(exampleRespChoice)),...
    'Marker','none','LineStyle','--','MarkerSize',8,'MarkerFaceColor','k',...
    'MarkerEdgeColor','none','LineWidth',1,'Color','k','capsize',0)
for c = 1:length(contrasts)
    pp = plot(contrasts(c),nanmedian(exampleRespChoice(c,2,:)),...
        'MarkerSize',8,'Marker','o','MarkerFaceColor',plotColors(c,:),'MarkerEdgeColor','k');
    uistack(pp,'top');
end
prettyPlot(gca)
xlim([-1.1 1.1])
ylabel('\DeltaNorm. activity')
xlabel('Ipsi contrast (%)   Contra contrast (%)')
xticks([-1 -.5 0 .5 1])
xticklabels({'-100', '-50', '0', '50', '100'})
set(gca, 'XDir','reverse')
title('Contra – ipsi choices')
ylim([.025 .065])

subplot(1,2,2)
hold on
errorbar(contrasts, nanmedian(exampleRespBlock(:,1,:),3),nanstd(exampleRespBlock(:,1,:),[],3)/sqrt(length(exampleRespBlock)),...
    'Marker','none','LineStyle','-','MarkerSize',8,'MarkerFaceColor','k',...
    'MarkerEdgeColor','none','LineWidth',1,'Color','k','capsize',0)
for c = 1:length(contrasts)
    pp = plot(contrasts(c),nanmedian(exampleRespBlock(c,1,:)),...
        'MarkerSize',8,'Marker','o','MarkerFaceColor',plotColors(c,:),'MarkerEdgeColor','k');
    uistack(pp,'top');
end
hold on
errorbar(contrasts, nanmedian(exampleRespBlock(:,2,:),3),nanstd(exampleRespBlock(:,2,:),[],3)/sqrt(length(exampleRespBlock)),...
    'Marker','none','LineStyle','--','MarkerSize',8,'MarkerFaceColor','k',...
    'MarkerEdgeColor','none','LineWidth',1,'Color','k','capsize',0)
for c = 1:length(contrasts)
    pp = plot(contrasts(c),nanmedian(exampleRespBlock(c,2,:)),...
        'MarkerSize',8,'Marker','o','MarkerFaceColor',plotColors(c,:),'MarkerEdgeColor','k');
    uistack(pp,'top');
end
prettyPlot(gca)
xlim([-1.1 1.1])
ylabel('Norm. activity')
xlabel('Ipsi contrast (%)   Contra contrast (%)')
xticks([-1 -.5 0 .5 1])
xticklabels({'-100', '-50', '0', '50', '100'})
set(gca, 'XDir','reverse')
title('Contra – ipsi choices')
ylim([.025 .065])


%%
plotColors = flipud([.25 0 0; .5 0 0;  1 0 0;  .8 .45 .45; .75 .75 .75; .6 .8 1; 0 .4 1; 0 0 1; 0 0 .5]);
figure;
subplot(1,2,1)
hold on
line([-1.1 1.1],[0 0],'Color',[.5 .5 .5],'LineStyle','--')
errorbar(contrasts, nanmean(medianContrastChoiceDiff,1),nanstd(medianContrastChoiceDiff,[],1),...
    'Marker','none','LineStyle','-','MarkerSize',8,'MarkerFaceColor','k','MarkerEdgeColor','none','LineWidth',1,'Color','k','capsize',0)
for c = 1:length(contrasts)
    xs = contrasts(c)+(rand(1,length(medianContrastChoiceDiff))-.5)*.08;
    ss = scatter(xs, medianContrastChoiceDiff(:,c),...
        'MarkerEdgeColor','none','MarkerFaceColor',plotColors(c,:),'MarkerFaceAlpha',.3);
    uistack(ss,'bottom');
    pp = plot(contrasts(c),nanmean(medianContrastChoiceDiff(:,c)),...
        'MarkerSize',8,'Marker','o','MarkerFaceColor',plotColors(c,:),'MarkerEdgeColor','k');
    uistack(pp,'top');
end
prettyPlot(gca)
xlim([-1.1 1.1])
ylabel('\DeltaNorm. activity')
xlabel('Ipsi contrast (%)   Contra contrast (%)')
xticks([-1 -.5 0 .5 1])
xticklabels({'-100', '-50', '0', '50', '100'})
set(gca, 'XDir','reverse')
title('Contra – ipsi choices')
ylim([-.06 .06])

subplot(1,2,2)
hold on
line([-1.1 1.1],[0 0],'Color',[.5 .5 .5],'LineStyle','--')
errorbar(contrasts, nanmean(medianContrastBlockDiff,1),nanstd(medianContrastBlockDiff,[],1),...
    'Marker','none','LineStyle','-','MarkerSize',8,'MarkerFaceColor','k','MarkerEdgeColor','none','LineWidth',1,'Color','k','capsize',0)
for c = 1:length(contrasts)
    xs = contrasts(c)+(rand(1,length(medianContrastBlockDiff))-.5)*.08;
    ss = scatter(xs, medianContrastBlockDiff(:,c),...
        'MarkerEdgeColor','none','MarkerFaceColor',plotColors(c,:),'MarkerFaceAlpha',.3);
    uistack(ss,'bottom');
    pp = plot(contrasts(c),nanmean(medianContrastBlockDiff(:,c)),...
        'MarkerSize',8,'Marker','o','MarkerFaceColor',plotColors(c,:),'MarkerEdgeColor','k');
    uistack(pp,'top');
end
prettyPlot(gca)
xlim([-1.1 1.1])
% ylabel('\DeltaNorm. activity')
xlabel('Ipsi contrast (%)   Contra contrast (%)')
xticks([-1 -.5 0 .5 1])
xticklabels({'-100', '-50', '0', '50', '100'})
set(gca, 'XDir','reverse')
title('Contra – ipsi blocks')
ylim([-.06 .06])

    %% split into quintiles

    [~, sidx] = sort(nanmean(stimResps(side{1},:)) - nanmean(stimResps(side{3},:)),'descend');
    nc = size(baselineResps,2);
    pctChunks = 5;
    pctLength = floor(nc/pctChunks);

    for pc = 1:pctChunks
        if pc == pctChunks
            whichNeurons{pc} = sidx((pc-1)*pctLength+1:end);
        else
            whichNeurons{pc} = sidx((pc-1)*pctLength+1:pc*pctLength);
        end
    end

    %%

    %%
    colors = [0 .4 1; .6 .8 1; .75 .75 .75; .8 .45 .45; 1 0 0];
    subplot(5,7,m)
   
    for pwn = 1:length(whichNeurons)
    clear sx sy
        for t = 1:41
            sx(t,:) = nanmean(stimIdx_mC(t,whichNeurons{pwn}) - stimIdx_mI(t,whichNeurons{pwn}));
            sy(t,:) = nanmean(stimIdx_bC(t,whichNeurons{pwn}) - stimIdx_bI(t,whichNeurons{pwn}));
        end
        sx = smooth(sx);
        sy = smooth(sy);
        plot((sx),(sy),'Color',colors(pwn,:),'LineWidth',1)
        hold on;
        plot(sx(1),sy(1),'o','MarkerFaceColor','none','MarkerEdgeColor',colors(pwn,:),'LineWidth',1)
        plot(sx(21),sy(21),'ko','MarkerFaceColor',colors(pwn,:))
        plot(sx(41),sy(41),'o','MarkerEdgeColor','none','MarkerFaceColor',colors(pwn,:))
        
        ma(pwn) = max(abs([sx; sy]));
    end
    
    maxax = max(ma);
     line([-maxax maxax],[0 0],'LineStyle','--','Color',[0.5 0.5 0.5])
    hold on
    line([0 0],[-maxax maxax],'LineStyle','--','Color',[0.5 0.5 0.5])
    try
        xlim([-maxax maxax])
        ylim([-maxax maxax])
    catch
    end
    prettyPlot(gca);
    axis square
    title(expRef,'Interpreter','none')
    
%     clearvars -except mouseList expList hemList