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
    
    lDirMean = nanmean(relMovResps(trialTypes.intVar.all.contrast_direction{5,1},whichCells));
    rDirMean = nanmean(relMovResps(trialTypes.intVar.all.contrast_direction{5,2},whichCells));
    
    lSideMean = nanmean(relStimResps(trialTypes.singleVar.side{1},whichCells));
    rSideMean = nanmean(relStimResps(trialTypes.singleVar.side{3},whichCells));
    
    cFeedMean = nanmean(relRewResps(trialTypes.intVar.all.contrast_outcome{5,1},whichCells));
    iFeedMean = nanmean(relRewResps(trialTypes.intVar.all.contrast_outcome{5,2},whichCells));
    
    hiMean = nanmean(pmovResps(correct_highSide,whichCells));
    loMean = nanmean(pmovResps(correct_lowSide,whichCells));
    
    impMean = nanmean(relMovResps(earlyTrials,whichCells));
    patMean = nanmean(relMovResps(lateTrials,whichCells));
    
    blMean = nanmean(pmovResps(trialTypes.singleVar.block{1},whichCells));
    brMean = nanmean(pmovResps(trialTypes.singleVar.block{2},whichCells));
    
     if expInfo.hemisphere > 0
        dir = (lDirMean - rDirMean);
        side = (lSideMean - rSideMean);
        block = (blMean - brMean);
    elseif expInfo.hemisphere < 0
        dir = (rDirMean - lDirMean);
        side = (rSideMean - lSideMean);
        block = (brMean - blMean);
    end
    feed = (cFeedMean - iFeedMean);
    value = (hiMean - loMean);
    timing = (impMean - patMean);
    
%     if expInfo.hemisphere < 0
%         dir = (lDirMean - rDirMean) ./ (lDirMean + rDirMean +eps);
%         side = (lSideMean - rSideMean) ./ (lSideMean + rSideMean +eps);
%         block = (blMean - brMean) ./ (blMean + brMean + eps);
%     elseif expInfo.hemisphere > 0
%         dir = (rDirMean - lDirMean) ./ (rDirMean + lDirMean +eps);
%         side = (rSideMean - lSideMean) ./ (rSideMean + lSideMean +eps);
%         block = (brMean - blMean) ./ (brMean + blMean + eps);
%     end
%     feed = (cFeedMean - iFeedMean) ./ (cFeedMean + iFeedMean +eps);
%     hilo = (hiMean - loMean) ./ (hiMean + loMean +eps);
%     timing = (impMean - patMean) ./ (impMean + patMean +eps);

%%
lim = max(max(abs(bigArray)))*1.1;
labels = {'Stimulus' 'Choice' 'Outcome' 'Block' 'Value'};
bigArray = [side; dir; feed; block; value];
figure;
set(gcf,'position',figpos)

for r = 1:size(bigArray,1)
    for c = 1:size(bigArray,1)
        subplot(size(bigArray,1),size(bigArray,1),sub2ind([size(bigArray,1) size(bigArray,1)],c,r))
        if r == c
            histogram(bigArray(c,:),linspace(-lim,lim,50),'FaceColor','k','LineStyle','none')
            prettyPlot(gca);
            xlim([-lim lim])
            title(labels{r})
            axis square
            cc(r,c) = 1;
        else 
            hold on;                 
            line([0 0],[-lim lim],'LineStyle',':','Color',[.5 .5 .5])
            line([-lim lim],[0 0],'LineStyle',':','Color',[.5 .5 .5])
            scatter(bigArray(r,:),bigArray(c,:),'k.');
            prettyPlot(gca);
            ylim([-lim lim])
            xlim([-lim lim])
            xlabel(labels{r});
            ylabel(labels{c});
            cc(r,c) = round(corr2(bigArray(r,:),bigArray(c,:))*1000)/1000;
            text(-lim*.9,lim*.9,strcat({'{\it r = }'},num2str(cc(r,c))))
            axis square
            if r ~= size(cc,1) || c ~= 1
            set(gca, 'XTickLabels', {''})
            set(gca, 'YTickLabels', {''})
            end
        end
    end
end

% for r = 1:size(bigArray,1)
%     for c = 1:size(bigArray,1)
%         if r == c
%             cc(r,c,m) = 1;
%         else 
%             cc(r,c,m) = round(corr2(bigArray(r,:),bigArray(c,:))*1000)/1000;
%         end
%     end
% end
% 
            
%%    
    
figure;
set(gcf,'position',figpos)

    
  
for r = 1:size(bigArray,1)
    for c = 1:size(bigArray,1)
        subplot(size(bigArray,1),size(bigArray,1),sub2ind([size(bigArray,1) size(bigArray,1)],c,r))
        if r == c
            prettyPlot(gca);
            xlim([-1 1])
            axis square
            axis off
        else 
            hold on;                 
            line([-2 2],[0 0],'LineStyle',':','Color',[.25 .25 .25])
%             histogram(cc(r,c,:),linspace(-.4,.4,15))
            violinPlotter(0, {squeeze(cc(r,c,:))},.3,.025,.003,'k','on')
%             beeswarm([ones(1,length(mouseList))'],[squeeze(cc(r,c,:))'],'sort_style','hex','dot_size',1,'overlay_style','sd');
            prettyPlot(gca);
            ylim([-.6 .6])
            xlim([-1 1])
            set(gca,'view',[90 90])
            mcc = round(nanmean(cc(r,c,:))*1000)/1000;
            scc = round(nanstd(cc(r,c,:))/sqrt(size(cc,3))*1000)/1000;
            [~,p]=ttest(cc(r,c,:));
            rp = round(p*1000)/1000;
            text(.9,-.55,strcat({'{\it r = }'},num2str(mcc),{'\pm'},num2str(scc)))
            if p > .05
                text(-.9,.5,'\it ns')
            elseif p <= 0.05 && p > 0.01
                text(-.9,.5,strcat({'{\it p = }'},num2str(p, '%0.2f')),'Color','r')
            elseif p <= 0.01 
                text(-.9,.5,strcat({'{\it p = }'},num2str(p, '%0.2e')),'Color','r')
            end
                ax = gca;
            ax.XColor = 'w';
            title(strcat(labels{c},{' vs. '},labels{r}))
            if r == 5
            ylabel('Corr. coefficient')
            end
        end
    end
end
    
    
    
    
    
    
    
    
    
    
    
%     
%     meanColor = [1 0 0];
%     ms = 4;
%     mm = 50;
%     ticksmax = .3;
%     % close all;
%     dirIndex = dir;
%     sideIndex = side;
%     rewIdx = feed;
%     timeIdx = timing;
%     hiloIdx = value;
%     [~,sortDirIdx] = sort(dirIndex);
%     [~,sortRewIdx] = sort(rewIdx);
%     
%     figure(1);
%     subplot(subDim,subDim,m)
%     lim = max([abs(sideIndex), abs(dirIndex)])*1.1;
%     line([0 0],[-lim lim],'LineStyle',':','Color',[.5 .5 .5])
%     hold on;
%     line([-lim lim],[0 0],'LineStyle',':','Color',[.5 .5 .5])
%     scatter(dirIndex,sideIndex,ms,...
%         'MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha', .7);
%     xlim([-lim lim]);
%     ylim([-lim lim]);
%     xticks([-ticksmax 0 ticksmax])
%     yticks([-ticksmax 0 ticksmax])
%     axis square
%     box off
%     set(gca,'tickdir','out')
%     ylabel('stim selectivity')
%     xlabel('choice selectivity')

% %% WORKBENCH
% figure;
% set(gcf,'position', [1000 1050 948 290])
% subplot(1,4,1)
% lim = max([abs(sideIndex), abs(dirIndex)])*1.1;
% line([0 0],[-lim lim],'LineStyle',':','Color',[.5 .5 .5])
% hold on;
% line([-lim lim],[0 0],'LineStyle',':','Color',[.5 .5 .5])
% scatter(dirIndex,sideIndex,ms,...
%     'MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha', .7);
% xlim([-lim lim]);
% ylim([-lim lim]);
% xticks([-.2 0 .2])
% yticks([-.2 0 .2])
% axis square
% box off
% set(gca,'tickdir','out')
% ylabel('stim selectivity')
% xlabel('choice selectivity')
% 
% subplot(1,4,2)
% lim = max([abs(rewIdx), abs(dirIndex)])*1.1;
% line([0 0],[-lim lim],'LineStyle',':','Color',[.5 .5 .5])
% hold on;
% line([-lim lim],[0 0],'LineStyle',':','Color',[.5 .5 .5])
% scatter(dirIndex,rewIdx,ms,...
%     'MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha', .7);
% xlim([-lim lim]);
% ylim([-lim lim]);
% xticks([-.2 0 .2])
% yticks([-.2 0 .2])
% axis square
% box off
% set(gca,'tickdir','out')
% ylabel('reward selectivity')
% xlabel('choice selectivity')
% 
% subplot(1,4,3)
% lim = max([abs(sideIndex), abs(rewIdx)])*1.1;
% line([0 0],[-lim lim],'LineStyle',':','Color',[.5 .5 .5])
% hold on;
% line([-lim lim],[0 0],'LineStyle',':','Color',[.5 .5 .5])
% scatter(sideIndex,rewIdx,ms,...
%     'MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha', .7);
% xlim([-lim lim]);
% ylim([-lim lim]);
% xticks([-.2 0 .2])
% yticks([-.2 0 .2])
% axis square
% box off
% set(gca,'tickdir','out')
% ylabel('reward selectivity')
% xlabel('stim selectivity')
% 
% subplot(1,4,4)
% lim = max([abs(rewIdx), abs(timeIdx)])*1.1;
% line([0 0],[-lim lim],'LineStyle',':','Color',[.5 .5 .5])
% hold on;
% line([-lim lim],[0 0],'LineStyle',':','Color',[.5 .5 .5])
% scatter(rewIdx,timeIdx,ms,...
%     'MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha', .7);
% xlim([-lim lim]);
% ylim([-lim lim]);
% xticks([-.2 0 .2])
% yticks([-.2 0 .2])
% axis square
% box off
% set(gca,'tickdir','out')
% ylabel('impulsivity selectivity')
% xlabel('reward selectivity')
% 
% fit{m}.dir_side = polyfit(dirIndex,sideIndex,1);
% fit{m}.dir_reward = polyfit(dirIndex,rewIdx,1);
% fit{m}.side_reward = polyfit(sideIndex,rewIdx,1);
% fit{m}.imp_reward = polyfit(timeIdx,rewIdx,1);
% 
% correlation{m}.dir_side = corr(dirIndex',sideIndex');
% correlation{m}.dir_reward = corr(dirIndex',rewIdx');
% correlation{m}.side_reward = corr(sideIndex',rewIdx');
% correlation{m}.imp_reward = corr(timeIdx',rewIdx');
%%
clearvars -except mouseList expList hemList subDim fit correlation cc
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
