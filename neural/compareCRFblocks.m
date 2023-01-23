for m = 1%:length(mouseList)
    mouseName = char(mouseList{m});
    expDate = char(expList{m}{1});
    expNum = expList{m}{2};
    hemisphere = hemList(m);
    [expRef, ~] = data.constructExpRef(mouseName,expDate,expNum);
    fprintf('fetching %s...\n',expRef)
    [behavioralData, expInfo, neuralData] = data.loadDataset(mouseName, expDate, expNum);
    [baselineResps, stimResps, pmovResps, movResps, rewResps, preCueResps] = getEpochResps(neuralData.eta);
    
    %% SET PARAMS
    trimLength = 10;
    nt = length(behavioralData.eventTimes(1).daqTime);
    np = 1000;
    nc = size(baselineResps,2);
    pctChunks = 5;
    pctLength = floor(nc/pctChunks);
    
    [whichTrialsLogical, whichTrials] = selectCondition(expInfo, getUniqueContrasts(expInfo), behavioralData, ...
        initTrialConditions('responseType','all','movementTime','late','specificRTs',[0.8 2]));
    
     %% EXTRACT TRIAL VARIABLES AND GENERATE PSEUDOBLOCKS
    
    % correct for hemisphere
    if hemisphere < 0
        trueStimuli = expInfo.block.events.contrastValues(1:nt);
        trueCorrectChoice = expInfo.block.events.correctResponseValues(1:nt);
        trueBlocks = expInfo.block.events.highRewardSideValues(1:nt);
        trueActualChoice = expInfo.block.events.responseValues(1:nt);
    else
        trueStimuli = -expInfo.block.events.contrastValues(1:nt);
        trueCorrectChoice = -expInfo.block.events.correctResponseValues(1:nt);
        trueBlocks = -expInfo.block.events.highRewardSideValues(1:nt);
        trueActualChoice = -expInfo.block.events.responseValues(1:nt);
    end
    trueFeedback = expInfo.block.events.feedbackValues(1:nt);
    
    % remove the first trials after the block switch by overwriting nans
    idxTrim = false(1,nt);
    switchPoints = find(diff(expInfo.block.events.highRewardSideValues(1:nt)) ~= 0);
    if numel(switchPoints) == 1
        idxTrim(trimLength+1:switchPoints) = true;
        idxTrim(switchPoints+trimLength+1:nt) = true;
    else
        for s = 1:length(switchPoints)
            if s == 1
                idxTrim(trimLength+1:switchPoints(s)) = true;
            elseif s < length(switchPoints)
                idxTrim(switchPoints(s-1)+trimLength+1:switchPoints(s)) = true;
            elseif s == length(switchPoints)
                idxTrim(switchPoints(s-1)+trimLength+1:switchPoints(s)) = true;
                idxTrim(switchPoints(s)+trimLength+1:nt) = true;
            end
        end
    end
    trueBlocks(~idxTrim) = nan;
    
%     % assign the 0% stimuli as either 'left' or 'right' depending on where
%     % the reward was preassigned (contra or ipsi)
%     trueStimuli(trueStimuli == 0) = eps;
%     trueStimuli(abs(trueStimuli) < .05) = ...
%         trueStimuli(abs(trueStimuli) < .05).* trueCorrectChoice(abs(trueStimuli) < .05);

    % low rewards are possible on sign-mismatched block and stimulus
    % high rewards are possible on sign-matched block and stimulus
    % 2 = high, 1 = low, 0 = error/none, nan = omitted trim trials
    trueValue(trueBlocks.*sign(trueStimuli) == -1) = 1;
    trueValue(trueBlocks.*sign(trueStimuli) == 1) = 2;
    trueValue(trueFeedback==0) = 0;
    trueValue(~idxTrim) = nan;

    % generate pseudoblocks
    blockStart = 'fixed';
    for p = 1:np
        if strcmp(blockStart,'fixed')
            firstSide = trueBlocks(trimLength+1);
        elseif strcmp(blockStart,'rand')
            firstSide = randsample([-1, 1],1,true);
        end
        b=nan(1,nt);
        switches = cumsum(125+randi(100,1,20));
        for s = 1:length(switches)
            if s == 1
                b((1+trimLength):switches(s)-1) = firstSide;
            elseif mod(s,2) == 1
                b((switches(s-1)+trimLength):switches(s)-1) = firstSide;
            elseif mod(s,2) == 0
                b((switches(s-1)+trimLength):switches(s)-1) = -firstSide;
            end
        end
        pseudoBlocks(:,p) = b(1:nt);
    end

    % assign new high/low trials based on the pseudoblocks + true stim
    for p = 1:np
        pseudoValue(pseudoBlocks(:,p).*sign(trueStimuli)' == -1,p) = 1;
        pseudoValue(pseudoBlocks(:,p).*sign(trueStimuli)' == 1,p) = 2;
        pseudoValue(trueFeedback==0,p) = 0;
        pseudoValue(isnan(pseudoBlocks(:,p)),p) = nan;
    end
    
    %% compute true blockDiff per neuron
    tic
    for c = 1:length(contrasts)
        plotTrials = intersect(find(trueStimuli == contrasts(c)),whichTrials);
        scbt = intersect(find(trueBlocks > 0),plotTrials);
        sibt = intersect(find(trueBlocks < 0),plotTrials);

%             contrastResps(:,:,c) = squeeze(nanmean(neuralData.eta.alignedResps{1}(plotTrials,:,:),1));
        bd(c,:) = nanmean(stimResps(scbt,:),1) - nanmean(stimResps(sibt,:),1);

    end
    trueBlockDiff = sum(bd,1)';
    
    for p = 1:np
        for c = 1:length(contrasts)
            plotTrials = intersect(find(trueStimuli == contrasts(c)),whichTrials);
            scbt = intersect(find(pseudoBlocks(:,p) > 0),plotTrials);
            sibt = intersect(find(pseudoBlocks(:,p) < 0),plotTrials);
            pbd(c,:,p) = nanmean(stimResps(scbt,:),1) - nanmean(stimResps(sibt,:),1);
        end
    end
    pseudoBlockDiff = squeeze(sum(pbd,1));
    toc
    %% compute stimulus selectivity index & sort neurons into ntiles
%     allResps(:,:,1) = baselineResps;
%     allResps(:,:,2) = stimResps;
%     allResps(:,:,3) = movResps;
%     allResps(:,:,4) = rewResps;
%     
    allResps(:,:,1) = stimResps;
    
%     cst = intersect(find(trueActualChoice > 0),whichTrials);
%     ist = intersect(find(trueActualChoice < 0),whichTrials);
    
    cst = intersect(find(trueStimuli > 0),whichTrials);
    ist = intersect(find(trueStimuli < 0),whichTrials);
    
    %%
   
    figure;
    set(gcf,'position',[1222 118 1722 1486])
    for r = 1:size(allResps,3)

        whichResps = allResps(:,:,r);
        contraResp = nanmean(stimResps(cst,:),1);
        ipsiResp = nanmean(stimResps(ist,:),1);

        stimSelectivityIdx = (contraResp - ipsiResp) ./ (contraResp + ipsiResp + eps);
        [~, sortIdx] = sort(stimSelectivityIdx);    
        for pc = 1:pctChunks
            if pc == pctChunks
                whichNeurons{pc} = sortIdx((pc-1)*pctLength+1:end);
            else
                whichNeurons{pc} = sortIdx((pc-1)*pctLength+1:pc*pctLength);
            end
        end

        whichNeurons = fliplr(whichNeurons);

        %% calculate crf (true)
        contrasts = getUniqueContrasts(expInfo);
        for cc = 1:nc  
            for c = 1:length(contrasts)
                plotTrials = find((trueStimuli == contrasts(c)) .* (whichTrialsLogical));
                scbt = find((trueStimuli == contrasts(c)) .* (trueBlocks > 0) .* (whichTrialsLogical));
                sibt = find((trueStimuli == contrasts(c)) .* (trueBlocks < 0) .* (whichTrialsLogical));

                pcrfcm(c) = nanmean(whichResps(scbt,cc),1);
                pcrfcs(c) = nanstd(whichResps(plotTrials,cc),1)/sqrt(length(scbt));
                pcrfim(c) = nanmean(whichResps(sibt,cc),1);
                pcrfis(c) = nanstd(whichResps(plotTrials,cc),1)/sqrt(length(sibt));

                dpcrf(cc,c) = pcrfcm(c) - pcrfim(c);

            end
        end

        for cc = 1:nc
            diffFits(cc,:) = polyfit(contrasts(~isnan(dpcrf(cc,:))), dpcrf(cc,~isnan(dpcrf(cc,:))),1);
        end

        meanBlockDiff = nanmean(dpcrf,2);
        
        
        
%         %% calculate crf (pseudo)
%         for cc = 1:nc  
%             for c = 1:length(contrasts)
%                 plotTrials = find((trueStimuli == contrasts(c)) .* (whichTrialsLogical));
%                 for p = 1:np
%                     scbt = find((trueStimuli == contrasts(c)) .* (pseudoBlocks(:,p)' > 0) .* (whichTrialsLogical));
%                     sibt = find((trueStimuli == contrasts(c)) .* (pseudoBlocks(:,p)' < 0) .* (whichTrialsLogical));
%                     pcrfcmp = nanmean(whichResps(scbt,cc),1);
%                     pcrfcsp = nanstd(whichResps(plotTrials,cc),1)/sqrt(length(scbt));
%                     pcrfimp = nanmean(whichResps(sibt,cc),1);
%                     pcrfisp = nanstd(whichResps(plotTrials,cc),1)/sqrt(length(sibt));
% 
%                     dpcrfp(cc,c,p) = pcrfcmp - pcrfimp;
% 
%                 end
%             end
%         end
% 
%         meanBlockDiff_pseudo = squeeze(nanmean(dpcrfp,2));
%         %%
%         %movmean
%         [~,sidx] = sort(stimSelectivityIdx);
%         movmeanDiff_true = movmean(meanBlockDiff(sidx),50);
%         movmeanDiff_pseudo = movmean(meanBlockDiff_pseudo(sidx,:),50,1);
%         
%         figure;
%         subplot(5,5,1)
%         maxy = 1.1*max(abs(meanBlockDiff));
%         maxx = 1.1*max(abs(stimSelectivityIdx));
%         hold on;
%         line([-maxx maxx],[0 0],'LineStyle','--','Color',[.5 .5 .5]);
%         line([0 0],[-maxy maxy],'LineStyle','--','Color',[.5 .5 .5]);
%         scatter(stimSelectivityIdx,meanBlockDiff,30,'MarkerFaceAlpha',.2,'MarkerEdgeColor','none','MarkerFaceColor','k','Marker','o');
%         xlim([-maxx maxx]);
%         ylim([-maxy maxy]);
%         plot(stimSelectivityIdx(sidx),movmeanDiff_true,'g')
%         xlabel('Stimulus selectivity index')
%         ylabel('Mean block difference')
% %         str = {char(strcat({'slope: '},num2str(allSlope_true(m)))),char(strcat({'int: '},num2str(allInt_true(m)))),char(strcat({'rho: '},num2str(allRho_true(m))))};
% %         text(-maxx*.9,maxy*.8,str)
%         box off
%         set(gca,'tickdir','out')
%         title('True')
%         
%         for p = 1:24
%             subplot(5,5,p+1)
%             maxy = 1.1*max(abs(meanBlockDiff_pseudo(:,p)));
%             maxx = 1.1*max(abs(stimSelectivityIdx));
%             hold on;
%             line([-maxx maxx],[0 0],'LineStyle','--','Color',[.5 .5 .5]);
%             line([0 0],[-maxy maxy],'LineStyle','--','Color',[.5 .5 .5]);
%             scatter(stimSelectivityIdx,meanBlockDiff_pseudo(:,p),30,'MarkerFaceAlpha',.2,'MarkerEdgeColor','none','MarkerFaceColor','k','Marker','o');
%             xlim([-maxx maxx]);
%             ylim([-maxy maxy]);
%             plot(stimSelectivityIdx(sidx),movmeanDiff_pseudo(:,p),'r')
%             xlabel('Stimulus selectivity index')
%             ylabel('Mean block difference')
%         end
 
        %% relationship between stim selectivity and mean block difference

        % true
        pf = polyfit(stimSelectivityIdx,meanBlockDiff,1);
        px = linspace(-1,1,10);
        py = px*pf(1) + pf(2);
        [rho,~] = corr(stimSelectivityIdx',meanBlockDiff,'Type','Spearman');

        allRho_true(m) = rho;
        allSlope_true(m) = pf(1);
        allInt_true(m) = pf(2);

        %pseudo
        % true
        for p = 1:np
            pf = polyfit(stimSelectivityIdx,meanBlockDiff_pseudo(:,p),1);
            [rho,~] = corr(stimSelectivityIdx',meanBlockDiff_pseudo(:,p),'Type','Spearman');

            allRho_pseudo(m,p) = rho;
            allSlope_pseudo(m,p) = pf(1);
            allInt_pseudo(m,p) = pf(2);
        end
    
        rank.slope = tiedrank(allSlope_pseudo(m,:));
        rank.int = tiedrank(allInt_pseudo(m,:));
        rank.rho = tiedrank(allRho_pseudo(m,:));

        [~, minIdx] = min(abs(allSlope_pseudo(m,:)-allSlope_true(m)));
        sigrank.slope = rank.slope(minIdx)/np;
        [~, minIdx] = min(abs(allInt_pseudo(m,:)-allInt_true(m)));
        sigrank.int = rank.int(minIdx)/np;
        [~, minIdx] = min(abs(allRho_pseudo(m,:)-allRho_true(m)));
        sigrank.rho = rank.rho(minIdx)/np;

        %% plot scatter of stim selectivity index and mean diffCRF between blocks

        
        subplot(size(allResps,3),4,(r-1)*size(allResps,3)+1)
        maxy = 1.1*max(abs(meanBlockDiff));
        maxx = 1.1*max(abs(stimSelectivityIdx));
        hold on;
        line([-maxx maxx],[0 0],'LineStyle','--','Color',[.5 .5 .5]);
        line([0 0],[-maxy maxy],'LineStyle','--','Color',[.5 .5 .5]);
        scatter(stimSelectivityIdx,meanBlockDiff,30,'MarkerFaceAlpha',.2,'MarkerEdgeColor','none','MarkerFaceColor','k','Marker','o');
        xlim([-maxx maxx]);
        ylim([-maxy maxy]);
        plot(px,py,'r')
        xlabel('Stimulus selectivity index')
        ylabel('Mean block difference')
        str = {char(strcat({'slope: '},num2str(allSlope_true(m)))),char(strcat({'int: '},num2str(allInt_true(m)))),char(strcat({'rho: '},num2str(allRho_true(m))))};
        text(-maxx*.9,maxy*.8,str)
        box off
        set(gca,'tickdir','out')
        title('Stim selectivity vs mean block difference')

        subplot(size(allResps,3),4,(r-1)*size(allResps,3)+2)
        h = histogram(allSlope_pseudo(m,:),41,'FaceColor',[.5 .5 .5]);
        hold on;
        line([allSlope_true(m) allSlope_true(m)],[0 max(h.Values)],'LineStyle','-','LineWidth',2,'Color','k');
        while 10 < 100
            if sigrank.slope >= 0.999 || sigrank.slope < 0.001
                text(allSlope_true(m),max(h.Values*1.05),'***','HorizontalAlignment','center')
                break
            elseif sigrank.slope >= 0.990 || sigrank.slope < 0.010
                text(allISlope_true(m),max(h.Values*1.05),'**','HorizontalAlignment','center')
                break
            elseif sigrank.slope >= 0.975 || sigrank.slope < 0.025
                text(allSlope_true(m),max(h.Values*1.05),'*','HorizontalAlignment','center')
                break
            else
                text(allSlope_true(m),max(h.Values*1.05),'ns','HorizontalAlignment','center')
                break
            end
        end
        box off
        ylim([0 max(h.Values)*1.1])
        set(gca,'tickdir','out')
        ylabel('Frequency')
        xlabel('Slope value')
        title('True slope vs pseudosessions')

        subplot(size(allResps,3),4,(r-1)*size(allResps,3)+3)
        h = histogram(allInt_pseudo(m,:),41,'FaceColor',[.5 .5 .5]);
        hold on;
        line([allInt_true(m) allInt_true(m)],[0 max(h.Values)],'LineStyle','-','LineWidth',2,'Color','k');
        while 10 < 100
            if sigrank.int >= 0.999 || sigrank.int < 0.001
                text(allInt_true(m),max(h.Values*1.05),'***','HorizontalAlignment','center')
                break
            elseif sigrank.int >= 0.990 || sigrank.int < 0.010
                text(allInt_true(m),max(h.Values*1.05),'**','HorizontalAlignment','center')
                break
            elseif sigrank.int >= 0.975 || sigrank.int < 0.025
                text(allInt_true(m),max(h.Values*1.05),'*','HorizontalAlignment','center')
                break
            else
                text(allInt_true(m),max(h.Values*1.05),'ns','HorizontalAlignment','center')
                break
            end
        end
        ylim([0 max(h.Values)*1.1])    
        box off
        set(gca,'tickdir','out')
        ylabel('Frequency')
        xlabel('Int value')
        title('True int vs pseudosessions')

        subplot(size(allResps,3),4,(r-1)*size(allResps,3)+4)
        h = histogram(allRho_pseudo(m,:),41,'FaceColor',[.5 .5 .5]);
        hold on;
        line([allRho_true(m) allRho_true(m)],[0 max(h.Values)],'LineStyle','-','LineWidth',2,'Color','k');
        while 10 < 100
            if sigrank.rho >= 0.999 || sigrank.rho < 0.001
                text(allRho_true(m),max(h.Values*1.05),'***','HorizontalAlignment','center')
                break
            elseif sigrank.rho >= 0.990 || sigrank.rho < 0.010
                text(allRho_true(m),max(h.Values*1.05),'**','HorizontalAlignment','center')
                break
            elseif sigrank.rho >= 0.975 || sigrank.rho < 0.025
                text(allRho_true(m),max(h.Values*1.05),'*','HorizontalAlignment','center')
                break
            else
                text(allRho_true(m),max(h.Values*1.05),'ns','HorizontalAlignment','center')
                break
            end
        end
        ylim([0 max(h.Values)*1.1])
        box off
        set(gca,'tickdir','out')
        ylabel('Frequency')
        xlabel('Rho value')
        title('True rho vs pseudosessions')
    end
    
    %% plot
    
    plotColors = [.25 0 0; .5 0 0;  1 0 0;  .8 .45 .45; .75 .75 .75; .6 .8 1; 0 .4 1; 0 0 1; 0 0 .5];
    
    figure;
    set(gcf,'position',[80 1360 1320 212])
    hold on;
    for c = 1:length(contrasts)
        plotTrials = intersect(find(trueStimuli == contrasts(c)),whichTrials);
        scbt = intersect(intersect(find(trueStimuli == contrasts(c)), find(trueBlocks > 0)),whichTrials);
        sibt = intersect(intersect(find(trueStimuli == contrasts(c)), find(trueBlocks < 0)),whichTrials);
        
        subplot(1,5,1)
        hold on;
        plot(neuralData.eta.eventWindow,nanmean(neuralData.eta.alignedResps{1}(plotTrials,:,sortIdx(end)),1),'LineWidth',2,'Color',plotColors(c,:));
        maxy1(c) = max(max(nanmean(neuralData.eta.alignedResps{1}(plotTrials,15:36,sortIdx(end)),1)));
        miny1(c) = min(min(nanmean(neuralData.eta.alignedResps{1}(plotTrials,15:36,sortIdx(end)),1)));
        
        subplot(1,5,2)
        hold on;
        pcrfm(c) = nanmean(stimResps(plotTrials,sortIdx(end)),1);
        pcrfs(c) = nanstd(stimResps(plotTrials,sortIdx(end)),1)/sqrt(length(plotTrials));
        plot(contrasts(c),pcrfm(c),'o','MarkerFaceColor',plotColors(c,:),'MarkerSize',8,'MarkerEdgeColor','none')
        maxy2(c) = max(pcrfm(c)+pcrfs(c));
        miny2(c) = min(pcrfm(c)-pcrfs(c));
        
        subplot(1,5,3)
        hold on;
        pcrfcm(c) = nanmean(stimResps(scbt,sortIdx(end)),1);
        pcrfcs(c) = nanstd(stimResps(plotTrials,sortIdx(end)),1)/sqrt(length(scbt));
        pcrfim(c) = nanmean(stimResps(sibt,sortIdx(end)),1);
        pcrfis(c) = nanstd(stimResps(plotTrials,sortIdx(end)),1)/sqrt(length(sibt));
        plot(contrasts(c),pcrfcm(c),'o','MarkerFaceColor',plotColors(c,:),'MarkerSize',8,'MarkerEdgeColor','none')
        plot(contrasts(c),pcrfim(c),'o','MarkerFaceColor',plotColors(c,:),'MarkerSize',8,'MarkerEdgeColor','none')
        maxy3(c) = max([pcrfcm(c)+pcrfcs(c) pcrfim(c)+pcrfis(c)]);
        miny3(c) = min([pcrfcm(c)-pcrfcs(c) pcrfim(c)-pcrfis(c)]);
        
        subplot(1,5,4)
        hold on;
        plot(contrasts(c),pcrfcm(c)-pcrfim(c),'o','MarkerFaceColor',plotColors(c,:),'MarkerSize',8,'MarkerEdgeColor','none')
        maxy4(c) = max(abs(pcrfcm(c)-pcrfim(c)));
    end
    
    
    subplot(1,5,1)
    hold on;
    ll = line([0 0],[0 5],'LineStyle','--','Color',[.5 .5 .5]);
    uistack(ll,'bottom') 
    xlim([-.5 1.5])
    ylim([min(miny1)-.1*abs(min(miny1)) max(maxy1)+.1*abs(max(maxy1))])
    set(gca,'tickdir','out')
    xlabel('Time from stim on (s)')
    title('Cell response')
    
    subplot(1,5,2)
    hold on;
    pp = plot(contrasts,pcrfm,'k');
    uistack(pp,'bottom')
    for c = 1:length(contrasts)
        ls = line([contrasts(c) contrasts(c)],[pcrfm(c)-pcrfs(c) pcrfm(c)+pcrfs(c)],'Color','k');
        uistack(ls,'bottom')
    end
    xlim([-1.1 1.1])
    ylim([min(miny2)-.1*abs(min(miny2)) max(maxy2)+.1*abs(max(maxy2))])
    set(gca,'tickdir','out')
    xlabel('Contrast')
    title('Cell CRF, all trials')
    
    subplot(1,5,3)
    hold on;
    pp = plot(contrasts,pcrfcm,'k');
    uistack(pp,'bottom')
    pp = plot(contrasts,pcrfim,'k--');
    uistack(pp,'bottom')
    for c = 1:length(contrasts)
        lsc = line([contrasts(c) contrasts(c)],[pcrfcm(c)-pcrfcs(c) pcrfcm(c)+pcrfcs(c)],'Color','k');
        uistack(lsc,'bottom')
        lsc = line([contrasts(c) contrasts(c)],[pcrfim(c)-pcrfis(c) pcrfim(c)+pcrfis(c)],'Color','k');
        uistack(lsc,'bottom')
    end
    xlim([-1.1 1.1])
    ylim([min(miny3)-.1*abs(min(miny3)) max(maxy3)+.1*abs(max(maxy3))])
    set(gca,'tickdir','out')
    xlabel('Contrast')
    title('Cell CRF, between blocks')
    
    subplot(1,5,4)
    hold on;
    pp = plot(contrasts,pcrfcm-pcrfim,'k');
    uistack(pp,'bottom')
    xlim([-1.1 1.1])
    ylim([-max(maxy4)*1.1 max(maxy4)*1.1])
    set(gca,'tickdir','out')
    xlabel('Contrast')
    title('Cell CRF, difference')
    l0 = line([-1.1 1.1],[0 0],'LineStyle','--','Color',[.5 .5 .5]);
    uistack(l0,'bottom') 
    
    subplot(1,5,5)
    scatter(diffFits(sortIdx(end),1),diffFits(sortIdx(end),2),30,'MarkerFaceAlpha',1,'MarkerEdgeColor','none','MarkerFaceColor','k','Marker','o');
    l1 = line([-1 1],[0 0],'LineStyle','--','Color',[.5 .5 .5]);
    uistack(l1,'bottom')
    l2 = line([0 0],[-1 1],'LineStyle','--','Color',[.5 .5 .5]);
    uistack(l2,'bottom')
    xlabel('slope')
    ylabel('int')
    xlim([-max(max(abs(diffFits))) max(max(abs(diffFits)))])
    ylim([-max(max(abs(diffFits))) max(max(abs(diffFits)))])
    title('Fit params')
    
    printfig(gcf,char(strcat(expRef,{' bCRFs - top contra cell'})))
    
    %% plot means by quintile
    subw = 5;
    figure;
    set(gcf,'position',[80 250 1320 1360])
    nw = length(whichNeurons);
    for ww = 1:nw
        for c = 1:length(contrasts)
            plotTrials = intersect(find(trueStimuli == contrasts(c)),whichTrials);
            scbt = intersect(intersect(find(trueStimuli == contrasts(c)), find(trueBlocks > 0)),whichTrials);
            sibt = intersect(intersect(find(trueStimuli == contrasts(c)), find(trueBlocks < 0)),whichTrials);

            subplot(nw,subw,(ww-1)*subw+1)
            hold on;
            plot(neuralData.eta.eventWindow,nanmean(nanmean(neuralData.eta.alignedResps{1}(plotTrials,:,whichNeurons{ww}),1),3),'LineWidth',2,'Color',plotColors(c,:));
            maxy1(ww,c) = max(max(nanmean(nanmean(neuralData.eta.alignedResps{1}(plotTrials,:,whichNeurons{ww}),1),3)));
            miny1(ww,c) = min(min(nanmean(nanmean(neuralData.eta.alignedResps{1}(plotTrials,:,whichNeurons{ww}),1),3)));
            
            subplot(nw,subw,(ww-1)*subw+2)
            hold on;
            pcrfm(c) = nanmean(nanmean(stimResps(plotTrials,whichNeurons{ww}),1));
            pcrfs(c) = nanstd(nanmean(stimResps(plotTrials,whichNeurons{ww}),1))/sqrt(length(plotTrials));
            plot(contrasts(c),pcrfm(c),'o','MarkerFaceColor',plotColors(c,:),'MarkerSize',8,'MarkerEdgeColor','none')
            maxy2(ww,c) = max(pcrfm(c)+pcrfs(c));
            miny2(ww,c) = min(pcrfm(c)-pcrfs(c));

            subplot(nw,subw,(ww-1)*subw+3)
            hold on;
            pcrfcm(c) = nanmean(nanmean(stimResps(scbt,whichNeurons{ww}),1));
            pcrfcs(c) = nanstd(nanmean(stimResps(scbt,whichNeurons{ww}),1))/sqrt(length(whichNeurons{ww}));
            pcrfim(c) = nanmean(nanmean(stimResps(sibt,whichNeurons{ww}),1));
            pcrfis(c) = nanstd(nanmean(stimResps(sibt,whichNeurons{ww}),1))/sqrt(length(whichNeurons{ww}));
            plot(contrasts(c),pcrfcm(c),'o','MarkerFaceColor',plotColors(c,:),'MarkerSize',8,'MarkerEdgeColor','none')
            plot(contrasts(c),pcrfim(c),'o','MarkerFaceColor',plotColors(c,:),'MarkerSize',8,'MarkerEdgeColor','none')
            maxy3(ww,c) = max([pcrfcm(c)+pcrfcs(c) pcrfim(c)+pcrfis(c)]);
            miny3(ww,c) = min([pcrfcm(c)-pcrfcs(c) pcrfim(c)-pcrfis(c)]);
            
            subplot(nw,subw,(ww-1)*subw+4)
            hold on;
            plot(contrasts(c),pcrfcm(c)-pcrfim(c),'o','MarkerFaceColor',plotColors(c,:),'MarkerSize',8,'MarkerEdgeColor','none')
            maxy4(ww,c) = max(abs(pcrfcm(c)-pcrfim(c)));            

        end
        
       

        subplot(nw,subw,(ww-1)*subw+1)
        hold on;
        ll = line([0 0],[0 5],'LineStyle','--','Color',[.5 .5 .5]);
        uistack(ll,'bottom') 
        xlim([-.5 1.5])
%         ylim([bigMin1-.1*bigMin1 bigMax1+.1*bigMax1])
        set(gca,'tickdir','out')
        xlabel('Time from stim on (s)')
        title('Cell response')

        subplot(nw,subw,(ww-1)*subw+2)
        hold on;
        pp = plot(contrasts,pcrfm,'k');
        uistack(pp,'bottom')
        for c = 1:length(contrasts)
            ls = line([contrasts(c) contrasts(c)],[pcrfm(c)-pcrfs(c) pcrfm(c)+pcrfs(c)],'Color','k');
            uistack(ls,'bottom')
        end
        xlim([-1.1 1.1])
%         ylim([bigMin2-.1*bigMin2 bigMax2+.1*bigMax2])
        set(gca,'tickdir','out')
        xlabel('Contrast')
        title('Cell CRF, all trials')

        subplot(nw,subw,(ww-1)*subw+3)
        hold on;
        pp = plot(contrasts,pcrfcm,'k');
        uistack(pp,'bottom')
        pp = plot(contrasts,pcrfim,'k--');
        uistack(pp,'bottom')
        for c = 1:length(contrasts)
            lsc = line([contrasts(c) contrasts(c)],[pcrfcm(c)-pcrfcs(c) pcrfcm(c)+pcrfcs(c)],'Color','k');
            uistack(lsc,'bottom')
            lsc = line([contrasts(c) contrasts(c)],[pcrfim(c)-pcrfis(c) pcrfim(c)+pcrfis(c)],'Color','k');
            uistack(lsc,'bottom')
        end
        xlim([-1.1 1.1])
%         ylim([bigMin3-.1*bigMin3 bigMax3+.1*bigMax3])
        set(gca,'tickdir','out')
        xlabel('Contrast')
        title('Cell CRF, between blocks')

        subplot(nw,subw,(ww-1)*subw+4)
        hold on;
        pp = plot(contrasts,pcrfcm-pcrfim,'k');
        uistack(pp,'bottom')
        l0 = line([-1.1 1.1],[0 0],'LineStyle','--','Color',[.5 .5 .5]);
        uistack(l0,'bottom') 
        xlim([-1.1 1.1])
%         ylim([bigMin4 bigMax4])
        set(gca,'tickdir','out')
        xlabel('Contrast')
        title('Cell CRF, difference')
        
        subplot(nw,subw,(ww-1)*subw+5)
        scatter(diffFits(whichNeurons{ww},1),diffFits(whichNeurons{ww},2),30,'MarkerFaceAlpha',.2,'MarkerEdgeColor','none','MarkerFaceColor','k','Marker','o');
        l1 = line([-1 1],[0 0],'LineStyle','--','Color',[.5 .5 .5]);
        uistack(l1,'bottom')
        l2 = line([0 0],[-1 1],'LineStyle','--','Color',[.5 .5 .5]);
        uistack(l2,'bottom')
        xlabel('slope')
        ylabel('int')
%         xlim([-max(max(abs(diffFits))) max(max(abs(diffFits)))])
%         ylim([-max(max(abs(diffFits))) max(max(abs(diffFits)))])
        title('Fit params')
    end
    
    bigMax1 = max(max(maxy1));
    bigMin1 = min(min(miny1));
    bigMax2 = max(max(maxy2));
    bigMin2 = min(min(miny2));
    bigMax3 = max(max(maxy3));
    bigMin3 = min(min(miny3));
    bigMax4 = max(max(maxy4));
    bigMin4 = -max(max(maxy4));
        
    for ww = 1:nw
        subplot(nw,subw,(ww-1)*subw+1)
        xlim([-.5 1.5])
        ylim([bigMin1-.1*bigMin1 bigMax1+.1*bigMax1])
        
        subplot(nw,subw,(ww-1)*subw+2)
        xlim([-1.1 1.1])
        ylim([bigMin2-.1*bigMin2 bigMax2+.1*bigMax2])
        
        subplot(nw,subw,(ww-1)*subw+3)
        xlim([-1.1 1.1])
        ylim([bigMin3-.1*bigMin3 bigMax3+.1*bigMax3])
        
        subplot(nw,subw,(ww-1)*subw+4)
        xlim([-1.1 1.1])
        ylim([bigMin4 bigMax4])
        
        subplot(nw,subw,(ww-1)*subw+5)
        xlim([-max(max(abs(diffFits))) max(max(abs(diffFits)))])
        ylim([-max(max(abs(diffFits))) max(max(abs(diffFits)))])
    end
        
%     printfig(gcf,char(strcat(expRef,{' bCRFs - cell quintiles'})))
    
    
    %%
    printfig(gcf,char(strcat(expRef,{' bCRFs test stats (4 epochs)'})))
    close all
    clearvars -except mouseList expList hemList m
    
end
    