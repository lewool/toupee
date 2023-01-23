subDim = ceil(sqrt(length(mouseList)));
figure;

for m = 1:length(mouseList)
    
    % fetch exp data
    mouseName = char(mouseList{m});
    expDate = char(expList{m}{1});
    expNum = expList{m}{2};
    hemisphere = hemList(m);
    [expRef, ~] = data.constructExpRef(mouseName,expDate,expNum);
    fprintf('fetching %s...\n',expRef)
    [behavioralData, expInfo, neuralData] = data.loadDataset(mouseName, expDate, expNum);
    %% choose which cells to include
    
    [baselineResps, stimResps, pmovResps, movResps, rewResps, preCueResps] = getEpochResps(neuralData.eta);
    whichCells = 'all'; %choose from 'pLabels' array
    if strcmp(whichCells, 'all')
        plotCells = (1:size(neuralData.eta.alignedResps{1},3))';
    elseif strcmp(whichCells,'contraStim')
        if hemisphere > 0
            plotCells = find(neuralData.stats.pValues(:,strcmp(neuralData.stats.labels,'leftStim')) < 0.05);
        else
            plotCells = find(neuralData.stats.pValues(:,strcmp(neuralData.stats.labels,'rightStim')) < 0.05);
        end
    else
        plotCells = find(neuralData.stats.pValues(:,strcmp(neuralData.stats.labels,whichCells)) < 0.05);
    end
    nc = length(plotCells);
    
    whichResps = baselineResps;
    
    %% get blocks (real, fake, pseudo)
    
    % extract block trials
    nt = length(behavioralData.eventTimes(1).daqTime);
%     contrasts = getUniqueContrasts(expInfo);
%     [~, lbTrials] = selectCondition(expInfo, contrasts, behavioralData, ...
%         initTrialConditions('highRewardSide','left','repeatType','all'));
%     [~, rbTrials] = selectCondition(expInfo, contrasts, behavioralData, ...
%         initTrialConditions('highRewardSide','right','repeatType','all'));
    
    rbTrials = [];
    lbTrials = [];
    trimLength = 0;
    newBlocks = [1 1+find((diff(expInfo.block.events.blockSwitchesValues)))];
    for b = 1:length(newBlocks)
        if expInfo.block.events.highRewardSideValues(newBlocks(b)+trimLength) == 1
            try
                rbTrials = [rbTrials (newBlocks(b)+trimLength):newBlocks(b+1)];
            catch
                rbTrials = [rbTrials (newBlocks(b)+trimLength):nt];
            end
        elseif expInfo.block.events.highRewardSideValues(newBlocks(b)+trimLength) == -1
            try
                lbTrials = [lbTrials (newBlocks(b)+trimLength):newBlocks(b+1)];
            catch
                lbTrials = [lbTrials (newBlocks(b)+trimLength):nt];
            end
        end
    end
    
    if hemisphere > 0
        bC_trials = lbTrials;
        bI_trials = rbTrials;
    else
        bC_trials = rbTrials;
        bI_trials = lbTrials;
    end
                    
    np = 1000;
    pseudoSessions{m} = nan(np,nt);
    blockStart = 'fixed';
    
     % generate pseudo sessions
    for p = 1:np
        if strcmp(blockStart,'fixed')
            firstSide = expInfo.block.paramsValues(1).firstHighSide;
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
        
%         b = randsample([-1, 1],nt,true);
        
        pseudoSessions{m}(p,:) = b(1:nt);
        if hemisphere > 0
            bC_trials_pseudo{p} = find(pseudoSessions{m}(p,:) < 0);
            bI_trials_pseudo{p} = find(pseudoSessions{m}(p,:) > 0);
        else
            bC_trials_pseudo{p} = find(pseudoSessions{m}(p,:) > 0);
            bI_trials_pseudo{p} = find(pseudoSessions{m}(p,:) < 0);
        end
    end
    
    %% compute block statistic for each cell
    
    % 1. using the real block
    
    bC_resps = nanmean(whichResps(bC_trials,:),1)';
    bI_resps = nanmean(whichResps(bI_trials,:),1)';
    blockResp = abs(bC_resps - bI_resps);

    % 2. for each pseudo session
    bC_resps_pseudo = nan(nc,np);
    bI_resps_pseudo = nan(nc,np);
    blockResp_pseudo = nan(nc,np);
    
    for p = 1:np
        bC_resps_pseudo(:,p) = nanmean(whichResps(bC_trials_pseudo{p},:),1)';
        bI_resps_pseudo(:,p) = nanmean(whichResps(bI_trials_pseudo{p},:),1)';
        blockResp_pseudo(:,p) = abs(bC_resps_pseudo(:,p) - bI_resps_pseudo(:,p));
    end
    
    %% sum so there is a single pop statistic
    
    trueStatistic = sum(blockResp);
    pseudoStatistics = sum(blockResp_pseudo,1);

    UB = prctile(pseudoStatistics,95);
    if trueStatistic > UB
        popSig = 1;
    else 
        popSig = 0;
    end    
    
    %% plot
    set(gcf,'position',[32 80 2560 1556]);

    mouseName = char(mouseList{m});
    expDate = char(expList{m}{1});
    expNum = expList{m}{2};
    hemisphere = hemList(m);
    [expRef, ~] = data.constructExpRef(mouseName,expDate,expNum);
    subplot(subDim,subDim,m)
    %     histogram(pValue{m},linspace(0,1,21))
    h = histogram(pseudoStatistics,41,'FaceColor',[.5 .5 .5]);
    hold on;
    maxy = max(h.Values);
    if popSig == 1
        line([trueStatistic trueStatistic],[0 maxy],'LineStyle','-','Color','r','LineWidth',2);
    else
        line([trueStatistic trueStatistic],[0 maxy],'LineStyle','-','Color','k','LineWidth',2);
    end
    box off
    
    set(gca,'tickdir','out')
    if m == 1
        xlabel('Session test statistic, \Sigma_c abs(f_{c,contra} – f_{c,ipsi})')
        ylabel('Frequency')
    end
    title(strcat(expRef),'Interpreter','none')
    
    clearvars -except mouseList expList hemList subDim m
    
end