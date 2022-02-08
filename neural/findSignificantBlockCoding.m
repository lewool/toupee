subDim = ceil(sqrt(length(mouseList)));
figure;
for m = 1:length(mouseList)
    
    mouseName = char(mouseList{m});
    expDate = char(expList{m}{1});
    expNum = expList{m}{2};
    hemisphere = hemList(m);
    [expRef, ~] = data.constructExpRef(mouseName,expDate,expNum);
    fprintf('fetching %s...\n',expRef)
    [behavioralData, expInfo, neuralData] = data.loadDataset(mouseName, expDate, expNum);
%     [neuralData] = getSignificantActivity(expInfo, behavioralData, neuralData,0);
    contrasts = getUniqueContrasts(expInfo);
%     allTrials = getTrialTypes(expInfo, behavioralData, 'all');
    [~, lbTrials] = selectCondition(expInfo, contrasts, behavioralData, ...
        initTrialConditions('highRewardSide','left','repeatType','random'));
    [~, rbTrials] = selectCondition(expInfo, contrasts, behavioralData, ...
        initTrialConditions('highRewardSide','right','repeatType','random'));
    [baselineResps, stimResps, pmovResps, movResps, rewResps, preCueResps] = getEpochResps(neuralData.eta);
    % [neuralData] = getSignificantActivity(expInfo, behavioralData, neuralData,0);
    
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
 
    whichResps = baselineResps;
    
    if hemisphere > 0

        bC_trials = lbTrials;
        bI_trials = rbTrials;

    else
        bC_trials = rbTrials;
        bI_trials = lbTrials;

    end

    bC_resps = nanmean(whichResps(bC_trials,plotCells),1)';
    bI_resps = nanmean(whichResps(bI_trials,plotCells),1)';
    blockResp = bC_resps - bI_resps;


    %% generate pseudosessions
    nt = length(behavioralData.eventTimes(1).daqTime);
    for p = 1:10000
        b=zeros(1,nt);
        switches = cumsum(125+randi(100,1,10));
        for s = 1:length(switches)
            if s == 1
                b(1:switches(s)-1) = -1;
            elseif mod(s,2) == 1
                b(switches(s-1):switches(s)-1) = -1;
            elseif mod(s,2) == 0
                b(switches(s-1):switches(s)-1) = 1;
            end
        end
        b = b(1:nt);
        flip = randsample([-1, 1],1,true);
%         if b(1) ~= expInfo.block.paramsValues(1).firstHighSide
%             flip = -1;
%         else
%             flip = 1;
%         end
        b = flip*b;
        bC_trials_pseudo{p} = find(b < 0);
        bI_trials_pseudo{p} = find(b > 0);
    end

    for p = 1:10000
    bC_resps_pseudo(:,p) = nanmean(whichResps(bC_trials_pseudo{p},plotCells),1)';
    bI_resps_pseudo(:,p) = nanmean(whichResps(bI_trials_pseudo{p},plotCells),1)';
    blockResp_pseudo(:,p) = bC_resps_pseudo(:,p) - bI_resps_pseudo(:,p);
    end

    blockSig = zeros(1,length(blockResp));
    for c = 1:length(blockResp)
        UB = prctile(blockResp_pseudo(c,:),97.5);
        LB = prctile(blockResp_pseudo(c,:),2.5);
        if blockResp(c) < LB ||  blockResp(c) > UB
            blockSig(c) = 1;
        end
        [~, idx] = min(abs(sort(blockResp_pseudo(c,:)) - blockResp(c)));
        if idx > 500
            pValue{m}(c) = 2*((1000-idx)/1000);
        else
            pValue{m}(c) = 2*idx/1000;
        end
        sigrank{m}(c) = idx/1000;
    end

    propBS(m) = sum(blockSig)/length(blockResp);

    %% plot 'block value' of every cell for real vs pseudo sessions
set(gcf,'position',[32 80 2560 1556]);

mouseName = char(mouseList{m});
expDate = char(expList{m}{1});
expNum = expList{m}{2};
hemisphere = hemList(m);
[expRef, ~] = data.constructExpRef(mouseName,expDate,expNum);
subplot(subDim,subDim,m)
%     histogram(pValue{m},linspace(0,1,21))
h = histogram(sigrank{m}*100,linspace(0,100,41),'FaceColor',[.5 .5 .5]);
hold on
maxy = max(h.Values);
LB = fill([0; 0; 2.5; 2.5],[0; h.Values(1); h.Values(1); 0],'r','LineStyle','none','FaceAlpha',.5);
LB = fill([97.5; 97.5; 100; 100],[0; h.Values(end); h.Values(end); 0],'r','LineStyle','none','FaceAlpha',.5);

ylim([0 maxy*1.2])
box off
set(gca,'tickdir','out')
if m ~= 1
set(gca, 'XTickLabels', {})
text(2,maxy*1.1,strcat(num2str(round(propBS(m)*100)),'%'),'Color',[.5 0 0])
else
xlabel('Percentile')
ylabel('No. neurons')
text(2,maxy*1.1,strcat(num2str(round(propBS(m)*100)),'% p < 0.05'),'Color',[.5 0 0])
end
title(strcat(expRef),'Interpreter','none')
%     legend(strcat({'significant = '},num2str(round(propBS(m)*100)),'%'),'Location','ne');
%     legend boxoff

    
%     subplot(subDim,subDim,m)
%     histogram(blockResp_pseudo,...
%         'normalization','pdf',...
%         'EdgeColor','none',...
%         'FaceAlpha',.5)
%     hold on;
%     histogram(blockResp,...
%         'normalization','pdf',...
%         'EdgeColor','r',...
%         'LineWidth',1.5,...
%         'DisplayStyle','stairs')
%     box off
%     set(gca,'tickdir','out')
%     if m == subDim*(subDim-1)
%         xlabel('Neural "block value"')
%         ylabel('pdf')
%     end
%     legend('pseudo',strcat('true (',num2str(round(propBS(m)*100)),'%)'),'Location','nw');
%     legend boxoff
%     title(expRef,'Interpreter','none')
%%
    clearvars -except mouseList expList hemList pValue propBS subDim sigrank
    
end
