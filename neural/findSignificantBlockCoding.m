subDim = ceil(sqrt(length(mouseList)));
figure;

sigrank = cell(1,length(mouseList));
pValue = cell(1,length(mouseList));
    
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
    trueBlock = expInfo.block.events.highRewardSideValues(1:nt);
    contrasts = getUniqueContrasts(expInfo);
    [~, lbTrials] = selectCondition(expInfo, contrasts, behavioralData, ...
        initTrialConditions('highRewardSide','left','repeatType','all'));
    [~, rbTrials] = selectCondition(expInfo, contrasts, behavioralData, ...
        initTrialConditions('highRewardSide','right','repeatType','all'));
    if hemisphere > 0
        bC_trials = lbTrials;
        bI_trials = rbTrials;
    else
        bC_trials = rbTrials;
        bI_trials = lbTrials;
    end
    
    np = 1000;
    nf = 1000;
    pseudoSessions{m} = nan(np,nt);
    fakeBlocks{m} = nan(np,nt);
    blockStart = 'fixed';
    
    % generate fake blocks
    for f = 1:nf
        if strcmp(blockStart,'fixed')
            firstSide = expInfo.block.paramsValues(1).firstHighSide;
        elseif strcmp(blockStart,'rand')
            firstSide = randsample([-1, 1],1,true);
        end
        b=nan(1,nt);
        switches = cumsum(125+randi(100,1,20));
        for s = 1:length(switches)
            if s == 1
                b(1:switches(s)-1) = firstSide;
            elseif mod(s,2) == 1
                b(switches(s-1):switches(s)-1) = firstSide;
            elseif mod(s,2) == 0
                b(switches(s-1):switches(s)-1) = -firstSide;
            end
        end
        fakeBlocks{m}(f,:) = b(1:nt);
        if hemisphere > 0
            bC_trials_fake{f} = find(fakeBlocks{m}(f,:) < 0);
            bI_trials_fake{f} = find(fakeBlocks{m}(f,:) > 0);
        else
            bC_trials_fake{f} = find(fakeBlocks{m}(f,:) > 0);
            bI_trials_fake{f} = find(fakeBlocks{m}(f,:) < 0);
        end
    end
    
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
                b(1:switches(s)-1) = firstSide;
            elseif mod(s,2) == 1
                b(switches(s-1):switches(s)-1) = firstSide;
            elseif mod(s,2) == 0
                b(switches(s-1):switches(s)-1) = -firstSide;
            end
        end
        pseudoSessions{m}(p,:) = b(1:nt);
        if hemisphere > 0
            bC_trials_pseudo{p} = find(pseudoSessions{m}(p,:) < 0);
            bI_trials_pseudo{p} = find(pseudoSessions{m}(p,:) > 0);
        else
            bC_trials_pseudo{p} = find(pseudoSessions{m}(p,:) > 0);
            bI_trials_pseudo{p} = find(pseudoSessions{m}(p,:) < 0);
        end
    end
    
%     figure;
%     set(gcf,'position',[740 228 560 1030]);
%     subplot(51,1,1)
%     plot(trueBlock,'k','LineWidth',2)
%     axis off
%     set(gca, 'XTickLabels', {})
%     set(gca, 'YTickLabels', {})
%     title('True vs. "fake true" vs. pseudo blocks')
%     
%     for p = 1:25
%         subplot(51,1,p+1)
%         plot(fakeBlocks{m}(p,:),'Color',[1 0 .7])
%         axis off
%         set(gca, 'XTickLabels', {})
%         set(gca, 'YTickLabels', {})
%     end
%     
%     for p = 26:50
%         subplot(51,1,p+1)
%         plot(pseudoSessions{m}(p,:),'Color',[0 .7 1])
%         axis off
%         set(gca, 'XTickLabels', {})
%         set(gca, 'YTickLabels', {})
%     end
    
    %% compute block statistic for each cell
    
    % 1. using the real block
    
    bC_resps = nanmean(whichResps(bC_trials,:),1)';
    bI_resps = nanmean(whichResps(bI_trials,:),1)';
    blockResp = bC_resps - bI_resps;
    
    % 2. for each pseudo session
    bC_resps_pseudo = nan(nc,np);
    bI_resps_pseudo = nan(nc,np);
    blockResp_pseudo = nan(nc,np);
    
    for p = 1:np
        bC_resps_pseudo(:,p) = nanmean(whichResps(bC_trials_pseudo{p},:),1)';
        bI_resps_pseudo(:,p) = nanmean(whichResps(bI_trials_pseudo{p},:),1)';
        blockResp_pseudo(:,p) = bC_resps_pseudo(:,p) - bI_resps_pseudo(:,p);
    end
    
    %3. for each fake block
    bC_resps_fake = nan(nc,nf);
    bI_resps_fake = nan(nc,nf);
    blockResp_fake = nan(nc,nf);
    
    for f = 1:nf
        bC_resps_fake(:,f) = nanmean(whichResps(bC_trials_fake{f},:),1)';
        bI_resps_fake(:,f) = nanmean(whichResps(bI_trials_fake{f},:),1)';
        blockResp_fake(:,f) = bC_resps_fake(:,f) - bI_resps_fake(:,f);
    end
    
%     sigIdx_true = (bC_resps - bI_resps)./(bC_resps + bI_resps);
%     sigAmp_true = bC_resps + bI_resps;
%     propContraCells(m) = sum(sigIdx_true>0)/numel(sigIdx_true);
% 
%     sigIdx_pseudo = (bC_resps_pseudo - bI_resps_pseudo)./(bC_resps_pseudo + bI_resps_pseudo);
%     sigAmp_pseudo = (bC_resps_pseudo + bI_resps_pseudo);
%     propContraCells_pseudo(m,:) = sum(sigIdx_pseudo>0,1)/size(sigIdx_pseudo,1);


%      for c = 1:size(whichResps,2)
%         r = corrcoef(whichResps(:,c),blockStructure');
%         blockCorr(c) = r(2,1);
%      end
    
    %% 
    
    blockSig.real = nan(1,nc);
    sigrank{m}.real = nan(1,nc);
    pValue{m} = nan(1,nc);
    
    for c = 1:nc
        UB = prctile(blockResp_pseudo(c,:),97.5);
        LB = prctile(blockResp_pseudo(c,:),2.5);
        if blockResp(c) < LB ||  blockResp(c) > UB
            blockSig.real(c) = 1;
        else 
            blockSig.real(c) = 0;
        end
        sigrank{m}.real(c) = sum(blockResp_pseudo(c,:) < blockResp(c))/np;
        for f = 1:nf
            if blockResp_fake(c,f) < LB ||  blockResp_fake(c,f) > UB
                blockSig.fake(f,c) = 1;
            else 
                blockSig.fake(f,c) = 0;
            end
            sigrank{m}.fake(f,c) = sum(blockResp_pseudo(c,:) < blockResp_fake(c,f))/np;
        end
    end

    propBS.real{m} = sum(blockSig.real)/nc;
    propBS.fake{m} = sum(blockSig.fake,2)./nc;
    
    propSigRank(m) = sum(propBS.fake{m} < propBS.real{m})/nf;
    
%     for p = 1:1000
%         for c = 1:size(whichResps,2)
%             r = corrcoef(whichResps(:,c),pseudoSessions{m}(p,:));
%             blockCorr_pseudo(c,p) = r(2,1);
%         end
%     end
%     
%     blockSig_corr = zeros(1,length(blockResp));
%     for c = 1:length(blockResp)
%         UB_corr = prctile(blockCorr_pseudo(c,:),97.5);
%         LB_corr = prctile(blockCorr_pseudo(c,:),2.5);
%         if blockCorr(c) < LB_corr ||  blockCorr(c) > UB_corr
%             blockSig_corr(c) = 1;
%         end
%         [~, idx_corr] = min(abs(sort(blockCorr_pseudo(c,:)) - blockCorr(c)));
%         sigrank_corr{m}(c) = idx_corr/1000;
%     end
%     propBS_corr(m) = sum(blockSig_corr)/length(blockCorr);



%%

% exampleSigCells = randsample(1:nc,20);
% exampleFakes = randsample(1:nf,20);
% 
% figure;
% set(gcf,'position',[20 10 560 1200])
% for p = 1:20
%     subplot(10,2,p)
%     hold on;
% set(gca, 'YTickLabels', {})
% h = histogram(blockResp_pseudo(exampleSigCells(p),:),linspace(-.1,.1,101));
% % w = histogram(blockResp_fake(exampleSigCells(p),:),linspace(-.1,.1,101));
% if blockSig.real(exampleSigCells(p))
%     line([blockResp(exampleSigCells(p)) blockResp(exampleSigCells(p))],[0 180],'LineWidth',2,'LineStyle','-','Color','r');
% else
%     line([blockResp(exampleSigCells(p)) blockResp(exampleSigCells(p))],[0 180],'LineWidth',2,'LineStyle','-','Color','k');
% end
% box off
% set(gca,'tickdir','out')
% xlim([-.06 .06]);
% end
% 
% figure;
% set(gcf,'position',[20 10 560 1200])
% 
% for p = 1:length(exampleFakes)
%     if p ~= 1
%         subplot(10,2,p)
%         plotHist = sigrank{m}.fake(exampleFakes(p-1),:);
%         plotProp = propBS.fake{m}(exampleFakes(p-1),:);
% %         set(gca, 'XTickLabels', {})
%     else
%         subplot(10,2,1)
%         plotHist = sigrank{m}.real;
%         plotProp = propBS.real{m};
% %         xlabel('Percentile')
%         ylabel('No. neurons')
%         title('Real block')
%     end
%     set(gca, 'YTickLabels', {})
% 
%     hold on;
%     h = histogram(plotHist*100,linspace(0,100,41),'FaceColor',[.5 .5 .5]);
%     hold on
%     maxy = max(h.Values);
%     LB = fill([0; 0; 2.5; 2.5],[0; h.Values(1); h.Values(1); 0],'r','LineStyle','none','FaceAlpha',.5);
%     UB = fill([97.5; 97.5; 100; 100],[0; h.Values(end); h.Values(end); 0],'r','LineStyle','none','FaceAlpha',.5);
% %         text(2,maxy*1.1,strcat(num2str(round(plotProp*100)),'% p < 0.05'),'Color',[.5 0 0])
%         text(2,maxy*1.1,strcat(num2str(round(plotProp*100)),'%'),'Color',[.5 0 0])
% 
%     ylim([0 maxy*1.2])
%     box off
%     set(gca,'tickdir','out')
% end
% 
% figure;
% h = histogram(sigrank{m}.fake*100,linspace(0,100,41),'FaceColor',[.5 .5 .5],'normalization','probability');
% hold on;
% maxy = max(h.Values);
% LB = fill([0; 0; 2.5; 2.5],[0; h.Values(1); h.Values(1); 0],'r','LineStyle','none','FaceAlpha',.5);
% UB = fill([97.5; 97.5; 100; 100],[0; h.Values(end); h.Values(end); 0],'r','LineStyle','none','FaceAlpha',.5);
% plotProp = (h.Values(1)+h.Values(end))/sum(h.Values);
% text(2,maxy*1,strcat(num2str(round(plotProp*100)),'%  p < 0.05'),'Color',[.5 0 0])
% box off
% set(gca,'tickdir','out')
% xlabel('Percentile')
% ylabel('Prop. neurons')
% title('Significant neurons over all fake blocks')
% 
% figure;
% h = histogram(propBS.fake{1}*100,41,'FaceColor',[1 0 1]);
% hold on;
% maxy = max(h.Values);
% line([propBS.real{m}*100 propBS.real{m}*100],[0 maxy],'Color',[0 0 0],'LineWidth',2)
% box off
% set(gca,'tickdir','out')
% xlabel('Percent sig. neurons')
% ylabel('No. sessions')
% title('Real block vs fake blocks')

    %% plot 'block value' of every cell for real vs pseudo sessions
%     set(gcf,'position',[32 80 2560 1556]);
% 
%     mouseName = char(mouseList{m});
%     expDate = char(expList{m}{1});
%     expNum = expList{m}{2};
%     hemisphere = hemList(m);
%     [expRef, ~] = data.constructExpRef(mouseName,expDate,expNum);
%     subplot(subDim,subDim,m)
%     %     histogram(pValue{m},linspace(0,1,21))
%     h = histogram((sigrank{m}.real)*100,linspace(0,100,41),'FaceColor',[.5 .5 .5]);
%     hold on
%     maxy = max(h.Values);
%     LB = fill([0; 0; 2.5; 2.5],[0; h.Values(1); h.Values(1); 0],'r','LineStyle','none','FaceAlpha',.5);
%     UB = fill([97.5; 97.5; 100; 100],[0; h.Values(end); h.Values(end); 0],'r','LineStyle','none','FaceAlpha',.5);
% 
%     ylim([0 maxy*1.2])
%     box off
%     set(gca,'tickdir','out')
%     if m ~= 1
%     set(gca, 'XTickLabels', {})
%     text(2,maxy*1.1,strcat(num2str(round(propBS.real{m}*100)),'%'),'Color',[.5 0 0])
%     else
%     xlabel('Percentile')
%     ylabel('No. neurons')
%     text(2,maxy*1.1,strcat(num2str(round(propBS.real{m}*100)),'% p < 0.05'),'Color',[.5 0 0])
%     end
%     title(strcat(expRef),'Interpreter','none')
    
%      set(gcf,'position',[32 80 2560 1556]);
% 
%     mouseName = char(mouseList{m});
%     expDate = char(expList{m}{1});
%     expNum = expList{m}{2};
%     hemisphere = hemList(m);
%     [expRef, ~] = data.constructExpRef(mouseName,expDate,expNum);
%     subplot(subDim,subDim,m)
%     UB = prctile(propContraCells_pseudo(m,:),97.5);
%     LB = prctile(propContraCells_pseudo(m,:),2.5);
%     %     histogram(pValue{m},linspace(0,1,21))
%     h = histogram(propContraCells_pseudo(m,:)*100,linspace(10,90,41),'FaceColor',[.5 .5 .5]);
%     hold on
%     maxy = max(h.Values);
%     if propContraCells(m) < LB || propContraCells(m) > UB
%         line([propContraCells(m)*100 propContraCells(m)*100],[0 maxy],'LineStyle','-','Color',[1 0 0],'LineWidth',2);
%     else
%         line([propContraCells(m)*100 propContraCells(m)*100],[0 maxy],'LineStyle','-','Color',[0 0 0],'LineWidth',2);
%     end
%     ylim([0 maxy*1.2])
%     box off
%     set(gca,'tickdir','out')
%     if m ~= 1
%     set(gca, 'XTickLabels', {})
%     else
%     xlabel('Contra block–preferring cells (%)')
%     ylabel('No. neurons')
%     end
%     title(strcat(expRef),'Interpreter','none')
    
%     figure;
    set(gcf,'position',[32 80 2560 1556]);
%     for m = 1:length(mouseList)
%     mouseName = char(mouseList{m});
%     expDate = char(expList{m}{1});
%     expNum = expList{m}{2};
%     hemisphere = hemList(m);
%     [expRef, ~] = data.constructExpRef(mouseName,expDate,expNum);
    subplot(subDim,subDim,m)
    h = histogram(propBS.fake{m}*100,41,'FaceColor',[1 0 1]);
    hold on;
    maxy = max(h.Values);
    line([propBS.real{m}*100 propBS.real{m}*100],[0 maxy],'Color',[0 0 0],'LineWidth',2)
    if propSigRank(m)*100 >= 95
        text(propBS.real{m}*100,maxy*1.05,'*','HorizontalAlignment','center')
    elseif propSigRank(m)*100 >= 99
        text(propBS.real{m}*100,maxy*1.05,'**','HorizontalAlignment','center')
    elseif propSigRank(m)*100 >= 99.9
        text(propBS.real{m}*100,maxy*1.05,'***','HorizontalAlignment','center')
    else
        text(propBS.real{m}*100,maxy*1.05,'ns','HorizontalAlignment','center')
    end
    box off
    set(gca,'tickdir','out')
    if m == 1
        xlabel('Percent sig. neurons')
        ylabel('No. sessions')
    else
        set(gca, 'XTickLabels', {})
    end
%     legend(strcat({'significant = '},num2str(round(propBS(m)*100)),'%'),'Location','ne');
%     legend boxoff
    title(strcat(expRef),'Interpreter','none')
%     end
    
    
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

% figure;
% subplot(5,10,1)
% line([0 .25],[0 .25],'LineStyle','--','Color',[.5 .5 .5]);
% hold on;
% p1 = scatter(bC_resps(blockSig.real==0,1),bI_resps(blockSig.real==0,1),'.','MarkerEdgeColor',[.5 .5 .5]);
% hold on;
% p2 = scatter(bC_resps(blockSig.real==1,1),bI_resps(blockSig.real==1,1),'r.');
% xlim([0 .25])
% ylim([0 .25])
% set(gca,'tickdir','out')
% legend([p1 p2],'not sig','sig (true)','location','northwest')
% legend boxoff
% axis square
% 
% for f = 2:50
%     subplot(5,10,f)
%     line([0 .25],[0 .25],'LineStyle','--','Color',[.5 .5 .5]);
%     hold on;
%     p1 = scatter(bC_resps_fake(blockSig.fake(f,:)==0,f),bI_resps_fake(blockSig.fake(f,:)==0,f),'.','MarkerEdgeColor',[.5 .5 .5]);
%     hold on;
%     p2 = scatter(bC_resps_fake(blockSig.fake(f,:)==1,f),bI_resps_fake(blockSig.fake(f,:)==1,f),'b.');
%     p3 = scatter(bC_resps_fake(blockSig.fake(f,:)==1 & blockSig.real==1,f),...
%         bI_resps_fake(blockSig.fake(f,:)==1 & blockSig.real==1,f),'.','MarkerEdgeColor',[1 0 1]);
%     xlim([0 .25])
%     ylim([0 .25])
%     set(gca,'tickdir','out')
%     if f == 2
%         legend([p1 p2 p3],'not sig','sig (fake)','sig (true & fake)','location','northwest');
%         legend boxoff
%     end
%     axis square
% end
% 
% sigIdx = (bC_resps - bI_resps)./(bC_resps + bI_resps);
% sigAmp = bC_resps + bI_resps;
% 
% figure;
% subplot(5,10,1)
% hold on;
% p1 = scatter(sigIdx(blockSig.real==0,1),sigAmp(blockSig.real==0,1),'.','MarkerEdgeColor',[.5 .5 .5]);
% hold on;
% p2 = scatter(sigIdx(blockSig.real==1,1),sigAmp(blockSig.real==1,1),'r.');
% 
% set(gca,'tickdir','out')
% legend([p1 p2],'not sig','sig (true)','location','northwest')
% legend boxoff
% axis square
% 
% for f = 2:50
%     sigIdx = (bC_resps_fake(:,f) - bI_resps_fake(:,f))./(bC_resps_fake(:,f) + bI_resps_fake(:,f));
% sigAmp = bC_resps_fake(:,f) + bI_resps_fake(:,f);
% 
%     subplot(5,10,f)
%     hold on;
%     p1 = scatter(sigIdx(blockSig.fake(f,:)==0,1),sigAmp(blockSig.fake(f,:)==0,1),'.','MarkerEdgeColor',[.5 .5 .5]);
%     hold on;
%     p2 = scatter(sigIdx(blockSig.fake(f,:)==1,1),sigAmp(blockSig.fake(f,:)==1,1),'b.');
%     p3 = scatter(sigIdx(blockSig.fake(f,:)==1 & blockSig.real==1,1),...
%         sigAmp(blockSig.fake(f,:)==1 & blockSig.real==1,1),'.','MarkerEdgeColor',[1 0 1]);
% 
%     set(gca,'tickdir','out')
%     if f == 2
%         legend([p1 p2 p3],'not sig','sig (fake)','sig (true & fake)','location','northwest');
%         legend boxoff
%     end
%     axis square
% end
%%
    clearvars -except mouseList expList hemList pValue propBS subDim sigrank pseudoSessions propContraCells propContraCells_pseudo propSigRank
    
end
