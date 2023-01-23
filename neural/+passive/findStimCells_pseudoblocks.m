if contains(expInfo.block.expDef,'Reward')
    sessionType = 'task';
elseif contains(expInfo.block.expDef,'passive')
    sessionType = 'passive';
else
    error('Unknown session type')
end

alignedResps = neuralData.eta.alignedResps{1};
eventTimes = behavioralData.eventTimes;
eventWindow = neuralData.eta.eventWindow;

Fs = 0.1;
event = 'stimulusOnTimes';
stim_alignedTraces = neuralData.eta.alignedResps{1};
stim_eventWindow = eventWindow;

%designate a baseline window
stim_eventIdx = find(stim_eventWindow == 0);
stim_preTime = [-0.5 0] / Fs;
baselineIdx = stim_eventIdx + stim_preTime(1) : stim_eventIdx - 1;

%compute the mean baseline activity per cell, per trial (trials x neurons)
baselineResps = squeeze(mean(stim_alignedTraces(:,baselineIdx,:),2));

%designate a peristimulus window
stimTime = [0 0.5] / Fs;
stimIdx = stim_eventIdx + stimTime(1) :stim_eventIdx + stimTime(2);

%compute the mean peristimulus activity per cell, per trial (trials x neurons)
stimResps = squeeze(mean(stim_alignedTraces(:,stimIdx,:),2));

whichResps = stimResps - baselineResps;

%% FETCH TRIALS

np = 2000;
nc = size(neuralData.cellResps,2);
nt = numel(expInfo.block.events.endTrialValues);

if strcmp(sessionType,'task')
    contrasts = getUniqueContrasts(expInfo);

    [~, patientTrials] = selectCondition(expInfo, getUniqueContrasts(expInfo), behavioralData,  ...
        initTrialConditions('movementTime','late'));

    trialContrasts = expInfo.block.events.contrastValues;

    [~, leftTrials] = selectCondition(expInfo, contrasts(contrasts<0), behavioralData,  ...
        initTrialConditions('movementTime','late'));
    [~, rightTrials] = selectCondition(expInfo, contrasts(contrasts>0), behavioralData,  ...
        initTrialConditions('movementTime','late'));

    trueResp_left = nanmean(whichResps(leftTrials,:),1);
    trueResp_right = nanmean(whichResps(rightTrials,:),1);

    for p = 1:np
        pseudoContrasts(p,:) = randsample(patientTrials,length(patientTrials));
    end

    for p = 1:np
        leftTrials_pseudo = pseudoContrasts(trialContrasts(pseudoContrasts(p,:)) < 0);
        rightTrials_pseudo = pseudoContrasts(trialContrasts(pseudoContrasts(p,:)) > 0);
        pseudoResps_left(p,:) = nanmean(whichResps(leftTrials_pseudo,:),1);
        pseudoResps_right(p,:) = nanmean(whichResps(rightTrials_pseudo,:),1);
    end
    
elseif strcmp(sessionType,'passive')

    trialContrasts = expInfo.block.events.contrastValues(1:nt);
    trialAzimuths = expInfo.block.events.azimuthValues(1:nt);
    passiveContrasts = sign(trialAzimuths).*trialContrasts;
    passiveContrasts(trialAzimuths == 0) = nan;

    leftTrials = find(passiveContrasts < 0);
    rightTrials = find(passiveContrasts > 0);

    trueResp_left = nanmean(whichResps(leftTrials,:),1);
    trueResp_right = nanmean(whichResps(rightTrials,:),1);
    patientTrials = 1:nt;

    for p = 1:np
        pseudoContrasts(p,:) = randsample(patientTrials,length(patientTrials));
    end

    for p = 1:np
        leftTrials_pseudo = pseudoContrasts(passiveContrasts(pseudoContrasts(p,:)) < 0);
        rightTrials_pseudo = pseudoContrasts(passiveContrasts(pseudoContrasts(p,:)) > 0);
        pseudoResps_left(p,:) = nanmean(whichResps(leftTrials_pseudo,:),1);
        pseudoResps_right(p,:) = nanmean(whichResps(rightTrials_pseudo,:),1);
    end
end


%% COMPUTE RESPONSES

if hemisphere < 0
    trueResp = trueResp_right - trueResp_left;
    pseudoResps = pseudoResps_right - pseudoResps_left;
else
    trueResp = trueResp_left - trueResp_right;
    pseudoResps = pseudoResps_left - pseudoResps_right;
end

%% STATS
stimSig = nan(1,nc);
sigrank = nan(1,nc);

for c = 1:nc
    UB = prctile(pseudoResps(:,c),95);
    if trueResp(c) > UB
        stimSig(c) = 1;
    else 
        stimSig(c) = 0;
    end
    [~, idx] = min(abs(sort(pseudoResps(:,c)) - trueResp(c)));
    sigrank(c) = idx/np;
end

%% PLOT

alpha = 0.01;

figure;
set(gcf,'position',[ 390    1050    1180    290]);
subplot(1,2,1)
nbins = 100/(alpha*100);
h = histogram(sigrank*100,nbins+1,'FaceColor',[.5 .5 .5],'normalization','probability');
hold on;
maxy = max(h.Values);
UB = fill([h.BinEdges(end-1); h.BinEdges(end-1); h.BinEdges(end); h.BinEdges(end)],[0; h.Values(end); h.Values(end); 0],'r','LineStyle','none','FaceAlpha',.5);
plotProp = (h.Values(end))/sum(h.Values);
text(2,maxy*1,strcat(num2str(round(plotProp*100)),{'%  p < '},num2str(alpha)),'Color',[.5 0 0]);
box off
xlim([0 100]);
set(gca,'tickdir','out')
xlabel('Percentile')
ylabel('Prop. neurons')
title('Significant neurons')


colors = {[0 0 .5],[0 0 1],[0 .4 1],[.6 .8 1],[.75 .75 .75],[.8 .7 .7],[.8 .45 .45],[1 0 0],[.5 0 0]};
pval = 1-alpha;
plotCells = find(sigrank > pval);
maxy = max(trueResp(plotCells));
subplot(1,2,2)
hold on;
xlim([-0.5 1.5])

contrasts = [-1 -.5 -.12 -.05 0 .05 .12 .5 1];

set(gca,'tickdir','out')
for c = 1:length(contrasts)
    try
    [~, whichTrials] = selectCondition(expInfo, contrasts(c), behavioralData,  ...
        initTrialConditions('movementTime','late'));
    catch
        whichTrials = find(passiveContrasts == contrasts(c));
    end

    p = plot(...
        eventWindow',...
        nanmean(squeeze(nanmean(neuralData.eta.alignedResps{1}(whichTrials,:,plotCells),1)),2),...
        'LineWidth',2,'Color',colors{c});
    y(c) = max(p.YData);
end

propSig = length(plotCells)/nc;
% text(-.45,max(y)*1.1,strcat(num2str(round(propSig*100)),{'%  p < '},num2str(1-pval)),'Color',[0 0 0]);
ylim([0 max(y)*1.1])
line([0 0],[.0 max(y)*1.1],'LineStyle','--','Color','k');
xlabel('Time from stimulus (s)')
ylabel('Mean activity')
title(strcat({'Stim response of '}, num2str(100*(1-alpha)), 'th percentile'))


%% WORKBENCH

% %% separate screens
% 
% contrastThresh = -10;
% 
% np = 2000;
% nc = size(neuralData.cellResps,2);
% nt = numel(expInfo.block.events.endTrialValues);
% 
% trialContrasts = expInfo.block.events.contrastValues(1:nt);
% trialAzimuths = expInfo.block.events.azimuthValues(1:nt);
% 
% whichTrialsLeft = find(trialAzimuths == -90 & trialContrasts > contrastThresh);
% whichTrialsCenter = find(trialAzimuths == 0 & trialContrasts > contrastThresh);
% whichTrialsRight = find(trialAzimuths == 90 & trialContrasts > contrastThresh);
% 
% leftStimResp = nanmean(whichResps(whichTrialsLeft,:),1);
% centerStimResp = nanmean(whichResps(whichTrialsCenter,:),1);
% rightStimResp = nanmean(whichResps(whichTrialsRight,:),1);
% 
% for p = 1:np
%     pseudoContrasts(p,:) = randsample(trialContrasts,nt);
%     pseudoAzimuths(p,:) = randsample(trialAzimuths,nt);
% end
% 
% for p = 1:np
%     
%     whichTrialsLeft_pseudo = find(pseudoAzimuths(p,:) == -90 & pseudoContrasts(p,:) > contrastThresh);
%     whichTrialsCenter_pseudo = find(pseudoAzimuths(p,:) == 0 & pseudoContrasts(p,:) > contrastThresh);
%     whichTrialsRight_pseudo = find(pseudoAzimuths(p,:) == 90 & pseudoContrasts(p,:) > contrastThresh);
% 
%     leftStimResp_pseudo(p,:) = nanmean(whichResps(whichTrialsLeft_pseudo,:),1);
%     centerStimResp_pseudo(p,:) = nanmean(whichResps(whichTrialsCenter_pseudo,:),1);
%     rightStimResp_pseudo(p,:) = nanmean(whichResps(whichTrialsRight_pseudo,:),1);
% end
% 
% stimSig_left = nan(1,nc);
% sigrank_left = nan(1,nc);
% 
% for c = 1:nc
%     UB = prctile(leftStimResp_pseudo(:,c),97.5);
%     LB = prctile(leftStimResp_pseudo(:,c),2.5);
%     if leftStimResp(c) < LB ||  leftStimResp(c) > UB
%         stimSig_left(c) = 1;
%     else 
%         stimSig_left(c) = 0;
%     end
%     [~, idx] = min(abs(sort(leftStimResp_pseudo(:,c)) - leftStimResp(c)));
%     sigrank_left(c) = idx/np;
% end
% 
% stimSig_center = nan(1,nc);
% sigrank_center = nan(1,nc);
% 
% for c = 1:nc
%     UB = prctile(centerStimResp_pseudo(:,c),97.5);
%     LB = prctile(centerStimResp_pseudo(:,c),2.5);
%     if centerStimResp(c) < LB ||  centerStimResp(c) > UB
%         stimSig_center(c) = 1;
%     else 
%         stimSig_center(c) = 0;
%     end
%     [~, idx] = min(abs(sort(centerStimResp_pseudo(:,c)) - centerStimResp(c)));
%     sigrank_center(c) = idx/np;
% end
% 
% stimSig_right = nan(1,nc);
% sigrank_right = nan(1,nc);
% 
% for c = 1:nc
%     UB = prctile(rightStimResp_pseudo(:,c),97.5);
%     LB = prctile(rightStimResp_pseudo(:,c),2.5);
%     if rightStimResp(c) < LB ||  rightStimResp(c) > UB
%         stimSig_right(c) = 1;
%     else 
%         stimSig_right(c) = 0;
%     end
%     [~, idx] = min(abs(sort(rightStimResp_pseudo(:,c)) - rightStimResp(c)));
%     sigrank_right(c) = idx/np;
% end
% 
% % figure;
% % set(gcf,'position',[981  1062    1566    268])
% % subplot(1,3,1)
% % h = histogram(sigrank_left*100,41,'FaceColor',[.5 .5 .5],'normalization','probability');
% % hold on;
% % maxy = max(h.Values);
% % LB = fill([0; 0; 2.5; 2.5],[0; h.Values(1); h.Values(1); 0],'r','LineStyle','none','FaceAlpha',.5);
% % UB = fill([97.5; 97.5; 100; 100],[0; h.Values(end); h.Values(end); 0],'r','LineStyle','none','FaceAlpha',.5);
% % plotProp = (h.Values(1)+h.Values(end))/sum(h.Values);
% % text(2,maxy*1,strcat(num2str(round(plotProp*100)),'%  p < 0.05'),'Color',[.5 0 0])
% % box off
% % xlim([0 100]);
% % set(gca,'tickdir','out')
% % xlabel('Percentile')
% % ylabel('Prop. neurons')
% % title('Left-stim-sensitive')
% % 
% % subplot(1,3,2)
% % h = histogram(sigrank_center*100,41,'FaceColor',[.5 .5 .5],'normalization','probability');
% % hold on;
% % maxy = max(h.Values);
% % LB = fill([0; 0; 2.5; 2.5],[0; h.Values(1); h.Values(1); 0],'r','LineStyle','none','FaceAlpha',.5);
% % UB = fill([97.5; 97.5; 100; 100],[0; h.Values(end); h.Values(end); 0],'r','LineStyle','none','FaceAlpha',.5);
% % plotProp = (h.Values(1)+h.Values(end))/sum(h.Values);
% % text(2,maxy*1,strcat(num2str(round(plotProp*100)),'%  p < 0.05'),'Color',[.5 0 0])
% % box off
% % xlim([0 100]);
% % set(gca,'tickdir','out')
% % xlabel('Percentile')
% % ylabel('Prop. neurons')
% % title('Center-stim-sensitive')
% % 
% % subplot(1,3,3)
% % h = histogram(sigrank_right*100,41,'FaceColor',[.5 .5 .5],'normalization','probability');
% % hold on;
% % maxy = max(h.Values);
% % LB = fill([0; 0; 2.5; 2.5],[0; h.Values(1); h.Values(1); 0],'r','LineStyle','none','FaceAlpha',.5);
% % UB = fill([97.5; 97.5; 100; 100],[0; h.Values(end); h.Values(end); 0],'r','LineStyle','none','FaceAlpha',.5);
% % plotProp = (h.Values(1)+h.Values(end))/sum(h.Values);
% % text(2,maxy*1,strcat(num2str(round(plotProp*100)),'%  p < 0.05'),'Color',[.5 0 0])
% % box off
% % xlim([0 100]);
% % set(gca,'tickdir','out')
% % xlabel('Percentile')
% % ylabel('Prop. neurons')
% % title('Right-stim-sensitive')
% 
% figure;
% set(gcf,'position',[981  1062    1566    268])
% subplot(1,3,1)
% h = histogram(sigrank_left*100,41,'FaceColor',[.5 .5 .5],'normalization','probability');
% hold on;
% maxy(1) = max(h.Values);
% LB = fill([95; 95; 97.5; 97.5],[0; h.Values(end-1); h.Values(end-1); 0],'r','LineStyle','none','FaceAlpha',.5);
% UB = fill([97.5; 97.5; 100; 100],[0; h.Values(end); h.Values(end); 0],'r','LineStyle','none','FaceAlpha',.5);
% plotProp = (h.Values(end-1)+h.Values(end))/sum(h.Values);
% t(1) = text(2,maxy(1)*1,strcat(num2str(round(plotProp*100)),'%  p < 0.05'),'Color',[.5 0 0]);
% box off
% xlim([0 100]);
% set(gca,'tickdir','out')
% xlabel('Percentile')
% ylabel('Prop. neurons')
% title('Left screen')
% 
% subplot(1,3,2)
% h = histogram(sigrank_center*100,41,'FaceColor',[.5 .5 .5],'normalization','probability');
% hold on;
% maxy(2) = max(h.Values);
% LB = fill([95; 95; 97.5; 97.5],[0; h.Values(end-1); h.Values(end-1); 0],'r','LineStyle','none','FaceAlpha',.5);
% UB = fill([97.5; 97.5; 100; 100],[0; h.Values(end); h.Values(end); 0],'r','LineStyle','none','FaceAlpha',.5);
% plotProp = (h.Values(end-1)+h.Values(end))/sum(h.Values);
% t(2) = text(2,maxy(2)*1,strcat(num2str(round(plotProp*100)),'%  p < 0.05'),'Color',[.5 0 0]);
% box off
% xlim([0 100]);
% set(gca,'tickdir','out')
% xlabel('Percentile')
% ylabel('Prop. neurons')
% title('Center screen')
% 
% subplot(1,3,3)
% h = histogram(sigrank_right*100,41,'FaceColor',[.5 .5 .5],'normalization','probability');
% hold on;
% maxy(3) = max(h.Values);
% LB = fill([95; 95; 97.5; 97.5],[0; h.Values(end-1); h.Values(end-1); 0],'r','LineStyle','none','FaceAlpha',.5);
% UB = fill([97.5; 97.5; 100; 100],[0; h.Values(end); h.Values(end); 0],'r','LineStyle','none','FaceAlpha',.5);
% plotProp = (h.Values(end-1)+h.Values(end))/sum(h.Values);
% t(3) = text(2,maxy(3)*1,strcat(num2str(round(plotProp*100)),'%  p < 0.05'),'Color',[.5 0 0]);
% box off
% xlim([0 100]);
% set(gca,'tickdir','out')
% xlabel('Percentile')
% ylabel('Prop. neurons')
% title('Right screen')
% 
% for s = 1:3
%     subplot(1,3,s)
%     ylim([0 max(maxy)*1.05])
%     t(s).Position = [2 max(maxy)*.95 0];
% end
%     
% 
% 
% 
% %% mean plot
% 
% whichScreen = [-90 0 90];
% colors = cell(1,3);
% colors{1} = {[0 0 .25],[0 0 .5],[0 0 1],[0 .4 1],[.6 .8 1],[.75 .75 .75]};
% colors{2} = {[0 0 0],[.15 .15 .15 ],[.3 .3 .3],[.45 .45 .45],[.6 .6 .6],[.75 .75 .75]};
% colors{3} = {[.25 0 0],[.5 0 0 ],[1 0 0],[.8 .45 .45],[.8 .7 .7],[.75 .75 .75]};
% contrasts = getUniqueContrasts(expInfo);
% 
% figure;
% set(gcf, 'position',[1256 440 828 784]);
% maxy = .30;
% pval = 0.999;
% 
% plotCells = find(sigrank_left > pval);
% for s = 1:length(whichScreen)
%     subplot(3,3,s)
%     line([0 0],[.0 maxy],'LineStyle','--','Color','k');
%     plotCol = fliplr(colors{s});
%     hold on;
%     ylim([0 maxy])
%     xlim([-0.5 1.5])
%     if s == 1
%         ylabel('norm. resp.')
%         text(1,maxy*.95,strcat(sprintf('%0.2f',length(plotCells)/nc*100),'%'))
%     end
%     set(gca,'tickdir','out')
%     for c = 1:length(contrasts)
%         whichTrials = find(...
%             expInfo.block.events.azimuthValues == whichScreen(s) ...
%             & ...
%             expInfo.block.events.contrastValues == contrasts(c));
%         plot(...
%             eventWindow,...
%             nanmean(squeeze(nanmean(neuralData.eta.alignedResps{1}(whichTrials,:,plotCells),1)),2),...
%             'LineWidth',2,...
%             'Color',plotCol{c})
%     end
%     
% end
% 
% plotCells = find(sigrank_center > pval);
% for s = 1:length(whichScreen)
%     subplot(3,3,s+3)
%     line([0 0],[.0 maxy],'LineStyle','--','Color','k');
%     plotCol = fliplr(colors{s});
%     hold on;
%     ylim([0 maxy])
%     xlim([-0.5 1.5])
%     if s == 1
%         ylabel('norm. resp.')
%     end
%     if s == 2
%         text(1,maxy*.95,strcat(sprintf('%0.2f',length(plotCells)/nc*100),'%'))
%     end
%     set(gca,'tickdir','out')
%     for c = 1:length(contrasts)
%         whichTrials = find(...
%             expInfo.block.events.azimuthValues == whichScreen(s) ...
%             & ...
%             expInfo.block.events.contrastValues == contrasts(c));
%         plot(...
%             eventWindow,...
%             nanmean(squeeze(nanmean(neuralData.eta.alignedResps{1}(whichTrials,:,plotCells),1)),2),...
%             'LineWidth',2,...
%             'Color',plotCol{c})
%     end
% end
% 
% plotCells = find(sigrank_right > pval);
% for s = 1:length(whichScreen)
%     subplot(3,3,s+6)
%     line([0 0],[.0 maxy],'LineStyle','--','Color','k');
%     plotCol = fliplr(colors{s});
%     hold on;
%     ylim([0 maxy])
%     xlim([-0.5 1.5])
%     xlabel('time (s)')
%     if s == 1
%         ylabel('norm. resp.')
%     end
%     if s == 3
%         text(1,maxy*.95,strcat(sprintf('%0.2f',length(plotCells)/nc*100),'%'))
%     end
%     set(gca,'tickdir','out')
%     for c = 1:length(contrasts)
%         whichTrials = find(...
%             expInfo.block.events.azimuthValues == whichScreen(s) ...
%             & ...
%             expInfo.block.events.contrastValues == contrasts(c));
%         plot(...
%             eventWindow,...
%             nanmean(squeeze(nanmean(neuralData.eta.alignedResps{1}(whichTrials,:,plotCells),1)),2),...
%             'LineWidth',2,...
%             'Color',plotCol{c})
%     end
% end