for m = 1:length(mouseList)
    
    mouseName = char(mouseList{m});
    expDate = char(expList{m}{1});
    expNum = expList{m}{2};
    hemisphere = hemList(m);
    [expRef, ~] = data.constructExpRef(mouseName,expDate,expNum);
    fprintf('Loading %s...\n',expRef)
    [behavioralData, expInfo] = data.loadBehavioralDataset(mouseName, expDate, expNum);
    

%%

    et = behavioralData;
    contrasts = getUniqueContrasts(expInfo);
    nt = length(et.eventTimes(1).daqTime);

    trialCorrectChoice = expInfo.block.events.correctResponseValues(1:nt);
    trueStimuli = expInfo.block.events.contrastValues(1:nt);
    sidedStimuli = expInfo.block.events.contrastValues(1:nt);
    sidedStimuli(sidedStimuli == 0) = eps;
    sidedStimuli(abs(sidedStimuli) < .05) = ...
        sidedStimuli(abs(sidedStimuli) < .05).* ...
        trialCorrectChoice(abs(sidedStimuli) < .05);
    trueChoices = et.wheelMoves.epochs(5).moveDir;
    trueBlock = expInfo.block.events.highRewardSideValues(1:nt);
    trueFeedback = zeros(1,nt);
    trueFeedback(trueChoices .* trialCorrectChoice > 0) = 1;
    trueFeedback(trueChoices .* trialCorrectChoice < 0) = 0;
    trueValue(sidedStimuli .* trueBlock > 0) = 1;
    trueValue(sidedStimuli .* trueBlock < 0) = 0;
    
    [impTrials, ~] = selectCondition(expInfo, contrasts, et, ...
        initTrialConditions('movementTime','early'));

    maxVels = behavioralData.wheelMoves.epochs(5).peakVel;
    RTs = behavioralData.wheelMoves.epochs(5).onsetTimes - behavioralData.eventTimes(1).daqTime;
    trials = intersect(...
        find((RTs < prctile(RTs,97.5)) & (RTs > prctile(RTs,2.5))),...
        find(~isnan(trueChoices)));
    
    fitData{m,1} = impTrials(trials)';
    fitData{m,2} = maxVels(trials)';
    fitData{m,3} = RTs(trials)';
    fitData{m,4} = sidedStimuli(trials)';
    fitData{m,5} = trueChoices(trials)';
    fitData{m,6} = trueBlock(trials)';
    fitData{m,7} = trueFeedback(trials)';
    fitData{m,8} = trueValue(trials)';

    clearvars -except mouseList expList hemList fitData
    
end

%% compute wheel/RT/imp indices for each mouse, per contrast/block

ml = cellfun(@char,mouseList,'UniformOutput',false);

for u = 1:length(ml)
    fitArray = cell2mat(fitData(u,:));
    
    fitArray(:,2) = abs(fitArray(:,2));
    fitArray(:,4) = abs(fitArray(:,4));

    contrasts = unique(fitArray(:,4));
    
    mmRT_imp = median(fitArray(fitArray(:,1) == 1,3));
    mmRT_pat = median(fitArray(fitArray(:,1) == 0,3));
    
    for c = 1:length(contrasts)
        hiTrials = ...
            (fitArray(:,4) == contrasts(c)) .* ...
            (fitArray(:,8) == 1);
        loTrials = ...
            (fitArray(:,4) == contrasts(c)) .* ...
            (fitArray(:,8) == 0);
        
        impTrials = ...
            (fitArray(:,4) == contrasts(c)) .* ...
            (fitArray(:,1) == 1);
        patTrials = ...
            (fitArray(:,4) == contrasts(c)) .* ...
            (fitArray(:,1) == 0);
        
        hiImpTrials = hiTrials .* impTrials;
        hiPatTrials = hiTrials .* patTrials;
        loImpTrials = loTrials .* impTrials;
        loPatTrials = loTrials .* patTrials;
        
        meanVel_impVal(u,c,1) = mean(fitArray(logical(hiImpTrials),2));
        meanVel_impVal(u,c,2) = mean(fitArray(logical(loImpTrials),2));
        
        meanVel_patVal(u,c,1) = mean(fitArray(logical(hiPatTrials),2));
        meanVel_patVal(u,c,2) = mean(fitArray(logical(loPatTrials),2));
        
        meanRT_impVal(u,c,1) = median(fitArray(logical(hiImpTrials),3));
        meanRT_impVal(u,c,2) = median(fitArray(logical(loImpTrials),3));
        
        meanRT_patVal(u,c,1) = median(fitArray(logical(hiPatTrials),3));
        meanRT_patVal(u,c,2) = median(fitArray(logical(loPatTrials),3));
            
        meanVel_val(u,c,1) = median(fitArray(logical(hiTrials),2));
        meanVel_val(u,c,2) = median(fitArray(logical(loTrials),2));
        
        propImp_val(u,c,1) = mean(fitArray(logical(hiTrials),1));
        propImp_val(u,c,2) = mean(fitArray(logical(loTrials),1));
        
        meanRT_val(u,c,1) = median(fitArray(logical(hiTrials),3));
        meanRT_val(u,c,2) = median(fitArray(logical(loTrials),3));
        
        meanVel_imp(u,c,1) = median(fitArray(logical(impTrials),2));
        meanVel_imp(u,c,2) = median(fitArray(logical(patTrials),2));
                
        meanRT_imp(u,c,1) = median(fitArray(logical(impTrials),3)) - mmRT_imp;
        meanRT_imp(u,c,2) = median(fitArray(logical(patTrials),3)) -mmRT_pat;
              
    end

end

%% average across mice

um = unique(ml);

for u = 1:length(um)
    whichSessions = strcmp(ml,um{u});
    gmVel_val(u,:,:) = mean(meanVel_val(whichSessions,:,:),1);
    gmImp_val(u,:,:) = mean(propImp_val(whichSessions,:,:),1);
    gmRT_val(u,:,:) = mean(meanRT_val(whichSessions,:,:),1);
    gmVel_imp(u,:,:) = mean(meanVel_imp(whichSessions,:,:),1);
    gmRT_imp(u,:,:) = mean(meanRT_imp(whichSessions,:,:),1);
end

%% plot - mean across mice

figure;
plotContrasts = contrasts*100;
subplot(1,4,1)
errorbar(plotContrasts,mean(gmVel_val(:,:,1)),std(gmVel_val(:,:,1))/sqrt(length(um)),...
    'capsize',0,'Marker','o','MarkerSize',10,'MarkerFaceColor',[0 0 1],...
    'MarkerEdgeColor','w','Color',[0 0 1],'LineWidth',1)
hold on
errorbar(plotContrasts,mean(gmVel_val(:,:,2)),std(gmVel_val(:,:,2))/sqrt(length(um)),...
    'capsize',0,'Marker','o','MarkerSize',10,'MarkerFaceColor',[ 0 .5 1],...
    'MarkerEdgeColor','w','Color',[ 0 .5 1],'LineWidth',1)
prettyPlot(gca)
xlim([-5 105])

% subplot(1,4,2)
% errorbar(contrasts,mean(gmImp_val(:,:,1)),std(gmImp_val(:,:,1))/sqrt(length(um)),...
%     'capsize',0,'Marker','o','MarkerSize',10,'MarkerFaceColor',[0 0 1],...
%     'MarkerEdgeColor','w','Color',[0 0 1],'LineWidth',1)
% hold on
% errorbar(contrasts,mean(gmImp_val(:,:,2)),std(gmImp_val(:,:,2))/sqrt(length(um)),...
%     'capsize',0,'Marker','o','MarkerSize',10,'MarkerFaceColor',[ 0 .5 1],...
%     'MarkerEdgeColor','w','Color',[ 0 .5 1],'LineWidth',1)
% prettyPlot(gca)
% xlim([-.05 1.05])

subplot(1,4,2)
errorbar(plotContrasts,mean(gmRT_val(:,:,1)),std(gmRT_val(:,:,1))/sqrt(length(um)),...
    'capsize',0,'Marker','o','MarkerSize',10,'MarkerFaceColor',[0 0 1],...
    'MarkerEdgeColor','w','Color',[0 0 1],'LineWidth',1)
hold on
errorbar(plotContrasts,mean(gmRT_val(:,:,2)),std(gmRT_val(:,:,2))/sqrt(length(um)),...
    'capsize',0,'Marker','o','MarkerSize',10,'MarkerFaceColor',[ 0 .5 1],...
    'MarkerEdgeColor','w','Color',[ 0 .5 1],'LineWidth',1)
prettyPlot(gca)
xlim([-5 105])

subplot(1,4,3)
errorbar(plotContrasts,mean(gmVel_imp(:,:,1)),std(gmVel_imp(:,:,1))/sqrt(length(um)),...
    'capsize',0,'Marker','o','MarkerSize',10,'MarkerFaceColor',[.6 .6 .6],...
    'MarkerEdgeColor','w','Color',[.6 .6 .6],'LineWidth',1)
hold on
errorbar(plotContrasts,mean(gmVel_imp(:,:,2)),std(gmVel_imp(:,:,2))/sqrt(length(um)),...
    'capsize',0,'Marker','o','MarkerSize',10,'MarkerFaceColor',[0 0 0],...
    'MarkerEdgeColor','w','Color',[.0 0 0],'LineWidth',1)
prettyPlot(gca)
xlim([-5 105])

subplot(1,4,4)
errorbar(plotContrasts,mean(gmRT_imp(:,:,1)),std(gmRT_imp(:,:,1))/sqrt(length(um)),...
    'capsize',0,'Marker','o','MarkerSize',10,'MarkerFaceColor',[.6 .6 .6],...
    'MarkerEdgeColor','w','Color',[.6 .6 .6],'LineWidth',1)
hold on
errorbar(plotContrasts,mean(gmRT_imp(:,:,2)),std(gmRT_imp(:,:,2))/sqrt(length(um)),...
    'capsize',0,'Marker','o','MarkerSize',10,'MarkerFaceColor',[0 0 0],...
    'MarkerEdgeColor','w','Color',[.0 0 0],'LineWidth',1)
prettyPlot(gca)
xlim([-5 105])



%% plot - mean across sessions

figure;
plotContrasts = contrasts * 100;
marker = 'none';

subplot(1,5,1)
errorbar(contrasts,mean(propImp_val(:,:,1)),std(propImp_val(:,:,1))/sqrt(length(ml)),...
    'capsize',0,'Marker',marker,'MarkerSize',10,'MarkerFaceColor',[0 0 1],...
    'MarkerEdgeColor','w','Color',[0 0 1],'LineWidth',1)
hold on
errorbar(contrasts,mean(propImp_val(:,:,2)),std(propImp_val(:,:,2))/sqrt(length(ml)),...
    'capsize',0,'Marker',marker,'MarkerSize',10,'MarkerFaceColor',[ 0 .5 1],...
    'MarkerEdgeColor','w','Color',[ 0 .5 1],'LineWidth',1)
prettyPlot(gca)
xlim([-.05 1.05])
xlabel('Contrast (%)')
ylabel('Prop. impulsive trials')
legend('High value','Low value','location','nw')
legend('boxoff')

subplot(1,5,2)
errorbar(plotContrasts,mean(meanVel_impVal(:,:,1)),std(meanVel_impVal(:,:,1))/sqrt(length(ml)),...
    'capsize',0,'Marker',marker,'MarkerSize',10,'MarkerFaceColor',[0 0 1],...
    'MarkerEdgeColor','w','Color',[0 0 1],'LineWidth',1)
hold on
errorbar(plotContrasts,mean(meanVel_impVal(:,:,2)),std(meanVel_impVal(:,:,2))/sqrt(length(ml)),...
    'capsize',0,'Marker',marker,'MarkerSize',10,'MarkerFaceColor',[ 0 .5 1],...
    'MarkerEdgeColor','w','Color',[ 0 .5 1],'LineWidth',1)
prettyPlot(gca)
xlim([-5 105])
xlabel('Contrast (%)')
ylabel('Response time (s)')
legend('High value','Low value','location','nw')
legend('boxoff')

subplot(1,5,3)
errorbar(plotContrasts,mean(meanVel_patVal(:,:,1)),std(meanVel_patVal(:,:,1))/sqrt(length(ml)),...
    'capsize',0,'Marker',marker,'MarkerSize',10,'MarkerFaceColor',[0 0 1],...
    'MarkerEdgeColor','w','Color',[0 0 1],'LineWidth',1)
hold on
errorbar(plotContrasts,mean(meanVel_patVal(:,:,2)),std(meanVel_patVal(:,:,2))/sqrt(length(ml)),...
    'capsize',0,'Marker',marker,'MarkerSize',10,'MarkerFaceColor',[ 0 .5 1],...
    'MarkerEdgeColor','w','Color',[ 0 .5 1],'LineWidth',1)
prettyPlot(gca)
xlim([-5 105])
xlabel('Contrast (%)')
ylabel('Response time (s)')
legend('High value','Low value','location','nw')
legend('boxoff')

subplot(1,5,4)
errorbar(plotContrasts,mean(meanRT_impVal(:,:,1)),std(meanRT_impVal(:,:,1))/sqrt(length(ml)),...
    'capsize',0,'Marker',marker,'MarkerSize',10,'MarkerFaceColor',[0 0 1],...
    'MarkerEdgeColor','w','Color',[0 0 1],'LineWidth',1)
hold on
errorbar(plotContrasts,mean(meanRT_impVal(:,:,2)),std(meanRT_impVal(:,:,2))/sqrt(length(ml)),...
    'capsize',0,'Marker',marker,'MarkerSize',10,'MarkerFaceColor',[ 0 .5 1],...
    'MarkerEdgeColor','w','Color',[ 0 .5 1],'LineWidth',1)
prettyPlot(gca)
xlim([-5 105])
xlabel('Contrast (%)')
ylabel('Response time (s)')
legend('High value','Low value','location','ne')
legend('boxoff')

subplot(1,5,5)
errorbar(plotContrasts,mean(meanRT_patVal(:,:,1)),std(meanRT_patVal(:,:,1))/sqrt(length(ml)),...
    'capsize',0,'Marker',marker,'MarkerSize',10,'MarkerFaceColor',[0 0 1],...
    'MarkerEdgeColor','w','Color',[0 0 1],'LineWidth',1)
hold on
errorbar(plotContrasts,mean(meanRT_patVal(:,:,2)),std(meanRT_patVal(:,:,2))/sqrt(length(ml)),...
    'capsize',0,'Marker',marker,'MarkerSize',10,'MarkerFaceColor',[ 0 .5 1],...
    'MarkerEdgeColor','w','Color',[ 0 .5 1],'LineWidth',1)
prettyPlot(gca)
xlim([-5 105])
xlabel('Contrast (%)')
ylabel('Response time (s)')
legend('High value','Low value','location','ne')
legend('boxoff')


% subplot(1,4,1)
% errorbar(plotContrasts,mean(meanVel_val(:,:,1)),std(meanVel_val(:,:,1))/sqrt(length(ml)),...
%     'capsize',0,'Marker','o','MarkerSize',10,'MarkerFaceColor',[0 0 1],...
%     'MarkerEdgeColor','w','Color',[0 0 1],'LineWidth',1)
% hold on
% errorbar(plotContrasts,mean(meanVel_val(:,:,2)),std(meanVel_val(:,:,2))/sqrt(length(ml)),...
%     'capsize',0,'Marker','o','MarkerSize',10,'MarkerFaceColor',[ 0 .5 1],...
%     'MarkerEdgeColor','w','Color',[ 0 .5 1],'LineWidth',1)
% prettyPlot(gca)
% xlim([-5 105])
% xlabel('Contrast (%)')
% ylabel('Wheel velocity (mm/s)')
% legend('High value','Low value','location','nw')
% legend('boxoff')
% 
% % subplot(1,4,2)
% % errorbar(contrasts,mean(propImp_val(:,:,1)),std(propImp_val(:,:,1))/sqrt(length(ml)),...
% %     'capsize',0,'Marker','o','MarkerSize',10,'MarkerFaceColor',[0 0 1],...
% %     'MarkerEdgeColor','w','Color',[0 0 1],'LineWidth',1)
% % hold on
% % errorbar(contrasts,mean(propImp_val(:,:,2)),std(propImp_val(:,:,2))/sqrt(length(ml)),...
% %     'capsize',0,'Marker','o','MarkerSize',10,'MarkerFaceColor',[ 0 .5 1],...
% %     'MarkerEdgeColor','w','Color',[ 0 .5 1],'LineWidth',1)
% % prettyPlot(gca)
% % xlim([-.05 1.05])
% % xlabel('Contrast (%)')
% % ylabel('Prop. impulsive trials')
% % legend('High value','Low value','location','nw')
% % legend('boxoff')
% 
% subplot(1,4,2)
% errorbar(plotContrasts,mean(meanRT_val(:,:,1)),std(meanRT_val(:,:,1))/sqrt(length(ml)),...
%     'capsize',0,'Marker','o','MarkerSize',10,'MarkerFaceColor',[0 0 1],...
%     'MarkerEdgeColor','w','Color',[0 0 1],'LineWidth',1)
% hold on
% errorbar(plotContrasts,mean(meanRT_val(:,:,2)),std(meanRT_val(:,:,2))/sqrt(length(ml)),...
%     'capsize',0,'Marker','o','MarkerSize',10,'MarkerFaceColor',[ 0 .5 1],...
%     'MarkerEdgeColor','w','Color',[ 0 .5 1],'LineWidth',1)
% prettyPlot(gca)
% xlim([-5 105])
% xlabel('Contrast (%)')
% ylabel('Response time (s)')
% legend('High value','Low value','location','ne')
% legend('boxoff')
% 
% subplot(1,4,3)
% errorbar(plotContrasts,mean(meanVel_imp(:,:,1)),std(meanVel_imp(:,:,1))/sqrt(length(ml)),...
%     'capsize',0,'Marker','o','MarkerSize',10,'MarkerFaceColor',[.6 .6 .6],...
%     'MarkerEdgeColor','w','Color',[.6 .6 .6],'LineWidth',1)
% hold on
% errorbar(plotContrasts,mean(meanVel_imp(:,:,2)),std(meanVel_imp(:,:,2))/sqrt(length(ml)),...
%     'capsize',0,'Marker','o','MarkerSize',10,'MarkerFaceColor',[0 0 0],...
%     'MarkerEdgeColor','w','Color',[.0 0 0],'LineWidth',1)
% prettyPlot(gca)
% xlim([-5 105])
% xlabel('Contrast (%)')
% ylabel('Wheel velocity (mm/s)')
% legend('Impulsive','Patient','location','nw')
% legend('boxoff')
% 
% subplot(1,4,4)
% errorbar(plotContrasts,mean(meanRT_imp(:,:,1)),std(meanRT_imp(:,:,1))/sqrt(length(ml)),...
%     'capsize',0,'Marker','o','MarkerSize',10,'MarkerFaceColor',[.6 .6 .6],...
%     'MarkerEdgeColor','w','Color',[.6 .6 .6],'LineWidth',1)
% hold on
% errorbar(plotContrasts,mean(meanRT_imp(:,:,2)),std(meanRT_imp(:,:,2))/sqrt(length(ml)),...
%     'capsize',0,'Marker','o','MarkerSize',10,'MarkerFaceColor',[0 0 0],...
%     'MarkerEdgeColor','w','Color',[.0 0 0],'LineWidth',1)
% prettyPlot(gca)
% xlim([-5 105])
% xlabel('Contrast (%)')
% ylabel('Mean-subtracted RT (s)')
% legend('Impulsive','Patient','location','ne')
% legend('boxoff')
% 
% 
% 

    
%% scatterplot of wheel vs RT, dots = trials, labeled by contrast or value

fitArray = cell2mat(fitData);
figure;

for c = 1:length(contrasts)   
    trials = abs(fitArray(:,4)) == contrasts(c);
    [RTdens,RTx, ~, ~] = ksdensity(fitArray(logical(trials),3),'bandwidth',.09);
    [wheeldens,wheelx, ~, ~] = ksdensity(abs(fitArray(logical(trials),2)),'bandwidth',3);

    subplot(2,1,1)
    plot(RTx,RTdens,'color',1 - [c c c]/5,'LineWidth',2)
    hold on
    xlim([-0 2])
    prettyPlot(gca)
    
    subplot(2,1,2)
    plot(wheelx,wheeldens,'color',1 - [c c c]/5,'LineWidth',2)
    hold on
    xlim([-0 300])
    prettyPlot(gca)
%     scatter(fitArray(logical(trials),3),fitArray(logical(trials),2),'MarkerEdgeColor','none','MarkerFaceColor',1 - [c c c]/5)
%     hold on
end



% 














