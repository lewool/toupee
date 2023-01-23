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

    trialTypes = getTrialTypes(expInfo, behavioralData, 'early');
    for b = 1:2        
        for c = 1:length(contrasts)
            if c < 5
                whichTrials = trialTypes.intVar.all.contrast_direction_block{c,1,b};
                maxVels.impulsive(m,c,b) = -nanmedian(behavioralData.wheelMoves.epochs(5).peakVel(whichTrials));
                RTs.impulsive(m,c,b) = nanmedian(behavioralData.wheelMoves.epochs(5).onsetTimes(whichTrials) - behavioralData.eventTimes(1).daqTime(whichTrials));
            elseif c > 5
                whichTrials = trialTypes.intVar.all.contrast_direction_block{c,2,b};
                maxVels.impulsive(m,c+1,b) = -nanmedian(behavioralData.wheelMoves.epochs(5).peakVel(whichTrials));
                RTs.impulsive(m,c+1,b) = nanmedian(behavioralData.wheelMoves.epochs(5).onsetTimes(whichTrials) - behavioralData.eventTimes(1).daqTime(whichTrials));
            else
                wA = trialTypes.intVar.all.contrast_direction_block{c,1,b};
                wB = trialTypes.intVar.all.contrast_direction_block{c,2,b};
                maxVels.impulsive(m,c,b) = -nanmedian(behavioralData.wheelMoves.epochs(5).peakVel(wA));
                RTs.impulsive(m,c,b) = nanmedian(behavioralData.wheelMoves.epochs(5).onsetTimes(wA) - behavioralData.eventTimes(1).daqTime(wA));
                maxVels.impulsive(m,c+1,b) = -nanmedian(behavioralData.wheelMoves.epochs(5).peakVel(wB));
                RTs.impulsive(m,c+1,b) = nanmedian(behavioralData.wheelMoves.epochs(5).onsetTimes(wB) - behavioralData.eventTimes(1).daqTime(wB));
            end
        end
    end

    trialTypes = getTrialTypes(expInfo, behavioralData, 'late');
    for b = 1:2
        for c = 1:length(contrasts)
            if c < 5
                whichTrials = trialTypes.intVar.all.contrast_direction_block{c,1,b};
                maxVels.patient(m,c,b) = -nanmedian(behavioralData.wheelMoves.epochs(5).peakVel(whichTrials));
                RTs.patient(m,c,b) = nanmedian(behavioralData.wheelMoves.epochs(5).onsetTimes(whichTrials) - behavioralData.eventTimes(2).daqTime(whichTrials));
            elseif c > 5
                whichTrials = trialTypes.intVar.all.contrast_direction_block{c,2,b};
                maxVels.patient(m,c+1,b) = -nanmedian(behavioralData.wheelMoves.epochs(5).peakVel(whichTrials));
                RTs.patient(m,c+1,b) = nanmedian(behavioralData.wheelMoves.epochs(5).onsetTimes(whichTrials) - behavioralData.eventTimes(2).daqTime(whichTrials));
            elseif c == 5
                wA = trialTypes.intVar.all.contrast_direction_block{c,1,b};
                wB = trialTypes.intVar.all.contrast_direction_block{c,2,b};
                maxVels.patient(m,c,b) = -nanmedian(behavioralData.wheelMoves.epochs(5).peakVel(wA));
                RTs.patient(m,c,b) = nanmedian(behavioralData.wheelMoves.epochs(5).onsetTimes(wA) - behavioralData.eventTimes(2).daqTime(wA));
                maxVels.patient(m,c+1,b) = -nanmedian(behavioralData.wheelMoves.epochs(5).peakVel(wB));
                RTs.patient(m,c+1,b) = nanmedian(behavioralData.wheelMoves.epochs(5).onsetTimes(wB) - behavioralData.eventTimes(2).daqTime(wB));
            end
        end
    end
    
end

%%
color1 = [0.1 0.7 0.1]; 
color2 = [1 .6 0];
ml = cellfun(@char,mouseList,'UniformOutput',false);
um = unique(ml);
timing = 1;
for u = 1:length(um)
    
    meanMaxVels(u,:,1,1) = nanmean(maxVels.impulsive(strcmp(ml,um{u}),:,1),1);
    meanMaxVels(u,:,1,2) = nanmean(maxVels.patient(strcmp(ml,um{u}),:,1),1);
    meanMaxVels(u,:,2,1) = nanmean(maxVels.impulsive(strcmp(ml,um{u}),:,2),1);
    meanMaxVels(u,:,2,2) = nanmean(maxVels.patient(strcmp(ml,um{u}),:,2),1);
    meanRTs(u,:,1,1) = nanmean(RTs.impulsive(strcmp(ml,um{u}),:,1),1);
    meanRTs(u,:,1,2) = nanmean(RTs.patient(strcmp(ml,um{u}),:,1),1);
    meanRTs(u,:,2,1) = nanmean(RTs.impulsive(strcmp(ml,um{u}),:,2),1);
    meanRTs(u,:,2,2) = nanmean(RTs.patient(strcmp(ml,um{u}),:,2),1);
end
pc = [-120 -70 -32 -25 -20 20 25 32 70 120];

figure;
subplot(1,4,1)
range = [1:5];
errorbar(pc(range),nanmean(meanMaxVels(:,range,1,1),1),nanstd(meanMaxVels(:,range,1,1),[],1)/sqrt(size(meanMaxVels(:,range,1,1),1)),...
        'capsize',0,'Marker','o','MarkerSize',10,'MarkerFaceColor',color1,'MarkerEdgeColor','w','Color',color1,'LineWidth',1,'LineStyle',':')
    hold on;
errorbar(pc(range),nanmean(meanMaxVels(:,range,2,1),1),nanstd(meanMaxVels(:,range,2,1),[],1)/sqrt(size(meanMaxVels(:,range,2,1),1)),...
        'capsize',0,'Marker','o','MarkerSize',10,'MarkerFaceColor',color2,'MarkerEdgeColor','w','Color',color2,'LineWidth',1,'LineStyle','-')
    
range = [6:10];
errorbar(pc(range),nanmean(meanMaxVels(:,range,1,1),1),nanstd(meanMaxVels(:,range,1,1),[],1)/sqrt(size(meanMaxVels(:,range,1,1),1)),...
        'capsize',0,'Marker','o','MarkerSize',10,'MarkerFaceColor',color1,'MarkerEdgeColor','w','Color',color1,'LineWidth',1,'LineStyle',':')
hold on;
errorbar(pc(range),nanmean(meanMaxVels(:,range,2,1),1),nanstd(meanMaxVels(:,range,2,1),[],1)/sqrt(size(meanMaxVels(:,range,2,1),1)),...
        'capsize',0,'Marker','o','MarkerSize',10,'MarkerFaceColor',color2,'MarkerEdgeColor','w','Color',color2,'LineWidth',1,'LineStyle','-')
xlim([-130 130])
ylim([-120 120])
prettyPlot(gca)
xticks([-120 -70 -20 20 70 120])
xticklabels({'100' '50' '0' '0' '50' '100'})
xlabel('Left contrast           Right contrast')
ylabel('Wheel velocity (mm/s)')
legend('Left','Right','location','nw')
legend('boxoff')

subplot(1,4,2)

range = [1:5];
errorbar(pc(range),nanmean(meanMaxVels(:,range,1,2),1),nanstd(meanMaxVels(:,range,1,2),[],1)/sqrt(size(meanMaxVels(:,range,1,2),1)),...
        'capsize',0,'Marker','o','MarkerSize',10,'MarkerFaceColor',color1,'MarkerEdgeColor','w','Color',color1,'LineWidth',1,'LineStyle',':')
    hold on;
errorbar(pc(range),nanmean(meanMaxVels(:,range,2,2),1),nanstd(meanMaxVels(:,range,2,2),[],1)/sqrt(size(meanMaxVels(:,range,2,2),1)),...
        'capsize',0,'Marker','o','MarkerSize',10,'MarkerFaceColor',color2,'MarkerEdgeColor','w','Color',color2,'LineWidth',1,'LineStyle','-')
    
range = [6:10];
errorbar(pc(range),nanmean(meanMaxVels(:,range,1,2),1),nanstd(meanMaxVels(:,range,1,2),[],1)/sqrt(size(meanMaxVels(:,range,1,2),1)),...
        'capsize',0,'Marker','o','MarkerSize',10,'MarkerFaceColor',color1,'MarkerEdgeColor','w','Color',color1,'LineWidth',1,'LineStyle',':')
hold on;
errorbar(pc(range),nanmean(meanMaxVels(:,range,2,2),1),nanstd(meanMaxVels(:,range,2,2),[],1)/sqrt(size(meanMaxVels(:,range,2,2),1)),...
        'capsize',0,'Marker','o','MarkerSize',10,'MarkerFaceColor',color2,'MarkerEdgeColor','w','Color',color2,'LineWidth',1,'LineStyle','-')
xlim([-130 130])
ylim([-120 120])
prettyPlot(gca)
xticks([-120 -70 -20 20 70 120])
xticklabels({'100' '50' '0' '0' '50' '100'})
xlabel('Left contrast           Right contrast')
ylabel('Wheel velocity (mm/s)')
legend('Left','Right','location','nw')
legend('boxoff')

subplot(1,4,3)
range = [1:5];
errorbar(pc(range),nanmedian(meanRTs(:,range,1,1),1),nanstd(meanRTs(:,range,1,1),[],1)/sqrt(size(meanRTs(:,range,1,1),1)),...
        'capsize',0,'Marker','o','MarkerSize',10,'MarkerFaceColor',color1,'MarkerEdgeColor','w','Color',color1,'LineWidth',1,'LineStyle',':')
hold on;
errorbar(pc(range),nanmedian(meanRTs(:,range,2,1),1),nanstd(meanRTs(:,range,2,1),[],1)/sqrt(size(meanRTs(:,range,2,1),1)),...
        'capsize',0,'Marker','o','MarkerSize',10,'MarkerFaceColor',color2,'MarkerEdgeColor','w','Color',color2,'LineWidth',1,'LineStyle','-')
    
range = [6:10];
errorbar(pc(range),nanmedian(meanRTs(:,range,1,1),1),nanstd(meanRTs(:,range,1,1),[],1)/sqrt(size(meanRTs(:,range,1,1),1)),...
        'capsize',0,'Marker','o','MarkerSize',10,'MarkerFaceColor',color1,'MarkerEdgeColor','w','Color',color1,'LineWidth',1,'LineStyle',':')
hold on;
errorbar(pc(range),nanmedian(meanRTs(:,range,2,1),1),nanstd(meanRTs(:,range,2,1),[],1)/sqrt(size(meanRTs(:,range,2,1),1)),...
        'capsize',0,'Marker','o','MarkerSize',10,'MarkerFaceColor',color2,'MarkerEdgeColor','w','Color',color2,'LineWidth',1,'LineStyle','-')
prettyPlot(gca)
xlim([-130 130])
xticks([-120 -70 -20 20 70 120])
xticklabels({'100' '50' '0' '0' '50' '100'})
xlabel('Left contrast           Right contrast')
ylabel('Response time (s)')

subplot(1,4,4)
range = [1:5];
errorbar(pc(range),nanmedian(meanRTs(:,range,1,2),1),nanstd(meanRTs(:,range,1,2),[],1)/sqrt(size(meanRTs(:,range,1,2),1)),...
        'capsize',0,'Marker','o','MarkerSize',10,'MarkerFaceColor',color1,'MarkerEdgeColor','w','Color',color1,'LineWidth',1,'LineStyle',':')
hold on;
errorbar(pc(range),nanmedian(meanRTs(:,range,2,2),1),nanstd(meanRTs(:,range,2,2),[],1)/sqrt(size(meanRTs(:,range,2,2),1)),...
        'capsize',0,'Marker','o','MarkerSize',10,'MarkerFaceColor',color2,'MarkerEdgeColor','w','Color',color2,'LineWidth',1,'LineStyle','-')
    
range = [6:10];
errorbar(pc(range),nanmedian(meanRTs(:,range,1,2),1),nanstd(meanRTs(:,range,1,2),[],1)/sqrt(size(meanRTs(:,range,1,2),1)),...
        'capsize',0,'Marker','o','MarkerSize',10,'MarkerFaceColor',color1,'MarkerEdgeColor','w','Color',color1,'LineWidth',1,'LineStyle',':')
hold on;
errorbar(pc(range),nanmedian(meanRTs(:,range,2,2),1),nanstd(meanRTs(:,range,2,2),[],1)/sqrt(size(meanRTs(:,range,2,2),1)),...
        'capsize',0,'Marker','o','MarkerSize',10,'MarkerFaceColor',color2,'MarkerEdgeColor','w','Color',color2,'LineWidth',1,'LineStyle','-')
prettyPlot(gca)
xlim([-130 130])
xticks([-120 -70 -20 20 70 120])
xticklabels({'100' '50' '0' '0' '50' '100'})
xlabel('Left contrast           Right contrast')
ylabel('Response time (s)')

%%
figure
for s = 1:2
    for m = 1:length(mouseList)
        subplot(1,2,s)
        range = [1:5];
        ll = plot(pc(range),maxVels_imp(m,range,s),'LineStyle','-','Color',[.6 .6 .6]);
        ll.Color(4) = .3;
        hold on
        range = [6:10];
        ll = plot(pc(range),maxVels_imp(m,range,s),'LineStyle','-','Color',[.6 .6 .6]);
        ll.Color(4) = .3;
        hold on    

    end

xlim([-130 130])
ylim([-220 220])
prettyPlot(gca)
xticks([-120 -70 -20 20 70 120])
xticklabels({'100' '50' '0' '0' '50' '100'})
xlabel('Left contrast           Right contrast')
ylabel('Wheel velocity (mm/s)')
end