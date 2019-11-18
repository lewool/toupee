%% overview

% this script plots ETAs aligned to any event you want, and compares them
% between two conditions of your choosing (high/low reward, t-1 left or 
% right, correct/incorrect, etc.) Should work with most lewieWorld_B2AFC
% expDefs

% 2018-11-20: LEW created (plotCRF_stimulusHistory.m)
% 2019-05-20: Changed name; modified to more easily change condition comparison
% and adapt to varied contrast conditions across subjects 
% 2019-07-02: moved to git repository 'toupee'; extracted cell type chooser
% to separate function call

%% load experiment details

% close all;
% clear all;
clearvars -except allCRFs allCRFs_L allCRFs_R

expInfo.mouseName = 'LEW008';
expInfo.expDate = '2019-01-29';
expInfo.expNum = 1;
expInfo.expSeries = [1];

%% load data
expInfo = data.loadExpData(expInfo);

%% get event timings and wheel trajectories

[eventTimes, wheelTrajectories] = getEventTimes(expInfo, {'stimulusOnTimes' 'interactiveOnTimes' 'stimulusOffTimes'});

 %% load traces
[allFcell, ops] = loadCellData(expInfo);

%% align calcium traces to the event you want

% cut the trace into trial-by-trial traces, aligned to a particular event
event = 'prestimulusQuiescenceEndTimes';
% event = 'stimulusOnTimes';
% event = 'rewardOnTimes';
[alignedTraces, eventWindow] = getAlignedTraces(expInfo, allFcell, eventTimes, event);

% [alignedTraces, alignedSpikes, eventWindow] = alignSpikes(mouseName, expDate, expNum, expSeries, allFcell, eventTimes, ops, event);

%% select cells with the properties you want
cd('C:\Users\Wool\Documents\GitHub\toupee')
plotAll = chooseCellType('movleft', expInfo, allFcell, eventTimes);
%% initialize some data values
block = expInfo.block;

Fs = 15;

%event window
eventIdx = find(eventWindow == 0);

%stimulus response index
crfTime = 0.2;
vPreIdx = eventIdx - 1;
vPeriIdx = eventIdx + ceil(crfTime*Fs);

% contrasts to plot
contrasts = unique(block.events.contrastValues);
% contrasts = unique(sign(contrasts));
dirs = [-1 1];

repeatType = {'all','all'};
movementDir = {'all','all'};
movementTime = {'all','all'};
highRewardSide = {'left','right'};
responseType = {'all','all'};
rewardOutcome = {'all','all'};
pastStimulus = {'all','all'};
pastMovementDir = {'all','all'};
pastResponseType = {'all','all'};


%% initialize some plot values
zeroIdx = find(contrasts == 0);
walkup = length(contrasts) - zeroIdx;
walkback = zeroIdx - 1;
% allColors = [.25 0 0;.5 0 0 ;1 0 0;.8 .45 .45;.75 .75 .75;.55 .55 .55;.35 .35 .35;.15 .15 .15;0 0 0];
allColors = [0 0 .5; 0 0 1; 0 .4 1; .6 .8 1; .75 .75 .75; .8 .45 .45; 1 0 0; .5 0 0; .25 0 0];

zeroGray = find(allColors(:,1) == .75);
colors = allColors(zeroGray-walkback:zeroGray + walkup,:);
dirColors = [0 .4 1; 1 0 0];

% condition comparison
title1 = 'high contra';
title2 = 'high ipsi';
legend1 = 'contra';
legend2 = 'ipsi';


%% grand mean plot

figure;hold on
set(gcf,'Position',[114   400   1080   315])

allRTs = eventTimes(7).daqTime - eventTimes(1).daqTime;
shortRTs = find(allRTs >= 0 & allRTs <= 5);
beforeSwitchTrials = [126:175,301:350];
afterSwitchTrials = [176:225,351:400];

subplot(1,3,1);cla
subplot(1,3,2);cla
subplot(1,3,3);cla
clear grandMean grandSEM_UB grandSEM_LB maxCalcium minCalcium periMean_1 preMean_1 periUCI_1 periLCI_1
clear grandMean grandSEM_UB grandSEM_LB maxCalcium minCalcium periMean_2 preMean_2 periUCI_2 periLCI_2

for c = 1:length(contrasts)

    % select all the correct trials where the left side had a
    % higher reward
    [~, condIdx_1] = selectCondition(block, contrasts(c), eventTimes, ...
        repeatType{1}, movementDir{1}, movementTime{1}, highRewardSide{1}, responseType{1}, rewardOutcome{1}, ...
        pastStimulus{1}, pastMovementDir{1}, pastResponseType{1});
%     condIdx_1 = intersect(afterSwitchTrials, condIdx_1);
    
    grandValues = [];
    for iCell = 1:length(plotAll)
        rawValues = alignedTraces{plotAll(iCell,1)}.eventSpikes(condIdx_1,:,plotAll(iCell,2));
        normValues = nanmean(rawValues,1);% - baselineValues)/baselineValues;
        grandValues = [grandValues; normValues];
    end
    
    grandMean = nanmean(grandValues);
    grandSEM_UB = grandMean + nanstd(grandValues)/sqrt(length(grandValues));
    grandSEM_LB = grandMean - nanstd(grandValues)/sqrt(length(grandValues));

    % retrieve calcium
    maxCalcium(c,1) = max(grandSEM_UB);
    minCalcium(c,1) = min(grandSEM_LB);
    periMean_1(c) = grandMean(vPeriIdx);
    preMean_1(c) = 0;
    periUCI_1(c) = grandSEM_UB(vPeriIdx);
    periLCI_1(c) = grandSEM_LB(vPeriIdx);

    % plot calcium
    subplot(1,3,1);hold on
    plot1 = plot(eventWindow,grandMean);
    plot1ci = fill([eventWindow';flipud(eventWindow')],[grandSEM_LB';flipud(grandSEM_UB')],colors(c,:), 'LineStyle', 'none');
    alpha(0.2);
    set(plot1, 'LineStyle', '-', 'LineWidth',1.5,'Color',colors(c,:));
    title(title1);
    xlim([-.4 0.8]);
    box off
    ylabel('z-scored spikes')
    xlabel('time (s)')
    hold on;
    clear grandMean grandSEM_UB grandSEM_LB

    % select all the correct trials where the right side had a
    % higher reward
    [~, condIdx_2] = selectCondition(block, contrasts(c), eventTimes, ...
        repeatType{2}, movementDir{2}, movementTime{2}, highRewardSide{2}, responseType{2}, rewardOutcome{2}, ...
        pastStimulus{2}, pastMovementDir{2}, pastResponseType{2});
%     condIdx_2 = intersect(afterSwitchTrials, condIdx_2);
    % retrieve calcium
    grandValues = [];
    for iCell = 1:length(plotAll)
        rawValues = alignedTraces{plotAll(iCell,1)}.eventSpikes(condIdx_2,:,plotAll(iCell,2));
        normValues = nanmean(rawValues,1);% - baselineValues)/baselineValues;
        grandValues = [grandValues; normValues];
    end
    
    grandMean = nanmean(grandValues);
    grandSEM_UB = grandMean + nanstd(grandValues)/sqrt(length(grandValues));
    grandSEM_LB = grandMean - nanstd(grandValues)/sqrt(length(grandValues));

    % retrieve calcium
    maxCalcium(c,2) = max(grandSEM_UB);
    minCalcium(c,2) = min(grandSEM_LB);
    periMean_2(c) = grandMean(vPeriIdx);
    preMean_2(c) = 0;
    periUCI_2(c) = grandSEM_UB(vPeriIdx);
    periLCI_2(c) = grandSEM_LB(vPeriIdx);

    % plot calcium
    subplot(1,3,2);hold on
    plot2 = plot(eventWindow,grandMean);
    plot2ci = fill([eventWindow';flipud(eventWindow')],[grandSEM_LB';flipud(grandSEM_UB')],colors(c,:), 'LineStyle', 'none');
    alpha(0.2);
    set(plot2, 'LineStyle', '-', 'LineWidth',1.5,'Color',colors(c,:));
    title(title2);
    xlim([-.4 0.8]);
    box off
    ylabel('z-scored spikes')
    xlabel('time (s)')
    hold on;
    clear grandMean grandSEM_UB grandSEM_LB

    % plot calcium CRF (markers)
    subplot(1,3,3);hold on
    crf1c = plot(contrasts(c),periMean_1(c)-preMean_1(c),'ko');
    crf2c = plot(contrasts(c),periMean_2(c)-preMean_2(c),'ko');
    crf1cerror = line([contrasts(c) contrasts(c)],[periUCI_1(c)-preMean_1(c) periLCI_1(c)-preMean_1(c)]); 
    crf2cerror = line([contrasts(c) contrasts(c)],[periUCI_2(c)-preMean_2(c) periLCI_2(c)-preMean_2(c)]); 
    set([crf1c,crf2c],'LineStyle', 'none', 'LineWidth',1,'MarkerEdgeColor','none','MarkerFaceColor',colors(c,:),'MarkerSize',7)
    set([crf1cerror,crf2cerror],'LineStyle', '-', 'LineWidth',1,'Marker','none','Color',colors(c,:))
    %set(gca,'XDir','reverse')

end

%format subplot limits
subplot(1,3,1);
ylim([min([min(min(minCalcium)) 0]) max([max(max(maxCalcium)) .3])]);
crfLine1 = line([eventWindow(eventIdx) eventWindow(eventIdx)],[-.2 15]);
uistack(crfLine1,'bottom');
set(crfLine1,'LineStyle', '--', 'LineWidth',1,'Marker','none','Color',[.5 .5 .5]);
legend off
ax1 = gca;
ax1.TickDir = 'out';

subplot(1,3,2);
ylim([min([min(min(minCalcium)) 0]) max([max(max(maxCalcium)) .3])]);
crfLine2 = line([eventWindow(eventIdx) eventWindow(eventIdx)],[-.2 15]);
uistack(crfLine2,'bottom');
set(crfLine2,'LineStyle', '--', 'LineWidth',1,'Marker','none','Color',[.5 .5 .5]);
legend off
ax2 = gca;
ax2.TickDir = 'out';

% plot calcium CRF (lines)
subplot(1,3,3);hold on
crf1p = plot(contrasts,periMean_1-preMean_1,'k-');
crf2p = plot(contrasts,periMean_2-preMean_2,'k--');
uistack(crf1p,'bottom');
uistack(crf2p,'bottom');
title(strcat({'CRF (n = '},num2str(length(plotAll)),')'));
xlim([min(contrasts)-0.05 max(contrasts)+0.05]);
ylim([min([min(periLCI_1-preMean_1) min(periLCI_2-preMean_2)]) max([max(periUCI_1+preMean_1) max(periUCI_2+preMean_2) .2])]);
ylabel('z-scored spikes')
xlabel('contrast (%)')
legend([crf1p crf2p],{legend1,legend2},'Location','northeast');
xticks([-1:.25:1])
xticklabels({'100', '75', '50', '25','0','25','50','75','100'})
ax3 = gca;
ax3.TickDir = 'out';

%% cell-by-cell browser    

fig = figure(201);
hold on
set(fig,'Position',[114   400   1080   630])

allRTs = eventTimes(7).daqTime - eventTimes(1).daqTime;
shortRTs = find(allRTs <= 2);

k = 1;
max_k = length(plotAll);

while k <= max_k
        subplot(2,3,1);cla
        subplot(2,3,2);cla
        subplot(2,3,3);cla
        subplot(2,3,4);cla
        subplot(2,3,5);cla
        subplot(2,3,6);cla
        for c = 1:length(contrasts)
%             if contrasts(c) == 0 
%                 movementDir = {'cw','ccw'};
%             else
%                 movementDir = {'all','all'};
%             end
            % select all the correct trials where the left side had a
            % higher reward
            [~, condIdx_1] = selectCondition(block, contrasts(c), eventTimes, ...
                repeatType{1}, movementDir{1}, movementTime{1}, highRewardSide{1}, responseType{1}, rewardOutcome{1}, ...
                pastStimulus{1}, pastMovementDir{1}, pastResponseType{1});
            condIdx_1 = intersect(shortRTs, condIdx_1);

            % retrieve calcium
            meanCellCalcium = nanmean((alignedTraces{plotAll(k,1)}.eventSpikes(condIdx_1,:,plotAll(k,2))),1);
            semCellCalcium = nanstd((alignedTraces{plotAll(k,1)}.eventSpikes(condIdx_1,:,plotAll(k,2))))/sqrt(length(condIdx_1));
            upperCICalcium = meanCellCalcium + semCellCalcium;
            lowerCICalcium = meanCellCalcium - semCellCalcium;
            maxCalcium(c,1) = max(upperCICalcium);
            minCalcium(c,1) = min(lowerCICalcium);
            preMean_1(c) = meanCellCalcium(vPreIdx);
            periMean_1(c) = meanCellCalcium(vPeriIdx);
            periUCI_1(c) = upperCICalcium(vPeriIdx);
            periLCI_1(c) = lowerCICalcium(vPeriIdx);
            
            % plot calcium
            subplot(2,3,1);hold on
            plot1 = plot(eventWindow,meanCellCalcium);
            plot1ci = fill([eventWindow';flipud(eventWindow')],[lowerCICalcium';flipud(upperCICalcium')],colors(c,:), 'LineStyle', 'none');
            alpha(0.2);
            set(plot1, 'LineStyle', '-', 'LineWidth',1.5,'Color',colors(c,:));
            title(strcat({'Plane '},num2str(plotAll(k,1)),{', Cell '},num2str(plotAll(k,2)),{': '},title1));
            xlim([-.4 .8]);
            box off
            ylabel('z-scored spikes')
            xlabel('time (s)')
            hold on;
            clear meanCellCalcium baselineCalcium normCellCalcium

            % select all the correct trials where the right side had a
            % higher reward
            [~, condIdx_2] = selectCondition(block, contrasts(c), eventTimes, ...
                repeatType{2}, movementDir{2}, movementTime{2}, highRewardSide{2}, responseType{2}, rewardOutcome{2}, ...
                pastStimulus{2}, pastMovementDir{2}, pastResponseType{2});
            condIdx_2 = intersect(shortRTs, condIdx_2);
            
            % retrieve calcium
            meanCellCalcium = nanmean((alignedTraces{plotAll(k,1)}.eventSpikes(condIdx_2,:,plotAll(k,2))),1);
            semCellCalcium = nanstd((alignedTraces{plotAll(k,1)}.eventSpikes(condIdx_2,:,plotAll(k,2))))/sqrt(length(condIdx_2));
            upperCICalcium = meanCellCalcium + semCellCalcium;
            lowerCICalcium = meanCellCalcium - semCellCalcium;
            maxCalcium(c,2) = max(upperCICalcium);
            minCalcium(c,2) = min(lowerCICalcium);
            preMean_2(c) = meanCellCalcium(vPreIdx);
            periMean_2(c) = meanCellCalcium(vPeriIdx);
            periUCI_2(c) = upperCICalcium(vPeriIdx);
            periLCI_2(c) = lowerCICalcium(vPeriIdx);
            
            % plot calcium
            subplot(2,3,2);hold on
            plot2 = plot(eventWindow,meanCellCalcium);
            plot2ci = fill([eventWindow';flipud(eventWindow')],[lowerCICalcium';flipud(upperCICalcium')],colors(c,:), 'LineStyle', 'none');
            alpha(0.2);
            set(plot2, 'LineStyle', '-', 'LineWidth',1.5,'Color',colors(c,:));
            title(strcat({'Plane '},num2str(plotAll(k,1)),{', Cell '},num2str(plotAll(k,2)),{': '},title2));
            xlim([-.4 .8]);
            box off
            ylabel('z-scored spikes')
            xlabel('time (s)')
            hold on;
            clear meanCellCalcium baselineCalcium normCellCalcium
                        
            % plot calcium CRF (markers)
            subplot(2,3,3);hold on
            crf1c = plot(contrasts(c),periMean_1(c)-preMean_1(c),'ko');
            crf2c = plot(contrasts(c),periMean_2(c)-preMean_2(c),'ko');
            crf1cerror = line([contrasts(c) contrasts(c)],[periUCI_1(c)-preMean_1(c) periLCI_1(c)-preMean_1(c)]); 
            crf2cerror = line([contrasts(c) contrasts(c)],[periUCI_2(c)-preMean_2(c) periLCI_2(c)-preMean_2(c)]); 
            set([crf1c,crf2c],'LineStyle', 'none', 'LineWidth',1,'MarkerEdgeColor','none','MarkerFaceColor',colors(c,:),'MarkerSize',7)
            set([crf1cerror,crf2cerror],'LineStyle', '-', 'LineWidth',1,'Marker','none','Color',colors(c,:))
            %set(gca,'XDir','reverse')
                        
        end
        
    %format subplot limits
    subplot(2,3,1);
    ylim([min([min(min(minCalcium)) 0]) max([max(max(maxCalcium)) .4])]);
    crfLine1 = line([eventWindow(eventIdx) eventWindow(eventIdx)],[-1 15]);
    uistack(crfLine1,'bottom');
    set(crfLine1,'LineStyle', '--', 'LineWidth',1,'Marker','none','Color',[.5 .5 .5]);
    legend off
    ax1 = gca;
    ax1.TickDir = 'out';

    subplot(2,3,2);
    ylim([min([min(min(minCalcium)) 0]) max([max(max(maxCalcium)) .4])]);
    crfLine2 = line([eventWindow(eventIdx) eventWindow(eventIdx)],[-1 15]);
    uistack(crfLine2,'bottom');
    set(crfLine2,'LineStyle', '--', 'LineWidth',1,'Marker','none','Color',[.5 .5 .5]);
    legend off
    ax2 = gca;
    ax2.TickDir = 'out';

    % plot calcium CRF (lines)
    subplot(2,3,3);hold on
    crf1p = plot(contrasts,periMean_1-preMean_1,'k-');
    crf2p = plot(contrasts,periMean_2-preMean_2,'k--');
    uistack(crf1p,'bottom');
    uistack(crf2p,'bottom');
    title(strcat({'Plane '},num2str(plotAll(k,1)),{', Cell '},num2str(plotAll(k,2)),{': CRF'}));
    xlim([min(contrasts) max(contrasts)]*1.05);
    ylim([min([min(periLCI_1-preMean_1) min(periLCI_2-preMean_2)]) max([max(periUCI_1-preMean_1) max(periUCI_2-preMean_2) 1])]);
    ylabel('z-scored spikes')
    xlabel('contrast (%)')
    xticks([-1:.2:1])
    xticklabels({'100', '80', '60', '40', '20','0','20','40','60','80','100'})
    ax3 = gca;
    ax3.TickDir = 'out';
    
    for s = 1:length(dirs)
            movementDir_dirs = {'cw','ccw'};

            % select all the correct trials where the left side had a
            % higher reward
            [~, condIdx_1] = selectCondition(block, contrasts, eventTimes, ...
                repeatType{1}, movementDir_dirs{s}, movementTime{1}, highRewardSide{1}, 'correct', rewardOutcome{1}, ...
                pastStimulus{1}, pastMovementDir{1}, pastResponseType{1});

            % retrieve calcium
            meanCellCalcium_sign = nanmean((alignedTraces{plotAll(k,1)}.eventSpikes(condIdx_1,:,plotAll(k,2))),1);
            semCellCalcium_sign = nanstd((alignedTraces{plotAll(k,1)}.eventSpikes(condIdx_1,:,plotAll(k,2))))/sqrt(length(condIdx_1));
            upperCICalcium_sign = meanCellCalcium_sign + semCellCalcium_sign;
            lowerCICalcium_sign = meanCellCalcium_sign - semCellCalcium_sign;
            maxCalcium_sign(s,1) = max(upperCICalcium_sign);
            minCalcium_sign(s,1) = min(lowerCICalcium_sign);
            preMean_1_sign(s) = meanCellCalcium_sign(vPreIdx);
            periMean_1_sign(s) = meanCellCalcium_sign(vPeriIdx);
            periUCI_1_sign(s) = upperCICalcium_sign(vPeriIdx);
            periLCI_1_sign(s) = lowerCICalcium_sign(vPeriIdx);
            
            % plot calcium
            subplot(2,3,4);hold on
            plot1 = plot(eventWindow,meanCellCalcium_sign);
            plot1ci = fill([eventWindow';flipud(eventWindow')],[lowerCICalcium_sign';flipud(upperCICalcium_sign')],dirColors(s,:), 'LineStyle', 'none');
            alpha(0.2);
            set(plot1, 'LineStyle', '-', 'LineWidth',1.5,'Color',dirColors(s,:));
            title(strcat({'Plane '},num2str(plotAll(k,1)),{', Cell '},num2str(plotAll(k,2)),{': '},title1));
            xlim([-.4 .8]);
            box off
            ylabel('z-scored spikes')
            xlabel('time (s)')
            hold on;
            clear meanCellCalcium baselineCalcium normCellCalcium

            % select all the correct trials where the right side had a
            % higher reward
            [~, condIdx_2] = selectCondition(block, contrasts, eventTimes, ...
                repeatType{2}, movementDir_dirs{s}, movementTime{2}, highRewardSide{2}, 'correct', rewardOutcome{2}, ...
                pastStimulus{2}, pastMovementDir{2}, pastResponseType{2});
            
            % retrieve calcium
            meanCellCalcium_sign = nanmean((alignedTraces{plotAll(k,1)}.eventSpikes(condIdx_2,:,plotAll(k,2))),1);
            semCellCalcium_sign = nanstd((alignedTraces{plotAll(k,1)}.eventSpikes(condIdx_2,:,plotAll(k,2))))/sqrt(length(condIdx_2));
            upperCICalcium_sign = meanCellCalcium_sign + semCellCalcium_sign;
            lowerCICalcium_sign = meanCellCalcium_sign - semCellCalcium_sign;
            maxCalcium_sign(s,2) = max(upperCICalcium_sign);
            minCalcium_sign(s,2) = min(lowerCICalcium_sign);
            preMean_2_sign(s) = meanCellCalcium_sign(vPreIdx);
            periMean_2_sign(s) = meanCellCalcium_sign(vPeriIdx);
            periUCI_2_sign(s) = upperCICalcium_sign(vPeriIdx);
            periLCI_2_sign(s) = lowerCICalcium_sign(vPeriIdx);
            
            % plot calcium
            subplot(2,3,5);hold on
            plot2 = plot(eventWindow,meanCellCalcium_sign);
            plot2ci = fill([eventWindow';flipud(eventWindow')],[lowerCICalcium_sign';flipud(upperCICalcium_sign')],dirColors(s,:), 'LineStyle', 'none');
            alpha(0.2);
            set(plot2, 'LineStyle', '-', 'LineWidth',1.5,'Color',dirColors(s,:));
            title(strcat({'Plane '},num2str(plotAll(k,1)),{', Cell '},num2str(plotAll(k,2)),{': '},title2));
            xlim([-.4 .8]);
            box off
            ylabel('z-scored spikes')
            xlabel('time (s)')
            hold on;
            clear meanCellCalcium baselineCalcium normCellCalcium
                        
            % plot calcium CRF (markers)
            subplot(2,3,6);hold on
            crf1c = plot(dirs(s),periMean_1_sign(s)-preMean_1_sign(s),'ko');
            crf2c = plot(dirs(s),periMean_2_sign(s)-preMean_2_sign(s),'ko');
            crf1cerror = line([dirs(s) dirs(s)],[periUCI_1_sign(s)-preMean_1_sign(s) periLCI_1_sign(s)-preMean_1_sign(s)]); 
            crf2cerror = line([dirs(s) dirs(s)],[periUCI_2_sign(s)-preMean_2_sign(s) periLCI_2_sign(s)-preMean_2_sign(s)]); 
            set([crf1c,crf2c],'LineStyle', 'none', 'LineWidth',1,'MarkerEdgeColor','none','MarkerFaceColor',dirColors(s,:),'MarkerSize',7)
            set([crf1cerror,crf2cerror],'LineStyle', '-', 'LineWidth',1,'Marker','none','Color',dirColors(s,:))
            %set(gca,'XDir','reverse')
                        
    end
        
    %format subplot limits
        subplot(2,3,4);
        ylim([min([min(min(minCalcium_sign)) 0]) max([max(max(maxCalcium_sign)) .4])]);
        crfLine1 = line([eventWindow(eventIdx) eventWindow(eventIdx)],[-1 15]);
        uistack(crfLine1,'bottom');
        set(crfLine1,'LineStyle', '--', 'LineWidth',1,'Marker','none','Color',[.5 .5 .5]);
        legend off
        ax1 = gca;
        ax1.TickDir = 'out';

            subplot(2,3,5);
            ylim([min([min(min(minCalcium_sign)) 0]) max([max(max(maxCalcium_sign)) .4])]);
            crfLine2 = line([eventWindow(eventIdx) eventWindow(eventIdx)],[-1 15]);
            uistack(crfLine2,'bottom');
            set(crfLine2,'LineStyle', '--', 'LineWidth',1,'Marker','none','Color',[.5 .5 .5]);
            legend off
            ax2 = gca;
            ax2.TickDir = 'out';

        % plot calcium CRF (lines)
        subplot(2,3,6);hold on
        crf1p = plot(dirs,periMean_1_sign-preMean_1_sign,'k-');
        crf2p = plot(dirs,periMean_2_sign-preMean_2_sign,'k--');
        uistack(crf1p,'bottom');
        uistack(crf2p,'bottom');
        title(strcat({'Plane '},num2str(plotAll(k,1)),{', Cell '},num2str(plotAll(k,2)),{': CRF'}));
        xlim([min(dirs) max(dirs)]*2);
        ylim([min([min(periLCI_1_sign-preMean_1_sign) min(periLCI_2_sign-preMean_2_sign)])*1.1 max([max(periUCI_1_sign-preMean_1_sign) max(periUCI_2_sign-preMean_2_sign) 1])]);
        ylabel('z-scored spikes')
        xlabel('choice')
        xticks([-1:1])
        xticklabels({'left (cw)', '','right (ccw)'})
        ax3 = gca;
        ax3.TickDir = 'out';


    was_a_key = waitforbuttonpress;
    if was_a_key && strcmp(get(fig, 'CurrentKey'), 'leftarrow')
      k = max(1, k - 1);
    elseif was_a_key && strcmp(get(fig, 'CurrentKey'), 'rightarrow')
      k = min(max_k, k + 1);
    end
 
end

%%
plotBlockScatter(plotAll, block, eventTimes, eventWindow, alignedTraces)