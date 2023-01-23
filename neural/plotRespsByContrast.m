function plotRespsByContrast
% this function generates an n x m plot of event-triggered average
% responses for a single cell or a subset of cells (mean of means), 
% where n = number of moveDirs (usually 2) and m = number of contrast 
% conditions (usually 9)
%% make a list of mice/experiments you want to analyze

mouseList = { ...
    {'LEW031'}...
    };

expList = { ...
    {'2020-02-13',1,[1]}};

%% load all the experiments into expInfo

expInfo = initExpInfo(mouseList,expList);


%% process the usual data
% the script knows to loop over all the experiments you listed above
% this will take a while but the command line will print progress

[expInfo, neuralData, behavioralData] = processExperiment(expInfo);
eyeData = getEyeData(expInfo);
[eyeData] = alignFace(expInfo, eyeData, behavioralData);

%% extract some variables

contrasts = getUniqueContrasts(expInfo);
eventWindow = neuralData(1).eta.eventWindow;

%% initialize the plot

zeroIdx = find(contrasts == 0);
walkup = length(contrasts) - zeroIdx;
walkback = zeroIdx - 1;
% allColors = [.25 0 0;.5 0 0 ;1 0 0;.8 .45 .45;.75 .75 .75;.55 .55 .55;.35 .35 .35;.15 .15 .15;0 0 0];
allColors = [0 0 .5; 0 0 1; 0 .4 1; .6 .8 1; .75 .75 .75; .8 .45 .45; 1 0 0; .5 0 0; .25 0 0];

zeroGray = find(allColors(:,1) == .75);
colors = allColors(zeroGray-walkback:zeroGray + walkup,:);
dirColors = [0 .4 1; 1 0 0];

%% plot responses to contrast x block

% whichNeuron = 41;
% whichNeuron = find(neuralData.stats.bfcH(:,strcmp(neuralData.stats.labels,'stim')) > 0);
whichNeuron = 1:length(neuralData(1).stats.bfcH);
cellResps = nanmean(neuralData.eta.alignedResps{1}(:,:,whichNeuron),3);

trialTypes = getTrialTypes(expInfo, behavioralData, 'late');
whichTrials = trialTypes.intVar.all.contrast_block;

[baselineResps, stimResps, pmovResps, movResps, rewResps, preCueResps] = getEpochResps(neuralData.eta);

for c = 1:size(whichTrials,1)
    hL = whichTrials{c,1};
    hR = whichTrials{c,2};
    
    R_hL(c,:) = nanmean(cellResps(hL,:),1);
    R_hR(c,:) = nanmean(cellResps(hR,:),1);
    
    CRF_hL(c) = nanmean(nanmean(stimResps(hL,whichNeuron),1));
    CRF_hR(c) = nanmean(nanmean(stimResps(hR,whichNeuron),1));
    
end   

figure;
set(gcf,'position',[120 560 1580 210])
for sp = 1:length(contrasts)*2
    subplot(2,length(contrasts),sp)
    hold on
    box off
    set(gca,'tickdir','out')
    if sp == 1
        ylabel('left choices')
    elseif sp == length(contrasts)+1
        ylabel('right choices')
    end
end

maxY = max(max([R_hL; R_hR]));
minY = min(min([R_hL; R_hR]));

for c = 1:size(whichTrials,1)
    subplot(1,length(contrasts)+2,c)
    plot(eventWindow, R_hL(c,:), 'Color',colors(c,:),'LineStyle','-','LineWidth',2);
    hold on
    plot(eventWindow, R_hR(c,:), 'Color',colors(c,:),'LineStyle',':','LineWidth',2);
    ylim([minY maxY]);
    xlim([-.5 1.5]);
    box off
    set(gca,'tickdir','out')
end

subplot(1,length(contrasts)+2,length(contrasts)+1:length(contrasts)+2)
 plot(contrasts,CRF_hL,'Color','k','LineStyle','-','LineWidth',1.5);
 hold on;
 plot(contrasts,CRF_hR,'Color','k','LineStyle',':','LineWidth',1.5);
 
for c = 1:length(contrasts)
    plot(contrasts(c),CRF_hL(c), 'ko', 'MarkerFaceColor',colors(c,:),'MarkerEdgeColor','w');
    plot(contrasts(c),CRF_hR(c), 'ko', 'MarkerFaceColor',colors(c,:),'MarkerEdgeColor','w');
end

box off
set(gca,'tickdir','out')
xlim([-1.1 1.1]);

%% plot responses to contrast x movement x block

whichNeuron = 623;
% whichNeuron = find(neuralData.stats.bfcH(:,strcmp(neuralData.stats.labels,'leftMov')) > 0);
% whichNeuron = 1:length(neuralData(1).stats.bfcH);
cellResps = nanmean(neuralData.eta.alignedResps{2}(:,:,whichNeuron),3);

trialTypes = getTrialTypes(expInfo, behavioralData, 'late');
whichTrials = trialTypes.intVar.all.contrast_direction_block;

for c = 1:size(whichTrials,1)
    mLhL = whichTrials{c,1,1};
    mRhL = whichTrials{c,2,1};
    mLhR = whichTrials{c,1,2};
    mRhR = whichTrials{c,2,2};
    
    R_mLhL(c,:) = nanmean(cellResps(mLhL,:),1);
    R_mRhL(c,:) = nanmean(cellResps(mRhL,:),1);
    R_mLhR(c,:) = nanmean(cellResps(mLhR,:),1);
    R_mRhR(c,:) = nanmean(cellResps(mRhR,:),1);
    
end   

figure;
set(gcf,'position',[120 560 1580 420])
for sp = 1:length(contrasts)*2
    subplot(2,length(contrasts),sp)
    hold on
    box off
    set(gca,'tickdir','out')
    if sp == 1
        ylabel('left choices')
    elseif sp == length(contrasts)+1
        ylabel('right choices')
    end
end

maxY = max(max([R_mLhL; R_mRhL; R_mLhR; R_mRhR]));
minY = min(min([R_mLhL; R_mRhL; R_mLhR; R_mRhR]));

for c = 1:size(whichTrials,1)
    subplot(2,length(contrasts),c)
    plot(eventWindow, R_mLhL(c,:), 'Color',colors(c,:),'LineStyle','-','LineWidth',2);
    hold on
    plot(eventWindow, R_mLhR(c,:), 'Color',colors(c,:),'LineStyle',':','LineWidth',2);
    ylim([minY maxY]);
    xlim([-1 1]);
    subplot(2,length(contrasts),c+length(contrasts))
    plot(eventWindow, R_mRhL(c,:), 'Color',colors(c,:),'LineStyle','-','LineWidth',2);
    hold on
    plot(eventWindow, R_mRhR(c,:), 'Color',colors(c,:),'LineStyle',':','LineWidth',2);
    ylim([minY maxY]);
    xlim([-1 1]);
end
% 
% for c = 3:size(whichTrials,1)-2
%     subplot(1,length(contrasts)-4,c-2)
%     plot(eventWindow, R_mLhL(c,:), 'Color',[0 .4 1],'LineStyle','-','LineWidth',2);
%     hold on
%     plot(eventWindow, R_mLhR(c,:), 'Color',[0 .4 1],'LineStyle',':','LineWidth',2);
%     plot(eventWindow, R_mRhL(c,:), 'Color','r','LineStyle','-','LineWidth',2);
%     plot(eventWindow, R_mRhR(c,:), 'Color','r','LineStyle',':','LineWidth',2);
%     line([0 0],[minY maxY],'Color',[.5 .5 .5],'LineStyle','--')
%     ylim([minY maxY]);
%     xlim([-1 1]);
%     prettyPlot(gca)
%     title(num2str(contrasts(c)))
% end
%     
%% plot responses to contrast x movement
plotCells = 1:length(neuralData(1).stats.bfcH);
fig = figure;
set(gcf,'position',[120 560 1580 180])

for sp = 1:length(contrasts)
    subplot(1,length(contrasts),sp)
    hold on
    box off
    set(gca,'tickdir','out')
end

if ~exist('k') == 1
    k = 1;
end

max_k = length(plotCells);

while k <= max_k
    
   
    
% whichNeuron = plotCells(k);
% whichNeuron = find(neuralData.stats.bfcH(:,strcmp(neuralData.stats.labels,'stim')) > 0);
whichNeuron = 1:length(neuralData(1).stats.hTest);
% whichNeuron = plotCells;
% whichNeuron = 41;
% cellResps = nanmedian(neuralData(1).eta.alignedResps{1}(:,:,whichNeuron),3);

% for iT = 1:size(neuralData.eta.alignedResps{ 2},1)
%     relResps(iT,:,:) = neuralData.eta.alignedResps{1}(iT,:,whichNeuron)-reshape(repmat(baselineResps(iT,whichNeuron),41,1),[1 41 sum(whichNeuron)]);
% end
% cellResps = nanmean(relResps,3);
trialTypes = getTrialTypes(expInfo(1), behavioralData(1), 'late');
whichTrials = trialTypes.intVar.all.contrast_direction;
% whichTrials = trialTypes.intVar.cb2D.contrast_direction;

clear R_mL R_mR vel_mR vel_mL

 %clear subplots
    for c = 1:size(whichTrials,1)
    subplot(1,length(contrasts),c)
    cla;
    end

for c = 1:size(whichTrials,1)
    mL = whichTrials{c,1};
    mR = whichTrials{c,2};
    
    R_mL(c,:) = nanmean(nanmean(neuralData(1).eta.alignedResps{1}(mL,:,whichNeuron),1),3);
    R_mR(c,:) = nanmean(nanmean(neuralData(1).eta.alignedResps{1}(mR,:,whichNeuron),1),3);
    vel_mL(c,:) = nanmean(behavioralData(1).wheelMoves.epochs(5).peakVel(mL));
    vel_mR(c,:) = nanmean(behavioralData(1).wheelMoves.epochs(5).peakVel(mR));
    
end   



maxY = max(max([R_mL; R_mR]))*1.1;
minY = min(min([R_mL; R_mR]));
maxY = max(max([vel_mL; vel_mR]))*1.1;
minY = min(min([vel_mL; vel_mR]));

for c = 1:size(whichTrials,1)
    subplot(1,length(contrasts),c)
    line([0 0],[-.25 10],'Color',[.5 .5 .5],'LineStyle',':');
    plot(eventWindow, vel_mL(c,:), 'Color',colors(c,:),'LineStyle','-','LineWidth',2);
    hold on
    plot(eventWindow, vel_mR(c,:), 'Color',colors(c,:),'LineStyle',':','LineWidth',2);
    ylim([minY maxY]);
    xlim([-2 2]);
    title(num2str(contrasts(c)*100),'Color',colors(c,:));
    if c == ceil(length(contrasts)/2)
        xlabel('Time – stimOn (s)')
    end
    if c == 1
        ylabel('Activity')
    end
end


% for c = 1:size(whichTrials,1)
%     mL = whichTrials{c,1};
%     mR = whichTrials{c,2};
%     
%     B_mL(c,:) = nanmean(baselineResps(mL,whichNeuron),1);
%     B_mR(c,:) = nanmean(baselineResps(mR,whichNeuron),1);    
% end   

was_a_key = waitforbuttonpress;
    
    if was_a_key && strcmp(get(fig, 'CurrentKey'), 'leftarrow')
      k = max(1, k - 1);
    elseif was_a_key && strcmp(get(fig, 'CurrentKey'), 'rightarrow')
      k = min(max_k, k + 1);
    elseif was_a_key && strcmp(get(fig, 'CurrentKey'), 'return')
        disp(strcat({'k = '},num2str(k)))
        figName = strcat('contrast_dir',expInfo.mouseName,'_',expInfo.expDate,'_cell_',num2str(plotCells(k)));
        printfig(gcf, figName)
        break
    elseif was_a_key && strcmp(get(fig, 'CurrentKey'), 'escape')
        disp(strcat({'k = '},num2str(k)))
        figName = strcat('contrast_dir',expInfo.mouseName,'_',expInfo.expDate,'_cell_',num2str(plotCells(k)));
        close(fig)
        break
    end
end

%% plot responses to contrast x movement x previous choice

% whichNeuron = 9;
% whichNeuron = neuralData(1).stats.bfcH(:,5) == 1;
whichNeuron = 1:length(neuralData(1).stats.bfcH);
cellResps = nanmean(neuralData.eta.alignedResps{2}(:,:,whichNeuron),3);

trialTypes = getTrialTypes(expInfo, behavioralData, 'late');
whichTrials = trialTypes.intVar.all.contrast_direction_block;

for c = 1:size(whichTrials,1)
    [~, mLpL] = selectCondition(expInfo, contrasts(c), behavioralData, ...
        initTrialConditions('movementDir','cw','movementTime','late','pastMovementDir','cw'));
    [~, mRpL] = selectCondition(expInfo, contrasts(c), behavioralData, ...
        initTrialConditions('movementDir','ccw','movementTime','late','pastMovementDir','cw'));
    [~, mLpR] = selectCondition(expInfo, contrasts(c), behavioralData, ...
        initTrialConditions('movementDir','cw','movementTime','late','pastMovementDir','ccw'));
    [~, mRpR] = selectCondition(expInfo, contrasts(c), behavioralData, ...
        initTrialConditions('movementDir','ccw','movementTime','late','pastMovementDir','ccw'));
    
    R_mLpL(c,:) = nanmean(cellResps(mLpL,:),1);
    R_mRpL(c,:) = nanmean(cellResps(mRpL,:),1);
    R_mLpR(c,:) = nanmean(cellResps(mLpR,:),1);
    R_mRpR(c,:) = nanmean(cellResps(mRpR,:),1);
    
end   

figure;
set(gcf,'position',[120 560 1580 420])
for sp = 1:length(contrasts)*2
    subplot(2,length(contrasts),sp)
    hold on
    box off
    set(gca,'tickdir','out')
    if sp == 1
        ylabel('left choices')
    elseif sp == length(contrasts)+1
        ylabel('right choices')
    end
end

maxY = max(max([R_mLpL; R_mRpL; R_mLpR; R_mRpR]));
minY = min(min([R_mLpL; R_mRpL; R_mLpR; R_mRpR]));

for c = 1:size(whichTrials,1)
    subplot(2,length(contrasts),c)
    plot(eventWindow, R_mLpL(c,:), 'Color',colors(c,:),'LineStyle','-','LineWidth',2);
    hold on
    plot(eventWindow, R_mLpR(c,:), 'Color',colors(c,:),'LineStyle',':','LineWidth',2);
    ylim([minY maxY]);
    xlim([-1 1]);
    subplot(2,length(contrasts),c+length(contrasts))
    plot(eventWindow, R_mRpL(c,:), 'Color',colors(c,:),'LineStyle','-','LineWidth',2);
    hold on
    plot(eventWindow, R_mRpR(c,:), 'Color',colors(c,:),'LineStyle',':','LineWidth',2);
    ylim([minY maxY]);
    xlim([-1 1]);
end

%% plot responses to contrast

whichNeuron = 768;
% whichNeuron = neuralData(1).stats.bfcH(:,2) == 1 ;
% whichNeuron = 1:length(neuralData(1).stats.bfcH);
cellResps = nanmean(neuralData(1).eta.alignedResps{1}(:,:,whichNeuron),3);
clear R_m vel_m
% for iT = 1:size(neuralData.eta.alignedResps{ 2},1)
%     relResps(iT,:,:) = neuralData.eta.alignedResps{1}(iT,:,whichNeuron)-reshape(repmat(baselineResps(iT,whichNeuron),41,1),[1 41 sum(whichNeuron)]);
% end
% cellResps = nanmean(relResps,3);
trialTypes = getTrialTypes(expInfo(1), behavioralData(1), 'late');
whichTrials = trialTypes.singleVar.contrast;

for c = 1:size(whichTrials,1)
    m = whichTrials{c,1};
    
    R_m(c,:) = nanmean(cellResps(m,:),1);
    vel_m(c,:) = nanmean(behavioralData(1).wheelMoves.epochs(5).peakVel(m));
    
end   

figure;
set(gcf,'position',[120 560 1580 180])
for sp = 1:length(contrasts)
    subplot(1,length(contrasts),sp)
    hold on
    box off
    set(gca,'tickdir','out')
end

maxY = max(max([R_m; R_m]))*1.1;
minY = min(min([R_m; R_m]))*1.1;

for c = 1:size(whichTrials,1)
    subplot(1,length(contrasts),c)
    line([0 0],[-.25 10],'Color',[.5 .5 .5],'LineStyle',':');
    plot(eventWindow, R_m(c,:), 'Color',colors(c,:),'LineStyle','-','LineWidth',2);
    hold on
    ylim([minY maxY]);
    xlim([-.5 2]);
    title(num2str(contrasts(c)*100),'Color',colors(c,:));
    if c == ceil(length(contrasts)/2)
        xlabel('Time – stimOn (s)')
    end
    if c == 1
        ylabel('Activity')
    end
    
end

%% plot responses to contrast x reward amount

whichNeuron = 773;
% whichNeuron = neuralData(1).stats.bfcH(:,8) == 1 ;
% whichNeuron = 1:length(neuralData(1).stats.bfcH);
ETA = 2;

cellResps = nanmean(neuralData(1).eta.alignedResps{ETA}(:,:,whichNeuron),3);

trialTypes = getTrialTypes(expInfo(1), behavioralData(1), 'late');
whichTrials = trialTypes.singleVar.contrast;
contrasts = getUniqueContrasts(expInfo);
clear R_m vel_m
for c = 1:size(whichTrials,1)
    for r = 1:3
        m = intersect(whichTrials{c,1},trialTypes.singleVar.reward{r});
        R_m(c,r,:) = nanmean(cellResps(m,:),1);
        vel_m(c,r,:) = nanmean(behavioralData(1).eta.alignedVels{ETA}(m,:));
    end
    
end   

figure;
set(gcf,'position',[120 560 1580 180])
for sp = 1:length(contrasts)
    subplot(1,length(contrasts),sp)
    hold on
    box off
    set(gca,'tickdir','out')
end

maxY = max(max(max([R_m])))*1.1;
minY = min(min(min([R_m])))*1.1;
rcol = [0 .5 0; .5 1 .5; .7 0.1 0];

for c = 1:size(whichTrials,1)
    subplot(1,length(contrasts),c)
    line([0 0],[-.25 10],'Color',[.5 .5 .5],'LineStyle',':');
    plot(eventWindow, squeeze(R_m(c,1,:)), 'Color',rcol(1,:),'LineStyle','-','LineWidth',2);
    hold on
    plot(eventWindow, squeeze(R_m(c,2,:)), 'Color',rcol(2,:),'LineStyle','-','LineWidth',2);
    plot(eventWindow, squeeze(R_m(c,3,:)), 'Color',rcol(3,:),'LineStyle','-','LineWidth',2);
    ylim([minY maxY]);
    xlim([-1 1]);
    title(num2str(contrasts(c)*100),'Color','k');
    if c == ceil(length(contrasts)/2)
        xlabel('Time – stimOn (s)')
    end
    if c == 1
        ylabel('Activity')
    end
    
end

figure;
set(gcf,'position',[120 560 1580 180])
for sp = 1:length(contrasts)
    subplot(1,length(contrasts),sp)
    hold on
    box off
    set(gca,'tickdir','out')
end

maxY = max(max(max([abs(vel_m)])))*1.1;
minY = min(min(min([abs(vel_m)])))*1.1;
rcol = [0 .5 0; .5 1 .5; .7 0.1 0];

for c = 1:size(whichTrials,1)
    subplot(1,length(contrasts),c)
    line([0 0],[-.25 100],'Color',[.5 .5 .5],'LineStyle',':');
    plot(behavioralData.eta.eventWindow, squeeze(abs(vel_m(c,1,:))), 'Color',rcol(1,:),'LineStyle','-','LineWidth',2);
    hold on
    plot(behavioralData.eta.eventWindow, squeeze(abs(vel_m(c,2,:))), 'Color',rcol(2,:),'LineStyle','-','LineWidth',2);
    plot(behavioralData.eta.eventWindow, squeeze(abs(vel_m(c,3,:))), 'Color',rcol(3,:),'LineStyle','-','LineWidth',2);
    ylim([minY maxY]);
    xlim([-1 1]);
    title(num2str(contrasts(c)*100),'Color','k');
    if c == ceil(length(contrasts)/2)
        xlabel('Time – stimOn (s)')
    end
    if c == 1
        ylabel('Activity')
    end
    
end

% for c = 1:size(whichTrials,1)
%     mL = whichTrials{c,1};
%     mR = whichTrials{c,2};
%     
%     B_mL(c,:) = nanmean(baselineResps(mL,whichNeuron),1);
%     B_mR(c,:) = nanmean(baselineResps(mR,whichNeuron),1);    
% end       