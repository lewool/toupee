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

%% plot responses to contrast x movement x block

% whichNeuron = 9;
whichNeuron = neuralData(1).stats.bfcH(:,5) == 1;
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

%% plot responses to contrast x movement

whichNeuron = 1;
whichNeuron = neuralData(1).stats.bfcH(:,6) == 1 ;
% whichNeuron = 1:length(neuralData(1).stats.bfcH);
cellResps = nanmean(neuralData(1).eta.alignedResps{1}(:,:,whichNeuron),3);

% for iT = 1:size(neuralData.eta.alignedResps{ 2},1)
%     relResps(iT,:,:) = neuralData.eta.alignedResps{1}(iT,:,whichNeuron)-reshape(repmat(baselineResps(iT,whichNeuron),41,1),[1 41 sum(whichNeuron)]);
% end
% cellResps = nanmean(relResps,3);
trialTypes = getTrialTypes(expInfo(1), behavioralData(1), 'late');
whichTrials = trialTypes.intVar.all.contrast_direction;
clear R_mL R_mR vel_mR vel_mL
for c = 1:size(whichTrials,1)
    mL = whichTrials{c,1};
    mR = whichTrials{c,2};
    
    R_mL(c,:) = nanmean(cellResps(mL,:),1);
    R_mR(c,:) = nanmean(cellResps(mR,:),1);
    vel_mL(c,:) = nanmean(behavioralData(1).wheelMoves.epochs(5).peakVel(mL));
    vel_mR(c,:) = nanmean(behavioralData(1).wheelMoves.epochs(5).peakVel(mR));
    
end   

figure;
set(gcf,'position',[120 560 1580 180])
for sp = 1:length(contrasts)
    subplot(1,length(contrasts),sp)
    hold on
    box off
    set(gca,'tickdir','out')
end

maxY = max(max([R_mL; R_mR]))*1.1;
minY = min(min([R_mL; R_mR]))*.9;

for c = 1:size(whichTrials,1)
    subplot(1,length(contrasts),c)
    line([0 0],[-.25 10],'Color',[.5 .5 .5],'LineStyle',':');
    plot(eventWindow, R_mL(c,:), 'Color',colors(c,:),'LineStyle','-','LineWidth',2);
    hold on
    plot(eventWindow, R_mR(c,:), 'Color',colors(c,:),'LineStyle',':','LineWidth',2);
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


% for c = 1:size(whichTrials,1)
%     mL = whichTrials{c,1};
%     mR = whichTrials{c,2};
%     
%     B_mL(c,:) = nanmean(baselineResps(mL,whichNeuron),1);
%     B_mR(c,:) = nanmean(baselineResps(mR,whichNeuron),1);    
% end   

%% plot responses to contrast x movement x previous choice

% whichNeuron = 9;
whichNeuron = neuralData(1).stats.bfcH(:,6) == 1;
% whichNeuron = 1:length(neuralData(1).stats.bfcH);
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

whichNeuron = 1;
% whichNeuron = neuralData(1).stats.bfcH(:,2) == 1 ;
% whichNeuron = 1:length(neuralData(1).stats.bfcH);
cellResps = nanmean(neuralData(1).eta.alignedResps{1}(:,:,whichNeuron),3);

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


% for c = 1:size(whichTrials,1)
%     mL = whichTrials{c,1};
%     mR = whichTrials{c,2};
%     
%     B_mL(c,:) = nanmean(baselineResps(mL,whichNeuron),1);
%     B_mR(c,:) = nanmean(baselineResps(mR,whichNeuron),1);    
% end       