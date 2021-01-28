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
    {'2020-02-03',1,[1]}};

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
eventWindow = neuralData.eta.eventWindow;

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
whichNeuron = 1:length(neuralData(1).stats.bfcH);
cellResps = nanmean(neuralData.eta.alignedResps{2}(:,:,whichNeuron),3);

trialTypes = getTrialTypes(expInfo, behavioralData, 'all');
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

    subplot(2,length(contrasts),c+length(contrasts))
    plot(eventWindow, R_mRhL(c,:), 'Color',colors(c,:),'LineStyle','-','LineWidth',2);
    hold on
    plot(eventWindow, R_mRhR(c,:), 'Color',colors(c,:),'LineStyle',':','LineWidth',2);
    ylim([minY maxY]);
end

%% plot responses to contrast x movement

% whichNeuron = 29;
whichNeuron = neuralData(1).stats.bfcH(:,6) == 1 ;
% whichNeuron = 1:length(neuralData(1).stats.bfcH);
cellResps = nanmean(neuralData.eta.alignedResps{2}(:,:,whichNeuron),3);

trialTypes = getTrialTypes(expInfo, behavioralData, 'late');
whichTrials = trialTypes.intVar.all.contrast_direction;

for c = 1:size(whichTrials,1)
    mL = whichTrials{c,1};
    mR = whichTrials{c,2};
    
    R_mL(c,:) = nanmean(cellResps(mL,:),1);
    R_mR(c,:) = nanmean(cellResps(mR,:),1);
    vel_mL(c,:) = nanmean(behavioralData.wheelMoves.epochs(5).peakVel(mL));
    vel_mR(c,:) = nanmean(behavioralData.wheelMoves.epochs(5).peakVel(mR));
    
end   

figure;
set(gcf,'position',[120 560 1580 180])
for sp = 1:length(contrasts)
    subplot(1,length(contrasts),sp)
    hold on
    box off
    set(gca,'tickdir','out')
end

maxY = max(max([R_mL; R_mR]));
minY = min(min([R_mL; R_mR]));

for c = 1:size(whichTrials,1)
    subplot(1,length(contrasts),c)
    plot(eventWindow, R_mL(c,:), 'Color',colors(c,:),'LineStyle','-','LineWidth',2);
    hold on
    plot(eventWindow, R_mR(c,:), 'Color',colors(c,:),'LineStyle',':','LineWidth',2);
    ylim([minY maxY]);
    title(num2str(contrasts(c)*100),'Color',colors(c,:));
    if c == ceil(length(contrasts)/2)
        xlabel('time from event onset')
    end
end

    