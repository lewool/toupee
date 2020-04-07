%% initialize experiment details
clear all;
% expInfo = initExpInfo({{'LEW032'}},{{'2020-02-13',1,[1]}});
expInfo = initExpInfo('LEW031');

%% load data & extract some variables to make this code still work - MATCHED

[expInfo, neuralData, behavioralData] = processExperiment(expInfo,'matched');
[neuralData] = alignResps(expInfo, neuralData, behavioralData);
[neuralData] = getSignificantActivity(expInfo, behavioralData, neuralData);
combinedNeuralData = combineNeuralData(expInfo, behavioralData, neuralData,'matched');

alignedResps = combinedNeuralData.matched.eta.alignedResps;
eventWindow = combinedNeuralData.matched.eta.eventWindow;
bfcH = combinedNeuralData.matched.stats.bfcH;
pLabels = combinedNeuralData.matched.stats.labels;
events = combinedNeuralData.matched.eta.events;

%% load data & extract some variables to make this code still work 

% [expInfo, neuralData, behavioralData] = processExperiment(expInfo);
% [neuralData] = alignResps(expInfo, neuralData, behavioralData);
% [neuralData] = getSignificantActivity(expInfo, behavioralData, neuralData);
% combinedNeuralData = combineNeuralData(expInfo, behavioralData, neuralData);
% 
% alignedResps = neuralData.eta.alignedResps;
% eventWindow = neuralData.eta.eventWindow;
% bfcH = neuralData.stats.bfcH;
% pLabels = neuralData.stats.labels;
% events = neuralData.eta.events;
% 
%% select cells with the properties you want
whichETA = 1;
Fs = 0.1;

%% get cel responses at a particular timepoint

%%%%%%%%%%%%%%%% compute baseline activity

% align traces to stim onset
event = 'stimulusOnTimes';
stim_alignedTraces = alignedResps{strcmp(events,event)};
stim_eventWindow = eventWindow;

%designate a baseline window
stim_eventIdx = find(stim_eventWindow == 0);
stim_preTime = [-0.5 0] / Fs;
baselineIdx = stim_eventIdx + stim_preTime(1) : stim_eventIdx;

%compute the mean baseline activity per cell, per trial (trials x neurons)
baselineResps = squeeze(mean(stim_alignedTraces(:,baselineIdx,:),2));

%%%% compute peristimulus activity

%designate a peristimulus window
stimTime = [0 0.3] / Fs;
stimIdx = stim_eventIdx + stimTime(1) :stim_eventIdx + stimTime(2);

%compute the mean peristimulus activity per cell, per trial (trials x neurons)
stimResps = squeeze(mean(stim_alignedTraces(:,stimIdx,:),2));

%%%%%%%%%%%%%%%% compute perimovement activity

% align traces to movement onset
event = 'prestimulusQuiescenceEndTimes';
mov_alignedTraces = alignedResps{strcmp(events,event)};
mov_eventWindow = eventWindow;

%designate a movement window
mov_eventIdx = find(mov_eventWindow == 0);
movTime = [-0.2 0.1] / Fs;
movIdx = mov_eventIdx + movTime(1) : mov_eventIdx + movTime(2);

%compute the mean perimovement activity per cell, per trial (trials x neurons)
movResps = squeeze(mean(mov_alignedTraces(:,movIdx,:),2));

%%%%%%%%%%%%%%%% compute premovement activity

%designate a movement window
pmov_eventIdx = find(mov_eventWindow == 0);
pmovTime = [-0.7 -0.1] / Fs;
pmovIdx = pmov_eventIdx + pmovTime(1) : pmov_eventIdx + pmovTime(2);

%compute the mean perimovement activity per cell, per trial (trials x neurons)
pmovResps = squeeze(mean(mov_alignedTraces(:,pmovIdx,:),2));

%%%%%%%%%%%%%%%% compute perireward activity

% align traces to movement onset
event = 'feedbackTimes';
rew_alignedTraces = alignedResps{strcmp(events,event)};
rew_eventWindow = eventWindow;

%designate a movement window
rew_eventIdx = find(rew_eventWindow == 0);
rewTime = [0 0.2] / Fs;
rewIdx = rew_eventIdx + rewTime(1) : rew_eventIdx + rewTime(2);

%compute the mean perireward activity per cell, per trial (trials x neurons)
rewResps = squeeze(mean(rew_alignedTraces(:,rewIdx,:),2));

%% set up trial conditions to compare

clear contrastConditions trialConditions labels condIdx
contrasts = getUniqueContrasts(expInfo);

%set up trial conditions for hi-L and hi-R blocks
trialConditions{1} = initTrialConditions('highRewardSide','left','movementDir','cw','movementTime','late');
trialConditions{2} = initTrialConditions('highRewardSide','left','movementDir','ccw','movementTime','late');
trialConditions{3} = initTrialConditions('highRewardSide','right','movementDir','cw','movementTime','late');
trialConditions{4} = initTrialConditions('highRewardSide','right','movementDir','ccw','movementTime','late');
trialConditions{5} = initTrialConditions('highRewardSide','left','movementTime','late');
trialConditions{6} = initTrialConditions('highRewardSide','right','movementTime','late');
trialConditions{7} = initTrialConditions('movementDir','cw','movementTime','late');
trialConditions{8} = initTrialConditions('movementDir','ccw','movementTime','late');

trialLabels{1} = 'bL_mL_';
trialLabels{2} = 'bL_mR_';
trialLabels{3} = 'bR_mL_';
trialLabels{4} = 'bR_mR_';
trialLabels{5} = 'bL_mAll_';
trialLabels{6} = 'bR_mAll_';
trialLabels{7} = 'bAll_mL_';
trialLabels{8} = 'bAll_mR_';

contrastConditions{1} = contrasts(contrasts<0);
contrastConditions{2} = contrasts(contrasts>0);
contrastConditions{3} = contrasts(contrasts~=0);
contrastLabels{1} = 'sL';
contrastLabels{2} = 'sR';
contrastLabels{3} = 'sAll';

testTrials = 1:2:size(alignedResps{whichETA},1);
trainTrials = 2:2:size(alignedResps{whichETA},1);

d = 1;
for c = 1:length(contrastConditions)
    for t = 1:length(trialConditions)
        [~, condIdx{d,:}.all] = selectCondition(expInfo, contrastConditions{c}, behavioralData, trialConditions{t});
        condIdx{d,:}.test = intersect(testTrials,condIdx{d}.all);
        condIdx{d,:}.train = intersect(trainTrials,condIdx{d}.all);
        labels{d,1} = strcat(trialLabels{t},contrastLabels{c});
        d = d+1;
    end
end

%%
close all
figure;
colors = [0.1 0.7 0.1; 1 .6 0; 0 .4 1; 1 0 0];

set(gcf,'position',[40 90 1230 900]);
respPeriod = stimResps;
plotCells = find(bfcH(:,strcmp(pLabels,'leftMov')) > 0);

% condVector = {'bR_mL_sL' 'bR_mR_sL' 'bR_mL_sR' 'bR_mR_sR'};
% color = colors(2,:);
% condVector = {'bL_mL_sL' 'bL_mR_sL' 'bL_mL_sR' 'bL_mR_sR'};
% color = colors(1,:);
condVector = {'bAll_mL_sL' 'bAll_mR_sL' 'bAll_mL_sR' 'bAll_mR_sR'};
color = [0 0 0];
mSize = 12;

X = [...
    mean(respPeriod(condIdx{strcmp(labels,condVector{1})}.all,plotCells),1);
    mean(respPeriod(condIdx{strcmp(labels,condVector{2})}.all,plotCells),1);
    mean(respPeriod(condIdx{strcmp(labels,condVector{3})}.all,plotCells),1);
    mean(respPeriod(condIdx{strcmp(labels,condVector{4})}.all,plotCells),1);
    ]';

axLim = 1.05 * [min(min(X)) max(max(X)) min(min(X)) max(max(X))];
for r = 1:4
    x1 = X(:,r);
    for c = 1:4
        x2 = X(:,c);
        subplot(4,4,sub2ind([4 4],r,c));
        hold on;
        if r == c
%             cla;
            h = histogram(x1,linspace(axLim(1),axLim(2),11),'Normalization','probability');
            set(h,'FaceColor',color,...
                'FaceAlpha',.2...
                );
            set(gca,'Color','none');
            box off;
            axis square;
            xlim([axLim(1) axLim(2)]);
        else
            p = scatter(x1,x2,mSize);
            l = line([axLim(1) axLim(2)],[axLim(1) axLim(2)]);
            set(p,...
                'MarkerEdgeColor','none',...
                'Marker','o',...
                'MarkerFaceColor',color,...
                'MarkerFaceAlpha',.2...
                );
            set(l,...
                'LineStyle', '--',...
                'LineWidth',1,...
                'Color',[.5 .5 .5]...
                );
            axis(axLim);
        end
        axis square
        
        if c == 4
            xlabel(condVector{r},'Interpreter','none')
        end
        
        if r == 1
            ylabel(condVector{c},'Interpreter','none')
        end
        
    end
end


condVectorY = {'bR_mL_sL' 'bR_mR_sL' 'bR_mL_sR' 'bR_mR_sR'};
condVectorX = {'bL_mL_sL' 'bL_mR_sL' 'bL_mL_sR' 'bL_mR_sR'};

mSize = 12;
color = [0 0 0];
X = [...
    mean(respPeriod(condIdx{strcmp(labels,condVectorX{1})}.all,plotCells),1);
    mean(respPeriod(condIdx{strcmp(labels,condVectorX{2})}.all,plotCells),1);
    mean(respPeriod(condIdx{strcmp(labels,condVectorX{3})}.all,plotCells),1);
    mean(respPeriod(condIdx{strcmp(labels,condVectorX{4})}.all,plotCells),1);
    ]';

Y = [...
    mean(respPeriod(condIdx{strcmp(labels,condVectorY{1})}.all,plotCells),1);
    mean(respPeriod(condIdx{strcmp(labels,condVectorY{2})}.all,plotCells),1);
    mean(respPeriod(condIdx{strcmp(labels,condVectorY{3})}.all,plotCells),1);
    mean(respPeriod(condIdx{strcmp(labels,condVectorY{4})}.all,plotCells),1);
    ]';

axLim = 1.05 * [min(min([X Y])) max(max([X Y])) min(min([X Y])) max(max([X Y]))];
figure;
set(gcf,'position',[1280 90 260 900]);

for s = 1:4
    subplot(4,1,s)
    x = X(:,s);
    y = Y(:,s);
    p = scatter(x,y,mSize);
        l = line([axLim(1) axLim(2)],[axLim(1) axLim(2)]);
    title(condVectorX{s}(end-4:end),'Interpreter','none');  
    xlabel('bL')
    ylabel('bR')
    set(p,...
        'MarkerEdgeColor','none',...
        'Marker','o',...
        'MarkerFaceColor',color,...
        'MarkerFaceAlpha',.2...
        );
    set(l,...
        'LineStyle', '--',...
        'LineWidth',1,...
        'Color',[.5 .5 .5]...
        );
    axis(axLim);
    axis square
end

%% set up trial conditions to compare - separate contrasts

clear contrastConditions trialConditions labels condIdx
contrasts = getUniqueContrasts(expInfo);

%set up trial conditions for hi-L and hi-R blocks
trialConditions{1} = initTrialConditions('highRewardSide','left','movementDir','cw','movementTime','late');
trialConditions{2} = initTrialConditions('highRewardSide','left','movementDir','ccw','movementTime','late');
trialConditions{3} = initTrialConditions('highRewardSide','right','movementDir','cw','movementTime','late');
trialConditions{4} = initTrialConditions('highRewardSide','right','movementDir','ccw','movementTime','late');
trialConditions{5} = initTrialConditions('highRewardSide','left','movementTime','late');
trialConditions{6} = initTrialConditions('highRewardSide','right','movementTime','late');
trialConditions{7} = initTrialConditions('movementDir','cw','movementTime','late');
trialConditions{8} = initTrialConditions('movementDir','ccw','movementTime','late');

trialLabels{1} = 'bL_mL_';
trialLabels{2} = 'bL_mR_';
trialLabels{3} = 'bR_mL_';
trialLabels{4} = 'bR_mR_';
trialLabels{5} = 'bL_mAll_';
trialLabels{6} = 'bR_mAll_';
trialLabels{7} = 'bAll_mL_';
trialLabels{8} = 'bAll_mR_';

for o = 1:length(contrasts)
    contrastConditions{o} = contrasts(o);
    if contrasts(o) < 0
        contrastLabels{o} = num2str(contrasts(o));
    elseif contrasts(o) > 0
        contrastLabels{o} = strcat('+',num2str(contrasts(o)));
    end
end

% contrastConditions{1} = contrasts(contrasts<0);
% contrastConditions{2} = contrasts(contrasts>0);
% contrastConditions{3} = contrasts(contrasts~=0);
% contrastLabels{1} = 'sL';
% contrastLabels{2} = 'sR';
% contrastLabels{3} = 'sAll';

testTrials = 1:2:size(alignedResps{whichETA},1);
trainTrials = 2:2:size(alignedResps{whichETA},1);

d = 1;
for c = 1:length(contrastConditions)
    for t = 1:length(trialConditions)
        [~, condIdx{d,:}.all] = selectCondition(expInfo, contrastConditions{c}, behavioralData, trialConditions{t});
        condIdx{d,:}.test = intersect(testTrials,condIdx{d}.all);
        condIdx{d,:}.train = intersect(trainTrials,condIdx{d}.all);
        labels{d,1} = strcat(trialLabels{t},contrastLabels{c});
        d = d+1;
    end
end

%%
close all
figure;
colors = [0.1 0.7 0.1; 1 .6 0; 0 .4 1; 1 0 0];

set(gcf,'position',[40 90 1230 900]);
respPeriod = stimResps;
plotCells = find(bfcH(:,strcmp(pLabels,'stim')) > 0);

% condVector = {'bR_mL_sL' 'bR_mR_sL' 'bR_mL_sR' 'bR_mR_sR'};
% color = colors(2,:);
% condVector = {'bL_mL_sL' 'bL_mR_sL' 'bL_mL_sR' 'bL_mR_sR'};
% color = colors(1,:);
condVector = {'bAll_mL_sL' 'bAll_mR_sL' 'bAll_mL_sR' 'bAll_mR_sR'};
color = [0 0 0];
mSize = 12;

for x = 1:length(condIdx)
    len(x,1) = length(condIdx{x}.all);
end

minCond = min(len(contains(labels,'bAll')));
minCond = 3;
bAll_mL_sL_idx = find(contains(labels,'-').*contains(labels,'bAll_mL'));
bAll_mL_sR_idx = find(contains(labels,'+').*contains(labels,'bAll_mL'));
bAll_mR_sL_idx = find(contains(labels,'-').*contains(labels,'bAll_mR'));
bAll_mR_sR_idx = find(contains(labels,'+').*contains(labels,'bAll_mR'));

bAll_mL_sL = []; bAll_mL_sR = []; bAll_mR_sL = []; bAll_mR_sR = [];
for i = 1:4
    bAll_mL_sL = [bAll_mL_sL, randsample(condIdx{bAll_mL_sL_idx(i)}.all,minCond)];
    bAll_mL_sR = [bAll_mL_sR, randsample(condIdx{bAll_mL_sR_idx(i)}.all,minCond)];
    bAll_mR_sL = [bAll_mR_sL, randsample(condIdx{bAll_mR_sL_idx(i)}.all,minCond)];
    bAll_mR_sR = [bAll_mR_sR, randsample(condIdx{bAll_mR_sR_idx(i)}.all,minCond)];
end

X = [...
    mean(respPeriod(bAll_mL_sL,plotCells),1);
    mean(respPeriod(bAll_mL_sR,plotCells),1);
    mean(respPeriod(bAll_mR_sL,plotCells),1);
    mean(respPeriod(bAll_mR_sR,plotCells),1);
    ]';

axLim = 1.05 * [min(min(X)) max(max(X)) min(min(X)) max(max(X))];
for r = 1:4
    x1 = X(:,r);
    for c = 1:4
        x2 = X(:,c);
        subplot(4,4,sub2ind([4 4],r,c));
        hold on;
        if r == c
%             cla;
            h = histogram(x1,linspace(axLim(1),axLim(2),11),'Normalization','probability');
            set(h,'FaceColor',color,...
                'FaceAlpha',.2...
                );
            set(gca,'Color','none');
            box off;
            axis square;
            xlim([axLim(1) axLim(2)]);
        else
            p = scatter(x1,x2,mSize);
            l = line([axLim(1) axLim(2)],[axLim(1) axLim(2)]);
            set(p,...
                'MarkerEdgeColor','none',...
                'Marker','o',...
                'MarkerFaceColor',color,...
                'MarkerFaceAlpha',.2...
                );
            set(l,...
                'LineStyle', '--',...
                'LineWidth',1,...
                'Color',[.5 .5 .5]...
                );
            axis(axLim);
        end
        axis square
        
        if c == 4
            xlabel(condVector{r},'Interpreter','none')
        end
        
        if r == 1
            ylabel(condVector{c},'Interpreter','none')
        end
        
    end
end

%%
condVectorY = {'bR_mL_sL' 'bR_mR_sL' 'bR_mL_sR' 'bR_mR_sR'};
condVectorX = {'bL_mL_sL' 'bL_mR_sL' 'bL_mL_sR' 'bL_mR_sR'};

mSize = 12;
color = [0 0 0];
X = [...
    mean(respPeriod(condIdx{strcmp(labels,condVectorX{1})}.all,plotCells),1);
    mean(respPeriod(condIdx{strcmp(labels,condVectorX{2})}.all,plotCells),1);
    mean(respPeriod(condIdx{strcmp(labels,condVectorX{3})}.all,plotCells),1);
    mean(respPeriod(condIdx{strcmp(labels,condVectorX{4})}.all,plotCells),1);
    ]';

Y = [...
    mean(respPeriod(condIdx{strcmp(labels,condVectorY{1})}.all,plotCells),1);
    mean(respPeriod(condIdx{strcmp(labels,condVectorY{2})}.all,plotCells),1);
    mean(respPeriod(condIdx{strcmp(labels,condVectorY{3})}.all,plotCells),1);
    mean(respPeriod(condIdx{strcmp(labels,condVectorY{4})}.all,plotCells),1);
    ]';

axLim = 1.05 * [min(min([X Y])) max(max([X Y])) min(min([X Y])) max(max([X Y]))];
figure;
set(gcf,'position',[1280 90 260 900]);

for s = 1:4
    subplot(4,1,s)
    x = X(:,s);
    y = Y(:,s);
    p = scatter(x,y,mSize);
        l = line([axLim(1) axLim(2)],[axLim(1) axLim(2)]);
    title(condVectorX{s}(end-4:end),'Interpreter','none');  
    xlabel('bL')
    ylabel('bR')
    set(p,...
        'MarkerEdgeColor','none',...
        'Marker','o',...
        'MarkerFaceColor',color,...
        'MarkerFaceAlpha',.2...
        );
    set(l,...
        'LineStyle', '--',...
        'LineWidth',1,...
        'Color',[.5 .5 .5]...
        );
    axis(axLim);
    axis square
end
