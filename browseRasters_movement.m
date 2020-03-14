function browseRasters_movement(cellLabel, response, expInfo, allFcell, eventTimes)
% 20 Nov 2019 LEW
% Plots four pseudoraster plots based on four trial conditions: left stim
% with left high reward, left stim with high right reward, right stim with
% high left reward, and right stim with high right reward
% cellLabel = the types of cells you want to look at; 'all', 'vis', 'movleft', or 'movright'
% response = filter trials by animal response; 'correct', 'incorrect', or 'all'
%% load data

block = expInfo.block;

%% align traces to both stimulus onset and movement onset (to be compared later)

%% align calcium traces to the event you want

% cut the trace into trial-by-trial traces, aligned to a particular event
events = {'stimulusOnTimes' 'prestimulusQuiescenceEndTimes' 'feedbackTimes'};
alignedResps = cell(1,length(events));
for e = 1:length(events)
    [alignedResps{e}, eventWindow] = alignResps(expInfo, cellResps, respTimes, eventTimes, events{e});
end

%% select cells with the properties you want
plotCells = chooseCellType('vis', expInfo, cellResps, respTimes, eventTimes, 0.1);

%% divide the trials up into 12 groups: stim (L, 0, R) x resp (L, R) x block (L, R) 
%(LLL, 0LL, RLL, LLR, 0LR, RLR, RRL, 0RL, RLL, RRR, 0RR, RLR)
contrasts = unique(expInfo.block.events.contrastValues);

%conditions for nonzero contrast
leftContrasts = contrasts(contrasts<0);
rightContrasts = contrasts(contrasts>0);
trialConditions_hiL_wentL = initTrialConditions('highRewardSide','left','movementDir','cw');
trialConditions_hiL_wentR = initTrialConditions('highRewardSide','left','movementDir','ccw');
trialConditions_hiR_wentL = initTrialConditions('highRewardSide','right','movementDir','cw');
trialConditions_hiR_wentR = initTrialConditions('highRewardSide','right','movementDir','ccw');

% conditions for zero contrast
zeroContrasts = 0;
trialConditions_hiL_wentL_correct = initTrialConditions('highRewardSide','left','responseType','correct','movementDir','cw');
trialConditions_hiL_wentR_correct = initTrialConditions('highRewardSide','left','responseType','correct','movementDir','ccw');
trialConditions_hiR_wentL_correct = initTrialConditions('highRewardSide','right','responseType','correct','movementDir','cw');
trialConditions_hiR_wentR_correct = initTrialConditions('highRewardSide','right','responseType','correct','movementDir','ccw');
trialConditions_hiL_wentL_incorrect = initTrialConditions('highRewardSide','left','responseType','incorrect','movementDir','cw');
trialConditions_hiL_wentR_incorrect = initTrialConditions('highRewardSide','left','responseType','incorrect','movementDir','ccw');
trialConditions_hiR_wentL_incorrect = initTrialConditions('highRewardSide','right','responseType','incorrect','movementDir','cw');
trialConditions_hiR_wentR_incorrect = initTrialConditions('highRewardSide','right','responseType','incorrect','movementDir','ccw');

%chose left during high left blocks
[~, condIdx.wentL.highL.stimL] = selectCondition(block, leftContrasts, eventTimes, trialConditions_hiL_wentL);
[~, condIdx.wentL.highL.stim0cor] = selectCondition(block, zeroContrasts, eventTimes, trialConditions_hiL_wentL_correct);
[~, condIdx.wentL.highL.stim0inc] = selectCondition(block, zeroContrasts, eventTimes, trialConditions_hiL_wentL_incorrect);
[~, condIdx.wentL.highL.stimR] = selectCondition(block, rightContrasts, eventTimes, trialConditions_hiL_wentL);

%chose left during high right blocks
[~, condIdx.wentL.highR.stimL] = selectCondition(block, leftContrasts, eventTimes, trialConditions_hiR_wentL);
[~, condIdx.wentL.highR.stim0cor] = selectCondition(block, zeroContrasts, eventTimes, trialConditions_hiR_wentL_correct);
[~, condIdx.wentL.highR.stim0inc] = selectCondition(block, zeroContrasts, eventTimes, trialConditions_hiR_wentL_incorrect);
[~, condIdx.wentL.highR.stimR] = selectCondition(block, rightContrasts, eventTimes, trialConditions_hiR_wentL);

%chose right during high left blocks
[~, condIdx.wentR.highL.stimR] = selectCondition(block, rightContrasts, eventTimes, trialConditions_hiL_wentR);
[~, condIdx.wentR.highL.stim0cor] = selectCondition(block, zeroContrasts, eventTimes, trialConditions_hiL_wentR_correct);
[~, condIdx.wentR.highL.stim0inc] = selectCondition(block, zeroContrasts, eventTimes, trialConditions_hiL_wentR_incorrect);
[~, condIdx.wentR.highL.stimL] = selectCondition(block, leftContrasts, eventTimes, trialConditions_hiL_wentR);

%chose right during high right blocks
[~, condIdx.wentR.highR.stimR] = selectCondition(block, rightContrasts, eventTimes, trialConditions_hiR_wentR);
[~, condIdx.wentR.highR.stim0cor] = selectCondition(block, zeroContrasts, eventTimes, trialConditions_hiR_wentR_correct);
[~, condIdx.wentR.highR.stim0inc] = selectCondition(block, zeroContrasts, eventTimes, trialConditions_hiR_wentR_incorrect);
[~, condIdx.wentR.highR.stimL] = selectCondition(block, leftContrasts, eventTimes, trialConditions_hiR_wentR);

%%

stimColors{1} = [0 .4 1; .4 .4 .4; .4 .4 .4; 1 0 0];
stimColors{2} = [0 .4 1; .4 .4 .4; .4 .4 .4; 1 0 0];
stimColors{3} = [1 0 0; .4 .4 .4; .4 .4 .4; 0 .4 1];
stimColors{4} = [1 0 0; .4 .4 .4; .4 .4 .4; 0 .4 1];

moveColors{1} = [0 0 1; 0 0 1; 0 0 1; 0 0 1];
moveColors{2} = [0 0 1; 0 0 1; 0 0 1; 0 0 1];
moveColors{3} = [.5 0 0; .5 0 0; .5 0 0; .5 0 0];
moveColors{4} = [.5 0 0; .5 0 0; .5 0 0; .5 0 0];

rewColors{1} = [1 1 0 0];
rewColors{2} = [1 1 0 0];
rewColors{3} = [1 1 0 0];
rewColors{4} = [1 1 0 0];

fig = figure;
hold on
set(fig,'Position',[50   300   2140   1000])

if ~exist('k') == 1
    k = 1;
end
max_k = length(plotCells);
while k <= max_k

    for s = 1:12
        subplot(3,4,s)
        cla;
    end

    for e = 1:3 

        event = events{e};
        % sort the trials based on stim relative to movement, keep this
        allSorts = cell(4,4);
        a=0;
        ww = fieldnames(condIdx);
        for w = 1:length(ww)
            hh = fieldnames(condIdx.(ww{w}));
            for h = 1:length(hh)
                ss = fieldnames(condIdx.(ww{w}).(hh{h}));
                a=a+1;
                for s = 1:length(ss)
                    relativeTimes{1} = eventTimes(7).daqTime(condIdx.(ww{w}).(hh{h}).(ss{s}))'-eventTimes(strcmp({eventTimes.event},event)).daqTime(condIdx.(ww{w}).(hh{h}).(ss{s}))';
                    relativeTimes{2} = eventTimes(1).daqTime(condIdx.(ww{w}).(hh{h}).(ss{s}))'-eventTimes(strcmp({eventTimes.event},event)).daqTime(condIdx.(ww{w}).(hh{h}).(ss{s}))';
                    relativeTimes{3} = eventTimes(1).daqTime(condIdx.(ww{w}).(hh{h}).(ss{s}))'-eventTimes(strcmp({eventTimes.event},event)).daqTime(condIdx.(ww{w}).(hh{h}).(ss{s}))';
                    [~,sortIdx] = sort(abs(relativeTimes{e}),'ascend');
                    allSorts{s,a} = condIdx.(ww{w}).(hh{h}).(ss{s})(sortIdx)';
                end
            end
        end

        for c = 1:4
            subplot(3,4,sub2ind([4 3],c,e));
            start = 0;
            for f = 1:4
                g = imagesc(eventWindow,start+1:(start+length(allSorts{f,c})),alignedResps{e}(allSorts{f,c},:,plotCells(k)));
                hold on;
                stimTimes = eventTimes(1).daqTime(allSorts{f,c})'-eventTimes(strcmp({eventTimes.event},event)).daqTime(allSorts{f,c})';
                movTimes = eventTimes(7).daqTime(allSorts{f,c})'-eventTimes(strcmp({eventTimes.event},event)).daqTime(allSorts{f,c})';
                rewTimes = eventTimes(5).daqTime(allSorts{f,c})'-eventTimes(strcmp({eventTimes.event},event)).daqTime(allSorts{f,c})';
                plot(stimTimes,start+1:(start+length(allSorts{f,c})),'.', 'MarkerEdgeColor',stimColors{c}(f,:));
                plot(movTimes,start+1:(start+length(allSorts{f,c})),'.', 'MarkerEdgeColor',moveColors{c}(f,:));
                pr = scatter(rewTimes,start+1:(start+length(allSorts{f,c})),'k.');
                pr.MarkerEdgeAlpha = rewColors{c}(f);
                colormap(flipud(gray));
                start = start+length(allSorts{f,c});
                ylim([1 length(vertcat(allSorts{:,c}))]);
                box off
                iMax(sub2ind([4 3],c,e)) = max(max(g.CData));
                iMin(sub2ind([4 3],c,e)) = min(min(g.CData));
                ax = gca;
                ax.TickDir = 'out';
                ax.YTick = [];
            end

        end
    end

%     for s = 1:12
%         subplot(3,4,s)
%         caxis([min(iMin) max(iMax)]);
%     end
    
    was_a_key = waitforbuttonpress;
    if was_a_key && strcmp(get(fig, 'CurrentKey'), 'leftarrow')
      k = max(1, k - 1);
    elseif was_a_key && strcmp(get(fig, 'CurrentKey'), 'rightarrow')
      k = min(max_k, k + 1);
    end
    
end
%% plot raster browser - equal size subplots, move-aligned
fig = figure;
hold on
set(fig,'Position',[290   100   740   210])
startIdx = find(eventWindow == -2);
endIdx = find(eventWindow == 2);

if ~exist('k') == 1
    k = 1;
end
max_k = length(plotCells);

% select the conditions
contrasts = unique(block.events.contrastValues);
if strcmp(response, 'correct') == 1
    trialConditions_hiL = initTrialConditions('highRewardSide','left','responseType','correct');
    trialConditions_hiR = initTrialConditions('highRewardSide','right','responseType','correct');
elseif strcmp(response, 'incorrect') == 1
    trialConditions_hiL = initTrialConditions('highRewardSide','left','responseType','incorrect');
    trialConditions_hiR = initTrialConditions('highRewardSide','right','responseType','incorrect');
elseif strcmp(response, 'all') == 1
    trialConditions_hiL = initTrialConditions('highRewardSide','left');
    trialConditions_hiR = initTrialConditions('highRewardSide','right');
end

[~, condIdx_highStimL_highL] = selectCondition(block, contrasts(contrasts<0), eventTimes, trialConditions_hiL);
[~, condIdx_highStimL_highR] = selectCondition(block, contrasts(contrasts<0), eventTimes, trialConditions_hiR);

[~, condIdx_highStimR_highL] = selectCondition(block, contrasts(contrasts>0), eventTimes, trialConditions_hiL);
[~, condIdx_highStimR_highR] = selectCondition(block, contrasts(contrasts>0), eventTimes, trialConditions_hiR);

highStimLHighL_axes = subplot(1,4,1);
highStimLHighR_axes = subplot(1,4,2);
highStimRHighL_axes = subplot(1,4,3);
highStimRHighR_axes = subplot(1,4,4);

while k <= max_k
    
    subplot(highStimLHighL_axes);cla;
    subplot(highStimLHighR_axes);cla;
    subplot(highStimRHighL_axes);cla;
    subplot(highStimRHighR_axes);cla; 

    subplot(highStimLHighL_axes);
    title(strcat('Left high (',num2str(length(condIdx_highStimL_highL)),')'));
    relativeMovTimes = zeros(length(condIdx_highStimL_highL),1)';
    relativeStimTimes = eventTimes(1).daqTime(condIdx_highStimL_highL)-eventTimes(7).daqTime(condIdx_highStimL_highL);
    [~,sortIdx] = sort(relativeStimTimes,'descend');
    f = imagesc(eventWindow(startIdx:endIdx), 1:length(condIdx_highStimL_highL),alignedResps{2}(condIdx_highStimL_highL(sortIdx),startIdx:endIdx,plotCells(k)));
    hold on;
    plot(relativeStimTimes(sortIdx),1:length(condIdx_highStimL_highL),'g.');
    etaLine = line([0 0],[1 length(condIdx_highStimL_highL)]);
    set(etaLine,'Color',[1 0 1]);
    colormap(flipud(gray));
    xlim([eventWindow(startIdx) eventWindow(endIdx)]);
    box off
    iMax(1) = max(max(f.CData));
    iMin(1) = min(min(f.CData));
    ax = gca;
    ax.TickDir = 'out';
    ax.YTick = [];

    subplot(highStimRHighL_axes);
    title(strcat('Right low (',num2str(length(condIdx_highStimR_highL)),')'));
    relativeMovTimes = zeros(length(condIdx_highStimR_highL),1)';
    relativeStimTimes = eventTimes(1).daqTime(condIdx_highStimR_highL)-eventTimes(7).daqTime(condIdx_highStimR_highL);
    [~,sortIdx] = sort(relativeStimTimes,'descend');
    f = imagesc(eventWindow(startIdx:endIdx), 1:length(condIdx_highStimR_highL),alignedResps{2}(condIdx_highStimR_highL(sortIdx),startIdx:endIdx,plotCells(k)));
    hold on;
    plot(relativeStimTimes(sortIdx),1:length(condIdx_highStimR_highL),'g.');
    etaLine = line([0 0],[1 length(condIdx_highStimR_highL)]);
    set(etaLine,'Color',[1 0 1]);
    colormap(flipud(gray));
    xlim([eventWindow(startIdx) eventWindow(endIdx)]);
    box off
    iMax(2) = max(max(f.CData));
    iMin(2) = min(min(f.CData));
    ax = gca;
    ax.TickDir = 'out';
    ax.YTick = [];

    subplot(highStimLHighR_axes);
    title(strcat('Left low (',num2str(length(condIdx_highStimL_highR)),')'));
    relativeMovTimes = zeros(length(condIdx_highStimL_highR),1)';
    relativeStimTimes = eventTimes(1).daqTime(condIdx_highStimL_highR)-eventTimes(7).daqTime(condIdx_highStimL_highR);
    [~,sortIdx] = sort(relativeStimTimes,'descend');
    f = imagesc(eventWindow(startIdx:endIdx), 1:length(condIdx_highStimL_highR),alignedResps{2}(condIdx_highStimL_highR(sortIdx),startIdx:endIdx,plotCells(k)));
    hold on;
    plot(relativeStimTimes(sortIdx),1:length(condIdx_highStimL_highR),'g.');
    etaLine = line([0 0],[1 length(condIdx_highStimL_highR)]);
    set(etaLine,'Color',[1 0 1]);
    colormap(flipud(gray));
    xlim([eventWindow(startIdx) eventWindow(endIdx)]);
    box off
    iMax(3) = max(max(f.CData));
    iMin(3) = min(min(f.CData));
    ax = gca;
    ax.TickDir = 'out';
    ax.YTick = [];
    
    subplot(highStimRHighR_axes);
    title(strcat('Right high (',num2str(length(condIdx_highStimR_highR)),')'));
    relativeMovTimes = zeros(length(condIdx_highStimR_highR),1)';
    relativeStimTimes = eventTimes(1).daqTime(condIdx_highStimR_highR)-eventTimes(7).daqTime(condIdx_highStimR_highR);
    [~,sortIdx] = sort(relativeStimTimes,'descend');
    f = imagesc(eventWindow(startIdx:endIdx), 1:length(condIdx_highStimR_highR),alignedResps{2}(condIdx_highStimR_highR(sortIdx),startIdx:endIdx,plotCells(k)));
    hold on;
    plot(relativeStimTimes(sortIdx),1:length(condIdx_highStimR_highR),'g.');
    etaLine = line([0 0],[1 length(condIdx_highStimR_highR)]);
    set(etaLine,'Color',[1 0 1]);
    colormap(flipud(gray));
    xlim([eventWindow(startIdx) eventWindow(endIdx)]);
    box off
    iMax(4) = max(max(f.CData));
    iMin(4) = min(min(f.CData));
    ax = gca;
    ax.TickDir = 'out';
    ax.YTick = [];
    
    subplot(highStimLHighL_axes);
    caxis([min(iMin) max(iMax)]);
    subplot(highStimLHighR_axes);
    caxis([min(iMin) max(iMax)]);
    subplot(highStimRHighL_axes);
    caxis([min(iMin) max(iMax)]);
    subplot(highStimRHighR_axes);
    caxis([min(iMin) max(iMax)]);
    
    was_a_key = waitforbuttonpress;
    if was_a_key && strcmp(get(fig, 'CurrentKey'), 'leftarrow')
      k = max(1, k - 1);
    elseif was_a_key && strcmp(get(fig, 'CurrentKey'), 'rightarrow')
      k = min(max_k, k + 1);
    end
end