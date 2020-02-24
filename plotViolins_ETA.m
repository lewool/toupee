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
% 2020-01-14 streamlining

%% load experiment details

expInfo = initExpInfo({{'LEW008'}},{{'2019-02-07',1,1}});

%% load data

expInfo = data.loadExpData(expInfo);

%% get event timings and wheel trajectories

[eventTimes, wheelTrajectories] = getEventTimes(expInfo, {'stimulusOnTimes' 'interactiveOnTimes' 'stimulusOffTimes'});

%% load traces

[allFcell, expInfo] = loadCellData(expInfo);
[cellResps, respTimes] = getCellResps(expInfo, allFcell);
cellResps = zscore(cellResps);

%% align calcium traces to the event you want

% cut the trace into trial-by-trial traces, aligned to a particular event
events = {'stimulusOnTimes' 'prestimulusQuiescenceEndTimes' 'feedbackTimes'};
alignedResps = cell(1,length(events));
for e = 1:length(events)
    [alignedResps{e}, eventWindow] = alignResps(expInfo, cellResps, respTimes, eventTimes, events{e});
end

%% select cells with the properties you want

plotCells = chooseCellType('movleft', expInfo, cellResps, respTimes, eventTimes, 0.1);

%% initialize some data values

contrasts = unique(expInfo.block.events.contrastValues);

%set up trial conditions for hi-L and hi-R blocks
trialConditions{1} = initTrialConditions('highRewardSide','left','responseType','correct','movementTime','late');
trialConditions{2} = initTrialConditions('highRewardSide','right','responseType','correct','movementTime','late');
trialConditions{3} = initTrialConditions('highRewardSide','left','responseType','incorrect','movementTime','late');
trialConditions{4} = initTrialConditions('highRewardSide','right','responseType','incorrect','movementTime','late');

% filter out the longest RTs
allRTs = eventTimes(7).daqTime - eventTimes(1).daqTime;
shortRTs = find(allRTs <= 50);

%initialize response arrays
s1 = length(trialConditions);
s2 = length(alignedResps);
s3 = length(contrasts);

meanResp = zeros(length(contrasts),length(eventWindow),s1,s2);
semResp = zeros(length(contrasts),length(eventWindow),s1,s2);
upperCI = zeros(length(contrasts),length(eventWindow),s1,s2);
lowerCI = zeros(length(contrasts),length(eventWindow),s1,s2);

%% initialize some plot values

% set color for each contrast
allColors = [0 0 .5; 0 0 1; 0 .4 1; .6 .8 1; .75 .75 .75; .8 .45 .45; 1 0 0; .5 0 0; .25 0 0];

Fs = 0.1;
zeroIdx = find(contrasts == 0);
walkup = length(contrasts) - zeroIdx;
walkback = zeroIdx - 1;
zeroGray = find(allColors(:,1) == .75);
colors = allColors(zeroGray-walkback:zeroGray + walkup,:);
bandwidth = 0.5;
width = 0.002;
timepoint = 0.3;
eventIdx = find(eventWindow == 0);
timeIdx = eventIdx + round(timepoint/Fs);

%initialize subplot titles
titles = {'stimOn' 'moveOn' 'rewardOn'};

%initialize figure shape & duration
m = 1;
max_m = length(plotCells);


%% cell-by-cell browser    

fig = figure(201);
hold on
set(fig,'Position',[114   400   1080   630])

while m <= max_m
    %clear subplots
    for s = 1:s1*s2
        subplot(s1,s2,s);
        cla;
    end
    k = plotCells(m);
     
    %for each ETA matrix
    for r = 1:length(svdResps)
        
        %for each trial reward condition
        for tcon = 1:length(trialConditions)

            for c = 1:length(contrasts)
                % select trials based on high-reward side x contrast
                [~, condIdx] = selectCondition(expInfo.block, contrasts(c), eventTimes, trialConditions{tcon});
                condIdx = intersect(shortRTs, condIdx);
                
                condY{c,tcon,r} = alignedResps{r}(condIdx,timeIdx,m);
                condX{c,tcon,r} = (rand(length(condIdx),1)-.5)*.02+contrasts(c);
                try
                    maxY(c,tcon,r) = max(condY{c,tcon,r});
                catch
                    maxY(c,tcon,r) = NaN;
                end
                try
                    minY(c,tcon,r) = min(condY{c,tcon,r});
                catch
                    minY(c,tcon,r) = NaN;
                end
                
                %compute mean+sem responses
                meanResp(c,tcon,r) = nanmean(condY{c,tcon,r});
                semResp(c,tcon,r) = nanstd(condY{c,tcon,r})/sqrt(length(condY{c,tcon,r}));
                upperCI(c,tcon,r) = meanResp(c,tcon,r) + semResp(c,tcon,r);
                lowerCI(c,tcon,r) = meanResp(c,tcon,r) - semResp(c,tcon,r);
    
                % plot responses
                figure(201);
                subplot(s1,s2,sub2ind([s2 s1],r,tcon));
                hold on;
                pResp = scatter(condX{c,tcon,r},condY{c,tcon,r});
                pResp.MarkerEdgeAlpha = .5;
                set(pResp, 'Marker', '.','MarkerEdgeColor',colors(c,:));
                
                % violin plot
                if ~isempty(condY{c,tcon,r})
                    [density,value] = ksdensity(condY{c,tcon,r},'Bandwidth', bandwidth);
    %                 density = density(value >= min(data) & value <= max(data));
    %                 value = value(value >= min(data) & value <= max(data));
    %                 value(1) = min(data);
    %                 value(end) = max(data);
                    w = width*length(condY{c,tcon,r});
                    pViolin = fill([density*w+contrasts(c) -density(end:-1:1)*w+contrasts(c)], ...
                     [value value(end:-1:1)], colors(c,:), 'LineStyle', '-','EdgeColor',colors(c,:),'EdgeAlpha',0.2,'FaceAlpha',0.2);
                    uistack(pViolin,'bottom');
                else
                end
                
                % plot mean
                pMean = plot(contrasts(c), meanResp(c,tcon,r),'ko');
                set(pMean, 'LineStyle','-','Color', colors(c,:),'MarkerSize',5);
              
                
                if tcon == 1
                    title(titles{r});
                end
                
                if r == 1
                    ylabel('U(:,1)');
                end
                
                xlim([-1.2 1.2]);
                box off
                xlabel('time (s)')
                hold on;
                                
            end
        
        end
        
    end
    
    for s = 1:s1*s2
        subplot(s1,s2,s);
        ylim([-1 max(max(max(maxY)))*1.1 ]);
%         ylim([0 .5]);
        ax3 = gca;
        ax3.TickDir = 'out';
    end
    
    %press arrow keys to browse through cells
    was_a_key = waitforbuttonpress;
    if was_a_key && strcmp(get(fig, 'CurrentKey'), 'leftarrow')
      m = max(1, m - 1);
    elseif was_a_key && strcmp(get(fig, 'CurrentKey'), 'rightarrow')
      m = min(max_m, m + 1);
    end
    
end