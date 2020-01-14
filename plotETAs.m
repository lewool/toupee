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

expInfo = initExpInfo({{'LEW008'}},{{'2019-01-29',1,[1]}});

%% load data

expInfo = data.loadExpData(expInfo);

%% get event timings and wheel trajectories

[eventTimes, wheelTrajectories] = getEventTimes(expInfo, {'stimulusOnTimes' 'interactiveOnTimes' 'stimulusOffTimes'});

 %% load traces
[allFcell, expInfo] = loadCellData(expInfo);

%% align calcium traces to the event you want

% cut the trace into trial-by-trial traces, aligned to a particular event
events = {'stimulusOnTimes' 'prestimulusQuiescenceEndTimes' 'rewardOnTimes'};
for e = 1:length(events)
    [alignedResps{e}, eventWindow] = alignResps(expInfo, cellResps, respTimes, eventTimes, events{e});
end

%% select cells with the properties you want

plotCells = chooseCellType('vis', expInfo, cellResps, respTimes, eventTimes, 0.1);

%% initialize some data values

%set up trial conditions for hi-L and hi-R blocks
trialConditions{1} = initTrialConditions('highRewardSide','left','responseType','all');
trialConditions{2} = initTrialConditions('highRewardSide','right','responseType','all');

% filter out the longest RTs
allRTs = eventTimes(7).daqTime - eventTimes(1).daqTime;
shortRTs = find(allRTs <= 2);

%initialize response arrays
meanResp = zeros(length(contrasts),length(eventWindow),s1,s2);
semResp = zeros(length(contrasts),length(eventWindow),s1,s2);
upperCI = zeros(length(contrasts),length(eventWindow),s1,s2);
lowerCI = zeros(length(contrasts),length(eventWindow),s1,s2);

%% initialize some plot values

% set color for each contrast
allColors = [0 0 .5; 0 0 1; 0 .4 1; .6 .8 1; .75 .75 .75; .8 .45 .45; 1 0 0; .5 0 0; .25 0 0];

contrasts = unique(expInfo.block.events.contrastValues);
zeroIdx = find(contrasts == 0);
walkup = length(contrasts) - zeroIdx;
walkback = zeroIdx - 1;
zeroGray = find(allColors(:,1) == .75);
colors = allColors(zeroGray-walkback:zeroGray + walkup,:);

%initialize subplot titles
titles = {'stimOn' 'moveOn' 'rewardOn'};

%initialize figure shape & duration
m = 1;
max_m = length(plotCells);
s1 = length(trialConditions);
s2 = length(alignedResps);

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
    for r = 1:length(alignedResps)
        
        %for each trial reward condition
        for tcon = 1:length(trialConditions)

            for c = 1:length(contrasts)
                % select trials based on high-reward side x contrast
                [~, condIdx] = selectCondition(expInfo.block, contrasts(c), eventTimes, trialConditions{tcon});
                condIdx = intersect(shortRTs, condIdx);
                
                %compute meam+sem responses
                meanResp(c,:,tcon,r) = nanmean(alignedResps{r}(condIdx,:,k),1);
                semResp(c,:,tcon,r) = nanstd(alignedResps{r}(condIdx,:,k)/sqrt(length(condIdx)));
                upperCI(c,:,tcon,r) = meanResp(c,:,tcon,r) + semResp(c,:,tcon,r);
                lowerCI(c,:,tcon,r) = meanResp(c,:,tcon,r) - semResp(c,:,tcon,r);
    
                % plot responses
                subplot(s1,s2,sub2ind([s2 s1],r,tcon));
                hold on;
                plotResp = plot(eventWindow,meanResp(c,:,tcon,r));
                plot1ci = fill([eventWindow';flipud(eventWindow')],[lowerCI(c,:,tcon,r)';flipud(upperCI(c,:,tcon,r)')],colors(c,:), 'LineStyle', 'none');
                alpha(0.2);
                set(plotResp, 'LineStyle', '-', 'LineWidth',1.5,'Color',colors(c,:));
                eventOn = line([0 0],[-1 15]);
                uistack(eventOn,'bottom');
                set(eventOn,'LineStyle', '--', 'LineWidth',1,'Marker','none','Color',[.5 .5 .5]);
                
                if tcon == 1
                    title(titles{r});
                end
                
                if r == 1
                    ylabel('z-scored spikes');
                end
                
                xlim([-.5 1.5]);
                box off
                xlabel('time (s)')
                hold on;
            end
        
        end
        
    end
    
    for s = 1:s1*s2
        subplot(s1,s2,s);
        ylim([ min(reshape(lowerCI,[prod(size(lowerCI)) 1]))*1.1 max(reshape(upperCI,[prod(size(upperCI)) 1]))*1.1 ]);
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


