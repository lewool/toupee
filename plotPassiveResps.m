%% overview

% this script plots responses to stimulus onsets in a contrast-dependent
% manner. Goes with lewieWorld_passive expDef but probably would work for
% any expDef with a stimOnset tag. Can omit trials based on movement
% artifacts (threshold is set in 'getEventTimes_passive.m') 

%2019-05-20: LEW created
%2020-01-22: updated for new expInfo struct, added to toupee

%% load experiment details

expInfo = initExpInfo({{'LEW025'}},{{'2020-01-22',2,2}});

%% load data

expInfo = data.loadExpData(expInfo);

%% get event timings and wheel trajectories
cd('C:\Users\Wool\Documents\MATLAB\expPipeline\expPipeline')
[eventTimes, wheelTrajectories] = getEventTimes_passive(expInfo, {'stimulusOnTimes' 'stimulusOffTimes'});
cd('C:\Users\Wool\Documents\GitHub\toupee');
%% times when animal was still around stim onset

isStill = [wheelTrajectories(:).goodQuiescenceFlag];
stillIdx = find(isStill > 0);
stimTimes = eventTimes(1).daqTime;
stimTimes_still = eventTimes(1).daqTime(stillIdx);

%% load traces

[allFcell, expInfo] = loadCellData(expInfo);
[cellResps, respTimes] = getCellResps(expInfo, allFcell);
cellResps = zscore(cellResps);

%% align calcium traces to the event you want

% cut the trace into trial-by-trial traces, aligned to a particular event
[alignedResps, eventWindow] = alignResps(expInfo, cellResps, respTimes, eventTimes, 'stimulusOnTimes');

%% select cells with the properties you want

plotCells = chooseCellType('vis_passive', expInfo, cellResps, respTimes, eventTimes, 0.1);

%% initialize some figure values
trialConditions = initTrialConditions;
whichScreen = [-90 0 90];
colors = cell(1,3);
colors{1} = {[0 0 .25],[0 0 .5],[0 0 1],[0 .4 1],[.6 .8 1],[.75 .75 .75]};
colors{2} = {[0 0 0],[.15 .15 .15 ],[.3 .3 .3],[.45 .45 .45],[.6 .6 .6],[.75 .75 .75]};
colors{3} = {[.25 0 0],[.5 0 0 ],[1 0 0],[.8 .45 .45],[.8 .7 .7],[.75 .75 .75]};

%initialize subplot titles
titles = {'left' 'center' 'right'};

m = 1;
max_m = length(plotCells);
contrasts = unique(expInfo.block.events.contrastValues);
%% cell by cell browser

fig = figure(201);
hold on

while m <= max_m
    %clear subplots
    for s = 1:3
        subplot(1,3,s);
        cla;
    end
    k = plotCells(m);
     
    %for each ETA matrix
    for sc = 1:length(whichScreen)
        screenIdx = find(expInfo.block.events.azimuthValues == whichScreen(sc));
        plotColors = fliplr(colors{sc});
        
        for c = 1:length(contrasts)
                % select trials based on high-reward side x contrast
                [~, condIdx] = selectCondition(expInfo.block, contrasts(c), eventTimes, trialConditions);
                condIdx = intersect(screenIdx, condIdx);
                
                %compute mean+sem responses
                meanResp(c,:,sc) = nanmean(alignedResps(condIdx,:,k),1);
                semResp(c,:,sc) = nanstd(alignedResps(condIdx,:,k)/sqrt(length(condIdx)));
                upperCI(c,:,sc) = meanResp(c,:,sc) + semResp(c,:,sc);
                lowerCI(c,:,sc) = meanResp(c,:,sc) - semResp(c,:,sc);
                
                % plot responses
                subplot(1,3,sc);
                hold on;
                plotResp = plot(eventWindow,meanResp(c,:,sc));
                plot1ci = fill([eventWindow';flipud(eventWindow')],[lowerCI(c,:,sc)';flipud(upperCI(c,:,sc)')],plotColors{c}, 'LineStyle', 'none');
                alpha(0.2);
                set(plotResp, 'LineStyle', '-', 'LineWidth',1.5,'Color',plotColors{c});
                eventOn = line([0 0],[-1 15]);
                uistack(eventOn,'bottom');
                set(eventOn,'LineStyle', '--', 'LineWidth',1,'Marker','none','Color',[.5 .5 .5]);
                
                title(titles{sc});
                
                if sc == 1
                    ylabel('z-scored spikes');
                end
                
                xlim([-.5 1.5]);
                box off
                xlabel('time (s)')
                hold on;
        end
    end
    
    for s = 1:3
        subplot(1,3,s);
        ylim([ min(reshape(lowerCI,[numel(lowerCI) 1]))*1.1 max(reshape(upperCI,[numel(upperCI) 1]))*1.1 ]);
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
                
                
                
                
                
