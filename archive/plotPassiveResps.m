%% overview

% this script plots responses to stimulus onsets in a contrast-dependent
% manner. Goes with lewieWorld_passive expDef but probably would work for
% any expDef with a stimOnset tag. Can omit trials based on movement
% artifacts (threshold is set in 'getEventTimes_passive.m') 

%2019-05-20: LEW created
%2020-01-22: updated for new expInfo struct, added to toupee

%% load experiment details

exx = initExpInfo({{'LEW046'}},{{'2021-06-17',1,1}});

%% load data

exx = data.loadExpData(exx);

%% get event timings and wheel trajectories
cd('C:\Users\Wool\Documents\MATLAB\expPipeline\expPipeline')
[eventTimes, wheelTrajectories] = getEventTimes_passive(exx, {'stimulusOnTimes' 'stimulusOffTimes'});
cd('C:\Users\Wool\Documents\GitHub\toupee');
%% times when animal was still around stim onset

isStill = [wheelTrajectories(:).goodQuiescenceFlag];
stillIdx = find(isStill > 0);
stimTimes = eventTimes(1).daqTime;
stimTimes_still = eventTimes(1).daqTime(stillIdx);

%% load traces

[allFcell, exx] = loadCellData(exx);
[cellResps, respTimes] = getCellResps(exx, allFcell);
cellResps = zscore(cellResps);

%% align calcium traces to the event you want

% cut the trace into trial-by-trial traces, aligned to a particular event
[alignedResps, eventWindow] = alignResps(exx, cellResps, respTimes, eventTimes, 'stimulusOnTimes');

%% select cells with the properties you want

plotCells = chooseCellType('vis_passive', exx, cellResps, respTimes, eventTimes, 0.1);

%% initialize some figure values

exx = expInfo(1);
alignedResps = neuralData(1).eta.alignedResps{1};
eventTimes = behavioralData(1).eventTimes;
eventWindow = neuralData(1).eta.eventWindow;

trialConditions = initTrialConditions;
whichScreen = [-90 0 90];
colors = cell(1,3);
colors{1} = {[0 0 .25],[0 0 .5],[0 0 1],[0 .4 1],[.6 .8 1],[.75 .75 .75]};
colors{2} = {[0 0 0],[.15 .15 .15 ],[.3 .3 .3],[.45 .45 .45],[.6 .6 .6],[.75 .75 .75]};
colors{3} = {[.25 0 0],[.5 0 0 ],[1 0 0],[.8 .45 .45],[.8 .7 .7],[.75 .75 .75]};

%initialize subplot titles
titles = {'left' 'center' 'right'};

m = 1;
max_m = size(alignedResps,3);
contrasts = unique(exx.block.events.contrastValues);
%% cell by cell browser

fig = figure(201);
set(gcf,'position',[1760 1230 760 218])
hold on

m=306;
if ~exist('m') == 1
    m = 1;
end
while m <= max_m
    %clear subplots
    for s = 1:3
        subplot(1,3,s);
        cla;
    end
    
    k = m;
     
    %for each ETA matrix
    for sc = 1:length(whichScreen)
        screenIdx = find(exx.block.events.azimuthValues == whichScreen(sc));
        plotColors = fliplr(colors{sc});
        
        for c = 1:length(contrasts)
                % select trials based on high-reward side x contrast
                [~, condIdx] = selectCondition(exx, contrasts(c), behavioralData(1), trialConditions);
                condIdx = intersect(screenIdx, condIdx);
                
                %compute mean+sem responses
                meanResp(c,:,sc) = nanmean(alignedResps(condIdx,:,k),1);
                semResp(c,:,sc) = nanstd(alignedResps(condIdx,:,k)/sqrt(length(condIdx)));
                upperCI(c,:,sc) = meanResp(c,:,sc) + semResp(c,:,sc);
                lowerCI(c,:,sc) = meanResp(c,:,sc) - semResp(c,:,sc);
                
                % plot responses
                subplot(1,3,sc);
                hold on;
                plotResp = plot(eventWindow,smooth(meanResp(c,:,sc),3));
%                 plot1ci = fill([eventWindow';flipud(eventWindow')],[lowerCI(c,:,sc)';flipud(upperCI(c,:,sc)')],plotColors{c}, 'LineStyle', 'none');
                alpha(0.2);
                set(plotResp, 'LineStyle', '-', 'LineWidth',2,'Color',plotColors{c});
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
        ylim([-0.3 2])
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
                
                
                
                
                
