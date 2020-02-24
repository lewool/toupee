function plotViolins(alignedResps, eventWindow, trialConditions, expInfo, eventTimes, plotCells, bandwidth, width)

%% overview



%% initialize some data values

contrasts = unique(expInfo.block.events.contrastValues);

% filter out the longest RTs
allRTs = eventTimes(7).daqTime - eventTimes(1).daqTime;
shortRTs = find(allRTs <= 5);

%initialize response arrays
s1 = length(trialConditions);
s2 = length(alignedResps);
s3 = length(contrasts);

meanResp = zeros(s3,s1,s2);
semResp = zeros(s3,s1,s2);
upperCI = zeros(s3,s1,s2);
lowerCI = zeros(s3,s1,s2);

%% initialize some plot values

% set color for each contrast
allColors = [0 0 .5; 0 0 1; 0 .4 1; .6 .8 1; .75 .75 .75; .8 .45 .45; 1 0 0; .5 0 0; .25 0 0];


zeroIdx = find(contrasts == 0);
walkup = length(contrasts) - zeroIdx;
walkback = zeroIdx - 1;
zeroGray = find(allColors(:,1) == .75);
colors = allColors(zeroGray-walkback:zeroGray + walkup,:);

if nargin < 7
    bandwidth = 0.01;
    width = 0.00005;
end

Fs = 0.1;
timepoint = 0.3;
eventIdx = find(eventWindow == 0);
timeIdx = eventIdx + round(timepoint/Fs);

%initialize figure shape & duration
m = 1;
max_m = length(plotCells);


%% cell-by-cell browser    

fig = figure(201);
hold on
set(fig,'Position',[124   300   1480   1030])
minY = nan(s3,s1,s2);
maxY = nan(s3,s1,s2);

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
                
                condY{c,tcon,r} = alignedResps{r}(condIdx,timeIdx,k);
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
                meanResp(c,tcon,r) = nanmedian(condY{c,tcon,r});
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
              
                
%                 if tcon == 1
%                     title(titles{r});
%                 end
                
                if r == 1
                    ylabel('z-scored spikes');
                end
                
                xlim([-max(contrasts)*1.1 max(contrasts)*1.1]);
                box off
                xlabel('contrast')
                hold on;
                                
            end
        
        end
        
    end
    
    for s = 1:s1*s2
        subplot(s1,s2,s);
        ylim([min(min(min(minY)))*3 max(max(max(maxY)))*1.1 ]);
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

% %% population response
% 
% fig = figure(203);
% hold on
% set(fig,'Position',[114   400   1080   630])
% 
% figure(203);
% for s = 1:s1*s2
%     subplot(s1,s2,s);
%     cla;
% end
% 
% figure(303);
% for s = 1:s3*s2
%     subplot(s2,s3,s);
%     cla;
% end
%     
% %for each ETA matrix
% for r = 1:length(svdResps)
% 
%     %for each trial reward condition
%     for tcon = 1:length(trialConditions)
% 
%         for c = 1:length(contrasts)
%             % select trials based on high-reward side x contrast
%             [~, condIdx] = selectCondition(expInfo.block, contrasts(c), eventTimes, trialConditions{tcon});
%             condIdx = intersect(shortRTs, condIdx);
% 
%             %compute mean+sem responses
%             meanResp(c,:,tcon,r) = nanmean(squeeze(nanmean(svdResps{r}(condIdx,:,plotCells),1))',1);
%             semResp(c,:,tcon,r) = nanstd(squeeze(nanmean(svdResps{r}(condIdx,:,plotCells),1))')/sqrt(length(plotCells));
%             upperCI(c,:,tcon,r) = meanResp(c,:,tcon,r) + semResp(c,:,tcon,r);
%             lowerCI(c,:,tcon,r) = meanResp(c,:,tcon,r) - semResp(c,:,tcon,r);
%             
%             % plot responses
%             figure(203);
%             subplot(s1,s2,sub2ind([s2 s1],r,tcon));
%             hold on;
%             plotResp = plot(eventWindow,meanResp(c,:,tcon,r));
%             plot1ci = fill([eventWindow';flipud(eventWindow')],[lowerCI(c,:,tcon,r)';flipud(upperCI(c,:,tcon,r)')],colors(c,:), 'LineStyle', 'none');
%             alpha(0.2);
%             set(plotResp, 'LineStyle', '-', 'LineWidth',1.5,'Color',colors(c,:));
%             eventOn = line([0 0],[-1 15]);
%             uistack(eventOn,'bottom');
%             set(eventOn,'LineStyle', '--', 'LineWidth',1,'Marker','none','Color',[.5 .5 .5]);
%             
%             if tcon == 1
%                 title(titles{r});
%             end
%                 
%             if r == 1
%                 ylabel('z-scored spikes');
%             end
% 
%             xlim([-.5 1.5]);
%             box off
%             xlabel('time (s)')
%             hold on;
%             
%             % plot responses
%             figure(303);
%             subplot(s2,s3,sub2ind([s3 s2],c,r));
%             hold on;
%             plotResp = plot(eventWindow,meanResp(c,:,tcon,r));
%             plot1ci = fill([eventWindow';flipud(eventWindow')],[lowerCI(c,:,tcon,r)';flipud(upperCI(c,:,tcon,r)')],colors(c,:), 'LineStyle', 'none');
%             alpha(0.2);
%             if tcon == 1
%                 set(plotResp, 'LineStyle', '-', 'LineWidth',1.5,'Color',colors(c,:));
%             else
%                 set(plotResp, 'LineStyle', '--', 'LineWidth',1.5,'Color',colors(c,:));
%             end
%             eventOn = line([0 0],[-1 15]);
%             uistack(eventOn,'bottom');
%             set(eventOn,'LineStyle', '--', 'LineWidth',1,'Marker','none','Color',[.5 .5 .5]);
% 
%             if tcon == 1
%                 title(titles{r});
%             end
% 
%             if c == 1
%                 ylabel('z-scored spikes');
%             end
% 
%             xlim([-.5 1.5]);
%             box off
%             if r == 3
%                 xlabel('time (s)')
%             end
%             hold on;
% 
%         end
%         
%     end
%         
% end
% 
% figure(203);    
% for s = 1:s1*s2
%     subplot(s1,s2,s);
%     ylim([ min(reshape(lowerCI,[numel(lowerCI) 1]))*1.1 max(reshape(upperCI,[numel(upperCI) 1]))*1.1 ]);
%     ax3 = gca;
%     ax3.TickDir = 'out';
% end
% 
% figure(303);    
% for s = 1:s3*s2
%     subplot(s2,s3,s);
%     ylim([ min(reshape(lowerCI,[numel(lowerCI) 1]))*1.1 max(reshape(upperCI,[numel(upperCI) 1]))*1.1 ]);
%     ax3 = gca;
%     ax3.TickDir = 'out';
% end