function plotBlockScatter(plotAll, block, eventTimes, eventWindow, alignedTraces)
% this plots each cell's preference for a (correct) left or (correct) right
% turn (this is equivalent to left contrast or right contrast, save for the
% handful of correct 0% trials) vs preference for high-left or high-right
% block

contrasts = unique(block.events.contrastValues);

%event window
eventIdx = find(eventWindow == 0);

% condition comparison
title1 = 'high reward contra';
title2 = 'high reward ipsi';
legend1 = 'high contra';
legend2 = 'high ipsi';

for iCell = 1:length(plotAll)
    
    % select all the trials with LEFT turns under CONTRA high reward (LC)
    [~, condIdx_LC] = selectCondition(block, contrasts, eventTimes, ...
        'all', 'cw', 'all', 'left', 'correct', 'all', ...
        'all', 'all', 'all');

    trials_LC = alignedTraces{plotAll(iCell,1)}.eventSpikes(condIdx_LC,:,plotAll(iCell,2));

    % select all the trials with LEFT turns under IPSI high reward (LI)
    [~, condIdx_LI] = selectCondition(block, contrasts, eventTimes, ...
    'all', 'cw', 'all', 'right', 'correct', 'all', ...
    'all', 'all', 'all');

    trials_LI = alignedTraces{plotAll(iCell,1)}.eventSpikes(condIdx_LI,:,plotAll(iCell,2));
    
    % select all the trials with RIGHT turns under CONTRA high reward (RC)
    [~, condIdx_RC] = selectCondition(block, contrasts, eventTimes, ...
    'all', 'ccw', 'all', 'left', 'correct', 'all', ...
    'all', 'all', 'all');

    trials_RC = alignedTraces{plotAll(iCell,1)}.eventSpikes(condIdx_RC,:,plotAll(iCell,2));
    
    % select all the trials with RIGHT turns under IPSI high reward (RI)
    [~, condIdx_RI] = selectCondition(block, contrasts, eventTimes, ...
    'all', 'ccw', 'all', 'right', 'correct', 'all', ...
    'all', 'all', 'all');

    trials_RI = alignedTraces{plotAll(iCell,1)}.eventSpikes(condIdx_RI,:,plotAll(iCell,2));
    
    odd_turnLeft = [trials_LC(1:2:end,:); trials_LI(1:2:end,:)];
    odd_turnRight = [trials_RC(1:2:end,:); trials_RI(1:2:end,:)];
    odd_highContra = [trials_LC(1:2:end,:); trials_RC(1:2:end,:)];
    odd_highIpsi = [trials_LI(1:2:end,:); trials_RI(1:2:end,:)];
    
    even_turnLeft = [trials_LC(2:2:end,:); trials_LI(2:2:end,:)];
    even_turnRight = [trials_RC(2:2:end,:); trials_RI(2:2:end,:)];
    even_highContra = [trials_LC(2:2:end,:); trials_RC(2:2:end,:)];
    even_highIpsi = [trials_LI(2:2:end,:); trials_RI(2:2:end,:)];
    
    odd_mean_turnLeft = mean(odd_turnLeft,1) - min(mean(odd_turnLeft,1));
    odd_mean_turnRight = mean(odd_turnRight,1) - min(mean(odd_turnRight,1));
    odd_mean_highContra = mean(odd_highContra,1) - min(mean(odd_highContra,1));
    odd_mean_highIpsi = mean(odd_highIpsi,1) - min(mean(odd_highIpsi,1));
    
    even_mean_turnLeft = mean(odd_turnLeft,1) - min(mean(odd_turnLeft,1));
    even_mean_turnRight = mean(odd_turnRight,1) - min(mean(odd_turnRight,1));
    even_mean_highContra = mean(odd_highContra,1) - min(mean(odd_highContra,1));
    even_mean_highIpsi = mean(odd_highIpsi,1) - min(mean(odd_highIpsi,1));

    
    [allMax, allIdx] = max([odd_mean_turnLeft(eventIdx:end); odd_mean_turnRight(eventIdx:end); odd_mean_highIpsi(eventIdx:end); odd_mean_highContra(eventIdx:end)],[],2);
    [oneMax, oneIdx] = max(allMax);
    ymaxIdx = allIdx(oneIdx) + eventIdx-1;
    ym(iCell) = ymaxIdx;
    
    turnPref(iCell) = (even_mean_turnRight(ymaxIdx) - even_mean_turnLeft(ymaxIdx))./(even_mean_turnRight(ymaxIdx) + even_mean_turnLeft(ymaxIdx));
    blockPref(iCell) = (even_mean_highIpsi(ymaxIdx) - even_mean_highContra(ymaxIdx))./(even_mean_highIpsi(ymaxIdx) + even_mean_highContra(ymaxIdx));

end

if sum(blockPref) == 0
    blockPref = rand(1,length(plotAll))-.5;
end
    
epochIdx = eventIdx:3:length(eventWindow);
epochTimes = eventWindow(eventIdx:3:length(eventWindow))*1000;

cm = BlueWhiteRed(30);
cm = cm(eventIdx:end,:);
cm(1:eventIdx-1,:) = zeros(eventIdx-1,3).*.75;

figure;
hold on
lines(1) = line([-1 1],[0 0]);
lines(2) = line([0 0],[-1 1]);
uistack(lines,'bottom');
set(lines,'LineStyle', '--', 'LineWidth',1,'Marker','none','Color',[.5 .5 .5]);
ax1 = gca;
ax1.TickDir = 'out';
axis square


for iCell = 1:length(plotAll)
       pp = plot(turnPref(iCell),blockPref(iCell),'ko');
       set(pp,'LineStyle', 'none', 'LineWidth',.5,'MarkerEdgeColor','k','MarkerFaceColor',cm(ym(iCell),:),'MarkerSize',5);
end

% ylabel('high contra rewards                     high ipsi rewards')
% xlabel('left (cw) turns                             right (ccw) turns')
set(gcf,'renderer','Painters');



figure;
hold on
for e = 2:length(epochIdx)
    
subplot(2,ceil(length(epochIdx(1:end-1))/2),e-1);
hold on;
lines(1) = line([-1 1],[0 0]);
lines(2) = line([0 0],[-1 1]);
uistack(lines,'bottom');
set(lines,'LineStyle', '--', 'LineWidth',1,'Marker','none','Color',[.5 .5 .5]);
ax1 = gca;
ax1.TickDir = 'out';
axis square
title(strcat({'peak < '},num2str(epochTimes(e)),{' ms'}));


for iCell = 1:length(plotAll)
    if ym(iCell) <= epochIdx(e) && ym(iCell) > epochIdx(e-1)
       pp = plot(turnPref(iCell),blockPref(iCell),'ko');
       set(pp,'LineStyle', 'none', 'LineWidth',.5,'MarkerEdgeColor','k','MarkerFaceColor',cm(ym(iCell),:),'MarkerSize',5);
    end
end

% ylabel('high contra rewards                     high ipsi rewards')
% xlabel('left (cw) turns                             right (ccw) turns')

end
xticks([-1:1]);
xticklabels({'-1', '0','1'});
yticks([-1:1]);
yticklabels({'-1', '0','1'});

set(gcf,'renderer','Painters');

for r = 1:10
    subplot(2,5,r)
    xticks([-1:1]);
xticklabels({'-1', '0','1'});
yticks([-1:1]);
yticklabels({'-1', '0','1'});
end