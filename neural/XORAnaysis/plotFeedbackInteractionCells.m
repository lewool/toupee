rasters = sum(cellfun(@length,trialTypes.intVar.all.side_direction_outcome),1);
totalRaster = sum(sum(cellfun(@length,trialTypes.intVar.all.side_direction_outcome),1),3);

reshaped = [...
    trialTypes.intVar.all.side_direction_outcome(:,1,1) trialTypes.intVar.all.side_direction_outcome(:,2,2); 
    trialTypes.intVar.all.side_direction_outcome(:,1,2) trialTypes.intVar.all.side_direction_outcome(:,2,1)];
totalRaster = sum(cellfun(@length,reshaped));
cellno = 1392;

yl = [-.05 1.2];
ca = [0 1.5];

figure;
hold on;
set(gcf,'position',[453 402 260 1100])
ETA = 3;
psth = 50;
pcolors = [0 .4 1; .75 .75 .75; 1 0 0; 0 .4 1; .75 .75 .75; 1 0 0];
rcolors = [0 .4 1; 0 0 0; 1 0 0; 0 .4 1; 0 0 0; 1 0 0];

for r = 1:size(reshaped,1)
    if r < 4
        ls = '-';
    else
        ls = ':';
    end
    if ~isempty(reshaped{r,1})
        subplot(max(totalRaster+psth),1,[1 psth])
        sem = std(neuralData.eta.alignedResps{ETA}(reshaped{r,1},:,cellno),[],1)/...
            sqrt(length(reshaped{r,1}));
        plotSignal(...
            neuralData.eta.eventWindow,...
            mean(neuralData.eta.alignedResps{ETA}(reshaped{r,1},:,cellno),1),...
            mean(neuralData.eta.alignedResps{ETA}(reshaped{r,1},:,cellno),1) + sem,...
            mean(neuralData.eta.alignedResps{ETA}(reshaped{r,1},:,cellno),1) - sem,...
            pcolors(r,:),ls)
        line([0 0],yl,'LineStyle',':','Color','k')
        box off
        prettyPlot(gca)
        ylim(yl)
        set(gca, 'XTickLabels', {''})
        title('Chose left')
        if r ~= 1
            rasStart = psth + sum(cellfun(@length,reshaped(1:r-1,1))) + 1;
        else
            rasStart = psth + 1;
        end
        rasEnd = psth + sum(cellfun(@length,reshaped(1:r,1)));
        ax = subplot(max(totalRaster+psth),1,[rasStart rasEnd]);
        imagesc(neuralData.eta.eventWindow, 1:length(reshaped{r,1}),neuralData.eta.alignedResps{ETA}(reshaped{r,1},:,cellno))
        caxis(ca);
        line([0 0],[1 length(reshaped{r,1})],'LineStyle',':','Color','k')
        ax.YColor = 'w';
        box off
        colormap(ax,rasterColor(rcolors(r,:),100))
        set(gca,'TickDir','out','TickLength', [.00 .00])
        if r ~= size(reshaped,1)
            set(gca, 'XTickLabels', {''})
        end
        if r == size(reshaped,1)
            xlabel('Time from feedback (s)')
        end
    end
end
printfig(gcf,char(strcat({'LEW031 2020-02-03 cell'},num2str(cellno),{' chose left'})));

figure;
hold on;
set(gcf,'position',[753 402 260 1100])
for r = 1:size(reshaped,1)
    if r > 3
        ls = '-';
    else
        ls = ':';
    end
    if ~isempty(reshaped{r,2})
        subplot(max(totalRaster+psth),1,[1 psth])
        sem = std(neuralData.eta.alignedResps{ETA}(reshaped{r,2},:,cellno),[],1)/...
            sqrt(length(reshaped{r,1}));
        plotSignal(...
            neuralData.eta.eventWindow,...
            mean(neuralData.eta.alignedResps{ETA}(reshaped{r,2},:,cellno),1),...
            mean(neuralData.eta.alignedResps{ETA}(reshaped{r,2},:,cellno),1) + sem,...
            mean(neuralData.eta.alignedResps{ETA}(reshaped{r,2},:,cellno),1) - sem,...
            pcolors(r,:),ls)
        line([0 0],yl,'LineStyle',':','Color','k')
        box off
        prettyPlot(gca)
        ylim(yl)
        set(gca, 'XTickLabels', {''})
        title('Chose right')
        if r ~= 1
            rasStart = psth + sum(cellfun(@length,reshaped(1:r-1,2))) + 1;
        else
            rasStart = psth + 1;
        end
        rasEnd = psth + sum(cellfun(@length,reshaped(1:r,2)));
        ax = subplot(max(totalRaster+psth),1,[rasStart rasEnd]);
        imagesc(neuralData.eta.eventWindow, 1:length(reshaped{r,2}),neuralData.eta.alignedResps{ETA}(reshaped{r,2},:,cellno))
        caxis(ca);
        line([0 0],[1 length(reshaped{r,2})],'LineStyle',':','Color','k')
        ax.YColor = 'w';
        box off
%         if r ~= size(reshaped,1)
%             axis off
%         end
        colormap(ax,rasterColor(rcolors(r,:),100))
        set(gca,'TickDir','out','TickLength', [.00 .00])
        if r ~= size(reshaped,1)
            set(gca, 'XTickLabels', {''})
        end
        if r == size(reshaped,1)
            xlabel('Time from feedback (s)')
        end
%         prettyPlot(gca)
    end
end
printfig(gcf,char(strcat({'LEW031 2020-02-03 cell'},num2str(cellno),{' chose right'})));

%% zero contrast only

rasters = sum(cellfun(@length,trialTypes.intVar.all.side_direction_outcome(2,:,:)),1);
totalRaster = sum(sum(cellfun(@length,trialTypes.intVar.all.side_direction_outcome(2,:,:)),1),3);

reshaped = [...
    trialTypes.intVar.all.side_direction_outcome(2,1,1); trialTypes.intVar.all.side_direction_outcome(2,1,2); 
    trialTypes.intVar.all.side_direction_outcome(2,2,1); trialTypes.intVar.all.side_direction_outcome(2,2,2)];
totalRaster = sum(cellfun(@length,reshaped));
ri = randi(length(neuralData.eta.alignedResps{ETA}));
cellno = 1329;

yl = [-.05 1.05];
ca = [0 1.2];

figure;
hold on;
set(gcf,'position',[453 402 260 600])
ETA = 3;
psth = 50;
pcolors = [.1 .7 1; .75 0 1; .1 .7 .1; .75 0 0];
rcolors = [.1 .7 1; .75 0 1; .1 .7 .1; .75 0 0];

for r = 1:size(reshaped,1)
    if r < 3
        ls = '-';
    else
        ls = ':';
    end
    if ~isempty(reshaped{r,1})
        subplot(max(totalRaster+psth),1,[1 psth])
        sem = std(neuralData.eta.alignedResps{ETA}(reshaped{r,1},:,cellno),[],1)/...
            sqrt(length(reshaped{r,1}));
        plotSignal(...
            neuralData.eta.eventWindow,...
            mean(neuralData.eta.alignedResps{ETA}(reshaped{r,1},:,cellno),1),...
            mean(neuralData.eta.alignedResps{ETA}(reshaped{r,1},:,cellno),1) + sem,...
            mean(neuralData.eta.alignedResps{ETA}(reshaped{r,1},:,cellno),1) - sem,...
            pcolors(r,:),ls)
        line([0 0],yl,'LineStyle',':','Color','k')
        box off
        prettyPlot(gca)
        ylim(yl)
        set(gca, 'XTickLabels', {''})
        title(num2str(cellno));
        if r ~= 1
            rasStart = psth + sum(cellfun(@length,reshaped(1:r-1,1))) + 1;
        else
            rasStart = psth + 1;
        end
        rasEnd = psth + sum(cellfun(@length,reshaped(1:r,1)));
        ax = subplot(max(totalRaster+psth),1,[rasStart rasEnd]);
        imagesc(neuralData.eta.eventWindow, 1:length(reshaped{r,1}),neuralData.eta.alignedResps{ETA}(reshaped{r,1},:,cellno))
        caxis(ca);
        line([0 0],[1 length(reshaped{r,1})],'LineStyle',':','Color','k')
        ax.YColor = 'w';
        box off
        colormap(ax,rasterColor(rcolors(r,:),100))
        set(gca,'TickDir','out','TickLength', [.00 .00])
        if r ~= size(reshaped,1)
            set(gca, 'XTickLabels', {''})
        end
        if r == size(reshaped,1)
            xlabel('Time from feedback (s)')
        end
    end
end
% printfig(gcf,char(strcat({'LEW031 2020-02-03 cell'},num2str(cellno))));

%%

contrasts = getUniqueContrasts(expInfo);

[~, ll] = selectCondition(expInfo, contrasts(contrasts<0), behavioralData, ...
    initTrialConditions('movementTime','late','movementDir','cw'));
[~, lr] = selectCondition(expInfo, contrasts(contrasts<0), behavioralData, ...
    initTrialConditions('movementTime','late','movementDir','ccw'));

[~, zlc] = selectCondition(expInfo, 0, behavioralData, ...
    initTrialConditions('movementTime','late','movementDir','cw','responseType','correct'));
[~, zli] = selectCondition(expInfo, 0, behavioralData, ...
    initTrialConditions('movementTime','late','movementDir','cw','responseType','incorrect'));
[~, zrc] = selectCondition(expInfo, 0, behavioralData, ...
    initTrialConditions('movementTime','late','movementDir','ccw','responseType','correct'));
[~, zri] = selectCondition(expInfo, 0, behavioralData, ...
    initTrialConditions('movementTime','late','movementDir','ccw','responseType','incorrect'));

[~, rl] = selectCondition(expInfo, contrasts(contrasts>0), behavioralData, ...
    initTrialConditions('movementTime','late','movementDir','cw'));
[~, rr] = selectCondition(expInfo, contrasts(contrasts>0), behavioralData, ...
    initTrialConditions('movementTime','late','movementDir','ccw'));
%%
colors = [0 .4 1;.7 .9 1; .9 .7 .7;1 0 0];

figure;
hold on

cellno = 115;

subplot(1,2,1)
hold on;
%stim was on the left, mouse detected it left
plot(neuralData.eta.eventWindow, ...
    mean(neuralData.eta.alignedResps{4}(ll,:,cellno),1),...
    'LineWidth',2,'Color',colors(1,:),'LineStyle','-')
plot(neuralData.eta.eventWindow, ...
    mean(neuralData.eta.alignedResps{4}(zlc,:,cellno),1),...
    'LineWidth',2,'Color',colors(2,:),'LineStyle','-')
%stim was on the right, mouse detected it
plot(neuralData.eta.eventWindow, ...
    mean(neuralData.eta.alignedResps{4}(rr,:,cellno),1),...
    'LineWidth',2,'Color',colors(4,:),'LineStyle','-')
plot(neuralData.eta.eventWindow, ...
    mean(neuralData.eta.alignedResps{4}(zrc,:,cellno),1),...
    'LineWidth',2,'Color',colors(3,:),'LineStyle','-')
prettyPlot(gca);

subplot(1,2,2)
hold on;
%stim was on the left, mouse didn't detect it
plot(neuralData.eta.eventWindow, ...
    mean(neuralData.eta.alignedResps{4}(lr,:,cellno),1),...
    'LineWidth',2,'Color',colors(1,:),'LineStyle',':')
plot(neuralData.eta.eventWindow, ...
    mean(neuralData.eta.alignedResps{4}(zri,:,cellno),1),...
    'LineWidth',2,'Color',colors(2,:),'LineStyle',':')
%stim was on the right, mouse didn't detect it
plot(neuralData.eta.eventWindow, ...
    mean(neuralData.eta.alignedResps{4}(rl,:,cellno),1),...
    'LineWidth',2,'Color',colors(4,:),'LineStyle',':')
plot(neuralData.eta.eventWindow, ...
    mean(neuralData.eta.alignedResps{4}(zli,:,cellno),1),...
    'LineWidth',2,'Color',colors(3,:),'LineStyle',':')
prettyPlot(gca);

%% psths only
colors = [0 .4 1;.75 .75 .75; .75 .75 .75;1 0 0];

figure;
hold on

cellno = 1957;

subplot(1,2,1)
hold on;
title('Chose "left"')
line([0 0],[0 1],'LineStyle','--','Color','k');
%MOVE LEFT
plot(neuralData.eta.eventWindow, ...
    mean(neuralData.eta.alignedResps{4}(ll,:,cellno),1),...
    'LineWidth',2,'Color',colors(1,:),'LineStyle','-')
plot(neuralData.eta.eventWindow, ...
    mean(neuralData.eta.alignedResps{4}(zlc,:,cellno),1),...
    'LineWidth',2,'Color',colors(2,:),'LineStyle','-')
plot(neuralData.eta.eventWindow, ...
    mean(neuralData.eta.alignedResps{4}(zli,:,cellno),1),...
    'LineWidth',2,'Color',colors(3,:),'LineStyle',':')
plot(neuralData.eta.eventWindow, ...
    mean(neuralData.eta.alignedResps{4}(rl,:,cellno),1),...
    'LineWidth',2,'Color',colors(4,:),'LineStyle',':')
prettyPlot(gca);
ax = gca;
mx(1) = max(ax.YLim);
legend('','left stim (rewarded)', 'zero stim (rewarded)', 'zero stim (unrewarded)', 'right stim (unrewarded)','Location','nw')
legend boxoff;
xlabel('Time from go cue (s)')

subplot(1,2,2)
hold on;
line([0 0],[0 1],'LineStyle','--','Color','k');
%MOVE RIGHT
title('Chose "right"')
plot(neuralData.eta.eventWindow, ...
    mean(neuralData.eta.alignedResps{4}(lr,:,cellno),1),...
    'LineWidth',2,'Color',colors(1,:),'LineStyle',':')
plot(neuralData.eta.eventWindow, ...
    mean(neuralData.eta.alignedResps{4}(zri,:,cellno),1),...
    'LineWidth',2,'Color',colors(2,:),'LineStyle',':')
plot(neuralData.eta.eventWindow, ...
    mean(neuralData.eta.alignedResps{4}(zrc,:,cellno),1),...
    'LineWidth',2,'Color',colors(3,:),'LineStyle','-')
plot(neuralData.eta.eventWindow, ...
    mean(neuralData.eta.alignedResps{4}(rr,:,cellno),1),...
    'LineWidth',2,'Color',colors(4,:),'LineStyle','-')
prettyPlot(gca);
ax = gca;
mx(2) = max(ax.YLim);
legend('','left stim (unrewarded)', 'zero stim (unrewarded)', 'zero stim (rewarded)', 'right stim (rewarded)','Location','nw')
legend boxoff;
xlabel('Time from go cue')

for s=1:2
    subplot(1,2,s)
    ylim([0 max(mx)])
end

set(gcf,'position',figpos)