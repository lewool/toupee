function figName = rasterBrowser_earlyLate(expInfo, behavioralData, neuralData, whichCells, k)

% pickTrials = {'side_direction', 'side_direction', 'outcome_direction'};
% contrastOverride = 'contrast_direction';
% trialStruct = trialTypes.intVar.cb2D;
%% initialize experiment details

alignedResps = neuralData.eta.alignedResps;
eventWindow = neuralData.eta.eventWindow;
bfcH = neuralData.stats.bfcH;
pLabels = neuralData.stats.labels;
et = behavioralData.eventTimes;
wm = behavioralData.wheelMoves;


%% choose cells

% whichCells = 'leftStim'; %choose from 'pLabels' array
if strcmp(whichCells, 'all')
    plotCells = 1:size(alignedResps{1},3);
else
    plotCells = find(bfcH(:,strcmp(pLabels,whichCells)) > 0);
end

%%
contrasts = getUniqueContrasts(expInfo);

for iA = 1:3
    for c = 1:length(contrasts)
    [~, trialLists{iA}{c,1}] = selectCondition(expInfo, contrasts(c), behavioralData, ...
        initTrialConditions('preStimMovement','all','movementDir','cw','movementTime','late','specificRTs',[0 3]));
    [~, trialLists{iA}{c,2}] = selectCondition(expInfo, contrasts(c), behavioralData, ...
        initTrialConditions('preStimMovement','all','movementDir','ccw','movementTime','late','specificRTs',[0 3]));
    end
end

blockID = expInfo.block.events.highRewardSideValues;

%% set up the plots
zeroIdx = find(contrasts == 0);
walkup = length(contrasts) - zeroIdx;
walkback = zeroIdx - 1;
% allColors = [.25 0 0;.5 0 0 ;1 0 0;.8 .45 .45;.75 .75 .75;.55 .55 .55;.35 .35 .35;.15 .15 .15;0 0 0];
allColors = [0 0 .5; 0 0 1; 0 .4 1; .6 .8 1; .75 .75 .75; .8 .45 .45; 1 0 0; .5 0 0; .25 0 0];

zeroGray = find(allColors(:,1) == .75);
colors = allColors(zeroGray-walkback:zeroGray + walkup,:);
cgain = .5;

fig = figure;
set(fig, 'Position', [80 380 1080 660]);

rowLength = 4;
colLength = length(contrasts);
for c = 1:length(contrasts)
    for iA = 1:2
        
        %'chose left' columns
        pIdx1 = sub2ind([rowLength colLength],iA*2-1,c);
        subplot(colLength,rowLength,pIdx1)
        if c == colLength
            set(gca,'ytick',[])
        else
            set(gca,'xtick',[],'ytick',[])
        end
        
        %'chose right' columns
        pIdx2 = sub2ind([rowLength colLength],iA*2,c);
        subplot(colLength,rowLength,pIdx2)
        if c == colLength
            set(gca,'ytick',[])
        else
            set(gca,'xtick',[],'ytick',[])
        end
        
%         %'summary' columns
%         pIdx2 = sub2ind([rowLength colLength],5,c);
%         subplot(colLength,rowLength,pIdx2)
%         if c == colLength
%             set(gca,'ytick',[])
%         else
%             set(gca,'xtick',[],'ytick',[])
%         end
    end
end

ha = get(gcf,'children');
sp1 = get(ha(1),'position');
col4pos = sp1(1)-.03;
sp3 = get(ha(3),'position');
col2pos = sp3(1)-.03;
%
%
% fig = figure;
% set(fig, 'Position', [80 250 870 660]);
% hold on;

if ~exist('k') == 1
    k = 1;
end

max_k = length(plotCells);

while k <= max_k

iMin = [];
iMax = [];

% rowLength = 5;
% colLength = length(contrasts);

for a = 1:2
    for c = 1:length(contrasts)
        % extract the relevant trials for that condition
        choseLeftTrials = trialLists{a}{c,1};
        choseRightTrials = trialLists{a}{c,2};
        numLeftTrials = size(choseLeftTrials,2);
        numRightTrials = size(choseRightTrials,2);
        
        pIdx1 = sub2ind([rowLength colLength],a*2-1,c);
        ax1 = subplot(colLength,rowLength,pIdx1);
        f1 = imagesc(eventWindow,1:numLeftTrials,alignedResps{a}(choseLeftTrials,:,plotCells(k)));
        colormap(ax1,rasterColor([0 0 0]));
        hold on;
        line([0 0],[0 100],'LineStyle','-','Color',[0 0 0],'linewidth',1);
        for b = 1:numLeftTrials
            if blockID(choseLeftTrials(b)) < 0
                if a > 1
                    plot(1,b,'.','Color',[.1 .7 .1]);
                else
                    plot(1.5,b,'.','Color',[.1 .7 .1]);
                end
            elseif blockID(choseLeftTrials(b)) > 0
                if a > 1
                    plot(1,b,'.','Color',[1 .6 0]);
                else
                    plot(1.5,b,'.','Color',[1 .6 0]);
                end
            end
        end
        if c == colLength
            set(gca,'xtick',[-1 0 1],'ytick',[])
        else
            set(gca,'xtick',[],'ytick',[])
        end
        if a == 1
            ylabel(num2str(contrasts(c)))
        end
        try
            iMax(end+1) = max(max(f1.CData));
            iMin(end+1) = min(min(f1.CData));
        catch
            iMax(end+1) = NaN;
            iMin(end+1) = NaN;
        end
        if a > 1
            xlim([-1.5 1]);
        else
            xlim([-.5 1.5]);
        end
        
        pIdx2 = sub2ind([rowLength colLength],a*2,c);
        ax2 = subplot(colLength,rowLength,pIdx2);
        f2 = imagesc(eventWindow,1:numRightTrials,alignedResps{a}(choseRightTrials,:,plotCells(k)));
        colormap(ax2,rasterColor([0 0 0]));
        hold on;
        line([0 0],[0 100],'LineStyle','-','Color',[0 0 0],'linewidth',1);
        for b = 1:numRightTrials
            if blockID(choseRightTrials(b)) < 0
                if a > 1
                    plot(1,b,'.','Color',[.1 .7 .1]);
                else
                    plot(1.5,b,'.','Color',[.1 .7 .1]);
                end
            elseif blockID(choseRightTrials(b)) > 0
                if a > 1
                    plot(1,b,'.','Color',[1 .6 0]);
                else
                    plot(1.5,b,'.','Color',[1 .6 0]);
                end
            end
        end
        if c == colLength
            set(gca,'xtick',[-1 0 1],'ytick',[])
        else
            set(gca,'xtick',[],'ytick',[])
        end
        try
            iMax(end+1) = max(max(f2.CData));
            iMin(end+1) = min(min(f2.CData));
        catch
            iMax(end+1) = NaN;
            iMin(end+1) = NaN;
        end
        if a > 1
            xlim([-1.5 1]);
        else
            xlim([-.5 1.5]);
        end
    end
end

for a = 1:2
    for c = 1:length(contrasts)
        pIdx1 = sub2ind([rowLength colLength],a*2-1,c);
        ax1 = subplot(colLength,rowLength,pIdx1);
        if a == 1
        ax1.YColor = colors(c,:);
%         ax1.LineWidth = 1.5;
        end
        
        caxis([min(iMin) max(iMax)*cgain]);
        if c == length(contrasts)
            xlabel('Time (s)')
        end
        if c == 1
            title('chose left')
        end
        box off
        set(gca,'tickdir','out')
        
        pIdx2 = sub2ind([rowLength colLength],a*2,c);
        ax2 = subplot(colLength,rowLength,pIdx2);
%         ax2.YColor = colors(c,:);
%         ax2.LineWidth = 1.5;
        caxis([min(iMin) max(iMax)*cgain]);
        if c == length(contrasts)
            xlabel('Time (s)')
        end
        if c == 1
            title('chose right')
        end
        box off
        set(gca,'tickdir','out')
    end
end

 was_a_key = waitforbuttonpress;
    
    if was_a_key && strcmp(get(fig, 'CurrentKey'), 'leftarrow')
      k = max(1, k - 1);
    elseif was_a_key && strcmp(get(fig, 'CurrentKey'), 'rightarrow')
      k = min(max_k, k + 1);
    elseif was_a_key && strcmp(get(fig, 'CurrentKey'), 'return')
        for h = 1:4:33
            sp = get(ha(h),'position');
            sp(1) = col4pos;
            set(ha(h),'position',sp);
        end

        for h = 3:4:35
            sp = get(ha(h),'position');
            sp(1) = col2pos;
            set(ha(h),'position',sp);
        end
        disp(strcat({'k = '},num2str(k)))
        figName = strcat('contrast_dir_block_',expInfo.mouseName,'_',expInfo.expDate,'_cell_',num2str(plotCells(k)));
        printfig(gcf, figName)
        break
    elseif was_a_key && strcmp(get(fig, 'CurrentKey'), 'escape')
        disp(strcat({'k = '},num2str(k)))
        figName = strcat('contrast_dir_block_',expInfo.mouseName,'_',expInfo.expDate,'_cell_',num2str(plotCells(k)));
        close(fig)
        break
    end
end      
            
            
            
            
            
            
            
            
            