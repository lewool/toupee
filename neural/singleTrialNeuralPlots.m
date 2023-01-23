% [expInfo, neuralData, behavioralData] = processExperiment(expInfo);
% neuralData = getCellPCs(neuralData);
% [neuralData] = alignPCs(expInfo, neuralData, behavioralData);
eyeData = getEyeData(expInfo);
[eyeData] = alignFace(expInfo, eyeData, behavioralData);

%%

cd('C:\Users\Wool\Documents\GitHub\rastermap\matlab')
csv = neuralData.cellResps';
[isort1, isort2, Sm] = mapTmap(csv);

%%

eventWindow = neuralData.eta.eventWindow;
contrasts = getUniqueContrasts(expInfo);

[~, whichTrials] = selectCondition(...
    expInfo, contrasts, behavioralData, initTrialConditions('movementTime','late','movementDir','all'));
evenTrials = whichTrials(2:2:end);
oddTrials = whichTrials(1:2:end);


useCV = 0;

if useCV == 1
    sortTrials = evenTrials;
    plotTrials = oddTrials;
else
    sortTrials = whichTrials;
    plotTrials = whichTrials;
end

% sort cells by mean firing time
ETA = 1;
meanPSTH = squeeze(nanmean(neuralData.eta.alignedResps{ETA}(sortTrials,21:28,:),1))';
% minPSTH = meanPSTH - min(meanPSTH,[],2);
normPSTH = meanPSTH./sum(meanPSTH,2);

countingVector = linspace(0,.7,8);
for c = 1:length(normPSTH)
    meanFiringTime(c,1) = dot(normPSTH(c,:),countingVector);
end

cellRange = 1:length(normPSTH);
% cellRange = 1000:1200;
setminy = 0.02;
setmaxy = 0.2;
% [~,sortIdx] = sort(meanFiringTime,'ascend');
sortIdx = isort1;

% plot all cell rasters per trial

fig100 = figure;
plotLength = 15;
vidHeight = 3;
pcColors = customColormap([.5 0 .4],[1 .5 0], 6);
close(fig100)
tracelabels = {...
    '';
    '\it neural';
    '\it pupil';
    '\it whisk';
    '\it lick';
    '\it paw';
    '\it wheel'};
    
for ETA = [1 2 3]
    
    for a = 1:size(eyeData.eta.alignedFrames{ETA},1)
        alignedFrames_int{ETA}(a,:) = interp1(eyeData.eta.eventWindow, eyeData.eta.alignedFrames{ETA}(a,:), neuralData.eta.eventWindow);
    end
    
    if ETA == 1
        fig1 = figure;
        xLimit = [-.5 2];
        set(fig1,'position',[160 -225 330 1550])
        set(fig1,'color','w');
    elseif ETA == 2
        fig2 = figure;
        xLimit = [-2 2];
        set(fig2,'position',[495 -225 330 1550])
        set(fig2,'color','w');
    elseif ETA == 3
        fig4 = figure;
        xLimit = [-2 2];
        set(fig4,'position',[830 -225 330 1550])
        set(fig4,'color','w');
    end
    
    
    trialContrasts = expInfo.block.events.contrastValues(whichTrials);
    trialChoices = behavioralData.wheelMoves.epochs(5).moveDir(whichTrials);
    trialRTs = behavioralData.wheelMoves.epochs(5).onsetTimes(whichTrials) - behavioralData.eventTimes(1).daqTime(whichTrials);
    cueTime = behavioralData.eventTimes(2).daqTime(whichTrials)  - behavioralData.eventTimes(1).daqTime(whichTrials);
    feedbackTime = behavioralData.eventTimes(5).daqTime(whichTrials) - behavioralData.eventTimes(1).daqTime(whichTrials);

    meanPupilActivity = squeeze(nanmean(eyeData.eta.alignedFace{ETA}(plotTrials,:,1),1))';
    stdPupilActivity = squeeze(nanstd(eyeData.eta.alignedFace{ETA}(plotTrials,:,1),[],1))';
    semPupilActivity = stdPupilActivity./sqrt(length(plotTrials));

    meanWhiskActivity = squeeze(nanmean(eyeData.eta.alignedFace{ETA}(plotTrials,:,2),1))';
    stdWhiskActivity = squeeze(nanstd(eyeData.eta.alignedFace{ETA}(plotTrials,:,2),[],1))';
    semWhiskActivity = stdWhiskActivity./sqrt(length(plotTrials));

    meanLickActivity = squeeze(nanmean(eyeData.eta.alignedFace{ETA}(plotTrials,:,3),1))';
    stdLickActivity = squeeze(nanstd(eyeData.eta.alignedFace{ETA}(plotTrials,:,3),[],1))';
    semLickActivity = stdLickActivity./sqrt(length(plotTrials));

    meanPawActivity = squeeze(nanmean(eyeData.eta.alignedFace{ETA}(plotTrials,:,4),1))';
    stdPawActivity = squeeze(nanstd(eyeData.eta.alignedFace{ETA}(plotTrials,:,4),[],1))';
    semPawActivity = stdPawActivity./sqrt(length(plotTrials));

    [alignedWheel, wheelWindow] = getAlignedWheelTraces(expInfo, behavioralData, ETA);
    meanWheelVel = nanmean(abs(alignedWheel(plotTrials,:)));
    stdWheelVel = nanstd(alignedWheel(plotTrials,:));
    semWheelVel = stdWheelVel./sqrt(length(plotTrials));
%     meanPCs = squeeze(nanmean(neuralData.eta.alignedPCs{ETA}(plotTrials,:,1:6),1))';
%     stdPCs = squeeze(std(neuralData.eta.alignedPCs{ETA}(plotTrials,:,1:6),1))';


    
    for t = 1:length(plotTrials)
        try
            tmp1 = read(eyeData.veye,[alignedFrames_int{ETA}(plotTrials(t),21)]);
            tmp2 = read(eyeData.veye,[alignedFrames_int{ETA}(plotTrials(t),21)+1]);
            meanFrame(:,:,t) = read(eyeData.veye,[alignedFrames_int{ETA}(plotTrials(t),21)]);
        catch
            tmp1 = nan;
            tmp2 = nan;
            meanFrame(:,:,t) = nan;
        end
        
        meanMove(:,:, t) = abs(tmp2 - tmp1);
    end
    
    ax(1) =subplot(plotLength,1,[1 vidHeight]);
    tmap = mean(meanMove,3);
    fmap = mean(meanFrame,3);
    red = cat(3, ones(size(tmap)), zeros(size(tmap)), zeros(size(tmap)))*125;
    imshow(fmap(:,:,:,1));
    caxis([0 125]);
    hold on;
    h = imshow(red);
    hold off
    set(h, 'AlphaData', tmap/15)
    axis off
    title('Mean(trial condition)')
    
    ax(2) = subplot(plotLength,1,[vidHeight+1:plotLength-6]);
    meanTrialActivity = squeeze(nanmean(neuralData.eta.alignedResps{ETA}(plotTrials,:,sortIdx),1))';
    meanPopActivity = squeeze(nanmean(neuralData.eta.alignedResps{ETA}(plotTrials,:,sortIdx),3))';
    imagesc(eventWindow,1:length(normPSTH),meanTrialActivity);
    colormap(flipud(gray));
    hold on;
    caxis([0.05 .5])

    ax(3) = subplot(plotLength,1,[plotLength-5]);
    plotPSTHs(neuralData.eta.eventWindow, mean(meanPopActivity,2)', std(meanPopActivity,[],2)',[0 0 0],'-')

    ax(4) = subplot(plotLength,1,[plotLength-4]);
    plotPSTHs(eyeData.eta.eventWindow, meanPupilActivity', stdPupilActivity',[.9 0 0],'-')

    ax(5) = subplot(plotLength,1,[plotLength-3]);
    plotPSTHs(eyeData.eta.eventWindow, meanWhiskActivity', stdWhiskActivity',[1 .5 0],'-')

    ax(6) = subplot(plotLength,1,[plotLength-2]);
    plotPSTHs(eyeData.eta.eventWindow, meanLickActivity', stdLickActivity',[.25 .75 0],'-')

    ax(7) = subplot(plotLength,1,[plotLength-1]);
    plotPSTHs(eyeData.eta.eventWindow, meanPawActivity', stdPawActivity',[0 .25 .75],'-')

    ax(8) = subplot(plotLength,1,[plotLength]);
    plotPSTHs(wheelWindow, (meanWheelVel), stdWheelVel,[.25 0 .75],'-')  
    
    for a = 2:length(ax)
        subplot(ax(a))
        if a == 8
            miny = -30; maxy = 30;
            if ETA == 1
                xlabel('Time from stimulus onset')
            elseif ETA == 2
                xlabel('Time from move onset')
            elseif ETA == 3
                xlabel('Time from feedback onset')
            end
        elseif a == 3
            miny = .05; maxy = .15;
            set(gca, 'XTickLabels', {})
        elseif a == 4
            miny = -1; maxy = 1.5;
            set(gca, 'XTickLabels', {})
        elseif a == 7 || a == 5 || a == 6
            miny = -1; maxy = 3;
            set(gca, 'XTickLabels', {})
        elseif a == 2
            miny = 1; maxy = size(meanTrialActivity,1);
            set(gca, 'XTickLabels', {})
            set(gca, 'YTickLabels', {})
        end
        
        text(xLimit(1)+.05,maxy,tracelabels{a-1})
        box off
        set(gca,'tickdir','out')
        xlim(xLimit);
        ylim([miny maxy]);
        hold on;
        line([0 0],[miny maxy],'LineStyle','--','Color','k');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig = figure;
set(fig,'position',[1165 -225 330 1550])
set(fig,'color','w');

% t = 349;
t=1;
max_t = length(whichTrials);
max_i = size(alignedFrames_int{1},2);
min_i = 16;
i = 16;
ihop = 1;
h1 = line;
h2 = line;
h3 = line;
h4 = line;
h5 = line;
h6 = line;
h7 = line;
h8 = line;
h9 = line;
h10 = line;
h11 = line;
h12 = line;
h13 = line;
h14 = line;

ax2(1) = subplot(plotLength,1,[1 vidHeight]);
ax2(2) = subplot(plotLength,1,[vidHeight+1:plotLength-6]);
ax2(3) = subplot(plotLength,1,[plotLength-5]);
ax2(4) = subplot(plotLength,1,[plotLength-4]);
ax2(5) = subplot(plotLength,1,[plotLength-3]);
ax2(6) = subplot(plotLength,1,[plotLength-2]);
ax2(7) = subplot(plotLength,1,[plotLength-1]);
ax2(8) = subplot(plotLength,1,[plotLength]);

while t <= max_t && i <= max_i
    
    for a = 1:length(ax2)
        subplot(ax2(a))
        cla;
    end
        
    subplot(ax2(1));
    tmp = read(eyeData.veye,[alignedFrames_int{1}(whichTrials(t),i)]);
    tmp1 = read(eyeData.veye,[alignedFrames_int{1}(whichTrials(t),i)-1]);
    tmp2 = read(eyeData.veye,[alignedFrames_int{1}(whichTrials(t),i)]);
    tmap = abs(tmp2(:,:,:,1)-tmp1(:,:,:,1));
    red = cat(3, ones(size(tmap)), zeros(size(tmap)), zeros(size(tmap)));
    imshow(tmp(:,:,:,1));
    hold on;
    h = imshow(red);
    hold off
    set(h, 'AlphaData', tmap*2)
    caxis([0 100]);
    axis off
%     axis equal
    title(strcat({'Trial: '},num2str(whichTrials(t)),{'; Contrast: '},num2str(trialContrasts(t)),{'; Choice: '},num2str(trialChoices(t))))
    
    faceActivity = squeeze(eyeData.eta.alignedFace{1}(whichTrials(t),:,:))';
    for f = 1:size(faceActivity,1)
        faceDot(f,:) = interp1(eyeData.eta.eventWindow, smooth(faceActivity(f,:)), neuralData.eta.eventWindow);
    end
    wheelActivity = behavioralData.wheelMoves.traces.pos{whichTrials(t)};
    wheelTime = behavioralData.wheelMoves.traces.time{whichTrials(t)} - behavioralData.eventTimes(1).daqTime(whichTrials(t));
    wheelDot = interp1(wheelTime,wheelActivity, neuralData.eta.eventWindow);
    trialActivity = squeeze(neuralData.eta.alignedResps{1}(whichTrials(t),:,sortIdx))';
%     pcActivity = squeeze(neuralData.eta.alignedPCs{1}(whichTrials(t),:,1:6))';
    
    subplot(ax2(2));
    imagesc(eventWindow,1:length(cellRange),trialActivity(cellRange,:));
    colormap(flipud(gray));
    hold on;
    caxis([-.0 .9])
    ylabel('Cells')
    h1 = line(...
        [neuralData.eta.eventWindow(i) neuralData.eta.eventWindow(i)],...
        [0 length(cellRange)],...
        'LineWidth',1.5,'LineStyle','-','Color','r');
    
    subplot(ax2(3));
    plotPSTHs(neuralData.eta.eventWindow, mean(trialActivity(cellRange,:)), 0,[0 0 0],'-')
    hold on;
    h2 = plot(...
        [neuralData.eta.eventWindow(i) neuralData.eta.eventWindow(i)],...
        [mean(trialActivity(cellRange,i)) mean(trialActivity(cellRange,i))],...
        'Marker','o','MarkerFaceColor','k','MarkerEdgeColor','none');
    
    subplot(ax2(4));
    plotPSTHs(eyeData.eta.eventWindow, smooth(faceActivity(1,:))', zeros(size(faceActivity(1,:))),[.9 0 0],'-')
    hold on;
    h3 = plot(...
        [neuralData.eta.eventWindow(i) neuralData.eta.eventWindow(i)],...
        [faceDot(1,i) faceDot(1,i)],...
        'Marker','o','MarkerFaceColor',[.9 0 0],'MarkerEdgeColor','none');
    
    subplot(ax2(5));
    plotPSTHs(eyeData.eta.eventWindow, smooth(faceActivity(2,:))', zeros(size(faceActivity(2,:))),[1 .5 0],'-')
    hold on;
    h4 = plot(...
        [neuralData.eta.eventWindow(i) neuralData.eta.eventWindow(i)],...
        [faceDot(2,i) faceDot(2,i)],...
        'Marker','o','MarkerFaceColor',[1 0.5 0],'MarkerEdgeColor','none');
    
    subplot(ax2(6));
    plotPSTHs(eyeData.eta.eventWindow, smooth(faceActivity(3,:))', zeros(size(faceActivity(3,:))),[.25 .75 0],'-')
    hold on;
    h5 = plot(...
        [neuralData.eta.eventWindow(i) neuralData.eta.eventWindow(i)],...
        [faceDot(3,i) faceDot(3,i)],...
        'Marker','o','MarkerFaceColor',[.25 .75 0],'MarkerEdgeColor','none');
    
    subplot(ax2(7));
    plotPSTHs(eyeData.eta.eventWindow, smooth(faceActivity(4,:))', zeros(size(faceActivity(4,:))),[0 .25 .75],'-')
    hold on;
    h6 = plot(...
        [neuralData.eta.eventWindow(i) neuralData.eta.eventWindow(i)],...
        [faceDot(4,i) faceDot(4,i)],...
        'Marker','o','MarkerFaceColor',[0 .25 .75],'MarkerEdgeColor','none');
    
    subplot(ax2(8));
    plotPSTHs(wheelTime, (wheelActivity), zeros(size(wheelActivity)),[.25 0 .75],'-')
    hold on;
    text(-.45,max(40),'\it wheel')
    h7 = plot(...
        [neuralData.eta.eventWindow(i) neuralData.eta.eventWindow(i)],...
        [wheelDot(i) wheelDot(i)],...
        'Marker','o','MarkerFaceColor',[.25 0 .75],'MarkerEdgeColor','none');
    
    for a = 2:length(ax2)
        subplot(ax2(a))
        if a == 2
            miny2 = 1; maxy2 = size(trialActivity(cellRange,:),1);
            set(gca, 'XTickLabels', {});
        elseif a == 3
            miny2 = setminy; maxy2 = setmaxy;
            set(gca, 'XTickLabels', {})
        elseif a == 8
            miny2 = -40; maxy2 = 40;
            xlabel('Time from stimulus onset')
        else
            miny2 = min(faceActivity(a-3,:)); maxy2 = max(faceActivity(a-3,:));
            set(gca, 'XTickLabels', {})
        end
        text(-.45,maxy2,tracelabels{a-1})
        line([0 0],[miny2 maxy2],'LineStyle','--','Color','k');
        line([trialRTs(t) trialRTs(t)],[miny2 maxy2],'LineStyle','--','Color','k');
        line([cueTime(t) cueTime(t)],[miny2 maxy2],'LineStyle','--','Color','k');
        line([feedbackTime(t) feedbackTime(t)],[miny2 maxy2],'LineStyle','--','Color','k');
        box off
        set(gca,'tickdir','out')
        xlim([-.5 2]);
        ylim([miny2 maxy2]);
    end
        
    was_a_key = waitforbuttonpress;
    
    if was_a_key && strcmp(get(fig, 'CurrentKey'), 'leftarrow')
        i = max(min_i, i - ihop);
    elseif was_a_key && strcmp(get(fig, 'CurrentKey'), 'rightarrow')
        i = min(max_i, i + ihop);
%         printfig(fig,strcat('vid',num2str(i)));
    elseif was_a_key && strcmp(get(fig, 'CurrentKey'), 'space')
        i = min_i;
    elseif was_a_key && strcmp(get(fig, 'CurrentKey'), 'uparrow')
        t = max(1, t - 1);
    elseif was_a_key && strcmp(get(fig, 'CurrentKey'), 'downarrow')
        t = min(max_t, t + 1);
    elseif was_a_key && strcmp(get(fig, 'CurrentKey'), 'escape')
        close(fig1);
        close(fig2);
        close(fig4);
        close(fig);
        break
    elseif was_a_key && strcmp(get(fig, 'CurrentKey'), 'return')
%         break
        printfig(fig,strcat('vid_trial',num2str(t),'_',num2str(i)));
    end
end    

%%

whichCells = 'patient'; %choose from 'pLabels' array
if strcmp(whichCells, 'all')
    plotCells = (1:size(neuralData.eta.alignedResps{1},3))';
else
    plotCells = find(neuralData.stats.pValues(:,strcmp(neuralData.stats.labels,whichCells)) < 0.01);
end

for p = 1:length(plotCells)
    idx(p) = find(isort1 == plotCells(p));
end

if strcmp(whichCells,'leftStim') || strcmp(whichCells,'rightStim')
    overlay = [0 1 .2];
elseif strcmp(whichCells,'leftMov')
    if hemisphere > 0
        overlay = [0 .4 1];
    else
        overlay = [1 0 0];
    end
elseif strcmp(whichCells,'rightMov')
    if hemisphere < 0
        overlay = [0 .4 1];
    else
        overlay = [1 0 0];
    end
elseif strcmp(whichCells,'value')
    overlay = [1 0 1];
elseif strcmp(whichCells,'hit')
    overlay = [0 .7 0];
elseif strcmp(whichCells,'miss')
    overlay = [.5 0 0];
elseif strcmp(whichCells,'impulsive')
    overlay = [1 .5 0];
elseif strcmp(whichCells,'patient')
    overlay = [0 .5 .5];
end
    

    
if hemisphere > 0 
    con = contrasts;
else
    con = fliplr(contrasts);
end

fig2 = figure;
set(fig2,'position',[106 476 2408 1150])
for c = 1: length(con)
    [~, whichTrials] = selectCondition(expInfo, con(c), behavioralData, initTrialConditions('movementTime','late','responseType','correct'));
    meanTrialActivity = flipud(squeeze(nanmean(neuralData.eta.alignedResps{1}(whichTrials,:,isort1),1))');
    subplot(1,length(con),c)
    imagesc(eventWindow,1:size(neuralData.eta.alignedResps{1},3),smoothdata(meanTrialActivity));
    colormap(flipud(gray));
    hold on;
    line([0 0],[0 size(neuralData.eta.alignedResps{1},3)],'LineStyle','--','Color','k');
    xlim([-.5 2]);
    caxis([-0 .5])
    if c == 1
        ylabel('Cells')
        xlabel('Time from stimulus onset')
    else axis off
    end
    
    box off 
    set(gca,'tickdir','out')
    title(num2str(contrasts(c)))
    
    for p = 1:length(plotCells)
        lh = line([-.5 2],[idx(p) idx(p)],'Color',overlay);
        lh.Color(4) = .3;
    end
end


%%
for p = 1:length(plotCells)
lh = line([-.5 2],[idx(p) idx(p)],'Color','g');
lh.Color(4) = .1;
end
%%


contrasts = getUniqueContrasts(expInfo);
zeroIdx = find(contrasts == 0);
walkup = length(contrasts) - zeroIdx;
walkback = zeroIdx - 1;
allColors = [0 0 .5; 0 0 1; 0 .4 1; .6 .8 1; .75 .75 .75; .8 .45 .45; 1 0 0; .5 0 0; .25 0 0];
zeroGray = find(allColors(:,1) == .75);
contrastColors = allColors(zeroGray-walkback:zeroGray + walkup,:);

fig3 = figure;
set(fig3,'position',[106 476 1140 316])
for c = 1: length(contrasts)
    for s = 1:3
    [~, whichTrials] = selectCondition(expInfo, contrasts(c), behavioralData, initTrialConditions('movementTime','late','movementDir','ccw'));
    plotPC1 = squeeze(neuralData.eta.alignedPCs{s}(whichTrials,:,1))';
    plotPC2 = squeeze(neuralData.eta.alignedPCs{s}(whichTrials,:,2))';
    subplot(1,3,s);
%     plot(mean(plotPC1,2),mean(plotPC2,2),'LineWidth',2,'Color',contrastColors(c,:))
    plot(eventWindow,mean(plotPC1,2),'LineWidth',2,'Color',contrastColors(c,:))
    hold on;

    
    box off 
    set(gca,'tickdir','out')
    title(num2str(contrasts(c)))
    end
end

%%

stimEpochAll = neuralData.eta.alignedResps{1}(:,21:28,:);
stimCat = nan(1957,909*8);
for c = 1:size(stimEpochAll,3)
    for t = 1:size(stimEpochAll,1)
        counter = (t-1)*8+1;
        stimCat(c,counter:counter+7) = stimEpochAll(t,:,c);
    end
end

%% WORKBENCH

subplot(plotLength,2,[2*plotLength-10]);
plotPSTHs(neuralData.eta.eventWindow,meanPCs(1,:), stdPCs(1,:),pcColors(1,:),'-')
miny = min(meanPCs(1,:)-stdPCs(1,:));
maxy = max(meanPCs(1,:)+stdPCs(1,:));
box off
set(gca,'tickdir','out')
xlim(xLimit);
ylim([miny maxy]);
hold on;
line([0 0],[miny maxy],'LineStyle','--','Color','k');
set(gca, 'XTickLabels', {})
set(gca, 'YTickLabels', {})
text(xLimit(1)+.05,maxy,'\it PC1')

subplot(plotLength,plotWidth,[plotWidth*plotLength-10]);
cla;
plot(neuralData.eta.eventWindow,pcActivity(1,:),'LineWidth',2,'Color',pcColors(1,:))
box off
set(gca,'tickdir','out')
xlim([-.5 2]);
ylim([min(pcActivity(1,:)) max(pcActivity(1,:))]);
hold on;
line([0 0],[min(pcActivity(1,:)) max(pcActivity(1,:))],'LineStyle','--','Color','k');
line([trialRTs(t) trialRTs(t)],[min(pcActivity(1,:)) max(pcActivity(1,:))],'LineStyle','--','Color','k');
line([cueTime(t) cueTime(t)],[min(pcActivity(1,:)) max(pcActivity(1,:))],'LineStyle','--','Color','k');
line([feedbackTime(t) feedbackTime(t)],[min(pcActivity(1,:)) max(pcActivity(1,:))],'LineStyle','--','Color','k');
set(gca, 'XTickLabels', {})
set(gca, 'YTickLabels', {})
text(-.45,max(pcActivity(1,:)),'\it PC1')
delete(h9);
h9 = plot(...
    [neuralData.eta.eventWindow(i) neuralData.eta.eventWindow(i)],...
    [pcActivity(1,i) pcActivity(1,i)],...
    'Marker','o','MarkerFaceColor',pcColors(1,:),'MarkerEdgeColor','none');
    
