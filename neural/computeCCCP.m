%% compute CP and BP
try
    parpool();
catch
end

for m = 1:length(mouseList)
    subpop = 'all';
    ETA = 1;

    mouseName = char(mouseList{m});
    expDate = char(expList{m}{1});
    expNum = expList{m}{2};
    [expRef, ~] = data.constructExpRef(mouseName,expDate,expNum);
    fprintf('Loading %s...',expRef)
    [behavioralData, expInfo, neuralData] = data.loadDataset(mouseName, expDate, expNum);
    expInfo.hemisphere = hemList(m);
    
    trialTypes = getTrialTypes(expInfo, behavioralData, 'late');
    alignedResps = neuralData.eta.alignedResps;
    
    if strcmp(subpop,'contraMov') && expInfo.hemisphere > 0
        plotCells = getWhichCells('leftMov',neuralData);
    elseif strcmp(subpop,'contraMov') && expInfo.hemisphere < 0
        plotCells = getWhichCells('rightMov',neuralData);
    elseif strcmp(subpop,'contraStim') && expInfo.hemisphere > 0
        plotCells = getWhichCells('leftStim',neuralData);
    elseif strcmp(subpop,'contraStim') && expInfo.hemisphere < 0
        plotCells = getWhichCells('rightStim',neuralData);
    elseif strcmp(subpop,'ipsiMov') && expInfo.hemisphere > 0
        plotCells = getWhichCells('rightMov',neuralData);
    elseif strcmp(subpop,'ipsiMov') && expInfo.hemisphere < 0
        plotCells = getWhichCells('leftMov',neuralData);
    elseif strcmp(subpop,'ipsiStim') && expInfo.hemisphere > 0
        plotCells = getWhichCells('rightStim',neuralData);
    elseif strcmp(subpop,'ipsiStim') && expInfo.hemisphere < 0
        plotCells = getWhichCells('leftStim',neuralData);
    elseif strcmp(subpop,'all')
        plotCells = getWhichCells('all',neuralData);
    end

    resps = squeeze(alignedResps{ETA}(:,:,:));
    timeRange = 15:2:35;
    
    tic
    [cp, ~, cpSig, bp, ~, bpSig] = getCP(...
        plotCells, resps, timeRange, ...
        behavioralData, expInfo);
    toc
%     [cp, ~, cp_sig] = getCP(plotCells, resps, timeRange, trialTypes.intVar.all.contrast_direction);
%     [bp, ~, bp_sig] = getCP(plotCells, resps, timeRange, trialTypes.intVar.all.contrast_block);
    
    pcp = sum(cpSig)/length(plotCells);
    pbp = sum(bpSig)/length(plotCells);
    pcp_all(m,:) = pcp;
    pbp_all(m,:) = pbp;
    
    clearvars -except mouseList expList hemList m pcp_all pbp_all eventWindow
end

%% plot

figure;
set(gcf,'position',[465   252   1200   407]);
    
timeRange = 15:2:35;
eventWindow = -2:.1:2;
eventWindow = eventWindow(timeRange);

subplot(1,2,1)
hold on;
line([0 0],[0 1],'LineStyle',':','Color',[.5 .5 .5])
line([0.8 0.8],[0 1],'LineStyle',':','Color',[.5 .5 .5])
line([-.5 1.5],[0.05 0.05],'LineStyle',':','Color',[.5 .5 .5])
xlim([-.4 1.4])
ylim([0 .25]);
yticks([0 .2 .4 .6 .8 1])
set(gca,'tickdir','out')
xlabel('Time – stimOn (s)')
ylabel('Proportion of significant cells')

for m = 1:length(mouseList)
    plot(eventWindow,pcp_all(m,:),'LineWidth',1,'Color',[.85 .85 .85])
end
plot(eventWindow, nanmean(pcp_all,1),'Color',[0 0 0],'LineWidth',2)
% plotSignal(eventWindow, nanmean(pcp_all,1),prctile(pcp_all,97.5),prctile(pcp_all,2.5),[0 0 0],'-')
title('Choice probability')

subplot(1,2,2)
hold on;
line([0 0],[0 1],'LineStyle',':','Color',[.5 .5 .5])
line([0.8 0.8],[0 1],'LineStyle',':','Color',[.5 .5 .5])
line([-.5 1.5],[0.05 0.05],'LineStyle',':','Color',[.5 .5 .5])
xlim([-.4 1.4])
ylim([0 .25]);
yticks([0 .2 .4 .6 .8 1])
set(gca,'tickdir','out')
xlabel('Time – stimOn (s)')
ylabel('Proportion of significant cells')

for m = 1:length(mouseList)
    plot(eventWindow,pbp_all(m,:),'LineWidth',1,'Color',[.85 .85 .85])
end
plot(eventWindow, nanmean(pbp_all,1),'Color',[0 0 0],'LineWidth',2)
% plotSignal(eventWindow, nanmean(pbp_all,1),prctile(pbp_all,97.5),prctile(pbp_all,2.5),[0 0 0],'-')
title('Block probability')

%% save data

ccCP.stimAligned = struct('expRef',[],'cp',[],'bp',[],'vp',[]);

ccCP.stimAligned(2).expRef = expRef;




