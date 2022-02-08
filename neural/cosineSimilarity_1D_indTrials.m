function cosineSimilarity_1D_indTrials(expInfo, behavioralData, neuralData, whichDim, pickTrials)

h =  findobj('type','figure');
n = length(h);

for iX = 1:length(expInfo)

%% set up some plot stuff 

% contrasts = getUniqueContrasts(expInfo(iX));
% zeroIdx = find(contrasts == 0);
% walkup = length(contrasts) - zeroIdx;
% walkback = zeroIdx - 1;
% allColors = [0 0 .5; 0 0 1; 0 .4 1; .6 .8 1; .75 .75 .75; .8 .45 .45; 1 0 0; .5 0 0; .25 0 0];
% 
% zeroGray = find(allColors(:,1) == .75);
% colors = allColors(zeroGray-walkback:zeroGray + walkup,:);

ms = 5;
for o = 1:2
%% extract trial IDs
%organize the trials into conditions
trialTypes = getTrialTypes(expInfo(iX), behavioralData(iX), 'late');
earlyTrialTypes = getTrialTypes(expInfo(iX), behavioralData(iX), 'early');

switch pickTrials
    case 'contrast'
        whichTrials = [];
        for c = 1:length(trialTypes.intVar.all.contrast_outcome)
            whichTrials{c} = trialTypes.intVar.all.contrast_direction{c,o};
        end
        
        contrasts = getUniqueContrasts(expInfo(iX));
        zeroIdx = find(contrasts == 0);
        walkup = length(contrasts) - zeroIdx;
        walkback = zeroIdx - 1;
        allColors = [0 0 .5; 0 0 1; 0 .4 1; .6 .8 1; .75 .75 .75; .8 .45 .45; 1 0 0; .5 0 0; .25 0 0];
        zeroGray = find(allColors(:,1) == .75);
        colors = allColors(zeroGray-walkback:zeroGray + walkup,:);
        
    case 'binned contrast'
        whichTrials = [];
        for b = 1:4
            whichTrials{1} = cat(2,trialTypes.intVar.all.contrast_outcome{1,o},trialTypes.intVar.all.contrast_outcome{2,o});
            whichTrials{2} = cat(2,trialTypes.intVar.all.contrast_outcome{3,o},trialTypes.intVar.all.contrast_outcome{4,o});
            whichTrials{3} = cat(2,trialTypes.intVar.all.contrast_outcome{6,o},trialTypes.intVar.all.contrast_outcome{7,o});
            whichTrials{4} = cat(2,trialTypes.intVar.all.contrast_outcome{8,o},trialTypes.intVar.all.contrast_outcome{9,o});
        end
        
        contrasts = getUniqueContrasts(expInfo(iX));
        colors = [0 .4 1; .6 .8 1; .8 .45 .45; 1 0 0];
        
    case 'low contrast'
        whichTrials = [];
        for b = 1:4
            whichTrials{1} = cat(2,trialTypes.intVar.all.contrast_outcome{3,o},trialTypes.intVar.all.contrast_outcome{4,o});
            whichTrials{2} = cat(2,trialTypes.intVar.all.contrast_outcome{6,o},trialTypes.intVar.all.contrast_outcome{7,o});
        end
        
        contrasts = getUniqueContrasts(expInfo(iX));
        colors = [0 .4 1; 1 0 0];
       
        
    case 'early contrast'
        whichTrials = [];
        for c = 1:length(earlyTrialTypes.intVar.all.contrast_outcome)
            whichTrials{c} = earlyTrialTypes.intVar.all.contrast_outcome{c,1};
        end
        
        contrasts = getUniqueContrasts(expInfo(iX));
        zeroIdx = find(contrasts == 0);
        walkup = length(contrasts) - zeroIdx;
        walkback = zeroIdx - 1;
        allColors = [0 0 .5; 0 0 1; 0 .4 1; .6 .8 1; .75 .75 .75; .8 .45 .45; 1 0 0; .5 0 0; .25 0 0];
        zeroGray = find(allColors(:,1) == .75);
        colors = allColors(zeroGray-walkback:zeroGray + walkup,:);

    case 'velocity'
        whichTrials = [];
        [velTrialIdx, ~] = getWheelVelPercentiles(expInfo(1), behavioralData(1),'late');
        for v = 1:length(velTrialIdx)
            whichTrials{v} = intersect(velTrialIdx{v},trialTypes.singleVar.outcome{1});
        end
        
        velColorMap = [linspace(.75, 0, length(whichTrials))', linspace(0, .75, length(whichTrials))', linspace(.75, 0, length(whichTrials))'];
        colors = velColorMap;
        
    case 'block'
        whichTrials{1} = intersect([trialTypes.singleVar.block{1}],trialTypes.singleVar.outcome{1});
        whichTrials{2} = intersect([trialTypes.singleVar.block{2}],trialTypes.singleVar.outcome{1});
        colors = [0 0 0.5; 0.5 0 0];
        
    case 'side'
        whichTrials{1} = trialTypes.intVar.all.side_outcome{1,o};
        whichTrials{2} = trialTypes.intVar.all.side_outcome{3,o};
        colors = [0 0.4 1; 1 0 0];
        
    case 'early side'
        whichTrials{1} = earlyTrialTypes.intVar.all.side_outcome{1,o};
        whichTrials{2} = earlyTrialTypes.intVar.all.side_outcome{3,o};
        colors = [0 0.4 1; 1 0 0];
        
    case 'outcome'
        whichTrials{1} = trialTypes.singleVar.outcome{1};
        whichTrials{2} = trialTypes.singleVar.outcome{2};
        colors = [0 0.5 0; 0.5 0 0];
end


%% set up mean vectors (aka 'dimensions')

%get epochs
[baselineResps, stimResps, pmovResps, movResps, rewResps, preCueResps] = getEpochResps(neuralData(iX).eta);

eventWindow = neuralData(iX).eta.eventWindow;

%define test and train sets (cv)
testTrials = 1:2:size(baselineResps,1);
trainTrials = 2:2:size(baselineResps,1);

%initialize mean vectors
mean_dim1 = nan(1,length(baselineResps));

%extract dims
dim1 = whichDim{1};

switch dim1
    case 'stimulus'
    %extract trials that go into a stimulus dimension
    dim1a_idx = [trialTypes.singleVar.contrast{1} trialTypes.singleVar.contrast{2}];
    dim1b_idx = [trialTypes.singleVar.contrast{8} trialTypes.singleVar.contrast{9}];
    mean_dim1a = nanmean(stimResps(intersect(dim1a_idx,trainTrials),:),1);
    mean_dim1b = nanmean(stimResps(intersect(dim1b_idx,trainTrials),:),1);
    
    case 'choice'
    %extract trials that go into a choice dimension
    dim1a_idx = cat(2,trialTypes.singleVar.direction{1});
    dim1b_idx = cat(2,trialTypes.singleVar.direction{2});
%     dim1a_idx = cat(2,trialTypes.intVar.cb2D.contrast_direction{:,1});
%     dim1b_idx = cat(2,trialTypes.intVar.cb2D.contrast_direction{:,2});
    mean_dim1a = nanmean(movResps(intersect(dim1a_idx,trainTrials),:),1);
    mean_dim1b = nanmean(movResps(intersect(dim1b_idx,trainTrials),:),1);
    
    case 'premove choice'
    %extract trials that go into a choice dimension
    dim1a_idx = cat(2,trialTypes.singleVar.direction{1});
    dim1b_idx = cat(2,trialTypes.singleVar.direction{2});
%     dim1a_idx = cat(2,trialTypes.intVar.cb2D.contrast_direction{:,1});
%     dim1b_idx = cat(2,trialTypes.intVar.cb2D.contrast_direction{:,2});
    mean_dim1a = nanmean(preCueResps(intersect(dim1a_idx,trainTrials),:),1);
    mean_dim1b = nanmean(preCueResps(intersect(dim1b_idx,trainTrials),:),1);
    
    case 'premove stimulus'
    %extract trials that go into a choice dimension
    dim1a_idx = [trialTypes.singleVar.contrast{1} trialTypes.singleVar.contrast{2}];
    dim1b_idx = [trialTypes.singleVar.contrast{8} trialTypes.singleVar.contrast{9}];
    mean_dim1a = nanmean(preCueResps(intersect(dim1a_idx,trainTrials),:),1);
    mean_dim1b = nanmean(preCueResps(intersect(dim1b_idx,trainTrials),:),1);

    case 'outcome'
    %extract trials that go into an outcome dimension
    dim1a_idx = cat(2,trialTypes.singleVar.outcome{2});
    dim1b_idx = cat(2,trialTypes.singleVar.outcome{1});
    mean_dim1a = nanmean(rewResps(intersect(dim1a_idx,trainTrials),:),1);
    mean_dim1b = nanmean(rewResps(intersect(dim1b_idx,trainTrials),:),1);

    case 'block'
    %extract trials that go into a block dimension
    dim1a_idx = intersect([trialTypes.singleVar.block{1}],trialTypes.singleVar.outcome{1});
    dim1b_idx = intersect([trialTypes.singleVar.block{2}],trialTypes.singleVar.outcome{1});
%     dim1a_idx = trialTypes.intVar.cb3D.block{1};
%     dim1b_idx = trialTypes.intVar.cb3D.block{2};
    mean_dim1a = nanmean(rewResps(intersect(dim1a_idx,trainTrials),:),1);
    mean_dim1b = nanmean(rewResps(intersect(dim1b_idx,trainTrials),:),1);
end

%% compute dot products

whichETA = 2; %use stimulus-aligned


for iT = 1:size(eventWindow,2)
    whichResps = neuralData(iX).eta.alignedResps{whichETA}(:,iT,:);

    fullNorm = 1; % 1 if just angle, 0 if angle and relative magnitude
    for iTrial = testTrials
        if fullNorm == 1
            tNorm = norm(whichResps(iTrial,:));
        else
            tNorm = 1;
        end
        dotProd_dim1a(iTrial,iT) = dot(whichResps(iTrial,:),mean_dim1a)/(tNorm*norm(mean_dim1a));
        dotProd_dim1b(iTrial,iT) = dot(whichResps(iTrial,:),mean_dim1b)/(tNorm*norm(mean_dim1b));
    end

end

dotProd_dim1a(trainTrials,:) = NaN;
dotProd_dim1b(trainTrials,:) = NaN;


%% plot all trials separately
figure(n+1);
if expInfo.hemisphere > 0
    for c = 2:length(contrasts)-1
        subplot(3,length(contrasts),sub2ind([length(contrasts) 3],c,o))
        try
        plot(eventWindow,...
            (dotProd_dim1b(whichTrials{c},:)- dotProd_dim1a(whichTrials{c},:))',...
            'Color',colors(c,:))
        maxY(c) = max(max(dotProd_dim1b(whichTrials{c},:)- dotProd_dim1a(whichTrials{c},:)));
        minY(c) = min(min(dotProd_dim1b(whichTrials{c},:)- dotProd_dim1a(whichTrials{c},:)));
        end
    end
elseif expInfo.hemisphere < 0
    for c = 1:length(contrasts)
        subplot(3,length(contrasts),sub2ind([length(contrasts) 3],c,o))
        plot(...
            (dotProd_dim1b(whichTrials{end-c+1},:)- dotProd_dim1a(whichTrials{end-c+1},:))',...
            'Color',colors(c,:),...
            'LineWidth',1)
        maxY(c) = max(max(dotProd_dim1b(whichTrials{end-c+1},:)- dotProd_dim1a(whichTrials{end-c+1},:)));
        minY(c) = min(min(dotProd_dim1b(whichTrials{end-c+1},:)- dotProd_dim1a(whichTrials{end-c+1},:)));
    end
end

lim = max(abs([maxY minY]));
for c = 2:length(contrasts)-1
    subplot(3,length(contrasts),sub2ind([length(contrasts) 3],c,o))
    line([-2 2],[0 0],'LineStyle','--','Color',[.5 .5 .5]);
    line([0 0],[-lim lim],'LineStyle','--','Color',[.5 .5 .5]);
    ylim([-lim lim]);
    xlim([-1 2])
end

%% split the trials dotProds into two groups based on trial identity

eventWindow = neuralData(iX).eta.eventWindow;

for iT = 1:length(eventWindow)
    for c = 1:length(whichTrials)
        if expInfo(iX).hemisphere < 0 && ~strcmp(dim1,'outcome')
            meanTrials_dim1(c,iT,iX) = nanmedian(dotProd_dim1a(whichTrials{c},iT) - dotProd_dim1b(whichTrials{c},iT));
        else
            meanTrials_dim1(c,iT,iX) = nanmedian(dotProd_dim1b(whichTrials{c},iT) - dotProd_dim1a(whichTrials{c},iT));
        end
    end
end

if expInfo(iX).hemisphere < 0
    meanTrials_dim1(:,:,iX) = flipud(meanTrials_dim1(:,:,iX));
end
%%
figure(n+iX);
hold on
if o == 1
    for c = 2:length(contrasts)-1
        subplot(3,length(contrasts),sub2ind([length(contrasts) 3],c,3))
        hold on
        try
        a = plot(eventWindow,...
            smooth(meanTrials_dim1(c,:,iX),1),...
            'Color',colors(c,:),...
            'LineStyle','-');
        set(a,'LineWidth',2);
        maxYm(c) = max(max(meanTrials_dim1));
        minYm(c) = min(min(meanTrials_dim1));
        end
    end
else
    for c = 2:length(contrasts)-1
        subplot(3,length(contrasts),sub2ind([length(contrasts) 3],c,3))
        hold on
        try
        a = plot(eventWindow,...
            smooth(meanTrials_dim1(c,:,iX),1),...
            'Color',colors(c,:),...
            'LineStyle',':');
        set(a,'LineWidth',2);
        maxYm(c) = max(max(meanTrials_dim1));
        minYm(c) = min(min(meanTrials_dim1));
        end
    end
end

limm = max(abs([maxYm minYm]));
for c = 2:length(contrasts)-1
    subplot(3,length(contrasts),sub2ind([length(contrasts) 3],c,3))
    line([-2 2],[0 0],'LineStyle','--','Color',[.5 .5 .5]);
    line([0 0],[-lim lim],'LineStyle','--','Color',[.5 .5 .5]);
    ylim([-lim lim]);
    xlim([-1 2])
end

end
end

%% correct for hemisphere
% 

% 
% %% plot trajectories
% tbi = 11:41;
% c = ceil(sqrt(length(expInfo)));
% r = ceil(length(expInfo)/c);
% figure(n+2);
% hold on
% 
%     subplot(r,c,iX)
%     maxXY = max(max([max(max(meanTrials_dim1))]));
%     minXY = min(min([min(min(meanTrials_dim1))]));
%     lim = max([abs(maxXY) abs(minXY)])*1.05;
%     line([-1 2],[0 0],'LineStyle','--','Color',[.5 .5 .5])
%     line([0 0],[-lim lim],'LineStyle','--','Color',[.5 .5 .5])
%     if whichETA == 1
%         xlim([-0.5 2]);
%         xlabel('time from stim on')
%     else
%         xlim([-1 1.5]);
%         xlabel('time from move on')
%     end
%     ylim([-lim lim]);
%     ylabel('cosine similarity')
%     hold on
%     box off
%         
%     for it = length(tbi)
%         for c = 1:length(whichTrials)
%             plot(eventWindow(tbi(1:it)), smooth(meanTrials_dim1(c,tbi(1:it),iX),3),'k-','Color',colors(c,:),'LineWidth',2)
%             hold on
%             box off
%             set(gca,'tickdir','out')
%         end
%     end

% %% plot epochs
% 
% epIdx{1} = 15:20;
% epIdx{2} = 22:25;
% epIdx{3} = 28:36;
% epIdx{4} = 37:41;
% 
% figure(n+2);
% set(gcf,'position', [1000 920 1120 420])
%     for e = 1:length(epIdx)
%         subplot(1,length(epIdx),e)
%             for c = 1:length(whichTrials)
%                 meanTrialsEpoch_dim1(c,iX) = nanmean(meanTrials_dim1(c,epIdx{e},iX));
%                 meanTrialsEpoch_dim2(c,iX) = nanmean(meanTrials_dim2(c,epIdx{e},iX));
%             end
%         
% 
%         maxXY(e,iX) = max([max(max(meanTrialsEpoch_dim1)) max(max(meanTrialsEpoch_dim2))]);
%         minXY(e,iX) = min([min(min(meanTrialsEpoch_dim1)) min(min(meanTrialsEpoch_dim2))]);
% 
%         for c = 1:length(whichTrials)
%             plot(meanTrialsEpoch_dim1(c,iX),meanTrialsEpoch_dim2(c,iX),'ko','MarkerSize',ms,'MarkerEdgeColor','none','MarkerFaceColor',colors(c,:),'LineStyle','none')
%             hold on
%         end
%     end
% end
% 
% for e = 1:length(epIdx)
%     lim = max([abs(max(maxXY)) abs(min(minXY))])*1.01;
%     subplot(1,length(epIdx),e)
%     h1 = line([-lim lim],[0 0],'LineStyle','--','Color',[.5 .5 .5]);
%     h2 = line([0 0],[-lim lim],'LineStyle','--','Color',[.5 .5 .5]);
%     uistack([h1 h2],'bottom');
%     xlim([-lim lim]);
%     ylim([-lim lim]);
%     axis square
%     box off
%     set(gca,'tickdir','out')
%     xlabel(strcat(dim1,{' '},'dimension'))
%     ylabel(strcat(dim2,{' '},'dimension'))
%     if e == 1
%         title('baseline')
%     elseif e == 2
%         title('stimulus')
%     elseif e == 3
%         title('movement')
%     elseif e == 4
%         title('feedback')
%     end
end
