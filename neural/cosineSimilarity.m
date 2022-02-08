function cosineSimilarity(expInfo, behavioralData, neuralData, whichDims, pickTrials)

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

ms = 8;

%% extract trial IDs
o=1;
%organize the trials into conditions
trialTypes = getTrialTypes(expInfo(iX), behavioralData(iX), 'late');

switch pickTrials
    case 'contrast'
        whichTrials = [];
        for c = 1:length(trialTypes.intVar.all.contrast_outcome)
            whichTrials{c} = trialTypes.intVar.all.contrast_outcome{c,o};
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
            whichTrials{v} = intersect(velTrialIdx{v},trialTypes.singleVar.outcome{o});
        end
        
        velColorMap = [linspace(.75, 0, length(whichTrials))', linspace(0, .75, length(whichTrials))', linspace(.75, 0, length(whichTrials))'];
        colors = velColorMap;
        
    case 'block'
        whichTrials{1} = intersect([trialTypes.singleVar.block{1}],trialTypes.singleVar.outcome{o});
        whichTrials{2} = intersect([trialTypes.singleVar.block{2}],trialTypes.singleVar.outcome{o});
        colors = [0 0 0.5; 0.5 0 0];
        
    case 'side'
        whichTrials{1} = trialTypes.intVar.all.side_outcome{1,o};
        whichTrials{2} = trialTypes.intVar.all.side_outcome{3,o};
        colors = [0 0.4 1; 1 0 0];
        
    case 'outcome'
        whichTrials{1} = trialTypes.singleVar.outcome{1};
        whichTrials{2} = trialTypes.singleVar.outcome{2};
        colors = [0 0.5 0; 0.5 0 0];
    
    case 'low contrast'
        whichTrials = [];
        for b = 1:4
            whichTrials{1} = cat(2,trialTypes.intVar.all.contrast_outcome{3,o},trialTypes.intVar.all.contrast_outcome{4,o});
            whichTrials{2} = cat(2,trialTypes.intVar.all.contrast_outcome{6,o},trialTypes.intVar.all.contrast_outcome{7,o});
        end
        
        contrasts = getUniqueContrasts(expInfo(iX));
        colors = [.6 .8 1; .8 .45 .45];
        
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
end


%% set up mean vectors (aka 'dimensions')

%get epochs
[baselineResps, stimResps, pmovResps, movResps, rewResps] = getEpochResps(neuralData(iX).eta);

%define test and train sets (cv)
testTrials = 1:2:size(baselineResps,1);
trainTrials = 2:2:size(baselineResps,1);

%initialize mean vectors
mean_dim1 = nan(1,length(baselineResps));
mean_dim2 = nan(1,length(baselineResps));

%extract dims
dim1 = whichDims{1};
dim2 = whichDims{2};

switch dim1
    case 'stimulus'
    %extract trials that go into a stimulus dimension
    dim1a_idx = [trialTypes.singleVar.contrast{1} trialTypes.singleVar.contrast{2}];
    dim1b_idx = [trialTypes.singleVar.contrast{8} trialTypes.singleVar.contrast{9}];
    mean_dim1a = nanmean(stimResps(intersect(dim1a_idx,trainTrials),:),1);
    mean_dim1b = nanmean(stimResps(intersect(dim1b_idx,trainTrials),:),1);
    if expInfo(iX).hemisphere < 0
        mean_dim1 = mean_dim1a - mean_dim1b;
    elseif expInfo(iX).hemisphere > 0
        mean_dim1 = mean_dim1b - mean_dim1a;
    end
    
    case 'choice'
    %extract trials that go into a choice dimension
    dim1a_idx = cat(2,trialTypes.intVar.cb2D.contrast_direction{:,1});
    dim1b_idx = cat(2,trialTypes.intVar.cb2D.contrast_direction{:,2});
    mean_dim1a = nanmean(movResps(intersect(dim1a_idx,trainTrials),:),1);
    mean_dim1b = nanmean(movResps(intersect(dim1b_idx,trainTrials),:),1);
    if expInfo(iX).hemisphere < 0
        mean_dim1 = mean_dim1a - mean_dim1b;
    elseif expInfo(iX).hemisphere > 0
        mean_dim1 = mean_dim1b - mean_dim1a;
    end
    
    case 'prev choice'
    %extract trials that go into a previous choice dimension
    dim1a_idx = cat(2,trialTypes.intVar.cb2D.contrast_direction{:,1}) + 1;
    dim1b_idx = cat(2,trialTypes.intVar.cb2D.contrast_direction{:,2}) + 1;
    mean_dim1a = nanmean(movResps(intersect(dim1a_idx,trainTrials),:),1);
    mean_dim1b = nanmean(movResps(intersect(dim1b_idx,trainTrials),:),1);
    if expInfo(iX).hemisphere < 0
        mean_dim1 = mean_dim1a - mean_dim1b;
    elseif expInfo(iX).hemisphere > 0
        mean_dim1 = mean_dim1b - mean_dim1a;
    end

    case 'outcome'
    %extract trials that go into an outcome dimension
    dim1a_idx = cat(2,trialTypes.singleVar.outcome{2});
    dim1b_idx = cat(2,trialTypes.singleVar.outcome{1});
    mean_dim1a = nanmean(rewResps(intersect(dim1a_idx,trainTrials),:),1);
    mean_dim1b = nanmean(rewResps(intersect(dim1b_idx,trainTrials),:),1);
    mean_dim1 = mean_dim1b - mean_dim1a;

    case 'block'
    %extract trials that go into a block dimension
%     dim1a_idx = intersect([trialTypes.singleVar.block{1}],trialTypes.singleVar.outcome{1});
%     dim1b_idx = intersect([trialTypes.singleVar.block{2}],trialTypes.singleVar.outcome{1});
    dim1a_idx = trialTypes.intVar.cb3D.block{1};
    dim1b_idx = trialTypes.intVar.cb3D.block{2};
    mean_dim1a = nanmean(rewResps(intersect(dim1a_idx,trainTrials),:),1);
    mean_dim1b = nanmean(rewResps(intersect(dim1b_idx,trainTrials),:),1);
    if expInfo(iX).hemisphere < 0
        mean_dim1 = mean_dim1a - mean_dim1b;
    elseif expInfo(iX).hemisphere > 0
        mean_dim1 = mean_dim1b - mean_dim1a;
    end
end

switch dim2
    case 'stimulus'
    %extract trials that go into a stimulus dimension
    dim2a_idx = [trialTypes.singleVar.contrast{1} trialTypes.singleVar.contrast{2}];
    dim2b_idx = [trialTypes.singleVar.contrast{8} trialTypes.singleVar.contrast{9}];
    mean_dim2a = nanmean(stimResps(intersect(dim2a_idx,trainTrials),:),1);
    mean_dim2b = nanmean(stimResps(intersect(dim2b_idx,trainTrials),:),1);
    if expInfo(iX).hemisphere < 0
        mean_dim2 = mean_dim2a - mean_dim2b;
    elseif expInfo(iX).hemisphere > 0
        mean_dim2 = mean_dim2b - mean_dim2a;
    end

    case 'choice'
    %extract trials that go into a choice dimension
    dim2a_idx = cat(2,trialTypes.intVar.cb2D.contrast_direction{:,1});
    dim2b_idx = cat(2,trialTypes.intVar.cb2D.contrast_direction{:,2});
    mean_dim2a = nanmean(movResps(intersect(dim2a_idx,trainTrials),:),1);
    mean_dim2b = nanmean(movResps(intersect(dim2b_idx,trainTrials),:),1);
    if expInfo(iX).hemisphere < 0
        mean_dim2 = mean_dim2a - mean_dim2b;
    elseif expInfo(iX).hemisphere > 0
        mean_dim2 = mean_dim2b - mean_dim2a;
    end
    
    case 'prev choice'
    %extract trials that go into a previous choice dimension
    dim1a_idx = cat(2,trialTypes.intVar.cb2D.contrast_direction{:,1}) + 1;
    dim1b_idx = cat(2,trialTypes.intVar.cb2D.contrast_direction{:,2}) + 1;
    mean_dim2a = nanmean(movResps(intersect(dim1a_idx,trainTrials),:),1);
    mean_dim2b = nanmean(movResps(intersect(dim1b_idx,trainTrials),:),1);
    if expInfo(iX).hemisphere < 0
        mean_dim2 = mean_dim2a - mean_dim2b;
    elseif expInfo(iX).hemisphere > 0
        mean_dim2 = mean_dim2b - mean_dim2a;
    end
    
    case 'outcome'
    %extract trials that go into an outcome dimension
    dim2a_idx = cat(2,trialTypes.singleVar.outcome{2});
    dim2b_idx = cat(2,trialTypes.singleVar.outcome{1});
    mean_dim2a = nanmean(rewResps(intersect(dim2a_idx,trainTrials),:),1);
    mean_dim2b = nanmean(rewResps(intersect(dim2b_idx,trainTrials),:),1);
    mean_dim2 = mean_dim2b - mean_dim2a;

    case 'block'
    %extract trials that go into a block dimension
%     dim2a_idx = intersect([trialTypes.singleVar.block{1}],trialTypes.singleVar.outcome{1});
%     dim2b_idx = intersect([trialTypes.singleVar.block{2}],trialTypes.singleVar.outcome{1});
    dim2a_idx = trialTypes.intVar.cb3D.block{1};
    dim2b_idx = trialTypes.intVar.cb3D.block{2};
    mean_dim2a = nanmean(rewResps(intersect(dim2a_idx,trainTrials),:),1);
    mean_dim2b = nanmean(rewResps(intersect(dim2b_idx,trainTrials),:),1);
    if expInfo(iX).hemisphere < 0
        mean_dim2 = mean_dim2a - mean_dim2b;
    elseif expInfo(iX).hemisphere > 0
        mean_dim2 = mean_dim2b - mean_dim2a;
    end
end

orthHist(iX) = acosd(dot(mean_dim1, mean_dim2)/(norm(mean_dim1)*norm(mean_dim2)));
%% compute dot products


eventWindow = neuralData(iX).eta.eventWindow;
for a = 1:3
    for iT = 1:size(eventWindow,2)
        whichResps = neuralData(iX).eta.alignedResps{a}(:,iT,:);

        fullNorm = 1; % 1 if just angle, 0 if angle and relative magnitude
        for iTrial = testTrials
            if fullNorm == 1
                tNorm = norm(whichResps(iTrial,:));
            else
                tNorm = 1;
            end
            dotProd_dim1(iTrial,iT,a) = dot(whichResps(iTrial,:),mean_dim1)/(tNorm*norm(mean_dim1));
            dotProd_dim2(iTrial,iT,a) = dot(whichResps(iTrial,:),mean_dim2)/(tNorm*norm(mean_dim2));
        end
    end
end

dotProd_dim1(trainTrials,:,:) = NaN;
dotProd_dim2(trainTrials,:,:) = NaN;

%% split the trials dotProds into two groups based on trial identity

eventWindow = neuralData(iX).eta.eventWindow;

for a = 1:3
for iT = 1:length(eventWindow)
    for c = 1:length(whichTrials)
        meanTrials_dim1(c,iT,a,iX) = nanmean(dotProd_dim1(whichTrials{c},iT,a));
        meanTrials_dim2(c,iT,a,iX) = nanmean(dotProd_dim2(whichTrials{c},iT,a));
    end
end
end


%% correct for hemisphere

if expInfo(iX).hemisphere < 0
    for a = 1:3
        meanTrials_dim1(:,:,a,iX) = flipud(meanTrials_dim1(:,:,a,iX));
        meanTrials_dim2(:,:,a,iX) = flipud(meanTrials_dim2(:,:,a,iX));
    end
end

%% plot trajectories

whichETA = 1; %use stimulus-aligned
tbi = 15:35;
c = ceil(sqrt(length(expInfo)));
r = ceil(length(expInfo)/c);
figure(n+1);
hold on

    subplot(r,c,iX)
    maxXY = max(max([max(max(meanTrials_dim1(:,:,whichETA))) max(max(meanTrials_dim2(:,:,whichETA)))]));
    minXY = min(min([min(min(meanTrials_dim1(:,:,whichETA))) min(min(meanTrials_dim2(:,:,whichETA)))]));
    lim = max([abs(maxXY) abs(minXY)])*1.01;
    line([-lim lim],[0 0],'LineStyle','--','Color',[.5 .5 .5])
    line([0 0],[-lim lim],'LineStyle','--','Color',[.5 .5 .5])
    xlim([-lim lim]);
    ylim([-lim lim]);
    xlabel(dim1)
    ylabel(dim2)
    hold on
    box off
    axis square
        
    for it = length(tbi)
        for c = 1:length(whichTrials)
            plot(smooth(meanTrials_dim1(c,tbi(1:it),whichETA,iX)),smooth(meanTrials_dim2(c,tbi(1:it),whichETA,iX)),'k-','Color',colors(c,:),'LineWidth',2)
            hold on
            box off
            axis square
            set(gca,'tickdir','out')
        end
    end

%% plot epochs

%baseline, stim-aligned
epIdx{1} = 15:20;
%stim, stim-aligned
epIdx{2} = 21:26;
%move, move-aligned
epIdx{3} = 20:25;
%outcome, rew-aligned
epIdx{4} = 21:26;

figure(n+2);
set(gcf,'position', [1000 620 1120 420])
    for e = 1:length(epIdx)
        if e == 1 || e == 2
            whichETA = 1;
        elseif e == 3
            whichETA = 2;
        elseif e == 4
            whichETA = 3;
        end
        subplot(1,length(epIdx),e)
            for c = 1:length(whichTrials)
                meanTrialsEpoch_dim1(c,iX) = nanmean(meanTrials_dim1(c,epIdx{e},whichETA,iX));
                meanTrialsEpoch_dim2(c,iX) = nanmean(meanTrials_dim2(c,epIdx{e},whichETA,iX));
            end
        

        maxXY(e,iX) = max([max(max(meanTrialsEpoch_dim1)) max(max(meanTrialsEpoch_dim2))]);
        minXY(e,iX) = min([min(min(meanTrialsEpoch_dim1)) min(min(meanTrialsEpoch_dim2))]);

        for c = 1:length(whichTrials)
            plot(meanTrialsEpoch_dim1(c,iX),meanTrialsEpoch_dim2(c,iX),'ko','MarkerSize',ms,'MarkerEdgeColor','none','MarkerFaceColor',colors(c,:),'LineStyle','none')
            hold on
        end
    end
end
%%
for e = 1:length(epIdx)
    lim = max([abs(max(maxXY)) abs(min(minXY))])*1.01;
    subplot(1,length(epIdx),e)
    h1 = line([-lim lim],[0 0],'LineStyle','--','Color',[.5 .5 .5]);
    h2 = line([0 0],[-lim lim],'LineStyle','--','Color',[.5 .5 .5]);
    uistack([h1 h2],'bottom');
    xlim([-lim lim]);
    ylim([-lim lim]);
    axis square
    box off
    set(gca,'tickdir','out')
    xlabel(strcat(dim1,{' '},'dimension'))
    ylabel(strcat(dim2,{' '},'dimension'))
    if e == 1
        title('baseline')
    elseif e == 2
        title('stimulus')
    elseif e == 3
        title('movement')
    elseif e == 4
        title('feedback')
    end
end

%% histogram of dimension orthogonality
figure;
histogram(orthHist,[45:10:135]);

