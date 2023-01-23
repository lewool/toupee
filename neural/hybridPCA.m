trialTypes = getTrialTypes(expInfo(1), behavioralData(1), 'late');

% pick conditions and reshape into a 1 x nc list
conditions = reshape(...
    trialTypes.intVar.all.contrast_direction,...
    [1 numel(trialTypes.intVar.all.contrast_direction)]);

%compute mean matrix (neurons x timepoints) for each condition 
for c = 1:length(conditions)
    meanMatrix(:,:,c) = squeeze(mean(neuralData.eta.alignedResps{1}(conditions{c},:,:),1))';
end

%concatenate matrices into a neurons x (timepoints x conditions)
cat_meanTrials = reshape(meanMatrix,[size(meanMatrix,1) size(meanMatrix,2)*size(meanMatrix,3)]);

%determine mean/std in order to scale/center the meanMatrix
scaleMatrix = [mean(cat_meanTrials,2) std(cat_meanTrials,[],2)];
scaled_meanTrials = (cat_meanTrials - scaleMatrix(:,1)) ./ scaleMatrix(:,2);

%reshape the trial matrix (neurons x timepoints) into a long matrix
%(neurons x (timepoints x numtrials))
for t = 1:size(neuralData.eta.alignedResps{1},1)
    trialMatrix(:,:,t) = squeeze(neuralData.eta.alignedResps{1}(t,:,:))';
end

cat_allTrials = reshape(trialMatrix,[size(trialMatrix,1) size(trialMatrix,2)*size(trialMatrix,3)]);

%scale the trial matrix by the same transformation used on the mean matrix
scaled_allTrials = (cat_allTrials - scaleMatrix(:,1)) ./ scaleMatrix(:,2);

%pca on mean matrix (take top ten PCs)
[coeff,score] = pca(scaled_meanTrials',10);

%
% projection = reshape(scaled_meanTrials'*score.M, ...
%     size(meanMatrix,2),...
%     size(meanMatrix,3),...
%     size(score.M,2));

projection = reshape(scaled_allTrials'*score.M, ...
    size(trialMatrix,2),...
    size(trialMatrix,3),...
    size(score.M,2));

%% plot
contrasts = getUniqueContrasts(expInfo);
zeroIdx = find(contrasts == 0);
walkup = length(contrasts) - zeroIdx;
walkback = zeroIdx - 1;
% allColors = [.25 0 0;.5 0 0 ;1 0 0;.8 .45 .45;.75 .75 .75;.55 .55 .55;.35 .35 .35;.15 .15 .15;0 0 0];
allColors = [0 0 .5; 0 0 1; 0 .4 1; .6 .8 1; .75 .75 .75; .8 .45 .45; 1 0 0; .5 0 0; .25 0 0];

zeroGray = find(allColors(:,1) == .75);
colors = [allColors(zeroGray-walkback:zeroGray + walkup,:); allColors(zeroGray-walkback:zeroGray + walkup,:)];
figure;
PC = 3;
for pc = [3 4 5 6 7 12 13 14 15 16]
% plotCond = trialTypes.singleVar.contrast{pc};
% plotCond = trialTypes.intVar.all.contrast_direction{pc,2};
plotCond = conditions{pc};

meanSignal = mean(projection(:,plotCond,PC),2)';
semSignal = std(projection(:,plotCond,PC),[],2)/sqrt(length(plotCond))';
plotSignal(neuralData.eta.eventWindow,meanSignal,meanSignal+semSignal,meanSignal-semSignal,colors(pc,:),'-')
end

%%
figure;
for pc = [3 4 5 6 7 12 13 14 15 16]
plotCond = conditions{pc};
hold on;
meanSignal1 = mean(projection(:,plotCond,3),2)';
meanSignal2 = mean(projection(:,plotCond,5),2)';

plot(smooth(meanSignal1),smooth(meanSignal2),'Color',colors(pc,:),'LineWidth',2)
end

%%
[baselineResps, stimResps, pmovResps, movResps, rewResps, preCueResps] = getEpochResps(neuralData.eta);

whichResps =stimResps;

%build discrimination axis 1
matrixA = movResps(trialTypes.intVar.all.contrast_direction{5,1},:)';
matrixB = movResps(trialTypes.intVar.all.contrast_direction{5,2},:)';
signal1 = mean(matrixA,2) - mean(matrixB,2);

%build discrimination axis 2
matrixC = stimResps(trialTypes.singleVar.contrast{1},:)';
matrixD = stimResps(trialTypes.singleVar.contrast{9},:)';
signal2 = mean(matrixC,2) - mean(matrixD,2);

%test trials
testA = whichResps(trialTypes.intVar.all.contrast_direction{4,1},:)';
testB = whichResps(trialTypes.intVar.all.contrast_direction{4,2},:)';


matrixE = [zscore(matrixA) zscore(matrixB) zscore(matrixC) zscore(matrixD)];
[coeff, score] = pca(matrixE,1);
e1 = matrixE*score.M;

signal = mean(matrixA,2) - mean(matrixB,2);

%rad2deg(acos(dot(e1,signal)/(norm(e1)*norm(signal))))

figure;
subplot(1,2,1)
hold on
xlim([-10 10])
ylim([-10 10])
scatter(testA'*signal1,testA'*signal2,'filled')
scatter(testB'*signal1,testB'*signal2,'filled')
prettyPlot(gca);
xlabel('signal1 dimension')
ylabel('signal2 dimension')

subplot(1,2,2)
scatter3(testA'*signal1,testA'*signal2,testA'*e1,'filled')
hold on
scatter3(testB'*signal1,testB'*signal2,testB'*e1,'filled')
xlim([-10 10])
ylim([-10 10])
zlim([-500 500])

%% SDT
close all
for t = 1:41
    testA = trialTypes.intVar.all.contrast_direction{2,1};
    testB = trialTypes.intVar.all.contrast_direction{2,2};

    popA = [projection(t,testA,1)' projection(t,testA,2)' projection(t,testA,3)' projection(t,testA,4)' projection(t,testA,5)'];
    popB = [projection(t,testB,1)' projection(t,testB,2)' projection(t,testB,3)' projection(t,testB,4)' projection(t,testB,5)'];
    
    for pc = 1:5
         [~,pv(t,pc),k(t,pc)]=kstest2(popA(:,pc),popB(:,pc));
    end
end

figure;
xlim([-.5 2])
ylim([0 1])
hold on
line([-.5 2],[0 0],'LineStyle','-','Color',[.5 .5 .5])
line([.81 .81],[-1.5 1.5],'LineStyle',':','Color',[.5 .5 .5])
line([0 0],[-1.5 1.5],'LineStyle',':','Color',[.5 .5 .5])
for pc = 1:5
    plot(neuralData.eta.eventWindow,smooth(k(:,pc)))
end
prettyPlot(gca)

%%

% close all
for t = 1:41
    testA = squeeze(neuralData.eta.alignedResps{1}(trialTypes.intVar.all.contrast_direction{3,1},t,:))';
    testB = squeeze(neuralData.eta.alignedResps{1}(trialTypes.intVar.all.contrast_direction{7,2},t,:))';

    popA = [testA'*signal1,testA'*signal2,testA'*e1];
    popB = [testB'*signal1,testB'*signal2,testB'*e1];
    
    [~,pv(t,1),k(t,1)]=kstest2(popA(:,1),popB(:,1));
    [~,pv(t,2),k(t,2)]=kstest2(popA(:,2),popB(:,2));
    [~,pv(t,3),k(t,3)]=kstest2(popA(:,3),popB(:,3));
        
    dprime(t,:) = (mean(popA) - mean(popB)) ./ sqrt(((std(popA).^2)+(std(popA).^2))*0.5);
end

figure;
xlim([-.5 2])
ylim([0 1])
hold on
line([-.5 2],[0 0],'LineStyle','-','Color',[.5 .5 .5])
line([.81 .81],[-1.5 1.5],'LineStyle',':','Color',[.5 .5 .5])
line([0 0],[-1.5 1.5],'LineStyle',':','Color',[.5 .5 .5])
plot(neuralData.eta.eventWindow,smooth(k(:,1)))
plot(neuralData.eta.eventWindow,smooth(k(:,2)))
plot(neuralData.eta.eventWindow,smooth(k(:,3)))
prettyPlot(gca)

