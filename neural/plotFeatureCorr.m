contrasts = getUniqueContrasts(expInfo);
trialTypes = getTrialTypes(expInfo, behavioralData, 'late');

%% choose cells
[baselineResps, stimResps, pmovResps, movResps, rewResps] = getEpochResps(neuralData.eta);

whichCells = 'all'; %choose from 'pLabels' array
if strcmp(whichCells, 'all')
    plotCells = 1:size(neuralData.eta.alignedResps{1},3);
else
    plotCells = find(neuralData.stats.bfcH(:,strcmp(neuralData.stats.labels,whichCells)) > 0);
end

%% discard early trials
contrasts = getUniqueContrasts(expInfo);
[~, whichTrials] = selectCondition(expInfo, contrasts, behavioralData, initTrialConditions('movementTime','late'));

%whichTrials = whichTrials(whichTrials ~= 1);
%%
nt = length(expInfo.block.events.endTrialTimes);
Y_choice = behavioralData.wheelMoves.epochs(5).moveDir(whichTrials)';
Y_stimulus = expInfo.block.events.contrastValues(whichTrials)';
Y_block = expInfo.block.events.highRewardSideValues(whichTrials)';
Y_block = expInfo.block.events.responseValues(whichTrials-1)';

X = baselineResps(whichTrials,plotCells);

stim_r = zeros(1,size(X,2));
choice_r = zeros(1,size(X,2));
block_r = zeros(1,size(X,2));
stim_p = zeros(1,size(X,2));
choice_p = zeros(1,size(X,2));
block_p = zeros(1,size(X,2));

for iCell = 1:size(X,2)
    [r_stim, p_stim] = corrcoef(X(:,iCell),Y_stimulus);
    [r_choice, p_choice] = corrcoef(X(:,iCell),Y_choice);
    [r_block, p_block] = corrcoef(X(:,iCell),Y_block);
    stim_r(iCell) = min(unique(r_stim));
    stim_p(iCell) = min(unique(p_stim));
    choice_r(iCell) = min(unique(r_choice));
    choice_p(iCell) = min(unique(p_choice));
    block_r(iCell) = min(unique(r_block));
    block_p(iCell) = min(unique(p_block));
end

stim_p(stim_p > 0.05) = 0.05;
stim_p(stim_p < 0.0002) = 0.0002;
choice_p(choice_p > 0.05) = 0.05;
choice_p(choice_p < 0.0002) = 0.0002;
block_p(block_p > 0.05) = 0.05;
block_p(block_p < 0.0002) = 0.0002;

% %%
% 
% if p1 > 0.05
%     p1 = 0.05;
% end
% if p2 > 0.05
%     p2 = 0.05;
% end
% pColor = I(ceil(-log10(p1)*150-195), ceil(-log10(p2)*150-195),:)/255;
% 
% 
% %%
%%
figure;
set(gcf,'position',[276 632 1215 332]);
hold on
subplot(1,3,1)
hold on
line([0 0],[-.7 0.7],'LineStyle',':','Color',[.5 .5 .5])
line([-.7 0.7],[0 0],'LineStyle',':','Color',[.5 .5 .5])
line([-.7 0.7],[-.7 0.7],'LineStyle',':','Color',[.5 .5 .5])
line([-.7 0.7],[.7 -0.7],'LineStyle',':','Color',[.5 .5 .5])
xlim([-.5 .5]);
ylim([-.5 .5]);
axis square
box off
set(gca,'tickdir','out')

for iS = 1:length(stim_p)
    pColor(iS,:) = I(ceil(-log10(stim_p(iS))*100-130), ceil(-log10(choice_p(iS))*100-130),:)/255;
end
for iP = 1:length(stim_p)
scatter(stim_r(iP),choice_r(iP),25,'MarkerEdgeColor','none','MarkerFaceColor', pColor(iP,:),'MarkerFaceAlpha', 0.3)
end

subplot(1,3,2)
hold on
line([0 0],[-.7 0.7],'LineStyle',':','Color',[.5 .5 .5])
line([-.7 0.7],[0 0],'LineStyle',':','Color',[.5 .5 .5])
line([-.7 0.7],[-.7 0.7],'LineStyle',':','Color',[.5 .5 .5])
line([-.7 0.7],[.7 -0.7],'LineStyle',':','Color',[.5 .5 .5])
for iS = 1:length(stim_p)
    pColor(iS,:) = I(ceil(-log10(choice_p(iS))*100-130), ceil(-log10(block_p(iS))*100-130),:)/255;
end
for iP = 1:length(stim_p)
    scatter(choice_r(iP),block_r(iP),25,'MarkerEdgeColor','none','MarkerFaceColor', pColor(iP,:),'MarkerFaceAlpha', 0.3)
end
xlim([-.5 .5]);
ylim([-.5 .5]);
axis square
box off
set(gca,'tickdir','out')

subplot(1,3,3)
hold on
line([0 0],[-.7 0.7],'LineStyle',':','Color',[.5 .5 .5])
line([-.7 0.7],[0 0],'LineStyle',':','Color',[.5 .5 .5])
line([-.7 0.7],[-.7 0.7],'LineStyle',':','Color',[.5 .5 .5])
line([-.7 0.7],[.7 -0.7],'LineStyle',':','Color',[.5 .5 .5])
for iS = 1:length(stim_p)
    pColor(iS,:) = I(ceil(-log10(stim_p(iS))*100-130), ceil(-log10(block_p(iS))*100-130),:)/255;
end

for iP = 1:length(stim_p)
scatter(stim_r(iP),block_r(iP),25,'MarkerEdgeColor','none','MarkerFaceColor', pColor(iP,:),'MarkerFaceAlpha', 0.3)
end
xlim([-.5 .5]);
ylim([-.5 .5]);
axis square
box off
set(gca,'tickdir','out')
%%
figure
R=[0 1;
   .5 1];
B=[1 1
   .5 0];
G=[0 0
   .5 0];
R = interp2(R,8);
G = interp2(G,8);
B = interp2(B,8);
I = flipud((255*cat(3,R,G,B)));
image(uint8(I))


for iS = 1:length(stim_p)
    pColor(iS,:) = I(ceil(-log10(stim_p(iS))*100-130), ceil(-log10(choice_p(iS))*100-130),:)/255;
end


%%
for b = 1:14
[baselineResps, stimResps, pmovResps, movResps, rewResps] = getEpochResps(neuralData(b).eta);

Y_vel = abs(behavioralData(b).wheelMoves.epochs(5).peakVel);
Y_vel(isnan(Y_vel)) = [];
X = movResps;
X(isnan(X(:,1)),:) = [];

vel_r = zeros(1,size(X,2));
vel_p = zeros(1,size(X,2));

for iCell = 1:size(X,2)
    [r_vel, p_vel] = corrcoef(X(:,iCell),Y_vel');
    vel_r(iCell) = min(unique(r_vel));
    vel_p(iCell) = min(unique(p_vel));
end

clear r_vel
vel_pctl = prctile(zscore((Y_vel)),[10 20 30 40 50 60 70 80 90 100]);
for v = 2:9
    
    r_vel(v,1) = mean(mean(X(zscore((Y_vel))>=vel_pctl(v-1) & zscore((Y_vel))<=vel_pctl(v),vel_p<0.00001)));
end
% r_vel(1,1) = mean(mean(X(zscore(Y_vel)<=vel_pctl(1),vel_p<0.001)));
hold on
plot(vel_pctl(2:9),r_vel(2:9))
end