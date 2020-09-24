clear all
expInfo = initExpInfo({{'LEW031'}},{{'2020-02-03',1,[1]}});
[expInfo, neuralData, behavioralData] = processExperiment(expInfo);
alignedResps = neuralData.eta.alignedResps;
pLabels = neuralData.stats.labels;
bfcH = neuralData.stats.bfcH;
trialTypes = getTrialTypes(expInfo, behavioralData, 'late');
[baselineResps, stimResps, pmovResps, movResps, rewResps] = getEpochResps(neuralData.eta);
eventWindow = neuralData.eta.eventWindow;

%% choose cells

whichCells = 'advanceDirection';
% whichCells = 'leftStim'; %choose from 'pLabels' array
if strcmp(whichCells, 'all')
    plotCells = 1:size(alignedResps{1},3);
elseif strcmp(whichCells, 'advanceDirection')
    plotCells = find(bfcH(:,strcmp(pLabels,'advanceMov')) > 0 & bfcH(:,strcmp(pLabels,'rightMov')) > 0);
else
    plotCells = find(bfcH(:,strcmp(pLabels,whichCells)) > 0);
end

ETA = 1;
color = [0 0 0];

%%
clear shuffTrials mr2_shuff mr1_shuff mr4_shuff mr3_shuff mr6_shuff mr5_shuff
trials{1} = trialTypes.intVar.cb3D.side{1};
trials{2} = trialTypes.intVar.cb3D.side{3};
trials{3} = trialTypes.intVar.cb3D.direction{1};
trials{4} = trialTypes.intVar.cb3D.direction{2};
trials{5} = trialTypes.intVar.cb3D.block{1};
trials{6} = trialTypes.intVar.cb3D.block{2};

resp = stimResps(:,plotCells);% - baselineResps(:,plotCells);
mr1 = mean(resp(trials{1},:));
mr2 = mean(resp(trials{2},:));
mr3 = mean(resp(trials{3},:));
mr4 = mean(resp(trials{4},:));
mr5 = mean(resp(trials{5},:));
mr6 = mean(resp(trials{6},:));

validTrials = sort(cat(2,trialTypes.singleVar.contrast{:}));
numTrials = size(baselineResps,1);
scatter_text = 'stimulus epoch';

for s = 1:length(plotCells)
    
    shuffTrials{1}(:,s) = randsample(validTrials,length(trials{1}));
    shuffTrials{2}(:,s) = randsample(validTrials,length(trials{2}));
    shuffTrials{3}(:,s) = randsample(validTrials,length(trials{3}));
    shuffTrials{4}(:,s) = randsample(validTrials,length(trials{4}));
    shuffTrials{5}(:,s) = randsample(validTrials,length(trials{5}));
    shuffTrials{6}(:,s) = randsample(validTrials,length(trials{6}));
end

for c = 1:length(plotCells)
    mr1_shuff(c) = mean(resp(shuffTrials{1}(:,c),c));
    mr2_shuff(c) = mean(resp(shuffTrials{2}(:,c),c));
    mr3_shuff(c) = mean(resp(shuffTrials{3}(:,c),c));
    mr4_shuff(c) = mean(resp(shuffTrials{4}(:,c),c));
    mr5_shuff(c) = mean(resp(shuffTrials{5}(:,c),c));
    mr6_shuff(c) = mean(resp(shuffTrials{6}(:,c),c));
end

figure();
set(gcf,'position',[68         337        1428         644]);
subplot(4,3,[1 4 7])
title('Stimulus (right versus left)')
axLim = max([max(mr1) max(mr2)]);
text(axLim*.9,axLim*.1,scatter_text,'FontAngle','italic','HorizontalAlignment','right');
% figure;hold on
line([0 axLim],[0 axLim],'LineStyle','-','Color','k');
hold on;
plsc = scatter(mr2,mr1,35);
set(plsc,'MarkerFaceColor',[1 .5 0],'MarkerEdgeColor','none','MarkerFaceAlpha', 0.5)
axis square
ylim([0 axLim])
xlim([0 axLim])
xlabel('Mean activity (right)')
ylabel('Mean activity (left)')
% xticks([0 axLim*.5 axLim])
% yticks([0 axLim*.5 axLim])
set(gcf,'renderer','Painters');
ax = gca;
ax.TickDir = 'out';

[a,~] = signrank(mr2-mr1,mr2_shuff-mr1_shuff,'tail','both');
disp(a);
subplot(4,3,[10])

[f_shuffle, x_shuffle] = ksdensity(mr2_shuff-mr1_shuff);
area(x_shuffle,f_shuffle,'FaceColor',[0.75 0.75 0.75],'FaceAlpha',.7,'EdgeColor','none')
% plot_fshuffle = plot(x_shuffle,f_shuffle,'LineWidth',2,'Color',[0.75 0.75 0.75]);
hold on;
[f_data, x_data] = ksdensity(mr2-mr1);
% area(x_data,f_data,'FaceColor',[0 0 0],'FaceAlpha',.4,'EdgeColor','none')
axlim_x = max([abs(x_data) abs(x_shuffle)]);
axlim_y = max([abs(f_data) abs(f_shuffle)]);
plot_fdata = plot(x_data,f_data,'LineWidth',2,'Color',[1 .5 0]);
line([0 0],[0 max([f_shuffle f_data])],'LineStyle','-','Color','k')
ax = gca;
ax.YAxis.Visible = 'off';
ax.TickDir = 'out';
% set(gca,'xtick',[0 0])
% xlabel('Right – Left')
xlim([-axlim_x axlim_x])
ylim([0 max([f_shuffle f_data])])
% xticklabels({''})
% axis square
% view([45 90])
set(gca, 'Color', 'None')
xlabel('Mean right – mean left')
box off
text(.9*axlim_x,.25*axlim_y,strcat({'p = '},num2str(a)),'HorizontalAlignment','right')

subplot(4,3,[2 5 8])
title('Direction (right versus left)')
axLim = max([max(mr3) max(mr4)]);
text(axLim*.9,axLim*.1,scatter_text,'FontAngle','italic','HorizontalAlignment','right');
% figure;hold on
line([0 axLim],[0 axLim],'LineStyle','-','Color','k');
hold on;
plsc = scatter(mr4,mr3,35);
set(plsc,'MarkerFaceColor',[1 0 1],'MarkerEdgeColor','none','MarkerFaceAlpha', 0.5)
axis square
ylim([0 axLim])
xlim([0 axLim])
xlabel('Mean activity (right)')
% xticks([0 axLim*.5 axLim])
% yticks([0 axLim*.5 axLim])
set(gcf,'renderer','Painters');
ax = gca;
ax.TickDir = 'out';

[a,~] = signrank(mr4-mr3,mr4_shuff-mr3_shuff,'tail','both');
disp(a);
subplot(4,3,[11])

[f_shuffle, x_shuffle] = ksdensity(mr4_shuff-mr3_shuff);
area(x_shuffle,f_shuffle,'FaceColor',[0.75 0.75 0.75],'FaceAlpha',.7,'EdgeColor','none')
% plot_fshuffle = plot(x_shuffle,f_shuffle,'LineWidth',2,'Color',[0.75 0.75 0.75]);
hold on;
[f_data, x_data] = ksdensity(mr4-mr3);
% area(x_data,f_data,'FaceColor',[0 0 0],'FaceAlpha',.4,'EdgeColor','none')
axlim_x = max([abs(x_data) abs(x_shuffle)]);
axlim_y = max([abs(f_data) abs(f_shuffle)]);
plot_fdata = plot(x_data,f_data,'LineWidth',2,'Color',[1 0 1]);
line([0 0],[0 max([f_shuffle f_data])],'LineStyle','-','Color','k')
ax = gca;
ax.YAxis.Visible = 'off';
ax.TickDir = 'out';
% set(gca,'xtick',[0 0])
% xlabel('Right – Left')
xlim([-axlim_x axlim_x])
ylim([0 max([f_shuffle f_data])])
% xticklabels({''})
% axis square
% view([45 90])
set(gca, 'Color', 'None')
xlabel('Mean right – mean left')
box off
text(.9*axlim_x,.25*axlim_y,strcat({'p = '},num2str(a)),'HorizontalAlignment','right')

subplot(4,3,[3 6 9])
title('Block (right versus left)')
axLim = max([max(mr5) max(mr6)]);
text(axLim*.9,axLim*.1,scatter_text,'FontAngle','italic','HorizontalAlignment','right');
% figure;hold on
line([0 axLim],[0 axLim],'LineStyle','-','Color','k');
hold on;
plsc = scatter(mr6,mr5,35);
set(plsc,'MarkerFaceColor',[0 .8 .7],'MarkerEdgeColor','none','MarkerFaceAlpha', 0.5)
axis square
ylim([0 axLim])
xlim([0 axLim])
xlabel('Mean activity (right)')
% xticks([0 axLim*.5 axLim])
% yticks([0 axLim*.5 axLim])
set(gcf,'renderer','Painters');
ax = gca;
ax.TickDir = 'out';

[a,~] = signrank(mr6-mr5,mr6_shuff-mr5_shuff,'tail','both');
disp(a);
subplot(4,3,[12])

[f_shuffle, x_shuffle] = ksdensity(mr6_shuff-mr5_shuff);
area(x_shuffle,f_shuffle,'FaceColor',[0.75 0.75 0.75],'FaceAlpha',.7,'EdgeColor','none')
% plot_fshuffle = plot(x_shuffle,f_shuffle,'LineWidth',2,'Color',[0.75 0.75 0.75]);
hold on;
[f_data, x_data] = ksdensity(mr6-mr5);
% area(x_data,f_data,'FaceColor',[0 0 0],'FaceAlpha',.4,'EdgeColor','none')
axlim_x = max([abs(x_data) abs(x_shuffle)]);
axlim_y = max([abs(f_data) abs(f_shuffle)]);
plot_fdata = plot(x_data,f_data,'LineWidth',2,'Color',[0 .8 .7]);
line([0 0],[0 max([f_shuffle f_data])],'LineStyle','-','Color','k')
ax = gca;
ax.YAxis.Visible = 'off';
ax.TickDir = 'out';
% set(gca,'xtick',[0 0])
% xlabel('Right – Left')
xlim([-axlim_x axlim_x])
ylim([0 max([f_shuffle f_data])])
% xticklabels({''})
% axis square
% view([45 90])
set(gca, 'Color', 'None')
xlabel('Mean right – mean left')
box off
text(.9*axlim_x,.25*axlim_y,strcat({'p = '},num2str(a)),'HorizontalAlignment','right')

%%
clear shuffTrials mr1_shuff_entireTrial mr2_shuff_entireTrial

trials{1} = trialTypes.intVar.cb3D.side{1};
trials{2} = trialTypes.intVar.cb3D.side{3};
trials{3} = trialTypes.intVar.cb3D.direction{1};
trials{4} = trialTypes.intVar.cb3D.direction{2};
trials{5} = trialTypes.intVar.cb3D.block{1};
trials{6} = trialTypes.intVar.cb3D.block{2};

for s = 1:length(plotCells)
    
    shuffTrials{1}(:,s) = randsample(validTrials,length(trials{1}));
    shuffTrials{2}(:,s) = randsample(validTrials,length(trials{2}));
    shuffTrials{3}(:,s) = randsample(validTrials,length(trials{3}));
    shuffTrials{4}(:,s) = randsample(validTrials,length(trials{4}));
    shuffTrials{5}(:,s) = randsample(validTrials,length(trials{5}));
    shuffTrials{6}(:,s) = randsample(validTrials,length(trials{6}));
end

mr1_entireTrial = squeeze(mean(alignedResps{ETA}(trials{1},:,plotCells)));
mr2_entireTrial = squeeze(mean(alignedResps{ETA}(trials{2},:,plotCells)));
mr3_entireTrial = squeeze(mean(alignedResps{ETA}(trials{3},:,plotCells)));
mr4_entireTrial = squeeze(mean(alignedResps{ETA}(trials{4},:,plotCells)));
mr5_entireTrial = squeeze(mean(alignedResps{ETA}(trials{5},:,plotCells)));
mr6_entireTrial = squeeze(mean(alignedResps{ETA}(trials{6},:,plotCells)));

for c = 1:length(plotCells)
    mr1_shuff_entireTrial(:,c) = mean(alignedResps{ETA}(shuffTrials{1}(:,c),:,c));
    mr2_shuff_entireTrial(:,c) = mean(alignedResps{ETA}(shuffTrials{2}(:,c),:,c));
    mr3_shuff_entireTrial(:,c) = mean(alignedResps{ETA}(shuffTrials{3}(:,c),:,c));
    mr4_shuff_entireTrial(:,c) = mean(alignedResps{ETA}(shuffTrials{4}(:,c),:,c));
    mr5_shuff_entireTrial(:,c) = mean(alignedResps{ETA}(shuffTrials{5}(:,c),:,c));
    mr6_shuff_entireTrial(:,c) = mean(alignedResps{ETA}(shuffTrials{6}(:,c),:,c));
end


for ti = 1:length(eventWindow)
    mu_s1(ti) = mean(mr2_entireTrial(ti,:)-mr1_entireTrial(ti,:));
    mu_n1(ti) = mean(mr2_shuff_entireTrial(ti,:)-mr1_shuff_entireTrial(ti,:));
    sigma_s1(ti) = std(mr2_entireTrial(ti,:)-mr1_entireTrial(ti,:));
    sigma_n1(ti) = std(mr2_shuff_entireTrial(ti,:)-mr1_shuff_entireTrial(ti,:));

    dprime1(ti) = (mu_s1(ti) - mu_n1(ti))/(sqrt((.5*(sigma_s1(ti)^2 + sigma_n1(ti)^2))));
    
    mu_s2(ti) = mean(mr4_entireTrial(ti,:)-mr3_entireTrial(ti,:));
    mu_n2(ti) = mean(mr4_shuff_entireTrial(ti,:)-mr3_shuff_entireTrial(ti,:));
    sigma_s2(ti) = std(mr4_entireTrial(ti,:)-mr3_entireTrial(ti,:));
    sigma_n2(ti) = std(mr4_shuff_entireTrial(ti,:)-mr3_shuff_entireTrial(ti,:));

    dprime2(ti) = (mu_s2(ti) - mu_n2(ti))/(sqrt((.5*(sigma_s2(ti)^2 + sigma_n2(ti)^2))));

    mu_s3(ti) = mean(mr6_entireTrial(ti,:)-mr5_entireTrial(ti,:));
    mu_n3(ti) = mean(mr6_shuff_entireTrial(ti,:)-mr5_shuff_entireTrial(ti,:));
    sigma_s3(ti) = std(mr6_entireTrial(ti,:)-mr5_entireTrial(ti,:));
    sigma_n3(ti) = std(mr6_shuff_entireTrial(ti,:)-mr5_shuff_entireTrial(ti,:));

    dprime3(ti) = (mu_s3(ti) - mu_n3(ti))/(sqrt((.5*(sigma_s3(ti)^2 + sigma_n3(ti)^2))));
end

% for ti = 1:length(eventWindow)
%     mu_s(ti) = mean(mr2_entireTrial(ti,:));
%     mu_n(ti) = mean(mr1_entireTrial(ti,:));
%     sigma_s(ti) = std(mr2_entireTrial(ti,:));
%     sigma_n(ti) = std(mr1_entireTrial(ti,:));
% 
%     dprime(ti) = (mu_s(ti) - mu_n(ti))/(sqrt((.5*(sigma_s(ti)^2 + sigma_n(ti)^2))));
% end

figure(2);
set(gcf,'position',[56   192   674   735])
hold on;
subplot(3,2,1)
line([0 0],[-.05 .05],'LineStyle',':','Color',[.5 .5 .5])
hold on;
line([-2 2],[0 0],'LineStyle',':','Color',[.5 .5 .5])
pdiff = plot(eventWindow, smooth(mean(mr2_entireTrial,2)) - smooth(mean(mr1_entireTrial,2)),'Color',[1 .5 0],'LineWidth',2);

pdiff.Color(4) = 1;
box off
xlim([-1 2])
ylim([-.02 0.02])
ylabel("Norm. activity")
title('R – L difference')
box off
ax = gca;
ax.TickDir = 'out';
text(1.8,-.015,'stimulus','HorizontalAlignment','right')

subplot(3,2,2)
line([0 0],[-1 1],'LineStyle',':','Color',[.5 .5 .5])
hold on;
line([-2 2],[0 0],'LineStyle',':','Color',[.5 .5 .5])
pdisc = plot(eventWindow,smooth(dprime1,5),'Color',[1 .5 0],'LineWidth',2);
pdisc.Color(4) = 1;
ylim([-1 1])
xlim([-1 2])
ylabel("d'")
title('R – L discriminability')
box off
ax = gca;
ax.TickDir = 'out';
text(1.8,-.75,'stimulus','HorizontalAlignment','right')

subplot(3,2,3)
line([0 0],[-.05 .05],'LineStyle',':','Color',[.5 .5 .5])
hold on;
line([-2 2],[0 0],'LineStyle',':','Color',[.5 .5 .5])
pdiff = plot(eventWindow, smooth(mean(mr4_entireTrial,2)) - smooth(mean(mr3_entireTrial,2)),'Color',[1 0 1],'LineWidth',2);
pdiff.Color(4) = 1;
box off
xlim([-1 2])
ylim([-.02 0.02])
ylabel("Norm. activity")
box off
ax = gca;
ax.TickDir = 'out';
text(1.8,-.015,'direction','HorizontalAlignment','right')

subplot(3,2,4)
line([0 0],[-1 1],'LineStyle',':','Color',[.5 .5 .5])
hold on;
line([-2 2],[0 0],'LineStyle',':','Color',[.5 .5 .5])
pdisc = plot(eventWindow,smooth(dprime2,5),'Color',[1 0 1],'LineWidth',2);
pdisc.Color(4) = 1;
ylim([-1 1])
xlim([-1 2])
ylabel("d'")
box off
ax = gca;
ax.TickDir = 'out';
text(1.8,-.75,'direction','HorizontalAlignment','right')

subplot(3,2,5)
line([0 0],[-.05 .05],'LineStyle',':','Color',[.5 .5 .5])
hold on;
line([-2 2],[0 0],'LineStyle',':','Color',[.5 .5 .5])
pdiff = plot(eventWindow, smooth(mean(mr6_entireTrial,2)) - smooth(mean(mr5_entireTrial,2)),'Color',[0 .8 .7],'LineWidth',2);
pdiff.Color(4) = 1;
box off
xlim([-1 2])
ylim([-.02 0.02])
xlabel('Time from stimulus onset (s)')
ylabel("Norm. activity")
box off
ax = gca;
ax.TickDir = 'out';
text(1.8,-.015,'block','HorizontalAlignment','right')

subplot(3,2,6)
line([0 0],[-1 1],'LineStyle',':','Color',[.5 .5 .5])
hold on;
line([-2 2],[0 0],'LineStyle',':','Color',[.5 .5 .5])
pdisc = plot(eventWindow,smooth(dprime3,5),'Color',[0 .8 .7],'LineWidth',2);
pdisc.Color(4) = 1;
ylim([-1 1])
xlim([-1 2])
xlabel('Time from stimulus onset (s)')
ylabel("d'")
box off
ax = gca;
ax.TickDir = 'out';
text(1.8,-.75,'block','HorizontalAlignment','right')

%%
if ETA == 1
    range = 21:26;
elseif ETA == 2
    range = 11:16;
end
clear mr1_shuff_entireTrial mr2_shuff_entireTrial mr3_shuff_entireTrial mr4_shuff_entireTrial shuffTrials

trials{1} = trialTypes.intVar.cb3D.side{1};
trials{2} = trialTypes.intVar.cb3D.side{3};
trials{3} = trialTypes.intVar.cb3D.direction{1};
trials{4} = trialTypes.intVar.cb3D.direction{2};
trials{5} = trialTypes.intVar.cb3D.block{1};
trials{6} = trialTypes.intVar.cb3D.block{2};

for s = 1:length(plotCells)
    
shuffTrials{1}(:,s) = randsample(validTrials,length(trials{1}));
shuffTrials{2}(:,s) = randsample(validTrials,length(trials{2}));
shuffTrials{3}(:,s) = randsample(validTrials,length(trials{3}));
shuffTrials{4}(:,s) = randsample(validTrials,length(trials{4}));
shuffTrials{5}(:,s) = randsample(validTrials,length(trials{5}));
shuffTrials{6}(:,s) = randsample(validTrials,length(trials{6}));
end

mr1_entireTrial = squeeze(mean(alignedResps{ETA}(trials{1},:,plotCells)));
mr2_entireTrial = squeeze(mean(alignedResps{ETA}(trials{2},:,plotCells)));
mr3_entireTrial = squeeze(mean(alignedResps{ETA}(trials{3},:,plotCells)));
mr4_entireTrial = squeeze(mean(alignedResps{ETA}(trials{4},:,plotCells)));
mr5_entireTrial = squeeze(mean(alignedResps{ETA}(trials{5},:,plotCells)));
mr6_entireTrial = squeeze(mean(alignedResps{ETA}(trials{6},:,plotCells)));

for c = 1:length(plotCells)
    mr1_shuff_entireTrial(:,c) = mean(alignedResps{ETA}(shuffTrials{1}(:,c),:,c));
    mr2_shuff_entireTrial(:,c) = mean(alignedResps{ETA}(shuffTrials{2}(:,c),:,c));
    mr3_shuff_entireTrial(:,c) = mean(alignedResps{ETA}(shuffTrials{3}(:,c),:,c));
    mr4_shuff_entireTrial(:,c) = mean(alignedResps{ETA}(shuffTrials{4}(:,c),:,c));
    mr5_shuff_entireTrial(:,c) = mean(alignedResps{ETA}(shuffTrials{5}(:,c),:,c));
    mr6_shuff_entireTrial(:,c) = mean(alignedResps{ETA}(shuffTrials{6}(:,c),:,c));
end


for ti = 1:length(eventWindow)
    mu_s1(ti) = mean(mr2_entireTrial(ti,:)-mr1_entireTrial(ti,:));
    mu_n1(ti) = mean(mr2_shuff_entireTrial(ti,:)-mr1_shuff_entireTrial(ti,:));
    sigma_s1(ti) = std(mr2_entireTrial(ti,:)-mr1_entireTrial(ti,:));
    sigma_n1(ti) = std(mr2_shuff_entireTrial(ti,:)-mr1_shuff_entireTrial(ti,:));

    dprime1(ti) = (mu_s1(ti) - mu_n1(ti))/(sqrt((.5*(sigma_s1(ti)^2 + sigma_n1(ti)^2))));
    
    mu_s2(ti) = mean(mr4_entireTrial(ti,:)-mr3_entireTrial(ti,:));
    mu_n2(ti) = mean(mr4_shuff_entireTrial(ti,:)-mr3_shuff_entireTrial(ti,:));
    sigma_s2(ti) = std(mr4_entireTrial(ti,:)-mr3_entireTrial(ti,:));
    sigma_n2(ti) = std(mr4_shuff_entireTrial(ti,:)-mr3_shuff_entireTrial(ti,:));

    dprime2(ti) = (mu_s2(ti) - mu_n2(ti))/(sqrt((.5*(sigma_s2(ti)^2 + sigma_n2(ti)^2))));

    mu_s3(ti) = mean(mr6_entireTrial(ti,:)-mr5_entireTrial(ti,:));
    mu_n3(ti) = mean(mr6_shuff_entireTrial(ti,:)-mr5_shuff_entireTrial(ti,:));
    sigma_s3(ti) = std(mr6_entireTrial(ti,:)-mr5_entireTrial(ti,:));
    sigma_n3(ti) = std(mr6_shuff_entireTrial(ti,:)-mr5_shuff_entireTrial(ti,:));

    dprime3(ti) = (mu_s3(ti) - mu_n3(ti))/(sqrt((.5*(sigma_s3(ti)^2 + sigma_n3(ti)^2))));
end

% for ti = 1:length(eventWindow)
%     mu_s1(ti) = mean(mr2_entireTrial(ti,:));
%     mu_n1(ti) = mean(mr1_entireTrial(ti,:));
%     sigma_s1(ti) = std(mr2_entireTrial(ti,:));
%     sigma_n1(ti) = std(mr1_entireTrial(ti,:));
% 
%     dprime1(ti) = (mu_s1(ti) - mu_n1(ti))/(sqrt((.5*(sigma_s1(ti)^2 + sigma_n1(ti)^2))));
%     
%     mu_s2(ti) = mean(mr4_entireTrial(ti,:));
%     mu_n2(ti) = mean(mr3_entireTrial(ti,:));
%     sigma_s2(ti) = std(mr4_entireTrial(ti,:));
%     sigma_n2(ti) = std(mr3_entireTrial(ti,:));
% 
%     dprime2(ti) = (mu_s2(ti) - mu_n2(ti))/(sqrt((.5*(sigma_s2(ti)^2 + sigma_n2(ti)^2))));
%     
%     mu_s3(ti) = mean(mr5_entireTrial(ti,:));
%     mu_n3(ti) = mean(mr6_entireTrial(ti,:));
%     sigma_s3(ti) = std(mr5_entireTrial(ti,:));
%     sigma_n3(ti) = std(mr6_entireTrial(ti,:));
% 
%     dprime3(ti) = (mu_s3(ti) - mu_n3(ti))/(sqrt((.5*(sigma_s3(ti)^2 + sigma_n3(ti)^2))));
% end

figure(3);
set(gcf,'position',[224   184   675   700]);
subplot(2,2,1)
hold on;
line([0 0],[-1 1],'LineStyle',':','Color',[.5 .5 .5])
line([-1 1],[0 0],'LineStyle',':','Color',[.5 .5 .5])
line([-1 -1],[-1 1],'LineStyle','-','LineWidth',2,'Color',[1 0 1])
line([-1 1],[-1 -1],'LineStyle','-','LineWidth',2,'Color',[1 .5 0])
axis square
ptraj = plot(smooth(dprime1(11:end)),smooth(dprime2(11:end)),'Color',color,'LineWidth',1.5);
% ptraj.Color(4) = .5;
alpha(.3)
xlim([-1 1])
ylim([-1 1])
yticks([ -1 -0.5 0 .5 1])
ylabel("d' (move R – move L)")
xlabel("d' (stim R – stim L)")
title('Full trial (-500–2000ms)')
ax = gca;
ax.TickDir = 'out';



subplot(2,2,2);
hold on;
line([0 0],[-1 1],'LineStyle',':','Color',[.5 .5 .5])
line([-1 1],[0 0],'LineStyle',':','Color',[.5 .5 .5])
line([-1 -1],[-1 1],'LineStyle','-','LineWidth',2,'Color',[1 0 1])
line([-1 1],[-1 -1],'LineStyle','-','LineWidth',2,'Color',[1 .5 0])
axis square
scatter(mean(dprime1(range)),mean(dprime2(range)),'o','MarkerFaceColor',color,'MarkerEdgeColor','none','MarkerFaceAlpha',1)
xlim([-1 1])
ylim([-1 1])
yticks([ -1 -0.5 0 .5 1])
xlabel("d' (stim R – stim L)")
title('Stim ON (0–300ms)')
ax = gca;
ax.TickDir = 'out';

subplot(2,2,3)
hold on;
line([0 0],[-1 1],'LineStyle',':','Color',[.5 .5 .5])
line([-1 1],[0 0],'LineStyle',':','Color',[.5 .5 .5])
line([-1 -1],[-1 1],'LineStyle','-','LineWidth',2,'Color',[0 .8 .7])
line([-1 1],[-1 -1],'LineStyle','-','LineWidth',2,'Color',[1 .5 0])
axis square
ptraj = plot(smooth(dprime1(11:end)),smooth(dprime3(11:end)),'Color',color,'LineWidth',1.5);
% ptraj.Color(4) = .5;
alpha(.3)
xlim([-1 1])
ylim([-1 1])
yticks([ -1 -0.5 0 .5 1])
ylabel("d' (block R – block L)")
xlabel("d' (stim R – stim L)")
ax = gca;
ax.TickDir = 'out';


subplot(2,2,4);
hold on;
line([0 0],[-1 1],'LineStyle',':','Color',[.5 .5 .5])
line([-1 1],[0 0],'LineStyle',':','Color',[.5 .5 .5])
line([-1 -1],[-1 1],'LineStyle','-','LineWidth',2,'Color',[0 .8 .7])
line([-1 1],[-1 -1],'LineStyle','-','LineWidth',2,'Color',[1 .5 0])
axis square
scatter(mean(dprime1(range)),mean(dprime3(range)),'o','MarkerFaceColor',color,'MarkerEdgeColor','none','MarkerFaceAlpha',1)
xlim([-1 1])
ylim([-1 1])
yticks([ -1 -0.5 0 .5 1])
xlabel("d' (stim R – stim L)")
ax = gca;
ax.TickDir = 'out';



%%
if ETA == 1
    range = 21:26;
elseif ETA == 2
    range = 11:16;
end
clear mr1_shuff_entireTrial mr2_shuff_entireTrial mr3_shuff_entireTrial mr4_shuff_entireTrial shuffTrials
trajLim = [-.02 .02];
meanLim = [-.01 .01];

trials{1} = trialTypes.intVar.cb3D.side{1};
trials{2} = trialTypes.intVar.cb3D.side{3};
trials{3} = trialTypes.intVar.cb3D.direction{1};
trials{4} = trialTypes.intVar.cb3D.direction{2};
trials{5} = trialTypes.intVar.cb3D.block{1};
trials{6} = trialTypes.intVar.cb3D.block{2};

mr1_entireTrial = squeeze(nanmean(alignedResps{ETA}(trials{1},:,plotCells)));
mr2_entireTrial = squeeze(nanmean(alignedResps{ETA}(trials{2},:,plotCells)));
mr3_entireTrial = squeeze(nanmean(alignedResps{ETA}(trials{3},:,plotCells)));
mr4_entireTrial = squeeze(nanmean(alignedResps{ETA}(trials{4},:,plotCells)));
mr5_entireTrial = squeeze(nanmean(alignedResps{ETA}(trials{5},:,plotCells)));
mr6_entireTrial = squeeze(nanmean(alignedResps{ETA}(trials{6},:,plotCells)));

mdiff1 = nanmean(mr2_entireTrial,2) - nanmean(mr1_entireTrial,2);
mdiff2 = nanmean(mr4_entireTrial,2) - nanmean(mr3_entireTrial,2);
mdiff3 = nanmean(mr6_entireTrial,2) - nanmean(mr5_entireTrial,2);

figure();
set(gcf,'position',[224   184   675   700]);
subplot(2,2,1)
hold on;
line([0 0],[-.02 .02],'LineStyle',':','Color',[.5 .5 .5])
line([-.02 .02],[0 0],'LineStyle',':','Color',[.5 .5 .5])
axis square
ptraj = plot(smooth(mdiff1(11:end)),smooth(mdiff2(11:end)),'Color',color,'LineWidth',1.5);
ptraj.Color(4) = .5;
alpha(.3)
xlim(trajLim)
ylim(trajLim)
% yticks([ -1 -0.5 0 .5 1])
ylabel("d' (move R – move L)")
xlabel("d' (stim R – stim L)")
title('Full trial (-500–2000ms)')


subplot(2,2,2);
hold on;
line([0 0],[-.02 .02],'LineStyle',':','Color',[.5 .5 .5])
line([-.02 .02],[0 0],'LineStyle',':','Color',[.5 .5 .5])
axis square
scatter(mean(mdiff1(range)),mean(mdiff2(range)),'o','MarkerFaceColor',color,'MarkerEdgeColor','none','MarkerFaceAlpha',.5)
xlim(meanLim)
ylim(meanLim)
% yticks([ -1 -0.5 0 .5 1])
xlabel("d' (stim R – stim L)")
title('Stim ON (0–300ms)')

subplot(2,2,3)
hold on;
line([0 0],[-.02 .02],'LineStyle',':','Color',[.5 .5 .5])
line([-.02 .02],[0 0],'LineStyle',':','Color',[.5 .5 .5])
axis square
ptraj = plot(smooth(mdiff1(11:end)),smooth(mdiff3(11:end)),'Color',color,'LineWidth',1.5);
ptraj.Color(4) = .5;
alpha(.3)
xlim(trajLim)
ylim(trajLim)
% yticks([ -1 -0.5 0 .5 1])
ylabel("d' (block R – block L)")
xlabel("d' (stim R – stim L)")
title('Full trial (-500–2000ms)')


subplot(2,2,4);
hold on;
line([0 0],[-.02 .02],'LineStyle',':','Color',[.5 .5 .5])
line([-.02 .02],[0 0],'LineStyle',':','Color',[.5 .5 .5])
axis square
scatter(mean(mdiff1(range)),mean(mdiff3(range)),'o','MarkerFaceColor',color,'MarkerEdgeColor','none','MarkerFaceAlpha',.5)
xlim(meanLim)
ylim(meanLim)
% yticks([ -1 -0.5 0 .5 1])
xlabel("d' (stim R – stim L)")
title('Stim ON (0–300ms)')

