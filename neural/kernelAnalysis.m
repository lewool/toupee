function kernelAnalysis
%%
for m = 7:length(mouseList)
    
    %load data
    mouseName = char(mouseList{m});
    expDate = char(expList{m}{1});
    expNum = expList{m}{2};
    [expRef, ~] = data.constructExpRef(mouseName,expDate,expNum);
    fprintf('Loading %s...\n',expRef)
    [behavioralData, expInfo, neuralData] = data.loadDataset(mouseName, expDate, expNum);
    expInfo.hemisphere = hemList(m);
    
    %select features
    features = {'stimulus' 'value' 'block' 'action' 'choice' 'outcome'};

    %fit model
    [kernelFunctions, predictedActivity, optR, maxEV, cellFeatureStrength] = ...
        kernelFitting(behavioralData, expInfo, neuralData, features);

    
    %extract kernels
    [predictors, windows] = getPredictors(...
        expInfo, ...
        behavioralData, ...
        features, ...
        0.1);
    
    [featureList, fitKernels] = unpackKernels(kernelFunctions, predictors, windows);

    %save

    kernelAnalysis.features = featureList;
    kernelAnalysis.predictors = predictors;
    kernelAnalysis.windows = windows;
    kernelAnalysis.fitKernels = fitKernels;
    kernelAnalysis.kernelFunctions = kernelFunctions;
    kernelAnalysis.predictedActivity = predictedActivity;
    kernelAnalysis.optR = optR;
    kernelAnalysis.maxEV = maxEV;
    kernelAnalysis.cellFeatureStrength = cellFeatureStrength;
    
    [kernelAnalysis] = alignKernels(expInfo, kernelAnalysis, neuralData, behavioralData);
        
    sessDir = fullfile('G:\Workspaces',mouseName,expDate,num2str(expNum));
    cd(sessDir)
    save('kernelAnalysis_new.mat','kernelAnalysis', '-v7.3')

    clearvars -except mouseList expList hemList fovList
end

%% plot kernels
clear colors
cellno = 669;
Fs = 0.1;

for f = 1:length(kernelAnalysis.features)
    colors(f,:) = kernelAnalysis.predictors.(matlab.lang.makeValidName(kernelAnalysis.features{f})).color;
end

limit = max(abs(kernelAnalysis.kernelFunctions(:,cellno)))*1.1;
figure;
set(gcf,'position',[1165 950 1220 270])

subplot(1,4,1)
hold on;
ps = find(contains(kernelAnalysis.features,'stimulus'));
for p = 1:length(ps)
    plot(kernelAnalysis.windows.value*Fs,kernelAnalysis.fitKernels{cellno,ps(p)},'Color',colors(ps(p),:),'LineWidth',2,'LineStyle','-')
end
line([0 0],[-limit limit],'Color',[.5 .5 .5],'LineStyle','--')
line([-.5 .9],[0 0],'Color',[.5 .5 .5],'LineStyle','--')
xlim([-.5 .9])
ylim([-limit limit])
prettyPlot(gca)
title('Stimulus')

subplot(1,4,2)
hold on;
pv = find(contains(kernelAnalysis.features,'value'));
pb = find(contains(kernelAnalysis.features,'block'));
plot(kernelAnalysis.windows.value*Fs,kernelAnalysis.fitKernels{cellno,pv},'Color',colors(pv,:),'LineWidth',2,'LineStyle',':')
plot(kernelAnalysis.windows.value*Fs,kernelAnalysis.fitKernels{cellno,pb},'Color',colors(pb,:),'LineWidth',2,'LineStyle','-')
line([0 0],[-limit limit],'Color',[.5 .5 .5],'LineStyle','--')
line([-.5 .9],[0 0],'Color',[.5 .5 .5],'LineStyle','--')
xlim([-.5 .9])
ylim([-limit limit])
prettyPlot(gca)
title('Value & block')

subplot(1,4,3)
hold on;
pa = find(contains(kernelAnalysis.features,'action'));
pc = find(contains(kernelAnalysis.features,'choice'));
plot(kernelAnalysis.windows.value*Fs,kernelAnalysis.fitKernels{cellno,pa},'Color',colors(pa,:),'LineWidth',2,'LineStyle','-')
plot(kernelAnalysis.windows.value*Fs,kernelAnalysis.fitKernels{cellno,pc},'Color',colors(pc,:),'LineWidth',2,'LineStyle',':')
line([0 0],[-limit limit],'Color',[.5 .5 .5],'LineStyle','--')
line([-.5 .9],[0 0],'Color',[.5 .5 .5],'LineStyle','--')
xlim([-.5 .9])
ylim([-limit limit])
prettyPlot(gca)
title('Action & choice')

subplot(1,4,4)
hold on
po = find(contains(kernelAnalysis.features,'outcome'));
for p = 1:length(po)
    plot(kernelAnalysis.windows.value*Fs,kernelAnalysis.fitKernels{cellno,po(p)},'Color',colors(po(p),:),'LineWidth',2,'LineStyle','-')
end
line([0 0],[-limit limit],'Color',[.5 .5 .5],'LineStyle','--')
line([-.5 .9],[0 0],'Color',[.5 .5 .5],'LineStyle','--')
xlim([-.5 .9])
ylim([-limit limit])
prettyPlot(gca)
title('Outcome')

% subplot(1,5,5)
% hold on
% pvel = find(contains(kernelAnalysis.features,'velocity'));
% for p = 1:length(pvel)
%     plot(kernelAnalysis.windows.action*Fs,kernelAnalysis.fitKernels{cellno,pvel(p)},'Color',colors(pvel(p),:),'LineWidth',2,'LineStyle','-')
% end
% line([0 0],[-limit limit],'Color',[.5 .5 .5],'LineStyle','--')
% line([-.5 .9],[0 0],'Color',[.5 .5 .5],'LineStyle','--')
% xlim([-.5 .9])
% ylim([-limit limit])
% prettyPlot(gca)
% title('Velocity')

% plot model vs data

% fetch trial types
et = behavioralData;
contrasts = getUniqueContrasts(expInfo);
nt = length(et.eventTimes(1).daqTime);
[whichTrials, ~] = selectCondition(expInfo, contrasts, behavioralData, initTrialConditions('movementTime','late'));

trialCorrectChoice = expInfo.block.events.correctResponseValues(1:nt);
trueStimuli = expInfo.block.events.contrastValues(1:nt);
trueStimuli(trueStimuli == 0) = eps;
trueStimuli(abs(trueStimuli) < .05) = ...
    trueStimuli(abs(trueStimuli) < .05).* trialCorrectChoice(abs(trueStimuli) < .05);

trueChoices = et.wheelMoves.epochs(5).moveDir;
trueBlock = expInfo.block.events.highRewardSideValues(1:nt);
trueValue(trueBlock .* sign(trueStimuli) > 0) = 2;
trueValue(trueBlock .* sign(trueStimuli) < 0) = 1;

hls = logical((trueStimuli < -.12) .* whichTrials);
lls = logical((trueStimuli >= -.12 & trueStimuli < -eps) .* whichTrials);
zs = logical((trueStimuli > -.005 & trueStimuli < .005) .* whichTrials);
lrs = logical((trueStimuli <= .12 & trueStimuli > eps) .* whichTrials);
hrs = logical((trueStimuli > .12) .* whichTrials);
ls = logical((trueStimuli < -.005) .* whichTrials);
rs = logical((trueStimuli > .005) .* whichTrials);

lc = logical((trueChoices < 0) .* whichTrials);
rc = logical((trueChoices > 0) .* whichTrials);

hv = logical((trueValue == 2) .* whichTrials);
lv = logical((trueValue == 1) .* whichTrials);
lb = logical((trueBlock == 1) .* whichTrials);
rb = logical((trueBlock == -1) .* whichTrials);

[co, ~] = selectCondition(expInfo, contrasts, et, initTrialConditions('responseType','correct','movementTime','late'));
[io, ~] = selectCondition(expInfo, contrasts, et, initTrialConditions('responseType','incorrect','movementTime','late'));

co = logical(co);
io = logical(io);

% fetch responses and model output
r_hls = neuralData.eta.alignedResps{1}(hls,:,cellno);
r_lls = neuralData.eta.alignedResps{1}(lls,:,cellno);
r_zs = neuralData.eta.alignedResps{1}(zs,:,cellno);
r_lrs = neuralData.eta.alignedResps{1}(lrs,:,cellno);
r_hrs = neuralData.eta.alignedResps{1}(hrs,:,cellno);
r_ls = neuralData.eta.alignedResps{1}(ls,:,cellno);
r_rs = neuralData.eta.alignedResps{1}(rs,:,cellno);
r_lc = neuralData.eta.alignedResps{2}(lc,:,cellno);
r_rc = neuralData.eta.alignedResps{2}(rc,:,cellno);
r_hv = neuralData.eta.alignedResps{1}(hv,:,cellno);
r_lv = neuralData.eta.alignedResps{1}(lv,:,cellno);
r_co = neuralData.eta.alignedResps{3}(co,:,cellno);
r_io = neuralData.eta.alignedResps{3}(io,:,cellno);
r_lb = neuralData.eta.alignedResps{1}(lb,:,cellno);
r_rb = neuralData.eta.alignedResps{1}(rb,:,cellno);

m_hls = kernelAnalysis.eta.alignedResps{1}(hls,:,cellno);
m_lls = kernelAnalysis.eta.alignedResps{1}(lls,:,cellno);
m_zs = kernelAnalysis.eta.alignedResps{1}(zs,:,cellno);
m_lrs = kernelAnalysis.eta.alignedResps{1}(lrs,:,cellno);
m_hrs = kernelAnalysis.eta.alignedResps{1}(hrs,:,cellno);
m_ls = kernelAnalysis.eta.alignedResps{1}(ls,:,cellno);
m_rs = kernelAnalysis.eta.alignedResps{1}(rs,:,cellno);
m_lc = kernelAnalysis.eta.alignedResps{2}(lc,:,cellno);
m_rc = kernelAnalysis.eta.alignedResps{2}(rc,:,cellno);
m_hv = kernelAnalysis.eta.alignedResps{1}(hv,:,cellno);
m_lv = kernelAnalysis.eta.alignedResps{1}(lv,:,cellno);
m_co = kernelAnalysis.eta.alignedResps{3}(co,:,cellno);
m_io = kernelAnalysis.eta.alignedResps{3}(io,:,cellno);
m_lb = kernelAnalysis.eta.alignedResps{1}(lb,:,cellno);
m_rb = kernelAnalysis.eta.alignedResps{1}(rb,:,cellno);

% psth vs model
figure;
set(gcf,'position',[1165 950 1220 270])
nlim = 1.1*max(abs([nanmean(r_ls) nanmean(r_rs) nanmean(r_hv) nanmean(r_lv) nanmean(r_rc) nanmean(r_lc) nanmean(r_co) nanmean(r_io)]));
subplot(1,5,1)
hold on
plotSignal(neuralData.eta.eventWindow,nanmean(r_ls,1),nanmean(r_ls,1)+nanstd(r_ls)/sqrt(size(r_ls,1)),nanmean(r_ls,1)-nanstd(r_ls)/sqrt(size(r_ls,1)),[0 .4 1],'-');
plotSignal(neuralData.eta.eventWindow,nanmean(r_rs,1),nanmean(r_rs,1)+nanstd(r_rs)/sqrt(size(r_rs,1)),nanmean(r_rs,1)-nanstd(r_rs)/sqrt(size(r_rs,1)),[1 0 0],'-');
plot(neuralData.eta.eventWindow,nanmean(m_ls,1),'Color',[0 .4 1],'LineWidth',2,'LineStyle','-')
plot(neuralData.eta.eventWindow,nanmean(m_rs,1),'Color',[1 0 0],'LineWidth',2,'LineStyle','-')
line([0 0],[-nlim nlim],'Color',[.5 .5 .5],'LineStyle','--')
xlim([-.5 .9])
ylim([-.05 nlim])
xlabel('Time from stimulus (s)')
prettyPlot(gca)

subplot(1,5,2)
hold on
plotSignal(neuralData.eta.eventWindow,nanmean(r_lb,1),nanmean(r_lb,1)+nanstd(r_lb)/sqrt(size(r_lb,1)),nanmean(r_lb,1)-nanstd(r_lb)/sqrt(size(r_lb,1)),[0.1 0.7 0.1],'-');
plotSignal(neuralData.eta.eventWindow,nanmean(r_rb,1),nanmean(r_rb,1)+nanstd(r_rb)/sqrt(size(r_rb,1)),nanmean(r_rb,1)-nanstd(r_rb)/sqrt(size(r_rb,1)),[1 .6 0],'-');
plot(neuralData.eta.eventWindow,nanmean(m_lb,1),'Color',[0.1 0.7 0.1],'LineWidth',2,'LineStyle','-')
plot(neuralData.eta.eventWindow,nanmean(m_rb,1),'Color',[1 .6 0],'LineWidth',2,'LineStyle','-')
line([0 0],[-nlim nlim],'Color',[.5 .5 .5],'LineStyle','--')
xlim([-.5 .9])
ylim([-.05 nlim])
xlabel('Time from stimulus (s)')
prettyPlot(gca)

subplot(1,5,3)
hold on
plotSignal(neuralData.eta.eventWindow,nanmean(r_hv,1),nanmean(r_hv,1)+nanstd(r_hv)/sqrt(size(r_hv,1)),nanmean(r_hv,1)-nanstd(r_hv)/sqrt(size(r_hv,1)),[0 0 1],'-');
plotSignal(neuralData.eta.eventWindow,nanmean(r_lv,1),nanmean(r_lv,1)+nanstd(r_lv)/sqrt(size(r_lv,1)),nanmean(r_lv,1)-nanstd(r_lv)/sqrt(size(r_lv,1)),[0 .5 1],'-');
plot(neuralData.eta.eventWindow,nanmean(m_hv,1),'Color',[0 0 1],'LineWidth',2,'LineStyle','-')
plot(neuralData.eta.eventWindow,nanmean(m_lv,1),'Color',[0 .5 1],'LineWidth',2,'LineStyle','-')
line([0 0],[-nlim nlim],'Color',[.5 .5 .5],'LineStyle','--')
xlim([-.5 .9])
ylim([-.05 nlim])
xlabel('Time from stimulus (s)')
prettyPlot(gca)

subplot(1,5,4)
hold on
plotSignal(neuralData.eta.eventWindow,nanmean(r_lc,1),nanmean(r_lc,1)+nanstd(r_lc)/sqrt(size(r_lc,1)),nanmean(r_lc,1)-nanstd(r_lc)/sqrt(size(r_lc,1)),[.5 0 1],'-');
plotSignal(neuralData.eta.eventWindow,nanmean(r_rc,1),nanmean(r_rc,1)+nanstd(r_rc)/sqrt(size(r_rc,1)),nanmean(r_rc,1)-nanstd(r_rc)/sqrt(size(r_rc,1)),[1 0 .5],'-');
plot(neuralData.eta.eventWindow,nanmean(m_lc,1),'Color',[.5 0 1],'LineWidth',2,'LineStyle','-')
plot(neuralData.eta.eventWindow,nanmean(m_rc,1),'Color',[1 0 .5],'LineWidth',2,'LineStyle','-')
line([0 0],[-nlim nlim],'Color',[.5 .5 .5],'LineStyle','--')
xlim([-.5 .9])
ylim([-.05 nlim])
xlabel('Time from move (s)')
prettyPlot(gca)

subplot(1,5,5)
hold on
plotSignal(neuralData.eta.eventWindow,nanmean(r_co,1),nanmean(r_co,1)+nanstd(r_co)/sqrt(size(r_co,1)),nanmean(r_co,1)-nanstd(r_co)/sqrt(size(r_co,1)),[0 .5 0],'-');
plotSignal(neuralData.eta.eventWindow,nanmean(r_io,1),nanmean(r_io,1)+nanstd(r_io)/sqrt(size(r_io,1)),nanmean(r_io,1)-nanstd(r_io)/sqrt(size(r_io,1)),[.7 .1 0],'-');
plot(neuralData.eta.eventWindow,nanmean(m_co,1),'Color',[.0 .5 0],'LineWidth',2,'LineStyle','-')
plot(neuralData.eta.eventWindow,nanmean(m_io,1),'Color',[.7 .1 0],'LineWidth',2,'LineStyle','-')
line([0 0],[-nlim nlim],'Color',[.5 .5 .5],'LineStyle','--')
xlim([-.5 .9])
ylim([-.05 nlim])
xlabel('Time from feedback (s)')
prettyPlot(gca)


% % psth 1st vs 2nd half
% figure;
% set(gcf,'position',[1165 950 1220 270])
% nlim = 1.1*max(abs([...
%     nanmean(r_ls(1:round(length(r_ls)/2),:),1) ...
%     nanmean(r_ls(round(length(r_ls)/2)+1:end,:),1) ...
%     nanmean(r_rs(1:round(length(r_rs)/2),:),1), ...
%     nanmean(r_rs(round(length(r_rs)/2)+1:end,:),1) ...
%     nanmean(r_hv(1:round(length(r_hv)/2),:),1) ...
%     nanmean(r_hv(round(length(r_hv)/2)+1:end,:),1) ...
%     nanmean(r_lv(1:round(length(r_lv)/2),:),1) ...
%     nanmean(r_lv(round(length(r_lv)/2)+1:end,:),1) ...
%     nanmean(r_lc(1:round(length(r_lc)/2),:),1) ...
%     nanmean(r_lc(round(length(r_lc)/2)+1:end,:),1) ...
%     nanmean(r_rc(1:round(length(r_rc)/2),:),1) ...
%     nanmean(r_rc(round(length(r_rc)/2)+1:end,:),1) ...
%     nanmean(r_co(1:round(length(r_co)/2),:),1) ...
%     nanmean(r_co(round(length(r_co)/2)+1:end,:),1) ...
%     nanmean(r_io(1:round(length(r_io)/2),:),1) ...
%     nanmean(r_io(round(length(r_io)/2)+1:end,:),1)]));
% subplot(1,4,1)
% hold on
% plotSignal(neuralData.eta.eventWindow,nanmean(r_ls(1:round(length(r_ls)/2),:),1),nanmean(r_ls(1:round(length(r_ls)/2),:),1)+nanstd(r_ls(1:round(length(r_ls)/2),:))/sqrt(size(r_ls(1:round(length(r_ls)/2),:),1)),nanmean(r_ls(1:round(length(r_ls)/2),:),1)-nanstd(r_ls(1:round(length(r_ls)/2),:))/sqrt(size(r_ls(1:round(length(r_ls)/2),:),1)),[0 .4 1],'-');
% plotSignal(neuralData.eta.eventWindow,nanmean(r_ls(round(length(r_ls)/2)+1:end,:),1),nanmean(r_ls(round(length(r_ls)/2)+1:end,:),1)+nanstd(r_ls(round(length(r_ls)/2)+1:end,:))/sqrt(size(r_ls(round(length(r_ls)/2)+1:end,:),1)),nanmean(r_ls(round(length(r_ls)/2)+1:end,:),1)-nanstd(r_ls(round(length(r_ls)/2)+1:end,:))/sqrt(size(r_ls(round(length(r_ls)/2)+1:end,:),1)),[0 .4 1],':');
% plotSignal(neuralData.eta.eventWindow,nanmean(r_rs(1:round(length(r_rs)/2),:),1),nanmean(r_rs(1:round(length(r_rs)/2),:),1)+nanstd(r_rs(1:round(length(r_rs)/2),:))/sqrt(size(r_rs(1:round(length(r_rs)/2),:),1)),nanmean(r_rs(1:round(length(r_rs)/2),:),1)-nanstd(r_rs(1:round(length(r_rs)/2),:))/sqrt(size(r_rs(1:round(length(r_rs)/2),:),1)),[1 0 0],'-');
% plotSignal(neuralData.eta.eventWindow,nanmean(r_rs(round(length(r_rs)/2)+1:end,:),1),nanmean(r_rs(round(length(r_rs)/2)+1:end,:),1)+nanstd(r_rs(round(length(r_rs)/2)+1:end,:))/sqrt(size(r_rs(round(length(r_rs)/2)+1:end,:),1)),nanmean(r_rs(round(length(r_rs)/2)+1:end,:),1)-nanstd(r_rs(round(length(r_rs)/2)+1:end,:))/sqrt(size(r_rs(round(length(r_rs)/2)+1:end,:),1)),[1 0 0],':');
% line([0 0],[-nlim nlim],'Color',[.5 .5 .5],'LineStyle','--')
% xlim([-.5 .9])
% ylim([-.05 nlim])
% xlabel('Time from stimulus (s)')
% prettyPlot(gca)
% 
% subplot(1,4,2)
% hold on
% plotSignal(neuralData.eta.eventWindow,nanmean(r_hv(1:round(length(r_hv)/2),:),1),nanmean(r_hv(1:round(length(r_hv)/2),:),1)+nanstd(r_hv(1:round(length(r_hv)/2),:))/sqrt(size(r_hv(1:round(length(r_hv)/2),:),1)),nanmean(r_hv(1:round(length(r_hv)/2),:),1)-nanstd(r_hv(1:round(length(r_hv)/2),:))/sqrt(size(r_hv(1:round(length(r_hv)/2),:),1)),[0 0 1],'-');
% plotSignal(neuralData.eta.eventWindow,nanmean(r_hv(round(length(r_hv)/2)+1:end,:),1),nanmean(r_hv(round(length(r_hv)/2)+1:end,:),1)+nanstd(r_hv(round(length(r_hv)/2)+1:end,:))/sqrt(size(r_hv(round(length(r_hv)/2)+1:end,:),1)),nanmean(r_hv(round(length(r_hv)/2)+1:end,:),1)-nanstd(r_hv(round(length(r_hv)/2)+1:end,:))/sqrt(size(r_hv(round(length(r_hv)/2)+1:end,:),1)),[0 0 1],':');
% plotSignal(neuralData.eta.eventWindow,nanmean(r_lv(1:round(length(r_lv)/2),:),1),nanmean(r_lv(1:round(length(r_lv)/2),:),1)+nanstd(r_lv(1:round(length(r_lv)/2),:))/sqrt(size(r_lv(1:round(length(r_lv)/2),:),1)),nanmean(r_lv(1:round(length(r_lv)/2),:),1)-nanstd(r_lv(1:round(length(r_lv)/2),:))/sqrt(size(r_lv(1:round(length(r_lv)/2),:),1)),[0 .5 1],'-');
% plotSignal(neuralData.eta.eventWindow,nanmean(r_lv(round(length(r_lv)/2)+1:end,:),1),nanmean(r_lv(round(length(r_lv)/2)+1:end,:),1)+nanstd(r_lv(round(length(r_lv)/2)+1:end,:))/sqrt(size(r_lv(round(length(r_lv)/2)+1:end,:),1)),nanmean(r_lv(round(length(r_lv)/2)+1:end,:),1)-nanstd(r_lv(round(length(r_lv)/2)+1:end,:))/sqrt(size(r_lv(round(length(r_lv)/2)+1:end,:),1)),[0 .5 1],':');
% line([0 0],[-nlim nlim],'Color',[.5 .5 .5],'LineStyle','--')
% xlim([-.5 .9])
% ylim([-.05 nlim])
% xlabel('Time from stimulus (s)')
% prettyPlot(gca)
% 
% subplot(1,4,3)
% hold on
% plotSignal(neuralData.eta.eventWindow,nanmean(r_lc(1:round(length(r_lc)/2),:),1),nanmean(r_lc(1:round(length(r_lc)/2),:),1)+nanstd(r_lc(1:round(length(r_lc)/2),:))/sqrt(size(r_lc(1:round(length(r_lc)/2),:),1)),nanmean(r_lc(1:round(length(r_lc)/2),:),1)-nanstd(r_lc(1:round(length(r_lc)/2),:))/sqrt(size(r_lc(1:round(length(r_lc)/2),:),1)),[.5 0 1],'-');
% plotSignal(neuralData.eta.eventWindow,nanmean(r_lc(round(length(r_lc)/2)+1:end,:),1),nanmean(r_lc(round(length(r_lc)/2)+1:end,:),1)+nanstd(r_lc(round(length(r_lc)/2)+1:end,:))/sqrt(size(r_lc(round(length(r_lc)/2)+1:end,:),1)),nanmean(r_lc(round(length(r_lc)/2)+1:end,:),1)-nanstd(r_lc(round(length(r_lc)/2)+1:end,:))/sqrt(size(r_lc(round(length(r_lc)/2)+1:end,:),1)),[.5 0 1],':');
% plotSignal(neuralData.eta.eventWindow,nanmean(r_rc(1:round(length(r_rc)/2),:),1),nanmean(r_rc(1:round(length(r_rc)/2),:),1)+nanstd(r_rc(1:round(length(r_rc)/2),:))/sqrt(size(r_rc(1:round(length(r_rc)/2),:),1)),nanmean(r_rc(1:round(length(r_rc)/2),:),1)-nanstd(r_rc(1:round(length(r_rc)/2),:))/sqrt(size(r_rc(1:round(length(r_rc)/2),:),1)),[1 0 .5],'-');
% plotSignal(neuralData.eta.eventWindow,nanmean(r_rc(round(length(r_rc)/2)+1:end,:),1),nanmean(r_rc(round(length(r_rc)/2)+1:end,:),1)+nanstd(r_rc(round(length(r_rc)/2)+1:end,:))/sqrt(size(r_rc(round(length(r_rc)/2)+1:end,:),1)),nanmean(r_rc(round(length(r_rc)/2)+1:end,:),1)-nanstd(r_rc(round(length(r_rc)/2)+1:end,:))/sqrt(size(r_rc(round(length(r_rc)/2)+1:end,:),1)),[1 0 .5],':');
% line([0 0],[-nlim nlim],'Color',[.5 .5 .5],'LineStyle','--')
% xlim([-.5 .9])
% ylim([-.05 nlim])
% xlabel('Time from move (s)')
% prettyPlot(gca)
% 
% subplot(1,4,4)
% hold on
% plotSignal(neuralData.eta.eventWindow,nanmean(r_co(1:round(length(r_co)/2),:),1),nanmean(r_co(1:round(length(r_co)/2),:),1)+nanstd(r_co(1:round(length(r_co)/2),:))/sqrt(size(r_co(1:round(length(r_co)/2),:),1)),nanmean(r_co(1:round(length(r_co)/2),:),1)-nanstd(r_co(1:round(length(r_co)/2),:))/sqrt(size(r_co(1:round(length(r_co)/2),:),1)),[0 .5 0],'-');
% plotSignal(neuralData.eta.eventWindow,nanmean(r_co(round(length(r_co)/2)+1:end,:),1),nanmean(r_co(round(length(r_co)/2)+1:end,:),1)+nanstd(r_co(round(length(r_co)/2)+1:end,:))/sqrt(size(r_co(round(length(r_co)/2)+1:end,:),1)),nanmean(r_co(round(length(r_co)/2)+1:end,:),1)-nanstd(r_co(round(length(r_co)/2)+1:end,:))/sqrt(size(r_co(round(length(r_co)/2)+1:end,:),1)),[0 .5 0],':');
% plotSignal(neuralData.eta.eventWindow,nanmean(r_io(1:round(length(r_io)/2),:),1),nanmean(r_io(1:round(length(r_io)/2),:),1)+nanstd(r_io(1:round(length(r_io)/2),:))/sqrt(size(r_io(1:round(length(r_io)/2),:),1)),nanmean(r_io(1:round(length(r_io)/2),:),1)-nanstd(r_io(1:round(length(r_io)/2),:))/sqrt(size(r_io(1:round(length(r_io)/2),:),1)),[.7 .1 0],'-');
% plotSignal(neuralData.eta.eventWindow,nanmean(r_io(round(length(r_io)/2)+1:end,:),1),nanmean(r_io(round(length(r_io)/2)+1:end,:),1)+nanstd(r_io(round(length(r_io)/2)+1:end,:))/sqrt(size(r_io(round(length(r_io)/2)+1:end,:),1)),nanmean(r_io(round(length(r_io)/2)+1:end,:),1)-nanstd(r_io(round(length(r_io)/2)+1:end,:))/sqrt(size(r_io(round(length(r_io)/2)+1:end,:),1)),[.7 .1 0],':');
% line([0 0],[-nlim nlim],'Color',[.5 .5 .5],'LineStyle','--')
% xlim([-.5 .9])
% ylim([-.05 nlim])
% xlabel('Time from feedback (s)')
% prettyPlot(gca)
% 

%%

[~, earlyTrials] = selectCondition(expInfo, contrasts, behavioralData, initTrialConditions('repeatType', 'random', 'movementTime','early'));
[~, lateTrials] = selectCondition(expInfo, contrasts, behavioralData, initTrialConditions('repeatType', 'random', 'movementTime','late'));

figure;
for tn = 1:96
    subplot(8,12,tn)
    hold on;
%     tn = 2;
    plot(neuralData.eta.eventWindow,neuralData.eta.alignedResps{1}(earlyTrials(tn),:,3)','LineStyle','-','Color','k','LineWidth',1)
    plot(neuralData.eta.eventWindow,kernelAnalysis.eta.alignedResps{1}(earlyTrials(tn),:,3)','LineStyle',':','Color',[0 0 0],'LineWidth',1)
    line([0 0],[-.2 1.2],'Color',[.5 .5 .5],'LineStyle','--')
    line([0.8 0.8],[-.2 1.2],'Color',[.5 .5 .5],'LineStyle','--')
    prettyPlot(gca)
    ylim([-.2 1.2])
    xlim([-.5 2])
    
end

%%
figure;
hold on;
for m = 1:length(mouseList)
    
    %load data
    mouseName = char(mouseList{m});
    expDate = char(expList{m}{1});
    expNum = expList{m}{2};
    [expRef, ~] = data.constructExpRef(mouseName,expDate,expNum);
    fprintf('Loading %s...\n',expRef)  
    sessDir = fullfile('G:\Workspaces',mouseName,expDate,num2str(expNum));
    cd(sessDir)
%     load('dataset.mat')
    load('kernelAnalysis.mat')

    line([m m],100*[prctile(kernelAnalysis.maxEV,10) prctile(kernelAnalysis.maxEV,90)],'Color','k','LineWidth',2)
    line([m m],100*[prctile(kernelAnalysis.maxEV,25) prctile(kernelAnalysis.maxEV,75)],'Color','k','LineWidth',8)
    line([m-.1 m+.1],100*[prctile(kernelAnalysis.maxEV,50) prctile(kernelAnalysis.maxEV,50)],'Color','k','LineWidth',2)
    
end

%% collect all sessions
allStrengths = {};
for m = 1:length(mouseList)
    
    %load data
    mouseName = char(mouseList{m});
    expDate = char(expList{m}{1});
    expNum = expList{m}{2};
    [expRef, ~] = data.constructExpRef(mouseName,expDate,expNum);
    fprintf('Loading %s...\n',expRef)  
    sessDir = fullfile('G:\Workspaces',mouseName,expDate,num2str(expNum));
    cd(sessDir)
%     load('dataset.mat')
    load('kernelAnalysis_new.mat')
    allMaxEV{m,1} = kernelAnalysis.maxEV;
    propFeatures(m,:) = sum(kernelAnalysis.cellFeatureStrength(kernelAnalysis.maxEV > .01,:) > .01)./(sum(kernelAnalysis.maxEV > .01))*100;
    allStrengths{m,1} = kernelAnalysis.cellFeatureStrength(kernelAnalysis.maxEV > .01,:) > .01;
end

%% bar chart
[~,i] = sort(mean(propFeatures),'descend');
lbls = {'Stimulus','Block','Value','Action','Choice','Outcome'};

figure;
cm = parula(6);
b = bar(mean(propFeatures(:,i)),'FaceColor','flat');
xtips1 = b(1).XEndPoints;
ytips1 = mean(propFeatures(:,i)) + std(propFeatures(:,i))/sqrt(length(propFeatures(:,i)));
vals = string(round(b(1).YData,2));
text(xtips1,ytips1,vals,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')
for j=1:size(propFeatures,2)
    b.CData(j,:) = cm(j,:);
end
hold on
errorbar(mean(propFeatures(:,i)),std(propFeatures(:,i))/sqrt(length(propFeatures(:,i))),'LineStyle','none','Color','k')
prettyPlot(gca)
set(gca, 'XTickLabels', lbls(i))
xlabel('Task feature')
ylabel('Percentage of cells')

%% venn diagram
%http://bioinformatics.psb.ugent.be/webtools/Venn/
as = cat(1,allStrengths{:});
as = as(:,[1 4 5 6]);

for f = 1:size(as,2)
    list{f} = find(as(:,f) == 1);
end

%% leave-one-out plotting
%select features
features = {'stimulus' 'value' 'block' 'action' 'choice' 'outcome'};

% design matrix
[predictors, windows] = getPredictors(...
    expInfo, ...
    behavioralData, ...
    features, ...
    0.1);

% set up regression matrices
X = makeToeplitz(neuralData.respTimes, predictors, windows);
LOO_kernelFunctions = kernelAnalysis.kernelFunctions;
LOO_kernelFunctions(1:15,:) = 0;
predictedActivity = X*LOO_kernelFunctions;

Fs = 0.1;
timeBefore = 2;
timeAfter = 2;
events = {'stimulusOnTimes' 'firstMoveTimes' 'feedbackTimes' 'interactiveOnTimes'};

nt = numel(expInfo.block.events.endTrialTimes);

for ev = 1:length(events)

    % get event time window for each trial
    eventWindow = -timeBefore:Fs:timeAfter;
    if strcmp(events(ev), 'firstMoveTimes')
        firstMoveTimes = min([behavioralData.wheelMoves.epochs(2).onsetTimes; behavioralData.wheelMoves.epochs(3).onsetTimes]);
        periEventTimes = bsxfun(@plus,firstMoveTimes',eventWindow);
    else
        periEventTimes = bsxfun(@plus,behavioralData.eventTimes(strcmp({behavioralData.eventTimes.event},events{ev})).daqTime',eventWindow);
    end
    periEventTimes = periEventTimes(1:nt,:);

    %initialize alignedResp cell
    loo_alignedResps{ev} = zeros(nt, length(eventWindow), size(predictedActivity,2));

    %grab cell responses associated with the event time windows 
    %(size is nTrials x windowLength x nCells)
    for iCell = 1:size(kernelAnalysis.predictedActivity,2)        
        loo_alignedResps{ev}(:,:,iCell) = interp1(neuralData.respTimes,predictedActivity(:,iCell),periEventTimes,'previous');
    end
end


