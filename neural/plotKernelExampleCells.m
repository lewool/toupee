%% remove features

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
LOO_kernelFunctions(136:165,:) = 0;
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

%% chose which model to plot

% whichAligned = kernelAnalysis.eta.alignedResps; linst = '-';
whichAligned = loo_alignedResps; linst = '--';

%% plot full model

cellno = 58;
Fs = 0.1;

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

m_hls = whichAligned{1}(hls,:,cellno);
m_lls = whichAligned{1}(lls,:,cellno);
m_zs = whichAligned{1}(zs,:,cellno);
m_lrs = whichAligned{1}(lrs,:,cellno);
m_hrs = whichAligned{1}(hrs,:,cellno);
m_ls = whichAligned{1}(ls,:,cellno);
m_rs = whichAligned{1}(rs,:,cellno);
m_lc = whichAligned{2}(lc,:,cellno);
m_rc = whichAligned{2}(rc,:,cellno);
m_hv = whichAligned{1}(hv,:,cellno);
m_lv = whichAligned{1}(lv,:,cellno);
m_co = whichAligned{3}(co,:,cellno);
m_io = whichAligned{3}(io,:,cellno);
m_lb = whichAligned{1}(lb,:,cellno);
m_rb = whichAligned{1}(rb,:,cellno);

colors = {[0 .4 1], [1 0 0]; ...
          [.1 .7 .1], [1 .6 0]; ...
          [0 0 1], [0 .5 1]; ...
          [.5 0 1], [1 0 .5]; ...
          [0 .5 0], [.7 .1 0]; ...
          };

whichR = {r_ls, r_rs; ...
          r_lb, r_rb; ...
          r_hv, r_lv; ...
          r_lc, r_rc; ...
          r_co, r_io; ...
          };
      
whichM = {m_ls, m_rs; ...
          m_lb, m_rb; ...
          m_hv, m_lv; ...
          m_lc, m_rc; ...
          m_co, m_io; ...
          };

whichETA = {'stimulus' 'stimulus' 'stimulus' 'movement' 'feedback'};
figure;   
set(gcf,'position',[1165 950 1220 270])
nlim = 1.1*max(abs([nanmean(r_ls) nanmean(r_rs) nanmean(r_hv) nanmean(r_lv) nanmean(r_rc) nanmean(r_lc) nanmean(r_co) nanmean(r_io)]));
nlim = .7;
for s = 1:length(whichR)
    subplot(1,5,s)
    hold on
    plotSignal(neuralData.eta.eventWindow,nanmean(whichR{s,1},1),nanmean(whichR{s,1},1)+nanstd(whichR{s,1})/sqrt(size(whichR{s,1},1)),nanmean(whichR{s,1},1)-nanstd(whichR{s,1})/sqrt(size(whichR{s,1},1)),colors{s,1},'-');
    plotSignal(neuralData.eta.eventWindow,nanmean(whichR{s,2},1),nanmean(whichR{s,2},1)+nanstd(whichR{s,2})/sqrt(size(whichR{s,2},1)),nanmean(whichR{s,2},1)-nanstd(whichR{s,2})/sqrt(size(whichR{s,2},1)),colors{s,2},'-');
    plot(neuralData.eta.eventWindow,nanmean(whichM{s,1},1),'Color',colors{s,1},'LineWidth',2,'LineStyle',linst)
    plot(neuralData.eta.eventWindow,nanmean(whichM{s,2},1),'Color',colors{s,2},'LineWidth',2,'LineStyle',linst)
    line([0 0],[-nlim nlim],'Color',[.5 .5 .5],'LineStyle','--')
    xlim([-.5 .9])
    ylim([-.05 nlim])
    xlabel(strcat({'Time from '},whichETA{s},{' (s)'}))
    prettyPlot(gca)
end

%% remove kernels

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

%% plot kernel visual timeline
figure;
ew=kernelAnalysis.windows.stimulus*.1;
yl = [-.15 .7];
cellno = 623;
subplot(1,5,1)
sk = reshape(kernelAnalysis.kernelFunctions(1:75,cellno),[15 5]);
plot(ew,mean(sk(:,1:2),2),'Color',[0 .4 1]);
hold on;
plot(ew,mean(sk(:,3),2),'Color',[.6 .6 .6]);
plot(ew,mean(sk(:,4:5),2),'Color',[1 0 0]);
line([0 0],yl,'Color',[.5 .5 .5],'LineStyle','--')
subplot(1,5,2)
plot(ew,kernelAnalysis.kernelFunctions(76:90,cellno),'Color',[0 0 0]);
hold on;
line([0 0],yl,'Color',[.5 .5 .5],'LineStyle','--')
subplot(1,5,3)
plot(ew,kernelAnalysis.kernelFunctions(91:105,cellno),'Color',[0 0 0]);
hold on;
line([0 0],yl,'Color',[.5 .5 .5],'LineStyle','--')
subplot(1,5,4)
plot(ew,kernelAnalysis.kernelFunctions(106:120,cellno),'Color',[0 0 0]);
hold on;
plot(ew,kernelAnalysis.kernelFunctions(121:135,cellno),'Color',[0 0 0],'linestyle',':');
line([0 0],yl,'Color',[.5 .5 .5],'LineStyle','--')
subplot(1,5,5)
plot(ew,reshape(kernelAnalysis.kernelFunctions(136:165,cellno),[15 2]));
hold on;
line([0 0],yl,'Color',[.5 .5 .5],'LineStyle','--')
for s = 1:5
    subplot(1,5,s)
    ylim(yl)
    prettyPlot(gca)
end

%% new timeline 
% figure;plot(kernelAnalysis.eta.eventWindow,mean(whichAligned{1}(1,:,cellno),1))
% hold on
% plot(kernelAnalysis.eta.eventWindow,mean(neuralData.eta.alignedResps{1}(1,:,cellno),1))
n = 623;

stimFunctions = X(:,1:125)*kernelAnalysis.kernelFunctions(1:125,n);
valueFunctions = X(:,126:150)*kernelAnalysis.kernelFunctions(126:150,n);
blockFunctions = X(:,151:175)*kernelAnalysis.kernelFunctions(151:175,n);
actionFunctions = X(:,176:190)*kernelAnalysis.kernelFunctions(176:190,n);
choiceFunctions = X(:,191:205)*kernelAnalysis.kernelFunctions(191:205,n);
feedbackFunctions = X(:,206:235)*kernelAnalysis.kernelFunctions(206:235,n);



Fs = 0.1;
timeBefore = 2;
timeAfter = 3;
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

    %grab cell responses associated with the event time windows 
    %(size is nTrials x windowLength x nCells)
    for iCell = 1:size(feedbackFunctions,2)        
        alignedPreds_stim{ev}(:,:,iCell) = interp1(neuralData.respTimes,stimFunctions(:,iCell),periEventTimes,'previous');
        alignedPreds_block{ev}(:,:,iCell) = interp1(neuralData.respTimes,blockFunctions(:,iCell),periEventTimes,'previous');
        alignedPreds_value{ev}(:,:,iCell) = interp1(neuralData.respTimes,valueFunctions(:,iCell),periEventTimes,'previous');
        alignedPreds_action{ev}(:,:,iCell) = interp1(neuralData.respTimes,actionFunctions(:,iCell),periEventTimes,'previous');
        alignedPreds_choice{ev}(:,:,iCell) = interp1(neuralData.respTimes,choiceFunctions(:,iCell),periEventTimes,'previous');
        alignedPreds_feedback{ev}(:,:,iCell) = interp1(neuralData.respTimes,feedbackFunctions(:,iCell),periEventTimes,'previous');
    end
end

%%

figure;
hold on;
subplot(1,3,1)
hold on
line([0 0],[-.4 1.5],'Color',[.5 .5 .5],'LineStyle','--')
line([1 1],[-.4 1.5],'Color',[.5 .5 .5],'LineStyle','--')

plot(neuralData.eta.eventWindow,mean(alignedPreds_stim{1}(:,:,1),1))
plot(neuralData.eta.eventWindow,mean(alignedPreds_block{1}(:,:,1),1))
plot(neuralData.eta.eventWindow,mean(alignedPreds_value{1}(:,:,1),1))
plot(neuralData.eta.eventWindow, mean(kernelAnalysis.eta.alignedResps{1}(:,:,n),1))
plot(neuralData.eta.eventWindow,mean(neuralData.eta.alignedResps{1}(:,:,n),1))
ylim([-.01 .3])
xlim([-.5 2])
subplot(1,3,2)
hold on
line([0 0],[-.4 1.5],'Color',[.5 .5 .5],'LineStyle','--')
line([1 1],[-.4 1.5],'Color',[.5 .5 .5],'LineStyle','--')

plot(neuralData.eta.eventWindow,nanmean(alignedPreds_action{2}(:,:,1),1))
plot(neuralData.eta.eventWindow,nanmean(alignedPreds_choice{2}(:,:,1),1))
ylim([-.01 .3])
xlim([-.5 2])
subplot(1,3,3)
hold on
line([0 0],[-.4 1.5],'Color',[.5 .5 .5],'LineStyle','--')
line([1 1],[-.4 1.5],'Color',[.5 .5 .5],'LineStyle','--')

plot(neuralData.eta.eventWindow,nanmean(alignedPreds_feedback{3}(:,:,1),1))
ylim([-.01 .3])
xlim([-.5 2])
