function kernelHoldout(behavioralData, expInfo, neuralData, features)

features = {'stimulus' 'value' 'block' 'action' 'choice' 'outcome'};

%% setup  
% design matrix
[predictors, windows] = getPredictors(...
    expInfo, ...
    behavioralData, ...
    features, ...
    0.1);

% set up regression matrices
X = makeToeplitz(neuralData.respTimes, predictors, windows);
Y = neuralData.cellResps;

%% RRR prep
% factorize to find basis functions for RRR
[B,~,~] = kFactorize(X, Y);

%% fit full kernels for each neuron
%for each neuron, determine the optimal rank of XB to use to find its
%kernels, then compute how much variance is explained by this model
cellList = [768 661 520 623 50];
cellList = [390 1827 520 623 1617 1392 363 153 781];

nFold = 5;
nComp = size(X,2);

% initialize variables
maxEV = zeros(size(cellList,2),1);
optR = zeros(size(cellList,2),1);
kernelFunctions = zeros(nComp,size(cellList,2));

for c = 1:length(cellList)
    
    %full-rank fit (nComp = [])
    [~, weightsK, intercept, explVarAll] = rrFit(X,Y(:,cellList(c)),B,[],nFold,0);

    %determine the rank with the highest EV
    [maxEV(c,:), optR(c,:)] = max(explVarAll);

    %generate the kernels based on this rank
    kernelFunctions(:,c) = B(:,1:optR(c,:))*weightsK(:,1:optR(c,:))';
    intercepts(c) = intercept;
end

predictedActivity = X*kernelFunctions + intercepts;

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
    alignedPreds{ev} = zeros(nt, length(eventWindow), size(predictedActivity,2));

    %grab cell responses associated with the event time windows 
    %(size is nTrials x windowLength x nCells)
    for iCell = 1:size(predictedActivity,2)        
        alignedPreds{ev}(:,:,iCell) = interp1(neuralData.respTimes,predictedActivity(:,iCell),periEventTimes,'previous');
    end
end


%% fit reduced kernels for each neuron

for f = 1:length(features)
    
    clear kernelFunctions_without;
    
    %determine the 'feature of interest' to test
    foi = features(f);

    %remove the foi from the feature set
    otherFeatures = features;
    otherFeatures(f) = [];

    %generate a predictor matrix without the foi
    [predictors_without, windows_without] = getPredictors(...
        expInfo, ...
        behavioralData, ...
        otherFeatures, ...
        0.1);

    %predictor matrices
    X_without = makeToeplitz(neuralData.respTimes, predictors_without, windows_without);

    %factorize
    [B_without,~,~] = kFactorize(X_without, Y);
    
    for c = 1:length(cellList)
        
        %fit without
        [~, weightsK_without, intercept, ~] = rrFit(X_without, Y(:,cellList(c)), B_without, optR(c), nFold, 0);

        %generate the kernels based on this rank
        kernelFunctions_without(:,c) = B_without(:,1:optR(c,:))*weightsK_without(:,1:optR(c,:))';
        intercepts(c) = intercept;

    end
    
    predictedActivity_without{f} = X_without*kernelFunctions_without + intercepts;

end

% align kernel responses

Fs = 0.1;
timeBefore = 2;
timeAfter = 2;
events = {'stimulusOnTimes' 'firstMoveTimes' 'feedbackTimes' 'interactiveOnTimes'};

nt = numel(expInfo.block.events.endTrialTimes);
for f = 1:length(features)
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
        alignedPreds_wo{f,ev} = zeros(nt, length(eventWindow), size(predictedActivity_without{f},2));

        %grab cell responses associated with the event time windows 
        %(size is nTrials x windowLength x nCells)
        for iCell = 1:size(predictedActivity_without{f},2)        
            alignedPreds_wo{f,ev}(:,:,iCell) = interp1(neuralData.respTimes,predictedActivity_without{f}(:,iCell),periEventTimes,'previous');
        end
    end
end

%% plot full model

cell = 6;
cellno = cellList(cell);

Fs = 0.1;
whichAligned = kernelAnalysis.eta.alignedResps; linst = '-';
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
lb = logical((trueBlock == -1) .* whichTrials);
rb = logical((trueBlock == 1) .* whichTrials);

[co, ~] = selectCondition(expInfo, contrasts, et, initTrialConditions('responseType','correct','movementTime','late'));
[io, ~] = selectCondition(expInfo, contrasts, et, initTrialConditions('responseType','incorrect','movementTime','late'));

co = logical(co);
io = logical(io);

% fetch responses and model output
% r_hls = neuralData.eta.alignedResps{1}(hls,:,cellno);
% r_lls = neuralData.eta.alignedResps{1}(lls,:,cellno);
% r_zs = neuralData.eta.alignedResps{1}(zs,:,cellno);
% r_lrs = neuralData.eta.alignedResps{1}(lrs,:,cellno);
% r_hrs = neuralData.eta.alignedResps{1}(hrs,:,cellno);
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

% m_hls = whichAligned{1}(hls,:,cellno);
% m_lls = whichAligned{1}(lls,:,cellno);
% m_zs = whichAligned{1}(zs,:,cellno);
% m_lrs = whichAligned{1}(lrs,:,cellno);
% m_hrs = whichAligned{1}(hrs,:,cellno);
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


% rm_hls = alignedKernelResps_wo{1}(hls,:,cellno);
% rm_lls = alignedKernelResps_wo{1}(lls,:,cellno);
% rm_zs = alignedKernelResps_wo{1}(zs,:,cellno);
% rm_lrs = alignedKernelResps_wo{1}(lrs,:,cellno);
% rm_hrs = alignedKernelResps_wo{1}(hrs,:,cellno);
rm_ls = alignedPreds_wo{1,1}(ls,:,cell);
rm_rs = alignedPreds_wo{1,1}(rs,:,cell);
rm_lc = alignedPreds_wo{5,2}(lc,:,cell);
rm_rc = alignedPreds_wo{5,2}(rc,:,cell);
rm_hv = alignedPreds_wo{2,1}(hv,:,cell);
rm_lv = alignedPreds_wo{2,1}(lv,:,cell);
rm_co = alignedPreds_wo{6,3}(co,:,cell);
rm_io = alignedPreds_wo{6,3}(io,:,cell);
rm_lb = alignedPreds_wo{3,1}(lb,:,cell);
rm_rb = alignedPreds_wo{3,1}(rb,:,cell);


colors = {[0 .4 1], [1 0 0]; ...
          [.5 0 1], [1 0 .5]; ...
          [0 .5 0], [.7 .1 0]; ...
          [.1 .7 .1], [1 .6 0]; ...
          [.8 0 1], [0 .6 .1]; ...
          };

whichR = {r_ls, r_rs; ...
          r_lc, r_rc; ...
          r_co, r_io; ...
          r_lb, r_rb; ...
          r_hv, r_lv; ...
          };
      
whichM = {m_ls, m_rs; ...
          m_lc, m_rc; ...
          m_co, m_io; ...
          m_lb, m_rb; ...
          m_hv, m_lv; ...
          };

whichRM = {rm_ls, rm_rs; ...
          rm_lc, rm_rc; ...
          rm_co, rm_io; ...
          rm_lb, rm_rb; ...
          rm_hv, rm_lv; ...
          };

whichETA = {'stimulus' 'movement' 'feedback' 'stimulus' 'stimulus'};
figure;   
set(gcf,'position',[1165 950 980 160])
nlim = 1.1*max(abs([nanmean(r_ls) nanmean(r_rs) nanmean(r_hv) nanmean(r_lv) nanmean(r_rc) nanmean(r_lc) nanmean(r_co) nanmean(r_io)]));
% nlim = .3;
for s = 1:length(whichR)
    subplot(1,5,s)
    hold on
    plotSignal(neuralData.eta.eventWindow,nanmean(whichR{s,1},1),nanmean(whichR{s,1},1)+nanstd(whichR{s,1})/sqrt(size(whichR{s,1},1)),nanmean(whichR{s,1},1)-nanstd(whichR{s,1})/sqrt(size(whichR{s,1},1)),colors{s,1},'-');
    plotSignal(neuralData.eta.eventWindow,nanmean(whichR{s,2},1),nanmean(whichR{s,2},1)+nanstd(whichR{s,2})/sqrt(size(whichR{s,2},1)),nanmean(whichR{s,2},1)-nanstd(whichR{s,2})/sqrt(size(whichR{s,2},1)),colors{s,2},'-');
    plot(neuralData.eta.eventWindow,nanmean(whichM{s,1},1),'Color',colors{s,1},'LineWidth',2,'LineStyle','-')
    plot(neuralData.eta.eventWindow,nanmean(whichM{s,2},1),'Color',colors{s,2},'LineWidth',2,'LineStyle','-')
    plot(neuralData.eta.eventWindow,nanmean(whichRM{s,1},1),'Color',colors{s,1},'LineWidth',2,'LineStyle',':')
    plot(neuralData.eta.eventWindow,nanmean(whichRM{s,2},1),'Color',colors{s,2},'LineWidth',2,'LineStyle',':')
    line([0 0],[-nlim nlim],'Color',[.5 .5 .5],'LineStyle','--')
    xlim([-.5 .9])
    ylim([-.05 nlim])
    xlabel(strcat({'Time from '},whichETA{s},{' (s)'}))
    prettyPlot(gca,.04)
end


