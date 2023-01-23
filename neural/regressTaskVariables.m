    mouseName = char(mouseList{m});
    expDate = char(expList{m}{1});
    expNum = expList{m}{2};
    hemisphere = hemList(m);
    [expRef, ~] = data.constructExpRef(mouseName,expDate,expNum);
    fprintf('fetching %s...\n',expRef)
    [behavioralData, expInfo, neuralData] = data.loadDataset(mouseName, expDate, expNum);
    [baselineResps, stimResps, pmovResps, movResps, rewResps, preCueResps] = getEpochResps(neuralData.eta);
    
    [~, whichTrials] = selectCondition(expInfo, getUniqueContrasts(expInfo), behavioralData, ...
        initTrialConditions('responseType','correct','movementTime','all','specificRTs',[0.5 Inf]));

    %% identify high and low trials based on the true block and stimulus
    nt = length(behavioralData.eventTimes(1).daqTime);
    if hemisphere < 0
        trialStimuli = expInfo.block.events.contrastValues(1:nt);
        trialCorrectChoice = expInfo.block.events.correctResponseValues(1:nt);
        trialBlockStructure = expInfo.block.events.highRewardSideValues(1:nt);
        trialActualChoice = expInfo.block.events.responseValues(1:nt);
    else
        trialStimuli = -expInfo.block.events.contrastValues(1:nt);
        trialCorrectChoice = -expInfo.block.events.correctResponseValues(1:nt);
        trialBlockStructure = -expInfo.block.events.highRewardSideValues(1:nt);
        trialActualChoice = -expInfo.block.events.responseValues(1:nt);
    end

    trialFeedback = expInfo.block.events.feedbackValues(1:nt);
    % assign the 0% stimuli as either 'left' or 'right' depending on the
    % preassigned correct choice (not the mouse's choice)
    trialStimuli(trialStimuli == 0) = eps;
    trialStimuli(abs(trialStimuli) < .05) = ...
        trialStimuli(abs(trialStimuli) < .05).* trialCorrectChoice(abs(trialStimuli) < .05);

    % low rewards are possible on sign-mismatched block and stimulus
    % high rewards are possible on sign-matched block and stimulus
    % 1 = high, 0 = low
    trueRewards(trialBlockStructure.*sign(trialStimuli) == -1) = 1;
    trueRewards(trialBlockStructure.*sign(trialStimuli) == 1) = 2;
    trueRewards(trialFeedback==0) = 0;
    
%% build X matrix of predictors (one value per trial)

X = [...
    trialStimuli' ...
    trialActualChoice' ...
    trialBlockStructure' ...
    trueRewards'] ;

X = X(whichTrials,:);

% gram-schmidt
[m,n] = size(X);
Q = zeros(m,n);
R = zeros(n,n);

for j = 1:n
    v = X(:,j);
    for i = 1:j-1
        R(i,j) = Q(:,i)'*X(:,j);
        v = v-R(i,j)*Q(:,i);
    end
    R(j,j) = norm(v);
    Q(:,j) = v/R(j,j);
end

%%
%set up some fitting options
options.alpha = 0;
options.nlambda = 20;
options.standardize = 'true';
family = 'gaussian';

try parpool();
catch
end
whichCell = 7; 
for t = 1:length(neuralData.eta.eventWindow)
    Y = neuralData.eta.alignedResps{1}(whichTrials,t,whichCell);    
    fit{t} = cvglmnet(X,Y,family,options,'deviance',10,[],true);
    bestIdx = find(fit{t}.lambda == fit{t}.lambda_min);
    weights(t,:) = fit{t}.glmnet_fit.beta(:,bestIdx);
end
figure;plot(neuralData.eta.eventWindow,weights,'LineWidth',2)
box off
set(gca,'tickdir','out');
xlim([-0.5 2]);
% ylim([-1 1]);

    
%% build X matrix of predictors (one value per time bin)

%stimulus
stimTrace = zeros(1,length(expInfo.Timeline.rawDAQTimestamps));
for t = 1:length(behavioralData.eventTimes(1).daqTime)
    [~, onset] = min(abs(expInfo.Timeline.rawDAQTimestamps - behavioralData.eventTimes(1).daqTime(t)));
    [~, offset] = min(abs(expInfo.Timeline.rawDAQTimestamps - behavioralData.eventTimes(3).daqTime(t)));
% 	onset = find(expInfo.Timeline.rawDAQTimestamps == behavioralData.eventTimes(1).daqTime(t));
%     offset = find(expInfo.Timeline.rawDAQTimestamps == behavioralData.eventTimes(3).daqTime(t));
    stimTrace(onset:offset) = trialStimuli(t);
end
stimTrace = interp1(expInfo.Timeline.rawDAQTimestamps,stimTrace,neuralData.respTimes);
stimTrace(isnan(stimTrace)) = 0;

%choice
choiceTrace = zeros(1,length(expInfo.Timeline.rawDAQTimestamps));
for t = 1:length(behavioralData.eventTimes(1).daqTime)
    [~, onset] = min(abs(expInfo.Timeline.rawDAQTimestamps - behavioralData.wheelMoves.epochs(5).onsetTimes(t) - 0.2));
    [~, offset] = min(abs(expInfo.Timeline.rawDAQTimestamps - behavioralData.wheelMoves.epochs(5).offsetTimes(t)));
    choiceTrace(onset:offset) = trialActualChoice(t);
end
choiceTrace = interp1(expInfo.Timeline.rawDAQTimestamps,choiceTrace,neuralData.respTimes);
choiceTrace(isnan(choiceTrace)) = 0;

%feedback
feedbackTrace = zeros(1,length(expInfo.Timeline.rawDAQTimestamps));
for t = 1:length(behavioralData.eventTimes(1).daqTime)
    [~, onset] = min(abs(expInfo.Timeline.rawDAQTimestamps - behavioralData.eventTimes(5).daqTime(t)));
    [~, offset] = min(abs(expInfo.Timeline.rawDAQTimestamps - (behavioralData.eventTimes(5).daqTime(t)+0.5)));
    feedbackTrace(onset:offset) = trueRewards(t);
end
feedbackTrace = interp1(expInfo.Timeline.rawDAQTimestamps,feedbackTrace,neuralData.respTimes);
feedbackTrace(isnan(feedbackTrace)) = 0;


%block
blockTrace = zeros(1,length(expInfo.Timeline.rawDAQTimestamps));
newBlocks = [0 find(diff(expInfo.block.events.highRewardSideValues)~=0)];
for b = 1:length(newBlocks)
    [~, onset] = min(abs(expInfo.Timeline.rawDAQTimestamps - behavioralData.eventTimes(5).daqTime(newBlocks(b)+1)));
    try
        [~, offset] = min(abs(expInfo.Timeline.rawDAQTimestamps - behavioralData.eventTimes(5).daqTime(newBlocks(b+1)+1)));
    catch
        [~, offset] = min(abs(expInfo.Timeline.rawDAQTimestamps - behavioralData.eventTimes(5).daqTime(end)));
    end
    blockTrace(onset:offset) = expInfo.block.events.highRewardSideValues(newBlocks(b)+1);
end
blockTrace = interp1(expInfo.Timeline.rawDAQTimestamps,blockTrace,neuralData.respTimes);
blockTrace(isnan(blockTrace)) = 0;
    
%% build predictor matrix

X = [...
    stimTrace' ...
    choiceTrace' ...
    feedbackTrace' ...
    blockTrace'] ;

% 
% % gram-schmidt
% [m,n] = size(X);
% Q = zeros(m,n);
% R = zeros(n,n);
% 
% for j = 1:n
%     v = X(:,j);
%     for i = 1:j-1
%         R(i,j) = Q(:,i)'*X(:,j);
%         v = v-R(i,j)*Q(:,i);
%     end
%     R(j,j) = norm(v);
%     Q(:,j) = v/R(j,j);
% end
ntp = length(X);
cvntp = round(ntp/3);
X_train = X(1:2*cvntp,:);
X_test = X(2*cvntp+1:end,:);
Y_train = Y(1:2*cvntp);
Y_test = Y(2*cvntp+1:end);

%%
%set up some fitting options
options.alpha = 0; %ridge
options.nlambda = 20;
options.standardize = 'true';
family = 'gaussian';

try parpool();
catch
end

whichCell = 4;
X = [...
    stimTrace' ...
    choiceTrace' ...
    feedbackTrace' ...
    blockTrace'] ;
Y = neuralData.cellResps(:,whichCell);
X_train = X(1:2*cvntp,:);
X_test = X(2*cvntp+1:end,:);
Y_train = Y(1:2*cvntp);
Y_test = Y(2*cvntp+1:end);
fit = cvglmnet(X_train,Y_train,family,options,'deviance',10,[],true);
Y_hat = cvglmnetPredict(fit,X_test,[],'response');
ssRes = sum((Y_test - Y_hat).^2);
ssTot = sum((Y_test - mean(Y_test)).^2);
r2 = 1 - ssRes/ssTot;

X = [...
    choiceTrace' ...
    feedbackTrace' ...
    blockTrace'] ;
Y = neuralData.cellResps(:,whichCell);
X_train = X(1:2*cvntp,:);
X_test = X(2*cvntp+1:end,:);
Y_train = Y(1:2*cvntp);
Y_test = Y(2*cvntp+1:end);
fit = cvglmnet(X_train,Y_train,family,options,'deviance',10,[],true);
Y_hat = cvglmnetPredict(fit,X_test,[],'response');
ssRes = sum((Y_test - Y_hat).^2);
ssTot = sum((Y_test - mean(Y_test)).^2);
r2 = 1 - ssRes/ssTot;



%%















whichCell = 1; 
for c = 1:size(neuralData.cellResps,2)
    Y = neuralData.cellResps(:,c);    
    fit = cvglmnet(X,Y,family,options,'deviance',10,[],true);
    bestIdx = find(fit.lambda == fit.lambda_min);
    weights(c,:) = fit.glmnet_fit.beta(:,bestIdx);
end

%which cells' fits are improved by adding stimulus?
X = [...
    choiceTrace' ...
    feedbackTrace' ...
    blockTrace'] ;

for c = 1:size(neuralData.cellResps,2)
    Y = neuralData.cellResps(:,c);    
    fit = cvglmnet(X,Y,family,options,'deviance',10,[],true);
    bestIdx = find(fit.lambda == fit.lambda_min);
    weights(c,:) = fit.glmnet_fit.beta(:,bestIdx);
end

%which cells' fits are improved by adding choice?
X = [...
    stimulusTrace' ...
    feedbackTrace' ...
    blockTrace'] ;

for c = 1:size(neuralData.cellResps,2)
    Y = neuralData.cellResps(:,c);    
    fit = cvglmnet(X,Y,family,options,'deviance',10,[],true);
    bestIdx = find(fit.lambda == fit.lambda_min);
    weights(c,:) = fit.glmnet_fit.beta(:,bestIdx);
end

%which cells' fits are improved by adding feedback?
X = [...
    stimulusTrace' ...
    choiceTrace' ...
    blockTrace'] ;

for c = 1:size(neuralData.cellResps,2)
    Y = neuralData.cellResps(:,c);    
    fit = cvglmnet(X,Y,family,options,'deviance',10,[],true);
    bestIdx = find(fit.lambda == fit.lambda_min);
    weights(c,:) = fit.glmnet_fit.beta(:,bestIdx);
end

%which cells' fits are improved by adding block?
X = [...
    stimulusTrace' ...
    choiceTrace' ...
    feedbackTrace'] ;

for c = 1:size(neuralData.cellResps,2)
    Y = neuralData.cellResps(:,c);    
    fit = cvglmnet(X,Y,family,options,'deviance',10,[],true);
    bestIdx = find(fit.lambda == fit.lambda_min);
    weights(c,:) = fit.glmnet_fit.beta(:,bestIdx);
end




%%
figure;
subplot(1,3,1)
hold on;
maxxy = max(max(abs(weights(:,1:3))));
xlim([-maxxy maxxy]);
ylim([-maxxy maxxy]);
line([-maxxy maxxy],[-0 0],'LineStyle','--','Color',[.5 .5 .5]);
line([-0 0],[-maxxy maxxy],'LineStyle','--','Color',[.5 .5 .5]);
scatter(weights(:,1),weights(:,2),20,'MarkerFaceColor','k','MarkerEdgeColor','none','MarkerFaceAlpha',0.5)
box off
set(gca,'tickdir','out');
axis square

subplot(1,3,2)
hold on;
maxxy = max(max(abs(weights(:,1:3))));
xlim([-maxxy maxxy]);
ylim([-maxxy maxxy]);
line([-maxxy maxxy],[-0 0],'LineStyle','--','Color',[.5 .5 .5]);
line([-0 0],[-maxxy maxxy],'LineStyle','--','Color',[.5 .5 .5]);
scatter(weights(:,2),weights(:,3),20,'MarkerFaceColor','k','MarkerEdgeColor','none','MarkerFaceAlpha',0.5)
box off
set(gca,'tickdir','out');
axis square

subplot(1,3,3)
hold on;
maxxy = max(max(abs(weights(:,1:3))));
xlim([-maxxy maxxy]);
ylim([-maxxy maxxy]);
line([-maxxy maxxy],[-0 0],'LineStyle','--','Color',[.5 .5 .5]);
line([-0 0],[-maxxy maxxy],'LineStyle','--','Color',[.5 .5 .5]);
scatter(weights(:,1),weights(:,4),20,'MarkerFaceColor','k','MarkerEdgeColor','none','MarkerFaceAlpha',0.5)
box off
set(gca,'tickdir','out');
axis square

 