function behaviorDecoder
% Predict choice, RT or wheel from task features

%%
significance=[];
Ytype = 'choice';
Xtype = 'block';
%

for m = 1:length(mouseList)
    
    mouseName = char(mouseList{m});
    expDate = char(expList{m}{1});
    expNum = expList{m}{2};
    hemisphere = hemList(m);
    [expRef, ~] = data.constructExpRef(mouseName,expDate,expNum);
    fprintf('Loading %s...\n',expRef)
    [behavioralData, expInfo] = data.loadBehavioralDataset(mouseName, expDate, expNum);
    
    nt = length(behavioralData.eventTimes(1).daqTime);
    
    %% get predictors
    
    trueStimuli = expInfo.block.events.contrastValues(1:nt);
    trialCorrectChoice = expInfo.block.events.correctResponseValues(1:nt);
    trueBlocks = expInfo.block.events.highRewardSideValues(1:nt);
    trueChoices = behavioralData.wheelMoves.epochs(5).moveDir;
%     trueChoices = expInfo.block.events.responseValues(1:nt);
    allFeedback = double(expInfo.block.events.feedbackValues(1:nt));
    feedbackIdx = find(trueStimuli == 0);
    trueFeedback = allFeedback;

    % assign the 0% stimuli as either 'left' or 'right' depending on the
    % preassigned correct choice (not the mouse's choice)
    trueStimuli(trueStimuli == 0) = eps;
    trueStimuli(abs(trueStimuli) < .05) = ...
        trueStimuli(abs(trueStimuli) < .05).* trialCorrectChoice(abs(trueStimuli) < .05);

    trueSide = trueStimuli > 0;

    % low rewards are possible on sign-mismatched block and stimulus
    % high rewards are possible on sign-matched block and stimulus
    % 1 = high, 0 = low
    trueValue(trueBlocks.*sign(trueStimuli) == -1) = 0;
    trueValue(trueBlocks.*sign(trueStimuli) == 1) = 1;

    % fetch RT and maxVel
    [trueImp, ~] = selectCondition(expInfo, getUniqueContrasts(expInfo), behavioralData, ...
        initTrialConditions('movementTime','early'));
    trueVel = abs(behavioralData.wheelMoves.epochs(5).peakVel(1:nt));
    trueRT = behavioralData.wheelMoves.epochs(5).onsetTimes(1:nt) - behavioralData.eventTimes(1).daqTime(1:nt);
    
    wt = 1:nt;
    wt(isnan(trueRT)) = [];
    
    trueStimuli = trueStimuli(wt);
    trueChoices = trueChoices(wt);
    trueBlocks = trueBlocks(wt);
    trueValue = trueValue(wt);
    trueVel = trueVel(wt);
    trueRT = trueRT(wt);
    trueImp = trueImp(wt);
    trueFeedback = trueFeedback(wt);
    
    ntt = length(trueChoices);

    %% generate pseudo X vectors
    np = 100;

    %generate randomly shuffled stims
    pseudoStimuli = nan(np,ntt);
    for p = 1:np
        pseudoStimuli(p,:) = randsample(trueStimuli,length(trueStimuli));
    end


    %generate linear-shifted choices
    trimLength = 20;
    for l = 1:trimLength*2+1
        ss = l;
        es = length(trueChoices) - (trimLength*2-l+1);
        shifts(l,:) = trueChoices(ss:es);
    end
    pseudoChoices = shifts([1:trimLength,trimLength+2:trimLength*2+1],:);
    
    %generate linear-shifted outcomes
    trimLength = 20;
    for l = 1:trimLength*2+1
        ss = l;
        es = length(trueFeedback) - (trimLength*2-l+1);
        fshifts(l,:) = trueFeedback(ss:es);
    end
    pseudoFeedback = fshifts([1:trimLength,trimLength+2:trimLength*2+1],:);
    

    %generate linear-shifted RTs
    trimLength = 20;
    for l = 1:trimLength*2+1
        ss = l;
        es = length(trueRT) - (trimLength*2-l+1);
        rt_shifts(l,:) = trueRT(ss:es);
    end
    pseudoRT = rt_shifts([1:trimLength,trimLength+2:trimLength*2+1],:);

    %generate linear-shifted maxVels
    trimLength = 20;
    for l = 1:trimLength*2+1
        ss = l;
        es = length(trueVel) - (trimLength*2-l+1);
        vel_shifts(l,:) = trueVel(ss:es);
    end
    pseudoVel = vel_shifts([1:trimLength,trimLength+2:trimLength*2+1],:);

    % generate pseudoblocks
    pseudoBlocks = nan(np,ntt);
    pseudoValue = nan(np,ntt);
    trimLength = 0;
    blockStart = 'fixed';
    for p = 1:np
        if strcmp(blockStart,'fixed')
            firstSide = trueBlocks(trimLength+1);
        elseif strcmp(blockStart,'rand')
            firstSide = randsample([-1, 1],1,true);
        end
        b=nan(1,ntt);
        switches = cumsum(125+randi(100,1,20));
        for s = 1:length(switches)
            if s == 1
                b((1+trimLength):switches(s)-1) = firstSide;
            elseif mod(s,2) == 1
                b((switches(s-1)+trimLength):switches(s)-1) = firstSide;
            elseif mod(s,2) == 0
                b((switches(s-1)+trimLength):switches(s)-1) = -firstSide;
            end
        end
        pseudoBlocks(p,:) = b(1:ntt);
    end

    % generate pseudo high/low trials based on the pseudoblocks + true stim
    for p = 1:np
        pseudoValue(p,pseudoBlocks(p,:).*sign(trueStimuli) == -1) = 0;
        pseudoValue(p,pseudoBlocks(p,:).*sign(trueStimuli) == 1) = 1;
    end
    %% set up Y outcome vector
    Y=[];
    if strcmp(Xtype,'choice') || strcmp(Xtype,'outcome')
        testTrials = [1:3:ntt-trimLength*2]+trimLength;
        trainTrials = [1:ntt-trimLength*2]+trimLength;
        trainTrials(ismember(trainTrials,testTrials)) = [];
        switch Ytype           
            case 'choice'
                Y.train = trueChoices(trainTrials)==1;
                Y.test = trueChoices(testTrials)==1;
                family = 'binomial';
            case 'imp'
                Y.train = trueImp(trainTrials)==1;
                Y.test = trueImp(testTrials)==1;
                family = 'binomial';
            case 'RT'
                Y.train = trueRT(trainTrials);
                Y.test = trueRT(testTrials);
                family = 'gaussian';
            case 'wheel'
                Y.train = trueVel(trainTrials);
                Y.test = trueVel(testTrials);
                family = 'gaussian';
        end
    else
        testTrials = 1:3:ntt;
        trainTrials = 1:ntt;
        trainTrials(ismember(trainTrials,testTrials)) = [];

        switch Ytype           
            case 'choice'
                Y.train = trueChoices(trainTrials)==1;
                Y.test = trueChoices(testTrials)==1;
                family = 'binomial';
            case 'imp'
                Y.train = trueImp(trainTrials)==1;
                Y.test = trueImp(testTrials)==1;
                family = 'binomial';
            case 'RT'
                Y.train = trueRT(trainTrials);
                Y.test = trueRT(testTrials);
                family = 'gaussian';
            case 'wheel'
                Y.train = trueVel(trainTrials);
                Y.test = trueVel(testTrials);
                family = 'gaussian';
        end
    end
           
    
    %% set up X predictor vector and pseudos
    
    switch Xtype
        case 'contrast'
            testTrials = 1:3:ntt;
            trainTrials = 1:ntt;
            trainTrials(ismember(trainTrials,testTrials)) = [];
            X.true.train = trueStimuli(trainTrials);
            X.true.test = trueStimuli(testTrials);
            X.pseudo.train = pseudoStimuli(:,trainTrials);
            X.pseudo.test = pseudoStimuli(:,testTrials);
        case 'block'
            testTrials = 1:3:ntt;
            trainTrials = 1:ntt;
            trainTrials(ismember(trainTrials,testTrials)) = [];
            X.true.train = trueBlocks(trainTrials);
            X.true.test = trueBlocks(testTrials);
            X.pseudo.train = pseudoBlocks(:,trainTrials);
            X.pseudo.test = pseudoBlocks(:,testTrials);
        case 'value'
            testTrials = 1:3:ntt;
            trainTrials = 1:ntt;
            trainTrials(ismember(trainTrials,testTrials)) = [];
            X.true.train = trueValue(trainTrials);
            X.true.test = trueValue(testTrials);
            X.pseudo.train = pseudoValue(:,trainTrials);
            X.pseudo.test = pseudoValue(:,testTrials);
        case 'outcome'
            trimLength=20;
            trueFeedback = trueFeedback(trimLength+1:ntt-trimLength);
            testTrials = 1:3:length(trueFeedback);
            trainTrials = 1:length(trueFeedback);
            trainTrials(ismember(trainTrials,testTrials)) = [];
            X.true.train = trueFeedback(trainTrials);
            X.true.test = trueFeedback(testTrials);
            X.pseudo.train = pseudoFeedback(:,trainTrials);
            X.pseudo.test = pseudoFeedback(:,testTrials);
        case 'choice'
            trimLength=20;
            trueChoices = trueChoices(trimLength+1:ntt-trimLength);
            testTrials = 1:3:length(trueChoices);
            trainTrials = 1:length(trueChoices);
            trainTrials(ismember(trainTrials,testTrials)) = [];
            X.true.train = trueChoices(trainTrials);
            X.true.test = trueChoices(testTrials);
            X.pseudo.train = pseudoChoices(:,trainTrials);
            X.pseudo.test = pseudoChoices(:,testTrials);
    end
    
    %% fitting 
    
    %set up some fitting options
    if strcmp(family,'binomial')
        dist = 'binomial';
        link = 'logit';
    else
        dist = 'normal';
        link = 'identity';
    end
    
    b = glmfit(X.true.train',Y.train',dist,'link',link);
    Y.true.hat = glmval(b,X.true.test,link);
    
    %goodness of fit
    %(gof = MSE for continuous (stimulus) or gof = LL for binomial)
    gof(m) = nan(1,1);
    if strcmp(family,'gaussian') %|| strcmp(family,'binomial')
        gof(m) = sum((Y.test' - (Y.true.hat)).^2)/length(Y.test);
    else
        term1 = Y.test'.*log(Y.true.hat);
        term2 = (1-Y.test').*log(1-Y.true.hat);
        ll = sum(term1 + term2);%/size(Y(t).true.hat,2);
        gof(m) = exp(ll/size(Y.true.hat',2));
    end
  
    
    %% fit pseudos
    np = size(X.pseudo.train,1);

    %fit
    for p = 1:np
        b = glmfit(X.pseudo.train(p,:)',Y.train',dist,'link',link);
        Y.pseudo.hat(p,:) = glmval(b,X.pseudo.test(p,:),link);
%         fit.pseudo{p} = cvglmnet(X.pseudo.train(p,:)',Y.train,family,options,'deviance',5,[],true);
%         Y.pseudo.hat(p,:) = cvglmnetPredict(fit.pseudo{p},X.pseudo.test(p,:)','lambda_min','response')';
    end

    % goodness of fit
    gof_pseudo(m,:) = nan(np,1);
    if strcmp(family,'gaussian')
        for p = 1:np
            gof_pseudo(m,p) = sum((Y.test - (Y.pseudo.hat(p,:))).^2)/length(Y.test);
        end
    else
        for p = 1:np
            term1 = Y.test.*log(Y.pseudo.hat(p,:));
            term2 = (1-Y.test).*log(1-Y.pseudo.hat(p,:));
            ll(p) = sum(term1 + term2);
            gof_pseudo(m,p) = exp(ll(p)/size(Y.pseudo.hat,2));
        end
    end

    
    %% test for significance
    
    significance(m) = gof(m) < prctile(gof_pseudo(m,:),2.5) || gof(m) > prctile(gof_pseudo(m,:),97.5);
    
    %% clear and restart for next session
    
    clearvars -except mouseList expList hemList gof gof_pseudo Xtype Ytype significance
end

%% plot true model accuracy against dist. of pseudosessions
subDims = ceil(sqrt(length(mouseList)));
figure;
set(gcf,'position',[32 80 2054 1250]);
plotIdx = [1:100,102:200];

for sp = 1:length(mouseList)
    mouseName = char(mouseList{sp});
    expDate = char(expList{sp}{1});
    expNum = expList{sp}{2};
    [expRef, ~] = data.constructExpRef(mouseName,expDate,expNum);
    subplot(3,4,sp)
    h = histogram(gof_pseudo(plotIdx,sp),'FaceColor',[.5 .5 .5]);
%     h = histogram(accuracy_pseudo(plotIdx,sp),'FaceColor',[.5 .5 .5]);
    maxy = max(h.Values);
    hold on
    if gof(sp) > prctile(gof_pseudo(:,sp),95)
        line([gof(sp) gof(sp)],[0 maxy],'LineWidth',2,'LineStyle','-','Color','r');
%         line([accuracy_true(sp) accuracy_true(sp)],[0 maxy],'LineWidth',2,'LineStyle','-','Color','r');
    else
        line([gof(sp) gof(sp)],[0 maxy],'LineWidth',2,'LineStyle','-','Color','k');
%         line([accuracy_true(sp) accuracy_true(sp)],[0 maxy],'LineWidth',2,'LineStyle','-','Color','k');
    end
    box off
    set(gca,'tickdir','out')
    title(expRef,'Interpreter','none')
    if sp < 31
        xlabel('Normalized likelihood')
        ylabel('\it n')
    end
%     xlim([.5 1])
    ylim([0 maxy*1.05])
%     if strcmp(fitFlag_true{sp},'smallest')
%         text(.51,maxy,strcat('poor fit (',num2str(100*1-(length(find(strcmp(fitFlag_pseudo(:,sp),'smallest')))/np)),'% good)'))
%     else
%         text(.51,maxy,strcat('good fit (',num2str(100*1-(length(find(strcmp(fitFlag_pseudo(:,sp),'smallest')))/np)),'% good)'))
%     end
end



%% plot model accuracy (as significance percentile) against block bias for each session
for m = 1:length(mouseList)
    
    mouseName = char(mouseList{m});
    expDate = char(expList{m}{1});
    expNum = expList{m}{2};
    hemisphere = hemList(m);
    [expRef, ~] = data.constructExpRef(mouseName,expDate,expNum);
    fprintf('Loading %s...',expRef)
    [behavioralData, expInfo, neuralData] = data.loadDataset(mouseName, expDate, expNum);
    
    if hemisphere < 0
        [~, iBlockTrials0] = selectCondition(expInfo, 0, behavioralData, ...
            initTrialConditions('highRewardSide','left'));
        [~, cBlockTrials0] = selectCondition(expInfo, 0, behavioralData, ...
            initTrialConditions('highRewardSide','right'));
        meanContraChoice = mean(expInfo.block.events.responseValues(cBlockTrials0) == 1);
        meanIpsiChoice = mean(expInfo.block.events.responseValues(iBlockTrials0) == 1);
    else
        [~, cBlockTrials0] = selectCondition(expInfo, 0, behavioralData, ...
            initTrialConditions('highRewardSide','left'));
        [~, iBlockTrials0] = selectCondition(expInfo, 0, behavioralData, ...
            initTrialConditions('highRewardSide','right'));
        meanContraChoice = mean(expInfo.block.events.responseValues(cBlockTrials0) == -1);
        meanIpsiChoice = mean(expInfo.block.events.responseValues(iBlockTrials0) == -1);
    end
    
    blockBias(m) = meanContraChoice - meanIpsiChoice;
    [~, pctVal(m)] = min(abs(prctile(accuracy_pseudo(:,m),[1:100]) - accuracy_true(m)));
end

%%
pctSig = [12 7 4 5 4 3 5 6 10 9 9 7 8 11 6 6 13 4 7 3 11 8 4 10 3 8 5 11 4 4 10 18 20 19]
xx = .6;
figure;
hold on
mouseColors = lines(5);
line([0 0],[0 100],'LineStyle','--','Color',[.5 .5 .5]);
line([-xx xx],[5 5],'LineStyle','--','Color',[.5 .5 .5]);
fill([-xx; -xx; xx; xx],[0; 2.5; 2.5; 0],'k','LineStyle','none','FaceAlpha',.2);

fill([-xx; -xx; xx; xx],[100; 97.5; 97.5; 100],'k','LineStyle','none','FaceAlpha',.2);

for m = 1:length(mouseList)
    if strcmp(mouseList{m},'LEW031')
        color = mouseColors(1,:);
    elseif strcmp(mouseList{m},'LEW032')
        color = mouseColors(2,:);
    elseif strcmp(mouseList{m},'LEW038')
        color = mouseColors(3,:);
    elseif strcmp(mouseList{m},'LEW046')
        color = mouseColors(4,:);
    elseif strcmp(mouseList{m},'LEW047')
        color = mouseColors(5,:);
    end

    scatter(blockBias(m),pctSig(m),60,'MarkerFaceColor',color,'MarkerEdgeColor','none');
    
end
xlim([-xx xx]);
box off    
set(gca,'tickdir','out') 
xlabel('Block bias')
ylabel('Model accuracy (significance percentile)')
axis square