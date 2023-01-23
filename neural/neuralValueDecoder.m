
%% Predict the block ID (1,0) from neural activity
% This takes a 1D vector of neural activity at a particular epoch (X) and uses
% this in lasso logistic regression to try and predict the value of each
% trial (which depends on stimulus and block)

% After determining the accuracy of the true model, it repeats the 
% fitting for np (typically = 1000) pseudogenerated blocks (same X, pseudo Y) 
% and reports the accuracy of those pseudomodels.

% Finally, this plots the true model accuracy against the distribution of
% pseudomodel accuracies to test whether it is significantly better than
% chance

% 'good fit' means the lamba selected by cvglmnet for the true model 
% falls safely within the range and not at the edge of the range of preselected lambdas. 
% 'P% good' reports the % of pseudo models with a safe lambda (defined as above)

parpool();

%%
for m = 1:length(mouseList)
    
    mouseName = char(mouseList{m});
    expDate = char(expList{m}{1});
    expNum = expList{m}{2};
    hemisphere = hemList(m);
    [expRef, ~] = data.constructExpRef(mouseName,expDate,expNum);
    fprintf('Loading %s...',expRef)
    [behavioralData, expInfo, neuralData] = data.loadDataset(mouseName, expDate, expNum);
    [baselineResps, stimResps, pmovResps, movResps, rewResps, preCueResps] = getEpochResps(neuralData.eta);
    
    nt = length(behavioralData.eventTimes(1).daqTime);
    
    %% identify high and low trials based on the true block and stimulus
    trialStimuli = expInfo.block.events.contrastValues(1:nt);
    trialCorrectChoice = expInfo.block.events.correctResponseValues(1:nt);
    trialBlockStructure = expInfo.block.events.highRewardSideValues(1:nt);
    trialActualChoice = expInfo.block.events.responseValues(1:nt);

    % assign the 0% stimuli as either 'left' or 'right' depending on the
    % preassigned correct choice (not the mouse's choice)
    trialStimuli(trialStimuli == 0) = eps;
    trialStimuli(abs(trialStimuli) < .05) = ...
        trialStimuli(abs(trialStimuli) < .05).* trialCorrectChoice(abs(trialStimuli) < .05);

    % low rewards are possible on sign-mismatched block and stimulus
    % high rewards are possible on sign-matched block and stimulus
    % 1 = high, 0 = low
    trueRewards(trialBlockStructure.*sign(trialStimuli) == -1) = 0;
    trueRewards(trialBlockStructure.*sign(trialStimuli) == 1) = 1;
    
    %% generate pseudoblocks
    trimLength = 0;
    np = 1000;
    blockStart = 'fixed';
    for p = 1:np
        if strcmp(blockStart,'fixed')
            firstSide = expInfo.block.paramsValues(1).firstHighSide;
        elseif strcmp(blockStart,'rand')
            firstSide = randsample([-1, 1],1,true);
        end
        b=nan(1,nt);
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
        y_pseudo(:,p) = b(1:nt);
    end
    
    % assign new high/low trials based on the pseudoblocks + true stim
    for p = 1:np
        pseudoRewards(y_pseudo(:,p).*sign(trialStimuli)' == -1,p) = 0;
        pseudoRewards(y_pseudo(:,p).*sign(trialStimuli)' == 1,p) = 1;
    end
    
    %%

    %set up predictor matrix
    whichEpoch = 'stim';
    if strcmp(whichEpoch,'stim')
        [~, whichTrials] = selectCondition(expInfo, getUniqueContrasts(expInfo), behavioralData, ...
        initTrialConditions('responseType','all','movementTime','all','specificRTs',[0.5 Inf]));
        X = stimResps(whichTrials,:);
        Y = trueRewards(whichTrials);
        Y_pseudo = pseudoRewards(whichTrials,:);
    elseif strcmp(whichEpoch,'baseline')
        [~, whichTrials] = selectCondition(expInfo, getUniqueContrasts(expInfo), behavioralData, ...
        initTrialConditions('responseType','all'));
        X = baselineResps(whichTrials,:);
        Y = trueRewards(whichTrials);
        Y_pseudo = pseudoRewards(whichTrials,:);
    elseif strcmp(whichEpoch,'move')
        [~, whichTrials] = selectCondition(expInfo, getUniqueContrasts(expInfo), behavioralData, ...
        initTrialConditions('responseType','all'));
        X = movResps(whichTrials,:);
        Y = trueRewards(whichTrials);
        Y_pseudo = pseudoRewards(whichTrials,:);
        [i, ~]= ind2sub(size(X),find(isnan(X)));
        X(unique(i),:) = [];
        Y(unique(i)) = [];
        Y_pseudo(unique(i),:) = [];
    end
    
    %% fitting 
    
    %set up some fitting options
    options.alpha = 0.9;
    options.nlambda = 20;
    options.standardize = 'false';
    family = 'binomial';
    

    %fit true block
    fprintf('fitting...')
    fitTrue = cvglmnet(X,Y,family, options,'deviance',5,[],true);
    bestLambda_true(m) = fitTrue.lambda_min;
    bestIdx = find(fitTrue.glmnet_fit.lambda == fitTrue.lambda_min);
    if bestLambda_true(m) == min(fitTrue.lambda)
        fitFlag_true{m} = 'smallest';
    elseif bestLambda_true(m) == max(fitTrue.lambda)
        fitFlag_true{m} = 'largest';
    else
        fitFlag_true{m} = 'inRange';
    end
    Y_hat = cvglmnetPredict(fitTrue, X,[],'response');
    accuracy_true(m) = sum((Y_hat>0.5) == Y')/length(Y);
    binDev_true(m) = fitTrue.cvm(bestIdx);
    dev_true(m) = fitTrue.glmnet_fit.dev(bestIdx);
    
    %% fitting 
    
    %fit pseudoblocks
    prog = 0;
    fprintf(1,'now fitting perms: %3d%%\n',prog);

    for p = 1:np
        fitPseudo = cvglmnet(X,Y_pseudo(:,p),family, options,'deviance',5,[],true);
        bestLambda_pseudo(p,m) = fitPseudo.lambda_min;
        bestIdx(p) = find(fitPseudo.glmnet_fit.lambda == fitPseudo.lambda_min);
        if bestLambda_pseudo(p,m) == min(fitPseudo.lambda)
            fitFlag_pseudo{p,m} = 'smallest';
        elseif bestLambda_pseudo(p,m) == max(fitPseudo.lambda)
            fitFlag_pseudo{p,m} = 'largest';
        else
            fitFlag_pseudo{p,m} = 'inRange';
        end
        Y_hat = cvglmnetPredict(fitPseudo, X,[],'response');
        accuracy_pseudo(p,m) = sum((Y_hat>0.5) == Y_pseudo(:,p))/length(Y_pseudo(:,p));
        binDev_pseudo(p,m) = fitPseudo.cvm(bestIdx(p));
        dev_pseudo(p,m) = fitPseudo.glmnet_fit.dev(bestIdx(p));
        prog = ( 100*(p/np) );
        fprintf(1,'\b\b\b\b%3.0f%%',prog);
    end
    fprintf('\n');

    clearvars -except mouseList expList hemList accuracy_true acc_pseudo accuracy_pseudo binDev_pseudo binDev_true propWeights blockBias fitFlag_true fitFlag_pseudo
end

%% plot true model accuracy against dist. of pseudosessions
subDims = ceil(sqrt(length(mouseList)));
figure;
set(gcf,'position',[32 80 2054 1250]);
for sp = 1:length(mouseList)
    mouseName = char(mouseList{sp});
    expDate = char(expList{sp}{1});
    expNum = expList{sp}{2};
    [expRef, ~] = data.constructExpRef(mouseName,expDate,expNum);
    subplot(subDims,subDims,sp)
    h = histogram(accuracy_pseudo(:,sp),linspace(0.5,1,41),'FaceColor',[.5 .5 .5]);
    maxy = max(h.Values);
    hold on
    if accuracy_true(sp) > prctile(accuracy_pseudo(:,sp),95)
        line([accuracy_true(sp) accuracy_true(sp)],[0 maxy],'LineWidth',2,'LineStyle','-','Color','r');
    else
        line([accuracy_true(sp) accuracy_true(sp)],[0 maxy],'LineWidth',2,'LineStyle','-','Color','k');
    end
    box off
    set(gca,'tickdir','out')
    title(expRef,'Interpreter','none')
    if sp == 31
        xlabel('Model accuracy')
        ylabel('\it n')
    end
    xlim([.5 1])
    ylim([0 maxy*1.05])
    if strcmp(fitFlag_true{sp},'smallest')
        text(.51,maxy,strcat('poor fit (',num2str(100*1-(length(find(strcmp(fitFlag_pseudo(:,sp),'smallest')))/np)),'% good)'))
    else
        text(.51,maxy,strcat('good fit (',num2str(100*1-(length(find(strcmp(fitFlag_pseudo(:,sp),'smallest')))/np)),'% good)'))
    end
end

%%

subDims = ceil(sqrt(length(mouseList)));
figure;
set(gcf,'position',[32 80 2054 1250]);
for sp = 1:length(mouseList)
    mouseName = char(mouseList{sp});
    expDate = char(expList{sp}{1});
    expNum = expList{sp}{2};
    [expRef, ~] = data.constructExpRef(mouseName,expDate,expNum);
    subplot(subDims,subDims,sp)
     h = histogram(binDev_pseudo(:,sp),41,'FaceColor',[.5 .5 .5]);
    maxy = max(h.Values);
    hold on
    if binDev_true(sp) < prctile(binDev_pseudo(:,sp),5)
        line([binDev_true(sp) binDev_true(sp)],[0 maxy],'LineWidth',2,'LineStyle','-','Color','r');
    else
        line([binDev_true(sp) binDev_true(sp)],[0 maxy],'LineWidth',2,'LineStyle','-','Color','k');
    end
    box off
    set(gca,'tickdir','out')
    title(expRef,'Interpreter','none')
    if sp == 31
        xlabel('CV binomial deviance')
        ylabel('\it n')
    end
    ylim([0 maxy*1.05])
    if strcmp(fitFlag_true{sp},'smallest')
        text(.51,maxy,strcat('poor fit (',num2str(100*1-(length(find(strcmp(fitFlag_pseudo(:,sp),'smallest')))/np)),'% good)'))
    else
        text(.51,maxy,strcat('good fit (',num2str(100*1-(length(find(strcmp(fitFlag_pseudo(:,sp),'smallest')))/np)),'% good)'))
    end
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