
%% Predict RT from neural activity
% This takes a 1D vector of neural activity at a particular epoch (X) and uses
% this in lasso logistic regression to try and predict whether the upcoming 
% animal move is impulsive or patient (Y). 

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
    
    %% get RTs
    
    allRTs = behavioralData.wheelMoves.epochs(5).onsetTimes(1:nt) - behavioralData.eventTimes(1).daqTime(1:nt);
    filterTrials = find((allRTs > prctile(allRTs,2.5)) & (allRTs < prctile(allRTs,97.5)));
    filteredRTs = allRTs(filterTrials);
    
    %% set up Y outcome vector and pseudos
    
    % generate linear-shifted RTs
    trimLength = 20;
    apertureIdx = trimLength+1:length(filteredRTs)-trimLength;
    for l = 1:trimLength*2+1
        ss = l;
        es = length(filteredRTs) - (trimLength*2-l+1);
        shifts(l,:) = filteredRTs(ss:es);
    end
    trueRT = shifts(trimLength+1,:);
    pseudoRT = shifts([1:trimLength,trimLength+2:trimLength*2+1],:);
    
    %% set up predictor matrix and choose which observations to use for X and Y
    
    allIdx = filterTrials(apertureIdx);
    cvp = cvpartition(size(allIdx,2),'HoldOut',round(length(allIdx)/3));
    testIdx = find(cvp.test);
    trainIdx = find(~cvp.test);
    
    X.all = baselineResps(allIdx,:);
    X.train = baselineResps(trainIdx,:);
    X.test = baselineResps(testIdx,:);
    
    Y.true.all = trueRT;
    Y.true.train = trueRT(trainIdx);
    Y.true.test = trueRT(testIdx);
    Y.pseudo.all = pseudoRT;
    Y.pseudo.train = pseudoRT(:,trainIdx);
    Y.pseudo.test = pseudoRT(:,testIdx);
    
    %% fitting 

    %set up some fitting options
    options.alpha = 0;
    options.nlambda = 20;
    options.standardize = 'false';
    family = 'gaussian';
    
    fitTrue = cvglmnet(X.train,Y.true.train',family, options,'deviance',5,[],true);
    Y_hat = cvglmnetPredict(fitTrue, X.test,'lambda_min','link');
    
    gof(m) = sum((Y.true.test - Y_hat').^2)/length(Y.true.test);    
    
    %% fit pseudoblocks
    prog = 0;
    fprintf(1,'now fitting perms: %3d%%\n',prog);
    
    np = size(Y.pseudo.train,1);
    cvTag = 0;
    for p = 1:np
        if cvTag == 1
            nFold = 5;
            cvp = cvpartition(size(X.train,1),'Kfold',nFold);
            T = 1:size(X.train,1);
        
            for i = 1:nFold

                %split timepoints for cv
                trainIdx = ismember(T, T(cvp.training(i)));
                testIdx = ismember(T, T(cvp.test(i)));

                Y_train = Y.pseudo.train(p,trainIdx);

                % fit X to Y using the training set
                fitPseudo = cvglmnet(X(trainIdx,:),Y_train,family, options,'deviance',5,[],true);

                % evaluate at the desired lambda value
                Y_hat_pseudo(testIdx,p) = cvglmnetPredict(fitPseudo, X(testIdx,:),'lambda_min','response');

            end
        else
            fitPseudo = cvglmnet(X.train,Y.pseudo.train(p,:),family, options,'deviance',5,[],true);
            Y_hat_pseudo(p,:) = cvglmnetPredict(fitPseudo, X.test,'lambda_min','response');
        end
        
        prog = ( 100*(p/np) );
        fprintf(1,'\b\b\b\b%3.0f%%',prog);
        
        gof_pseudo(p,m) = sum((Y.pseudo.test(p,:) - Y_hat_pseudo(p,:)).^2)/length(Y.pseudo.test(p,:));
    
    end
    fprintf('\n');
    
    %% clear and restart for next session
    
    clearvars -except mouseList expList hemList gof gof_pseudo accuracy_true accuracy_pseudo
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