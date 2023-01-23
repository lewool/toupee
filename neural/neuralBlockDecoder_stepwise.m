
%% Predict the block ID (1,0) from neural activity
% This takes a 1D vector of neural activity at a particular epoch and uses
% this in lasso logistic regression to try and predict the block that the 
% animal is currently in. 

% After determining the accuracy of the real model, it repeats the 
% training/fitting for a set number of pseudosessions and reports the accuracy 
% of those pseudomodels.

% Finally, this plots the real model accuracy against the distribution of
% pseudomodel accuracies to test whether it is significantly better (or
% worse?)


for m = 1:length(mouseList)
    
    mouseName = char(mouseList{m});
    expDate = char(expList{m}{1});
    expNum = expList{m}{2};
    hemisphere = hemList(m);
    [expRef, ~] = data.constructExpRef(mouseName,expDate,expNum);
    fprintf('Loading %s...',expRef)
    [behavioralData, expInfo, neuralData] = data.loadDataset(mouseName, expDate, expNum);
    [baselineResps, stimResps, pmovResps, movResps, rewResps, preCueResps] = getEpochResps(neuralData.eta);
    
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
    
    nt = length(behavioralData.eventTimes(1).daqTime);
    halfway = floor(nt/2);
    
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
        Y_pseudo(:,p) = b(1:nt);
    end
    
    %set up binomial Y outcome vector (block L vs R)
    if hemisphere > 0
        Y_all = double(expInfo.block.events.highRewardSideValues(1:nt) == -1)';
        Y_pseudo = Y_pseudo == -1;
    else
        Y_all = double(expInfo.block.events.highRewardSideValues(1:nt) == 1)';
        Y_pseudo = Y_pseudo == 1;
    end
    
   
    Y=Y_all;

    %retrieve the activity of each cell on every trial, at timebin t
    X = baselineResps;
    
    fprintf('fitting...')
    mdl =  stepwiseglm(X,Y,'constant','Distribution','binomial','upper','linear','nsteps',2);
    Y_hat = glmval([stats.intercept;b(finalmodel)],X(:,finalmodel),'logit');
    accuracy_true(m) = sum((Y_hat>0.5) == Y)/length(Y);

    prog = 0;
    fprintf(1,'now fitting perms: %3d%%\n',prog);

    for p = 1:1000
        clear b se pval finalmodel stats
        [b, se, pval, finalmodel, stats] = stepwisefit(X,Y_pseudo(:,p),'MaxIter',20);
        Y_hat = glmval([stats.intercept;b(finalmodel)],X(:,finalmodel),'logit');
        accuracy_pseudo(p,m) = sum((Y_hat>0.5) == Y_pseudo(:,p))/length(Y_pseudo(:,p));
        prog = ( 100*(p/1000) );
        fprintf(1,'\b\b\b\b%3.0f%%',prog);
    end
    fprintf('\n');

    clearvars -except mouseList expList hemList accuracy_true acc_pseudo accuracy_pseudo propWeights blockBias
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
    histogram(accuracy_pseudo(:,sp),linspace(0.5,1,41),'FaceColor',[.5 .5 .5])
    hold on
    if accuracy_true(sp) > prctile(accuracy_pseudo(:,sp),95)
        line([accuracy_true(sp) accuracy_true(sp)],[0 160],'LineWidth',2,'LineStyle','-','Color','r');
    else
        line([accuracy_true(sp) accuracy_true(sp)],[0 160],'LineWidth',2,'LineStyle','-','Color','k');
    end
    box off
    set(gca,'tickdir','out')
    title(expRef,'Interpreter','none')
    if sp == 31
        xlabel('Model accuracy')
        ylabel('\it n')
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