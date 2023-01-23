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
    
    %% get early vs late trials
    [impIdx, impTrials] = selectCondition(expInfo, getUniqueContrasts(expInfo), behavioralData, ...
        initTrialConditions('responseType','all','movementTime','early','specificRTs',[0.1 .8]));
    [~, patTrials] = selectCondition(expInfo, getUniqueContrasts(expInfo), behavioralData, ...
        initTrialConditions('responseType','all','movementTime','late','specificRTs',[0.8 3]));
    
    whichTrials = false(1,nt);
    whichTrials(impTrials) = true;
    whichTrials(patTrials) = true;
    whichTrialsList = sort([impTrials patTrials]);
    
    %% set up binomial Y outcome vector (impulsive vs nonimpulsive)
    
    yFull = impIdx(whichTrials)';
    
      %% set up predictor matrix and choose which observations to use for X and Y
    X = baselineResps(whichTrials,:);
    Y = yFull;
    
    %% fitting 
    
    %set up some fitting options
    options.alpha = 1;
    options.nlambda = 20;
    options.standardize = 'false';
    family = 'binomial';
    
    cvTag = 1;
    
    if cvTag == 1
        nFold = 5;
        cvp = cvpartition(size(X,1),'Kfold',nFold);
        T = 1:size(X,1);

        for i = 1:nFold

            %split timepoints for cv
            trainIdx = ismember(T, T(cvp.training(i)));
            testIdx = ismember(T, T(cvp.test(i)));

            Y_train = Y(trainIdx,:);

            % fit X to Y using the training set
            fitTrue = cvglmnet(X(trainIdx,:),Y_train,family, options,'deviance',5,[],true);

            % evaluate at the desired lambda value
            Y_hat(testIdx) = cvglmnetPredict(fitTrue, X(testIdx,:),'lambda_min','response');

        end
    else
        fitTrue = cvglmnet(X,Y,family, options,'deviance',5,[],true);
        Y_hat = cvglmnetPredict(fitTrue, X,'lambda_min','response');
    end

    term1 = Y.*log(Y_hat)';
    term2 = (1-Y).*log(1-Y_hat)';
    ll = sum(term1 + term2);%/size(Y(t).true.hat,2);
    gof(m) = exp(ll/size(Y,1));
    accuracy_true(m) = sum((Y_hat>0.5) == Y')/length(Y);
    
    predictions(m).Y = Y;
    predictions(m).Y_hat = Y_hat';
    predictions(m).nl = gof(m);
    predictions(m).accuracy = accuracy_true(m);
    
    %% compute 'impulsivity index' for each cell
    
    numer = mean(baselineResps(impTrials,:)) - mean(baselineResps(patTrials,:));
    denom = mean(baselineResps(impTrials,:)) + mean(baselineResps(patTrials,:));
    impSelIdx = numer./denom;
    [~, sortIdx] = sort(impSelIdx,'descend');
    
    %%
    
    figure;
    set(gcf,'position',[198 506 1734 1120])
    subplot(52,1,[3:48])
    imagesc(X(:,sortIdx)')
%     plot(X(:,sortIdx)')
    colormap(flipud(gray));
    caxis([-0 .4])
    axis off
    
    subplot(52,1,[2])
    imagesc(Y')
%     colormap(flipud(gray));
    caxis([-0 1])
    axis off
    
    subplot(52,1,[1])
    plot(smooth(Y_hat))
    colormap(flipud(gray));
    xlim([1 length(Y_hat)])
    ylim([.4 .6])
    prettyPlot(gca)
    axis off
    
    subplot(52,1,[49])
    imagesc(Y')
%     colormap(flipud(gray));
    caxis([-0 1])
    axis off
    
    subplot(52,1,[50])
    imagesc(Y_hat)
%     colormap(flipud(gray));
    caxis([0 1])
    prettyPlot(gca)
    set(gca, 'YTickLabels', {''})
    xlabel('Trials')

%     printfig(gcf,char(strcat(expRef,{' imp pseudocolor'})))
%     close all
%%
    cm = customColormap([1 .5 0],[0 .6 .6]);
    trialSnapshot = 1:400;
    selectNeurons = 100;
    plotNeuronsA = sortIdx(1:selectNeurons);
    plotNeuronsB = sortIdx(end-selectNeurons+1:end);
    plotNeurons = [plotNeuronsA plotNeuronsB];
    
    psthA = mean(X(trialSnapshot,plotNeuronsA),2)';
    psthB = mean(X(trialSnapshot,plotNeuronsB),2)';
    figure;
    set(gcf,'position',[198 506 934 820]);
    subplot(42,1,[1 2])
%     plot(smooth(Y_hat(trialSnapshot),3));
    plot(smooth(psthA,2))
    xlim([1 length(trialSnapshot)])
    prettyPlot(gca)
    axis off
    
    ax2 = subplot(42,1,[3 4]);
    imagesc(Y(trialSnapshot)');
    colormap(ax2,flipud(cm));
    caxis([-0 1])
    axis off
    
    ax3 = subplot(42,1,[5:20]);
    imagesc(X(trialSnapshot,plotNeuronsA)');
    colormap(ax3,flipud(gray));
    caxis([-0 .4])
    axis off
    
    ax4 = subplot(42,1,[23:38]);
    imagesc(X(trialSnapshot,plotNeuronsB)');
    colormap(ax4,flipud(gray));
    caxis([-0 .4])
    axis off
    
    subplot(42,1,[39 40])
    imagesc(Y(trialSnapshot)');
    caxis([-0 1])
    axis off
    
    subplot(42,1,[41 42])
%     plot(smooth(Y_hat(trialSnapshot),3));
    plot(smooth(psthB,2))
    xlim([1 length(trialSnapshot)])
    prettyPlot(gca)
    axis off
    %%
    
    cm = customColormap([1 .5 0],[0 .6 .6]);
    trialSnapshot = 100:200;
    selectNeurons = 250;
    cax = [0 .3];
    plotNeuronsA = sortIdx(1:selectNeurons);
    plotNeuronsB = sortIdx(end-selectNeurons+1:end);
    plotNeurons = [plotNeuronsA plotNeuronsB];
    
    psthA = mean(X(trialSnapshot,plotNeuronsA),2)';
    psthB = mean(X(trialSnapshot,plotNeuronsB),2)';
    figure;
    set(gcf,'position',[198 506 934 820]);
    
    
    
    ax3 = subplot(42,1,[1:16]);
    imagesc(X(trialSnapshot,plotNeuronsB)');
    colormap(ax3,flipud(gray));
    caxis(cax)
    axis off
    
    subplot(42,1,[17 18])
%     plot(smooth(Y_hat(trialSnapshot),3));
    plot(smooth(psthB,1),'k','LineWidth',1)
    xlim([1 length(trialSnapshot)])
    prettyPlot(gca)
    axis off
    
    ax2 = subplot(42,1,[19 20]);
    imagesc(Y(trialSnapshot)');
    colormap(ax2,flipud(cm));
    caxis(cax)
    axis off
    
    subplot(42,1,[21 22])
%     plot(smooth(Y_hat(trialSnapshot),3));
    plot(smooth(zscore(psthA),1),'k','LineWidth',1)
    xlim([1 length(trialSnapshot)])
    prettyPlot(gca)
    axis off
    
    ax4 = subplot(42,1,[23:38]);
    imagesc(X(trialSnapshot,plotNeuronsA)');
    colormap(ax4,flipud(gray));
    caxis(cax)
    axis off
    
    %%
    cm = customColormap([1 .5 0],[0 .6 .6],2);
    figure;
    scatter(psthB,psthA,30,cm(Y(trialSnapshot)'+1,:),'filled')
prettyPlot(gca)
xlim([0 .18])
ylim([0 .18])
axis square
line([0 .18],[0 .18],'Color',[.5 .5 .5],'LineStyle','--');
xlabel('Positively modulated')
ylabel('Negatively modulated')
 
    
    

     %% clear and restart for next session
    
    clearvars -except mouseList expList hemList m predictions
end