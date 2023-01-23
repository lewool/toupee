
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
    [expRef, ~] = data.constructExpRef(mouseName,expDate,expNum);
    fprintf('Loading %s...',expRef)
    [behavioralData, expInfo, neuralData] = data.loadDataset(mouseName, expDate, expNum);
    expInfo.hemisphere = hemList(m);
    [baselineResps, stimResps, pmovResps, movResps, rewResps, preCueResps] = getEpochResps(neuralData.eta);
    bigArminMatrix = arminizeBehavior(expInfo);
    [FitResult, NLLData] = ArminQCModel(expInfo, bigArminMatrix);
    QC_true = FitResult.QC;
    if expInfo.hemisphere > 0
        QContra_true = FitResult.QL;
        QIpsi_true = FitResult.QR;
    else
        QContra_true = FitResult.QR;
        QIpsi_true = FitResult.QL;
    end
    
    nt = length(behavioralData.eventTimes(1).daqTime);
    
    trimLength = 0;
    np = 1000;
    blockStart = 'fixed';
    for p = 1:np
        if strcmp(blockStart,'fixed')
            if expInfo.block.paramsValues(1).firstHighSide == 1
                firstSide = 2;
                secondSide = 1;
            elseif expInfo.block.paramsValues(1).firstHighSide == -1
                firstSide = 1;
                secondSide = 2;
            end
        elseif strcmp(blockStart,'rand')
            firstSide = randsample([1, 2],1,true);
        end
        b=nan(1,nt);
        switches = cumsum(125+randi(100,1,20));
        for s = 1:length(switches)
            if s == 1
                b((1+trimLength):switches(s)-1) = firstSide;
            elseif mod(s,2) == 1
                b((switches(s-1)+trimLength):switches(s)-1) = firstSide;
            elseif mod(s,2) == 0
                b((switches(s-1)+trimLength):switches(s)-1) = secondSide;
            end
        end
        Y_pseudo(:,p) = b(1:nt);
    end
    
    for p = 1:10%:np
        bigArminMatrix(:,8) = Y_pseudo(:,p);
        [FitResult, NLLData] = ArminQCModel(expInfo, bigArminMatrix);
        QC_pseudo(:,p) = FitResult.QC;
        if expInfo.hemisphere > 0
            QContra_pseudo(:,p) = FitResult.QL;
            QIpsi_pseudo(:,p) = FitResult.QR;
        else
            QContra_pseudo(:,p) = FitResult.QR;
            QIpsi_pseudo(:,p) = FitResult.QL;
        end
    end
    
%     for p = 1:1000
%         b=zeros(1,nt);
%         switches = cumsum(125+randi(100,1,10));
%         for s = 1:length(switches)
%             if s == 1
%                 b(1:switches(s)-1) = -1;
%             elseif mod(s,2) == 1
%                 b(switches(s-1):switches(s)-1) = -1;
%             elseif mod(s,2) == 0
%                 b(switches(s-1):switches(s)-1) = 1;
%             end
%         end
%         b = b(1:nt);
%         flip = randsample([-1, 1],1,true);
%         b = flip*b;
%         Y_pseudo(:,p) = b;
%     end
    
    %set up binomial Y outcome vector (block L vs R)
    if hemisphere > 0
        Y_all = double(expInfo.block.events.highRewardSideValues(1:nt) == -1)';
        Y_pseudo = Y_pseudo == -1;
    else
        Y_all = double(expInfo.block.events.highRewardSideValues(1:nt) == 1)';
        Y_pseudo = Y_pseudo == 1;
    end
   
    Y = QC_true;

     %set up some lambdas
%     lambdas = [0, 1e-10, 1e-8, 1e-6, 1e-4, 1e-2, 1];
    lambdas = 1e-2;

    %retrieve the activity of each cell on every trial, at timebin t
%     X = squeeze(neuralData.eta.alignedResps{2}(lateTrials,timerange(t),:));
    X = stimResps;
        
    fprintf('fitting...')
    Y = QContra_true;
    [B,FitInfo] = lassoglm(X,Y,'normal','Lambda',lambdas);
    Y_hat = glmval([FitInfo.Intercept;B],X,'identity');
    propWeights(m,:) = sum(B>0)/size(X,1);
    sset = sum((Y - Y_hat).^2);
    [~, bestIdx] = max(sset);
    SSE_true(m) = sset(bestIdx);

    lambda_match = lambdas(bestIdx);
    prog = 0;
    fprintf(1,'now fitting perms: %3d%%\n',prog);

    for p = 1:10%00
        Y = QContra_pseudo(:,p);
        [B,FitInfo] = lassoglm(X,Y,'normal','Lambda',lambdas);
        Y_hat = glmval([FitInfo.Intercept;B],X,'identity');
        ssep = sum((Y - Y_hat).^2);
        [~, bestIdx] = max(ssep);
        SSE_pseudo(p,m) = ssep(bestIdx);
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
    if bestSSE_true(sp) > prctile(accuracy_pseudo(:,sp),95)
        line([bestSSE_true(sp) bestSSE_true(sp)],[0 160],'LineWidth',2,'LineStyle','-','Color','r');
    else
        line([bestSSE_true(sp) bestSSE_true(sp)],[0 160],'LineWidth',2,'LineStyle','-','Color','k');
    end
    box off
    set(gca,'tickdir','out')
    title(expRef,'Interpreter','none')
    if sp == 31
        xlabel('Model accuracy')
        ylabel('\it n')
    end
    xlim([.5 1])
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
    [~, pctVal(m)] = min(abs(prctile(accuracy_pseudo(:,m),[1:100]) - bestSSE_true(m)));
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

[B,FitInfo] = lassoglm(X,Y,'normal','Lambda',lambdas);
    Y_hat = glmval([FitInfo.Intercept;B],X,'identity');