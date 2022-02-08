for m = 1:6
    
    mouseName = char(mouseList{1});
    expDate = char(expList{m}{1});
    expNum = expList{m}{2};
    
    [expRef, ~] = data.constructExpRef(mouseName,expDate,expNum);
    fprintf('Loading %s...',expRef)
    [sessions(m).behavioralData, sessions(m).expInfo, sessions(m).neuralData] = data.loadDataset(mouseName, expDate, expNum);
end
   
for s = 1:length(sessions)
    nt(s) = length(sessions(s).behavioralData.eventTimes(1).daqTime);
end
mint = min(nt);
halfway = floor(mint/2);

%%
for s = 1:length(sessions)
    [baselineResps, stimResps, pmovResps, movResps, rewResps, preCueResps] = getEpochResps(sessions(s).neuralData.eta);
    X{s} = movResps(1:mint,:);
    
    %split into test/train for simple CV
    Xfirst{s} = X{s}(1:halfway,:);
    Xlast{s} = X{s}(halfway+1:end,:);
    
    hemisphere = 1;%hemList(s);
    
     %set up binomial Y outcome vector (block L vs R)
    if hemisphere > 0
        Y_all = double((sessions(s).expInfo.block.events.highRewardSideValues == -1)');
    else
        Y_all = double((sessions(s).expInfo.block.events.highRewardSideValues == 1)');
    end

     %set up binomial Y outcome vector (choice L vs R)
%     if hemisphere > 0
%         Y_all = double((sessions(s).behavioralData.wheelMoves.epochs(5).moveDir == -1)');
%     else
%         Y_all = double((sessions(s).behavioralData.wheelMoves.epochs(5).moveDir == 1)');
%     end
    
    Y(:,s)=Y_all(1:mint);
     
end

%split into test/train for simple CV
Yfirst = Y(1:halfway,:);
Ylast = Y(halfway+1:end,:);

%set up some lambdas
lambdas = [0, 1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1];
% lambdas = 1e-3;

for x = 1:length(sessions)
    for y = 1:length(sessions) 
        [B,FitInfo] = lassoglm(Xfirst{x},Yfirst(:,y),'binomial','Lambda',lambdas,'Alpha',1,'Standardize',false);
        Y_hat = glmval([FitInfo.Intercept;B],Xlast{x},'logit');
        allAcc = sum((Y_hat>0.5) == Ylast(:,y))/length(Ylast(:,y));
        [~, bestInd] = max(allAcc);
        accuracy(x,y) = allAcc(bestInd);
    end
end

trueAccuracy = sum(diag(accuracy))/length(sessions);

shuffles = perms(1:length(sessions));
shuffles = shuffles(1:end-1,:);

for s = 1:length(shuffles)
    shuffAcc = accuracy(:,shuffles(s,:));
    shuffAccuracy(s) = sum(diag(shuffAcc))/length(sessions);
end

%%
figure;
histogram(shuffAccuracy);
hold on
line([trueAccuracy trueAccuracy],[1 100],'Color','r','LineWidth',2)
box off
set(gca,'tickdir','out')
ylabel('No. sessions')
xlabel('Mean accuracy (pseudo vs true)')
% title('LEW031: block prediction, move epoch (six sessions)')
