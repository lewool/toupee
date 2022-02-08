%% check for short sessions and truncate the rest
sessionIDs = 1:length(expInfo);

for ex = 1:length(expInfo)
    numTrials(ex) = length(expInfo(ex).block.events.endTrialValues);
end

minLength = median(numTrials)-2*mad(numTrials);
sessionIDs = sessionIDs(numTrials > minLength);
numTrials = numTrials(numTrials > minLength);
truncLength = min(numTrials);

%% generate the permutation matrix for pseudo-pairs

numPerms = 5000;
for iP = 1:numPerms
    onePerm = zeros(length(sessionIDs)-1,length(sessionIDs));
    for iS = 1:length(sessionIDs)
        pseudoPairs = sessionIDs(sessionIDs ~= sessionIDs(iS));
        onePerm(:,iS) = pseudoPairs(randperm(length(pseudoPairs),length(pseudoPairs)))';
    end
    allPerms{iP,:} = onePerm;
end

% each row randomly assigns a pseudo-pair for each session in 'sessionIDs'
permList = cat(1,allPerms{:});

%% compute metastatistic across true pairs

for iS = 1:length(sessionIDs)
    X = mean(eyeData(sessionIDs(iS)).eta.alignedFace{1}(1:truncLength,91:101,2),2);
    Y = behavioralData(sessionIDs(iS)).wheelMoves.epochs(5).onsetTimes(1:truncLength)' ...
        - behavioralData(sessionIDs(iS)).eventTimes(1).daqTime(1:truncLength)';
    X(find(isnan(Y))) = [];
    Y(find(isnan(Y))) = [];
    maxRT = nanmedian(Y) + 3*mad(Y);
    X(Y>maxRT) = [];
    Y(Y>maxRT) = [];
    
    mdl = fitlm(X,Y);
    r(iS) = mdl.Rsquared.Ordinary;
%     r(iS) = min(unique(corrcoef(X, Y)));
end

metaR_true = sum(r);

for iP = 1:numPerms
    for iS = 1:length(sessionIDs)
        X = mean(eyeData(sessionIDs(iS)).eta.alignedFace{1}(1:truncLength,91:101,2),2);
        Y = behavioralData(permList(iP,iS)).wheelMoves.epochs(5).onsetTimes(1:truncLength)' ...
            - behavioralData(permList(iP,iS)).eventTimes(1).daqTime(1:truncLength)';
        X(find(isnan(Y))) = [];
        Y(find(isnan(Y))) = [];
        maxRT = nanmedian(Y) + 3*mad(Y);
        X(Y>maxRT) = [];
        Y(Y>maxRT) = [];
        
        mdl = fitlm(X,Y);
        r(iS) = mdl.Rsquared.Ordinary;
%         r(iS) = min(unique(corrcoef(X, Y)));
    end
    
    metaR_pseudo(iP) = sum(r);
end

