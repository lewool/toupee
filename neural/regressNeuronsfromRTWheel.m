
for m = 1:length(mouseList)
    
    %load data
    mouseName = char(mouseList{m});
    expDate = char(expList{m}{1});
    expNum = expList{m}{2};
    [expRef, ~] = data.constructExpRef(mouseName,expDate,expNum);
    fprintf('Loading %s...\n',expRef)
    [behavioralData, expInfo, neuralData] = data.loadDataset(mouseName, expDate, expNum);
    
    contrasts = getUniqueContrasts(expInfo);
    nt = length(behavioralData.eventTimes(1).daqTime);

    [impTrials, ~] = selectCondition(expInfo, contrasts, behavioralData, ...
        initTrialConditions('movementTime','early'));
    
    trueChoices = behavioralData.wheelMoves.epochs(5).moveDir;
    maxVels = behavioralData.wheelMoves.epochs(5).peakVel;
    RTs = behavioralData.wheelMoves.epochs(5).onsetTimes - behavioralData.eventTimes(1).daqTime;
    trials = intersect(...
        find((RTs < prctile(RTs,97.5)) & (RTs > prctile(RTs,2.5))),...
        find(~isnan(trueChoices)));


    %% get neural responses

    [baselineResps, stimResps, pmovResps, movResps, rewResps, preCueResps] = getEpochResps(neuralData.eta);


    %% regress RT
    X = [ones(length(RTs(trials)),1) RTs(trials)'];
    Y = baselineResps(trials,:);

    trimLength = 20;
    for c = 1:length(Y)
        clear B R stats
        trueIdx = trimLength+1;
        Y_trunc = Y(trimLength+1:end - trimLength,c);
        for l = 1:trimLength*2+1
            ss = l;
            es = length(X) - (trimLength*2-l+1);
            [B(l,:), ~, ~, ~, stats(l,:)] = regress(Y_trunc,X(ss:es,:));
        end
        trueB_RT{m}(c,:) = B(trueIdx,:);
        trueR2_RT{m}(c) = stats(trueIdx,1);
        pseudoR2_RT = stats([1:trimLength,trimLength+2:trimLength*2+1],1);
        [~,h] = ranksum(trueR2_RT{m}(c),pseudoR2_RT,'method','exact','tail','right');
        hTest{m,1}(c,1) = double(h);
    end
    
    
    %% regress max vel

    X = [ones(length(maxVels(trials)),1) abs(maxVels(trials))'];
    Y = movResps(trials,:);

    trimLength = 20;
    for c = 1:length(Y)
        clear B R stats
        trueIdx = trimLength+1;
        Y_trunc = Y(trimLength+1:end - trimLength,c);
        for l = 1:trimLength*2+1
            ss = l;
            es = length(X) - (trimLength*2-l+1);
            [B(l,:), ~, ~, ~, stats(l,:)] = regress(Y_trunc,X(ss:es,:));
        end
        
        trueB_wheel{m}(c,:) = B(trueIdx,:);
        trueR2_wheel{m}(c) = stats(trueIdx,1);
        pseudoR2_wheel = stats([1:trimLength,trimLength+2:trimLength*2+1],1);
        [~,h] = ranksum(trueR2_wheel{m}(c),pseudoR2_wheel,'method','exact','tail','right');
        hTest{m,1}(c,2) = double(h);
    end    
    
    %% clean up
    clearvars -except mouseList expList hTest trueB_wheel trueB_RT trueR2_wheel trueR2_RT
end

%%
trueB_RT_all = [];
trueB_wheel_all = [];
for m = 1:length(mouseList)
    trueB_RT_all = cat(1,trueB_RT_all,trueB_RT{m}(:,2));
    trueB_wheel_all = cat(1,trueB_wheel_all,trueB_wheel{m}(:,2));
end

%% WORKBENCH

%% regress instantaneous wheel velocity
% 
%     velTrace = cat(2, behavioralData.wheelMoves.traces.pos{:});
%     timeTrace = cat(2, behavioralData.wheelMoves.traces.time{:});
% 
%     velTrace_int = interp1(timeTrace, velTrace, neuralData.respTimes);
% 
%     Y = neuralData.cellResps(~isnan(velTrace_int),:);
%     X = [ones(length(velTrace_int(~isnan(velTrace_int))),1) abs(velTrace_int(~isnan(velTrace_int))')];
% 
%     trimLength = 20;
%     for c = 1:length(Y)
%         clear B R stats
%         trueIdx = trimLength+1;
%         Y_trunc = Y(trimLength+1:end - trimLength,c);
%         for l = 1:trimLength*2+1
%             ss = l;
%             es = length(X) - (trimLength*2-l+1);
%             [B(l,:), ~, R(l,:), ~, stats(l,:)] = regress(Y_trunc,X(ss:es,:));
%         end
% 
%         trueVal = stats(trueIdx,1);
%         pseudoVals = stats([1:trimLength,trimLength+2:trimLength*2+1],1);
%         [~,h] = ranksum(trueVal,pseudoVals,'method','exact','tail','right');
%         hTest(c,3) = double(h);
%     end
