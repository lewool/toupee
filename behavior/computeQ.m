function [Q, Q_pseudo] = computeQ(behavioralData, expInfo)

nt = numel(expInfo.block.events.endTrialTimes);
contrasts = getUniqueContrasts(expInfo);
[~, ppl, ppr, ~, ~] = getPsychometric(expInfo, behavioralData, 1:nt, contrasts);

pp = mean([ppl; ppr]);

trueStimuli = expInfo.block.events.contrastValues(1:nt);
trialCorrectChoice = expInfo.block.events.correctResponseValues(1:nt);
trueBlocks = expInfo.block.events.highRewardSideValues(1:nt);
trueChoices = behavioralData.wheelMoves.epochs(5).moveDir(1:nt);

firstMoveFeedback = (trialCorrectChoice .* trueChoices == 1);


% assign the 0% stimuli as either 'left' or 'right' depending on the
% preassigned correct choice (not the mouse's choice)
trueStimuli(trueStimuli == 0) = eps;
trueStimuli(abs(trueStimuli) < .05) = ...
    trueStimuli(abs(trueStimuli) < .05).* trialCorrectChoice(abs(trueStimuli) < .05);

sidedStimuli(trueStimuli < 0,1) = abs(trueStimuli(trueStimuli < 0));
sidedStimuli(trueStimuli > 0,2) = abs(trueStimuli(trueStimuli > 0));
sidedStimuli(sidedStimuli == 0) = eps;

uniqueStim = unique(trueStimuli);
for c = 1:length(uniqueStim)
    hitRate(c) = mean(firstMoveFeedback(trueStimuli == uniqueStim(c)));
end

uniqueAbsStim = fliplr(unique(abs(trueStimuli)));
hitMatrix = reshape(hitRate,[5,2]);
hitMatrix(:,2) = flipud(hitMatrix(:,2));
        
% for c = 1:length(uniqueAbsStim)
%     hitSeries(sidedStimuli(:,1) == uniqueAbsStim(c),1) = hitMatrix(c,1);
%     hitSeries(sidedStimuli(:,2) == uniqueAbsStim(c),2) = hitMatrix(c,2);
% end

trueStimuli = expInfo.block.events.contrastValues(1:nt);
for c = 1:length(contrasts)
    hitSeries(trueStimuli == contrasts(c),1) = 1 - pp(c);
    hitSeries(trueStimuli == contrasts(c),2) = pp(c);
end


% low rewards are possible on sign-mismatched block and stimulus
% high rewards are possible on sign-matched block and stimulus
% 1 = high, 0 = low
sidedValue(trueBlocks == -1,1) = 2;
sidedValue(trueBlocks == 1,1) = 1;
sidedValue(trueBlocks == -1,2) = 1;
sidedValue(trueBlocks == 1,2) = 2;

Qseries = sidedValue.*hitSeries;

Q.Qc = bsxfun(@max,Qseries(:,1),Qseries(:,2));
Q.Qtotal = sum(Qseries,2);
Q.Qrel = Qseries(:,2) - Qseries(:,1);

%% generate pseudo-Q series for stats

np = 1000;
pseudoBlocks = nan(np,nt);
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
            b(1:switches(s)-1) = firstSide;
        elseif mod(s,2) == 1
            b(switches(s-1):switches(s)-1) = firstSide;
        elseif mod(s,2) == 0
            b(switches(s-1):switches(s)-1) = -firstSide;
        end
    end
    pseudoBlocks(p,:) = b(1:nt);
end

for p = 1:np
    sidedValue_pseudo(pseudoBlocks(p,:) == -1,1,p) = 2;
    sidedValue_pseudo(pseudoBlocks(p,:) == 1,1,p) = 1;
    sidedValue_pseudo(pseudoBlocks(p,:) == -1,2,p) = 1;
    sidedValue_pseudo(pseudoBlocks(p,:) == 1,2,p) = 2;
end

Qseries_pseudo = sidedValue_pseudo.*hitSeries;

Q_pseudo.Qc = squeeze(bsxfun(@max,Qseries_pseudo(:,1,:),Qseries_pseudo(:,2,:)));
Q_pseudo.Qtotal = squeeze(sum(Qseries_pseudo,2));
Q_pseudo.Qrel = squeeze(Qseries_pseudo(:,2,:) - Qseries_pseudo(:,1,:));



%% plot
% 
% figure;
% 
% subplot(3,1,1)
% plot(Qc);
% plot(smooth(Qc,'moving',40));
% prettyPlot(gca)
% xlim([1 nt])
% title('Q_{chosen} = max(Q_c,Q_i)');
% 
% subplot(3,1,2)
% plot(Qtotal);
% plot(smooth(Qtotal,'moving',40));
% prettyPlot(gca)
% xlim([1 nt])
% title('Q_{total} = Q_c + Q_i');
% 
% subplot(3,1,3)
% plot(Qrel);
% plot(smooth(Qrel,'moving',40));
% prettyPlot(gca)
% xlim([1 nt])
% title('Q_{relative} = Q_c – Q_i');
% 
