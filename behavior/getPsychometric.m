function [cc, ppl, ppr, ppt, pcil, pcir, pcit] = getPsychometric(expInfo, behavioralData, trialList, cc)

block = expInfo.block;
et = behavioralData.eventTimes;
wm = behavioralData.wheelMoves;

allContrastValues = block.events.contrastValues(trialList);
allRewardSideValues = block.events.highRewardSideValues(trialList);
allRepeatValues = block.events.repeatNumValues(trialList);
allChoiceValues = wm.epochs(5).moveDir(trialList);
allFeedbackValues = block.events.feedbackValues(trialList);

ppl = nan(size(cc));
nnl = nan(size(cc));
ppr = nan(size(cc));
nnr = nan(size(cc));
ppt = nan(size(cc));
nnt = nan(size(cc));
alpha = 0.1; % 90% CIs - optional argument? 


% LEFT BLOCK
for cont = 1:length(cc)
    contIdxs = [];

    %determine which trials had contrast = cc(cont)
    contIdxs = (allRewardSideValues < 0) & (allContrastValues == cc(cont)) & (allRepeatValues == 1);

    %number of trials with contrast = cc(cont)
    nnl(cont) = sum(contIdxs); 
    try
        ppl(cont) = sum(allChoiceValues(contIdxs) == 1)/sum(contIdxs);
    catch
        if nnl(cont) == 0
             ppl(cont) = 0;
        end
    end 
end
    
% compute CIs
[~, pcil] = binofit(round(ppl.*nnl), nnl, alpha); % get CIs of the binomial distribution    
    
% RIGHT BLOCK
for cont = 1:length(cc)
    contIdxs = [];

    %determine which trials had contrast = cc(cont)
    contIdxs = (allRewardSideValues > 0) & (allContrastValues == cc(cont)) & (allRepeatValues == 1);

    %number of trials with contrast = cc(cont)
    nnr(cont) = sum(contIdxs); 
    try
        ppr(cont) = sum(allChoiceValues(contIdxs) == 1)/sum(contIdxs);
    catch
        if nnr(cont) == 0
             ppr(cont) = 0;
        end
    end 
end
    
% compute CIs
[~, pcir] = binofit(round(ppr.*nnr), nnr, alpha); % get CIs of the binomial distribution

% BOTH BLOCKS
for cont = 1:length(cc)
    contIdxs = [];

    %determine which trials had contrast = cc(cont)
    contIdxs = (allContrastValues == cc(cont)) & (allRepeatValues == 1);

    %number of trials with contrast = cc(cont)
    nnt(cont) = sum(contIdxs); 
    try
        ppt(cont) = sum(allChoiceValues(contIdxs) == 1)/sum(contIdxs);
    catch
        if nnt(cont) == 0
             ppt(cont) = 0;
        end
    end 
end
    
% compute CIs
[~, pcit] = binofit(round(ppt.*nnt), nnt, alpha); % get CIs of the binomial distribution

