function [correct_highSide, correct_lowSide, incorrect] = getHighLoRewardTrials(expInfo, behavioralData)

contrasts = getUniqueContrasts(expInfo);
[~, mL_bL_correct] = selectCondition(expInfo, contrasts, behavioralData, initTrialConditions(...
    'movementTime','late',...
    'responseType','correct',...
    'highRewardSide','left',...
    'movementDir','cw'));

[~, mL_bL_incorrect] = selectCondition(expInfo, contrasts, behavioralData, initTrialConditions(...
    'movementTime','late',...
    'responseType','incorrect',...
    'highRewardSide','left',...
    'movementDir','cw'));

[~, mR_bL_correct] = selectCondition(expInfo, contrasts, behavioralData, initTrialConditions(...
    'movementTime','late',...
    'responseType','correct',...
    'highRewardSide','left',...
    'movementDir','ccw'));

[~, mR_bL_incorrect] = selectCondition(expInfo, contrasts, behavioralData, initTrialConditions(...
    'movementTime','late',...
    'responseType','incorrect',...
    'highRewardSide','left',...
    'movementDir','ccw'));

[~, mL_bR_correct] = selectCondition(expInfo, contrasts, behavioralData, initTrialConditions(...
    'movementTime','late',...
    'responseType','correct',...
    'highRewardSide','right',...
    'movementDir','cw'));

[~, mL_bR_incorrect] = selectCondition(expInfo, contrasts, behavioralData, initTrialConditions(...
    'movementTime','late',...
    'responseType','incorrect',...
    'highRewardSide','right',...
    'movementDir','cw'));

[~, mR_bR_correct] = selectCondition(expInfo, contrasts, behavioralData, initTrialConditions(...
    'movementTime','late',...
    'responseType','correct',...
    'highRewardSide','right',...
    'movementDir','ccw'));

[~, mR_bR_incorrect] = selectCondition(expInfo, contrasts, behavioralData, initTrialConditions(...
    'movementTime','late',...
    'responseType','incorrect',...
    'highRewardSide','right',...
    'movementDir','ccw'));

%%

correct_highSide = [mL_bL_correct mR_bR_correct];
correct_lowSide = [mR_bL_correct mL_bR_correct];
incorrect_highSide = [mL_bL_incorrect mR_bR_incorrect];
incorrect_lowSide = [mR_bL_incorrect mL_bR_incorrect];
incorrect = [incorrect_highSide incorrect_lowSide];

%%

% ep = movResps;
% corHi = nanmean(ep(correct_highSide,:),1);
% corLo = nanmean(ep(correct_lowSide,:),1);
% incorHi = nanmean(ep(incorrect_highSide,:),1);
% incorLo = nanmean(ep(incorrect_lowSide,:),1);
% incor = nanmean(ep(incorrect,:),1);


%%

% labels = [ones(1,1277) ones(1,1277)*2 ones(1,1277)*3];
% resps = [corHi corLo incor];
% figure;
% beeswarm(labels',resps','sort_style','hex','dot_size',.02,'overlay_style','sd');
