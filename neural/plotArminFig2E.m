for c = 1:5
    [~, hconds1] = selectCondition(expInfo, contrasts(c), behavioralData, initTrialConditions('movementTime','all','responseType','correct','highRewardSide','left'));
    [~, hconds2] = selectCondition(expInfo, -contrasts(c), behavioralData, initTrialConditions('movementTime','all','responseType','correct','highRewardSide','right'));
    highConds{c,1} = [hconds1 hconds2];
    
    [~, lconds1] = selectCondition(expInfo, contrasts(c), behavioralData, initTrialConditions('movementTime','all','responseType','correct','highRewardSide','right'));
    [~, lconds2] = selectCondition(expInfo, -contrasts(c), behavioralData, initTrialConditions('movementTime','all','responseType','correct','highRewardSide','left'));
    lowConds{c,1} = [lconds1 lconds2];
    
    [~, iconds1] = selectCondition(expInfo, contrasts(c), behavioralData, initTrialConditions('movementTime','all','responseType','incorrect','highRewardSide','left'));
    [~, iconds2] = selectCondition(expInfo, -contrasts(c), behavioralData, initTrialConditions('movementTime','all','responseType','incorrect','highRewardSide','right'));
    incConds{c,1} = [iconds1 iconds2];
end

%%

% plotCells = getWhichCells('all',neuralData);

figure;
grays = [0 .2 .4 .6 .8];
subplot(1,3,1)
hold on
for pc = 1:length(highConds)
    ph = plot(eventWindow,(nanmean(nanmean(bigAlignedResps{2}(highConds{pc},:,:),3))));
    set(ph,'LineWidth',1,'Color', grays(pc)*ones(1,3));
end
xlim([-1 1]);
ylim([.01 .04]);
line([0 0],[0.0 .07],'LineStyle',':','Color',[.5 .5 .5])

subplot(1,3,2)
hold on
for pc = 1:length(lowConds)
    pl = plot(eventWindow,(nanmean(nanmean(bigAlignedResps{2}(lowConds{pc},:,:),3))));
    set(pl,'LineWidth',1,'Color', grays(pc)*ones(1,3));
end
xlim([-1 1]);
ylim([.01 .04]);
line([0 0],[0.0 .07],'LineStyle',':','Color',[.5 .5 .5])

subplot(1,3,3)
hold on
for pc = 1:length(incConds)
    pl = plot(eventWindow,(nanmean(nanmean(bigAlignedResps{2}(incConds{pc},:,:),3))));
    set(pl,'LineWidth',1,'Color', grays(pc)*ones(1,3));
end
xlim([-1 1]);
ylim([.01 .04]);
line([0 0],[0.0 .07],'LineStyle',':','Color',[.5 .5 .5])

set(gcf,'renderer','Painters');

