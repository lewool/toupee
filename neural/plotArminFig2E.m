for iX = 1:length(expInfo)
    
    contrasts = getUniqueContrasts(expInfo(iX));
    for c = 1:5
        [~, hconds1] = selectCondition(expInfo(iX), contrasts(c), behavioralData(iX), initTrialConditions('movementTime','all','responseType','correct','highRewardSide','left'));
        [~, hconds2] = selectCondition(expInfo(iX), -contrasts(c), behavioralData(iX), initTrialConditions('movementTime','all','responseType','correct','highRewardSide','right'));
        highConds{c,1} = [hconds1 hconds2];

        [~, lconds1] = selectCondition(expInfo(iX), contrasts(c), behavioralData(iX), initTrialConditions('movementTime','all','responseType','correct','highRewardSide','right'));
        [~, lconds2] = selectCondition(expInfo(iX), -contrasts(c), behavioralData(iX), initTrialConditions('movementTime','all','responseType','correct','highRewardSide','left'));
        lowConds{c,1} = [lconds1 lconds2];

        [~, iconds1] = selectCondition(expInfo(iX), contrasts(c), behavioralData(iX), initTrialConditions('movementTime','all','responseType','incorrect','highRewardSide','left'));
        [~, iconds2] = selectCondition(expInfo(iX), -contrasts(c), behavioralData(iX), initTrialConditions('movementTime','all','responseType','incorrect','highRewardSide','right'));
        incConds{c,1} = [iconds1 iconds2];
    end

    for c = 1:5
        highRew(c,:,iX) = nanmean(nanmean(neuralData(1).eta.alignedResps{2}(highConds{c},:,:),3));
        lowRew(c,:,iX) = nanmean(nanmean(neuralData(1).eta.alignedResps{2}(lowConds{c},:,:),3));
        noRew(c,:,iX) = nanmean(nanmean(neuralData(1).eta.alignedResps{2}(incConds{c},:,:),3));
    end
end
%%

mean_highRew = flipud(nanmean(highRew,3));
mean_lowRew = flipud(nanmean(lowRew,3));
mean_noRew = flipud(nanmean(noRew,3));

% plotCells = getWhichCells('all',neuralData);

figure;
grays = [.8 .6 .4 .2 .0];
for iC = 1:size(highRew,1)
    subplot(1,3,1)
    hold on
    line([0 0],[-.06 .14],'LineStyle','--','LineWidth',1,'Color',[.5 .5 .5])
    ph = plot(eventWindow,smooth(mean_highRew(iC,:)));
    set(ph,'LineWidth',2,'Color', grays(iC)*ones(1,3));
    xlim([-1 1])
    axis off
    
    subplot(1,3,2)
    hold on
    line([0 0],[-.06 .14],'LineStyle','--','LineWidth',1,'Color',[.5 .5 .5])
    ph = plot(eventWindow,smooth(mean_lowRew(iC,:)));
    set(ph,'LineWidth',2,'Color', grays(iC)*ones(1,3));
    xlim([-1 1])
    axis off
    
    subplot(1,3,3)
    hold on
    line([0 0],[-.06 .14],'LineStyle','--','LineWidth',1,'Color',[.5 .5 .5])
    pn = plot(eventWindow,smooth(mean_noRew(iC,:)));
    set(pn,'LineWidth',2,'Color', grays(iC)*ones(1,3));
    xlim([-1 1])
    axis off
end
%%
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

