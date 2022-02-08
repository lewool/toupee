for m = 1:length(mouseList)
    
    mouseName = char(mouseList{m});
    expDate = char(expList{m}{1});
    expNum = expList{m}{2};
    [expRef, ~] = data.constructExpRef(mouseName,expDate,expNum);
    fprintf('fetching %s...\n',expRef)
    [behavioralData, expInfo, neuralData] = data.loadDataset(mouseName, expDate, expNum);
    [neuralData] = getSignificantActivity(expInfo, behavioralData, neuralData,0);
    eventWindow = neuralData.eta.eventWindow;

    patience = {'early' 'late'};
    etas = [1 2];

    for p = 1:length(patience)
        for e = 1:length(etas)

            contrasts = getUniqueContrasts(expInfo);
            for c = 1:5
                if c == 5
                    [~, hconds01] = selectCondition(expInfo, contrasts(c), behavioralData, initTrialConditions('movementTime',patience{p},'movementDir','cw','responseType','correct','highRewardSide','left'));
                    [~, hconds02] = selectCondition(expInfo, contrasts(c), behavioralData, initTrialConditions('movementTime',patience{p},'movementDir','ccw','responseType','correct','highRewardSide','right'));
                    highConds{c,1} = [hconds01 hconds02];

                    [~, lconds01] = selectCondition(expInfo, contrasts(c), behavioralData, initTrialConditions('movementTime',patience{p},'movementDir','cw','responseType','correct','highRewardSide','right'));
                    [~, lconds02] = selectCondition(expInfo, contrasts(c), behavioralData, initTrialConditions('movementTime',patience{p},'movementDir','ccw','responseType','correct','highRewardSide','left'));
                    lowConds{c,1} = [lconds01 lconds02];

                    [~, iconds1] = selectCondition(expInfo, contrasts(c), behavioralData, initTrialConditions('movementTime',patience{p},'responseType','incorrect','highRewardSide','left'));
                    [~, iconds2] = selectCondition(expInfo, -contrasts(c), behavioralData, initTrialConditions('movementTime',patience{p},'responseType','incorrect','highRewardSide','right'));
                    incConds{c,1} = [iconds1 iconds2];
                else
                    [~, hconds1] = selectCondition(expInfo, contrasts(c), behavioralData, initTrialConditions('movementTime',patience{p},'responseType','correct','highRewardSide','left'));
                    [~, hconds2] = selectCondition(expInfo, -contrasts(c), behavioralData, initTrialConditions('movementTime',patience{p},'responseType','correct','highRewardSide','right'));
                    highConds{c,1} = [hconds1 hconds2];

                    [~, lconds1] = selectCondition(expInfo, contrasts(c), behavioralData, initTrialConditions('movementTime',patience{p},'responseType','correct','highRewardSide','right'));
                    [~, lconds2] = selectCondition(expInfo, -contrasts(c), behavioralData, initTrialConditions('movementTime',patience{p},'responseType','correct','highRewardSide','left'));
                    lowConds{c,1} = [lconds1 lconds2];

                    [~, iconds1] = selectCondition(expInfo, contrasts(c), behavioralData, initTrialConditions('movementTime',patience{p},'responseType','incorrect','highRewardSide','left'));
                    [~, iconds2] = selectCondition(expInfo, -contrasts(c), behavioralData, initTrialConditions('movementTime',patience{p},'responseType','incorrect','highRewardSide','right'));
                    incConds{c,1} = [iconds1 iconds2];
                end
            end

            for c = 1:5
                highRew(c,:,1) = nanmean(nanmean(neuralData.eta.alignedResps{e}(highConds{c},:,:),3),1);
                lowRew(c,:,1) = nanmean(nanmean(neuralData.eta.alignedResps{e}(lowConds{c},:,:),3),1);
                noRew(c,:,1) = nanmean(nanmean(neuralData.eta.alignedResps{e}(incConds{c},:,:),3),1);
            end
            %
            eventWindow = neuralData(1).eta.eventWindow;

            mean_highRew = flipud(nanmean(highRew,3));
            mean_lowRew = flipud(nanmean(lowRew,3));
            mean_noRew = flipud(nanmean(noRew,3));

            % plotCells = getWhichCells('all',neuralData);
            ymax = .65;
            ymin = -.15;
            figure;
            set(gcf,'position',[1000 1400 800 230])
            grays = [.8 .6 .4 .2 .0];
            for iC = 1:size(highRew,1)
                subplot(1,3,1)
                hold on
                line([0 0],[ymin ymax],'LineStyle','--','LineWidth',1,'Color',[.5 .5 .5])
                ph = plot(eventWindow,(mean_highRew(iC,:)));
                set(ph,'LineWidth',2,'Color', grays(iC)*ones(1,3));
                xlim([-1 1])
                ylim([ymin ymax])
                title('high')
                if p == 1
                    ylabel('IMPATIENT')
                else
                    ylabel('PATIENT')
                end
                if e == 1
                    xlabel('time from stim')
                else
                    xlabel('time from move')
                end
            %     axis off

                subplot(1,3,2)
                hold on
                line([0 0],[ymin ymax],'LineStyle','--','LineWidth',1,'Color',[.5 .5 .5])
                ph = plot(eventWindow,(mean_lowRew(iC,:)));
                set(ph,'LineWidth',2,'Color', grays(iC)*ones(1,3));
                xlim([-1 1])
                ylim([ymin ymax])
                title('low')
                if e == 1
                    xlabel('time from stim')
                else
                    xlabel('time from move')
                end
            %     axis off

                subplot(1,3,3)
                hold on
                line([0 0],[ymin ymax],'LineStyle','--','LineWidth',1,'Color',[.5 .5 .5])
                pn = plot(eventWindow,(mean_noRew(iC,:)));
                set(pn,'LineWidth',2,'Color', grays(iC)*ones(1,3));
                xlim([-1 1])
                ylim([ymin ymax])
                title('none')
                if e == 1
                    xlabel('time from stim')
                else
                    xlabel('time from move')
                end
            %     axis off
            end
            if e == 1
                if p == 1
                    printfig(gcf,strcat(mouseName,'_',expDate,'_','stimOn_impatient'));
                else
                    printfig(gcf,strcat(mouseName,'_',expDate,'_','stimOn_patient'));
                end
            else
                if p == 1
                    printfig(gcf,strcat(mouseName,'_',expDate,'_','moveOn_impatient'));
                else
                    printfig(gcf,strcat(mouseName,'_',expDate,'_','moveOn_patient'));
                end
            end
        end
    end
    close all
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

