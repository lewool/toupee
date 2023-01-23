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
    etas = [2];
       whichCells = 'all'; %choose from 'pLabels' array
    if strcmp(whichCells, 'all')
        plotCells = (1:size(neuralData.eta.alignedResps{1},3))';
    elseif strcmp(whichCells,'contraStim')
        if hemisphere > 0
            plotCells = find(neuralData.stats.pValues(:,strcmp(neuralData.stats.labels,'leftStim')) < 0.05);
        else
            plotCells = find(neuralData.stats.pValues(:,strcmp(neuralData.stats.labels,'rightStim')) < 0.05);
        end
    else
        plotCells = find(neuralData.stats.pValues(:,strcmp(neuralData.stats.labels,whichCells)) < 0.05);
    end
    nc = length(plotCells);
    [alignedWheel, wheelWindow] = getAlignedWheelTraces(expInfo, behavioralData, 2);
    
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
                    [~, iconds2] = selectCondition(expInfo, contrasts(c), behavioralData, initTrialConditions('movementTime',patience{p},'responseType','incorrect','highRewardSide','right'));
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
                highRew(c,:,1) = nanmean(nanmean(neuralData.eta.alignedResps{etas(e)}(highConds{c},:,plotCells),3),1);
                lowRew(c,:,1) = nanmean(nanmean(neuralData.eta.alignedResps{etas(e)}(lowConds{c},:,plotCells),3),1);
                noRew(c,:,1) = nanmean(nanmean(neuralData.eta.alignedResps{etas(e)}(incConds{c},:,plotCells),3),1);
                
                highRewWheel(c,:,1) = nanmean(abs(alignedWheel(highConds{c},:)),1);
                lowRewWheel(c,:,1) = nanmean(abs(alignedWheel(lowConds{c},:)),1);
                noRewWheel(c,:,1) = nanmean(abs(alignedWheel(incConds{c},:)),1);
            end
            %
            eventWindow = neuralData(1).eta.eventWindow;

            mean_highRew = flipud(nanmean(highRew,3));
            mean_lowRew = flipud(nanmean(lowRew,3));
            mean_noRew = flipud(nanmean(noRew,3));
            
            mean_highRewWheel = flipud(nanmean(highRewWheel,3));
            mean_lowRewWheel = flipud(nanmean(lowRewWheel,3));
            mean_noRewWheel = flipud(nanmean(noRewWheel,3));
            
            contResps(m).(matlab.lang.makeValidName(patience{p})).highRew = mean_highRew;
            contResps(m).(matlab.lang.makeValidName(patience{p})).lowRew = mean_lowRew;
            contResps(m).(matlab.lang.makeValidName(patience{p})).noRew = mean_noRew;
            
            wheelPos(m).(matlab.lang.makeValidName(patience{p})).highRew = mean_highRewWheel;
            wheelPos(m).(matlab.lang.makeValidName(patience{p})).lowRew = mean_lowRewWheel;
            wheelPos(m).(matlab.lang.makeValidName(patience{p})).noRew = mean_noRewWheel;

        end
    end
    clearvars -except mouseList expList hemList fovList contResps wheelPos eventWindow wheelWindow
end
%%
whichStruct = 'wheel';
if strcmp(whichStruct,'wheel')
    whichData = wheelPos;
    ew = wheelWindow(2:end);
else
    whichData = contResps;
    ew = eventWindow;
end

for m = 1:length(mouseList)
    siteList{m} = strcat(mouseList{m},'_',num2str(fovList(m)));
end
try
    
    ml = siteList;
    mlu = unique(ml);
catch
    ml = cellfun(@char,siteList,'UniformOutput',false);
    mlu = unique(ml);
end

clear meh mel men mlh mll mln allEarlyHigh allEarlyLow allEarlyNone allLateHigh allLateLow allLateNone
for m = 1:length(mlu)
    mn = mlu{m};

    range = find(strcmp(ml,mn));
    if isempty(range)
        range = 1:length(ml);
    end
    
    for r = 1:length(range)
        meh(:,:,r) = whichData(range(r)).early.highRew;
        mel(:,:,r) = whichData(range(r)).early.lowRew;
        men(:,:,r) = whichData(range(r)).early.noRew;
        mlh(:,:,r) = whichData(range(r)).late.highRew;
        mll(:,:,r) = whichData(range(r)).late.lowRew;
        mln(:,:,r) = whichData(range(r)).late.noRew;
    end
    
    if strcmp(whichStruct,'wheel')
        allEarlyHigh(:,:,m) = (diff(nanmean(meh,3),[],2));
        allEarlyLow(:,:,m) = (diff(nanmean(mel,3),[],2));
        allEarlyNone(:,:,m) = (diff(nanmean(men,3),[],2));
        allLateHigh(:,:,m) = (diff(nanmean(mlh,3),[],2));
        allLateLow(:,:,m) = (diff(nanmean(mll,3),[],2));
        allLateNone(:,:,m) = (diff(nanmean(mln,3),[],2));
    else
        allEarlyHigh(:,:,m) = nanmean(meh,3);
        allEarlyLow(:,:,m) = nanmean(mel,3);
        allEarlyNone(:,:,m) = nanmean(men,3);
        allLateHigh(:,:,m) = nanmean(mlh,3);
        allLateLow(:,:,m) = nanmean(mll,3);
        allLateNone(:,:,m) = nanmean(mln,3);
    end

end

%%

figure;
set(gcf,'position',[1000 1400 800 230])
grays = [.8 .6 .4 .2 .0];
yl = [-0.01 .18];

subplot(2,3,1)
for c = 1:size(allEarlyHigh,1)
    if strcmp(whichStruct,'wheel')
        plot(ew,smooth(nanmean(allLateHigh(c,:,:),3),20),'Color',grays(c)*[1 1 1],'LineWidth',1.5)
    else
        plot(ew,nanmean(allLateHigh(c,:,:),3),'Color',grays(c)*[1 1 1],'LineWidth',1.5)
    end
    hold on
end
prettyPlot(gca);
xlim([-1 1])
ylim(yl)
line([0 0],yl,'Color',[.5 .5 .5],'LineStyle',':')

subplot(2,3,2)
for c = 1:size(allEarlyHigh,1)
    if strcmp(whichStruct,'wheel')
        plot(ew,smooth(nanmean(allLateLow(c,:,:),3),20),'Color',grays(c)*[1 1 1],'LineWidth',1.5)
    else
        plot(ew,nanmean(allLateLow(c,:,:),3),'Color',grays(c)*[1 1 1],'LineWidth',1.5)
    end
    hold on
end
prettyPlot(gca);
xlim([-1 1])
ylim(yl)
line([0 0],yl,'Color',[.5 .5 .5],'LineStyle',':')

subplot(2,3,3)
for c = 1:size(allEarlyHigh,1)
    if strcmp(whichStruct,'wheel')
        plot(ew,smooth(nanmean(allLateNone(c,:,:),3),20),'Color',grays(c)*[1 1 1],'LineWidth',1.5)
    else
        plot(ew,nanmean(allLateNone(c,:,:),3),'Color',grays(c)*[1 1 1],'LineWidth',1.5)
    end
    hold on
end
prettyPlot(gca);
xlim([-1 1])
ylim(yl)
line([0 0],yl,'Color',[.5 .5 .5],'LineStyle',':')

yl = [-0.01 .15];

subplot(2,3,4)
for c = 1:size(allEarlyHigh,1)
    if strcmp(whichStruct,'wheel')
        plot(ew,smooth(nanmean(allEarlyHigh(c,:,:),3),20),'Color',grays(c)*[1 1 1],'LineWidth',1.5)
    else
        plot(ew,nanmean(allEarlyHigh(c,:,:),3),'Color',grays(c)*[1 1 1],'LineWidth',1.5)
    end
    hold on
end
prettyPlot(gca);
xlim([-1 1])
ylim(yl)
line([0 0],yl,'Color',[.5 .5 .5],'LineStyle',':')

subplot(2,3,5)
for c = 1:size(allEarlyHigh,1)
     if strcmp(whichStruct,'wheel')
        plot(ew,smooth(nanmean(allEarlyLow(c,:,:),3),20),'Color',grays(c)*[1 1 1],'LineWidth',1.5)
    else
        plot(ew,nanmean(allEarlyLow(c,:,:),3),'Color',grays(c)*[1 1 1],'LineWidth',1.5)
    end
    hold on
end
prettyPlot(gca);
xlim([-1 1])
ylim(yl)
line([0 0],yl,'Color',[.5 .5 .5],'LineStyle',':')

subplot(2,3,6)
for c = 1:size(allEarlyHigh,1)
     if strcmp(whichStruct,'wheel')
        plot(ew,smooth(nanmean(allEarlyNone(c,:,:),3),20),'Color',grays(c)*[1 1 1],'LineWidth',1.5)
    else
        plot(ew,nanmean(allEarlyNone(c,:,:),3),'Color',grays(c)*[1 1 1],'LineWidth',1.5)
    end
    hold on
end
prettyPlot(gca);
xlim([-1 1])
ylim(yl)
line([0 0],yl,'Color',[.5 .5 .5],'LineStyle',':')

for s = 1:6
    subplot(2,3,s)
    xlabel('Time from movement (s)')
    if s == 1
        title('High')
        legend('0', '5', '12', '50', '100%','Location','nw')
        legend boxoff
    elseif s == 2
        title('Low')
    elseif s == 3
        title('None')
    end
end

%%

%
figure;
set(gcf,'position',[1145 1275 1670 350])
for tp = 1:size(ppl,1)
    subplot(1,size(ppl,1),tp)
    hold on;
    line([-105 105],[.5 .5],'Color',[.5 .5 .5],'LineStyle',':');
    line([0 0],[-.05 1.05],'Color',[.5 .5 .5],'LineStyle',':')
    errorbar(pc,nanmean(ppr(tp,:,range),3),nanstd(ppr(tp,:,range),[],3)/sqrt(length(range)),...
        'capsize',0,'Marker','o','MarkerSize',10,'MarkerFaceColor',[1 .6 0],'MarkerEdgeColor','w','Color',[1 .6 0],'LineWidth',1)
    hold on
    % plot(contrastSet{1},mean(pplmat(1:14,:),1),'Color',[0.1 0.7 0.1],'LineWidth',1.5)
    errorbar(pc,nanmean(ppl(tp,:,range),3),nanstd(ppl(tp,:,range),[],3)/sqrt(length(range)),...
        'capsize',0,'Marker','o','MarkerSize',10,'MarkerFaceColor',[0.1 0.7 0.1],'MarkerEdgeColor','w','Color',[0.1 0.7 0.1],'LineWidth',1)

    prettyPlot(gca)
    xlim([-105 105])
    ylim([-.05 1.05])
    xticks([-100 -50 0 50 100])
    yticks([0 .5 1])
    xlabel('Contrast (%)')
    ylabel('Proportion of right choices')
    title(strcat(mn,{' '},'mean performance'))
    prettyPlot(gca)

end

%             % plotCells = getWhichCells('all',neuralData);
%             ymax = .14;
%             ymin = .0;
%             figure;
%             set(gcf,'position',[1000 1400 800 230])
%             grays = [.8 .6 .4 .2 .0];
%             for iC = 1:size(highRew,1)
%                 subplot(1,3,1)
%                 hold on
%                 line([0 0],[ymin ymax],'LineStyle','--','LineWidth',1,'Color',[.5 .5 .5])
%                 ph = plot(eventWindow,(mean_highRew(iC,:)));
%                 set(ph,'LineWidth',2,'Color', grays(iC)*ones(1,3));
%                 xlim([-1 1])
%                 ylim([ymin ymax])
%                 title('high')
%                 if p == 1
%                     ylabel('IMPATIENT')
%                 else
%                     ylabel('PATIENT')
%                 end
%                 if etas(e) == 1
%                     xlabel('time from stim')
%                 elseif etas(e) == 2
%                     xlabel('time from move')
%                     elseif etas(e) == 3
%                     xlabel('time from feedback')
%                     elseif etas(e) == 4
%                     xlabel('time from cue')
%                 end
%             %     axis off
% 
%                 subplot(1,3,2)
%                 hold on
%                 line([0 0],[ymin ymax],'LineStyle','--','LineWidth',1,'Color',[.5 .5 .5])
%                 ph = plot(eventWindow,(mean_lowRew(iC,:)));
%                 set(ph,'LineWidth',2,'Color', grays(iC)*ones(1,3));
%                 xlim([-1 1])
%                 ylim([ymin ymax])
%                 title('low')
%                 if etas(e) == 1
%                     xlabel('time from stim')
%                 elseif etas(e) == 2
%                     xlabel('time from move')
%                     elseif etas(e) == 3
%                     xlabel('time from feedback')
%                     elseif etas(e) == 4
%                     xlabel('time from cue')
%                 end
%             %     axis off
% 
%                 subplot(1,3,3)
%                 hold on
%                 line([0 0],[ymin ymax],'LineStyle','--','LineWidth',1,'Color',[.5 .5 .5])
%                 pn = plot(eventWindow,(mean_noRew(iC,:)));
%                 set(pn,'LineWidth',2,'Color', grays(iC)*ones(1,3));
%                 xlim([-1 1])
%                 ylim([ymin ymax])
%                 title('none')
%                 if etas(e) == 1
%                     xlabel('time from stim')
%                 elseif etas(e) == 2
%                     xlabel('time from move')
%                     elseif etas(e) == 3
%                     xlabel('time from feedback')
%                     elseif etas(e) == 4
%                     xlabel('time from cue')
%                 end
%             %     axis off
%             end
% %             if etas(e) == 1
% %                 if p == 1
% %                     printfig(gcf,strcat(mouseName,'_',expDate,'_','stimOn_impatient'));
% %                 else
% %                     printfig(gcf,strcat(mouseName,'_',expDate,'_','stimOn_patient'));
% %                 end
% %             else
% %                 if p == 1
% %                     printfig(gcf,strcat(mouseName,'_',expDate,'_','moveOn_impatient'));
% %                 else
% %                     printfig(gcf,strcat(mouseName,'_',expDate,'_','moveOn_patient'));
% %                 end
% %             end

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

