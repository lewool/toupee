clearvars -except mouseList expList hemList
trialTypes = 'all';
contrasts = [-1 -.5 -.12 -.05 0 .05 .12 .5 1];

for m = 1:length(mouseList)
   
    mouseName = char(mouseList{m});
    expDate = char(expList{m}{1});
    expNum = expList{m}{2};
    [expRef, ~] = data.constructExpRef(mouseName,expDate,expNum);
    fprintf('Loading %s...\n',expRef)
    [behavioralData, expInfo] = data.loadBehavioralDataset(mouseName, expDate, expNum);
    expInfo.hemisphere = hemList(m);
    
    nt = numel(expInfo.block.events.endTrialTimes);
    stimOnTimes = behavioralData.eventTimes(1).daqTime;
    moveOnTimes = behavioralData.wheelMoves.epochs(5).onsetTimes;    
    RTs = moveOnTimes - stimOnTimes;
    trials = intersect(...
            find((RTs < prctile(RTs,97.5)) & (RTs > prctile(RTs,2.5))),...
            find(~isnan(RTs)));    
    RTs = RTs(trials);
    vels = abs(behavioralData.wheelMoves.epochs(5).peakVel(trials));
    
    rewardSizeValues = expInfo.block.events.rewardSizeValues;
    prevHigh = find(rewardSizeValues == 2.8) + 1;
    prevHigh(prevHigh>nt) = [];
    prevLow = find(rewardSizeValues == 1.4) + 1;
    prevLow(prevLow>nt) = [];
    
    [~, lateTrials] = selectCondition(expInfo,getUniqueContrasts(expInfo),behavioralData,initTrialConditions('movementTime','late'));
    [~, earlyTrials] = selectCondition(expInfo,getUniqueContrasts(expInfo),behavioralData,initTrialConditions('movementTime','early'));
    
    [~, prevCorrect] = selectCondition(expInfo,getUniqueContrasts(expInfo),behavioralData,initTrialConditions('pastResponseType','correct'));
    [~, prevIncorrect] = selectCondition(expInfo,getUniqueContrasts(expInfo),behavioralData,initTrialConditions('pastResponseType','incorrect'));
    
    [~, prevEasy] = selectCondition(expInfo,[-.1 -.5 .5 1],behavioralData,initTrialConditions('responseType','correct'));
    prevEasy = prevEasy+1;
    prevEasy(prevEasy>nt) = [];
    [~, prevDifficult] = selectCondition(expInfo,[-.05 0 .05],behavioralData,initTrialConditions('responseType','correct'));
    prevDifficult = prevDifficult+1;
    prevDifficult(prevDifficult>nt) = [];
    
    [~, prevLeftChoice] = selectCondition(expInfo,getUniqueContrasts(expInfo),behavioralData,initTrialConditions('pastMovementDir','cw'));
    [~, prevRightChoice] = selectCondition(expInfo,getUniqueContrasts(expInfo),behavioralData,initTrialConditions('pastMovementDir','ccw'));
    
    [~, prevLeftStim] = selectCondition(expInfo,contrasts(contrasts<0),behavioralData,initTrialConditions('pastMovementDir','all'));
    [~, prevRightStim] = selectCondition(expInfo,contrasts(contrasts>0),behavioralData,initTrialConditions('pastMovementDir','all'));
    prevLeftStim = prevLeftStim+1;
    prevLeftStim(prevLeftStim>nt) = [];
    prevRightStim = prevRightStim+1;
    prevRightStim(prevRightStim>nt) = [];
    
    whichTrials = intersect(lateTrials,trials);
    allnt(m) = nt;


    %%
    
    switch trialTypes
        case 'RT'
            trialList{1} = find(RTs < prctile(RTs,25));
%             trialList{2} = find((RTs >= prctile(RTs,33.3)) & (RTs < prctile(RTs,66.6)));
            trialList{2} = find((RTs >= prctile(RTs,75)));
%             trialList{1} = find(RTs < prctile(RTs,50));
%             trialList{2} = find(RTs >= prctile(RTs,50));
            lg = {'Fastest RTs' 'Slowest RTs'};
            ttl = 'Reaction times';
%             colors = [0 .25 0; 0 .5 0; 0 .75 0];
        case 'sessionTime'
            ss = floor(length(trials)/2);
            trialList{1} = 1:ss;
            trialList{2} = ss+1:2*ss;
            lg = {'First half' 'Second half'};
            ttl = 'Time in session';
        case 'impulsive'
            ss = floor(length(trials)/2);
            trialList{1} = earlyTrials;
            trialList{2} = lateTrials;
            lg = {'Impulsive' 'Patient'};
            ttl = 'Response relative to cue';
        case 'velocity'
            trialList{1} = find(vels < prctile(vels,25));
%             trialList{2} = find((vels >= prctile(vels,33.3)) & (vels < prctile(vels,66.6)));
            trialList{2} = find((vels >= prctile(vels,75)));
%             trialList{1} = find(vels < prctile(vels,50));
%             trialList{2} = find(vels >= prctile(vels,50));
            lg = {'Slowest wheel' 'Fastest wheel'};
            ttl = 'Wheel speed';
%             colors = [.25 0 .8; .25 .4 .8; .25 .8 .8];
        case 'rewardHistory'
            trialList{1} = prevCorrect;
            trialList{2} = prevIncorrect;
            lg = {'Prev. correct' 'Prev. incorrect'};
            ttl = 'Outcome history';
%             colors = [0 .5 0; .7 .1 0];
        case 'valueHistory'
            trialList{1} = prevHigh;
            trialList{2} = prevLow;
            lg = {'Prev. high' 'Prev. low'};
            ttl = 'Value history';
        case 'choiceHistory'
            trialList{1} = prevLeftChoice;
            trialList{2} = prevRightChoice;
            lg = {'Prev. left' 'Prev. right'};
            ttl = 'Choice history';
        case 'stimHistory'
            trialList{1} = prevLeftStim;
            trialList{2} = prevRightStim;
            lg = {'Prev. left' 'Prev. right'};
            ttl = 'Stim side history';
        case 'diffHistory'
            trialList{1} = prevEasy;
            trialList{2} = prevDifficult;
            lg = {'Prev. easy' 'Prev. hard'};
            ttl = 'Difficulty history';
        case 'all'
            trialList{1} = whichTrials;
            lg = {'All trials'};
            ttl = 'All trials';
    end     

    for t = 1:length(trialList)
        [cc, ppl, ppr, ppt, ~, ~, ~] = getPsychometric(expInfo, behavioralData, trialList{t},contrasts);           
        psychos(m).cc = cc;
        psychos(m).ppL(t,:) = ppl;
        psychos(m).ppR(t,:) = ppr;
        psychos(m).ppT(t,:) = ppt;
    end

    clearvars -except mouseList expList hemList m psychos cc trialList trialTypes contrasts ttl lg allnt

end

%%
lapseIdx = [1 2 8 9];
biasIdx = 4:6;

for m = 1:length(mouseList)
    ppl(:,:,m) = psychos(m).ppL;
    ppr(:,:,m) = psychos(m).ppR;
    ppt(:,:,m) = psychos(m).ppT;
end

for m = 1:length(mouseList)
    deltaBias(m,:) = nanmean(ppr(:,biasIdx,m) - ppl(:,biasIdx,m),2);
    deltaLapse(m,:) = nanmean(ppr(:,lapseIdx,m) - ppl(:,lapseIdx,m),2);
end

umn = unique(mouseList);
for u = 1:length(umn)   
    mn = umn{u};
    range = find(strcmp(mouseList,mn));
    meanDB(u,:) = nanmean(deltaBias(range,:),1);
    meanDL(u,:) = nanmean(deltaLapse(range,:),1);
    semDB(u,:) = nanstd(deltaBias(range,:))/sqrt(length(range));
    semDL(u,:) = nanstd(deltaLapse(range,:))/sqrt(length(range));
end

%% compare, split over blocks
mn = 'all';

range = find(strcmp(mouseList,mn));
if isempty(range)
    range = 1:length(mouseList);
end

colors = [0.1 0.7 0.1; 1 .6 0];
pc = cc*100;

%
% figure;
set(gcf,'position',[1145 1275 1670 350])
for tp = 1:size(ppl,1)
    subplot(1,size(ppl,1),tp)
    hold on;
    line([-105 105],[.5 .5],'Color',[.5 .5 .5],'LineStyle',':');
    line([0 0],[-.05 1.05],'Color',[.5 .5 .5],'LineStyle',':')
    errorbar(pc,nanmean(ppr(tp,:,range),3),nanstd(ppr(tp,:,range),[],3)/sqrt(length(range)),...
        'capsize',0,'Marker','o','MarkerSize',10,'MarkerFaceColor',[1 .6 0],'MarkerEdgeColor','w','Color',[1 .6 0],'LineWidth',2)
    hold on
    % plot(contrastSet{1},mean(pplmat(1:14,:),1),'Color',[0.1 0.7 0.1],'LineWidth',1.5)
    errorbar(pc,nanmean(ppl(tp,:,range),3),nanstd(ppl(tp,:,range),[],3)/sqrt(length(range)),...
        'capsize',0,'Marker','o','MarkerSize',10,'MarkerFaceColor',[0.1 0.7 0.1],'MarkerEdgeColor','w','Color',[0.1 0.7 0.1],'LineWidth',2)

    prettyPlot(gca)
    xlim([-110 110])
    ylim([-.05 1.05])
%     xticks([-100 -50 0 50 100])
    yticks([0 .5 1])
    xlabel('Contrast (%)')
    ylabel('Proportion of right choices')
    title(strcat(mn,{' '},'mean performance'))
    prettyPlot(gca)

end

%% compare (no block split)
mn = 'all';

range = find(strcmp(mouseList,mn));
if isempty(range)
    range = 1:length(mouseList);
end

colors = [0 .4 1; 1 0 0];
pc = cc*100;

figure;
% set(gcf,'position',[1145 1275 670 350])
for tp = 1:size(ppt,1)
    hold on;
    line([-105 105],[.5 .5],'Color',[.5 .5 .5],'LineStyle',':');
    line([0 0],[-.05 1.05],'Color',[.5 .5 .5],'LineStyle',':')
    errorbar(pc,nanmean(ppt(tp,:,range),3),nanstd(ppt(tp,:,range),[],3)/sqrt(length(range)),...
        'capsize',0,'Marker','o','MarkerSize',10,'MarkerFaceColor',colors(tp,:),'MarkerEdgeColor','w','Color',colors(tp,:),'LineWidth',2)
    hold on
%     % plot(contrastSet{1},mean(pplmat(1:14,:),1),'Color',[0.1 0.7 0.1],'LineWidth',1.5)
%     errorbar(pc,nanmean(ppt(tp,:,range),3),nanstd(ppt(tp,:,range),[],3)/sqrt(length(range)),...
%         'capsize',0,'Marker','o','MarkerSize',10,'MarkerFaceColor',colors(tp,:),'MarkerEdgeColor','w','Color',colors(tp,:),'LineWidth',2)

    prettyPlot(gca)
    xlim([-110 110])
    ylim([-.05 1.05])
    xticks([-100 -50 0 50 100])
    yticks([0 .5 1])
    xlabel('Contrast (%)')
    ylabel('Proportion of right choices')
    title(strcat(mn,{' '},'mean performance'))
    prettyPlot(gca)

end

%%
figure;
alp = 0.5;
for u = 1:length(umn)   
    mn = umn{u};
    
    range = find(strcmp(mouseList,mn));
    if isempty(range)
        range = 1:length(mouseList);
    end
    meanPPL(u,:) = nanmean(ppl(:,:,range),3);
    meanPPR(u,:) = nanmean(ppr(:,:,range),3);

    colors = [0.1 0.7 0.1; 1 .6 0];
    pc = cc*100;
    %
    for tp = 1:size(ppl,1)
        subplot(1,size(ppl,1),tp)
        hold on;
        line([-105 105],[.5 .5],'Color',[.5 .5 .5],'LineStyle',':');
        line([0 0],[-.05 1.05],'Color',[.5 .5 .5],'LineStyle',':')
        ln1 = plot(pc,nanmean(ppr(tp,:,range),3));
        ln1.Color = [1,.6,0,alp];
        hold on
        ln2 = plot(pc,nanmean(ppl(tp,:,range),3));
        ln2.Color = [0.1,.7,0.1,alp];
        prettyPlot(gca)
        xlim([-110 110])
        ylim([-.05 1.05])
        xticks([-100 -50 -12 0 12 50 100])
        set(gca, 'XTickLabels', {'-100' '-50' '-12' '0' '12' '50' '100'})
        xlabel('Left contrast (%)                          Right contrast (%)')

        yticks([0 .5 1])
        
        ylabel('Proportion of right choices')
        prettyPlot(gca)

    end
end

errorbar(pc,nanmean(meanPPR,1),nanstd(meanPPR,[],1)/sqrt(size(meanPPR,1)),...
        'capsize',0,'Marker','o','MarkerSize',10,'MarkerFaceColor',[1 .6 0],'MarkerEdgeColor','w','Color',[1 .6 0],'LineWidth',2)
errorbar(pc,nanmean(meanPPL,1),nanstd(meanPPL,[],1)/sqrt(size(meanPPL,1)),...
        'capsize',0,'Marker','o','MarkerSize',10,'MarkerFaceColor',[0.1 0.7 0.1],'MarkerEdgeColor','w','Color',[0.1 0.7 0.1],'LineWidth',2)

%%

fl= polyfit(deltaBias,deltaLapse,1);
xfit = linspace(-.3,.5,2);
yfit = fl(1)*xfit + fl(2);
colors = [0 .447 .741; ...
          .85 .325 .098; ...
          .929 .694 .125; ...
          .494 .184 .556; ...
          .466 .674 .188; ...
          .301 .745 .933];
figure;
for tp = 1:size(psychos(m).ppL,1)
    subplot(1,size(psychos(m).ppL,1),tp)
    lim1 = -.3;
    lim2 = .5;
    hold on;
    line([0 0],[-1 1],'Color',[.25 .25 .25],'LineStyle',':')
    line([-1 1],[0 0],'Color',[.25 .25 .25],'LineStyle',':')
    line([xfit(1) xfit(end)],[yfit(1) yfit(end)],'Color',[.5 .5 .5],'LineStyle','--')
%     scatter(deltaBias(:,tp),deltaLapse(:,tp),200,[.75 .75 .75],'.')
    prettyPlot(gca)
    for u = 1:length(umn)
        mn = umn{u};
        range = find(strcmp(mouseList,mn));
        scatter(deltaBias(range,tp),deltaLapse(range,tp),20,[.75 .75 .75],'o','MarkerFaceColor',colors(u,:),'MarkerEdgeColor','none','MarkerFaceAlpha',.3)
        errorbar(meanDB(u,tp),meanDL(u,tp),semDL(u,tp),semDL(u,tp),semDB(u,tp),semDB(u,tp),'Marker','.','MarkerEdgeColor',colors(u,:),'MarkerSize',20,'capsize',0,'color',colors(u,:))
        hold on
    end
    xlim([lim1 lim2])
    ylim([lim1 lim2])
    xticks([lim1:.2:lim2])
    yticks([lim1:.2:lim2])
    axis square
    xlabel('\Delta Bias')
    ylabel('\Delta Lapse')
    
end

%%
colors = [0 0 0; .5 .5 .5; .8 .8 .8];
% pc = sqrt(abs(cc)).*sign(cc);
pc = cc*100;

for m = 1:length(mouseList)
    ppl(:,:,m) = psychos(m).ppL;
    ppr(:,:,m) = psychos(m).ppR;
end

clear meanDeltaPsycho
umn = unique(mouseList);
for u = 1:length(umn)   
    for t = 1:length(trialList)
        mn = umn{u};
        range = find(strcmp(mouseList,mn));
        meanDeltaPsycho(t,:,u) = nanmean(ppr(t,:,range),3) - nanmean(ppl(t,:,range),3);
    end
end

% for m = 1:length(mouseList)
%     for t = 1:length(trialList)
%         meanDeltaPsycho(t,:,m) = nanmean(ppr(t,:,m),3) - nanmean(ppl(t,:,m),3);
%     end
% end

figure;
% set(gcf,'position',[1256 1298 425 327]);
hold on;
for u = 1:length(umn)
    ln = plot(pc,meanDeltaPsycho(:,:,u),'k','LineWidth',.5);
    ln.Color(4) = .2;
end

for t = 1:length(trialList)
    errorbar(pc,nanmean(meanDeltaPsycho(t,:,:),3),nanstd(meanDeltaPsycho(t,:,:),[],3)/sqrt(length(umn)),...
        'capsize',0,'Marker','o','MarkerSize',10,'MarkerFaceColor',colors(t,:),'MarkerEdgeColor','w','Color',colors(t,:),'LineWidth',2,'LineStyle','-')
    hold on
end
line([-100 100],[0 0],'Color',[.25 .25 .25],'LineStyle',':')
xlim([-107 107])
prettyPlot(gca)
xticks([-100 -50 -12 0 12 50 100])
set(gca, 'XTickLabels', {'-100' '-50' '-12' '0' '12' '50' '100'})
xlabel('Left contrast (%)                          Right contrast (%)')
ylabel('\Delta Proportion right choices')
% title(ttl)
% legend(lg,'location','nw')
% legend boxoff

%% delta ANOVA

dp = squeeze(ppr - ppl)';
abs_dp = [(dp(:,1:4) + fliplr(dp(:,6:9))/2) dp(:,5)];
anmat = reshape(abs_dp,numel(abs_dp),1);
anmat(:,2) = [ones(size(abs_dp,1),1);...
    .5*ones(size(abs_dp,1),1);...
    .12*ones(size(abs_dp,1),1);...
    .05*ones(size(abs_dp,1),1);...
    0*ones(size(abs_dp,1),1)];

anovan(anmat(:,1),anmat(:,2))
%%

mn = 'all';

range = find(strcmp(mouseList,mn));
if isempty(range)
    range = 1:length(mouseList);
end
ppA = (squeeze(ppl(1,:,:)) + squeeze(ppr(1,:,:)))/2;
ppB = (squeeze(ppl(2,:,:)) + squeeze(ppr(2,:,:)))/2;

colors = [0.1 0.7 0.1; 1 .6 0];
pc = cc;

% figure;
set(gcf,'position',[1145 1275 1670 350])
hold on;
line([-105 105],[.5 .5],'Color',[.5 .5 .5],'LineStyle',':');
line([0 0],[-.05 1.05],'Color',[.5 .5 .5],'LineStyle',':')
errorbar(pc,nanmean(ppA(:,range),2),nanstd(ppA(:,range),[],2)/sqrt(length(range)),...
    'capsize',0,'Marker','o','MarkerSize',10,'MarkerFaceColor',[1 .6 0],'MarkerEdgeColor','w','Color',[1 .6 0],'LineWidth',2)
hold on
% plot(contrastSet{1},mean(pplmat(1:14,:),1),'Color',[0.1 0.7 0.1],'LineWidth',1.5)
errorbar(pc,nanmean(ppB(:,range),2),nanstd(ppB(:,range),[],2)/sqrt(length(range)),...
    'capsize',0,'Marker','o','MarkerSize',10,'MarkerFaceColor',[0.1 0.7 0.1],'MarkerEdgeColor','w','Color',[0.1 0.7 0.1],'LineWidth',2)

xlim([-1.05 1.05])
prettyPlot(gca)
xticks([-1 -.5 0 .5 1])
set(gca, 'XTickLabels', {'-100' '-50' '0' '50' '100'})
xlabel('Contrast (%)')
ylabel('Chose right')
title(ttl)
legend(lg,'location','nw')
legend boxoff
