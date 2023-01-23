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
    
    [~, prevLeft] = selectCondition(expInfo,getUniqueContrasts(expInfo),behavioralData,initTrialConditions('pastMovementDir','cw'));
    [~, prevRight] = selectCondition(expInfo,getUniqueContrasts(expInfo),behavioralData,initTrialConditions('pastMovementDir','ccw'));
    
    
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
            trialList{1} = prevLeft;
            trialList{2} = prevRight;
            lg = {'Prev. left' 'Prev. right'};
            ttl = 'Choice history';
        case 'diffHistory'
            trialList{1} = prevEasy;
            trialList{2} = prevDifficult;
            lg = {'Prev. easy' 'Prev. hard'};
            ttl = 'Difficulty history';
        case 'all'
            trialList{1} = trials;
            lg = {'All trials'};
            ttl = 'All trials';
    end     

    for t = 1:length(trialList)
        [cc, ppl, ppr, ~, ~] = getPsychometric(expInfo, behavioralData, trialList{t},contrasts);           
        psychos(m).cc = cc;
        psychos(m).ppL(t,:) = ppl;
        psychos(m).ppR(t,:) = ppr;
    end

    clearvars -except mouseList expList hemList m psychos cc trialList trialTypes contrasts ttl lg

end

%%
lapseIdx = [1 2 8 9];
biasIdx = 4:6;

for m = 1:length(mouseList)
    ppl(:,:,m) = psychos(m).ppL;
    ppr(:,:,m) = psychos(m).ppR;
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

%%
mn = 'all';

range = find(strcmp(mouseList,mn));
if isempty(range)
    range = 1:length(mouseList);
end
meanPPL = nanmean(ppl(:,:,range),3);
meanPPR = nanmean(ppr(:,:,range),3);

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
    meanPPL = nanmean(ppl(:,:,range),3);
    meanPPR = nanmean(ppr(:,:,range),3);

    colors = [0.1 0.7 0.1; 1 .6 0];
    pc = cc*100;

    %
    set(gcf,'position',[1145 1275 1670 350])
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
        xticks([-100 -50 0 50 100])
        yticks([0 .5 1])
        xlabel('Contrast (%)')
        ylabel('Proportion of right choices')
        title(strcat(mn,{' '},'mean performance'))
        prettyPlot(gca)

    end
end
%%
figure;
for tp = 1:size(psychos(m).ppL,1)
    subplot(1,size(psychos(m).ppL,1),tp)
    lim1 = -.4;
    lim2 = .7;
    hold on;
    line([0 0],[-1 1],'Color',[.25 .25 .25],'LineStyle','--')
    line([-1 1],[0 0],'Color',[.25 .25 .25],'LineStyle','--')

    scatter(deltaBias(:,tp),deltaLapse(:,tp),200,[.75 .75 .75],'.')
    prettyPlot(gca)
    for u = 1:length(umn)
        errorbar(meanDB(u,tp),meanDL(u,tp),semDL(u,tp),semDL(u,tp),semDB(u,tp),semDB(u,tp),'Marker','.','MarkerSize',20,'capsize',0)
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
set(gcf,'position',[1256 1298 425 327]);
for t = 1:length(trialList)
    errorbar(cc,nanmean(meanDeltaPsycho(t,:,:),3),nanstd(meanDeltaPsycho(t,:,:),[],3)/sqrt(length(umn)),...
        'capsize',0,'Marker','o','MarkerSize',10,'MarkerFaceColor',colors(t,:),'MarkerEdgeColor','w','Color',colors(t,:),'LineWidth',1,'LineStyle',':')
    hold on
end
xlim([-1.05 1.05])
prettyPlot(gca)
xticks([-1 -.5 0 .5 1])
set(gca, 'XTickLabels', {'-100' '-50' '0' '50' '100'})
xlabel('Contrast (%)')
ylabel('\Delta Chose right')
title(ttl)
legend(lg,'location','nw')
legend boxoff

