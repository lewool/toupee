figure;
for m = 1:21
     %load data
    mouseName = char(mouseList{m});
    expDate = char(expList{m}{1});
    expNum = expList{m}{2};
    [expRef, ~] = data.constructExpRef(mouseName,expDate,expNum);
    sessDir = fullfile('G:\Workspaces',mouseName,expDate,num2str(expNum));
    cd(sessDir)
    load('KernelAnalysis.mat')
    subplot(3,7,m)
    bar(sum(kernelAnalysis.cellFeatureStrength(kernelAnalysis.maxEV > .01,:) > .01)./(sum(kernelAnalysis.maxEV > .01))*100)
    prettyPlot(gca)
    title(expRef)
    
end

%%

for m=1:21
    subplot(3,7,m);
    mouseName = char(mouseList{m});
    expDate = char(expList{m}{1});
    expNum = expList{m}{2};
    [expRef, ~] = data.constructExpRef(mouseName,expDate,expNum);
    title(expRef,'Interpreter','none')
end

%%

figure;
for tn = 1:96
    subplot(8,12,tn)
    yyaxis left
    plot(neuralData.eta.eventWindow, neuralData.eta.alignedResps{1}(earlyTrials(tn),:,3)','LineWidth',2)
    hold on;
    plot(neuralData.eta.eventWindow, kernelAnalysis.eta.alignedResps{1}(earlyTrials(tn),:,3)','LineStyle',':','LineWidth',2)
    line([0 0],[-.2 1.2],'Color',[.5 .5 .5],'LineStyle','--')
    line([0.8 0.8],[-.2 1.2],'Color',[.5 .5 .5],'LineStyle','--')
    ylim([-.2 1.2])
    yyaxis right
    [~,i] = max(abs(alignedWheel(tn,:)));
    if alignedWheel(tn,i)<0
        plot(plotWindow, -alignedWheel(tn,:))
    else
        plot(plotWindow, alignedWheel(tn,:))
    end
    xlim([-.5 2])
    ylim([-50 300])
    prettyPlot(gca)
end

%%
predictY = 'choiceIndStim';

ew = -2:.1:2;
timeRange = 6:2:36;
ew = ew(timeRange);
pdist = zeros(1,16);
figure;
set(gcf,'position',[135 41 2666 1608]);
for t=1:36
    subplot(6,6,t)
    mouseName = char(mouseList{t});
    expDate = char(expList{t}{1});
    expNum = expList{t}{2};
    [expRef, ~] = data.constructExpRef(mouseName,expDate,expNum);
    try
    ptest = predictions(t).(matlab.lang.makeValidName(predictY)).gocue.true > prctile(predictions(t).(matlab.lang.makeValidName(predictY)).gocue.pseudo,97.5) ...
        | predictions(t).(matlab.lang.makeValidName(predictY)).gocue.true < prctile(predictions(t).(matlab.lang.makeValidName(predictY)).gocue.pseudo,2.5);
    pdist = pdist + ptest;
    plot(ew, predictions(t).(matlab.lang.makeValidName(predictY)).gocue.pseudo,'Color',[.65 .65 .65]);
    hold on
    plot(ew, predictions(t).(matlab.lang.makeValidName(predictY)).gocue.true,'Color','k','LineWidth',3)
    line([0 0],[0 1],'Color',[.5 .5 .5],'LineStyle','--')
    line([-1 -1],[0 1],'Color',[.5 .5 .5],'LineStyle','--')
    line([-1.5 1.5],[0.5 0.5],'Color','k','LineStyle','-')
    plot(ew,ptest-.05,'Marker','*','LineStyle','none','Color','k')
    catch
    end
    ylim([0.2 1])
    xlim([-1.5 1.5])
    prettyPlot(gca)
    title(expRef,'Interpreter','none')
end

%
for b = 1:16
    binomtest(b) = myBinomTest(pdist(b),36,.05,'two');
end
allPseudo = [];
for t=1:36
    allTrue(t,:) = predictions(t).(matlab.lang.makeValidName(predictY)).gocue.true;
    for p = 1:length(predictions(t).(matlab.lang.makeValidName(predictY)).gocue.pseudo)
        allPseudo(t,p,:) = predictions(t).(matlab.lang.makeValidName(predictY)).gocue.pseudo(p,:);
    end
end

meanTrue = mean(allTrue,1);
meanPseudo = squeeze(mean(allPseudo,1));

figure;
set(gcf,'position',[996 920 639 312])
hold on;
plot(ew,meanPseudo,'Color',[.65 .65 .65]);
plot(ew,meanTrue,'Color','k','LineWidth',3);
line([0 0],[0 1],'Color',[.5 .5 .5],'LineStyle','--')
line([-1 -1],[0 1],'Color',[.5 .5 .5],'LineStyle','--')
plot(ew,(binomtest < .05)-.276,'Marker','*','LineStyle','none','Color','k')
ylim([0.475 .725])
xlim([-1.5 1.5])
prettyPlot(gca)
xlabel('Time from go cue (s)')
ylabel('Normalized log-likelihood')
% for t=1:16
%     if binomtest(t) <= 0.05 && binomtest(t) > 0.01
%         text(ew(t),0.724,'*','HorizontalAlignment','center')
%     elseif binomtest(t) <= 0.01 && binomtest(t) > 0.001
%         text(ew(t),0.724,'**','HorizontalAlignment','center')
%     elseif binomtest(t) <= 0.001
%         text(ew(t),0.724,'***','HorizontalAlignment','center')
%     end
% end
%
figure;
set(gcf,'position',[996 600 639 312])
bar(ew, pdist/36,'FaceColor','k')
prettyPlot(gca)
line([0 0],[0 1],'Color',[.5 .5 .5],'LineStyle','--')
line([-1 -1],[0 1],'Color',[.5 .5 .5],'LineStyle','--')
hold on
ylabel('Proportion of sessions')
xlim([-1.6 1.6])
ylim([0 1])
xlabel('Time from go cue (s)')
% plot(ew,(binomtest<.05)-.01,'Marker','*','LineStyle','none','Color','k')
line([-1.6 1.6],[0.1389 .1389],'Color','r','LineStyle','-')

%%
predictY = 'stimIndChoice';

ew = -2:.1:2;
timeRange = 6:2:36;
ew = ew(timeRange);
pdist = zeros(1,16);
figure;
set(gcf,'position',[135 41 2666 1608]);
for t=1:36
    subplot(6,6,t)
    mouseName = char(mouseList{t});
    expDate = char(expList{t}{1});
    expNum = expList{t}{2};
    [expRef, ~] = data.constructExpRef(mouseName,expDate,expNum);
    try
    ptest = mutualInformation(t).(matlab.lang.makeValidName(predictY)).gocue.true > prctile(mutualInformation(t).(matlab.lang.makeValidName(predictY)).gocue.pseudo,97.5) ...
        | mutualInformation(t).(matlab.lang.makeValidName(predictY)).gocue.true < prctile(mutualInformation(t).(matlab.lang.makeValidName(predictY)).gocue.pseudo,2.5);
    pdist = pdist + ptest;
    plot(ew, mutualInformation(t).(matlab.lang.makeValidName(predictY)).gocue.pseudo,'Color',[.65 .65 .65]);
    hold on
    plot(ew, mutualInformation(t).(matlab.lang.makeValidName(predictY)).gocue.true,'Color','k','LineWidth',3)
    line([0 0],[-.1 1],'Color',[.5 .5 .5],'LineStyle','--')
    line([-1 -1],[-.1 1],'Color',[.5 .5 .5],'LineStyle','--')
%     line([-1.5 1.5],[0.5 0.5],'Color','k','LineStyle','-')
    plot(ew,ptest-.15,'Marker','*','LineStyle','none','Color','k')
    catch
    end
    ylim([-.1 1])
    xlim([-1.5 1.5])
    prettyPlot(gca)
    title(expRef,'Interpreter','none')
end

%
for b = 1:16
    binomtest(b) = myBinomTest(pdist(b),36,.05,'two');
end
allPseudo = [];
for t=1:36
    allTrue(t,:) = mutualInformation(t).(matlab.lang.makeValidName(predictY)).gocue.true;
    for p = 1:length(mutualInformation(t).(matlab.lang.makeValidName(predictY)).gocue.pseudo)
        allPseudo(t,p,:) = mutualInformation(t).(matlab.lang.makeValidName(predictY)).gocue.pseudo(p,:);
    end
end

meanTrue = mean(allTrue,1);
meanPseudo = squeeze(mean(allPseudo,1));

figure;
set(gcf,'position',[996 920 639 312])
hold on;
plot(ew,meanPseudo,'Color',[.65 .65 .65]);
plot(ew,meanTrue,'Color','k','LineWidth',3);
line([0 0],[0 1],'Color',[.5 .5 .5],'LineStyle','--')
line([-1 -1],[0 1],'Color',[.5 .5 .5],'LineStyle','--')
plot(ew,(binomtest < .05)-.276,'Marker','*','LineStyle','none','Color','k')
ylim([-.1 .725])
xlim([-1.5 1.5])
prettyPlot(gca)
xlabel('Time from go cue (s)')
ylabel('Mutual information')
% for t=1:16
%     if binomtest(t) <= 0.05 && binomtest(t) > 0.01
%         text(ew(t),0.724,'*','HorizontalAlignment','center')
%     elseif binomtest(t) <= 0.01 && binomtest(t) > 0.001
%         text(ew(t),0.724,'**','HorizontalAlignment','center')
%     elseif binomtest(t) <= 0.001
%         text(ew(t),0.724,'***','HorizontalAlignment','center')
%     end
% end
%
figure;
set(gcf,'position',[996 600 639 312])
bar(ew, pdist/36,'FaceColor','k')
prettyPlot(gca)
line([0 0],[0 1],'Color',[.5 .5 .5],'LineStyle','--')
line([-1 -1],[0 1],'Color',[.5 .5 .5],'LineStyle','--')
hold on
ylabel('Proportion of sessions')
xlim([-1.6 1.6])
ylim([0 1])
xlabel('Time from go cue (s)')
% plot(ew,(binomtest<.05)-.01,'Marker','*','LineStyle','none','Color','k')
line([-1.6 1.6],[0.1389 .1389],'Color','r','LineStyle','-')



