for m = 1:length(mouseList)
    
    mouseName = char(mouseList{m});
    expDate = char(expList{m}{1});
    expNum = expList{m}{2};
    hemisphere = hemList(m);
    [expRef, ~] = data.constructExpRef(mouseName,expDate,expNum);
    fprintf('Loading %s...\n',expRef)
    [behavioralData, expInfo] = data.loadBehavioralDataset(mouseName, expDate, expNum);
    

%%

    et = behavioralData;
    contrasts = getUniqueContrasts(expInfo);
    nt = length(et.eventTimes(1).daqTime);

    trialCorrectChoice = expInfo.block.events.correctResponseValues(1:nt);
%     trueStimuli = expInfo.block.events.contrastValues(1:nt);
    sidedStimuli = expInfo.block.events.contrastValues(1:nt);
    sidedStimuli(sidedStimuli == 0) = eps;
    sidedStimuli(abs(sidedStimuli) < .05) = ...
        sidedStimuli(abs(sidedStimuli) < .05).* ...
        trialCorrectChoice(abs(sidedStimuli) < .05);
    trueChoices = et.wheelMoves.epochs(5).moveDir;
    trueBlock = expInfo.block.events.highRewardSideValues(1:nt);
    trueFeedback = zeros(1,nt);
    trueFeedback(trueChoices .* trialCorrectChoice > 0) = 1;
    trueFeedback(trueChoices .* trialCorrectChoice < 0) = 0;
    trueValue(sidedStimuli .* trueBlock > 0) = 1;
    trueValue(sidedStimuli .* trueBlock < 0) = 0;
    
    [impTrials, ~] = selectCondition(expInfo, contrasts, et, ...
        initTrialConditions('movementTime','early'));

    maxVels = behavioralData.wheelMoves.epochs(5).peakVel;
    RTs = behavioralData.wheelMoves.epochs(5).onsetTimes - behavioralData.eventTimes(1).daqTime;
    trials = intersect(...
        find((RTs < prctile(RTs,97.5)) & (RTs > prctile(RTs,2.5))),...
        find(~isnan(trueChoices)));
    
    fitData{m,1} = impTrials(trials)';
    fitData{m,2} = maxVels(trials)';
    fitData{m,3} = RTs(trials)';
    fitData{m,4} = sidedStimuli(trials)';
    fitData{m,5} = trueChoices(trials)';
    fitData{m,6} = trueBlock(trials)';
    fitData{m,7} = trueFeedback(trials)';
    fitData{m,8} = trueValue(trials)';


    clearvars -except mouseList expList hemList fitData
    
end

%%
fitArray = cell2mat(fitData);
contrasts = [1 .5 .12 .05 eps];
values = [0 1];
blocks = [-1 1];
feed = [0 1];
choices = [-1 1];
colors = [0 0 0; .2 .2 .2; .4 .4 .4; .6 .6 .6; .8 .8 .8];
vcolors = [0 .5 1; 0 0 1];
bcolors = [.1 .7 0; 1 .6 0];
ccolors = [0 .4 1; 1 0 0];


figure;
subplot(6,2,1);
h = histogram(abs(fitArray(:,3)),linspace(0,2,41),'FaceColor',[.9 .9 .9]);
ylabel('No. trials')
hold on
line([.8 .8],[0 max(h.Values)*1.1],'Color',[.5 .5 .5],'LineStyle',':')
xlim([0 2])
prettyPlot(gca)
ylim([0 max(h.Values)])

subplot(6,2,2);
histogram(abs(fitArray(:,2)),linspace(0,300,41),'FaceColor',[.9 .9 .9])
hold on
xlim([0 300])
prettyPlot(gca)

subplot(6,2,3);
hold on
for c = 1:length(contrasts)
    [fout0,xout,~,~] = ksdensity(fitArray(abs(fitArray(:,4)) == contrasts(c),3),[0:.05:2],'function','pdf');
    fm(c) = max(fout0);
    plot(xout,fout0,'color',colors(c,:),'LineWidth',2)
end
line([.8 .8],[0 max(fm)*1.1],'Color',[.5 .5 .5],'LineStyle',':')
legend({'100%' '50' '12' '5' '0'},'location','ne')
legend('boxoff')
yticks([])
ylabel('PDF')
ylim([0 max(fm)*1.1])
prettyPlot(gca)

subplot(6,2,4);
hold on
for c = 1:length(contrasts)
    [fout0,xout,~,~] = ksdensity(abs(fitArray(abs(fitArray(:,4)) == contrasts(c),2)),[0:.1:300],'function','pdf');
    plot(xout,fout0,'color',colors(c,:),'LineWidth',2)
end
prettyPlot(gca)
yticks([])
xticks([0 100 200 300])

subplot(6,2,5);
hold on
for c = 1:length(values)
    [fout0,xout,~,~] = ksdensity(fitArray(abs(fitArray(:,8)) == values(c),3),[0:.05:2],'function','pdf');
    plot(xout,fout0,'color',vcolors(c,:),'LineWidth',2)
    fm(c) = max(fout0);
end
prettyPlot(gca)
line([.8 .8],[0 max(fm)*1.1],'Color',[.5 .5 .5],'LineStyle',':')

xlabel('Reaction time (s)')
legend({'Low value' 'High value'},'location','ne')
legend('boxoff')
yticks([])
ylabel('PDF')
ylim([0 max(fm)*1.1])

subplot(6,2,6);
hold on
for c = 1:length(values)
    [fout0,xout,~,~] = ksdensity(abs(fitArray(abs(fitArray(:,8)) == values(c),2)),[0:.1:300],'function','pdf');
    plot(xout,fout0,'color',vcolors(c,:),'LineWidth',2)
end
prettyPlot(gca)
xlabel('Peak wheel velocity (mm/s)')
yticks([])
xticks([0 100 200 300])

subplot(6,2,7);
hold on
for c = 1:length(blocks)
    [fout0,xout,~,~] = ksdensity(fitArray(fitArray(:,6) == blocks(c),3),[0:.05:2],'function','pdf');
    plot(xout,fout0,'color',bcolors(c,:),'LineWidth',2)
    fm(c) = max(fout0);
end
prettyPlot(gca)
line([.8 .8],[0 max(fm)*1.1],'Color',[.5 .5 .5],'LineStyle',':')

xlabel('Reaction time (s)')
legend({'Left block' 'Right block'},'location','ne')
legend('boxoff')
yticks([])
ylabel('PDF')
ylim([0 max(fm)*1.1])

subplot(6,2,8);
hold on
for c = 1:length(blocks)
    [fout0,xout,~,~] = ksdensity(abs(fitArray(fitArray(:,6) == blocks(c),2)),[0:.1:300],'function','pdf');
    plot(xout,fout0,'color',bcolors(c,:),'LineWidth',2)
end
prettyPlot(gca)
xlabel('Peak wheel velocity (mm/s)')
yticks([])
xticks([0 100 200 300])

subplot(6,2,9);
hold on
for c = 1:length(choices)
    [fout0,xout,~,~] = ksdensity(fitArray(fitArray(:,5) == choices(c),3),[0:.05:2],'function','pdf');
    plot(xout,fout0,'color',ccolors(c,:),'LineWidth',2)
    fm(c) = max(fout0);
end
prettyPlot(gca)
line([.8 .8],[0 max(fm)*1.1],'Color',[.5 .5 .5],'LineStyle',':')

xlabel('Reaction time (s)')
legend({'Left choice' 'Right choice'},'location','ne')
legend('boxoff')
yticks([])
ylabel('PDF')
ylim([0 max(fm)*1.1])

subplot(6,2,10);
hold on
for c = 1:length(choices)
    [fout0,xout,~,~] = ksdensity(abs(fitArray(fitArray(:,5) == choices(c),2)),[0:.1:300],'function','pdf');
    plot(xout,fout0,'color',ccolors(c,:),'LineWidth',2)
end
prettyPlot(gca)
xlabel('Peak wheel velocity (mm/s)')
yticks([])
xticks([0 100 200 300])

subplot(6,2,11);
hold on
for c = 1:length(feed)
    [fout0,xout,~,~] = ksdensity(fitArray(fitArray(:,7) == feed(c),3),[0:.05:2],'function','pdf');
    plot(xout,fout0,'color',ccolors(c,:),'LineWidth',2)
    fm(c) = max(fout0);
end
prettyPlot(gca)
line([.8 .8],[0 max(fm)*1.1],'Color',[.5 .5 .5],'LineStyle',':')

xlabel('Reaction time (s)')
legend({'Left choice' 'Right choice'},'location','ne')
legend('boxoff')
yticks([])
ylabel('PDF')
ylim([0 max(fm)*1.1])

subplot(6,2,12);
hold on
for c = 1:length(blocks)
    [fout0,xout,~,~] = ksdensity(abs(fitArray(fitArray(:,7) == feed(c),2)),[0:.1:300],'function','pdf');
    plot(xout,fout0,'color',ccolors(c,:),'LineWidth',2)
end
prettyPlot(gca)
xlabel('Peak wheel velocity (mm/s)')
yticks([])
xticks([0 100 200 300])
