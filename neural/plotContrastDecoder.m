figure;
hold on;
cons = {'con0' 'con0_05' 'con0_12' 'con0_5' 'con1'};
colors = [.8 .8 .8; .6 .6 .6; .4 .4 .4; .2 .2 .2; 0 0 0];

for m = 1:length(mutualInformation)
    mouseName = char(mouseList{m});
    expDate = char(expList{m}{1});
    expNum = expList{m}{2};
    [expRef, ~] = data.constructExpRef(mouseName,expDate,expNum);
    
%     subplot(5,5,m);
    hold on
    for c = 1:length(fieldnames(mutualInformation(m).side.gocue.true))
        plp = plot(ew,...
            mutualInformation(m).side.gocue.pseudo.(matlab.lang.makeValidName(cons{c}))',...
            'Color',[.5 .5 .5],...
            'LineWidth',0.5);
        hold on 
        plot(ew,...
            mutualInformation(m).side.gocue.true.(matlab.lang.makeValidName(cons{c}))',...
            'Color',colors(c,:),...
            'LineWidth',2);
        uistack(plp,'bottom')
        
        line([-1 -1],[-.4 1.5],'Color',[.5 .5 .5],'LineStyle','--')
        line([0 0],[-.4 1.5],'Color',[.5 .5 .5],'LineStyle','--')
        ylim([-.4 1])
        xlim([ew(1) ew(end)])
        prettyPlot(gca)
        xlabel('Time from go cue (s)')
        ylabel('Mutual information')
        title(expRef,'Interpreter','none')
    end
end

%% plot
figure;

hold on;
% mean prediction across all sessions
for c = 1:length(cons)
    clear pb
    for m = 1:length(fieldnames(mutualInformation(m).side.gocue.true))
        pb(m,:) = mutualInformation(m).side.gocue.true.(matlab.lang.makeValidName(cons{c}));
    end
pb = squeeze(pb);
hold on;
plotSignal(ew,median(pb,1),median(pb,1)+std(pb,1)/sqrt(26),median(pb,1)-std(pb,1)/sqrt(26),colors(c,:),'-')
line([0 0],[-.4 1],'Color',[.5 .5 .5],'LineStyle',':')
line([-1 -1],[-.4 1],'Color',[.5 .5 .5],'LineStyle',':')
line([-1.5 3],[0 0],'Color',[.5 .5 .5],'LineStyle',':')
xlim([-1.5 2.8])
prettyPlot(gca)
xlabel('Time from go cue')

ylabel('Mutual Information')
end