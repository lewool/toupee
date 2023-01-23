contrastSet = cell(1,length(mouseList));
ppLeft = cell(1,length(mouseList));
ppRight = cell(1,length(mouseList));
pciLeft = cell(1,length(mouseList));
pciRight = cell(1,length(mouseList));

for m = 1:length(mouseList)
    fprintf('%s \n',num2str(m))
    [cc, ppl, ppr, pcil, pcir] = plotPsychometric({mouseList{m}}, {expList{m}});
    contrastSet{m} = cc;
    ppLeft{m} = ppl;
    ppRight{m} = ppr;
    pciLeft{m} = pcil;
    pciRight{m} = pcir;
    close all
end
%%
pplmat = cell2mat(ppLeft');
pprmat = cell2mat(ppRight');
ppDeltaLapse = (mean(pprmat(:,1),2) - mean(pplmat(:,1),2)) + (mean(pprmat(:,end),2) - mean(pplmat(:,end),2));
% pplLapse = (pplmat(:,1)) + (1 - pplmat(:,end));
% pprLapse = (pprmat(:,1)) + (1 - pprmat(:,end));

ppBias = mean(pprmat(:,5),2) - mean(pplmat(:,5),2);
% ppDeltaLapse = pprLapse - pplLapse;

%%
figure;
umn = unique(mouseList);
for u = 1:length(umn)   
    mn = umn{u};
    range = find(strcmp(mouseList,mn));
    meanDB(u) = mean(pprmat(range,5) - pplmat(range,5),1);
    meanDL(u) = mean(pprmat(range,1) - pplmat(range,1) + pprmat(range,end) - pplmat(range,end),1);
    semDB(u) = std(pprmat(range,5) - pplmat(range,5),1)/sqrt(length(range));
    semDL(u) = std(pprmat(range,1) - pplmat(range,1) + pprmat(range,end) - pplmat(range,end),1)/sqrt(length(range));
    subplot(2,3,u)
    hold on
    line([-105 105],[0.5 0.5],'Color',[.5 .5 .5],'LineStyle',':')
    line([0 0],[-.05 1.05],'Color',[.5 .5 .5],'LineStyle',':')
    % plot(contrastSet{1},mean(pprmat(1:14,:),1),'Color',[; 1 .6 0],'LineWidth',1.5)
    errorbar(contrastSet{1},mean(pprmat(range,:),1),std(pprmat(range,:))/sqrt(length(range)),...
        'capsize',0,'Marker','o','MarkerSize',10,'MarkerFaceColor',[1 .6 0],'MarkerEdgeColor','w','Color',[1 .6 0],'LineWidth',1)
    hold on
    % plot(contrastSet{1},mean(pplmat(1:14,:),1),'Color',[0.1 0.7 0.1],'LineWidth',1.5)
    errorbar(contrastSet{1},mean(pplmat(range,:),1),std(pprmat(range,:))/sqrt(length(range)),...
        'capsize',0,'Marker','o','MarkerSize',10,'MarkerFaceColor',[0.1 0.7 0.1],'MarkerEdgeColor','w','Color',[0.1 0.7 0.1],'LineWidth',1)

    prettyPlot(gca)
    xlim([-105 105])
    ylim([-.05 1.05])
    xticks([-100 -50 0 50 100])
    yticks([0 .5 1])
    xlabel('Contrast (%)')
    ylabel('Proportion of right choices')
    title(strcat(mn,{' '},'mean performance'))
end

%%
figure
lim1 = -.2;
lim2 = .8;
hold on;
line([0 0],[-1 1],'Color',[.25 .25 .25],'LineStyle','--')
line([-1 1],[0 0],'Color',[.25 .25 .25],'LineStyle','--')
scatter(ppBias,ppDeltaLapse,100,[.75 .75 .75],'.')
for u = 1:length(umn)
    
    errorbar(meanDB(u),meanDL(u),meanDL(u)-semDL(u),meanDL(u)+semDL(u),meanDB(u)-semDB(u),meanDB(u)+semDB(u),'Marker','.','MarkerSize',20,'capsize',0)
    hold on
end
prettyPlot(gca)
xlim([lim1 lim2])
ylim([lim1 lim2])
xticks([lim1:.2:lim2])
yticks([lim1:.2:lim2])
axis square
xlabel('\Delta Bias')
ylabel('\Delta Lapse')
