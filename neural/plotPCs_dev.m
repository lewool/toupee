numPCs = 3;
allPCs = getPCs(alignedResps{2});

%%
trials{1} = trialTypes.intVar.cb2D.side_direction{3,1};
trials{2} = trialTypes.intVar.cb2D.side_direction{3,2};

eventWindow_int = linspace(eventWindow(1),eventWindow(end),100);

for t = 1:length(trials)
    for p = 1:numPCs
        pcResps_mean(t,p,:) = interp1(eventWindow, nanmean(allPCs(trials{t},:,p)), eventWindow_int);
        pcResps_sem(t,p,:) = interp1(eventWindow, nanstd(allPCs(trials{t},:,p))/sqrt(length(trials{t})), eventWindow_int);
    end
end

%%
figure;hold on
color = [0 .4 1];
color = [.75 .75 .75];
color = [1 0 0];
for pp = 1:3
subplot(1,3,pp)
hold on;
plotPSTHs(eventWindow_int, pcResps_mean(1,pp,:), pcResps_sem(1,pp,:), color,'-')
plotPSTHs(eventWindow_int, pcResps_mean(2,pp,:), pcResps_sem(2,pp,:), color,':')
line([0 0],[0 3],'LineStyle','--','Color',[0.5 0.5 0.5]);
ylim([0 1.4])

end


%%
figure;hold on
xlim([0 max(max(max(pcResps_mean)))])
ylim([0 .8])
zlim([0 .7])
set(gca,'view',[-90 0])

for pli = 1:length(eventWindow_int)-1
    plot1 = plot3(squeeze(pcResps_mean(1,1,pli:pli+1)),squeeze(pcResps_mean(1,2,pli:pli+1)),squeeze(pcResps_mean(1,3,pli:pli+1)),'Color','b','LineWidth',1);
    hold on
    plot2 = plot3(squeeze(pcResps_mean(2,1,pli:pli+1)),squeeze(pcResps_mean(2,2,pli:pli+1)),squeeze(pcResps_mean(2,3,pli:pli+1)),'Color','r','LineWidth',1);
    hold on
    if pli == 50
        plot3(squeeze(pcResps_mean(1,1,pli)),squeeze(pcResps_mean(1,2,pli)),squeeze(pcResps_mean(1,3,pli)),'ko');
        plot3(squeeze(pcResps_mean(2,1,pli)),squeeze(pcResps_mean(2,2,pli)),squeeze(pcResps_mean(2,3,pli)),'ko');
    end
    pause(.2)
end


