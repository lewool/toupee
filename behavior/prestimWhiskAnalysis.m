function prestimWhiskAnalysis(eyeData,earlyTrialsWhisk,lateTrialsWhisk)

%datavectors are generated in loopOverTrials.m
%time window 120-0 ms aligned to stim

x = 1000*(eyeData(1).eta.eventWindow(:,95:101));
yEarly = earlyTrialsWhisk(:,95:101);
yLate = lateTrialsWhisk(:,95:101);

%%
%calculate the fit for all early Trials and plot whisking + linear fit
for irow = 1:length(yEarly(1,1))
    coefficientsE{irow} = polyfit(x, yEarly(irow,:), 1);
    xFit{irow} = linspace(min(x), max(x), 1000);
    yEarlyFit{irow} = polyval(coefficientsE{irow} , xFit{irow});
    hold on;
    plot(x, yEarly(irow,:),'Color', [255/255, 128/255, 0/255],'LineWidth', 2)
    plot(xFit{irow}, yEarlyFit{irow},'--','Color', [255/255, 128/255, 0/255], 'LineWidth', 1)
    %gtext('b = 0.191')
    hold on
end


%calculate the fit for all late Trials and plot whisking + linear fit
for irow = 1:length(yLate(1,1))
    coefficientsL{irow} = polyfit(x, yLate(irow,:), 1);
    xFit{irow} = linspace(min(x), max(x), 1000);
    yLateFit{irow} = polyval(coefficientsL{irow} , xFit{irow});
    hold on;
    plot(x, yLate(irow,:),'Color', [0/255, 153/255, 153/255],'LineWidth', 2)
    plot(xFit{irow}, yLateFit{irow},'--',...
        'Color', [0/255, 153/255, 153/255], 'LineWidth', 1)
    xlabel('Time (ms) aligned to stimulus')
    ylabel('Whisking')
    title('LEW031')
    grid off
    %legend('Early (s1)', 'Early fit', 'Late (s1)', 'Late fit')
end

%%
%perform paired ttest between slopes
 slopes = zeros(length(yEarly(:,1)),2);
for sesh = 1:length(yEarly(:,1))
    slopes(sesh,1)= coefficientsE{sesh}(1);
    slopes(sesh,2)= coefficientsL{sesh}(1);
end
[h,p]= ttest(slopes(:,1),slopes(:,2))

%plot a paired line graph for the slopes of early and late trials 
trials = categorical(cellstr('Early'));
trials(end+1) = categorical(cellstr('Late'));
for i = 1:length(slopes(:,:))
    plot(trials,slopes(i,:), 'LineWidth', 3, 'Color',[160/255, 160/255, 160/255])
    ylabel('Slope')
    hold on 
end
plot(trials,mean(slopes),'--','LineWidth', 2, 'Color','black')

%% create data for quantified mean paired line plot 
yEarly = earlyTrialsWhisk(:,95:101);
yLate = lateTrialsWhisk(:,95:101);

%calculate the session mean for all early & late trials 
meanWhisk = zeros(length(yLate(:,1)),2);
for sesh = 1:length(yEarly(:,1))
    meanWhisk(sesh,1)= mean(yEarly(sesh,:));
    meanWhisk(sesh,2)= mean(yLate(sesh,:));
end

%perform statistical test
[h,p]= ttest(meanWhisk(:,1),meanWhisk(:,2))

%plot a paired line graph for the means of early and late trials 
trials = categorical(cellstr('Early'));
trials(end+1) = categorical(cellstr('Late'));
for i = 1:length(meanWhisk(:,:))
    plot(trials,meanWhisk(i,:), 'LineWidth', 3, 'Color',[160/255, 160/255, 160/255])
    ylabel('Mean pre-stim whisking')
    hold on 
end
plot(trials,mean(meanWhisk),'--','LineWidth', 2, 'Color','black')
