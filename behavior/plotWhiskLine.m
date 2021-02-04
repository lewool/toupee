%datavectors are generated in loopOverSessions.m
%y is a mean of the whisking in all sessions for the chosen animal
%choose the time window 0.2-0 s aligned to stim
%yEarly = mean(earlySessionsWhisk(:,91:101));
yEarly = earlySessionsWhisk(:,95:101);
yLate = lateSessionsWhisk(:,95:101);
x = 1000*(eyeData(iX).eta.eventWindow(:,95:101));


%calculate the fit for early sessions and plot whisking + linear fit
for irow = 1:length(yEarly(:,1))
    coefficientsE{irow} = polyfit(x, yEarly(irow,:), 1);
    xFit{irow} = linspace(min(x), max(x), 1000);
    yEarlyFit{irow} = polyval(coefficientsE{irow} , xFit{irow});
    hold on;
    plot(x, yEarly(irow,:),'Color', [0, 0.4470, 0.7410],'LineWidth', 1)
    plot(xFit{irow}, yEarlyFit{irow},'--','Color', [0, 0.75, 0.75], 'LineWidth', 1)
    %gtext('b = 0.191')
    hold on
end


%calculate the fit for late sessions and plot  whisking + linear fit
for irow = 1:length(yLate(:,1))
    coefficientsL{irow} = polyfit(x, yLate(irow,:), 1);
    xFit{irow} = linspace(min(x), max(x), 1000);
    yLateFit{irow} = polyval(coefficientsL{irow} , xFit{irow});
    hold on;
    plot(x, yLate(irow,:),'Color',[0.9290, 0.6940, 0.1250],'LineWidth', 1)
    plot(xFit{irow}, yLateFit{irow},'--','Color',[0.8500 0.3250 0.0980], 'LineWidth', 1)
    xlabel('Time (ms) aligned to stimulus')
    ylabel('Whisking')
    title('LEW031')
    grid off
end

%%
    %perform paired ttest
 slopeE = [];
 slopeL = [];
for cell = 1:5
    slopeE(end+1)= coefficientsE{cell}(1);
    slopeL(end+1)= coefficientsL{cell}(1);
end
[h,p]= ttest(slopeE,slopeL)
