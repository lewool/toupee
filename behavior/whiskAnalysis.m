%datavectors are generated in loopOverSessions.m
%choose the time window 0.2-0 s aligned to stim
%yEarly = mean(earlySessionsWhisk(:,91:101));
yAll = (eyeData(iX).eta.alignedFace{1}(:,95:101,2));
yEarly = earlySessionsWhisk(:,95:101);
yLate = lateSessionsWhisk(:,95:101);
x = 1000*(eyeData(iX).eta.eventWindow(:,95:101));


%calculate the fit for all early sessions and plot whisking + linear fit
for irow = 1:length(yEarly(:,1))
    coefficientsE{irow} = polyfit(x, yEarly(irow,:), 1);
    xFit{irow} = linspace(min(x), max(x), 1000);
    yEarlyFit{irow} = polyval(coefficientsE{irow} , xFit{irow});
    hold on;
    plot(x, yEarly(irow,:),'Color', [0, 0.4470, 0.7410],'LineWidth', 2)
    plot(xFit{irow}, yEarlyFit{irow},'--','Color', [0, 0.75, 0.75], 'LineWidth', 1)
    %gtext('b = 0.191')
    hold on
end


%calculate the fit for all late sessions and plot whisking + linear fit
for irow = 1:length(yLate(:,1))
    coefficientsL{irow} = polyfit(x, yLate(irow,:), 1);
    xFit{irow} = linspace(min(x), max(x), 1000);
    yLateFit{irow} = polyval(coefficientsL{irow} , xFit{irow});
    hold on;
    plot(x, yLate(irow,:),'Color', [77/255, 137/255, 99/255],'LineWidth', 2)
    plot(xFit{irow}, yLateFit{irow},'--',...
        'Color', [105/255, 165/255, 131/255], 'LineWidth', 1)
    xlabel('Time (ms) aligned to stimulus')
    ylabel('Whisking')
    title('LEW031')
    grid off
    legend('Early (s1)', 'Early fit', 'Late (s1)', 'Late fit')
end

%calculate the fit of whisking for all sessions to use for trial indexing 
for irow = 1:length(yAll(:,1))
    coefficientsA{irow} = polyfit(x, yAll(irow,:), 1);
    xFit{irow} = linspace(min(x), max(x), 1000);
    yAllFit{irow} = polyval(coefficientsA{irow} , xFit{irow});
end

%%
    %perform paired ttest
 slopes = zeros(length(yLate(:,1)),2);
for cell = 1:5
    slopes(cell,1)= coefficientsE{cell}(1);
    slopes(cell,2)= coefficientsL{cell}(1);
end
[h,p]= ttest(slopes(:,1),slopes(:,2))

%%
%plot a paired line graph for the slopes of early and late trials 

%boxplot(slopes,'Labels',{'Early', 'Late'},'PlotStyle', 'compact')
%plot([1:2],slopes)
trials = categorical(cellstr('Early'));
trials(end+1) = categorical(cellstr('Late'));
for i = 1:length(slopes(:,:))
    plot(trials,slopes(i,:), 'LineWidth', 3, 'Color',[160/255, 160/255, 160/255])
    ylabel('Slope')
    hold on 
end
plot(trials,mean(slopes),'--','LineWidth', 2, 'Color','black')


%%
%filtering type 1: whisk vs non-whisk
%filter only trials where whisking is significantly INCREASING (pos slope) pre-stim 
%store the trials with increase in whisking in whiskTrials
lm_slopes=(zeros(length(yAll),1));
whiskTrials=[];
nowhiskTrials=[];
for itrial = 1:length(yAll)
    lm = fitlm(x, yAll(itrial,:)); 
    lm_slopes(itrial)=lm.Coefficients.Estimate(2);
    if lm.Coefficients.pValue(2)<0.05
        if lm.Coefficients.Estimate(2)>0
            whiskTrials(end+1) = itrial;
        else
            nowhiskTrials(end+1) = itrial;
        end
    else
        nowhiskTrials(end+1) = itrial;
    end
end
%add whisk and non whisk trials to 1 matrix for  raster plotting 

whichTrials{1,1} = whiskTrials;
whichTrials{1,2} = nowhiskTrials;


%%
[baselineResps, stimResps, pmovResps, movResps, rewResps] = getEpochResps(neuralData.eta);
%c1=baselineResps(whiskTrials,1);
whiskCells=[];
for icell = 1:length(neuralData.eta.alignedResps{1,3})
    [p,h]=ranksum(baselineResps(whiskTrials,icell),baselineResps(nowhiskTrials,icell));
    if h==1
        whiskCells(end+1)= icell;
    end
end
        




