function plotWheelTrajectories(expInfo, behavioralData, blockTag)
    
    if nargin < 3
        blockTag = 'none';
    end
        
    if strcmp(blockTag,'none')

        clear trialConditions labels condIdx
        contrasts = getUniqueContrasts(expInfo(ex));
        allContrasts = getAllContrasts(expInfo(ex));

        % %set up trial conditions for hi-L and hi-R blocks
        trialConditions{1} = initTrialConditions('movementDir','cw','movementTime','early');
        trialConditions{2} = initTrialConditions('movementDir','ccw','movementTime','early');
        trialConditions{3} = initTrialConditions('movementDir','cw','movementTime','late');
        trialConditions{4} = initTrialConditions('movementDir','ccw','movementTime','late');

        trialLabels{1} = 'leftEarly';
        trialLabels{2} = 'rightEarly';
        trialLabels{3} = 'leftLate';
        trialLabels{4} = 'rightLate';

    else

        

        % %set up trial conditions for hi-L and hi-R blocks
        trialConditions{1} = initTrialConditions('highRewardSide','left','movementDir','cw','movementTime','early');
        trialConditions{2} = initTrialConditions('highRewardSide','left','movementDir','ccw','movementTime','early');
        trialConditions{3} = initTrialConditions('highRewardSide','left','movementDir','cw','movementTime','late');
        trialConditions{4} = initTrialConditions('highRewardSide','left','movementDir','ccw','movementTime','late');
        trialConditions{5} = initTrialConditions('highRewardSide','right','movementDir','cw','movementTime','early');
        trialConditions{6} = initTrialConditions('highRewardSide','right','movementDir','ccw','movementTime','early');
        trialConditions{7} = initTrialConditions('highRewardSide','right','movementDir','cw','movementTime','late');
        trialConditions{8} = initTrialConditions('highRewardSide','right','movementDir','ccw','movementTime','late');


        trialLabels{1} = 'leftEarly_hiL';
        trialLabels{2} = 'rightEarly_hiL';
        trialLabels{3} = 'leftLate_hiL';
        trialLabels{4} = 'rightLate_hiL';
        trialLabels{5} = 'leftEarly_hiR';
        trialLabels{6} = 'rightEarly_hiR';
        trialLabels{7} = 'leftLate_hiR';
        trialLabels{8} = 'rightLate_hiR';
    end
    
    for ex = 1:size(expInfo,2)
        clear labels condIdx
        contrasts = getUniqueContrasts(expInfo(ex));

        d = 1;
        for t = 1:length(trialConditions)
            [~, condIdx{d,:}.all] = selectCondition(expInfo(ex), contrasts, behavioralData(ex), trialConditions{t});
            labels{d,1} = strcat(trialLabels{t});
            d = d+1;
        end

        % 
        
        for iTrial = 1:length(behavioralData(ex).eventTimes(1).daqTime)
            rawTraceTimes = behavioralData(ex).wheelMoves.traces.time{iTrial};
            rawTraceValues = behavioralData(ex).wheelMoves.traces.pos{iTrial};

            eventTime = interp1(rawTraceTimes,rawTraceTimes,behavioralData(ex).eventTimes(1).daqTime(iTrial)-.5,'nearest','extrap');
            feedbackTime = interp1(rawTraceTimes,rawTraceTimes,behavioralData(ex).eventTimes(5).daqTime(iTrial),'nearest','extrap');
            eventIdx = find(rawTraceTimes == eventTime);
            feedbackIdx = find(rawTraceTimes == feedbackTime);
            traceTime = rawTraceTimes - eventTime;
            traceValue = rawTraceValues;


            startIdx = find(traceTime == -.5);
            endIdx = find(traceTime == 1.5);
            timeSeries = eventIdx:eventIdx+2000;

            truncValues(iTrial,:) = traceValue(timeSeries);
            

        end
        
        tLen = length(trialLabels);
        for t = 1:tLen   
            testIdx = condIdx{t}.all;
            meanTraces(t,:,ex) = nanmean(truncValues(testIdx,:));
        end
        
        clear truncValues

    end

    % plot
    plotTime = linspace(-.5,1.5,2001);
    figure;hold on;
    set(gcf,'position',[150 580 920 400]);

    for t = 1:tLen
        if contains(trialLabels{t},'hiL')
            color = [.1 .7 .1];
        elseif contains(trialLabels{t},'hiR')
            color = [1 .6 0];
        else
            color = [0 0 0];
        end
        if contains(trialLabels{t},'leftEarly')
            subplot(2,2,1)
            hold on;
        elseif contains(trialLabels{t},'leftLate')
            subplot(2,2,3)
            hold on;
        elseif contains(trialLabels{t},'rightEarly')
            subplot(2,2,2)
            hold on;
        elseif contains(trialLabels{t},'rightLate')
            subplot(2,2,4)
            hold on;
        end
        meanTrace = nanmean(meanTraces(t,:,:),3)';
        semTrace = nanstd(meanTraces(t,:,:),[],3)'/sqrt(size(meanTraces,3));
        plot(plotTime,meanTrace,'MarkerFaceColor',color,'LineStyle','-','LineWidth',1,'color',color,'MarkerEdgeColor','none');
        plotCI = fill([plotTime';flipud(plotTime')],[(meanTrace-semTrace);flipud(meanTrace+semTrace)],color, 'LineStyle', 'none');
        alpha(0.2); 
    end
    
    titles = {'choose left (early)' 'choose right (early)' 'choose left (late)' 'choose right (late)'};
    for t = 1:4
        ax1 = subplot(2,2,t);
        ax1.TickDir = 'out';
        xlabel('Time (s)');
        ylabel('Position (mm)')
        title(titles{t})
        line([0 0],[-50 50],'linestyle','--','color',[.5 .5 .5])
        axis([-0.5 1.5 -50 50])
    end

   

