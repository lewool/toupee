function plotVelocities(expInfo, behavioralData, blockTag)
    
    if nargin < 3
        blockTag = 'none';
    end

    for ex = 1:size(expInfo,2)
        
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

            clear trialConditions labels condIdx
            contrasts = getUniqueContrasts(expInfo(ex));
            allContrasts = getAllContrasts(expInfo(ex));

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

        d = 1;
        for c = 1:length(contrasts)
            for t = 1:length(trialConditions)
                [~, condIdx{d,:}.all] = selectCondition(expInfo(ex), contrasts(c), behavioralData(ex), trialConditions{t});
                labels{d,1} = strcat(trialLabels{t},num2str(contrasts(c)));
                d = d+1;
            end
        end

        % 
        tLen = length(trialLabels);
        for t = 1:tLen   
            conds = find(contains(labels,trialLabels{t}));
            for iCond = 1:length(conds)
                testIdx = condIdx{conds(iCond)}.all;
                stimOnsets = behavioralData(ex).eventTimes(1).daqTime(testIdx);  
                firstMoves = behavioralData(ex).wheelMoves.epochs(5).onsetTimes(testIdx);
                firstMoveVels = abs(behavioralData(ex).wheelMoves.epochs(5).peakVel(testIdx));
                responseTimes = firstMoves - stimOnsets;
                firstMoveVels(responseTimes > 3) = NaN;
                meanVels(t,iCond,ex) = mean(firstMoveVels);
                semVels(t,iCond,ex) = std(firstMoveVels)/sqrt(length(firstMoveVels));
            end

        end

    end

    % plot

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
        if contains(trialLabels{t},'left')
            subplot(1,2,1)
            hold on;
        else
            subplot(1,2,2)
            hold on;
        end
        
        errorbar(contrasts,nanmean(meanVels(t,:,:),3),nanstd(meanVels(t,:,:),[],3)/sqrt(size(meanVels,3)),'o','MarkerFaceColor',color,'LineStyle','-','color',color,'MarkerEdgeColor','none','capsize',0);
        
    end

    ax1 = subplot(1,2,1);
    ax1.Color = [0 .4 1];
    ax1.Color(4) = 0.1;
    ax1.TickDir = 'out';
%     text(-1,1.8,'left choices')
    xticklabels({'-100' '-50' '0' '50' '100'})
    xlabel('Left (%)          Right (%)');
    ylabel('Velocity (mm/s)')
    line([0 0],[0 25],'linestyle','--','color',[.5 .5 .5])
    axis([-1.2 1.2 2 25])

    ax2 = subplot(1,2,2);
    ax2.Color = [1 0 0];
    ax2.Color(4) = 0.1;
    ax2.TickDir = 'out';
%     text(-1,1.8,'right choices')
    xticklabels({'-100' '-50' '0' '50' '100'})
    xlabel('Left (%)           Right (%)');
    ylabel('Velocity (mm/s)')
    line([0 0],[0 25],'linestyle','--','color',[.5 .5 .5])
    axis([-1.2 1.2 2 25])

