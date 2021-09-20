function plotMoveTypeHistory(expInfo, behavioralData, blockTag)
    
    if nargin < 3
        blockTag = 'none';
    end
    
    if strcmp(blockTag,'none')

        clear trialConditions labels condIdx
%         contrasts = getUniqueContrasts(expInfo(ex));
%         allContrasts = getAllContrasts(expInfo(ex));

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

        %set up trial conditions for hi-L and hi-R blocks
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
    
    % plot proportion of early vs late as a function of trial number
%     allSessionLengths = [];
%     for ex = 1:size(expInfo,2)
%         sessionLength = length(behavioralData(ex).wheelMoves.epochs(3).onsetTimes);
%         allSessionLengths = [allSessionLengths sessionLength];
%     end
%     
%     binEdges = linspace(1,max(allSessionLengths),21);
%     
%     
%     for ex = 1:size(expInfo,2)
%         earlyTrials = find(~isnan(behavioralData(ex).wheelMoves.epochs(2).onsetTimes));
%         lateTrials = find(~isnan(behavioralData(ex).wheelMoves.epochs(3).onsetTimes) & isnan(behavioralData(ex).wheelMoves.epochs(2).onsetTimes));
%         earlyHist = histcounts(earlyTrials,binEdges); 
%         lateHist = histcounts(lateTrials,binEdges);
%         sumHist = earlyHist + lateHist;
%         earlyHists(ex,:) = earlyHist;
%         lateHists(ex,:) = lateHist;
%     end
%     
%     meanEarly = nanmean(earlyHists);
%     semEarly = std(earlyHists)/sqrt(sum(earlyHists>0));
%     meanLate = nanmean(lateHists);
%     semLate = std(lateHists)/sqrt(sum(lateHists>0));
%     meanSum = meanEarly + meanLate;
%     meanEarly_norm = meanEarly./meanSum;
%     semEarly_norm = semEarly./meanSum;
%     meanLate_norm = meanLate./meanSum;
%     semLate_norm = semLate./meanSum;
%     
%     figure;
%     hold on
% %     plot(binEdges(1:end-1),meanLate_norm,'c','linewidth',2);
% %     plotCI = fill([binEdges(1:end-1)';flipud(binEdges(1:end-1)')],[(meanLate_norm'-semLate_norm');flipud(meanLate_norm'+semLate_norm')],'c', 'LineStyle', 'none');
% %     alpha(0.2); 
%     plot(binEdges(1:end-1),meanEarly_norm,'m','linewidth',2);
%     plotCI = fill([binEdges(1:end-1)';flipud(binEdges(1:end-1)')],[(meanEarly_norm'-semEarly_norm');flipud(meanEarly_norm'+semEarly_norm')],'m', 'LineStyle', 'none');
%     alpha(0.2); 
%     
%     ax1 = gca;
%     ax1.TickDir = 'out';
%     xlabel('Trial number');
%     ylabel('Proportion choice type')
%     axis([1 700 0 1])
    
    % compare proportion of early vs late across trial conditions

    for ex = 1:size(expInfo,2)
        clear labels condIdx
        contrasts = getUniqueContrasts(expInfo(ex));

        d = 1;
        for c = 1:length(contrasts)
        for t = 1:length(trialConditions)
            [~, condIdx{d,:}.all] = selectCondition(expInfo(ex), contrasts(c), behavioralData(ex), trialConditions{t});
            labels{d,1} = strcat(trialLabels{t},'_',num2str(contrasts(c)));
            d = d+1;
        end
        end
        
        tLen = length(trialLabels);
        for t = 1:tLen   
            conds = find(contains(labels,trialLabels{t}));
            for iCond = 1:length(conds)
                testIdx = condIdx{conds(iCond)}.all;
                numMoves(t,iCond,ex) = length(testIdx);
            end
        end
    end
    
    meanMoves = squeeze(nanmean(numMoves,3));
    semMoves = nanstd(numMoves,[],3)/sqrt(size(numMoves,3));
    
    normMean(1,:) = meanMoves(1,:)./(meanMoves(1,:)+meanMoves(3,:));
    normSem(1,:) = semMoves(1,:)./(meanMoves(1,:)+meanMoves(3,:));
    
    normMean(3,:) = meanMoves(2,:)./(meanMoves(2,:)+meanMoves(4,:));
    normSem(3,:) = semMoves(2,:)./(meanMoves(2,:)+meanMoves(4,:));
    
    normMean(2,:) = meanMoves(5,:)./(meanMoves(5,:)+meanMoves(7,:));
    normSem(2,:) = semMoves(5,:)./(meanMoves(5,:)+meanMoves(7,:));
    
    normMean(4,:) = meanMoves(6,:)./(meanMoves(6,:)+meanMoves(8,:));
    normSem(4,:) = semMoves(6,:)./(meanMoves(6,:)+meanMoves(8,:));
    
    figure;
    for tt = 1:size(normMean,1)
        if mod(tt,2) == 1
            color = [.1 .7 .1];
        else
            color = [1 .6 0];
        end
        if tt < 3
            subplot(1,2,1);
            hold on;
        else
            subplot(1,2,2);
            hold on;
        end
        errorbar(contrasts,normMean(tt,:),normSem(tt,:),'ko', 'Color', color,'MarkerFaceColor',color,'MarkerEdgeColor','none', 'MarkerSize', 6, ...
            'LineWidth', .5,'LineStyle','-','capsize',0)
    end
        
    ax1 = subplot(1,2,1);
    ax1.Color = [0 .4 1];
    ax1.Color(4) = 0.1;
    ax1.TickDir = 'out';
%     text(-1,1.8,'left choices')
    xticklabels({'-100' '-50' '0' '50' '100'})
    xlabel('Left (%)          Right (%)');
    ylabel('Proportion early choices')
    line([0 0],[0 25],'linestyle','--','color',[.5 .5 .5])
    axis([-1.2 1.2 0 1])

    ax2 = subplot(1,2,2);
    ax2.Color = [1 0 0];
    ax2.Color(4) = 0.1;
    ax2.TickDir = 'out';
%     text(-1,1.8,'right choices')
    xticklabels({'-100' '-50' '0' '50' '100'})
    xlabel('Left (%)           Right (%)');
    ylabel('Proportion early choices')
    line([0 0],[0 25],'linestyle','--','color',[.5 .5 .5])
    axis([-1.2 1.2 0 1])
        
    
        
        