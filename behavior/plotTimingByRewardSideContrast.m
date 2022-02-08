function plotTimingByRewardSideContrast(expInfo, behavioralData)

    %set up trial conditions for hi-L and hi-R blocks
    trialConditions{1} = initTrialConditions('highRewardSide','left','movementTime','early');
    trialConditions{2} = initTrialConditions('highRewardSide','left','movementTime','late');
    trialConditions{3} = initTrialConditions('highRewardSide','right','movementTime','early');
    trialConditions{4} = initTrialConditions('highRewardSide','right','movementTime','late');


    trialLabels{1} = 'early_hiL';
    trialLabels{2} = 'late_hiL';
    trialLabels{3} = 'early_hiR';
    trialLabels{4} = 'late_hiR';

    for ex = 1:size(expInfo,2)
        
        nt = numel(expInfo(ex).block.events.endTrialValues);

        %absolute visual contrast
        V = [abs(expInfo(ex).block.events.contrastValues(1:nt))];

        % high reward side?
        H = [sign(expInfo(ex).block.events.contrastValues(1:nt) .*...
            expInfo(ex).block.events.highRewardSideValues(1:nt))];

        %outcome
        Y = [behavioralData(ex).wheelMoves.epochs(2).isMoving]';
        
        X = V'.*H';
        
        uniqueX = unique(X);
        for u = 1:length(uniqueX)
            ybar(u,ex) = mean(Y(X == uniqueX(u)));
%             semY(u) = std(Y(X == uniqueX(u)))/sqrt(length(Y(X == uniqueX(u))));
        end
        
        
        
    end
    
    meanY = mean(ybar,2);
    semY = std(ybar,[],2)/sqrt(size(ybar,2));
   
    %%
    figure;
    lineColor = 'k';
    errorbar(uniqueX, meanY, semY, ...
        'k', 'Color', lineColor,'MarkerFaceColor',lineColor,'MarkerEdgeColor','w', 'Marker','o','MarkerSize', 8, ...
        'LineWidth', 1,'capsize',0)
    hold on;
    box off
    
     xlim([-1.2 1.2])
    set(gca, 'XTick', [-1, -.5 -.25, 0, .25, .50, 1.00])
    set(gca, 'XTickLabels', {'100', '50', '25', '0', '25', '50', '100'})
    xlabel('Contrast (%), low rew.              Contrast (%), high rew.')
    ylabel('Prop. impulsive moves')
    set(gca,'tickdir','out')