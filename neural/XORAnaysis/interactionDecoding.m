try
    parpool();
catch
end

figure;

for m = 1:length(mouseList)
    
    %load data
    mouseName = char(mouseList{m});
    expDate = char(expList{m}{1});
    expNum = expList{m}{2};
    [expRef, ~] = data.constructExpRef(mouseName,expDate,expNum);
    fprintf('Loading %s...\n',expRef)
    [behavioralData, expInfo, neuralData] = data.loadDataset(mouseName, expDate, expNum);
    expInfo.hemisphere = hemList(m);
    
    [baselineResps, stimResps, pmovResps, movResps, rewResps, preCueResps, postCueResps] = getEpochResps(neuralData.eta);

    %% set up X and Y

    nt = length(behavioralData.eventTimes(1).daqTime);

    % extract stimulus, choice, feedback, value, and block values
    trueStimuli = expInfo.block.events.contrastValues(1:nt);
    trialCorrectChoice = expInfo.block.events.correctResponseValues(1:nt);

    trueChoices = expInfo.block.events.responseValues(1:nt);
    trueFeedback = -1 + 2*double(expInfo.block.events.feedbackValues(1:nt));

    trueStimEPS = trueStimuli;
    trueStimEPS(trueStimEPS == 0) = eps;
    trueStimEPS(abs(trueStimEPS) < .05) = ...
        trueStimEPS(abs(trueStimEPS) < .05).* trialCorrectChoice(abs(trueStimEPS) < .05);
    trueSide = sign(trueStimEPS);

    labels = [];
    for n = 1:nt
        if trueChoices(n) == -1 && trueFeedback(n) == 1
            labels(n) = 1; %left correct
        elseif trueChoices(n) == -1 && trueFeedback(n) == -1
            labels(n) = 2; %left incorrect
        elseif trueChoices(n) == 1 && trueFeedback(n) == -1
            labels(n) = 3; %right incorrect
        elseif trueChoices(n) == 1 && trueFeedback(n) == 1
            labels(n) = 4; %right correct
        end
    end
    
    [~, whichTrials] = selectCondition(expInfo, getUniqueContrasts(expInfo), behavioralData, ...
        initTrialConditions('repeatType','random','movementTime','late','specificRTs',[.8 3]));
    
    trainTrials = find(trueStimuli ~= 0);
    testTrials = find(trueStimuli == 0);
    
    X_train = rewResps(intersect(whichTrials,trainTrials),:);
    X_test = rewResps(intersect(whichTrials,testTrials),:);

    Y_train = ...
        [trueSide(intersect(whichTrials,trainTrials))' ...
         trueChoices(intersect(whichTrials,trainTrials))' ...
         trueFeedback(intersect(whichTrials,trainTrials))'];
    Y_test = ...
        [trueSide(intersect(whichTrials,testTrials))' ...
         trueChoices(intersect(whichTrials,testTrials))' ...
         trueFeedback(intersect(whichTrials,testTrials))'];

    %keep only 0% labels
    labels = labels(intersect(whichTrials,testTrials))';
    
    %% set up X and Y for linear shift test
    trimLength = 20;
    np = 2*trimLength;
    for l = 1:trimLength*2+1
        ss = l;
        es = length(Y_train) - (trimLength*2-l+1);
        Y_train_shifted(:,:,l) = Y_train(ss:es,:);
    end
    Y_train_shifted(:,:,trimLength+1) = [];
    X_train_shifted = X_train(trimLength+1:end-trimLength,:);
    
    %% fit true
    
    fit_type = 'linear';
    fit_type = 'logistic';
    
    if strcmp(fit_type,'linear')
        options.alpha = 1;
        options.nlambda = 20;
        options.standardize = 'true';
        nf = 3;
        family = 'gaussian';

        for d = 1:size(Y_train,2)
            fit{d} = cvglmnet(X_train, Y_train(:,d), family, options, 'deviance',nf,[],true);
            Y_pred(:,d) = cvglmnetPredict(fit{d}, X_test,'lambda_min');
        end
    elseif strcmp(fit_type,'logistic')
        Y_train = (Y_train+1)/2;
        options.alpha = 1;
        options.nlambda = 20;
        options.standardize = 'true';
        nf = 3;
        family = 'binomial';

        for d = 1:size(Y_train,2)
            fit{d} = cvglmnet(X_train, Y_train(:,d), family, options, 'deviance',nf,[],true);
            Y_pred(:,d) = cvglmnetPredict(fit{d}, X_test,'lambda_min');
        end
    end
    
    % linear shift
    for l = 1:trimLength*2+1
        ss = l;
        es = length(Y_pred) - (trimLength*2-l+1);
        Y_pred_shifted(:,:,l) = Y_pred(ss:es,:);
    end
    Y_pred_shifted(:,:,trimLength+1) = [];
    Y_test_shifted = Y_test(trimLength+1:end-trimLength,:);
    
    all_Y_pred{m} = Y_pred;
    all_Y_test{m} = Y_test;
    all_labels{m} = labels;
    
    %% fit pseudo
    
    if strcmp(fit_type,'linear')
        options.alpha = 1;
        options.nlambda = 20;
        options.standardize = 'true';
        nf = 3;
        family = 'gaussian';
        prog = 0;
        fprintf(1,'fitting pseudos: %3d%%\n',prog);
        for p = 1:np
            for d = 1:size(Y_train,2)
                fit{d,p} = cvglmnet(X_train_shifted, Y_train_shifted(:,d,p), family, options, 'deviance',nf,[],true);
                Y_pred_pseudo(:,d,p) = cvglmnetPredict(fit{d,p}, X_test,'lambda_min');
            end
            prog = (100*(p/np));
            fprintf(1,'\b\b\b\b%3.0f%%',prog);
        end
        fprintf('\n');       
    elseif strcmp(fit_type,'logistic')
        options.alpha = 1;
        options.nlambda = 20;
        options.standardize = 'true';
        nf = 3;
        family = 'binomial';
        prog = 0;
        fprintf(1,'fitting pseudos: %3d%%\n',prog);
        for p = 1:np
            for d = 1:size(Y_train,2)
                fit{d,p} = cvglmnet(X_train_shifted, Y_train_shifted(:,d,p), family, options, 'deviance',nf,[],true);
                Y_pred_pseudo(:,d,p) = cvglmnetPredict(fit{d,p}, X_test,'lambda_min');
            end
            prog = (100*(p/np));
            fprintf(1,'\b\b\b\b%3.0f%%',prog);
        end
    end

    %% plot

    for l = 1:length(unique(labels))
        pred_means(l,:) = mean(Y_pred(labels==l,:),1);
        pred_stds(l,:) = std(Y_pred(labels==l,:),[],1);
    end

    subplot(6,6,m)
    scatter3(Y_pred(labels==1,1), Y_pred(labels==1,2), Y_pred(labels==1,3),...
        'MarkerFaceColor',[.1 .7 .1],'MarkerFaceAlpha',1,'MarkerEdgeColor','none')
    hold on;
    scatter3(Y_pred(labels==2,1), Y_pred(labels==2,2), Y_pred(labels==2,3),...
        'MarkerFaceColor',[.75 0 0],'MarkerFaceAlpha',1,'MarkerEdgeColor','none')
    scatter3(Y_pred(labels==3,1), Y_pred(labels==3,2), Y_pred(labels==3,3),...
        'MarkerFaceColor','w','MarkerEdgeAlpha',1,'MarkerEdgeColor',[.75 0 0])
    scatter3(Y_pred(labels==4,1), Y_pred(labels==4,2), Y_pred(labels==4,3),...
        'MarkerFaceColor','w','MarkerEdgeAlpha',1,'MarkerEdgeColor',[.1 .7 .1])

    ax_lim(1) = min(min(Y_pred))*1.1;
    ax_lim(2) = max(max(Y_pred))*1.1;
    xlim([ax_lim(1) ax_lim(2)])
    ylim([ax_lim(1) ax_lim(2)])
    zlim([ax_lim(1) ax_lim(2)])
    
    %tetrahedron
    line_pairs = nchoosek([1:max(labels)],2);
    for f = 1:size(line_pairs,1)
        line(...
            [pred_means(line_pairs(f,1),1) pred_means(line_pairs(f,2),1)],...
            [pred_means(line_pairs(f,1),2) pred_means(line_pairs(f,2),2)],...
            [pred_means(line_pairs(f,1),3) pred_means(line_pairs(f,2),3)],'Color','k')
    end
    
    %shadow
    for f = 1:size(line_pairs,1)
        line(...
            [pred_means(line_pairs(f,1),1) pred_means(line_pairs(f,2),1)],...
            [pred_means(line_pairs(f,1),2) pred_means(line_pairs(f,2),2)],...
            [ax_lim(1) ax_lim(1)],'Color',[.5 .5 .5])
    end
    
    %lines to floor
    for f = 1:size(pred_means,1)
        line(...
            [pred_means(f,1) pred_means(f,1)],...
            [pred_means(f,2) pred_means(f,2)],...
            [pred_means(f,3) ax_lim(1)],'Color',[.5 .5 .5])
    end

    xlabel('stim')
    ylabel('choice')
    zlabel('outcome')
    view(126,23)
    
    %% fetch summary stats
    
    for f = 1:max(labels)
        label_accuracy(m,f,:) = ...
            sum(Y_test(labels==f,:) == ...
            sign(Y_pred(labels==f,:)))/length(Y_pred(labels==f,:));
    end
    
    all_pred_means(:,:,m) = pred_means;
    total_accuracy(m,:) = sum(Y_test == sign(Y_pred))/size(Y_pred,1);
    
    %% null stats

%     for p = 1:np
%         shift_accuracy(m,p,:) = sum(Y_test_shifted == sign(Y_pred_shifted(:,:,p)))/length(Y_pred_shifted);
%     end
    
    for p = 1:np
        shift_accuracy(m,p,:) = sum(Y_test == sign(Y_pred_pseudo(:,:,p)))/size(Y_pred_pseudo,1);
    end    
    
    %% reset
    clearvars -except mouseList expList hemList m all_pred_means label_accuracy total_accuracy shift_accuracy all_Y_pred all_Y_test all_labels
end

%% plot example session
mm=1;
figure;
scatter3(...
    all_Y_pred{mm}(all_labels{mm}==1,1), ...
    all_Y_pred{mm}(all_labels{mm}==1,2),...
    all_Y_pred{mm}(all_labels{mm}==1,3),...
    'MarkerFaceColor',[.1 .7 .1],'MarkerFaceAlpha',1,'MarkerEdgeColor','none')
hold on;
scatter3(...
    all_Y_pred{mm}(all_labels{mm}==2,1), ...
    all_Y_pred{mm}(all_labels{mm}==2,2),...
    all_Y_pred{mm}(all_labels{mm}==2,3),...
    'MarkerFaceColor',[.75 0 0],'MarkerFaceAlpha',1,'MarkerEdgeColor','none')
scatter3(...
    all_Y_pred{mm}(all_labels{mm}==3,1), ...
    all_Y_pred{mm}(all_labels{mm}==3,2),...
    all_Y_pred{mm}(all_labels{mm}==3,3),...
    'MarkerFaceColor','w','MarkerEdgeAlpha',1,'MarkerEdgeColor',[.75 0 0])
scatter3(...
    all_Y_pred{mm}(all_labels{mm}==4,1), ...
    all_Y_pred{mm}(all_labels{mm}==4,2),...
    all_Y_pred{mm}(all_labels{mm}==4,3),...
    'MarkerFaceColor','w','MarkerEdgeAlpha',1,'MarkerEdgeColor',[.1 .7 .1])

ax_lim(1) = min(min(min(all_Y_pred{mm})))*1.1;
ax_lim(2) = max(max(max(all_Y_pred{mm})))*1.1;
xlim([ax_lim(1) ax_lim(2)])
ylim([ax_lim(1) ax_lim(2)])
zlim([ax_lim(1) ax_lim(2)])

%tetrahedron
for i = 1:4
    for d = 1:3
        grand_means(i,d) = mean(all_Y_pred{mm}(all_labels{mm}==i,d));
    end
end
line_pairs = nchoosek([1:size(grand_means,1)],2);
for f = 1:size(line_pairs,1)
    line(...
        [grand_means(line_pairs(f,1),1) grand_means(line_pairs(f,2),1)],...
        [grand_means(line_pairs(f,1),2) grand_means(line_pairs(f,2),2)],...
        [grand_means(line_pairs(f,1),3) grand_means(line_pairs(f,2),3)],'Color','k')
end

%shadow
for f = 1:size(line_pairs,1)
    line(...
        [grand_means(line_pairs(f,1),1) grand_means(line_pairs(f,2),1)],...
        [grand_means(line_pairs(f,1),2) grand_means(line_pairs(f,2),2)],...
        [ax_lim(1) ax_lim(1)],'Color',[.5 .5 .5])
end

%lines to floor
for f = 1:size(grand_means,1)
    line(...
        [grand_means(f,1) grand_means(f,1)],...
        [grand_means(f,2) grand_means(f,2)],...
        [grand_means(f,3) ax_lim(1)],'Color',[.5 .5 .5])
end

xlabel('stim')
ylabel('choice')
zlabel('outcome')
view(126,23)

%% plot means

figure;
scatter3(all_pred_means(1,1,:), all_pred_means(1,2,:), all_pred_means(1,3,:),...
    'MarkerFaceColor',[.1 .7 .1],'MarkerFaceAlpha',1,'MarkerEdgeColor','none')
hold on;
scatter3(all_pred_means(2,1,:), all_pred_means(2,2,:), all_pred_means(2,3,:),...
    'MarkerFaceColor',[.75 0 0],'MarkerFaceAlpha',1,'MarkerEdgeColor','none')
scatter3(all_pred_means(3,1,:), all_pred_means(3,2,:), all_pred_means(3,3,:),...
    'MarkerFaceColor','w','MarkerEdgeAlpha',1,'MarkerEdgeColor',[.75 0 0])
scatter3(all_pred_means(4,1,:), all_pred_means(4,2,:), all_pred_means(4,3,:),...
    'MarkerFaceColor','w','MarkerEdgeAlpha',1,'MarkerEdgeColor',[.1 .7 .1])

ax_lim(1) = min(min(min(all_pred_means)))*1.1;
ax_lim(2) = max(max(max(all_pred_means)))*1.1;
xlim([ax_lim(1) ax_lim(2)])
ylim([ax_lim(1) ax_lim(2)])
zlim([ax_lim(1) ax_lim(2)])

%tetrahedron
grand_means = mean(all_pred_means,3);
line_pairs = nchoosek([1:size(all_pred_means,1)],2);
for f = 1:size(line_pairs,1)
    line(...
        [grand_means(line_pairs(f,1),1) grand_means(line_pairs(f,2),1)],...
        [grand_means(line_pairs(f,1),2) grand_means(line_pairs(f,2),2)],...
        [grand_means(line_pairs(f,1),3) grand_means(line_pairs(f,2),3)],'Color','k')
end

%shadow
for f = 1:size(line_pairs,1)
    line(...
        [grand_means(line_pairs(f,1),1) grand_means(line_pairs(f,2),1)],...
        [grand_means(line_pairs(f,1),2) grand_means(line_pairs(f,2),2)],...
        [ax_lim(1) ax_lim(1)],'Color',[.5 .5 .5])
end

%lines to floor
for f = 1:size(grand_means,1)
    line(...
        [grand_means(f,1) grand_means(f,1)],...
        [grand_means(f,2) grand_means(f,2)],...
        [grand_means(f,3) ax_lim(1)],'Color',[.5 .5 .5])
end

xlabel('stim')
ylabel('choice')
zlabel('outcome')
view(126,23)

%% plot accuracy

mean_total_accuracy = squeeze(mean(total_accuracy));
sem_total_accuracy = squeeze(std(total_accuracy))/sqrt(length(total_accuracy));
figure;
b = bar([1:3],mean_total_accuracy);
hold on
errorbar([1:3],mean_total_accuracy,sem_total_accuracy,'Color','k','LineStyle','none','capsize',0)
b.FaceColor = 'flat';
b.CData(1,:) = [.78 .62 1];
b.CData(2,:) = [.26 .40 1];
b.CData(3,:) = [.24 .15 .66];
line([0 4],[.5 .5],'Color',[.6 .6 .6],'LineStyle',':')
xlim([0 4])
xticklabels({'XOR' 'Choice' 'Outcome'})
xlabel('Decoder')
ylabel('Accuracy')
prettyPlot(gca)

%% plot accuracy against pseudo

colors = [.78 .62 1; .26 .40 1; .24 .15 .66];
figure;
hold on;
for s = 1:3
    subplot(3,1,s)
    h = histogram(squeeze(mean(shift_accuracy(:,:,s),1)),linspace(0,1,101),'FaceColor',[.6 .6 .6]);
    mx = max(h.Values);
    line([mean_total_accuracy(s) mean_total_accuracy(s)],[0 mx],'Color',colors(s,:),'LineWidth',2)
    xlim([0 1])
    prettyPlot(gca)
end

%% plot accuracy against pseudo - violin
figure;
hold on;
for s = 1:3
    [f,xi] = ksdensity(squeeze(mean(shift_accuracy(:,:,s),1)));
    fill(10*(f*(xi(end-1)-xi(end)))+s,xi,'k','LineStyle','none')
    alpha(.2);
    fill(10*(f*(xi(end)-xi(end-1)))+s,xi,'k','LineStyle','none')
    alpha(.2);
    errorbar(s,mean_total_accuracy(s),sem_total_accuracy(s),...
        'Marker','o','Color',colors(s,:),...
        'LineWidth',2,...
        'MarkerFaceColor',colors(s,:),...
        'MarkerEdgeColor','w',...
        'MarkerSize',10,...
        'capsize',0)
end
ylim([0.3 1])
xlim([.5 3.5])
prettyPlot(gca)
xticks(1:3)
xticklabels({'XOR' 'Choice' 'Outcome'})
xlabel('Decoder')
ylabel('Accuracy')

%% plot accuracy against pseudo - shift order
figure;
dl = {'XOR' 'Choice' 'Outcome'};
hold on;
for s = 1:3
    subplot(3,1,s)
    errorbar(0,mean(total_accuracy(:,s),1),...
        std(total_accuracy(:,s),1)/sqrt(size(total_accuracy,1)),...
        'Marker','o','Color',colors(s,:),...
        'LineWidth',1.5,...
        'LineStyle','none',...
        'MarkerFaceColor',colors(s,:),...
        'MarkerEdgeColor','w',...
        'MarkerSize',6,...
        'capsize',0);
    hold on;
    errorbar([-20:-1],mean(shift_accuracy(:,1:20,s),1),...
        std(shift_accuracy(:,1:20,s),1)/sqrt(size(shift_accuracy,1)),...
        'Marker','none','Color',[.6 .6 .6],...
        'LineWidth',1.5,...
        'LineStyle','none',...
        'MarkerFaceColor',[.6 .6 .6],...
        'MarkerEdgeColor','w',...
        'MarkerSize',10,...
        'capsize',0);
    errorbar([1:20],mean(shift_accuracy(:,21:40,s),1),...
        std(shift_accuracy(:,21:40,s),1)/sqrt(size(shift_accuracy,1)),...
        'Marker','none','Color',[.6 .6 .6],...
        'LineWidth',1.5,...
        'LineStyle','none',...
        'MarkerFaceColor',[.6 .6 .6],...
        'MarkerEdgeColor','w',...
        'MarkerSize',10,...
        'capsize',0)
    
    ylim([0.45 1])
    xlim([-21 21])
    xticks([-20 -10 0 10 20])
    prettyPlot(gca)
    ylabel('Accuracy')
    title(dl{s});
end
    xlabel('Shifts')

%% plot accuracy by label
mean_label_accuracy = squeeze(mean(label_accuracy));
sem_label_accuracy = squeeze(std(label_accuracy)/sqrt(length(label_accuracy)));

figure;
b = bar([1:3],mean_label_accuracy);
b(1).FaceColor = [.1 .7 .1]; b(1).EdgeColor = 'none';
b(2).FaceColor = [.75 0 0]; b(2).EdgeColor = 'none';
b(3).FaceColor = 'w'; b(3).EdgeColor = [.75 0 0]; b(3).LineWidth = 1.25;
b(4).FaceColor = 'w'; b(4).EdgeColor = [.1 .7 .1]; b(4).LineWidth = 1.25;

% Get the x coordinate of the bars
[ngroups,nbars] = size(mean_label_accuracy');
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end
hold on
for f = 1:size(mean_label_accuracy,1)
    errorbar(x(f,:),mean_label_accuracy(f,:),sem_label_accuracy(f,:),'Color','k','LineStyle','none','capsize',0)
end

line([0.5 4.5],[.5 .5],'Color',[.6 .6 .6],'LineStyle',':')
xticklabels({'XOR' 'Choice' 'Outcome'})
xlabel('Decoder')
ylabel('Accuracy')
legend({'Contra rewarded' 'contra unrewarded' 'Ipsi unrewarded' 'Ipsi rewarded'},'location','ne')
legend boxoff

prettyPlot(gca)
%% WORKBENCH
% ns = 1.5;
% x=c1.mean(1); y=c1.mean(2); z=-2; rx=ns*c1.std(1); ry=ns*c1.std(2);
% hold on
% th = 0:pi/50:2*pi;
% x_circle = rx * cos(th) + x;
% y_circle = ry * sin(th) + y;
% z_circle = 0 * sin(th) + z;
% circles = plot3(x_circle, y_circle,z_circle);
% 
% x=c2.mean(1); y=c2.mean(2); z=-2; rx=ns*c2.std(1); ry=ns*c2.std(2);
% hold on
% th = 0:pi/50:2*pi;
% x_circle = rx * cos(th) + x;
% y_circle = ry * sin(th) + y;
% z_circle = 0 * sin(th) + z;
% circles = plot3(x_circle, y_circle,z_circle);
% 
% x=c3.mean(1); y=c3.mean(2); z=-2; rx=ns*c3.std(1); ry=ns*c3.std(2);
% hold on
% th = 0:pi/50:2*pi;
% x_circle = rx * cos(th) + x;
% y_circle = ry * sin(th) + y;
% z_circle = 0 * sin(th) + z;
% circles = plot3(x_circle, y_circle,z_circle);
% 
% x=c4.mean(1); y=c4.mean(2); z=-2; rx=ns*c4.std(1); ry=ns*c4.std(2);
% hold on
% th = 0:pi/50:2*pi;
% x_circle = rx * cos(th) + x;
% y_circle = ry * sin(th) + y;
% z_circle = 0 * sin(th) + z;
% circles = plot3(x_circle, y_circle,z_circle);

  %% fit pseudo
%     
%     fit_type = 'linear';
% %     fit_type = 'logistic';
%     
%     if strcmp(fit_type,'linear')
%         options.alpha = 1;
%         options.nlambda = 20;
%         options.standardize = 'true';
%         nf = 3;
%         family = 'gaussian';
%         prog = 0;
%         fprintf(1,'fitting pseudos: %3d%%\n',prog);
%         for p = 1:np
%             for d = 1:size(Y_train,2)
%                 fit{d,p} = cvglmnet(X_train_shifted, Y_train_shifted(:,d,p), family, options, 'deviance',nf,[],true);
%                 Y_pred_pseudo(:,d,p) = cvglmnetPredict(fit{d,p}, X_test,'lambda_min');
%             end
%             prog = (100*(p/np));
%             fprintf(1,'\b\b\b\b%3.0f%%',prog);
%         end
%         fprintf('\n');
%         
%     elseif strcmp(fit_type,'logistic')
%         Y_train = (Y_train+1)/2;
%         options.alpha = 1;
%         options.nlambda = 20;
%         options.standardize = 'true';
%         nf = 3;
%         family = 'binomial';
%         prog = 0;
%         fprintf(1,'fitting pseudos: %3d%%\n',prog);
%         for p = 1:np
%             for d = 1:size(Y_train,2)
%                 fit{d,p} = cvglmnet(X_train_shifted, Y_train_shifted(:,d,p), family, options, 'deviance',nf,[],false);
%                 Y_pred(:,d) = cvglmnetPredict(fit{d,p}, X_test,'lambda_min');
%             end
%             prog = (100*(p/np));
%             fprintf(1,'\b\b\b\b%3.0f%%',prog);
%         end
%     end
