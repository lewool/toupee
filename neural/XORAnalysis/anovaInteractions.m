figure;
for m = 1:length(mouseList)
    mouseName = char(mouseList{m});
    expDate = char(expList{m}{1});
    expNum = expList{m}{2};  
    [expRef, ~] = data.constructExpRef(mouseName,expDate,expNum);
    
    fprintf('Loading %s...\n',expRef)
    [behavioralData, expInfo, neuralData] = data.loadDataset(mouseName, expDate, expNum);
    nt = length(behavioralData.eventTimes(1).daqTime);
    
    % extract stimulus, choice, feedback, value, and block values
    trueStimuli = expInfo.block.events.contrastValues(1:nt);
    trialCorrectChoice = expInfo.block.events.correctResponseValues(1:nt);
    trueBlocks = expInfo.block.events.highRewardSideValues(1:nt);
    trueChoices = expInfo.block.events.responseValues(1:nt);
    trueFeedback = double(expInfo.block.events.feedbackValues(1:nt));

    [baselineResps, stimResps, pmovResps, movResps, rewResps, preCueResps] = getEpochResps(neuralData.eta);
    trialTypes = getTrialTypes(expInfo,behavioralData,'late');

    for c = 1:size(rewResps,2)
        p{m}(c,:) = anovan(rewResps(trialTypes.singleVar.contrast{5, 1},c),{trueChoices(trialTypes.singleVar.contrast{5, 1})' trueFeedback(trialTypes.singleVar.contrast{5, 1})'},'model','interaction','display','off');
    end
    
%     propIntCells(m) = sum((p{m}(:,3) < .05))/size(rewResps,2);
    propIntCells(m) = sum((p{m}(:,3) < .05))/sum((p{m}(:,1) < .05));
    
    subplot(6,6,m)
    histogram(p{m}(:,3),linspace(0,1,41),'normalization','probability')
    prettyPlot(gca)
    clearvars -except mouseList expList hemList m propIntCells p
end

%%
figure;
for m = 1:length(mouseList)
    subplot(6,6,m)
    histogram(p{m}(p{m}(:,1) < .05,3),linspace(0,1,21),'normalization','probability') %choice
%     histogram(p{m}(p{m}(:,2) < .05,3),linspace(0,1,21),'normalization','probability') %feedback
%     histogram(p{m}(p{m}(:,1) < .05 & p{m}(:,2) < .05,3),linspace(0,1,21),'normalization','probability') % both
    prettyPlot(gca)

    mouseName = char(mouseList{m});
    expDate = char(expList{m}{1});
    expNum = expList{m}{2};  
    [expRef, ~] = data.constructExpRef(mouseName,expDate,expNum);
    subplot(6,6,m)
    title(expRef, 'interpreter','none')
end
set(gcf,'position',figpos)
%% venn
all_p = cat(1,p{:});
%main-effect choice neurons
me_c = find(all_p(:,1) < 0.05);
me_f = find(all_p(:,2) < 0.05);
int_cf = find(all_p(:,3) < 0.05);

setListData = {me_c, me_f, int_cf};
h = vennEulerDiagram(setListData, 'drawProportional', true, 'SetLabels', ...
    [strcat("Choice (",num2str(length(me_c)),")"); ...
     strcat("Outcome (",num2str(length(me_f)),")"); ...
     strcat("Interaction (",num2str(length(int_cf)),")")]);
 
 n_me_c_me_f = length(intersect(me_f,me_c));
 n_me_c_int = length(intersect(int_cf,me_c));
 n_me_f_int = length(intersect(int_cf,me_f));
 n_allint = length(intersect(intersect(me_f,me_c),int_cf));
 n_no_sig = find(all_p(:,1) >= 0.05 & all_p(:,2) >= 0.05 & all_p(:,3) >= 0.05)
 

