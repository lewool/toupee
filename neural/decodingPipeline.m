for m = 1:length(mouseList)
    
    mouseName = char(mouseList{m});
    expDate = char(expList{m}{1});
    expNum = expList{m}{2};
    [expRef, ~] = data.constructExpRef(mouseName,expDate,expNum);
    fprintf('Loading %s...',expRef)
    [behavioralData, expInfo, neuralData] = data.loadDataset(mouseName, expDate, expNum);
    [neuralData] = alignResps(expInfo, neuralData, behavioralData, 5);
    expInfo.hemisphere = hemList(m);
    fprintf('done\n')
    sessDir = fullfile('G:\Workspaces',mouseName,expDate,num2str(expNum));
    cd(sessDir);
    if ~exist('decoderAnalysis')
        fprintf('Loading decodingAnalysis.mat...');
        try
            load(fullfile(sessDir,'decodingAnalysis.mat'));
            fprintf('done\n');
        catch
            fprintf('no decoder file (yet); will generate one after this run\n')
        end
    end
    
    eta = 'gocue';    
    features = {'choice' 'side' 'feedback' 'block' 'value' ...
                'choiceGivenStim' 'sideGivenChoice' ...
                'choice_contrast' 'side_contrast' ...
                'choice_0' 'side_0'};
            
    for f = 1:length(features)
        predictY = features{f};
        fprintf('Predicting %s...',predictY);        
        if exist('decoderAnalysis') && ...
            isfield(decoderAnalysis, predictY) && ...
            isfield(decoderAnalysis.(matlab.lang.makeValidName(predictY)), eta)
            fprintf('Analysis already done; skipping!\n')
        else
            if strcmp(eta,'stim')
                ETA = 1;
            elseif strcmp(eta,'gocue')
                ETA = 4;
            elseif strcmp(eta,'move')
                ETA = 2;
            elseif strcmp(eta,'feedback')
                ETA = 3;
            end

            if strcmp(predictY,'side_contrast') || strcmp(predictY,'choice_contrast')  || strcmp(predictY,'feedback_contrast')
                [fits, mutual_info, X, Y, D] = neuralDecoder_contrasts(expInfo, behavioralData, neuralData, predictY, ETA);
            elseif strcmp(predictY,'side_0') || strcmp(predictY,'choice_0')
                [fits, mutual_info, X, Y, D] = neuralDecoder_contrasts_easyModel(expInfo, behavioralData, neuralData, predictY, ETA);
            else
                [fits, mutual_info, X, Y, D] = neuralDecoder(expInfo, behavioralData, neuralData, predictY, ETA);
            end

            decoderAnalysis.(matlab.lang.makeValidName(predictY)).(matlab.lang.makeValidName(eta)).fits = fits;
            decoderAnalysis.(matlab.lang.makeValidName(predictY)).(matlab.lang.makeValidName(eta)).mutual_info = mutual_info;
            fprintf('done\n');
        end
    end
    fprintf('Saving...')
    save('decodingAnalysis.mat','decoderAnalysis', '-v7.3')
    fprintf('done\n')
    
    clearvars -except mouseList expList hemList m
    
end

%% plot

plotY = 'choiceIndStim';

for m=1:length(mouseList)
prcts = prctile(predictions(m).(matlab.lang.makeValidName(plotY)).gocue.pseudo,0:.05:100);

for t = 1:size(prcts,2)
    [~, ff] = min(abs((prcts(:,t) - predictions(m).(matlab.lang.makeValidName(plotY)).gocue.true(t))));
    ptile(m,t) = ((ff/length(prcts) > 0.975) | (ff/length(prcts) < 0.025));
end
end

% example session

eventWindow = -2:.1:2;
timeRange = 6:2:36;
ew = eventWindow(timeRange);
figure;
set(gcf,'position',[214 1286 1214 282]);
yl1 = [.49 .72];
yl2 = [.49 .72];
% yl1 = [5.9 11];
% yl2 = [5.9 11];
ex = 22;
if strcmp(plotY,'stimulus')
    pt = 1./predictions(ex).(matlab.lang.makeValidName(plotY)).gocue.true;
    pp = 1./predictions(ex).(matlab.lang.makeValidName(plotY)).gocue.pseudo;
else
    pt = predictions(ex).(matlab.lang.makeValidName(plotY)).gocue.true;
    pp = predictions(ex).(matlab.lang.makeValidName(plotY)).gocue.pseudo;
end
subplot(1,3,1)
hold on;
pp_UB = prctile(pp,97.5);
pp_LB = prctile(pp,2.5);
plotSignal(ew,pt,pp_UB,pp_LB,[0 0 0],'-');
line([0 0],[yl1],'Color',[.5 .5 .5],'LineStyle',':')
line([-.8 -.8],[yl1],'Color',[.5 .5 .5],'LineStyle',':')
xlim([-1.5 1.5]);
prettyPlot(gca)
xlabel('Time from go cue')
if strcmp(plotY,'stimulus')
    ylabel('1/MSE')
else
ylabel('Norm. log-likelihood')
end
title('Decoding (1 example session)')
ylim(yl1)
% plot number of significant timebins
legend boxoff
legend('true','95% CI pseudo','Location','nw');
propp = sum(ptile)/size(ptile,1)*100;

subplot(1,3,2)
hold on;

bar(ew,propp,'k','LineWidth',2)
line([0 0],[min(propp)*.9 max(propp)*1.05],'Color',[.5 .5 .5],'LineStyle',':')
line([-.8 -.8],[min(propp)*.9 max(propp)*1.05],'Color',[.5 .5 .5],'LineStyle',':')

prettyPlot(gca)
xlabel('Time from go cue')
ylabel('Percent')
title('Sessions with significance')
ylim([min(propp)*.9 max(propp)*1.05])


% mean prediction across all sessions
if strcmp(plotY,'stimulus')
    for m = 1:length(mouseList)
        pb(m,:) = 1./predictions(m).(matlab.lang.makeValidName(plotY)).gocue.true;
    end
else
    for m = 1:length(mouseList)
    pb(m,:) = predictions(m).(matlab.lang.makeValidName(plotY)).gocue.true;
    end
end


subplot(1,3,3)
% plot(ew,pb,'Color',[.5 .5 .5]);
hold on;
% plot(ew, mean(pb,1),'LineWidth',2)
    plotSignal(ew,mean(pb,1),mean(pb,1)+std(pb,1)/sqrt(36),mean(pb,1)-std(pb,1)/sqrt(36),[0 .3 .9],'-')
line([0 0],[yl2],'Color',[.5 .5 .5],'LineStyle',':')
line([-.8 -.8],[yl2],'Color',[.5 .5 .5],'LineStyle',':')
prettyPlot(gca)
xlabel('Time from go cue')
if strcmp(plotY,'stimulus')
    ylabel('1/MSE')
else
ylabel('Norm. log-likelihood')
end
% title('Mean \pm SEM decoding (all sessions)')
% for t = 1:size(prcts,2)
%     p(t) = myBinomTest(sum(ptile(:,t)),length(ptile(:,t)),0.05);
%     if p(t) < 0.05 && p(t) >= 0.01
%         text(ew(t), max(yl2)*.9,'*')
%     elseif p(t) < 0.01
%         text(ew(t), max(yl2)*.9,'**')
%     end
% end 
xlim([-1.5 1.5])
ylim([yl2])
set(gcf,'position',figpos)