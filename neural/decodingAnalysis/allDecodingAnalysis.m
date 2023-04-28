%% collect data across files

for m = 1:length(mouseList)
    
    mouseName = char(mouseList{m});
    expDate = char(expList{m}{1});
    expNum = expList{m}{2};  
    [expRef, ~] = data.constructExpRef(mouseName,expDate,expNum);
    fprintf('Loading %s...\n',expRef)
    sessDir = fullfile('G:\Workspaces',mouseName,expDate,num2str(expNum));
    cd(sessDir);
    load(fullfile(sessDir,'decodingAnalysis.mat'));       
    
    eta = 'gocue';    
    features = {'choice' 'side' 'feedback' 'block' 'value' ...
                'choiceGivenStim' 'sideGivenChoice' ...
                'choice_contrast' 'side_contrast'};
    
    conList = fieldnames(decoderAnalysis.side_contrast.gocue.mutual_info.true);
    for f = 1:length(features)
        if contains(features{f},'_contrast') || contains(features{f},'_0')
%             conList = fieldnames(decoderAnalysis.(matlab.lang.makeValidName(features{f})).gocue.mutual_info.true);
            try
                for c = 1:length(conList)
                    con_ttmp(c,:) =  decoderAnalysis.(matlab.lang.makeValidName(features{f})).gocue.mutual_info.true.(matlab.lang.makeValidName(conList{c}));
                    con_ptmp(c,:,:) =  decoderAnalysis.(matlab.lang.makeValidName(features{f})).gocue.mutual_info.pseudo.(matlab.lang.makeValidName(conList{c}));
                end
            catch
                con_ttmp = nan(5,23);
                con_ptmp = nan(5,40,23);
            end
            allDecoding.(matlab.lang.makeValidName(features{f})).true(m,:,:) = con_ttmp;
            allDecoding.(matlab.lang.makeValidName(features{f})).pseudo(m,:,:,:) = con_ptmp;
        else
            true_tmp = decoderAnalysis.(matlab.lang.makeValidName(features{f})).gocue.mutual_info.true;
            pseudo_tmp = decoderAnalysis.(matlab.lang.makeValidName(features{f})).gocue.mutual_info.pseudo;
            allDecoding.(matlab.lang.makeValidName(features{f})).true(m,:) = true_tmp;
            allDecoding.(matlab.lang.makeValidName(features{f})).pseudo(m,:,:) = pseudo_tmp;  
        end
    end
end
            
%% test for significance

for f = 1:length(features)
    if contains(features{f},'_contrast')
        for m = 1:length(mouseList)
            for c = 1:length(conList)
                tt = squeeze(allDecoding.(matlab.lang.makeValidName(features{f})).true(m,c,:))';
                pt = squeeze(allDecoding.(matlab.lang.makeValidName(features{f})).pseudo(m,c,:,:));
                allDecoding.(matlab.lang.makeValidName(features{f})).sig(m,c,:) = tt > prctile(pt,95);
            end
        end
    else
        for m = 1:length(mouseList)
            tt = allDecoding.(matlab.lang.makeValidName(features{f})).true(m,:);
            pt = squeeze(allDecoding.(matlab.lang.makeValidName(features{f})).pseudo(m,:,:));
            allDecoding.(matlab.lang.makeValidName(features{f})).sig(m,:) = tt > prctile(pt,95);
        end
    end
end

%% separate sessions
figure;
for m = 1:length(mouseList)
    subplot(6,6,m)
    hold on;
    plot(ew, squeeze(allDecoding.choiceGivenStim.pseudo(m,:,:)),'Color',[.6 .6 .6])
    plot(ew, allDecoding.choiceGivenStim.true(m,:),'Color','k','LineWidth',2)
    xlim([-1.5 2.9])
    ylim([-.05 1])
    prettyPlot(gca)
end
%% total significance - binomial test


for f = 1:length(features)
    if contains(features{f},'_contrast')
        for c = 1:length(conList)
            s = squeeze(sum(allDecoding.(matlab.lang.makeValidName(features{f})).sig(:,c,:)))';
            p(c,:) = myBinomTest(s,length(mouseList),.05,'one');
        end
    else
        s = squeeze(sum(allDecoding.(matlab.lang.makeValidName(features{f})).sig));
        p = myBinomTest(s,length(mouseList),.05,'one');    
    end
    p1 = double(p<0.05);
    p1(p1==0)=nan;
end

%%
heights = [.85 .9 .95 1 1.05];
figure;
hold on;
for f = 1:length(features)
    subplot(3,3,f)
    title(features{f},'interpreter','none')
    ylabel('Mutual information')
    xlabel('Time from go cue (s)')
    line([0 0],[-.4 1.5],'Color',[.5 .5 .5],'LineStyle','--')
    line([-1 -1],[-.4 1.5],'Color',[.5 .5 .5],'LineStyle','--')
    hold on
    if contains(features{f},'_contrast')
        for c = 1:length(conList)
            psm = squeeze(nanmean(allDecoding.(matlab.lang.makeValidName(features{f})).true(:,c,:),1))';
            pss = squeeze(nanstd(allDecoding.(matlab.lang.makeValidName(features{f})).true(:,c,:),[],1)/sqrt(length(allDecoding.choice_contrast.true(:,c,:))))';
            plot(ew, reshape(squeeze(nanmean(allDecoding.(matlab.lang.makeValidName(features{f})).pseudo,1)),[length(conList)*40,length(ew)]),'Color',[.6 .6 .6])
            errorbar(ew,psm,pss,'capsize',0,'Marker','o','MarkerSize',8,'MarkerFaceColor',colors(c,:),'MarkerEdgeColor','w','Color',colors(c,:),'LineStyle','-','LineWidth',1)
        end
        
        for c = 1:length(conList)
            s = squeeze(sum(allDecoding.(matlab.lang.makeValidName(features{f})).sig(:,c,:)))';
            p(c,:) = myBinomTest(s,length(mouseList),.05,'one');
        end
        p1 = double(p<0.05);
        p1(p1==0)=nan;
        for c = 1:length(conList)
            plot(ew,p1(c,:).*heights(c),'*','Color',colors(c,:))
        end
    else
        psm = squeeze(nanmean(allDecoding.(matlab.lang.makeValidName(features{f})).true,1));
        pss = squeeze(nanstd(allDecoding.(matlab.lang.makeValidName(features{f})).true,[],1)/sqrt(length(allDecoding.(matlab.lang.makeValidName(features{f})).true(:,c,:))));
        plot(ew, squeeze(mean(allDecoding.(matlab.lang.makeValidName(features{f})).pseudo,1)),'Color',[.6 .6 .6])
        errorbar(ew,psm,pss,'capsize',0,'Marker','o','MarkerSize',8,'MarkerFaceColor','k','MarkerEdgeColor','w','Color','k','LineStyle','-','LineWidth',1)
        
        s = squeeze(sum(allDecoding.(matlab.lang.makeValidName(features{f})).sig));
        p = myBinomTest(s,length(mouseList),.05,'one');    
        p1 = double(p<0.05);
        p1(p1==0)=nan;
        plot(ew,p1,'k*')
    end
    xlim([-1.5 2.9])
    ylim([-.15 1.1])
    prettyPlot(gca)
end

%% test for significance & plot
colors = [0 0 0; .6 .6 .6; .4 .4 .4; .2 .2 .2; 0 0 0];
ew = -1.5:.2:3;
features = {'feedback' 'block' 'value' 'choice' 'side' ...
            '' '' '' 'choiceGivenStim' 'sideGivenChoice' ...
            '' '' '' 'choice_contrast' 'side_contrast' ...
            };

titles =   {'Feedback' 'Block' 'Value' 'Choice' 'Visual stimulus side' ...
            '' '' '' 'Choice, given stimulus' 'Stimulus side, given choice'  ...
            '' '' '' 'Choice, given stimulus (0%)' 'Stimulus side, given choice (0%)' ...
            };
        
heights = [.85 .9 .95 1 1.05];
figure;
hold on;
for f = 1:length(features)
    subplot(3,5,f)
    if ~isempty(features{f})
        title(titles{f},'interpreter','none')
        ylabel('Mutual information')
        xlabel('Time from go cue (s)')
        line([0 0],[-.4 1.5],'Color',[.5 .5 .5],'LineStyle','--')
        line([-1 -1],[-.4 1.5],'Color',[.5 .5 .5],'LineStyle','--')
        hold on
        if contains(features{f},'_contrast')
            for c = 1%:length(conList)
                psm = squeeze(nanmean(allDecoding.(matlab.lang.makeValidName(features{f})).true(:,c,:),1))';
                pss = squeeze(nanstd(allDecoding.(matlab.lang.makeValidName(features{f})).true(:,c,:),[],1)/sqrt(length(allDecoding.choice_contrast.true(:,c,:))))';
                pseudoFill = squeeze(nanmean(allDecoding.(matlab.lang.makeValidName(features{f})).pseudo(:,1,:,:),1));
                plotSignal(ew,nan(1,23),prctile(pseudoFill,2.5),prctile(pseudoFill,97.5),[.2 .2 .2],'-');
%                 plot(ew, reshape(squeeze(nanmean(allDecoding.(matlab.lang.makeValidName(features{f})).pseudo,1)),[length(conList)*40,length(ew)]),'Color',[.6 .6 .6])
                errorbar(ew,psm,pss,'capsize',0,'Marker','o','MarkerSize',7,'MarkerFaceColor',colors(c,:),'MarkerEdgeColor','w','Color',colors(c,:),'LineStyle','-','LineWidth',1)
            end

            for c = 1%:length(conList)
                pt = squeeze(nanmean(allDecoding.(matlab.lang.makeValidName(features{f})).true(:,c,:),1))';
                pp = squeeze(nanmean(allDecoding.(matlab.lang.makeValidName(features{f})).pseudo(:,c,:,:),1));
                p(c,:) = double(pt > max(pp(2:39,:)));
            end
            p(p==0)=nan;
            for c = 1%:length(conList)
                plot(ew,p(c,:).*heights(c),'*','Color',colors(c,:),'MarkerSize',5)
            end
        else
            psm = squeeze(nanmean(allDecoding.(matlab.lang.makeValidName(features{f})).true,1));
            pss = squeeze(nanstd(allDecoding.(matlab.lang.makeValidName(features{f})).true,[],1)/sqrt(length(allDecoding.(matlab.lang.makeValidName(features{f})).true(:,c,:))));
            pseudoFill = squeeze(nanmean(allDecoding.(matlab.lang.makeValidName(features{f})).pseudo,1));
            plotSignal(ew,nan(1,23),prctile(pseudoFill,2.5),prctile(pseudoFill,97.5),[.2 .2 .2],'-');
%             plot(ew, squeeze(mean(allDecoding.(matlab.lang.makeValidName(features{f})).pseudo,1)),'Color',[.6 .6 .6])
            errorbar(ew,psm,pss,'capsize',0,'Marker','o','MarkerSize',7,'MarkerFaceColor','k','MarkerEdgeColor','w','Color','k','LineStyle','-','LineWidth',1)

            pt = squeeze(nanmean(allDecoding.(matlab.lang.makeValidName(features{f})).true,1));
            pp = squeeze(nanmean(allDecoding.(matlab.lang.makeValidName(features{f})).pseudo,1));
            p = double(pt > max(pp(2:39,:)));
            p(p==0)=nan;
            plot(ew,p,'k*','MarkerSize',5)
        end
        xlim([-1.8 3])
        ylim([-.15 1.1])
        prettyPlot(gca)
    end
end



























%%
revList = fliplr(conList);
figure;
hold on;
for f = 1:length(features)
    subplot(3,3,f)
    title(features{f},'interpreter','none')
    ylabel('Proportion of sessions')
    xlabel('Time from go cue (s)')
    line([0 0],[-.4 1.5],'Color',[.5 .5 .5],'LineStyle','--')
    line([-1 -1],[-.4 1.5],'Color',[.5 .5 .5],'LineStyle','--')
    line([-1.6 3],[5/36 5/36],'Color',[.5 .5 .5],'LineStyle',':')
    line([-1.6 3],[0.05 0.05],'Color',[.5 .5 .5],'LineStyle',':')
    hold on
    if contains(features{f},'_contrast')
        for c = 1:length(conList)
            line([-1.6 3],[0 0],'Color','k','LineStyle','-')
            plot(ew,...
                squeeze(sum(allDecoding.(matlab.lang.makeValidName(features{f})).sig(:,c,:)))/length(allDecoding.(matlab.lang.makeValidName(features{f})).true),...
                'Color',colors(c,:),'Marker','_','LineStyle','none','LineWIdth',2)
%             bar(ew,...
%                 squeeze(sum(allDecoding.(matlab.lang.makeValidName(features{f})).sig(:,c,:)))/length(allDecoding.(matlab.lang.makeValidName(features{f})).true),...
%                 'FaceColor',colors(c,:),'FaceAlpha',.5)
            
        end
        
    else
        bar(ew,...
            sum(allDecoding.(matlab.lang.makeValidName(features{f})).sig)/length(allDecoding.(matlab.lang.makeValidName(features{f})).true),...
            'FaceColor','k')
    end
    xlim([-1.6 3])
    ylim([-.1 1.05])
    prettyPlot(gca)
end