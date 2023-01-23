% expInfo = initExpInfo({{'LEW032'}},{{'2020-02-13',1,[1]}});
% expInfo = data.loadExpData(expInfo);

%% 1. generate 'true' block

close all; clearvars -except expInfo; clc;

blockType = 'real';
nt = 1000;
blockStart = 'rand';

if strcmp(blockStart,'fixed')
    firstSide = -1;
elseif strcmp(blockStart,'rand')
    firstSide = randsample([-1, 1],1,true);
end

if strcmp(blockType, 'real')
    nt = numel(expInfo.block.events.endTrialValues);
    trueBlock = expInfo.block.events.highRewardSideValues(1:nt);
    lbTrials = find(trueBlock<0);
    rbTrials =find(trueBlock>0);
    bC_trials = lbTrials;
    bI_trials = rbTrials;
    
elseif strcmp(blockType, 'task')
    switches = cumsum(125+randi(100,1,20));
    for s = 1:length(switches)
        if s == 1
            trueBlock(1:switches(s)-1) = firstSide;
        elseif mod(s,2) == 1
            trueBlock(switches(s-1):switches(s)-1) = firstSide;
        elseif mod(s,2) == 0
            trueBlock(switches(s-1):switches(s)-1) = -firstSide;
        end
    end
    trueBlock = trueBlock(1:nt);

    lbTrials = find(trueBlock<0);
    rbTrials =find(trueBlock>0);
    bC_trials = lbTrials;
    bI_trials = rbTrials;
    
elseif strcmp(blockType,'poisson')
            
    switches = find(binornd(1,0.004,[nt+2000 1]));
    for s = 1:length(switches)
        if s == 1
            b(1:switches(s)-1) = firstSide;
        elseif mod(s,2) == 1
            b(switches(s-1):switches(s)-1) = firstSide;
        elseif mod(s,2) == 0
            b(switches(s-1):switches(s)-1) = -firstSide;
        end
    end
    b=b(1:nt);
    trueBlock = b;
    lbTrials = find(trueBlock<0);
    rbTrials =find(trueBlock>0);  
    bC_trials = lbTrials;
    bI_trials = rbTrials;
    
elseif strcmp(blockType,'rand')
    firstSide = -1;
    trueBlock = randsample([-1, 1],nt-1,true);
    trueBlock = [firstSide trueBlock];
    lbTrials = find(trueBlock<0);
    rbTrials =find(trueBlock>0);
    bC_trials = lbTrials;
    bI_trials = rbTrials;
end

fig = figure;
set(gcf,'position',[20 120 560 1200])
subplot(50,1,1)
hold on;
plot(trueBlock,'r')
axis off
set(gca, 'XTickLabels', {})
set(gca, 'YTickLabels', {})
title('"True" block vs. example pseudoblocks')


%% 2. generate fake drift neurons
nc = 2000;
fakeNeurons = poissrnd(5,nt,nc);
for f = 1:size(fakeNeurons,2)
    fakeDriftNeurons(:,f) = fakeNeurons(:,f) + linspace(randi(10),randi(10),nt)';
end
whichResps = fakeDriftNeurons;

fig2 = figure;
set(gcf,'position',[600 120 560 1200])
for c = 1:20
    subplot(10,2,c)
    plot(whichResps(:,c))
    box off
    set(gca,'tickdir','out')
    ylim([0 25])
    xlim([1 nt])
end

%% 3. compute mean diff between blocks for each neuron

bC_resps = nanmean(whichResps(bC_trials,:),1)';
bI_resps = nanmean(whichResps(bI_trials,:),1)';
blockResp = bC_resps - bI_resps;

%% 4. generate pseudosessions
np = 10000;
pseudoSessions = nan(np,nt);
bC_trials_pseudo = cell(1,np);
bI_trials_pseudo = cell(1,np);

if strcmp(blockType,'task') || strcmp(blockType,'real')
    for p = 1:np
        if strcmp(blockStart,'fixed')
            firstSide = -1;
        elseif strcmp(blockStart,'rand')
            firstSide = randsample([-1, 1],1,true);
        end
        b=nan(1,nt);
        switches = cumsum(125+randi(100,1,20));
        for s = 1:length(switches)
            if s == 1
                b(1:switches(s)-1) = firstSide;
            elseif mod(s,2) == 1
                b(switches(s-1):switches(s)-1) = firstSide;
            elseif mod(s,2) == 0
                b(switches(s-1):switches(s)-1) = -firstSide;
            end
        end
        pseudoSessions(p,:) = b(1:nt);
        bC_trials_pseudo{p} = find(pseudoSessions(p,:) < 0);
        bI_trials_pseudo{p} = find(pseudoSessions(p,:) > 0);
    end
    
elseif strcmp(blockType,'poisson')
    for p = 1:np
        if strcmp(blockStart,'fixed')
            firstSide = -1;
        elseif strcmp(blockStart,'rand')
            firstSide = randsample([-1, 1],1,true);
        end
        b=nan(1,nt);
        switches = find(binornd(1,0.004,[nt+2000 1]));
        for s = 1:length(switches)
            if s == 1
                b(1:switches(s)-1) = firstSide;
            elseif mod(s,2) == 1
                b(switches(s-1):switches(s)-1) = firstSide;
            elseif mod(s,2) == 0
                b(switches(s-1):switches(s)-1) = -firstSide;
            end
        end
        pseudoSessions(p,:) = b(1:nt);
        bC_trials_pseudo{p} = find(pseudoSessions(p,:) < 0);
        bI_trials_pseudo{p} = find(pseudoSessions(p,:) > 0);
    end
elseif strcmp(blockType,'rand')
    for p = 1:np
        if strcmp(blockStart,'fixed')
            firstSide = -1;
        elseif strcmp(blockStart,'rand')
            firstSide = randsample([-1, 1],1,true);
        end
        b=nan(1,nt);
        b = randsample([-1, 1],nt-1,true);
        b = [firstSide b];
        pseudoSessions(p,:) = b(1:nt);
        bC_trials_pseudo{p} = find(pseudoSessions(p,:) < 0);
        bI_trials_pseudo{p} = find(pseudoSessions(p,:) > 0);
    end
end

figure(fig);
for p = 2:50
    subplot(50,1,p)
    plot(pseudoSessions(p,:))
    axis off
    set(gca, 'XTickLabels', {})
    set(gca, 'YTickLabels', {})
end

%% 5. compute mean diff between each block for each neuron

bC_resps_pseudo = nan(nc,np);
bI_resps_pseudo = nan(nc,np);
blockResp_pseudo = nan(nc,np);

for p = 1:np
    bC_resps_pseudo(:,p) = nanmean(whichResps(bC_trials_pseudo{p},:),1)';
    bI_resps_pseudo(:,p) = nanmean(whichResps(bI_trials_pseudo{p},:),1)';
    blockResp_pseudo(:,p) = bC_resps_pseudo(:,p) - bI_resps_pseudo(:,p);
end

%% compute significance distribution

blockSig = nan(1,nc);
sigrank = nan(1,nc);
pValue = nan(1,nc);

for c = 1:nc
    UB = prctile(blockResp_pseudo(c,:),97.5);
    LB = prctile(blockResp_pseudo(c,:),2.5);
    if blockResp(c) < LB ||  blockResp(c) > UB
        blockSig(c) = 1;
    else 
        blockSig(c) = 0;
    end
    brp = [];
    brp = blockResp_pseudo(c,~isnan(blockResp_pseudo(c,:)));
    [~, idx] = min(abs(sort(brp) - blockResp(c)));
    if idx > 500
        pValue(c) = 2*((length(brp)-idx)/length(brp));
    else
        pValue(c) = 2*idx/length(brp);
    end
    sigrank(c) = idx/length(brp);
end

propBS = sum(blockSig)/length(blockResp);

figure;
set(gcf,'position',[1780 120 560 420])
h = histogram(sigrank*100,linspace(0,100,41),'FaceColor',[.5 .5 .5]);
hold on
maxy = max(h.Values);
LB = fill([0; 0; 2.5; 2.5],[0; h.Values(1); h.Values(1); 0],'r','LineStyle','none','FaceAlpha',.5);
LB = fill([97.5; 97.5; 100; 100],[0; h.Values(end); h.Values(end); 0],'r','LineStyle','none','FaceAlpha',.5);

ylim([0 maxy*1.2])
box off
set(gca,'tickdir','out')
xlabel('Percentile')
ylabel('No. neurons')
text(2,maxy*1.1,strcat(num2str(round(propBS*100)),'% p < 0.05'),'Color',[.5 0 0])

figure;
set(gcf,'position',[1180 120 560 1200])
for s = 1:20
subplot(10,2,s)
histogram(blockResp_pseudo(s,:),linspace(-5,5,41));
hold on
if blockSig(s)
    line([blockResp(s) blockResp(s)],[0 180],'LineStyle','-','Color','r');
else
    line([blockResp(s) blockResp(s)],[0 180],'LineStyle','-','Color','k');
end
text(-5,200,num2str(round(sigrank(s)*100)))
box off
axis off
end


% %%
% pseudoSessions = nan(np,nt);
% bC_trials_pseudo = cell(1,np);
% bI_trials_pseudo = cell(1,np);
% 
% for p = 1:np
%     firstSide = -1;
%     switches = find(binornd(1,0.004,[nt+500 1]));
%     for s = 1:length(switches)
%         if s == 1
%             b(1:switches(s)-1) = firstSide;
%         elseif mod(s,2) == 1
%             b(switches(s-1):switches(s)-1) = firstSide;
%         elseif mod(s,2) == 0
%             b(switches(s-1):switches(s)-1) = -firstSide;
%         end
%     end
%     pseudoSessions(p,:) = b(1:nt);
%     bC_trials_pseudo{p} = find(pseudoSessions(p,:) < 0);
%     bI_trials_pseudo{p} = find(pseudoSessions(p,:) > 0);
% end
% 
%  
% for p = 1:np
%     bC_resps_pseudo(:,p) = nanmean(whichResps(bC_trials_pseudo{p},:),1)';
%     bI_resps_pseudo(:,p) = nanmean(whichResps(bI_trials_pseudo{p},:),1)';
%     blockResp_pseudo(:,p) = bC_resps_pseudo(:,p) - bI_resps_pseudo(:,p);
% end
% 
% 
% 
% 
% 
% 
% 
% 
% %%
% tic
% np = 1000;
% pseudoSessions{m} = nan(np,length(plotCells),nt);
% 
% for c = 1:length(plotCells)
% for p = 1:length(np)
%     randSwitches = cumsum(125+randi(100,1,10));
%     totalSwitches = 0;
%     firstSide = randsample([-1, 1],1,true);
%         for t = 1:nt
%             if any(t == randSwitches)
%                 totalSwitches = totalSwitches + 1;
%             end
%             stayOrSwitch = -((mod(totalSwitches,2)*2)-1);
%             pseudoSessions{m}(p,c,t) = stayOrSwitch * firstSide;
%         end
% end  
% end
%     
%     
% pseudoSessions = nan(np,nc,nt);
% 
% for c = 1:nc
%     for p = 1:np
% %         firstSide = randsample([-1, 1],1,true);
%         firstSide = -1;
%         switches = find(binornd(1,0.004,[1500 1]));
%         for s = 1:length(switches)
%             if s == 1
%                 b(1:switches(s)-1) = firstSide;
%             elseif mod(s,2) == 1
%                 b(switches(s-1):switches(s)-1) = firstSide;
%             elseif mod(s,2) == 0
%                 b(switches(s-1):switches(s)-1) = -firstSide;
%             end
%         end
%         pseudoSessions(p,c,:) = b(1:nt);
%     end
% end
% 
% bC_trials_pseudo = cell(nc,np);
% bI_trials_pseudo = cell(nc,np);
% for c = 1:nc
%     for p = 1:1000
%         bC_trials_pseudo{c,p} = find(pseudoSessions(p,c,:) < 0);
%         bI_trials_pseudo{c,p} = find(pseudoSessions(p,c,:) > 0);
%     end
% end   
% 
% for c = 1:nc
%     for p = 1:np
%         bC_resps_pseudo(c,p) = nanmean(whichResps(bC_trials_pseudo{p},c),1)';
%         bI_resps_pseudo(c,p) = nanmean(whichResps(bI_trials_pseudo{p},c),1)';
%         blockResp_pseudo(c,p) = bC_resps_pseudo(c,p) - bI_resps_pseudo(c,p);
%     end
% end